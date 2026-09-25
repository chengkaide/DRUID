"""
druid.depth.windows —— 滑动窗口切分与逐窗口还原
=============================================

把一个剥蚀段切成一串有重叠的窗口，逐窗口还原出比值与 τ（归一化深度）。

为什么要重叠
------------
窗口宽 win=4 s、步长 step=1 s，意味着相邻窗口有 3 s 是共用数据的。
重叠换来的是**更高的时间分辨率**：虽然每个窗口的时间精度仍是 4 s，
但窗口间隔只有 1 s，能看清 1 s 尺度内的变化。
代价是相邻窗口的点不再独立（不能直接用点数当自由度），
这一点在 domains 模块的 MSWD 判据里要留意。

参数怎么选
----------
win 太小 → 窗口里只剩几个 cycle（4 s 才 13 个），各通道计数统计都变差；
           先撑不住的是计数本来就少的 204Pb（4 s 窗口净计数中位约 37）
win 太大 → 把核-边过渡带糊在一起，看不出结构
step 太小 → 相邻点几乎重复，图会很"滑"，但看不出更多信息
本批实测的合理取值：win = 4 s，step = 1 s。
cycle 周期实测 0.31871 s，所以 4 s 窗口只有 **12~13 个 cycle**
（`剖面窗口` 表的 `n_cycles` 列就是这个数，示例批次中位 13）。
重算：`python tools/window_scale_facts.py`
"""
from __future__ import annotations

from typing import List, Sequence

import numpy as np
import pandas as pd

from ..reduction.ratios import reduce_interval
from ..reduction.trace import window_mask


# ─────────────────────────────────────────────────────────────────────────────
# 一、窗口边界
# ─────────────────────────────────────────────────────────────────────────────
def window_edges(t0: float, t1: float, win: float, step: float) -> np.ndarray:
    """
    生成一串窗口的 [起, 止] 时刻，形状 (n_win, 2)。

    参数
    ----
    t0, t1 : 剥蚀积分窗口的起止时刻 (s)
    win    : 窗口宽度 (s)
    step   : 步长 (s)

    边界情况
    --------
    若整段比一个窗口还短，就退化为"一个覆盖全段的窗口"，
    而不是返回空集 —— 因为调用方通常希望在图上至少看到一个点。
    """
    if t1 - t0 <= win:
        return np.array([[t0, t1]])
    # np.arange 的上界加 1e-9：浮点误差可能让最后一个本应落在范围内的
    # 起点刚好差一点点被排除
    starts = np.arange(t0, t1 - win + 1e-9, step)
    return np.column_stack([starts, starts + win])


# ─────────────────────────────────────────────────────────────────────────────
# 二、逐窗口还原
# ─────────────────────────────────────────────────────────────────────────────
def window_profile(tr, win: float, step: float, fix=None,
                   n_sigma: float = 2.0) -> pd.DataFrame:
    """
    对一个测点做滑动窗口还原，返回逐窗口的比值表。

    参数
    ----
    tr       : Tra 对象
    win/step : 窗口宽度与步长 (s)
    fix      : (f206, sk) —— 强制使用整段的普通铅比例，见 ratios.reduce_interval。
               深度剖面**必须**传这个值，否则每个短窗口独立判 204 会全是噪声。
    n_sigma  : 204Pb 显著性门槛（独立模式下才生效）

    返回
    ----
    DataFrame，每行一个窗口，包含 ratios.reduce_interval 的全部字段，
    外加：
        t0, t1   : 窗口起止时刻
        t_mid    : 窗口中点时刻（画图用横轴）
        tau      : **归一化深度**，0 = 坑口，1 = 坑底
                   τ = (t_mid − t0_剥蚀) / 剥蚀时长
    """
    rows: List[dict] = []
    span = max(tr.t1 - tr.t0, 1e-9)      # 剥蚀总时长，防除零

    for a, b in window_edges(tr.t0, tr.t1, win, step):
        # 逐窗口调用同一套还原逻辑 —— 与整段还原共用代码，保证一致
        r = reduce_interval(tr.net, window_mask(tr.t, a, b),
                            n_sigma=n_sigma, fix=fix)
        if r is None:
            # 窗口太短 / 点数不足 → 跳过（行数为 0 时调用方会检查 empty）
            continue
        # 补上窗口坐标信息
        r.update(t0=a, t1=b,
                 t_mid=0.5 * (a + b),
                 tau=(0.5 * (a + b) - tr.t0) / span)
        rows.append(r)
    return pd.DataFrame(rows)


def window_sums(tr, win: float, step: float,
                masses: Sequence[int] = (206, 207, 208, 232, 238)):
    """
    统计每个窗口内各通道的**净计数积分和**（不是均值！）。

    用途
    ----
    计算"整段比值 + F(τ) 逐窗口校正"的加权合成：
        R68 = Σ_w (206Pb_w · F68(τ_w)) / Σ_w 238U_w
    这里用的是**计数求和**而非"比值的平均"，因为前者才等价于
    "把所有 cycle 合成一次积分"，符合计数统计的自然规律。

    返回
    ----
    (tau, dict_of_arrays)
        tau           : 每个窗口的归一化深度
        dict_of_arrays: {质量数: 该质量数在各窗口内的积分和}
    """
    S = {m: [] for m in masses}
    tau = []
    for a, b in window_edges(tr.t0, tr.t1, win, step):
        m = window_mask(tr.t, a, b)
        for k in masses:
            S[k].append(float(tr.net[k][m].sum()))
        tau.append((0.5 * (a + b) - tr.t0) / max(tr.t1 - tr.t0, 1e-9))
    return np.asarray(tau), {k: np.asarray(v) for k, v in S.items()}
