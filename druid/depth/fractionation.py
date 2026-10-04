"""
druid.depth.fractionation —— 深度依赖的分馏校正因子 F(τ)
======================================================

核心公式只有一行：

    F(τ) = R_ref / R_measured(τ)

通俗说：**标样在同一深度处"偏了多少"，就认为样品在同一深度处也偏了多少**。

为什么必须"同深度"
------------------
Down-hole 分馏是坑深的函数，不是常数：
    τ = 0（坑口）  → 气溶胶粒径细、传输效率高、Pb 相对富集
    τ = 1（坑底）  → 坑深增大、气溶胶变粗、传输效率下降，Pb/U 逐渐降低
本批次 206Pb/238U 从坑口到坑底漂移约 **20%**，这是一个巨大的系统误差，
如果只用"整段的平均值"（即 F_bulk），年龄剖面会明显倾斜，
甚至被误判成"两期年龄"。

夹逼（bracketing）的含义
------------------------
    标样 — 样品 — 样品 — 样品 — 标样
     ↑________________________↑
              一对夹逼标样

取时间上一个在前、一个在后的两个（或更多）标样，
在每个 τ 上取它们的**平均实测比值**，用来构造样品的 F(τ)。
这样仪器随时间的线性漂移（ Detektor 灵敏度缓慢下降等）也被一阶抵消。
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from ..core.constants import L238
from ..core.geochronology import age68, age76
from ..core.statistics import weighted_mean


# ─────────────────────────────────────────────────────────────────────────────
# 一、构造 F(τ)
# ─────────────────────────────────────────────────────────────────────────────
def bracket_F(tau: np.ndarray,
              profiles: Sequence[pd.DataFrame],
              key: str,
              ref: float) -> np.ndarray:
    """
    由若干夹逼标样的窗口剖面构造出目标深度网格上的 F(τ)。

    参数
    ----
    tau      : 目标深度网格（0~1），通常是样品窗口的 τ
    profiles : 标样的 window_profile DataFrame 列表（一般 1~2 个）
    key      : 用哪个比值构造，"R68" 或 "R76"
    ref      : 标样的参考比值

    返回
    ----
    F(τ) 数组，形状与 tau 相同

    细节
    ----
    · np.interp 做线性插值：标样的 τ 网格与样品的 τ 网格一般不会重合，
      需要先插值到同一个网格上才能逐点相除。
    · 多个标样时在**每个 τ 上取算术平均**：
      等价于假设两个标样对真值的偏离是对称的（线性漂移假设）。
    · meas ≤ 0 的位置返回 nan 而不是 inf：
      负比值没有物理意义，用 nan 标记为无效，让上层决定是否丢弃。

    ⚠⚠ 两个静默失效点（2026-10-04 修）
    ----------------------------------
    ① **NaN 传染**：原先用 `np.mean(vals, axis=0)`，于是**一个**标样在某个 τ
       上是 nan，该 τ 的均值就是 nan —— 一个坏标样窗口废掉一整段深度。
       改成 `np.nanmean`，并要求**至少一个**标样在该点有效。
    ② **外插被夹成常数**：`np.interp` 在查询点超出 `xp` 范围时**钳到端点值**、
       **不外推**。而 τ=0 是坑口、τ=1 是坑底，恰好是漂移最大的两端
       （模块文档说 206Pb/238U 从坑口到坑底漂移约 20%）；样品的 τ 网格由
       **它自己的**剥蚀时长归一化，标样的由**标样自己的**归一化，
       两者时长不同是常态 ⇒ **标样覆盖不到的深度段会拿到端点值**，
       误差方向单一、随深度系统性变化，正好造出「剖面两端弯曲」的假结构。
       现在把那两段显式置为 nan（**不外插**），让上层决定，而不是给一个
       看起来正常、其实是夹出来的数。
       实测影响：示例批次 **31/1677** 个窗口（1.8%）的 Age68 变 NaN。
    """
    # 把每个标样的实测比值插值到目标 τ 网格上
    vals = [np.interp(tau, df["tau"].to_numpy(), df[key].to_numpy())
            for df in profiles]
    arr = np.asarray(vals, float)                  # shape = (n_standard, n_tau)

    # ① NaN 传染：原来 np.mean 会让「一个标样坏 ⇒ 整段变 nan」。
    with np.errstate(invalid="ignore"):
        n_valid = np.sum(np.isfinite(arr), axis=0)
        meas = np.nanmean(np.where(np.isfinite(arr), arr, np.nan), axis=0)
    # 一个有效标样都没有的点保持 nan
    meas = np.where(n_valid > 0, meas, np.nan)

    # ② 外插钳位：标样 τ 覆盖不到的两端不外推，显式标为无效。
    #    覆盖范围取**各标样的交集**（保守：任一标样没覆盖就不算覆盖）。
    if profiles:
        lo = max(float(df["tau"].to_numpy().min()) for df in profiles)
        hi = min(float(df["tau"].to_numpy().max()) for df in profiles)
        outside = (np.asarray(tau, float) < lo) | (np.asarray(tau, float) > hi)
        meas = np.where(outside, np.nan, meas)

    with np.errstate(divide="ignore", invalid="ignore"):
        # 先造一个"安全的分母"：meas ≤ 0 的位置用 1.0 占位，避免除零，
        # 再用 np.where 把这些位置的结果替换成 nan
        safe = np.where(meas > 0, meas, 1.0)
        return np.where(meas > 0, ref / safe, np.nan)


# ─────────────────────────────────────────────────────────────────────────────
# 二、把窗口比值转成窗口年龄
# ─────────────────────────────────────────────────────────────────────────────
def profile_ages(prof: pd.DataFrame, F: np.ndarray, sigma_ext: float):
    """
    窗口比值 + F(τ) → 窗口年龄与 1σ。

    参数
    ----
    prof      : 某测点的 window_profile 表
    F         : 与 prof 行数相同的分馏因子数组
    sigma_ext : 外部重现性（相对值，如 0.007）

    返回
    ----
    (age, sigma_age)，单位 Ma

    不确定度怎么合成
    ----------------
    两个来源，正交相加（平方和再开方）：
        ① 内部误差 σ68 · F          —— 本窗口内的计数/闪烁噪声
        ② 外部误差 σ_ext · R68_corr  —— 逐点之间的重现性差异
    注意两者都是**绝对量**（比值的标准差），可以平方相加。

    比值 σ → 年龄 σ 的换算
    -----------------------
        age = ln(1+R)/λ ⇒ d(age)/dR = 1 / (λ·(1+R))
    所以 σ_age = σ_R / (λ·(1+R))，最后 /1e6 把年换成 Ma。
    """
    R = prof["R68"].to_numpy() * F                       # 校正后的 206Pb/238U
    s = np.hypot(prof["s68"].to_numpy() * F,             # 内部误差，随 F 一起缩放
                 sigma_ext * R)                          # 外部误差（相对量×比值）
    return age68(R), s / (L238 * (1.0 + R)) / 1e6


def profile_ages_76(prof: pd.DataFrame, F76: np.ndarray, sigma_ext76: float):
    """
    窗口的 207Pb/206Pb 年龄（较老锆石才用得上）。

    ⚠ 用法上的重要提醒
    -------------------
    207Pb/206Pb 是 **Pb 同位素内部的比值**，两个 Pb 同位素的分馏行为几乎一样，
    所以它理论上**不随坑深漂移**，也不应该套 F(τ)。
    计数上也压不住这种校正：4 s 窗口的 207Pb 净计数中位约 9×10³，
    R76 的计数统计已在 1% 量级，而能指望的 Pb-Pb 漂移远小于它 ——
    套 F(τ) 只会把标样的噪声**注入**到样品里，让剖面抖得没法看。

    所以调用此函数时，通常传的 F76 是全 1 的数组（即不校正），
    或者由上一层显式决定。这里保留参数位是为了保持接口对称。

    ⚠ 关于调用方（2026-10-02 核对）
    ------------------------------
    本函数**当前没有任何调用方** —— `workflow.py` 只 import 了 `bracket_F` 与
    `profile_ages`，`depth/__init__.py` 把它当作公开 API 导出，但全仓搜不到使用者
    （对应的原因见 `进度与交接.md`：`剖面窗口` 表刻意不补 `Age76/Age76_1s` 两列）。
    所以"改了它会不会动到结果"的答案是**动不到** —— 这也意味着**端到端基线
    （`check_example_batch.py`）对它无效**，它唯一的守护是
    `tests/test_depth.py::test_profile_ages_76_matches_the_scalar_definition`。
    要动这个函数，先看那条测试，别指望基线会替你发现。

    实现
    ----
    原先这里是逐点循环调三次**标量** `age76`，而 `age76` 本身已经向量化
    （内部是 80 次数组二分，逐点调用等于把同一段循环跑 137 遍）。
    2026-10-02 改为整体向量化，并在 137 个随机点上与旧写法逐位比对：
    `age` 与 `sigma` 的**最大差都是 0.0**（逐位相等，不是"接近"）。
    """
    R = prof["R76"].to_numpy() * F76
    s = np.hypot(prof["s76"].to_numpy() * F76, sigma_ext76 * R)
    # σ 用"上下各偏 s 再取半差"而不是解析式：207/206 的年龄方程是隐式的
    # （age76 内部是二分法），没有解析导数，差分是最稳的近似。
    return age76(R), np.abs(age76(R + s) - age76(R - s)) / 2.0


# ─────────────────────────────────────────────────────────────────────────────
# 三、诊断工具：F(τ) vs F_bulk 的对照实验
# ─────────────────────────────────────────────────────────────────────────────
def compare_fractionation(prof: pd.DataFrame,
                          profiles: Sequence[pd.DataFrame],
                          ref68: float,
                          sigma_ext: float) -> dict:
    """
    在同一个点上对比两种分馏校正模型（最好挑一个**已知均一**的锆石来做）。

        F(τ)    逐窗口取标样同一深度的分馏因子   ← 正确做法
        F_bulk  整段用一个常数分馏因子           ← 常见错误

    判据
    ----
    均一样品本该给出水平的年龄剖面。若模型错了，
    "末窗口年龄 − 首窗口年龄"（drift）会明显不为 0，MSWD 也会偏大。
    """
    tau = prof["tau"].to_numpy()
    F_tau = bracket_F(tau, profiles, "R68", ref68)
    F_bulk = np.full_like(F_tau, float(np.nanmean(F_tau)))   # 常数版本

    out = {}
    for tag, F in (("F(tau)", F_tau), ("F_bulk", F_bulk)):
        a, sa = profile_ages(prof, F, sigma_ext)
        mu, se, mswd, k = weighted_mean(a, sa)
        out[tag] = dict(ages=a, sig=sa, mean=mu, se=se, mswd=mswd, n=k,
                        first=float(a[0]), last=float(a[-1]),
                        drift=float(a[-1] - a[0]))
    return out
