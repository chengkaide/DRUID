"""
druid.reduction.trace —— 测点容器 Tra 与装载流程
===============================================

Tra（TRansient Acquisition，瞬时采集）是一个测点的"全部家当"：
原始信号、净信号、时间轴、剥蚀区间的各种时刻坐标、样品名与角色。

为什么要做成一个 dataclass 而不是到处传 tuple
---------------------------------------------
原先的版本返回 (t, raw, t0, t1, ...) 这种长元组，调用方必须记位置顺序；
一旦某天要在中间插一个字段，所有调用点全炸。
用 dataclass 后字段名就是契约，加字段只要给个默认值即可。

一个测点对象里为什么要存 raw 和 net 两份
----------------------------------------
· raw 是"事实"，任何时候都能回溯检查（比如复查空白取得合不合理）；
· net 是"计算用版本"，绝大多数算法都基于它。
两者都不大（几千个 float），没必要为了省内存牺牲可回溯性。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

import numpy as np

from ..core.constants import MASSES_NEEDED
from ..core.deadtime import correct_channels
from ..io.qtegra import read_qtegra
from .ablation import find_ablation, subtract_blank


@dataclass
class Tra:
    """一个测点（一次激光剥蚀）的完整数据。"""

    # ── 身份信息 ──
    idx: int                 # 在批次中的序号（从 0 开始）
    sample: str              # 样品名，如 "91500" / "YL-46-1"
    role: str                # 角色：primary_std / secondary_std / unknown / glass
    path: str                # 源 CSV 路径（排查问题时要用）

    # ── 时间轴与信号 ──
    t: np.ndarray                            # 时间轴 (s)
    raw: Dict[int, np.ndarray] = field(default_factory=dict)   # 原始 cps
    net: Dict[int, np.ndarray] = field(default_factory=dict)   # 扣空白后 cps
    blank: Dict[int, float] = field(default_factory=dict)      # 各通道扣掉的空白值

    # ── 时刻坐标（单位 s） ──
    t_ab0: float = 0.0       # 剥蚀起点（自动识别得到）
    t_ab1: float = 0.0       # 剥蚀终点
    t0: float = 0.0          # **积分窗口**起点（已裁掉 trim 秒）
    t1: float = 0.0          # 积分窗口终点

    # ── 元信息 ──
    au_deadtime_ns: float = 0.0   # 本次装载实际使用的死时间 (ns)，0 表示未校正
    note: str = ""                # 备注（比如"因窗口过短未裁边"）

    @property
    def name(self) -> str:
        """不含路径和扩展名的测点名，用于图件标题和文件名。"""
        return Path(self.path).stem

    @property
    def label(self) -> str:
        """人类可读标签：序号 + 样品名。"""
        return f"{self.idx + 1:02d} {self.sample}"

    def duration(self) -> float:
        """有效积分窗口长度（秒）。"""
        return float(self.t1 - self.t0)


def load_spot(idx: int,
              path,
              sample: str,
              role: str,
              trim: float = 1.5,
              deadtime_ns: float = 0.0,
              blank_dur: float = 15.0) -> Tra:
    """
    从磁盘上的一个 Qtegra CSV 装载出一个 Tra。

    处理顺序（**不能调换**）
    -----------------------
        读文件 → 死时间校正 → 识别剥蚀区间 → 扣气体空白 → 裁剪瞬态边缘

    ① 为什么死时间校正在扣空白之前：
       空白段也是计数，同样被探测器"漏计"，必须先还原真计数率再谈本底。

    ② 为什么最后才裁边：
       裁剪要基于已经识别好的剥蚀区间，而识别本身要用完整信号。

    参数
    ----
    idx         : 批次内序号（0 起）
    path        : CSV 路径
    sample      : 样品名
    role        : 角色（见 io.sequence.sample_role）
    trim        : 两端各裁掉多少秒。激光开/关的瞬间会有一个信号尖峰与
                  快速爬升段，此时气溶胶粒径分布不稳定，分馏行为与稳定
                  剥蚀阶段不同，必须去掉。
    deadtime_ns : 探测器死时间（纳秒）。0 或负数表示不做校正（默认）。
    blank_dur   : 气体空白取样时长（秒）

    返回
    ----
    Tra 对象

    异常
    ----
    KeyError   : 缺少 U-Pb 定年必需的通道
    ValueError : 文件本身读不出来（由 read_qtegra 抛出）
    """
    path = Path(path)

    # ── ① 读原始文件 ──
    _title, t, raw, _cols = read_qtegra(path)

    # ── ② 通道完整性检查 ──
    # 少了任何一个关键质量数，后面的比值就无从算起，早点报错比中途 NaN 好
    missing = [m for m in MASSES_NEEDED if m not in raw]
    if missing:
        raise KeyError(f"缺少通道 {missing}")

    # ── ③ 死时间校正（可选）──
    if deadtime_ns and deadtime_ns > 0:
        raw = correct_channels(raw, deadtime_ns)

    # ── ④ 定位剥蚀区间 ──
    # 用 238U 通道：信号最强，"激光是否开火"的形态最清晰
    ta0, ta1 = find_ablation(t, raw[238])

    # ── ⑤ 扣气体空白 ──
    net, blank = subtract_blank(t, raw, ta0, blank_dur)

    # ── ⑥ 裁剪瞬态边缘 ──
    t0, t1 = ta0 + trim, ta1 - trim
    note = ""
    if t1 - t0 < 3.0:
        # 极短的剥蚀（可能是信号太弱被误判），此时再裁边就什么都不剩了，
        # 退回使用完整的剥蚀区间，并在备注里留下痕迹便于复核
        t0, t1 = ta0, ta1
        note = "剥蚀窗口过短，未裁边"

    return Tra(idx=idx, sample=sample, role=role, path=str(path),
               t=t, raw=raw, net=net, blank=blank,
               t_ab0=ta0, t_ab1=ta1, t0=t0, t1=t1,
               au_deadtime_ns=deadtime_ns, note=note)


def window_mask(t, a, b):
    """
    构造 [a, b] 时间窗口的布尔掩码 —— 全流程里最不起眼却最常用的一个函数。

    单独抽出来是为了：所有地方（整段积分、滑动窗口、短窗口 Jackknife）用**同一个**
    闭开区间约定，避免某处写 >=、另一处写 > 这种隐蔽的不一致。
    """
    return (t >= a) & (t <= b)
