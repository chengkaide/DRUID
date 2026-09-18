"""
druid.core.deadtime —— 电子倍增器死时间校正（**可选**，默认不启用）
================================================================

为什么会有这件事
----------------
在计数模式下，每一个被记录到的脉冲之后，探测器会有一段"失明"时间 τ
（dead time，典型 10~50 ns），此间到达的离子 **不会被记录**。
于是：

    计数率越高 → 漏计的越多 → 实测值系统性偏低

非瘫痪型（non-paralyzable）模型下的修正公式：

    I_true = I_obs / (1 − I_obs · τ)

其中 I_obs 单位是"计数/秒"，τ 单位是"秒"，
两者的乘积就是"每秒漏计的时间比例"。

量级感（τ = 14.76 ns，本批次由承保拟合得到）：
    I_obs = 1e5 cps  → 修正 +0.15%
    I_obs = 1e6 cps  → 修正 +1.50%   ← 样品锆石就是这个量级
    I_obs = 3e6 cps  → 修正 +4.63%

为什么关键是 **比值**
---------------------
如果所有通道同比例被压低，比值完全不受影响——但实际上各通道计数率差了几个数量级
（238U 可达 1e6 cps，而 207Pb 只有几千 cps），压低的比例完全不同，
**206Pb/238U 这种比值就会被扭曲**，而且随 down-hole 分馏导致的信号强度变化而变化。

本项目为什么默认不走这条路
---------------------------
1. τ 的具体数值依赖探测器高压老化状态，厂商标称值常与实际有出入；
2. 同一批次里用主标（91500，~1e5 cps）去归一化样品（~1e6 cps）时，
   死时间效应本身就是"基体/强度失配"的一部分；
3. 保守做法是用**计数率相匹配的二级标样做 QC 二次校准**（见 workflow 模块），
   这不需要事先知道 τ。
因此这里保留完整的實現，通过命令行开关 `--deadtime-ns` 决定是否启用，
便于做对照实验。
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def correct_rate(cps: np.ndarray, tau_ns: float) -> np.ndarray:
    """
    把观测计数率 (cps) 校正为真计数率。

    参数
    ----
    cps     : 观测计数率数组，单位 counts per second
    tau_ns  : 探测器死时间，单位 **纳秒**

    返回
    ----
    校正后的计数率数组（形状同输入）

    安全性
    ------
    分母 (1 − I·τ) 在 I 极大（≥ 1/τ）时会穿过 0 变负，公式失效。
    这里用 np.maximum 把分母下限钳在 0.5，即最多允许 2 倍修正：
    超过这个范围说明数据或 τ 有问题，钳位比给出天文数字更有用。
    """
    if tau_ns is None or tau_ns <= 0:
        # τ ≤ 0 表示不需要校正，原样返回（不要返回引用，避免调用方意外改到原数组）
        return np.asarray(cps, float).copy()

    tau_s = float(tau_ns) * 1e-9        # 纳秒 → 秒，与 cps 相乘才是无量纲比例
    obs = np.asarray(cps, float)
    denom = 1.0 - obs * tau_s           # 未被"失明"占据的时间比例
    denom = np.maximum(denom, 0.5)      # 防御：极端情况下公式失去物理意义
    return obs / denom


def correct_channels(raw: Dict[int, np.ndarray], tau_ns: float) \
        -> Dict[int, np.ndarray]:
    """
    对一个测点的**全部通道**批量施加死时间校正。

    参数
    ----
    raw    : {质量数: cps 数组}
    tau_ns : 死时间（纳秒）；≤ 0 时返回原字典的浅拷贝

    注意
    ----
    死时间校正必须在 **扣气体空白之前** 做：
    空白也是计数，同样受死时间压缩，先后顺序不能颠倒。
    """
    if tau_ns is None or tau_ns <= 0:
        return dict(raw)
    return {m: correct_rate(v, tau_ns) for m, v in raw.items()}
