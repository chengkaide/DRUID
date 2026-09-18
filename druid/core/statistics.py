"""
druid.core.statistics —— 同位素年代学里反复用到的几个统计量
========================================================

为什么要单独一个文件
--------------------
"加权平均 + MSWD" 在函数资料的加权、年龄域的合并、标样 QC 里各写了一遍
（原先 _legacy 版本就是这样重复了三处），一旦其中一处公式写错，
结果会微妙地不一致。集中定义后三处共享同一份实现。

MSWD（Mean Square of Weighted Deviates，加权偏差均方）
------------------------------------------------------
    MSWD = Σ wᵢ(xᵢ − x̄)² / (n − 1)，其中权重 wᵢ = 1/σᵢ²

怎么读这个数：
    MSWD ≈ 1        → 数据散布程度与各自的不确定度相符，一致年龄没问题
    MSWD ≫ 1        → 散布远大于不确定度：要么有地质意义上的多期年龄，
                      要么 σ 被低估了（常见于只算了计数统计、没算外部重现性）
    MSWD ≪ 1        → 散布远小于不确定度：σ 被高估了
"""
from __future__ import annotations

import numpy as np


def weighted_mean(x, s):
    """
    反比方差加权平均，同时给出 MSWD。

    参数
    ----
    x : 观测值数组
    s : 每个观测值的 **1σ** 标准误数组（不是 2σ！传 2σ 会让 MSWD 缩小 4 倍）

    返回
    ----
    (均值, 均值的1σ标准误, MSWD, 有效点数)
        均值的1σ = 1 / sqrt(Σwᵢ)

    边界处理
    --------
    · 自动剔除 NaN/Inf 和非正的不确定度（不合法点不应让整个平均失败）；
    · n = 0 时返回全 NaN 而不是抛异常，调用方可以用 np.isfinite 判断；
    · n = 1 时 MSWD 的分母为 0，返回 NaN（MSWD 在单点上没有定义）。
    """
    x = np.asarray(x, float)
    s = np.asarray(s, float)
    # 只保留"值有限 + 不确定度有限且 > 0"的点
    keep = np.isfinite(x) & np.isfinite(s) & (s > 0)
    x, s = x[keep], s[keep]
    n = x.size
    if n == 0:
        return np.nan, np.nan, np.nan, 0

    w = 1.0 / s ** 2                          # 权重 = 方差的倒数
    mu = float((w * x).sum() / w.sum())       # 加权均值
    se = float(1.0 / np.sqrt(w.sum()))        # 加权均值的标准误
    mswd = float((w * (x - mu) ** 2).sum() / (n - 1)) if n > 1 else np.nan
    return mu, se, mswd, n


def robust_mask(v, k=3.0):
    """
    用 **中位数 + MAD** 做稳健离群点判别，返回保留掩码。

    为什么不用"均值 ± 3σ"
    ----------------------
    均值和标准差本身就会被离群点带偏（掩盖效应 masking）。
    中位数和中位绝对偏差（MAD）有 50% 的崩溃点：
    即使一半的数据是离群点，它们依然稳健。

    MAD → σ 的换算系数 1.4826
    -------------------------
    对正态分布，MAD = 0.6745σ，所以 σ ≈ MAD / 0.6745 = 1.4826·MAD。

    返回
    ----
    (保留掩码 bool 数组, 中位数)
    若 MAD = 0（数据严重离散或几乎全等），无从判别，一律保留。
    """
    v = np.asarray(v, float)
    med = float(np.median(v))
    # MAD = median(|x − median|)，再乘 1.4826 换算成正态等效 σ
    mad = 1.4826 * float(np.median(np.abs(v - med)))
    if mad <= 0:
        return np.ones(v.size, bool), med
    return np.abs(v - med) <= k * mad, med


def external_scatter(values, sigmas, lo=0.003, hi=0.05):
    """
    由监控标样的实测散度估计 **外部重现性**（相对，1σ）。

    背景
    ----
    单点不确定度若只由 jackknife（内部回放误差）给出，会系统性偏小：
    它反映的是"本窗口内各 cycle 之间的随机涨落"，
    抓不到"逐点之间的仪器漂移/基体差异/剥蚀条件差异"。

    于是把观测到的总散度按方差分解：
        σ²_观测 = σ²_内部均值 + σ²_外部
        ⇒ σ_外部 = sqrt( max(σ²_观测 − mean(σ²_内部), 0) )

    返回
    ----
    外部重现性 σ_ext（相对值，如 0.007 表示 0.7%）

    为什么要 clip 到 [lo, hi]
    --------------------------
    · 下限：防止 QC 点太少、偶然高度一致时算出 0，导致后续真实不确定度崩塌；
    · 上限：防止某次 QC 严重失败（比如打到包裹体）时给出荒谬的外部误差，
            把全部样品的不确定度撑大到没有分辨力。
    """
    v = np.asarray(values, float)
    s = np.asarray(sigmas, float)
    ok = np.isfinite(v) & np.isfinite(s)
    v, s = v[ok], s[ok]
    if v.size < 2 or np.mean(v) == 0:
        return float(lo)
    # 观测到的相对散度（用样本标准差）
    obs = (np.std(v, ddof=1) / np.mean(v)) ** 2
    # 内部不确定度的平均贡献（同样用相对量）
    internal = float(np.mean((s / v) ** 2))
    # 方差相减后开方，并做上下限保护
    return float(np.clip(np.sqrt(max(obs - internal, 0.0)), lo, hi))
