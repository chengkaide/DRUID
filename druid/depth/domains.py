"""
druid.depth.domains —— 年龄域分割（核 / 边 / 过渡带识别）
======================================================

拿到一条"年龄随坑深变化"的曲线后，要回答：**这条曲线几段？边界在哪？**

三个环节
--------
    ① segment        BIC 引导的二叉递归分割：切成几段，段界在哪
    ② refine_domains 边界精修：把紧邻边界的"混合窗口"剥掉
    ③ summarize      汇总成表，并把中间的短域标成"过渡带"

为什么还需要 ②
--------------
核→边不是一个数学上的阶跃，而是一个**过渡带**（本批数据约 1 s 宽，
对应窗口宽度的一部分）落在过渡带上的窗口里混合了两个年龄域的物质，
它报出的数值是两域的加权平均，**不代表任何一期地质事件**。

如果把这样的窗口算进某一个年龄域：
    · 该域的加权平均年龄被拉向另一个域；
    · MSWD 显著大于 1（因为混合点与均值不相容）。
所以判据很自然：**一个真正的年龄域必须与"常数年龄模型"统计相容**。
不满足就从"离均值更远"的那一侧剥掉一个窗口，反复迭代。
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd

from ..core.statistics import chi2_sf, weighted_mean


# ─────────────────────────────────────────────────────────────────────────────
# 〇、"常数年龄模型"的相容判据（全模块唯一定义）
# ─────────────────────────────────────────────────────────────────────────────
def mswd_acceptance(n_win: int, alpha: float = 2.0) -> float:
    """
    MSWD 的相容上限：低于它才认为"这一段只有一个年龄"。

        MSWD ≤ 1 + α·√(2/(n−1))        （默认 α = 2，约 2σ）

    为什么必须只有一份定义
    ----------------------
    这个式子原先在 `refine_domains`（默认参数 α=2.0）与 `summarize_segments`
    （硬编码 `1.0 + 2.0*sqrt(...)`）里各写了一遍。两处一旦不一致，就会出现
    "精修时认为这个域没问题、汇总时又把它标成过渡带"这种自相矛盾的结果，
    而且**不报错**——只是某一个域莫名其妙被丢掉。现在两处都调这里。

    ⚠ 这个判据用的是**窗口数** n，而窗口是重叠的（win/step > 1，
    实测 ρ = 0.75），真实独立观测数比 n 小约一倍多。所以上限偏**宽松**：
    方向是漏杀混合窗口，不会错杀真年龄域。详见 AGENTS.md「已知地雷」。
    """
    return 1.0 + alpha * np.sqrt(2.0 / max(int(n_win) - 1, 1))


# ─────────────────────────────────────────────────────────────────────────────
# 一、BIC 引导的二叉分割
# ─────────────────────────────────────────────────────────────────────────────
def _sse(a: np.ndarray, w: np.ndarray) -> float:
    """加权残差平方和（RSS）：用 w 加权后到加权均值的偏离总量。"""
    m = (w * a).sum() / w.sum()
    return float((w * (a - m) ** 2).sum())


def _best_split(a: np.ndarray, w: np.ndarray):
    """
    穷举所有可能的单一分割点，返回使"两段各自 RSS 之和"最小的那一个。

    为什么用穷举而不是先找最大落差
    ------------------------------
    单点噪声会让"最大落差"落在随机噪声最大的地方；
    而 RSS 是整段累积的统计量，对单个噪声点不敏感，稳健得多。
    窗口数通常只有 20~40 个，穷举 O(n²) 完全可以接受。
    """
    n = a.size
    best = None
    for k in range(1, n):
        wl, wr = w[:k].sum(), w[k:].sum()
        if wl <= 0 or wr <= 0:          # 权重为 0 的段无法定义均值，跳过
            continue
        ml = (w[:k] * a[:k]).sum() / wl          # 左段加权均值
        mr = (w[k:] * a[k:]).sum() / wr          # 右段加权均值
        s = float((w[:k] * (a[:k] - ml) ** 2).sum()
                  + (w[k:] * (a[k:] - mr) ** 2).sum())
        if best is None or s < best[0]:
            best = (s, k)
    return best


def segment(age, sig, n_min: int = 3, max_segs: int = 3) -> List[Tuple[int, int]]:
    """
    BIC 引导的二叉递归分割。

    参数
    ----
    age     : 窗口年龄数组
    sig     : 对应的 1σ 数组
    n_min   : 每个域最少要包含多少个窗口（少于就是不值得分的东西）
    max_segs: 最多分几段（锆石通常核+边两段，偶尔有个幔，3 够用）

    为什么用 BIC 而不是"看 P 值"
    ---------------------------
    BIC = n·ln(RSS/n) + k·ln(n)
         ↑ 拟合优度            ↑ 参数个数惩罚
    它自动惩罚模型复杂度：多一段虽然必然让 RSS 变小，
    但要多付现金 3·ln(n) 的复杂度代价。**避免了无限细分**。
    每多一个分割点，代价是 2 个额外参数（两段各一个均值）。

    返回
    ----
    [(lo, hi), ...] 元素为 (起始窗口下标, 结束下标+1)，按下标排序。
    若样本太少则返回 []（调用方视为"整个剖面是一个域"）。
    """
    a = np.asarray(age, float)
    s = np.asarray(sig, float)
    ok = np.isfinite(a) & np.isfinite(s) & (s > 0)
    a, s = a[ok], s[ok]
    if a.size < 2 * n_min:
        return []
    w = 1.0 / s ** 2                       # 权重
    out: List[Tuple[int, int]] = []
    _recurse(a, w, 0, a.size, max_segs - 1, n_min, out)
    out.sort(key=lambda x: x[0])           # 按窗口顺序排列
    return out


def _recurse(a, w, lo, hi, depth, n_min, out):
    """递归体的实现，不直接对外暴露。"""
    n = hi - lo
    # 停止条件 1：样本太少；停止条件 2：已达最大段数
    if n < 2 * n_min or depth <= 0:
        out.append((lo, hi))
        return

    seg_a, seg_w = a[lo:hi], w[lo:hi]
    s1 = _sse(seg_a, seg_w)                       # 单段模型的 RSS
    bs = _best_split(seg_a, seg_w)
    if bs is None:
        out.append((lo, hi))
        return
    s2, k = bs                                    # 两段的合计 RSS + 分割点偏移
    # 停止条件 3：任何一侧段长不足 n_min
    if k < n_min or n - k < n_min:
        out.append((lo, hi))
        return

    # 两个模型的 BIC。+1 / +3 是各自的参数个数
    bic1 = n * np.log(max(s1, 1e-12) / n) + 1 * np.log(n)
    bic2 = n * np.log(max(s2, 1e-12) / n) + 3 * np.log(n)
    # BIC 越小越好；两段的 BIC 没变小 → 保留单段
    if bic2 >= bic1:
        out.append((lo, hi))
        return

    # 两边各自继续递归
    _recurse(a, w, lo, lo + k, depth - 1, n_min, out)
    _recurse(a, w, lo + k, hi, depth - 1, n_min, out)


# ─────────────────────────────────────────────────────────────────────────────
# 二、边界精修（去混合窗口）
# ─────────────────────────────────────────────────────────────────────────────
def refine_domains(prof: pd.DataFrame,
                   segs: Sequence[Tuple[int, int]],
                   n_min: int = 4,
                   alpha: float = 2.0,
                   max_iter: int = 30) -> List[Tuple[int, int]]:
    """
    逐个年龄域检查"是否与常数年龄模型相容"，不相容就从两侧剥窗口。

    相容判据
    --------
    MSWD 的理论期望是 1，其抽样标准差约为 sqrt(2/(n−1))。
    这里要求：
        MSWD ≤ 1 + α·sqrt(2/(n−1))      （默认 α = 2，约 2σ）
    超过就说明域内有系统性的结构（大概率是混合了边界另一侧的物质）。

    剥哪一侧
    --------
    比较"左边界外第一个点"和"右边界外第一个点"谁离本域均值更远，
    剥掉离得更远的那个方向的一个窗口 —— 因为它更可能是"另一域的物质"。

    参数里的两道"安全阀"
    --------------------
    · n_min：剥到只剩 n_min 个窗口就停止（剥光了就没意义了）
    · max_iter：防止极端情况下来回震荡永不收敛
    """
    a = prof["age68"].to_numpy()
    s = prof["s_age68"].to_numpy()
    n = a.size
    segs = [list(x) for x in segs]        # 转成可变 list 便于原地修改

    for _ in range(max_iter):
        changed = False
        for j, (lo, hi) in enumerate(segs):
            if hi - lo <= n_min:          # 已经太短，不再剥
                continue
            mu, _, mswd, k = weighted_mean(a[lo:hi], s[lo:hi])
            if not np.isfinite(mswd):     # 只有 1 个点时 MSWD 无定义
                continue
            crit = mswd_acceptance(k, alpha)
            if mswd <= crit:              # 通过相容性检验 → 本域 OK
                continue
            # 看看左右边界外侧各有一个候选窗口
            dl = abs(a[lo - 1] - mu) if lo > 0 else -1.0
            dr = abs(a[hi] - mu) if hi < n else -1.0
            if dl < 0 and dr < 0:         # 两侧都没有可剥的了
                continue
            if dl >= dr and lo + 1 < hi:  # 左侧更远 → 左边界右移一格
                segs[j][0] = lo + 1
                changed = True
            elif hi - 1 > lo:             # 右侧更远 → 右边界左移一格
                segs[j][1] = hi - 1
                changed = True
        if not changed:                   # 一轮下来没有任何变化 → 收敛
            break
    # 过滤掉被剥空的域；保留长度 ≥ 1 的
    return [tuple(x) for x in segs if x[1] - x[0] >= 1]


# ─────────────────────────────────────────────────────────────────────────────
# 三、合并"伪分域"（残余漂移导致的假差别）
# ─────────────────────────────────────────────────────────────────────────────
def merge_close(prof: pd.DataFrame,
                segs: Sequence[Tuple[int, int]],
                n_sigma: float = 3.0,
                min_frac: float = 0.05) -> List[Tuple[int, int]]:
    """
    把相邻但"差别不大"的域合并回去。

    为什么要这一步
    --------------
    本批次存在未完全校正的计数率非线性，会在**本来均一**的颗粒上
    造成轻微的年龄剖面倾斜。BIC 对这种微弱但系统性的倾斜很敏感，
    会把一颗均一锆石切成两三个"伪域"。

    两条件同时满足才算真分域：
        · 统计显著：|Δ年龄| > n_sigma × 合成σ
        · 地质显著：|Δ年龄| > min_frac × 年龄（默认 5%）

    为什么还要"地质显著"这一条
    --------------------------
    两个锆石域相差 0.5%、尽管在统计上可能很显著（因为 σ 很小），
    但它对地质解释毫无意义 —— 那可能只是同一个岩浆房里
    结晶时间相差一二十万年的差别。**统计显著 ≠ 地质显著**。
    """
    if len(segs) <= 1:
        return list(segs)
    a = prof["age68"].to_numpy()
    s = prof["s_age68"].to_numpy()

    out = [list(segs[0])]
    for lo, hi in segs[1:]:
        m1, s1, _, _ = weighted_mean(a[out[-1][0]:out[-1][1]], s[out[-1][0]:out[-1][1]])
        m2, s2, _, _ = weighted_mean(a[lo:hi], s[lo:hi])
        # 阈值取两条判据里 **更宽松** 的那个（max），即两条都要跨过
        crit = max(n_sigma * np.hypot(s1, s2), min_frac * abs(m1))
        if abs(m1 - m2) < crit:
            out[-1][1] = hi        # 差别不够 → 并入前一段
        else:
            out.append([lo, hi])
    return [tuple(x) for x in out]


# ─────────────────────────────────────────────────────────────────────────────
# 四、汇总
# ─────────────────────────────────────────────────────────────────────────────
def summarize_segments(prof: pd.DataFrame,
                       segs: Sequence[Tuple[int, int]],
                       mixed_max_frac: float = 0.20) -> pd.DataFrame:
    """
    把分割结果汇总成一张"年龄域表"。

    判定某个域是"过渡带（混合）"还是"真年龄域"的依据：
        · 它是中间段（既非首段也非末段）**且**
          (窗口数 ≤ 全剖面的 mixed_max_frac  **或**  MSWD 显著 > 1)
        · 窗口数 < 3 时无论如何都不算独立年龄域（自由度为 0，不可靠）

    为什么中间的短段默认是混合带
    ----------------------------
    核→边的几何关系决定了：过渡带只可能出现在两个端元之间，
    不可能出现在最浅或最深的位置。中间一个只有 2~3 个窗口的短段，
    几乎必然是跨越边界的混合信号，而不是第三期地质事件。

    返回
    ----
    DataFrame，列：
        domain   域名（D1/D2…，混合带为 "—"）
        i0/i1    窗口下标范围
        tau0/tau1 归一化深度范围
        n_win    窗口数
        age_Ma   加权平均年龄
        se_1sig  加权平均的 1σ
        mswd     加权偏差均方
        mswd_prob 卡方上尾概率 —— "该域只有一个年龄"这一假设下，
                 MSWD 至少这么大的概率。读法与 R 端 ADEPT 的
                 `MSWD probability` 完全一致，两边可比。
        flag     判定依据与结果
        ThU/U_cps 该域的平均 Th/U 与 U 信号（判岩性用）
    """
    rows = []
    a = prof["age68"].to_numpy()
    s = prof["s_age68"].to_numpy()
    n_tot = a.size
    nseg = len(segs)

    for j, (lo, hi) in enumerate(segs):
        mu, se, mswd, k = weighted_mean(a[lo:hi], s[lo:hi])
        crit = mswd_acceptance(k)
        interior = 0 < j < nseg - 1                      # 是否为中间段
        mixed = (interior and (k <= mixed_max_frac * n_tot or mswd > crit)) or (k < 3)
        # ADEPT 用 pf(mswd, k-1, Inf, lower.tail=FALSE)，等价于
        # chi2_sf(mswd*(k-1), k-1) —— **要乘 (k-1)**，见 core.statistics.chi2_sf。
        prob = chi2_sf(mswd * (k - 1), k - 1) if k >= 2 else float("nan")
        rows.append(dict(
            domain=f"D{j + 1}",
            i0=lo, i1=hi,
            tau0=float(prof["tau"].iloc[lo]),
            tau1=float(prof["tau"].iloc[hi - 1]),
            n_win=k, age_Ma=mu, se_1sig=se, mswd=mswd, mswd_prob=prob,
            flag="mixed/过渡带" if mixed else "age domain",
            ThU=float(prof["ThU"].iloc[lo:hi].mean()),
            U_cps=float(prof["U_cps"].iloc[lo:hi].mean()),
        ))

    df = pd.DataFrame(rows)
    # 剔除混合带之后重新编号：D1、D2 应该是按地质顺序排列的年龄域
    keep = df["flag"] == "age domain"
    df.loc[keep, "domain"] = [f"D{i + 1}" for i in range(int(keep.sum()))]
    df.loc[~keep, "domain"] = "—"
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 五、不分域：把整段当作一个域
# ─────────────────────────────────────────────────────────────────────────────
def whole_spot_stats(prof: pd.DataFrame, alpha: float = 2.0) -> dict:
    """
    **不分域**的整段年龄 —— 刻意不做 BIC 分割、不做边界精修、不做伪域合并，
    把全部窗口当作一个域，用与域级**完全相同的口径**（反比方差加权平均 +
    χ² 上尾）算出一个年龄。

    它回答的不是"这个颗粒几岁"，而是：

        "如果我不做任何分域，这个测点会报出什么？"

    有了它，`剖面窗口`（逐窗口）→ 本函数（整段一个数）→「深度剖面域」
    （分域后的各个数）就是一条完整可对照的链：**相邻两级之差，正好是那一级
    处理步骤的全部贡献**，中间不掺任何口径差异。这也是它对"跨测点横向对比"
    特别有用的原因 —— 全部测点走同一条路，不再有的点多一行、有的点没有行
    （`深度剖面域` 只收多域点，均一点在里面是没有行的）。

    ⚠ 读法（这条必须传到下游）
    --------------------------
    `mswd` 明显大于 1 时，这个年龄**不是任何一期地质事件的年龄**，而是几个
    年龄按权重的混合值。实测示例批次（窗口 4 s / 步长 1 s）48 个样品测点里
    只有 5 个（10%）的整段 MSWD ≤ 1.5，中位 4.4，最大的到 196：

        mswd_ok 为真 → 整段确实只有一个年龄，这个数可以直接用；
        mswd_ok 为假 → 只可用于横向对比，**不可用于定年**。

    另有一个次要成因：窗口重叠让相邻窗口不独立，会把 MSWD 系统性压小一点。
    所以这个判据偏宽松 —— 真报出 MSWD 大，那就一定更大。

    参数
    ----
    prof  : window_profile 结果，需含 age68 / s_age68（即已经过 F(τ) 校正）
    alpha : 相容判据的倍数，含义同 refine_domains，默认 2

    返回
    ----
    dict，字段：
        n_win     参与加权平均的窗口数（weighted_mean 已剔除 NaN 与非正 σ）
        tau0/tau1 全部窗口的归一化深度范围
        age_Ma    整段加权平均年龄
        se_1sig   该均值的 1σ（**未**按 MSWD 膨胀，与域级表同口径、可直接相减）
        mswd      加权偏差均方
        mswd_prob 卡方上尾概率（与 ADEPT 的 MSWD probability 同口径）
        mswd_crit 相容上限 = mswd_acceptance(n_win, alpha)
        mswd_ok   整段是否与"单一常数年龄"相容
        drift_Ma  末窗口 − 首窗口的年龄差，残余倾斜的直接度量
    """
    a = prof["age68"].to_numpy(float)
    s = prof["s_age68"].to_numpy(float)
    mu, se, mswd, k = weighted_mean(a, s)
    crit = mswd_acceptance(k, alpha)
    ok = bool(np.isfinite(mswd) and mswd <= crit)
    # ADEPT 用 pf(mswd, k-1, Inf, lower.tail=FALSE)，等价于 chi2_sf(mswd*(k-1), k-1)
    # —— **要乘 (k-1)**，与 summarize_segments 同一个口径。
    prob = chi2_sf(mswd * (k - 1), k - 1) if k >= 2 else float("nan")
    tau = prof["tau"].to_numpy(float)
    # 漂移只用**参与加权平均的那批窗口**的首末，这样才能与 mswd 自洽
    # （若首窗口的年龄是 NaN，直接取 a[0] 会得到一个假的 nan 漂移）
    keep = np.isfinite(a) & np.isfinite(s) & (s > 0)
    av = a[keep]
    return dict(
        n_win=k,
        tau0=float(tau[0]) if tau.size else float("nan"),
        tau1=float(tau[-1]) if tau.size else float("nan"),
        age_Ma=mu, se_1sig=se,
        mswd=mswd, mswd_prob=prob, mswd_crit=crit, mswd_ok=ok,
        drift_Ma=float(av[-1] - av[0]) if av.size >= 2 else float("nan"),
    )
