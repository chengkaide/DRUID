"""
druid 的深度剖面自检 —— 分域（`depth/domains.py`）与分馏（`depth/fractionation.py`）：

    python tests/test_depth.py
    python -m pytest tests/            （装了 pytest 时）

这两块之前**一个测试都没有**，而它们是"从剖面解出几个年龄"这一步的全部逻辑。

分域这条链是四步流水，每一步都可能把结果静默改掉：

    segment            BIC 二叉分割      → 下标区间
    refine_domains     剥边界混合窗口    → 下标区间
    merge_close        合并伪分域        → 下标区间
    summarize_segments 定"真域 / 过渡带" → 域表

本文件按这个顺序逐步钉住。重点在两处**判据**（它们决定了报出去几个年龄）：
· `mswd_acceptance`：全模块唯一的"相容上限"，两个地方各写过一遍（现在只此一处）；
· `merge_close` 的**两条**判据（统计显著 **且** 地质显著）——只满足一条不算分域。

fractionation 那边最重要的一条：`profile_ages` 的 σ 是**两路正交合成**
（内部 σ68·F 与外部 σ_ext·R68），漏掉任何一路都会让误差棒偏小。
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                              # noqa: E402
import pandas as pd                                             # noqa: E402

from druid.core.constants import L238                           # noqa: E402
from druid.core.geochronology import age76                      # noqa: E402
from druid.depth.domains import (merge_close, mswd_acceptance,   # noqa: E402
                                 refine_domains, segment,
                                 summarize_segments, whole_spot_stats)
from druid.depth.fractionation import (bracket_F, profile_ages,  # noqa: E402
                                       profile_ages_76)


def _prof(ages, sig=5.0):
    """最小可用剖面：只要 age68 / s_age68 / tau（分域链路读的就是这三列）。"""
    n = len(ages)
    return pd.DataFrame(dict(
        age68=np.array(ages, float), s_age68=np.full(n, float(sig)),
        tau=np.linspace(0.0, 1.0, n), ThU=np.full(n, 0.5),
        U_cps=np.full(n, 1.0e6)))


# ═════════════════════════════════════════════════════════════════════════════
# 一、相容判据：全模块唯一定义
# ═════════════════════════════════════════════════════════════════════════════
def test_mswd_acceptance_formula():
    """
    `MSWD ≤ 1 + α·√(2/(n−1))`，α 默认 2。

    两个容易忽略的细节：
    · n 用的是**窗口数**，而窗口是重叠的 ⇒ 真实独立观测数更少 ⇒ 上限偏**宽松**
      （方向是漏杀混合点，不会错杀真域）；
    · n ≤ 2 时 `max(n−1, 1)` 把分母钉在 1，所以 n=1 与 n=2 得上限相同 ——
      不能让 n=1 走到除以 0。
    """
    assert abs(mswd_acceptance(5, 2.0) - (1.0 + 2.0 * math.sqrt(2.0 / 4.0))) < 1e-15
    assert abs(mswd_acceptance(40) - (1.0 + 2.0 * math.sqrt(2.0 / 39.0))) < 1e-15
    assert mswd_acceptance(1) == mswd_acceptance(2), "n≤2 都应落在 max(n−1,1)=1 上"
    assert abs(mswd_acceptance(1) - (1.0 + 2.0 * math.sqrt(2.0))) < 1e-15
    # α 要真的参与运算（防"传了没用"）
    assert mswd_acceptance(10, 3.0) > mswd_acceptance(10, 1.0)
    # 上限随 n 单调收紧，且恒 > 1
    vals = [mswd_acceptance(n) for n in range(3, 60)]
    assert all(b < a for a, b in zip(vals, vals[1:]))
    assert all(v > 1.0 for v in vals)


# ═════════════════════════════════════════════════════════════════════════════
# 二、BIC 二叉分割
# ═════════════════════════════════════════════════════════════════════════════
def test_segment_finds_a_clean_two_level_step():
    """一个干净的两级台阶 → 恰好两段，边界落在台阶上。"""
    segs = segment(np.array([500.0] * 15 + [300.0] * 15), np.full(30, 5.0),
                   n_min=3, max_segs=3)
    assert segs == [(0, 15), (15, 30)], segs


def test_segment_leaves_uniform_data_alone():
    """均一剖面**不许**被切 —— BIC 的复杂度罚项就是为这件事存在的。"""
    assert segment(np.full(30, 400.0), np.full(30, 5.0)) == [(0, 30)]


def test_segment_respects_minimum_and_maximum():
    """
    两道闸门：样本少于 `2·n_min` 直接返回 []（调用方视作"整段一个域"）；
    每一段的长度都不许低于 n_min。
    """
    assert segment(np.array([400.0, 410.0]), np.array([5.0, 5.0]), n_min=3) == []
    segs = segment(np.array([600.0] * 10 + [400.0] * 10 + [200.0] * 10),
                   np.full(30, 5.0), n_min=3, max_segs=3)
    assert segs == [(0, 10), (10, 20), (20, 30)], segs
    assert all(hi - lo >= 3 for lo, hi in segs)
    # n_min 调大到一段装不下 → 只能少切
    assert len(segment(np.array([600.0] * 10 + [400.0] * 10 + [200.0] * 10),
                       np.full(30, 5.0), n_min=12, max_segs=3)) == 1


def test_segment_recursion_depth_is_what_limits_the_count():
    """
    `max_segs` 是通过**递归深度**生效的，不是"切完再合并"。
    所以同一批三级台阶：深度给 1 时只切一刀，而且切在**第一次 BIC 最优**的位置
    （这里是 20，不是 10）—— 这个细节决定了"只允许两域"时切在哪里。
    """
    a = np.array([600.0] * 10 + [400.0] * 10 + [200.0] * 10)
    two = segment(a, np.full(30, 5.0), n_min=3, max_segs=2)
    three = segment(a, np.full(30, 5.0), n_min=3, max_segs=3)
    assert len(two) == 2 and len(three) == 3
    assert two == [(0, 20), (20, 30)], two
    assert three == [(0, 10), (10, 20), (20, 30)], three


def test_segment_indices_are_into_the_original_array():
    """
    ★ **D-1 修掉的那条口径（2026-10-02）：`segment` 返回的是原数组下标。**

    它内部会把非有限年龄 / 非正 σ 的窗口滤掉，但**返回前映射回原下标**。
    这条契约必须钉住，因为调用方 `refine_domains(prof, segs)` 用的是
    **未过滤**的 prof —— 两处口径一错位，只要剖面里有 NaN，此后所有下标
    就会**整体左移且不报错**（旧实现正是如此：8 个点、第 2 个 NaN → `[(0, 7)]`）。

    实测全库 26314 个窗口**零 NaN** ⇒ 修它在本库上是**恒等映射、基线不变**；
    但"碰巧没有 NaN"不是保证（`bracket_F` 在标样比值 ≤ 0 处会给 NaN）。
    """
    ages = np.array([400.0, np.nan, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0])
    segs = segment(ages, np.full(8, 5.0))
    # 滤后有 7 个点（原下标 0, 2…7）⇒ 覆盖范围是原数组的 [0, 8)
    assert segs == [(0, 8)], segs
    assert segs[0][1] == ages.size, "下标必须是原数组的，不是滤后长度"

    # 被滤掉的点夹在两段之间：它不属于任何一段，但**不会**让后一段的下标左移
    ages2 = np.array([600.0] * 4 + [np.nan] + [200.0] * 4)
    segs2 = segment(ages2, np.full(9, 5.0), n_min=3)
    assert segs2 == [(0, 4), (5, 9)], segs2
    assert segs2[-1][1] == ages2.size, "末段必须覆盖到原数组末尾"

    # 全 NaN / 全零 σ 也不能炸，只是返回 []
    assert segment(np.full(8, np.nan), np.full(8, 5.0)) == []
    assert segment(np.full(8, 400.0), np.zeros(8)) == []


def test_segment_and_refine_domains_share_one_index_convention():
    """
    D-1 的**端到端回归**：`segment` 的输出要能直接喂给 `refine_domains`，
    即使剖面里有 NaN。改坏任何一侧的下标口径，这条都会红。
    """
    n = 12
    prof = _prof([500.0] * 6 + [np.nan] + [300.0] * 5)
    segs = segment(prof["age68"].to_numpy(), prof["s_age68"].to_numpy(), n_min=3)
    assert segs, "两个明显的台阶应当被切出来"
    assert all(0 <= lo < hi <= n for lo, hi in segs), ("segment 越界", segs, n)
    assert segs[-1][1] == n, "末段必须覆盖到 prof 的最后一行"
    # 两者同一套下标 ⇒ refine_domains 不会越界、也不会静默错位
    refined = refine_domains(prof, segs)
    assert refined, refined
    assert all(0 <= lo < hi <= n for lo, hi in refined), ("refine 越界", refined, n)
    assert refined[-1][1] <= n


# ═════════════════════════════════════════════════════════════════════════════
# 三、边界精修
# ═════════════════════════════════════════════════════════════════════════════
def test_refine_domains_peels_mixed_windows():
    """
    段内混进 2 个异类窗口 ⇒ MSWD 远超上限 ⇒ 从**离均值更远的那一侧**逐个剥掉，
    剥到相容为止（本例正好剥掉那 2 个）。

    注意剥的方向只可能朝"本段之外"：贴着数组两端的段没有可剥的邻居，
    所以 refine 对"整段一个域"是**无能为力**的（这正是它后面还要 merge_close 的原因）。
    """
    p = _prof([500.0] * 12 + [300.0] * 12)
    out = refine_domains(p, [(0, 14), (14, 24)], n_min=4)
    assert out == [(0, 12), (14, 24)], out


def test_refine_domains_leaves_a_clean_split_alone():
    """本来就干净的切分必须**原样返回** —— 精修不是"总会改动"，它是条件触发。"""
    p = _prof([500.0] * 12 + [300.0] * 12)
    assert refine_domains(p, [(0, 12), (12, 24)], n_min=4) == [(0, 12), (12, 24)]


def test_refine_domains_cannot_peel_a_whole_spot_segment():
    """整段一个域、且段内有异类 → MSWD 很大但**剥不动**（两侧都没邻居）。"""
    p = _prof([500.0] * 10 + [300.0] * 2 + [500.0] * 10)
    assert refine_domains(p, [(0, 22)], n_min=4) == [(0, 22)]


def test_refine_domains_respects_n_min():
    """`n_min` 是安全阀：段长已经 ≤ n_min 就不再剥，免得把域剥没了。"""
    p = _prof([500.0] * 12 + [300.0] * 12)
    assert refine_domains(p, [(0, 14), (14, 24)], n_min=20) == [(0, 14), (14, 24)]


# ═════════════════════════════════════════════════════════════════════════════
# 四、合并伪分域（两条判据同时满足才算真分域）
# ═════════════════════════════════════════════════════════════════════════════
def test_merge_close_merges_statistically_insignificant_differences():
    """
    相差 10 Ma / 500 Ma = 2%，统计上可能很显著（σ 很小），但**地质上没意义** ⇒
    必须合并。`min_frac × 年龄`（默认 5% = 25 Ma）这一条就是为它设的。
    """
    p = _prof([500.0] * 15 + [490.0] * 15)
    assert merge_close(p, [(0, 15), (15, 30)]) == [(0, 30)]


def test_merge_close_keeps_real_differences():
    """相差 40% ⇒ 两条判据都过 ⇒ 保留两个域。"""
    p = _prof([500.0] * 15 + [300.0] * 15)
    assert merge_close(p, [(0, 15), (15, 30)]) == [(0, 15), (15, 30)]


def test_merge_close_needs_both_criteria():
    """
    ★ 两条判据是**与**关系，阈值取更宽松的那个（max）。
    所以：
      · 若把 min_frac 调到 0，唯一的门槛就是统计显著（3σ 合成）——
        这时 10 Ma 的差别（3·hypot(1.29,1.29) ≈ 5.5 Ma）就该被**保留**；
      · 反过来把 min_frac 调到 0.5（250 Ma），40% 的差别反而会被**合并**。
    这条同时证明两个参数都真的参与运算。
    """
    p_small = _prof([500.0] * 15 + [490.0] * 15)
    assert merge_close(p_small, [(0, 15), (15, 30)], min_frac=0.0) == [(0, 15), (15, 30)]
    p_big = _prof([500.0] * 15 + [300.0] * 15)
    assert merge_close(p_big, [(0, 15), (15, 30)], min_frac=0.5) == [(0, 30)]


def test_merge_close_single_segment_is_a_noop():
    p = _prof([500.0] * 30)
    assert merge_close(p, [(0, 30)]) == [(0, 30)]


def test_merge_close_cascades_left():
    """
    合并是**向左串**的：三段 a/b/c，b 与前段相近就并进前段，
    于是 c 比较的对象变成"a+b 的合并段"而不是 b。这个顺序会影响结果，
    必须钉住（否则"总是合并相邻两段"的直觉写法会给出不同答案）。
    """
    ages = [500.0] * 10 + [495.0] * 10 + [300.0] * 10
    out = merge_close(_prof(ages), [(0, 10), (10, 20), (20, 30)])
    assert out == [(0, 20), (20, 30)], out


# ═════════════════════════════════════════════════════════════════════════════
# 五、汇总：真域 vs 过渡带
# ═════════════════════════════════════════════════════════════════════════════
def test_summarize_flags_interior_short_segments_as_mixed():
    """
    中间的短段（既非首段也非末段，且窗口数少于全剖面 20%）判为过渡带。

    为什么中间段一律可疑：核→边的几何决定了混合带只可能夹在两个端元之间，
    中间一个只有 2~3 个窗口的短段几乎必然是"跨边界"的混合信号。
    """
    p = _prof([600.0] * 10 + [400.0] * 3 + [200.0] * 10)
    df = summarize_segments(p, [(0, 10), (10, 13), (13, 23)])
    assert list(df["flag"]) == ["age domain", "mixed/过渡带", "age domain"]


def test_summarize_keeps_both_end_segments_even_if_short():
    """首段 / 末段再短也**不算**过渡带（判据里那个 `interior` 的作用）。"""
    p = _prof([600.0] * 3 + [500.0] * 10 + [200.0] * 10)
    df = summarize_segments(p, [(0, 3), (3, 13), (13, 23)])
    assert list(df["flag"]) == ["age domain", "age domain", "age domain"]


def test_summarize_treats_any_segment_shorter_than_three_as_mixed():
    """窗口数 < 3 时**无论在哪**都不算独立年龄域 —— k=2 的自由度只有 1，撑不住。"""
    p = _prof([600.0] * 10 + [400.0] * 2 + [200.0] * 11)
    df = summarize_segments(p, [(0, 10), (10, 12), (12, 23)])
    assert list(df["flag"]) == ["age domain", "mixed/过渡带", "age domain"]


def test_summarize_renumbers_only_real_domains():
    """
    域号按**地质顺序**重排（D1、D2…），过渡带拿 "—"。
    所以"表中最大编号 = 真域数"，下游按这个数报"多域(n)"。
    """
    p = _prof([600.0] * 10 + [400.0] * 3 + [200.0] * 10)
    df = summarize_segments(p, [(0, 10), (10, 13), (13, 23)])
    assert list(df["domain"]) == ["D1", "—", "D2"]
    assert int((df["flag"] == "age domain").sum()) == 2


def test_summarize_probability_matches_the_adept_convention():
    """
    `mswd_prob = chi2_sf(mswd·(k−1), k−1)` —— **要乘 (k−1)**。
    与 R 端 ADEPT 的 `pf(mswd, k-1, Inf, lower.tail=FALSE)` 同口径，两边可比。

    这条专门防"忘了乘自由度"：漏乘会让概率值系统性偏小（把好域判成坏域）。
    """
    from druid.core.statistics import chi2_sf
    p = _prof([600.0] * 10)                       # 常数年龄 → MSWD = 0
    df = summarize_segments(p, [(0, 10)])
    row = df.iloc[0]
    assert abs(row["mswd"]) < 1e-12
    assert abs(row["mswd_prob"] - 1.0) < 1e-12    # chi2_sf(0, 9) = 1
    assert abs(row["se_1sig"] - 5.0 / math.sqrt(10)) < 1e-12

    # 造一个 MSWD 显著 > 1 的段，逐位对回 chi2_sf(mswd*(k-1), k-1)
    q = _prof([500.0, 400.0, 600.0, 450.0, 550.0])
    df2 = summarize_segments(q, [(0, 5)])
    r2 = df2.iloc[0]
    assert abs(r2["mswd_prob"] - chi2_sf(r2["mswd"] * 4, 4)) < 1e-15


# ═════════════════════════════════════════════════════════════════════════════
# 六、不分域（整段）年龄
# ═════════════════════════════════════════════════════════════════════════════
def test_whole_spot_stats_reports_incompatibility():
    """
    整段一个域的口径：不分割、不精修、不合并，只为给出"如果我不分域会报什么"。
    两个内部混合的异类窗口就足以让 MSWD 爆表 ⇒ `mswd_ok` 必须是 False。

    这条的用途写在返回值的读法里：`mswd_ok` 为假时，这个年龄**不是任何一期
    地质事件的年龄**，只能用于横向对比。
    """
    p = _prof([500.0] * 10 + [300.0] * 2 + [500.0] * 10)
    w = whole_spot_stats(p)
    assert w["n_win"] == 22
    assert w["mswd"] > 100.0
    assert w["mswd_ok"] is False
    assert w["mswd_prob"] == 0.0
    assert abs(w["mswd_crit"] - mswd_acceptance(22)) < 1e-15
    assert abs(w["age_Ma"] - 481.8181818) < 1e-6   # 被两个 300 Ma 拉低
    assert w["drift_Ma"] == 0.0                    # 首末都是 500 → 无倾斜

    # 均一剖面则相反
    ok = whole_spot_stats(_prof([450.0] * 20))
    assert ok["mswd_ok"] is True and abs(ok["mswd"]) < 1e-12


# ═════════════════════════════════════════════════════════════════════════════
# 七、分馏因子 F(τ)
# ═════════════════════════════════════════════════════════════════════════════
def test_bracket_F_is_reference_over_measured():
    """
    `F(τ) = R_ref / R_measured(τ)` —— 标样在同一深度偏了多少，样品就补多少。
    标样比值随深度**下降**时，F 必须随之**上升**（这正是分馏被校正的方向）。
    """
    std = pd.DataFrame(dict(tau=[0.0, 0.5, 1.0], R68=[0.1700, 0.1680, 0.1660]))
    F = bracket_F(np.array([0.0, 0.25, 1.0]), [std], "R68", 0.16780)
    assert abs(F[0] - 0.16780 / 0.1700) < 1e-12
    assert abs(F[1] - 0.16780 / 0.1690) < 1e-12     # 线性插值后的中点
    assert abs(F[2] - 0.16780 / 0.1660) < 1e-12
    assert all(b > a for a, b in zip(F, F[1:])), F


def test_bracket_F_returns_nan_for_nonpositive_measured():
    """
    标样比值 ≤ 0 没有物理意义（负比值）⇒ 该 τ 处的 F 必须是 **nan 而不是 inf**，
    让上层按"这个窗口无效"处理。返回 inf 会把年龄变成 0/NaN 更难追。

    ⚠ 注意：只有当**目标 τ 网格正好落在那个坏点附近**时才会命中 ——
    插值会把它"绕过去"。所以造样例时目标网格要盖住它。
    """
    bad = pd.DataFrame(dict(tau=[0.0, 0.5, 1.0], R68=[0.1700, -0.01, 0.1660]))
    F = bracket_F(np.array([0.0, 0.5, 1.0]), [bad], "R68", 0.16780)
    assert np.isnan(F[1])
    assert np.isfinite(F[0]) and np.isfinite(F[2])


def test_bracket_F_averages_multiple_standards_per_tau():
    """多个夹逼标样在**每个 τ 上取算术平均**（一阶抵消仪器漂移）。"""
    std = pd.DataFrame(dict(tau=[0.0, 0.5, 1.0], R68=[0.1700, 0.1680, 0.1660]))
    assert abs(float(bracket_F(np.array([0.5]), [std], "R68", 0.5)[0])
               - 0.5 / 0.1680) < 1e-12
    two = bracket_F(np.array([0.5]), [std, std], "R68", 0.5)
    assert abs(float(two[0]) - 0.5 / 0.1680) < 1e-12     # 两份一样的 → 平均不变
    half = bracket_F(np.array([0.5]), [std], "R68", 0.084)
    assert abs(float(half[0]) * 2.0 - float(bracket_F(np.array([0.5]), [std],
                                                      "R68", 0.168)[0])) < 1e-12


# ═════════════════════════════════════════════════════════════════════════════
# 八、窗口年龄与 σ 合成
# ═════════════════════════════════════════════════════════════════════════════
def _one_window(r68=0.06, s68=0.001, r76=0.055, s76=0.003):
    return pd.DataFrame(dict(R68=[r68], s68=[s68], R76=[r76], s76=[s76]))


def test_profile_ages_sigma_is_two_orthogonal_sources():
    """
    σ 的合成规则（本文件最容易被漏掉的一处）：

        σ_age = hypot(σ68·F, σ_ext·R68) / (λ238·(1+R68·F)) / 1e6

    两路来源：① 本窗口的计数/闪烁噪声（随 F 一起缩放）；
              ② 逐点之间的外部重现性（相对量 × 校正后的比值）。
    漏掉②会让误差棒系统性偏小 —— 而"误差棒刚好小一点"是看不出来的。

    这里逐位对回手算值，并且验证 F ≠ 1 时两路**都跟着 F 缩放**。
    """
    p = _one_window()
    F = np.array([1.0])
    a, sa = profile_ages(p, F, 0.007)
    R = 0.06 * 1.0
    want_a = math.log1p(R) / L238 / 1e6
    want_s = math.hypot(0.001 * 1.0, 0.007 * R) / (L238 * (1.0 + R)) / 1e6
    assert abs(float(a[0]) - want_a) < 1e-12
    assert abs(float(sa[0]) - want_s) < 1e-15, (sa[0], want_s)
    # 去掉外部项会明显变小 —— 证明第二路真的在起作用
    assert want_s > 0.001 / (L238 * (1.0 + R)) / 1e6
    # F 变化时年龄与 σ 同步变化
    a2, sa2 = profile_ages(p, np.array([1.05]), 0.007)
    assert float(a2[0]) > float(a[0])
    assert float(sa2[0]) != float(sa[0])


def test_profile_ages_76_matches_the_scalar_definition():
    """
    ↑ 207Pb/206Pb 版本的年龄与 σ：

        age = age76(R76·F76)
        σ   = |age76(R+s) − age76(R−s)| / 2 ，  s = hypot(s76·F76, σ_ext76·R76)

    为什么 σ 要用"上下各偏 s 再取半差"而不是解析式：207/206 的年龄方程是**隐式**的
    （`age76` 内部是二分法），没有解析导数，用差分是最稳的近似。

    ⚠ 这条也顺便钉住了一件事：**复合 s 里含外部项**（σ_ext76·R76）。
    只拿 s76 去差分得到的 σ 会偏小（本例 122.2 vs 正确的 123.2）。
    """
    p = _one_window()
    ag, sa = profile_ages_76(p, np.array([1.0]), 0.007)
    R = 0.055
    assert abs(float(ag[0]) - float(age76(R))) < 1e-12
    s = math.hypot(0.003 * 1.0, 0.007 * R)
    want = abs(age76(R + s) - age76(R - s)) / 2.0
    assert abs(float(sa[0]) - want) < 1e-9, (sa[0], want)
    assert float(sa[0]) > abs(age76(R + 0.003) - age76(R - 0.003)) / 2.0

    # F76 ≠ 1 时年龄跟着走（注意 F76 会先乘到比值上）
    ag2, _ = profile_ages_76(p, np.array([1.05]), 0.007)
    assert abs(float(ag2[0]) - float(age76(R * 1.05))) < 1e-12


def _run_standalone() -> int:
    """见 tests/_selftest.py —— 让这个文件不装 pytest 也能直接跑。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _selftest
    return _selftest.run(globals())


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
