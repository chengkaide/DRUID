"""
druid 的心脏：`reduction.ratios.reduce_interval` 契约自检 —— 不依赖 pytest：

    python tests/test_reduction.py
    python -m pytest tests/            （装了 pytest 时）

为什么这个文件最重要
--------------------
`reduce_interval` 把一段净信号还原成五个比值 + σ + ρ，全包所有年龄都从这里出来。
它有两处"**错得很安静**"的地方：

① **204Pb 显著性检验**。真实数据里 204Pb 只有几~几十 cps，而 204Hg 本底几百 cps，
   204Pb 是"两个大数相减"。不加检验一律扣普通铅的话，会把噪声当铅扣掉 ——
   206Pb 影响不大，**207Pb 被扣的比例大得多** ⇒ 207Pb/206Pb 系统性偏低 ⇒
   年龄大幅偏年轻。这条判据是"数据能不能用"的分水岭，所以必须钉住
   "不显著就不扣"这件事。
② **jackknife 的 σ 与 ρ**。σ 若写成 `Σ(x−x̄)²/n`（少了 (n−1)/n 修正）会偏小，
   症状是"误差棒刚好小一点点"，看不出问题；ρ 若算错，协和椭圆会画歪。

策略：用**合成信号**（自己造 net 字典）精确控制"有没有普通铅、噪声多大"，
所以每条断言测的都是一个明确契约，不依赖真实批次。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                              # noqa: E402

from druid.core.constants import HG202_204, MASSES_NEEDED       # noqa: E402
from druid.reduction.ratios import reduce_interval              # noqa: E402


def _synth(n=20, seed=1, excess204=0.0, noise=0.0,
           r68=0.06, r76=0.055, i238=1.0e6, i202=100.0):
    """
    造一个测点的净信号字典。

    excess204 : 叠加在"Hg 扣净后"的**真实普通铅** 204Pb 强度上（cps）。
                0 表示只有 Hg，没有普通铅。
    noise     : 每个通道加的同方差高斯噪声（cps），用来让 jackknife 有东西可估。

    注意 202 通道是**按几何关系**给的：204 里的 Hg 份额恰好是 202/4.35，
    所以扣完之后 i204_raw 正好等于 excess204（再叠加噪声）。
    """
    rng = np.random.default_rng(seed)

    def noisy(v):
        base = np.full(n, float(v))
        return base + (rng.normal(0.0, noise, n) if noise > 0 else 0.0)

    i206 = i238 * r68
    i207 = i206 * r76
    return {238: noisy(i238), 206: noisy(i206), 207: noisy(i207),
            208: noisy(i206 * 0.02), 232: noisy(i206 * 0.04),
            202: noisy(i202), 204: noisy(i202 / HG202_204 + excess204)}


def _mask(n, k=None):
    """前 k 个时间点属于本窗口（缺省全选）。"""
    m = np.zeros(n, bool)
    m[:n if k is None else k] = True
    return m


def test_channel_set_matches_constants():
    """通道全集只认 `constants.MASSES_NEEDED` 这一处定义（不许本文件再写一份）。"""
    assert MASSES_NEEDED == (202, 204, 206, 207, 208, 232, 238)
    assert 235 not in MASSES_NEEDED, "235U 不测，207/235 靠 238/235 换算"


def test_missing_channel_raises_with_the_offender_named():
    """
    少任何一个通道都必须**当场报错**，而且要说清是哪一个。
    静默算错年龄是本文件全篇要防的头号事故。
    """
    net = _synth()
    del net[207]
    try:
        reduce_interval(net, _mask(20))
    except KeyError as e:
        assert "207" in str(e), str(e)
    else:
        raise AssertionError("缺通道居然没报错")


def test_too_few_points_returns_none():
    """
    窗口内有效点数 < 5 → 返回 None（不是空 dict、不是抛异常）。

    为什么是 5：jackknife 要留掉一个点还有 ≥4 个，delete-1 的方差估计才勉强
    可用。这是与调用方（滑窗）之间的契约 —— 它按 `is None` 决定这个窗口要不要。
    """
    net = _synth(n=20)
    for k in (0, 1, 4):
        assert reduce_interval(net, _mask(20, k)) is None, k
    assert reduce_interval(net, _mask(20, 5)) is not None, "正好 5 个点应当可用"


def test_no_common_lead_skips_the_correction():
    """
    ★ 本文件最重要的一条：**204Pb 不显著时不许扣普通铅**。

    判据是 `i204 > n_sigma · σ(i204)`。没有显著 204Pb 时（这里连噪声都没加，
    所以 σ = 0、必然不显著），必须：
        · 不做扣除（f206 = 0、sk 全 0、i204 = 0）；
        · R68 / R76 就是**原始比值**，一个数都不改。
    """
    net = _synth(n=20, excess204=0.0, noise=0.0)
    d = reduce_interval(net, _mask(20))
    assert d["i204_significant"] is False
    assert d["i204"] == 0.0
    assert d["f206"] == 0.0
    assert d["sk"] == (0.0, 0.0, 0.0)
    assert abs(d["R68"] - (net[206].mean() / net[238].mean())) < 1e-15
    assert abs(d["R76"] - (net[207].mean() / net[206].mean())) < 1e-15


def test_significant_common_lead_is_subtracted():
    """
    204Pb 显著时：扣掉、给出 f206 > 0，并改变两个比值。

    ⚠ **方向别记反**。这里造的是 460 Ma 量级的年轻样品（R76 ≈ 0.055），
    而普通铅自己的 `207/206 = sk[1]/sk[0] ≈ 15.60/18.14 ≈ 0.860` ——
    **远高于**样品。从混合物里扣掉一个"高比值组分"，残差比值必然**下降**。
    同时 207Pb 被扣掉的**相对份额**也大得多（本例 2.2% vs 206Pb 的 0.14%，
    因为 207Pb 本来就少），两个理由指向同一个方向。

    所以这里不写死"变大/变小"，而是**按扣除组分相对样品偏高还是偏低**判方向：
    这样换成老锆石（R76 接近甚至超过 0.86）时断言依然成立。
    """
    net = _synth(n=20, seed=1, excess204=5.0, noise=1.0)
    d = reduce_interval(net, _mask(20))
    raw = reduce_interval(net, _mask(20), n_sigma=1e9)      # 门槛抬到天上 = 不扣

    assert d["i204_significant"] is True
    assert d["i204"] > 0.0
    assert d["f206"] > 0.0
    assert raw["f206"] == 0.0

    sk76 = d["sk"][1] / d["sk"][0]                          # 普通铅的 207/206
    assert abs(sk76 - 0.860) < 0.01, sk76
    if sk76 > raw["R76"]:
        assert d["R76"] < raw["R76"], (d["R76"], raw["R76"])
    else:
        assert d["R76"] > raw["R76"], (d["R76"], raw["R76"])

    # 206Pb 被扣掉的**相对份额**远小于 207Pb（这才是不做显著性检验的代价所在）
    #   frac6 = c6/i206 = 1 − R68_corr/R68_raw
    #   R76_corr/R76_raw = (1−frac7)/(1−frac6)
    frac6 = 1.0 - d["R68"] / raw["R68"]
    frac7 = 1.0 - (d["R76"] / raw["R76"]) * (1.0 - frac6)
    assert 0.0 < frac6 < frac7, (frac6, frac7)


def test_significance_threshold_is_a_real_gate():
    """
    n_sigma 是个**真的闸门**，不是摆设：同一个信号，门槛 2σ 判"有铅"、
    门槛 50σ 判"没铅"，两条路给出的结果必须分别落到两个分支上。

    这条专门防"阈值传下去了但没参与判断"（把 `sig204` 写成常量、
    或者忘了把 `n_sigma` 用上）。
    """
    net = _synth(n=20, seed=1, excess204=5.0, noise=1.0)
    low = reduce_interval(net, _mask(20), n_sigma=2.0)
    high = reduce_interval(net, _mask(20), n_sigma=50.0)
    assert low["i204_significant"] is True and low["f206"] > 0
    assert high["i204_significant"] is False and high["f206"] == 0.0


def test_jackknife_sigma_scales_as_one_over_sqrt_n():
    """
    delete-1 jackknife 的 σ 是"均值的标准误"，所以 **σ ∝ 1/√n**。
    噪声固定、把窗口拉长 4 倍，σ 必须掉到约 1/2。

    这条能同时抓住：漏乘修正因子 (n−1)/n、以及把"求和"写成"求和÷n"
    之类会让 σ 不随 n 变的错法。
    """
    s68 = {}
    for n in (20, 80, 320):
        d = reduce_interval(_synth(n, seed=7, noise=1.0), _mask(n))
        assert d["s68"] > 0 and d["s76"] > 0
        s68[n] = d["s68"]
    for small, big in ((20, 80), (80, 320)):
        ratio = s68[big] / s68[small]
        assert 0.35 < ratio < 0.7, (small, big, ratio)   # 理论 0.5，留抽样余地


def test_rho_is_always_a_valid_correlation():
    """ρ 必须落在 [-1, 1]（数值误差可能让它越界一点点，代码里已 clip）。"""
    for seed in range(1, 8):
        d = reduce_interval(_synth(20, seed=seed, noise=1.0), _mask(20))
        assert -1.0 <= d["rho"] <= 1.0, (seed, d["rho"])


def test_noiseless_input_gives_a_degenerate_rho_not_a_meaningful_one():
    """
    ⚠ 一个**必须知道**的数值陷阱：无噪声时 ρ 会算成 ±1。

    因为这时候两个 σ 都只是"浮点残差"（实测 ~1e-17，来自 `sum() − x[j]`
    的舍入），并不是 0，于是 `s68 > 0 and s76 > 0` 成立、`cov/(σ68·σ76)`
    给出 ±1。**ρ=1 是"分母趋近于 0"造成的，不是"两个比值完美相关"的证据。**
    判读真实数据的 ρ 之前，先确认 σ 不是浮点残差。
    """
    net = _synth(n=20, excess204=0.0, noise=0.0)
    d = reduce_interval(net, _mask(20))
    assert d["s68"] < 1e-15 and d["s76"] < 1e-15, "无噪声时 σ 应为浮点残差"
    assert abs(abs(d["rho"]) - 1.0) < 1e-12, d["rho"]


def test_fix_mode_forces_the_bulk_common_lead():
    """
    fix 模式（深度剖面用）：整段算好的 f206 被"借用"到每个短窗口上。

        c206 = f206_fix · i206(本窗口)   ⇒   c204 = c206 / (206/204)_common

    所以回填的 i204 必须等于 `f206_fix · mean(206) / sk[0]`。
    另外：fix 模式下**不做显著性检验**（判据是整段给的），
    所以 `i204_significant` 必须是 None，不能报 True/False 骗人。
    """
    net = _synth(n=20, seed=1, excess204=5.0, noise=1.0)
    sk = (18.137097204708997, 15.600267822093336, 38.11384303117978)
    d = reduce_interval(net, _mask(20), fix=(0.01, sk))
    assert d["i204_significant"] is None
    assert abs(d["i204"] - 0.01 * net[206].mean() / sk[0]) < 1e-9
    assert abs(d["f206"] - 0.01) < 1e-12, "f206 必须就是整段那个值"
    assert d["sk"] == sk


def test_fix_mode_with_zero_f206_does_not_correct():
    """整段没检测到普通铅时，短窗口也不扣（f206_fix = 0 走"不扣"分支）。"""
    net = _synth(n=20, seed=1, excess204=5.0, noise=1.0)
    d = reduce_interval(net, _mask(20), fix=(0.0, (0.0, 0.0, 0.0)))
    assert d["i204_significant"] is None
    assert d["i204"] == 0.0
    assert d["f206"] == 0.0
    assert d["sk"] == (0.0, 0.0, 0.0)
    # 与"完全不扣"的独立模式结果一致（除了 i204_significant 的语义不同）
    raw = reduce_interval(net, _mask(20), n_sigma=1e9)
    assert abs(d["R68"] - raw["R68"]) < 1e-15
    assert abs(d["R76"] - raw["R76"]) < 1e-15


def test_passthrough_fields():
    """质控要用的三个直通量：U_cps / n_cycles / ThU / R82。"""
    net = _synth(n=20, seed=5, noise=1.0)
    d = reduce_interval(net, _mask(20, 13))
    assert d["n_cycles"] == 13
    assert abs(d["U_cps"] - net[238][:13].mean()) < 1e-9
    assert abs(d["ThU"] - net[232][:13].mean() / net[238][:13].mean()) < 1e-12
    assert abs(d["R82"] - net[208][:13].mean() / net[232][:13].mean()) < 1e-12


def _run_standalone() -> int:
    """见 tests/_selftest.py —— 让这个文件不装 pytest 也能直接跑。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _selftest
    return _selftest.run(globals())


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
