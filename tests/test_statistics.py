"""
druid 的统计内核自检 —— 不依赖 pytest，直接跑也认：

    python tests/test_statistics.py
    python -m pytest tests/            （装了 pytest 时）

为什么这个文件必须存在
----------------------
`core.statistics.chi2_sf()` 是手写的不完全 Gamma 函数，用来给 MSWD 配一个
与 R 端 ADEPT 同口径的卡方上尾概率。这类数值代码**错了不会报错**，
只会让"这个坪自不自洽"的判断悄悄反过来，所以必须钉死。

下面 EXPECTED 里的参考值全部来自 R：

    pf(mswd, df, Inf, lower.tail = FALSE)

即 ADEPT 内部用的那一句。逐位对照过 153 组 (df, mswd) 网格，
最大绝对差 3.2e-15（机器精度）。这里固化成少量代表点。
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                              # noqa: E402

from druid.core.statistics import chi2_sf, robust_mask, weighted_mean   # noqa: E402


# (df, mswd, P(χ²_df > mswd·df))，参考值取自 R 的 pf(mswd, df, Inf, lower.tail=FALSE)
EXPECTED = [
    (1, 0.500, 0.47950012218695298),
    (2, 1.000, 0.36787944117144200),
    (3, 2.000, 0.11161022509471299),
    (5, 1.000, 0.41588018699550799),
    (8, 1.500, 0.15120388277664801),
    (12, 1.000, 0.44567964136461102),
    (20, 1.200, 0.24239216167051200),
    (40, 1.000, 0.47025726683924002),
    (99, 1.500, 0.00095164399668730),
]


def test_chi2_sf_matches_r():
    """与 R 的 pf(mswd, df, Inf, lower.tail=FALSE) 逐位一致。"""
    for df, mswd, ref in EXPECTED:
        got = chi2_sf(mswd * df, df)      # ← 乘以 df 这一步是口径的关键
        assert abs(got - ref) < 1e-14, (df, mswd, got, ref)


def test_chi2_sf_df2_has_a_closed_form():
    """
    自由度为 2 时卡方上尾概率有解析解 P = exp(−χ²/2)。
    用它当独立于参考值的第二道检查 —— 万一整条链路都被同一个错误带偏，
    这一条能抓住。
    """
    for chi2 in (0.5, 1.0, 2.5, 7.0, 20.0):
        assert abs(chi2_sf(chi2, 2) - math.exp(-chi2 / 2.0)) < 1e-15, chi2


def test_chi2_sf_boundaries():
    assert chi2_sf(0.0, 5) == 1.0                  # P(χ² > 0) = 1
    assert chi2_sf(float("inf"), 5) == 0.0         # 上尾为 0，不能退化成 nan
    assert math.isnan(chi2_sf(float("nan"), 5))
    assert math.isnan(chi2_sf(1.0, 0))             # df 必须为正
    assert math.isnan(chi2_sf(1.0, -3))
    # 概率必须落在 [0, 1]，且随 χ² 单调不增
    xs = [0.01, 0.5, 1, 5, 20, 100, 1000]
    ps = [chi2_sf(x, 7) for x in xs]
    assert all(0.0 <= p <= 1.0 for p in ps)
    assert all(a >= b for a, b in zip(ps, ps[1:])), ps


def test_chi2_sf_is_accurate_in_the_far_tail():
    """
    尾部最容易被 1−P 的相消误差毁掉。这三个值直接取自 R：
        pchisq(200, 10, lower.tail = FALSE)  = 1.6139305336977309e-37
        pchisq(50,   4, lower.tail = FALSE)  = 3.6108654048906463e-10
        pchisq(1000, 40, lower.tail = FALSE) = 1.161138236365771e-183
    相对误差都要在 1e-12 以内。（顺带在 R 里核对过
    pchisq(x, df, FALSE) 与 pf(x/df, df, Inf, FALSE) 严格相等，
    所以用 pf 口径的 ADEPT 与本函数可比。）
    """
    for chi2, df, ref in ((200.0, 10, 1.6139305336977309e-37),
                          (50.0, 4, 3.6108654048906463e-10),
                          (1000.0, 40, 1.161138236365771e-183)):
        got = chi2_sf(chi2, df)
        assert abs(got / ref - 1.0) < 1e-12, (chi2, df, got, ref)


def test_weighted_mean_unchanged():
    """加权平均的老行为不许变 —— 深度剖面域与标样 QC 都依赖它。"""
    x = np.array([10.0, 12.0, 11.0])
    s = np.array([1.0, 2.0, 1.0])
    mu, se, mswd, n = weighted_mean(x, s)
    w = 1.0 / s ** 2
    assert n == 3
    assert abs(mu - (w * x).sum() / w.sum()) < 1e-14
    assert abs(se - 1.0 / math.sqrt(w.sum())) < 1e-14
    assert abs(mswd - (w * (x - mu) ** 2).sum() / 2.0) < 1e-14

    # 等权时退化为算术平均 —— 与 ADEPT 的 "σ 为常数则加权均值 = 算术均值" 同理
    mu2, _, _, _ = weighted_mean(x, np.full(3, 0.5))
    assert abs(mu2 - x.mean()) < 1e-14

    # 非法点被剔除而不是让整段失败
    mu3, _, _, n3 = weighted_mean([1.0, np.nan, 3.0], [1.0, 1.0, 0.0])
    assert n3 == 1 and mu3 == 1.0

    # n = 0 与 n = 1
    assert weighted_mean([], [])[3] == 0
    assert math.isnan(weighted_mean([5.0], [1.0])[2])


def test_robust_mask_unchanged():
    """中位数 + MAD 判离群，以及 MAD = 0 时的兜底行为。"""
    # 20 个 1.0 加一个明显偏高的点，MAD 非零 → 偏高点被剔除
    v = np.array([1.0, 1.0, 1.05, 0.95, 1.02, 0.98, 1.01, 0.99,
                  1.03, 0.97, 1.04, 0.96, 1.06, 0.94, 1.07, 0.93,
                  1.08, 0.92, 1.09, 0.91, 50.0])
    keep, med = robust_mask(v)
    assert not keep[-1], "50.0 应被判为离群"
    assert keep.sum() == 20

    # MAD = 0（数据几乎全等）时无从判别 → 一律保留，绝不返回空掩码
    keep2, _ = robust_mask(np.ones(5))
    assert keep2.all()
    keep3, _ = robust_mask(np.array([1.0] * 20 + [50.0]))
    assert keep3.all(), "MAD=0 时按文档是一律保留（宁可漏杀，不可错杀）"


def _run_standalone() -> int:
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  ok    {name}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} 通过")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
