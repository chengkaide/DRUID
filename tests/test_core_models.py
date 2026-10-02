"""
druid 的两个"物理模型"自检 —— 不依赖 pytest，直接跑也认：

    python tests/test_core_models.py
    python -m pytest tests/            （装了 pytest 时）

被这两条链路坑过的方式
----------------------
**Stacey–Kramers 普通铅模型**（`core/common_lead.py`）：
它给出的是"给定年龄 t 时，普通铅的 206/204、207/204、208/204 应该是多少"。
错法很隐蔽 —— 三个系数的初始值（9.307 / 10.294 / 29.476）与两个阶段的 μ
写错一个，结果仍然"看起来像一组合理的 Pb 同位素比值"，只有拿已知样品
（SK 模型自己公布的表）去对才会发现。这里钉住的是**结构性质**：
三个比值随年龄单调、端点钳在 [0, 4570] Ma、数组输入不退化。

**死时间**（`core/deadtime.py`）：
公式只有一行 `I_true = I_obs / (1 − I_obs·τ)`，但它有一个必须存在的
**夹具**：分母在 I ≥ 1/τ 处穿过 0，若不钳位就会给出天文数字或负数。
这里钉住公式、钳位上限（2 倍）、以及 τ ≤ 0 时"原样返回且不返回引用"。
最后一条最容易被忽略：返回引用的话，调用方就地改动会**静默污染**输入。
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                              # noqa: E402

from druid.core.common_lead import stacey_kramers               # noqa: E402
from druid.core.deadtime import correct_channels, correct_rate   # noqa: E402

# 冻结锚点：模型在几个年龄上的取值。改了系数或分阶段方式，这里当场红。
SK_AT = {
    0.0: (18.702613, 15.630722, 38.821480),
    1000.0: (17.068203, 15.512174, 36.845455),
    3000.0: (12.930579, 14.344963, 32.587654),
}
T0 = 4570.0                                     # 模型的时间上界（地球年龄）


def test_stacey_kramers_anchors():
    """三个比值在代表年龄上的字面值 —— 防"悄悄改了模型系数"。"""
    for t, (a6, a7, a8) in SK_AT.items():
        p6, p7, p8 = stacey_kramers(t)
        assert abs(float(p6) - a6) < 5e-7, (t, p6, a6)
        assert abs(float(p7) - a7) < 5e-7, (t, p7, a7)
        assert abs(float(p8) - a8) < 5e-7, (t, p8, a8)


def test_stacey_kramers_is_monotonic_in_age():
    """
    年龄越大 → 圈闭时铅演化得越少 → 三个比值都必须**单调下降**。

    这条比数值锚点更耐改：将来若换成别的普通铅模型（只要仍是同一族演化模型），
    锚点会红但单调性仍应成立。
    """
    ts = [0.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0, T0]
    p6 = [float(stacey_kramers(t)[0]) for t in ts]
    p7 = [float(stacey_kramers(t)[1]) for t in ts]
    p8 = [float(stacey_kramers(t)[2]) for t in ts]
    for name, seq in (("206/204", p6), ("207/204", p7), ("208/204", p8)):
        assert all(b < a for a, b in zip(seq, seq[1:])), (name, seq)
    # 三个比值的量级必须依次递增（208 由 Th 贡献，最大）
    assert p6[0] < p8[0] and p7[0] < p8[0]


def test_stacey_kramers_clamps_outside_the_model_domain():
    """
    模型只定义在 [0, 4570] Ma 上。超出范围**钳位而不是外推** ——
    外推会让 exp() 溢出、或者给出"比原始铅还原始"的荒谬比值。
    钳位后超界输入必须与端点**完全相等**（不是近似）。
    """
    lo = tuple(float(x) for x in stacey_kramers(0.0))
    hi = tuple(float(x) for x in stacey_kramers(T0))
    assert tuple(float(x) for x in stacey_kramers(-50.0)) == lo
    assert tuple(float(x) for x in stacey_kramers(9000.0)) == hi
    assert tuple(float(x) for x in stacey_kramers(T0 + 1e9)) == hi


def test_stacey_kramers_is_vectorised():
    """数组进 → 三个同形状数组出；标量进 → 标量。"""
    arr = stacey_kramers(np.array([100.0, 1000.0, 3000.0]))
    assert len(arr) == 3
    assert np.shape(arr[0]) == (3,) and np.shape(arr[1]) == (3,)
    for i, t in enumerate((100.0, 1000.0, 3000.0)):
        assert abs(float(arr[0][i]) - float(stacey_kramers(t)[0])) < 1e-12
    assert np.ndim(stacey_kramers(500.0)[0]) == 0


def test_deadtime_formula():
    """
    公式直查：I_true = I_obs / (1 − I_obs·τ)，τ 以秒参与运算（入参是纳秒）。

    顺带把文档里那句量级感钉住：14.7645 ns 下 1e6 cps 的修正 = +1.50%。
    这个数是"样品锆石工作点"的修正幅度，改公式会让它偏移。
    """
    tau_ns = 14.7645
    obs = 1.0e6
    want = obs / (1.0 - obs * tau_ns * 1e-9)
    got = float(correct_rate(np.array([obs]), tau_ns)[0])
    assert abs(got - want) < 1e-9, (got, want)
    assert abs((got / obs - 1.0) - 0.0150) < 1e-3


def test_deadtime_is_identity_when_disabled():
    """
    τ ≤ 0（含 None）＝ 不校正，原样返回 —— 而且**必须是拷贝**。

    返回引用是最阴的一类 bug：调用方就地改一下"校正后"的数组，
    原数组跟着变，而调用方以为自己动的是副本。
    """
    for tau in (None, 0.0, -1.0):
        r0 = np.array([1.0, 2.0, 3.0])
        out = correct_rate(r0, tau)
        assert np.array_equal(out, r0)
        out[0] = 99.0
        assert r0[0] == 1.0, f"tau={tau} 时返回了引用，调用方改动污染了输入"
    # 通道批量版本同理：τ ≤ 0 走浅拷贝分支，也不该改到原字典
    raw = {238: np.array([1.0, 2.0])}
    cc = correct_channels(raw, 0.0)
    assert cc[238] is raw[238]                  # 浅拷贝：值对象共享，但字典是新的
    assert set(cc) == set(raw)


def test_deadtime_clamps_the_correction_at_two_times():
    """
    分母 `1 − I·τ` 在 I ≥ 1/τ 处穿过 0 —— 公式在那里失去物理意义。
    夹具把它钳在 0.5，即**最多允许 2 倍修正**。

    没有夹具会怎样：1e12 cps、τ = 14.7645 ns 时 `I·τ ≈ 1.5e4`，
    分母是个大的**负数**，结果会是 **−6.8e7**（比输入还小、还是负的）。
    这种数一旦进了比值，会发现成 "负年龄" 而不是报错。有夹具则恒 ≤ 2×。
    """
    huge = np.array([1e12, 1e15, 1e30])
    out = correct_rate(huge, 14.7645)
    assert np.allclose(out / huge, 2.0), out / huge
    assert np.isfinite(out).all()
    assert (out >= huge).all()                  # 校正只会抬高，不会压低


def test_deadtime_corrects_every_channel():
    """批量版逐通道施加同一个 τ，且不改原字典的内容。"""
    raw = {238: np.array([1e6]), 207: np.array([5e3]), 204: np.array([20.0])}
    out = correct_channels(raw, 14.7645)
    assert set(out) == set(raw)
    for k in raw:
        assert np.allclose(out[k], correct_rate(raw[k], 14.7645))
    # 高强度通道被抬得更多 —— 这正是它扭曲比值、必须逐通道做校正的原因
    up238 = out[238][0] / raw[238][0] - 1.0
    up204 = out[204][0] / raw[204][0] - 1.0
    assert up238 > up204 > 0.0


def _run_standalone() -> int:
    """见 tests/_selftest.py —— 让这个文件不装 pytest 也能直接跑。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _selftest
    return _selftest.run(globals())


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
