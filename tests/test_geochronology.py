"""
druid 的年龄 ⇄ 比值换算自检 —— 不依赖 pytest，直接跑也认：

    python tests/test_geochronology.py
    python -m pytest tests/            （装了 pytest 时）

为什么这个文件必须存在
----------------------
`core/geochronology.py` 是全包的**度量衡**：所有年龄、所有比值、所有不确定度
最后都要落到这四个函数上。它错的方式非常"安静"：

  · 把 λ 的指数写成 `exp(λt)` 而不是 `expm1(λt)` —— 大年龄上差到 1% 看不出来；
  · `r75_from` 里 `238U/235U` 乘成除 —— 207/235 直接差 138 倍，但只在
    老锆石上才看得出来，年轻样品看起来"还行"；
  · `age68` 忘了 clip —— 空白扣除后净信号为负的窗口变成 NaN，而 NaN 会在
    加权平均里**传染整段**，症状是"这个点莫名其妙没有年龄"。

上面任何一条都不会抛异常。所以这里用三条互相独立的判据把它钉住：
① **往返**（正演 → 反演，两条独立实现）；② **交叉恒等式**（`r75_from`
   必须复现 `expm1(λ235·t)`）；③ **边界行为**（钳位、饱和、标量/数组形状）。
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                              # noqa: E402

from druid.core.constants import L235, L238, U238_U235          # noqa: E402
from druid.core.geochronology import (age68, age75, age76, r68_of_age,   # noqa: E402
                                      r75_from, r76_of_age)

# 覆盖年轻锆石（100 Ma）、本批主战场（460 Ma）、老核（1~3 Ga）
AGES = [1.0, 100.0, 460.0, 1000.0, 2000.0, 3000.0]

# 回归锚点：正演比值的字面值。这不是"照抄代码"，而是把**当前口径**固化下来 ——
# 换了衰变常数或改了 expm1 写法，这里会当场红，逼你确认是不是有意的。
R68_1000MA = 0.16780392747297115
R76_1000MA = 0.0725322627457843


def test_r68_age68_round_trip():
    """正演 r68_of_age → 反演 age68，必须回到原值。"""
    for t in AGES:
        back = float(age68(r68_of_age(t)))
        assert abs(back - t) <= 1e-9 * max(t, 1.0), (t, back)


def test_r76_age76_round_trip():
    """
    r76_of_age → age76（二分法）同样往返。

    容差给 1e-6 relative：二分法在 [1, 5000] 上跑 80 次，区间已收到 5e-21 Ma，
    所以误差其实完全由正演那一步的浮点决定，不该有可见偏差。
    """
    for t in AGES:
        back = age76(r76_of_age(t))
        assert abs(back - t) <= 1e-6 * max(t, 1.0), (t, back)


def test_forward_ratios_match_frozen_anchors():
    """正演比值的口径锚点（防"悄悄换了衰变常数"）。"""
    assert abs(float(r68_of_age(1000.0)) / R68_1000MA - 1.0) < 1e-12
    assert abs(float(r76_of_age(1000.0)) / R76_1000MA - 1.0) < 1e-12


def test_r75_from_reproduces_the_direct_ratio():
    """
    交叉恒等式 —— 本文件里最硬的一条。

        (207/235) = (207/206) · (206/238) · (238/235)

    左边由 `expm1(λ235·t)` 直接算，右边由 `r75_from` 从另外两个比值拼出来。
    两条路只共用 `238U/235U` 一个常数，所以能同时抓住"乘除搞反"（差 138 倍）
    和"漏乘丰度比"（差 1.9 万倍）这两类错法。
    """
    for t in AGES:
        got = float(r75_from(r68_of_age(t), r76_of_age(t)))
        want = math.expm1(L235 * t * 1e6)
        assert abs(got / want - 1.0) < 1e-14, (t, got, want)


def test_age75_is_the_sister_of_age68():
    """age75 只是把 λ238 换成 λ235 —— 用同一条 460 Ma 的比值钉住它俩别写串。"""
    t = 460.0
    r68 = float(r68_of_age(t))
    r75 = float(r75_from(r68, float(r76_of_age(t))))
    assert abs(float(age68(r68)) - t) < 1e-9
    assert abs(float(age75(r75)) - t) < 1e-9
    # 同一个年龄下 207/235 必须比 206/238 大（235U 衰变快得多）
    assert r75 > r68 * 5.0


def test_age76_is_monotonic_and_bounded():
    """
    二分法的收敛前提是 f(t) 单调 —— 这里从外部确认一遍，
    并钉住它的**饱和行为**：区间是 [1, 5000] Ma，落在区间外的比值会被钳在端点上。

    ⚠ 取样点必须**高于 r76 的下限** 0.0460662（见下一条测试）：
    t→0 时 207Pb/206Pb 趋于 λ235/(λ238·238U/235U) 而不是 0，
    所以 0.01 / 0.03 这种值全都被钳在 1.0 Ma，拿它们测单调性只会测到钳位。
    """
    rs = [0.047, 0.05, 0.06, 0.08, 0.15, 0.3, 0.6]
    ages = [age76(r) for r in rs]
    assert all(b > a for a, b in zip(ages, ages[1:])), ages
    assert all(1.0 <= a <= 5000.0 for a in ages), ages

    # 比值小到 1 Ma 以下 → 钳在下界；比值超出 5000 Ma 的正演值 → 钳在上界
    assert age76(1e-9) == 1.0
    assert age76(-0.3) == 1.0, "负比值没有物理意义，应落在下界而不是抛异常"
    assert age76(1.0) == 5000.0


def test_r76_has_a_floor_which_is_not_zero():
    """
    `207Pb/206Pb` 的**下限不是 0**，而是 t→0 时的极限：

        lim(t→0) (207/206) = λ235 / (λ238 · 238U/235U) = 0.0460662

    因为两个 Pb 同位素都是按 (exp(λt)−1)/U 累积的，比值里 λ 的比留了下来，
    t 约掉了。这是**物理上的**，不是数值缺陷。做窗口级质量筛查时若拿
    "R76 ≈ 0"当异常判据，就会把一堆正常的老年龄窗口判成异常。
    """
    floor = L235 / (L238 * U238_U235)
    assert abs(floor - 0.04606619605024173) < 1e-16
    # 钳位（下界 1 Ma）只在比值落在下限以下时发生
    assert age76(floor * 0.99) == 1.0
    assert age76(floor * 1.01) > 1.0


def test_age68_clips_instead_of_producing_nan():
    """
    `age68` 的核心承诺：净信号为负（比值 < 0）时**不许产生 NaN**。

    为什么这么要紧：NaN 一旦进入 `weighted_mean` / `np.nanmean` 之外的任何
    聚合，整段就废了，而症状只是"某个点没有年龄" —— 极难追。
    所以这里直接断言"finite"，而不是断言某个具体值。
    """
    for r in (-0.5, -0.999999, -2.0, -100.0):
        v = float(age68(r))
        assert math.isfinite(v), (r, v)
        assert v < 0.0, (r, v)          # 退化成负年龄，不是 0，也不是 nan


def test_age68_is_vectorised_and_keeps_scalar_shape():
    """标量进标量出、数组进数组出 —— 调用方（深度剖面）依赖这个约定。"""
    arr = age68(np.array([0.05, 0.06, 0.07]))
    assert isinstance(arr, np.ndarray) and arr.shape == (3,)
    assert abs(float(arr[1]) - float(age68(0.06))) < 1e-12
    assert isinstance(float(age68(0.06)), float)
    assert np.isfinite(age68(np.array([-0.1, 0.06]))).all()


def test_r68_of_age_is_linear_at_the_origin():
    """
    `expm1` 的用处在这里体现：t→0 时比值必须趋于 0，且**偏离线性的量恰好是
    泰勒展开的第二项** ——

        expm1(x)/x = 1 + x/2 + x²/6 + …      （x = λ238·t）

    只断言"约等于线性"是没意义的（任何实现都约等于）。断言残差**精确等于 x/2**
    才能同时抓住两类错法：写成 `exp(x) − 1`（大 x 下相消）、或者误用 `log1p`。
    """
    assert float(r68_of_age(0.0)) == 0.0
    tiny = 1e-3                                    # 1 ka，比值约 1.55e-7
    x = L238 * tiny * 1e6
    rel = float(r68_of_age(tiny)) / x - 1.0
    assert abs(rel - x / 2.0) < 1e-13, (rel, x / 2.0)


def test_r76_of_age_at_zero_is_nan_by_construction():
    """
    `r76_of_age(0)` 是 0/0 —— 这是**构造性的**，不是 bug：
    207Pb/206Pb 在零年龄上没有定义（两条链都没累积出铅）。
    钉住它是为了说明"看到 nan 别慌"，同时保证它是 nan 而不是 0
    （返回 0 会被下游当成一个合法的、极年轻的年龄）。
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        v = r76_of_age(0.0)
    assert math.isnan(float(v))


def test_u238_u235_is_the_only_abundance_constant_in_play():
    """
    `r75_from` 与 `r76_of_age` 允许传入自定义丰度比 —— 两条路必须同步跟着变。
    这条抓的是"某一处把参数忘了传下去"。
    """
    u = 137.0                                       # 故意换一个值
    t = 1000.0
    got = float(r75_from(r68_of_age(t), r76_of_age(t, u238_u235=u), u238_u235=u))
    want = math.expm1(L235 * t * 1e6)
    assert abs(got / want - 1.0) < 1e-13, (got, want)
    assert U238_U235 > 100.0                        # 默认值仍在合理范围


def _run_standalone() -> int:
    """见 tests/_selftest.py —— 让这个文件不装 pytest 也能直接跑。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _selftest
    return _selftest.run(globals())


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
