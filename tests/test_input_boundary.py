"""
输入边界自检 —— Qtegra CSV 解析、剥蚀区间定位、气体空白、窗口切分。

为什么用合成数据
----------------
真实批次是未发表的分析数据，不进这个仓库（`.gitignore` 里有 `结果/`，
而且那 98 个文件是从历史里剔出来的）。但这一层的契约完全可以用几条合成记录
表达清楚，跑得也快 —— 这正是它适合放进 CI 的原因。

为什么这一层最值得测
--------------------
上面的函数全是"静默出错"型：读错一列、把 washout 拖尾算进积分、
窗口边界差一个步长，都不会抛异常，只会让年龄悄悄偏一点。等发现时，
那个偏差已经在论文里了。
"""
from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np                                              # noqa: E402

import _selftest                                                # noqa: E402
from druid.depth.windows import window_edges                    # noqa: E402
from druid.io.qtegra import read_qtegra                         # noqa: E402
from druid.reduction.ablation import (blank_mean,               # noqa: E402
                                      find_ablation, subtract_blank)

# ── 一份形态与 Qtegra 输出一致的记录 ────────────────────────────────────────
# 前 4 行元数据 + 表头 + **单位行** + 3 行数据。
# 单位行（第二个字段是 "dwell time=0.01;..."）是最容易出事的地方：
# 谁要是用 pandas.read_csv 加 skiprows，早晚会因为不同批次元数据行数不同而错位。
CSV_TEXT = "\n".join([
    "MY SAMPLE:03/01/2022 06:50:39 AM;",
    "Software:Name=Qtegra;Version=2.8.3170.309;File Version=1;",
    "RF Generator:RF Plasma Lit Readback=1;Plasma Power Readback=1548.61;",
    "Pulse Counting:Threshold=2500000;",
    "",
    "Time,29Si,232Th,204Pb,206Pb,207Pb,208Pb,238U,",
    ",dwell time=0.01;xcal factor=67325.46529,dwell time=0.01,dwell time=0.01,"
    "dwell time=0.01,dwell time=0.01,dwell time=0.01,dwell time=0.01,",
    "0.01234,58436.27,0.5,300.0,1000.0,200.0,500.0,200000.0,",
    "0.33125,60043.86,0.6,320.0,1100.0,210.0,510.0,210000.0,",
    "0.65002,54016.45,0.4,290.0,900.0,190.0,490.0,190000.0,",
]) + "\n"


def _write_csv(text: str, newline: str = "\r\n") -> Path:
    """写一个临时 CSV。默认 CRLF（Qtegra 导出就是 CRLF），LF 要显式指定。"""
    tmp = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False,
                                      encoding="utf-8", newline=newline)
    tmp.write(text)
    tmp.close()
    return Path(tmp.name)


def _raises(exc, fn, *a, **kw) -> None:
    """标准库版的 assertRaises —— 不用 pytest 的 raises。"""
    try:
        fn(*a, **kw)
    except exc:
        return
    except Exception as e:                                      # noqa: BLE001
        raise AssertionError(f"应抛 {exc.__name__}，实际抛 {type(e).__name__}: {e}")
    raise AssertionError(f"应抛 {exc.__name__}，但没有抛")


# ═════════════════════════════════════════════════════════════════════════════
# 一、Qtegra CSV 解析
# ═════════════════════════════════════════════════════════════════════════════
def test_read_qtegra_title_time_and_masses():
    p = _write_csv(CSV_TEXT)
    title, t, data, columns = read_qtegra(p)

    assert title == "MY SAMPLE", title
    assert columns[0] == "Time", columns[:2]
    assert t.size == 3, f"应只有 3 行数据，实际 {t.size}"
    assert abs(t[0] - 0.01234) < 1e-12
    assert abs(t[-1] - 0.65002) < 1e-12

    # 只保留 U-Pb 体系用得到的质量数；29Si 这类微量元素不进 data
    assert set(data) == {232, 204, 206, 207, 208, 238}, sorted(data)
    assert 29 not in data
    assert abs(data[238][1] - 210000.0) < 1e-9
    assert abs(data[204][2] - 290.0) < 1e-9


def test_read_qtegra_skips_the_unit_row():
    """
    单位行含 "dwell time=0.01;xcal factor=..."，float() 解析必然失败 →
    被自然跳过。这条断言是防止有人"优化"成固定 skiprows：
    一旦不同批次元数据行数变了，就会把单位行或表头当数据读进来。
    """
    p = _write_csv(CSV_TEXT)
    _, t, _, _ = read_qtegra(p)
    assert t.size == 3           # 4 行数据形态（表头+单位行+3 行数值）里只有 3 行是数据
    assert not np.any(np.isnan(t))


def test_read_qtegra_tolerates_lf_line_endings():
    """CRLF 是常态，但别在 LF 上崩掉（有人会用编辑器另存一次）。"""
    p = _write_csv(CSV_TEXT, newline="")        # newline="" → 原样写出 LF
    title, t, data, _ = read_qtegra(p)
    assert title == "MY SAMPLE" and t.size == 3 and 238 in data


def test_read_qtegra_rejects_malformed_input():
    _raises(ValueError, read_qtegra, _write_csv(""))
    _raises(ValueError, read_qtegra,
            _write_csv("NO SAMPLE:foo;\nMetadata:only=1;\n"))
    # 有表头但一行数值都没有 → 不能返回空数组让下游去猜
    _raises(ValueError, read_qtegra,
            _write_csv("S:1;\n\nTime,238U,\n,dwell time=0.01,\n"))


# ═════════════════════════════════════════════════════════════════════════════
# 二、剥蚀区间定位
# ═════════════════════════════════════════════════════════════════════════════
def _synthetic_spot(blank=50.0, plateau=1.0e6, tail=180.0, n=180,
                    dt=0.33, seed=0):
    """
    造一条形态真实的 238U 记录：空白 → 剥蚀平台 → washout 拖尾。

    尾部取 180 cps 是刻意的：它**高于**第一级阈值（bg + k·σ ≈ 107 cps），
    但远低于第二级阈值（bg + 3% × 峰高 ≈ 3 万 cps）。只有两级阈值都用上，
    才能把拖尾砍掉 —— 这正是 find_ablation 文档里说的那件事。
    """
    t = np.arange(n) * dt
    u = np.empty(n)
    u[:60] = blank
    u[60:120] = plateau
    u[120:] = tail
    u = u + np.random.default_rng(seed).normal(0.0, 2.0, n)
    return t, u


def test_find_ablation_cuts_the_washout_tail():
    t, u = _synthetic_spot()
    t0, t1 = find_ablation(t, u)

    assert 19.0 < t0 < 21.0, f"起点应在平台前缘（约 19.8 s），实际 {t0}"
    assert 39.0 < t1 < 40.5, f"终点应在平台后缘（约 39.6 s），实际 {t1}"
    # 关键：拖尾（180 cps > thr1）必须被砍掉。若哪天有人把 thr2 删掉，
    # 这里会一直拉到记录末尾，也就是 59 s 附近。
    assert t1 < t[-1] - 10.0, f"washout 拖尾没被砍掉，t1={t1}，记录末={t[-1]}"


def test_find_ablation_handles_a_degenerate_series():
    """信号全程平坦（点没打上 / 走空文件）时不许抛异常，也不许返回 NaN。"""
    t = np.arange(60) * 0.5
    u = np.full(60, 40.0)
    t0, t1 = find_ablation(t, u)
    assert math.isfinite(t0) and math.isfinite(t1)
    assert t0 <= t1


# ═════════════════════════════════════════════════════════════════════════════
# 三、气体空白
# ═════════════════════════════════════════════════════════════════════════════
def test_blank_mean_takes_the_segment_just_before_ablation():
    """
    气体空白的主要成分是载气里的 Hg，会缓慢漂移，所以必须取**紧挨剥蚀之前**
    的那一段，而不是文件开头。这里让文件开头本底 500、剥蚀前降到 60：
    正确实现给 60，误改成"取开头"就会给 500。
    """
    t = np.arange(0.0, 40.0, 0.5)
    v = np.full(t.size, 500.0)
    v[t >= 18.0] = 60.0
    b = blank_mean(t, v, t_ab0=35.0, dur=15.0)      # 窗口 = [19.5, 34.5] s
    assert abs(b - 60.0) < 1e-9, b


def test_blank_mean_gives_up_rather_than_guessing():
    """可用点太少时返回 0（宁可不扣空白，也不用一个离谱的值）。"""
    t = np.array([0.0, 100.0])
    v = np.array([500.0, 900.0])
    assert blank_mean(t, v, t_ab0=0.2, dur=15.0) == 0.0


def test_subtract_blank_keeps_negative_values():
    """
    净信号允许为负 —— 这是统计学事实（信号等于本底时相减随机波动）。
    一旦在这里 clip 到 0，就会引入系统性正偏差。这条断言把它钉住。
    """
    t = np.arange(0.0, 40.0, 0.5)
    u = np.full(t.size, 100.0)
    u[t >= 30.0] = 20.0                      # 剥蚀期间反而低于本底
    net, blank = subtract_blank(t, {238: u}, t_ab0=30.0)

    assert abs(blank[238] - 100.0) < 1e-9, blank[238]
    assert net[238][t >= 30.0].max() < 0.0, "负值被 clip 掉了"


# ═════════════════════════════════════════════════════════════════════════════
# 四、窗口切分
# ═════════════════════════════════════════════════════════════════════════════
def test_window_edges_overlap_by_win_minus_step():
    e = window_edges(10.0, 20.0, win=4.0, step=1.0)
    assert e.shape == (7, 2), e.shape
    assert np.allclose(e[0], [10.0, 14.0])
    assert np.allclose(e[-1], [16.0, 20.0])
    assert np.allclose(e[:, 1] - e[:, 0], 4.0)                 # 每窗宽 4 s
    assert np.allclose(e[1:, 0] - e[:-1, 0], 1.0)              # 步长 1 s
    # 相邻窗口共享 3 s 数据 → 它们不独立。这就是 MSWD 判据里自由度要打折的根源，
    # domains.py 用 k_eff ≈ k/2.5 近似，别把它当成独立点数。
    assert np.allclose(np.diff(e[:, 0]), 1.0)


def test_window_edges_span_shorter_than_one_window():
    """整段比一个窗口还短时退化为单窗口，而不是返回空集（图上至少要有一个点）。"""
    e = window_edges(10.0, 12.0, win=4.0, step=1.0)
    assert e.shape == (1, 2), e.shape
    assert np.allclose(e[0], [10.0, 12.0])


def test_window_edges_never_overrun_the_ablation_span():
    """
    窗口不许越过剥蚀段的终点（越过去就把 washout 混进来了）。

    顺带钉住一个**容易误判成 bug** 的事实：固定步长的滑动窗口留不下一个
    完整的窗口时，末尾那段就空着 —— 最后一个窗口的终点可能距 t1 还有
    最多一个步长（win=4 / step=1 时就是不到 1 s，约 35 s 剖面的 2%）。
    这是标准做法（否则最后会多出一个宽度不等的窗口，τ 间距不再均匀），
    所以这里断言的是"残差 < 步长"，不是"正好贴到 t1"。
    """
    for t1 in (20.0, 20.5, 33.3, 47.77):
        e = window_edges(10.0, t1, win=4.0, step=1.0)
        assert e[-1, 1] <= t1 + 1e-9, (t1, e[-1])
        assert t1 - e[-1, 1] < 1.0, (t1, e[-1])


if __name__ == "__main__":
    raise SystemExit(_selftest.run(globals()))
