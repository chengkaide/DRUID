"""
druid.qc 自检 —— 不依赖 pytest，直接跑也认：

    python tests/test_qc.py
    python -m pytest tests/test_qc.py

为什么这个文件必须存在
----------------------
`druid/qc.py` 里每一条检查都对应一个**对外结论**：这批数据能不能定年、
误差棒能不能用来论证一致性、交给 ADEPT 会掉多少点。这些结论一旦悄悄反转
（比如把"σext 是假设值"判成 pass），下游会拿着假前提去解释真实年龄，
而且**任何数值回归测试都抓不到** —— 数值一个都没变。

所以这里逐条钉住"输入长这样 → 必须给这个 level"。用的是
`tests/_fixtures.py` 手搓的假 BatchResult（只模拟形状，不模拟数值）。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _fixtures import keyed, make_result                       # noqa: E402
from druid.qc import (                                          # noqa: E402
    FAIL, INFO, LEVEL_ORDER, PASS, WARN,
    Check, QCThresholds, assess_batch, counts, jsonable, summary_lines,
    verdict, worst, _NO_LINE_START, _display_width, _wrap_cjk,
)


# ═════════════════════════════════════════════════════════════════════════════
# 序列化：两个真会踩到的坑
# ═════════════════════════════════════════════════════════════════════════════
def test_jsonable_converts_numpy_and_nan():
    """
    numpy 标量不可被 json 序列化（`np.int64` 不是 `int`、`np.bool_` 不是 `bool`），
    而 NaN/Inf 不是合法 JSON 字面量。两个错都**只在真有数据的路径上出现**，
    拿小样本手测很容易漏。
    """
    out = jsonable({
        "i": np.int64(7), "b": np.bool_(True), "f": np.float64(1.5),
        "nan": float("nan"), "inf": float("inf"), "ninf": float("-inf"),
        "p": Path("a/b"), "lst": [np.int64(1), float("nan")],
        "nested": {"x": np.int32(3)},
    })
    assert out["i"] == 7 and isinstance(out["i"], int) and not isinstance(out["i"], bool)
    assert out["b"] is True
    assert isinstance(out["f"], float) and abs(out["f"] - 1.5) < 1e-12
    assert out["nan"] is None and out["inf"] is None and out["ninf"] is None
    # ⚠ 不要写成 == "a/b"：Windows 上 Path 的字符串是 `a\b`。这条测试要在
    #    ubuntu 与 windows 两个 CI job 上都过，所以平台相关的东西一律现算。
    assert out["p"] == str(Path("a/b"))
    assert out["lst"] == [1, None]
    assert out["nested"]["x"] == 3


def test_jsonable_output_survives_strict_json():
    """
    转换后的东西必须能被 `allow_nan=False` 的 json.dumps 接受 ——
    这正是写交接文件时用的模式。Python 自己能读回 `NaN`，
    别的语言的解析器一律报错，所以严格模式是必须的。
    """
    payload = jsonable({"a": np.float64("nan"), "b": np.int64(2),
                        "c": [np.bool_(False)], "d": {"e": float("inf")}})
    text = json.dumps(payload, allow_nan=False)
    back = json.loads(text)
    assert back == {"a": None, "b": 2, "c": [False], "d": {"e": None}}


# ═════════════════════════════════════════════════════════════════════════════
# 等级与聚合
# ═════════════════════════════════════════════════════════════════════════════
def test_level_ordering_and_worst():
    assert LEVEL_ORDER[FAIL] < LEVEL_ORDER[WARN] < LEVEL_ORDER[INFO] < LEVEL_ORDER[PASS]
    mk = lambda lv: Check(key="k", level=lv, title="t")     # noqa: E731
    assert worst([mk(PASS), mk(WARN), mk(PASS)]) == WARN
    assert worst([mk(PASS), mk(FAIL), mk(WARN)]) == FAIL
    assert worst([]) == PASS, "空列表必须当作 pass，不能抛"


def test_counts_always_has_all_four_levels():
    """四种等级都要有键，下游才不用 .get —— 少一个键会让下游悄悄走错分支。"""
    c = counts([Check(key="k", level=WARN, title="t")])
    assert set(c) == {FAIL, WARN, INFO, PASS}
    assert c[WARN] == 1 and c[FAIL] == 0


# ═════════════════════════════════════════════════════════════════════════════
# 决定性结论
# ═════════════════════════════════════════════════════════════════════════════
def test_uncalibrated_batch_is_fail():
    """没有主标 → 年龄不能用于定年。这条必须是 fail，不能降级成 warn。"""
    r = make_result(mode="未校准", standards=[], calibrate_column=False)
    k = keyed(assess_batch(r))
    assert k["calibration.mode"].level == FAIL
    assert k["calibration.primary_spots"].level == FAIL
    assert k["calibration.primary_spots"].data["n_primary"] == 0


def test_calibrated_batch_with_primary_passes():
    r = make_result()
    k = keyed(assess_batch(r))
    assert k["calibration.mode"].level == PASS
    assert k["calibration.primary_spots"].level == PASS


def test_assumed_sigma_ext_is_warned():
    """
    没有监控标样 → σext 退回写死的经验值。这时样品年龄的误差棒里有一块是假设，
    必须报出来；报成 pass 会让下游拿误差棒去论证"两个年龄一致"。
    """
    r = make_result(standards=[("91500", "主标", 5, 1062.4, 1059.75, 0.51)],
                    sd68=0.007, sd76=0.0025)
    k = keyed(assess_batch(r))
    chk = k["uncertainty.sigma_ext_source"]
    assert chk.level == WARN
    assert chk.data["source"] == "assumed"
    assert chk.data["n_secondary"] == 0


def test_measured_sigma_ext_is_pass():
    r = make_result()                     # 默认带 Ple 3 个点
    k = keyed(assess_batch(r))
    chk = k["uncertainty.sigma_ext_source"]
    assert chk.level == PASS and chk.data["source"] == "measured"


def test_sigma_ext_bound_is_flagged():
    """
    `external_scatter()` 在监控标样不足 2 个时返回**下限**、散度算爆时返回**上限**。
    取到边界说明这个数不是测出来的，只是缺口的形状。
    """
    from druid.core.statistics import EXTERNAL_SCATTER_HI, EXTERNAL_SCATTER_LO

    r = make_result(sd76=EXTERNAL_SCATTER_HI)          # 207/206 撞上限
    k = keyed(assess_batch(r))
    chk = k["uncertainty.sigma_ext_bound"]
    assert chk.level == WARN
    assert any("上限" in b for b in chk.data["bound"])

    r2 = make_result(sd68=EXTERNAL_SCATTER_LO)         # 206/238 撞下限
    k2 = keyed(assess_batch(r2))
    assert k2["uncertainty.sigma_ext_bound"].level == WARN

    r3 = make_result(sd68=0.011, sd76=0.036)           # 两个都不在边界
    assert keyed(assess_batch(r3))["uncertainty.sigma_ext_bound"].level == PASS


# ═════════════════════════════════════════════════════════════════════════════
# 参考值自洽性：把上游未决的口径变成一条带数值的检查
# ═════════════════════════════════════════════════════════════════════════════
def test_reference_self_consistency_flags_the_91500_anchor():
    """
    91500 同时给了文献比值 R68=0.17928 与公认年龄 1062.4 Ma，两者不自洽：
    R68 反算回年龄是 1063.04 Ma，差 **+0.0602%**。后果是全部样品年龄按同一个
    因子整体偏老，且该标样 QC「偏差%」带一个固定地板。

    这条检查就是把那个数字钉在这里 —— 以后谁把 R68 改成 0.17917（或把
    age_Ma 改成 1063.04），这条会当场变 pass，提醒"全批年龄会整体平移"。
    """
    r = make_result()
    k = keyed(assess_batch(r))
    chk = k["reference.self_consistency"]
    assert chk.level == WARN
    assert "91500" in chk.data["standards"]
    d = chk.data["standards"]["91500"]
    assert abs(d["R68"] - 0.17928) < 1e-9
    assert abs(d["age_Ma"] - 1062.4) < 1e-9
    assert 0.058 < d["diff_pct"] < 0.062, d["diff_pct"]

    # 而本身自洽的标样（Plesovice 的 R68 是由年龄反算来的）不该被误报
    assert abs(d["diff_pct"]) > 0.02, "91500 就应当超出阈值"


# ═════════════════════════════════════════════════════════════════════════════
# ADEPT 交接：本模块存在的主要理由
# ═════════════════════════════════════════════════════════════════════════════
def test_adept_dropout_classifies_old_core_and_short_ablation():
    """
    两类会被 ADEPT 静默丢掉的测点，必须分得开 —— 因为**修法完全不同**：
    老核要靠补 Age76/Age76_1s 两列救回，短采只能重采。
    """
    r = make_result(unknown=[
        dict(sample="GOOD", age=450.0, windows=[450.0] * 30),        # 正常
        dict(sample="OLD", age=1100.0, windows=[1100.0] * 30),       # 整点老核
        dict(sample="SHORT", age=455.0, windows=[455.0] * 8),        # 采集窗口不足
    ])
    chk = keyed(assess_batch(r))["handoff.adept_dropout"]
    d = chk.data
    assert d["n_dropped"] == 2
    got = {x["sample"]: x["reason"] for x in d["dropped"]}
    assert got == {"OLD": "old_core", "SHORT": "short_ablation"}
    assert d["n_old_core"] == 1 and d["n_short_ablation"] == 1
    # 掉了 2/3 = 66.7% > 10% -> 必须升成 fail，不能停在 warn
    assert d["frac"] > 0.10
    assert chk.level == FAIL


def test_adept_dropout_high_fraction_is_fail():
    """掉点占比 >10% 时不能只 warn —— 大批测点连管线都没进，结论口径会变。"""
    r = make_result(unknown=[dict(sample=f"S{i}", age=1100.0, windows=[1100.0] * 30)
                             for i in range(10)])
    chk = keyed(assess_batch(r))["handoff.adept_dropout"]
    assert chk.level == FAIL
    assert chk.data["n_dropped"] == 10


def test_adept_dropout_reads_windows_not_the_results_table():
    """
    判据必须来自 `剖面窗口` 表（ADEPT 真正拿到的东西），不能来自结果表。
    构造一个"结果表里年龄很正常、但窗口全是老核"的测点：
    只有读 windows 才判得出来。
    """
    r = make_result(unknown=[dict(sample="SNEAKY", age=450.0, windows=[1050.0] * 30)])
    chk = keyed(assess_batch(r))["handoff.adept_dropout"]
    assert chk.data["n_dropped"] == 1
    assert chk.data["dropped"][0]["n_windows"] == 30
    assert chk.data["dropped"][0]["n_usable"] == 0


def test_clean_batch_has_no_adept_dropout():
    r = make_result(unknown=[dict(sample=f"S{i}", age=450.0, windows=[450.0] * 30)
                             for i in range(10)])
    chk = keyed(assess_batch(r))["handoff.adept_dropout"]
    assert chk.level == PASS and chk.data["n_dropped"] == 0


def test_missing_windows_table_is_reported_not_crashed():
    """
    没做深度分析（--no-depth）时不得抛异常。而且**三个 key 都要在** ——
    下游按 key 取值时不该再处理"这个 key 不存在"，那会把"缺表"和"写错 key"
    混成同一类错误。判不了就报 warn，绝不报 pass（会把"没查"说成"没问题"）。
    """
    r = make_result()
    r.windows = r.windows.iloc[0:0]
    k = keyed(assess_batch(r))
    for key in ("handoff.windows", "handoff.adept_dropout", "handoff.old_core_windows"):
        assert key in k, f"缺表时也必须保留 {key} 这个 key"
        assert k[key].level == WARN, key
    assert k["handoff.adept_dropout"].data["n_dropped"] is None, "判不了就不能报 0，0 会被读成'没问题'"


# ═════════════════════════════════════════════════════════════════════════════
# 数据卫生
# ═════════════════════════════════════════════════════════════════════════════
def test_name_hygiene_catches_the_double_separator_that_splits_a_rock():
    """
    真实批次里出现过：同一批 8 个测点写成 `X--2-9`…`-16`，其余写成 `X-2-17`…。

    这些名字**不会和别人撞名**（末尾编号本来就不同），所以"查重名"抓不到它们。
    真正的后果在**归组层级**：去掉末尾测点编号后 `X-2-9` 属于 `X-2`、
    `X--2-9` 属于 `X` —— 同一个岩样被算成两个。所以这条检查必须同时报
    "非规范写法"和"归一化前后的岩样个数差"。
    """
    r = make_result(unknown=[dict(sample="XL-46-2-1", age=450.0),
                             dict(sample="XL-46-2-2", age=451.0),
                             dict(sample="XL-46--2-9", age=452.0)])
    chk = keyed(assess_batch(r))["data.name_hygiene"]
    assert chk.level == WARN
    assert chk.data["n_variant"] == 1
    assert chk.data["variant"] == ["XL-46--2-9"]
    # 归一化前有 {XL-46-2, XL-46} 两个岩样，归一化后只剩 {XL-46-2}
    assert chk.data["n_rocks_raw"] == 2
    assert chk.data["n_rocks_canonical"] == 1
    assert chk.data["n_group_split"] == 1


def test_name_hygiene_not_fooled_by_merely_similar_names():
    """真正不同的样品不能被当成"非规范写法" —— 归组过头比不归组更糟。"""
    r = make_result(unknown=[dict(sample="XL-46-2-1", age=450.0),
                             dict(sample="XL-46-2-9", age=451.0)])
    chk = keyed(assess_batch(r))["data.name_hygiene"]
    assert chk.level == PASS
    assert chk.data["n_variant"] == 0 and chk.data["n_group_split"] == 0


def test_clean_names_pass():
    r = make_result(unknown=[dict(sample="S01", age=450.0), dict(sample="S02", age=451.0)])
    chk = keyed(assess_batch(r))["data.name_hygiene"]
    assert chk.level == PASS
    assert chk.data["n_group_split"] == 0


# ═════════════════════════════════════════════════════════════════════════════
# 契约层面的不变量
# ═════════════════════════════════════════════════════════════════════════════
def test_check_keys_are_ascii_unique_and_stable():
    """
    下游按 key 分支，所以 key 必须：① 全 ASCII（下游可能是别的语言写的）；
    ② 不重复；③ 是下面这份清单。**改这份清单等于改对外契约**，
    改的时候必须同时升 druid.handoff 的 schema 版本。
    """
    k = keyed(assess_batch(make_result()))
    for key in k:
        assert key.isascii(), f"key 里出现非 ASCII：{key}"
    assert len(k) == len(assess_batch(make_result())), "key 有重复"

    expected = {
        "calibration.mode", "calibration.primary_spots", "calibration.secondary_spots",
        "calibration.primary_rejected", "samples.present",
        "standards.primary_bias", "standards.primary_mswd",
        "standards.secondary_bias", "standards.secondary_mswd",
        "reference.self_consistency",
        "uncertainty.sigma_ext_source", "uncertainty.sigma_ext_bound",
        "uncertainty.count_rate_match", "uncertainty.secondary_correction",
        "samples.concordance", "samples.common_lead", "data.name_hygiene",
        "handoff.windows", "handoff.adept_dropout", "handoff.old_core_windows",
        "samples.multi_domain", "data.skipped_files", "samples.method_difference",
    }
    missing = expected - set(k)
    extra = set(k) - expected
    assert not missing, f"少了检查项（下游在等它）：{sorted(missing)}"
    assert not extra, f"多了未登记的检查项（请一并更新本清单与 schema 版本）：{sorted(extra)}"


def test_checks_are_sorted_by_severity():
    """顺序必须固定（先 level 后 key），否则同一份输入两次导出会 diff 出一堆噪声。"""
    checks = assess_batch(make_result(mode="未校准",
                                      standards=[("91500", "主标", 5, 1062.4, 1100.0, 3.0)]))
    levels = [LEVEL_ORDER[c.level] for c in checks]
    assert levels == sorted(levels), levels
    keys = [c.key for c in checks]
    for i in range(1, len(checks)):
        if checks[i].level == checks[i - 1].level:
            assert keys[i - 1] < keys[i], (keys[i - 1], keys[i])


def test_assess_batch_works_without_a_config_object():
    """
    下游（agent / 脚本）很可能只拿到 BatchResult 而没有 BatchConfig。
    这时标样名要从 info 里取，阈值用默认 —— 不能抛。
    """
    checks = assess_batch(make_result(), None)
    k = keyed(checks)
    assert k["standards.primary_bias"].data["name"] == "91500"
    assert k["standards.secondary_bias"].data["name"] == "Ple"
    assert len(checks) > 0


def test_thresholds_are_overridable_and_change_only_the_verdict():
    """
    阈值只影响**判词**，不影响任何数值。把 91500 的偏差门槛收紧到 0.1%，
    原本 pass 的那条应当变 warn，而 observed（实测偏差）必须一字不变。
    """
    r = make_result(standards=[("91500", "主标", 5, 1062.4, 1059.75, 0.51),
                               ("Ple", "监控标样", 3, 337.13, 343.24, 1.00)])
    base = keyed(assess_batch(r))["standards.primary_bias"]
    tight = keyed(assess_batch(
        r, None, QCThresholds(primary_bias_warn_pct=0.1)))["standards.primary_bias"]
    assert base.level == PASS
    assert tight.level != PASS
    assert base.observed == tight.observed
    assert abs(base.data["bias_pct"] - tight.data["bias_pct"]) < 1e-12


# ═════════════════════════════════════════════════════════════════════════════
# 给人看的排版
# ═════════════════════════════════════════════════════════════════════════════
#: 造一个能同时产出 fail / warn / info / pass 的批次：
#: 1200 Ma 的测点会被 ADEPT 丢掉（fail），91500 参考值不自洽（warn），
#: 其余照常（info/pass）。
def _mixed_checks():
    r = make_result(
        unknown=[dict(sample="S01", age=440.0, conc=99.0),
                 dict(sample="S01", age=1200.0, conc=99.0),
                 dict(sample="S02", age=445.0, conc=99.0, struct="多域(2)", f206=3.5)],
        standards=[("91500", "主标", 5, 1062.4, 1059.75, 0.51),
                   ("Ple", "监控标样", 3, 337.13, 343.24, 1.00)])
    return assess_batch(r)


def test_verdict_agrees_with_the_checks():
    """
    `verdict()` 是网页界面与交接 JSON 共用的那一份结论。它必须与逐条检查
    一致 —— 否则会出现「headline 说能用、level 是 fail」这种最难被发现的分歧。
    """
    checks = _mixed_checks()
    v = verdict(checks)
    assert v["level"] == worst(checks)
    assert v["counts"] == counts(checks)
    for c in checks:
        if c.level in (FAIL, WARN):
            assert c.title in v["headline"], f"{c.key} 是 {c.level} 却没进 headline"
        else:
            assert c.title not in v["headline"], f"{c.key} 是 {c.level} 却进了 headline"


def test_summary_lines_puts_fail_first_and_collapses_pass():
    """
    排版的两条约定：**最严重的排在最前**；info/pass 只报条数不逐条列
    （逐条铺满控制台会把 fail 挤到屏幕外面）。
    """
    lines = summary_lines(_mixed_checks())
    joined = "\n".join(lines)
    assert lines[0].strip().startswith("质控结论：")
    # 判词与条数都要在第一行里，只读这一行就知道能不能用
    assert "项不通过" in lines[0] and "项通过" in lines[0]
    # fail 的条目出现在 warn 的条目之前
    assert joined.index("[不通过]") < joined.index("[有前提]")
    # info/pass 的标题不该逐条出现，只该有一个汇总行
    assert "；未逐条列出" in joined


def test_summary_lines_does_not_point_at_a_file_that_was_not_written():
    """
    `--no-json` 时清单没有写盘，末行就不能再说「完整清单见 handoff.json」——
    用户照着去找会以为是自己弄丢了。这个小地方专门钉一条，因为它只在
    "不写盘"的那条分支上出错，平时跑默认路径永远看不到。
    """
    checks = _mixed_checks()
    assert "handoff.json" in "\n".join(summary_lines(checks))
    off = "\n".join(summary_lines(checks, detail_ref=None))
    assert "handoff.json" not in off
    assert "本次未写出清单" in off


def test_wrap_keeps_identifiers_whole():
    """
    折行不能把 `Age76_1s`、`0.000%` 这类记号从中间切开 —— 逐字符折行会，
    而技术报告里被切开的记号会让人认不出来是哪个量。
    """
    text = "甲" * 10 + " Age76/Age76_1s 两列"
    lines = _wrap_cjk(text, width=30)
    assert any("Age76/Age76_1s" in ln for ln in lines), lines
    # 内容一个字符都不能丢
    assert "".join(ln.strip() for ln in lines).replace(" ", "") == text.replace(" ", "")


def test_wrap_never_starts_a_line_with_punctuation():
    """
    中文排版的禁则：折下来的行不能以 `。，、／）` 开头，否则读者会以为
    句子在那里结束了。这里对**真实的排版输出**做检查，而不是只测辅助函数 ——
    禁则曾经被同一循环里的"硬切"逻辑撤销过一次，只有端到端才看得出来。
    """
    lines = summary_lines(_mixed_checks())
    for ln in lines:
        body = ln.strip()
        assert not body or body[0] not in _NO_LINE_START, f"行首是禁则字符：{ln!r}"
        # 80 列控制台要放得下（全角按 2 列算）
        assert _display_width(ln) <= 80, f"太宽（{_display_width(ln)} 列）：{ln!r}"


def test_wrap_hard_splits_only_when_a_token_is_wider_than_the_line():
    """
    单词元比整行还宽（例如一条长 URL、一个长样品名）时只能硬切 ——
    硬切是最后手段，但不能切出死循环，也不能丢字符。
    """
    lines = _wrap_cjk("x" * 50, width=20)
    assert len(lines) > 1
    assert all(_display_width(ln) <= 20 for ln in lines)
    assert "".join(ln.strip() for ln in lines) == "x" * 50


if __name__ == "__main__":
    import _selftest                                              # noqa: E402
    raise SystemExit(_selftest.run(globals()))
