"""
druid.io.handoff 自检 —— 不依赖 pytest，直接跑也认：

    python tests/test_handoff.py
    python -m pytest tests/test_handoff.py

为什么这个文件必须存在
----------------------
交接 JSON 是**给程序用的契约**：下游（自动化质控/解释流程、ADEPT 封装）按字段名取值。
它和 Excel 不一样 —— Excel 少一列，人会看见；JSON 少一个键，下游要么崩、
要么悄悄走进 `if key in data` 的 else 分支，拿一个默认值继续算。

所以这里钉三件事：

1. **键集合**。改字段名会当场红，不会悄悄漏给下游。
2. **两个序列化坑**：NaN 不是合法 JSON；numpy 标量不可序列化。
   这两个错只在真有数据的路径上出现，手测很容易漏。
3. **几个刻意的"不给"**：不给每个样品的加权平均年龄、不把 `_1s` 说成 2σ。
   这类决定不像 bug 那样会红，只能靠测试钉住，否则半年后会被人"顺手补上"。
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _fixtures import make_result                                    # noqa: E402
from druid.io.handoff import (                                       # noqa: E402
    ADEPT_INPUT, SCHEMA, build_payload, default_path, export_handoff,
)

#: 顶层键集合。**增删都要同步这里**，否则 `test_top_level_keys_are_stable`
#: 会红。但注意：`druid.handoff` 的 schema 版本**只在做破坏性改动时才升**
#: （改字段名 / 改语义 / 改类型），增字段属于兼容扩展，不升 —— 见 handoff.py。
TOP_KEYS = {
    "schema", "generated_at", "tool", "batch", "verdict", "config", "counts",
    "conventions", "calibration", "standards", "samples", "whole_spot",
    "checks", "handoff", "intentionally_omitted",
}

#: 每个样品允许出现的键。**故意没有加权平均年龄** —— 见模块 docstring。
SAMPLE_KEYS = {
    "name", "n_spots", "age68", "age68_pooled_diagnostic", "concordance",
    "f206_pct_median", "u238_cps_median", "th_u_median",
    "depth_structures", "n_multi_domain",
}

AGE68_STAT_KEYS = {"n", "median", "q05", "q25", "q75", "q95", "min", "max"}


def _strict_load(text: str):
    """
    严格解析：`parse_constant` 只在遇到 `NaN` / `Infinity` / `-Infinity` 这三种
    **非法 JSON 字面量**时被调用（Python 的 json 默认容忍它们）。
    所以这里让它抛 —— 能抛就说明文件里有非法字面量。
    """
    def boom(tok):
        raise AssertionError(f"JSON 里出现非法字面量：{tok}")

    return json.loads(text, parse_constant=boom)


# ═════════════════════════════════════════════════════════════════════════════
# 契约的骨架
# ═════════════════════════════════════════════════════════════════════════════
def test_top_level_keys_are_stable():
    p = build_payload(None, make_result(), version="2.3.0")
    assert set(p) == TOP_KEYS, sorted(set(p) ^ TOP_KEYS)
    assert p["schema"] == SCHEMA
    assert p["tool"] == {"name": "druid", "version": "2.3.0"}


def test_verdict_block_agrees_with_the_checks():
    p = build_payload(None, make_result(), version="2.3.0")
    levels = [c["level"] for c in p["checks"]]
    n = p["verdict"]["counts"]
    assert n["fail"] + n["warn"] + n["info"] + n["pass"] == len(levels)
    assert n["pass"] == levels.count("pass")
    # headline 只摘 fail/warn，pass 的标题不该出现在里面
    for c in p["checks"]:
        if c["level"] == "pass":
            assert c["title"] not in p["verdict"]["headline"]


def test_check_entries_carry_the_full_shape():
    p = build_payload(None, make_result(), version="2.3.0")
    for c in p["checks"]:
        assert set(c) == {"key", "level", "title", "observed", "criterion",
                          "detail", "data"}, set(c)
        assert c["key"].isascii(), c["key"]
        assert c["level"] in ("fail", "warn", "info", "pass")


def test_conventions_declare_the_sign_conventions():
    """
    交接文件不只要能取到数，还要能取到"这个数是什么意思"。
    `_1s` 是 1σ 还是 2σ，拿错 MSWD 会差 4 倍 —— 这条必须写在文件里，
    不能只写在文档里靠人记。
    """
    p = build_payload(None, make_result(), version="2.3.0")
    c = p["conventions"]
    assert c["age68_1s_sigma"] == 1
    assert c["age_unit"] == "Ma"
    assert c["age_column_used"] in ("年龄206_238", "年龄206_238_QC校正")
    assert "反比方差加权" in c["wmean_definition"]


# ═════════════════════════════════════════════════════════════════════════════
# 两个序列化坑
# ═════════════════════════════════════════════════════════════════════════════
def test_payload_is_strict_json_without_nan_literals():
    p = build_payload(None, make_result(), version="2.3.0")
    text = json.dumps(p, ensure_ascii=False, allow_nan=False)
    _strict_load(text)                       # 出现 NaN/Infinity 会当场 AssertionError
    assert "NaN" not in text and "Infinity" not in text


def test_nan_values_become_null_not_the_nan_token():
    """
    NaN 必须变 `null`。Python 自己能把裸 `NaN` 读回来，所以本地手测会"通过"；
    别的语言的解析器一律报错 —— 这个坑只在跨语言时才现形。
    """
    r = make_result(unknown=[
        dict(sample="SOK", age=450.0),
        dict(sample="SNAN", age=float("nan"), u=float("nan"), conc=float("nan")),
    ])
    p = build_payload(None, r, version="2.3.0")
    text = json.dumps(p, ensure_ascii=False, allow_nan=False)
    back = _strict_load(text)
    nan_sample = next(s for s in back["samples"] if s["name"] == "SNAN")
    assert nan_sample["age68"]["n"] == 0
    assert nan_sample["age68"]["median"] is None
    assert nan_sample["u238_cps_median"] is None


def test_numpy_scalars_are_converted_to_native_types():
    """
    `np.int64` 不是 `int`、`np.bool_` 不是 `bool`，json 直接抛 TypeError。
    转换之后必须能被 json 序列化，且**能再读回来**（说明是真原生类型）。
    """
    p = build_payload(None, make_result(unknown=[
        dict(sample=f"S{i}", age=450.0 + i) for i in range(6)]), version="2.3.0")
    text = json.dumps(p, ensure_ascii=False, allow_nan=False)     # 不抛就说明转干净了
    back = _strict_load(text)
    s = back["samples"][0]
    assert isinstance(s["n_spots"], int) and not isinstance(s["n_spots"], bool)
    assert isinstance(s["age68"]["n"], int)
    assert isinstance(s["age68"]["median"], float)
    assert isinstance(back["verdict"]["counts"]["pass"], int)


# ═════════════════════════════════════════════════════════════════════════════
# 刻意的"不给"：这类决定只能靠测试钉
# ═════════════════════════════════════════════════════════════════════════════
def test_no_per_sample_weighted_mean_age():
    """
    一个样品里各测点的散布通常远大于各自误差（实测池化 MSWD 可达 10³~10⁴），
    加权平均只能靠剔掉大部分测点来达标 —— 得到的数看着精确、其实是被挑出来的子集。

    所以这里**只给描述统计**，并附上"要达标得剔掉几个测点"当证据。
    以后若有人想"顺手补一个岩样加权平均年龄"，这条会红。
    """
    p = build_payload(None, make_result(unknown=[
        dict(sample="SAME", age=450.0 + 3 * i, s1=1.0) for i in range(8)]),
        version="2.3.0")
    for s in p["samples"]:
        assert set(s) <= SAMPLE_KEYS, sorted(set(s) - SAMPLE_KEYS)
        assert set(s["age68"]) == AGE68_STAT_KEYS, sorted(s["age68"])
        # 描述统计里不许出现任何"均值"字段
        assert not any("mean" in k for k in s["age68"]), s["age68"]
    assert "per_sample_weighted_mean_age" in p["intentionally_omitted"]


def test_pooled_mswd_is_reported_as_a_diagnostic_with_the_drop_count():
    """
    pool 诊断量要给两个东西：MSWD 本身，以及"把 MSWD 压到阈值要剔掉多少个测点"。
    只给 MSWD 会被误当成"这个样品均一"；后半句才是它到底有多勉强的证据。
    """
    p = build_payload(None, make_result(unknown=[
        dict(sample="SAME", age=450.0 + 30 * i, s1=1.0) for i in range(8)]),
        version="2.3.0")
    d = p["samples"][0]["age68_pooled_diagnostic"]
    assert set(d) == {"n_spots", "mswd", "n_dropped_for_threshold", "threshold", "prob"}
    assert d["n_spots"] == 8
    assert d["mswd"] is not None and d["mswd"] > 2.5
    assert d["n_dropped_for_threshold"] > 0


# ═════════════════════════════════════════════════════════════════════════════
# ADEPT 交接块
# ═════════════════════════════════════════════════════════════════════════════
def test_adept_block_pins_the_call_parameters():
    """
    四个调用参数都必须**显式**给：ADEPT 刻意不猜数据是否已预处理。
    `age68_1s_sigma = 1` 更是硬约束 —— 拿 2σ 当 1σ 喂进去，MSWD 会差 4 倍。
    """
    p = build_payload(None, make_result(), version="2.3.0")
    a = p["handoff"]["adept"]
    assert a["sheet"] == ADEPT_INPUT["sheet"] == "剖面窗口"
    assert a["columns"] == ["Analysis", "Time", "Age68", "Age68_1s"]
    assert a["age68_1s_sigma"] == 1
    assert a["call"] == {"lower_ablation_time": 0, "upper_ablation_time": 1000000,
                         "smooth": "none", "calibration_uncertainty": 0}


def test_adept_block_agrees_with_the_dropout_check():
    """
    交接块里的掉点信息必须与质控检查项一致 —— 两处各算一遍迟早会对不上，
    所以交接块是**读**检查项拿的。这条测试就是钉住这一点。
    """
    r = make_result(unknown=[dict(sample="GOOD", age=450.0, windows=[450.0] * 30),
                             dict(sample="OLD", age=1150.0, windows=[1150.0] * 30)])
    checks = build_payload(None, r, version="2.3.0")["checks"]
    chk = next(c for c in checks if c["key"] == "handoff.adept_dropout")
    a = build_payload(None, r, version="2.3.0")["handoff"]["adept"]
    assert a["will_drop_n"] == chk["data"]["n_dropped"] == 1
    assert a["will_drop_spots"] == chk["data"]["dropped"]
    assert a["will_drop_spots"][0]["sample"] == "OLD"


# ═════════════════════════════════════════════════════════════════════════════
# 路径与写盘
# ═════════════════════════════════════════════════════════════════════════════
def test_whole_spot_block_lists_every_spot_with_a_compatibility_flag():
    """
    每个样品测点都要有一行 —— 这正是它与 `samples`（岩样级描述统计）
    以及「深度剖面域」表（只收多域点）的区别。`mswd_compatible` 把
    "这个整段平均值能不能用"变成可机读的布尔量，下游不必自己猜阈值。
    """
    p = build_payload(None, make_result(unknown=[
        dict(sample="A", age=450.0, whole_ok=True),
        dict(sample="B", age=1150.0, struct="多域(2)", whole_ok=False,
             whole_mswd=50.0),
    ]), version="2.4.0")
    w = p["whole_spot"]
    assert w["n_spots"] == 2
    assert w["n_compatible"] == 1 and abs(w["compatible_frac"] - 0.5) < 1e-12
    assert [d["sample"] for d in w["per_spot"]] == ["A", "B"]
    a, b = w["per_spot"]
    assert a["mswd_compatible"] is True and b["mswd_compatible"] is False
    assert b["n_domains"] == 2
    assert isinstance(a["age_ma"], float) and a["age_ma"] == 450.0
    # 键名一律 ASCII（下游可能是任何语言写的），中文只出现在说明性字符串里
    for d in w["per_spot"]:
        assert all(k.isascii() for k in d), sorted(d)
    # numpy 标量与 NaN 两个序列化坑对每个块都适用，这里也要过一遍
    _strict_load(json.dumps(p, ensure_ascii=False, allow_nan=False))


def test_whole_spot_block_is_empty_when_the_table_is_absent():
    """没有这张表就给空 dict —— 不要顺手编一份出来（编的那份没人知道是假的）。"""
    p = build_payload(None, make_result(whole=[]), version="2.4.0")
    assert p["whole_spot"] == {}


def test_default_path_sits_next_to_the_excel():
    class _C:
        out_excel = Path("D:/out/BATCH_U-Pb结果.xlsx")

    assert default_path(_C).name == "BATCH_U-Pb结果.handoff.json"
    assert default_path(_C).parent == Path("D:/out")

    class _C2:
        out_excel = Path("D:/out/noext")

    assert default_path(_C2).name == "noext.handoff.json"


def test_export_writes_utf8_and_leaves_no_tmp_behind():
    """
    原子写：先写 `.tmp` 再 `os.replace`。下游可能在**轮询**这个文件，
    读到写了一半的 JSON 会解析失败 —— 所以必须确保换上去的永远是完整内容，
    而且不留残留的 `.tmp`（否则下次轮询可能读到旧的那份）。
    """
    r = make_result()
    with tempfile.TemporaryDirectory() as td:
        target = Path(td) / "sub" / "B.handoff.json"
        p = export_handoff(None, r, version="2.3.0", path=target)
        assert p == target and p.exists()
        text = p.read_text(encoding="utf-8")
        assert "剖" in text, "中文应原样写出（ensure_ascii=False）"
        _strict_load(text)
        leftovers = list((Path(td) / "sub").glob("*.tmp"))
        assert not leftovers, leftovers


def test_build_payload_does_not_touch_the_disk():
    """`build_payload` 必须能在内存里跑完 —— 测试与检查都靠它。"""
    with tempfile.TemporaryDirectory() as td:
        before = set(Path(td).iterdir())
        build_payload(None, make_result(), version="2.3.0")
        assert set(Path(td).iterdir()) == before


def test_payload_works_with_a_real_batch_config_shape():
    """
    传真的 BatchConfig 时，`config` 块要把私有字段（下划线开头）挡在外面，
    并把嵌套的阈值对象摊平 —— 否则下游会拿到 `_out_dir` 这种内部路径。
    """
    from druid.workflow import BatchConfig

    cfg = BatchConfig(data_dir="examples/EX2022A")
    p = build_payload(cfg, make_result(), version="2.3.0")
    assert "data_dir" in p["config"]
    assert not [k for k in p["config"] if k.startswith("_")]
    assert isinstance(p["config"]["qc_thresholds"], dict)
    assert p["config"]["qc_thresholds"]["min_usable_windows"] == 10


if __name__ == "__main__":
    import _selftest                                              # noqa: E402
    raise SystemExit(_selftest.run(globals()))
