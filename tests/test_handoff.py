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

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _fixtures import FakeResult, make_result                           # noqa: E402
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
    "name", "n_spots", "age68", "age68_pooled_diagnostic", "age68_robustness",
    "concordance",
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
# 稳健性：剔掉最偏离的 k% 之后中位怎么动
# ═════════════════════════════════════════════════════════════════════════════
def test_robustness_reports_the_median_shift_at_each_trim_level():
    """
    这一块是论文那句"n 大时误差被压小"的仓库版：让人能机读地检查
    "这个中位背后有多少测点在拉偏"，而不是只拿到一个看起来很稳的数。

    判据是**可复算**的：给一组人为构造的数据（一个明显的离群点 + 一堆
    紧密聚在一起的），剔掉 20% 之后离群点必然 gone、中位必须回到聚堆中心。
    关键：**中位本身几乎不动**（中位数对单个离群点本就稳健），
    动的是"剔了多少、还剩多少"这两个数。

    ⚠ 基准值是 **445.5** 而不是 450：`age68_robustness` 与 `age68` 一样
    走 `age68_column()`，在已校准的假数据上那是 `年龄206_238_QC校正`
    （= 年龄 × 0.99）。这里必须用真值断言，否则改一次口径就会静默错。
    """
    ages = [450.0] * 9 + [800.0]                       # 10 个点，1 个离群
    p = build_payload(None, make_result(unknown=[
        dict(sample="SAME", age=a, s1=2.0) for a in ages]), version="2.3.0")
    rb = p["samples"][0]["age68_robustness"]

    base = 450.0 * 0.99                                # 走的是 QC 校正列
    assert rb["n_input"] == 10
    assert abs(rb["median_full"] - base) < 1e-9
    lv = {v["frac"]: v for v in rb["levels"]}
    assert set(lv) == {0.05, 0.10, 0.20}
    # 20% ⇒ 剔 2 个；离群点(792)必在其中，中位仍是聚堆中心
    assert lv[0.20]["n_dropped"] == 2
    assert lv[0.20]["n_remaining"] == 8
    assert abs(lv[0.20]["median_Ma"] - base) < 1e-9
    assert abs(lv[0.20]["shift_Ma"]) < 1e-9, "中位对单个离群点应当不动"
    # shift 与 shift_pct 必须自洽（下游常常只要百分比那个）
    assert abs(lv[0.20]["shift_pct"] - lv[0.20]["shift_Ma"] / base * 100) < 1e-9


def test_robustness_never_reports_a_sigma_change():
    """
    ★ 刻意**不报** σ 的变化，这是本块最容易出事的地方。

    剔除量是 |age − 中位| ⁄ σ ⇒ **离群点天然就是 σ 大的点**（示例批次实测
    corr(|age−中位|, σ) = 0.822；被剔的 10 个点平均 σ 11.96 Ma，
    留下的 38 个只有 6.68 Ma）。所以"剔得越多、平均 1σ 越小"是**定义
    带来的**，不是测量变精了。

    一旦把 σ 的变化一并报出去，读者必然会读成"精度提高了" ——
    那正是这块数据最容易被误读的地方。所以测试在这里把"不许报"钉死：
    任何一层里出现 σ / sigma / se 之类的键，本条就红。
    """
    ages = [450.0] * 9 + [800.0]
    p = build_payload(None, make_result(unknown=[
        dict(sample="SAME", age=a, s1=1.0 if a < 500 else 30.0) for a in ages]),
        version="2.3.0")
    rb = p["samples"][0]["age68_robustness"]
    blob = json.dumps(rb, ensure_ascii=False)
    for bad in ("sigma_Ma", "mean_sigma", "s1_median", "sigma_change", "se_"):
        assert bad not in blob, f"稳健性块不许报 σ 的变化，却出现了 {bad}"
    # 但必须把这件事**写明**，否则读者仍会自己去算
    assert "不代表测量变精" in rb["sigma_trend_note"]


def test_robustness_survives_a_duplicated_index():
    """
    ★ 这条是**同轮真踩到的坑**：`x` 与 `sig` 必须落在**同一个 0..n-1 坐标系**里。

    实测崩过，报 `index 144 is out of bounds for axis 0 with size 144`，
    症状是整批崩在 `build_payload` 里，而不是安静地算错。

    ⚠⚠ 触发它的**不是** NaN，而是**索引重复**：写成
    `sig.loc[x.index]` 之后，若 `x.index` 里有重复值，返回的 Series 会
    **按重复次数膨胀**。实测：x 30 行、索引 0..9 各重复 3 次 ⇒
    `sig.loc[x.index]` 变成 **90 行**，而 `cur` 里的位置最大只有 29 ⇒ 越界。
    第一版守护只造了 NaN，索引唯一 ⇒ `loc` 与 `reindex` 结果完全一样，
    **守护抓不到**，负向测试当场放行。这条就是把那次的缺口补上。
    """
    rows = [dict(sample="SAME", age=450.0 + i, s1=2.0) for i in range(10)]
    r = make_result(unknown=rows)
    res = r.results
    # 造重复索引：**同一行**复制三次（索引 0 会出现 3 次）。
    # ⚠ 用 `iloc[[0,1,2]]` 那种"三行叠起来"造不出重复 —— 它们的索引是
    #   0/1/2，本来就不同；只有复制**同一行**才会重复。
    u = res[res["类型"] == "样品"]
    dup = pd.concat([u.iloc[[0]], u.iloc[[0]], u.iloc[[0]]])
    assert dup.index.duplicated().any(), "这个测试的前提就是索引重复"
    # 只留样品行，避免标样行混进来改变形状
    dup = dup[dup["类型"] == "样品"]
    rb = build_payload(None, FakeResult(dup, pd.DataFrame(), dict(dup.iloc[0]),
                                        pd.DataFrame(), pd.DataFrame()),
                       version="2.3.0")["samples"][0]["age68_robustness"]
    assert rb["n_input"] == 3
    for v in rb["levels"]:
        assert v["n_remaining"] <= rb["n_input"], "剩余点数不该超过输入点数"


def test_robustness_counts_only_the_rows_it_actually_used():
    """
    索引重复时**不能**把同一行数两遍：`n_input` 必须是真正参与计算的点数。

    上一条钉的是"不崩"，这条钉的是"**不虚报**" —— 崩溃修好了但把 3 行
    报成 9 行，同样是错的，而且更安静。
    """
    rows = [dict(sample="SAME", age=450.0 + i, s1=2.0) for i in range(6)]
    r = make_result(unknown=rows)
    res = r.results
    u = res[res["类型"] == "样品"]
    dup = pd.concat([u.iloc[[0]], u.iloc[[0]]])
    rb = build_payload(None, FakeResult(dup, pd.DataFrame(), dict(dup.iloc[0]),
                                        pd.DataFrame(), pd.DataFrame()),
                       version="2.3.0")["samples"][0]["age68_robustness"]
    assert rb["n_input"] == 2, f"索引重复不该让 n 虚增，实得 {rb['n_input']}"


def test_robustness_says_nothing_on_a_single_spot():
    """
    n = 1 时"稳健性"是个没有意义的词：任何 k% 都剔不到点。
    这时必须给 `levels: []`，而不是给一串 shift = 0 让人以为"很稳健"。

    示例批次正好就是这个情形（48 个岩样各 1 个测点），所以真实批次
    跑出来的 48 条 `age68_robustness` 全是空的 —— 这是**如实**，
    不是没算出来。
    """
    p = build_payload(None, make_result(unknown=[
        dict(sample="A", age=450.0, s1=3.0), dict(sample="B", age=460.0, s1=3.0)]),
        version="2.3.0")
    for s in p["samples"]:
        rb = s["age68_robustness"]
        assert rb["levels"] == [], "单点样品不该有 levels"
        assert rb["n_input"] == 1
        assert rb["median_full"] is not None, "但仍然要给出那个中位本身"


def test_robustness_keeps_at_least_three_spots():
    """
    k% × n 向下取整、且**至少留 3 个点**。小样本上"至少留 3"会盖过 k%
    （n=4、k=20% ⇒ 本想剔 0 个；n=5、k=40% ⇒ 本想剔 2 个但只能剔 1 个）——
    少于 3 个点的中位没有任何意义。

    钉的是 `n_remaining >= 3`，不是"k% 一定剔得到 k 个"。
    """
    p = build_payload(None, make_result(unknown=[
        dict(sample="SAME", age=450.0 + 40 * i, s1=1.0) for i in range(5)]),
        version="2.3.0")
    rb = p["samples"][0]["age68_robustness"]
    for v in rb["levels"]:
        assert v["n_remaining"] >= 3, v
        assert v["n_dropped"] <= rb["n_input"] - 3, v


def test_robustness_is_registered_in_the_sample_key_whitelist():
    """
    `SAMPLE_KEYS` 是**白名单**：样品块里出现未登记的键，测试就红。
    新增字段必须登记 —— 否则下游会拿到一个没人认识、也没人保证的字段。
    """
    assert "age68_robustness" in SAMPLE_KEYS
    p = build_payload(None, make_result(unknown=[
        dict(sample="SAME", age=450.0, s1=2.0)]), version="2.3.0")
    for s in p["samples"]:
        assert set(s) <= SAMPLE_KEYS, sorted(set(s) - SAMPLE_KEYS)


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
