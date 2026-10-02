"""
druid 的序列表自检 —— 角色判定与样品名规范化（`io/sequence.py`）：

    python tests/test_io_sequence.py
    python -m pytest tests/            （装了 pytest 时）

为什么这两件小事值得一个测试文件
--------------------------------
角色判定是整个流程的**入口分流**。它一旦判错，后果不是"某处报错"，
而是**整批数据静默走错分支**：

  · 主标没被认出来 → 归一化因子曲线为空 → 最后在空数组上插值崩掉；
  · 玻璃标样（SRM 612/610）被当成样品 → 拿一个 U 加标、Pb 为现代普通铅的
    玻璃去"定年"，报出一个荒谬年龄，还会污染整批的统计；
  · 监控标样没被认出来 → 第 ⑩ 步二次校正拿不到基体匹配的依据。

而这一层的输入是**用户手写的 Excel**，脏数据是常态：`91500` 在 xls 里
读出来是浮点 `91500.0`；样品名前后带空格；大小写不一致（`Ple` / `PLE`）。
所以判定的**第一步必须是规范化**，顺序不能反 —— 本文件最重要的那条测试
（`test_the_float_trap_requires_normalize_first`）就是钉这个顺序。
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from druid.core.constants import (ROLE_GLASS, ROLE_LABEL_CN, ROLE_PRIMARY,   # noqa: E402
                                  ROLE_SECONDARY, ROLE_UNKNOWN, ROLE_VOID,
                                  VOID_NAMES)
from druid.io.sequence import normalize_name, read_sequence, sample_role      # noqa: E402


# ═════════════════════════════════════════════════════════════════════════════
# 一、角色常量
# ═════════════════════════════════════════════════════════════════════════════
def test_five_roles_and_their_labels():
    """
    角色共五种，且每个都要有中文标签（`ROLE_LABEL_CN` 用于写进结果表与网页）。
    标签表缺一项的症状是 **KeyError 在写报表时才炸**，离出错点多很远。

    第五种（`void`）是 2026-10-02 加的**实验者显式作废标记**：它不参与任何
    「样品」统计，但必须有标签，否则结果表的 `类型` 列会写英文原值。
    """
    roles = {ROLE_PRIMARY, ROLE_SECONDARY, ROLE_GLASS, ROLE_VOID, ROLE_UNKNOWN}
    assert len(roles) == 5
    assert roles == {"primary_std", "secondary_std", "glass", "void", "unknown"}
    assert set(ROLE_LABEL_CN) == roles, "标签表必须覆盖全部五种角色"
    assert ROLE_LABEL_CN[ROLE_VOID] == "作废"


def test_void_markers_are_matched_by_whole_name_only():
    """
    ★ 作废标记（2026-10-02）：序列表里**实验者自己写的作废标记**。

    实测某批 4 个测点的样品名**就叫 `wrong`**（年龄 998 / 990 / 836 Ma，
    而该批群体约 155 Ma —— 明显是打偏或打到了 91500 标样上）。这类点若被
    当成普通样品，会静默污染整批的中位年龄与多域率统计。

    判据是**整名相等**、不做子串匹配 —— 否则 `wrong-1`、`skip-2`
    这类真样品名会被误伤；`x`、`bad` 这种太宽泛的词也刻意不收。
    """
    for n in ("wrong", "WRONG", " wrong ", "void", "skip", "废弃", "作废", "无效"):
        assert sample_role(n) == ROLE_VOID, n
    # 整名相等，因此这些**不是**作废点
    for n in ("wrong-1", "wrong1", "skip-2", "voided", "废品", "x", "bad"):
        assert sample_role(n) == ROLE_UNKNOWN, n
    # 不能抢标样 / 玻璃的名分
    assert sample_role("91500") == ROLE_PRIMARY
    assert sample_role("SRM 612") == ROLE_GLASS
    # 词表要能被外部读出来（下游据此决定"哪些名字算作废"）
    assert "wrong" in VOID_NAMES
    assert all(v == v.strip().lower() for v in VOID_NAMES), "词表必须已规范化"


# ═════════════════════════════════════════════════════════════════════════════
# 二、规范化
# ═════════════════════════════════════════════════════════════════════════════
def test_normalize_name_handles_the_three_kinds_of_excel_dirt():
    """
    规范化处理三类脏数据（模块 docstring 里承诺过的三条）：
      ① 浮点整数 `91500.0` → `"91500"`
      ② 首尾空格 `" Ple "`   → `"Ple"`
      ③ 空值 / None / 纯空白 → `""`（调用方按空串判无效行）
    """
    assert normalize_name(91500.0) == "91500"
    assert normalize_name(91500) == "91500"
    assert normalize_name("  Ple ") == "Ple"
    assert normalize_name(None) == ""
    assert normalize_name("   ") == ""
    assert normalize_name("") == ""
    # 非整数浮点**不许**被取整（612.5 不能被悄悄变成 613）
    assert normalize_name(612.5) == "612.5"
    # 字符串里本来就带小数点时不动它 —— 这是 CSV 路径的行为，见下一条
    assert normalize_name("91500.0") == "91500.0"


def test_the_float_trap_requires_normalize_first():
    """
    ★ 本文件最重要的一条：**`sample_role` 自己不处理浮点陷阱**。

    `sample_role("91500.0")` 返回的是 `unknown`（别名表里只有 `"91500"`），
    必须先过 `normalize_name(91500.0) == "91500"` 才能认出来。

    这不是缺陷，是**分层契约**：`_read_rows` 负责规范化、`sample_role`
    只负责判等。但如果哪天有人跳过规范化直接调 `sample_role`，
    症状就是"主标认不出来、整批失去归一化" —— 而且不报错。
    """
    assert sample_role("91500.0") == ROLE_UNKNOWN, "翻转这条说明分层契约变了"
    assert sample_role(normalize_name(91500.0)) == ROLE_PRIMARY
    # 大小写与空格 sample_role 自己是处理的（str.lower + strip）
    assert sample_role(" 91500 ") == ROLE_PRIMARY
    assert sample_role("PLE") == ROLE_SECONDARY


# ═════════════════════════════════════════════════════════════════════════════
# 三、默认别名
# ═════════════════════════════════════════════════════════════════════════════
def test_default_aliases():
    """别名表：主标 3 种写法、监控标样 5 种写法（含带变音符的 plešovice）。"""
    for n in ("91500", "91500 zircon", "91500z"):
        assert sample_role(n) == ROLE_PRIMARY, n
    for n in ("Ple", "PLE", "Plesovice", "pl", "plešovice", "Plešovice"):
        assert sample_role(n) == ROLE_SECONDARY, n


def test_glass_keywords_are_caught_by_substring():
    """
    玻璃标样走**子串**匹配（srm / nist / 612 / 610），所以
    "SRM 612"、"NIST610"、以及被规范化成 `"612"` 的浮点 `612.0` 都能命中。

    为什么玻璃必须单独拦：NIST SRM 612/610 的 U 是人为加标的（~40 ppm）、
    Pb 是现代普通铅，给它"定年"会得到荒谬结果，混进归一化会污染整批。
    """
    assert sample_role("SRM 612") == ROLE_GLASS
    assert sample_role("srm612") == ROLE_GLASS
    assert sample_role("NIST610") == ROLE_GLASS
    assert sample_role(normalize_name(612.0)) == ROLE_GLASS
    assert sample_role(normalize_name(610.0)) == ROLE_GLASS


def test_custom_names_are_recognized():
    """
    用户把自己实际用的标样名传进来时，这些名字要被认出来
    —— 否则它们会被当成"样品"混进定年统计。
    """
    assert sample_role("GJ1", primary="GJ1") == ROLE_PRIMARY
    assert sample_role("gj1", primary="GJ1") == ROLE_PRIMARY, "参数比对也要大小写不敏感"
    assert sample_role("Temora", secondary="Temora") == ROLE_SECONDARY


def test_custom_names_do_not_disable_the_default_aliases():
    """
    ⚠ **实测行为，与两处文档不符**：别名表是**无条件**生效的。

    `sample_role` 的注释写着别名是"不指定参数时的兜底"、
    `read_sequence` 的 docstring 写着自定义名"会**覆盖**默认的 91500 / Ple 判定"
    —— 但代码里别名那两条 `if` 不带任何条件，所以即使指定了
    `primary="GJ1"`，`"91500"` 仍被判为主标（`"Ple"` 同理）。

    这里有两条出路：改文档，或把别名改成"仅当参数为默认值时才启用"。
    后者会改变自定义标样批次的行为（本仓 25 批全用 91500/Ple，所以数值上中性），
    属于需要人工拍板的事 —— 在此之前先把**现状**钉住。
    """
    assert sample_role("91500", primary="GJ1") == ROLE_PRIMARY
    assert sample_role("Ple", secondary="Temora") == ROLE_SECONDARY
    # 自定义名与默认别名可以同时成立（一批里出现两个"主标"）
    assert sample_role("GJ1", primary="GJ1") == ROLE_PRIMARY


def test_unknown_is_the_safe_fallback():
    """认不出来的一律落到 unknown（真正的定年对象），不猜、不报错。"""
    for n in ("S01", "SAMPLE-A1", "", "X", "样品甲"):
        assert sample_role(n) == ROLE_UNKNOWN, n


# ═════════════════════════════════════════════════════════════════════════════
# 四、端到端：读一个真实格式的序列文件
# ═════════════════════════════════════════════════════════════════════════════
def test_read_sequence_normalizes_then_assigns_roles():
    """
    两列 CSV（文件名 / 样品名），无表头 —— 与 `_read_rows` 的 csv 分支一致。

    这里验证的是**顺序**：规范化在角色判定之前发生，所以样例名里的空格被去掉。
    返回表必须正好四列（file / sample / role / order），`order` 从 1 起。
    """
    rows = ["EX2022A_1,91500", "EX2022A_2, Ple ", "EX2022A_3,SRM 612", "EX2022A_4,S01"]
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "T_LIST.csv"
        p.write_text("\n".join(rows) + "\n", encoding="utf-8")
        df = read_sequence(p)

    assert list(df.columns) == ["file", "sample", "role", "order"]
    assert list(df["sample"]) == ["91500", "Ple", "SRM 612", "S01"]
    assert list(df["role"]) == [ROLE_PRIMARY, ROLE_SECONDARY, ROLE_GLASS, ROLE_UNKNOWN]
    assert list(df["order"]) == [1, 2, 3, 4]


def test_csv_path_depends_on_pandas_dtype_inference():
    """
    ⚠ 一条**脆弱点**，值得知道但不必现在改：csv 路径下 `91500.0` 能不能被救回来，
    取决于 **pandas 把那一列推断成什么 dtype**，而 dtype 又取决于**同列里还有什么**。

    实测（同两份文件，只差一行）：

      · 只有一行 `EX_1,91500.0`   → pandas 推断为**浮点** → `normalize_name`
        拿到 91500.0 → 救回 `"91500"` → **primary_std**；
      · 再加一行 `EX_2,Ple`       → 推断为**字符串** → `normalize_name` 拿到
        字符串 `"91500.0"`（`isinstance(..., float)` 不成立，不动它）
        → **unknown**。

    也就是说，同一份序列内容**存成 csv 还是 xls、甚至同一列里多一行别的样品**，
    主标可能一个认得出、一个认不出 —— 而且都不报错。
    本仓的序列表是 xls/xlsx（xlrd 一律给浮点），走的是被救回来的那条路，
    所以生产路径不受影响；但"将来支持 csv 序列"时必须一起处理。
    """
    with tempfile.TemporaryDirectory() as d:
        # ① 整列都是数字 → pandas 推成浮点 → 被救回
        p1 = Path(d) / "ONLY_NUM_LIST.csv"
        p1.write_text("EX_1,91500.0\n", encoding="utf-8")
        df1 = read_sequence(p1)
        assert list(df1["sample"]) == ["91500"]
        assert list(df1["role"]) == [ROLE_PRIMARY]

        # ② 同列里混进非数字 → pandas 推成字符串 → 救不回来
        p2 = Path(d) / "MIXED_LIST.csv"
        p2.write_text("EX_1,91500.0\nEX_2,Ple\n", encoding="utf-8")
        df2 = read_sequence(p2)
        assert list(df2["sample"]) == ["91500.0", "Ple"]
        assert list(df2["role"]) == [ROLE_UNKNOWN, ROLE_SECONDARY]


def _run_standalone() -> int:
    """见 tests/_selftest.py —— 让这个文件不装 pytest 也能直接跑。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _selftest
    return _selftest.run(globals())


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
