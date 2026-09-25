"""
druid.core 参考值预设（`REFERENCE_PRESETS`）自检 —— 不依赖 pytest：
==============================

    python tests/test_reference_presets.py
    python -m pytest tests/test_reference_presets.py

为什么这个文件必须存在
----------------------
`BatchConfig.ref_preset` 是一个**会改变全部年龄**的开关：它决定标样的参考值
是取"ID-TIMS 实测比值"还是"由年龄反算"。这类开关有两个必须钉死的性质：

1. **`repo` 档必须逐位等于 STANDARDS。** 2026-09-25 起默认档改为 horstwood2016，
   于是 `repo` 成了「复现 2026-09-25 之前结果」的唯一入口 —— 它必须永远等于第一版，
   否则历史结果再也复现不出来。
2. **选了档必须真的换。** 否则用户以为换了口径、实际没换（比报错更糟）。

顺带钉住：预设只改它自己列出的标样、每条必须带出处、拼错必须报错。
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import _selftest                                                    # noqa: E402
from druid.core.constants import REFERENCE_PRESETS, STANDARDS       # noqa: E402
from druid.core.geochronology import age68, r68_of_age, r76_of_age  # noqa: E402
from druid.core.references import std_age, std_ref                  # noqa: E402
from druid.workflow import BatchConfig                              # noqa: E402


# ═════════════════════════════════════════════════════════════════════════════
# 1. 默认档：不选就一点都不能变
# ═════════════════════════════════════════════════════════════════════════════
def test_default_preset_is_horstwood2016_and_repo_table_stays_empty():
    """默认档 2026-09-25 起是 horstwood2016；`repo` 必须仍是**空表**（可复现历史）。

    `repo` 是"按第一版口径复现旧结果"的唯一入口 —— 它一旦非空，就说明
    有人往历史档里塞了新数，那条路就不再可信了。
    """
    assert BatchConfig(data_dir=".").ref_preset == "horstwood2016"
    assert REFERENCE_PRESETS["repo"] == {}
    # 默认档要真的生效（不能只是字段上写着）
    assert std_ref("91500") == std_ref("91500", "horstwood2016")
    assert std_age("91500") == std_age("91500", "horstwood2016") == 1063.51


def test_repo_preset_equals_standards_bit_for_bit():
    """选 repo（或空串）时，取值必须与直接读 STANDARDS 逐位相同。

    ⚠ 这里**不能**拿 `std_ref(name)`（默认档）去比 —— 2026-09-25 起默认档已是
    horstwood2016，两者对 91500/GJ1/Plešovice/MudTank 本来就不相等。
    """
    for name in STANDARDS:
        assert std_ref(name, "repo") == std_ref(name, "")
        assert std_age(name, "repo") == std_age(name, "")
    # 直接把 91500 的**第一版**现值钉在这里：动了它会红。
    assert std_ref("91500", "repo") == (0.17928, 0.07494)
    assert std_age("91500", "repo") == 1062.4


# ═════════════════════════════════════════════════════════════════════════════
# 2. 各档的数值：逐字对齐来源
# ═════════════════════════════════════════════════════════════════════════════
def test_horstwood2016_matches_table_s2():
    """horstwood2016 的值必须逐字等于 GGR 补充材料 Table S2 块 1（CA-ID-TIMS）。"""
    assert std_ref("91500", "horstwood2016") == (0.179365, 0.074941)
    assert std_age("91500", "horstwood2016") == 1063.51
    assert std_ref("GJ1", "horstwood2016") == (0.09786, 0.060139)
    assert std_age("GJ1", "horstwood2016") == 601.86
    assert std_ref("Plesovice", "horstwood2016") == (0.053694, 0.053244)
    assert std_age("Plesovice", "horstwood2016") == 337.16
    assert std_ref("MudTank", "horstwood2016") == (0.120188, 0.063802)
    assert std_age("MudTank", "horstwood2016") == 731.65


def test_wiedenbeck1995_matches_the_paper():
    """Wiedenbeck 1995 的 91500：R68=0.17917 / R76=0.07488（ID-TIMS 原始文献）。"""
    assert std_ref("91500", "wiedenbeck1995") == (0.17917, 0.07488)
    assert std_age("91500", "wiedenbeck1995") == 1062.4


def test_age_only_presets_are_back_calculated():
    """R68/R76 为 None 的条目要由年龄反算，且与直接调 r68_of_age/r76_of_age 一致。"""
    for preset, age in (("self-consistent", 1062.4), ("isoclock", 1062.0)):
        r68, r76 = std_ref("91500", preset)
        assert abs(r68 - float(r68_of_age(age))) < 1e-15
        assert abs(r76 - float(r76_of_age(age))) < 1e-15
        assert std_age("91500", preset) == age


def test_self_consistent_preset_actually_removes_the_known_inconsistency():
    """选 self-consistent 后，91500 的比值与年龄必须自洽（差 ~0）。

    这正是"用年龄反算"这条路线的全部卖点：比值与年龄不可能互相打架。
    默认档（repo）在这一条上是 **+0.0602%** —— 差别就在这里，不是纠错。
    """
    r68, _ = std_ref("91500", "self-consistent")
    a = float(age68(r68))
    assert abs(a / std_age("91500", "self-consistent") - 1) < 1e-9
    # 顺手确认**第一版**口径确实不自洽，否则这个测试就失去意义了
    # （必须显式 "repo"：默认档 horstwood2016 本身是自洽的）
    r68_repo, _ = std_ref("91500", "repo")
    assert abs(float(age68(r68_repo)) / std_age("91500", "repo") - 1) > 5e-4


# ═════════════════════════════════════════════════════════════════════════════
# 3. 契约：只改列出的、必须带出处、拼错必须报错
# ═════════════════════════════════════════════════════════════════════════════
def test_preset_only_touches_the_standards_it_lists():
    """预设没列出的标样必须不受影响 —— 这是「只改它列出的」契约。"""
    # Temora1 不在 horstwood2016 里
    assert std_ref("Temora1", "horstwood2016") == std_ref("Temora1")
    assert std_age("Temora1", "horstwood2016") == std_age("Temora1") == 416.78
    # weidenbeck1995 / isoclock 只列了 91500
    assert std_age("Plesovice", "isoclock") == 337.13
    assert std_age("Plesovice", "wiedenbeck1995") == 337.13


def test_every_preset_entry_has_provenance():
    """每条预设条目都要带 ref（出处）—— 数值可以不同，出处不能没有。"""
    for pname, table in REFERENCE_PRESETS.items():
        for sname, entry in table.items():
            assert sname in STANDARDS, f"{pname}/{sname} 不在 STANDARDS 里"
            assert set(entry) == {"R68", "R76", "age_Ma", "ref"}, \
                f"{pname}/{sname} 字段不对：{sorted(entry)}"
            assert entry["ref"], f"{pname}/{sname} 缺 ref 出处"


def test_unknown_preset_fails_loudly():
    """预设名拼错必须**报错**，不能静默退回默认（否则用户以为换了口径）。"""
    for fn in (std_ref, std_age):
        try:
            fn("91500", "no-such-preset")
        except KeyError:
            pass
        else:
            raise AssertionError(f"{fn.__name__} 没有对未知 preset 报错")
    # BatchConfig 也要在建对象时就挡住
    try:
        BatchConfig(data_dir=".", ref_preset="no-such-preset")
    except ValueError:
        pass
    else:
        raise AssertionError("BatchConfig 没有对未知 ref_preset 报错")


if __name__ == "__main__":
    raise SystemExit(_selftest.run(globals()))
