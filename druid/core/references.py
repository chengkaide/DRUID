"""
druid.core.references —— 按标准物质名取参考比值
=============================================

职责单一：把"标准物质叫什么名字"翻译成"它的参考比值是多少"。

为什么要绕一层
--------------
constants.STANDARDS 里，有些标样给了文献直接列出的比值（91500），
有些只给了公认年龄（Plesovice、GJ1…）。如果每次用都手写
`if R68 is None: 从年龄反算`，这段判断会散落各处。

集中到这里后：
    · 用到 ratio / 年龄的地方统一调用 std_ref() / std_age()，永远是同一套规则；
    · 换"参考值口径"（ID-TIMS 实测比值 ↔ 由年龄反算）只改这一个文件 ——
      候选口径写在 constants.REFERENCE_PRESETS，用 `preset` 参数选；
      默认档由 constants.DEFAULT_REF_PRESET 统一给出（2026-09-25 起为
      "horstwood2016"）。传 "repo" 可复现 2026-09-25 之前的第一版口径。
"""
from __future__ import annotations

from .constants import DEFAULT_REF_PRESET, REFERENCE_PRESETS, STANDARDS
from .geochronology import r68_of_age, r76_of_age


def _preset_table(preset=DEFAULT_REF_PRESET):
    """
    把一个预设名解析成"覆盖表"（标样名 → 参考值条目）。

    "repo"（或空）→ 返回空表，即不覆盖任何标样（＝第一版口径，可复现历史）。
    默认档取自 constants.DEFAULT_REF_PRESET（见函数签名的默认值）——
    **不在这里写死默认名**，免得与 BatchConfig 的默认值各说各话。
    未知名字直接报错：预设名拼错必须**早失败**，不能静默当成默认值
    （否则用户以为换了口径、实际没换，最坏）。
    """
    if not preset or preset == "repo":
        return {}
    table = REFERENCE_PRESETS.get(preset)
    if table is None:
        raise KeyError(
            f"未知参考值预设：{preset}；可用：{list(REFERENCE_PRESETS)}")
    return table


def _entry(name, preset=DEFAULT_REF_PRESET):
    """
    取一个标样的参考值条目（dict），可被预设覆盖。

    顺序：预设表里若有该标样 → 用它；否则回落到 STANDARDS。
    （预设只改它自己列出的标样，别的标样不受影响。）
    """
    table = _preset_table(preset)
    if name in table:
        return table[name]
    if name not in STANDARDS:
        raise KeyError(f"未知标准物质：{name}；可用：{list(STANDARDS)}")
    return STANDARDS[name]


def _std_ratio(name, key, preset=DEFAULT_REF_PRESET):
    """
    取某个标准物质的某一个参考比值。

    参数
    ----
    name   : 标准物质内部名，取值见 constants.STANDARDS 的键
    key    : "R68" 或 "R76"
    preset : 参考值预设名（见 constants.REFERENCE_PRESETS）；默认 "repo"

    规则
    ----
    1) 条目里直接给了比值 → 直接返回（精度最高，优先）；
    2) 没给 → 用该条目的公认年龄经衰变方程反算。
       ⚠ 反算隐含假设"该锆石 U-Pb 体系完全封闭"，对好标样这是合理的。
    """
    s = _entry(name, preset)
    if s[key] is not None:
        # 情形 1：文献直给
        return float(s[key])
    # 情形 2：由年龄反算
    return float(r68_of_age(s["age_Ma"])) if key == "R68" else float(r76_of_age(s["age_Ma"]))


def std_ref(name, preset=DEFAULT_REF_PRESET):
    """
    取某个标准物质的参考比值对。

    返回
    ----
    (R68_ref, R76_ref)
        R68_ref : 参考的 206Pb/238U
        R76_ref : 参考的 207Pb/206Pb

    这两个值是"外标归一化"的分母基准，校正因子定义为 F = 参考值 / 实测值。
    `preset` 见 `_std_ratio`；不传则用 constants.DEFAULT_REF_PRESET。
    """
    return _std_ratio(name, "R68", preset), _std_ratio(name, "R76", preset)


def std_age(name, preset=DEFAULT_REF_PRESET):
    """取标准物质的公认年龄 (Ma)；标样不存在则返回 nan。

    未知的 `preset` 名会**报错**（不吞成 nan）—— 与 `std_ref` 保持一致。
    对存在但未被该预设覆盖的标样，仍返回 STANDARDS 里的年龄。
    """
    table = _preset_table(preset)
    s = table.get(name) or STANDARDS.get(name)
    return float(s["age_Ma"]) if s else float("nan")


def std_alias(sample_name):
    """
    把序列文件里写的样品名映射到 STANDARDS 的键。

    序列（LIST.xls）里通常写简称："Ple"、"91500"、"GJ1"，
    而 STANDARDS 里用的是规范名："Plesovice"、"91500"、"GJ1"。
    """
    return {
        "91500": "91500",
        "ple": "Plesovice",
        "plesovice": "Plesovice",
        "pl": "Plesovice",
        "gj1": "GJ1",
        "gj-1": "GJ1",
        "temora": "Temora1",
        "temora1": "Temora1",
        "mudtank": "MudTank",
    }.get(str(sample_name).strip().lower(), None)
