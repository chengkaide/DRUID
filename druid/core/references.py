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
    · 用到 ratio 的地方统一调用 std_ref()，永远是同一套规则；
    · 以后加"按实验室自定义值覆盖"也只改这一个文件。
"""
from __future__ import annotations

from .constants import STANDARDS
from .geochronology import r68_of_age, r76_of_age


def _std_ratio(name, key):
    """
    取某个标准物质的某一个参考比值。

    参数
    ----
    name : 标准物质内部名，取值见 constants.STANDARDS 的键
    key  : "R68" 或 "R76"

    规则
    ----
    1) 文献直接给了比值 → 直接返回（精度最高，优先）；
    2) 没给 → 用公认年龄经衰变方程反算。
       ⚠ 反算隐含假设"该锆石 U-Pb 体系完全封闭"，对好标样这是合理的。
    """
    if name not in STANDARDS:
        raise KeyError(f"未知标准物质：{name}；可用：{list(STANDARDS)}")
    s = STANDARDS[name]
    if s[key] is not None:
        # 情形 1：文献直给
        return float(s[key])
    # 情形 2：由年龄反算
    return float(r68_of_age(s["age_Ma"])) if key == "R68" else float(r76_of_age(s["age_Ma"]))


def std_ref(name):
    """
    取某个标准物质的参考比值对。

    返回
    ----
    (R68_ref, R76_ref)
        R68_ref : 参考的 206Pb/238U
        R76_ref : 参考的 207Pb/206Pb

    这两个值是"外标归一化"的分母基准，校正因子定义为 F = 参考值 / 实测值。
    """
    return _std_ratio(name, "R68"), _std_ratio(name, "R76")


def std_age(name):
    """取标准物质的公认年龄 (Ma)；不存在则返回 nan。"""
    s = STANDARDS.get(name)
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
