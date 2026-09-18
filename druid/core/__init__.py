"""
druid.core —— 纯计算内核（无 I/O）

为方便上层书写，这里把最常用的符号一次性导出，
上层可以写 `from druid.core import age68, weighted_mean` 而不必记住具体子模块。
"""
from .constants import (
    HG202_204,
    L232,
    L235,
    L238,
    U238_U235,
    MASSES_NEEDED,
    SUM_MASSES,
    STANDARDS,
    DOMAIN_SPAN_COLORS,
    CJK_FONTS,
    ROLE_GLASS,
    ROLE_LABEL_CN,
    ROLE_PRIMARY,
    ROLE_SECONDARY,
    ROLE_UNKNOWN,
)
from .geochronology import age68, age75, age76, r68_of_age, r75_from, r76_of_age
from .common_lead import stacey_kramers
from .references import std_age, std_alias, std_ref
from .statistics import external_scatter, robust_mask, weighted_mean
from .deadtime import correct_channels, correct_rate

__all__ = [
    # 常量
    "L238", "L235", "L232", "U238_U235", "HG202_204",
    "MASSES_NEEDED", "SUM_MASSES", "STANDARDS",
    "DOMAIN_SPAN_COLORS", "CJK_FONTS",
    "ROLE_PRIMARY", "ROLE_SECONDARY", "ROLE_UNKNOWN", "ROLE_GLASS", "ROLE_LABEL_CN",
    # 年代学
    "r68_of_age", "r76_of_age", "age68", "age75", "age76", "r75_from",
    # 普通铅
    "stacey_kramers",
    # 参考值
    "std_ref", "std_age", "std_alias",
    # 统计
    "weighted_mean", "robust_mask", "external_scatter",
    # 死时间（可选）
    "correct_rate", "correct_channels",
]
