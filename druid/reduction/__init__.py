"""
druid.reduction —— 单点还原层

Tra 是唯一的数据结构契约；装载流程与算法分离，
使得"换仪器格式"和"改算法"两件事互不干扰。
"""
from .trace import Tra, load_spot, window_mask
from .ablation import find_ablation, blank_mean, subtract_blank
from .ratios import reduce_interval

__all__ = [
    "Tra", "load_spot", "window_mask",
    "find_ablation", "blank_mean", "subtract_blank",
    "reduce_interval",
]
