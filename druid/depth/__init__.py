"""
druid.depth —— 深度剖面层
"""
from .windows import window_edges, window_profile, window_sums
from .fractionation import (
    bracket_F,
    profile_ages,
    profile_ages_76,
    compare_fractionation,
)
from .domains import (
    segment,
    refine_domains,
    merge_close,
    summarize_segments,
)
from .figures import (
    make_depth_figure,
    save_depth_figure,
)

__all__ = [
    "window_edges", "window_profile", "window_sums",
    "bracket_F", "profile_ages", "profile_ages_76", "compare_fractionation",
    "segment", "refine_domains", "merge_close", "summarize_segments",
    "make_depth_figure", "save_depth_figure",
]
