"""
druid.qc —— 批处理质控判定：把"能不能用"变成可机读的检查项
==========================================================

为什么要有这一层
----------------
`run_batch()` 已经把每个数都算出来了，但**判断**一直只存在于三处：控制台上
打印的那几行、Excel 里给人看的表、以及分析者脑子里。三者都无法被程序消费。

下游要接的是"锆石数据质控与解释"这类自动化流程，它需要的是**结论 + 理由 +
数值**这样的三元组，而不是一段文字。所以这里把口径固定下来：

    assess_batch(result) -> [Check(key, level, title, observed, criterion, detail, data)]

设计约束
--------
1. **纯函数、无 I/O、不打印。** 进 DataFrames 与 info 字典，出 list[Check]。
   这样它可以被单测覆盖（`tests/test_qc.py`），也可以在没写盘的情况下调用。
   连给人看的排版也只到「返回字符串」为止（`summary_lines`），
   **打印这个动作留给调用方** —— 命令行和网页界面要往不同的地方写。
2. **只读取，不重算。** 所有数值都取自 `run_batch()` 的产出（results / qc /
   windows / info），这里绝不重新推导年龄。质检层一旦开始算数，就会出现
   "报告里的数和表里的数不一样"这种最难查的问题。
3. **阈值集中在一个 frozen dataclass**（`QCThresholds`），并且可以被
   `BatchConfig.qc_thresholds` 覆盖 —— AGENTS.md 的守则：魔法数字不许散落。
4. **每条检查都有稳定 key。** 下游按 key 分支，不按 title 匹配文字
   （文字会改，key 不会）。改 key 等于改对外契约，要当成破坏性变更对待。
5. **本文件全是中文，写字符串时不要用 ASCII 直引号。** 在中文散文里顺手敲一个
   `"`，会把字符串提前截断，症状是 `SyntaxError: invalid syntax. Perhaps you
   forgot a comma?` 指着一行看起来完全正常的代码（本文件第一次提交前就这么挂了）。
   中文引文一律用「」或者全角引号，英文才用 `\"`。

`level` 取四个值，含义固定
--------------------------
    fail  结论不成立。例如没有主标 → 年龄不可用于定年。
    warn  结论成立但附带前提，必须向下游交代。例如 σext 是假设值而非实测。
    info  陈述事实、不含判断。例如"48 个样品测点里 33 个检出多年龄域"。
    pass  显式通过。**不要因为"没问题"就省掉一条检查** —— 下游需要把
          "查过了、没问题"和"没查"区分开。

比数量更重要的是**检查项的挑选口径**：只保留"会让下游做出不同决定"的那些。
把一条检查加进来，成本不只是这十几行代码，还有以后每次都要读它。
"""
from __future__ import annotations

import math
import re
import unicodedata
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .core.constants import ROLE_LABEL_CN, ROLE_PRIMARY, ROLE_SECONDARY, ROLE_UNKNOWN
from .core.geochronology import age68
from .core.constants import DEFAULT_REF_PRESET
from .core.references import std_age, std_ref
from .core.statistics import EXTERNAL_SCATTER_HI, EXTERNAL_SCATTER_LO

# ── 等级与排序 ──
FAIL, WARN, INFO, PASS = "fail", "warn", "info", "pass"
#: 越小越严重。`worst()` 用它排序，报告里按这个顺序排。
LEVEL_ORDER = {FAIL: 0, WARN: 1, INFO: 2, PASS: 3}
LEVEL_LABEL_CN = {FAIL: "不通过", WARN: "有前提", INFO: "说明", PASS: "通过"}

#: 「本批年龄能用来做什么」的声明 id。**是契约**：下游按 id 分支，不要改字面量。
#:
#: 口径分两根轴 —— 「值」（年龄数值本身可不可信）与「不确定度」（误差棒能不能
#: 用来论证一致性）。两者**不共命运**，这是这个表存在的全部理由：
#: 没有监控标样的批次（如桂北），值照样能用，只有误差棒失去依据。
#: 合成一句话报出去，就把这个区别抹平了。
#:
#: 只保留"会让下游做出不同决定"的四项。加一项的成本不只是这里几行，
#: 还有以后每次都要读它。
CLAIM_IDS = {
    "relative_ordering": "批内按年龄排序、划分相对早晚",
    "relative_age_comparison": "同一批、同一口径下的相对年龄对比",
    "absolute_age_value": "把绝对年龄数值对外引用",
    "absolute_age_uncertainty": "用年龄的误差棒论证两个年龄一致",
}


# ═════════════════════════════════════════════════════════════════════════════
# 阈值
# ═════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class QCThresholds:
    """
    质控判定的全部门槛。**改这里等于改对外结论**，所以要改就改一处。

    默认值的出处分三类，注释里逐条写明，免得下一个人以为是随手拍的：

    * 标样偏差/MSWD 一类取行业惯例（偏差按 1%/2% 两级、MSWD 取 2.5）；
    * `min_usable_windows` 与 `old_age_cutoff_ma` **不是本工具的口径，是
      ADEPT 的口径**，照抄过来只为让下游"会不会被丢掉"的判断与 ADEPT 一致；
    * `count_rate_ratio_warn` 与 `f206_warn_pct` 是**本实验室的经验值**，
      换仪器应当重新标定。
    """
    # ── 标样准确度（相对偏差，%）──
    primary_bias_warn_pct: float = 1.0
    primary_bias_fail_pct: float = 2.0
    secondary_bias_warn_pct: float = 3.0
    secondary_bias_fail_pct: float = 5.0

    # ── 标样精度（加权平均的 MSWD）──
    # 2.5：MSWD 落在 1 附近才算"给出的 σ 是诚实的"。远大于 1 说明 σ 偏小
    # （只含计数统计、漏了外部重现性就是最常见的原因）。
    std_mswd_warn: float = 2.5

    # ── 参考值自洽性 ──
    # 标样若同时给了文献比值与公认年龄，两者应当互相吻合到远好于分析精度。
    # 0.02% 是"比任何一台仪器的重现性都小一个量级"的量级。
    reference_inconsistency_warn_pct: float = 0.02

    # ── 主标稳健剔除 ──
    rejected_primary_frac_warn: float = 0.25

    # ── 样品侧 ──
    concordance_lo: float = 90.0
    concordance_hi: float = 110.0
    concordance_warn_frac: float = 0.85      # 落进区间的占比低于此值 → warn
    concordance_fail_frac: float = 0.70      # 低于此值 → fail
    f206_warn_pct: float = 2.0               # 普通铅占比高于此值算"高"
    f206_warn_frac: float = 0.10             # 高于该值的测点占比超过此值 → warn

    # ── 不分域（整段）口径 ──
    # "分域后判为均一、但整段的 MSWD 已经拒绝了常数年龄模型"的测点占比。
    # 实测示例批次 10/15（67%）—— 这一类在「深度剖面域」表里**没有行**
    # （那张表只收多域点），所以不看「不分域年龄」表就永远发现不了。
    # 0.5：超过一半就值得停下来看一眼，而不是继续往下交。
    unresolved_structure_frac_warn: float = 0.50

    # ── 计数率匹配（基体效应）──
    # 5 倍：单点线性归一化在这个量级上还过得去，再大就会留下明显的非线性残差。
    count_rate_ratio_warn: float = 5.0

    # ── ADEPT 交接（照抄 ADEPT 的判据，不要按自己的喜好改）──
    #: ADEPT 判坪要求"可用年龄窗口 ≥ 10 个"，否则 ok=FALSE、静默丢弃。
    min_usable_windows: int = 10
    #: ADEPT 在 Age68 ≥ 1000 Ma 时会退回 207Pb/206Pb；DRUID 的 `剖面窗口`
    #: 表不带 207/206 列，所以这些窗口对 ADEPT 而言是不可用的。
    old_age_cutoff_ma: float = 1000.0

    # ── 对照方法差异 ──
    method_diff_warn_pct: float = 5.0        # simple 与 ftau 的相对差中位


# ═════════════════════════════════════════════════════════════════════════════
# 检查项
# ═════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class Check:
    """
    一条质控结论。

    字段分工（别混用，下游要按分工来的）：
        key        稳定标识，供程序分支。**改它等于改对外契约。**
        level      fail / warn / info / pass
        title      一句话结论，给人看
        observed   实测值（写成字符串，因为经常是"41/48"这种）
        criterion  判据（写成字符串，因为经常是区间）
        detail     为什么重要 / 该怎么办；只有在**确有需要**时才写
        data       机读的原始数值。下游要算东西就读这里，不要解析 observed。
    """
    key: str
    level: str
    title: str
    observed: str = ""
    criterion: str = ""
    detail: str = ""
    data: Dict[str, object] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, object]:
        """转成 `json.dumps` 能直接吃的字典。"""
        d = asdict(self)
        d["data"] = jsonable(d["data"])
        return d

    def line(self) -> str:
        """一行式摘要，供 CLI 打印。"""
        return (f"[{LEVEL_LABEL_CN.get(self.level, self.level)}] {self.title}"
                + (f" —— 实测 {self.observed}" if self.observed else ""))


def _mk(key, level, title, observed="", criterion="", detail="", **data) -> Check:
    return Check(key=key, level=level, title=title, observed=observed,
                 criterion=criterion, detail=detail, data=data)


def jsonable(obj):
    """
    把 numpy / pandas 的东西转成 `json.dumps` 能接受的原生类型。

    两个必须处理的坑，都真的会踩到：

    1. **NaN / Inf 不是合法 JSON。** `json.dumps(nan)` 默认写出 `NaN`，
       那不是一个合法 JSON 字面量 —— Python 自己能读回来，别的语言的解析器
       一律报错。写契约文件时要用 `allow_nan=False`，于是 NaN 必须变 `None`。
    2. **numpy 标量不可序列化。** `np.int64` 不是 `int`（`np.bool_` 也不是
       `bool`），`json.dumps` 直接抛 TypeError。而且这个错**只在真正有数据的
       路径上出现**，拿小样本手测很容易漏掉。
    """
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if obj is pd.NA:                         # pandas 的缺失值单例
        return None
    if isinstance(obj, np.generic):          # np.int64 / np.bool_ / np.float64 …
        obj = obj.item()                     # 先降成 Python 原生，再走下面
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [jsonable(v) for v in obj]
    return str(obj)


def worst(checks: Sequence[Check]) -> str:
    """最严重的等级；空列表返回 pass。"""
    if not checks:
        return PASS
    return min((c.level for c in checks), key=lambda lv: LEVEL_ORDER.get(lv, 99))


def counts(checks: Sequence[Check]) -> Dict[str, int]:
    """各等级的条数，四种等级都保证有键（下游不用 .get）。"""
    out = {FAIL: 0, WARN: 0, INFO: 0, PASS: 0}
    for c in checks:
        if c.level in out:
            out[c.level] += 1
    return out


def by_key(checks: Sequence[Check]) -> Dict[str, Check]:
    """{key: Check}，便于下游直接取某一条。"""
    return {c.key: c for c in checks}


def verdict(checks: Sequence[Check]) -> Dict[str, object]:
    """
    汇总成 `{level, counts, headline}` —— 交接 JSON 与网页界面共用这一份。

    为什么不各写一遍：这三样是**同一件事的三个视图**，分别实现在两处时，
    最容易出现的是 headline 说「已校准、可用」而 level 是 fail。
    下游拿 JSON、用户看网页，谁对不上都会让人不再信任这套判据。

    `headline` 把 fail/warn 的标题用「；」连起来，全过时给一句肯定的话 ——
    **不是空字符串**，因为前端拿去显示时空串会渲染成一个没有内容的提示条。
    """
    bad = [c for c in sorted(checks,
                             key=lambda c: (LEVEL_ORDER.get(c.level, 99), c.key))
           if c.level in (FAIL, WARN)]
    return {
        "level": worst(checks),
        "counts": counts(checks),
        "headline": ("全部质控项通过" if not bad
                     else "；".join(c.title for c in bad)),
    }


# ═════════════════════════════════════════════════════════════════════════════
# 排版（仍然不打印：返回字符串列表，谁调用谁负责输出）
# ═════════════════════════════════════════════════════════════════════════════
def _is_wide(ch: str) -> bool:
    """是不是全角字符（中日韩、全角标点）。终端里占 2 列。"""
    return unicodedata.east_asian_width(ch) in ("W", "F")


def _display_width(s: str) -> int:
    """终端显示宽度：全角占 2 列、其余占 1 列。

    按 `len()` 折行在中文段落上会折错 —— 一行 60 个汉字实际占 120 列，
    在 80 列的终端里会被自动再折一次，折出来的位置跟缩进完全对不上。
    """
    return sum(2 if _is_wide(ch) else 1 for ch in s)


def _tokens(text: str) -> List[str]:
    """
    切成折行的最小单位：**连续的半角字符算一个词元，每个汉字各算一个**。

    为什么不能逐字符折行：`Age76_1s`、`0.000%`、`EX2022A_01`
    这些都是不可分割的记号，从中间断开会让人认不出来（`Age7` / `6_1s`）。
    连续空白折成一个空格，避免输出右边缘参差。
    """
    out, buf = [], ""
    for ch in text:
        if ch.isspace():
            if buf:
                out.append(buf)
                buf = ""
            if out and out[-1] != " ":
                out.append(" ")
        elif _is_wide(ch):
            if buf:
                out.append(buf)
                buf = ""
            out.append(ch)
        else:
            buf += ch
    if buf:
        out.append(buf)
    return out


#: 不能出现在行首的字符（中文排版的「禁则」）。落到行首会让人以为句子在那里结束。
_NO_LINE_START = set("，。、；：？！）〕】》」』”’…‰℃%°·/")
#: 不能出现在行尾的字符（开括号/开引号留在行尾会跟后面的内容断开）。
_NO_LINE_END = set("（〔【《「『“‘")


def _wrap_cjk(text: str, width: int = 66, indent: str = "",
              indent_rest: Optional[str] = None) -> List[str]:
    """
    按**显示宽度**折行，折点是词元边界，并做中英文混排的禁则处理。

    `indent` 给首行、`indent_rest` 给后续行（不给则与首行相同）。
    两者分开是为了悬挂缩进：折下来的行往里缩两格，一眼能看出是同一句的延续，
    而不是下一条检查项。

    单词元超宽时（例如一条 URL）才会被硬切开 —— 硬切是最后手段，
    因为切点会落在记号中间。
    """
    if indent_rest is None:
        indent_rest = indent
    if not text:
        return []

    lines: List[str] = []
    cur, pre = "", indent
    for tok in _tokens(text):
        cand = cur + tok
        # 禁则允许的「故意超宽」：这一行会比 width 多占一两格。
        # 必须记下来 —— 否则下面的硬切会在同一次循环里把它切回去，
        # 看上去像禁则完全没生效（本函数第一版就是这么错的）。
        overflow_ok = False
        if cur and _display_width(cand) > width:
            if tok[:1] in _NO_LINE_START:
                # 标点宁可跟着上一行多占一格，也不落到行首
                cur, overflow_ok = cand, True
            else:
                while cur and cur[-1] in _NO_LINE_END:
                    tok = cur[-1] + tok
                    cur = cur[:-1]
                lines.append(pre + cur.rstrip())
                cur, pre = tok.lstrip(), indent_rest
        else:
            cur = cand
        # 单词元就超宽 —— 只能硬切，切完再把剩下的当新行继续
        while not overflow_ok and _display_width(cur) > width:
            cut = ""
            for ch in cur:
                if cut and _display_width(cut) + _display_width(ch) > width:
                    break
                cut += ch
            lines.append(pre + cut.rstrip())
            cur, pre = cur[len(cut):].lstrip(), indent_rest
    if cur:
        lines.append(pre + cur.rstrip())
    return lines


def summary_lines(checks: Sequence[Check], *, max_items: int = 8,
                  indent: str = "    ",
                  detail_ref: Optional[str] = "handoff.json") -> List[str]:
    """
    把检查项排成给人看的几行文字。**只排版，不打印** —— 输出由调用方决定。

    为什么放在 qc.py 而不是 CLI：命令行与网页界面都要打印同一套结论，各写一遍
    必然慢慢跑偏（`io.report.age68_column` 那次分家就是这么来的：同一个意思的
    片段在两处各有一份）。

    `detail_ref` 是「完整清单在哪儿」的文件名。**没写盘时必须传 None** ——
    否则会指着 handoff.json 说"见这里"，而用户当时用的是 `--no-json`，
    那个文件根本不存在。

    排版规则（每条都有理由，不是为了好看）
    --------------------------------------
    · 首行给汇总：最严重等级 + 各等级条数。**只读这一行也知道能不能用。**
    · `fail` / `warn` 逐条列出 —— 这些是「会让下游做出不同决定」的条目，
      带上判据，必要时带处置建议。
    · `info` / `pass` 只报条数、不逐条列。它们的作用是证明「查过了」，
      逐条铺满控制台只会把 fail 挤到屏幕外面去。
    · `max_items` 是每个等级的列出上限，超出只报「还有 N 条」；
      设上限而不是无脑全列，因为一个坏批次能有几十条同类 warn。
    """
    if not checks:
        return [indent + "质控：没有收到任何检查项（assess_batch 没跑？）"]

    n = counts(checks)
    ordered = sorted(checks,
                     key=lambda c: (LEVEL_ORDER.get(c.level, 99), c.key))
    out = [indent + f"质控结论：{LEVEL_LABEL_CN.get(worst(checks), '?')}"
                   f"（{n[FAIL]} 项不通过 / {n[WARN]} 项有前提 / "
                   f"{n[INFO]} 项说明 / {n[PASS]} 项通过）"]
    where = f"，见 {detail_ref}" if detail_ref else ""

    for lv in (FAIL, WARN):
        group = [c for c in ordered if c.level == lv]
        for c in group[:max_items]:
            head = f"[{LEVEL_LABEL_CN[lv]}] {c.title}"
            if c.observed:
                head += f" —— 实测 {c.observed}"
            out.extend(_wrap_cjk(head, indent=indent, indent_rest=indent + "  "))
            if c.criterion:
                out.extend(_wrap_cjk(f"判据：{c.criterion}",
                                     indent=indent + "  ", indent_rest=indent + "    "))
            if c.detail:
                out.extend(_wrap_cjk(f"→ {c.detail}",
                                     indent=indent + "  ", indent_rest=indent + "    "))
        if len(group) > max_items:
            out.append(indent + f"[{LEVEL_LABEL_CN[lv]}] …还有 "
                               f"{len(group) - max_items} 条{where}")

    tail = [f"{n[lv]} 项{LEVEL_LABEL_CN[lv]}" for lv in (INFO, PASS) if n[lv]]
    if tail:
        out.append(indent + "（" + "、".join(tail)
                   + "；未逐条列出" + (f"，完整清单见 {detail_ref}" if detail_ref
                                     else "，本次未写出清单") + "）")
    return out


# ═════════════════════════════════════════════════════════════════════════════
# 取值小工具（一律"读"，不"算"）
# ═════════════════════════════════════════════════════════════════════════════
def _finite(x) -> Optional[float]:
    """转成有限 float，否则 None。质控层到处要用，别每处 `math.isfinite`。"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _col(df, name):
    """取列；不存在返回 None（而不是抛 KeyError —— 缺列本身要被报成一条检查）。"""
    if df is None or getattr(df, "empty", True) or name not in df.columns:
        return None
    return df[name]


def _sub(df, mask):
    """按布尔掩码取子表，掩码长度对不上时返回空表而不是抛异常。"""
    try:
        return df[mask]
    except (IndexError, KeyError, ValueError):
        return df.iloc[0:0]


def _pct(x: float) -> str:
    return f"{x:.2f}%"


#: 样品名里连续出现的分隔符（`X--2-9`、`X -  2` 这类笔误）
_RUNS = re.compile(r"[\s\-_]{1,}")
#: 末尾的测点编号（`-9` / `_12`），去掉它得到"岩样名"
_TAIL_INDEX = re.compile(r"[-_]\d+$")


def _canonical_name(name: str) -> str:
    """
    把样品名收成规范形式：连续分隔符收成一个、去掉首尾分隔符与空白。

    只做**分隔符**层面的归一，不碰字母数字 —— 归一化过头会把真正不同的样品
    合并到一起，那是比原问题更糟的错误。
    """
    s = _RUNS.sub("-", str(name).strip())
    return s.strip("-_ ")


def _rock_of(name: str) -> str:
    """
    取"岩样名"：去掉样品名末尾的测点编号，再收掉尾部的分隔符。

    `X-2-17` → `X-2`；`X--2-9` 去掉末尾 `-9` 后剩 `X-` → `X`。
    **这里刻意不做归一化** —— 归一化的先后顺序正是要考察的东西：
    对原名取岩样得 `{X-2, X}`（同一个岩样被拆成两个），
    先归一化再取岩样才得 `{X-2}`。两者之差就是这个检查要报的量。
    """
    return _TAIL_INDEX.sub("", str(name).strip()).strip("-_ ")


def _window_stats(result, th: QCThresholds):
    """
    逐测点的窗口统计，只用 `剖面窗口` 表（即 ADEPT 会拿到的那张）。

    ⚠ **必须用 `windows` 而不是 `results`。** 这张表就是交接给 ADEPT 的全部内容，
    "ADEPT 会不会把这个测点丢掉"只能从这张表算出来；从 results 算会得出
    另一个数（results 里的年龄是整段值，与窗口级的可用性无关）。

    返回 [{analysis, sample, n_windows, n_usable}]，无窗口表时返回 None。
    """
    w = getattr(result, "windows", None)
    if w is None or getattr(w, "empty", True):
        return None
    if not {"Analysis", "Age68"}.issubset(set(w.columns)):
        return None
    has_sample = "Sample" in w.columns
    rows = []
    for name, g in w.groupby("Analysis", sort=False):
        n = int(len(g))
        n_usable = int((g["Age68"] < th.old_age_cutoff_ma).sum())
        rows.append(dict(
            analysis=str(name),
            sample=str(g["Sample"].iloc[0]) if has_sample else "",
            n_windows=n, n_usable=n_usable))
    return rows


# ═════════════════════════════════════════════════════════════════════════════
# 各类检查
# ═════════════════════════════════════════════════════════════════════════════
def _check_calibration(result, cfg, th, out: List[Check]) -> None:
    """① 能不能给出绝对年龄 —— 这一组不过，后面的年龄就不必看了。"""
    res = result.results
    roles = _col(res, "类型")

    mode_col = _col(res, "校准状态")
    mode = str(mode_col.iloc[0]) if mode_col is not None and len(mode_col) else "未知"
    calibrated = (mode == "已校准")
    out.append(_mk(
        "calibration.mode", PASS if calibrated else FAIL,
        "已做外标归一化" if calibrated else "未校准：比值是实测值，年龄不能用于定年",
        observed=mode,
        criterion="校准状态 = 已校准",
        detail="" if calibrated else
               "批次里没有找到主标样，或主标名没设对。此时 F=1、年龄仅作相对参考。"
               "在序列表里加入主标（如 91500），或把『主标名』设成实际使用的标样后重跑。",
        calibrated=calibrated))

    n_primary = int((roles == ROLE_LABEL_CN[ROLE_PRIMARY]).sum()) if roles is not None else 0
    n_secondary = int((roles == ROLE_LABEL_CN[ROLE_SECONDARY]).sum()) if roles is not None else 0
    n_unknown = int((roles == ROLE_LABEL_CN[ROLE_UNKNOWN]).sum()) if roles is not None else 0

    # 主标太少时归一化因子本身的误差会成为主导项，所以单独报一条。
    if n_primary == 0:
        lv, why = FAIL, "没有主标就无法构造归一化因子 F，绝对年龄无从谈起。"
    elif n_primary < 5:
        lv, why = WARN, ("主标点数偏少，F 因子本身的不确定度会主导结果。"
                         "夹逼归一化需要主标在序列里分布足够密。")
    else:
        lv, why = PASS, ""
    out.append(_mk(
        "calibration.primary_spots", lv, f"主标测点 {n_primary} 个",
        observed=str(n_primary), criterion="≥5（建议 ≥10）", detail=why,
        n_primary=n_primary))

    # 监控标样只影响 σext 与二次校正，缺了不致命 —— 但必须报，因为它把
    # "外部重现性"从实测值变成了假设值（见 uncertainty 那一条）。
    out.append(_mk(
        "calibration.secondary_spots", PASS if n_secondary else WARN,
        f"监控标样测点 {n_secondary} 个" if n_secondary else "没有监控标样",
        observed=str(n_secondary), criterion="≥1",
        detail="" if n_secondary else
               "没有第二标样就无法实测外部重现性，也不能做基体匹配的二次校正。",
        n_secondary=n_secondary))

    # 稳健统计把哪些主标剔掉了 —— 剔除本身正常（异常点该剔），但一次剔掉
    # 四分之一以上说明序列质量有问题。
    rejected = list(result.info.get("rejected_primary") or [])
    frac = len(rejected) / n_primary if n_primary else 0.0
    if not rejected:
        lv = PASS
    elif frac > th.rejected_primary_frac_warn:
        lv = WARN
    else:
        lv = INFO
    out.append(_mk(
        "calibration.primary_rejected", lv,
        f"稳健统计剔除 {len(rejected)} 个主标测点",
        observed=f"{len(rejected)}/{n_primary} = {frac:.0%}",
        criterion=f"剔除比例 ≤ {th.rejected_primary_frac_warn:.0%}",
        detail="剔除的是 204Pb 异常或 Hg 瞬时波动导致的离群点，属正常操作；"
               "比例过高则说明序列本身不稳，F 因子可能被少数点主导。",
        rejected=list(rejected), frac=float(frac)))

    out.append(_mk(
        "samples.present", PASS if n_unknown else FAIL,
        f"样品测点 {n_unknown} 个" if n_unknown else "本批次没有样品，只有标样",
        observed=str(n_unknown), criterion="≥1",
        detail="" if n_unknown else
               "纯标样/校准批次只产出 QC 表，不存在可解释的样品年龄。",
        n_unknown=n_unknown))


def _check_standards(result, cfg, th, out: List[Check]) -> None:
    """② 标样准确度与精度 —— 所有样品年龄的可信度都挂在这里。"""
    qc = getattr(result, "qc", None)

    def _row(name):
        if qc is None or getattr(qc, "empty", True) or "标样" not in qc.columns:
            return None
        g = qc[qc["标样"] == name]
        return g.iloc[0] if len(g) else None

    rows = {}
    for tag, name, wan, fal in (
            ("primary", cfg.primary, th.primary_bias_warn_pct, th.primary_bias_fail_pct),
            ("secondary", cfg.secondary, th.secondary_bias_warn_pct, th.secondary_bias_fail_pct)):
        q = _row(name)
        if q is None:
            out.append(_mk(
                f"standards.{tag}_bias", WARN, f"标样 {name} 没有 QC 结果",
                criterion=f"|偏差| < {wan}%",
                detail="该标样在结果表里没有测点，或全部被剔除。"))
            continue
        rows[tag] = q
        bias = _finite(q.get("偏差_pct"))
        n = int(q.get("点数", 0))
        if bias is None:
            lv, obs = WARN, "参考年龄缺失，无法比较"
        else:
            ab = abs(bias)
            lv = PASS if ab < wan else (WARN if ab < fal else FAIL)
            obs = f"{bias:+.2f}%  (n={n})"
        extra = ""
        if tag == "primary":
            extra = ("这是归一化锚点，它的偏差会**整体平移全部样品年龄**，"
                     "不是随机误差。")
        out.append(_mk(
            f"standards.{tag}_bias", lv, f"标样 {name} 加权平均与参考值之差",
            observed=obs, criterion=f"|偏差| < {wan}%（≥{fal}% 判 fail）",
            detail=extra,
            name=str(name), n=n, bias_pct=bias,
            measured_Ma=_finite(q.get("加权平均年龄_Ma")),
            reference_Ma=_finite(q.get("参考年龄_Ma")),
            s2_Ma=_finite(q.get("s2_Ma"))))

        mswd = _finite(q.get("MSWD"))
        out.append(_mk(
            f"standards.{tag}_mswd", PASS if (mswd is not None and mswd <= th.std_mswd_warn) else WARN,
            f"标样 {name} 加权平均的 MSWD",
            observed=("—" if mswd is None else f"{mswd:.2f}  (n={n})"),
            criterion=f"≤ {th.std_mswd_warn}",
            detail="" if (mswd is None or mswd <= th.std_mswd_warn) else
                   "MSWD 明显大于 1 说明测点散布超过各自给出的 σ —— 通常是 σ 偏小"
                   "（只含计数统计、漏了外部重现性），而不是标样真的不均一。"
                   "这时样品年龄的误差棒同样偏小。",
            name=str(name), mswd=mswd))

    # ── 参考值自洽性：标样同时给了文献比值与公认年龄时，两者必须互相吻合 ──
    # 这一条是本项目真实存在的一个问题，把它变成检查项而不是口头约定：
    # 91500 的 R68=0.17928 反算回年龄是 1063.04 Ma，而库里写的 age_Ma=1062.4，
    # 差 +0.0602%。后果是**全部样品年龄一致偏老 0.0602%**，
    # 而且 QC「偏差%」有一个 +0.0602% 的固定地板。
    # 参考值预设（见 workflow.BatchConfig.ref_preset）：选了非默认口径时，
    # 这条自洽性检查按**所选口径**判 —— 例如选了 horstwood2016，91500 的
    # 比值与年龄本来就自洽，这条就该变绿（这正是"换个口径"的可见后果之一）。
    preset = getattr(cfg, "ref_preset", DEFAULT_REF_PRESET)
    inc = {}
    for name in dict.fromkeys([cfg.primary, cfg.secondary]):
        try:
            r68, _ = std_ref(name, preset)
        except KeyError:
            continue
        a_ref = std_age(name, preset)
        a_from = _finite(age68(r68))
        if not (a_ref and math.isfinite(a_ref) and a_from):
            continue
        dpct = (a_from - a_ref) / a_ref * 100
        inc[name] = dict(R68=float(r68), age_from_R68_Ma=a_from,
                         age_Ma=float(a_ref), diff_pct=float(dpct))
    bad = {k: v for k, v in inc.items() if abs(v["diff_pct"]) > th.reference_inconsistency_warn_pct}
    if bad:
        bits = "; ".join(f"{k} 差 {v['diff_pct']:+.4f}%" for k, v in bad.items())
        out.append(_mk(
            "reference.self_consistency", WARN,
            "参考标样库里比值与年龄不自洽：全体年龄会被系统性平移",
            observed=bits,
            criterion=f"|反算年龄 − 公认年龄| ≤ {th.reference_inconsistency_warn_pct}%",
            detail="标样同时给了文献比值（R68/R76）与公认年龄时，前者反算回年龄应当"
                   "与后者一致。不一致的后果有两个：① 全部样品年龄按同一个因子整体"
                   "偏老/偏新，与文献做千分之几级别的对比时会出错；② 该标样的 QC"
                   "「偏差%」带一个固定地板，看起来像分析质量、实际是口径。"
                   "这是**上游口径未决**，不是本工具算错；要消掉它得先定改哪一头。",
            standards=inc))
    elif inc:
        out.append(_mk(
            "reference.self_consistency", PASS,
            "参考标样库的比值与年龄自洽",
            observed="; ".join(f"{k} 差 {v['diff_pct']:+.4f}%" for k, v in inc.items()),
            criterion=f"≤ {th.reference_inconsistency_warn_pct}%",
            standards=inc))

    # ── 计数率匹配：主标与样品的信号强度差一个数量级时，单一线性因子消不掉
    #    非线性残差，这正是"QC 二次校正"存在的原因。
    u = _col(result.results, "U238_cps")
    roles = _col(result.results, "类型")
    ratio = None
    u_prim_med = u_unk_med = None
    if u is not None and roles is not None:
        up = u[roles == ROLE_LABEL_CN[ROLE_PRIMARY]]
        uu = u[roles == ROLE_LABEL_CN[ROLE_UNKNOWN]]
        if len(up) and len(uu):
            u_prim_med = _finite(up.median())
            u_unk_med = _finite(uu.median())
            if u_prim_med and u_unk_med and u_prim_med > 0:
                ratio = u_unk_med / u_prim_med
    if ratio is None:
        lv = INFO
    elif ratio > th.count_rate_ratio_warn or ratio < 1.0 / th.count_rate_ratio_warn:
        lv = WARN
    else:
        lv = PASS
    out.append(_mk(
        "uncertainty.count_rate_match", lv,
        "标样与样品的 U 信号强度之比",
        observed=("—" if ratio is None else f"{ratio:.1f}×"),
        criterion=f"1/{th.count_rate_ratio_warn:g} ~ {th.count_rate_ratio_warn:g}",
        detail="" if lv == PASS else
               "主标与样品计数率差得越远，探测器/电子学的非线性残差越难被一个线性"
               "归一化因子消掉，表现为监控标样系统性偏老。基体匹配的二次校正只是"
               "治标；治本是缩小信号强度差距（换小束斑/降低能量）或施加可靠的死时间校正。",
        ratio=(None if ratio is None else float(ratio)),
        primary_median_cps=u_prim_med, unknown_median_cps=u_unk_med))


def _check_uncertainty(result, cfg, th, out: List[Check]) -> None:
    """③ 不确定度从哪来 —— 是测出来的还是假设的，必须说清楚。"""
    res = result.results
    roles = _col(res, "类型")
    n_secondary = int((roles == ROLE_LABEL_CN[ROLE_SECONDARY]).sum()) if roles is not None else 0
    sd68 = _finite(result.info.get("sd68"))
    sd76 = _finite(result.info.get("sd76"))

    forced = (cfg.sigma_ext68 is not None or cfg.sigma_ext76 is not None)
    if forced:
        source, lv = "forced", WARN
        detail = ("外部重现性由调用方在参数里指定（--sigma-ext68 / --sigma-ext76），"
                  "不是从本批数据测出来的 —— 换批次时这个值不会跟着变，"
                  "对比不同批次时要留意。")
    elif n_secondary:
        source, lv = "measured", PASS
        detail = ""
    else:
        source, lv = "assumed", WARN
        detail = ("本批次没有监控标样，外部重现性退回到写死的经验值"
                  "（206/238 用 0.70%、207/206 用 0.25%）。也就是说**样品年龄的"
                  "误差棒里有一块是假设而非实测的**，不能用它去论证"
                  "「两个年龄在误差内一致」。补做监控标样是唯一的根治办法。")

    out.append(_mk(
        "uncertainty.sigma_ext_source", lv,
        f"外部重现性来源：{'实测' if source == 'measured' else ('参数指定' if source == 'forced' else '写死的假设值')}",
        observed=f"206/238 {sd68:.2%}、207/206 {sd76:.2%}  (n_监控={n_secondary})"
                 if (sd68 is not None and sd76 is not None) else "—",
        criterion="优先用监控标样实测",
        detail=detail,
        source=source, n_secondary=n_secondary, sd68=sd68, sd76=sd76))

    # ── σext 被上下限截断：这时它表示的**不是一个测出来的量，而是一个缺口** ──
    # `external_scatter()` 在监控标样不足 2 个时返回下限、在散度算爆时返回上限。
    # 两者都会让下游以为"外部重现性已经考虑过了"，实际是没测出来。
    bound = []
    for label, v in (("206/238", sd68), ("207/206", sd76)):
        if v is None:
            continue
        if abs(v - EXTERNAL_SCATTER_LO) < 1e-12:
            bound.append(f"{label} 触到下限 {EXTERNAL_SCATTER_LO:.2%}")
        elif abs(v - EXTERNAL_SCATTER_HI) < 1e-12:
            bound.append(f"{label} 触到上限 {EXTERNAL_SCATTER_HI:.2%}")
    out.append(_mk(
        "uncertainty.sigma_ext_bound", WARN if bound else PASS,
        "外部重现性触到了上下限保护" if bound else "外部重现性未触边界",
        observed="; ".join(bound) if bound else "—",
        criterion=f"不取到 {EXTERNAL_SCATTER_LO:.2%} / {EXTERNAL_SCATTER_HI:.2%} 两端",
        detail="" if not bound else
               "触到边界说明这个数**不是测出来的**：下限来自"
               f"「监控标样不足 2 个」（函数直接返回 {EXTERNAL_SCATTER_LO:.2%}）"
               "或实测散度小于内部精度（方差相减后截到 0）；"
               f"上限来自散度算爆被截在 {EXTERNAL_SCATTER_HI:.2%}。"
               "用这种 σext 算出的误差棒只能当量级参考，不能拿来做"
               "「两个年龄在误差内一致」这类判断。",
        bound=bound, sd68=sd68, sd76=sd76,
        lo=EXTERNAL_SCATTER_LO, hi=EXTERNAL_SCATTER_HI))

    # 二次校正被施加过就一定要说 —— 它是一个乘性平移，会改变每个年龄的绝对值。
    corr = result.info.get("secondary_correction") or None
    if corr:
        kfac = _finite(corr.get("factor"))
        # 系数 <1 时年龄是被**乘小**的，所以带符号写成负数 —— 用 abs() 会写成
        # "+1.78%"，读者会以为年龄被抬高了。
        shift = (kfac - 1) * 100 if kfac else float("nan")
        out.append(_mk(
            "uncertainty.secondary_correction", WARN,
            "已施加监控标样二次校正（基体匹配）",
            observed=(f"系数 {kfac:.5f}，全部年龄平移 {shift:+.2f}%" if kfac else "—"),
            criterion="—（这是补救措施，不是常规步骤）",
            detail="二次校正用计数率与样品相当的监控标样做再校准，是行业通行的治标"
                   "办法。它把全部年龄乘了同一个系数，所以：① 与未做该校正的结果"
                   "**不可直接比数值**；② 它并没有消除非线性残差，只是把它压小。"
                   "报告里必须交代这一点。",
            factor=kfac, sample=str(corr.get("sample", "")),
            measured_Ma=_finite(corr.get("measured")),
            reference_Ma=_finite(corr.get("ref_age")), n=int(corr.get("n", 0) or 0)))
    else:
        out.append(_mk(
            "uncertainty.secondary_correction", PASS,
            "未施加二次校正",
            criterion="—", detail="", factor=None))


def _check_samples(result, cfg, th, out: List[Check]) -> None:
    """④ 样品侧的可解释前提：协和度、普通铅、样品名卫生。"""
    res = result.results
    roles = _col(res, "类型")
    if roles is None:
        return
    unk = _sub(res, roles == ROLE_LABEL_CN[ROLE_UNKNOWN])
    if getattr(unk, "empty", True):
        return

    conc = _col(unk, "协和度_pct")
    if conc is not None:
        c = conc.dropna()
        frac = _finite((c.between(th.concordance_lo, th.concordance_hi)).mean()) or 0.0
        med = _finite(c.median())
        lv = (PASS if frac >= th.concordance_warn_frac
              else (WARN if frac >= th.concordance_fail_frac else FAIL))
        out.append(_mk(
            "samples.concordance", lv,
            "样品协和度落在 %g–%g%% 的占比" % (th.concordance_lo, th.concordance_hi),
            observed=f"{frac:.0%}  (中位 {med:.1f}%，n={len(c)})" if med is not None else f"{frac:.0%}",
            criterion=f"≥ {th.concordance_warn_frac:.0%}（< {th.concordance_fail_frac:.0%} 判 fail）",
            detail="" if lv == PASS else
                   "协和度偏离 100% 说明 206/238 与 207/235 两个体系不自洽，"
                   "最常见的原因是普通铅未校正（偏老）或 Pb 丢失（偏新）。"
                   "这类测点的 206/238 年龄不能单独当定年结果用。",
            frac=float(frac), median_pct=med, n=int(len(c))))

    f206 = _col(unk, "f206_pct")
    if f206 is not None:
        f = f206.dropna()
        n_hi = int((f > th.f206_warn_pct).sum())
        frac = n_hi / len(f) if len(f) else 0.0
        lv = WARN if frac > th.f206_warn_frac else PASS
        out.append(_mk(
            "samples.common_lead", lv,
            f"普通铅占比高于 {th.f206_warn_pct:g}% 的样品测点",
            observed=f"{n_hi}/{len(f)} = {frac:.0%}  (中位 {_finite(f.median()):.3f}%)",
            criterion=f"占比 ≤ {th.f206_warn_frac:.0%}",
            detail="" if lv == PASS else
                   "f206（普通铅在总 206Pb 里的份额）偏高说明普通铅校正触发得多。"
                   "204Pb 只有几~几十 cps、而 204Hg 本底几百 cps，"
                   "校正的残余本身带噪声，会同时抬高年龄与它的不确定度。",
            n_high=n_hi, frac=float(frac), median_pct=_finite(f.median())))

    # ── 样品名卫生：不干净的名字会让"按样品分组"悄悄出错 ──
    # 真实批次里出现过：同一批里 8 个测点写成 `X--2-9`…`-16`，其余写成 `X-2-17`…。
    # 这些名字**不会撞名**，所以"查重名"抓不到它们；真正的后果是**归组层级错位** ——
    # 去掉末尾测点编号后，`X--2-9` 属于 `X`、`X-2-17` 属于 `X-2`，
    # 同一个岩样被拆成了两组。所以这里同时报"非规范写法"和"归组个数差"。
    names = [str(x) for x in _col(unk, "样品").tolist()] if _col(unk, "样品") is not None else []
    blank = [n for n in names if not n.strip()]
    spaced = [n for n in names if n != n.strip()]
    canon = {n: _canonical_name(n) for n in names}
    variant = [n for n in names if n and n != canon[n]]
    grouped = {}
    for n in names:
        grouped.setdefault(canon[n], []).append(n)
    collide = {k: sorted(set(v)) for k, v in grouped.items()
               if len(v) > 1 and len(set(v)) > 1}
    rocks_raw = {_rock_of(n) for n in names if n}
    rocks_canon = {_rock_of(canon[n]) for n in names if n}
    n_split = len(rocks_raw) - len(rocks_canon)
    n_problem = len(blank) + len(spaced) + len(variant)
    out.append(_mk(
        "data.name_hygiene", PASS if n_problem == 0 and n_split == 0 else WARN,
        "样品名规范" if n_problem == 0 and n_split == 0 else
        "样品名有非规范写法，按名归组会把同一个岩样拆成两组",
        observed=(f"空名 {len(blank)}、首尾空白 {len(spaced)}、非规范写法 {len(variant)}、"
                  f"归组个数差 {n_split}（{len(rocks_raw)} → {len(rocks_canon)}）"),
        criterion="全部为 0",
        detail="" if (n_problem == 0 and n_split == 0) else
               "典型写法是连续两个分隔符（`X--2-9`）。它**不会和别人撞名**，"
               "所以查重名抓不到；问题是去掉末尾测点编号后，它归到 `X`，"
               "而同批其余测点归到 `X-2` —— 同一个岩样被算成两个。"
               "「归组个数差」就是归一化前后岩样个数的差，非 0 说明这件事真的发生了。"
               "修法：归组前先把连续分隔符收成一个、去掉首尾分隔符；"
               "或者更省事——直接修序列表。",
        n_blank=len(blank), n_spaced=len(spaced), n_variant=len(variant),
        n_rocks_raw=len(rocks_raw), n_rocks_canonical=len(rocks_canon),
        n_group_split=int(n_split),
        blank=blank, spaced=spaced, variant=variant,
        collide={k: v for k, v in collide.items()}))


def _check_whole_spot(result, cfg, th, out: List[Check]) -> None:
    """
    **不分域（整段）口径**的两条检查项。

    为什么值得单独一组
    ------------------
    「深度剖面域」表只收多域点 —— 均一点在里面**没有行**。所以只看那张表，
    "分域之后仍然只报均一、而整段的 MSWD 其实已经拒绝了常数年龄模型"这类
    测点是看不见的。`overall`（不分域年龄表）把**全部**测点按同一口径列出来，
    下面两件事才变得可见：

        ① 有多少测点的整段年龄可以直接用（与常数年龄模型相容）；
        ② 有多少测点在分域之后仍带着未解决的结构 —— 它们给出的"年龄"
           要么是几个年龄的加权混合值，要么是被残余倾斜拉出来的数。

    ⚠ 本组是**只读**的：所有数值都取自 depth 层已经算好的 `overall` 表，
    这里一个数都不重算。质检层一旦开始自己算年龄，就会出现"报告里的数
    和表里的数不一样"，而且没人知道该信哪个。
    """
    ov = getattr(result, "overall", None)
    if ov is None or getattr(ov, "empty", True):
        # 缺数据既不能报 pass 也不能报 0 —— 两种都会被下游读成"没问题"。
        out.append(_mk(
            "samples.whole_spot", WARN, "不分域整段年龄：无数据",
            observed="结果对象里没有「不分域年龄」表",
            criterion="该表由 depth 层产出；缺失时本组检查项无法进行",
            detail="若用的是旧版结果对象、或本次运行跳过了深度分析，这张表"
                   "本来就不存在。**不要把这一条读成『通过』。**"))
        out.append(_mk(
            "samples.unresolved_structure", WARN, "分域未解决的结构：无数据",
            observed="结果对象里没有「不分域年龄」表",
            criterion="同上", detail="缺数据不等于没问题。"))
        return

    mswd_s, verdict_s = _col(ov, "MSWD"), _col(ov, "判定")
    if mswd_s is None or verdict_s is None:
        out.append(_mk(
            "samples.whole_spot", WARN, "不分域整段年龄：列缺失",
            observed="「不分域年龄」表里没有 MSWD / 判定 列",
            criterion="该表应按固定列名产出（见 io/report.export_batch）",
            detail="缺列说明产出方与质控层对表的定义已经不一致，"
                   "**同样不要读成通过**。"))
        out.append(_mk(
            "samples.unresolved_structure", WARN, "分域未解决的结构：列缺失",
            observed="同上", criterion="同上", detail="缺列不等于没问题。"))
        return

    mswd = pd.to_numeric(mswd_s, errors="coerce")
    is_const = verdict_s.astype(str) == "整段常数"
    n_total = int(len(ov))
    n_ok = int(is_const.sum())
    frac_ok = (n_ok / n_total) if n_total else float("nan")
    mswd_med = _finite(mswd.median())
    mswd_max = _finite(mswd.max())
    d_bulk = pd.to_numeric(_col(ov, "Δ整段_pct"), errors="coerce")
    age = pd.to_numeric(_col(ov, "年龄_Ma"), errors="coerce")
    out.append(_mk(
        "samples.whole_spot", INFO, "不分域整段年龄：跨测点的统一口径基准",
        observed=f"{n_ok} / {n_total} 个测点与『整段只有一个年龄』相容"
                 f"（{_pct(frac_ok * 100)}），中位 MSWD "
                 + (f"{mswd_med:.2f}" if mswd_med is not None else "—"),
        criterion="信息项：这条不判好坏 —— 它给出的是横向对比的**统一基准**，"
                  "以及有多少点的整段平均值本身有意义",
        detail="相容的测点，整段年龄可以直接用；不相容的是几个年龄的加权"
               "混合值，只能用于跨测点对比，**不能当定年结果**。"
               "逐点数值见「不分域年龄」表。",
        n_spots=n_total, n_compatible=n_ok, compatible_frac=frac_ok,
        mswd_median=mswd_med, mswd_max=mswd_max,
        age_median=_finite(age.median()),
        delta_vs_bulk_pct_median=_finite(d_bulk.median())))

    # ② 分域之后仍报均一、但整段已经非常数 —— 只看域表发现不了的那一类
    n_dom = pd.to_numeric(_col(ov, "域数"), errors="coerce")
    uni = (n_dom == 1) if n_dom is not None else pd.Series(False, index=ov.index)
    n_uni = int(uni.sum())
    bad = ov[uni & ~is_const]
    n_bad = int(len(bad))
    frac_bad = (n_bad / n_uni) if n_uni else float("nan")
    lv = WARN if (n_bad and frac_bad > th.unresolved_structure_frac_warn) else INFO
    labels = [f"{int(r['序号'])} {r['样品']}" for _, r in bad.head(8).iterrows()]
    out.append(_mk(
        "samples.unresolved_structure", lv,
        "分域后仍报『均一』、但整段已不是常数年龄的测点",
        observed=f"{n_bad} / {n_uni} 个均一测点的整段 MSWD 拒绝了常数模型"
                 f"（{_pct(frac_bad * 100)}）",
        criterion=f"占比 ≤ {_pct(th.unresolved_structure_frac_warn * 100)}",
        detail="这些点在「深度剖面域」表里**没有行**，只看域表是发现不了的。"
               "多为未完全校正的残余倾斜，也可能是 σext 偏小 —— 用"
               "「不分域年龄」表的 Δ年龄_pct 与逐点剖面图确认，"
               "别直接把它们的整段年龄当定年结果。",
        n_uniform=n_uni, n_unresolved=n_bad, frac=frac_bad, spots=labels))


def _check_handoff(result, cfg, th, out: List[Check]) -> None:
    """
    ⑤ 交接给 ADEPT 会发生什么 —— 这一组是本模块存在的主要理由。

    现实中"喂 ADEPT 判坪"最反直觉的一点：**输入表里有的测点不等于 ADEPT 会算的
    测点**。ADEPT 对每个测点先判"可用年龄窗口 ≥10 个"，不满足就返回 ok=FALSE，
    **既不报错也不出图**，只在段表里留一行年龄全 NA 的空行。而"可用"的定义是
    `Age68 < 1000 Ma`（否则它要退回 207Pb/206Pb），偏偏 DRUID 的 `剖面窗口`
    表不带 207/206 列 —— 于是老旧锆石（继承核）测点会被整点丢掉。

    这个信息以前只能靠人工对比"输入几个点 / 出来几张图"才发现。放在这里之后，
    下游在**喂之前**就能知道会掉多少、是哪些、为什么。
    """
    stats = _window_stats(result, th)
    if stats is None:
        # 契约上的决定：**三个 key 都要在**，哪怕判不了。
        # 下游按 key 取值时不应该还要处理"这个 key 不存在"这一种情况 ——
        # 那会把"缺表"和"写错 key"混成同一类错误。判不了就报 warn，不报 pass。
        why = ("批次可能用了 --no-depth，或所有测点都没能建立深度剖面。"
               "缺这张表时 ADEPT 也无法判坪 —— 它读的就是这张表。")
        out.append(_mk(
            "handoff.windows", WARN, "没有生成 `剖面窗口` 表，无法交接给 ADEPT",
            observed="—", criterion="存在逐窗口剖面表", detail=why))
        out.append(_mk(
            "handoff.adept_dropout", WARN, "无法判断 ADEPT 会静默丢掉哪些测点",
            observed="—", criterion="需要 `剖面窗口` 表", detail=why,
            n_spots=None, n_dropped=None, frac=None,
            n_old_core=None, n_short_ablation=None, dropped=[]))
        out.append(_mk(
            "handoff.old_core_windows", WARN, "无法统计老核窗口占比",
            observed="—", criterion="需要 `剖面窗口` 表", detail=why,
            n_old_windows=None, n_windows=None, frac=None))
        struct0 = _col(result.results, "深度结构")
        multi0 = int(struct0.astype(str).str.startswith("多域").sum()) if struct0 is not None else 0
        out.append(_mk(
            "samples.multi_domain", INFO,
            f"{multi0} 个样品测点检出多个年龄域（无窗口表，按测点数计）",
            observed=str(multi0), criterion="—（是事实，不是缺陷）",
            detail="一个测点内有多个年龄域时，"
                   "取哪一个作为该测点的年龄是地质判断，算法替不了。",
            n_multi_domain=multi0, n_spots=None))
        return

    n_spots = len(stats)
    n_win = sum(s["n_windows"] for s in stats)
    n_old_win = sum(s["n_windows"] - s["n_usable"] for s in stats)

    out.append(_mk(
        "handoff.windows", PASS,
        f"逐窗口剖面表已生成：{n_spots} 个测点 / {n_win} 个窗口",
        observed=f"{n_spots} 测点 · {n_win} 窗口",
        criterion="存在即可交接",
        n_spots=n_spots, n_windows=n_win))

    # 会被 ADEPT 整点丢掉的
    dropped = []
    for s in stats:
        if s["n_usable"] >= th.min_usable_windows:
            continue
        if s["n_windows"] < th.min_usable_windows:
            reason = "short_ablation"
            why = f"采集窗口只有 {s['n_windows']} 个（<{th.min_usable_windows}），本身数据不足"
        else:
            reason = "old_core"
            why = (f"{s['n_windows']} 个窗口里只有 {s['n_usable']} 个 <"
                   f"{th.old_age_cutoff_ma:g} Ma，其余是老核；"
                   f"输入表没有 207Pb/206Pb 列可回退")
        dropped.append(dict(analysis=s["analysis"], sample=s["sample"],
                            n_windows=s["n_windows"], n_usable=s["n_usable"],
                            reason=reason, why=why))

    n_drop = len(dropped)
    frac = n_drop / n_spots if n_spots else 0.0
    if n_drop == 0:
        lv = PASS
    elif frac > 0.10:
        lv = FAIL
    else:
        lv = WARN
    n_short = sum(1 for d in dropped if d["reason"] == "short_ablation")
    n_core = n_drop - n_short
    out.append(_mk(
        "handoff.adept_dropout", lv,
        ("ADEPT 不会报错，但会静默丢掉这些测点" if n_drop else
         "所有测点都能进入 ADEPT 的判坪流程"),
        observed=f"{n_drop}/{n_spots} = {frac:.1%}"
                 + (f"（老核 {n_core}、短采 {n_short}）" if n_drop else ""),
        criterion=f"占比 ≤ 10%（可用窗口 ≥ {th.min_usable_windows} 个）",
        detail="" if n_drop == 0 else
               "这些测点在 ADEPT 里会被判为 ok=FALSE：输出一行年龄全 NA 的空行、"
               "不画图、不报错。**统计成功率时若只数坪表，会以为它们是\"没有坪\"，"
               "其实它们连管线都没进。** 老核那一类要靠给 `剖面窗口` 表"
               "补 Age76/Age76_1s 两列才能救回来（需升版本、下游重跑）；"
               "短采那一类只能重采。",
        n_spots=n_spots, n_dropped=n_drop, frac=float(frac),
        n_old_core=n_core, n_short_ablation=n_short, dropped=dropped))

    # 老核窗口占全批的比例 —— 即使没到"丢掉整点"的程度，也说明样品里有继承核，
    # 这直接影响"一个测点是不是只有一个年龄域"的解释。
    frac_old = n_old_win / n_win if n_win else 0.0
    lv = WARN if frac_old > 0.05 else PASS
    out.append(_mk(
        "handoff.old_core_windows", lv,
        f"206Pb/238U 年龄 ≥{th.old_age_cutoff_ma:g} Ma 的窗口",
        observed=f"{n_old_win}/{n_win} = {frac_old:.1%}",
        criterion="占比 ≤ 5%",
        detail="" if lv == PASS else
               "这些窗口对 ADEPT 不可用（它要退回 207Pb/206Pb，而输入没这一列）。"
               "存在即说明样品含继承核/老核，解释时要按"
               "「核—边年龄谱系」读，不能直接取加权平均。",
        n_old_windows=n_old_win, n_windows=n_win, frac=float(frac_old)))

    # 多域占比 —— 不是问题，是解释前提。仍然显式报出来，免得下游以为"一个测点一个年龄"。
    struct = _col(result.results, "深度结构")
    multi = 0
    if struct is not None:
        multi = int(struct.astype(str).str.startswith("多域").sum())
    out.append(_mk(
        "samples.multi_domain", INFO,
        f"{multi}/{n_spots} 个样品测点检出多个年龄域",
        observed=f"{multi}/{n_spots}",
        criterion="—（是事实，不是缺陷）",
        detail="一个测点内有多个年龄域时，"
               "取哪一个作为该测点的年龄是地质判断，算法替不了。"
               "下游做加权平均前必须先定这个口径。",
        n_multi_domain=multi, n_spots=n_spots))


def _check_process(result, cfg, th, out: List[Check]) -> None:
    """⑥ 流程与对照 —— 跳过的文件、以及两套方法之间差多少。"""
    skipped = list(result.info.get("skipped") or [])
    out.append(_mk(
        "data.skipped_files", INFO if skipped else PASS,
        f"有 {len(skipped)} 个文件被跳过" if skipped else "没有文件被跳过",
        observed=str(len(skipped)), criterion="—",
        detail="跳过玻璃标样（SRM 612 等）是设计如此，它不参与 U-Pb。"
               "若跳过的不是玻璃标样，说明该测点装载阶段就失败了。",
        skipped=[list(map(str, s)) if isinstance(s, (list, tuple)) else str(s)
                 for s in skipped]))

    alt = next((c for c in ("年龄206_238_ftau", "年龄206_238_simple")
                if c in getattr(result.results, "columns", [])), None)
    roles = _col(result.results, "类型")
    if alt is None or roles is None:
        return
    unk = _sub(result.results, roles == ROLE_LABEL_CN[ROLE_UNKNOWN])
    a = _col(unk, "年龄206_238")
    b = _col(unk, alt)
    if a is None or b is None or getattr(unk, "empty", True):
        return
    m = (a.notna() & b.notna() & (a > 0))
    if not m.any():
        return
    rel = ((b[m] - a[m]) / a[m] * 100)
    med = _finite(rel.median())
    if med is None:
        return
    lv = WARN if abs(med) > th.method_diff_warn_pct else PASS
    out.append(_mk(
        "samples.method_difference", lv,
        f"两套整段比值方法（{cfg.bulk} 与另一套）的年龄差中位",
        observed=f"{med:+.2f}%  (n={int(m.sum())})",
        criterion=f"|中位差| ≤ {th.method_diff_warn_pct:g}%",
        detail="" if lv == PASS else
               "两套方法在同一批数据上的差别本身就是诊断信息：差得远说明"
               "down-hole 分馏严重，或存在别的结构性问题（例如计数率失配）。"
               "同一批数据换一套方法就得到不同年龄，说明结论对方法敏感，"
               "报告里要写清用的是哪一套。",
        method_alt=alt.split("_")[-1], median_pct=med))


def _check_claim_scope(result, cfg, th, out: List[Check]) -> None:
    """
    ⑦ **本批年龄「能用来做什么、不能用来做什么」** —— 把使用限制写成机读的。

    为什么需要这一条
    ----------------
    前面每条检查各管一个事实：校准了没有、σext 是实测还是假设。但下游真正要问的
    是另一句话 —— **这批年龄我能拿它干什么** —— 而这句话现在要靠人把几条 WARN
    在脑子里拼起来。`verdict.level` 也帮不上忙：`warn` 既可能是「标样偏差略大、
    不影响结论」，也可能是「误差棒里有一块根本是假设、结论要改口径」，
    两者同一个等级，下游按等级分支就会把后者当噪音放过。

    最典型的是**没有监控标样的批次**（桂北：样品仓限制放不下 Ple）：
    绝对年龄的**数值**照样可用（它是主标归一化出来的，值本身不含 σext），
    但**误差棒**里有一块是写死的假设值。说成一句「这批判 warn」就把区别抹平了，
    实际后果是有人拿误差棒去论证「两个年龄在误差内一致」。

    判据（两根轴，各自单一门槛，不做加权）
    --------------------------------------
        相对口径（排序 / 批内比较）   只要本批有样品测点就成立，与 σext 无关
        绝对年龄的**值**              要求 校准状态 = 已校准
        绝对年龄的**误差棒**          要求 σext 是**本批实测**、且未触保护上下限

    为什么这里**不判 fail**
    -----------------------
    本条是**派生**的：它一个数都不算，只把前面几条已经报出来的事实翻译成一张
    「用途清单」，并在 `reasons` 里指回**责任检查项的 key**（理由不在本条重述）。
    最严重的等级由原始条目承担 —— 未校准时 `calibration.mode` 本身就是 fail，
    这里再判一次只会让 `counts[fail]` 虚高、`headline` 把同一件事说两遍。
    """
    res = result.results
    roles = _col(res, "类型")
    n_unknown = int((roles == ROLE_LABEL_CN[ROLE_UNKNOWN]).sum()) if roles is not None else 0

    mode_col = _col(res, "校准状态")
    calibrated = bool(mode_col is not None and len(mode_col)
                      and str(mode_col.iloc[0]) == "已校准")

    # σext 的来源与是否触边界：与 `_check_uncertainty` 同一套判法。这里刻意各写
    # 一遍而不去引用那一条的结论 —— 两条检查服务的对象不同：那条说「这个数从哪
    # 来」，这条说「据此能主张什么」；耦合起来会让改动一处牵动另一处的措辞。
    n_sec = int((roles == ROLE_LABEL_CN[ROLE_SECONDARY]).sum()) if roles is not None else 0
    if cfg.sigma_ext68 is not None or cfg.sigma_ext76 is not None:
        src = "forced"
    elif n_sec:
        src = "measured"
    else:
        src = "assumed"

    sd68, sd76 = _finite(result.info.get("sd68")), _finite(result.info.get("sd76"))
    bound = []
    for label, v in (("206/238", sd68), ("207/206", sd76)):
        if v is not None and (abs(v - EXTERNAL_SCATTER_LO) < 1e-12
                              or abs(v - EXTERNAL_SCATTER_HI) < 1e-12):
            bound.append(label)

    # {claim id: None ＝ 可用；否则是**责任检查项的 key**}
    scope: Dict[str, Optional[str]] = {cid: None for cid in CLAIM_IDS}
    if not n_unknown:
        # 没有样品测点就不存在"可解释的年龄"，四项全部不成立。
        for cid in scope:
            scope[cid] = "samples.present"
    elif not calibrated:
        scope["absolute_age_value"] = "calibration.mode"
        scope["absolute_age_uncertainty"] = "calibration.mode"
    elif src != "measured":
        scope["absolute_age_uncertainty"] = "uncertainty.sigma_ext_source"
    elif bound:
        scope["absolute_age_uncertainty"] = "uncertainty.sigma_ext_bound"

    usable = [c for c in CLAIM_IDS if scope[c] is None]
    blocked = {c: r for c, r in scope.items() if r}
    listed = "、".join(f"「{CLAIM_IDS[c]}」" for c in blocked)
    # ⚠ `detail` 必须按**实际挡下它的原因**分派，不能写成「没有监控标样」那一个
    # 故事的通用说明 —— 未校准也走这条分支，而那句话在那时是假的：它跟 σext 无关，
    # 而且那种批次往往**有**监控标样。（本函数第一版就是这么写的，实跑出来才发现。）
    # `tests/test_qc.py` 的结构性测试会检查每个 reason 都在 `detail` 里露过面，
    # 所以这里漏一支会当场红，不会静默留一句错话。
    why = []
    if "samples.present" in blocked.values():
        why.append("本批没有样品测点，不存在可解释的年龄（见 `samples.present`）。")
    if "calibration.mode" in blocked.values():
        why.append("未校准时 F=1，年龄只作相对参考 —— 连绝对年龄的**数值**都不成立"
                   "（见 `calibration.mode`），误差棒更不必说。")
    if "uncertainty.sigma_ext_source" in blocked.values():
        why.append("没有监控标样时，绝对年龄的数值照样是主标归一化出来的（值本身不含"
                   " σext），但误差棒里有一块是写死的假设值 —— 于是用它论证『两个年龄"
                   "在误差内一致』没有依据（见 `uncertainty.sigma_ext_source`）。")
    if "uncertainty.sigma_ext_bound" in blocked.values():
        why.append("σext 触到了保护上下限，它表示的是一个缺口而不是一个测出来的量"
                   "（见 `uncertainty.sigma_ext_bound`）。")
    why.append("`reasons` 给出每一项限制的责任检查项 key，理由不在本条重述；"
               "下游按 key 分支时**只读 `data`**，不要解析 `observed`。")
    out.append(_mk(
        "uncertainty.claim_scope", PASS if not blocked else WARN,
        ("本批年龄的各项用途均成立" if not blocked
         else f"本批年龄有使用限制，不可用于{listed}"),
        observed=(f"可用 {len(usable)}/{len(CLAIM_IDS)} 项"
                  + (f"、不可用 {len(blocked)} 项" if blocked else "")),
        criterion="相对口径要求有样品测点；绝对年龄的『值』要求已校准；"
                  "绝对年龄的『误差棒』要求 σext 为本批实测且未触保护边界",
        detail="" if not blocked else "".join(why),
        usable_for=usable, not_usable_for=list(blocked), reasons=dict(blocked),
        labels=dict(CLAIM_IDS), n_claims=len(CLAIM_IDS),
        calibrated=calibrated, n_unknown=n_unknown,
        sigma_ext_source=src, sigma_ext_bound=bound))


# ═════════════════════════════════════════════════════════════════════════════
# 主入口
# ═════════════════════════════════════════════════════════════════════════════
def assess_batch(result, cfg=None, thresholds: Optional[QCThresholds] = None) -> List[Check]:
    """
    对一个批次的产出做全面质控，返回检查项列表。

    参数
    ----
    result     : `druid.workflow.BatchResult`（鸭子类型：只要有
                 results / qc / info / windows 四个属性即可）
    cfg        : `BatchConfig`，可选。给了才能知道"主标/监控标样叫什么"
                 以及阈值覆盖。缺省时从 `result.info` 里取标样名。
    thresholds : 覆盖阈值；缺省用 `cfg.qc_thresholds`，再缺省用 `QCThresholds()`。

    返回
    ----
    list[Check]，**顺序固定为"越严重越前"**（先按 level，再按 key），
    这样下游和报告都不需要自己排序。

    为什么返回列表而不是一个对象
    ----------------------------
    检查项会随工具演进增减。列表 + 稳定 key 让下游能"只关心它在意的几条"，
    而不是被 Schema 变更牵着走。
    """
    if thresholds is None:
        thresholds = getattr(cfg, "qc_thresholds", None) or QCThresholds()
    th = thresholds

    if cfg is None:
        class _Cfg:                                   # 只读的最小替身
            primary = str(result.info.get("primary", "91500"))
            secondary = str(result.info.get("secondary", "Ple"))
            sigma_ext68 = None
            sigma_ext76 = None
            bulk = str(result.info.get("bulk", "simple"))
            ref_preset = DEFAULT_REF_PRESET   # 没给 cfg 时按产品默认档判自洽性
        cfg = _Cfg()

    out: List[Check] = []
    _check_calibration(result, cfg, th, out)
    _check_standards(result, cfg, th, out)
    _check_uncertainty(result, cfg, th, out)
    _check_samples(result, cfg, th, out)
    _check_whole_spot(result, cfg, th, out)
    _check_handoff(result, cfg, th, out)
    _check_process(result, cfg, th, out)
    # 放在最后：它是对前面全部事实的**派生汇总**（用途清单），自己不判事实。
    _check_claim_scope(result, cfg, th, out)

    # 排序：先按严重程度，同级别按 key 字典序（保证同一份输入永远同一顺序，
    # 否则 diff 出来的报告会到处是无意义的行序变化）。
    out.sort(key=lambda c: (LEVEL_ORDER.get(c.level, 99), c.key))
    return out
