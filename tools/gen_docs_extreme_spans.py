#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools/gen_docs_extreme_spans.py —— 极端跨域测点的综合图与在线展示页
===================================================================

**它做什么**

    python tools/gen_docs_extreme_spans.py scan     # 重算全库域级协和度 → CSV
    python tools/gen_docs_extreme_spans.py build    # 出两份图 + 生成在线页

一次跑出**三张图**（两张是同一批数的本地/线上两版，一张两版共用）：

    · 线上脱敏版 —— 测点一律叫 `G01`…`G12`，不含任何真实批次/样品/测点号，
      就地注入 `docs/extreme-domain-spans.html`（页面里的每个小数都是这里算的）；
    · 本地实名版 —— `汇总库/_分析/extreme_spans_local.svg`（+ .png），
      标签带真实批次与样品号，并另写一份 `extreme_spans_labels.csv` 当对照表；
    · **协和度对照图** —— `concordia_local.svg`（+ .png）。它不含任何测点标签，
      所以线上/本地是同一份，直接注入页面第一节。

都不该手改：改完算法重跑本脚本即可。

**为什么要有这个图**

跨域跨度 > 500 Ma 的测点，是"一个剥蚀坑里同时解出 ~200 Ma 与 ~2000 Ma"
这类结果的集合。它们有两种完全不同的成因，处置相反：

    · 真继承核 / 捕获锆石 —— 值得写进结果；
    · 普通铅污染（抛光粉、环氧、镀碳层、表面吸附铅）—— 必须剔掉。

`206Pb/238U` 分不开这两者（普通铅让两个年龄**一起**变老，看起来仍"协和"），
要的是**域级 207Pb/235U 年龄**。算法见 `tools/domain_span_check.py`（本脚本复用它）。

**图为什么长这样（两个面板，同一批 12 行）**

左：每个域的 `206Pb/238U` 年龄 ±2σ（对数轴），行 = 测点，按跨度降序。
右：同一个域落到 `τ`-协和度平面上，**并垫上全库同 τ 位置的分位包络**。

右边那层灰底是关键，不是装饰：剖面两端的窗口少、`207Pb/206Pb` 计数统计差，
**协和度的自然散布在两端本来就大得多**。不垫这层基准，读者会把 145–170%
当成"明显的污染信号"，而它在最深处只是这个尺度里的普通尾巴。

**口径**

    · 只读 `汇总库/DRUID汇总.db` 与原始 CSV；产出只写到 `汇总库/_分析/`
      （该目录不入版本库）与 `docs/` 一个页面。
    · 跨度只算 `标记 = 'age domain'` 的域；参考值预设取 `批次清单.csv` 的
      `参考值口径` 列（与出结果表时一致）。
    · 协和度 = 域级 207Pb/235U 年龄 ÷ 206Pb/238U 年龄 × 100%。
    · `τ` 是**积分段的相对位置**（0 = 积分段起点、1 = 终点），不是"激光开/关"。

**改这个脚本时别再犯的坑**

1. SVG 根元素必须带 `class="fig"`，且线上版**不能**内嵌 `<style>` ——
   配色走页面 CSS 变量，深色主题才能跟着变；本地版则**必须**内嵌一套浅色
   `<style>`，否则单独打开会是"尺寸正确、内容全黑"。
2. 图例与图注放在 SVG 外面（HTML）。画进坐标系一定会压住某个角上的数据。
3. 线上版走 `test_no_identifying_strings_in_tracked_files`：标签、表头、
   图注一个真名都不能有，且**未跟踪的新文件也在扫描范围内**。
4. 分位数一律**先滤掉非有限值**再排序：本批实测有 NaN 与负协和度，
   混进去会让 `median` 静默给出一个自相矛盾的值（本轮踩过）。
"""
from __future__ import annotations

import argparse
import csv
import html
import math
import pathlib
import subprocess
import sys
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "tools"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from domain_span_check import (group_by_batch, read_spot_list,   # noqa: E402
                               recompute_batch)

DB_DEFAULT = ROOT / "汇总库" / "DRUID汇总.db"
LIB_CSV = ROOT / "汇总库" / "_分析" / "domain_concordance.csv"
OUT_DIR = ROOT / "汇总库" / "_分析"
PAGE = ROOT / "docs" / "extreme-domain-spans.html"
EDGE = pathlib.Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")

BEGIN = "<!-- ==== 极端跨域测点综合图：由 tools/gen_docs_extreme_spans.py 生成，请勿手改 ==== -->"
END = "<!-- ==== 生成区结束 ==== -->"
CONC_BEGIN = "<!-- ==== 协和度对照图：由 tools/gen_docs_extreme_spans.py 生成，请勿手改 ==== -->"
CONC_END = "<!-- ==== 生成区结束 ==== -->"

SPAN_DEFAULT = 500.0        # Ma，跨度门槛
NBIN = 10                   # 全库背景按 τ 分 10 箱
CONC_MAX = 200.0            # 右面板 x 轴上限（%），超出者夹到轴端
AGE_LO, AGE_HI = 150.0, 3000.0
AGE_TICKS = (200, 300, 500, 700, 1000, 1500, 2000, 2500)

# 协和度分箱：<105 / 105–125 / 125–150 / ≥150。四档，两套主题各给一组颜色。
CONC_EDGES = (105.0, 125.0, 150.0)

# 版面。面板宽度是照着"线上版整幅 ≈ 站点正文宽度（876 px）"定的 ——
# 否则 SVG 被等比缩小，11 px 的标注会掉到 9 px，等于白画。
# 行标签是**多行块**（G 码 / 真实标识 / 跨度），全在绘图区左边 ——
# 早先把"跨度"画在行内右上角，正好压住最深端那个点。
ROW_H = 40
LINE_H = 12
HEAD = 32
TAIL = 44
TOP = 8
PA_W = 460
GAP = 52
PB_W = 250
MR = 24


# ─────────────────────────────────────────────────────────────────────────────
# 一、数据
# ─────────────────────────────────────────────────────────────────────────────
def _pct(vals, p: float) -> float:
    """最近秩分位数。调用方保证 vals 非空且已滤掉非有限值。"""
    v = sorted(vals)
    k = min(len(v) - 1, max(0, int(round(p / 100.0 * (len(v) - 1)))))
    return float(v[k])


def _finite(x) -> bool:
    return isinstance(x, float) and math.isfinite(x)


def library_bands(csv_path: pathlib.Path):
    """
    全库（正式年龄域）按 τ 分箱的分位包络。返回 (bands, 总域数)。

    `τ` 用域中点。非有限值与超范围值都保留参与分位计算（它们是最深端真实
    存在的极端值），只在**画**的时候夹到轴端 —— 反过来做会把右尾削掉，
    恰好把要讨论的现象滤掉。
    """
    if not csv_path.exists():
        return [], 0
    bins: dict[int, list] = defaultdict(list)
    total = 0
    with open(csv_path, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r.get("flag") != "age domain":
                continue
            total += 1
            try:
                tm = 0.5 * (float(r["tau0"]) + float(r["tau1"]))
                c = float(r["concordance"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(c):
                continue
            k = min(NBIN - 1, max(0, int(tm * NBIN)))
            bins[k].append(c)
    out = []
    for k in range(NBIN):
        v = bins.get(k)
        if not v:
            continue
        out.append(dict(tau=(k + 0.5) / NBIN, n=len(v),
                        p05=_pct(v, 5), p25=_pct(v, 25), p50=_pct(v, 50),
                        p75=_pct(v, 75), p95=_pct(v, 95),
                        frac_hi=100.0 * sum(1 for x in v if x > 130.0) / len(v)))
    return out, total


def collect_groups(db: pathlib.Path, span: float):
    """挑出候选测点并逐批重算。返回 (spots 顺序表, {批次: {序号: 审计结果}})。"""
    spots = read_spot_list(db, span, include_transitional=False)
    audited: dict[str, dict] = {}
    for batch, group in group_by_batch(spots).items():
        folder = group[0]["文件夹"]
        data_dir = ROOT / "汇总库" / folder / "原始数据"
        preset = group[0].get("参考值口径") or "horstwood2016"
        if not data_dir.exists():
            print("[跳过] %s：找不到 %s" % (batch, data_dir), file=sys.stderr)
            continue
        print("[重算] %s（%d 点，预设 %s）" % (batch, len(group), preset), file=sys.stderr)
        try:
            audited[batch] = recompute_batch(data_dir, preset, {s["序号"] for s in group})
        except Exception as exc:                                       # noqa: BLE001
            print("[失败] %s：%s" % (batch, exc), file=sys.stderr)
    return spots, audited


def domain_rows(spots, audited):
    """
    摊平成"一行一个域"。只收正式年龄域，按 τ 升序。

    没解出域的测点照样保留（`domains` 空），否则图上会少一行、与页面上
    声称的个数对不上。
    """
    rows = []
    for i, s in enumerate(spots, 1):
        aud = audited.get(s["批次"], {}).get(s["序号"])
        ds = [d for d in (aud or {}).get("domains", []) if d["flag"] == "age domain"]
        ds.sort(key=lambda d: d["tau0"])
        rows.append(dict(code="G%02d" % i, 序号=s["序号"], 样品=s["样品"],
                         批次=s["批次"], 区组=s["区组"], span=float(s["span"]),
                         domains=ds,
                         note=(aud or {}).get("note", "")))
    return rows


def _clip(c: float) -> float:
    if not _finite(c):
        return 0.0
    return min(CONC_MAX, max(0.0, c))


# ─────────────────────────────────────────────────────────────────────────────
# 二、画图
# ─────────────────────────────────────────────────────────────────────────────
INLINE_STYLE = """<style>
text{font-family:"Segoe UI","Microsoft YaHei",sans-serif}
.mono{font-family:"Cascadia Mono",Consolas,monospace}
.fig-bg{fill:#FFFFFF;stroke:#D3D1C7}
.fig-frame{fill:none;stroke:#C4C0B3}
.fig-grid{stroke:#D3D1C7;stroke-width:.7;opacity:.75}
.fig-rowsep{stroke:#E8E6DF;stroke-width:.9}
.fig-tick{fill:#888780;font-size:11px}
.fig-axis{fill:#5F5E5A;font-size:11.5px}
.fig-title{fill:#2C2C2A;font-size:12.5px;font-weight:600}
.fig-lab{fill:#2C2C2A;font-size:12px;font-weight:600}
.fig-lab2{fill:#888780;font-size:10.5px}
.fig-sub{fill:#888780;font-size:11px}
.c1{stroke:#0F6E56;fill:#0F6E56}
.c2{stroke:#854F0B;fill:#854F0B}
.c3{stroke:#C2410C;fill:#C2410C}
.c4{stroke:#A32D2D;fill:#A32D2D}
.fig-ebar{stroke-width:2.6;stroke-linecap:round}
.fig-dot{stroke:none}
.fig-halo{stroke:#FFFFFF;stroke-width:1}
.bg-band{fill:#8B8983;opacity:.22;stroke:none}
.bg-band2{fill:#5F5E5A;opacity:.30;stroke:none}
.bg-med{stroke:#2C2C2A;stroke-width:1.6;fill:none}
.bg-ref{stroke:#A32D2D;stroke-width:1;stroke-dasharray:5 4;opacity:.75}
.cc-band{fill:#8B8983;opacity:.16;stroke:none}
.cc-ray{stroke:#888780;stroke-width:1;stroke-dasharray:4 3}
.cc-ray2{fill:none;stroke:#888780;stroke-width:1;stroke-dasharray:4 3;opacity:.75}
.cc-curve{fill:none;stroke:#A32D2D;stroke-width:2}
.cc-line{stroke:#A32D2D;stroke-width:2}
.cc-bg{fill:#8B8983;opacity:.30;stroke:none}
.cc-lab{fill:#5F5E5A;font-size:10.5px;font-family:Consolas,monospace;
        stroke:#FFFFFF;stroke-width:2.6;paint-order:stroke}
</style>"""


def _conc_class(c: float) -> str:
    if not _finite(c):
        return "c4"
    if c < CONC_EDGES[0]:
        return "c1"
    if c < CONC_EDGES[1]:
        return "c2"
    if c < CONC_EDGES[2]:
        return "c3"
    return "c4"


def _f(v, spec="%.1f", dash="--"):
    return dash if not _finite(v) else spec % v


def svg_figure(rows, bands, labels, inline_style: bool) -> tuple[str, int, int]:
    """
    拼出整张 SVG。`labels[i]` 是第 i 行的**多行**标签（第 1 行是 G 码加粗，
    其余为小字）。返回 (svg, 宽, 高)。
    """
    ml = max((max((len(s) for s in block), default=1) for block in labels), default=8)
    ML = int(9.0 * ml) + 18                      # 主标签 mono 12px ≈ 7.2px/字
    W = ML + PA_W + GAP + PB_W + MR
    n = len(rows)
    PT = TOP + HEAD
    PH = n * ROW_H
    PB_ = PT + PH
    H = PB_ + TAIL + 10           # 末行与轴标题之下再留一点白

    def xa(a: float) -> float:
        lo, hi = math.log10(AGE_LO), math.log10(AGE_HI)
        return ML + (math.log10(max(a, 1e-9)) - lo) / (hi - lo) * PA_W

    def xb(c: float) -> float:
        return ML + PA_W + GAP + _clip(c) / CONC_MAX * PB_W

    def yb(t: float) -> float:
        return PT + max(0.0, min(1.0, t)) * PH

    P = ['<svg class="fig" %sviewBox="0 0 %d %d" role="img" '
         'aria-labelledby="figSpanT figSpanD" aria-describedby="figSpanD">'
         % ('width="%d" height="%d" ' % (W, H) if inline_style else '', W, H)]
    P.append('<title id="figSpanT">跨域跨度异常测点的域级年龄与域级协和度</title>')
    P.append('<desc id="figSpanD">左：每个域的 206Pb/238U 年龄与 2σ 误差，'
             '按测点分行、年龄取对数轴。右：同一个域落在 τ-域级协和度平面上；'
             '灰色包络是全库全部正式年龄域在同一 τ 位置的分位区间（中位、'
             '25–75% 与 95%），红线是 100% 协和线。</desc>')
    if inline_style:
        P.append(INLINE_STYLE)
    P.append('<rect class="fig-bg" x="0" y="0" width="%d" height="%d" rx="12"/>' % (W, H))

    # 列标题
    P.append('<text class="fig-title" x="%d" y="%d">年龄 206Pb/238U（Ma，对数轴，横杠 = 2σ）</text>'
             % (ML, TOP + 14))
    P.append('<text class="fig-title" x="%d" y="%d">域级协和度 207Pb/235U ÷ 206Pb/238U（%%）</text>'
             % (ML + PA_W + GAP, TOP + 14))
    P.append('<text class="fig-sub" x="%d" y="%d" text-anchor="end">'
             '灰底 = 全库同 τ 位置的分位包络</text>' % (ML + PA_W + GAP + PB_W, TOP + 27))
    # 右面板的纵轴名：它是 τ（0 = 积分段起点 → 1 = 终点），光有 0.0~1.0 的刻度
    # 读者不知道那是什么。放不下旋转的轴标题（两面板之间只有 52 px），
    # 就把 "τ" 摆在刻度列的正上方。
    P.append('<text class="fig-sub" x="%d" y="%d" text-anchor="end">τ</text>'
             % (ML + PA_W + GAP - 8, TOP + 27))

    # ── 左面板 ──
    for t in AGE_TICKS:
        x = xa(t)
        P.append('<line class="fig-grid" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>' % (x, PT, x, PB_))
        P.append('<text class="fig-tick" x="%.1f" y="%d" text-anchor="middle">%d</text>'
                 % (x, PB_ + 16, t))
    for i in range(1, n):
        y = PT + i * ROW_H
        P.append('<line class="fig-rowsep" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (ML, y, ML + PA_W, y))
    P.append('<text class="fig-axis" x="%.1f" y="%d" text-anchor="middle">'
             '206Pb/238U 年龄（Ma）</text>' % (ML + PA_W / 2, PB_ + 36))

    for i, r in enumerate(rows):
        yc = PT + i * ROW_H + ROW_H / 2
        block = labels[i]
        y0 = yc - (len(block) - 1) * LINE_H / 2.0 + 4
        for j, s in enumerate(block):
            if not s:
                continue
            cls = "fig-lab mono" if j == 0 else "fig-lab2 mono"
            P.append('<text class="%s" x="%d" y="%.1f" text-anchor="end">%s</text>'
                     % (cls, ML - 14, y0 + j * LINE_H, html.escape(s)))
        for d in r["domains"]:
            cls = _conc_class(d["concordance"])
            lo, hi = d["age68"] - d["s2_68"], d["age68"] + d["s2_68"]
            P.append('<line class="fig-ebar %s" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                     % (cls, xa(lo), yc, xa(hi), yc))
            P.append('<circle class="fig-dot %s" cx="%.1f" cy="%.1f" r="3.8"/>'
                     % (cls, xa(d["age68"]), yc))
            P.append('<text class="fig-sub mono" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                     % (xa(d["age68"]), yc + 15, d["domain"]))

    # ── 右面板：全库背景（判读基准，必须看得清） ──
    # 两层：P5–P95 浅、P25–P75 深。只画 25–75 会细成一条线 ——
    # 中段那一箱 IQR 只有 98–102，而它要表达的是"这里几乎不散"，对比才能读出来。
    bh = ROW_H * 0.46
    for bd in bands:
        y = yb(bd["tau"])
        P.append('<rect class="bg-band" x="%.1f" y="%.1f" width="%.1f" height="%.1f"/>'
                 % (xb(bd["p05"]), y - bh / 2, max(0.6, xb(bd["p95"]) - xb(bd["p05"])), bh))
        P.append('<rect class="bg-band2" x="%.1f" y="%.1f" width="%.1f" height="%.1f"/>'
                 % (xb(bd["p25"]), y - bh / 2, max(0.6, xb(bd["p75"]) - xb(bd["p25"])), bh))
        P.append('<line class="bg-med" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (xb(bd["p50"]), y - bh / 2, xb(bd["p50"]), y + bh / 2))
    x100 = xb(100.0)
    P.append('<line class="bg-ref" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>' % (x100, PT, x100, PB_))

    for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        y = yb(t)
        P.append('<line class="fig-grid" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (ML + PA_W + GAP, y, ML + PA_W + GAP + PB_W, y))
        P.append('<text class="fig-tick" x="%d" y="%.1f" text-anchor="end">%.1f</text>'
                 % (ML + PA_W + GAP - 8, y + 4, t))
    for c in (0, 50, 100, 150, 200):
        P.append('<text class="fig-tick" x="%.1f" y="%d" text-anchor="middle">%d</text>'
                 % (xb(c), PB_ + 16, c))
    P.append('<text class="fig-axis" x="%.1f" y="%d" text-anchor="middle">协和度（%%）</text>'
             % (ML + PA_W + GAP + PB_W / 2, PB_ + 36))

    # ── 右面板：本图这 12 行（同一个域，同一个颜色） ──
    for r in rows:
        for d in r["domains"]:
            tm = 0.5 * (d["tau0"] + d["tau1"])
            P.append('<circle class="fig-dot fig-halo %s" cx="%.1f" cy="%.1f" r="4.4"/>'
                     % (_conc_class(d["concordance"]), xb(d["concordance"]), yb(tm)))

    P.append('<rect class="fig-frame" x="%d" y="%d" width="%d" height="%d"/>' % (ML, PT, PA_W, PH))
    P.append('<rect class="fig-frame" x="%d" y="%d" width="%d" height="%d"/>'
             % (ML + PA_W + GAP, PT, PB_W, PH))
    P.append('</svg>')
    return "\n".join(P), W, H


# ─────────────────────────────────────────────────────────────────────────────
# 三、HTML 片段
# ─────────────────────────────────────────────────────────────────────────────
def legend_html() -> str:
    return """<ul class="figlegend">
<li class="key"><b>怎么读</b>
<span class="k"><i class="sw c1"></i>&lt; 105%　与协和线相容</span>
<span class="k"><i class="sw c2"></i>105–125%</span>
<span class="k"><i class="sw c3"></i>125–150%</span>
<span class="k"><i class="sw c4"></i>≥ 150%　207Pb 明显过量</span>
<span class="k">右图灰底 = 全库同 τ 位置的分位包络（深框 25–75%、浅框 5–95%、竖线 = 中位）</span>
<span class="k"><i class="sw ref"></i>红色虚线 = 100% 协和线</span></li>
</ul>"""


def bench_html(rows, bands, total):
    """读图基准：把"两端本来就更散"这件事用全库数字说清楚。"""
    deep = [b for b in bands if b["tau"] >= 0.85]
    mid = [b for b in bands if 0.35 <= b["tau"] < 0.65]
    lo = min(bands, key=lambda b: b["tau"]) if bands else None
    p95_deep = max((b["p95"] for b in deep), default=float("nan"))
    frac_deep = max((b["frac_hi"] for b in deep), default=float("nan"))
    med_mid = (sum(b["p50"] for b in mid) / len(mid)) if mid else float("nan")
    p95_mid = max((b["p95"] for b in mid), default=float("nan"))
    n_dom = sum(len(r["domains"]) for r in rows)
    n_multi = sum(1 for r in rows if len(r["domains"]) >= 2)
    over = sum(1 for r in rows for d in r["domains"] if d["concordance"] > 130.0)
    C = []
    C.append('<div class="kpi">')
    C.append('<div><b>%d</b><span>跨度 &gt; %g Ma 的测点（只算正式年龄域）</span></div>'
             % (len(rows), SPAN_DEFAULT))
    C.append('<div><b>%.0f</b><span>其中解出 ≥2 个域的测点</span></div>' % n_multi)
    C.append('<div><b>%d / %d</b><span>域级协和度 &gt; 130%% 的域 / 全部域</span></div>' % (over, n_dom))
    C.append('<div><b>%.0f%%</b><span>全库最深端（τ≥0.85）域级协和度的 95 分位</span></div>' % p95_deep)
    C.append('</div>')
    C.append('<div class="note">')
    C.append('<p><b>右图的灰底不是装饰，是判读基准。</b>全库 %d 个正式年龄域，'
             '按域中点 τ 分 %d 箱算分位：中段（τ 0.35–0.65）中位协和度 %.1f%%、'
             '95 分位只 %.0f%%；而最深端（τ≥0.85）的 95 分位是 <b>%.0f%%</b>'
             '、超过 130%% 的占到 <b>%.1f%%</b>。'
             '剖面两端窗口少、207Pb/206Pb 的计数统计差，协和度的自然散布本来就在那里最大 ——'
             '所以"某个域报出 145%%"这件事，<b>要跟同 τ 位置比，不能跟 100%% 比</b>。</p>'
             % (total, NBIN, med_mid, p95_mid, p95_deep, frac_deep))
    if lo is not None:
        C.append('<p>最浅端（τ≈%.2f，n=%d）的中位只有 %.1f%%，看着很乖，'
                 '但它的 95 分位到了 <b>%.0f%%</b>、超过 130%% 的占 %.1f%% —— '
                 '与最深端同一个量级。<b>两端都不稳，只是成因不同</b>：'
                 '最浅端紧挨积分段起点，那里 207Pb 的净信号可能低到 0 附近，'
                 '比值在分母很小时被放大，所以既出很大的数、也出负数；'
                 '最深端则是窗口少、207Pb/206Pb 的计数统计差。'
                 '中段（τ 0.35–0.65）才谈得上"接近 100%% = 真正协和"。</p>'
                 % (lo["tau"], lo["n"], lo["p50"], lo["p95"], lo["frac_hi"]))
    C.append('</div>')
    return "\n".join(C)


def table_html(rows, show_real: bool) -> str:
    """逐点明细表。`show_real=False` 时只出现 G 码，不出现任何真名。"""
    T = ['<table>']
    T.append('<thead><tr>'
             '<th>测点</th><th>跨度 (Ma)</th><th>域</th><th>τ 区间</th><th>窗数</th>'
             '<th>206Pb/238U (Ma)</th><th>2σ</th><th>207Pb/235U (Ma)</th><th>2σ</th>'
             '<th>域级协和度</th><th>域内 MSWD</th><th>Th/U</th></tr></thead><tbody>')
    for r in rows:
        ds = r["domains"]
        n = max(1, len(ds))
        if show_real:
            name = '%s<br><span class="dim">%s · %s#%d</span>' % (
                r["code"], html.escape(r["批次"]), html.escape(r["样品"]), r["序号"])
        else:
            name = r["code"]
        first = ('<td class="n" rowspan="%d">%s</td>'
                 '<td class="n" rowspan="%d">%.0f</td>' % (n, name, n, r["span"]))
        if not ds:
            T.append('<tr>%s<td colspan="10">该点没有解出正式年龄域（%s）</td></tr>'
                     % (first, html.escape(r["note"] or "未说明")))
            continue
        for j, d in enumerate(ds):
            tds = first if j == 0 else ""
            cls = "n hi" if d["concordance"] >= 150.0 else "n"
            T.append('<tr>%s<td class="n">%s</td><td class="n">%.2f – %.2f</td>'
                     '<td class="n">%d</td><td class="n">%.1f</td><td class="n">%.1f</td>'
                     '<td class="n">%.1f</td><td class="n">%.1f</td>'
                     '<td class="%s">%.1f%%</td><td class="n">%.2f</td>'
                     '<td class="n">%.3f</td></tr>'
                     % (tds, d["domain"], d["tau0"], d["tau1"], d["n_win"],
                        d["age68"], d["s2_68"], d["age75"], d["s2_75"],
                        cls, d["concordance"], d["mswd"], d["th_u"]))
    T.append('</tbody></table>')
    return "\n".join(T)


def caption_html(rows, bands) -> str:
    d = [b for b in bands if b["tau"] >= 0.85]
    up = sum(1 for r in rows if len(r["domains"]) >= 2
             and r["domains"][0]["age68"] < r["domains"][-1]["age68"])
    dn = sum(1 for r in rows if len(r["domains"]) >= 2
             and r["domains"][0]["age68"] > r["domains"][-1]["age68"])
    n_dom = sum(len(r["domains"]) for r in rows)
    three = sum(1 for r in rows if len(r["domains"]) >= 3)
    high = [x for r in rows for x in r["domains"] if x["concordance"] >= 150.0]
    hi = len(high)

    def _end(x):
        return 0.5 * (x["tau0"] + x["tau1"]) < 0.15 or 0.5 * (x["tau0"] + x["tau1"]) > 0.8

    n_end = sum(1 for x in high if _end(x))
    if not high:
        where = "这一批里没有域级协和度 ≥ 150% 的域"
    elif n_end == hi:
        where = "而它们<b>全部落在剖面的两端</b>（域中点 τ &lt; 0.15 或 &gt; 0.8）"
    else:
        where = ("其中 <b>%d / %d 个落在剖面的两端</b>（域中点 τ &lt; 0.15 或 &gt; 0.8）"
                 % (n_end, hi))
    p95_deep = max((b["p95"] for b in d), default=float("nan"))
    return ('<figcaption>%d 个测点、共 %d 个正式年龄域。左图每行一个测点（按跨度降序），'
            '圆点是一个域，横杠是它的 2σ。右图把同一个域放回 τ 轴：'
            '<b>%d 个域里有 %d 个域级协和度 ≥ 150%%</b>，%s —— '
            '而全库同位置（τ≥0.85）的 95 分位就有 %.0f%%。也就是说这种值'
            '在本剖面的两端并不罕见，<b>要跟同 τ 位置比，不能跟 100%% 比</b>。'
            '域序<b>沿坑深递减</b>的 %d 个、<b>递增</b>的 %d 个，'
            '另有 %d 个测点解出三个域且不单调 —— 两种方向同时存在，'
            '说明至少一部分是真实的空间结构，不是单一系统偏差。'
            '两端的域窗数普遍只有 3–8 个，MSWD 也大 —— 明细见下表。'
            '图上没有任何批次号、样品号或测点序号；测点一律叫 G01–G%02d。</figcaption>'
            % (len(rows), n_dom, n_dom, hi, where, p95_deep, dn, up, three, len(rows)))


def row_labels(rows, show_real: bool) -> list:
    """行标签块。脱敏版只有 G 码 + 跨度；实名版中间插一行真实批次与样品号。"""
    out = []
    for r in rows:
        blk = [r["code"]]
        if show_real:
            blk.append("%s %s#%d" % (r["批次"], r["样品"], r["序号"]))
        blk.append("跨度 %.0f Ma" % r["span"])
        out.append(blk)
    return out


def build_fragment(rows, bands, total, show_real: bool) -> str:
    labels = row_labels(rows, show_real)
    svg, _w, _h = svg_figure(rows, bands, labels, inline_style=False)
    F = ['<figure>', svg, caption_html(rows, bands), '</figure>']
    F.append(legend_html())
    F.append(bench_html(rows, bands, total))
    F.append('<h3>逐点明细</h3>')
    F.append(table_html(rows, show_real))
    return "\n".join(F)


def concordia_legend(n_lib: int) -> str:
    return """<ul class="figlegend">
<li class="key"><b>怎么读</b>
<span class="k"><i class="sw c1"></i>&lt; 105%%</span>
<span class="k"><i class="sw c2"></i>105–125%%</span>
<span class="k"><i class="sw c3"></i>125–150%%</span>
<span class="k"><i class="sw c4"></i>≥ 150%%</span>
<span class="k">虚线 = 等协和度线 50/70/80/90/110/120/150/200%%（左图因取对数而互相平行）</span>
<span class="k">灰带 = 常用判据 90–110%%</span>
<span class="k">灰点 = 全库 %d 个正式年龄域（两图各按自己的轴范围裁剪）</span>
<span class="k">彩点 = 本文这 12 个测点解出的域</span></li>
</ul>""" % n_lib


def concordia_fragment(rows, lib) -> str:
    svg, _w, _h = svg_concordia(rows, lib, inline_style=False)
    return "\n".join(['<figure>', svg, concordia_caption(rows, lib), '</figure>',
                      concordia_legend(len(lib))])


def inject(fragment: str, page: pathlib.Path,
           begin: str = BEGIN, end: str = END) -> None:
    text = page.read_text(encoding="utf-8")
    i = text.find(begin)
    if i == -1:
        raise SystemExit("%s 里找不到生成区起点 %r" % (page, begin))
    j = text.find(end, i + len(begin))          # 必须从起点之后找：
    if j == -1:                                 # 页面里有两个区，END 是同一个串
        raise SystemExit("%s 里找不到生成区终点" % page)
    page.write_text(text[:i] + begin + "\n" + fragment + "\n" + end
                    + text[j + len(end):], encoding="utf-8")


# ─────────────────────────────────────────────────────────────────────────────
# 三之二、协和度对照图：把"等协和度"画出来
# ─────────────────────────────────────────────────────────────────────────────
# 关键一句话：**协和度就是斜率**。
#   在「年龄–年龄」空间里，横轴 = 206/238 年龄、纵轴 = 207/235 年龄，
#   于是 协和度 = 纵轴 ÷ 横轴 = 从原点出发的斜率 ⇒ 等协和度是一条**射线族**，
#   协和线（100%）就是 45° 线。这比在比值空间里读曲线直观得多。
# 右图给出同一批点在**比值空间**（教科书上的 Wetherill 谐和线）里的样子，
#   那里的等协和度**不是**直线（两条衰变链的 λ 不同），但同样从原点发散 ——
#   两张图放一起，"协和度 = 相对原点的斜率"这件事就没法误解了。
CONC_RAYS = (50.0, 70.0, 80.0, 90.0, 110.0, 120.0, 150.0, 200.0)
CONC_BAND = (90.0, 110.0)          # 常用判据带
AGE_AX_LO, AGE_AX_HI = 100.0, 3500.0
R68_AX_HI, R75_AX_HI = 0.66, 18.5   # 0.66 ≈ 3140 Ma：把 3000 Ma 那根刻度留在轴内
# 两条衰变链的常数（与 druid.core.geochronology 同值）。**只在这里写一次**：
# 四个调用点各写一遍的话，改了一处忘了另一处，图上的协和线与数据点就不再自洽。
LAM68 = 1.55125e-10                 # 238U，yr⁻¹
LAM35 = 9.8485e-10                  # 235U，yr⁻¹


def library_points(csv_path: pathlib.Path):
    """全库正式年龄域的 (206/238 年龄, 207/235 年龄)，用作两图的背景云。"""
    pts = []
    if not csv_path.exists():
        return pts
    with open(csv_path, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r.get("flag") != "age domain":
                continue
            try:
                a68, a75 = float(r["age68"]), float(r["age75"])
            except (TypeError, ValueError):
                continue
            if math.isfinite(a68) and math.isfinite(a75) and a68 > 0 and a75 > 0:
                pts.append((a68, a75))
    return pts


def _to_panel(panel, a68, a75, xa, ya, xb, yb):
    """
    把一个 (206/238 年龄, 207/235 年龄) 投到指定面板；越界返回 None。

    `panel="age"` 直接画年龄；`panel="ratio"` 先按两条衰变链换算成比值 ——
    λ₂₃₈ = 1.55125e-10、λ₂₃₅ = 9.8485e-10（与 druid.core.geochronology 同值）。
    用的是本模块顶部的 `LAM68` / `LAM35` —— **不 import druid**，这个脚本
    只读数据库和一张 csv，不该牵扯包的运行时状态。
    """
    if panel == "age":
        if not (AGE_AX_LO <= a68 <= AGE_AX_HI and AGE_AX_LO <= a75 <= AGE_AX_HI):
            return None
        return xa(a68), ya(a75)
    r68 = math.expm1(LAM68 * a68 * 1e6)
    r75 = math.expm1(LAM35 * a75 * 1e6)
    if not (0 < r68 <= R68_AX_HI and 0 < r75 <= R75_AX_HI):
        return None
    return xb(r68), yb(r75)


def _cascade(items, sep, key, lower_is_first):
    """
    沿同一条边把标签错开：同边标签坐标太近就沿边平移，直到两两相隔 ≥ sep。

    只平移几像素，标签仍然贴着自己那条线；比"把标签挪到画布另一头"诚实得多。
    `items` 是 [(sort_coord, ...), …]（会被原地改 sort_coord）。
    """
    items.sort(key=key, reverse=not lower_is_first)
    for i in range(1, len(items)):
        prev = key(items[i - 1])
        if abs(key(items[i]) - prev) < sep:
            items[i][0] = prev + sep if lower_is_first else prev - sep
    return items


def svg_concordia(rows, lib, inline_style: bool) -> tuple[str, int, int]:
    """
    两面板，回答同一个问题："协和度这个数，落到图上是什么样子。"

    · 左 = **年龄–年龄**（两轴都取对数）：协和度 = 纵值 ÷ 横值。两轴都取对数后，
      等协和度线是**一族斜率同为 1 的平行线** —— 100% 那条就是 45° 协和线，
      离它越远 = 偏离 100% 越多。
    · 右 = **比值空间**（教科书上的 Wetherill 谐和线，两轴都线性）：协和线是
      曲线；等协和度线从原点发散，但**既不是直线、也不平行**。

    ⚠ 左图**不要**写成"从原点发散的射线" —— 那句话是在**线性**坐标下成立的，
    本图是对数轴，画出来是平行线。写的话要跟着画走。
    """
    MLA, PW = 60, 370
    GAPX = 74                  # 两面板之间要放下右图的 y 轴标题 + 刻度数字
    MLB = MLA + PW + GAPX
    MR = 84                  # 右侧留一列「对应年龄」小字
    W = MLB + PW + MR
    TOP, TT = 8, 46          # TT：面板标题行（标题 + 两行小字）
    PT = TOP + TT
    PH = PW                  # 与 PW 相等 ⇒ 对数坐标下 45° 线真的是 45°
    PB_ = PT + PH
    H = PB_ + 48
    n_lib = len(lib)
    n_obs = sum(len(r["domains"]) for r in rows)

    lo, hi = math.log10(AGE_AX_LO), math.log10(AGE_AX_HI)

    def xa(a):
        return MLA + (math.log10(max(a, 1e-9)) - lo) / (hi - lo) * PW

    def ya(a):
        return PB_ - (math.log10(max(a, 1e-9)) - lo) / (hi - lo) * PH

    def xb(r):
        return MLB + min(R68_AX_HI, max(0.0, r)) / R68_AX_HI * PW

    def yb(r):
        return PB_ - min(R75_AX_HI, max(0.0, r)) / R75_AX_HI * PH

    # 两个面板各自真正画得出的灰点数。**不要写死** —— 轴范围一改，
    # 手写的"两图都是 1480 个"立刻变成假话（这一版就踩到：左图只画得出 1401）。
    n_draw = {"age": 0, "ratio": 0}
    for a68, a75 in lib:
        for panel in ("age", "ratio"):
            if _to_panel(panel, a68, a75, xa, ya, xb, yb) is not None:
                n_draw[panel] += 1

    P = ['<svg class="fig" %sviewBox="0 0 %d %d" role="img" '
         'aria-labelledby="figConcT figConcD" aria-describedby="figConcD">'
         % ('width="%d" height="%d" ' % (W, H) if inline_style else '', W, H)]
    P.append('<title id="figConcT">协和度代表什么：等协和度线族与实测点</title>')
    P.append('<desc id="figConcD">左图横轴为 206Pb/238U 年龄、纵轴为 207Pb/235U 年龄，'
             '两者都取对数，于是协和度等于纵值除以横值；100%% 的协和线是 45 度线，'
             '虚线是 50 到 200%% 的等协和度线，因两轴同取对数而互相平行。'
             '右图是同一批点在 207Pb/235U 对 206Pb/238U 的比值空间'
             '（教科书上的 Wetherill 谐和线），协和线是一条曲线，'
             '等协和度线从原点发散。两图的灰点都取自全库 %d 个正式年龄域，'
             '各自按轴范围裁剪（左图 %d 个、右图 %d 个）；'
             '彩点是本文讨论的 %d 个域。</desc>'
             % (n_lib, n_draw["age"], n_draw["ratio"], n_obs))
    if inline_style:
        P.append(INLINE_STYLE)
    P.append('<rect class="fig-bg" x="0" y="0" width="%d" height="%d" rx="12"/>' % (W, H))

    P.append('<text class="fig-title" x="%d" y="%d">'
             '左：年龄–年龄（两轴都取对数）—— 协和度 = 纵轴 ÷ 横轴</text>'
             % (MLA, TOP + 14))
    P.append('<text class="fig-sub" x="%d" y="%d">粗红线 = 协和线 100%%；'
             '虚线 = 等协和度线</text>' % (MLA, TOP + 27))
    P.append('<text class="fig-sub" x="%d" y="%d">两轴同取对数 ⇒ 这些虚线斜率同为 1、'
             '互相平行；离红线越远 = 偏离 100%% 越多</text>' % (MLA, TOP + 39))
    P.append('<text class="fig-title" x="%d" y="%d">'
             '右：比值空间（Wetherill 谐和线，两轴都线性）</text>' % (MLB, TOP + 14))
    P.append('<text class="fig-sub" x="%d" y="%d">粗红线 = 协和线；'
             '虚线 = 等协和度线</text>' % (MLB, TOP + 27))
    P.append('<text class="fig-sub" x="%d" y="%d">等协和度线从原点发散，'
             '但不是直线、也不平行</text>' % (MLB, TOP + 39))

    # ── 左：年龄–年龄 ──
    for t in (100, 200, 300, 500, 700, 1000, 1500, 2000, 3000):
        P.append('<line class="fig-grid" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>'
                 % (xa(t), PT, xa(t), PB_))
        P.append('<line class="fig-grid" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (MLA, ya(t), MLA + PW, ya(t)))
        P.append('<text class="fig-tick" x="%.1f" y="%d" text-anchor="middle">%d</text>'
                 % (xa(t), PB_ + 15, t))
        P.append('<text class="fig-tick" x="%d" y="%.1f" text-anchor="end">%d</text>'
                 % (MLA - 7, ya(t) + 4, t))
    P.append('<text class="fig-axis" x="%.1f" y="%d" text-anchor="middle">'
             '206Pb/238U 年龄（Ma，对数轴）</text>' % (MLA + PW / 2, PB_ + 40))
    P.append('<text class="fig-axis" transform="translate(14,%.1f) rotate(-90)" '
             'text-anchor="middle">207Pb/235U 年龄（Ma，对数轴）</text>' % (PT + PH / 2))

    # 90–110% 带：夹在两条射线之间，画成一个多边形
    poly = [(xa(AGE_AX_LO), ya(AGE_AX_LO * CONC_BAND[0] / 100)),
            (xa(AGE_AX_HI), ya(AGE_AX_HI * CONC_BAND[0] / 100)),
            (xa(AGE_AX_HI), ya(AGE_AX_HI * CONC_BAND[1] / 100)),
            (xa(AGE_AX_LO), ya(AGE_AX_LO * CONC_BAND[1] / 100))]
    P.append('<polygon class="cc-band" points="%s"/>'
             % " ".join("%.1f,%.1f" % p for p in poly))
    # 射线标签：放在**射线出图框的位置**（顶边或右边），同边的沿边错开。
    tops, rights = [], []                      # [x, y, 文本]（错开只改第一个分量）
    for c in CONC_RAYS:
        f = c / 100.0
        P.append('<line class="cc-ray" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (xa(AGE_AX_LO), ya(AGE_AX_LO * f), MLA + PW, ya(AGE_AX_HI * f)))
        if f > 1.0:                            # 从顶边出去
            tops.append([xa(AGE_AX_HI * 100 / c), PT + 12, "%g%%" % c])
        elif f < 1.0:                          # 从右边出去
            rights.append([ya(AGE_AX_HI * f) - 4, MLA + PW - 4, "%g%%" % c])
    _cascade(tops, 30.0, lambda r: r[0], lower_is_first=False)
    # 右边这组竖着排，会和顶边那排（横着排）在右上角撞在一起（90% 压住 110%，踩过）。
    # 所以先**按出框位置排好**（y 越小 = 协和度越高），再给第 i 条一个下界：
    # 既不高于顶行之下，也不比上一条近。⚠ 别把顺序交给"键相等时的稳定排序"——
    # 那样排出来的先后看不出来，实测就把 70/80/90 的顺序搞反了。
    rights.sort(key=lambda r: r[0])
    for i, r in enumerate(rights):
        r[0] = max(r[0], PT + 28.0 + i * 16.0)
    _cascade(rights, 16.0, lambda r: r[0], lower_is_first=True)
    for x, y, s in tops:
        P.append('<text class="cc-lab" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (max(x, MLA + 16), y, s))
    for y, x, s in rights:
        P.append('<text class="cc-lab" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x, y, s))
    P.append('<line class="cc-line" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (xa(AGE_AX_LO), ya(AGE_AX_LO), MLA + PW, ya(AGE_AX_HI)))

    # ── 右：比值空间 ──
    # 右侧再加一列「对应年龄」：比值轴上的数字单看没有直觉，标上年龄才知道
    # 它意味着什么（比值与年龄一一对应）。
    P.append('<text class="fig-sub" x="%d" y="%d">对应年龄</text>'
             % (MLB + PW + 8, PT - 6))
    # 刻度选点：比值轴是**线性**的，等年龄间隔在图上并不等距（小年龄全挤在左下角，
    # 200/500 Ma 的刻度和数字会叠在一起）。所以这里既挑间隔得开的一串年龄，
    # 又按**像素间距**再筛一遍：谁跟上一个挨得比 MINSEP 还近就丢掉。
    ticks, last_x, last_y = [], -1e9, 1e9
    for t in (500, 1000, 1500, 2000, 2500, 3000):
        r68 = math.expm1(LAM68 * t * 1e6)
        r75 = math.expm1(LAM35 * t * 1e6)
        if r68 > R68_AX_HI or r75 > R75_AX_HI:
            continue                  # 刻度落在轴外就不画，别贴在边上骗人
        px, py = xb(r68), yb(r75)
        if px - last_x < 46.0 or last_y - py < 15.0:
            continue
        ticks.append((t, r68, r75, px, py))
        last_x, last_y = px, py
    for t, r68, r75, px, py in ticks:
        P.append('<line class="fig-grid" x1="%.1f" y1="%d" x2="%.1f" y2="%d"/>'
                 % (px, PT, px, PB_))
        P.append('<line class="fig-grid" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (MLB, py, MLB + PW, py))
        P.append('<text class="fig-tick" x="%.1f" y="%d" text-anchor="middle">%.2f</text>'
                 % (px, PB_ + 15, r68))
        P.append('<text class="fig-tick" x="%d" y="%.1f" text-anchor="end">%.2f</text>'
                 % (MLB - 7, py + 4, r75))
        P.append('<text class="fig-sub" x="%.1f" y="%d" text-anchor="middle">%d Ma</text>'
                 % (px, PB_ + 28, t))
        P.append('<text class="fig-sub" x="%d" y="%.1f">%d Ma</text>'
                 % (MLB + PW + 8, py + 4, t))
    P.append('<text class="fig-axis" x="%.1f" y="%d" text-anchor="middle">'
             '206Pb/238U 比值（线性轴；小字 = 对应年龄）</text>'
             % (MLB + PW / 2, PB_ + 40))
    P.append('<text class="fig-axis" transform="translate(%d,%.1f) rotate(-90)" '
             'text-anchor="middle">207Pb/235U 比值（线性轴）</text>'
             % (MLB - 46, PT + PH / 2))

    # 协和曲线 + 等协和度曲线
    tt = [AGE_AX_LO * (AGE_AX_HI / AGE_AX_LO) ** (i / 240.0) for i in range(241)]
    pts = [(xb(math.expm1(LAM68 * t * 1e6)), yb(math.expm1(LAM35 * t * 1e6)))
           for t in tt]
    P.append('<polyline class="cc-curve" points="%s"/>'
             % " ".join("%.1f,%.1f" % p for p in pts))
    r_tops, r_rights = [], []
    for c in CONC_RAYS:
        f = c / 100.0
        raw, stop = [], None
        for t in tt:
            r68 = math.expm1(LAM68 * t * 1e6)
            r75 = math.expm1(LAM35 * f * t * 1e6)
            if r68 > R68_AX_HI:
                stop = "right"
                break
            if r75 > R75_AX_HI:
                stop = "top"
                break
            raw.append((xb(r68), yb(r75)))
        if len(raw) < 2:
            continue
        P.append('<polyline class="cc-ray2" points="%s"/>'
                 % " ".join("%.1f,%.1f" % p for p in raw))
        # 标签放在**射线真正出框的那个点**上。⚠ 不要用"最后一点的比值是否
        # 接近轴上限"来猜出框边：采样是有步长的（这 241 个点每步约 1.5%），
        # 最后一点可能离轴上限还差一整步，于是**所有**射线都被判成"从右边
        # 出去"，标签全堆在右上角压成一团（踩过）。这里直接解析求交点。
        if stop == "top":
            te = math.log1p(R75_AX_HI) / (LAM35 * f) / 1e6      # 使 r75 = 轴上限
            r_tops.append([xb(math.expm1(LAM68 * te * 1e6)), PT + 12, "%g%%" % c])
        else:                                                    # 撞的是 r68 上限
            te = math.log1p(R68_AX_HI) / LAM68 / 1e6
            r_rights.append([yb(math.expm1(LAM35 * f * te * 1e6)) - 4,
                             MLB + PW - 4, "%g%%" % c])
    _cascade(r_tops, 30.0, lambda r: r[0], lower_is_first=False)
    r_rights.sort(key=lambda r: r[0])            # 同上：顺序要自己排，不靠稳定排序
    for i, r in enumerate(r_rights):
        r[0] = max(r[0], PT + 28.0 + i * 16.0)
    _cascade(r_rights, 16.0, lambda r: r[0], lower_is_first=True)
    for x, y, s in r_tops:
        P.append('<text class="cc-lab" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (max(x, MLB + 16), y, s))
    for y, x, s in r_rights:
        P.append('<text class="cc-lab" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x, y, s))

    # ── 两图共用的数据云 ──
    for panel in ("age", "ratio"):
        for a68, a75 in lib:
            got = _to_panel(panel, a68, a75, xa, ya, xb, yb)
            if got is None:
                continue
            P.append('<circle class="cc-bg" cx="%.1f" cy="%.1f" r="1.5"/>' % got)
        for r in rows:
            for d in r["domains"]:
                got = _to_panel(panel, d["age68"], d["age75"], xa, ya, xb, yb)
                if got is None:
                    continue
                P.append('<circle class="fig-dot fig-halo %s" cx="%.1f" cy="%.1f" r="4"/>'
                         % ((_conc_class(d["concordance"]),) + got))

    P.append('<rect class="fig-frame" x="%d" y="%d" width="%d" height="%d"/>' % (MLA, PT, PW, PH))
    P.append('<rect class="fig-frame" x="%d" y="%d" width="%d" height="%d"/>' % (MLB, PT, PW, PH))
    P.append('</svg>')
    return "\n".join(P), W, H


def concordia_caption(rows, lib) -> str:
    """图注 + 一句"到底代表什么"的算例，数字全部现算。"""
    n_band = sum(1 for a68, a75 in lib if CONC_BAND[0] <= a75 / a68 * 100 <= CONC_BAND[1])
    # 左图轴内画得出的灰点数（越界的点在 svg_concordia 里被跳过）
    n_in_age = sum(1 for a68, a75 in lib
                   if AGE_AX_LO <= a68 <= AGE_AX_HI and AGE_AX_LO <= a75 <= AGE_AX_HI)
    obs = [d for r in rows for d in r["domains"]]
    n_hi = sum(1 for d in obs if d["concordance"] >= 150.0)
    lo1, hi1 = CONC_BAND
    return ('<figcaption><b>协和度就是一个比值</b>：'
            '<code>age(207Pb/235U) ÷ age(206Pb/238U)</code>。'
            '左图横轴取 206Pb/238U 年龄、纵轴取 207Pb/235U 年龄，'
            '于是这个比值就是"纵值 ÷ 横值"这一个数。'
            '100%% 是那条 45° 的<b>协和线</b>（两个衰变体系给出同一个年龄）；'
            '90–110%% 是它两侧的灰带。落在红线<b>上方</b>说明 207Pb 侧偏老'
            '（普通铅未扣净、或 <sup>206</sup>PbH<sup>+</sup> 干扰）；'
            '<b>下方</b>说明铅丢失或 206Pb 过量。'
            '⚠ 左图两轴都取了对数（否则 100 Ma 与 3000 Ma 没法画进同一张图），'
            '所以等协和度线在图上是一族<b>斜率同为 1 的平行线</b>、而不是扇形 ——'
            '那是坐标轴的功劳，不是协和度的性质变了；'
            '离红线越远，就代表偏离 100%% 越多。'
            '右图是同一批点在比值空间里的样子（教科书上的 Wetherill 谐和线）：'
            '协和线是一条曲线，等协和度线在这里<b>从原点发散</b>，'
            '但既不是直线、也不互相平行（两条衰变链的 λ 不同）。'
            '<b>全库 %d 个正式年龄域里有 %.0f%% 落在 90–110%% 的带内</b>'
            '（左图只画得出其中落在 100–3500 Ma 框内的 %d 个）。'
            '彩标的 %d 个域里有 %d 个 ≥ 150%%，全部在带外上方。'
            '算例：206/238 给 460 Ma 时，207/235 给 %.0f Ma 就正好是 %g%%，'
            '给 %.0f Ma 就是 %g%%。</figcaption>'
            % (len(lib), 100.0 * n_band / len(lib) if lib else float("nan"),
               n_in_age, len(obs), n_hi, 460 * lo1 / 100, lo1, 460 * hi1 / 100, hi1))


# ─────────────────────────────────────────────────────────────────────────────
# 四、本地实名版
# ─────────────────────────────────────────────────────────────────────────────

def write_local(rows, bands, total) -> list:
    svg, W, H = svg_figure(rows, bands, row_labels(rows, True), inline_style=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = OUT_DIR / "extreme_spans_local.svg"
    svg_path.write_text(svg, encoding="utf-8")

    csv_path = OUT_DIR / "extreme_spans_labels.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["G码", "批次", "区组", "序号", "样品", "跨度Ma", "域", "tau0", "tau1",
                    "窗数", "年龄206_238_Ma", "2σ", "年龄207_235_Ma", "2σ",
                    "域级协和度_pct", "域内MSWD", "ThU"])
        for r in rows:
            for d in r["domains"]:
                w.writerow([r["code"], r["批次"], r["区组"], r["序号"], r["样品"],
                            "%.1f" % r["span"], d["domain"], "%.4f" % d["tau0"],
                            "%.4f" % d["tau1"], d["n_win"], "%.2f" % d["age68"],
                            "%.2f" % d["s2_68"], "%.2f" % d["age75"], "%.2f" % d["s2_75"],
                            "%.2f" % d["concordance"], "%.3f" % d["mswd"], "%.4f" % d["th_u"]])

    png_path = OUT_DIR / "extreme_spans_local.png"
    ok = _render_png(svg, png_path, W, H)
    return [svg_path, csv_path] + ([png_path] if ok else [])


def _render_png(svg: str, out: pathlib.Path, w: int, h: int) -> bool:
    """
    用无头 Edge 把 SVG 截成 PNG（2 倍尺寸）。Edge 不在就跳过 —— 这不是失败。

    三处细节都是实测踩出来的：
      · `--user-data-dir` 指到临时目录：不指就用默认 profile，
        本机正开着 Edge 时无头实例直接退 21；
      · **不要用 `--force-device-scale-factor`** —— 它与 `--window-size` 的
        换算关系在 headless=new 下不是 1:1，实测截出来四周一片空白。
        改用"把 SVG 的 CSS 尺寸设成 2 倍 + 窗口也开 2 倍"：图是矢量，
        放大不损失清晰度。
      · ★ **`--window-size` ≠ 视口尺寸**：headless=new 会扣掉窗口边框与
        工具栏（本机实测宽少 16 px、高少 88 px）。照窗口尺寸设，图的**底边
        和右边会被整条裁掉** —— 第一版就是这么丢掉 x 轴刻度与轴标题的
        （症状：轴标题明明写在 SVG 里，图上就是没有）。
        所以这里改成：**窗口开得比图大一圈，再用 Pillow 按 SVG 的
        精确像素尺寸从左上角裁齐**。这样输出尺寸不再依赖浏览器版本。
    """
    if not EDGE.exists():
        print("[提示] 找不到 Edge，跳过 PNG：%s" % EDGE, file=sys.stderr)
        return False
    sc = 2
    want_w, want_h = w * sc, h * sc
    pad = 240                      # 留够窗口 chrome 的余量（实测 16 / 88）
    html_path = OUT_DIR / "_local_preview.html"
    html_path.write_text(
        "<!DOCTYPE html><html><head><meta charset='utf-8'><style>"
        "html,body{margin:0;padding:0;background:#fff}"
        "svg{display:block;width:%dpx;height:%dpx}"
        "</style></head><body>%s</body></html>" % (want_w, want_h, svg),
        encoding="utf-8")
    cmd = [str(EDGE), "--headless=new", "--disable-gpu", "--hide-scrollbars",
           "--user-data-dir=%s" % (OUT_DIR / "_edge_tmp"),
           "--screenshot=%s" % out,
           "--window-size=%d,%d" % (want_w + pad, want_h + pad),
           html_path.as_uri()]
    try:
        r = subprocess.run(cmd, capture_output=True, timeout=120)
    except Exception as exc:                                           # noqa: BLE001
        print("[提示] PNG 截图失败：%s" % exc, file=sys.stderr)
        return False
    if r.returncode != 0 or not out.exists():
        print("[提示] PNG 截图返回 %s（SVG 已生成，不影响使用）" % r.returncode,
              file=sys.stderr)
        return False
    try:
        from PIL import Image
        with Image.open(out) as im:
            got = im.size
            if got[0] < want_w or got[1] < want_h:
                print("[提示] 截图只有 %dx%d，比图小的 %dx%d 还小，不做裁切"
                      % (got[0], got[1], want_w, want_h), file=sys.stderr)
                return True
            im.crop((0, 0, want_w, want_h)).save(out)
    except Exception as exc:                                           # noqa: BLE001
        print("[提示] PNG 裁切失败：%s" % exc, file=sys.stderr)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# 五、全库背景重算
# ─────────────────────────────────────────────────────────────────────────────
def listing_map() -> dict:
    """`批次清单.csv` → {批次: {"文件夹": …, "参考值口径": …}}。文件夹是反斜杠拼的。"""
    out: dict[str, dict] = {}
    p = ROOT / "汇总库" / "批次清单.csv"
    if not p.exists():
        return out
    with open(p, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            out[r["批次"]] = dict(文件夹=(r.get("文件夹") or "").replace("\\", "/"),
                                  参考值口径=(r.get("参考值口径") or "").strip())
    return out


def cmd_scan(db: pathlib.Path, out: pathlib.Path) -> int:
    """全库（全部样品测点）逐域算域级协和度 → CSV。这是右图灰底的来源。"""
    import sqlite3
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    batches = [r["批次"] for r in con.execute("select 批次 from 批次")]
    want: dict[str, set] = defaultdict(set)
    for r in con.execute("select 批次, 序号 from 分析"):
        want[r["批次"]].add(r["序号"])
    con.close()

    meta = listing_map()
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for batch in batches:
        info = meta.get(batch, {})
        preset = info.get("参考值口径") or "horstwood2016"
        data_dir = ROOT / "汇总库" / info.get("文件夹", batch) / "原始数据"
        if not data_dir.exists():
            print("[跳过] %s：找不到 %s" % (batch, data_dir), file=sys.stderr)
            continue
        try:
            res = recompute_batch(data_dir, preset, want[batch])
        except Exception as exc:                                       # noqa: BLE001
            print("[失败] %s：%s" % (batch, exc), file=sys.stderr)
            continue
        for idx, a in res.items():
            for d in a["domains"]:
                rows.append(dict(批次=batch, 序号=idx, **d))
        print("[扫描] %s → %d 点" % (batch, len(res)), file=sys.stderr)

    if not rows:
        print("没有扫到任何域。", file=sys.stderr)
        return 1
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("写出 %d 行 → %s" % (len(rows), out))
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# 六、命令行
# ─────────────────────────────────────────────────────────────────────────────
def cmd_build(db: pathlib.Path, span: float, page: pathlib.Path,
              local_only: bool) -> int:
    if not LIB_CSV.exists():
        print("缺少全库背景 %s —— 先跑：python %s scan"
              % (LIB_CSV, pathlib.Path(__file__).name), file=sys.stderr)
        return 2
    bands, total = library_bands(LIB_CSV)
    spots, audited = collect_groups(db, span)
    rows = domain_rows(spots, audited)

    print()
    print("=" * 108)
    print("跨度 > %g Ma 的测点复核（只算正式年龄域；协和度 = 域级 207/235 ÷ 206/238）"
          % span)
    print("=" * 108)
    for r in rows:
        print("%s  %-16s %-3d 跨度 %7.1f Ma  域数 %d" % (
            r["code"], "%s#%d" % (r["样品"], r["序号"]), r["序号"], r["span"],
            len(r["domains"])))
        for d in r["domains"]:
            print("      %-3s τ %.2f–%.2f  %2d 窗  %8.1f ±%-5.1f Ma  "
                  "207/235 %8.1f ±%-5.1f  协和度 %7.1f%%  MSWD %5.2f  Th/U %.3f"
                  % (d["domain"], d["tau0"], d["tau1"], d["n_win"], d["age68"],
                     d["s2_68"], d["age75"], d["s2_75"], d["concordance"],
                     d["mswd"], d["th_u"]))
    print("-" * 108)
    print("全库背景：%d 个正式年龄域，%d 箱" % (total, len(bands)))
    for b in bands:
        print("   τ≈%.2f  n=%4d  中位 %7.1f  P25 %7.1f  P75 %7.1f  P95 %7.1f"
              "  >130%% 占 %4.1f%%" % (b["tau"], b["n"], b["p50"], b["p25"],
                                       b["p75"], b["p95"], b["frac_hi"]))

    written = write_local(rows, bands, total)
    print("\n本地实名版：")
    for p in written:
        print("   %s" % p)

    if not local_only:
        frag = build_fragment(rows, bands, total, show_real=False)
        if not page.exists():
            print("找不到页面 %s（先建页面骨架）" % page, file=sys.stderr)
            return 2
        inject(frag, page)
        print("\n线上脱敏版：已注入 %s" % page)

        # 协和度对照图：**两份图完全一样**（图中没有测点标签），
        # 所以本地只多存一份 SVG/PNG，线上直接注入同一份。
        lib = library_points(LIB_CSV)
        conc_svg, cw, ch = svg_concordia(rows, lib, inline_style=True)
        (OUT_DIR / "concordia_local.svg").write_text(conc_svg, encoding="utf-8")
        if _render_png(conc_svg, OUT_DIR / "concordia_local.png", cw, ch):
            print("   协和度对照图：%s" % (OUT_DIR / "concordia_local.png"))
        inject(concordia_fragment(rows, lib), page, CONC_BEGIN, CONC_END)
        print("   协和度对照图已注入（%d 个背景域 + %d 个彩点）"
              % (len(lib), sum(len(r["domains"]) for r in rows)))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[2],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    p1 = sub.add_parser("scan", help="重算全库域级协和度 → CSV（右图灰底）")
    p1.add_argument("--db", default=str(DB_DEFAULT))
    p1.add_argument("--out", default=str(LIB_CSV))
    p2 = sub.add_parser("build", help="出两份图 + 生成在线页")
    p2.add_argument("--db", default=str(DB_DEFAULT))
    p2.add_argument("--span", type=float, default=SPAN_DEFAULT)
    p2.add_argument("--page", default=str(PAGE))
    p2.add_argument("--local-only", action="store_true", help="只出本地实名版，不动页面")
    args = ap.parse_args(argv)

    if args.cmd == "scan":
        return cmd_scan(pathlib.Path(args.db), pathlib.Path(args.out))
    if args.cmd == "build":
        return cmd_build(pathlib.Path(args.db), args.span,
                         pathlib.Path(args.page), args.local_only)
    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
