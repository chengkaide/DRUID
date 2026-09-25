# -*- coding: utf-8 -*-
"""
重新生成落地页 `docs/index.html` 里那张「分域前 / 分域后」对比图，并就地注入。

    python tools/gen_docs_split_figure.py

画的是示例批次 `examples/EX2022A` 里三个**真实测点**，对应三种剖面形态：

    ① S43  均一      —— 整段与「只有一个年龄」相容。域表里根本没有它，不需要分域
    ② S34  核边分明  —— 两段平台差 45.6 Ma，一个窗口都没被剥掉，一分为二正好
    ③ S18  复杂变化  —— 三段平台 + 一片 8 个窗口的过渡带，只拆得动一部分

为什么要有这个脚本
------------------
图上每一个点、每一条线都是**真实数据**：坐标取自「剖面窗口」表，域均值取自
「深度剖面域」表，灰虚线取自「不分域年龄」表。既然是算出来的，就必须能从仓库里
的数据重新算一遍 —— 否则它和手抄的数字没有区别，改了算法也没人知道它过期了。

只用**公开输出**（`run_batch` 返回的四张表），不碰内部结构：
这样图上的点与文档里引用的数字保证来自同一次运行。

`tests/check_example_batch.py` 会反过来核对这张图里的数字（按测点分块抓取），
所以算法一变、图没跟上，检查当场就会红。

踩过的坑（改这个脚本时别再犯）
------------------------------
1. SVG 根元素必须带 `class="fig"`。页面里所有配色规则都写成 `.fig .xxx`，
   漏了这个 class 会得到一张"尺寸正确、内容完全空白"的图。
2. 图例必须放在坐标系**外面**（这里输出成 HTML 放在图下方）——
   画进坐标系里一定会压住某个角上的数据（实测盖掉过右上角的三个高年龄窗口）。
3. `tau` 是窗口**起点**，而域表的 `tau` 写成"首个窗口起点-末个窗口起点"，
   所以归域要按半个步长做容差；直接用区间过滤会漏掉每个域的最后一个窗口。
4. 页面配色全部走 CSS 变量，深浅两套主题共用同一份 SVG；
   所以这里一个颜色都不写死，只给 class。
"""
from __future__ import annotations

import argparse
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from druid.workflow import BatchConfig, run_batch   # noqa: E402

BATCH = ROOT / "examples" / "EX2022A"
HTML = ROOT / "docs" / "index.html"

#: 三个示例测点。`label` 会出现在图上与图例里，改它要同步改回归守护的正则。
CASES = (
    {"spot": "S43", "num": "①", "label": "均一", "head": "整段就是答案"},
    {"spot": "S34", "num": "②", "label": "核边分明", "head": "一分为二正好"},
    {"spot": "S18", "num": "③", "label": "复杂变化", "head": "只拆得动一部分"},
)

# ── 版面 ──────────────────────────────────────────────────────────────────
W = 900
ML, MR = 78, 28
PW = W - ML - MR
PH = 150          # 单个面板的绘图区高度
HEAD = 26         # 面板标题行
TAIL = 22         # 面板自己的 x 轴刻度行
BLOCK = HEAD + PH + TAIL
GAP = 26
TOP = 14
BOTTOM = 30
H = TOP + len(CASES) * BLOCK + (len(CASES) - 1) * GAP + BOTTOM

BEGIN = "<!-- ==== 分域前后对比图：由 tools/gen_docs_split_figure.py 生成，请勿手改 ==== -->"
LEGACY = "<!-- 这张图由 build_tmp/"          # 迁移用，兼容最早的标记


def nice_ticks(lo: float, hi: float, want: int = 4) -> list[float]:
    """取"好看的整数"刻度，而不是把 y0~y1 等分 —— 等分会给出 417/441/465 这种数。"""
    span = hi - lo
    mag = 10 ** math.floor(math.log10(span / want))
    for m in (1, 2, 2.5, 5, 10):
        s = m * mag
        if span / s <= want + 1:
            break
    v = math.ceil(lo / s) * s
    out = []
    while v <= hi + 1e-9:
        out.append(round(v, 6))
        v += s
    return out


def collect(res, spot: str) -> dict:
    """把一个测点的窗口、域、整段口径全部抽出来，并核对窗口数守恒。"""
    win = res.windows[res.windows["Sample"] == spot].reset_index(drop=True)
    dom = res.domains[res.domains["样品"] == spot].sort_values("tau").reset_index(drop=True)
    ov = res.overall[res.overall["样品"] == spot]
    if not len(win) or not len(ov):
        raise SystemExit(f"示例批次里找不到 {spot} 的剖面窗口或整段口径，图没法画")
    ov = ov.iloc[0]

    tau = win["Tau"].to_numpy(float)
    age = win["Age68"].to_numpy(float)
    sig = win["Age68_1s"].to_numpy(float)
    step = float(np.median(np.diff(tau)))
    tol = step * 0.55

    # 逐窗口归域：域表的 tau 上界是**末个窗口的起点**，所以要带半个步长的容差
    belong = np.full(tau.size, "", dtype=object)
    bands = []
    for _, r in dom.iterrows():
        a, b = (float(x) for x in str(r["tau"]).split("-"))
        belong[(tau >= a - tol) & (tau <= b + tol)] = str(r["域"])
        bands.append({"name": str(r["域"]), "a": a, "b": min(b + step, 1.0),
                      "age": float(r["年龄_Ma"]), "n": int(r["n_win"]),
                      "mswd": float(r["MSWD"]), "s2": float(r["s2_Ma"])})

    n_win = int(ov["n_win"])
    n_off = int((belong == "").sum())
    for bd in bands:
        got = int((belong == bd["name"]).sum())
        assert got == bd["n"], f"{spot}/{bd['name']} 窗口数对不上：算得 {got}，表里 {bd['n']}"
    if not bands:
        # 均一测点：域表里没有行，全部窗口都属于同一个年龄，不存在"被剥掉"这回事
        assert n_off == n_win, f"{spot} 域表里没有行，却有窗口落进了某个域"
        belong[:] = "D1"
        n_off = 0
    else:
        assert n_off == n_win - sum(bd["n"] for bd in bands), f"{spot} 未归域窗口数对不上"

    return {
        "spot": spot, "tau": tau, "age": age, "sig": sig, "belong": belong,
        "bands": bands, "n_win": n_win, "n_off": n_off,
        "whole": float(ov["年龄_Ma"]), "whole_mswd": float(ov["MSWD"]),
        "crit": float(ov["相容上限"]), "whole_s2": float(ov["s2_Ma"]),
        "bulk": float(ov["整段积分年龄_Ma"]),
        "main_age": float(ov["主域年龄_Ma"]), "main_dom": str(ov["主域"]),
        "d_age": float(ov["Δ年龄_Ma"]), "d_pct": float(ov["Δ年龄_pct"]),
        "verdict": str(ov["判定"]), "structure": str(ov["深度结构"]),
    }


def panel(i: int, case: dict, d: dict) -> list[str]:
    """画第 i 个面板，返回 SVG 片段。"""
    top = TOP + i * (BLOCK + GAP)
    pt = top + HEAD                     # 绘图区上边
    pb = pt + PH                        # 绘图区下边

    lo = float(np.nanmin(d["age"] - d["sig"]))
    hi = float(np.nanmax(d["age"] + d["sig"]))
    pad = (hi - lo) * 0.10
    y0v, y1v = lo - pad, hi + pad
    ticks = nice_ticks(y0v, y1v, want=4)

    def X(t: float) -> float:
        return ML + (t - 0.0) / 1.0 * PW

    def Y(a: float) -> float:
        return pt + (y1v - a) / (y1v - y0v) * PH

    n_dom = len(d["bands"])
    if n_dom == 0:
        verdict = f"整段 MSWD {d['whole_mswd']:.2f} / 上限 {d['crit']:.2f} → 整段常数"
    else:
        verdict = f"整段 MSWD {d['whole_mswd']:.2f} / 上限 {d['crit']:.2f} → 整段非常数"

    P: list[str] = [f'<g class="case" data-spot="{d["spot"]}">']
    # 面板标题行：左边是形态，右边是整段判据
    P.append(f'<text class="fig-cap" x="{ML}" y="{pt-11:.0f}">'
             f'<tspan class="fig-num">{case["num"]}</tspan>  {d["spot"]} · {case["label"]}'
             f'<tspan class="fig-cap-sub">　{case["head"]}</tspan></text>')
    P.append(f'<text class="fig-sub" x="{ML+PW}" y="{pt-11:.0f}" text-anchor="end">{verdict}</text>')

    for a in ticks:
        P.append(f'<line class="fig-grid" x1="{ML}" y1="{Y(a):.1f}" x2="{ML+PW}" y2="{Y(a):.1f}"/>')
        P.append(f'<text class="fig-tick" x="{ML-12}" y="{Y(a)+4:.1f}" text-anchor="end">{a:.0f}</text>')
    if i == 1:                          # 纵轴标题只画一次，免得三行重复
        P.append(f'<text class="fig-axis" transform="translate(22,{pt+PH/2:.0f}) rotate(-90)" '
                 f'text-anchor="middle">年龄 (Ma)</text>')

    for bd in d["bands"]:               # 域底纹画在数据点下面
        P.append(f'<rect class="band-{bd["name"]}" x="{X(bd["a"]):.1f}" y="{pt}" '
                 f'width="{X(bd["b"])-X(bd["a"]):.1f}" height="{PH}"/>')

    for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        P.append(f'<line class="fig-grid" x1="{X(t):.1f}" y1="{pt}" x2="{X(t):.1f}" y2="{pb}"/>')
        P.append(f'<text class="fig-tick" x="{X(t):.1f}" y="{pb+17}" text-anchor="middle">'
                 f'{t:g}</text>')

    P.append(f'<clipPath id="cpCase{i}"><rect x="{ML}" y="{pt}" width="{PW}" height="{PH}"/></clipPath>')
    P.append(f'<g clip-path="url(#cpCase{i})">')
    for t, a, s, b in zip(d["tau"], d["age"], d["sig"], d["belong"]):
        if not np.isfinite(a):
            continue
        x, y = X(t), Y(a)
        cls = f'dot-{b}' if b else "fig-dot-off"
        P.append(f'<line class="{"fig-ebar" if b else "fig-ebar-off"}" '
                 f'x1="{x:.1f}" y1="{Y(a-s):.1f}" x2="{x:.1f}" y2="{Y(a+s):.1f}"/>')
        P.append(f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="3.4"/>')
    P.append('</g>')

    # 灰虚线＝整段不分域口径，三个面板都画：①里它就是答案，②③里它不是任何一个域
    P.append(f'<line class="fig-whole" x1="{ML}" y1="{Y(d["whole"]):.1f}" x2="{ML+PW}" '
             f'y2="{Y(d["whole"]):.1f}" stroke-dasharray="7 4"/>')
    for bd in d["bands"]:
        ym = Y(bd["age"])
        P.append(f'<line class="solid-{bd["name"]}" x1="{X(bd["a"]):.1f}" y1="{ym:.1f}" '
                 f'x2="{X(bd["b"]):.1f}" y2="{ym:.1f}"/>')
        P.append(f'<text class="txt-{bd["name"]}" x="{X((bd["a"]+bd["b"])/2):.1f}" y="{ym-9:.1f}" '
                 f'text-anchor="middle">{bd["name"]} {bd["age"]:.1f}</text>')
    if d["n_off"]:
        edges = [bd for bd in d["bands"]]
        mid = X((edges[0]["b"] + edges[-1]["a"]) / 2) if n_dom > 1 else X(0.5)
        P.append(f'<text class="fig-off-note" x="{mid:.1f}" y="{pt+17}" text-anchor="middle">'
                 f'{d["n_off"]} 个混合/过渡窗口（未进入任何年龄域）</text>')
    P.append(f'<rect class="fig-frame" x="{ML}" y="{pt}" width="{PW}" height="{PH}"/>')
    if i == len(CASES) - 1:
        P.append(f'<text class="fig-axis" x="{ML+PW/2:.0f}" y="{H-10}" text-anchor="middle">'
                 f'剥蚀时间 τ（0 = 激光开，1 = 激光关）　→ 坑深方向</text>')
    P.append('</g>')
    return P


def build() -> tuple[str, list[dict]]:
    res = run_batch(BatchConfig(data_dir=str(BATCH), plot=False))
    data = []
    for case in CASES:
        d = collect(res, case["spot"])
        d["num"], d["label"], d["head"] = case["num"], case["label"], case["head"]
        data.append(d)

    P: list[str] = ['<div class="figscroll">']
    P.append(f'<svg class="fig" viewBox="0 0 {W} {H}" role="img" '
             f'aria-labelledby="figSplitT figSplitD">')
    P.append('<title id="figSplitT">同一个示例批次里的三种剖面形态：均一 / 核边分明 / 复杂变化</title>')
    P.append('<desc id="figSplitD">三个真实测点，每个点是一个 4 秒滑窗，误差棒是 1σ。'
             '灰虚线是整段不分域的加权平均：在 ① 里它就是答案，在 ②③ 里它不对应任何一个窗口。'
             '彩色实线是两个或三个年龄域的均值，空心圈是没通过域内相容性检验、'
             '不进入任何年龄域的窗口。</desc>')
    P.append(f'<rect x="0" y="0" width="{W}" height="{H}" rx="14" class="fig-bg"/>')

    rows = []
    for i, (case, d) in enumerate(zip(CASES, data)):
        P.extend(panel(i, case, d))
        rows.append({"case": case, **{k: d[k] for k in
                    ("spot", "n_win", "n_off", "bands", "whole", "whole_mswd", "crit",
                     "bulk", "main_age", "main_dom", "d_age", "d_pct", "verdict",
                     "structure", "whole_s2")}})
    P.append('</svg>')
    P.append('</div>')
    P.append('<p class="figscroll-hint">图较宽，可左右滑动查看完整的坐标轴</p>')

    # ── 图例（HTML）：每行一个测点，做成回归守护能按行抓取的形状 ──
    P.append('<ul class="figlegend">')
    P.append('<li class="key"><b>怎么读</b>'
             '<span class="k"><i class="sw-whole"></i>整段不分域：把全部窗口当成唯一一个域</span>'
             '<span class="k"><i class="sw-d1"></i><i class="sw-d2"></i><i class="sw-d3"></i>'
             '各年龄域均值（D1 / D2 / D3）</span>'
             '<span class="k"><i class="sw-off"></i>没通过域内相容性检验、不进入任何域</span>'
             '<span class="k">Δ 一律是「整段 − 主域」</span></li>')
    for case, d in zip(CASES, data):
        tag = f'{d["spot"]} {case["label"]}'
        if not d["bands"]:
            em = (f'共 {d["n_win"]} 个滑窗。整段 MSWD {d["whole_mswd"]:.2f} &lt; 判据上限 '
                  f'{d["crit"]:.2f} → <b>整段就是常数</b>。域表里根本没有这个测点：'
                  f'没有域可拆，灰虚线就是答案。整段积分（另一种算法）给 {d["bulk"]:.1f} Ma，'
                  f'与窗口加权差 {abs(d["bulk"]-d["whole"])/d["whole"]*100:.2f}%。')
        else:
            segs = " · ".join(f'{b["name"]} {b["age"]:.1f} Ma（{b["n"]} 窗，MSWD {b["mswd"]:.2f}）'
                              for b in d["bands"])
            tail = (f'{d["n_off"]} 个窗口没通过域内相容性检验，被分域步骤剥掉。'
                    if d["n_off"] else '没有一个窗口被剥掉。')
            em = (f'共 {d["n_win"]} 个滑窗。{len(d["bands"])} 个域：{segs}。{tail}'
                  f'整段 MSWD {d["whole_mswd"]:.2f} &gt; 上限 {d["crit"]:.2f} → 整段非常数，'
                  f'整段 {d["whole"]:.1f} − 主域 {d["main_age"]:.1f} = '
                  f'<b>{d["d_age"]:+.2f} Ma（{d["d_pct"]:+.2f}%）</b>。')
        P.append(f'<li class="case"><b>{case["num"]} {tag} · 整段 {d["whole"]:.1f} Ma</b>'
                 f'<em>{em}</em></li>')
    P.append('</ul>')
    P.append(figcaption(data))
    return "\n".join(P), rows


def figcaption(data: list[dict]) -> str:
    """图注。数字同样全部现算 —— 它和图上那三张图是同一批数。"""
    c = []
    c.append('<figcaption>三个测点都取自示例批次 <code>examples/EX2022A</code>，'
             '上面的数字可以自行复现。')
    for d in data:
        if not d["bands"]:
            c.append(f'<b>{d["spot"]}</b> 整段 MSWD {d["whole_mswd"]:.2f} 低于判据上限 '
                     f'{d["crit"]:.2f}，<b>整段就是常数</b>：域表里根本没有它，'
                     f'所以这里没有彩色实线，灰虚线就是答案。'
                     f'整段积分（另一种算法）给 {d["bulk"]:.1f} Ma，'
                     f'与窗口加权差 {abs(d["bulk"]-d["whole"])/d["whole"]*100:.2f}%。')
        else:
            ages = " / ".join(f'{b["age"]:.1f}' for b in d["bands"])
            aa = d["bands"]
            spread = max(b["age"] for b in aa) - min(b["age"] for b in aa)
            up = all(aa[i]["age"] < aa[i+1]["age"] for i in range(len(aa)-1))
            dn = all(aa[i]["age"] > aa[i+1]["age"] for i in range(len(aa)-1))
            if up:
                trend = (f'沿坑深年龄<b>递增</b>：浅部 {aa[0]["age"]:.1f} → '
                         f'深部 {aa[-1]["age"]:.1f} Ma，是"核老边新"的常序')
            elif dn:
                trend = (f'沿坑深年龄<b>递减</b>：浅部 {aa[0]["age"]:.1f} → '
                         f'深部 {aa[-1]["age"]:.1f} Ma')
            else:
                trend = ('沿坑深年龄<b>不单调</b>（' +
                         " → ".join(f'{b["age"]:.1f}' for b in aa) +
                         ' Ma）—— 不能用一个"核"一个"边"解释')
            if d["n_off"]:
                tail = (f'中间 <b>{d["n_off"]} 个空心点</b>没通过域内相容性检验，'
                        f'被分域步骤剥掉了')
            else:
                tail = '<b>一个窗口都没被剥掉</b>'
            extra = ""
            if len(aa) >= 3 and abs(aa[-1]["age"] - aa[0]["age"]) < 5:
                extra = (f'　注意 {aa[0]["name"]} 与 {aa[-1]["name"]} 的年龄几乎相同'
                         f'（{aa[0]["age"]:.1f} 与 {aa[-1]["age"]:.1f} Ma）却是两个域 ——'
                         f'「域」是<b>深度上连续的一段</b>，不是"年龄相同的所有窗口"。')
            c.append(f'<b>{d["spot"]}</b> 分为 <b>{len(aa)} 段</b>平台 {ages} Ma，'
                     f'域均值最大相差 <b>{spread:.1f} Ma（{spread/aa[0]["age"]*100:.1f}%）</b>；'
                     f'{trend}，{tail}。整段不分域给 {d["whole"]:.1f} Ma，'
                     f'它<b>不等于任何一段</b>。{extra}')
    flat = "、".join(d["num"] for d in data if not d["bands"])
    spl = "、".join(d["num"] for d in data if d["bands"])
    if flat and spl:
        c.append(f'{flat} 的灰虚线<b>就是年龄</b>（整段只有一个域，分域无事可做）；'
                 f'{spl} 的灰虚线只是在回答"不做分域会得到什么"，<b>不是任何一个年龄</b>。')
    c.append('</figcaption>')
    return "\n      ".join(c)


def inject(frag: str, html_path: pathlib.Path) -> None:
    """把生成的片段就地替换进页面。

    区间从 BEGIN 标记开始，到 `</ul>`（图例）为止；如果紧跟其后还有一个
    `<figcaption>`，一并吃掉 —— 图注里的数字也是现算的，留着旧的会自相矛盾。
    """
    html = html_path.read_text(encoding="utf-8")
    i = html.find(BEGIN)
    if i == -1:
        i = html.find(LEGACY)
        if i == -1:
            raise SystemExit(f"{html_path} 里找不到图的标记（既不是新版也不是旧版）")
    tail = html[i:]
    k = tail.find('<ul class="figlegend">')
    if k == -1:
        j = i + tail.index("</svg>") + len("</svg>")          # 更早的版本没有图例
    else:
        j = i + tail.index("</ul>", k) + len("</ul>")
        ahead = html[j:]
        s = len(ahead) - len(ahead.lstrip())                  # 跳过空白后紧跟图注？
        if ahead[s:].startswith("<figcaption>"):
            j += s + ahead[s:].index("</figcaption>") + len("</figcaption>")
    assert html[i:j].count("<svg") == 1, "待替换区间里有多个 svg"
    assert html[i:j].count("<figcaption") <= 1, "待替换区间里有多个图注"
    body = BEGIN + "\n" + "\n".join("      " + ln for ln in frag.splitlines())
    html_path.write_text(html[:i] + body + html[j:], encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--html", default=str(HTML), help="注入目标（草稿时可以指向临时副本）")
    ap.add_argument("--stdout", action="store_true", help="只把片段打到标准输出，不写文件")
    args = ap.parse_args()

    frag, rows = build()
    if args.stdout:
        print(frag)
        return 0
    inject(frag, pathlib.Path(args.html))

    print(f"示例批次 {BATCH.name}，三个剖面形态：")
    for r in rows:
        doms = " | ".join(f'{b["name"]} {b["age"]:.1f}(n={b["n"]}, MSWD {b["mswd"]:.2f})'
                          for b in r["bands"]) or "（域表里没有这个测点）"
        print(f'  {r["spot"]} {r["case"]["label"]:<6s} 窗口 {r["n_win"]:>2d} 未归域 {r["n_off"]:>2d} | '
              f'整段 {r["whole"]:.4f} MSWD {r["whole_mswd"]:.2f}/{r["crit"]:.2f} | {doms}')
    print(f"已写入 {args.html}（{len(frag)} 字节的图 + 图例）")
    print("别忘了跑 `python tests/check_example_batch.py`：它会核对这张图里的数字。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
