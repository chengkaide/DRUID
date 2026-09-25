# -*- coding: utf-8 -*-
"""
重新生成落地页 `docs/index.html` 里那张「分域前 / 分域后」对比图，并就地注入。

    python tools/gen_docs_split_figure.py

为什么要有这个脚本
------------------
那张图上的每一个点都是**真实数据**：示例批次 `examples/EX2022A` 里 S01 测点的
35 个滑窗，坐标来自「剖面窗口」表，两条域均值来自「深度剖面域」表，灰虚线来自
「不分域年龄」表。既然是算出来的，就必须能从仓库里的数据重新算一遍 ——
否则它和手抄的数字没有区别，改了算法也没人知道它过期了。

只用**公开输出**（run_batch 返回的四张表），不碰内部结构：
这样图上的点与文档里引用的数字保证来自同一次运行。

`tests/check_example_batch.py` 会反过来核对这张图里的数字，
所以算法一变、图没跟上，检查当场就会红。

踩过的两个坑（改这个脚本时别再犯）
----------------------------------
1. SVG 根元素必须带 `class="fig"`。页面里所有配色规则都写成 `.fig .xxx`，
   漏了这个 class 会得到一张"尺寸正确、内容完全空白"的图。
   图例也必须放在坐标系**外面**（这里输出成 HTML 放在图下方）——
   放图内右上角时，正好压住了 D2 那 3 个高年龄窗口。
2. 页面上的配色全部走 CSS 变量，深浅两套主题共用同一份 SVG；
   所以这里一个颜色都不写死，只给 class。
"""
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from druid.workflow import BatchConfig, run_batch   # noqa: E402

SPOT = "S01"
BATCH = ROOT / "examples" / "EX2022A"
HTML = ROOT / "docs" / "index.html"

W, H = 900, 372
ML, MR, MT, MB = 88, 30, 40, 66
PW, PH = W - ML - MR, H - MT - MB

BEGIN = "<!-- ==== 分域前后对比图：由 tools/gen_docs_split_figure.py 生成，请勿手改 ==== -->"
LEGACY = "<!-- 这张图由 build_tmp/_gen_split_svg.py"     # 迁移用，兼容旧标记


def nice_ticks(lo, hi, want=5):
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


def build() -> tuple[str, dict]:
    res = run_batch(BatchConfig(data_dir=str(BATCH), plot=False))
    win = res.windows[res.windows["Sample"] == SPOT].reset_index(drop=True)
    dom = res.domains[res.domains["样品"] == SPOT].reset_index(drop=True)
    ov = res.overall[res.overall["样品"] == SPOT].iloc[0]
    if not len(win) or not len(dom):
        raise SystemExit(f"示例批次里找不到 {SPOT} 的剖面窗口或年龄域，图没法画")

    tau = win["Tau"].to_numpy(float)
    age = win["Age68"].to_numpy(float)
    sig = win["Age68_1s"].to_numpy(float)
    step = float(np.median(np.diff(tau)))
    tol = step * 0.55

    # ── 真实数值：全部取自公开输出，一个都不手写 ──
    whole_age = float(ov["年龄_Ma"])
    whole_mswd = float(ov["MSWD"])
    whole_crit = float(ov["相容上限"])
    main_age = float(ov["主域年龄_Ma"])
    d_age = float(ov["Δ年龄_Ma"])
    d_pct = float(ov["Δ年龄_pct"])
    bulk_age = float(ov["整段积分年龄_Ma"])
    n_win = int(ov["n_win"])

    # ── 逐窗口归域。tau 是窗口起点，域表的 tau 写成「首个窗口起点-末个窗口起点」，
    #    所以按半个步长做容差（直接用区间过滤会漏掉最后一个窗口）。 ──
    belong = np.full(tau.size, "", dtype=object)
    bands = []
    for _, r in dom.iterrows():
        a, b = [float(x) for x in str(r["tau"]).split("-")]
        belong[(tau >= a - tol) & (tau <= b + tol)] = str(r["域"])
        bands.append((str(r["域"]), a, min(b + step, 1.0), float(r["年龄_Ma"]), int(r["n_win"])))
    n_off = int((belong == "").sum())
    assert int((belong == "D1").sum()) == int(dom["n_win"].iloc[0]), "D1 窗口数对不上"
    assert int((belong == "D2").sum()) == int(dom["n_win"].iloc[1]), "D2 窗口数对不上"
    assert n_off == n_win - int(dom["n_win"].sum()), "未归域窗口数对不上"

    lo, hi = float(np.nanmin(age - sig)), float(np.nanmax(age + sig))
    pad = (hi - lo) * 0.09
    y0, y1 = lo - pad, hi + pad
    yticks, x0, x1 = nice_ticks(y0, y1), 0.0, 1.0

    def X(t):
        return ML + (t - x0) / (x1 - x0) * PW

    def Y(a):
        return MT + (y1 - a) / (y1 - y0) * PH

    P = []
    P.append('<div class="figscroll">')
    P.append(f'<svg class="fig" viewBox="0 0 {W} {H}" role="img" '
             f'aria-labelledby="figSplitT figSplitD">')
    P.append(f'<title id="figSplitT">{SPOT} 测点：整段一个数，与分域后的阶梯</title>')
    P.append(f'<desc id="figSplitD">同一个剥蚀坑的 {n_win} 个滑窗年龄。灰虚线是整段加权平均出的 '
             f'{whole_age:.1f} Ma —— 它不对应任何一个窗口；两条实线是两个年龄域的均值。'
             f'中间 {n_off} 个窗口没通过域内相容性检验，被分域步骤剥掉了。</desc>')
    P.append(f'<clipPath id="cpSplit"><rect x="{ML}" y="{MT}" width="{PW}" height="{PH}"/></clipPath>')
    P.append(f'<rect x="0" y="0" width="{W}" height="{H}" rx="14" class="fig-bg"/>')

    for a in yticks:
        P.append(f'<line class="fig-grid" x1="{ML}" y1="{Y(a):.1f}" x2="{ML+PW}" y2="{Y(a):.1f}"/>')
        P.append(f'<text class="fig-tick" x="{ML-12}" y="{Y(a)+4:.1f}" '
                 f'text-anchor="end">{a:.0f}</text>')
    P.append(f'<text class="fig-axis" transform="translate(24,{MT+PH/2:.0f}) rotate(-90)" '
             f'text-anchor="middle">年龄 (Ma)</text>')

    for name, a, b, _, _ in bands:                       # 域底纹画在数据点下面
        P.append(f'<rect class="band-{name}" x="{X(a):.1f}" y="{MT}" '
                 f'width="{X(b)-X(a):.1f}" height="{PH}"/>')
    for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        P.append(f'<line class="fig-grid" x1="{X(t):.1f}" y1="{MT}" x2="{X(t):.1f}" y2="{MT+PH}"/>')
        P.append(f'<text class="fig-tick" x="{X(t):.1f}" y="{MT+PH+20}" '
                 f'text-anchor="middle">{t:g}</text>')

    P.append('<g clip-path="url(#cpSplit)">')
    for t, a, s, b in zip(tau, age, sig, belong):
        if not np.isfinite(a):
            continue
        x, y = X(t), Y(a)
        P.append(f'<line class="{"fig-ebar" if b else "fig-ebar-off"}" '
                 f'x1="{x:.1f}" y1="{Y(a-s):.1f}" x2="{x:.1f}" y2="{Y(a+s):.1f}"/>')
        P.append(f'<circle class="{"dot-"+b if b else "fig-dot-off"}" '
                 f'cx="{x:.1f}" cy="{y:.1f}" r="3.4"/>')
    P.append('</g>')

    P.append(f'<line class="fig-whole" x1="{ML}" y1="{Y(whole_age):.1f}" x2="{ML+PW}" '
             f'y2="{Y(whole_age):.1f}" stroke-dasharray="7 4"/>')
    for name, a, b, mage, _ in bands:
        ym = Y(mage)
        P.append(f'<line class="solid-{name}" x1="{X(a):.1f}" y1="{ym:.1f}" '
                 f'x2="{X(b):.1f}" y2="{ym:.1f}"/>')
        P.append(f'<text class="txt-{name}" x="{X((a+b)/2):.1f}" y="{ym-10:.1f}" '
                 f'text-anchor="middle">{name} {mage:.1f}</text>')
    if n_off:
        mid = X((bands[0][2] + bands[1][1]) / 2)
        P.append(f'<text class="fig-off-note" x="{mid:.1f}" y="{MT+18}" text-anchor="middle">'
                 f'{n_off} 个混合/过渡窗口（未进入任何年龄域）</text>')

    P.append(f'<rect class="fig-frame" x="{ML}" y="{MT}" width="{PW}" height="{PH}"/>')
    P.append(f'<text class="fig-axis" x="{ML+PW/2:.0f}" y="{H-14}" text-anchor="middle">'
             f'剥蚀时间 τ（0 = 激光开，1 = 激光关）　→ 坑深方向</text>')
    P.append('</svg>')
    P.append('</div>')
    P.append('<p class="figscroll-hint">图较宽，可左右滑动查看完整的坐标轴</p>')

    # 图例放图下方（HTML）：画进坐标系里一定会压住某个角上的数据
    P.append('<ul class="figlegend">')
    P.append(f'<li><i class="sw-whole"></i><span class="lg"><b>整段不分域 {whole_age:.1f} Ma</b>'
             f'<em>MSWD {whole_mswd:.2f} > 判据上限 {whole_crit:.2f} → 整段非常数</em></span></li>')
    P.append(f'<li><i class="sw-d1"></i><span class="lg"><b>D1 域均值 {main_age:.1f} Ma</b>'
             f'<em>{int(dom["n_win"].iloc[0])} 个窗口，MSWD {float(dom["MSWD"].iloc[0]):.2f}'
             f'</em></span></li>')
    P.append(f'<li><i class="sw-d2"></i><span class="lg"><b>D2 域均值 {bands[1][3]:.1f} Ma</b>'
             f'<em>{int(dom["n_win"].iloc[1])} 个窗口，MSWD {float(dom["MSWD"].iloc[1]):.2f}'
             f'</em></span></li>')
    P.append(f'<li><i class="sw-off"></i><span class="lg"><b>未进入任何年龄域：{n_off} 个窗口</b>'
             f'<em>域内相容性检验未通过，被分域步骤剥掉</em></span></li>')
    P.append(f'<li class="wide">逐点差额：整段 {whole_age:.1f} − 主域 {main_age:.1f} = '
             f'<b>{d_age:+.2f} Ma（{d_pct:+.2f}%）</b>，这就是「分域」这一步的贡献；'
             f'整段积分是另一种算法，本例给 {bulk_age:.1f} Ma。</li>')
    P.append('</ul>')

    nums = {"spot": SPOT, "n_win": n_win, "n_off": n_off,
            "whole_age": round(whole_age, 1), "main_age": round(main_age, 1),
            "d2_age": round(bands[1][3], 1)}
    return "\n".join(P), nums


def inject(frag: str) -> None:
    html = HTML.read_text(encoding="utf-8")
    i = html.find(BEGIN)
    if i == -1:
        i = html.find(LEGACY)
        if i == -1:
            raise SystemExit(f"{HTML} 里找不到图的标记（既不是新版也不是旧版）")
    tail = html[i:]
    k = tail.find('<ul class="figlegend">')
    if k == -1:
        j = i + tail.index("</svg>") + len("</svg>")          # 更早的版本没有图例
    else:
        j = i + tail.index("</ul>", k) + len("</ul>")
    assert html[i:j].count("<svg") == 1, "待替换区间里有多个 svg"
    body = BEGIN + "\n" + "\n".join("      " + ln for ln in frag.splitlines())
    HTML.write_text(html[:i] + body + html[j:], encoding="utf-8")


def main() -> int:
    frag, nums = build()
    inject(frag)
    print(f"示例批次 {BATCH.name} / {SPOT}：窗口 {nums['n_win']} 个，"
          f"其中 {nums['n_off']} 个未进入年龄域")
    print(f"  整段 {nums['whole_age']} Ma · D1 {nums['main_age']} Ma · D2 {nums['d2_age']} Ma")
    print(f"已写入 {HTML.relative_to(ROOT)}（{len(frag)} 字节的图 + 图例）")
    print("别忘了跑 `python tests/check_example_batch.py`：它会核对这张图里的数字。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
