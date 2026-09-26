# -*- coding: utf-8 -*-
"""
重新生成落地页 `docs/index.html` 首屏的「一眼看到什么」区块，并就地注入。

    python tools/gen_docs_preview.py

它出三样东西，全部来自示例批次 `examples/EX2022A` 的**真实输出**：

  1. 一张**逐深度年龄剖面**：某个真实测点的窗口点（含 1σ 误差棒）、各年龄域的
     均值线，以及整段不分域口径 —— 最后那个用右侧竖直箭头标出"不分域会差多少"。
     哪个测点见下面 `SPOT` 那段注释（判据是"能不能一眼看懂分域"，不是"落差最大"）。
  2. 三条核心质控判据的状态条（`standards.primary_bias` /
     `standards.secondary_bias` / `samples.concordance`），
     值、判据、结论**直接取自 `qc.assess_batch()`**，不另抄一份。
  3. 「标样QC」表的全部行（`res.qc`）—— 这张表页面上别处没有展示过。

为什么要有这个脚本
------------------
首屏是访客 5 秒内建立信任的地方，所以它必须放**真实输出**而不是示意图。
真实输出就必须能重新算出来 —— 否则它跟手抄的数字没区别，改了算法也没人知道它过期了。
`tests/check_example_batch.py` 会从注入区间里把这些数字抓回来跟实跑比，
抓不到或对不上都算失败。

和 `tools/gen_docs_split_figure.py` 的分工
------------------------------------------
那个脚本画**三种剖面形态**（均一 / 核边分明 / 复杂变化）放在 §02；
这个脚本画**首屏那一张**，用另一个测点，顺带给质控与标样QC 一个位置。
两者共用 `collect()`（窗口归域、两类域行的区分、窗口数守恒的断言都在那里）。

踩过的坑（改这个脚本时别再犯）
------------------------------
1. SVG 根元素必须带 `class="fig"`，否则配色规则全部落空 → 一张空白图。
2. `clipPath` 的 id 必须全页唯一（三联图用 `cpCase0..2`，这里用 `cpPreview`）。
3. y 轴范围要把**整段线与各域线**一起框进来，否则右侧那个 Δ 标注会跑到画框外面。
4. 域表的行分 `age domain` 与 `mixed/过渡带` 两类，只有前者能画成均值实线 ——
   过渡带的 `域` 列是破折号，拿它拼 CSS class 会得到一条默认黑色的线。
   这件事已经收进 `collect()` 了，这里只管用。
5. 质控层给的字符串里带 `<` `>`（判据）和 ASCII 负号（偏差），进 HTML 前要转义、
   负号要换成排版用的 −，否则页面上会出现裸标签和两种长的减号。
6. 描述测点分域情况的那句话曾经是个**写死的常量**（"两个年龄域，中间一段过渡带"）——
   换测点它不会跟着变，图注就开始撒谎。现在由 `spot_label()` 从 `collect()`
   的结果推出来。往图注里加任何描述性文字之前，先问一句"它是算出来的吗"。
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from druid.workflow import BatchConfig, run_batch        # noqa: E402
from druid.qc import assess_batch                        # noqa: E402
import gen_docs_split_figure as gds                      # noqa: E402

BATCH = ROOT / "examples" / "EX2022A"
HTML = ROOT / "docs" / "index.html"

#: 首屏那个测点。挑的是"一眼就能看懂分域"的那一类：
#:   ① 全部窗口都归进了年龄域 —— 图上没有一片空心点要解释；
#:   ② 每个域的 MSWD 都远在限内 —— 拆出来的每一段本身都是常数；
#:   ③ 整段 MSWD 又明显超限 —— 不分域就是错的。
#: 示例批次里同时满足这三条的有 S04 / S10 / S20；取其中域内 MSWD 最漂亮、
#: 三个域窗口数最均衡的 S10。不含 §02 三联图用掉的 S43 / S34 / S18。
#: 换测点：`python tools/gen_docs_preview.py --spot S20`（图注文字会跟着改）。
SPOT = "S10"

_NUM_CN = {1: "一", 2: "两", 3: "三", 4: "四", 5: "五"}


def spot_label(d: dict) -> str:
    """按**实际**分域结果生成一句测点描述。

    以前这里是一个写死的常量（"两个年龄域，中间一段过渡带"）—— 换个测点它不会
    跟着变，图注就会开始撒谎。现在从 `collect()` 的结果推出来。
    """
    n = len(d["bands"])
    s = f'{_NUM_CN.get(n, str(n))}个年龄域'
    if d["n_mixed"]:
        s += "，中间一段过渡带"
    n_bad = d["n_hollow"] - d["n_mixed"]
    if n_bad:
        s += f'，另有 {n_bad} 个窗口没通过相容性检验'
    return s


#: 首屏状态条上放哪三条检查项。key 是契约，改它要同步改回归守护。
MINI_KEYS = ("standards.primary_bias", "standards.secondary_bias", "samples.concordance")
MINI_LABELS = {
    "standards.primary_bias": "主标 91500 偏差",
    "standards.secondary_bias": "监控标样 Ple 偏差",
    "samples.concordance": "样品协和度占比",
}
LEVEL_TXT = {"pass": "通过", "warn": "注意", "fail": "不通过", "info": "提示"}

BEGIN = "<!-- ==== 首屏预览：由 tools/gen_docs_preview.py 生成，请勿手改 ==== -->"
END = "<!-- ==== 首屏预览结束 ==== -->"

# ── 版面 ──────────────────────────────────────────────────────────────────
W = 860
ML, MR = 76, 26
PW = W - ML - MR
PH = 196
HEAD, TAIL, TOP, BOTTOM = 30, 26, 14, 30
H = TOP + HEAD + PH + TAIL + BOTTOM


def esc(s) -> str:
    """把质控层的字符串放进 HTML：先转义 & < >，再把 ASCII 负号换成 −。"""
    s = " ".join(str(s).split())
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return s.replace("-", "\u2212")


def panel(d: dict) -> list[str]:
    """画首屏那张单面板剖面图，返回 SVG 片段。"""
    pt = TOP + HEAD
    pb = pt + PH

    lo = float(np.nanmin(d["age"] - d["sig"]))
    hi = float(np.nanmax(d["age"] + d["sig"]))
    pad = (hi - lo) * 0.08
    # 整段线与各域线都得落在视野里 —— 否则右侧的 Δ 标注会画到画框外
    ys = [lo - pad, hi + pad, d["whole"] - pad, d["whole"] + pad]
    for bd in d["bands"]:
        ys += [bd["age"] - pad, bd["age"] + pad]
    y0v, y1v = min(ys), max(ys)
    ticks = gds.nice_ticks(y0v, y1v, want=5)

    def X(t: float) -> float:
        return ML + (t - 0.0) / 1.0 * PW

    def Y(a: float) -> float:
        return pt + (y1v - a) / (y1v - y0v) * PH

    P: list[str] = [f'<g class="case" data-spot="{d["spot"]}">']
    P.append(f'<text class="fig-cap" x="{ML}" y="{pt-11:.0f}">'
             f'{d["spot"]} · {spot_label(d)}'
             f'<tspan class="fig-cap-sub">　{d["n_win"]} 个 4 秒滑窗</tspan></text>')
    P.append(f'<text class="fig-sub" x="{ML+PW}" y="{pt-11:.0f}" text-anchor="end">'
             f'整段 MSWD {d["whole_mswd"]:.2f} / 上限 {d["crit"]:.2f} → 整段非常数</text>')

    for a in ticks:
        P.append(f'<line class="fig-grid" x1="{ML}" y1="{Y(a):.1f}" x2="{ML+PW}" y2="{Y(a):.1f}"/>')
        P.append(f'<text class="fig-tick" x="{ML-12}" y="{Y(a)+4:.1f}" text-anchor="end">'
                 f'{a:.0f}</text>')
    P.append(f'<text class="fig-axis" transform="translate(22,{pt+PH/2:.0f}) rotate(-90)" '
             f'text-anchor="middle">年龄 (Ma)</text>')

    for bd in d["bands"]:
        P.append(f'<rect class="band-{bd["name"]}" x="{X(bd["a"]):.1f}" y="{pt}" '
                 f'width="{X(bd["b"])-X(bd["a"]):.1f}" height="{PH}"/>')

    for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        P.append(f'<line class="fig-grid" x1="{X(t):.1f}" y1="{pt}" x2="{X(t):.1f}" y2="{pb}"/>')
        P.append(f'<text class="fig-tick" x="{X(t):.1f}" y="{pb+17}" text-anchor="middle">'
                 f'{t:g}</text>')

    P.append(f'<clipPath id="cpPreview"><rect x="{ML}" y="{pt}" width="{PW}" height="{PH}"/></clipPath>')
    P.append('<g clip-path="url(#cpPreview)">')
    for t, a, s, b in zip(d["tau"], d["age"], d["sig"], d["belong"]):
        if not np.isfinite(a):
            continue
        x, y = X(t), Y(a)
        cls = f'dot-{b}' if b else "fig-dot-off"
        P.append(f'<line class="{"fig-ebar" if b else "fig-ebar-off"}" '
                 f'x1="{x:.1f}" y1="{Y(a-s):.1f}" x2="{x:.1f}" y2="{Y(a+s):.1f}"/>')
        P.append(f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="3.4"/>')
    P.append('</g>')

    # 灰虚线＝整段不分域口径；彩色实线＝各**年龄域**的均值
    P.append(f'<line class="fig-whole" x1="{ML}" y1="{Y(d["whole"]):.1f}" x2="{ML+PW}" '
             f'y2="{Y(d["whole"]):.1f}" stroke-dasharray="7 4"/>')
    P.append(f'<text class="fig-whole-t" x="{ML+4}" y="{Y(d["whole"])-7:.1f}">'
             f'整段不分域 {d["whole"]:.1f}</text>')
    for bd in d["bands"]:
        ym = Y(bd["age"])
        P.append(f'<line class="solid-{bd["name"]}" x1="{X(bd["a"]):.1f}" y1="{ym:.1f}" '
                 f'x2="{X(bd["b"]):.1f}" y2="{ym:.1f}"/>')
        P.append(f'<text class="txt-{bd["name"]}" x="{X((bd["a"]+bd["b"])/2):.1f}" y="{ym-9:.1f}" '
                 f'text-anchor="middle">{bd["name"]} {bd["age"]:.1f}</text>')

    # 右侧竖直箭头：把"分域这一步"的价值标成一个数字（整段 − 主域）
    main = next(bd for bd in d["bands"] if bd["name"] == d["main_dom"])
    xa = ML + PW - 7
    ya, yb = Y(d["whole"]), Y(main["age"])
    P.append(f'<line class="fig-delta" x1="{xa}" y1="{ya:.1f}" x2="{xa}" y2="{yb:.1f}"/>')
    for yy in (ya, yb):
        P.append(f'<line class="fig-delta" x1="{xa-6}" y1="{yy:.1f}" x2="{xa}" y2="{yy:.1f}"/>')
    delta_txt = esc(f"Δ {d['d_age']:+.1f} Ma")
    P.append(f'<text class="fig-delta-t" x="{xa-11:.1f}" y="{(ya+yb)/2+4:.1f}" '
             f'text-anchor="end">{delta_txt}</text>')

    if d["n_hollow"]:
        extra = (f'（{d["n_mixed"]} 个在过渡带，{d["n_hollow"]-d["n_mixed"]} 个没通过相容性检验）'
                 if d["n_mixed"] else '（没进入任何年龄域）')
        P.append(f'<text class="fig-off-note" x="{ML+PW/2:.1f}" y="{pt+17}" '
                 f'text-anchor="middle">{d["n_hollow"]} 个空心点{extra}</text>')
    P.append(f'<rect class="fig-frame" x="{ML}" y="{pt}" width="{PW}" height="{PH}"/>')
    P.append(f'<text class="fig-axis" x="{ML+PW/2:.0f}" y="{H-12}" text-anchor="middle">'
             f'剥蚀时间 τ（0 = 激光开，1 = 激光关）　→ 坑深方向</text>')
    P.append('</g>')
    return P


def profile_figure(d: dict) -> str:
    P = ['<div class="figscroll">']
    P.append(f'<svg class="fig" viewBox="0 0 {W} {H}" role="img" '
             f'aria-labelledby="figPrevT figPrevD" aria-describedby="figPrevD">')
    P.append('<title id="figPrevT">示例批次 EX2022A 的一个真实测点：'
             '逐深度年龄剖面与它的两个年龄域</title>')
    P.append('<desc id="figPrevD">横轴是剥蚀时间（对应坑深），纵轴是年龄。'
             '每个点是一个 4 秒滑窗，误差棒是 1σ。灰虚线是整段不分域的加权平均，'
             '两条彩色实线是两个年龄域的均值，空心点是没进入任何年龄域的窗口。'
             '右侧竖直箭头标出整段口径与主域之差。'
             '数据取自示例批次的公开结果表，可用仓库里的命令重新算出同样的数字。</desc>')
    P.append(f'<rect x="0" y="0" width="{W}" height="{H}" rx="14" class="fig-bg"/>')
    P.extend(panel(d))
    P.append('</svg>')
    P.append('</div>')
    # 窄屏下这张图会横向滚动（`.figscroll` 在 760px 以下给 .fig 一个最小宽度），
    # 得有人告诉用户"右边还有" —— 三联图那边也是这么处理的。
    P.append('<p class="figscroll-hint">图较宽，可左右滑动查看完整的坐标轴</p>')
    return "\n".join("        " + ln for ln in "\n".join(P).splitlines())


def mini_cards(checks: dict) -> str:
    """三条核心判据的状态条。值 / 判据 / 结论全部取自质控层本身。"""
    P = ['<div class="mini">']
    for i, key in enumerate(MINI_KEYS):
        c = checks[key]
        P.append(f'  <div class="mcard {c.level}" data-key="{key}">')
        P.append(f'    <div class="mtop"><span class="mn">{"①②③"[i]}</span>'
                 f'<span class="mname">{MINI_LABELS[key]}</span>'
                 f'<span class="pill {c.level}">{LEVEL_TXT.get(c.level, c.level)}</span></div>')
        P.append(f'    <div class="mval">{esc(c.observed)}</div>')
        P.append(f'    <div class="mcri"><code class="mkey">{key}</code> · '
                 f'{esc(c.criterion)}</div>')
        P.append('  </div>')
    P.append('</div>')
    return "\n    ".join(P)


def qc_table(qc) -> str:
    """「标样QC」表 —— 示例批次只有两行（主标 / 监控标样），照实排。"""
    P = ['<table class="qctab">']
    P.append('<caption>标样QC 表（结果.xlsx 的第 2 张表）</caption>')
    P.append('<thead><tr><th>标样</th><th>点数</th><th>参考年龄</th><th>加权平均</th>'
             '<th>2σ</th><th>MSWD</th><th>偏差</th></tr></thead>')
    P.append('<tbody>')
    for _, r in qc.iterrows():
        bias = float(r["偏差_pct"])
        cls = "bias-neg" if bias < 0 else "bias-pos"
        # 负号只替换**这一格**的显示文本 —— 整行一起换会把 class="bias-neg" 也改掉
        bt = f'{bias:+.2f}%'.replace("-", "\u2212")
        P.append(f'  <tr><td>{esc(r["标样"])}</td><td>{int(r["点数"])}</td>'
                 f'<td>{float(r["参考年龄_Ma"]):.2f}</td>'
                 f'<td><b>{float(r["加权平均年龄_Ma"]):.2f}</b></td>'
                 f'<td>{float(r["s2_Ma"]):.2f}</td>'
                 f'<td>{float(r["MSWD"]):.2f}</td>'
                 f'<td class="{cls}">{bt}</td></tr>')
    P.append('</tbody></table>')
    P.append('<p class="qcfoot">年龄单位 Ma。偏差 = 加权平均 − 参考值'
             '（参考值取 Horstwood et al. 2016 表 S2；参考值口径 / 归一化 / 二次校正'
             '三件事都写在「运行参数」表里）。</p>')
    return "\n      ".join(P)


def legend_line(d: dict) -> str:
    d_age, d_pct = d["d_age"], d["d_pct"]
    segs = " · ".join(f'{bd["name"]} {bd["age"]:.1f} Ma（{bd["n"]} 窗，MSWD {bd["mswd"]:.2f}）'
                      for bd in d["bands"])
    parts = []
    if d["n_mixed"]:
        parts.append(f'中间的过渡带另有 {d["n_mixed"]} 个窗口，它不是一个年龄，'
                     f'所以不画均值线')
    if d["n_off"] - d["n_mixed"]:
        parts.append(f'{d["n_off"]-d["n_mixed"]} 个窗口没通过域内相容性检验')
    tail = "；".join(parts) + "。" if parts else "没有一个窗口被剥掉。"
    return (f'<b>{d["spot"]} {spot_label(d)}</b>：共 {d["n_win"]} 个滑窗，'
            f'{len(d["bands"])} 个域 —— {segs}。{tail}整段不分域给 {d["whole"]:.1f} Ma，'
            f'离主域 {d["main_age"]:.1f} Ma 差 <b>{esc(f"{d_age:+.1f} Ma"
                                                  f"（{d_pct:+.2f}%）")}</b>'
            f'（Δ = 整段 − 主域）—— 这就是"分域"这一步买到的东西。')


def build(spot: str) -> tuple[str, dict]:
    cfg = BatchConfig(data_dir=str(BATCH), plot=False)
    res = run_batch(cfg)
    d = gds.collect(res, spot)
    checks = {c.key: c for c in assess_batch(res, cfg)}
    missing = [k for k in MINI_KEYS if k not in checks]
    if missing:
        raise SystemExit(f"质控层里没有这些检查项：{missing}（key 改过？）")
    if len(d["bands"]) < 2:
        raise SystemExit(f"{spot} 只有 {len(d['bands'])} 个年龄域 —— 首屏这张图需要两个")
    n_spot = len(res.overall)
    n_multi = int((res.overall["域数"] > 1).sum())

    P = ['<div class="pv">']
    P.append('  <figure class="pvfig">')
    P.append('    ' + profile_figure(d))
    P.append('    <figcaption>这就是<b>一个测点</b>给出的定量结论'
             f'（整批 {n_spot} 个样品测点各出一张，其中 {n_multi} 个检出多年龄域）。'
             '上面这张图是示例批次 <code>examples/EX2022A</code> 的<b>真实输出</b>，'
             '不是示意图 —— 仓库里的命令跑一遍就能得到同样的数字。'
             f'{legend_line(d)}</figcaption>')
    P.append('  </figure>')
    P.append('  <div class="pvside">')
    P.append('    <p class="pvhead">先看这三个数　<span>判据原文，不是我们另拍的</span></p>')
    P.append('    ' + mini_cards(checks).replace("\n", "\n    "))
    P.append('    <p class="pvmore">三条都过了才往下看年龄。'
             '完整的 26 条检查项在<a href="#gate">第 05 节</a>，'
             '也可以用结果旁边的 <code>.handoff.json</code> 给程序读。</p>')
    P.append('    ' + qc_table(res.qc))
    P.append('  </div>')
    P.append('</div>')

    info = {
        "spot": d["spot"], "n_win": d["n_win"], "n_off": d["n_off"],
        "n_mixed": d["n_mixed"], "n_hollow": d["n_hollow"],
        "whole": d["whole"], "whole_mswd": d["whole_mswd"], "crit": d["crit"],
        "d_age": d["d_age"], "d_pct": d["d_pct"], "main_dom": d["main_dom"],
        "main_age": d["main_age"],
        "bands": [(bd["name"], bd["age"], bd["n"], bd["mswd"]) for bd in d["bands"]],
        "qc": [(str(r["标样"]), int(r["点数"]), float(r["参考年龄_Ma"]),
                float(r["加权平均年龄_Ma"]), float(r["s2_Ma"]), float(r["MSWD"]),
                float(r["偏差_pct"])) for _, r in res.qc.iterrows()],
        "mini": [(k, checks[k].level, checks[k].observed, checks[k].criterion)
                 for k in MINI_KEYS],
    }
    return "\n".join(P), info


def inject(frag: str, html_path: pathlib.Path) -> None:
    """把生成的片段就地替换进页面（含首尾标记之间的一切）。"""
    html = html_path.read_text(encoding="utf-8")
    i = html.find(BEGIN)
    j = html.find(END)
    if i == -1 or j == -1 or j < i:
        raise SystemExit(f"{html_path} 里找不到首屏预览的标记：{BEGIN} … {END}")
    body = BEGIN + "\n" + "\n".join("      " + ln for ln in frag.splitlines()) + "\n" + END
    html_path.write_text(html[:i] + body + html[j + len(END):], encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--html", default=str(HTML), help="注入目标（草稿时可以指向临时副本）")
    ap.add_argument("--spot", default=SPOT,
                    help=f"首屏那个测点（默认 {SPOT}）。换它之前先看本模块顶部那段挑选标准")
    ap.add_argument("--stdout", action="store_true", help="只把片段打到标准输出，不写文件")
    args = ap.parse_args()

    frag, info = build(args.spot)
    if args.stdout:
        print(frag)
        return 0
    inject(frag, pathlib.Path(args.html))

    print(f"示例批次 {BATCH.name}，首屏测点 {info['spot']}：")
    print(f"  窗口 {info['n_win']} 个：年龄域内 {info['n_win']-info['n_hollow']}，"
          f"过渡带 {info['n_mixed']}，没有任何域 {info['n_off']}；"
          f"整段 {info['whole']:.4f} Ma（MSWD {info['whole_mswd']:.2f}/{info['crit']:.2f}）")
    for name, age, n, mswd in info["bands"]:
        print(f"    {name} {age:.4f} Ma（n={n}，MSWD {mswd:.2f}）")
    print(f"  整段 − 主域 {info['main_dom']} = {info['d_age']:+.4f} Ma（{info['d_pct']:+.2f}%）")
    for name, n, ref, avg, s2, mswd, bias in info["qc"]:
        print(f"  {name}: n={n} 参考 {ref:.2f} 实测 {avg:.4f} 2σ {s2:.2f} "
              f"MSWD {mswd:.2f} 偏差 {bias:+.4f}%")
    for key, lv, obs, cri in info["mini"]:
        print(f"  [{lv:>4s}] {key}: {obs}  ({cri})")
    print(f"已写入 {args.html}（{len(frag)} 字节）")
    print("别忘了跑 `python tests/check_example_batch.py`：它会从这段里把数字抓回去核对。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
