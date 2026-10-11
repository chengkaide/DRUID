#!/usr/bin/env python
"""
tools/gen_docs_segmented.py —— 生成《一段一段地读》那页的全部数字、表格与图（只读）
================================================================================

它做什么
--------
把模板 `tools/segmented_reading.page.html` 里所有 `{{占位符}}` 与三张图，全部
**现算**后注入，写出 `docs/segmented-reading.html`：

    · 扫描示例批次 `examples/EX2022A` 的每个样品测点（口径复用
      `tools/spot_concordia.py` 的 `prepare()`，一个数都不另算）；
    · 汇总五把"域间分辨力"尺子、域级协和度偏离、协和度沿深度的单调性，
      以及"把域当独立观测"会伪造出多少精度；
    · 画三张图（matplotlib → base64，页面因此仍是自包含单文件）；
    · 把结果写进页面，并断言页面上不再有未替换的 `{{...}}`。

⚠ **模板与产物是两个文件，这是有意的。** 第一版把页面本身当模板（就地替换占位符），
   于是生成过一次之后占位符全没了 —— 再跑一遍是空操作，`--check` 永远说"一致"。
   那种检查比没有更糟：它会替一个已经过期的页面背书。现在模板独立存放，
   `--check` 是"重算一遍 ≠ 已提交的页面"，真的会红。

⚠ **域级 1σ 一律乘 `synth_sigma_audit.SIGMA_INFLATION`**（合成实测的乐观倍数，
   表里另留"用自报 σ"的对照列）。这是本页与工具默认输出的唯一差别：
   比值、年龄、判定一概不动。

为什么写成生成器而不是手写页面
------------------------------
仓库守则：**文档里的数字要么是脚本算出来的，要么别写**。手抄的数字会在下一次
改算法时静默过期 —— 页面不报错，只会开始骗人。所以这一页的每一个小数都从这里
出来；`--check` 可以随时重算一遍并与已提交的页面比对（不一致就非零退出）。

用法
----
    python tools/gen_docs_segmented.py            # 重算并按模板写出页面
    python tools/gen_docs_segmented.py --check    # 只比对，不写文件

⚠ 本脚本**只读**：不写结果表、不改任何数值输出、不碰 `结果/`。
"""
from __future__ import annotations

import argparse
import base64
import importlib.util
import io
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")                                    # 批处理环境无显示器
import matplotlib.pyplot as plt                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "spot_concordia", Path(__file__).resolve().parent / "spot_concordia.py")
sc = importlib.util.module_from_spec(_spec)
sys.modules["spot_concordia"] = sc
_spec.loader.exec_module(sc)

_aspec = importlib.util.spec_from_file_location(
    "synth_sigma_audit", Path(__file__).resolve().parent / "synth_sigma_audit.py")
AUDIT = importlib.util.module_from_spec(_aspec)
sys.modules["synth_sigma_audit"] = AUDIT
_aspec.loader.exec_module(AUDIT)

# 域级 1σ 的校正因子 —— **唯一出处是 tools/synth_sigma_audit.py**
# （那里有实测：条件⑤ 全套效应下 pull 标准差 1.81、median|pull|/0.6745 1.90 …）
K = float(AUDIT.SIGMA_INFLATION)

PAGE = ROOT / "docs" / "segmented-reading.html"
TEMPLATE = ROOT / "tools" / "segmented_reading.page.html"
BATCH = ROOT / "examples" / "EX2022A"

# 图里的字必须与页面同一套配色（浅色主题下看着才不突兀）
INK = "#2C2C2A"
INK2 = "#5F5E5A"
GRID = "#D3D1C7"
C_AGE = "#185FA5"
C_PLANE = "#A32D2D"
C_CONC = "#7A3FA5"
C_CHEM = "#0F6E56"
C_WARN = "#C0392B"


# ─────────────────────────────────────────────────────────────────────────────
# 一、扫描（批次级步骤只算一次）
# ─────────────────────────────────────────────────────────────────────────────
def build_session():
    """跑一次批次级前置步骤，并把 `spot_concordia` 的模块级函数顶掉成常量。"""
    cfg = sc.BatchConfig(data_dir=str(BATCH), plot=False, verbose=False)
    spots, seq, skipped = sc.load_batch(cfg)
    brack = sc.build_bracketing(spots, cfg)
    ref68, ref76 = sc.std_ref(sc.std_alias(cfg.primary) or cfg.primary,
                              cfg.ref_preset)
    variants, (b68, b76) = sc.compute_bulk_ratios(spots, brack, ref68, cfg)
    pidx, F68, F76, cal = sc.calibrate_primary(spots, b68, b76, ref68, ref76, cfg)
    sd68, sd76 = sc.estimate_external(spots, pidx, F68, F76, b68, b76, cfg)
    tbl = sc.build_results(spots, pidx, F68, F76, b68, b76, sd68, sd76, cfg)

    sc.load_batch = lambda c: (spots, seq, skipped)
    sc.build_bracketing = lambda s, c: brack
    sc.std_ref = lambda *a, **k: (ref68, ref76)
    sc.compute_bulk_ratios = lambda *a, **k: (variants, (b68, b76))
    sc.calibrate_primary = lambda *a, **k: (pidx, F68, F76, cal)
    sc.estimate_external = lambda *a, **k: (sd68, sd76)
    sc.build_results = lambda *a, **k: tbl
    return dict(cfg=cfg, sd68=float(sd68), sd76=float(sd76),
                names=sorted({t.sample for t, _ in spots
                              if t.role == sc.ROLE_UNKNOWN}))


def _pair_stats(P, a, b, sd68, sd76, k=K):
    """相邻两域的分辨力：五把尺子各给一个等效 σ 数（与工具自己的打印一致）。

    `k` 把域级 1σ 放大成"合成实测的诚实值"（`synth_sigma_audit.SIGMA_INFLATION`）：
    三把与 σ 有关的尺子（z_age / z_conc / z_plane）按 1/k 缩小；
    两把化学尺子（z_ThU / z_U）用的是**窗口内部散布**，严格说需要它自己的合成实验，
    这里一并按同一个 k 缩放，只为让五把尺子可比 —— 见页面里的说明。
    """
    thu = P["prof"]["ThU"].to_numpy(float)
    ucp = P["prof"]["U_cps"].to_numpy(float)
    d_age = abs(a["age68"] - b["age68"])
    z_age = d_age / np.hypot(a["age68_2s"] / 2.0, b["age68_2s"] / 2.0) / k
    z_conc = (abs(a["concordance"] - b["concordance"]) / 100.0
              / np.hypot(sc._sigma_ln_concordance(a, sd68, sd76),
                         sc._sigma_ln_concordance(b, sd68, sd76)) / k)
    # 马氏距离：σ 放大 k ⇒ 距离缩小 k（在协方差里放大，而不是事后除）
    cov = np.array([[
        (k * a["s68"]) ** 2 + (k * b["s68"]) ** 2,
        0.5 * (a["rho"] * k * a["s68"] * k * a["s75"]
               + b["rho"] * k * b["s68"] * k * b["s75"])],
        [0.5 * (a["rho"] * k * a["s68"] * k * a["s75"]
                + b["rho"] * k * b["s68"] * k * b["s75"]),
         (k * a["s75"]) ** 2 + (k * b["s75"]) ** 2]])
    dv = np.array([b["r68"] - a["r68"], b["r75"] - a["r75"]])
    try:
        z_plane = sc._z_equiv(float(np.sqrt(dv @ np.linalg.inv(cov) @ dv)))
    except np.linalg.LinAlgError:
        z_plane = float("nan")
    z_age0 = z_age * k
    z_plane0 = sc._z_equiv(float(np.sqrt(dv @ np.linalg.inv(cov / k ** 2) @ dv))) \
        if np.isfinite(z_plane) else float("nan")
    return dict(z_age=float(z_age), z_conc=float(z_conc), z_plane=float(z_plane),
                z_age0=float(z_age0), z_conc0=float(z_conc * k),
                z_plane0=float(z_plane0),
                z_ThU=float(sc._chem_z(thu, a, b) / k),
                z_U=float(sc._chem_z(ucp, a, b) / k),
                z_ThU0=float(sc._chem_z(thu, a, b)),
                z_U0=float(sc._chem_z(ucp, a, b)),
                d_age=float(d_age))


def scan(sess):
    """逐样品取回：域、域对分辨力、域级协和度、逐窗口 207/235 的相对 1σ。"""
    cfg, sd68, sd76 = sess["cfg"], sess["sd68"], sess["sd76"]
    recs = []
    for nm in sess["names"]:
        try:
            P = sc.prepare(cfg, nm)
        except SystemExit:
            continue
        doms = P["doms"]
        with np.errstate(divide="ignore", invalid="ignore"):
            rel75 = P["s75"] / P["r75"]
        rec = dict(sample=nm, tag=str(P["tag"]), n_win=int(len(P["prof"])),
                   win_rel75=float(np.nanmedian(rel75)),
                   ws_age=float(P["ws"].get("age_Ma", np.nan)),
                   ws_mswd=float(P["ws"].get("mswd", np.nan)),
                   ws_ok=bool(P["ws"].get("mswd_ok", False)),
                   bulk=({k: float(v) for k, v in P["bulk"].items()}
                         if P["bulk"] else None),
                   doms=[dict(name=d["name"], n_win=d["n_win"], age68=d["age68"],
                              age68_2s=d["age68_2s"], mswd=d["mswd"],
                              age75=d["age75"], concordance=d["concordance"],
                              conc_sigma=d["concordance"]
                              * sc._sigma_ln_concordance(d, sd68, sd76),
                              th_u=d["th_u"]) for d in doms],
                   pairs=[])
        rec["pairs"] = [_pair_stats(P, a, b, sd68, sd76)
                        for a, b in zip(doms, doms[1:])]
        recs.append(rec)
    return recs


# ─────────────────────────────────────────────────────────────────────────────
# 二、汇总：占位符的值 + 四张表
# ─────────────────────────────────────────────────────────────────────────────
def _pct(x):
    return f"{x * 100:.0f}"


def metrics(recs):
    V = {}
    doms = [(r["sample"], d) for r in recs for d in r["doms"]]
    pairs = [p for r in recs for p in r["pairs"]]
    multi = [r for r in recs if len(r["doms"]) >= 2]
    multi3 = [r for r in recs if len(r["doms"]) >= 3]

    V["N_SPOTS"] = str(len(recs))
    V["N_DOM"] = str(len(doms))
    V["N_MULTI"] = str(len(multi))
    V["N_UNI"] = str(len(recs) - len(multi))
    V["MULTI_PCT"] = _pct(len(multi) / len(recs))
    V["N_PAIR"] = str(len(pairs))

    # ── 五把尺子（z 已按合成实测的乐观倍数 k 校正；同时留下列未校正的对照）──
    rulers = [("年龄轴 <code>z_age</code>", "z_age", C_AGE),
              ("协和平面 <code>z_plane</code>", "z_plane", C_PLANE),
              ("协和度轴 <code>z_conc</code>", "z_conc", C_CONC),
              ("化学 <code>z_ThU</code>", "z_ThU", C_CHEM),
              ("化学 <code>z_U</code>（238U 计数率）", "z_U", C_CHEM)]
    rows = []
    stat = {}
    for label, key, _c in rulers:
        v = np.array([p[key] for p in pairs], float)
        v = v[np.isfinite(v)]
        stat[key] = v
        rows.append((label, key))
    for key in ("z_age", "z_plane", "z_conc", "z_ThU", "z_U"):
        v0 = np.array([p[key + "0"] for p in pairs], float)
        stat[key + "0"] = v0[np.isfinite(v0)]
    V["K"] = f"{K:g}"
    V["Z_AGE_MED"] = f"{np.median(stat['z_age']):.2f}"
    V["Z_PLANE_MED"] = f"{np.median(stat['z_plane']):.2f}"
    V["Z_CONC_MED"] = f"{np.median(stat['z_conc']):.2f}"
    V["Z_U_MED"] = f"{np.median(stat['z_U']):.2f}"
    V["Z_AGE_MED0"] = f"{np.median(stat['z_age0']):.2f}"
    V["Z_CONC_MED0"] = f"{np.median(stat['z_conc0']):.2f}"
    V["Z_CONC_HIT"] = _pct((stat["z_conc"] >= 3).mean())
    V["Z_CONC_HIT0"] = _pct((stat["z_conc0"] >= 3).mean())
    V["Z_U_HIT"] = _pct((stat["z_U"] >= 3).mean())
    V["Z_AGE_HIT"] = _pct((stat["z_age"] >= 3).mean())
    ok = (np.isfinite(stat["z_age"]) & np.isfinite(stat["z_plane"])
          & np.isfinite(stat["z_conc"]))
    V["PLANE_OVER_AGE"] = f"{np.median(stat['z_plane'][ok] / stat['z_age'][ok]):.2f}"
    V["PLANE_WINS"] = str(int(((stat["z_age"] < 3) & (stat["z_plane"] >= 3) & ok).sum()))
    V["AGE_WINS"] = str(int(((stat["z_age"] >= 3) & (stat["z_plane"] < 3) & ok).sum()))

    #   左两列 = 本页展示值（三把 σ 尺子按 1/K 校正；化学两把也按 1/K 缩放，
    #   见 `_pair_stats` 的说明）；右两列 = 工具自报口径，留作对照。
    t1 = ['<div class="tscroll"><table>',
          "<tr><th>尺子</th><th>中位 z<br><span class='dim'>σ 校正后</span></th>"
          "<th>≥ 3 的比例<br><span class='dim'>σ 校正后</span></th>"
          "<th>中位 z<br><span class='dim'>用自报 σ</span></th>"
          "<th>≥ 3 的比例<br><span class='dim'>用自报 σ</span></th></tr>"]
    for label, key in rows:
        v, v0 = stat[key], stat[key + "0"]
        med, hit = float(np.median(v)), float((v >= 3).mean())
        m0, h0 = float(np.median(v0)), float((v0 >= 3).mean())
        cls = "n hi" if "协和度轴" in label else "n"
        t1.append(f"<tr><td>{label}</td><td class='{cls}'>{med:.2f}</td>"
                  f"<td class='n'>{hit*100:.0f}%</td>"
                  f"<td class='n dim'>{m0:.2f}</td><td class='n dim'>{h0*100:.0f}%</td></tr>")
    t1.append("</table></div>")
    V["T1"] = "\n".join(t1)

    # ── 域级协和度偏离 ──
    #   `z` = 用工具自报的 σ（未校正，留作对照）；`z / K` = 本页展示的校正后值。
    dev = []
    for s, d in doms:
        if not np.isfinite(d["conc_sigma"]) or d["conc_sigma"] <= 0:
            continue
        z = (d["concordance"] - 100.0) / d["conc_sigma"]
        dev.append((abs(z), s, d, z, z / K))
    zdev = np.array([x[3] for x in dev], float)
    zdev_k = np.array([x[4] for x in dev], float)
    V["DEV_GE2"] = str(int((np.abs(zdev) >= 2).sum()))
    V["DEV_GE3"] = str(int((np.abs(zdev) >= 3).sum()))
    V["DEV_MED"] = f"{np.median(np.abs(zdev)):.2f}"
    V["DEV_GE2_K"] = str(int((np.abs(zdev_k) >= 2).sum()))
    V["DEV_GE3_K"] = str(int((np.abs(zdev_k) >= 3).sum()))
    V["DEV_MED_K"] = f"{np.median(np.abs(zdev_k)):.2f}"
    V["CONC_POS"] = str(int((zdev > 0).sum()))
    V["CONC_NEG"] = str(int((zdev < 0).sum()))
    # 极端个例（按校正后的 |z| 排序取前三）—— 页面正文与图 1(b) 的标注用同一份名单
    #   （图 1(b) 用的是同一个"除以 K"的口径，两处独立算、结论必然一致）
    V["EXTREME"] = "、".join(
        f"{s} {d['name']}（{d['concordance']:.1f}%，z={zk:+.1f}）"
        for _a, s, d, _z, zk in sorted(dev, key=lambda x: -abs(x[4]))[:3])
    t2 = ['<div class="tscroll"><table>',
          "<tr><th>测点</th><th>域</th><th>窗数</th><th>206Pb/238U 年龄</th>"
          "<th>域级协和度</th><th>其 1σ<br><span class='dim'>σ 校正后</span></th>"
          "<th>z<br><span class='dim'>σ 校正后</span></th>"
          "<th>z<br><span class='dim'>用自报 σ</span></th>"
          "<th>域内 MSWD</th></tr>"]
    for _a, s, d, z, zk in sorted(dev, key=lambda x: -abs(x[4]))[:10]:
        cls = " hi" if abs(zk) >= 3 else ""
        t2.append(f"<tr><td class='n'>{s}</td><td class='n'>{d['name']}</td>"
                  f"<td class='n'>{d['n_win']}</td>"
                  f"<td class='n'>{d['age68']:.1f} Ma</td>"
                  f"<td class='n'>{d['concordance']:.1f}%</td>"
                  f"<td class='n'>±{d['conc_sigma'] * K:.1f}%</td>"
                  f"<td class='n{cls}'>{zk:+.1f}</td>"
                  f"<td class='n dim'>{z:+.1f}</td>"
                  f"<td class='n'>{d['mswd']:.2f}</td></tr>")
    t2.append("</table></div>")
    V["T2"] = "\n".join(t2)

    # ── 协和度沿深度的单调性 ──
    mono, nonmono = [], []
    for r in multi3:
        c = [d["concordance"] for d in r["doms"]]
        inc = all(b > a for a, b in zip(c, c[1:]))
        dec = all(b < a for a, b in zip(c, c[1:]))
        (mono if (inc or dec) else nonmono).append(
            (r["sample"], " → ".join(f"{x:.1f}%" for x in c), len(c)))
    V["MULTI3_TOTAL"] = str(len(multi3))
    V["MONO_N"] = str(len(mono))
    V["NONMONO_N"] = str(len(nonmono))
    t4 = ['<div class="two">',
          '<div class="card t"><h4>单调（%d 个）—— 与单条混合轨迹相容</h4>'
          "<ul>" % len(mono)]
    for s, seq, k in mono:
        t4.append(f"<li><code>{s}</code>（{k} 域）{seq}</li>")
    t4.append("</ul></div>")
    t4.append('<div class="card p"><h4>非单调（%d 个）—— 单条轨迹造不出</h4>'
              "<ul>" % len(nonmono))
    for s, seq, k in nonmono:
        t4.append(f"<li><code>{s}</code>（{k} 域）{seq}</li>")
    t4.append("</ul></div></div>")
    V["T4"] = "\n".join(t4)

    # ── 整段与常数模型 ──
    mswd = np.array([r["ws_mswd"] for r in recs], float)
    mswd = mswd[np.isfinite(mswd)]
    V["WS_MSWD_MED"] = f"{np.median(mswd):.2f}"
    V["WS_OK"] = str(sum(1 for r in recs if r["ws_ok"]))
    bc = np.array([r["bulk"]["concordance"] for r in recs if r["bulk"]], float)
    V["BULK_MED"] = f"{np.median(bc):.2f}"
    V["BULK_OFF"] = str(int(((bc < 90) | (bc > 110)).sum()))

    # ── 把域当独立观测：伪精度 ──
    #   ★ 表里展示的 1σ 一律乘 K（与全页一致），但"看起来缩小多少倍"与 χ²/dof
    #     都是比值，两边同时放大 ⇒ 不受影响。
    pool = []
    for r in recs:
        d = r["doms"]
        if len(d) < 2:
            continue
        a = np.array([x["age68"] for x in d], float)
        s1 = np.array([x["age68_2s"] / 2.0 for x in d], float)
        w = 1.0 / s1 ** 2
        mu = float(np.sum(w * a) / np.sum(w))
        se = float(np.sqrt(1.0 / np.sum(w)))
        chi2 = float(np.sum((a - mu) ** 2 * w))
        pool.append((r["sample"], len(a), float(np.median(s1)), se,
                     s1.min(), chi2 / (len(a) - 1), float(a.max() - a.min())))
    med_se = np.array([x[2] for x in pool], float)
    pool_se = np.array([x[3] for x in pool], float)
    V["DOM_SE"] = f"{np.median(med_se) * K:.2f}"
    V["POOL_SE"] = f"{np.median(pool_se) * K:.2f}"
    V["POOL_GAIN"] = f"{np.median(med_se / pool_se):.1f}"
    V["POOL_CHI2_MED"] = f"{np.median([x[5] for x in pool]):.1f}"
    t3 = ['<div class="tscroll"><table>',
          "<tr><th>测点</th><th>域数</th>"
          "<th>单域 1σ 中位<br><span class='dim'>σ 校正后</span></th>"
          "<th>池化后 1σ<br><span class='dim'>σ 校正后</span></th>"
          "<th>看起来缩小</th><th>χ²/dof</th><th>域间跨度</th></tr>"]
    for s, k, s1, se, smin, chi2, span in sorted(pool, key=lambda x: -x[5])[:8]:
        t3.append(f"<tr><td class='n'>{s}</td><td class='n'>{k}</td>"
                  f"<td class='n'>{s1 * K:.2f} Ma</td><td class='n'>{se * K:.2f} Ma</td>"
                  f"<td class='n'>{s1/se:.1f}×</td>"
                  f"<td class='n hi'>{chi2:.0f}</td>"
                  f"<td class='n'>{span:.0f} Ma</td></tr>")
    t3.append("</table></div>")
    V["T3"] = "\n".join(t3)

    # ── 207/235 的误差（**同一口径比**：都是相对 1σ）──
    wr = np.array([r["win_rel75"] for r in recs], float)
    wr = wr[np.isfinite(wr)]
    V["S75_WIN_PCT"] = f"{np.median(wr) * 100:.1f}"
    both = [r for r in recs if r["bulk"]]
    lg = np.array([r["bulk"]["s75_legacy_ma"] for r in both], float)
    st = np.array([r["bulk"]["s75_strict_ma"] for r in both], float)
    a75 = np.array([r["bulk"]["age75"] for r in both], float)
    good = np.isfinite(a75) & (a75 > 0)
    V["S75_LEGACY"] = f"{np.median(lg):.2f}"
    V["S75_STRICT"] = f"{np.median(st):.2f}"
    V["S75_RATIO"] = f"{np.median(st) / np.median(lg):.2f}"
    V["S75_SPOT_REL"] = f"{np.median(st[good] / a75[good]) * 100:.1f}"
    V["S75_SPOT_REL_LEGACY"] = f"{np.median(lg[good] / a75[good]) * 100:.1f}"
    V["S75_WIN_OVER_SPOT"] = (f"{np.median(wr) / np.median(st[good] / a75[good]):.1f}"
                              if good.any() else "--")

    # ── σ 审计（唯一出处：tools/synth_sigma_audit.py）──
    M = AUDIT.MEASURED
    V["K_DATE"] = str(M["date"])
    V["K_SEEDS"] = str(M["seeds"])
    V["K_FULL_N"] = str(M["full"]["n"])
    V["K_FULL_PULL_STD"] = f"{M['full']['pull_std']:.2f}"
    V["K_FULL_COVER1"] = f"{M['full']['cover1'] * 100:.1f}"
    V["K_FULL_MED_ABS"] = f"{M['full']['med_abs']:.2f}"
    V["K_FULL_SIG"] = f"{M['full']['sig_med']:.2f}"
    V["K_FULL_ERR"] = f"{M['full']['err_med']:.2f}"
    V["K_COUNT_N"] = str(M["counting"]["n"])
    V["K_COUNT_PULL_STD"] = f"{M['counting']['pull_std']:.2f}"
    V["K_COUNT_COVER1"] = f"{M['counting']['cover1'] * 100:.1f}"
    V["K_COUNT_MED_ABS"] = f"{M['counting']['med_abs']:.2f}"
    return V


# ─────────────────────────────────────────────────────────────────────────────
# 三、三张图
# ─────────────────────────────────────────────────────────────────────────────
def _png_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=118, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _style():
    sc._style()
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["figure.facecolor"] = "white"


def _fig_frame(b64, caption, alt):
    return (f'<figure><img src="data:image/png;base64,{b64}" alt="{alt}">'
            f"<figcaption>{caption}</figcaption></figure>")


def _inflate(P, k):
    """
    把图上出现的域级 1σ 一律乘 k 倍 —— 这是本页与工具默认输出的唯一差别。

    `tools/synth_sigma_audit.py` 实测：条件⑤（全套效应）下报出的 1σ 偏乐观约 k 倍
    （pull 标准差 1.81、median|pull|/0.6745 1.90）。图与表都按 k 校正后展示，
    同时把未校正的对照值留在表里 —— 免得读者以为工具本来就报这么宽的误差。
    比值、年龄、判定一概不动，只动 σ。
    """
    for key in ("s68", "s75", "s76"):
        if P.get(key) is not None:
            P[key] = np.asarray(P[key], float) * k
    P["prof"] = P["prof"].copy()
    if "s_age68" in P["prof"]:
        P["prof"]["s_age68"] = P["prof"]["s_age68"].to_numpy(float) * k
    for d in P["doms"]:
        for key in ("s68", "s75", "s76"):
            d[key] = float(d[key]) * k
        # 图上那句 "D1 440±23 Ma" 用的是这个 2σ，得跟着一起放大，
        # 否则同一条域均值线旁边会同时出现两个口径的误差。
        d["age68_2s"] = float(d["age68_2s"]) * k
    if P.get("bulk"):
        for key in ("s68", "s75", "s76"):
            P["bulk"][key] = float(P["bulk"][key]) * k
    return P


def figure_increment(recs, V):
    """图 1：全批次四面板。"""
    _style()
    doms = [(r["sample"], d) for r in recs for d in r["doms"]]
    pairs = [p for r in recs for p in r["pairs"]]
    multi3 = [r for r in recs if len(r["doms"]) >= 3]

    fig = plt.figure(figsize=(14.2, 9.4))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.24)

    # (a) 五把尺子
    ax = fig.add_subplot(gs[0, 0])
    items = [("年龄轴\nz_age", "z_age", C_AGE), ("协和平面\nz_plane", "z_plane", C_PLANE),
             ("协和度轴\nz_conc", "z_conc", C_CONC), ("化学\nz_ThU", "z_ThU", C_CHEM),
             ("化学\nz_U", "z_U", C_CHEM)]
    vals = [np.clip(np.asarray([p[k] for p in pairs], float), 0, 14) for _l, k, _c in items]
    labs = [l for l, _k, _c in items]
    cols = [c for _l, _k, c in items]
    bp = ax.boxplot(vals, widths=0.55, patch_artist=True, showfliers=False)
    for b, c in zip(bp["boxes"], cols):
        b.set_facecolor(c); b.set_alpha(0.28); b.set_edgecolor(c)
    for m in bp["medians"]:
        m.set_color(INK); m.set_linewidth(1.6)
    rng = np.random.default_rng(7)
    for i, (v, c) in enumerate(zip(vals, cols), 1):
        j = i + rng.uniform(-0.16, 0.16, len(v))
        ax.plot(j, v, "o", ms=3.0, mfc=c, mec="white", mew=0.4, alpha=0.65, ls="none")
    ax.axhline(3.0, color=C_WARN, lw=1.1, ls="--")
    ax.annotate("3σ = 真分开", xy=(0.02, 0.235), xycoords="axes fraction",
                fontsize=8.5, color=C_WARN)
    ax.set_xticklabels(labs, fontsize=8.6)
    ax.set_ylabel("相邻域之间：等效 σ 数", fontsize=10)
    ax.set_ylim(0, 14.6)
    ax.set_title("(a) 分域本身，分辨力并没有提高\n"
                 f"协和平面 ≈ 年龄轴（中位 {V['PLANE_OVER_AGE']}×）；"
                 f"协和度轴几乎判不出（中位 {V['Z_CONC_MED']}）", fontsize=10.5, loc="left")
    ax.grid(alpha=0.2, axis="y", lw=0.5)

    # (b) 域级协和度（误差棒按 K 放大；着色与标注用同一个口径的 z）
    ax = fig.add_subplot(gs[0, 1])
    xs, ys, es, cs, zs, ns = [], [], [], [], [], []
    for s, d in doms:
        if not (np.isfinite(d["conc_sigma"]) and d["conc_sigma"] > 0):
            continue
        e = d["conc_sigma"] * K                       # ← 与表里同一口径
        z = (d["concordance"] - 100.0) / e
        xs.append(d["age68"]); ys.append(d["concordance"]); es.append(e)
        zs.append(z); ns.append(f"{s} {d['name']}")
        cs.append(C_WARN if abs(z) >= 3 else ("#E8A33D" if abs(z) >= 2 else "#7F7F7F"))
    xs, ys, es = np.array(xs), np.array(ys), np.array(es)
    zs = np.array(zs, float)
    ax.axhline(100.0, color="#333333", lw=1.2, ls="--", zorder=1)
    ax.fill_between([320, 1000], 100 - 2 * np.median(es), 100 + 2 * np.median(es),
                    color="#333333", alpha=0.08, zorder=0,
                    label=f"100% ± 2×(域 1σ 中位 {np.median(es):.1f}%)")
    ax.errorbar(xs, ys, yerr=es, fmt="none", ecolor="#BBBBBB", elinewidth=0.8, zorder=2)
    ax.scatter(xs, ys, s=34, c=cs, edgecolors="white", linewidths=0.6, zorder=3)
    # 只点三个最极端的（与页面正文 {{EXTREME}} 同一份名单：按校正后的 |z| 取前三）
    for i in np.argsort(-np.abs(zs))[:3]:
        ax.annotate(ns[i], xy=(xs[i], ys[i]), xytext=(6, 5),
                    textcoords="offset points", fontsize=8.2, color=C_WARN)
    ax.set_xlabel("该域的 206Pb/238U 年龄 (Ma)", fontsize=9.5)
    ax.set_ylabel("该域的协和度 (%)", fontsize=9.5)
    ax.set_ylim(20, 185)
    ax.set_title("(b) 谐和图的真正增量：单个域的协和度\n"
                 f"{V['N_DOM']} 个域里 |z|≥2 有 {V['DEV_GE2_K']} 个、"
                 f"|z|≥3 有 {V['DEV_GE3_K']} 个（σ 已校正）—— 坪图这一维完全看不见",
                 fontsize=10.5, loc="left")
    ax.grid(alpha=0.2, lw=0.5)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    # (c) 协和度沿深度
    ax = fig.add_subplot(gs[1, 0])
    k = 0
    for r in multi3:
        c = [d["concordance"] for d in r["doms"]]
        inc = all(b > a for a, b in zip(c, c[1:]))
        dec = all(b < a for a, b in zip(c, c[1:]))
        k += 1
        mono = inc or dec
        ax.plot(range(1, len(c) + 1), c, "-o", ms=5,
                color=C_PLANE if mono else C_AGE,
                lw=2.0 if mono else 1.3, alpha=0.9 if mono else 0.75,
                zorder=4 if mono else 3)
        ax.annotate(r["sample"], xy=(len(c), c[-1]),
                    xytext=(5, [2.4, -2.8, 5.2, -5.6][k % 4]),
                    textcoords="offset points", fontsize=7.6,
                    color=C_PLANE if mono else C_AGE, va="center")
    ax.axhline(100.0, color="#333333", lw=1.2, ls="--", zorder=1)
    ax.plot([], [], "-o", color=C_PLANE, lw=2.0,
            label=f"单调（{V['MONO_N']} 个）—— 与单一混合轨迹相容")
    ax.plot([], [], "-o", color=C_AGE, lw=1.3,
            label=f"非单调（{V['NONMONO_N']} 个）—— 混合轨迹造不出")
    ax.set_xlabel("域序号（浅 → 深）", fontsize=9.5)
    ax.set_ylabel("域均值协和度 (%)", fontsize=9.5)
    ax.set_xticks(range(1, 5))
    ax.set_xlim(0.7, 4.9)
    ax.set_title("(c) 3 域及以上的测点：协和度沿深度的走向\n"
                 "普通铅/单次 Pb 丢失只会让它单调漂移 —— 非单调的那一半是真结构",
                 fontsize=10.5, loc="left")
    ax.grid(alpha=0.2, lw=0.5)
    ax.legend(fontsize=8.2, loc="best", framealpha=0.9)

    # (d) 整段协和度
    ax = fig.add_subplot(gs[1, 1])
    for r in recs:
        if not r["bulk"]:
            continue
        a, c = r["bulk"]["age68"], r["bulk"]["concordance"]
        bad = (c < 90 or c > 110)
        ax.plot(a, c, "o", ms=7 if bad else 5, mfc=C_WARN if bad else "#7F7F7F",
                mec="white", mew=0.7, alpha=0.95 if bad else 0.6, ls="none")
        if bad:
            ax.annotate(f"{r['sample']} {c:.1f}%", xy=(a, c), xytext=(6, 4),
                        textcoords="offset points", fontsize=8.0, color=C_WARN)
    ax.axhline(100.0, color="#333333", lw=1.2, ls="--")
    ax.set_xlim(320, 900)
    ax.set_xlabel("整段积分 206Pb/238U 年龄 (Ma)", fontsize=9.5)
    ax.set_ylabel("整段积分协和度 (%)", fontsize=9.5)
    ax.set_title(f"(d) 整段（结果表）：中位 {V['BULK_MED']}%，"
                 f"但 {V['BULK_OFF']} 个测点偏离 100% 超 10%\n"
                 "这些点的坪图看着正常，年龄的协和前提却已经不成立",
                 fontsize=10.5, loc="left")
    ax.grid(alpha=0.2, lw=0.5)

    fig.suptitle("一段一段地读：它多给出了什么、没多给出什么", fontsize=13.5)
    return _png_b64(fig)


def figure_contrast(sess, V):
    """图 2：S10（非单调）与 S06（单调）的 Wetherill + 协和度–深度对照。"""
    _style()
    picks = [("S10", "非单调（老–新–老）：单条混合轨迹解释不了"),
             ("S06", "单调升：与单一混合轨迹相容")]
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.2),
                             gridspec_kw=dict(hspace=0.42, wspace=0.22))
    for row, (name, note) in enumerate(picks):
        P = _inflate(sc.prepare(sess["cfg"], name), K)
        ax = axes[row, 0]
        cx, cy = P["r75"], P["r68"]
        cxx, cyy, _t = sc._curve_wetherill(sc._tmax_ma(P["prof"]["age68"]))
        ax.plot(cxx, cyy, color="#333333", lw=1.3, zorder=1, label="协和线")
        for tt in (300, 400, 500, 600):
            xv, yv = np.expm1(sc.L235 * tt * 1e6), np.expm1(sc.L238 * tt * 1e6)
            ax.annotate(f"{tt}", xy=(xv, yv), xytext=(3, -9),
                        textcoords="offset points", fontsize=7.6, color="#7A7A7A")
        used = np.zeros(len(cx), bool)
        for j, d in enumerate(P["doms"]):
            m = np.zeros(len(cx), bool); m[d["i0"]:d["i1"]] = True
            used |= m
            col = sc.DOMAIN_SPAN_COLORS[(2 * j) % len(sc.DOMAIN_SPAN_COLORS)]
            ax.plot(cx[m], cy[m], "o", ms=5.0, mfc=col, mec="white", mew=0.6,
                    ls="none", zorder=4,
                    label=f"{d['name']} {d['age68']:.0f} Ma · 协和度 {d['concordance']:.1f}%")
            ax.errorbar(d["r75"], d["r68"], xerr=2 * d["s75"], yerr=2 * d["s68"],
                        fmt="D", ms=7, mfc=col, mec="#222222", mew=0.8, ecolor=col,
                        elinewidth=1.1, capsize=3, zorder=6)
        if (~used).any():
            ax.plot(cx[~used], cy[~used], "o", ms=4.0, mfc="none", mec="#8A8A8A",
                    mew=0.9, ls="none", zorder=2, label="过渡带窗口")
        if P["bulk"]:
            b = P["bulk"]
            ax.errorbar(b["r75"], b["r68"], xerr=2 * b["s75"], yerr=2 * b["s68"],
                        fmt="*", ms=16, mfc="#FFD24C", mec="#5A4A00", mew=1.0,
                        ecolor="#5A4A00", zorder=7,
                        label=f"整段 {b['concordance']:.1f}%")
        ax.set_xlim(float(np.nanmin(cx)) * 0.94, float(np.nanmax(cx)) * 1.06)
        ax.set_ylim(float(np.nanmin(cy)) * 0.94, float(np.nanmax(cy)) * 1.06)
        ax.set_xlabel("207Pb/235U", fontsize=9.5)
        ax.set_ylabel("206Pb/238U", fontsize=9.5)
        ax.set_title(f"{name} —— Wetherill 协和图（域点 + 整段点）", fontsize=10.2, loc="left")
        ax.grid(alpha=0.22, lw=0.5)
        ax.legend(fontsize=7.4, loc="best", framealpha=0.9)

        ax = axes[row, 1]
        t = P["t"]
        a68 = P["prof"]["age68"].to_numpy(float)
        sa68 = P["prof"]["s_age68"].to_numpy(float)
        a75 = sc.age75(P["r75"])
        sa75 = np.abs(sc.age75(P["r75"] + P["s75"])
                      - sc.age75(P["r75"] - P["s75"])) / 2.0
        with np.errstate(divide="ignore", invalid="ignore"):
            conc = a75 / a68 * 100.0
            rel = np.hypot(sa75 / a75, sa68 / a68)
        for j, (lo, hi) in enumerate(P["segs"]):
            ax.axvspan(t[lo], t[max(hi - 1, lo)], alpha=0.20, zorder=0,
                       color=sc.DOMAIN_SPAN_COLORS[j % len(sc.DOMAIN_SPAN_COLORS)])
        ax.axhline(100.0, color="#333333", lw=1.2, ls="--", zorder=1, label="100%（协和）")
        ax.errorbar(t, conc, yerr=2 * conc * rel, fmt="o-", ms=3.4, lw=1.0,
                    color=C_CONC, ecolor="#C9A3E0", capsize=2, zorder=3,
                    label="逐窗口协和度 (±2σ)")
        for j, d in enumerate(P["doms"]):
            col = sc.DOMAIN_SPAN_COLORS[(2 * j) % len(sc.DOMAIN_SPAN_COLORS)]
            ax.hlines(d["concordance"], t[d["i0"]], t[max(d["i1"] - 1, d["i0"])],
                      color=col, lw=2.2, zorder=4, label=f"{d['name']} {d['concordance']:.1f}%")
        ax.set_xlabel("剥蚀时间 (s) —— 由浅到深", fontsize=9.5)
        ax.set_ylabel("协和度 (%)", fontsize=9.5)
        ax.set_title(f"{name} —— 协和度沿深度：{note}", fontsize=10.2, loc="left")
        ax.grid(alpha=0.22, lw=0.5)
        ax.legend(fontsize=7.6, loc="best", framealpha=0.9, ncol=2)
    fig.suptitle("协和度沿深度的走向：同一批窗口，多出来的一条可判读的轴", fontsize=13)
    return _png_b64(fig)


def figure_extreme(sess, spot, V):
    """图 3：一个测点的三栏 —— 坪图（看着正常）／协和平面／协和度沿深度。"""
    _style()
    P = _inflate(sc.prepare(sess["cfg"], spot), K)
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.2))
    sc._panel_age(axes[0], P)
    axes[0].set_xlabel("剥蚀时间 (s) —— 由浅到深", fontsize=9.5)
    sc._panel_concordia(axes[1], P, "wetherill", 2.0)
    sc._panel_concordance(axes[2], P, 2.0)
    fig.suptitle(f"{spot}：坪图上它平得像教科书，协和度却整段偏离"
                 f"（整段 {P['bulk']['concordance']:.1f}%；误差棒已按 ×{K:g} 校正）",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _png_b64(fig)


# ─────────────────────────────────────────────────────────────────────────────
# 四、注入
# ─────────────────────────────────────────────────────────────────────────────
def render(recs, figs_model, V):
    V = dict(V)
    V["FIG1"] = _fig_frame(
        figs_model["fig1"],
        "<b>图 1</b>　示例批次全览（所有误差棒已按实测倍数放大）。<b>(a)</b> 相邻域对的分辨力："
        "五把尺子各给一个等效 σ 数（箱体=四分位、点=每个域对、虚线=3σ）。<b>(b)</b> "
        "每个正式年龄域的协和度与其自身 1σ —— 灰点为 |z|&lt;2、橙点 2–3、红点 ≥3，"
        "标注的是校正后最极端的三个域。<b>(c)</b> 三域及以上的测点，域均值协和度沿深度："
        "红线单调、蓝线非单调。<b>(d)</b> 整段积分协和度与测点年龄，红点为偏离 100% 超 10% 者。",
        "分域之后多得到什么：五把尺子的分辨力、域级协和度、协和度沿深度的走向、整段协和度")
    V["FIG2"] = _fig_frame(
        figs_model["fig2"],
        "<b>图 2</b>　同一批窗口的两条读法。<b>左</b>：域点投在 Wetherill 协和图上，"
        "菱形为域均值（±2σ），星号为整段积分。<b>右</b>：协和度沿深度，紫点为逐窗口"
        "（带 ±2σ 误差棒）、横线为域均值。上排这个测点的协和度是「高—低—高」，"
        "下排是一路单调升 —— 后者与「某一条混合轨迹」相容，前者只能解释成真实的结构。",
        "两个测点的 Wetherill 协和分布与协和度沿深度对照")
    V["FIG3"] = _fig_frame(
        figs_model["fig3"],
        "<b>图 3</b>　坪图看不见的那一半。<b>左</b>坪图：逐窗口年龄带 ±2σ（已校正），"
        "整段加权后与常数模型勉强相容，看不出任何异常；<b>中</b>协和平面：域均值点离协和线很远，"
        "域级协和度只有几十个百分点；<b>右</b>协和度沿深度：从约 30% 一路飘到约 140%，"
        "是系统性的，不是噪声。⚠ 这类测点有两种可能（真混入普通铅 / <code>204Pb</code> 校正过度），"
        "本工具分不开，需要回看 204 显著性与 CL。",
        "坪图看似正常而协和度整段偏离的测点")

    text = TEMPLATE.read_text(encoding="utf-8")
    for key, val in V.items():
        if not isinstance(val, str):
            raise SystemExit(f"占位符 {key} 的值不是字符串：{type(val).__name__}")
        text = text.replace("{{" + key + "}}", val)
    left = re.findall(r"\{\{[A-Za-z0-9_]+\}\}", text)
    if left:
        raise SystemExit(f"页面里还有没替换的占位符：{sorted(set(left))}")
    return text, V


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="生成《一段一段地读》页面的数字与图（只读）。")
    ap.add_argument("--check", action="store_true",
                    help="只重算并与已提交的页面比对，不写文件")
    ap.add_argument("--extreme-spot", default="S02",
                    help="图 3 用的测点（默认 S02：坪图正常、域级协和度崩掉的那个）")
    args = ap.parse_args(argv)

    print(f"扫描 {BATCH.name} …", flush=True)
    sess = build_session()
    recs = scan(sess)
    V = metrics(recs)
    print(f"  {len(recs)} 个样品测点；σ 校正因子 k = {K:g}；画图 …", flush=True)
    figs = dict(fig1=figure_increment(recs, V),
                fig2=figure_contrast(sess, V),
                fig3=figure_extreme(sess, args.extreme_spot, V))
    text, V = render(recs, figs, V)

    if args.check:
        cur = PAGE.read_text(encoding="utf-8")
        if cur != text:
            print("⚠ 页面与现算值不一致（把 --check 去掉可重写）。")
            for k in ("N_DOM", "Z_CONC_MED", "DEV_GE3_K", "DEV_GE3",
                      "MONO_N", "POOL_GAIN"):
                print(f"   现算 {k} = {V[k]}")
            return 1
        print(f"✓ 页面与现算值一致（由 {TEMPLATE.name} + 现算值生成）")
        return 0

    PAGE.write_text(text, encoding="utf-8")
    print(f"已按模板写出 {PAGE.relative_to(ROOT)}（{len(text) / 1024:.0f} KB）")
    for k in ("N_SPOTS", "N_DOM", "N_MULTI", "Z_AGE_MED", "Z_PLANE_MED",
              "Z_CONC_MED", "Z_U_MED", "DEV_GE3_K", "DEV_GE3", "DEV_GE2_K",
              "MONO_N", "NONMONO_N", "WS_MSWD_MED", "WS_OK", "POOL_GAIN",
              "POOL_CHI2_MED", "DOM_SE", "POOL_SE", "S75_RATIO",
              "S75_WIN_PCT", "S75_SPOT_REL", "S75_WIN_OVER_SPOT", "K"):
        print(f"   {k} = {V[k]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
