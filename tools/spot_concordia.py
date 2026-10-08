#!/usr/bin/env python
"""
tools/spot_concordia.py —— 单点的「逐深度坪图 + 协和分布图」（只读）
====================================================================

**它回答什么**

一个测点沿深度切了几十个窗口后，通常只画一张坪图（age68 vs 时间）。
坪图只回答"年龄沿深度变不变"；它回答不了另一半问题：

    · 这些窗口年龄是不是**同一个协和体系**？（有没有普通铅、有没有 Pb 丢失）
    · 两个域分开，是"真的两期生长"，还是"一个体系被污染搅开"？

把**同一批窗口**再投到协和图上，这两条轴就同时有了。所以本脚本一次给三样：

    ① 坪图        逐窗口 age68 vs 剥蚀时间（τ 由浅到深），带域底色与各域加权平均
    ② 协和图      同一批窗口的 (207Pb/235U, 206Pb/238U) 投到 Wetherill 协和图上，
                  另给一张 Tera–Wasserburg（207Pb/206Pb vs 238U/206Pb）——
                  后者是判普通铅的标准平面
    ③ 协和度–深度 age(207Pb/235U) / age(206Pb/238U) 沿深度的变化 —— 三张里最直读：
                  平坦且贴 100% = 整段一个协和体系；系统性偏高 = 过量 207Pb
                  （普通铅没扣干净 / 混进老核）；偏低 = Pb 丢失

⚠ 逐窗口的 207/235 精度很低（4 s 窗口相对 1σ 中位约 7%），**2σ 椭圆比整个
  坐标轴还宽**。所以逐窗口默认只画点，椭圆要用 `--window-ellipse` 显式打开；
  逐窗口的精度信息由 ③ 栏的误差棒与域均值点的十字棒承担。这不是偷懒 ——
  把一屏都盖住的椭圆画上去，读者只会得到"看起来什么都说不清"的印象。

**口径：一个数都不另算**

全部调用流水线自己的函数（`load_batch` → `build_bracketing` →
`compute_bulk_ratios` → `calibrate_primary` → `estimate_external` →
`analyse_depth_spot`），与出正式结果表时**同一套代码、同一套参考值口径**。
本脚本里只有"换算成图上的坐标"与画图。

    · 206/238 ：逐窗口 R68 · F68(τ)，F68 由夹逼标样在同一 τ 处给出
    · 207/206 ：**不套 F(τ)**。Pb 同位素之间的分馏几乎相同，套 F 只会把
                标样噪声注进样品（见 `depth/fractionation.py::profile_ages_76`）
    · 207/235 ：= (206/238) × (207/206) × 238U/235U，**不是直接测的**
                ⇒ 相对 1σ 必须带协方差项（`relative_sigma_product`，
                σ 与 ρ 同基：对角线用合成相对 1σ，交叉项只用内部分量配内部 ρ）

**误差椭圆的相关系数**

Wetherill 平面画的是 (207/235, 206/238)，它不是两个独立测量，ρ 不能取 0：

    ρ(206/238, 207/235) = (c68² + a68·a76·ρ) / (c68 · c75)
    ρ(238/206, 207/206) = −a68·a76·ρ / (c68 · c76)        （TW 平面）

c ＝ 合成相对 1σ，a ＝ 对应的内部分量，ρ ＝ `reduce_interval` 给的内部分基
相关系数（实测中位 −0.12，负号是构造性的）。推导见 `判读细目.md` §二.2。
**域均值点**画成十字误差棒而不是椭圆：均值点的 ρ 需要另行假设，
用十字棒如实一些（脚本里注明的近似只有一处：域均值椭圆的 ρ 取该域
窗口 ρ 的中位，仅当 `--mean-ellipse` 打开时才会用到）。

**用法**

    python tools/spot_concordia.py --dir <批次目录> --spot <样品名>
    python tools/spot_concordia.py --dir examples/EX2022A --spot S10 --out <目录>
    python tools/spot_concordia.py --dir <批次目录> --list        # 只列样品名

输出（写到 `--out`，默认「批次目录的上一级/spot_concordia」）：
    <序号>_<样品>_单点总览.png   一页答案：左坪图（年龄/U/Th-U 三栏）+ 右协和三栏
    <序号>_<样品>_谐和分布.png   只要协和三栏（放大版，便于单图放论文）
    <序号>_<样品>_深度剖面.png   官方样式的深度剖面图（调用流水线自己的画图函数）
    <序号>_<样品>_单点总览.pdf   上述拼成一份 PDF

⚠ 本脚本**只读**：不写结果表、不改任何数值输出、不碰 `结果/`。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")                                   # 批处理环境无显示器
import matplotlib.pyplot as plt                          # noqa: E402
from matplotlib.patches import Ellipse                    # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from druid.core.constants import (CJK_FONTS, L235, L238,   # noqa: E402
                                  ROLE_UNKNOWN, U238_U235)
from druid.core.geochronology import age68, age75          # noqa: E402
from druid.core.references import std_alias, std_ref       # noqa: E402
from druid.core.statistics import weighted_mean            # noqa: E402
from druid.depth.figures import (DOMAIN_SPAN_COLORS,        # noqa: E402
                                 save_depth_figure)
from druid.depth.fractionation import bracket_F            # noqa: E402
from druid.workflow import (BatchConfig, analyse_depth_spot,  # noqa: E402
                            build_bracketing, build_results,
                            calibrate_primary, compute_bulk_ratios,
                            estimate_external, load_batch,
                            whole_spot_stats)


# ─────────────────────────────────────────────────────────────────────────────
# 〇、全部落在同一个坐标系里的换算
# ─────────────────────────────────────────────────────────────────────────────
def _overlap_step(cfg: BatchConfig):
    """窗口重叠校正要的 step/win；与 `workflow._overlap_ratio` 同一判据。"""
    if not cfg.overlap_correct or cfg.win <= 0.0:
        return None
    return float(cfg.step) / float(cfg.win)


def _rho_wetherill(c68, a68, c76, a76, rho):
    """ρ(206/238, 207/235)，并顺带给出 207/235 的相对 1σ。见模块 docstring。

    与 `core.statistics.relative_sigma_product` **同一个式子**（那个函数是标量版、
    一次只能算一个点；这里要整条窗口数组，所以按同一公式写成向量式，
    公式本身不另立）。
    """
    r = np.clip(np.asarray(rho, float), -1.0, 1.0)
    c75 = np.sqrt(np.maximum(c68 ** 2 + c76 ** 2 + 2.0 * a68 * a76 * r, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = (c68 ** 2 + a68 * a76 * r) / (c68 * c75)
    return c75, np.clip(rr, -1.0, 1.0)


def _rho_tw(a68, a76, rho, c68, c76):
    """ρ(238/206, 207/206)：只是 206/238 取倒数，协方差变号。"""
    with np.errstate(divide="ignore", invalid="ignore"):
        r = -(a68 * a76 * rho) / (c68 * c76)
    return np.clip(r, -1.0, 1.0)


def _curve_wetherill(t_hi_ma: float, n: int = 700):
    """协和线上各时刻的 (207Pb/235U, 206Pb/238U)。用 expm1 保住小量精度。"""
    t = np.linspace(0.0, t_hi_ma * 1e6, n)               # 年
    return np.expm1(L235 * t), np.expm1(L238 * t), t / 1e6


def _curve_tw(t_hi_ma: float, n: int = 700):
    """协和线上各时刻的 (238U/206Pb, 207Pb/206Pb)。"""
    r75, r68, t_ma = _curve_wetherill(t_hi_ma, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 / r68, r75 / (U238_U235 * r68), t_ma


def _ellipse(ax, x, y, sx, sy, rho, nsigma, **kw):
    """按协方差矩阵画误差椭圆（nsigma σ）。任何非有限量都直接跳过。"""
    if not all(np.isfinite([x, y, sx, sy])) or sx <= 0 or sy <= 0:
        return None
    rho = float(np.clip(rho, -0.999, 0.999))
    cov = np.array([[sx * sx, rho * sx * sy], [rho * sx * sy, sy * sy]])
    val, vec = np.linalg.eigh(cov)
    order = np.argsort(val)[::-1]
    val, vec = val[order], vec[:, order]
    ang = float(np.degrees(np.arctan2(vec[1, 0], vec[0, 0])))
    w, h = 2.0 * nsigma * np.sqrt(np.maximum(val, 0.0))
    e = Ellipse((x, y), width=w, height=h, angle=ang, fill=False, **kw)
    ax.add_patch(e)
    return e


# ─────────────────────────────────────────────────────────────────────────────
# 一、把一个测点算成图上要用的数组
# ─────────────────────────────────────────────────────────────────────────────
def prepare(cfg: BatchConfig, spot: str):
    """
    重跑批次前置步骤，取回目标测点的逐窗口数据与域划分。

    返回 dict：
        prof / segs / summ / tag       —— 与 `analyse_depth_spot` 返回的同一批对象
        tau, t                                 窗口坐标
        r68, s68, r76, s76, r75, s75          校正后的比值与 1σ（206/238 已套 F(τ)）
        rho_w, rho_tw                          两个平面各自的相关系数
        doms                                   每个正式年龄域的均值点（比值空间）
        bulk                                   结果表里的"整段积分"点
        ws                                     whole_spot_stats（不分域窗口加权）
    """
    spots, _seq, _skipped = load_batch(cfg)
    brack_for = build_bracketing(spots, cfg)
    ref68, ref76 = std_ref(std_alias(cfg.primary) or cfg.primary, cfg.ref_preset)

    _variants, (b68, b76) = compute_bulk_ratios(spots, brack_for, ref68, cfg)
    pidx, F68, F76, _cal = calibrate_primary(spots, b68, b76, ref68, ref76, cfg)
    sd68, sd76 = estimate_external(spots, pidx, F68, F76, b68, b76, cfg)

    hit = None
    for k, (tr, r) in enumerate(spots):
        if tr.role == ROLE_UNKNOWN and tr.sample == spot:
            hit = (k, tr, r)
            break
    if hit is None:
        names = sorted({tr.sample for tr, _ in spots if tr.role == ROLE_UNKNOWN})
        raise SystemExit(f"批次里没有样品 『{spot}』。可用样品名：{', '.join(names)}")

    k, tr, r = hit
    brack = brack_for(k)
    got = analyse_depth_spot(tr, r, brack, ref68, sd68, cfg)
    if got is None:
        raise SystemExit(f"测点 {spot} 无法做深度剖面（窗口为空或没有可用夹逼标样）。")
    prof, segs, summ, tag = got

    tau = prof["tau"].to_numpy(float)
    # 206/238：套 F(τ)。与 analyse_depth_spot 内部逐位同一条路径 ——
    # 下面这行 assert 就是钉这一点的（差一位就说明口径被改过）。
    F = bracket_F(tau, brack, "R68", ref68)
    r68 = prof["R68"].to_numpy(float) * F
    s68 = np.hypot(prof["s68"].to_numpy(float) * F, sd68 * r68)
    if not np.array_equal(age68(r68), prof["age68"].to_numpy(float)):
        raise AssertionError("逐窗口年龄与流水线不一致：F(τ) 口径被改动过？")

    # 207/206：**不**套 F(τ)（理由见模块 docstring）
    r76 = prof["R76"].to_numpy(float)
    s76 = np.hypot(prof["s76"].to_numpy(float), sd76 * r76)

    # 207/235：两个比值的乘积，相对 1σ 带协方差项（严格传播）
    rho = prof["rho"].to_numpy(float)
    c68 = s68 / r68
    c76 = s76 / r76
    a68 = prof["s68"].to_numpy(float) * F / r68            # 内部分量
    a76 = prof["s76"].to_numpy(float) / r76
    c75, rho_w = _rho_wetherill(c68, a68, c76, a76, rho)
    r75 = r76 * r68 * U238_U235
    s75 = c75 * r75
    rho_tw = _rho_tw(a68, a76, rho, c68, c76)

    # 各正式年龄域的均值点（比值空间加权，不是把年龄平均后再反算）
    sow = _overlap_step(cfg)
    doms = []
    for _, d in summ.iterrows():
        if str(d["flag"]) != "age domain":
            continue
        i0, i1 = int(d["i0"]), int(d["i1"])
        mu68, se68, mswd68, n68 = weighted_mean(r68[i0:i1], s68[i0:i1], sow)
        mu76, se76, _mswd76, _n76 = weighted_mean(r76[i0:i1], s76[i0:i1], sow)
        mu75 = mu76 * mu68 * U238_U235
        se75 = mu75 * np.hypot(se68 / mu68, se76 / mu76)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.hypot(se68 / mu68, se76 / mu76)
            rho_m = np.nanmedian(rho_w[i0:i1]) if rel > 0 else 0.0
        doms.append(dict(
            name=str(d["domain"]), i0=i0, i1=i1, n_win=int(d["n_win"]),
            tau0=float(d["tau0"]), tau1=float(d["tau1"]),
            r68=mu68, s68=se68, r76=mu76, s76=se76, r75=mu75, s75=se75,
            rho=float(rho_m),
            age68=float(d["age_Ma"]), age68_2s=2.0 * float(d["se_1sig"]),
            mswd=float(d["mswd"]), mswd_prob=float(d["mswd_prob"]),
            age75=float(age75(np.array([mu75]))[0]),
            concordance=float(age75(np.array([mu75]))[0]) / float(d["age_Ma"]) * 100.0,
            th_u=float(d["ThU"]),
        ))

    # 结果表里的"整段积分"点：与域级同一套 F，只差"不切窗口"
    tbl = build_results(spots, pidx, F68, F76, b68, b76, sd68, sd76, cfg)
    row = tbl[tbl["序号"] == int(tr.idx) + 1]
    bulk = None
    if len(row):
        rr = row.iloc[0]
        c68b = float(rr["s68_pct"]) / 100.0
        c76b = float(rr["s76_pct"]) / 100.0
        # 结果表的 1σ 是 hypot(内部, 外部)，反解内部分量（同 `build_results` 公式）
        a68b = np.sqrt(max(c68b ** 2 - float(sd68) ** 2, 0.0))
        a76b = np.sqrt(max(c76b ** 2 - float(sd76) ** 2, 0.0))
        rhob = float(rr["rho_68_76"])
        c75b, rho_wb = _rho_wetherill(c68b, a68b, c76b, a76b, rhob)
        r75b = float(rr["Pb207_235U"])
        s75b = c75b * r75b                       # 比值空间的 1σ（画图要用这个）
        bulk = dict(
            r68=float(rr["Pb206_238U"]), s68=c68b * float(rr["Pb206_238U"]),
            r76=float(rr["Pb207_206Pb"]), s76=c76b * float(rr["Pb207_206Pb"]),
            r75=r75b, s75=s75b,
            rho_w=float(rho_wb),
            rho_tw=float(_rho_tw(a68b, a76b, rhob, c68b, c76b)),
            age68=float(rr["年龄206_238"]), age75=float(rr["年龄207_235"]),
            age76=float(rr["年龄207_206"]),
            # 旧口径（只搬 206/238 那一项）与新口径的比 —— 就是那个 2.69 倍。
            # ⚠ 两者都换成 **Ma** 再打印：比值空间的 σ 与年龄空间的 σ 差着
            #   λ·(1+R) 这个因子，混起来会得到一个荒谬的比值（第一版就踩了）。
            s75_legacy_ma=float(rr["s75_2sig"]) / 2.0,
            s75_strict_ma=float(abs(age75(r75b + s75b) - age75(r75b - s75b)) / 2.0),
            concordance=float(rr["协和度_pct"]),
        )

    ws = whole_spot_stats(prof, step_over_win=sow)
    return dict(prof=prof, segs=segs, summ=summ, tag=tag, tau=tau, t=prof["t_mid"].to_numpy(float),
                r68=r68, s68=s68, r76=r76, s76=s76, r75=r75, s75=s75,
                rho_w=rho_w, rho_tw=rho_tw, doms=doms, bulk=bulk, ws=ws,
                label=f"{int(tr.idx) + 1:02d} {tr.sample}", sample=tr.sample,
                win=float(cfg.win), step=float(cfg.step),
                ref68=ref68, sd68=float(sd68), sd76=float(sd76))


# ─────────────────────────────────────────────────────────────────────────────
# 二、画图
# ─────────────────────────────────────────────────────────────────────────────
def _style():
    plt.rcParams["font.sans-serif"] = CJK_FONTS
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.facecolor"] = "white"
    return plt


def _panel_age(ax, P, annotate=True):
    """坪图主栏：逐窗口 age68（±2σ）+ 域底色 + 各域加权平均线。"""
    t, segs = P["t"], P["segs"]
    if len(segs) > 1:
        for j, (lo, hi) in enumerate(segs):
            ax.axvspan(t[lo], t[max(hi - 1, lo)], alpha=0.30, zorder=0,
                       color=DOMAIN_SPAN_COLORS[j % len(DOMAIN_SPAN_COLORS)])
    ax.errorbar(t, P["prof"]["age68"], yerr=2 * P["prof"]["s_age68"], fmt="o-",
                ms=3.4, lw=1.1, color="#185FA5", ecolor="#9FC3E8", capsize=2,
                zorder=3, label="206Pb/238U 逐窗口 (±2σ)")
    for j, d in enumerate(P["doms"]):
        ax.hlines(d["age68"], t[d["i0"]], t[max(d["i1"] - 1, d["i0"])],
                  color="#A32D2D", lw=1.8, zorder=4,
                  label=None if j else "域加权平均")
        if annotate:
            ax.annotate(f"{d['name']} {d['age68']:.0f}±{d['age68_2s']:.0f} Ma",
                        xy=(0.5 * (t[d["i0"]] + t[max(d["i1"] - 1, d["i0"])]), d["age68"]),
                        xytext=(0, 4.0), textcoords="offset points",
                        fontsize=8.5, color="#A32D2D", ha="center", va="bottom", zorder=5)
    ws = P["ws"]
    if np.isfinite(ws.get("age_Ma", np.nan)):
        ax.axhline(float(ws["age_Ma"]), color="#4A4A4A", lw=1.4, ls=(0, (6, 3)),
                   zorder=3.5,
                   label="不分域窗口加权 %.0f Ma（MSWD %.3g）"
                         % (float(ws["age_Ma"]), float(ws["mswd"])))
    ax.set_ylabel("年龄 (Ma)", fontsize=10)
    ax.legend(fontsize=8.2, loc="best", framealpha=0.9)


def _panel_concordia(ax, P, plane: str, nsigma: float, mean_ellipse=False,
                     window_ellipse=False):
    """协和分布图。plane = 'wetherill' | 'tw'。"""
    if plane == "wetherill":
        cx, cy = P["r75"], P["r68"]
        sx, sy, sr = P["s75"], P["s68"], P["rho_w"]
        bx, by, bsx, bsy, br = (P["bulk"]["r75"], P["bulk"]["r68"],
                                P["bulk"]["s75"], P["bulk"]["s68"],
                                P["bulk"]["rho_w"]) if P["bulk"] else (None,) * 5
        curve = _curve_wetherill(_tmax_ma(P["prof"]["age68"]))
        xlab, ylab = "207Pb/235U", "206Pb/238U"
        title = "Wetherill 协和图（逐深度窗口）"
        # 协和线上的等时刻度：标出这些年龄的位置，便于直接读图上落在哪
        ticks = [100, 200, 300, 500, 700, 1000, 1500]
        for tt in ticks:
            xv = np.expm1(L235 * tt * 1e6)
            yv = np.expm1(L238 * tt * 1e6)
            if xv <= 1.6 * np.nanmax(cx) and yv <= 2.0 * np.nanmax(cy):
                ax.annotate(f"{tt}", xy=(xv, yv), xytext=(3, -8),
                            textcoords="offset points", fontsize=7.5, color="#7A7A7A")
    else:
        cx, cy = 1.0 / P["r68"], P["r76"]
        sx = P["s68"] / P["r68"] ** 2            # σ(1/R68) = σ68 / R68²
        sy, sr = P["s76"], P["rho_tw"]
        bx, by, bsx, bsy, br = (1.0 / P["bulk"]["r68"], P["bulk"]["r76"],
                                P["bulk"]["s68"] / P["bulk"]["r68"] ** 2,
                                P["bulk"]["s76"], P["bulk"]["rho_tw"]) \
            if P["bulk"] else (None,) * 5
        curve = _curve_tw(_tmax_ma(P["prof"]["age68"]))
        xlab, ylab = "238U/206Pb", "207Pb/206Pb"
        title = "Tera–Wasserburg 图（共用一个 206Pb，判普通铅的标准平面）"

    cxx, cyy, ctt = curve
    ax.plot(cxx, cyy, color="#333333", lw=1.3, zorder=1, label="协和线")

    for j, d in enumerate(P["doms"]):
        m = np.zeros(len(cx), bool)
        m[d["i0"]:d["i1"]] = True
        col = DOMAIN_SPAN_COLORS[(2 * j) % len(DOMAIN_SPAN_COLORS)]
        dx = (d["r75"], d["r68"]) if plane == "wetherill" else (1.0 / d["r68"], d["r76"])
        dsx = d["s75"] if plane == "wetherill" else d["s68"] / d["r68"] ** 2
        dsy = d["s68"] if plane == "wetherill" else d["s76"]
        if mean_ellipse:
            _ellipse(ax, dx[0], dx[1], dsx, dsy, d["rho"], nsigma,
                     color=col, lw=1.8, ls="--", alpha=0.95, zorder=5)
        ax.errorbar(dx[0], dx[1], xerr=nsigma * dsx, yerr=nsigma * dsy,
                    fmt="D", ms=7, mfc=col, mec="#222222", mew=0.8,
                    ecolor=col, elinewidth=1.2, capsize=3, zorder=6,
                    label=f"{d['name']}  {d['age68']:.0f} Ma（{d['n_win']} 窗）")
        # 该域的逐窗口点；椭圆默认**不画** —— 4 s 窗口的 207/235 相对 1σ
        # 中位 7%，2σ 椭圆比整个坐标轴还宽，全画出来会把图糊成一片。
        # 逐窗口的精度信息由右下的「协和度–深度」栏与域均值十字棒承担。
        if window_ellipse:
            for i in np.flatnonzero(m):
                _ellipse(ax, cx[i], cy[i], sx[i], sy[i], sr[i], nsigma,
                         color=col, lw=0.7, alpha=0.5, zorder=3)
        ax.plot(cx[m], cy[m], "o", ms=4.6, mfc=col, mec="white", mew=0.6,
                ls="none", zorder=4)

    # 不属于任何正式年龄域的窗口（过渡带 / 未通过相容性检验）画成空心点
    used = np.zeros(len(cx), bool)
    for d in P["doms"]:
        used[d["i0"]:d["i1"]] = True
    if (~used).any():
        ax.plot(cx[~used], cy[~used], "o", ms=4.0, mfc="none", mec="#8A8A8A",
                mew=0.9, ls="none", zorder=2, label="未归入年龄域（过渡带）")

    if P["bulk"] is not None:
        ax.errorbar(bx, by, xerr=nsigma * bsx, yerr=nsigma * bsy, fmt="*",
                    ms=15, mfc="#FFD24C", mec="#5A4A00", mew=1.0,
                    ecolor="#5A4A00", elinewidth=1.4, capsize=3, zorder=7,
                    label="整段积分（结果表）")

    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_title(title, fontsize=10.5)
    ax.grid(alpha=0.22, lw=0.5)
    ax.legend(fontsize=7.8, loc="best", framealpha=0.9)
    _autoscale_with_curve(ax, cxx, cyy, [cx] + ([bx] if P["bulk"] else []),
                          [cy] + ([by] if P["bulk"] else []))


def _tmax_ma(ages):
    a = np.asarray(ages, float)
    a = a[np.isfinite(a)]
    return float(np.nanmax(a)) * 1.35 if a.size else 1000.0


def _autoscale_with_curve(ax, cx, cy, xs_list, ys_list):
    """把数据与协和线一起纳入视野 —— 协和线只在数据量级附近才有意义。"""
    xs = np.concatenate([np.atleast_1d(np.asarray(v, float)).ravel()
                         for v in xs_list if np.size(v)])
    ys = np.concatenate([np.atleast_1d(np.asarray(v, float)).ravel()
                         for v in ys_list if np.size(v)])
    xs = xs[np.isfinite(xs)]
    ys = ys[np.isfinite(ys)]
    if not xs.size or not ys.size:
        return
    x0, x1 = float(np.min(xs)), float(np.max(xs))
    y0, y1 = float(np.min(ys)), float(np.max(ys))
    dx, dy = max(x1 - x0, 1e-9), max(y1 - y0, 1e-9)
    xlo, xhi = x0 - 0.18 * dx, x1 + 0.18 * dx
    ylo, yhi = y0 - 0.18 * dy, y1 + 0.18 * dy
    cx = np.asarray(cx, float)
    cy = np.asarray(cy, float)
    m = (cx >= xlo) & (cx <= xhi) & np.isfinite(cy)
    if m.any():
        ylo = min(ylo, float(np.nanmin(cy[m])) - 0.10 * dy)
        yhi = max(yhi, float(np.nanmax(cy[m])) + 0.10 * dy)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)


def _panel_concordance(ax, P, nsigma):
    """
    协和度沿深度：age(207Pb/235U) / age(206Pb/238U) × 100%。

    这一栏是三个协和视角里**最直读**的一个：
        · 平坦且贴 100%  → 整段是一个协和体系
        · 系统性偏高     → 有过量 207Pb（普通铅没扣干净 / 混进老核）
        · 系统性偏低     → 有 Pb 丢失（或 206Pb 被表面污染抬高）
    """
    t = P["t"]
    if len(P["segs"]) > 1:
        for j, (lo, hi) in enumerate(P["segs"]):
            ax.axvspan(t[lo], t[max(hi - 1, lo)], alpha=0.30, zorder=0,
                       color=DOMAIN_SPAN_COLORS[j % len(DOMAIN_SPAN_COLORS)])
    a68 = P["prof"]["age68"].to_numpy(float)
    sa68 = P["prof"]["s_age68"].to_numpy(float)
    a75 = age75(P["r75"])
    sa75 = np.abs(age75(P["r75"] + P["s75"]) - age75(P["r75"] - P["s75"])) / 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        conc = a75 / a68 * 100.0
        rel = np.hypot(sa75 / a75, sa68 / a68)
    ax.axhline(100.0, color="#333333", lw=1.2, ls="--", zorder=1,
               label="100%（协和）")
    ax.errorbar(t, conc, yerr=nsigma * conc * rel, fmt="o-", ms=3.4, lw=1.0,
                color="#7A3FA5", ecolor="#C9A3E0", capsize=2, zorder=3,
                label=f"逐窗口协和度 (mean ±{nsigma:g}σ)")
    for j, d in enumerate(P["doms"]):
        col = DOMAIN_SPAN_COLORS[(2 * j) % len(DOMAIN_SPAN_COLORS)]
        ax.hlines(d["concordance"], t[d["i0"]], t[max(d["i1"] - 1, d["i0"])],
                  color=col, lw=1.8, zorder=4,
                  label=None if j else "域均值协和度")
    ax.set_ylabel("协和度 (%)", fontsize=9.5)
    ax.set_xlabel("剥蚀时间 (s) —— 由浅到深", fontsize=9.5)
    ax.grid(alpha=0.22, lw=0.5)
    ax.legend(fontsize=7.6, loc="best", framealpha=0.9)


def draw_overview(P, out_png: Path, nsigma: float, dpi: int = 150,
                  mean_ellipse=False, window_ellipse=False):
    """一页答案：左（年龄 / 238U / Th-U 三栏）+ 右（Wetherill / TW / 协和度）。"""
    _style()
    fig = plt.figure(figsize=(15.5, 10.6))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.19)
    left = outer[0].subgridspec(3, 1, height_ratios=[2.6, 1.0, 1.0], hspace=0.10)
    right = outer[1].subgridspec(3, 1, hspace=0.52)

    ax0 = fig.add_subplot(left[0])
    ax1 = fig.add_subplot(left[1], sharex=ax0)
    ax2 = fig.add_subplot(left[2], sharex=ax0)
    _panel_age(ax0, P)
    ax0.set_title(f"{P['label']} —— 逐深度窗口 [{P['tag']}]   "
                  f"窗口 {P['win']:g}s / 步长 {P['step']:g}s", fontsize=11.5)
    prof = P["prof"]
    ax1.plot(P["t"], prof["U_cps"], color="#534AB7", lw=1.3)
    ax1.set_ylabel("238U (cps)", fontsize=9.5)
    ax2.plot(P["t"], prof["ThU"], color="#0F6E56", lw=1.3)
    ax2.set_ylabel("Th/U", fontsize=9.5)
    ax2.set_xlabel("剥蚀时间 (s) —— 由浅到深", fontsize=9.5)

    _panel_concordia(fig.add_subplot(right[0]), P, "wetherill", nsigma,
                     mean_ellipse, window_ellipse)
    _panel_concordia(fig.add_subplot(right[1]), P, "tw", nsigma,
                     mean_ellipse, window_ellipse)
    _panel_concordance(fig.add_subplot(right[2]), P, nsigma)

    fig.suptitle(f"{P['label']} —— 单点逐深度：坪图 + 协和分布（椭圆/误差棒 = {nsigma:g}σ）",
                 fontsize=13)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    return fig, out_png


def draw_concordia_only(P, out_png: Path, nsigma: float, dpi: int = 160,
                        mean_ellipse=False, window_ellipse=False):
    """只要协和三图（放大版）。"""
    _style()
    fig, axes = plt.subplots(1, 3, figsize=(19.0, 5.6))
    _panel_concordia(axes[0], P, "wetherill", nsigma, mean_ellipse, window_ellipse)
    _panel_concordia(axes[1], P, "tw", nsigma, mean_ellipse, window_ellipse)
    _panel_concordance(axes[2], P, nsigma)
    fig.suptitle(f"{P['label']} 逐深度窗口的协和分布（误差棒/椭圆 = {nsigma:g}σ）",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    return fig, out_png


# ─────────────────────────────────────────────────────────────────────────────
# 三、控制台表
# ─────────────────────────────────────────────────────────────────────────────
def print_report(P):
    """把图上的数同时打成文字表 —— 图上读不准，表里能核。"""
    print()
    print("=" * 100)
    print(f"单点逐深度复核：{P['label']}   结构 [{P['tag']}]   "
          f"窗口 {len(P['prof'])} 个   σext(206/238) = {P['sd68'] * 100:.2f}% / "
          f"σext(207/206) = {P['sd76'] * 100:.2f}%")
    print("=" * 100)
    print("%-5s %-11s %4s %8s %10s %8s %10s %9s %7s" %
          ("域", "τ 区间", "窗数", "age68", "±2σ", "MSWD", "age207/235", "协和度", "Th/U"))
    print("-" * 100)
    for d in P["doms"]:
        print("%-5s %-11s %4d %8.1f %10.1f %8.2f %10.1f %8.1f%% %7.3f" %
              (d["name"], "%.2f-%.2f" % (d["tau0"], d["tau1"]), d["n_win"],
               d["age68"], d["age68_2s"], d["mswd"], d["age75"], d["concordance"],
               d["th_u"]))
    ws = P["ws"]
    if np.isfinite(ws.get("age_Ma", np.nan)):
        print("-" * 100)
        print("整段不分域（窗口加权）  age68 = %.1f ± %.1f Ma (2σ)   MSWD = %.2f   %s"
              % (float(ws["age_Ma"]), 2 * float(ws["se_1sig"]), float(ws["mswd"]),
                 "与常数相容" if ws.get("mswd_ok") else "与常数不相容"))
    b = P["bulk"]
    if b:
        print("整段积分（结果表）      age68 = %.1f Ma   age75 = %.1f Ma   "
              "age76 = %.1f Ma   协和度 = %.1f%%"
              % (b["age68"], b["age75"], b["age76"], b["concordance"]))
        print("                       207/235 的 1σ：旧口径 %.2f Ma → 严格传播 %.2f Ma"
              "（%.2f×，含 ρ 交叉项）"
              % (b["s75_legacy_ma"], b["s75_strict_ma"],
                 b["s75_strict_ma"] / b["s75_legacy_ma"]
                 if b["s75_legacy_ma"] else float("nan")))
    print("-" * 100)
    print("读法：协和图上同一测点的窗口若沿一条斜率 1 的线分布 → 单个协和体系；")
    print("      向 207/206 高值一侧散开 → 有普通铅（或老核混入）；向低值一侧 → Pb 丢失。")
    print("      椭圆只描述该窗口的计数统计，不包含窗口之间的系统项。")


# ─────────────────────────────────────────────────────────────────────────────
# 四、命令行
# ─────────────────────────────────────────────────────────────────────────────
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="单点的逐深度窗口：坪图 + 协和分布图（只读，不改任何数据）。")
    ap.add_argument("--dir", default=str(ROOT / "examples" / "EX2022A"),
                    help="批次目录（含 *_LIST.csv 与逐点 CSV）")
    ap.add_argument("--spot", default="S10", help="样品名（如 S10）")
    ap.add_argument("--out", default=None, help="输出目录；默认「批次目录上一级/spot_concordia」")
    ap.add_argument("--win", type=float, default=None, help="窗口宽度 (s)，默认 4")
    ap.add_argument("--step", type=float, default=None, help="窗口步长 (s)，默认 1")
    ap.add_argument("--nsigma", type=float, default=2.0, help="椭圆/误差棒的 σ 倍数，默认 2")
    ap.add_argument("--preset", default=None,
                    help="参考值口径（默认与流水线一致：constants.DEFAULT_REF_PRESET）")
    ap.add_argument("--list", action="store_true", help="只列出批次里的样品名后退出")
    ap.add_argument("--no-pdf", action="store_true", help="不写 PDF")
    ap.add_argument("--mean-ellipse", action="store_true",
                    help="域均值点也画椭圆（ρ 取该域窗口 ρ 的中位，属近似）")
    ap.add_argument("--window-ellipse", action="store_true",
                    help="逐窗口也画误差椭圆（默认不画：4 s 窗口的 207/235 "
                         "2σ 椭圆比坐标轴还宽，会把图糊死）")
    args = ap.parse_args(argv)

    kw = dict(data_dir=str(args.dir), plot=False, verbose=False)
    if args.win is not None:
        kw["win"] = args.win
    if args.step is not None:
        kw["step"] = args.step
    if args.preset:
        kw["ref_preset"] = args.preset
    cfg = BatchConfig(**kw)

    if args.list:
        spots, _s, _k = load_batch(cfg)
        print("样品名：", ", ".join(sorted({t.sample for t, _ in spots
                                        if t.role == ROLE_UNKNOWN})))
        return 0

    P = prepare(cfg, args.spot)
    out_dir = Path(args.out) if args.out else Path(args.dir).resolve().parent / "spot_concordia"
    stem = f"{P['label'].split()[0]}_{P['sample']}"

    fig1, p1 = draw_overview(P, out_dir / f"{stem}_单点总览.png", args.nsigma,
                             mean_ellipse=args.mean_ellipse,
                             window_ellipse=args.window_ellipse)
    fig2, p2 = draw_concordia_only(P, out_dir / f"{stem}_谐和分布.png", args.nsigma,
                                   mean_ellipse=args.mean_ellipse,
                                   window_ellipse=args.window_ellipse)
    p3 = save_depth_figure(P["prof"], P["segs"],
                           f"{P['label']}   深度剖面 [{P['tag']}]   "
                           f"窗口 {cfg.win:g}s / 步长 {cfg.step:g}s",
                           out_dir / f"{stem}_深度剖面.png", summ=P["summ"],
                           with_207=True, whole=P["ws"])

    if not args.no_pdf:
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(out_dir / f"{stem}_单点总览.pdf") as pdf:
            pdf.savefig(fig1, bbox_inches="tight")
            pdf.savefig(fig2, bbox_inches="tight")

    plt.close(fig1)
    plt.close(fig2)
    print_report(P)
    print()
    print("输出：")
    for p in ([p1, p2, p3] + ([] if args.no_pdf else [out_dir / f"{stem}_单点总览.pdf"])):
        print("   ", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
