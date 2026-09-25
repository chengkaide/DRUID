"""
druid.depth.figures —— 深度剖面图
===============================

三栏共享横轴的一张图，从上到下：

    ┌──────────────────────────────────────────────────────┐
    │  年龄剖面（206Pb/238U ±2σ）                            │
    │  ＋ 207Pb/206Pb 叠加（可选）                           │
    │  ＋ 年龄域底色 ＋ 各域加权平均年龄的红线标注            │
    ├──────────────────────────────────────────────────────┤
    │  238U 信号强度 —— 激光打到不同环带的直接指示            │
    ├──────────────────────────────────────────────────────┤
    │  Th/U —— 判别核/边的关键地球化学指标                    │
    └──────────────────────────────────────────────────────┘
                剥蚀时间 (s)  ——  由浅到深

为什么要把信号强度和 Th/U 放在一起看
------------------------------------
年龄曲线自己不会说话。要判断"这里是不是真的有一个边"，
必须有独立的证据：
    · U 含量突变 → 激光进入了一个成分不同的环带
    · Th/U 突变 → 岩浆成因（Th/U > 0.5）↔ 变质成因（Th/U < 0.1）的分野
三者在同一时刻改变，**才算站得住**。
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..core.constants import CJK_FONTS, DOMAIN_SPAN_COLORS


def _setup_style():
    """
    统一配置 matplotlib：Agg 后端 + 中文字体 + 负号显示。

    为什么要用 Agg
    --------------
    Agg 是纯计算后端，不需要图形界面。批处理时会在没有显示器的
    环境里跑（远程/后台），用默认后端可能在第一个 show() 时就卡住。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = CJK_FONTS   # 中文字体候选
    plt.rcParams["axes.unicode_minus"] = False    # 否则负号会渲染成方块
    return plt


def make_depth_figure(prof: pd.DataFrame,
                      segs: Sequence[Tuple[int, int]],
                      title: str,
                      summ: Optional[pd.DataFrame] = None,
                      with_207: bool = True):
    """
    绘制一个测点的深度剖面图。

    参数
    ----
    prof      : window_profile 结果，需含 age68 / s_age68 / U_cps / ThU / t_mid
    segs      : 年龄域下标范围列表
    title     : 图标题
    summ      : summarize_segments 的输出，用于画每个域的平均年龄红线
    with_207  : 是否叠加 207Pb/206Pb（年轻锆石叠加只会增加混乱，可关掉）

    返回
    ----
    matplotlib Figure（**不负责 close**，由调用方决定何时释放）
    """
    plt = _setup_style()

    # 三栏共享横轴；高度比 2.6:1:1 让年龄曲线占主要视觉空间
    fig, ax = plt.subplots(
        3, 1, figsize=(9, 8.6), sharex=True,
        gridspec_kw=dict(height_ratios=[2.6, 1, 1], hspace=0.07))

    t = prof["t_mid"].to_numpy()
    # 只有确实存在有限值才叠加 207 曲线
    has76 = with_207 and ("age76" in prof.columns) \
        and bool(np.isfinite(pd.to_numeric(prof["age76"], errors="coerce")).any())

    # ── 年龄域底色：三栏全部打上，便于跨栏对齐观察 ──
    if len(segs) > 1:
        for k in (0, 1, 2):
            for j, (lo, hi) in enumerate(segs):
                x0 = prof["t_mid"].iloc[lo]
                x1 = prof["t_mid"].iloc[hi - 1]
                ax[k].axvspan(x0, x1,
                              color=DOMAIN_SPAN_COLORS[j % len(DOMAIN_SPAN_COLORS)],
                              alpha=0.75, zorder=0)

    # ── 第一栏：年龄 ──
    # errorbar 的 yerr 用的是 **2σ**，理由：出版物惯例是给 95% 置信区间，
    # 而且 1σ 的误差棒在这张图上小到看不见，容易让人误以为精度极高。
    ax[0].errorbar(t, prof["age68"], yerr=prof["s_age68"].to_numpy() * 2,
                   fmt="o-", ms=3.6, lw=1.2, color="#185FA5", ecolor="#9FC3E8",
                   capsize=2, zorder=3, label="206Pb/238U (±2σ)")
    if has76:
        y76 = pd.to_numeric(prof["age76"], errors="coerce").to_numpy()
        e76 = pd.to_numeric(prof.get("s_age76"), errors="coerce").to_numpy() \
            if "s_age76" in prof.columns else np.zeros_like(y76)
        ax[0].errorbar(t, y76, yerr=np.nan_to_num(e76) * 2,
                       fmt="s--", ms=3, lw=1.0, color="#993C1D", ecolor="#F0997B",
                       capsize=2, alpha=0.75, zorder=2, label="207Pb/206Pb (±2σ)")

    # ── 每个年龄域的加权平均年龄：红线 + 文字标注 ──
    if summ is not None and len(summ):
        # 文字抬离红线的距离用**点**而不是数据单位 —— 同一个偏移量画在
        # 450 Ma 的点和 3000 Ma 的点上，视觉效果差着十几倍。
        # 4 pt 在 fontsize=9 下约 0.6 个字高：够把文字抬离红线，
        # 又不至于顶到相邻域的标注上去。
        # （原先直接 va="bottom" 贴在 age_Ma 上，实测红线从字底穿过去，
        #   几个域的标注看着像被划掉。）
        t_lo = float(prof["t_mid"].iloc[0])
        t_hi = float(prof["t_mid"].iloc[-1])
        span = max(t_hi - t_lo, 1e-9)
        for _, s in summ.iterrows():
            if str(s["flag"]).startswith("mixed"):
                continue                       # 过渡带不画，它没有地质意义
            lo, hi = int(s["i0"]), int(s["i1"]) - 1
            x0, x1 = prof["t_mid"].iloc[lo], prof["t_mid"].iloc[hi]
            ax[0].hlines(s["age_Ma"], x0, x1, color="#A32D2D", lw=1.8, zorder=4)
            txt = f"{s['domain']}: {s['age_Ma']:.0f}±{2 * s['se_1sig']:.0f} Ma"
            xc = 0.5 * (x0 + x1)
            # 贴着左右边界的域（最常见的是最后一个域）若还居中排版，
            # 文字会越过坐标轴画到轴外去。改成贴住该域的内侧端点、向里展开。
            if xc > t_lo + 0.80 * span:
                tx, ha = x1, "right"
            elif xc < t_lo + 0.20 * span:
                tx, ha = x0, "left"
            else:
                tx, ha = xc, "center"
            ax[0].annotate(txt, xy=(tx, s["age_Ma"]), xytext=(0, 4.0),
                           textcoords="offset points", fontsize=9,
                           color="#A32D2D", va="bottom", ha=ha, zorder=5)

    ax[0].set_ylabel("年龄 (Ma)", fontsize=11)
    ax[0].set_title(title, fontsize=12)
    ax[0].legend(fontsize=9, loc="best", framealpha=0.9)

    # ── 第二栏：238U 信号 ──
    ax[1].plot(t, prof["U_cps"], color="#534AB7", lw=1.4)
    ax[1].set_ylabel("238U (cps)", fontsize=11)

    # ── 第三栏：Th/U ──
    ax[2].plot(t, prof["ThU"], color="#0F6E56", lw=1.4)
    ax[2].set_ylabel("Th/U", fontsize=11)
    ax[2].set_xlabel("剥蚀时间 (s)  ——  由浅到深对应坑深", fontsize=11)

    for a_ in ax:
        a_.grid(alpha=0.25, lw=0.5)
        a_.tick_params(labelsize=9)
    return fig


def save_depth_figure(prof, segs, title, out_path, summ=None,
                      dpi: int = 140, with_207: bool = True,
                      pdf=None) -> Path:
    """
    画完直接存盘并释放内存。

    为什么必须显式 plt.close(fig)
    ----------------------------
    一个批次几十个测点，若不关 Figure，matplotlib 会把它们全留在内存里，
    几百张之后会触发 "More than 20 figures opened" 警告甚至 OOM。

    参数
    ----
    pdf : matplotlib.backends.backend_pdf.PdfPages 对象（可选）
          传入则同时把该图追加进多页 PDF。**推荐在循环里一边画一边写入 PDF
          并立即关闭**，避免几十个 Figure 同时驻留内存。
    """
    fig = make_depth_figure(prof, segs, title, summ, with_207=with_207)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    if pdf is not None:                      # 立刻写进 PDF，之后就可以释放
        pdf.savefig(fig, bbox_inches="tight")
    # 立即关闭句柄，防止批处理几十上百张图时内存堆积
    import matplotlib.pyplot as plt
    plt.close(fig)
    return out_path
