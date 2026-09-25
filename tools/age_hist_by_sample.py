# -*- coding: utf-8 -*-
"""
tools/age_hist_by_sample.py —— 按样品名分画幅，画年龄分布柱状图（只读）
=======================================================================

    python tools/age_hist_by_sample.py "G:\\_gx\\_out_h2016"

    # 也可以逐个文件指定，或换个输出目录 / 分组规则
    python tools/age_hist_by_sample.py a.xlsx b.xlsx --out-dir D:\\图 --bins 12

它回答的是**中位数回答不了的问题**：中位数只说"这批年龄的中心在哪"，
看不出一个样品的测点是挤成一团，还是分两群（继承核 + 新生边）、
还是有一条长拖尾。柱状图看得见。

口径（与 `tests/check_example_batch.py` 完全一致）
--------------------------------------------------
1. **只统计 `类型 == 样品` 的测点**，不含标样。标样是 91500（约 1060 Ma）
   与 Ple（约 343 Ma）这种极端值，混进来会把分布整体拽偏、把 5–95% 撑开
   —— 仓库历史上真出过这个事故（中位 458.1 被抬到 460.9），所以这条是硬规矩。
2. 年龄列：**优先 `年龄206_238_QC校正`，没有才退回 `年龄206_238`**。
   前者是再乘了监控标样基体匹配系数的值，才是对外报出的"样品年龄"；
   桂北那类**没有监控标样**的批次（样品仓放不下 Ple）只有后者，用它即可。
3. 样品名分层：默认把末尾的测点序号去掉（`18-GP-1-7` -> `18-GP-1`），
   于是"一个取样点下打了很多激光点"会合成一个画幅。
   命名规则不一样时用 `--group-regex` 改，或传空字符串表示不分组（一个测点一个画幅）。

两张图，用途不同，都出
----------------------
A. `..._自适应`：每个画幅的横轴各自贴合本样品。
   用来**看单样品的分布形态**。样品之间的年龄跨度常常很大（桂北 19 个样品
   从 70 到 1838 Ma），统一坐标会把大多数样品压成一根竖线。
B. `..._统一坐标`：全部画幅共用一条横轴。
   用来**看样品之间的相对位置**——哪些样品是一回事、哪些明显是另一回事。

每个画幅里画了什么
------------------
- 柱：年龄直方图（各画幅各自分 `--bins` 箱）
- 红虚线：中位年龄（标了数字）
- 底部短竖线（rug）：**每一个测点一根**。直方图会骗人 —— 柱子高度取决于分箱，
  同一批数据换个箱宽就换个长相；rug 是原始数据，n 数得出来。

关于离群点（重要）
------------------
自适应模式下，**离群点不参与决定横轴范围**。理由是实测的：`18-GP-18` 有 24 个点，
其中 1 个 450 Ma，若按 min–max 定范围，主群（150–200 Ma）会被压成画幅左边一根
竖线，形态完全看不出来。判据用"相对 5–95% 跨度"，不用 Tukey 外栅栏 ——
每样品只有 24–29 个点，IQR 小到会让外栅栏误触发（`18-SH-21` 曾因此被移走 4 个点、
只画了剩下 23 个，而它本身根本没有离群点）。

被移出范围的点**在画幅右上角标数量**（`↑1 点 >263`）：是"看得见"，不是"消失"。
统一坐标模式不截断 —— 那条横轴本来就由调用方给定。
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from druid.core.constants import CJK_FONTS            # noqa: E402
from druid.io.report import AGE68_COLUMN, AGE68_QC_COLUMN   # noqa: E402

import matplotlib                                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402

RESULT_SHEET = "结果"
DEFAULT_GLOB = "*_U-Pb结果.xlsx"

BAR_FACE = "#9EC9E8"     # 柱子填充（浅蓝）
BAR_EDGE = "#2B6CA3"     # 柱子描边
MED_COLOR = "#C0392B"    # 中位线（红）
RUG_COLOR = "#34495E"    # 单点刻度（深灰蓝）
NOTE_COLOR = "#7F8C8D"   # 离群点标注


def collect_xlsx(paths):
    """把位置参数（目录或文件）展开成 xlsx 列表。目录按 DEFAULT_GLOB 收集。"""
    out = []
    for p in paths:
        p = pathlib.Path(p)
        if p.is_dir():
            out += sorted(p.glob(DEFAULT_GLOB))
        elif p.suffix.lower() in (".xlsx", ".xls"):
            out.append(p)
        else:
            raise SystemExit(f"既不存在的文件、也不是目录：{p}")
    if not out:
        raise SystemExit("没找到任何结果表（目录下应有 %s）" % DEFAULT_GLOB)
    return sorted(set(out))


def load_groups(xlsx, group_regex):
    """
    读一批结果表 -> [{batch, spot, age, column}, ...]。

    每张表读 `结果` sheet，只取 `类型 == 样品` 的行，按样品名（去序号）分组。
    """
    groups = []
    pat = re.compile(group_regex) if group_regex else None
    for p in xlsx:
        try:
            d = pd.read_excel(p, sheet_name=RESULT_SHEET)
        except ValueError as e:
            print("  跳过 %s（没有 `%s` sheet：%s）" % (p.name, RESULT_SHEET, e))
            continue
        if "类型" not in d.columns or "样品" not in d.columns:
            print("  跳过 %s（缺 `类型` / `样品` 列）" % p.name)
            continue

        # 列口径：优先 QC 校正列，没有才退回未校正列
        col = AGE68_QC_COLUMN if AGE68_QC_COLUMN in d.columns else AGE68_COLUMN
        if col not in d.columns:
            print("  跳过 %s（既没有 %s 也没有 %s）" % (p.name, AGE68_QC_COLUMN, AGE68_COLUMN))
            continue

        s = d[d["类型"] == "样品"].copy()
        if s.empty:
            print("  跳过 %s（没有样品测点）" % p.name)
            continue
        names = s["样品"].astype(str)
        s["_组"] = names.str.replace(pat, "", regex=True) if pat else names

        batch = (p.name[:-len("_U-Pb结果.xlsx")]
                 if p.name.endswith("_U-Pb结果.xlsx") else p.stem)
        s["_批"] = batch
        # 分组键带上批次：不同批次完全可能有同名样品，只按样品名分会把它们合并
        for (b, g), sub in s.groupby(["_批", "_组"]):
            a = pd.to_numeric(sub[col], errors="coerce").dropna().values.astype(float)
            if len(a):
                groups.append(dict(batch=b, spot=str(g), age=a, column=col))
    return groups


def sort_key(name):
    """18-GP-2 排在 18-GP-10 前面（按前缀、再按编号）；别的命名退回字典序。"""
    m = re.match(r"^(.*?)[-_ ]?(\d+)$", name)
    if not m:
        return (name, 0)
    return (m.group(1), int(m.group(2)))


def draw(groups, share_x, xlim, out_png, out_pdf, title, ncol, bins):
    """画一张多画幅图。share_x=True 时所有画幅共用横轴（由 xlim 给定）。"""
    plt.rcParams["font.sans-serif"] = CJK_FONTS
    plt.rcParams["axes.unicode_minus"] = False

    n = len(groups)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.15 * ncol, 2.55 * nrow),
                             sharex=share_x)
    axes = np.atleast_1d(axes).ravel()
    n_moved = 0

    for ax, r in zip(axes, groups):
        a = r["age"]
        med = float(np.median(a))
        lo, hi = float(a.min()), float(a.max())
        n_out_hi = n_out_lo = 0

        if share_x and xlim is not None:
            x0, x1 = xlim
        else:
            # 离群点不参与定范围（理由见模块 docstring）。判据用相对 5–95% 跨度。
            q05, q95 = np.percentile(a, [5, 95])
            spread = float(q95 - q05)
            if spread > 0 and hi > q95 + 3 * spread:
                hi = q95 + 1.5 * spread
                n_out_hi = int((a > hi).sum())
            if spread > 0 and lo < q05 - 3 * spread:
                lo = q05 - 1.5 * spread
                n_out_lo = int((a < lo).sum())
            span = hi - lo
            pad = max(span * 0.08, 1.0)
            x0, x1 = (lo - pad, hi + pad) if span > 0 else (lo - 20, hi + 20)

        inside = a[(a >= x0) & (a <= x1)]
        n_moved += len(a) - len(inside)

        cnt, edges = np.histogram(inside, bins=bins, range=(x0, x1))
        ymax = max(cnt.max(), 1)
        assert cnt.sum() == len(inside), "直方图计数与入选点数不符"
        ax.bar(edges[:-1], cnt, width=np.diff(edges), align="edge",
               color=BAR_FACE, edgecolor=BAR_EDGE, linewidth=0.8, zorder=2)

        ax.axvline(med, color=MED_COLOR, ls="--", lw=1.4, zorder=4)
        pos = (med - x0) / (x1 - x0) if x1 > x0 else 0.5
        ax.annotate("中位 %.1f" % med, xy=(med, ymax * 1.07),
                    ha=("left" if pos < 0.18 else "right" if pos > 0.82 else "center"),
                    va="bottom", fontsize=8, color=MED_COLOR, zorder=5)

        ax.vlines(inside, -0.14 * ymax, 0, color=RUG_COLOR, lw=0.8,
                  zorder=3, clip_on=False)

        if n_out_hi or n_out_lo:
            txt = "、".join(t for t in (
                ("↑%d 点 >%.0f" % (n_out_hi, x1)) if n_out_hi else "",
                ("↓%d 点 <%.0f" % (n_out_lo, x0)) if n_out_lo else "") if t)
            ax.text(0.985, 0.97, txt, transform=ax.transAxes, ha="right",
                    va="top", fontsize=7.5, color=NOTE_COLOR, zorder=6)

        ax.set_xlim(x0, x1)
        ax.set_ylim(-0.16 * ymax, ymax * 1.45)
        ax.set_yticks(np.arange(0, ymax + 1, max(1, int(np.ceil(ymax / 4)))))
        ax.tick_params(labelsize=8)
        ax.grid(axis="y", ls=":", lw=0.5, color="#BBBBBB", zorder=1)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_title("%s\nn=%d，%.0f–%.0f Ma"
                     % (r.get("label", r["spot"]), len(a), a.min(), a.max()),
                     fontsize=9, linespacing=1.5)

    for ax in axes[n:]:
        ax.set_visible(False)
    for ax in axes[max(0, n - ncol):n]:
        ax.set_xlabel("206Pb/238U 年龄 (Ma)", fontsize=9)

    fig.suptitle(title, fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)
    return n_moved


def main() -> int:
    ap = argparse.ArgumentParser(
        description="按样品名分画幅，画年龄分布柱状图（只读，不改任何结果表）")
    ap.add_argument("paths", nargs="+", help="结果表所在目录，或一个个 xlsx 文件")
    ap.add_argument("--out-dir", default=None, help="输出目录（默认与第一个输入同目录）")
    ap.add_argument("--prefix", default="按样品", help="输出文件名前缀")
    ap.add_argument("--bins", type=int, default=10, help="每个画幅的直方图箱数（默认 10）")
    ap.add_argument("--ncol", type=int, default=5, help="画幅列数（默认 5）")
    ap.add_argument("--group-regex", default=r"-\d+$",
                    help=r"去掉这个正则匹配的后缀得到样品名（默认 '-\d+$'）；传空串表示不分组")
    ap.add_argument("--title", default=None, help="总标题（默认自动写）")
    args = ap.parse_args()

    xlsx = collect_xlsx(args.paths)
    groups = load_groups(xlsx, args.group_regex)
    if not groups:
        raise SystemExit("没有读到任何样品测点。")
    groups.sort(key=lambda r: (sort_key(r["spot"]), r["batch"]))
    # 同名样品分散在多个批次时，画幅标题得带批次，否则两张画幅重名、分不清谁是谁
    seen = {}
    for r in groups:
        seen[r["spot"]] = seen.get(r["spot"], 0) + 1
    for r in groups:
        r["label"] = (r["spot"] if seen[r["spot"]] == 1
                      else "%s·%s" % (r["batch"], r["spot"]))

    out_dir = pathlib.Path(args.out_dir) if args.out_dir else xlsx[0].parent
    out_dir.mkdir(parents=True, exist_ok=True)

    allage = np.concatenate([r["age"] for r in groups])
    col = groups[0]["column"]
    print("结果表 %d 份，样品 %d 个，样品测点合计 %d" % (len(xlsx), len(groups), len(allage)))
    print("年龄列：%s　|　范围 %.1f–%.1f Ma，合并中位 %.1f Ma"
          % (col, allage.min(), allage.max(), float(np.median(allage))))
    print()
    for r in groups:
        a = r["age"]
        print("  %-18s n=%2d  中位 %8.2f  5–95%% %7.2f ~ %7.2f"
              % (r["label"], len(a), np.median(a),
                 np.percentile(a, 5), np.percentile(a, 95)))
    print()

    note = ("口径：仅样品测点（不含标样）· 206Pb/238U · %s · 每一根不是柱高而是单点（底部刻度）"
            % col)
    moved = draw(groups, False, None,
                 out_dir / (args.prefix + "_年龄分布柱状图_自适应.png"),
                 out_dir / (args.prefix + "_年龄分布柱状图_自适应.pdf"),
                 (args.title or "按样品名分画幅 · 年龄分布柱状图（各画幅横轴自适应）")
                 + "\n共 %d 个样品 / %d 个样品测点　|　%s" % (len(groups), len(allage), note),
                 args.ncol, args.bins)
    lo = max(0.0, np.floor((allage.min() - 30) / 50) * 50)
    hi = np.ceil((allage.max() + 50) / 50) * 50
    draw(groups, True, (lo, hi),
         out_dir / (args.prefix + "_年龄分布柱状图_统一坐标.png"),
         out_dir / (args.prefix + "_年龄分布柱状图_统一坐标.pdf"),
         (args.title or "按样品名分画幅 · 年龄分布柱状图（横轴统一 %g–%g Ma，便于横向比较）"
          % (lo, hi))
         + "\n共 %d 个样品 / %d 个样品测点　|　%s" % (len(groups), len(allage), note),
         args.ncol, args.bins)

    print("写出到 %s：" % out_dir)
    for f in sorted(out_dir.glob(args.prefix + "_年龄分布柱状图_*")):
        print("  %s  %d 字节" % (f.name, f.stat().st_size))
    if moved:
        print("\n注：自适应版有 %d 个测点因离群被移出横轴范围，各画幅右上角已标出数量"
              "（统一坐标版不截断）。" % moved)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
