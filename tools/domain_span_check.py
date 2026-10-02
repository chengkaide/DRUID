#!/usr/bin/env python
"""
tools/domain_span_check.py —— 跨域跨度异常的复核（只读，可复算）
================================================================

**它回答什么**

分域流程偶尔会解出一个"老得离谱"的域：同一个测点里既有 ~200 Ma 的域，
又有 ~1900 Ma 的域，跨度上 1000 Ma。这类点的共同特征是**老域总在最浅的几个
窗口（τ 小）**，因此嫌疑对象只有两个，而它们的处置完全相反：

    · **真继承核 / 捕获锆石** —— 值得写进结果，是地质信息；
    · **普通铅污染**（抛光粉、环氧树脂、镀碳层、表面吸附铅）——
      必须剔掉，否则会污染整个样品的年龄谱。

只靠 206Pb/238U 分不开这两者：普通铅同时抬高 206Pb 与 207Pb，
两个年龄**一起变老**，看起来仍然"协和"。分开它们要的是
**域级的 207Pb/235U 年龄**：

    域级协和度 = age(207Pb/235U) / age(206Pb/238U) × 100%

    · ≈100%  → 两点在协和线上 → 真继承核（有分量）
    · ≫100%  → 有过量 207Pb   → 普通铅污染（该剔）

**为什么不能只看 `分析` 表的 `协和度_pct`**

那一列是**整点**的：整点把老域、过渡带、新域混在一起平均，
一个 480 Ma 的真年龄域配一个 2270 Ma 的污染域，整点协和度会被拉平到
看不出所以然。必须**按域**算 —— 而 `汇总库` 的 `窗口` 表只入了 R68，
没有 R76，所以这里回到原始 CSV 重跑逐窗口。

**为什么重跑而不是从库里推算**

207Pb/235U = (207Pb/206Pb) × (206Pb/238U) × 238U/235U，
其中 206Pb/238U 必须用**分馏校正后**的值（R68 × F(τ)）。
F(τ) 是标样逐深度的归一化因子，只存在于流水线内部，库里没有。
所以本脚本复用 `druid.workflow` 的前置步骤，口径与出结果表时**完全一致**。

**口径**

    · 只读 `汇总库/DRUID汇总.db`；原始数据也只读。
    · 参考值预设与汇总库一致（`批次清单.csv` 的 `参考值口径` 列）。
    · 跨度只算 `标记 = 'age domain'` 的域；`--include-transitional`
      可把 `mixed/过渡带` 也算进来（那会把口径放宽，见 `改进清单.md`）。

用法
----
    python tools/domain_span_check.py            # 打印复核表（默认）
    python tools/domain_span_check.py --span 300 # 换跨度门槛

综合图与在线展示页在 `tools/gen_docs_extreme_spans.py`（它复用本模块的
`recompute_batch` 与 `read_spot_list`）—— 这里只打印表，不写任何文件。

增删任何判据阈值前请先读 `AGENTS.md`：这里的数是**判读输出**，
不是训练/校准参数，改阈值不会改变年龄本身。
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from druid.core.constants import ROLE_UNKNOWN                    # noqa: E402
from druid.core.geochronology import age75, r75_from             # noqa: E402
from druid.core.references import std_alias, std_ref             # noqa: E402
from druid.core.statistics import weighted_mean                  # noqa: E402
from druid.depth.fractionation import bracket_F                  # noqa: E402
from druid.workflow import (BatchConfig, analyse_depth_spot,     # noqa: E402
                            build_bracketing, calibrate_primary,
                            compute_bulk_ratios, estimate_external,
                            load_batch)

DB_DEFAULT = ROOT / "汇总库" / "DRUID汇总.db"
LIST_DEFAULT = ROOT / "汇总库" / "批次清单.csv"


# ─────────────────────────────────────────────────────────────────────────────
# 一、挑出候选测点
# ─────────────────────────────────────────────────────────────────────────────
def read_spot_list(db: Path, span: float, include_transitional: bool):
    """
    返回跨度超过 `span` 的测点，按跨度降序。

    每条记录带上所属批次、区组、样品名与文件夹 —— 后两者是重跑时找
    原始数据用的。批次清单里的一层目录名（形如 `区组\\样品__批次号`）
    是**用反斜杠拼的**，在非 Windows 上 Path 会把它当成一个整体名字，
    所以这里显式拆开。
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    where = "" if include_transitional else "where 标记 = 'age domain'"
    rows = con.execute(f"""
        select 批次, 区组, 序号, 样品, count(*) n,
               min(年龄_Ma) lo, max(年龄_Ma) hi,
               max(年龄_Ma) - min(年龄_Ma) span
        from 域 {where}
        group by 批次, 序号 having n >= 2 and span > ?
        order by span desc
    """, (span,)).fetchall()
    out = [dict(r) for r in rows]

    # 文件夹与参考值口径都从批次清单取：前者用来找原始数据，后者必须与
    # 出结果表时**同一套参考值**，否则这里算出来的域级年龄与库里对不上。
    # 早期版本只取了文件夹，`参考值口径` 一直是 `.get(..., 默认)` ——
    # 25 批恰好都是 horstwood2016，所以症状是静默的。
    folder, preset = {}, {}
    if LIST_DEFAULT.exists():
        import csv
        with open(LIST_DEFAULT, encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                folder[r["批次"]] = r["文件夹"]
                preset[r["批次"]] = (r.get("参考值口径") or "").strip()
    for r in out:
        rel = folder.get(r["批次"], r["批次"])
        r["文件夹"] = rel.replace("\\", "/")
        r["参考值口径"] = preset.get(r["批次"]) or "horstwood2016"
    con.close()
    return out


def group_by_batch(spots):
    """{批次: [测点…]}，保持原来的跨度降序。"""
    out: dict[str, list] = {}
    for s in spots:
        out.setdefault(s["批次"], []).append(s)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 二、域级量：207Pb/235U 年龄与协和度
# ─────────────────────────────────────────────────────────────────────────────
def _sigma_75(r75, s75):
    """
    207Pb/235U 年龄的 1σ。与 fractionation.age76 同一套做法：
    比值方程无法解析反解，用两端求差代替微分（在 σ 很小时等价）。
    """
    a_hi = age75(r75 + s75)
    a_lo = age75(r75 - s75)
    return np.abs(a_hi - a_lo) / 2.0


def domain_audit(prof, summ, brack, ref68: float) -> list:
    """
    把逐窗口剖面按域汇总，补上**域级 207Pb/235U 年龄**与协和度。

    206Pb/238U 侧直接用 `summ` 里已经算好的 `age_Ma`（反比方差加权），
    不另起一套口径；207Pb/235U 侧用同样的权重（1/σ²）加权。
    """
    tau = prof["tau"].to_numpy()
    r68c = prof["R68"].to_numpy() * bracket_F(tau, brack, "R68", ref68)
    r76 = prof["R76"].to_numpy()
    # 207Pb/235U 的 1σ：由 R68、R76 的相对误差按平方和传播。
    # 这里**不加 rho 交叉项** —— 它只用于给加权平均定权重，
    # 而权重对中位数级别的差异不敏感（域内窗口数只有 3~20 个）。
    s75 = r75_from(r68c, r76) * np.hypot(
        prof["s68"].to_numpy() / prof["R68"].to_numpy(),
        prof["s76"].to_numpy() / r76)
    a75 = age75(r75_from(r68c, r76))
    s_a75 = _sigma_75(r75_from(r68c, r76), s75)

    out = []
    for _, d in summ.iterrows():
        i0, i1 = int(d["i0"]), int(d["i1"])
        mu75, se75, mswd75, k = weighted_mean(a75[i0:i1], s_a75[i0:i1])
        age68_dom = float(d["age_Ma"])
        out.append(dict(
            domain=str(d["domain"]),
            flag=str(d["flag"]),
            tau0=float(d["tau0"]), tau1=float(d["tau1"]),
            n_win=int(d["n_win"]),
            age68=age68_dom,
            s2_68=2.0 * float(d["se_1sig"]),
            age75=float(mu75),
            s2_75=2.0 * float(se75),
            mswd=float(d["mswd"]),
            mswd_prob=float(d["mswd_prob"]),
            mswd75=float(mswd75),
            concordance=(float(mu75) / age68_dom * 100.0) if age68_dom else float("nan"),
            th_u=float(d["ThU"]),
            u_cps=float(d["U_cps"]),
            k75=int(k),
        ))
    return out


def recompute_batch(data_dir: Path, ref_preset: str, want_idx: set):
    """
    重跑一个批次的**前置步骤**（装载 → 夹逼 → 归一化因子 → σext），
    然后只对 `want_idx`（1 起的序号）里的点做深度剖面。

    返回 {序号: {"domains": […], "f206": …, "时段": …, "prof": …}}
    """
    cfg = BatchConfig(data_dir=str(data_dir), ref_preset=ref_preset,
                      plot=False, verbose=False)
    spots, _seq, _skipped = load_batch(cfg)
    brack_for = build_bracketing(spots, cfg)
    ref68, ref76 = std_ref(std_alias(cfg.primary) or cfg.primary, cfg.ref_preset)
    _variants, (b68, b76) = compute_bulk_ratios(spots, brack_for, ref68, cfg)
    pidx, F68, F76, _cal = calibrate_primary(spots, b68, b76, ref68, ref76, cfg)
    sd68, sd76 = estimate_external(spots, pidx, F68, F76, b68, b76, cfg)

    out = {}
    for k, (tr, r) in enumerate(spots):
        if tr.role != ROLE_UNKNOWN:
            continue
        idx = int(tr.idx) + 1
        if idx not in want_idx:
            continue
        got = analyse_depth_spot(tr, r, brack_for(k), ref68, sd68, cfg)
        if got is None:
            out[idx] = {"domains": [], "note": "窗口为空或无可用夹逼标样"}
            continue
        prof, _segs, summ, tag = got
        out[idx] = dict(
            domains=domain_audit(prof, summ, brack_for(k), ref68),
            tag=tag,
            f206=float(r["f206"]),
            span_s=float(tr.t1 - tr.t0),
            n_win=len(prof),
            sd68=float(sd68), sd76=float(sd76),
        )
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 三、控制台表
# ─────────────────────────────────────────────────────────────────────────────
def _fmt(v, spec="%.1f", dash="--"):
    return dash if v is None or (isinstance(v, float) and not np.isfinite(v)) else spec % v


def print_table(spots, audited):
    """把复核结果打成一张可肉眼扫的表。"""
    print()
    print("=" * 112)
    print("跨域跨度异常的逐点复核（只算正式年龄域；协和度 = 域级 207/235 ÷ 206/238）")
    print("=" * 112)
    print("%-3s %-16s %6s %8s  %-24s" % ("#", "标签", "跨度", "τ 区间", "各域（年龄 / τ / n / MSWD / 协和度 / Th/U）"))
    print("-" * 112)
    for i, s in enumerate(spots, 1):
        aud = audited.get(s["批次"], {}).get(s["序号"])
        if aud is None:
            print("%-3d %-16s %8.1f   （该点未参与复核）" % (i, f"{s['样品']}#{s['序号']}", s["span"]))
            continue
        ds = [d for d in aud["domains"] if d["flag"] == "age domain"]
        if not ds:
            print("%-3d %-16s %8.1f   %s" % (i, f"{s['样品']}#{s['序号']}", s["span"],
                                             aud.get("note", "无正式年龄域")))
            continue
        head = "%-3d %-16s %8.1f  %s" % (
            i, f"{s['样品']}#{s['序号']}", s["span"],
            "%s~%s" % (_fmt(min(d["tau0"] for d in ds), "%.2f"),
                       _fmt(max(d["tau1"] for d in ds), "%.2f")))
        print(head)
        for d in ds:
            print("      D%-2s τ=%4.2f-%-4.2f n=%-2d MSWD=%5.2f  %8s±%-7s Ma"
                  "  协和度=%7s%%  Th/U=%.3f"
                  % (d["domain"], d["tau0"], d["tau1"], d["n_win"], d["mswd"],
                     _fmt(d["age68"], "%8.1f"), _fmt(d["s2_68"], "%.1f"),
                     _fmt(d["concordance"], "%7.2f"), d["th_u"]))
    print("-" * 112)


# ─────────────────────────────────────────────────────────────────────────────
# 四、命令行
# ─────────────────────────────────────────────────────────────────────────────
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="复核跨域跨度过大的测点：算域级 207/235 年龄与协和度，"
                    "分辨真继承核与普通铅污染。只读，不改任何数据。")
    ap.add_argument("--db", default=str(DB_DEFAULT), help="汇总库路径")
    ap.add_argument("--span", type=float, default=500.0,
                    help="跨度门槛 (Ma)，默认 500")
    ap.add_argument("--include-transitional", action="store_true",
                    help="把 mixed/过渡带也算进跨度（默认只算正式年龄域）")
    args = ap.parse_args(argv)

    db = Path(args.db)
    if not db.exists():
        print("汇总库不存在：%s\n先跑 `python 汇总库/建库.py 全部`。" % db, file=sys.stderr)
        return 2

    spots = read_spot_list(db, args.span, args.include_transitional)
    if not spots:
        print("没有跨度 > %g Ma 的测点。" % args.span)
        return 0

    audited: dict[str, dict] = {}
    for batch, group in group_by_batch(spots).items():
        folder = group[0]["文件夹"]
        data_dir = ROOT / "汇总库" / folder / "原始数据"
        preset = group[0].get("参考值口径", "horstwood2016")
        print("[复核] %s  (%d 个候选点)  %s" % (batch, len(group), data_dir.name),
              file=sys.stderr)
        if not data_dir.exists():
            print("   !! 原始数据不存在，跳过：%s" % data_dir, file=sys.stderr)
            continue
        want = {s["序号"] for s in group}
        try:
            audited[batch] = recompute_batch(data_dir, preset, want)
        except Exception as exc:                      # noqa: BLE001
            print("   !! 重跑失败：%s" % exc, file=sys.stderr)

    print_table(spots, audited)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
