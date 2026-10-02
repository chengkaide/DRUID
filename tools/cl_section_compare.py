"""
tools/cl_section_compare.py —— 把 CL 图像（空间）与深度剖面（时间）对上
=====================================================================

为什么需要这个工具
------------------
CL 图像给的是**空间**信息：颗粒内的核 / 幔 / 增生边、环带宽度，单位 μm。
LA-ICP-MS 深度剖面给的是**时间**信息：逐窗口的 τ 与年龄，τ 是归一化进度。
两者要对话，中间只有一座桥：

    z(μm，从积分起点算) = τ × 剥蚀总时长(s) × 剥蚀速率(μm/s)

这是一步除法，**完全可以定量**；不定量的是"桥墩"（速率）要从外面搬进来。

⚠ 三个必须先知道的坑
--------------------
1. **τ 的零点不是坑口。**
   τ = (t_mid − t0) / (t1 − t0)，而 t0 / t1 是**裁掉 trim 秒之后**的积分
   窗口端点（见 `druid/depth/windows.py` 的 `window_profile`）。默认
   trim = 1.5 s ⇒ **坑口那 1.5 s 的信息根本不在剖面里**。
   要换算成"从坑口起算的深度"，得把 trim 加回去：

       z_自坑口 = (trim + τ × 剥蚀总时长) × 速率

   直接把 τ 当成"从坑口算的深度比例"会系统性偏掉开头那一段。

2. **剥蚀速率无法从数据反算，只能外部给。**
   工具链里既没有测点坐标，也没有剥蚀坑深度（原始 CSV 只有计数与时间）。
   速率要么做完后实测（SEM / 白光干涉量坑深，除以积分时长），
   要么用本实验室的标称值。**本工具不接受"猜"的默认值** —— 不传
   `--rate` 就只打印 τ，不打印 μm。

3. **坑径必须远小于结构宽度，这个换算才成立。**
   LA 坑直径是几十 μm 量级；只有坑径 ≪ 环带宽度时，剥蚀才近似"沿一条
   线"穿过结构。若坑本身跨过了两个结构，这个测点的剖面天生就是混合的。
   此时 CL 的价值反而最大：**它能在做之前就把这种点挑出来**，
   而不是等深度剖面解出一个说不清的"中间年龄"。

CL 侧要提供什么
----------------
不用改任何数据，只要自己填一张表（`template` 子命令生成空模板）：

    测点, 颗粒编号, 颗粒内位置, CL结构数, 界面深度_um, 备注

· 颗粒编号：同一颗锆石上的多个测点填同一个号 —— 这是最强的自检
  （同一颗粒的深度剖面本该一致，不一致就该回头看 CL）。
· CL结构数：你在图上数出几个结构（核 / 幔 / 边）。
· 界面深度_um：如果量到了界面位置就填，工具会与 τ 换算出的深度并列。

三个子命令
----------
    # ① 打印工具侧的域结构 + τ / 深度换算（拿到 CL 图后照着核对的清单）
    python tools/cl_section_compare.py list --batch 示例批次 --rate 0.35

    # ② 生成 CL 判读空模板（把测点行铺好，你只填 CL 那几列）
    python tools/cl_section_compare.py template --batch 示例批次 --out CL判读.csv

    # ③ 读回填好的表，算"工具 vs CL"的一致 / 漏检 / 过拟合
    python tools/cl_section_compare.py compare --cl CL判读.csv

只读：只查询 `汇总库/DRUID汇总.db`，不写任何数据、不改任何输出口径。
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from pathlib import Path

DB_DEFAULT = Path(__file__).resolve().parent.parent / "汇总库" / "DRUID汇总.db"

# 与 cfg.trim 的默认值保持一致（见 druid/cli/reduce_batch.py 的 --trim）
TRIM_DEFAULT = 1.5

CL_COLUMNS = ["测点", "颗粒编号", "颗粒内位置", "CL结构数", "界面深度_um", "备注"]


# ─────────────────────────────────────────────────────────────────────────────
# 读库
# ─────────────────────────────────────────────────────────────────────────────
def con(db: Path) -> sqlite3.Connection:
    """以**只读**方式打开汇总库（mode=ro），避免误写。"""
    if not db.exists():
        raise SystemExit(f"找不到汇总库：{db}\n先跑 python 汇总库/建库.py 全部")
    return sqlite3.connect(f"file:{db}?mode=ro", uri=True)


def lo_hi(rng: str) -> tuple[float, float]:
    """
    把 `域.tau` 的字符串区间（形如 '0.05-0.46'）拆成 (起点, 终点)。

    ⚠ 这里有个真踩过的坑：`域` 表的 `tau` 是**字符串**，`窗口` 表的 `Tau`
      才是数字。直接 float() 会抛异常；外层若套裸 `except: pass`，会静默
      跳过全部行（症状是结果为空但脚本不报错）。
    """
    a, b = str(rng).split("-", 1)
    return float(a), float(b)


def spot_rows(c: sqlite3.Connection, batch: str | None) -> list[dict]:
    """
    把「每个样品测点的域结构 + 剥蚀窗口绝对时刻」拼成一条记录。

    为什么要把两种信息拼在一起：τ 是**相对**进度，单看 τ 不知道它落在
    第几秒；而 `分析.剥蚀窗口_s` 存着绝对时刻区间（形如 '33.4-72.1'），
    两者一拼就能把 τ 换算回绝对时间，进而换算成 μm 深度。
    """
    # 批次过滤统一在这里拼，避免两处 where 各写一遍写岔
    w_sql, w_args = ("where 批次 = ?", (batch,)) if batch else ("", ())
    w_type = ("and" if w_sql else "where") + " 类型='样品'"

    # 1) 剥蚀窗口：样品测点的绝对时刻区间与样品名
    win: dict[tuple[str, int], dict] = {}
    for b, g, no, f, sp, w in c.execute(
        f"select 批次,区组,序号,文件,样品,剥蚀窗口_s from 分析 {w_sql} {w_type}",
        w_args
    ):
        if not w or "-" not in str(w):
            continue
        t0, t1 = (float(x) for x in str(w).split("-", 1))
        win[(b, no)] = dict(区组=g, 文件=f, 样品=sp, t0=t0, t1=t1, span=t1 - t0)

    # 2) 域结构：只取真正定年的域（'age domain'），过渡带单独记
    doms: dict[tuple[str, int], list[dict]] = {}
    mixed: dict[tuple[str, int], int] = {}
    for b, no, name, mark, tau, age, s2, ms, nw in c.execute(
        f"select 批次,序号,域,标记,tau,年龄_Ma,s2_Ma,MSWD,n_win from 域 {w_sql}",
        w_args
    ):
        if (b, no) not in win:
            continue
        if mark == "age domain":
            a, z = lo_hi(tau)
            doms.setdefault((b, no), []).append(
                dict(域=name, tau0=a, tau1=z, 年龄=age, s2=s2, MSWD=ms, n_win=nw))
        else:
            mixed[(b, no)] = mixed.get((b, no), 0) + 1

    out = []
    for key, w in sorted(win.items()):
        ds = sorted(doms.get(key, []), key=lambda d: d["tau0"])
        out.append(dict(key=key, 域=ds, 过渡带数=mixed.get(key, 0), **w))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 换算
# ─────────────────────────────────────────────────────────────────────────────
def to_depth(tau: float, span: float, trim: float, rate: float | None
             ) -> tuple[float, float | None, float | None]:
    """
    τ → (绝对时刻 s, 自积分起点的深度 μm, 自坑口的深度 μm)。

    · 绝对时刻  : t = t0 + τ × span（t0 是**积分起点**，不含 trim）
    · 自积分起点: z = τ × span × rate
    · 自坑口    : z = (trim + τ × span) × rate   ← 多算了坑口那 trim 秒

    rate 为 None 时后两个返回 None —— 宁可空着，也不给一个猜的深度。
    """
    t_rel = tau * span
    if rate is None:
        return t_rel, None, None
    return t_rel, t_rel * rate, (trim + t_rel) * rate


def interfaces(ds: list[dict]) -> list[dict]:
    """
    相邻两个年龄域之间的界面位置。

    界面取"前一域的终点"与"后一域的起点"的**中点** —— 工具在这个间隙里
    没有给出确定年龄（那段被算作过渡带）。用中点只是给它一个代表位置，
    不代表界面精确落在这里。
    """
    out = []
    for i in range(len(ds) - 1):
        out.append(dict(after=ds[i]["域"], before=ds[i + 1]["域"],
                        tau=0.5 * (ds[i]["tau1"] + ds[i + 1]["tau0"])))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 子命令 ①：list —— 打印工具侧的清单
# ─────────────────────────────────────────────────────────────────────────────
def cmd_list(c: sqlite3.Connection, a: argparse.Namespace) -> int:
    rows = spot_rows(c, a.batch)
    if not rows:
        print("没有查到样品测点。检查 --batch 拼写，或先建库。")
        return 1

    rate = a.rate
    print("=" * 74)
    print("工具侧：深度剖面解出的年龄域（这就是拿到 CL 图后要对照的清单）")
    print("=" * 74)
    if rate is None:
        print("⚠ 未给 --rate，只打印 τ 与绝对时刻；μm 深度留空。")
        print("  剥蚀速率无法从数据反算，必须自己量（见文件头「坑 2」）。")
    else:
        print(f"剥蚀速率 {rate} μm/s，trim {a.trim} s"
              f"（τ=0 在积分起点，不是坑口；「自坑口」一列已把 trim 加回）")
    print()

    n_multi = 0
    for r in rows:
        ds = r["域"]
        if not ds:
            continue
        if len(ds) >= 2:
            n_multi += 1
        tag = f"{len(ds)} 域" if len(ds) >= 2 else "均一"
        print(f"── {r['文件']:<22} {r['样品']:<14} {tag}"
              f"   剥蚀 {r['t0']:.1f}–{r['t1']:.1f} s（{r['span']:.1f} s）"
              f"   过渡带 {r['过渡带数']} 段")
        for d in ds:
            t_rel0, z0, _ = to_depth(d["tau0"], r["span"], a.trim, rate)
            t_rel1, z1, _ = to_depth(d["tau1"], r["span"], a.trim, rate)
            depth = "" if z0 is None else f"  深度 {z0:6.1f}–{z1:6.1f} μm（自起点）"
            print(f"     {d['域']:<3} τ {d['tau0']:.2f}–{d['tau1']:.2f}"
                  f"   年龄 {d['年龄']:7.1f} ± {d['s2']:4.1f} Ma(2σ)"
                  f"   MSWD {d['MSWD']:5.2f}  n_win {d['n_win']:<3}"
                  f"  时刻 {t_rel0 + r['t0']:5.1f}–{t_rel1 + r['t0']:5.1f} s{depth}")
        for it in interfaces(ds):
            _, _, zg = to_depth(it["tau"], r["span"], a.trim, rate)
            tg = r["t0"] + it["tau"] * r["span"]
            ztxt = "" if zg is None else f"   自坑口 {zg:6.1f} μm"
            print(f"     ↑ 界面 {it['after']}→{it['before']}  τ≈{it['tau']:.2f}"
                  f"   时刻≈{tg:.1f} s{ztxt}")
        print()

    # 界面位置的分布 —— 顺带把「τ 到不了 1」这件事说清楚，免得误读成盲区
    allif = [i["tau"] for r in rows for i in interfaces(r["域"])]
    if allif:
        far = sum(1 for t in allif if t > 0.8)
        print(f"统计：{len(rows)} 个样品测点，{n_multi} 个解出 ≥2 个域"
              f"（{100 * n_multi / len(rows):.0f}%）；界面 {len(allif)} 个，"
              f"其中 τ>0.8 的 {far} 个。")
        print("口径提醒：**τ 的可达范围不是 0~1**，而是")
        print("    [win/(2·总时长), 1 − win/(2·总时长)]")
        print("  因为 τ 取的是滑窗**中点**，窗口中心永远到不了剥蚀段的两端。")
        print("  本批 win=4 s、总时长≈39 s ⇒ 端点 ≈ [0.05, 0.95]（全库实测")
        print("  0.0459~0.9494，与这个预测吻合）。这是几何必然，")
        print("  不要把「τ 到不了 1」误读成「分域在后段失灵」。")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# 子命令 ②：template —— 生成 CL 判读空模板
# ─────────────────────────────────────────────────────────────────────────────
def cmd_template(c: sqlite3.Connection, a: argparse.Namespace) -> int:
    rows = spot_rows(c, a.batch)
    if not rows:
        print("没有查到样品测点。")
        return 1
    out = Path(a.out)
    with out.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(CL_COLUMNS)
        for r in rows:
            ds = r["域"]
            if not ds:
                continue
            # 「工具说几个域」放在备注里，填表时能直接对照
            note = f"工具解出 {len(ds)} 个域：" + " / ".join(
                f"{d['域']} {d['年龄']:.0f}Ma" for d in ds)
            w.writerow([r["文件"], "", "", "", "", note])
    print(f"已写出模板：{out}")
    print(f"  {len(rows)} 行测点。你只需填「颗粒编号 / 颗粒内位置 / CL结构数 /"
          f" 界面深度_um」四列，然后跑：")
    print(f"  python tools/cl_section_compare.py compare --cl \"{out}\"")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# 子命令 ③：compare —— 工具 vs CL 的一致性
# ─────────────────────────────────────────────────────────────────────────────
def _as_int(v) -> int | None:
    try:
        s = str(v).strip()
        return int(float(s)) if s else None
    except (TypeError, ValueError):
        return None


def cmd_compare(c: sqlite3.Connection, a: argparse.Namespace) -> int:
    path = Path(a.cl)
    if not path.exists():
        print(f"找不到 CL 判读表：{path}")
        return 1
    with path.open(encoding="utf-8-sig") as f:
        cl = {r["测点"].strip(): r for r in csv.DictReader(f) if r.get("测点", "").strip()}

    # 按测点文件反查工具的解（文件名即可对上 —— 模板里铺的就是它）
    tool: dict[str, list[dict]] = {r["文件"]: r["域"] for r in spot_rows(c, None)}

    both_multi = both_one = 0
    miss: list[tuple[str, int, int]] = []
    over: list[tuple[str, int, int]] = []
    unknown = []
    for name, row in cl.items():
        ds = tool.get(name)
        if ds is None:
            unknown.append(name)
            continue
        n_tool = len(ds)
        n_cl = _as_int(row.get("CL结构数"))
        if n_cl is None:
            continue
        t_multi, c_multi = n_tool >= 2, n_cl >= 2
        if t_multi and c_multi:
            both_multi += 1
        elif not t_multi and not c_multi:
            both_one += 1
        elif c_multi and not t_multi:
            miss.append((name, n_tool, n_cl))          # CL 有结构，工具只给 1 个域
        else:
            over.append((name, n_tool, n_cl))          # 工具给多域，CL 只看到 1 个结构

    print("=" * 74)
    print("工具 vs CL：域数一致性")
    print("=" * 74)
    print(f"参与对照（两边都填了 CL结构数）：{both_multi + both_one + len(miss) + len(over)} 个测点")
    print(f"  两边都说多域      {both_multi}")
    print(f"  两边都说均一      {both_one}")
    print(f"  ★ 工具漏检（CL≥2，工具=1）  {len(miss)}")
    for n, t, cln in miss:
        print(f"       {n:<22} 工具 {t} 域 / CL {cln} 结构")
    print(f"  ★ 工具过拟合（工具≥2，CL=1）{len(over)}")
    for n, t, cln in over:
        print(f"       {n:<22} 工具 {t} 域 / CL {cln} 结构")
    if unknown:
        print(f"\n⚠ 有 {len(unknown)} 个测点名在库里查不到（名字拼写？）：{unknown[:6]}")

    total = both_multi + both_one + len(miss) + len(over)
    if total:
        agree = both_multi + both_one
        print(f"\n一致率 {(100.0 * agree / total):.1f}%（{agree}/{total}）")
        print("口径：只比「是不是多域」这一件事。域数具体是 2 还是 3 不再细分 ——")
        print("      CL 上数结构本身带主观性，比到个位数会把读数噪声当结论。")
    print("\n提醒：'一致性'不等于'正确'。两边同错（比如都把一个混合带当成了")
    print("      一个真实的域）时，这张表是看不出来的。它只负责回答")
    print("      「工具漏掉了 CL 上明明看得见的结构吗」。")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="把 CL 图像（空间）与深度剖面（时间）对上：τ × 剥蚀速率 = 深度",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("三个子命令")[-1])
    ap.add_argument("--db", default=str(DB_DEFAULT), help="汇总库路径")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("list", help="打印工具侧的域结构 + τ/深度换算")
    p1.add_argument("--batch", default=None, help="只看某个批次（默认全部）")
    p1.add_argument("--rate", type=float, default=None,
                    help="剥蚀速率 μm/s（不给则只打印 τ，不换算深度）")
    p1.add_argument("--trim", type=float, default=TRIM_DEFAULT,
                    help=f"裁边秒数（默认 {TRIM_DEFAULT}）")
    p1.set_defaults(fn=cmd_list)

    p2 = sub.add_parser("template", help="生成 CL 判读空模板")
    p2.add_argument("--batch", default=None)
    p2.add_argument("--out", default="CL判读.csv")
    p2.set_defaults(fn=cmd_template)

    p3 = sub.add_parser("compare", help="读入 CL 判读表，算一致性")
    p3.add_argument("--cl", required=True, help="CL 判读表 CSV")
    p3.set_defaults(fn=cmd_compare)

    a = ap.parse_args(argv)
    c = con(Path(a.db))
    try:
        return a.fn(c, a)
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
