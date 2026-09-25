# -*- coding: utf-8 -*-
"""
端到端数值回归 —— 拿仓库自带的示例批次 `examples/EX2022A/` 真跑一遍，
把关键数字钉死。

    python tests/check_example_batch.py

为什么它不叫 `test_*.py`
------------------------
这个检查要跑完整条流水线（85 个测点，约 1 分钟）。放进 `pytest` / `run_all.py`
的默认集合里，会让"改了一行注释想快速验一下"变成等一分钟。
所以它和 `check_wheel.py` 一样是**显式调用**的脚本：
日常改动跑 `tests/run_all.py`，改动碰到数值链路时再跑这个，CI 里每次都跑。

它钉的是什么
------------
1. **标样 QC**：91500 与 Ple 的加权平均年龄与偏差。这是整批数据可信度的根，
   标样不对，样品年龄就没有意义。
2. **样品年龄的分布**（中位与 5–95% 区间）。注意必须是**只统计样品**——
   这里曾经有过一个 bug：`DataFrame.get()` 返回整表的那一列，
   把 35 个标样测点混进了"样品年龄"，中位从 458.1 抬到 460.9，
   区间从 420~644 撑成 333~1044 Ma。
3. **行数与角色分布**：序列表有没有被正确解析（顺带覆盖 CSV 版序列表）。
4. **深度结构分布**：域判别的结果，改动 BIC / MSWD 判据会立刻反映到这里。
5. **不分域（整段）口径**：整段加权平均的年龄与"与常数模型相容"的测点数。
   这条不经过分域，是**独立于第 4 条**的另一条判据：两者一起变，说明问题在
   更上游（归一化、外部重现性）；只有一条变，说明问题就在那一层。
6. **落地页那张「分域前后」图**：它由 `tools/gen_docs_split_figure.py` 从本批次
   算出来并写进 `docs/index.html`，是**发布出去的东西**。算法一改、图忘了重新生成，
   页面上的数字就会与真实结果不符。这里把页面上的数字抓回来逐项比 ——
   于是"过期的图"变成一次确定的失败，而不是等读者发现。
   数字对不上时：重跑 `python tools/gen_docs_split_figure.py`，
   确认新数值是对的，再更新本文件顶部的常量。

改了算法就有数字变化是正常的 —— 那时候要**重新确认基线并更新这里的常量**，
而不是把容差放宽。容差放宽等于把这个检查废掉。
"""
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BATCH = ROOT / "examples" / "EX2022A"

# ── 基线（2026-09-19 由 examples/EX2022A 实跑得到，与真实批次逐格一致）──────
ROWS = {"结果": 83, "标样QC": 2, "深度剖面域": 91, "剖面窗口": 1677,
        "不分域年龄": 48}
ROLES = {"样品": 48, "主标": 21, "监控标样": 14}
QC = {                       # 标样 -> (加权平均年龄 Ma, 偏差 %)
    "91500": (1059.7543, -0.2490),
    "Ple": (343.2381, 1.8118),
}
AGE = {"中位": 458.0686, "5%": 420.1799, "95%": 643.5333}   # 只统计样品
CONC = {"中位": 101.3006, "90-110%占比": 0.9375}
STRUCT = {None: 35, "多域(2)": 18, "均一": 15, "多域(3)": 15}
# 不分域（整段）口径 —— 2.4.0 新增的那条路径：整段加权平均 + MSWD 相容判定。
# 它不经过 BIC 分域，所以动 `whole_spot_stats` / `mswd_acceptance` /
# `weighted_mean` 都会立刻反映到这里。首点锚（S01）是防"全员同步偏移"
# 这类从分布上看不出来的变化。
WHOLE = {"行数": 48, "整段常数": 5, "MSWD中位": 4.4376,
         "首点年龄": 459.0069, "首点主域": 453.8410}
# 落地页那张图上的数字（图上是大字，只到十分位）。
# 这些不是"另抄一份基线"，而是**从 docs/index.html 抓回来、与本次实跑比**——
# 抓不到才算失败，抓到了对不上也算失败。
DOC_FIG = {"整段": 459.0, "D1": 453.8, "D2": 494.3, "未归域": 18, "窗口数": 35}

TOL_AGE = 1e-3         # Ma。基线是从导出的 xlsx 里取的（导出时 round(4)），
                       # 内存里的完整精度与之可能差 1e-4 量级，所以留 1e-3。
TOL_PCT = 1e-4         # 百分点


def main() -> int:
    # 这里的输出含中文；Windows 上重定向时会踩编码坑（见 druid/console.py）
    from druid.console import ensure_utf8_streams
    ensure_utf8_streams()

    if not BATCH.is_dir():
        print(f"找不到示例批次 {BATCH}")
        print("（它应该在仓库里：examples/EX2022A/ —— 见 .gitignore 的说明）")
        return 1

    from druid import __version__
    from druid.workflow import BatchConfig, run_batch

    print(f"DRUID {__version__}  ·  示例批次 {BATCH.name}")
    print("-" * 62)
    t0 = time.time()
    cfg = BatchConfig(data_dir=str(BATCH), plot=False)
    result = run_batch(cfg)
    dt = time.time() - t0

    res, qc = result.results, result.qc
    dom, win = result.domains, result.windows
    print(f"跑完，用时 {dt:.1f} s")

    bad = []

    def check(label, got, want, tol=None):
        if tol is None:
            ok = got == want
        else:
            ok = abs(got - want) <= tol
        flag = "ok  " if ok else "FAIL"
        if not ok:
            bad.append(f"{label}: 得到 {got!r}，应为 {want!r}")
        print(f"  {flag} {label:34s} {got!r}"
              + (f"   (基线 {want!r})" if not ok else ""))

    print()
    print("=== 行数 ===")
    for name, want in ROWS.items():
        got = {"结果": len(res), "标样QC": len(qc), "深度剖面域": len(dom),
               "剖面窗口": len(win), "不分域年龄": len(result.overall)}[name]
        check(f"{name} 行数", got, want)

    print()
    print("=== 角色分布（序列表解析）===")
    for role, want in ROLES.items():
        check(f"{role} 测点数", int((res["类型"] == role).sum()), want)

    print()
    print("=== 标样 QC ===")
    for std, (age, dev) in QC.items():
        row = qc[qc["标样"] == std]
        if not len(row):
            bad.append(f"标样QC 里没有 {std}")
            print(f"  FAIL 缺少 {std}")
            continue
        check(f"{std} 加权平均年龄 (Ma)", round(float(row["加权平均年龄_Ma"].iloc[0]), 4),
              age, TOL_AGE)
        check(f"{std} 偏差 (%)", round(float(row["偏差_pct"].iloc[0]), 4),
              dev, TOL_PCT)

    print()
    print("=== 样品年龄（只统计样品，不含标样）===")
    unk = res[res["类型"] == "样品"]
    col = ("年龄206_238_QC校正"
           if "年龄206_238_QC校正" in res.columns else "年龄206_238")
    a = unk[col]
    check("中位 (Ma)", round(float(a.quantile(0.5)), 4), AGE["中位"], TOL_AGE)
    check("5% (Ma)", round(float(a.quantile(0.05)), 4), AGE["5%"], TOL_AGE)
    check("95% (Ma)", round(float(a.quantile(0.95)), 4), AGE["95%"], TOL_AGE)
    # 这一条专门盯住上面说的那个 bug：若统计口径退回整表，区间会立刻撑大
    if not (400 < float(a.quantile(0.05)) < 440):
        bad.append("5% 分位跑到 400~440 Ma 之外 —— 统计口径可能又混进标样了")
        print("  FAIL 5% 分位不在 400~440 Ma 内（标样是否混进来了？）")

    print()
    print("=== 协和度 ===")
    c = unk["协和度_pct"]
    check("中位 (%)", round(float(c.median()), 4), CONC["中位"], TOL_PCT)
    check("90–110% 占比", round(float(c.between(90, 110).mean()), 4),
          CONC["90-110%占比"], 1e-6)

    print()
    print("=== 深度结构分布 ===")
    st = res["深度结构"]
    # ⚠ 内存里"没有剖面"是**空字符串**，写成 Excel 再读回来才变成 NaN。
    # 只数 isna() 会得到 0 而基线是 35（这里踩过一次）。
    # 所以两种表示都要算上。
    blank = st.isna() | (st.astype(str).str.strip() == "")
    check("（单点/无剖面）", int(blank.sum()), STRUCT[None])
    vc = st.value_counts().to_dict()
    for k, want in STRUCT.items():
        if k is None:
            continue
        check(k, int(vc.get(k, 0)), want)

    print()
    print("=== 不分域（整段）口径 ===")
    ov = result.overall
    if getattr(ov, "empty", True):
        # 空表必须判失败：它是"深度分析根本没跑"的信号，
        # 而不是"这批数据没问题"。
        bad.append("不分域年龄表是空的 —— 深度分析没跑到？")
        print("  FAIL 不分域年龄表为空")
    else:
        check("行数（= 样品测点数）", len(ov), WHOLE["行数"])
        is_const = ov["判定"].astype(str) == "整段常数"
        check("整段常数测点数", int(is_const.sum()), WHOLE["整段常数"])
        check("MSWD 中位", round(float(ov["MSWD"].median()), 4),
              WHOLE["MSWD中位"], 1e-3)
        check("首点整段年龄 (Ma)", round(float(ov["年龄_Ma"].iloc[0]), 4),
              WHOLE["首点年龄"], TOL_AGE)
        check("首点主域年龄 (Ma)", round(float(ov["主域年龄_Ma"].iloc[0]), 4),
              WHOLE["首点主域"], TOL_AGE)

    print()
    print("=== 文档里的「分域前后」图 ===")
    doc = ROOT / "docs" / "index.html"
    if not getattr(ov, "empty", True) and doc.is_file() and "gen_docs_split_figure" in \
            doc.read_text(encoding="utf-8"):
        txt = doc.read_text(encoding="utf-8")

        def from_doc(pat, cast, label):
            m = re.search(pat, txt)
            if not m:
                bad.append(f"文档图里找不到「{label}」—— 图被改坏或没生成？")
                print(f"  FAIL 文档图里找不到「{label}」")
                return None
            return cast(m.group(1))

        first = ov["样品"].iloc[0]
        d0 = dom[dom["样品"] == first]
        got = {
            "整段": from_doc(r"整段不分域 ([\d.]+) Ma</b>", float, "整段不分域"),
            "D1": from_doc(r"D1 域均值 ([\d.]+) Ma</b>", float, "D1 域均值"),
            "D2": from_doc(r"D2 域均值 ([\d.]+) Ma</b>", float, "D2 域均值"),
            "未归域": from_doc(r"未进入任何年龄域：(\d+) 个窗口", int, "未归域窗口数"),
            "窗口数": from_doc(r"同一个剥蚀坑的 (\d+) 个滑窗", int, "窗口数"),
        }
        want = {
            "整段": round(float(ov["年龄_Ma"].iloc[0]), 1),
            "D1": round(float(ov["主域年龄_Ma"].iloc[0]), 1),
            "D2": round(float(d0["年龄_Ma"].iloc[1]), 1) if len(d0) > 1 else None,
            "未归域": int(ov["n_win"].iloc[0]) - int(d0["n_win"].sum()),
            "窗口数": int(ov["n_win"].iloc[0]),
        }
        for k, w in want.items():
            g = got.get(k)
            if g is None:
                continue          # 已经在 from_doc 里记过失败了
            if w is None:
                print(f"  skip  {k}：本批次的该测点只有一个域，图上没有 D2")
                continue
            check(f"图上 {k}", g, w)
        for k, v in DOC_FIG.items():
            check(f"图上 {k}（静态锚）", got.get(k), v)
        print(f"  （图取自测点 {first}；重生成：python tools/gen_docs_split_figure.py）")
    else:
        print("  skip  落地页里没有这张自动生成的图，或不分域表为空")

    print()
    if bad:
        print(f"✗ {len(bad)} 项与基线不符：")
        for b in bad:
            print(f"    {b}")
        print()
        print("确认是有意改动的话，请核对后更新本文件顶部的基线常量。")
        return 1
    print("✓ 全部与基线一致")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
