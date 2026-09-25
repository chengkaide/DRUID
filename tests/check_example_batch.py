# -*- coding: utf-8 -*-
"""
端到端数值回归 —— 拿仓库自带的示例批次 `examples/EX2022A/` 真跑一遍，
把关键数字钉死。

    python tests/check_example_batch.py

为什么它不叫 `test_*.py`
------------------------
这个检查要跑完整条流水线（85 个测点），2026-09-25 实测 **6.7 s**
（流水线本身 2.1 s，其余是「窗口尺度」那节又独立装载了一遍 85 个测点）。
放进 `pytest` / `run_all.py` 的默认集合里没必要 ——
"改了一行注释想快速验一下"不该顺带跑一条端到端。
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
6. **落地页那张三联图**（均一 / 核边分明 / 复杂变化）：它由
   `tools/gen_docs_split_figure.py` 从本批次算出三个真实测点并写进
   `docs/index.html`，是**发布出去的东西**。算法一改、图忘了重新生成，
   页面上的数字就会与真实结果不符。这里把页面上每个测点的数字都抓回来逐项比 ——
   于是"过期的图"变成一次确定的失败，而不是等读者发现。
   数字对不上时：重跑 `python tools/gen_docs_split_figure.py`，
   确认新数值是对的，再更新本文件顶部的常量。
7. **窗口尺度**：4 s 窗口里到底有多少个 cycle、每个通道有多少净计数
   （由 `tools/window_scale_facts.py` 现算）。文档里引用这些数字的地方有 8 处，
   而它们**曾经全是错的**：长期写着"窗口级 207Pb 只有约 130 个计数"当
   "207Pb/206Pb 不套 F(τ)"的第二条理由，可 130 的量级属于 **204Pb** ——
   错安了通道，在 8 个文件里活了好几轮。所以这里守两层：
   （a）现算的数字要与基线一致；（b）那三个错串**全仓禁用**，
   除非同一行里带着"曾经写/错安了通道"这类纠正标记（即那段警示本身）。

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
# 落地页那张三联图上的数字（图上是大字，只到十分位）。三个测点各代表一种形态：
# S43 均一（域表里没有它）、S34 核边分明（两段平台、无过渡带）、
# S18 复杂变化（三段平台 + 8 个窗口的过渡带）。
# 这些不是"另抄一份基线"，而是**从 docs/index.html 抓回来、与本次实跑比**——
# 抓不到才算失败，抓到了对不上也算失败。
# "积分"只在均一那个测点的图注里出现，其余为 None（不检查）。
DOC_FIG = {
    "S43": {"整段": 458.5, "窗口数": 35, "未归域": 0, "整段MSWD": 1.11, "上限": 1.49,
            "积分": 459.1, "域": {}},
    "S34": {"整段": 496.9, "窗口数": 35, "未归域": 0, "整段MSWD": 5.08, "上限": 1.49,
            "积分": None,
            "域": {"D1": (469.3, 13, 1.40), "D2": (514.9, 22, 0.88)}},
    "S18": {"整段": 458.6, "窗口数": 35, "未归域": 8, "整段MSWD": 3.39, "上限": 1.49,
            "积分": None,
            "域": {"D1": (471.6, 7, 0.15), "D2": (442.3, 12, 0.86),
                   "D3": (474.7, 8, 1.96)}},
}
# 窗口尺度（`tools/window_scale_facts.py` 现算）。文档与代码注释里引用的
# "12~13 个 cycle""207Pb 净计数中位约 9×10³""204Pb 净计数中位约 37"就是这几个数。
# 口径：Σ(net cps) × cycle 周期，全部窗口（不是"测点中位再取中位"）。
SCALE = {
    "cycle_s": 0.31871,
    "cycle_med": 13,          # 每窗口 cycle 数
    "cycle_range": (12, 13),
    "cnt207_k": 9,            # 207Pb 窗口净计数中位，单位千
    "cnt204_med": 37,         # 204Pb 窗口净计数中位（实测 36.6）
    "ratio_207_204_lo": 150,  # 两者的比值区间 —— 这条专门盯"错安通道"那类错误
    "ratio_207_204_hi": 400,
    "n_spot": 83,             # 玻璃标样不参与（与 load_batch 同口径）
}
#: 引用了"207Pb 净计数中位约 N×10³"这个说法的页面/文件。
#: 只要求发布出去的主页与接手必读的 AGENTS.md 必须带这个数 ——
#: 其余位置由"错串全仓禁用"兜住，免得换个说法就误报。
SCALE_QUOTERS = ("docs/index.html", "AGENTS.md")
#: 引用了"204Pb 净计数中位约 N"、且必须被核对的文件。
#: 这个数不是给读者看的，它是 A-11 那次"错安通道"的反证 ——
#: 130 的量级本来属于它，所以它必须一直挂在实测值上。
SCALE_QUOTERS_204 = ("druid/__init__.py", "druid/depth/windows.py",
                     "tests/check_example_batch.py",
                     "AGENTS.md", "README.md")
#: 出处里必须出现的纠正标记：带着它的行才允许出现旧错串
#: （"长期写着…"这类句式就是那段警示本身，不是"又把错数抄回来了"）。
EXEMPT_MARKERS = ("曾经写", "错安了通道", "STALE_CLAIMS", "长期写着")

TOL_AGE = 1e-3         # Ma。基线是从导出的 xlsx 里取的（导出时 round(4)），
                       # 内存里的完整精度与之可能差 1e-4 量级，所以留 1e-3。
TOL_PCT = 1e-4         # 百分点


def _viol(rel: str, line_no: int, what: str) -> str:
    """旧错串的违规描述。

    报错里直接写清怎么改：**引用旧数是允许的，只要同一行说明它是错的**
    （带 EXEMPT_MARKERS 里的标记）—— 因为"记下这个数曾经是错的"本身就是
    最有效的防复发手段，不该被禁掉。
    """
    return (f"{rel}:{line_no} 仍有「{what}」（要么改掉，要么在同一行"
            "写清它是错的，如“曾经写…/错安了通道”）")


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
    print("=== 文档里的「分域前后」三联图 ===")
    doc = ROOT / "docs" / "index.html"
    if not getattr(ov, "empty", True) and doc.is_file() and "gen_docs_split_figure" in \
            doc.read_text(encoding="utf-8"):
        txt = doc.read_text(encoding="utf-8")

        # 图例是"一个测点一行"（li.case）。必须先按行切开再抓 ——
        # 在一整页上直接 re.search，某个测点的数会串到另一个测点头上，
        # 而且串错了照样全部通过。
        rows = dict(re.findall(r'<li class="case"><b>[^<]*?(S\d\d)(.*?)</li>', txt, re.S))
        if not rows:
            bad.append("文档里找不到三联图的图例行 —— 图被改坏或没生成？")
            print("  FAIL 找不到图例行（li.case）")

        def from_doc(pat, cast, label, src):
            m = re.search(pat, src)
            if not m:
                bad.append(f"文档图里找不到「{label}」—— 图被改坏或没生成？")
                print(f"  FAIL 文档图里找不到「{label}」")
                return None
            return cast(m.group(1))

        scraped = {}
        for spot, anchor in DOC_FIG.items():
            row = rows.get(spot)
            if row is None:
                bad.append(f"文档图里没有测点 {spot} 的那一行")
                print(f"  FAIL 文档图里没有测点 {spot}")
                continue
            live = ov[ov["样品"] == spot]
            if not len(live):
                bad.append(f"实跑结果里没有测点 {spot}，而文档图里有")
                print(f"  FAIL 实跑里没有测点 {spot}")
                continue
            live = live.iloc[0]
            d0 = dom[dom["样品"] == spot].sort_values("tau")

            got = {
                "整段": from_doc(r"整段 ([\d.]+) Ma</b>", float, f"{spot} 整段", row),
                "窗口数": from_doc(r"共 (\d+) 个滑窗", int, f"{spot} 窗口数", row),
                "整段MSWD": from_doc(r"整段 MSWD ([\d.]+) &", float, f"{spot} 整段 MSWD", row),
                "上限": from_doc(r"上限 ([\d.]+) →", float, f"{spot} 判据上限", row),
            }
            scraped[spot] = got

            # ① 与本次实跑比
            check(f"图上 {spot} 整段 (Ma)", got["整段"], round(float(live["年龄_Ma"]), 1))
            check(f"图上 {spot} 窗口数", got["窗口数"], int(live["n_win"]))
            check(f"图上 {spot} 整段 MSWD", got["整段MSWD"], round(float(live["MSWD"]), 2))
            check(f"图上 {spot} 判据上限", got["上限"], round(float(live["相容上限"]), 2))

            # 均一测点在域表里没有行 —— 那是"没有域可拆"，不是"35 个窗口全被剥掉"，
            # 所以这里必须先判空，否则相减会得到一个很唬人的 35。
            n_off = 0 if not len(d0) else int(live["n_win"]) - int(d0["n_win"].sum())
            if not len(d0):
                # 均一测点：域表里没有行，图注里会另给一个"整段积分"的数
                got["积分"] = from_doc(r"整段积分（另一种算法）给 ([\d.]+) Ma", float,
                                      f"{spot} 整段积分", row)
                check(f"图上 {spot} 整段积分 (Ma)", got["积分"],
                      round(float(live["整段积分年龄_Ma"]), 1))
            if n_off:
                got["未归域"] = from_doc(r"(\d+) 个窗口没通过域内相容性检验", int,
                                        f"{spot} 未归域窗口数", row)
                check(f"图上 {spot} 未归域窗口", got["未归域"], n_off)
            else:
                got["未归域"] = 0
                # 均一测点没有"被剥掉的窗口"这回事，图上用的是另一句话
                key = "没有一个窗口被剥掉" if len(d0) else "域表里根本没有这个测点"
                if key in row:
                    print(f"  ok   图上 {spot} 未归域窗口                  0")
                else:
                    bad.append(f"{spot} 实跑是 0 个未归域窗口，图上却没说（应有「{key}」）")
                    print(f"  FAIL {spot} 未归域：实跑 0，图上没写「{key}」")

            # ② 逐域比（域号 / 年龄 / 窗口数 / MSWD 四项一起，少一项都算不符）
            got_dom = {n: (float(a), int(w), float(m)) for n, a, w, m in re.findall(
                r"(D\d) ([\d.]+) Ma（(\d+) 窗，MSWD ([\d.]+)）", row)}
            want_dom = {str(r["域"]): (round(float(r["年龄_Ma"]), 1), int(r["n_win"]),
                                       round(float(r["MSWD"]), 2)) for _, r in d0.iterrows()}
            if got_dom == want_dom:
                print(f"  ok   图上 {spot} 的 {len(want_dom)} 个域（年龄/窗口数/MSWD）")
            else:
                bad.append(f"文档图里 {spot} 的域与实跑不符：页面 {got_dom}，实跑 {want_dom}")
                print(f"  FAIL {spot} 域：页面 {got_dom} ≠ 实跑 {want_dom}")

        # ③ 与静态锚比：两份独立来源同时对得上，才说明图既没手改也没过期
        for spot, anchor in DOC_FIG.items():
            got = scraped.get(spot)
            if not got:
                continue
            for k, v in anchor.items():
                if k == "域" or v is None:
                    continue          # 域已在 ② 里按四项比过；None = 该测点没有这个数
                if k not in got:
                    bad.append(f"文档图里没抓到 {spot} 的「{k}」，而静态锚里有")
                    print(f"  FAIL 图上 {spot} 缺 {k}")
                    continue
                check(f"图上 {spot} {k}（静态锚）", got[k], v)
        print("  （三个测点各自代表一种形态；重生成：python tools/gen_docs_split_figure.py）")
    else:
        print("  skip  落地页里没有这张自动生成的图，或不分域表为空")

    print()
    print("=== 窗口尺度（一个 4 s 窗口有多大）===")
    # 数字现算，不抄；口径写在 tools/window_scale_facts.py 的 docstring 里。
    sys.path.insert(0, str(ROOT / "tools"))
    import window_scale_facts as wsf                       # noqa: E402

    sc = wsf.measure(BATCH)
    check("cycle 周期 (s)", round(sc["cycle_s"], 5), SCALE["cycle_s"], 1e-4)
    check("每窗口 cycle 数 中位", int(round(sc["n_cycles_med"])), SCALE["cycle_med"])
    check("每窗口 cycle 数 范围",
          (sc["n_cycles_min"], sc["n_cycles_max"]), tuple(SCALE["cycle_range"]))
    check("参与统计的测点数（玻璃标样不算）", sc["n_spot"], SCALE["n_spot"])

    cnt207 = sc["cnt"][207]["med"]
    cnt204 = sc["cnt"][204]["med"]
    check("207Pb 窗口净计数 中位（千）", int(round(cnt207 / 1e3)), SCALE["cnt207_k"])
    check("204Pb 窗口净计数 中位", int(round(cnt204)), SCALE["cnt204_med"])
    # 这一条才是关键：错安通道那次，就是把 204Pb 的量级写成了 207Pb 的。
    # 两个通道差了 2 个数量级以上，只要谁再抄错一次，这里立刻红。
    # 不断言具体比值（两个中位数相除会随小改动抖动），只卡数量级。
    ratio = cnt207 / cnt204
    if not (SCALE["ratio_207_204_lo"] < ratio < SCALE["ratio_207_204_hi"]):
        bad.append(f"207Pb 窗口计数（{cnt207:,.0f}）与 204Pb（{cnt204:,.0f}）之比 "
                   f"{ratio:.0f} 跑到 {SCALE['ratio_207_204_lo']}~"
                   f"{SCALE['ratio_207_204_hi']} 之外 —— 是不是又把通道写反了？")
        print(f"  FAIL 207Pb / 204Pb 窗口计数之比 {ratio:.0f} 不在合理区间")
    else:
        print(f"  ok   207Pb / 204Pb 窗口计数之比              {ratio:.0f}"
              f"   (期望 {SCALE['ratio_207_204_lo']}~{SCALE['ratio_207_204_hi']})")

    # 引用这个数字的地方：页面上的数必须落在现算值的 ±20% 内。
    for rel in SCALE_QUOTERS:
        path = ROOT / rel
        if not path.is_file():
            bad.append(f"找不到 {rel}（SCALE_QUOTERS 里登记了它）")
            print(f"  FAIL 找不到 {rel}")
            continue
        m = re.search(r"Pb 净计数中位约 (\d+)×10³", path.read_text(encoding="utf-8"))
        if not m:
            bad.append(f"{rel} 里找不到「207Pb 净计数中位约 N×10³」这个说法 —— "
                       "要么句子被改写了，要么数字被删了")
            print(f"  FAIL {rel} 里没有引用的窗口计数")
            continue
        quoted = int(m.group(1)) * 1000
        # 只要求"量级对得上"（±20%）：页面写的是"约"，钉死到个位会让
        # 改一句措辞就得动基线。真正的回归保护在上面那两条现算基线里。
        off = abs(quoted - cnt207) / cnt207
        if off > 0.2:
            bad.append(f"{rel} 写着 {quoted:.0f}，实测 {cnt207:.0f}"
                       f"（差 {off * 100:.0f}%，超出 ±20%）")
            print(f"  FAIL {rel} 引用的窗口计数与实测差得太远")
        else:
            print(f"  ok   {rel} 引用 {quoted:,.0f} vs 实测 {cnt207:,.0f}"
                  f"（差 {off * 100:.0f}%）")

    # 204Pb 那一侧的同类保护：引用"204Pb 净计数中位约 N"的地方，
    # 数必须落在现算值的 ±20% 内。和上面的比值检查一起把"错安通道"钉死。
    # 只认"净计数中位约 N"这一种写法：不跨行，且必须带"约"。
    # 放宽一点就会连本段自己的正则、以及"（实测 36.6）"这类旁注一起抓进来。
    rx204 = re.compile(r"204Pb[^\n]{0,20}?净计数中位约 (\d+)")
    for rel in SCALE_QUOTERS_204:
        path = ROOT / rel
        if not path.is_file():
            bad.append(f"找不到 {rel}（SCALE_QUOTERS_204 里登记了它）")
            print(f"  FAIL 找不到 {rel}")
            continue
        ms = rx204.findall(path.read_text(encoding="utf-8"))
        if not ms:
            bad.append(f"{rel} 里找不到「204Pb 净计数中位约 N」这个说法")
            print(f"  FAIL {rel} 里没有引用的 204Pb 计数")
            continue
        for s in ms:
            quoted = int(s)
            off = abs(quoted - cnt204) / cnt204
            if off > 0.2:
                bad.append(f"{rel} 把 204Pb 写成 {quoted}，实测 {cnt204:.1f}"
                           f"（差 {off * 100:.0f}%）—— 是不是又和 207Pb 换位了？")
                print(f"  FAIL {rel} 引用的 204Pb 计数与实测差得太远")
            else:
                print(f"  ok   {rel} 引用 204Pb {quoted} vs 实测 {cnt204:.1f}"
                      f"（差 {off * 100:.0f}%）")

    # 错串全仓禁用。带着"曾经写 / 错安了通道"标记的行是警示本身，豁免。
    # 用 set：同一条文本既被子串命中、又被正则命中时只报一次。
    hits = set()
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        if any(rel.startswith(x) for x in (".git/", "build_tmp/", "dist/", ".workbuddy/")):
            continue
        if path.suffix.lower() not in (".py", ".md", ".html", ".txt", ".toml"):
            continue
        if rel == "tools/window_scale_facts.py":
            continue          # 这个文件是那份禁用清单的定义处
        for i, line in enumerate(
                path.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            if any(mk in line for mk in EXEMPT_MARKERS):
                continue
            for pat in wsf.STALE_CLAIMS:
                if pat in line:
                    hits.add(_viol(rel, i, pat))
            # 子串只挡得住"照抄原文"，所以同一批谬误还要按正则再扫一遍：
            # 实测它曾经写成另一种样子（"25–35 个采样周期"）藏在另一页里。
            for rx in wsf.STALE_PATTERNS:
                mm = re.search(rx, line)
                if mm:
                    hits.add(_viol(rel, i, mm.group(0)))
    if hits:
        for h in hits:
            bad.append(h)
            print(f"  FAIL {h}")
    else:
        print(f"  ok   那 {len(wsf.STALE_CLAIMS)} 个旧错串全仓找不到"
              f"（重算：python tools/window_scale_facts.py）")

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
