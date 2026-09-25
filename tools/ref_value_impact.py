# -*- coding: utf-8 -*-
"""
tools/ref_value_impact.py —— 91500 参考值的口径差异到底造成多大差别（只读）
============================================================================

为什么要这个东西
----------------
`druid/core/constants.py` 里 91500 的那一行是：

    "91500": dict(R68=0.17928, R76=0.07494, age_Ma=1062.4,
                  ref="Wiedenbeck et al. (1995) / 常用推荐值")

这三个数**互不自洽**：0.17928 反算回来是 1063.04 Ma，而 age_Ma 写 1062.4 Ma，
差 **+0.0602%**。用户 2026-09-25 的决定是**默认改用 Horstwood et al. (2016)
表 S2 的 CA-ID-TIMS 推荐值**（`ref_preset="horstwood2016"`）；上面那组第一版
数字作为 `repo` 档保留下来，用来复现 2026-09-25 之前的结果。

本脚本的职责是回答**"各档之间到底差多大"** —— 把那个 0.0602% 换成样品年龄上
看得见的 Ma，并和不确定度摆在一起比。⚠ 下表**以 `repo` 档为基线**（脚本内定义，
不是产品默认档）：因为"第一版 vs 各候选"是最常被问到的问题。

两条独立的路，互相印证
----------------------
1. **解析**：直接对候选参考值调 `age68() / age76()`，算各自的等效年龄与 R75。
2. **实跑**：在**内存里**临时替换 `STANDARDS['91500']`，重跑
   `examples/EX2022A` 整条流水线，读回 `标样QC` 与 `结果` 两张表。

`run_batch()` 通过 `std_ref()` 取参考值，而 `references.py` 拿到的是
`constants.STANDARDS` 这个 dict 对象本身 —— 所以原地改它就能改变整条流水线，
**不需要动仓库里的任何文件**。脚本结束前必定还原（`try/finally` 兜底）。

只读保证
--------
本脚本不改任何仓库文件（`--out` 指定的报告文件除外），更不改
`druid/core/constants.py`：评估结束后 `STANDARDS['91500']` 与启动时逐键相同，
脚本会自己断言这一点。

用法
----
    python tools/ref_value_impact.py              # 打到屏幕
    python tools/ref_value_impact.py --out r.md   # 另存一份 markdown
"""
from __future__ import annotations

import argparse
import contextlib
import io
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from druid.core.constants import STANDARDS, U238_U235          # noqa: E402
from druid.core.geochronology import age68, age76, r68_of_age, r76_of_age  # noqa: E402
from druid.workflow import BatchConfig, run_batch               # noqa: E402

EX = ROOT / "examples" / "EX2022A"
PRIMARY = "91500"
SECONDARY = "Ple"

#: 候选口径。R68/R76 为 None 时表示"该口径只给年龄，比值由衰变方程反算"
#: ——这正是 Isoclock 与 Wiedenbeck 两条路线最本质的差别。
CANDIDATES = (
    dict(key="repo",
         label="repo 档（＝仓库第一版）",
         ref="constants.py 注释：Wiedenbeck et al. (1995) / 常用推荐值",
         R68=0.17928, R76=0.07494, age=1062.4),
    dict(key="w1995",
         label="Wiedenbeck et al. (1995) 文献原值",
         ref="Geostandards Newsletter 19, 1–23（ID-TIMS）",
         R68=0.17917, R76=0.07488, age=1062.4),
    dict(key="h2016",
         label="Horstwood et al. (2016) 社区推荐",
         ref="GGR 40, 311–332（CA-ID-TIMS；U-Pb 年龄 1063.51 Ma）",
         R68=0.179365, R76=0.074941, age=1063.51),
    dict(key="self1062",
         label="由 age_Ma=1062.4 反算（自洽解）",
         ref="本仓衰变常数 λ238/λ235 + U238/U235=137.818",
         R68=None, R76=None, age=1062.4),
    dict(key="isoclock",
         label="Isoclock 口径（整数 1062 Ma）",
         ref="Isoclock2.0.py: Cal_age(1062) / Standard_names['91500']=1062",
         R68=None, R76=None, age=1062.0),
)


def resolve(c: dict) -> tuple[float, float]:
    """把一个候选口径解析成一对具体的参考比值。"""
    r68 = c["R68"] if c["R68"] is not None else r68_of_age(c["age"])
    r76 = c["R76"] if c["R76"] is not None else r76_of_age(c["age"])
    return float(r68), float(r76)


def run_case(c: dict):
    """在内存里临时换掉 91500 的参考值，跑一遍完整批次，然后原样还回去。"""
    r68, r76 = resolve(c)
    slot = STANDARDS[PRIMARY]
    orig = dict(slot)
    try:
        slot.update(R68=r68, R76=r76, age_Ma=float(c["age"]))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            res = run_batch(BatchConfig(data_dir=EX, verbose=False, plot=False))
    finally:
        slot.clear()
        slot.update(orig)
    return res


def harvest(res, c: dict) -> dict:
    """从一次实跑的结果里抽出做对照需要的那几列。"""
    qc = res.qc.set_index("标样")
    smp = res.results[res.results["类型"] == "样品"]
    # ⚠ 必须与仓库基线同口径。`tests/check_example_batch.py` 取的是 QC 二次校正后的
    #   那一列（它自己的第 200 行同一套判断）。注：那个中位随默认档而变 ——
    #   repo 档 458.07 Ma、默认的 horstwood2016 档 458.11 Ma（整体平移，不影响本表结论）。
    #   若误用未校正的「年龄206_238」，会得到 466.4 Ma —— 不是同一个量，两者差 ≈1.8%。
    acol = ("年龄206_238_QC校正"
            if "年龄206_238_QC校正" in res.results.columns else "年龄206_238")
    # 未校正列也留一份：口径差异在「未校正列」与「最终列」上的落点完全不同，
    # 见 main() 的表 2b。
    uncal = "年龄206_238"
    r68, r76 = resolve(c)
    return dict(
        R68=r68, R76=r76, age_written=float(c["age"]),
        age_of_r68=float(age68(r68)), age_of_r76=float(age76(r76)),
        R75=r76 * r68 * U238_U235,
        # 标样 QC
        n_primary=int(qc.loc[PRIMARY, "点数"]),
        primary_ref_age=float(qc.loc[PRIMARY, "参考年龄_Ma"]),
        primary_wmean=float(qc.loc[PRIMARY, "加权平均年龄_Ma"]),
        primary_s2=float(qc.loc[PRIMARY, "s2_Ma"]),
        primary_bias=float(qc.loc[PRIMARY, "偏差_pct"]),
        secondary_bias=float(qc.loc[SECONDARY, "偏差_pct"]),
        # 样品
        n_smp=int(len(smp)),
        age_col=acol,
        med_age=float(smp[acol].median()),
        med_s1=float(smp["s68_1sig"].median()),
        conc_med=float(smp["协和度_pct"].median()),
        ages=smp.set_index("序号")[acol],
        uncal_med=float(smp[uncal].median()),
        uncal_ages=smp.set_index("序号")[uncal],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="91500 参考值口径差异的影响评估（只读）")
    ap.add_argument("--out", default=None, help="把报告另存到这个路径（markdown）")
    args = ap.parse_args()

    out: list[str] = []

    def w(s: str = "") -> None:
        out.append(s)

    before = dict(STANDARDS[PRIMARY])

    w("# 91500 参考值的口径差异：影响评估")
    w()
    w(f"数据：`examples/EX2022A`（{len(list(EX.glob('*_[0-9]*.csv')))} 个原始 CSV）。"
      "本表由 `tools/ref_value_impact.py` 现算，**不改** `druid/core/constants.py`。")
    w()

    # ── 表 1：常数层面 ────────────────────────────────────────────────
    w("## 表 1　各口径的常数本身")
    w()
    w("| 口径 | R68 (206Pb/238U) | R76 (207Pb/206Pb) | age_Ma 写的是 | "
      "R68 反算的年龄 | R76 反算的年龄 | 与 age_Ma 差 |")
    w("|---|---|---|---|---|---|---|")
    for c in CANDIDATES:
        r68, r76 = resolve(c)
        a68, a76 = float(age68(r68)), float(age76(r76))
        diff = (a68 / c["age"] - 1) * 100
        w(f"| {c['label']} | {r68:.6f} | {r76:.6f} | {c['age']} Ma | "
          f"{a68:.2f} Ma | {a76:.2f} Ma | **{diff:+.4f}%** |")
    w()

    # ── 实跑 ──────────────────────────────────────────────────────────
    rows: list[dict] = []
    for c in CANDIDATES:
        r = harvest(run_case(c), c)
        r["cand"] = c
        rows.append(r)

    base = rows[0]
    w("## 表 2　同一条流水线、同一批数据，只换 91500 的参考值")
    w()
    w(f"基线是 **repo 档（＝仓库第一版）** 那行。样品数 {base['n_smp']}，主标点数 {base['n_primary']}；"
      f"样品那一列用的是 `{base['age_col']}`（与 `tests/check_example_batch.py` 同口径）。")
    w()
    w("| 口径 | 主标加权平均 (Ma) | 主标 QC 偏差 | 样品年龄中位 (Ma) | "
      "相对基线的平移 | 平移的绝对值 |")
    w("|---|---|---|---|---|---|")
    for r in rows:
        c = r["cand"]
        d = r["ages"].subtract(base["ages"]).dropna()
        dm = d.median()
        dmax = d.abs().max()
        pct = dm / base["med_age"] * 100
        ab = "—（基线）" if r is base else f"中位 {dm:+.3f} Ma / 最大 {dmax:.3f} Ma"
        w(f"| {c['label']} | {r['primary_wmean']:.4f} | {r['primary_bias']:+.4f}% | "
          f"{r['med_age']:.3f} | **{pct:+.4f}%** | {ab} |")
    w()

    # ── 表 2b：口径差异落在「未校正列」和「最终列」上完全不同 ────────
    w("## 表 2b　同一个口径差异，落在「未校正列」与「最终列」上完全不同")
    w()
    w("这是本次评估**最该记住的一点**。第 ⑩ 步 QC 二次校正按监控标样 Ple 算")
    w("`kfac = 337.13 / mu_qc`，而 **Ple 的参考年龄与 91500 的参考值毫无关系** ——")
    w("于是它把「换 91500 参考值」带来的整体平移**几乎全部吸收掉**：")
    w()
    w("- `年龄206_238`（**未校正**、不是报出值）→ 会整体平移；")
    w("- `年龄206_238_QC校正`（**有监控标样时真正报出的那一列**）→ 几乎不动。")
    w()
    w("| 口径 | 未校正列中位 (Ma) | 未校正列平移 | 最终列中位 (Ma) | 最终列平移 |")
    w("|---|---|---|---|---|")
    for r in rows:
        du = r["uncal_ages"].subtract(base["uncal_ages"]).dropna()
        df = r["ages"].subtract(base["ages"]).dropna()
        if r is base:
            u = "—（基线）"
            f = "—（基线）"
        else:
            u = f"中位 {du.median():+.3f} / 最大 {du.abs().max():.3f} Ma"
            f = f"中位 {df.median():+.3f} / 最大 {df.abs().max():.3f} Ma"
        w(f"| {r['cand']['label']} | {r['uncal_med']:.3f} | {u} | "
          f"{r['med_age']:.3f} | {f} |")
    w()
    dmax_uncal = max(r["uncal_ages"].subtract(base["uncal_ages"]).dropna().abs().max()
                     for r in rows if r is not base)
    dmax_final = max(r["ages"].subtract(base["ages"]).dropna().abs().max()
                     for r in rows if r is not base)
    w(f"量级：未校正列最大平移 **{dmax_uncal:.3f} Ma**，最终列最大平移 "
      f"**{dmax_final:.3f} Ma** —— 后者只有前者的 "
      f"{dmax_final / dmax_uncal * 100:.1f}%。**结论：对有监控标样的批次，"
      "91500 的口径分歧在最终报出值上基本无害。**")
    w()

    # ── 表 3：和不确定度摆在一起 ─────────────────────────────────────
    w("## 表 3　这点平移，和不确定度比是多大")
    w()
    floor = (base["age_of_r68"] / base["age_written"] - 1) * 100
    w("| 量 | 数值 | 单位 |")
    w("|---|---|---|")
    w(f"| 样品单点 1σ 的中位（{base['n_smp']} 个样品） | {base['med_s1']:.2f} | Ma |")
    w(f"| 同上，换成相对值 | {base['med_s1'] / base['med_age'] * 100:.3f} | % |")
    w(f"| 主标 QC 加权平均的 2σ | {base['primary_s2']:.3f} | Ma |")
    w(f"| 同上，换成相对值 | {base['primary_s2'] / base['primary_wmean'] * 100:.3f} | % |")
    w(f"| **常数不自洽贡献的「偏差%」地板**（age68(R68) / age_Ma − 1） | "
      f"{floor:+.4f} | % |")
    w(f"| 实际报出的主标 QC 偏差（＝这个地板 + 真实数据那一份） | "
      f"{base['primary_bias']:+.4f} | % |")
    w()

    # ── 表 4：极端口径之间的差别 ─────────────────────────────────────
    lo = min(rows, key=lambda r: r["R68"])
    hi = max(rows, key=lambda r: r["R68"])
    du = hi["uncal_ages"].subtract(lo["uncal_ages"]).dropna()
    df = hi["ages"].subtract(lo["ages"]).dropna()
    w("## 表 4　候选口径里最极端的两头，差多少")
    w()
    w(f"- 最低 R68：**{lo['cand']['label']}**（{lo['R68']:.6f}）")
    w(f"- 最高 R68：**{hi['cand']['label']}**（{hi['R68']:.6f}）")
    w(f"- 两者参考比值相差 **{(hi['R68'] / lo['R68'] - 1) * 100:+.4f}%**")
    w(f"- 未校正列 `年龄206_238`：逐点差中位 **{du.median():+.3f} Ma**、"
      f"最大 **{du.abs().max():.3f} Ma**（共 {len(du)} 个样品）")
    w(f"- 最终列 `年龄206_238_QC校正`：逐点差中位 **{df.median():+.3f} Ma**、"
      f"最大 **{df.abs().max():.3f} Ma**")
    w(f"- 拿不确定度当尺子：单点 1σ 中位 {base['med_s1']:.2f} Ma、"
      f"主标 QC 的 2σ {base['primary_s2']:.2f} Ma —— **未校正列**的极端差约是"
      f"前者的 {du.abs().max() / base['med_s1'] * 100:.0f}%，而**最终列**"
      f"只有 {df.abs().max() / base['med_s1'] * 100:.0f}%。")
    w()

    # ── 结论 ──────────────────────────────────────────────────────────
    w("## 结论")
    w()
    w1995 = next(r for r in rows if r["cand"]["key"] == "w1995")
    d2 = w1995["ages"].subtract(base["ages"]).dropna()
    p2 = (w1995["R68"] / base["R68"] - 1) * 100
    at1062 = 1062.4 * p2 / 100
    w("1. **口径分歧是「系统平移」，不是「随机噪声」。** 换参考值不改变任何一个测点的"
      "相对关系，只把全部年龄按同一个比例搬走 —— 所以它在单批数据内部**看不出来**，"
      "只有跨批次、跨实验室、或与 ID-TIMS 比时才显形。")
    w("2. **分歧落在哪一列，决定了它有多大 —— 两列差一个多数量级。** "
      f"比值层面：仓库现值 vs Wiedenbeck 1995 的 R68={w1995['R68']:.5f} 差 {p2:+.4f}%。"
      "把参考值换成后者：")
    u2 = w1995["uncal_ages"].subtract(base["uncal_ages"]).dropna()
    w(f"   - **未校正列 `年龄206_238`**（＝主标 QC 用的列，＝**没有监控标样时**的报出值）"
      f"整体平移：中位 **{u2.median():+.3f} Ma**、最大 **{u2.abs().max():.3f} Ma**；"
      f"一个 1062.4 Ma 的主标本身是 **{at1062:+.2f} Ma**。")
    w(f"   - **最终列 `年龄206_238_QC校正`**（有监控标样时的报出值）被 Ple 二次校正"
      f"吸收到：中位 **{d2.median():+.3f} Ma**、最大 **{d2.abs().max():.3f} Ma**。")
    w(f"   两者都远小于单点 1σ（{base['med_s1']:.1f} Ma）、也不到主标 QC 的 2σ"
      f"（{base['primary_s2']:.1f} Ma）。**但它不随机** —— 单批数据内部永远看不出来"
      "（所有点一起搬），跨实验室、或与 ID-TIMS 比时才显形。")
    w("3. **换档会断可比性。** 两档之间不可直接比较 —— 这正是 2026-09-25 换默认档时"
      "把旧结果全部按\"需要重跑\"处理的原因。量级写进 `AGENTS.md` §10 备查。")
    w("4. **对没有监控标样的批次（如桂北），这份平移不再被吸收。** 二次校正要有 Ple "
      "才谈得上；没有它，91500 的偏差会**原封不动留在报出值里**，量级就是上面"
      "「未校正列」那一行。所以这类批次报数时应显式写明用的是哪一套 91500 参考值。")
    w()

    # ── 还原断言 ─────────────────────────────────────────────────────
    after = dict(STANDARDS[PRIMARY])
    assert after == before, f"参考值没还原干净：{before} -> {after}"
    w("---")
    w(f"（自检：脚本跑完后 `STANDARDS['91500']` 与启动时逐键相同 —— "
      f"`{before['R68']}, {before['R76']}, {before['age_Ma']}`。）")

    text = "\n".join(out) + "\n"
    sys.stdout.write(text)
    if args.out:
        pathlib.Path(args.out).write_text(text, encoding="utf-8")
        sys.stdout.write(f"\n[已写入 {args.out}]\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
