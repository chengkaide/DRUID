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

改了算法就有数字变化是正常的 —— 那时候要**重新确认基线并更新这里的常量**，
而不是把容差放宽。容差放宽等于把这个检查废掉。
"""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BATCH = ROOT / "examples" / "EX2022A"

# ── 基线（2026-09-19 由 examples/EX2022A 实跑得到，与真实批次逐格一致）──────
ROWS = {"结果": 83, "标样QC": 2, "深度剖面域": 91, "剖面窗口": 1677}
ROLES = {"样品": 48, "主标": 21, "监控标样": 14}
QC = {                       # 标样 -> (加权平均年龄 Ma, 偏差 %)
    "91500": (1059.7543, -0.2490),
    "Ple": (343.2381, 1.8118),
}
AGE = {"中位": 458.0686, "5%": 420.1799, "95%": 643.5333}   # 只统计样品
CONC = {"中位": 101.3006, "90-110%占比": 0.9375}
STRUCT = {None: 35, "多域(2)": 18, "均一": 15, "多域(3)": 15}

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
               "剖面窗口": len(win)}[name]
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
