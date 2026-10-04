# -*- coding: utf-8 -*-
"""
自实现的数值 vs scipy —— 把"我们不依赖 scipy"这件事变成**可核验**的。

    python tools/verify_against_scipy.py

为什么要有这个脚本
------------------
本包**刻意不依赖 scipy**（`AGENTS.md` §1：只有 numpy/pandas/matplotlib/openpyxl/xlrd）。
理由是"不依赖会变的第三方行为"—— 但这句话本身是个**承诺**，不是证据。
承诺要被人相信，就得能被别人自己验一遍。这个脚本就是那份证据。

它量三件事
----------
1. **卡方分布上尾概率**（`core.statistics.chi2_sf`）vs `scipy.stats.chi2.sf`
   —— 这是全包唯一一处"自己实现统计分布"的地方，也是最该被质疑的。
2. **年龄 ↔ 比值**（`age68` / `age76` / `r68_of_age` / `r76_of_age`）
   vs 闭式解与往返自洽。
3. **Stacey–Kramers 普通铅**（`core.common_lead.stacey_kramers`）
   vs 按教科书公式**独立重写**一遍的实现（不复用仓库任何常数）。

⚠ 装 scipy 不污染生产环境
------------------------
`upb` 环境**不许**装 scipy —— 装了就无法再验证"没有 scipy 也能跑"这条。
所以这个脚本要在**另一个环境**里跑：

    "C:/Users/<用户名>/.workbuddy/binaries/python/versions/3.13.12/python.exe" \
        -m venv G:/_druid_check/venv_scipy
    "G:/_druid_check/venv_scipy/Scripts/python.exe" -m pip install scipy
    "G:/_druid_check/venv_scipy/Scripts/python.exe" tools/verify_against_scipy.py

⚠ 两边 numpy 版本必须相同，否则"一致"只能说明是 numpy 变了。
   实测本仓库的结论建立在 **numpy 2.5.3 / scipy 1.18.1** 上。

口径（不写清就没法复现）
------------------------
· `chi2_sf` **只收标量**（内部 `float(x)`），所以逐点对比而不是整数组。
  这是接口事实，不是缺陷 —— 调用点都是标量。
· `stacey_kramers` 的参考实现**独立重写**：初始比 9.307 / 10.294 / 29.476、
  两阶段 ²³⁸U/²⁰⁴Pb = 7.19 / 9.74、Th/U κ=4.0、分界 3700 Ma（S&K 1975）。
  衰变常数取 **λ₂₃₂ = 4.9475e-11 /a**（²³²Th 半衰期 14.05 Gy ⇒ ln2/14.05e9
  = 4.93e-11）—— ⚠ 若这里写成 1.16e-10 会得到 208/204 差 30% 的假警报，
  那是**参考实现写错**，不是仓库错。第一版就踩过这个坑。
· `overlap_factor` 没有 scipy 对应实现，验的是「**趋近** 4」而不是「等于 4」：
  K→4 是 n→∞ 的极限，有限 n 的精确值是 3.95（n=100）、3.995（n=1000）。
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser(description="自实现数值 vs scipy 的逐项对比")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="相对误差容差（默认 1e-9，即双精度的 ~7 位有效）")
    args = ap.parse_args()

    try:
        import scipy
        from scipy import stats
    except ImportError:
        print("这个脚本需要在**装了 scipy 的隔离环境**里跑。\n"
              "生产环境（upb）刻意不装 scipy —— 装了就无法验证"
              "『没有 scipy 也能跑』。\n\n"
              "  python -m venv G:/_druid_check/venv_scipy\n"
              "  G:/_druid_check/venv_scipy/Scripts/python.exe -m pip install scipy\n"
              "  G:/_druid_check/venv_scipy/Scripts/python.exe "
              "tools/verify_against_scipy.py", file=sys.stderr)
        return 2

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from druid.core.common_lead import stacey_kramers
    from druid.core.geochronology import age68, age76, r68_of_age, r76_of_age
    from druid.core.statistics import chi2_sf, overlap_factor, weighted_mean

    print(f"环境：numpy {np.__version__}  scipy {scipy.__version__}")
    # numpy 版本必须与生产环境一致，否则"一致"只能说明是 numpy 变了。
    # 生产解释器的路径**不写死**（写死既含本机用户名、换台机器也必然失败）：
    # 用环境变量 DRUID_PY 传入，不设就跳过这项检查。
    prod = os.environ.get("DRUID_PY", "")
    if prod and Path(prod).exists():
        r = subprocess.run([prod, "-c", "import numpy;print(numpy.__version__)"],
                           capture_output=True, text=True)
        pver = r.stdout.strip()
        same = (pver == np.__version__)
        print(f"生产环境 numpy {pver or '(读取失败)'} —— "
              f"{'一致' if same else '★不一致，两边对比无效'}")
        if not same:
            print("⚠ numpy 版本不同，『一致』只能说明是 numpy 变了。"
                  "请在同一版本的 numpy 下重跑。")
            return 1
    else:
        print("（未设 DRUID_PY，跳过 numpy 版本核对；"
              "设上它才能确认两边 numpy 同版本）")
    print()

    over = []

    def cmp(name, ours, ref, unit="", tol=None):
        tol = args.tol if tol is None else tol
        ours = np.asarray(ours, float)
        ref = np.asarray(ref, float)
        rel = np.abs(ours - ref) / np.maximum(np.abs(ref), 1e-300)
        i = int(np.nanargmax(rel))
        worst = float(np.nanmax(rel))
        bad = worst > tol
        if bad:
            over.append((name, worst, tol))
        print(f"  {name:36s} 最大相对误差 {worst:.3e}  @ {unit}[{i}]"
              + (f"  <- 超过 {tol:.0e}" if bad else ""))

    # ── 1. 卡方上尾概率 ──
    print("【1】chi2_sf vs scipy.stats.chi2.sf（全包唯一自己实现的统计分布）")
    xs = np.array([1e-3, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 100.0, 1074.7])
    for df in (1, 2, 3, 5, 10, 30, 100, 1000):
        cmp(f"chi2_sf(x, df={df})",
            np.array([chi2_sf(v, df) for v in xs]), stats.chi2.sf(xs, df), "x=")
    xe = np.array([100.0, 500.0, 2000.0, 1e4, 1e6])
    cmp("chi2_sf 极端尾部 df=10",
        np.array([chi2_sf(v, 10) for v in xe]), stats.chi2.sf(xe, 10), "x=")
    xd = np.array([1e4, 1e5, 1e6])
    cmp("chi2_sf 大 df=1e5 极端",
        np.array([chi2_sf(v, 1e5) for v in xd]), stats.chi2.sf(xd, 1e5), "x=", tol=1e-4)

    # ── 2. 年龄 ↔ 比值 ──
    print("\n【2】年龄 ↔ 比值（可对照的独立参考只有 age68 的闭式解）")
    ages = np.array([1.0, 10.0, 100.0, 250.0, 460.0, 1000.0, 1500.0, 2500.0, 4000.0])
    r68 = np.array([0.005, 0.02, 0.05, 0.1, 0.17917, 0.3, 0.5, 0.8])
    lam238 = 1.55125e-10                                    # 1/a
    cmp("age68 vs 闭式解 ln(1+R)/lambda",
        age68(r68), np.log(1.0 + r68) / lam238 / 1e6, "R68=")
    cmp("age68 往返 age→R→age", age68(r68_of_age(ages)), ages, "age=")
    a76 = np.linspace(1, 4570, 40)
    cmp("age76 往返 全域 1–4570 Ma", age76(r76_of_age(a76)), a76, "age=")
    r76 = np.array([0.0461, 0.05, 0.0749, 0.1, 0.2, 0.5, 0.8456])
    cmp("age76 → R76 往返", r76_of_age(age76(r76)), r76, "R76=")

    # ── 3. Stacey–Kramers vs 独立重写 ──
    print("\n【3】stacey_kramers vs 教科书公式的独立实现")

    def sk_ref(t_ma):
        """S&K (1975) 两步模型，单位 Ma。**不复用仓库任何常数**。"""
        l238, l235, l232 = 1.55125e-10, 9.84850e-10, 4.94750e-11
        t = np.clip(np.asarray(t_ma, float), 0.0, 4570.0) * 1e6
        t0, t1, mu1, mu2, kap, u = 4570e6, 3700e6, 7.19, 9.74, 4.0, 137.818
        p6 = 9.307 + mu1 * (np.exp(l238 * t0) - np.exp(l238 * t1)) \
            + mu2 * (np.exp(l238 * t1) - np.exp(l238 * t))
        p7 = 10.294 + (mu1 / u) * (np.exp(l235 * t0) - np.exp(l235 * t1)) \
            + (mu2 / u) * (np.exp(l235 * t1) - np.exp(l235 * t))
        p8 = 29.476 + mu1 * kap * (np.exp(l232 * t0) - np.exp(l232 * t1)) \
            + mu2 * kap * (np.exp(l232 * t1) - np.exp(l232 * t))
        return np.stack([p6, p7, p8], axis=-1)

    ts = np.array([0.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0, 4570.0])
    ours = np.array([stacey_kramers(t) for t in ts])
    ref = sk_ref(ts)
    for i, nm in enumerate(("206/204", "207/204", "208/204")):
        cmp(f"stacey_kramers {nm}", ours[:, i], ref[:, i], "t=")

    # ── 4. MSWD → 概率 ──
    print("\n【4】weighted_mean 的 MSWD → 概率（自实现 chi2_sf vs scipy）")
    rng = np.random.default_rng(20261004)
    for n in (5, 10, 30, 100):
        x = rng.normal(450.0, 3.0, n)
        _mu, _se, mswd, k = weighted_mean(x, np.full(n, 3.0))
        obs = mswd * (k - 1)
        o, r = chi2_sf(obs, k - 1), stats.chi2.sf(obs, k - 1)
        print(f"  n={k:3d}  MSWD {mswd:.6f}  χ² {obs:.6f}  "
              f"自实现 p={o:.12e}  scipy p={r:.12e}  绝对差 {abs(o - r):.3e}")

    # ── 5. overlap_factor ──
    print("\n【5】overlap_factor（几何因子，scipy 无对应实现）")
    print("    K→4 是 n→∞ 的**极限**；有限 n 精确值略小，所以只验「趋近」")
    for n in (10, 100, 1000, 10000):
        k = overlap_factor(n, 0.25)
        print(f"  n={n:6d}  K = {k:.12f}  距 4 的差 {4.0 - k:.3e}")
    for stp in (0.0, 1.0):
        k = overlap_factor(100, stp)
        ok = "OK" if abs(k - 1) < 1e-12 else "MISMATCH"
        print(f"  step/win={stp:g}  K = {k:.12f}（应为 1）{ok}")

    # ── 6. 极端 df 下还剩多少有效位 ──
    print("\n【6】chi2_sf 在极端 df 下还剩多少有效位（相对 scipy）")
    for df in (1, 10, 100, 1e3, 1e4, 1e5):
        x = df * 1.5
        o, r = chi2_sf(x, df), stats.chi2.sf(x, df)
        if r == 0 or o == 0:
            print(f"  df={df:>8.0f}  p={o:.3e}  （两侧都下溢到 0，无法比较）")
            continue
        rel = abs(o - r) / r
        print(f"  df={df:>8.0f}  p={o:.6e}  相对误差 {rel:.3e}"
              f"  ~ {-np.log10(rel):.1f} 位有效")

    print("\n" + "=" * 68)
    if over:
        print("超过容差的项：")
        for nm, w, tol in over:
            print(f"  {nm}: {w:.3e}（容差 {tol:.0e}）")
        return 1
    print("★ 全部对比项在容差内 —— 自实现与 scipy / 闭式解一致")
    print("  结论：「不依赖 scipy」是可核验的，不是口号。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
