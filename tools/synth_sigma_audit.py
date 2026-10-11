#!/usr/bin/env python
"""
tools/synth_sigma_audit.py —— 量一件事：**流水线报出的域级 1σ 到底被低估了多少倍**（只读）
=========================================================================================

为什么需要它
------------
`结果表` / `深度剖面域` 里那些 1σ 是用"独立观测"的式子算的：窗口之间共用计数、
且校准漂移/分馏失配/普通铅这些**系统性**项并不在计数统计里。于是报出的误差棒
偏窄 —— 而这是**真数据上无法自查**的：真值不知道。只能造合成数据。

本脚本用仓库自己的合成机制（`tools/synth_profile.py`，一个数都不另算）：
造出一批**真值已知**的两域剖面 → 跑完整 `run_batch` → 取每个测点解出的两个
年龄域 (age, σ)，与喂进去的真值配对，算

    pull = (报回年龄 − 真值) / 报出的 1σ

σ 诚实 ⇒ pull 的标准差 ≈ 1、|pull| ≤ 1 的比例 ≈ 68%、median|pull| ≈ 0.6745。
实测（2026-10-11，条件⑤ 全套效应）得到 **median|pull|/0.6745 ≈ 1.94**。

    ★ `SIGMA_INFLATION` 就是那个倍数，**本文件是它唯一的出处** ——
      `tools/gen_docs_segmented.py` 直接 import 它去校正页面上的误差棒与 z。

两档对照（同一次运行）：
    · 条件① 「只有计数统计」 → pull 标准差约 1.2 ⇒ **计数统计这一层是诚实的**
    · 条件⑤ 「全套效应」     → pull 标准差约 1.8 ⇒ 差额来自额外效应

用法
----
    python tools/synth_sigma_audit.py            # 重新量一遍并打印
    python tools/synth_sigma_audit.py --check    # 与 SIGMA_INFLATION 比对，不一致非零退出
    python tools/synth_sigma_audit.py --seeds 12 # 加大样本（默认 6）

⚠ 只读：临时批次写到 G 盘临时目录、跑完自删（`synth_profile.fast_rmtree`）。
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "synth_profile", Path(__file__).resolve().parent / "synth_profile.py")
sp = importlib.util.module_from_spec(_spec)
sys.modules["synth_profile"] = sp
_spec.loader.exec_module(sp)

from druid.core.constants import DEFAULT_REF_PRESET          # noqa: E402

# ★ 唯一的出处：域级 1σ 应乘的倍数（由本脚本在条件⑤下实测）
#   三个互相独立的估计都给同一个量级：pull 标准差 1.81、median|pull|/0.6745 1.90、
#   err_med/(sig_med·0.6745) 1.94 ⇒ 取 1.9 作为校正因子。
SIGMA_INFLATION = 1.9
SIGMA_INFLATION_TOL = 0.35          # --check 的容差（N 小时 σ 本身就有 ±20~25% 的抖动）

# 参考运行的完整读数（2026-10-11，--seeds 6；`--check` 会连这些一起比对）
MEASURED = {
    "seeds": 6,
    "date": "2026-10-11",
    "counting": dict(name="① 只有计数统计", n=38, n_two=19, n_all=36,
                     pull_std=1.19, cover1=0.632, cover2=0.921, med_abs=0.68,
                     sig_med=2.67, err_med=1.73),
    "full": dict(name="⑤ 全套（＋普通铅 5%）", n=36, n_two=18, n_all=36,
                 pull_std=1.81, cover1=0.417, cover2=0.778, med_abs=1.28,
                 sig_med=2.68, err_med=3.51),
}


def _pulls(det, cond_name):
    """抽出某一条件下、判出两域的测点，逐个域算 pull = (报回 − 真值)/σ。"""
    sub = det[det["条件"] == cond_name]
    two = sub[sub["ndet"] == 2]
    pulls, errs, sigs = [], [], []
    for _, r in two.iterrows():
        ages, ses = list(r["ages"]), list(r["ses"])
        if len(ages) != 2 or len(ses) != 2:
            continue
        truth = sorted([float(r["top"]), float(r["top"]) - float(r["dt"])])
        for (a, s), t in zip(sorted(zip(ages, ses)), truth):
            pulls.append((a - t) / s)
            errs.append(a - t)
            sigs.append(s)
    return (np.array(pulls, float), np.array(errs, float), np.array(sigs, float),
            len(two), len(sub))


def audit(seeds=6, verbose=False):
    """跑两档条件，返回 dict：inflation / pull_std / cover1 / ...（全部实测）。"""
    conds = sp.base_conditions()
    pair = [conds[0], conds[4]]            # ① 只有计数统计 / ⑤ 全套效应
    out = Path(r"G:/_zd_run/synth_sigma") if Path("G:/").exists() \
        else Path.home() / "druid_synth_sigma"
    det = sp.run_all(dts=[20.0, 30.0], conds=pair, seeds=int(seeds), repeat=3,
                     tops=[450.0], betas=[0.5], u_cps=1.0e6, th_u=0.5,
                     ref_preset=DEFAULT_REF_PRESET, win=4.0, step=1.0,
                     out_root=out, keep=False, verbose=bool(verbose))
    res = {}
    for c in pair:
        pulls, errs, sigs, n2, n_all = _pulls(det, c["name"])
        res[c["name"]] = dict(
            n=int(len(pulls)), n_two=int(n2), n_all=int(n_all),
            pull_std=float(pulls.std(ddof=1)) if len(pulls) > 1 else float("nan"),
            cover1=float(np.mean(np.abs(pulls) <= 1)) if len(pulls) else float("nan"),
            cover2=float(np.mean(np.abs(pulls) <= 2)) if len(pulls) else float("nan"),
            med_abs=float(np.median(np.abs(pulls))) if len(pulls) else float("nan"),
            sig_med=float(np.median(sigs)) if len(sigs) else float("nan"),
            err_med=float(np.median(np.abs(errs))) if len(errs) else float("nan"),
        )
    full = res[pair[1]["name"]]
    res["inflation"] = (full["med_abs"] / 0.6745) if np.isfinite(full["med_abs"]) \
        else float("nan")
    return res


def _print(res):
    for name, r in res.items():
        if name == "inflation":
            continue
        print(f"\n条件：{name}")
        if not np.isfinite(r["pull_std"]):
            print("  没有可统计的域")
            continue
        print(f"  判出两域 {r['n_two']}/{r['n_all']} 个，可统计的域 {r['n']} 个")
        print(f"  |pull| ≤ 1 的比例 {r['cover1'] * 100:.1f}%（σ 诚实应 ≈ 68%）"
              f"；≤ 2 的比例 {r['cover2'] * 100:.1f}%（应 ≈ 95%）")
        print(f"  pull 标准差 {r['pull_std']:.2f}（应 ≈ 1）；"
              f"median|pull| {r['med_abs']:.2f}（应 ≈ 0.6745）")
        print(f"  报出 1σ 中位 {r['sig_med']:.2f} Ma；实际 |误差| 中位 {r['err_med']:.2f} Ma")
    print(f"\n★ 域级 1σ 的乐观倍数（条件⑤）：median|pull|/0.6745 = "
          f"{res['inflation']:.2f}")
    print(f"  本文件记的 SIGMA_INFLATION = {SIGMA_INFLATION}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="量域级 1σ 的乐观倍数（只读）。")
    ap.add_argument("--seeds", type=int, default=6, help="每个条件跑几个种子（默认 6）")
    ap.add_argument("--check", action="store_true",
                    help=f"与 SIGMA_INFLATION 比对（容差 ±{SIGMA_INFLATION_TOL}）")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    res = audit(seeds=args.seeds, verbose=args.verbose)
    _print(res)
    if args.check:
        bad = []
        k = res["inflation"]
        if not np.isfinite(k) or abs(k - SIGMA_INFLATION) > SIGMA_INFLATION_TOL:
            bad.append(f"乐观倍数 实测 {k:.2f} vs 记录 {SIGMA_INFLATION}")
        for key, ref in (("counting", MEASURED["counting"]), ("full", MEASURED["full"])):
            got = res[ref["name"]]
            for f in ("pull_std", "cover1", "med_abs"):
                if abs(got[f] - ref[f]) > 0.25:
                    bad.append(f"{key}.{f} 实测 {got[f]:.3f} vs 记录 {ref[f]:.3f}")
            if got["n"] != ref["n"]:
                bad.append(f"{key}.n 实测 {got['n']} vs 记录 {ref['n']}")
        if bad:
            print("\n✗ 与记录值不一致：")
            for b in bad:
                print("   ", b)
            print("   要么改常数（并说明为什么），要么查流水线。")
            return 1
        print(f"\n✓ 全部读数与记录值一致（乐观倍数 {k:.2f}，"
              f"记 {SIGMA_INFLATION} ±{SIGMA_INFLATION_TOL}）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
