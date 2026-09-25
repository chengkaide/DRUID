# -*- coding: utf-8 -*-
"""
量「一个 4 s 窗口到底有多大」—— cycle 数，以及各通道在一个窗口里的净计数。

    python tools/window_scale_facts.py

为什么要有这个脚本
------------------
文档与代码注释里有几处引用这些数字，比如：

    · 4 s 窗口含 **12~13 个 cycle**（不是"25~35 个"）
    · 窗口级 **207Pb 净计数中位约 9×10³**
    · 真正计数稀少的通道是 **204Pb**，窗口净计数中位约 37（5% 分位为负）

这些说法**曾经全是错的**。文档里长期写着"窗口级 207Pb 只有约 130 个计数"，
用来当"207Pb/206Pb 不套 F(τ)"的第二条理由 —— 但 130 这个量级属于 **204Pb**，
是**错安了通道**（207Pb 的计数从来没这么少）。它同时出现在 8 个文件里，
活了好几轮没人发现，因为**没有任何东西能把它重算一遍**。

所以这个脚本的职责就一个：让这些数字变成可重算的。`tests/check_example_batch.py`
会调用 `measure()`，与基线、与页面上的说法逐项比；那三个错串则全线禁用（见下）。

口径（必须写清，"计数"两个字不写口径可以差三倍）
------------------------------------------------
· 走工具自己的路径：`read_sequence` → `load_spot`（死时间 → 找剥蚀段 → 扣空白 →
  裁边）→ 窗口切分，与 `run_batch` 装载阶段同一套代码，不另搓一条。
· 窗口 = `BatchConfig` 的默认 `win=4.0 / step=1.0`，从配置里读，不写死。
· **"净计数" = Σ(net cps) × cycle 周期。** `Tra.net` 的单位是 **cps**，
  乘上 cycle 周期才是计数。直接对 cps 求和会得到一个没有物理含义的数
  （这正是上一轮笔记里把 206Pb 的 1.76×10⁴ 当成"207Pb 计数"记下来的由来）。
· 统计的是**全部窗口**（不是"每个测点先取中位、再在测点间取中位"）。
  两种口径差得不小，所以只保留一种，并在返回值里标明窗口总数。
· 窗口内的点数是整数（由 `剖面窗口` 表也有的 `n_cycles` 列给出），与 `fix` 参数无关；
  计数则用 `window_sums`（它也不看 `fix`）。所以这里不需要普通铅比例。

返回值
------
`measure()` 返回 dict：

    batch          示例批次目录
    win_s/step_s   窗口宽度与步长 (s)
    cycle_s        相邻 cycle 的时间差中位 (s)
    n_cycles_med   每窗口 cycle 数中位（int 化前的浮点中位）
    n_cycles_min   最少 / 最多
    n_cycles_max
    n_win          窗口总数
    n_spot         参与统计的测点数
    cnt             {质量数: {"med":…, "p05":…, "p95":…}}，单位是计数
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from druid.depth.windows import window_profile, window_sums   # noqa: E402
from druid.io import read_sequence                            # noqa: E402
from druid.io.sequence import ROLE_GLASS                      # noqa: E402
from druid.io.qtegra import read_qtegra                       # noqa: E402
from druid.reduction.trace import load_spot                   # noqa: E402
from druid.workflow import BatchConfig, spot_csv_path         # noqa: E402

BATCH = ROOT / "examples" / "EX2022A"

#: 关心的通道。207 是"不套 F(τ)"那条结论的主角，204 是真正计数稀少的那个。
MASSES = (204, 206, 207, 238)

#: **全仓禁用**的旧说法（纯子串）。它们要么数字错、要么错安了通道，
#: 一旦有人从旧文档/旧笔记里抄回来，`check_example_batch.py` 当场就会红。
STALE_CLAIMS = (
    "约 130 个计数",
    "百来个计数",
    "25~35 个 cycle",
    "25–35 个采样周期",
    "25-35 个采样周期",
)

#: 同一批谬误的**正则**形式。为什么光有子串不够 —— 本轮实测吃过一次亏：
#: "25~35 个 cycle" 在 `docs/from-first-principles.html` 里其实写作
#: "大约 25–35 个 **采样周期**"（连字符 + 另一个词），按子串搜**完全搜不到**，
#: 差点漏掉。子串只能挡住"照抄原文"那一种写法，换个说法就溜过去了。
STALE_PATTERNS = (
    # 每窗口 25~35 个 cycle / 采样周期（真值 12~13）
    r"25\s*[~\-–—～]\s*35\s*(?:个)?\s*(?:cycle|cycles|采样周期)",
    # "207Pb 只有约 130 个计数" 及其变体
    r"207\s*Pb[^。；\n]{0,14}?(?:只有|仅约|仅|才)\s*(?:约)?\s*130\s*个计数",
    # "207Pb 的计数本来就很少"——204 扣铅伤的是 207Pb/206Pb 的**放射成因份额**，
    # 不是计数（见 druid/__init__.py 第 1 条）
    r"207\s*Pb[^。；\n]{0,12}?计数(?:极少|很少|本来就少|本来就很少)",
)


def measure(batch: pathlib.Path | None = None) -> dict:
    """量一遍。见模块 docstring 的"返回值"。"""
    cfg = BatchConfig(data_dir=str(batch or BATCH))
    seq = read_sequence(cfg.list_file, cfg.primary, cfg.secondary)

    # ── cycle 周期 ────────────────────────────────────────────────────────
    # 窗口切分只认"行数"，周期必须单独量：原始 CSV 相邻行的时间差。
    # （示例批次 85 个文件的周期完全一致，所以取全局中位是安全的；
    #   换批次若不一致，这里的中位就只是一个代表性数字 —— 这也是为什么
    #   文档里该说"本批实测"。）
    dts = []
    for name in seq["file"]:
        path = spot_csv_path(cfg.data_dir, name)
        if not path.exists():
            continue
        _, t, _, _ = read_qtegra(path)
        if len(t) > 3:
            dts.append(np.diff(t))
    if not dts:
        raise RuntimeError(f"{cfg.data_dir} 里读不到任何时间轴")
    cycle_s = float(np.median(np.concatenate(dts)))

    # ── 逐测点：窗口点数 + 各通道净计数 ──────────────────────────────────
    ncyc, cnt = [], {m: [] for m in MASSES}
    n_spot = 0
    for _, row in seq.iterrows():
        path = spot_csv_path(cfg.data_dir, row["file"])
        if not path.exists() or row["role"] == ROLE_GLASS:
            continue                       # 与 load_batch 一致：玻璃标样不参与
        try:
            tr = load_spot(int(row["order"]) - 1, path, row["sample"], row["role"],
                           trim=cfg.trim, deadtime_ns=cfg.deadtime_ns,
                           blank_dur=cfg.blank_dur)
        except Exception:                  # noqa: BLE001  单点失败不影响整批
            continue

        prof = window_profile(tr, win=cfg.win, step=cfg.step)
        if prof is None or not len(prof):
            continue
        ncyc.append(prof["n_cycles"].to_numpy(float))

        _, S = window_sums(tr, cfg.win, cfg.step, masses=MASSES)
        for m in MASSES:
            # net 是 cps → 乘周期才是计数
            cnt[m].append(np.asarray(S[m], float) * cycle_s)
        n_spot += 1

    if not ncyc:
        raise RuntimeError(f"{cfg.data_dir} 里没有任何测点切出了窗口")

    nc = np.concatenate(ncyc)
    out = {
        "batch": str(pathlib.Path(cfg.data_dir).name),
        "win_s": float(cfg.win),
        "step_s": float(cfg.step),
        "cycle_s": cycle_s,
        "n_cycles_med": float(np.median(nc)),
        "n_cycles_min": int(nc.min()),
        "n_cycles_max": int(nc.max()),
        "n_win": int(nc.size),
        "n_spot": n_spot,
        "cnt": {},
    }
    for m in MASSES:
        v = np.concatenate(cnt[m])
        out["cnt"][m] = {"med": float(np.median(v)),
                         "p05": float(np.quantile(v, 0.05)),
                         "p95": float(np.quantile(v, 0.95))}
    return out


def main() -> int:
    from druid.console import ensure_utf8_streams
    ensure_utf8_streams()

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--batch", type=pathlib.Path, default=BATCH,
                    help=f"批次目录（默认 {BATCH}）")
    args = ap.parse_args()

    d = measure(args.batch)
    print(f"批次 {d['batch']}  ·  {d['n_spot']} 个测点 / {d['n_win']} 个窗口")
    print(f"窗口 {d['win_s']:g} s / 步长 {d['step_s']:g} s")
    print("-" * 60)
    print(f"cycle 周期            {d['cycle_s']:.5f} s")
    print(f"每窗口 cycle 数       中位 {d['n_cycles_med']:.0f}"
          f"（{d['n_cycles_min']} ~ {d['n_cycles_max']}）")
    for m in MASSES:
        c = d["cnt"][m]
        print(f"{m:>3}Pb 窗口净计数     中位 {c['med']:>10,.0f}"
              f"  5% {c['p05']:>10,.0f}  95% {c['p95']:>12,.0f}")
    print()
    print("（口径：Σ(net cps) × cycle 周期 —— 对 cps 直接求和得到的不是计数。）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
