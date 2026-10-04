# -*- coding: utf-8 -*-
"""
一个批次里到底有哪几个通道？哪些进了包、哪些在解析阶段就被丢掉了。

    python tools/trace_channels.py [--batch 目录]

为什么要有这个脚本
------------------
`druid/io/qtegra.py` 的 `MASS_RE` 只匹配 **8 个质量数**
（202/204/206/207/208/232/235/238）。其余通道 —— 29Si、49Ti、89Y、91Zr、93Nb、
**全部 14 个 REE**、178Hf、181Ta —— 在**解析阶段**就不进 `data` 字典了。

这不是疏漏，是设计：它们不参与 U-Pb 归一化，进来只会占内存。
但后果必须说清楚：**任何"用化学给年龄当旁证"的做法**
（Lim et al. 2024 的六步滤波、Ce/Ce*、(Sm−Nd)N、Th/U 台阶、Y/Ho…）
第一步都要先跨过这个事实 —— 知道丢了什么、还拿不拿得回来。

`read_qtegra()` 顺手把**原始列名**放在返回的 `columns` 里，所以不用另写解析器。
这个脚本把两件事量出来：

    1. 每个通道出现在**多少个测点**的表头里（覆盖度）；
    2. 它有没有进 `data`（即是否被 `MASS_RE` 留下）。

另一个数字也一并报出：表头里**有没有 235U**。它决定了 `207Pb/206Pb` 能不能
直接测出来 —— 本示例批次没有它，所以 R76 只能靠 `238U/235U = 137.818`
这个假定换算（这条结论在 `AGENTS.md` §7 有记录）。

口径（不写清就能差好几倍）
--------------------------
· 走 `read_qtegra()`，与装载阶段**同一个**解析器 —— 它按内容找 `Time` 表头，
  不硬编码行号，所以不同批次的元数据行数不同也不影响。不另搓一条 CSV 读取。
· 覆盖度按"该列名出现在该文件的表头里"计，**不看数值是否有效** ——
  这里问的是"采集时有没有这一路"，不是"这一路测得准不准"。
· 序列表（`*_LIST.csv`）不是采集文件，**必须排除**；判据取自
  `BatchConfig.list_file`（DRUID 按目录名找它），不靠文件名猜。
· 玻璃 / 标样 / 样品**都统计**，并分别报覆盖度：它们来自同一套采集模板，
  少了任一类都只是少一层旁证，没有理由替读者过滤。
· **`data` 里的 key 才是"进包"**，而不是"名字里有没有 8 个质量数之一" ——
  两者在本批恰好一致，但前者才是定义（`read_qtegra` 用 `setdefault` 只保留
  第一次读到的同名质量）。

返回值
------
`measure()` 返回 dict：

    batch           批次目录名
    n_file          参与统计的采集文件数
    list_file       识别出来的序列表文件名（已排除）
    n_col           表头里出现过的不同列名个数（含第 0 列的时间轴）
    n_channel       真正的通道数 = n_col − 1（时间轴不算通道）
    channels        [{name, coverage, n_file, mass, is_axis, kept}]，按表头顺序
    kept_names      进包的通道名（有序）
    dropped_names   没进包的通道名（有序）
    kept_masses     进包的质量数（升序）
    has_235         表头里有没有 235U
"""
from __future__ import annotations

import argparse
import collections
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from druid.io.qtegra import MASS_RE, read_qtegra          # noqa: E402
from druid.workflow import BatchConfig                    # noqa: E402

BATCH = ROOT / "examples" / "EX2022A"


def measure(batch: pathlib.Path | None = None) -> dict:
    """量一遍。见模块 docstring 的"返回值"。"""
    cfg = BatchConfig(data_dir=str(batch or BATCH))
    data_dir = pathlib.Path(cfg.data_dir)
    if not data_dir.is_dir():
        raise FileNotFoundError(f"批次目录不存在: {data_dir}")

    # ── 序列表必须排除：它不是采集文件，表头结构与数据文件不同 ──
    list_name = pathlib.Path(cfg.list_file).name
    files = [p for p in sorted(data_dir.glob("*.csv")) if p.name != list_name]
    if not files:
        raise RuntimeError(f"{data_dir} 里没有任何采集 CSV"
                           f"（已排除序列表 {list_name}）")

    cov: collections.Counter = collections.Counter()
    order: list[str] = []                  # 表头出现顺序（第一次出现的次序）
    kept: set[int] = set()

    for path in files:
        try:
            _, _, data, columns, _info = read_qtegra(path)
        except ValueError:
            # 单个文件解析失败不该让整批统计垮掉 —— 与 load_batch 的宽容度一致
            continue
        kept |= set(data.keys())
        for c in columns:
            if c not in cov:
                order.append(c)
            cov[c] += 1

    channels = []
    for i, c in enumerate(order):
        m = MASS_RE.search(c)
        mass = int(m.group(1)) if m else None
        # 第 0 列永远是时间轴（read_qtegra 的约定），它不是"通道" ——
        # 把它算进"被丢掉"会让丢失数从 21 变成 22，是个纯口径错误。
        is_axis = (i == 0)
        channels.append({
            "name": c,
            "coverage": cov[c],
            "n_file": len(files),
            "mass": mass,
            "is_axis": is_axis,
            # "进包"由 data 的 key 定义，不是由名字里有没有 8 个质量数之一
            "kept": (not is_axis) and mass is not None and mass in kept,
        })

    chans = [c for c in channels if not c["is_axis"]]
    return {
        "batch": data_dir.name,
        "n_file": len(files),
        "list_file": list_name,
        "n_col": len(order),
        "n_channel": len(chans),
        "channels": channels,
        "kept_names": [c["name"] for c in chans if c["kept"]],
        "dropped_names": [c["name"] for c in chans if not c["kept"]],
        "kept_masses": sorted(kept),
        "has_235": any(c["mass"] == 235 for c in channels),
    }


def main(argv: list[str] | None = None) -> int:
    from druid.console import ensure_utf8_streams
    ensure_utf8_streams()

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--batch", type=pathlib.Path, default=BATCH,
                    help=f"批次目录（默认 {BATCH}）")
    args = ap.parse_args(argv)

    d = measure(args.batch)
    n_kept, n_drop = len(d["kept_names"]), len(d["dropped_names"])
    print(f"批次 {d['batch']}  ·  {d['n_file']} 个采集文件"
          f"（已排除序列表 {d['list_file']}）")
    print(f"表头 {d['n_col']} 列 = 1 个时间轴 + {d['n_channel']} 个通道"
          f" → 进包 {n_kept} / 丢掉 {n_drop}")
    print(f"进包的质量数 {d['kept_masses']}"
          f"   表头里有 235U 吗：{'有' if d['has_235'] else '没有'}")
    print("-" * 62)
    print(f"{'通道':<10}{'覆盖':>6}{'占比':>8}   进包?")
    for c in d["channels"]:
        pct = 100.0 * c["coverage"] / d["n_file"]
        tag = "时间轴" if c["is_axis"] else ("进包" if c["kept"] else "丢掉")
        print(f"{c['name']:<10}{c['coverage']:>6}{pct:>7.1f}%   {tag}")
    print("-" * 62)
    print("被丢掉的通道：" + "、".join(d["dropped_names"]))
    print()
    print("（口径：进不进包由 read_qtegra 返回的 data 的 key 定义，"
          "不看名字里有没有 8 个质量数之一。）")
    print("（想拿回这些通道：`MASS_RE` 是唯一的关口 —— 它匹配 8 个质量数，"
          "本批表头命中 7 个；")
    print("  本项目的库侧（汇总库/建库.py 的 MASS）已全取，包侧仍是这一份。）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
