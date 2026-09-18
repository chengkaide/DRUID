"""
druid.cli.reduce_batch —— 批处理命令行入口
========================================

用法
----
    # 进入 druid 包的**上级**目录（即 ...\\UPb处理）后执行
    python -m druid.cli.reduce_batch ^
        --dir "G:\\1.云龙锡矿\\云龙锡矿锆石\\20220301CKDB" ^
        --out "G:\\1.云龙锡矿\\云龙锡矿锆石\\UPb处理\\结果\\20220301CKDB_U-Pb结果.xlsx" ^
        --plot

    # 做对照实验：换一种整段比值方法 / 施加死时间校正
    python -m druid.cli.reduce_batch --dir ... --bulk ftau
    python -m druid.cli.reduce_batch --dir ... --deadtime-ns 14.7645

运行 -h 可查看全部参数。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# 支持两种调用方式：
#   ① python -m druid.cli.reduce_batch（包上下文完整，包根目录已在 sys.path 上）
#   ② python druid/cli/reduce_batch.py（直接跑脚本，需要自己把包根目录塞进 sys.path）
#
# 判据用 __package__：只有"直接跑脚本"时它才是空/None。
# 这里刻意不用 try/except ImportError —— 那样会把 numpy 缺失之类的真实错误
# 一并吞掉，然后在另一条路径上重新抛出，报错位置与真正的原因对不上；
# 而且导入清单要写两遍，改一处漏一处。
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from druid import __version__
from druid.io.report import export_batch
from druid.workflow import BatchConfig, run_batch


def build_parser() -> argparse.ArgumentParser:
    """组装参数解析器。**每一个可调参数都要写清楚它影响哪一步**。"""
    ap = argparse.ArgumentParser(
        prog="druid.cli.reduce_batch",
        description="LA-ICP-MS 锆石 U-Pb 批处理：Qtegra CSV → 年龄表 + 深度剖面图",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # ── 输入输出 ──
    ap.add_argument("--dir", required=True,
                    help="批次目录，内含 XXXX_LIST.xls 与 XXXX_N.csv")
    ap.add_argument("--out", default=None,
                    help="结果 Excel 路径（默认写到 <批次目录>/<批次名>_U-Pb结果.xlsx）")
    ap.add_argument("--plot", action="store_true",
                    help="生成逐点深度剖面 PNG 并汇总成 PDF")
    ap.add_argument("--plot-dir", default=None, help="图件输出目录")
    ap.add_argument("--quiet", action="store_true", help="不打印中间过程")

    # ── 标样 ──
    ap.add_argument("--primary", default="91500", help="主标（用于归一化）")
    ap.add_argument("--secondary", default="Ple", help="监控标样（用于 QC 与 σext）")

    # ── 单点还原 ──
    ap.add_argument("--trim", type=float, default=1.5,
                    help="剥蚀段两端各裁掉多少秒（避开开关激光瞬态）")
    ap.add_argument("--blank-dur", type=float, default=15.0,
                    help="气体空白取样时长 (s)")
    ap.add_argument("--n-sigma-common-pb", type=float, default=2.0,
                    help="204Pb 显著性门槛，越大越不容易触发普通铅校正")
    ap.add_argument("--deadtime-ns", type=float, default=0.0,
                    help="探测器死时间 (ns)，0 表示不校正（默认）。"
                         "本实验室实测值记在 druid/deadtime_ns.txt，"
                         "要用请显式传 --deadtime-ns 14.7645")

    # ── 深度剖面 ──
    ap.add_argument("--win", type=float, default=4.0, help="滑动窗口宽度 (s)")
    ap.add_argument("--step", type=float, default=1.0, help="滑动步长 (s)")
    ap.add_argument("--no-depth", action="store_true", help="跳过深度剖面年龄域判别")

    # ── 整段比值方法 ──
    ap.add_argument("--bulk", choices=["simple", "ftau"], default="simple",
                    help="simple=整段单一校正因子（默认，本批次更准）；"
                         "ftau=逐窗口 F(τ) 深度校正后按计数合成")

    # ── 外部重现性（留 None 则由监控标样自动估计）──
    ap.add_argument("--sigma-ext68", type=float, default=None,
                    help="强制指定 206Pb/238U 外部重现性（相对值，如 0.007）")
    ap.add_argument("--sigma-ext76", type=float, default=None,
                    help="强制指定 207Pb/206Pb 外部重现性（相对值）")
    return ap


def _summary_lines(result) -> list:
    """把结果里的几个关键统计拼成几行文字，用于收尾打印。"""
    res = result.results
    lines = []
    unk = res[res["类型"] == "样品"]
    if len(unk):
        a = res.get("年龄206_238_QC校正", unk["年龄206_238"])
        q = a.quantile([0.05, 0.5, 0.95])
        lines.append(
            f"    样品年龄 206Pb/238U：中位 {q[0.5]:.1f} Ma，"
            f"5–95% 区间 {q[0.05]:.1f} ~ {q[0.95]:.1f} Ma (n={len(unk)})")
        conc = unk["协和度_pct"]
        lines.append(f"    协和度：中位 {conc.median():.1f}%，"
                     f"落在 90–110% 的占 {(conc.between(90, 110)).mean():.0%}")
    return lines


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    cfg = BatchConfig(
        data_dir=args.dir,
        out_excel=args.out,
        primary=args.primary,
        secondary=args.secondary,
        trim=args.trim,
        blank_dur=args.blank_dur,
        n_sigma_common_pb=args.n_sigma_common_pb,
        deadtime_ns=args.deadtime_ns,
        win=args.win,
        step=args.step,
        bulk=args.bulk,
        sigma_ext68=args.sigma_ext68,
        sigma_ext76=args.sigma_ext76,
        do_depth=not args.no_depth,
        plot=args.plot,
        plot_dir=args.plot_dir,
        verbose=not args.quiet,
    )

    print("=" * 72)
    print(f"druid v{__version__}  LA-ICP-MS 锆石 U-Pb 数据还原")
    print(f"批次目录：{cfg.data_dir}")
    print("=" * 72)

    # 真正干活：全部业务逻辑都在 workflow 里
    result = run_batch(cfg)

    # 落盘：与网页界面共用 io.report.export_batch，保证两边产出的表完全一致
    out_path = export_batch(cfg, result, version=__version__)

    print("\n" + "=" * 72)
    print(f"[7] 结果已写出：{out_path}")
    print(f"    结果 {len(result.results)} 行 / QC {len(result.qc)} 行 "
          f"/ 多域明细 {len(result.domains)} 行")
    for line in _summary_lines(result):
        print(line)
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
