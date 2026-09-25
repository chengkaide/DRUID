"""
druid.cli.reduce_batch —— 批处理命令行入口
========================================

用法
----
    # 进入包根目录（含 druid/ 的那一层）后执行
    python -m druid.cli.reduce_batch ^
        --dir "examples\\EX2022A" ^
        --out "结果\\EX2022A_U-Pb结果.xlsx" ^
        --plot

    # 做对照实验：换一种整段比值方法 / 施加死时间校正
    python -m druid.cli.reduce_batch --dir ... --bulk ftau
    python -m druid.cli.reduce_batch --dir ... --deadtime-ns 14.7645

运行 -h 可查看全部参数。

输出两个文件
------------
    <...>_U-Pb结果.xlsx             给人看的表（多 sheet）
    <...>_U-Pb结果.handoff.json     给程序看的质控结论与交接契约

第二个是下游（自动化质控 / 解释流程、ADEPT）的输入，schema 为
`druid.handoff/1`；不想要就加 `--no-json`。质控结论同时会在末尾打印出来，
所以人只看控制台也不会漏掉「这批数据能不能用」。

`examples/EX2022A/` 是仓库自带的示例批次（真实数据，样品名已匿名成 S01…S48，
标样原样保留），可以直接用来验证环境是否装好。
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
from druid.console import ensure_utf8_streams
from druid.core.constants import REFERENCE_PRESETS
from druid.io.handoff import export_handoff
from druid.io.report import age68_column, export_batch
from druid.qc import assess_batch, summary_lines
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

    # ── 质控交接 ──
    # 两个参数互斥：同一条命令里既说"写到哪儿"又说"别写"是自相矛盾，
    # 让它当场报错，别猜用户想要哪个。
    js = ap.add_mutually_exclusive_group()
    js.add_argument("--json", default=None, metavar="PATH",
                    help="质控交接 JSON 的输出路径。"
                         "默认写在结果 Excel 旁边：<结果名>.handoff.json")
    js.add_argument("--no-json", action="store_true",
                    help="不写交接 JSON（只打印质控结论）")

    # ── 标样 ──
    ap.add_argument("--primary", default="91500", help="主标（用于归一化）")
    ap.add_argument("--secondary", default="Ple", help="监控标样（用于 QC 与 σext）")
    ap.add_argument("--ref-preset", default="horstwood2016", choices=list(REFERENCE_PRESETS),
                    help="标样参考值的口径预设（见 core.constants.REFERENCE_PRESETS）。"
                         "默认 horstwood2016 = Horstwood et al. (2016) CA-ID-TIMS 推荐值；"
                         "repo = 仓库第一版旧口径（复现历史结果用）；"
                         "wiedenbeck1995 / self-consistent / isoclock 见源码注释。"
                         "换口径是按同一个因子整体平移年龄，不是纠错。")

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
        # ⚠ 这里曾经写成 `res.get("年龄206_238_QC校正", unk["年龄206_238"])`。
        # `DataFrame.get(key)` 返回的是**整表**的那一列，不是筛选后的子集 ——
        # 于是 35 个标样测点（91500 约 1060 Ma、Ple 约 343 Ma）混进了
        # "样品年龄"的统计：中位数从 458.1 抬到 460.9，
        # 5–95% 区间从 420~644 撑成 333~1044 Ma。
        # 打印出来的 n=48 又让人以为只统计了样品，所以这个错很难被看出来。
        # 先定列名，再从 unk 里取列，两步分开，不要再合成一步。
        # 列名怎么定统一在 io.report.age68_column()，别在这儿再写一遍。
        col = age68_column(res)
        a = unk[col]
        q = a.quantile([0.05, 0.5, 0.95])
        lines.append(
            f"    样品年龄 206Pb/238U：中位 {q[0.5]:.1f} Ma，"
            f"5–95% 区间 {q[0.05]:.1f} ~ {q[0.95]:.1f} Ma (n={len(unk)})")
        conc = unk["协和度_pct"]
        lines.append(f"    协和度：中位 {conc.median():.1f}%，"
                     f"落在 90–110% 的占 {(conc.between(90, 110)).mean():.0%}")
    # 不分域（整段）口径 —— 全部测点走同一条路，跨测点对比就靠它。
    # getattr 容错：旧版结果对象（或测试里伪造的）没有 overall 这个字段。
    ov = getattr(result, "overall", None)
    if ov is not None and not getattr(ov, "empty", True) and "判定" in ov.columns:
        n_ok = int((ov["判定"].astype(str) == "整段常数").sum())
        lines.append(
            f"    不分域整段年龄：{n_ok}/{len(ov)} 个测点与『整段只有一个年龄』"
            f"相容，中位 MSWD {ov['MSWD'].median():.2f}"
            f"（其余仅可用于横向对比，不可定年）")
    return lines


def main(argv=None) -> int:
    # 这个命令的输出全是中文。Windows 上把输出重定向到文件/管道时，
    # Python 会用 locale 编码（cp1252 之类）编码 stdout，一句 print 就抛
    # UnicodeEncodeError 把进程带走 —— 而且是在跑完批处理、准备打印结果的时候。
    # 在入口把流切到 UTF-8 一次解决，别在每个 print 上包 try。
    ensure_utf8_streams()

    args = build_parser().parse_args(argv)

    cfg = BatchConfig(
        data_dir=args.dir,
        out_excel=args.out,
        primary=args.primary,
        secondary=args.secondary,
        ref_preset=args.ref_preset,
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

    # 质控：把「这批数据能不能用」从控制台散文变成可机读的检查项，
    # 并写成交接 JSON 供下游（自动化质控 / 解释流程、ADEPT）读取。
    # 先算检查项再写盘：写盘函数吃 checks，避免同一次运行里算两遍
    # （算两遍就有一边被改坏而另一边没改的风险，早晚对不上）。
    checks = assess_batch(result, cfg)
    handoff_path = None
    if not args.no_json:
        handoff_path = export_handoff(cfg, result, checks=checks,
                                      version=__version__, path=args.json)

    print("\n" + "=" * 72)
    print(f"[7] 结果已写出：{out_path}")
    print(f"    结果 {len(result.results)} 行 / QC {len(result.qc)} 行 "
          f"/ 多域明细 {len(result.domains)} 行 "
          f"/ 不分域年龄 {len(getattr(result, 'overall', []))} 行")
    for line in _summary_lines(result):
        print(line)
    if handoff_path is not None:
        print(f"[8] 质控交接已写出：{handoff_path}")
    else:
        print("[8] 质控交接：按 --no-json 跳过，仅打印结论")
    # detail_ref 跟着实际情况走：没写盘就别指着一个不存在的文件说"见这里"。
    for line in summary_lines(checks,
                              detail_ref=(handoff_path.name
                                          if handoff_path is not None else None)):
        print(line)
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
