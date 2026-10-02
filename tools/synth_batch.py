# -*- coding: utf-8 -*-
"""
tools/synth_batch.py —— 合成一个「标样 + 样品」批次，走完整流水线（写盘、可复算）
================================================================================

为什么需要它
------------
现有两套验证证的都是「算得和上次一样」，**不是「算得对」**：

    · `tests/run_all.py`（自检 N 项）—— 测的是各函数自己的契约；
    · `tests/check_example_batch.py` —— 拿真实批次跑一遍，与写死的基线逐位比。

真数据的根本问题在于**没有真值**：分域年龄与边界、MSWD 的期望、2σ 的覆盖率、
F(τ) 有没有偏，都没有可以对照的答案。所以「算法在什么条件下会静默偏掉」
这类问题，只能靠**知道自己喂进去的是什么**的数据来问。

这个脚本造的就是那种数据：按 Qtegra 的导出格式写盘，走**完整**的
`druid.workflow.run_batch`，然后用「喂进去的年龄」去比对「报回来的年龄」。

三件事都不是从零造轮子
----------------------
1. **CSV 形态**：与 `tests/test_input_boundary.py` 的 `CSV_TEXT` 同构
   （标题行 / 元数据 / 表头 / 单位行 / 数据行；表头与数据行都带尾逗号）。
   ⚠ 那份模板只有 7 列、**没有 202Hg**，而 `core.constants.MASSES_NEEDED`
   要 202 ⇒ 直接拿它当输入会在 `load_spot` 里报「缺少通道 [202]」。这里补上。
2. **正演**：`core.geochronology.r68_of_age / r76_of_age`，不另写公式。
3. **参考值**：`core.references.std_ref`，不手抄比值。

噪声模型（刻意的、也是有限的）
------------------------------
逐行按 **Poisson 计数**采样：把 cps 当成「每 dwell 秒的计数 ÷ dwell」，
dwell 取 0.01 s（与真实采集一致）⇒ 一行的计数 = cps × 0.01。
这是默认**唯一**的噪声源；额外的闪烁噪声用 `--flicker` 打开。

⚠ **已知局限**（写在这里，免得被当成"已经验证过"）
    · 模型是**完全线性**的：没有死时间、没有检测器非线性、没有 down-hole 分馏、
      没有普通铅（净 204Pb 期望为 0）、没有剥蚀速率变化。
      因此它**不能**用来检验二次校正要解决的那个问题
      （主标与样品 238U 差一个数量级时的非线性残差）。
    · 时间剖面是**均一**的（每个测点只有一个真值年龄），所以它与
      「分域分辨力」「分域盲区」那类问题无关 —— 要让信号随深度变才行，
      那是下一步的事。
    · 只检验 **206Pb/238U** 这一个钟。207Pb/206Pb、208Pb/232Th 的年龄虽然也
      出现在结果表里（且是按同一批比值算的），本脚本不拿它们做判据。

用法
----
    python tools/synth_batch.py                        # 460 Ma × 10 个种子
    python tools/synth_batch.py --age 460 --verbose    # 看完整流水线日志
    python tools/synth_batch.py --age 30,155,460,1063  # 跨年龄闭合
    python tools/synth_batch.py --keep                 # 保留写下的批次目录
    python tools/synth_batch.py --nsigma-204 1e9       # 完全不扣普通铅（对照）

首跑量到了什么（2026-10-02，报告见 `G:/UPb归档/2026-10-02/合成数据验证_评估/`）
--------------------------------------------------------------------------------
· **闭合成立**：30 / 60 / 100 / 155 / 260 / 460 / 1063 Ma 七个真值，
  各 20 个种子，pull 均值 −0.40…+0.36、pull 标准差 0.84…1.07
  ⇒ 报出的 σ 是诚实的（不多不少）。
· **204 判据的实际门槛只有 1.15σ，不是 2σ**：假阳性率实测 11.9%（480 个测点）
  而名义 2.3%。原因是刀切 σ(i204) 只算了窗口内的散度（86 行），**漏掉了
  空白段那 47 行带进来的方差** ⇒ σ_真/σ_jk = 1.73 ⇒ 2/1.73 = 1.15σ，
  对应的理论假阳性率 12.4%，与实测吻合。
· **均一剖面被误判"多域"的比例强烈依赖年龄**：30 Ma 65%、60 Ma 35%、
  ≥100 Ma 0%（各 20 个测点）。留下的域差是 3.3–6.6σ（工具自洽），
  但出现率说明**域差的真实标准差约为报出 σ 的 3.6 倍**；把滑窗改成不重叠
  （win/step = 4/4）后这个误判率降到 0 —— 指向**重叠窗口不被当作独立观测**。

验收标志
--------
**正演 → 反演闭合**：报回的样品加权平均年龄与喂进去的真值之差，
应当落在它自己的不确定度之内。脚本打印 `pull = (报回 − 真值) / σ`，
并给出多种子下的 pull 分布：**均值应 ≈ 0**（无系统偏差）、
**标准差应 ≈ 1**（σ 报得不多不少）。
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from druid.console import ensure_utf8_streams                        # noqa: E402
from druid.core.constants import (DEFAULT_REF_PRESET, HG202_204,     # noqa: E402
                                  L232, MASSES_NEEDED, ROLE_LABEL_CN,
                                  ROLE_UNKNOWN)
from druid.core.geochronology import r68_of_age, r76_of_age          # noqa: E402
from druid.core.references import std_alias, std_ref                 # noqa: E402
from druid.core.statistics import weighted_mean                      # noqa: E402
from druid.workflow import BatchConfig, run_batch                    # noqa: E402


# ═════════════════════════════════════════════════════════════════════════════
# 一、采集形态（全部有实测出处，不是随手取的圆整数）
# ═════════════════════════════════════════════════════════════════════════════
DT = 0.31871          # 行间隔中位 (s)：一行 = 一轮全质量扫描，与激光无关
N_BLANK = 62          # 剥蚀前的行数（62 × DT ≈ 19.8 s，够 blank_dur=15 s 取窗）
N_PLATEAU = 95        # 剥蚀平台行数（≈30.3 s）
N_WASH = 82           # 冲洗拖尾行数
DWELL = 0.01          # 每个质量数的驻留时间 (s)
TAIL_FRAC = 0.015     # 拖尾 = 净峰高的 1.5%：低于 find_ablation 的第二级阈值
                      # （3%），所以它**必须**被砍掉；砍不掉就说明两级阈值坏了
N_ROWS = N_BLANK + N_PLATEAU + N_WASH      # = 239，末点 238 × DT ≈ 75.85 s

# 气体空白（cps）。204 的空白刻意取 202/4.35 ——
# 这样「净 204Pb」的期望正好是 0（没有普通铅），而噪声照旧是真实的：
# 「两个大数相减求小数」这件事照样会发生，只是期望值落在 0 上。
BG = {238: 300.0, 206: 20.0, 207: 5.0, 208: 5.0, 232: 10.0, 202: 300.0,
      204: 300.0 / HG202_204}

# 表头列顺序：Time + 七个 MASSES_NEEDED。⚠ 202Hg 必须在（见模块 docstring）。
HEADER = ["Time", "202Hg", "204Pb", "206Pb", "207Pb", "208Pb", "232Th", "238U"]
COL_MASS = [202, 204, 206, 207, 208, 232, 238]

PRIMARY = "91500"
SECONDARY = "Ple"


# ═════════════════════════════════════════════════════════════════════════════
# 二、正演：真值年龄 → 期望信号
# ═════════════════════════════════════════════════════════════════════════════
def forward_ratios(age_ma):
    """
    真值年龄 (Ma) → (R68, R76, R82)。

    R68 / R76 直接调 `core.geochronology` 的正演函数（**不另写公式**）。
    R82 = 208Pb/232Th 走同一个形式、换 λ232 —— 它不进本脚本的判据，
    只是让 208/232 两个通道有物理上说得过去的值，免得造出一份
    "Th 通道随便填"的数据误导后来的人。
    """
    return (float(r68_of_age(age_ma)), float(r76_of_age(age_ma)),
            float(np.expm1(L232 * age_ma * 1e6)))


def spot_levels(r68, r76, r82, u_cps, th_u):
    """
    造出三阶段（空白 / 剥蚀平台 / 冲洗拖尾）的**期望** cps 数组。

    拖尾取净峰高的 TAIL_FRAC —— 这是 `find_ablation` 那段注释里描述的
    真实形态（"关激光后 10 s 仍有 100~200 cps，远高于本底噪声"）。

    返回 {质量数: 长度 N_ROWS 的期望 cps}
    """
    plat = np.zeros(N_ROWS, dtype=bool)
    plat[N_BLANK:N_BLANK + N_PLATEAU] = True
    tail = ~plat
    tail[:N_BLANK] = False

    i238 = float(u_cps)
    i232 = th_u * i238
    peak = {238: i238, 206: r68 * i238, 207: r76 * r68 * i238,
            232: i232, 208: r82 * i232,
            202: BG[202], 204: BG[204]}          # Hg 与普通铅不随剥蚀变化

    levels = {}
    for m in MASSES_NEEDED:
        net = peak[m] - BG[m]
        v = np.full(N_ROWS, BG[m], dtype=float)
        v[plat] = BG[m] + net
        v[tail] = BG[m] + TAIL_FRAC * net
        levels[m] = v
    return levels


def _counts_to_cps(level, rng, flicker):
    """
    期望 cps → 实测 cps：先按 Poisson 抽计数，再 ÷ dwell 换回 cps。

    为什么要绕这一圈：真实文件里就是 cps，而 cps 的噪声来自**计数**。
    直接给 cps 加高斯噪声会让 σ 与信号强度的关系（∝1/√信号）消失，
    而那正是"弱峰比强峰更不准"的来源。
    """
    lam = np.maximum(np.asarray(level, float) * DWELL, 0.0)
    cps = rng.poisson(lam) / DWELL
    if flicker > 0.0:
        # 额外的乘性闪烁噪声（等离子体 1/f 噪声的粗模型），默认关闭。
        cps = cps * (1.0 + rng.normal(0.0, flicker, size=cps.shape))
    return cps


def spot_csv_text(title, levels, rng, flicker):
    """按 Qtegra 的导出格式渲染一个测点的 CSV 文本。"""
    t = np.arange(N_ROWS) * DT
    cols = [_counts_to_cps(levels[m], rng, flicker) for m in COL_MASS]

    lines = [
        f"{title}:03/01/2022 06:50:39 AM;",
        "Software:Name=Qtegra;Version=2.8.3170.309;File Version=1;",
        "RF Generator:RF Plasma Lit Readback=1;Plasma Power Readback=1548.61;",
        "Pulse Counting:Threshold=2500000;",
        "",
        ",".join(HEADER) + ",",
        # 单位行：第一个字段为空 ⇒ read_qtegra 的 float("") 抛异常后自然跳过。
        # 这正是 test_input_boundary 里钉住的那个契约，别改成固定 skiprows。
        "," + ",".join(f"dwell time={DWELL}" for _ in HEADER[1:]) + ",",
    ]
    for j in range(N_ROWS):
        lines.append(f"{t[j]:.5f}," + ",".join(f"{c[j]:.4f}" for c in cols) + ",")
    return "\n".join(lines) + "\n"


# ═════════════════════════════════════════════════════════════════════════════
# 三、造整批（写盘）
# ═════════════════════════════════════════════════════════════════════════════
def _rmtree(path):
    """
    只用 pathlib 逐个删（这台机器上删除被改道回收站、且批量删会撞保护，
    所以逐个来、失败不抛）。只删自己刚写下的东西。
    """
    p = Path(path)
    if not p.exists():
        return
    for f in sorted(p.rglob("*"), reverse=True):
        try:
            if f.is_file() or f.is_symlink():
                f.unlink()
            elif f.is_dir():
                f.rmdir()
        except OSError:
            pass
    try:
        p.rmdir()
    except OSError:
        pass


def make_batch(out_dir, batch, ages, repeat, u_cps, th_u, seed,
               ref_preset, flicker=0.0):
    """
    造一个完整批次目录：`<out_dir>/<batch>/` 下是 `<batch>_LIST.csv` + 每个测点一个 CSV。

    目录名就是批次名 —— DRUID 靠**目录名**找序列表
    （`BatchConfig.list_file` → `<目录名>_LIST.*`），改不得。

    返回 (batch_dir, plan)：
        plan 里每人一项 `dict(file, sample, role, true_age)`，
        `true_age` 为 None 表示这是标样（标样按参考值造，见下）。
    """
    d = Path(out_dir) / batch
    _rmtree(d)
    d.mkdir(parents=True)

    rng = np.random.default_rng(seed)
    ref68_p, ref76_p = std_ref(std_alias(PRIMARY) or PRIMARY, ref_preset)
    ref68_s, ref76_s = std_ref(std_alias(SECONDARY) or SECONDARY, ref_preset)

    # 标样按**参考值**造（不是按参考年龄反算）——
    # 这正是"外标归一化"的假设本身：主标就是参考值。于是 F = 参考/实测 ≈ 1，
    # 而样品报回的年龄 = 真值（在噪声范围内），这就是要检验的闭合。
    # ⚠ 若要检验"标样不在参考值上会怎样"，那是另一件事（外推配对），
    #   把标样的真值年龄改成别的即可，本函数的接口已经够用。
    # 排布：开头两个主标打底，之后**每个样品测点后面跟一对 [主标, 监控标样]**。
    #
    # 为什么标样要这么密：`workflow.estimate_external` 用监控标样的散度估 σext，
    # `apply_secondary_correction` 又用它算二次校正系数 kfac —— 两者都是**样本量
    # 敏感**的。只有 1 个监控标样时，kfac 完全由那一个点的噪声决定，会把整个
    # 样品年龄整体平移（实测过：单点时 kfac = 0.9909 ⇒ 样品偏 −0.9%），
    # 而 kfac 的这份不确定度**没有被传播进 s68**（只乘了个系数）。
    # 真实批次里标样是穿插在样品之间的（示例批次 48 个样品点配 21 主标 + 14 监控），
    # 这里按同样的密度造，闭合检验才有意义。
    plan = [dict(sample=PRIMARY, r68=ref68_p, r76=ref76_p, true_age=None),
            dict(sample=PRIMARY, r68=ref68_p, r76=ref76_p, true_age=None),
            dict(sample=SECONDARY, r68=ref68_s, r76=ref76_s, true_age=None)]
    for age in ages:
        r68, r76, r82 = forward_ratios(age)
        for _ in range(repeat):
            plan.append(dict(sample=f"SYN-{age:g}", r68=r68, r76=r76, r82=r82,
                             true_age=float(age)))
            plan.append(dict(sample=PRIMARY, r68=ref68_p, r76=ref76_p,
                             true_age=None))
            plan.append(dict(sample=SECONDARY, r68=ref68_s, r76=ref76_s,
                             true_age=None))

    entries = []
    for i, item in enumerate(plan):
        fname = f"{batch}_{i + 1}"
        r82 = item.get("r82")
        if r82 is None:
            # 标样的 208/232 也给一个说得过去的值（用参考年龄正演），
            # 免得造出"标样 Th 通道是垃圾值"的数据。
            r82 = forward_ratios(985.0)[2]
        levels = spot_levels(item["r68"], item["r76"], r82, u_cps, th_u)
        # CRLF：Qtegra 导出就是 CRLF（read_qtegra 两种都认，但真的更像真的）
        text = spot_csv_text(item["sample"], levels, rng, flicker)
        (d / f"{fname}.csv").write_bytes(text.replace("\n", "\r\n").encode("utf-8"))
        item["file"] = fname
        entries.append((fname, item["sample"]))

    # 序列表：两列、无表头。分隔符定死逗号 —— 样品名里有空格、
    # 有的行没有，交给 pandas 嗅探有被误判成空格分隔的风险（见 io/sequence.py）。
    (d / f"{batch}_LIST.csv").write_text(
        "".join(f"{f},{s}\n" for f, s in entries), encoding="utf-8")
    return d, plan


# ═════════════════════════════════════════════════════════════════════════════
# 四、跑完整流水线并读回结果
# ═════════════════════════════════════════════════════════════════════════════
def run_and_collect(batch_dir, cfg_kwargs):
    """跑一次 `run_batch`，把该看的量整理成普通 dict。"""
    cfg = BatchConfig(data_dir=str(batch_dir), **cfg_kwargs)
    br = run_batch(cfg)
    res = br.results
    samples = res[res["类型"] == ROLE_LABEL_CN[ROLE_UNKNOWN]]
    # 按**样品名**分组（不是把所有样品并在一起 —— 多年龄时那样会把几个
    # 真值不同的测点混成一堆，MSWD 直接飞到几万）。
    depth_by_name = {}
    for _, row in samples.iterrows():
        depth_by_name.setdefault(str(row["样品"]), []).append(
            str(row["深度结构"]))
    # 触发了普通铅扣除的测点数（f206 > 0 ⇔ `reduce_interval` 的 204 显著性检验通过）。
    # 合成数据的**真值**是"没有普通铅"（净 204Pb 期望为 0），所以这里的每一个
    # 都是**假阳性** —— 而假阳性会真的把 206Pb 扣掉一点、把年龄拉年轻。
    # 把它数出来，是为了让"这条判据有多容易误报"变成可观测的量，而不是靠猜。
    n_common = int((res["f206_pct"] > 0).sum())
    n_common_sample = int((samples["f206_pct"] > 0).sum())
    return dict(
        res=res, qc=br.qc, info=br.info, overall=br.overall,
        depth_by_name=depth_by_name,
        n_sample=int(len(samples)),
        n_common=n_common, n_spots=int(len(res)),
        n_common_sample=n_common_sample,
        primary_rsd=br.info.get("primary_rsd68"),
        sd68=br.info.get("sd68"),
    )


def sample_stats(res, name, true_age):
    """
    把**同一个真值年龄**的那几个样品测点合起来，返回加权平均与 pull。

    用 `年龄206_238_QC校正` 而不是 `年龄206_238`：前者才是工具对外报的
    "样品中的最终年龄"（见 workflow.apply_secondary_correction）。
    """
    g = res[res["样品"] == name]
    mu, se, mswd, n = weighted_mean(g["年龄206_238_QC校正"],
                                    g["s68_1sig_QC校正"])
    if not np.isfinite(mu) or not np.isfinite(se) or se <= 0:
        return None
    return dict(mean=float(mu), sigma=float(se), mswd=float(mswd), n=int(n),
                delta=float(mu - true_age),
                pct=float((mu - true_age) / true_age * 100.0),
                pull=float((mu - true_age) / se))


# ═════════════════════════════════════════════════════════════════════════════
# 五、命令行
# ═════════════════════════════════════════════════════════════════════════════
def _default_out() -> Path:
    """
    默认写哪儿：**G 盘优先**（这台机器 C 盘长期贴红线），没有 G 盘才退系统临时目录。
    可用环境变量 DRUID_SYNTH_DIR 覆盖。
    """
    import os
    env = os.environ.get("DRUID_SYNTH_DIR")
    if env:
        return Path(env)
    g = Path("G:/")
    if g.exists():
        return g / "_zd_run" / "synth_batch"
    return Path(tempfile.gettempdir()) / "druid_synth"


def _parse_ages(text):
    ages = [float(x) for x in str(text).replace("，", ",").split(",") if x.strip()]
    if not ages:
        raise SystemExit("--age 至少要给一个年龄")
    return ages


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="合成「标样 + 样品」批次，走完整 run_batch，检验正演→反演闭合。")
    ap.add_argument("--age", default="460",
                    help="样品真值年龄 (Ma)，可逗号分隔多个（默认 460）")
    ap.add_argument("--repeat", type=int, default=3,
                    help="每个年龄造几个测点（默认 3；每个样品后跟一对 [主标, 监控标样]）")
    ap.add_argument("--seeds", type=int, default=10,
                    help="跑多少个随机种子（默认 10；1 = 单次）")
    ap.add_argument("--seed", type=int, default=20261002, help="起始种子")
    ap.add_argument("--u-cps", type=float, default=1.0e6,
                    help="238U 剥蚀平台的计数率 (cps，默认 1e6)")
    ap.add_argument("--th-u", type=float, default=0.5, help="表观 Th/U（默认 0.5）")
    ap.add_argument("--flicker", type=float, default=0.0,
                    help="额外乘性闪烁噪声（相对，默认 0 = 只有计数噪声）")
    ap.add_argument("--nsigma-204", type=float, default=2.0,
                    help="204Pb 显著性门槛（σ 倍数，默认 2；调到 1e9 等于不扣普通铅）")
    ap.add_argument("--ref-preset", default=DEFAULT_REF_PRESET,
                    help=f"参考值档（默认 {DEFAULT_REF_PRESET}）")
    ap.add_argument("--batch", default="SYNTH", help="批次名（＝生成的目录名）")
    ap.add_argument("--out", default=None, help="写到哪个目录（默认见 --keep 说明）")
    ap.add_argument("--no-depth", action="store_true", help="关闭深度域判别")
    ap.add_argument("--keep", action="store_true",
                    help="保留写下的批次目录（默认跑完即删）")
    ap.add_argument("--verbose", action="store_true",
                    help="转发 run_batch 的完整中文进度日志")
    args = ap.parse_args(argv)

    ensure_utf8_streams()

    ages = _parse_ages(args.age)
    out_root = Path(args.out) if args.out else _default_out()
    out_root.mkdir(parents=True, exist_ok=True)

    cfg_kwargs = dict(ref_preset=args.ref_preset, do_depth=not args.no_depth,
                      verbose=bool(args.verbose), plot=False,
                      n_sigma_common_pb=args.nsigma_204)

    print(f"合成批次：{len(ages)} 个真值年龄 × repeat={args.repeat} × "
          f"{args.seeds} 个种子；238U ≈ {args.u_cps:.3g} cps，"
          f"同行 {N_ROWS} 行 / {N_ROWS * DT:.2f} s")
    print(f"输出目录：{out_root}（{'保留' if args.keep else '跑完即删'}）")
    print()

    hdr = (f"{'种子':>6}  {'真值Ma':>8}  {'报回Ma':>10}  {'1σ':>6}  "
           f"{'Δ Ma':>8}  {'Δ%':>7}  {'pull':>6}  {'MSWD':>5}  {'n':>2}  "
           f"{'多域':>4}")
    per_age = {a: [] for a in ages}
    multi_by_age = {a: [0, 0] for a in ages}      # [多域测点数, 该年龄测点数]
    qc_seen = []
    n_common_total = 0
    n_common_sample_total = 0
    n_sample_total = 0
    n_spots_total = 0
    first_dir = None

    for s in range(args.seeds):
        seed = args.seed + s
        # 每个种子一个独立子目录 —— 目录名必须仍然是批次名（DRUID 靠它找序列表）
        sub = out_root if args.seeds == 1 else out_root / f".seed_{s + 1}"
        d, plan = make_batch(sub, args.batch, ages, args.repeat,
                             args.u_cps, args.th_u, seed, args.ref_preset,
                             args.flicker)
        if first_dir is None:
            first_dir = d
        got = run_and_collect(d, cfg_kwargs)
        qc_seen.append(got["qc"])
        n_common_total += got["n_common"]
        n_common_sample_total += got["n_common_sample"]
        n_sample_total += got["n_sample"]
        n_spots_total += got["n_spots"]

        if s == 0:
            print("— 单次流水线的自检项 —")
            print(f"  主标 F 的逐点 RSD（唯一直的有分辨力的标样稳不稳的量）："
                  f"{_pct(got['primary_rsd'])}")
            print(f"  由监控标样散度估出的外部重现性 σext(206/238)："
                  f"{_pct(got['sd68'])}")
            print()
            print(hdr)
            print("-" * len(hdr))
        for age in ages:
            name = f"SYN-{age:g}"
            structs = got["depth_by_name"].get(name, [])
            multi_by_age[age][0] += sum(s0.startswith("多域") for s0 in structs)
            multi_by_age[age][1] += len(structs)
            st = sample_stats(got["res"], name, age)
            if st is None:
                print(f"{seed:>6}  {age:>8.1f}  （该年龄无有效测点）")
                continue
            per_age[age].append(st)
            print(f"{seed:>6}  {age:>8.1f}  {st['mean']:>10.2f}  {st['sigma']:>6.2f}  "
                  f"{st['delta']:>+8.2f}  {st['pct']:>+7.3f}  {st['pull']:>+6.2f}  "
                  f"{st['mswd']:>5.2f}  {st['n']:>2}  "
                  f"{('是' if any(s0.startswith('多域') for s0 in structs) else '否'):>4}")

    print()
    print("— 汇总（跨种子）—")
    print(f"{'真值Ma':>8}  {'平均报回':>10}  {'平均1σ':>7}  {'平均ΔMa':>9}  "
          f"{'平均Δ%':>8}  {'pull均值':>9}  {'pull标准差':>10}  {'判读':>6}")
    verdict_ok = True
    for age in ages:
        rows = per_age[age]
        if not rows:
            continue
        mean_age = float(np.mean([r["mean"] for r in rows]))
        mean_sig = float(np.mean([r["sigma"] for r in rows]))
        mean_d = float(np.mean([r["delta"] for r in rows]))
        mean_p = float(np.mean([r["pct"] for r in rows]))
        pulls = np.array([r["pull"] for r in rows])
        pm = float(pulls.mean())
        ps = float(pulls.std(ddof=1)) if pulls.size > 1 else float("nan")
        # 判据：均值不能显著偏离 0（3/√N 的宽容带）；标准差应 ≈ 1。
        # ⚠ 种子少时 σ 的估计本身就不稳（N=10 时 ±23%），所以带宽给得宽。
        good = abs(pm) < 3.0 / max(np.sqrt(pulls.size), 1.0)
        if np.isfinite(ps):
            good = good and (0.5 <= ps <= 1.8)
        verdict_ok &= bool(good)
        print(f"{age:>8.1f}  {mean_age:>10.2f}  {mean_sig:>7.2f}  {mean_d:>+9.3f}  "
              f"{mean_p:>+8.4f}  {pm:>+9.2f}  {ps:>10.2f}  "
              f"{'✓ 闭合' if good else '✗ 检查':>6}")

    # 标样 QC：把各种子的 QC 表按标样名合起来看偏差
    print()
    print("— 标样 QC（各种子合并）—")
    allqc = _concat(qc_seen)
    if allqc is not None and not allqc.empty:
        for name, grp in allqc.groupby("标样", sort=False):
            print(f"  {name:<10s} 参考 {grp['参考年龄_Ma'].iloc[0]:8.2f} Ma  "
                  f"偏差 中位 {grp['偏差_pct'].median():+.4f}%  "
                  f"范围 [{grp['偏差_pct'].min():+.4f}, {grp['偏差_pct'].max():+.4f}]%  "
                  f"MSWD 中位 {grp['MSWD'].median():.2f}  n={int(grp['点数'].iloc[0])}")

    print()
    print("— 深度结构（真剖面全程均一，本不该被判成多域）—")
    kn, kt = 0, 0
    for age in ages:
        m, t = multi_by_age[age]
        kn += m
        kt += t
        if t:
            print(f"  {age:>8.1f} Ma：{m} / {t}（{m / t:.1%}）")
    print(f"  合计：{kn} / {kt}（{kn / max(kt, 1):.1%}）")

    # 204 显著性检验的假阳性率。
    # 合成数据的真值是"没有普通铅"，所以这里每一个扣除都是误报 ——
    # 误报会真把 206Pb 扣掉一点（样品上实测过 −0.6% 的 R68）、年龄跟着变年轻。
    print()
    print("— 普通铅判据的假阳性率（真值：无普通铅）—")
    print(f"  全部测点 {n_spots_total} 个，触发扣除 {n_common_total} 个 "
          f"({n_common_total / max(n_spots_total, 1):.1%}；"
          f"理论上单侧 2σ 门槛 ≈ 2.3%)")
    print(f"  其中样品测点 {n_common_sample_total} / {n_sample_total} "
          f"({n_common_sample_total / max(n_sample_total, 1):.1%})")

    if args.keep:
        print()
        print(f"批次目录保留在：{first_dir}")
    else:
        _rmtree(out_root if args.seeds == 1 else out_root)
        print()
        print("已删除临时批次目录（要保留加 --keep）")

    print()
    print("结论：" + ("正演 → 反演闭合成立（报回年龄在自身 1σ 内、无系统偏差）"
                     if verdict_ok else "有年龄点未闭合，需要检查"))
    return 0 if verdict_ok else 1


def _pct(x):
    return "（未计算）" if x is None else f"{float(x) * 100:.3f}%"


def _concat(frames):
    import pandas as pd
    keep = [f for f in frames if f is not None and not f.empty]
    return pd.concat(keep, ignore_index=True) if keep else None


if __name__ == "__main__":
    raise SystemExit(main())
