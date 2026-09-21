"""
druid.workflow —— 批处理编排层
=============================

这一层**不包含任何新的算法**，只做一件事：
把 io / reduction / depth 提供的零件按正确的顺序组装成一条流水线。

为什么要单独一层
----------------
· core / reduction / depth **互不依赖**，可以单独换掉某一环；
· 算法的"是什么"和流程的"先做什么后做什么"是两种不同的知识，
  混在一起会让两者都难以修改；
· 命令行只是一个参数解析器，真正的业务逻辑在这里，
  将来要做 GUI 或 Jupyter 调用，直接复用 run_batch() 即可。

流水线十步
----------
    ① 读序列 LIST.xls                      io.sequence
    ② 逐点装载 + 整段还原                   reduction.trace / ratios
    ③ 主标建立深度剖面（供 F(τ) 使用）      depth.windows
    ④ 构造夹逼关系 bracketing
    ⑤ 计算整段比值（simple / ftau 两种）
    ⑥ 主标归一化因子 F（稳健剔除离群后插值）
    ⑦ 用监控标样散度估计外部重现性 σext
    ⑧ 逐点算年龄与不确定度
    ⑨ 深度剖面年龄域判别
    ⑩ QC 汇总 + 监控标样二次校正
"""
from __future__ import annotations

import sys
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from .console import ensure_utf8_streams
from .core.constants import (
    L238,
    ROLE_GLASS,
    ROLE_LABEL_CN,
    ROLE_PRIMARY,
    ROLE_SECONDARY,
    ROLE_UNKNOWN,
    STANDARDS,
    U238_U235,
)
from .core.geochronology import age68, age75, age76
from .core.references import std_age, std_alias, std_ref
from .core.statistics import external_scatter, robust_mask, weighted_mean
from .depth.domains import merge_close, refine_domains, segment, summarize_segments
from .depth.fractionation import bracket_F, profile_ages
from .depth.windows import window_profile, window_sums
from .io.sequence import read_sequence, sequence_summary, spot_csv_path
from .reduction.ratios import reduce_interval
from .reduction.trace import Tra, load_spot, window_mask


# ═════════════════════════════════════════════════════════════════════════════
# 配置对象
# ═════════════════════════════════════════════════════════════════════════════
@dataclass
class BatchConfig:
    """一次批处理需要的全部可调参数。**所有魔法数字都应该出现在这里**。"""

    # ── 输入输出 ──
    data_dir: str                       # 批次目录，内含 *_LIST.xls 与 *_N.csv
    out_excel: Optional[str] = None     # 结果 Excel 路径

    # ── 标样设置 ──
    primary: str = "91500"              # 主标名（用于归一化）
    secondary: str = "Ple"              # 监控标样名（用于 QC 与 σext 估计）

    # ── 单点还原 ──
    trim: float = 1.5                   # 剥蚀段两端各裁掉多少秒（避开开关激光瞬态）
    blank_dur: float = 15.0             # 气体空白取样时长 (s)
    n_sigma_common_pb: float = 2.0      # 204Pb 显著性门槛（σ 倍数）
    deadtime_ns: float = 0.0            # 探测器死时间 (ns)，0 = 不校正

    # ── 深度剖面 ──
    win: float = 4.0                    # 滑动窗口宽度 (s)
    step: float = 1.0                   # 滑动步长 (s)

    # ── 整段比值方法 ──
    #   "simple"：整段用单一校正因子（本批次实测更准，默认）
    #   "ftau"  ：逐窗口 F(τ) 深度校正后按计数合成（剖面更平，但整段偏差略大）
    bulk: str = "simple"

    # ── 不确定度 ──
    sigma_ext68: Optional[float] = None  # 强制指定外部重现性 (206Pb/238U，相对)
    sigma_ext76: Optional[float] = None  # 强制指定外部重现性 (207Pb/206Pb，相对)

    # ── 深度域判别 ──
    do_depth: bool = True                # 是否做逐点深度剖面年龄域判别
    merge_n_sigma: float = 3.0           # 合并伪分域的统计门槛
    merge_min_frac: float = 0.05         # 合并伪分域的地质门槛（占年龄的百分比）

    # ── 输出 ──
    plot: bool = False                   # 是否生成逐点深度剖面图
    plot_dir: Optional[str] = None       # 图件目录
    verbose: bool = True                 # 是否打印进度

    # ── 内部派生 ──
    _out_dir: Path = field(default=None, repr=False)

    def __post_init__(self):
        """把字符串路径统一转成 Path，并准备好输出目录。"""
        self.data_dir = Path(self.data_dir)
        base = self.data_dir.name                       # 例如 "EX2022A"
        if self.out_excel is None:
            # 默认和原始数据放在一起，便于"数据—结果"一一对应、不会丢
            self.out_excel = self.data_dir / f"{base}_U-Pb结果.xlsx"
        self.out_excel = Path(self.out_excel)
        self._out_dir = self.out_excel.parent
        if self.plot_dir is None:
            self.plot_dir = self._out_dir / f"{base}_深度剖面图"
        else:
            self.plot_dir = Path(self.plot_dir)

    @property
    def list_file(self) -> Path:
        """
        序列文件路径。按**优先级**依次找，第一个存在的就用它。

        为什么把 .csv / .tsv 排在 .xls 前面
        ----------------------------------
        序列表只有两列、几十行，用文本存比二进制好得多：能 diff、能在
        code review 里看懂、出了错能用编辑器直接修。仓库里的示例批次
        （`examples/EX2022A/`）就是这么存的 —— 谁 clone 下来都能看清
        "第几个测点是标样、第几个是样品"。

        两列都不带表头（第 1 行就是数据），这与 xls 版的约定一致。

        仍然兼容 .xls/.xlsx：仪器工作站导出的原始 LIST 就是这个格式，
        日常跑真实批次时不去动它。
        """
        d = self.data_dir
        for cand in (f"{d.name}_LIST.csv", f"{d.name}_list.csv",
                     f"{d.name}_LIST.tsv", f"{d.name}_list.tsv",
                     f"{d.name}_LIST.xls", f"{d.name}_list.xls",
                     f"{d.name}_LIST.xlsx", f"{d.name}_list.xlsx"):
            p = d / cand
            if p.exists():
                return p
        # 找不到就返回最常见的一种，让 FileNotFoundError 在后面自然抛出
        return d / f"{d.name}_LIST.xls"


# ═════════════════════════════════════════════════════════════════════════════
# 结果容器
# ═════════════════════════════════════════════════════════════════════════════
@dataclass
class BatchResult:
    """一次批处理的完整产出。"""
    results: pd.DataFrame        # 逐点结果表
    qc: pd.DataFrame             # 标样 QC 表
    domains: pd.DataFrame        # 多域剖面明细
    info: dict                   # 关键参数与中间量（便于复现与追溯）
    spots: List[Tuple[Tra, dict]] = field(default_factory=list)  # 原始中间数据
    # 逐窗口年龄剖面。列名按 ADEPT 的 Format 4 对齐
    # （Analysis / Time / Age68 / Age68_1s），可以直接写盘交给 R 端做
    # 加权平均与 MSWD，不需要转录或改列名。
    windows: pd.DataFrame = field(default_factory=pd.DataFrame)


def _log(cfg: BatchConfig, msg: str) -> None:
    """统一的进度打印：受 cfg.verbose 控制，并立即 flush（便于重定向到日志）。"""
    if cfg.verbose:
        print(msg)
        sys.stdout.flush()


# ═════════════════════════════════════════════════════════════════════════════
# ①─② 装载整批数据
# ═════════════════════════════════════════════════════════════════════════════
def load_batch(cfg: BatchConfig):
    """
    读序列 + 逐个 CSV 装载 + 整段还原。

    返回
    ----
    (spots, seq, skipped)
        spots   : [(Tra, ratio_dict), ...] 只含成功还原的点
        seq     : 原始序列表
        skipped : [(文件名, 失败原因), ...] 便于排查
    """
    seq = read_sequence(cfg.list_file, cfg.primary, cfg.secondary)
    _log(cfg, f"[1] 序列：{len(seq)} 个测点  →  {sequence_summary(seq)}")

    spots, skipped = [], []
    for i, row in seq.iterrows():
        fname = row["file"]
        # 用统一的 helper 拼路径：LIST 里登记的名字不带 .csv 后缀
        csv_path = spot_csv_path(cfg.data_dir, fname)

        if not csv_path.exists():                       # LIST 里有但实际文件缺失
            skipped.append((fname, "数据文件不存在"))
            continue
        if row["role"] == ROLE_GLASS:                   # 玻璃标样：不参与 U-Pb
            skipped.append((fname, "玻璃标样，不参与 U-Pb"))
            continue

        try:
            # 装载：读文件 → 死时间 → 找剥蚀段 → 扣空白 → 裁边
            tr = load_spot(int(row["order"]) - 1, csv_path, row["sample"], row["role"],
                           trim=cfg.trim, deadtime_ns=cfg.deadtime_ns,
                           blank_dur=cfg.blank_dur)
            # 整段还原：在 [t0, t1] 窗口内算比值与不确定度
            r = reduce_interval(tr.net, window_mask(tr.t, tr.t0, tr.t1),
                                n_sigma=cfg.n_sigma_common_pb)
            if r is None:
                raise ValueError("有效窗口内的点数不足")
        except Exception as e:                          # 单点失败不影响整批
            skipped.append((fname, f"{type(e).__name__}: {e}"))
            continue
        spots.append((tr, r))

    for fname, why in skipped:
        _log(cfg, f"    -- {fname} 跳过（{why}）")
    _log(cfg, f"    有效测点 {len(spots)} 个")
    return spots, seq, skipped


# ═════════════════════════════════════════════════════════════════════════════
# ③─④ 主标剖面与夹逼关系
# ═════════════════════════════════════════════════════════════════════════════
def build_bracketing(spots, cfg: BatchConfig):
    """
    为每个主标建立窗口剖面，并给出"某个序号应取哪几个主标来夹逼"的函数。

    为什么每个主标都要做剖面
    ------------------------
    不光深度分析需要 F(τ)，连"整段比值要不要用 F(τ)"这个选项也需要它
    （见 bulk_ratios_ftau）。主标的窗口剖面是整批处理的基础设施。

    返回
    ----
    brack_for(k) → 第 k 个点应该使用的夹逼标样剖面列表（最多两个：前后各一）

    只返回这一个函数：主标序号与剖面只在闭包里用，调用方拿到也用不上，
    以前返回三元组、调用方用 `_std_prof` 接住再丢掉，等于把内部状态暴露成接口。
    """
    std_prof, std_pos = [], []
    for k, (tr, r) in enumerate(spots):
        if tr.role != ROLE_PRIMARY:
            continue
        # fix=(f206, sk)：主标自己也按整段的普通铅比例处理短窗口，规则要一致
        std_prof.append(window_profile(tr, cfg.win, cfg.step,
                                       fix=(r["f206"], r["sk"]),
                                       n_sigma=cfg.n_sigma_common_pb))
        std_pos.append(k)

    def brack_for(k: int) -> list:
        """
        找到第 k 个点前后最近的两个主标（各一个），构成夹逼对。

        ⚠ 序号 k 是 spots 列表的下标，不是 LIST 里的原始序号；
          因为中间可能跳过玻璃标样或失败的点，两者并不一致。

        用二分定邻居，而不是列表推导筛一遍：这段逻辑对每个测点要跑两次
        （整段比值一次、深度剖面一次），原写法每次都要重建两个长度 n 的
        列表，整批下来是 O(n²)。
        """
        i = bisect_left(std_pos, k)      # 第一个 std_pos >= k
        j = bisect_right(std_pos, k)     # 第一个 std_pos > k
        out = []
        if i > 0:
            out.append(std_prof[i - 1])  # 前一个主标
        if j < len(std_pos):
            out.append(std_prof[j])      # 后一个主标
        return out

    return brack_for


# ═════════════════════════════════════════════════════════════════════════════
# ⑤ 两种"整段比值"
# ═════════════════════════════════════════════════════════════════════════════
def bulk_ratios_ftau(tr: Tra, r: dict, brack, ref68: float,
                     cfg: BatchConfig) -> Tuple[float, float]:
    """
    **ftau 模式**：整段比值 = 逐窗口用 F(τ) 校正后按**计数**合成。

        R68 = Σ_w (206Pb_w · F68(τ_w)) / Σ_w 238U_w

    先看分子为什么要扣除普通铅后才加权：
    普通铅不是随 U 衰变产生的，不该参与"深度分馏"的校正，
    所以先用整段的 f206 把普通铅剔掉，剩下的放射成因部分才乘 F(τ)。

    ⚠ 只有 **206Pb/238U 做 F(τ) 逐窗口校正**；**207Pb/206Pb 直接沿用整段值**。
      它是 Pb 同位素比值，几乎不随坑深分馏（随坑深漂移的是 Pb/U 这类
      元素对元素的比值），而窗口级 207Pb 计数极少（4 s 窗口仅约 130 个计数），
      套 F(τ) 只会注入噪声。详见 druid/__init__.py 的
      "两个针对本批数据的关键设计"。

    返回
    ----
    (R68, R76)。R76 恒等于 simple 模式的值（见上）；窗口为空或分母非正时
    R68 也回退到 simple 模式的值。
    """
    # 每个窗口内各通道的净计数积分和 + 对应的归一化深度 τ
    tau, S = window_sums(tr, cfg.win, cfg.step)
    if len(tau) == 0 or not brack:
        # 没有可用的夹逼标样时（例如批次里只有主标、或单点主标、或前后标样
        # 都被稳健剔除），逐窗口 F(τ) 校正无从构造。直接退回整段（simple）
        # 比值，这样"只有 91500 标样"的校准批次也能正常还原并输出 QC。
        return r["R68"], r["R76"]

    # 在这张 τ 网格上构造主标的分馏因子。
    # 只构造 F68 —— 207Pb/206Pb 不参与逐窗口校正，构造 F76 再丢掉是误导。
    F68 = bracket_F(tau, brack, "R68", ref68)

    # 逐窗口按各窗口 206 强度的比例分配普通铅
    c6 = r["f206"] * S[206]
    p206 = (S[206] - c6) * F68             # 扣除普通铅 + 深度校正后的放射成因 206Pb
    den68 = S[238].sum()                   # 分母：238U 总计数（无需分馏校正）
    if den68 <= 0:
        return r["R68"], r["R76"]
    return float(p206.sum() / den68), r["R76"]


def compute_bulk_ratios(spots, brack_for, ref68, cfg: BatchConfig):
    """
    算出全部测点的两套整段比值，并按 cfg.bulk 选出启用的那套。

    为什么要把两套都算出来
    ----------------------
    两套方案在同一批数据上的差别，本身就是重要的诊断信息：
    如果差别远大于它们各自的 σ，说明 down-hole 分馏很严重，
    或者存在别的结构性问题。结果表里会并列输出两种年龄便于对照。

    两套的 207Pb/206Pb 是同一个数组：ftau 只对 206Pb/238U 做逐窗口校正
    （见 bulk_ratios_ftau）。保留二元组是为了让两种方法取用方式一致。
    """
    raw68 = np.array([r["R68"] for _, r in spots])      # simple：来自整段 reduce_interval
    raw76 = np.array([r["R76"] for _, r in spots])

    ftau68 = np.array([
        bulk_ratios_ftau(tr, r, brack_for(k), ref68, cfg)[0]
        for k, (tr, r) in enumerate(spots)
    ])

    variants = {"simple": (raw68, raw76), "ftau": (ftau68, raw76)}
    if cfg.bulk not in variants:
        raise ValueError(f"未知的整段比值方法 {cfg.bulk}，可选：{list(variants)}")
    return variants, variants[cfg.bulk]


# ═════════════════════════════════════════════════════════════════════════════
# ⑥ 主标归一化因子 F
# ═════════════════════════════════════════════════════════════════════════════
def calibrate_primary(spots, b68, b76, ref68, ref76, cfg: BatchConfig):
    """
    由主标算出整条的归一化因子曲线 F(idx)。

        F = 参考比值 / 实测比值

    为什么 F 是"随序号变化的曲线"而不是一个常数
    -------------------------------------------
    仪器的灵敏度在一次数小时的测试中会缓慢漂移（锥口沉积、等离子体条件变化…）。
    用一个全程常数去归一化，等于假设这段时间仪器纹丝不动。
    正确做法是：把每个主标处的 F 算出来，再**沿序号线性插值**，
    让相邻两个主标之间的样品分享漂移的中间状态。

    为什么要做离群主标剔除
    ----------------------
    偶尔某个主标会打到包裹体（f206 异常偏大）或遇上 Hg 瞬时波动，
    这样的点会把 F 拉偏，进而影响它夹逼的一整段样品。
    用 MAD 稳健统计把它们找出来剔除（不删数据，只是不参与拟合 F）。
    """
    pmask = np.array([tr.role == ROLE_PRIMARY for tr, _ in spots])
    F68_all = ref68 / b68[pmask]
    F76_all = ref76 / b76[pmask]

    # 分别在 206/238 与 207/206 两个体系里找离群，两者都正常的才算合格
    ok68, _ = robust_mask(F68_all)
    ok76, _ = robust_mask(F76_all)
    ok = ok68 & ok76

    rejected = [int(spots[j][0].idx) + 1
                for j, flag in zip(np.flatnonzero(pmask), ok) if not flag]

    # 插值用的节点：合格主标在 spots 列表中的下标
    pidx = np.flatnonzero(pmask)[ok].astype(float)
    F68, F76 = F68_all[ok], F76_all[ok]

    _log(cfg, f"\n[2] 主标 {cfg.primary}   R68_ref={ref68:.6f}   R76_ref={ref76:.6f}")
    _log(cfg, f"    F68 = {F68.mean():.5f} ± {F68.std(ddof=1):.5f}"
              f"  (RSD {F68.std(ddof=1) / F68.mean():.2%}, n={len(F68)})")
    _log(cfg, f"    F76 = {F76.mean():.5f} ± {F76.std(ddof=1):.5f}"
              f"  (RSD {F76.std(ddof=1) / F76.mean():.2%})")
    _log(cfg, f"    整段比值方法：{cfg.bulk}"
              f"（simple = 单一因子；ftau = 逐窗口 F(τ) 深度校正）")
    if rejected:
        _log(cfg, f"    ⚠ 稳健统计剔除离群主标 序号 {rejected}"
                  f"（f206 异常或 Hg 瞬时波动，不参与归一化曲线拟合）")
    return pidx, F68, F76, dict(rejected=rejected, pmask=pmask)


def interp_F(k: int, pidx, F68, F76):
    """
    在序号 k 处插值出该点的归一化因子。超出范围自动取端点值。

    兜底：当批次里没有任何主标（pidx 为空）时，返回 (1.0, 1.0)，
    即"不做归一化"——上层此时已进入未校准诊断模式，比值直接取实测值。
    这样即便调用路径绕过了主流程的判定，也不会再出现
    `np.interp(k, [], [])` 在空数组上抛 `ValueError` 的崩溃。
    """
    if pidx.size == 0:
        return 1.0, 1.0
    return float(np.interp(k, pidx, F68)), float(np.interp(k, pidx, F76))


# ═════════════════════════════════════════════════════════════════════════════
# ⑦ 外部重现性
# ═════════════════════════════════════════════════════════════════════════════
def estimate_external(spots, pidx, F68, F76, b68, b76, cfg: BatchConfig):
    """
    用监控标样的实测散度，估计"内部精度之外"的那部分不确定度。

    两遍法
    ------
    第一遍：σext = 0，先把所有点校正一遍，看看监控标样的
            age68 / R68 到底散到什么程度；
    第二遍：用第一遍得到的 σext 重算，得到最终的不确定度。

    为什么必须用两遍：σext 本身依赖于校正后的比值，而校正是线性的，
    一次迭代即可达到足够精度（不会有明显偏差）。
    """
    qc_values = {"R68": [], "s68": [], "R76": [], "s76": []}
    for k, (tr, r) in enumerate(spots):
        if tr.sample != cfg.secondary:
            continue
        f68, f76 = interp_F(k, pidx, F68, F76)
        qc_values["R68"].append(b68[k] * f68)
        qc_values["s68"].append(r["s68"] * f68)
        qc_values["R76"].append(b76[k] * f76)
        qc_values["s76"].append(r["s76"] * f76)

    sd68 = cfg.sigma_ext68
    sd76 = cfg.sigma_ext76
    if qc_values["R68"]:
        sd68 = sd68 if sd68 is not None else external_scatter(
            qc_values["R68"], qc_values["s68"])
        sd76 = sd76 if sd76 is not None else external_scatter(
            qc_values["R76"], qc_values["s76"], lo=0.001, hi=0.05)
    # 兜底：万一序列里一个监控标样都没有
    sd68 = 0.007 if sd68 is None else float(sd68)
    sd76 = 0.0025 if sd76 is None else float(sd76)

    _log(cfg, f"\n[3] 外部重现性（由 {cfg.secondary} 散度扣除内部精度得到，1σ 相对）")
    _log(cfg, f"    206Pb/238U  {sd68:.2%}      207Pb/206Pb  {sd76:.2%}"
              f"   (n={len(qc_values['R68'])})")
    return sd68, sd76


# ═════════════════════════════════════════════════════════════════════════════
# ⑧ 逐点结果与年龄
# ═════════════════════════════════════════════════════════════════════════════
def build_results(spots, pidx, F68, F76, b68, b76, sd68, sd76, cfg: BatchConfig):
    """
    把所有校正过的比值转成最终结果表。

    每一行一个测点，包含：
        标识信息（序号/文件/样品/角色/剥蚀窗口）
        质控量（U 信号强度、Th/U、普通铅占比 f206）
        比值及其相对不确定度（206/238、207/206、207/235）
        三个年龄体系及其不确定度
        协和度  = 207Pb/235U 年龄 / 206Pb/238U 年龄 × 100%
    """
    rows = []
    for k, (tr, r) in enumerate(spots):
        f68, f76 = interp_F(k, pidx, F68, F76)

        # 校正后的比值
        R68 = b68[k] * f68
        R76 = b76[k] * f76
        # 不确定度：内部误差也随同一个 F 缩放；外部误差是相对量，乘校正后的比值
        s68 = float(np.hypot(r["s68"] * f68, sd68 * R68))
        s76 = float(np.hypot(r["s76"] * f76, sd76 * R76))

        # 207Pb/235U 由另外两个比值换算（不直接测）
        R75 = R76 * R68 * U238_U235

        a68 = float(age68(R68))
        # σ_age = σ_R / (λ·(1+R))，再 /1e6 换成 Ma。
        # λ 取自 core.constants，绝不在这里另写一个数：同一个衰变常数写两遍，
        # 改了一处忘了另一处，结果不报错、只是年龄不确定度悄悄错掉。
        sa68 = float(s68 / (L238 * (1 + R68)) / 1e6)
        a75 = float(age75(R75))
        a76 = float(age76(R76))
        # 207Pb/206Pb 的年龄方程没有解析导数，用有限差分近似 σ
        sa76 = float(abs(age76(R76 + s76) - age76(R76 - s76)) / 2.0)

        # 207Pb/235U 的年龄不确定度是**近似**，不是独立传播出来的：
        # 这里把 206Pb/238U 的相对不确定度原样搬到 207/235 上（两者共享同一个
        # 206Pb/238U 测量项），而不是从 R76 与 R68 的协方差严格传播。
        # 保持自 1.x 起的口径，以免已发表的数字发生无解释的变动；
        # 若要严格化，需要 ratios.reduce_interval 额外输出 R68-R76 的相关系数。
        sa75_approx = sa68 * (a75 / max(a68, 1e-9))

        rows.append(dict(
            序号=int(tr.idx) + 1,
            文件=Path(tr.path).stem,
            样品=tr.sample,
            类型=ROLE_LABEL_CN.get(tr.role, tr.role),
            剥蚀窗口_s=f"{tr.t0:.1f}-{tr.t1:.1f}",
            U238_cps=r["U_cps"],
            Th_U=r["ThU"],
            f206_pct=r["f206"] * 100,
            Pb206_238U=R68, s68_pct=s68 / R68 * 100,
            Pb207_206Pb=R76, s76_pct=s76 / R76 * 100,
            Pb207_235U=R75,
            年龄206_238=a68, s68_1sig=sa68, s68_2sig=2 * sa68,
            年龄207_235=a75, s75_2sig=2 * sa75_approx,
            年龄207_206=a76, s76_2sig=2 * sa76,
            协和度_pct=a75 / a68 * 100 if a68 > 0 else np.nan,
        ))
    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
# ⑨ 深度剖面年龄域判别
# ═════════════════════════════════════════════════════════════════════════════
def analyse_depth_spot(tr: Tra, r: dict, brack, ref68: float, sd68: float,
                       cfg: BatchConfig):
    """
    对单个样品点做深度剖面分析。

    返回
    ----
    (prof, segs, summ, tag)
        prof : 逐窗口还原结果（含 age68 / s_age68）
        segs : 年龄域下标范围
        summ : 域汇总表
        tag  : "均一" / "多域(n)"
    """
    # 逐窗口还原。这里必须传 fix=(f206, sk)：短窗口里 204 计数不足以独立判普通铅
    prof = window_profile(tr, cfg.win, cfg.step,
                          fix=(r["f206"], r["sk"]),
                          n_sigma=cfg.n_sigma_common_pb)
    if prof.empty or not brack:
        return None

    # ★ 关键：用"同一归一化深度"的分馏因子逐窗口校正 —— 不是整段单一因子！
    F = bracket_F(prof["tau"].to_numpy(), brack, "R68", ref68)
    ages, sigs = profile_ages(prof, F, sd68)
    prof["age68"], prof["s_age68"] = ages, sigs

    # BIC 二叉分割 → 边界精修（去混合窗口）→ 合并伪分域
    segs = segment(ages, sigs)
    segs = refine_domains(prof, segs) if segs else [(0, len(prof))]
    segs = merge_close(prof, segs, n_sigma=cfg.merge_n_sigma,
                       min_frac=cfg.merge_min_frac)
    summ = summarize_segments(prof, segs)

    nd = int((summ["flag"] == "age domain").sum())
    tag = "均一" if nd <= 1 else f"多域({nd})"
    return prof, segs, summ, tag


def run_depth_analysis(spots, brack_for, ref68, sd68, cfg: BatchConfig):
    """
    对全部样品点跑深度分析，顺便出图。

    返回
    ----
    (struct, domains, windows)
        struct  : 每个点的结构标签列表，如 "均一" / "多域(2)"，长度 = len(spots)
        domains : 多域点的各年龄域明细 DataFrame
        windows : 逐窗口年龄剖面，列名按 ADEPT 的 Format 4 对齐
    """
    struct: List[str] = []
    dom_rows: List[dict] = []
    win_rows: List[dict] = []
    plot_dir = cfg.plot_dir if cfg.plot else None
    pdf = None
    if cfg.plot:
        # 一边画一边写入多页 PDF 并立即关闭 Figure —— 不在内存里攒图
        from matplotlib.backends.backend_pdf import PdfPages
        plot_dir.mkdir(parents=True, exist_ok=True)
        pdf = PdfPages(plot_dir / f"{cfg.data_dir.name}_深度剖面_全部.pdf")

    for k, (tr, r) in enumerate(spots):
        # 只有未知样品需要做多域判别；标样理论上应该是均一的
        if tr.role != ROLE_UNKNOWN or not cfg.do_depth:
            struct.append("")
            continue
        try:
            out = analyse_depth_spot(tr, r, brack_for(k), ref68, sd68, cfg)
        except Exception as e:
            print(f"    !! {Path(tr.path).stem} 深度剖面失败: {e}")
            struct.append("")
            continue
        if out is None:
            struct.append("")
            continue

        prof, segs, summ, tag = out
        struct.append(tag)

        # 逐窗口剖面：全部窗口都收，不只是多域点 —— 均一的点同样需要一个
        # 带不确定度的坪年龄。列名在这里就对齐 ADEPT 的 Format 4，
        # 下游不需要任何转录或改名。
        for _, w in prof.iterrows():
            win_rows.append(dict(
                Analysis=f"{int(tr.idx) + 1:02d} {tr.sample}",
                Sample=tr.sample,
                Time=float(w["t_mid"]),
                Tau=float(w["tau"]),
                Age68=float(w["age68"]),
                Age68_1s=float(w["s_age68"]),
                R68=float(w["R68"]),
                s68=float(w["s68"]),
                ThU=float(w["ThU"]),
                U_cps=float(w["U_cps"]),
                f206_pct=float(w["f206"]) * 100,
                n_cycles=int(w["n_cycles"]),
            ))

        if cfg.plot:
            title = (f"{tr.label}   深度剖面 [{tag}]   "
                     f"窗口 {cfg.win}s / 步长 {cfg.step}s")
            from .depth.figures import save_depth_figure
            png = plot_dir / f"{int(tr.idx) + 1:02d}_{tr.sample}_{tag}.png"
            save_depth_figure(prof, segs, title, png, summ=summ,
                              with_207=False, pdf=pdf)

        # 多域的才写明细表，避免结果表里塞进一堆无信息的行
        nd = int((summ["flag"] == "age domain").sum())
        if nd > 1:
            for _, s in summ.iterrows():
                dom_rows.append(dict(
                    序号=int(tr.idx) + 1, 样品=tr.sample,
                    域=s["domain"], 标记=s["flag"],
                    tau=f"{s['tau0']:.2f}-{s['tau1']:.2f}",
                    年龄_Ma=s["age_Ma"], s2_Ma=2 * s["se_1sig"],
                    # MSWD 与它的卡方上尾概率。与 R 端 ADEPT 的
                    # `MSWD` / `MSWD probability` 同一口径 —— 两边都是
                    # 反比方差加权平均 + χ² 上尾，可以直接对照。
                    # （年龄_Ma 本身就是反比方差加权平均，不是算术平均。）
                    MSWD=s["mswd"], MSWD_概率=s["mswd_prob"],
                    Th_U=s["ThU"], n_win=int(s["n_win"])))

    if cfg.plot and pdf is not None:
        pdf.close()
        _log(cfg, f"    图件目录：{plot_dir}")
        _log(cfg, f"    汇总 PDF：{plot_dir / (cfg.data_dir.name + '_深度剖面_全部.pdf')}")

    # 防御：异常路径可能导致 struct 与结果行数不一致
    if len(struct) != len(spots):
        struct = struct[:len(spots)] + [""] * max(0, len(spots) - len(struct))
    return struct, pd.DataFrame(dom_rows), pd.DataFrame(win_rows)


# ═════════════════════════════════════════════════════════════════════════════
# ⑩ QC 汇总 + 监控标样二次校正
# ═════════════════════════════════════════════════════════════════════════════
def build_qc_table(res: pd.DataFrame, cfg: BatchConfig) -> pd.DataFrame:
    """对每个标样做加权平均，与其参考年龄比较，输出 QC 表。"""
    rows = []
    for name, grp in res.groupby("样品", sort=False):
        if grp["类型"].iloc[0] not in (ROLE_LABEL_CN[ROLE_PRIMARY],
                                       ROLE_LABEL_CN[ROLE_SECONDARY]):
            continue
        mu, se, mswd, n = weighted_mean(grp["年龄206_238"], grp["s68_1sig"])
        ref = std_age(std_alias(name)) if std_alias(name) else float("nan")
        rows.append(dict(
            标样=name, 点数=n, 参考年龄_Ma=ref,
            加权平均年龄_Ma=mu, s2_Ma=2 * se, MSWD=mswd,
            偏差_pct=(mu - ref) / ref * 100 if np.isfinite(ref) else np.nan))
    qc = pd.DataFrame(rows)
    _log(cfg, "\n[5] 标样 QC（206Pb/238U 加权平均）")
    for _, q in qc.iterrows():
        extra = (f"   参考 {q['参考年龄_Ma']:.2f}   偏差 {q['偏差_pct']:+.2f}%"
                 if np.isfinite(q["参考年龄_Ma"]) else "")
        _log(cfg, f"    {q['标样']:<10s} {q['加权平均年龄_Ma']:8.2f} ± {q['s2_Ma']:5.2f} Ma"
                  f" (2σ, n={int(q['点数'])}, MSWD={q['MSWD']:.2f}){extra}")
    return qc


def apply_secondary_correction(res: pd.DataFrame, cfg: BatchConfig):
    """
    用监控标样做**二次校正**（基体匹配标准化）。

    为什么需要它
    ------------
    本批次存在一个结构性问题：
        主标 91500 的 238U ≈ 1×10⁵ cps
        样品（示例样品）的 238U ≈ 1×10⁶ cps
    整整差一个数量级。在高计数率下，探测器/电子学的脉冲计数非线性
    无法被单一的线性归一化因子消除，表现为监控标样系统性偏老。

    用**计数率与样品相当**的监控标样做再校准，等价于"基体匹配"，
    是这种情况下行业内通行的补救办法。

    ⚠ 这是治标不治本。正式发表前应尽量：
        · 缩小标样与样品的信号强度差距（换小束斑 / 调低能量）；
        · 或施加可靠的死时间校正（见 core.deadtime）。

    做法
    ----
    系数 kfac = 参考年龄 / 监控标样实测加权平均年龄，
    把全部年龄乘以 kfac。所有 result 列的 σ 同步缩放。
    """
    alias = std_alias(cfg.secondary)
    ref_age = std_age(alias) if alias else float("nan")
    if not np.isfinite(ref_age) or cfg.secondary == cfg.primary:
        return res, None

    g = res[res["样品"] == cfg.secondary]
    if len(g) == 0:
        return res, None
    mu_qc, se_qc, mswd_qc, n_qc = weighted_mean(g["年龄206_238"], g["s68_1sig"])
    if not np.isfinite(mu_qc) or mu_qc <= 0:
        return res, None

    kfac = float(ref_age / mu_qc)
    res = res.copy()
    res["QC校正系数"] = kfac
    res["年龄206_238_QC校正"] = res["年龄206_238"] * kfac
    res["s68_2sig_QC校正"] = res["s68_2sig"] * kfac
    # 1σ 与 2σ 都按同一个系数缩放：kfac 是乘性因子，误差一同放大/缩小。
    # （s68_2sig 上一行已缩放，这里补 s68_1sig —— ADEPT 画加权平均图用的是 1σ。）
    res["s68_1sig_QC校正"] = res["s68_1sig"] * kfac

    _log(cfg, f"\n[6] QC 二次校正：{cfg.secondary} 实测 {mu_qc:.2f} Ma "
              f"(n={n_qc}, MSWD={mswd_qc:.2f})  vs 参考 {ref_age:.2f} Ma"
              f"  →  系数 {kfac:.5f}")
    _log(cfg, "    已写入列：年龄206_238_QC校正 / s68_2sig_QC校正 / s68_1sig_QC校正")
    _log(cfg, f"    成因：主标 {cfg.primary}(238U≈1e5 cps) 与样品(≈1e6 cps)"
              f"计数率相差一个数量级，单一归一化因子无法完全消除非线性残差。")
    return res, dict(sample=cfg.secondary, ref_age=ref_age, measured=mu_qc,
                     factor=kfac, n=n_qc, mswd=mswd_qc)


# ═════════════════════════════════════════════════════════════════════════════
# 对照方法列
# ═════════════════════════════════════════════════════════════════════════════
def _alternative_method_column(spots, variants, alt: str, ref68: float,
                               pmask, cfg: BatchConfig) -> List[float]:
    """
    用另一种整段比值方法算一列对照年龄。

    为什么值得单独一列
    ------------------
    两套方法在同一批数据上的差别本身就是诊断信息：若差别远大于它们各自的 σ，
    说明 down-hole 分馏严重，或存在别的结构性问题。所以默认那套之外再算一列，
    不必让用户重跑一遍去对照。

    三种算不出来的情形（整批无主标、主标在该方法下被稳健统计全部剔除）
    统一返回全 nan：宁可少一列，也不要让空数组进 robust_mask / np.interp
    触发无意义的警告或崩溃。

    ⚠ 本函数只做归一化，不做深度校正——"对照"的含义是比值方法不同，
      不是流程不同。与原实现逐位一致（同样的 robust_mask + np.interp 顺序）。
    """
    n = len(spots)
    a68b = variants[alt][0]
    if pmask.sum() == 0:
        _log(cfg, f"    （注：对照方法 {alt} 无有效主标，未生成对比列）")
        return [float("nan")] * n

    Fb = ref68 / a68b[pmask]
    okb = robust_mask(Fb)[0]
    pb = np.flatnonzero(pmask)[okb].astype(float)
    if len(pb) == 0:
        _log(cfg, f"    （注：对照方法 {alt} 无有效主标，未生成对比列）")
        return [float("nan")] * n

    return [float(age68(a68b[k] * float(np.interp(k, pb, Fb[okb]))))
            for k in range(n)]


# ═════════════════════════════════════════════════════════════════════════════
# 主入口
# ═════════════════════════════════════════════════════════════════════════════
def run_batch(cfg: BatchConfig) -> BatchResult:
    """
    执行完整的批处理流水线（上面十步全部在此串起来）。

    返回
    ----
    BatchResult，其中 results 是主结果表，qc 是标样质控表，
    domains 是多域剖面明细，info 记录了本次运行的全部关键量。
    """
    # 流水线全程用中文打印进度。Windows 上把输出重定向到文件时，
    # Python 会用 locale 编码（cp1252 之类）编码 stdout，第一句 print 就抛
    # UnicodeEncodeError —— 而且是在算完之后、准备写结果的时候。
    # 放在这里（而不是只在 CLI 入口）是因为作为库被调用时同样会踩到。
    ensure_utf8_streams()

    # ①② 装载
    spots, seq, skipped = load_batch(cfg)
    if not spots:
        raise RuntimeError(f"没有任何测点成功还原，检查数据目录：{cfg.data_dir}")

    # 只有主标/监控标样、没有未知样品的批次（纯校准/QC 批次）也允许处理：
    # 这类数据只产出标样 QC 表，深度域判别与二次校正自动跳过。
    n_unknown = sum(1 for tr, _ in spots if tr.role == ROLE_UNKNOWN)
    if n_unknown == 0:
        _log(cfg, "    ℹ 本批次没有未知样品，仅输出标样 QC（深度域判别、二次校正均跳过）。")

    # ── 主标判定：决定走"校准模式"还是"未校准诊断模式" ──
    # 只数主标：监控标样只用于 QC 与 σext，缺了不影响能否归一化。
    calibrated = any(tr.role == ROLE_PRIMARY for tr, _ in spots)

    if not calibrated:
        # ⚠ 批次里没有任何被识别为主标的测点（既没有 91500，
        # 也没在参数里把"主标名"设成实际使用的标样）。
        # 没有主标就构造不出归一化因子 F，绝对年龄无从谈起。
        # 这里**不崩溃**，而是退化为"未校准原始比值"模式：
        #   比值直接取实测值（F=1），年龄仅作相对参考、不可用于定年；
        #   深度剖面与二次校正因缺少标样自动跳过。
        # 同时给出醒目警告，提示用户补一个主标样后重跑才是有科学意义的结果。
        names = ", ".join(sorted({tr.sample for tr, _ in spots})) or "（空）"
        _log(cfg, "=" * 60)
        _log(cfg, "  ⚠ 警告：本批次未找到主标样（默认 91500）！")
        _log(cfg, "     已切换为【未校准 / 原始比值】诊断模式：")
        _log(cfg, "        · 比值 = 实测值，未做外标归一化；")
        _log(cfg, "        · 输出的年龄仅作相对参考，严禁直接用于定年；")
        _log(cfg, "        · 深度域判别与监控标样二次校正自动跳过。")
        _log(cfg, f"     批次中出现的样品名：{names}")
        _log(cfg, "     若要得到校准年龄：在序列表中加入主标样（如 91500），")
        _log(cfg, "     或在参数里把『主标名』设为你实际使用的标样名后重跑。")
        _log(cfg, "=" * 60)

        ref68 = ref76 = 1.0
        # 无夹逼关系；整段比值直接取实测值
        brack_for = lambda k: []
        b68 = np.array([r["R68"] for _, r in spots])
        b76 = np.array([r["R76"] for _, r in spots])
        variants = {"simple": (b68, b76), "ftau": (b68, b76)}
        pidx = np.array([], dtype=float)
        F68 = np.array([], dtype=float)
        F76 = np.array([], dtype=float)
        # 未校准时不确定度用默认外部重现性（仅作占位，量级仅供参考）
        sd68, sd76 = 0.007, 0.0025
        cal = dict(rejected=[], pmask=np.array([False] * len(spots)))
        uncal_note = ("未校准：批次无主标样，比值为实测值、年龄仅供参考，"
                      "不可用于定年。建议加入主标样后重跑。")
    else:
        try:
            ref68, ref76 = std_ref(std_alias(cfg.primary) or cfg.primary)
        except KeyError:
            # 用户把"主标名"设成了一个标准库里没有的标样。
            # 没法查到它的参考比值，无法做外标归一化——明确报错而不是事后给错年龄。
            known = ", ".join(sorted(STANDARDS))
            raise RuntimeError(
                f"主标样『{cfg.primary}』不在参考标准库内，无法查到参考比值做归一化。"
                f"可用主标：{known}。请在参数里改用已知标样名，"
                f"或先在 core.constants.STANDARDS 中添加该标样。")

        # ③④ 主标剖面 + 夹逼关系
        brack_for = build_bracketing(spots, cfg)

        # ⑤ 两套整段比值，选启用的一套
        variants, (b68, b76) = compute_bulk_ratios(spots, brack_for, ref68, cfg)

        # ⑥ 归一化因子曲线
        pidx, F68, F76, cal = calibrate_primary(spots, b68, b76, ref68, ref76, cfg)

        # ⑦ 外部重现性
        sd68, sd76 = estimate_external(spots, pidx, F68, F76, b68, b76, cfg)
        uncal_note = None

    # ⑧ 结果表
    res = build_results(spots, pidx, F68, F76, b68, b76, sd68, sd76, cfg)
    # 在校准状态列里明确标注，避免把未校准年龄误当定年结果
    res["校准状态"] = "未校准" if not calibrated else "已校准"

    # 另一套方法的结果并列输出，便于对照两者的系统性差异
    alt = "ftau" if cfg.bulk == "simple" else "simple"
    res[f"年龄206_238_{alt}"] = _alternative_method_column(
        spots, variants, alt, ref68, cal["pmask"], cfg)

    # ⑨ 深度剖面
    struct, domains, windows = run_depth_analysis(spots, brack_for, ref68, sd68, cfg)
    res["深度结构"] = struct
    ns = int(res["深度结构"].str.startswith("多域").sum())
    _log(cfg, f"\n[4] 深度剖面结构判别：{ns} / "
              f"{int((res['类型'] == ROLE_LABEL_CN[ROLE_UNKNOWN]).sum())} 个样品测点检出多年龄域")

    # ⑩ QC + 二次校正
    qc = build_qc_table(res, cfg)
    res, corr = apply_secondary_correction(res, cfg)

    info = dict(
        data_dir=str(cfg.data_dir),
        primary=cfg.primary, secondary=cfg.secondary,
        ref68=ref68, ref76=ref76,
        calibrated=calibrated,
        bulk=cfg.bulk, win=cfg.win, step=cfg.step, trim=cfg.trim,
        deadtime_ns=cfg.deadtime_ns,
        sd68=sd68, sd76=sd76,
        rejected_primary=cal["rejected"],
        skipped=skipped,
        secondary_correction=corr,
        warning=uncal_note,
    )
    return BatchResult(results=res, qc=qc, domains=domains, info=info,
                       spots=spots, windows=windows)
