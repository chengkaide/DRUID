# -*- coding: utf-8 -*-
"""
tools/synth_profile.py —— 两域深度剖面正演：年龄差 Δt → P(正确检出)
======================================================================

要回答的问题
------------
`tools/synth_batch.py` 造的是**均一**剖面（一个测点一个真年龄），所以它证明不了
任何与"分域"有关的事。这个脚本把它推到**一个测点内部有两个年龄域**：

    浅处（坑口）一个域，深处（坑底）另一个域，界在归一化深度 τ = β 处。

然后逐项加上真实数据里躲不掉的效应，扫描**两域的年龄差 Δt**，最后画

    P(正确检出两个域) —— Δt

这条曲线。它回答的是"**要多大年龄差才看得见**"，而不是"算法算得准不准"。

加进去的效应（每一步都可以单独关掉，见 `--` 参数）
--------------------------------------------------
    · Pb 计数噪声        —— 逐行 Poisson（把 cps 当「计数 ÷ dwell」），**始终在**
    · 1% 随机噪声        —— 逐通道乘性高斯（等离子体白噪声的粗模型）`--noise`
    · 2% 校准漂移        —— 全批次的慢漂移（Pb/U 分馏因子），`--drift`
    · 5% 深度分馏 (DHF)  —— 平台段内 Pb 随 τ 的线性变化，`--dhf`
    · 普通铅             —— 按 f206 = 206Pb_common/206Pb_total 加进样品，`--common-pb`

⚠ 一条**必须知道**的结构性事实（本脚本会把它的后果画出来）
--------------------------------------------------------
`depth/domains.py: merge_close` 判"这两个域是不是同一个"时，同时要求

    |Δ年龄| > 3σ        （统计显著）      **并且**
    |Δ年龄| > 5% × 年龄  （地质显著）

在 450 Ma 上第二条 = **22.5 Ma**。也就是说：**只要 Δt < 22.5 Ma，不论数据多干净，
工具都必然把它们合并成一个域** —— 这不是统计分辨力不够，是软件的政策。
所以本脚本默认同时跑一条 `--min-frac 0` 的对照曲线（只留 3σ 统计门槛），
两条曲线之差恰好就是**"算法政策"从"数据信息量"手里拿走的那一段**。

怎么判"正确检出"
----------------
一个样品测点算**正确检出**，要同时满足三件事：

    ① 结构标签是"多域(2)"（不是均一、也不是多域(3)）          —— 数对
    ② 两个域的年龄各自落在真值的 3σ 内（σ 取工具自己报的那一个）—— 年龄对
    ③ 域界位置与真值 β 相差 ≤ 0.25 τ                          —— 边界对

②用**工具自己报的 σ**做容差，是为了让"报回得准"这件事由工具自己背书：
它说 ±2 Ma 却偏了 10 Ma，那就不算对。另有 2 Ma 的容差地板（域 σ 极小时用）。

三条各自单独统计（`P(2域)` / `P(年龄对)` / `P(严格)`），因为它们的缺口
指向完全不同的原因。

不是从零造轮子
--------------
    · 采集形态、CSV 形态、噪声采样、写盘 —— 全部复用 `tools/synth_batch.py`
      （常量与 `spot_csv_text` / `_rmtree` 直接 import，不复制代码）。
    · 正演比值 —— `core.geochronology.r68_of_age / r76_of_age`。
    · 普通铅成分 —— `core.common_lead.stacey_kramers`（**工具自己那一份**）。
    · 参考值 —— `core.references.std_ref`。

⚠ 由此产生的**已知局限**（写在这里，免得被当成"已经验证过"）
    · 普通铅用的是**工具自己的模型** ⇒ 这是"模型匹配"检验：它检验的是
      **统计链条**（204 的显著性判定、扣除的传播、噪声下的稳定性），
      **不是**"换一个普通铅模型会怎样"。
    · DHF 只改 Pb 三个通道（206/207/208），U/Th 不动 ⇒ "206/238 随深度变化
      dhf"这句话是逐字精确的；但它不含质量依赖分馏（207/206 不变）。
    · 校准漂移是**整批平滑**的，而标样排得很密（每个样品后跟一对）⇒
      夹逼线性插值几乎能全部吸收。稀疏标样下的漂移是**另一件事**，本脚本答不了。
    · 域界是**行级的阶跃**（再被 4 s 窗天然平滑成 ~4 s 过渡）——不是真实的
      坑底混合几何，也不含"过渡带正好跨在窗口上"的随机性。
    · 只判 **206Pb/238U** 这一个钟。

用法
----
    python tools/synth_profile.py --pilot              # 1 个种子、少量 Δt，先量时间
    python tools/synth_profile.py                      # 默认全跑（12 种子 × 6 条件）
    python tools/synth_profile.py --dt 5,10,20,40,80   # 自定 Δt 网格
    python tools/synth_profile.py --seeds 24 --repeat 4
    python tools/synth_profile.py --plot-only path.csv # 只重画，不重算
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "tools") not in sys.path:
    sys.path.insert(0, str(ROOT / "tools"))

import synth_batch as SB                                             # noqa: E402
from druid.console import ensure_utf8_streams                        # noqa: E402
from druid.core.common_lead import stacey_kramers                    # noqa: E402
from druid.core.constants import (CJK_FONTS, DEFAULT_REF_PRESET,     # noqa: E402
                                  L232, MASSES_NEEDED)
from druid.core.geochronology import r68_of_age, r76_of_age          # noqa: E402
from druid.core.references import std_alias, std_ref                 # noqa: E402
from druid.workflow import BatchConfig, run_batch                    # noqa: E402


# 采集形态与空白全部沿用 synth_batch：两边的测点必须**同构**，
# 否则"分域"这一步的输入与之前那批均一剖面的输入不是一个东西，无法对照。
N_BLANK, N_PLATEAU, N_ROWS = SB.N_BLANK, SB.N_PLATEAU, SB.N_ROWS
TAIL_FRAC = SB.TAIL_FRAC
BG = SB.BG
PRIMARY, SECONDARY = SB.PRIMARY, SB.SECONDARY

TOP_DEFAULT = 450.0        # 浅处（坑口）那一域的年龄
BETA_DEFAULT = 0.5         # 域界所在的归一化深度


# ═════════════════════════════════════════════════════════════════════════════
# 一、正演：逐行（逐深度）的期望信号
# ═════════════════════════════════════════════════════════════════════════════
def domain_rows(top, bot, beta, n_plateau=N_PLATEAU):
    """
    逐行的真值比值 —— 前 round(β·n) 行是 `top` 年龄，之后是 `bot` 年龄。

    返回 (age_row, r68_row, r76_row, r82_row)，长度都是 n_plateau。
    `age_row` 是**逐行的真值年龄**，加普通铅时按行取成分要用它。
    """
    k = int(round(beta * n_plateau))
    k = max(1, min(n_plateau - 1, k))        # 两侧各至少留 1 行，否则不成"两域"
    idx = np.arange(n_plateau)
    age_row = np.where(idx < k, float(top), float(bot))

    def rows_of(fn):
        return np.where(idx < k, float(fn(top)), float(fn(bot)))

    return (age_row, rows_of(r68_of_age), rows_of(r76_of_age),
            rows_of(lambda a: float(np.expm1(L232 * a * 1e6))))


def profile_levels(r68_row, r76_row, r82_row, u_cps, th_u, dhf=0.0):
    """
    三阶段（空白 / 剥蚀平台 / 冲洗拖尾）的**期望** cps，平台段**逐行不同**。

    与 `synth_batch.spot_levels` 的唯一区别：平台段不是常数，而是一条
    "阶跃 + 线性分馏"的曲线。空白与拖尾的处理**完全照抄**（拖尾取平台
    **末行**净峰高的 TAIL_FRAC），这样 `find_ablation` 遇到的两级阈值
    形态与上一批一字不差。

    DHF 只乘 Pb 三个通道（206/207/208），U 与 Th 不动 —— 于是
    "206/238 在平台段内线性变化 dhf"这句话在下游是**逐字精确**的。
    φ(τ) = 1 + dhf·(τ − 0.5)，τ 从平台首行 0 到末行 1 ⇒ 首末相差 dhf。
    """
    tau_plateau = np.linspace(0.0, 1.0, N_PLATEAU)
    phi = 1.0 + dhf * (tau_plateau - 0.5)
    i238 = float(u_cps)
    i232 = th_u * i238

    peak = {
        238: np.full(N_PLATEAU, i238),
        232: np.full(N_PLATEAU, i232),
        206: np.asarray(r68_row, float) * i238 * phi,
        207: np.asarray(r76_row, float) * np.asarray(r68_row, float) * i238 * phi,
        208: np.asarray(r82_row, float) * i232 * phi,
        202: np.zeros(N_PLATEAU),        # Hg：不随剥蚀变化（净 0）
        204: np.zeros(N_PLATEAU),        # 普通铅：默认净 0
    }
    levels = {}
    for m in MASSES_NEEDED:
        net = np.zeros(N_ROWS, dtype=float)
        net[N_BLANK:N_BLANK + N_PLATEAU] = peak[m]
        net[N_BLANK + N_PLATEAU:] = TAIL_FRAC * float(peak[m][-1])
        levels[m] = BG[m] + net
    return levels


def add_common_pb(levels, f206, age_row):
    """
    往**样品**里加普通铅：f206 = 206Pb_common / 206Pb_total（按质量计）。

    做法（与 `reduce_interval` 的扣除是互逆的）：
        rad206 = 净 206 信号                （已经只含放射成因）
        com206 = f206/(1−f206) · rad206     ⇒ com206/(com206+rad206) = f206
        其余按 `stacey_kramers(该行真值年龄)` 的成分一起补进 204/207/208。

    ⚠ 用**工具自己的**普通铅模型 ⇒ 这是模型匹配检验，见模块 docstring 的局限。
    """
    if f206 <= 0.0:
        return levels
    # 把平台段的**逐行年龄**铺满整条曲线，成分才与 rad206 同长：
    # 空白段 rad206 = 0（成分取什么都不影响结果）；拖尾段用末行的年龄。
    age_full = np.empty(N_ROWS, dtype=float)
    age_full[:N_BLANK] = age_row[0]
    age_full[N_BLANK:N_BLANK + N_PLATEAU] = age_row
    age_full[N_BLANK + N_PLATEAU:] = age_row[-1]
    rad206 = levels[206] - BG[206]
    com206 = f206 / (1.0 - f206) * rad206
    p6, p7, p8 = stacey_kramers(age_full)         # 206/204, 207/204, 208/204
    com204 = com206 / p6
    levels[206] = levels[206] + com206
    levels[204] = levels[204] + com204
    levels[207] = levels[207] + com204 * p7
    levels[208] = levels[208] + com204 * p8
    return levels


def scale_pb(levels, g):
    """把 Pb 三个通道的**净**信号乘以 g（Pb/U 分馏/灵敏度的漂移），U 与 Th 不动。"""
    for m in (206, 207, 208):
        levels[m] = BG[m] + (levels[m] - BG[m]) * g
    return levels


# ═════════════════════════════════════════════════════════════════════════════
# 二、造整批（写盘）
# ═════════════════════════════════════════════════════════════════════════════
def fast_rmtree(path, root=None):
    """
    **真删除**（不进回收站），只用于自己刚写下的临时目录。

    为什么非要绕开 `pathlib`
    ------------------------
    这台机器上 WorkBuddy 的 `sitecustomize.py` 把 `pathlib.Path.unlink`、
    `os.remove`、`os.rmdir`、`shutil.rmtree` **全部改道 Windows 回收站**
    （`SHFileOperationW`）。实测同一批 84 个 17 KB 的文件：

        · `pathlib.unlink` / `os.remove`  →  34.0 s   （≈0.4 s/个）
        · 原生 `kernel32.DeleteFileW`      →  0.018 s  （同一条命令，1900 倍）

    而这个脚本一跑就是几百个批次、上万个测点文件 ⇒ 用垫片那条路，
    整个实验会从 4 分钟变成 1 小时以上，且**时间全花在把临时文件挪进回收站**。

    ⚠ 安全阀：给定 `root` 时，`path` 必须落在它之下，否则直接拒绝 ——
    免得调用方传错路径删到别的东西。
    """
    p = Path(path)
    if not p.exists():
        return
    if root is not None:
        try:
            p.resolve().relative_to(Path(root).resolve())
        except ValueError:
            raise ValueError(f"拒绝删除 {p}：不在允许的根 {root} 之下")
    if sys.platform == "win32":
        try:
            import ctypes
            k32 = ctypes.WinDLL("kernel32", use_last_error=True)
            for fn, arg in (("DeleteFileW", ctypes.c_wchar_p),
                            ("RemoveDirectoryW", ctypes.c_wchar_p)):
                getattr(k32, fn).argtypes = [arg]
                getattr(k32, fn).restype = ctypes.c_int
            # 深的先删（目录必须在它里面的东西之后）
            for f in sorted(p.rglob("*"), key=lambda x: len(str(x)), reverse=True):
                if f.is_dir():
                    k32.RemoveDirectoryW(str(f))
                else:
                    k32.DeleteFileW(str(f))
            k32.RemoveDirectoryW(str(p))
            return
        except Exception:
            pass
    SB._rmtree(p)          # 非 Windows / 原生调用失败 → 退回逐文件删（慢但安全）


def _rmdir_if_empty(p):
    """只删**空**目录（不递归）；非空就原样留着。"""
    try:
        d = Path(p)
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()
    except OSError:
        pass


def make_profile_batch(out_dir, batch, dts, top, beta, repeat, u_cps, th_u,
                       seed, ref_preset, flicker, drift, dhf, dhf_std, f206):
    """
    一个批次里**放齐全部 Δt**（共用同一批夹逼标样）—— 一趟流水线就能拿到
    整条曲线，比"每个 Δt 起一次批"快一个量级，而且各 Δt 处在**同一套**
    漂移/噪声条件下，可比。

    返回 (batch_dir, plan)；plan 里样品项带 `name` / `top` / `bot` / `dt`。
    """
    d = Path(out_dir) / batch
    fast_rmtree(d, root=out_dir)
    d.mkdir(parents=True)

    rng = np.random.default_rng(seed)
    ref68_p, ref76_p = std_ref(std_alias(PRIMARY) or PRIMARY, ref_preset)
    ref68_s, ref76_s = std_ref(std_alias(SECONDARY) or SECONDARY, ref_preset)
    r82_std = SB.forward_ratios(985.0)[2]
    dhf_std = dhf if dhf_std is None else dhf_std

    # 排布与 synth_batch 同构：开头两主标 + 一监控打底，之后**每个样品测点后面
    # 跟一对 [主标, 监控标样]**（标样密度是闭合检验的前提，见 synth_batch 的注释）。
    plan = [dict(kind="std", sample=PRIMARY),
            dict(kind="std", sample=PRIMARY),
            dict(kind="std", sample=SECONDARY)]
    for dt in dts:
        bot = float(top) - float(dt)
        for k in range(repeat):
            plan.append(dict(kind="sample", sample=f"SYN-{top:g}-{bot:g}_{k + 1}",
                             top=float(top), bot=bot, dt=float(dt)))
            plan.append(dict(kind="std", sample=PRIMARY))
            plan.append(dict(kind="std", sample=SECONDARY))

    n = len(plan)
    entries = []
    for i, item in enumerate(plan):
        # 校准漂移：在整批序列上**线性**扫过 ±drift/2（峰-峰 = drift）
        pos = 0.0 if n <= 1 else (2.0 * i / (n - 1) - 1.0)
        g = 1.0 + 0.5 * float(drift) * pos
        fname = f"{batch}_{i + 1}"

        if item["kind"] == "std":
            if item["sample"] == PRIMARY:
                r68, r76 = ref68_p, ref76_p
            else:
                r68, r76 = ref68_s, ref76_s
            levels = profile_levels(np.full(N_PLATEAU, r68),
                                    np.full(N_PLATEAU, r76),
                                    np.full(N_PLATEAU, r82_std),
                                    u_cps, th_u, dhf=dhf_std)
        else:
            age_row, r68r, r76r, r82r = domain_rows(item["top"], item["bot"], beta)
            levels = profile_levels(r68r, r76r, r82r, u_cps, th_u, dhf=dhf)
            levels = add_common_pb(levels, f206, age_row)

        levels = scale_pb(levels, g)
        text = SB.spot_csv_text(item["sample"], levels, rng, flicker)
        (d / f"{fname}.csv").write_bytes(text.replace("\n", "\r\n").encode("utf-8"))
        item["file"] = fname
        entries.append((fname, item["sample"]))

    (d / f"{batch}_LIST.csv").write_text(
        "".join(f"{f},{s}\n" for f, s in entries), encoding="utf-8")
    return d, plan


# ═════════════════════════════════════════════════════════════════════════════
# 三、判"正确检出"
# ═════════════════════════════════════════════════════════════════════════════
def _ndet(struct):
    """从「深度结构」标签取域数："" = 没做深度分析；"均一" = 1；"多域(n)" = n。"""
    if not isinstance(struct, str) or not struct:
        return 0
    if struct.startswith("多域"):
        try:
            return int(struct[struct.index("(") + 1:struct.index(")")])
        except (ValueError, IndexError):
            return 0
    return 1 if struct == "均一" else 0


def _tau_pair(text):
    """解析「深度剖面域」里的 tau 字符串 "0.05-0.45" → (0.05, 0.45)。"""
    parts = [p for p in str(text).split("-") if p.strip()]
    if len(parts) < 2:
        return (float("nan"), float("nan"))
    return (float(parts[0]), float(parts[-1]))


def verdict(struct, dom, true_ages, beta, nsigma=3.0, tol_floor=2.0,
            bnd_tol=0.25):
    """
    把一个测点的深度判别结果判成四个互斥结局之一。

    返回 dict：ndet / ok2 / okage / okbnd / okstrict / order_ok / bnd_err /
    ages / ses / 结局。
    """
    nd = _ndet(struct)
    real = dom[dom["标记"] == "age domain"] if (dom is not None and len(dom)) \
        else pd.DataFrame()
    ages = [float(x) for x in real["年龄_Ma"]] if len(real) else []
    ses = [float(x) / 2.0 for x in real["s2_Ma"]] if len(real) else []

    ok2 = (nd == 2)
    okage = False
    order_ok = False
    bnd_err = float("nan")
    if ok2 and len(ages) == 2 and len(ses) == 2:
        det = sorted(zip(ages, ses))
        for (a, s), t in zip(det, sorted(true_ages)):
            if abs(a - t) > max(nsigma * s, tol_floor):
                break
        else:
            okage = True
        tau = [_tau_pair(x) for x in real["tau"]]
        if len(tau) == 2 and np.isfinite(tau[0][1]) and np.isfinite(tau[1][0]):
            bnd_err = 0.5 * (tau[0][1] + tau[1][0]) - beta
        order_ok = ages[0] > ages[1]        # 核老边新：τ 小的那一域应更老
    okbnd = bool(np.isfinite(bnd_err) and abs(bnd_err) <= bnd_tol)
    okstrict = bool(okage and okbnd and order_ok)
    # 检出的两个域**之间**还剩多少年龄差。真值 Δt 与它的比 = "跨度保留比"：
    # <1 说明域被向中间拉（过渡带混合 + 窗口只覆盖 τ 的一部分），
    # 这是**与检出与否无关**的一个独立退化量，必须单独看。
    dt_det = abs(ages[0] - ages[1]) if (ok2 and len(ages) == 2) else float("nan")

    if nd == 0:
        tag = "无结果"
    elif nd == 1:
        tag = "漏检"
    elif nd == 2:
        tag = "正确" if okstrict else ("年龄偏" if not okage else "边界偏")
    else:
        tag = "过分割"
    return dict(ndet=nd, ok2=ok2, okage=bool(okage), okbnd=okbnd,
                okstrict=okstrict, order_ok=bool(order_ok),
                bnd_err=bnd_err, ages=ages, ses=ses, dt_det=dt_det, 结局=tag)


# ═════════════════════════════════════════════════════════════════════════════
# 四、条件（"逐项加入效应"的阶梯）
# ═════════════════════════════════════════════════════════════════════════════
def base_conditions():
    """
    阶梯（①②③④⑤⑥）：每一项都**只加不减**，相邻两条之差就是**那一项的代价**。
    随后三条是**对照**，各自回答一个单独的问题，不参与"阶梯"的读法：

        ⑤*  把 5% 地质门槛归零  → 检出上限里有多少是**软件政策**而非数据信息量
        ⑤†  深度分馏**不匹配**（标样 0%）→ F(τ) 这个乘法假设错了会怎样
        ⑤‡  普通铅降到 1%       → 204 逼近显著性门槛时是什么形态
    """
    e = dict(flicker=0.0, drift=0.0, dhf=0.0, dhf_std=None, f206=0.0)
    full = dict(e, flicker=0.01, drift=0.02, dhf=0.05, f206=0.05)
    return [
        dict(name="① 只有计数统计", min_frac=0.05, **e),
        dict(name="② ＋1% 随机噪声", min_frac=0.05, **dict(e, flicker=0.01)),
        dict(name="③ ＋2% 校准漂移", min_frac=0.05,
             **dict(e, flicker=0.01, drift=0.02)),
        dict(name="④ ＋5% 深度分馏", min_frac=0.05,
             **dict(e, flicker=0.01, drift=0.02, dhf=0.05)),
        dict(name="⑤ 全套（＋普通铅 5%）", min_frac=0.05, **full),
        dict(name="⑥ 全套＋普通铅 1%（温和）", min_frac=0.05,
             **dict(full, f206=0.01)),
        dict(name="⑤* 全套但去掉 5% 地质门槛", min_frac=0.0, **full),
        dict(name="⑤† 全套但深度分馏不匹配", min_frac=0.05,
             **dict(full, dhf_std=0.0)),
    ]


# 合并判据的默认值（与 BatchConfig 一致）。条件里不带这两个键时用它们。
MERGE_N_SIGMA_DEFAULT = 3.0


def roc_conditions(min_fracs=(0.0, 0.01, 0.03, 0.05),
                   n_sigmas=(2.0, 3.0, 4.0)):
    """
    取舍曲线用的条件集：把 `merge_close` 的两条腿（统计门槛 `n_sigma` ×
    地质门槛 `min_frac`）逐一组合，其余效应固定为"全套"。

    这样曲线上每一点对应**同一批数据**，差别只在合并判据 —— 才是一条干净的
    ROC：横轴是**假阳性率**（Δt=0，即完全均一的剖面被判成两域的比例），
    纵轴是**灵敏度**（真两域被正确检出的比例）。`min_frac=0` 意味着只留统计腿。
    """
    full = dict(flicker=0.01, drift=0.02, dhf=0.05, dhf_std=None, f206=0.05)
    return [dict(name=f"mf={mf:g}|ns={ns:g}", min_frac=float(mf),
                 merge_n_sigma=float(ns), **full)
            for mf in min_fracs for ns in n_sigmas]


# ═════════════════════════════════════════════════════════════════════════════
# 五、主循环
# ═════════════════════════════════════════════════════════════════════════════
def run_all(dts, conds, seeds, repeat, tops, betas, u_cps, th_u, ref_preset,
            win, step, out_root, keep, verbose, dt_rel=None):
    """
    在 (核年龄 top × 域界 beta × 条件 × 种子) 上循环，逐测点判读，返回长表。

    每个批次里放齐该 (top, beta) 的全部 Δt —— 一趟流水线出一条完整曲线，
    且各 Δt 共用**同一套**夹逼标样与**同一条**漂移实现，彼此可比。

    `dt_rel` 给的是**相对**网格（占核年龄的百分比）；给了它就忽略 `dts`，
    每个 top 各按自己的比例换算成绝对 Δt —— 只有这样才能回答"相对跨度相同的
    两个年龄，检出难度是否相同"（绝对 Δt 网格天然把老年份的 Δt/age 压得很小）。
    """
    rows = []
    for ti, top in enumerate(tops):
        top = float(top)
        dts_i = [float(r) / 100.0 * top for r in dt_rel] if dt_rel else list(dts)
        for bi, beta in enumerate(betas):
            beta = float(beta)
            for ci, c in enumerate(conds):
                for s in range(seeds):
                    seed = 20261002 + 1000 * ci + s + 7919 * bi
                    sub = out_root / f".t{ti}_b{bi}_c{ci}_s{s}"
                    d, plan = make_profile_batch(
                        sub, "SYNPROF", dts_i, top, beta, repeat, u_cps, th_u,
                        seed, ref_preset, c["flicker"], c["drift"], c["dhf"],
                        c["dhf_std"], c["f206"])
                    cfg = BatchConfig(data_dir=str(d), ref_preset=ref_preset,
                                      do_depth=True, verbose=bool(verbose),
                                      plot=False, win=win, step=step,
                                      merge_min_frac=c["min_frac"],
                                      merge_n_sigma=c.get("merge_n_sigma", 3.0))
                    br = run_batch(cfg)
                    struct_by = {str(r["样品"]): str(r["深度结构"])
                                 for _, r in br.results.iterrows()}
                    dom = br.domains
                    for item in plan:
                        if item["kind"] != "sample":
                            continue
                        name = item["sample"]
                        one = dom[dom["样品"] == name] \
                            if (dom is not None and len(dom)) else pd.DataFrame()
                        v = verdict(struct_by.get(name, ""), one,
                                    (item["top"], item["bot"]), beta)
                        v.update(条件=c["name"], dt=item["dt"], 种子=seed,
                                 样品=name, min_frac=c["min_frac"],
                                 merge_n_sigma=float(c.get("merge_n_sigma", 3.0)),
                                 top=top, beta=beta,
                                 dt_rel=(item["dt"] / top * 100.0) if top else
                                 float("nan"))
                        rows.append(v)
                    if not keep:
                        fast_rmtree(sub, root=out_root)
    return pd.DataFrame(rows)


def summarize(det: pd.DataFrame) -> pd.DataFrame:
    """逐（核年龄 × 域界 × 条件 × Δt）算概率，外加两个与检出无关的退化量。"""
    det = det.copy()
    for col, val in (("top", float("nan")), ("beta", float("nan"))):
        if col not in det.columns:
            det[col] = val
    out = []
    for (top, beta, cond, dt), g in det.groupby(["top", "beta", "条件", "dt"],
                                                sort=False):
        n = len(g)
        two = g[g["ndet"] == 2]
        ratio = (two["dt_det"] / float(dt)).replace([np.inf, -np.inf], np.nan).dropna()
        be = g["bnd_err"].dropna().abs()
        out.append(dict(top=float(top), beta=float(beta), 条件=cond, dt=float(dt),
                        dt_rel=(float(dt) / float(top) * 100.0) if top else
                        float("nan"), n=n,
                        P_2域=float((g["ndet"] == 2).mean()),
                        P_年龄=float(g["okage"].mean()),
                        P_边界=float(g["okbnd"].mean()),
                        P_严格=float(g["okstrict"].mean()),
                        P_漏检=float((g["ndet"] == 1).mean()),
                        P_过分割=float((g["ndet"] >= 3).mean()),
                        跨度保留=float(ratio.median()) if len(ratio) else float("nan"),
                        边界误差中位=float(be.median()) if len(be) else float("nan")))
    return pd.DataFrame(out)


# ═════════════════════════════════════════════════════════════════════════════
# 六、出图
# ═════════════════════════════════════════════════════════════════════════════
def _is_control(name):
    """对照曲线（不参与"阶梯"读法）：去门槛 / 分馏不匹配 / 温和普通铅。"""
    return any(k in str(name) for k in ("⑤*", "⑤†", "⑥"))


def make_figure(sum_df, out_png, top, beta, win, step, min_frac, ycol="P_严格"):
    """
    两张子图：
        A 阶梯 —— 逐项加入效应，相邻两条之差 = 那一项的代价
        B 对照 —— 三条对照各自回答一个单独的问题
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = CJK_FONTS
    plt.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), dpi=150, sharey=True)
    order = list(dict.fromkeys(map(str, sum_df["条件"])))
    panels = [("A 阶梯：逐项加入效应",
               [c for c in order if not _is_control(c)]),
              ("B 对照：各问一个问题",
               [c for c in order if _is_control(c)])]
    floor_ma = float(min_frac) * float(top)

    for ax, (title, names) in zip(axes, panels):
        for cond in names:
            g = sum_df[sum_df["条件"] == cond].sort_values("dt")
            if not len(g):
                continue
            ax.plot(g["dt"], g[ycol], marker="o", ms=5, lw=1.6,
                    label=f"{cond}（n={int(g['n'].sum())}）")
        if floor_ma > 0:
            ax.axvline(floor_ma, color="0.35", ls="--", lw=1.1)
            ax.annotate(f"5% 地质门槛\n{floor_ma:.1f} Ma", xy=(floor_ma, 0.04),
                        xytext=(floor_ma + 2, 0.14), fontsize=8, color="0.25",
                        arrowprops=dict(arrowstyle="->", color="0.5", lw=1))
        ax.set_xlabel("两个年龄域的年龄差 Δt (Ma)", fontsize=10.5)
        ax.set_title(title, fontsize=11)
        ax.set_ylim(-0.04, 1.05)
        ax.grid(alpha=0.30)
        ax.legend(fontsize=8, loc="upper left", framealpha=0.92)

    axes[0].set_ylabel("P（正确检出两个域）", fontsize=11)
    fig.suptitle(f"两域深度剖面的检出概率：核 {top:g} Ma，边 {top:g}−Δt，域界 τ={beta:g}"
                 f"　｜　窗口 {win:g} s / 步长 {step:g} s"
                 f"　｜　判据＝数对＋年龄在 3σ 内＋边界误差 ≤ 0.25 τ",
                 fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png)
    plt.close(fig)
    return out_png


def make_age_figure(sum_df, out_png, cond, win, step, ycol="P_严格"):
    """
    年龄扫描图：**同一个条件**、不同的核年龄，横轴分别是绝对 Δt 与相对 Δt/年龄。

    左边按 Ma 看"绝对差要多大才检得出"，右边按 % 看"标度是否唯一"——
    若年龄绝对值只有**标度**上的影响，右图应当塌缩成一条线；塌缩不了，
    说明还有与年龄有关的非标度效应（207Pb 计数、门槛的相对大小、灵敏度曲线）。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = CJK_FONTS
    plt.rcParams["axes.unicode_minus"] = False

    sub = sum_df[sum_df["条件"] == cond] if cond else sum_df
    tops = sorted(float(t) for t in sub["top"].dropna().unique())
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), dpi=150)
    for ax, xcol, xlabel in ((axes[0], "dt", "两个年龄域的年龄差 Δt (Ma)"),
                             (axes[1], "dt_rel", "相对年龄差 Δt / 核年龄 (%)")):
        for t in tops:
            g = sub[sub["top"] == t].sort_values(xcol)
            ax.plot(g[xcol], g[ycol], marker="o", ms=4.5, lw=1.6,
                    label=f"{t:g} Ma（n={int(g['n'].sum())}）")
        ax.set_xlabel(xlabel, fontsize=10.5)
        ax.set_ylim(-0.04, 1.05)
        ax.grid(alpha=0.30)
        ax.legend(fontsize=8.5, loc="upper left", framealpha=0.92)
    axes[0].set_ylabel("P（正确检出两个域）", fontsize=11)
    fig.suptitle(f"检出概率的年龄依赖：条件「{cond}」　｜　窗口 {win:g} s / 步长 {step:g} s\n"
                 f"判据＝数对＋两域年龄各在 max(3σ, 2 Ma) 内＋边界 |Δτ| ≤ 0.25"
                 f"（左＝绝对刻度，右＝相对刻度）", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out_png)
    plt.close(fig)
    return out_png


# ═════════════════════════════════════════════════════════════════════════════
# 七、命令行
# ═════════════════════════════════════════════════════════════════════════════
def _default_out() -> Path:
    import os
    env = os.environ.get("DRUID_SYNTH_DIR")
    if env:
        return Path(env)
    g = Path("G:/")
    if g.exists():
        return g / "_zd_run" / "synth_profile"
    return Path(tempfile.gettempdir()) / "druid_synth_profile"


def _parse_floats(text):
    return [float(x) for x in str(text).replace("，", ",").split(",") if x.strip()]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="两域深度剖面正演：扫描年龄差 Δt（可同时扫核年龄与域界位置），"
                    "画 P(正确检出)。--roc 可改跑合并判据的取舍网格。")
    ap.add_argument("--dt", default="5,10,15,20,25,30,40,60,100",
                    help="两域年龄差网格 (Ma)，逗号分隔")
    ap.add_argument("--dt-rel", default=None,
                    help="相对 Δt 网格（占核年龄的百分比，逗号分隔）。给了它就"
                         "忽略 --dt：每个核年龄各按自己的比例换算 —— 跨年龄比较"
                         "相对跨度时必须用它")
    ap.add_argument("--top", default=f"{TOP_DEFAULT:g}",
                    help="浅处（坑口）那一域的年龄 (Ma)，可给多个（逗号分隔）"
                         f"做年龄扫描（默认 {TOP_DEFAULT:g}）")
    ap.add_argument("--beta", default=f"{BETA_DEFAULT:g}",
                    help="域界所在的归一化深度，可给多个（逗号分隔）做厚薄扫描"
                         f"（默认 {BETA_DEFAULT:g}）")
    ap.add_argument("--seeds", type=int, default=12, help="每个条件跑几个批次 (默认 12)")
    ap.add_argument("--repeat", type=int, default=4,
                    help="每个 Δt 在一个批次里造几个测点 (默认 4)")
    ap.add_argument("--win", type=float, default=4.0, help="滑动窗口宽度 (s，默认 4)")
    ap.add_argument("--step", type=float, default=1.0, help="滑动步长 (s，默认 1)")
    ap.add_argument("--u-cps", type=float, default=1.0e6, help="238U 平台计数率 (cps)")
    ap.add_argument("--th-u", type=float, default=0.5, help="表观 Th/U (默认 0.5)")
    ap.add_argument("--ref-preset", default=DEFAULT_REF_PRESET, help="参考值档")
    ap.add_argument("--out", default=None, help="临时批次写哪儿（跑完即删）")
    ap.add_argument("--fig", default=None, help="图写到哪儿 (PNG)")
    ap.add_argument("--fig-age", default=None,
                    help="年龄扫描图写到哪儿 (PNG)：左＝绝对 Δt，右＝相对 Δt/年龄")
    ap.add_argument("--fig-age-cond", default=None,
                    help="年龄扫描图用哪条条件：名字子串或**序号**（从 1 起；"
                         "默认取测点最多的那条）")
    ap.add_argument("--csv", default=None, help="逐测点判读长表写到哪儿 (CSV)")
    ap.add_argument("--plot-only", default=None, metavar="CSV",
                    help="只读已有的逐测点 CSV 重画，不重算")
    ap.add_argument("--roc", action="store_true",
                    help="改跑合并判据的取舍网格（min_frac × merge_n_sigma），"
                         "其余效应固定为全套；Δt=0 那一列就是假阳性率")
    ap.add_argument("--conditions", default=None,
                    help="只跑这些条件：给名字子串（逗号分隔），或给**序号**"
                         "（从 1 起，如 1,5,7 —— 条件名是中文，序号更稳）")
    ap.add_argument("--pilot", action="store_true",
                    help="快速试跑：1 个种子 + 少量 Δt（用来先量时间）")
    ap.add_argument("--keep", action="store_true", help="保留写下的批次目录")
    ap.add_argument("--verbose", action="store_true", help="转发完整流水线日志")
    args = ap.parse_args(argv)

    ensure_utf8_streams()
    conds = roc_conditions() if args.roc else base_conditions()
    if args.conditions:
        keys = [k.strip() for k in args.conditions.split(",") if k.strip()]
        if keys and all(k.isdigit() for k in keys):
            # 纯数字 ⇒ 按条件表的**序号**（从 1 起）取。条件名是中文，从命令行
            # 传中文在 Windows 上有编码风险，序号是稳的。
            idx = [int(k) - 1 for k in keys]
            conds = [conds[i] for i in idx if 0 <= i < len(conds)]
        else:
            conds = [c for c in conds if any(k in c["name"] for k in keys)]
        if not conds:
            raise SystemExit("--conditions 没匹配到任何条件")

    tops = _parse_floats(args.top)
    betas = _parse_floats(args.beta)
    dt_rel = _parse_floats(args.dt_rel) if args.dt_rel else None

    if args.plot_only:
        det = pd.read_csv(args.plot_only)
        dts = sorted(det["dt"].unique())
    else:
        dts = _parse_floats(args.dt)
        seeds, repeat = args.seeds, args.repeat
        if args.pilot:
            dts, dt_rel, seeds, repeat = [20.0, 40.0], None, 1, 2
            tops, betas = tops[:1], betas[:1]
        n_dt = len(dt_rel) if dt_rel else len(dts)
        out_root = Path(args.out) if args.out else _default_out()
        out_root.mkdir(parents=True, exist_ok=True)
        print(f"两域深度剖面正演：核 {', '.join(f'{t:g}' for t in tops)} Ma，"
              f"边 = 核 − Δt，域界 τ = {', '.join(f'{b:g}' for b in betas)}")
        if dt_rel:
            print(f"Δt 网格（相对）：{', '.join(f'{r:g}%' for r in dt_rel)} × 核年龄")
        else:
            print(f"Δt 网格（绝对）：{', '.join(f'{d:g}' for d in dts)} Ma")
        print(f"条件 {len(conds)} × 种子 {seeds} × 每格 {repeat} 测点 × "
              f"核年龄 {len(tops)} × 域界 {len(betas)} × Δt {n_dt} "
              f"⇒ 约 {len(conds) * seeds * repeat * len(tops) * len(betas) * n_dt} "
              f"个样品测点")
        print(f"窗口 {args.win:g} s / 步长 {args.step:g} s；"
              f"临时目录 {out_root}（{'保留' if args.keep else '跑完即删'}）")
        print()
        for c in conds:
            print(f"  · {c['name']}  (min_frac={c['min_frac']:g}, "
                  f"ns={c.get('merge_n_sigma', MERGE_N_SIGMA_DEFAULT):g}, "
                  f"noise={c['flicker']:.0%}, drift={c['drift']:.0%}, "
                  f"dhf={c['dhf']:.0%}, f206={c['f206']:.0%})")
        print()
        det = run_all(dts, conds, seeds, repeat, tops, betas,
                      args.u_cps, args.th_u, args.ref_preset, args.win,
                      args.step, out_root, args.keep, args.verbose,
                      dt_rel=dt_rel)
        if not args.keep:
            # 各批次的子目录已在 run_all 里删掉，这里只用**非递归**的空目录删除
            # 收掉临时根 —— 绝不递归，免得把用户 `--out` 指向的目录连内容一起清掉。
            _rmdir_if_empty(out_root)
        if args.csv:
            Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
            det.to_csv(args.csv, index=False, encoding="utf-8-sig")
            print(f"逐测点长表：{args.csv}")

    if det is None or det.empty:
        print("没有拿到任何测点判读结果。")
        return 1

    sum_df = summarize(det)
    multi = (sum_df["top"].nunique(dropna=True) > 1
             or sum_df["beta"].nunique(dropna=True) > 1)
    print()
    print("— 汇总（P = 正确检出两个域的比例；跨度保留 = 检出两域的年龄差 / 真值 Δt）—")
    pre_h = f"{'核(Ma)':>7} {'β':>5} " if multi else ""
    print(f"{pre_h}{'条件':<24s} {'Δt(Ma)':>7} {'Δt/核%':>7} {'n':>4} "
          f"{'P(2域)':>7} {'P(年龄)':>8} {'P(严格)':>8} {'漏检':>6} "
          f"{'过分割':>7} {'跨度保留':>8} {'边界|Δτ|':>9}")
    for _, r in sum_df.iterrows():
        pre = f"{r['top']:>7.0f} {r['beta']:>5.2f} " if multi else ""
        print(f"{pre}{r['条件']:<24s} {r['dt']:>7.1f} {r['dt_rel']:>7.2f} "
              f"{int(r['n']):>4} "
              f"{r['P_2域']:>7.2f} {r['P_年龄']:>8.2f} {r['P_严格']:>8.2f} "
              f"{r['P_漏检']:>6.2f} {r['P_过分割']:>7.2f} "
              f"{r['跨度保留']:>8.2f} {r['边界误差中位']:>9.3f}")

    # 逐 (核 × 域界 × 条件) 的"边缘"：P 首次 ≥ 0.5 的位置
    print()
    print("— P(严格) 的检出边缘（首次 ≥ 0.5）—")
    keys = ["top", "beta", "条件"] if multi else ["条件"]
    for k, g in sum_df.groupby(keys, sort=False):
        g = g.sort_values("dt")
        hit = g[g["P_严格"] >= 0.5]
        if len(hit):
            r0 = hit.iloc[0]
            tail = (f"Δt ≈ {r0['dt']:.1f} Ma"
                    f"（Δt/核 = {r0['dt_rel']:.1f}%）")
        else:
            tail = "在网格内未达到 0.5"
        tag = (f"核 {k[0]:g} / β={k[1]:g} / {k[2]}" if multi
               else str(k[0] if isinstance(k, tuple) else k))
        print(f"  {tag:<44s} {tail}")

    single = (len(tops) == 1 and len(betas) == 1)
    fig = args.fig
    if fig is None and args.csv:
        fig = str(Path(args.csv).with_suffix(".png"))
    if fig:
        if not single:
            print()
            print("（--fig 只画单核单域界；多值请用 --fig-age）")
        else:
            Path(fig).parent.mkdir(parents=True, exist_ok=True)
            make_figure(sum_df, fig, tops[0], betas[0], args.win, args.step,
                        min_frac=0.05)
            print()
            print(f"图：{fig}")

    if args.fig_age:
        conds_seen = list(dict.fromkeys(map(str, sum_df["条件"])))
        if args.fig_age_cond:
            key = args.fig_age_cond.strip()
            if key.isdigit():
                i = int(key) - 1
                cond0 = conds_seen[i] if 0 <= i < len(conds_seen) else conds_seen[0]
            else:
                cand = [c for c in conds_seen if key in c]
                cond0 = cand[0] if cand else conds_seen[0]
        else:
            tot = sum_df.groupby("条件")["n"].sum()
            cond0 = str(tot.idxmax())
        Path(args.fig_age).parent.mkdir(parents=True, exist_ok=True)
        make_age_figure(sum_df, args.fig_age, cond0, args.win, args.step)
        print()
        print(f"年龄扫描图：{args.fig_age}（条件「{cond0}」）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
