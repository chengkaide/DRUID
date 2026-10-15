"""
tools/export_stages.py —— 把流水线每一步的中间数据落成机器读的文件（**只读**）
==============================================================================

为什么要有这个工具
------------------
DRUID 目前只有两处落盘出口：`io/report.py` 的 Excel 六表、`io/handoff.py` 的
`.handoff.json`。**十步里其余中间量全部只活在内存里** ——
`BatchResult.spots` 里的 `Tra`（逐 cycle 的 raw / net cps、每通道扣掉的空白、
四个时刻）、主标窗口剖面、夹逼关系、`pidx / F68 / F76` 归一化曲线、两套整段
比值里没启用的那一套 —— **进程一退就没了**。

后果是**没有任何一步可以"外接出去"**。最典型的是 Isoclock
（github.com/sndjgm/Isoclock）：它是**从原始质谱文件自己做全套还原**的工具
（导入 → 扣背景 → 剔离群 → 算普通铅 → 分馏校正 → 定年），所以能喂给它的既
不是比值也不是年龄，而是**逐 cycle 的原始信号** —— 而那一层恰好零出口。

这个工具把这条链补上：重跑一次批次前置步骤，把每一步的产物写成 CSV，
并配一份 `manifest.json` 记录「每个文件的字节数 + sha256 + 行数 + 生成参数」，
外接方据此可以声明"我从第 N 步接"，并校验手上这份是不是当前参数下的产物。

⚠ 三条边界（与仓库其它 `tools/` 脚本同一约定）
----------------------------------------------
1. **只读**：不改 `druid/` 一个字节、不碰任何数值输出，产物只写到 `--out` 目录。
   算法**全部**取自流水线自己的公开函数（`load_batch` / `build_bracketing` /
   `compute_bulk_ratios` / `calibrate_primary` / `estimate_external` /
   `window_profile`），这里**一个数都不另算**。
2. ⚠ **它不是 `run_batch` 的替代品**，也不重写 Excel / handoff：
   它自己重跑一遍前置步骤，只为了把那些"不出口的中间量"拿到手。
3. ⚠ `04_bracketing.csv` 的夹逼邻居是**照着 `workflow.build_bracketing` 的
   `bisect_left / bisect_right` 规则重写的**（`build_bracketing` 只吐出
   `brack_for` 闭包、不吐下标）。那条规则一改，这里要跟着改 ——
   `manifest.json` 的 `notes` 里记了这条。

产物（`--out` 目录）
--------------------
    01_blank_and_window.csv     每点一行：四个时刻 + 装载元信息 + 各通道扣掉的空白
    01_traces/<NN>_<样品>.csv   逐 cycle：时间轴 + 各通道 raw / net cps  ← Isoclock 接点
    02_interval.csv             每点一行：`reduce_interval` 的**全部 15 个字段**
                                （含从没进过任何表的 R82 与三个普通铅标志）
    03_primary_profile.csv      主标逐窗口剖面（F(τ) 标定的基础）
    04_bracketing.csv           每个测点用了哪几个主标做夹逼
    05_bulk_ratios.csv          simple 与 ftau 两套整段比值**都落盘**
    06_calibration.csv          pidx / F68 / F76 曲线节点 + 是否被稳健剔除
    09_windows_full.csv         逐窗口**全字段**（含 R76 / s76 / R82 / rho /
                                i204 / sk / 三个标志 + t0 / t1 / t_mid / tau）
    README.txt                  每个文件一行说明
    manifest.json               版本、参数、每个文件的 sha256 与行数

列名一律 **ASCII**（中文只出现在本文件的说明、README 与 manifest 的
`columns` 描述里）—— 与 `io/handoff.py` 同一条约定：下游可能是任何语言写的。

用法
----
    python tools/export_stages.py --dir examples/EX2022A
    python tools/export_stages.py --dir <批次目录> --out <输出目录>
    python tools/export_stages.py --dir <批次目录> --no-traces      # 不写逐 cycle 大表
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from druid import __version__                                # noqa: E402
from druid.console import ensure_utf8_streams                # noqa: E402
from druid.core.constants import ROLE_PRIMARY, ROLE_UNKNOWN  # noqa: E402
from druid.core.geochronology import age68                   # noqa: E402
from druid.core.references import std_alias, std_ref         # noqa: E402
from druid.depth.windows import window_profile               # noqa: E402
from druid.qc import jsonable                                # noqa: E402
from druid.workflow import (BatchConfig, build_bracketing,   # noqa: E402
                            calibrate_primary, compute_bulk_ratios,
                            estimate_external, load_batch)

#: 本工具的产物契约版本。**破坏性改动才升**（改文件名 / 列名 / 语义），
#: 加文件、加列不升 —— 与 `io/handoff.py` 的 SCHEMA 同一规则。
SCHEMA = "druid.stages/1"

#: 浮点写出用几位有效数字。**不四舍五入到 4 位**（那是给人看的 Excel 的做法）：
#: 这些文件是给机器做进一步计算的，截断会引入看不见的偏差。
_FLOAT = "%.12g"

# 文件名里不允许出现的字符（样品名可能有空格、斜杠、中文标点）
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


# ─────────────────────────────────────────────────────────────────────────────
# 〇、写盘小工具
# ─────────────────────────────────────────────────────────────────────────────
def _safe(name: str) -> str:
    """样品名 → 安全的文件名片段（`SRM 612` → `SRM_612`）。"""
    return _UNSAFE.sub("_", str(name)).strip("_") or "unnamed"


def _dump(df: pd.DataFrame, path: Path, files: list) -> Path:
    """
    写一个 CSV 并登记进 manifest 的文件清单。

    · `encoding="utf-8-sig"`：带 BOM 的 UTF-8。列名是 ASCII 时本不需要 BOM，
      但 `README.txt` 与上游 Excel 都是中文环境，带上 BOM 能让中文 Windows 的
      Excel 双击直接开对编码；机器读时 BOM 也会被 pandas / Python 自动吃掉。
    · 先写 `.tmp` 再 `os.replace`：与 `io/handoff.py` 同一做法 ——
      下游可能在**轮询**这些文件，读到写了一半的 CSV 会解析失败。
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False, encoding="utf-8-sig", float_format=_FLOAT)
    data = tmp.read_bytes()
    tmp.replace(path)
    files.append({
        "path": path.name if path.parent == files_root[0] else
                f"{path.parent.name}/{path.name}",
        "rows": int(len(df)),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "columns": [str(c) for c in df.columns],
    })
    return path


# 记录 `--out` 根目录，让 `_dump` 能算出相对路径
files_root: list = [None]


def _split_sk(df: pd.DataFrame) -> pd.DataFrame:
    """
    把 `window_profile` 的 `sk` 元组列拆成三个数值列。

    `reduce_interval` 返回的 `sk` 是 (206/204, 207/204, 208/204) 三元组，
    pandas 把它当 object 列、存进 CSV 就是字符串 `"(0.12, 4.5, 12.3)"` ——
    下游还得自己解元组字符串，那就不算"机器读"了。拆成三个 float 列才算。
    """
    if "sk" not in df.columns:
        return df
    sk = np.array([list(v) if isinstance(v, (tuple, list)) and len(v) == 3
                   else [np.nan] * 3 for v in df["sk"]], float)
    df = df.drop(columns=["sk"])
    df["sk_206_204"], df["sk_207_204"], df["sk_208_204"] = (
        sk[:, 0], sk[:, 1], sk[:, 2])
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 一、逐 cycle 原始信号（① ── ②a）
# ─────────────────────────────────────────────────────────────────────────────
def export_traces(spots, out: Path, files: list) -> pd.DataFrame:
    """
    每个测点一份 CSV：时间轴 + 各通道 raw / net。

    这是**整条链最底层**、也是唯一能交给 Isoclock 那一类工具的东西：
    它要的是仪器原始的时间分辨信号，后面所有 DRUID 算出来的东西它都会自己重算。

    ⚠ 通道数是**解析阶段就定死的**：`io/qtegra.py` 的 `MASS_RE` 只认
    `202|204|206|207|208|232|235|238`，而 Qtegra 导出的每个 CSV 有 29 列
    （28 通道 + 时间轴）⇒ **21 个通道连这里都进不来**（见 `改进清单.md` A-14 与
    `tools/trace_channels.py`）。要拿回它们得先改 `MASS_RE`。
    """
    meta = []
    trace_dir = out / "01_traces"
    for tr, r in spots:
        masses = sorted(set(tr.raw) | set(tr.net))
        cols = {"t": np.asarray(tr.t, float)}
        for m in masses:
            if m in tr.raw:
                cols[f"{m}_raw"] = np.asarray(tr.raw[m], float)
            if m in tr.net:
                cols[f"{m}_net"] = np.asarray(tr.net[m], float)
        df = pd.DataFrame(cols)
        _dump(df, trace_dir / f"{int(tr.idx) + 1:02d}_{_safe(tr.sample)}.csv", files)

        row = {
            "idx": int(tr.idx) + 1,
            "file": Path(tr.path).stem,
            "sample": tr.sample,
            "role": tr.role,
            "n_cycles": int(len(tr.t)),
            "t_first_s": float(tr.t[0]) if len(tr.t) else None,
            "t_last_s": float(tr.t[-1]) if len(tr.t) else None,
            # 四个时刻。⚠ 它们是**两对**，含义不同：
            #   t_ab0/t_ab1 = 自动识别出的剥蚀起止
            #   t0/t1       = 再裁掉 trim 秒之后的**积分窗口**端点（下游要的是这一对）
            "t_ab0_s": float(tr.t_ab0),
            "t_ab1_s": float(tr.t_ab1),
            "t0_s": float(tr.t0),
            "t1_s": float(tr.t1),
            "duration_s": float(tr.t1 - tr.t0),
            "deadtime_ns": float(tr.au_deadtime_ns),
            "note": str(tr.note),
            "f206": float(r["f206"]),
        }
        for m in masses:
            row[f"blank_{m}"] = float(tr.blank.get(m, float("nan")))
        meta.append(row)
    return pd.DataFrame(meta)


# ─────────────────────────────────────────────────────────────────────────────
# 二、逐点整段还原的全部字段（②b）
# ─────────────────────────────────────────────────────────────────────────────
def export_interval(spots) -> pd.DataFrame:
    """
    `reduce_interval` 返回 15 个键，结果表只带出其中一部分。**全部落盘。**

    最要紧的是三个**从没进过任何表**的量：
        R82             208Pb/232Th（算了，但四张表一张都不带）
        i204 / sk       估算的 204Pb 强度、当时用的普通铅比值
    另有三个普通铅标志（`i204_significant` / `common_lead_applied` /
    `sk_converged`）—— 它们在结果表里有，这里再带一遍是为了让这份文件**自足**。
    """
    rows = []
    for tr, r in spots:
        d = {"idx": int(tr.idx) + 1, "file": Path(tr.path).stem,
             "sample": tr.sample, "role": tr.role}
        for k in ("R68", "R76", "R82", "ThU", "f206", "s68", "s76", "rho",
                  "U_cps", "n_cycles", "i204"):
            d[k] = r.get(k)
        sk = r.get("sk")
        d["sk_206_204"], d["sk_207_204"], d["sk_208_204"] = (
            (tuple(sk) + (None, None, None))[:3] if sk is not None else (None, None, None))
        for k in ("i204_significant", "common_lead_applied", "sk_converged"):
            d[k] = r.get(k)
        rows.append(d)
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# 三、主标窗口剖面（③）与夹逼关系（④）
# ─────────────────────────────────────────────────────────────────────────────
def export_primary_profile(spots, cfg: BatchConfig) -> pd.DataFrame:
    """
    主标的逐窗口剖面 —— F(τ) 曲线的标定基础。

    ⚠ 这段是**照着 `workflow.build_bracketing` 的调用重跑的**（同一个
    `window_profile`、同一个 `fix=(f206, sk)`、同一个 `n_sigma`）。那一步算完
    只把剖面留在闭包里，外面拿不到 —— 所以要外接就只能重跑一次。
    两边一旦不一致，`03` 与流水线真实用的 F(τ) 就会对不上。
    """
    rows = []
    for tr, r in spots:
        if tr.role != ROLE_PRIMARY:
            continue
        prof = window_profile(tr, cfg.win, cfg.step,
                              fix=(r["f206"], r["sk"]),
                              n_sigma=cfg.n_sigma_common_pb)
        if prof.empty:
            continue
        prof = prof.copy()
        prof.insert(0, "idx", int(tr.idx) + 1)
        prof.insert(1, "sample", tr.sample)
        rows.append(_split_sk(prof))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def export_bracketing(spots) -> pd.DataFrame:
    """
    每个测点用哪两个主标做夹逼。

    ⚠ 下标是**照 `workflow.build_bracketing` 的 `bisect_left / bisect_right`
    规则重写的** —— `build_bracketing` 只返回 `brack_for` 闭包，不返回主标序号。
    见模块 docstring 的边界第 3 条。
    """
    from bisect import bisect_left, bisect_right
    std_pos = [k for k, (tr, _) in enumerate(spots) if tr.role == ROLE_PRIMARY]
    rows = []
    for k, (tr, _) in enumerate(spots):
        if tr.role != ROLE_UNKNOWN:
            continue
        i, j = bisect_left(std_pos, k), bisect_right(std_pos, k)
        nb = []
        if i > 0:
            nb.append(std_pos[i - 1])
        if j < len(std_pos):
            nb.append(std_pos[j])
        rows.append({
            "idx": int(tr.idx) + 1, "sample": tr.sample,
            "n_bracket": len(nb),
            "bracket_idx": "|".join(str(int(spots[m][0].idx) + 1) for m in nb),
            "bracket_samples": "|".join(spots[m][0].sample for m in nb),
        })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# 四、两套整段比值（⑤）与归一化曲线（⑥）
# ─────────────────────────────────────────────────────────────────────────────
def export_bulk_ratios(spots, variants) -> pd.DataFrame:
    """
    `workflow.compute_bulk_ratios` **两套都算了**，但只有启用的那一套进结果表，
    另一套只剩一列对照年龄。这里两套都落盘 —— 两者的差别本身就是诊断信息。

    ⚠ **`R68_ftau` 会有一大片 NaN，那是如实、不是算错**：`bracket_F` 自
    2026-10-04 起**不再外插**（外插夹成端点值会造出"剖面两端弯曲"的假结构），
    所以测点的 τ 只要越出它那两个夹逼主标剖面的覆盖范围，该点的整条 ftau
    路径就是 NaN。示例批次实测 **50 / 83 行**如此（`R68_simple` 一个都不缺）。
    ⚠ 同一个现象会出现在结果表的对照列 `年龄206_238_ftau` 上，而
    `workflow._alternative_method_column` 的注释**只预想了两种** NaN 情形
    （整批无主标、主标全被剔除），没提这第三种、也是实际发生的这一种 ——
    已记进 `改进清单.md` 的待办（未改动代码）。
    """
    raw68, raw76 = variants["simple"]
    ftau68, _ = variants["ftau"]
    rows = []
    for k, (tr, _) in enumerate(spots):
        d = {"idx": int(tr.idx) + 1, "sample": tr.sample, "role": tr.role,
             "R68_simple": float(raw68[k]), "R68_ftau": float(ftau68[k]),
             "R76": float(raw76[k])}
        # 207/206 两套共用同一个数组（ftau 只对 206/238 做逐窗口校正），
        # 所以只在 206/238 上报差，免得读者以为 207/206 也有两套。
        d["delta_R68_pct"] = ((d["R68_ftau"] - d["R68_simple"])
                              / d["R68_simple"] * 100.0
                              if d["R68_simple"] else None)
        rows.append(d)
    return pd.DataFrame(rows)


def export_calibration(spots, b68, b76, ref68, ref76, pidx, cal) -> pd.DataFrame:
    """
    归一化曲线 F 的**节点**（每个主标一个），以及它是否被稳健统计剔除。

    这张表回答的是"漂移曲线长什么样"以及"哪个主标没参与拟合" ——
    两者在 `info` 里只有一个 RSD 数字，看不出形状。
    """
    pmask = np.asarray(cal["pmask"], bool)
    all_idx = np.flatnonzero(pmask)
    if all_idx.size == 0:
        return pd.DataFrame(columns=["idx", "sample", "kept", "F68", "F76"])
    F68_all, F76_all = ref68 / b68[pmask], ref76 / b76[pmask]
    kept = set(int(x) for x in np.atleast_1d(pidx))
    rows = []
    for j, f68, f76 in zip(all_idx, F68_all, F76_all):
        tr = spots[int(j)][0]
        rows.append({
            "idx": int(tr.idx) + 1, "sample": tr.sample,
            "kept": bool(int(j) in kept),
            "spot_index": int(j),
            "F68": float(f68), "F76": float(f76),
        })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# 五、逐窗口全字段（⑨）
# ─────────────────────────────────────────────────────────────────────────────
def export_windows_full(spots, brack_for, ref68, sd68, cfg: BatchConfig) -> pd.DataFrame:
    """
    逐窗口剖面 + 逐窗口的 F(τ) 与合成 1σ —— **全部字段**。

    与结果 Excel 的「剖面窗口」表**不是替代关系**：那张表按 ADEPT 的 Format 4
    只挑 12 列（`Analysis / Sample / Time / Tau / Age68 / Age68_1s / R68 / s68 /
    ThU / U_cps / f206_pct / n_cycles`），**丢了** `R76 / s76 / R82 / rho /
    i204 / sk / 三个普通铅标志`，也丢了 `t0 / t1`。丢这些是有意的（那张表要喂
    ADEPT），所以**不往那里加列**，而是另出一张全字段的表。

    `Age68_1s_corrected` = 套了 F(τ) 并合成 σext 之后的 1σ，与流水线分域用的
    那个 σ 同一条路径（`analyse_depth_spot` 里 `profile_ages` 的输入）。
    """
    from druid.depth.fractionation import bracket_F
    from druid.depth.fractionation import profile_ages

    rows = []
    for k, (tr, r) in enumerate(spots):
        if tr.role != ROLE_UNKNOWN or not cfg.do_depth:
            continue
        prof = window_profile(tr, cfg.win, cfg.step,
                              fix=(r["f206"], r["sk"]),
                              n_sigma=cfg.n_sigma_common_pb)
        if prof.empty:
            continue
        prof = prof.copy()
        brack = brack_for(k)
        if brack:
            F = bracket_F(prof["tau"].to_numpy(float), brack, "R68", ref68)
            ages, sigs = profile_ages(prof, F, sd68)
            # 自检：与 `analyse_depth_spot` 内部同一条路径（先套 F(τ)、再算年龄）。
            # ⚠ 必须 `equal_nan=True`：cycle 不足的窗口本来就是 NaN，
            #   而 `np.array_equal` 对 NaN 一律判不等（`spot_concordia.py` 踩过同一个坑）。
            if not np.array_equal(
                    age68(prof["R68"].to_numpy(float) * F), ages, equal_nan=True):
                raise AssertionError(
                    "逐窗口年龄与流水线不一致：F(τ) 这条路径被改动过？")
            prof["F68_tau"] = F
            prof["age68_corrected"] = ages
            prof["age68_1s_corrected"] = sigs
        else:
            prof["F68_tau"] = np.nan
            prof["age68_corrected"] = np.nan
            prof["age68_1s_corrected"] = np.nan
        prof.insert(0, "idx", int(tr.idx) + 1)
        rows.append(_split_sk(prof))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# ─────────────────────────────────────────────────────────────────────────────
# 主入口
# ─────────────────────────────────────────────────────────────────────────────
README_LINES = [
    ("manifest.json", "本目录的清单：schema、工具版本、生成参数、每个文件的 sha256 与行数"),
    ("01_blank_and_window.csv", "每点一行：四个时刻（剥蚀起止 / 积分窗口）、装载元信息、各通道扣掉的空白"),
    ("01_traces/", "逐测点逐 cycle：时间轴 + 各通道 raw / net cps。**这是给外部还原软件（如 Isoclock）的接点**"),
    ("02_interval.csv", "每点一行：单点整段还原 reduce_interval 的全部字段（含 R82、i204、sk 与三个普通铅标志）"),
    ("03_primary_profile.csv", "主标逐窗口剖面 —— F(τ) 曲线的标定基础"),
    ("04_bracketing.csv", "每个测点用了哪两个主标做夹逼"),
    ("05_bulk_ratios.csv", "simple 与 ftau 两套整段 206Pb/238U 都落盘，外加 207Pb/206Pb"),
    ("06_calibration.csv", "归一化曲线 F 的逐主标节点，以及它是否被稳健统计剔除"),
    ("09_windows_full.csv", "逐窗口全字段（含 R76 / s76 / R82 / rho / i204 / sk 与 t0 / t1），比 Excel 的「剖面窗口」多一倍列"),
    ("README.txt", "本文件"),
]


def main(argv=None) -> int:
    ensure_utf8_streams()
    ap = argparse.ArgumentParser(
        description="把 DRUID 流水线每一步的中间数据落成机器读的文件（只读，"
                    "不改 druid/ 任何代码、不改任何数值输出）。")
    ap.add_argument("--dir", default=str(ROOT / "examples" / "EX2022A"),
                    help="批次目录（含 *_LIST.csv 与逐点 CSV）")
    ap.add_argument("--out", default=None,
                    help="输出目录；默认「批次目录上一级/<批次名>_stage」")
    ap.add_argument("--win", type=float, default=None, help="窗口宽度 (s)，默认与流水线一致")
    ap.add_argument("--step", type=float, default=None, help="窗口步长 (s)，默认与流水线一致")
    ap.add_argument("--preset", default=None,
                    help="参考值来源（默认与流水线一致：constants.DEFAULT_REF_PRESET）")
    ap.add_argument("--no-traces", action="store_true",
                    help="不写 01_traces/（逐 cycle 大表；整批可能几十 MB）")
    args = ap.parse_args(argv)

    kw = dict(data_dir=str(args.dir), plot=False, verbose=False)
    for key, val in (("win", args.win), ("step", args.step),
                     ("ref_preset", args.preset)):
        if val is not None:
            kw[key] = val
    cfg = BatchConfig(**kw)

    out = (Path(args.out) if args.out
           else Path(args.dir).resolve().parent / f"{Path(args.dir).name}_stage")
    out.mkdir(parents=True, exist_ok=True)
    files_root[0] = out
    files: list = []

    # ── 重跑前置步骤，把不出口的中间量拿到手 ──
    print(f"[1] 装载 {cfg.data_dir}")
    spots, _seq, skipped = load_batch(cfg)
    if not spots:
        raise SystemExit(f"没有任何测点成功还原：{cfg.data_dir}")

    primary = cfg.primary if any(tr.role == ROLE_PRIMARY for tr, _ in spots) else None
    if primary:
        ref68, ref76 = std_ref(std_alias(primary) or primary, cfg.ref_preset)
        brack_for = build_bracketing(spots, cfg)
        variants, (b68, b76) = compute_bulk_ratios(spots, brack_for, ref68, cfg)
        pidx, _F68, _F76, cal = calibrate_primary(spots, b68, b76, ref68, ref76, cfg)
        sd68, _sd76 = estimate_external(spots, pidx, _F68, _F76, b68, b76, cfg)
    else:
        print("    ⚠ 批次里没有主标，只导出 ①② 两层（③④⑤⑥ 需要主标）")
        ref68 = ref76 = float("nan")
        brack_for = lambda k: []                             # noqa: E731
        variants = {"simple": (np.array([]), np.array([])),
                    "ftau": (np.array([]), np.array([]))}
        b68 = b76 = np.array([])
        cal = {"pmask": np.array([False] * len(spots)), "pidx": np.array([])}

    # ── 逐层落盘 ──
    print("[2] 写 01 / 02 ……")
    meta = export_traces(spots, out, files) if not args.no_traces else None
    if meta is not None:
        _dump(meta, out / "01_blank_and_window.csv", files)
    _dump(export_interval(spots), out / "02_interval.csv", files)

    if primary:
        print("[3] 写 03 / 04 / 05 / 06 ……")
        prof = export_primary_profile(spots, cfg)
        if not prof.empty:
            _dump(prof, out / "03_primary_profile.csv", files)
        _dump(export_bracketing(spots), out / "04_bracketing.csv", files)
        _dump(export_bulk_ratios(spots, variants), out / "05_bulk_ratios.csv", files)
        calib = export_calibration(spots, b68, b76, ref68, ref76, pidx, cal)
        if not calib.empty:
            _dump(calib, out / "06_calibration.csv", files)

    print("[4] 写 09（逐窗口全字段）……")
    if cfg.do_depth:
        wf = export_windows_full(spots, brack_for, ref68, sd68, cfg)
        if not wf.empty:
            _dump(wf, out / "09_windows_full.csv", files)

    # ── README + manifest ──
    readme = ["DRUID 分级导出（只读工具 tools/export_stages.py 生成）", ""]
    readme += [f"  {a:<26s} {b}" for a, b in README_LINES]
    readme += [
        "",
        "列名一律 ASCII。浮点用 %.12g 写出（不四舍五入到 4 位）——这些文件是给",
        "机器做进一步计算的，截断会引入看不见的偏差。",
        "",
        "★ 想从第 N 步外接别的软件：看 manifest.json 的 files[].sha256 核对版本，",
        "  再按上面这张表挑文件。最底层、也是唯一能喂给外部还原软件（如 Isoclock）",
        "  的是 01_traces/ —— 那一层之外的东西都已经是 DRUID 加工过的了。",
    ]
    (out / "README.txt").write_text("\n".join(readme) + "\n", encoding="utf-8-sig")

    manifest = {
        "schema": SCHEMA,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "tool": {"name": "druid", "version": __version__},
        "generated_by": "tools/export_stages.py",
        "batch": {"name": Path(str(cfg.data_dir)).name, "data_dir": str(cfg.data_dir)},
        "config": {k: jsonable(v) for k, v in kw.items()},
        "counts": {
            "spots": int(len(spots)),
            "skipped": int(len(skipped)),
            # 逐窗口表的行数（写在哪儿由第 ⑨ 层决定，这里不重算）
            "windows": next((f["rows"] for f in files
                             if f["path"] == "09_windows_full.csv"), 0),
        },
        "steps": {a: b for a, b in README_LINES},
        "files": sorted(files, key=lambda d: d["path"]),
        "notes": [
            "本目录由 tools/export_stages.py 生成，**不改 druid/ 任何代码、"
            "不改任何数值输出**；算法全部取自流水线自己的公开函数。",
            "01_traces/ 的通道数由 io/qtegra.py 的 MASS_RE 决定（目前只认 "
            "202|204|206|207|208|232|235|238），原始 CSV 的另外 21 个通道"
            "连这里也进不来 —— 见 改进清单.md A-14 与 tools/trace_channels.py。",
            "04_bracketing.csv 的夹逼邻居是照 workflow.build_bracketing 的 "
            "bisect 规则重写的（那一步只返回闭包、不返回下标）；那条规则一改，"
            "这里要跟着改。",
            "03_primary_profile.csv 与 09_windows_full.csv 是**重跑** "
            "window_profile 得到的（该步的产物原本只留在闭包与局部变量里）；"
            "调用参数与流水线一致，且 09 带一条自检：重算年龄必须与套 F(τ) "
            "之后的年龄逐位相等，否则直接报错。",
            "05_bulk_ratios.csv 的 R68_ftau 会有一大片 NaN（示例批次实测 50/83）——"
            "bracket_F 自 2026-10-04 起不再外插，测点 τ 越出夹逼主标剖面的覆盖范围"
            "即为 NaN。这是如实，不是算错；R68_simple 一个都不缺。",
        ],
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8")

    print(f"\n输出目录：{out}")
    for f in manifest["files"]:
        print(f"    {f['path']:<34s} {f['rows']:>7d} 行  "
              f"{f['bytes'] / 1024:>9.1f} KiB")
    print(f"    {'manifest.json':<34s} {len(manifest['files']):>7d} 项登记")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
