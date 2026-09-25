"""
druid.io.handoff —— 给下游（自动化质控 / 解释流程、ADEPT）的 JSON 交接文件
==========================================================================

为什么要有这份文件
------------------
Excel 是给人看的：列名是中文（`年龄206_238_QC校正`）、一张表一个 sheet、
缺失值表示方式随引擎变化、空单元格读回来可能是 `None` 也可能是 `nan`。
用它当**程序之间的接口**，每一个坑都得下游自己踩一遍。

这份 JSON 是给程序看的稳定契约。它不放新数据 —— 每个数都能在 Excel 里找到 ——
它的价值在三件事：

1. **口径写明**：`_1s` 是 1σ 还是 2σ、年龄取的是哪一列、剥蚀窗口用的是哪个区间。
   这些以前只写在文档里，下游要靠读文档猜。
2. **前提写明**：外部重现性是实测还是假设、参考标样库自不自洽、
   是否施加了二次校正。定量结论的可信度取决于这些，不能丢。
3. **交接风险写明**：哪些测点交给 ADEPT 会被静默丢掉、为什么。

契约的稳定性
------------
`schema` 字段是版本号（`druid.handoff/<n>`），**只在做破坏性改动时才升**。
增字段不升版本 —— 下游必须容忍未知字段。改字段名、改语义、改类型才升。
`tests/test_handoff.py` 钉住键集合，所以改名会当场红，不会悄悄漏给下游。

两个刻意决定
------------
* **键名一律 ASCII，中文只出现在给人看的字符串里**（`title`/`detail`/`why`）。
  下游可能是任何语言写的，让它们处理中文键名是平白增加一整个编码类的问题。
* **不给每个样品的加权平均年龄。** 一个样品里各测点的散布通常远大于各自的误差
  （实测岩样级池化 MSWD 能到 10³~10⁴），加权平均只能靠剔掉大部分测点来达标，
  得到的数看着精确、实际是被挑选出来的子集。所以这里只给**描述统计**
  （中位 / 分位 / 范围），把"是否均一"交给 `mswd_pooled` 当**诊断量**，
  并显式附上"要压到 MSWD≤2.5 得剔掉多少个测点"作为证据。
  真出岩样年龄，筛选标准必须来自地质证据（CL 图像、Th/U、协和度），
  不能来自 MSWD。
"""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..core.constants import ROLE_LABEL_CN, ROLE_PRIMARY, ROLE_SECONDARY, ROLE_UNKNOWN
from ..core.statistics import weighted_mean
from ..qc import QCThresholds, assess_batch, by_key, jsonable, verdict
from .report import age68_column

#: 契约版本。破坏性改动才升；增字段不升。
SCHEMA = "druid.handoff/1"

#: 输入表从 Excel 里读出来时，ADEPT / 解释流程该按什么口径用。
#: 这几条是**跨工具契约**，不是本工具的实现细节，所以集中在这里而不是散在文档里。
ADEPT_INPUT = {
    "sheet": "剖面窗口",
    "columns": ["Analysis", "Time", "Age68", "Age68_1s"],
    "age68_1s_sigma": 1,          # 是 1σ，不是 2σ。拿错 MSWD 会差 4 倍。
    "call": {
        "lower_ablation_time": 0,       # DRUID 已切好剥蚀窗口，不要让 ADEPT 再猜
        "upper_ablation_time": 1000000,  # 上界放到最后一个窗口之后
        "smooth": "none",                # 已做 F(τ)，再平滑会抹平真实域边界
        "calibration_uncertainty": 0,    # σext 已含在 Age68_1s 里，再乘一次是重复计入
    },
    "note": "四个调用参数都必须显式给 —— ADEPT 刻意不猜数据是否已预处理。",
}


# ═════════════════════════════════════════════════════════════════════════════
# 小工具
# ═════════════════════════════════════════════════════════════════════════════
def _f(x) -> Optional[float]:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _stats(s: pd.Series) -> Dict[str, object]:
    """
    一个 Series 的描述统计。**只有描述统计，没有加权平均** —— 见模块 docstring。

    `n` 是有效值个数（去 NaN），不是 `len(s)`：两者不等时下游必须看得出来。
    """
    s = pd.to_numeric(s, errors="coerce").dropna()
    if not len(s):
        return {"n": 0, "median": None, "q05": None, "q25": None, "q75": None,
                "q95": None, "min": None, "max": None}
    q = s.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    return {
        "n": int(len(s)),
        "median": _f(q[0.5]), "q05": _f(q[0.05]), "q25": _f(q[0.25]),
        "q75": _f(q[0.75]), "q95": _f(q[0.95]),
        "min": _f(s.min()), "max": _f(s.max()),
    }


def _pooled_mswd(s: pd.Series, sig: pd.Series, threshold: float = 2.5,
                 max_drop: int = 200) -> Dict[str, object]:
    """
    池化 MSWD 及其"要达标得剔掉多少个测点"。

    这是**诊断量，不是年龄**。返回值里没有平均年龄 —— 故意不给，
    免得有人顺手把它当岩样年龄贴出去。

    做法：反复剔掉"离当前加权平均最远的那个（按它自己的 σ 归一化）"，
    直到 MSWD ≤ threshold 或点数不足。
    """
    x = pd.to_numeric(s, errors="coerce")
    e = pd.to_numeric(sig, errors="coerce")
    keep = (x.notna() & e.notna() & (e > 0)).to_numpy()
    idx = list(np.nonzero(keep)[0])
    n0 = len(idx)
    out = {"n_spots": n0, "mswd": None, "n_dropped_for_threshold": None,
           "threshold": float(threshold), "prob": None}
    if n0 < 2:
        return out
    try:
        mu, se, mswd, prob = weighted_mean(x.to_numpy()[keep], e.to_numpy()[keep])
    except Exception:                                   # noqa: BLE001 —— 诊断量不值得中断流程
        return out
    out["mswd"], out["prob"] = _f(mswd), _f(prob)

    cur = list(idx)
    dropped = 0
    while len(cur) >= 3 and dropped < max_drop:
        xv = x.to_numpy()[cur]
        ev = e.to_numpy()[cur]
        try:
            mu, se, mswd, _ = weighted_mean(xv, ev)
        except Exception:                               # noqa: BLE001
            break
        if _f(mswd) is None or mswd <= threshold:
            break
        z = np.abs((xv - mu) / ev)
        cur.pop(int(np.argmax(z)))
        dropped += 1
    out["n_dropped_for_threshold"] = int(dropped)
    return out


def _counts(series: pd.Series) -> Dict[str, int]:
    """{取值: 个数}，按个数降序、同数按名字升序 —— 固定顺序便于 diff。"""
    if series is None:
        return {}
    vc = series.astype(str).value_counts()
    return {str(k): int(v) for k, v in
            sorted(vc.items(), key=lambda kv: (-kv[1], str(kv[0])))}


def _config_dict(cfg) -> Dict[str, object]:
    """把 BatchConfig 摊平成 dict。私有字段（下划线开头）不外泄。"""
    if cfg is None:
        return {}
    out = {}
    for name in getattr(cfg, "__dataclass_fields__", {}):
        if name.startswith("_"):
            continue
        val = getattr(cfg, name, None)
        if hasattr(val, "__dataclass_fields__"):         # 嵌套的阈值对象
            val = {k: jsonable(getattr(val, k))
                   for k in val.__dataclass_fields__}
        elif isinstance(val, Path):
            val = str(val)
        out[name] = jsonable(val)
    return out


# ═════════════════════════════════════════════════════════════════════════════
# 各个块
# ═════════════════════════════════════════════════════════════════════════════
def _block_batch(cfg, result) -> Dict[str, object]:
    data_dir = Path(str(result.info.get("data_dir", ".")))
    return {
        "name": data_dir.name or "",
        "data_dir": str(data_dir),
        "list_file": str(getattr(cfg, "list_file", "") or ""),
    }


def _block_counts(result) -> Dict[str, object]:
    res = result.results
    roles = res["类型"] if "类型" in res.columns else pd.Series(dtype=object)
    struct = res["深度结构"] if "深度结构" in res.columns else pd.Series(dtype=object)
    return {
        "spots": int(len(res)),
        "by_role": _counts(roles),
        "by_depth_structure": {k or "(空)": v for k, v in _counts(struct).items()},
        "windows": int(len(result.windows)) if getattr(result, "windows", None) is not None else 0,
        "skipped_files": [list(map(str, s)) if isinstance(s, (list, tuple)) else str(s)
                          for s in (result.info.get("skipped") or [])],
    }


def _block_calibration(cfg, result, checks) -> Dict[str, object]:
    info = result.info
    res = result.results
    roles = res["类型"] if "类型" in res.columns else pd.Series(dtype=object)
    qc = getattr(result, "qc", None)

    def _std_row(name):
        if qc is None or getattr(qc, "empty", True) or "标样" not in qc.columns:
            return None
        g = qc[qc["标样"] == name]
        return g.iloc[0] if len(g) else None

    def _std(name, role_label):
        q = _std_row(name)
        d = {"name": str(name), "role": role_label,
             "n_spots": int((roles == ROLE_LABEL_CN[role_label]).sum())}
        if q is not None:
            d.update(
                n_used=int(q.get("点数", 0) or 0),
                reference_age_Ma=_f(q.get("参考年龄_Ma")),
                weighted_mean_age_Ma=_f(q.get("加权平均年龄_Ma")),
                s2_Ma=_f(q.get("s2_Ma")),
                mswd=_f(q.get("MSWD")),
                bias_pct=_f(q.get("偏差_pct")),
            )
        return d

    # σext 的来源 —— 下游判断"误差棒能不能用来论证一致性"全靠这一条
    k_src = by_key(checks).get("uncertainty.sigma_ext_source")
    return {
        "mode": str(res["校准状态"].iloc[0]) if "校准状态" in res.columns and len(res) else "未知",
        "primary": _std(info.get("primary", cfg.primary), ROLE_PRIMARY),
        "secondary": _std(info.get("secondary", cfg.secondary), ROLE_SECONDARY),
        "reference_ratios": {
            "R68": _f(info.get("ref68")), "R76": _f(info.get("ref76")),
        },
        "primary_rejected_indices": [int(x) for x in (info.get("rejected_primary") or [])],
        "secondary_correction": jsonable(info.get("secondary_correction")),
        "uncertainty": {
            "sigma_ext68": _f(info.get("sd68")),
            "sigma_ext76": _f(info.get("sd76")),
            "source": (k_src.data.get("source") if k_src else None),
            "n_secondary_spots": int(
                (roles == ROLE_LABEL_CN[ROLE_SECONDARY]).sum()) if len(roles) else 0,
        },
        "note": ("`source=measured` 表示由本批监控标样的实测散度估出；"
                 "`assumed` 表示没有监控标样、退回写死的经验值（206/238 0.70%、"
                 "207/206 0.25%）—— 后者不能用误差棒论证年龄一致性。"),
    }


def _block_standards(result) -> List[Dict[str, object]]:
    qc = getattr(result, "qc", None)
    if qc is None or getattr(qc, "empty", True):
        return []
    out = []
    for _, q in qc.iterrows():
        out.append({
            "name": str(q.get("标样", "")),
            "n": int(q.get("点数", 0) or 0),
            "reference_age_Ma": _f(q.get("参考年龄_Ma")),
            "weighted_mean_age_Ma": _f(q.get("加权平均年龄_Ma")),
            "s2_Ma": _f(q.get("s2_Ma")),
            "mswd": _f(q.get("MSWD")),
            "bias_pct": _f(q.get("偏差_pct")),
        })
    return out


def _block_samples(result, th: QCThresholds) -> List[Dict[str, object]]:
    res = result.results
    if "类型" not in res.columns:
        return []
    unk = res[res["类型"] == ROLE_LABEL_CN[ROLE_UNKNOWN]]
    if getattr(unk, "empty", True):
        return []
    col = age68_column(res)
    out = []
    for name, g in unk.groupby("样品", sort=True):
        d: Dict[str, object] = {"name": str(name), "n_spots": int(len(g))}
        d["age68"] = _stats(g[col]) if col in g.columns else {}
        d["age68_pooled_diagnostic"] = (
            _pooled_mswd(g[col], g["s68_1sig"], th.std_mswd_warn)
            if (col in g.columns and "s68_1sig" in g.columns) else {})
        if "协和度_pct" in g.columns:
            c = pd.to_numeric(g["协和度_pct"], errors="coerce").dropna()
            d["concordance"] = {
                "n": int(len(c)),
                "median_pct": _f(c.median()) if len(c) else None,
                "frac_in_90_110": (float(c.between(90, 110).mean()) if len(c) else None),
            }
        if "f206_pct" in g.columns:
            d["f206_pct_median"] = _f(pd.to_numeric(g["f206_pct"], errors="coerce").median())
        if "U238_cps" in g.columns:
            d["u238_cps_median"] = _f(pd.to_numeric(g["U238_cps"], errors="coerce").median())
        if "Th_U" in g.columns:
            d["th_u_median"] = _f(pd.to_numeric(g["Th_U"], errors="coerce").median())
        if "深度结构" in g.columns:
            strs = g["深度结构"].astype(str)
            d["depth_structures"] = _counts(strs)
            d["n_multi_domain"] = int(strs.str.startswith("多域").sum())
        out.append(d)
    return out


def _block_handoff(checks) -> Dict[str, object]:
    """
    交接块。**信息全部取自质控检查项**，不在这里重算 ——
    否则"检查项说会丢 29 个、交接块说丢 26 个"这种不一致迟早出现。
    """
    byk = by_key(checks)
    k_drop = byk.get("handoff.adept_dropout")
    k_win = byk.get("handoff.windows")
    d = k_drop.data if k_drop else {}
    w = k_win.data if k_win else {}
    return {
        "adept": {
            **ADEPT_INPUT,
            "n_spots": w.get("n_spots"),
            "n_windows": w.get("n_windows"),
            "will_drop_n": d.get("n_dropped", 0),
            "will_drop_frac": d.get("frac"),
            "will_drop_spots": jsonable(d.get("dropped") or []),
        }
    }


# ═════════════════════════════════════════════════════════════════════════════
# 主入口
# ═════════════════════════════════════════════════════════════════════════════
def default_path(cfg) -> Path:
    """默认与结果 Excel 并排同名：`<...>_U-Pb结果.xlsx` → `<...>_U-Pb结果.handoff.json`。"""
    out = Path(str(cfg.out_excel))
    stem = out.with_suffix("") if out.suffix else out
    return stem.parent / (stem.name + ".handoff.json")


def build_payload(cfg, result, checks: Optional[Sequence] = None,
                  version: str = "") -> Dict[str, object]:
    """
    组装交接字典。**不写盘**，方便测试与在内存里检查。

    参数
    ----
    cfg     : `BatchConfig`（可为 None，此时配置块为空 —— 但标样名会从 info 里补）
    result  : `BatchResult`（鸭子类型：results / qc / info / windows）
    checks  : 质控检查项；不给就现算（`druid.qc.assess_batch`）
    version : 工具版本号，写进 `tool.version`
    """
    if checks is None:
        checks = assess_batch(result, cfg)
    checks = list(checks)

    th = getattr(cfg, "qc_thresholds", None) or QCThresholds()

    # cfg 只给了 result 没给时，用一个最小替身，免得下面到处判 None
    if cfg is None:
        class _Cfg:
            primary = str(result.info.get("primary", "91500"))
            secondary = str(result.info.get("secondary", "Ple"))
            out_excel = Path(str(result.info.get("data_dir", "."))) / "out.xlsx"
            list_file = Path("")
        cfg = _Cfg()

    payload = {
        "schema": SCHEMA,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "tool": {"name": "druid", "version": version},
        "batch": _block_batch(cfg, result),
        # 结论块由 qc.verdict() 独家生产 —— 网页界面显示的也是它，
        # 两边不会出现"JSON 说 fail、页面说能用"。
        "verdict": verdict(checks),
        "config": _config_dict(cfg),
        "counts": _block_counts(result),
        # 口径声明。下游不必读文档就能知道这些数是"什么意思"。
        "conventions": {
            "age_column_used": age68_column(result.results),
            "age_unit": "Ma",
            "age68_1s_sigma": 1,
            "age68_2sig_is_double": True,
            "wmean_definition": "反比方差加权平均（与 ADEPT 同口径），不是算术平均",
            "mswd_probability": "chi2 upper tail，等价于 R 的 pchisq(x, df, lower.tail=FALSE)",
            "ablation_window": "已按 trim 裁掉两端瞬态，与 `剖面窗口` 表一致",
        },
        "calibration": _block_calibration(cfg, result, checks),
        "standards": _block_standards(result),
        "samples": _block_samples(result, th),
        "checks": [c.as_dict() for c in checks],
        "handoff": _block_handoff(checks),
        # 这个字段是给"拿到 JSON 就开始算"的流程看的：说清楚哪些量**故意没给**。
        "intentionally_omitted": {
            "per_sample_weighted_mean_age":
                "不给。样品内测点散布通常远大于各自误差（实测池化 MSWD 可达 10³~10⁴），"
                "加权平均只能靠剔除大部分测点达标，得到的数是选出来的子集。"
                "请用 age68 的描述统计 + age68_pooled_diagnostic 判断是否均一，"
                "筛选标准由地质证据（CL、Th/U、协和度）决定。",
        },
    }
    return payload


def export_handoff(cfg, result, checks: Optional[Sequence] = None,
                   version: str = "", path=None) -> Path:
    """
    写交接 JSON，返回实际路径。

    写盘细节
    --------
    · `ensure_ascii=False` + UTF-8：中文原样可读，出问题能直接用编辑器看。
    · `allow_nan=False`：**NaN 不是合法 JSON**。Python 自己能把 `NaN` 读回来，
      别的语言的解析器一律报错 —— 所以必须让它在写的时候就把问题暴露出来。
      `qc.jsonable` 已经把 NaN/Inf 转成 None。
    · 先写 `.tmp` 再 `os.replace`：下游可能在**轮询**这个文件，
      读到写了一半的 JSON 会解析失败。原子替换避免这一整类竞态。
    """
    payload = build_payload(cfg, result, checks=checks, version=version)
    if path is None:
        path = default_path(cfg)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, ensure_ascii=False, indent=2,
                      allow_nan=False, sort_keys=False)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
    return path
