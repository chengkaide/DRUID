"""
druid.webui.handlers —— Web 界面的业务实现
=========================================

这一层是「网页点一下」和「流水线跑起来」之间的翻译器。

为什么要和 server.py 分开
--------------------------
server.py 只关心 HTTP 协议细节（路由、状态码、静态文件），
handlers.py 只关心业务语义（这个目录是不是有效批次？该配什么参数？）。
分开之后：
    · 将来换成 Flask / FastAPI，server.py 重写即可，业务逻辑一行不用动；
    · 加 REST API / 命令行 / 定时任务，都能直接复用 handlers 里的函数；
    · 单元测试可以不启 HTTP 服务，直接调用这里的函数。
"""
from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from .. import __version__
from ..core.constants import ROLE_LABEL_CN, ROLE_PRIMARY, ROLE_SECONDARY, ROLE_UNKNOWN
from ..io.handoff import export_handoff
from ..io.report import age68_column, export_batch
from ..io.sequence import read_sequence, sequence_summary, spot_csv_path
from ..qc import assess_batch, summary_lines, verdict
from ..reduction.trace import load_spot
from ..workflow import BatchConfig, run_batch

# 一次 BatchConfig 里所有允许从网页端调节的参数白名单。
# 用白名单而不是直接用 request 里的键，是为了防止前端传脏参数直接进核心流程。
_TUNABLE = {
    "primary", "secondary",
    "trim", "blank_dur", "n_sigma_common_pb", "deadtime_ns",
    "win", "step", "bulk",
    "sigma_ext68", "sigma_ext76",
    "do_depth", "plot",
    "merge_n_sigma", "merge_min_frac",
}


# ═════════════════════════════════════════════════════════════════════════════
# 文件系统浏览
# ═════════════════════════════════════════════════════════════════════════════
def list_drives() -> List[str]:
    """
    列出可用的磁盘根目录。

    Windows 下返回所有存在的盘符（C:\\、G:\\ …）；
    类 Unix 直接返回 /。
    """
    if platform.system() == "Windows":
        drives = []
        for code in range(ord("A"), ord("Z") + 1):
            root = f"{chr(code)}:\\"
            if Path(root).exists():
                drives.append(root)
        return drives
    return ["/"]


def browse(path: Optional[str] = None) -> Dict[str, Any]:
    """
    浏览一个目录，返回它的子目录清单与统计信息。

    为什么不用浏览器的 <input type=file webkitdirectory>
    ----------------------------------------------------
    浏览器出于安全限制，不允许网页读取本地任意路径，只能由用户手动选择。
    而这里的典型场景是"数据已经躺在 G:\\某矿区\\某批次 里了"，
    让用户在网页里一层层点进去选，比让他切出去找路径再粘贴回来要顺手。

    返回
    ----
    {
      cwd, parent, drives, dirs:[{name,path,n_csv,has_list}],
      has_list, list_name, n_csv, is_valid_batch
    }
    """
    roots = list_drives()

    # 没给路径就返回根列表。
    # 注意 list_drives() 返回的是字符串，_probe 要的是 Path —— 必须转换，
    # 否则会拿到 "'str' object has no attribute 'name'"。
    if not path:
        return dict(cwd="", parent=None, drives=roots, dirs=[_probe(Path(p)) for p in roots],
                    has_list=False, n_csv=0, is_valid_batch=False)

    p = Path(path)
    if not p.exists():
        return dict(cwd=str(p), parent=str(p.parent) if p.parent != p else None,
                    drives=roots, dirs=[], has_list=False, n_csv=0,
                    is_valid_batch=False, error="路径不存在")

    # 只看子目录 —— 浏览的目的就是找批次文件夹
    subdirs = []
    try:
        for d in sorted(p.iterdir()):
            if d.is_dir() and not d.name.startswith((".", "$")):
                try:
                    subdirs.append(_probe(d))
                except Exception:
                    continue   # 无权限访问的子目录直接忽略，不能让它拖垮整个浏览
    except PermissionError:
        return dict(cwd=str(p), parent=str(p.parent) if p.parent != p else None,
                    drives=roots, dirs=[], has_list=False, n_csv=0,
                    is_valid_batch=False, error="没有访问权限")

    info = _probe(p)
    info.update(drives=roots, dirs=subdirs,
                parent=str(p.parent) if p.parent != p else None)
    return info


def _probe(d) -> Dict[str, Any]:
    """
    探测一个目录"像不像一个 U-Pb 批次目录"。

    判据很朴素：里面有 *_LIST.xls(x) 且有一定数量的 *_N.csv。
    这个探测要快 —— 浏览时会对每个子目录调用，所以只做文件名匹配，不读文件内容。
    """
    # 入参可能是 str（盘符），统一成 Path，避免后续属性访问炸掉
    d = Path(d)
    try:
        entries = list(d.iterdir())
    except Exception:
        return dict(name=d.name, path=str(d), n_csv=0, has_list=False)

    n_csv = sum(1 for e in entries if e.suffix.lower() == ".csv")
    list_name = None
    for e in entries:
        if e.suffix.lower() in (".xls", ".xlsx") and "_LIST" in e.name.upper():
            list_name = e.name
            break
    return dict(name=d.name, path=str(d), n_csv=n_csv,
                has_list=list_name is not None, list_name=list_name)


# ═════════════════════════════════════════════════════════════════════════════
# 批次探测（点选目录后，先把"这里有什么"告诉用户）
# ═════════════════════════════════════════════════════════════════════════════
def inspect_batch(data_dir: str) -> Dict[str, Any]:
    """
    检查一个批次目录，返回给人看的结构摘要 + 给机器用的合法性判据。

    这一步的意义在于**让用户点完目录立刻知道自己点对了没**，
    而不是等到跑了两分钟才在日志里看到 "LIST.xls 不存在"。
    """
    cfg = BatchConfig(data_dir=data_dir)
    res: Dict[str, Any] = dict(
        data_dir=str(cfg.data_dir),
        name=cfg.data_dir.name,
        list_file=str(cfg.list_file),
        list_exists=cfg.list_file.exists(),
        out_excel=str(cfg.out_excel),
        csv_count=0,
        valid=False,
        problems=[],
    )

    # CSV 数量：直接数文件名符合 <批次名>_N.csv 的
    res["csv_count"] = sum(
        1 for p in cfg.data_dir.glob("*.csv")
        if re.search(r"_\d+\.csv$", p.name, re.IGNORECASE))

    if not res["list_exists"]:
        res["problems"].append(f"找不到序列表 {cfg.list_file.name}"
                               f"（期望文件名 <目录名>_LIST.xls/.xlsx）")
    if res["csv_count"] == 0:
        res["problems"].append("目录下没有任何 *_N.csv 数据文件")

    # 序列表存在才去读，读的代价很低（只有几十行）
    if res["list_exists"]:
        try:
            seq = read_sequence(cfg.list_file)
            counts: Dict[str, int] = {}
            for role in seq["role"]:
                counts[role] = counts.get(role, 0) + 1
            res["sequence"] = dict(
                rows=int(len(seq)),
                summary=sequence_summary(seq),
                counts={ROLE_LABEL_CN.get(k, k): v for k, v in counts.items()},
                preview=seq.head(6).astype(str).to_dict("records"),
            )
            # 序列里登记的文件是否在磁盘上都找得到。
            # 必须用 spot_csv_path —— LIST 里的名字不带 .csv 后缀，
            # 直接拼路径会误报"全部文件缺失"。
            missing = [f for f in seq["file"]
                       if not spot_csv_path(cfg.data_dir, f).exists()]
            res["missing_files"] = [str(m) for m in missing[:10]]
            if missing:
                res["problems"].append(
                    f"序列表里登记但有 {len(missing)} 个文件在磁盘上找不到")

            # ── 角色计数 + 剥蚀时长离散度（用于给用户的提示，不阻断处理） ──
            n_primary = int((seq["role"] == ROLE_PRIMARY).sum())
            n_secondary = int((seq["role"] == ROLE_SECONDARY).sum())
            n_unknown = int((seq["role"] == ROLE_UNKNOWN).sum())
            durs = []
            for _, row in seq.iterrows():
                cp = spot_csv_path(cfg.data_dir, row["file"])
                if not cp.exists():
                    continue
                try:
                    tr = load_spot(int(row["order"]) - 1, str(cp),
                                   str(row["sample"]), str(row["role"]), trim=cfg.trim)
                    durs.append(tr.t1 - tr.t0)
                except Exception:
                    continue
        except Exception as e:                       # noqa: BLE001
            res["problems"].append(f"序列表读取失败：{type(e).__name__}: {e}")
            n_primary = n_secondary = n_unknown = 0
            durs = []

    res["valid"] = res["list_exists"] and res["csv_count"] > 0 and not res["problems"]

    # ── 非阻断提示（notices）：让用户在点"开始处理"之前就知道自己数据的特点 ──
    # 这些不是错误，不会挡住运行——即便没有标样，也会进入"未校准诊断模式"
    # 把数据还原出来；只是结果的含义需要按提示理解。
    notices: List[str] = []
    if res.get("sequence"):
        # 1) 没有主标 → 无法做外标归一化，进入未校准模式
        if n_primary == 0:
            notices.append(
                "⚠ 未找到主标样（默认 91500）。将以『未校准 / 原始比值』模式处理："
                "比值取实测值、年龄仅供参考、不可用于定年；深度域判别与二次校正自动跳过。"
                "请加入主标样，或在参数里把『主标名』设成你实际使用的标样名后重跑。")
        # 2) 没有监控标样 → 外部重现性用默认值，QC 与二次校正跳过
        if n_secondary == 0:
            notices.append(
                "未找到监控标样（默认 Ple）。将用默认外部重现性估计不确定度，"
                "监控标样 QC 与二次校正会跳过。")
        # 3) 只有标样、没有未知样品
        if n_unknown == 0:
            notices.append(
                "本批次只有标样（如 91500），没有未知样品；"
                "处理后将只输出标样 QC 表，深度域判别与二次校正自动跳过。")
        # 4) 剥蚀时长差异较大
        if len(durs) >= 2:
            lo, hi = min(durs), max(durs)
            if hi - lo > 8.0:
                notices.append(
                    f"各测点剥蚀时长差异较大（约 {lo:.0f}~{hi:.0f} s）。"
                    "已按每个测点自身的归一化深度 τ 分别处理，深度剖面与对照仍可靠；"
                    "若差异来自激光条件不稳定，建议复核剥蚀参数。")

    res["notices"] = notices
    # 向后兼容：保留单一 notice 字段（旧前端/旧逻辑可能仍读它）
    res["notice"] = notices[0] if notices else ""
    return res


# ═════════════════════════════════════════════════════════════════════════════
# 配置构造
# ═════════════════════════════════════════════════════════════════════════════
def make_config(payload: Dict[str, Any]) -> BatchConfig:
    """
    把前端传来的 JSON 字典翻译成 BatchConfig。

    白名单过滤：只接受 BatchConfig 里真实存在、且我们主动选择放开的字段。
    这样即使前端被改坏、或者有人手工伪造请求，也传不进脏参数把核心流程搞崩。
    """
    data_dir = str(payload["data_dir"])

    kwargs: Dict[str, Any] = {}
    for key in _TUNABLE:
        if key in payload and payload[key] is not None:
            kwargs[key] = payload[key]

    # 输出路径可以直接指定；缺省就用 BatchConfig 的默认（和数据放一起）
    if payload.get("out_excel"):
        kwargs["out_excel"] = str(payload["out_excel"])
    if payload.get("plot_dir"):
        kwargs["plot_dir"] = str(payload["plot_dir"])

    kwargs["verbose"] = True     # Web 端必须开，否则抓不到任何日志
    return BatchConfig(data_dir=data_dir, **kwargs)


# ═════════════════════════════════════════════════════════════════════════════
# 任务体：在工作线程里执行的那一段
# ═════════════════════════════════════════════════════════════════════════════
def run_job(cfg: BatchConfig, task) -> None:
    """
    在工作线程里被调用的完整任务体：跑流水线 → 写 Excel → 登记产物。

    注意这里**不能抛异常出去**（抛了就没人接），异常由 TaskManager 兜住。
    所有 print 都会被 TaskManager 重定向到任务日志里，前端实时可见。
    """
    print("=" * 68)
    print(f"druid v{__version__}   批次 {cfg.data_dir.name}")
    print(f"数据目录：{cfg.data_dir}")
    print(f"整段比值方法：{cfg.bulk}   204 检验门槛：{cfg.n_sigma_common_pb}σ   "
          f"窗宽/步长：{cfg.win}/{cfg.step} s")
    print("=" * 68)

    result = run_batch(cfg)
    out_path = export_batch(cfg, result, version=__version__)

    # 质控：与命令行共用同一套检查项与同一套排版（druid.qc），
    # 免得"网页说能用、命令行说不能用"这种最难查的不一致。
    try:
        checks = assess_batch(result, cfg)
        handoff_path = export_handoff(cfg, result, checks=checks,
                                      version=__version__)
        print("[8] 质控交接：%s" % handoff_path)
        for line in summary_lines(checks, detail_ref=handoff_path.name):
            print(line)
    except Exception as e:                       # noqa: BLE001
        # 交接文件写失败不该把整个任务判死 —— Excel 已经写好了。
        # 但**必须**让用户在日志里看到，否则下游拿不到文件还以为是没跑。
        checks, handoff_path = [], None
        print("（质控交接生成失败，结果表不受影响）%s: %s"
              % (type(e).__name__, e))

    # 登记产物，前端据此渲染下载清单。
    # 注意路径一律转 str ——  dict 里放 Path 对象没法 JSON 序列化。
    outputs = [dict(kind="excel", label="U-Pb 结果总表", path=str(out_path))]
    if handoff_path is not None:
        outputs.append(dict(kind="json", label="质控交接（给下游程序读）",
                            path=str(handoff_path)))

    # 如果生成了图件，把 PDF 汇总册也挂上去（单张 PNG 太多，只给册子）
    if cfg.plot:
        pdf = Path(cfg.plot_dir) / f"{cfg.data_dir.name}_深度剖面_全部.pdf"
        if pdf.exists():
            outputs.append(dict(kind="pdf", label="深度剖面图册", path=str(pdf)))

    task.outputs = outputs

    # 摘要：给前端在完成后直接显示几个关键数字，不用用户自己去翻 Excel
    try:
        res = result.results
        unk = res[res["类型"] == ROLE_LABEL_CN["unknown"]]
        summary: Dict[str, Any] = dict(
            n_points=int(len(res)),
            n_samples=int(len(unk)),
            qc=result.qc.astype(object).where(result.qc.notna(), None).to_dict("records"),
            secondary_correction=result.info.get("secondary_correction"),
        )
        # 质控结论给前端直接用，前端不必解析 JSON 才知道该显示什么颜色
        if checks:
            summary["verdict"] = verdict(checks)
        if len(unk):
            # ⚠ 与 CLI `_summary_lines` 同一段历史 bug 的副本：
            # `res.get("年龄206_238_QC校正", unk["年龄206_238"])` 在 QC 列存在时
            # 返回的是**整表**那一列，标样测点会混进"样品年龄"统计
            # （中位数被抬高、5–95% 区间被撑宽，详见 cli/reduce_batch.py 的注释）。
            # 先定列名，再从 unk 里取列，两步分开，不要再合成一步。
            # 列名怎么定统一在 io.report.age68_column()，别在这儿再写一遍。
            col = age68_column(res)
            a = unk[col]
            q = a.quantile([0.05, 0.5, 0.95])
            summary["age_median"] = round(float(q[0.5]), 1)
            summary["age_p5"] = round(float(q[0.05]), 1)
            summary["age_p95"] = round(float(q[0.95]), 1)
            conc = unk["协和度_pct"]
            summary["concordance_median"] = round(float(conc.median()), 1)
            summary["concordance_in_range"] = round(float(conc.between(90, 110).mean()) * 100)
            struct = unk["深度结构"]
            summary["multi_domain"] = int(struct.str.startswith("多域").sum())
        # domains 可能为空（关掉了深度判别）
        summary["n_domains"] = 0 if result.domains is None else int(len(result.domains))
        # 不分域（整段）口径：全部测点走同一条路。前端据此可以直接显示
        # "48 个测点里 5 个的整段年龄能直接用"，不必让用户去翻表。
        ov = getattr(result, "overall", None)
        if ov is not None and not getattr(ov, "empty", True) and "判定" in ov.columns:
            summary["n_whole_spot"] = int(len(ov))
            summary["whole_spot_compatible"] = int(
                (ov["判定"].astype(str) == "整段常数").sum())
        task.summary = summary
    except Exception as e:                           # noqa: BLE001
        # 摘要算失败不应该让整个任务失败 —— Excel 已经写好了
        print(f"（结果摘要生成失败，可忽略）{type(e).__name__}: {e}")

    print("\n" + "=" * 68)
    print(f"[7] 完成。结果文件：{out_path}")
    print("=" * 68)


# ═════════════════════════════════════════════════════════════════════════════
# 上传：把一些曲别针拽过来的文件拼成一个临时批次
# ═════════════════════════════════════════════════════════════════════════════
def save_uploads(files: List[Dict[str, str]], batch_name: str) -> Dict[str, Any]:
    """
    把前端 base64 传上来的文件写到一个临时批次目录里。

    用处：数据还在 U 盘/别的机器上，不想先手动拷到一个固定目录。
    写一个临时目录后就退化成正常的"批次目录"流程，后面的代码完全不用改。

    files: [{"name": "EX2022A_1.csv", "b64": "..."}]
    """
    import base64

    target = Path(tempfile.gettempdir()) / "upb_uploads" / batch_name
    if target.exists():
        shutil.rmtree(target, ignore_errors=True)
    target.mkdir(parents=True, exist_ok=True)

    saved = []
    for f in files:
        # 防文件名注入：只取最后一段文件名，丢掉任何路径成分
        name = Path(str(f.get("name", ""))).name
        if not name:
            continue
        data = base64.b64decode(f.get("b64", "") or "")
        (target / name).write_bytes(data)
        saved.append(name)

    return dict(dir=str(target), saved=saved,
                inspect=inspect_batch(str(target)) if saved else {})


# ═════════════════════════════════════════════════════════════════════════════
# 打开本机文件夹
# ═════════════════════════════════════════════════════════════════════════════
def reveal(path: str) -> Dict[str, Any]:
    """用系统文件管理器打开某个文件/目录所在的文件夹。"""
    p = Path(path)
    target = p.parent if p.is_file() else p
    if not target.exists():
        target = Path.cwd()

    try:
        if platform.system() == "Windows":
            # /select, 的形式可以在打开资源管理器的同时选中该文件，体验最好
            os.startfile(str(p if p.is_file() else target))  # noqa: S606
        elif platform.system() == "Darwin":
            subprocess.Popen(["open", str(target)])
        else:
            subprocess.Popen(["xdg-open", str(target)])
        return dict(ok=True, opened=str(target))
    except Exception as e:                           # noqa: BLE001
        return dict(ok=False, error=f"{type(e).__name__}: {e}")
