"""
测试用的假数据工厂 —— 名字以 `_` 开头，所以不会被 `run_all.py` 收集成测试文件。

为什么要有它
------------
质控层（`druid.qc`）与交接契约（`druid.io.handoff`）的测试必须能**精确构造**
几种异常情形：没有主标、没有监控标样、σext 触到上下限、有老核测点会被 ADEPT 丢掉、
样品名写成双横线……这些在示例批次里未必同时出现，而且真跑一遍要几秒。

所以这里手搓一个**鸭子类型的 BatchResult**（只需要 results / qc / info / windows
四个属性）。注意：它模拟的是**形状**，不是数值 —— 数值正确性由
`tests/check_example_batch.py` 用真实批次负责。两者分工别搞混。

⚠ 造这个对象时不要顺手把 draid 的算法重新实现一遍。一旦这里开始"照着印象算年龄"，
测试就会变成"验证我的印象"，而不是"验证代码"。
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


class FakeResult:
    """够 `qc.assess_batch` / `handoff.build_payload` 用的最小替身。"""

    def __init__(self, results, qc, info, windows):
        self.results = results
        self.qc = qc
        self.info = info
        self.windows = windows


def make_windows(plan: dict) -> pd.DataFrame:
    """
    plan : {Analysis 名: [该测点逐窗口的 Age68, ...]}

    刻意用 `Analysis` 的值直接把"属于哪个测点"编码进去，与 workflow 里
    `f"{序号:02d} {样品}"` 的写法一致。
    """
    rows = []
    for analysis, ages in plan.items():
        sample = str(analysis).split(" ", 1)[1] if " " in str(analysis) else str(analysis)
        for i, a in enumerate(ages):
            rows.append(dict(Analysis=str(analysis), Sample=sample, Time=float(30 + i),
                             Tau=float(0.05 + 0.02 * i), Age68=float(a),
                             Age68_1s=5.0, R68=0.06, s68=0.0008,
                             ThU=0.5, U_cps=1.0e6, f206_pct=0.0, n_cycles=13))
    return pd.DataFrame(rows)


def make_result(unknown=None, standards=None, mode="已校准",
                sd68=0.011, sd76=0.036, corr=None, info_extra=None,
                calibrate_column=True, rejected=None, skipped=None,
                bulk="simple"):
    """
    构造一个假 BatchResult。

    参数
    ----
    unknown : [{sample, age, conc, f206, u, struct, s1, ftau}]，每个 dict 一个样品测点。
              缺省字段有默认值。`struct` 用 "" / "均一" / "多域(2)"。
    standards : [(名字, 角色, 点数, 参考年龄, 实测加权平均, MSWD)]；给 None 用一套
              "一切正常"的默认值。
    mode : 结果表 `校准状态` 列的值（"已校准" / "未校准"）。
    sd68/sd76 : 写进 info 的外部重现性。
    corr : 二次校正 dict；None 表示未施加。
    calibrate_column : 是否生成 `年龄206_238_QC校正` 列。
    rejected : info["rejected_primary"]
    skipped : info["skipped"]
    """
    unk = unknown if unknown is not None else [
        dict(sample="S01", age=450.0), dict(sample="S02", age=455.0),
    ]
    rows, win_plan = [], {}
    for i, u in enumerate(unk, start=1):
        name = u["sample"]
        rows.append(dict(
            序号=i, 文件=f"B_{i}", 样品=name, 类型="样品",
            剥蚀窗口_s="30.0-60.0",
            U238_cps=u.get("u", 1.0e6), Th_U=u.get("thu", 0.5),
            f206_pct=u.get("f206", 0.0),
            Pb206_238U=0.06, s68_pct=0.8, Pb207_206Pb=0.055, s76_pct=3.0,
            Pb207_235U=0.55,
            年龄206_238=u["age"], s68_1sig=u.get("s1", 3.0),
            s68_2sig=2 * u.get("s1", 3.0),
            年龄207_235=u["age"], s75_2sig=6.0,
            年龄207_206=u["age"], s76_2sig=8.0,
            协和度_pct=u.get("conc", 101.0),
            校准状态=mode,
            年龄206_238_ftau=u.get("ftau", u["age"]),
            深度结构=u.get("struct", ""),
        ))
        win_plan[f"{i:02d} {name}"] = u.get("windows", [u["age"]] * 30)

    # 标样行：只需要让 `类型` 计数正确，其余列给占位值
    n_prim = n_sec = 0
    base = len(unk)
    if standards is None:
        standards = [("91500", "主标", 5, 1062.4, 1059.75, 0.51),
                     ("Ple", "监控标样", 3, 337.13, 343.24, 1.00)]
    for j, (sname, role, n, _ref, _meas, _mswd) in enumerate(standards):
        if role == "主标":
            n_prim += n
        elif role == "监控标样":
            n_sec += n
        for k in range(n):
            rows.append(dict(
                序号=base + j * 10 + k + 1, 文件=f"STD_{j}_{k}", 样品=sname, 类型=role,
                剥蚀窗口_s="30.0-60.0", U238_cps=1.0e5, Th_U=0.4, f206_pct=0.0,
                Pb206_238U=0.06, s68_pct=0.8, Pb207_206Pb=0.055, s76_pct=3.0,
                Pb207_235U=0.55, 年龄206_238=_meas, s68_1sig=4.0, s68_2sig=8.0,
                年龄207_235=_meas, s75_2sig=8.0, 年龄207_206=_meas, s76_2sig=10.0,
                协和度_pct=100.0, 校准状态=mode, 年龄206_238_ftau=_meas,
                深度结构=""))
    res = pd.DataFrame(rows)

    if calibrate_column and mode == "已校准":
        res["年龄206_238_QC校正"] = res["年龄206_238"] * 0.99
        res["s68_1sig_QC校正"] = res["s68_1sig"] * 0.99

    qc = pd.DataFrame([
        dict(标样=s, 点数=n, 参考年龄_Ma=ref, 加权平均年龄_Ma=meas,
             s2_Ma=2.0, MSWD=mswd,
             偏差_pct=(meas - ref) / ref * 100)
        for (s, _role, n, ref, meas, mswd) in standards
        if _role in ("主标", "监控标样")
    ])

    info = dict(
        data_dir=str(Path("C:/fake/BATCH")),
        primary="91500", secondary="Ple",
        ref68=0.17928, ref76=0.07494,
        calibrated=(mode == "已校准"),
        bulk=bulk, win=4.0, step=1.0, trim=1.5, deadtime_ns=0.0,
        sd68=sd68, sd76=sd76,
        rejected_primary=list(rejected or []),
        skipped=list(skipped or []),
        secondary_correction=corr,
        warning=None if mode == "已校准" else "未校准诊断模式",
        n_primary=n_prim, n_secondary=n_sec,
    )
    if info_extra:
        info.update(info_extra)

    return FakeResult(res, qc, info, make_windows(win_plan))


def keyed(checks) -> dict:
    """{key: Check} —— 测试里按 key 取某一条，不按文字匹配。"""
    return {c.key: c for c in checks}
