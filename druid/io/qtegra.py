"""
druid.io.qtegra —— Thermo Fisher Qtegra 时间分辨数据解析器
========================================================

实测文件结构（本项目 示例批次 EX2022A，iCAP RQ，Qtegra 2.8）
--------------------------------------------------------------------
    第 1 行      "SRM 612:03/01/2022 06:36:35 AM;"   ← 样品名:日期 时间;
    第 2~13 行   仪器元数据（Software / Configuration / STD / RF / Ion Optics …
                 Detector / Cooling / Power Supply / Gas Supply / Pulse Counting）
    第 14 行     表头：Time,29Si,49Ti,…,202Hg,204Pb,206Pb,207Pb,208Pb,232Th,238U,
    第 15 行     单位行：,dwell time=0.01;xcal factor=67325.46529,…   ← **必须跳过**
    第 16 行起   真正的数值数据

为什么不能直接用 pandas.read_csv
--------------------------------
pandas 默认把第 1 行当表头，而这里前 13 行是元数据、
且中间还夹着一行非数值的"单位行"。
最稳的做法是：手工找到 "Time" 表头所在行，然后**自己解析**，
遇到解析失败的行（单位行）自然跳过 —— 比写一堆 skiprows 参数可靠得多。

返回值单位
----------
所有通道都是 **计数率 cps（counts per second）**，不是原始计数。
Qtegra 已经按 dwell time 换算好了。
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

# 从表头字符串里抽质量数的正则。
# 例："204Pb" → 204 ；"238U" → 238 ；"49Ti" → 49
# ⚠ 只关心 U-Pb 定年要用到的这几个质量数（202/204/206/207/208/232/235/238），
#   微量元素（29Si、178Hf…）不处理 —— 这个分组就是"属于 U-Pb 体系"的判据。
MASS_RE = re.compile(r"(202|204|206|207|208|232|235|238)")


def read_qtegra(path) -> Tuple[str, np.ndarray, Dict[int, np.ndarray], list]:
    """
    解析单个 Qtegra CSV 文件。

    参数
    ----
    path : str 或 Path，CSV 文件路径

    返回
    ----
    (title, t, data, columns)
        title   : 第 1 行冒号前截取出的样品名（如 "SRM 612"）
        t       : 时间轴数组，单位秒 (s)
        data    : {质量数(int): cps 数组}
        columns : 原始列名列表（便于回溯，比如想知道第几个元素是什么）

    异常
    ----
    ValueError : 空文件 / 找不到 Time 表头 / 找不到任何有效数据行
    """
    path = Path(path)
    # errors="replace"：万一元数据里混了非 UTF-8 字节也不至于整个崩掉
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"空文件: {path}")

    # ── 1. 样品名 ──
    # 首行形如 "SRM 612:03/01/2022 06:36:35 AM;"
    # split(":")[0] 取冒号左边，再 strip 去掉空白
    title = lines[0].split(":")[0].strip()

    # ── 2. 定位 Time 表头所在行 ──
    # 不用固定行号：不同批次元数据行数可能不同，用内容匹配更稳
    header_idx = None
    for i, line in enumerate(lines):
        if line.lstrip().lower().startswith("time"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"未找到 Time 表头: {path}")

    # 列名：去空白、去末尾多余逗号再切分
    columns = [c.strip() for c in lines[header_idx].strip().rstrip(",").split(",")]

    # ── 3. 逐行读取数值 ──
    rows = []
    for line in lines[header_idx + 1:]:
        line = line.strip().rstrip(",")
        if not line:                       # 空行直接跳过
            continue
        try:
            # 单位行的第一个字段为空字符串 → float("") 抛 ValueError → 被跳过
            # 这就是"无需显式跳行"的技巧
            rows.append([float(x) for x in line.split(",")])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"无有效数据行: {path}")

    # ── 4. 对齐列数与行数 ──
    # 有些行末尾可能多出字段，取各行列数的最小值，保证能构成规整矩阵
    ncols = min(len(columns), min(len(r) for r in rows))
    arr = np.asarray([r[:ncols] for r in rows], dtype=float)
    columns = columns[:ncols]

    # ── 5. 拆分时间轴与各个质量通道 ──
    t = arr[:, 0]                          # 第 0 列永远是 Time
    data: Dict[int, np.ndarray] = {}
    for j, col in enumerate(columns):
        if j == 0:                         # 跳过时间列自身
            continue
        m = MASS_RE.search(col)
        if m:
            mass = int(m.group(1))
            # setdefault：若同名质量数出现了多次（罕见），保留第一次读到的
            data.setdefault(mass, arr[:, j])

    return title, t, data, columns
