"""
druid.io.report —— 结果报表写出（Excel）
======================================

为什么用 openpyxl 而不是 pandas 直接 to_excel
---------------------------------------------
pandas.to_excel 只能写数据，写完之后还得再打开文件去调列宽、冻结窗格、
条件着色。这里统一用一个 writer 处理，让每张表都自动具备：
    · 冻结首行（往下翻时表头一直在）
    · 列宽自适应 + 长中文表头不会被挤成一团
    · 表头加粗加底色
这样每次打开结果都是能直接看的样子，省掉一遍手工调整。
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

# ── 样品年龄该取哪一列 ──
# 结果表里有两套年龄：`年龄206_238`（只做外标归一化）与
# `年龄206_238_QC校正`（再乘了监控标样的基体匹配系数）。报告"样品年龄"时应当
# 优先用后者 —— 但**必须先判断该列在不在**，不能写成
# `res.get("年龄206_238_QC校正", unk["年龄206_238"])`：
# `DataFrame.get(key)` 拿到的是**整表**那一列，不是筛选后的子表，
# 于是 43 个标样测点会混进"样品年龄"统计（中位 458.1 被抬到 460.9、
# 5–95% 区间从 420~644 撑成 333~1044 Ma）。这个错在 CLI 与 webui 里
# **各自发生过一次**，所以现在只有这一处定义。
AGE68_COLUMN = "年龄206_238"
AGE68_QC_COLUMN = "年龄206_238_QC校正"


def age68_column(res: pd.DataFrame) -> str:
    """
    报告"样品年龄"时应当取哪一列。

    有二次校正列就用它（那是最终交付口径），否则退回未校正列。
    只返回**列名**，取值一律由调用方从自己的子表里取 —— 两步分开，别再合成一步。
    """
    return AGE68_QC_COLUMN if AGE68_QC_COLUMN in res.columns else AGE68_COLUMN


def _auto_width(series: pd.Series, header: str) -> int:
    """
    估算某一列应有的字符宽度。

    为什么中文要按 2 个字符宽度算
    -----------------------------
    中文字形在等宽/非等宽字体下视觉宽度约为拉丁字符的 2 倍，
    只按 len() 算的话，"年龄206_238_QC校正" 这种长表头会被压成 "###"。
    """
    try:
        body = series.astype(str).str.len().max()
    # 只兜"这一列量不出宽度"这一种情况：传进来的不是 Series（AttributeError）、
    # 或者列里的东西 astype/len 处理不了（TypeError / ValueError）。
    # ⚠ 不要写成裸 except Exception：列宽只是排版小事，一旦把 KeyError、
    # MemoryError 这类真问题也吞掉，坏的是整张表却查不出原因。
    except (AttributeError, TypeError, ValueError):
        body = 10
    body = 0 if (body is None or (isinstance(body, float) and math.isnan(body))) else float(body)
    head = len(str(header))
    # 粗略折算：中文字符出现就记 2 倍
    cn = sum(1 for ch in str(header) if "\u4e00" <= ch <= "\u9fff")
    width = max(body, head + cn) + 2
    return int(min(max(width, 9), 40))     # 夹到 [9, 40]，避免过宽或过窄


def write_excel(out_path, sheets: dict, freeze_header: bool = True) -> Path:
    """
    把若干 DataFrame 写成一个多 Sheet 的 Excel。

    参数
    ----
    out_path : 输出文件路径（.xlsx）
    sheets   : {sheet 名: DataFrame} 有序字典，顺序即 sheet 顺序
    freeze_header : 是否冻结首行

    返回
    ----
    实际写出的 Path（方便链式调用 / 打印）

    细节
    ----
    · 空的 DataFrame 也照写，保留表头，避免调用方拿到"少一个 sheet"的意外；
    · openpyxl 引擎：.xlsx 格式标准、兼容性好，不需要装 Excel。
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            df = df if df is not None else pd.DataFrame()
            df.to_excel(writer, sheet_name=str(name)[:31], index=False)

            ws = writer.sheets[str(name)[:31]]
            if freeze_header:
                # 冻结窗格设为 A2：即第 1 行（表头）保持可见，向下滚动时不动
                ws.freeze_panes = "A2"

            # 逐个在内存里的可视列设置宽度
            for idx, col in enumerate(df.columns, start=1):
                ws.column_dimensions[ws.cell(row=1, column=idx).column_letter].width = \
                    _auto_width(df[col] if len(df) else pd.Series(dtype=float), col)
    return out_path


def export_batch(cfg, result, version: str = "") -> Path:
    """
    把 run_batch 的结果落成标准四表 Excel。

    为什么要单独拎出来
    ------------------
    命令行（cli.reduce_batch）和网页界面（webui）都需要"跑完之后写出同样的表"。
    如果两边各写一遍，早晚会出现"命令行跑的表多一列、网页跑的表少一列"
    这种最让人困惑的不一致。落盘逻辑必须只有一份。

    这里刻意**没有 import workflow**：只要传入的对象带 results/qc/domains/info
    四个属性即可（鸭子类型），因此不存在循环依赖，写层的单元测试也好做。

    返回
    ----
    实际写出的 Excel 路径
    """
    sheets = {
        "结果": result.results.round(4),
        "标样QC": result.qc.round(4),
    }
    if result.domains is not None and not getattr(result.domains, "empty", True):
        sheets["深度剖面域"] = result.domains.round(3)
    # 逐窗口年龄剖面。列名按 ADEPT 的 Format 4 对齐
    # （Analysis / Time / Age68 / Age68_1s），这张表可以直接交给 R 端
    # 做加权平均与 MSWD。
    windows = getattr(result, "windows", None)
    if windows is not None and not getattr(windows, "empty", True):
        sheets["剖面窗口"] = windows.round(4)

    # 把本次运行的关键参数也写进 Excel，便于半年后回溯"当初用的是什么参数"
    info_rows = [(k, str(v)) for k, v in result.info.items() if k != "skipped"]
    skipped = result.info.get("skipped") or []
    for item in skipped:
        # 被跳过的文件要留痕，否则用户会奇怪"为什么少了几个点"
        info_rows.append(("跳过文件", str(item)))
    if version:
        info_rows.insert(0, ("版本", version))
    sheets["运行参数"] = pd.DataFrame(info_rows, columns=["参数", "值"])

    return write_excel(cfg.out_excel, sheets)
