"""
druid —— 激光剥蚀 ICP-MS（LA-ICP-MS）锆石 U-Pb 定年数据还原包
==========================================================

设计目标
--------
把一条完整的数据还原流程拆成 **职责单一、可单独测试、可单独复用** 的若干层：

    druid.core        物理/数学内核：衰变常数、标准物质参考值、年龄方程、
                    普通铅模型、统计量 —— 不依赖任何 I/O
    druid.io          输入/输出：Qtegra CSV 解析、测点序列表、Excel 报表写出
    druid.reduction   单点数据还原：剥蚀区间识别、气体空白、204Hg 干扰校正、
                    普通铅校正、比值与 jackknife 不确定度
    druid.depth       深度剖面/时间分辨：滑动窗口、F(τ) 逐深度分馏校正、
                    年龄域分割（核/边识别）、图件
    druid.workflow    编排层：把上面各层串成一条完整批处理流水线
    druid.cli         命令行入口：参数解析 + 调用编排层

一句话流程
----------
    Qtegra CSV
      → 自动识别剥蚀区间（扣掉未剥蚀的气体空白段与冲洗尾巴）
      → 扣气体空白
      → 204Hg 干扰校正 + 普通铅校正（带显著性检验！）
      → 主标 91500 夹逼归一化（含逐深度的分馏因子 F(τ)）
      → 外部重现性估计 + 年龄与不确定度
      → 逐点深度剖面年龄域判别
      → Excel 报表 + 剖面图件

两个针对本批数据的关键设计（务必保留，删掉结果会失真）
------------------------------------------------------
1. **204 普通铅校正做显著性检验**
   204Pb 通常只有几~几十 cps，而 204Hg 本底可达数百 cps，两者做差之后
   残余常常就是噪声。若 204Pb 不显著却强行扣除，等于把噪声当成普通铅扣掉，
   对计数本来就很少的 207Pb 是灾难性的。

2. **207Pb/206Pb 不做 F(τ) 逐窗口深度校正**
   它是 Pb 同位素比值，几乎不存在随坑深的元素分馏（随坑深漂移的是 Pb/U
   这种"元素对元素"比值）。窗口级 207Pb 计数极少（4 s 窗口仅约 130 个计数），
   套 F(τ) 只会注入噪声。

用法示例
--------
    # 命令行
    python -m druid.cli.reduce_batch --dir "D:/data/EX2022A" --out 结果.xlsx --plot

    # 作为库调用
    from druid.workflow import BatchConfig, run_batch
    result = run_batch(BatchConfig(data_dir="D:/data/EX2022A"))
"""
from __future__ import annotations

__version__ = "2.2.3"

# 版本说明
# --------
# 2.2.3 新增 docs/getting-started.html —— 面向"第一次拿到这个工具的人"的上手指南：
#   装环境四条命令、数据怎么摆、结果先看哪三个数、出错怎么定位、改代码的边界。
#   同时做了一批纯卫生修正（不改任何数值输出）：
#   · 删掉两个全仓无调用的死函数（io.qtegra.channels_of_interest、
#     depth.figures.multi_page_pdf）及其导出，连带删掉只被前者使用的 UA_MASSES；
#   · reduction.ratios 不再自带一份与 core.constants.MASSES_NEEDED 重复的通道元组；
#   · webui 的 HTTP `Server` 头不再硬编码 "druid-webui/2.0"，改读 __version__；
#   · 修正两处会把读者带偏的注释（QC 校正的 1σ/2σ、去掉扩展名为何不能用 rstrip）；
#   · io.report 的裸 except 收窄为 (AttributeError, TypeError, ValueError)。
#   新手指南里的每个数字都取自示例批次 EX2022A 的端到端回归基线，可自行复现。
#
# 2.2.2 新增 docs/from-first-principles.html —— 从物理原理讲到 MSWD 的完整说明。
#   纯文档改动，不动任何数值（端到端回归的基线值不变）。
#
# 2.2.1 修掉一个只在 Windows + 输出重定向时出现的崩溃：
#   本包的进度与摘要输出全是中文。把 stdout 重定向到文件/管道时，
#   Python 用 locale 编码（英文系统上是 cp1252）编码它，第一句
#   `print("...通过")` 就抛 UnicodeEncodeError 把进程带走 ——
#   而且是在跑完批处理、准备打印结果的时刻，最难受的时机。
#   现在 `druid/console.py` 在入口把输出流切成 UTF-8（errors="replace"）。
#   交互式控制台不一定会踩到（Python 对控制台用宽字符 API），
#   所以本地开发几乎遇不到 —— CI 的 windows job 是第一个撞上的。
#
# 2.0.0 起由原先的 _legacy/*.py 单文件脚本重构为分层子包，
# 旧脚本保留在 druid/_legacy/ 下仅供对照，不再参与运行。
#
# 2.1.0 起 `深度剖面域` 表增加 `MSWD` 与 `MSWD_概率` 两列。这两个数一直算得出来
# （`summarize_segments()` 内部就在调 `weighted_mean()`，它本来就返回 MSWD），
# 只是导出时被丢掉了。补上是为了让本工具与 R 端 ADEPT 的坪年龄口径能直接并列
# 对照 —— 两边都是"反比方差加权平均 + 卡方上尾概率"。只增加列，既有列与数值未动。
#
# 2.2.0 起：
#   · 序列表支持 CSV/TSV 文本格式（原来只认 .xls/.xlsx），优先于二进制格式查找。
#     两列几十行的东西用文本存才能 diff、才能在 review 里看懂。
#   · 随仓库附带示例批次 `examples/EX2022A/`（85 个测点，样品名匿名成 S01…S48，
#     标样保留），CI 用它做端到端数值回归 —— 合成数据覆盖不了"数值算得对"。
#   · 修掉 CLI 收尾摘要里的一处统计错误：`DataFrame.get()` 返回的是**整表**的那一列，
#     于是 35 个标样测点混进了"样品年龄"，中位 458.1 被抬成 460.9 Ma、
#     5–95% 区间从 420~644 撑成 333~1044 Ma。
#   · `_read_rows` 对 .csv/.tsv 改为显式指定分隔符，不再交给嗅探 ——
#     样品名里有空格（`SRM 612`）时嗅探可能误判，把样品名劈成两半。
