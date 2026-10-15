"""
tools/export_isoclock.py —— 把 DRUID 选好的积分区间写成 Isoclock 的 `Time_setting.npy`（**只读**）
================================================================================================

要解决的问题：**选区间不要做两遍**
----------------------------------
两个工具在这件事上正好相反：

* **DRUID 自动选**：`reduction/ablation.find_ablation`（两级阈值 + 最长连续段）定位剥蚀区间，
  再按 `trim` 裁出积分窗口，四个时刻 `t_ab0 / t_ab1 / t0 / t1` 在 `Tra` 上是现成的；
* **Isoclock 手选**：它**不做信号区间自动检测**。区间由使用者在界面上逐个样品点，
  存进 `Time_setting.npy`，`dataprocess()` 再按 `[b0, b1, s0, s1]` 切出背景段与样品段
  （`Isoclock2.0.py` 行 1346-1353、1370-1386）。

于是同一批 83 个点要在界面上手点 83 次。更要紧的是：**两边区间不一致时，
"算法差异"和"区间差异"就分不开了** —— 而逐个比对两套还原算法正是接 Isoclock 的目的。

这个工具是那座桥：把 DRUID 选好的区间**原样**翻写成 Isoclock 读得懂的文件，
同时把**这套区间本身**写成一份跟软件无关、带完整标注的机读产物 —— 见下一节。

**选一次区间，多处用**（这份产物的第二重用途）
------------------------------------------------
`Time_setting.npy` 是**只给 Isoclock 看的**（pickle、键是文件名、值只能放五个数）。
要"选一次区间做多种数据运算"，光有那个文件不够：别的软件看不懂 pickle，也看不出
"这个数到底是哪一头、闭区间还是半开、下标从 0 还是从 1 数"。

所以本工具另外写两份**与软件无关**的表，每个时刻都带名字、单位与口径：

    interval_selection.csv   一个测点一行：六个时刻（秒）+ 四个 cycle 下标 + 四个时间戳 + 诊断
    interval_cycles.csv      一个 cycle 一行：它落在哪一段（seg）+ 它的归一化深度 τ

有了这两份表，任何接得上的人（Python、R、Excel、别的定年软件）都能**按秒或按下标
复现同一段**，并且能立刻看出自己那一刀切在哪儿。`interval_cycles.csv` 的 `seg`
一列是把"时间点位置"标到**逐 cycle 分辨率**上的那一列 —— 想按别的分段方式重算
（比如换个 trim、多切几段做分域），直接按 `t_s` 重选即可，不需要再回 DRUID 跑一遍。

**时间点命名与口径（读上面两份表之前先看这一段）**
---------------------------------------------------
| 列名 | 对应 DRUID 内部量 | 单位 | 含义 |
|---|---|---|---|
| `ablation_start_s` | `Tra.t_ab0` | s | 剥蚀起点（`find_ablation` 自动识别） |
| `ablation_end_s` | `Tra.t_ab1` | s | 剥蚀终点 |
| `blank_start_s` | `ablation.blank_mean` 的取窗下界 | s | 气体空白取样窗起点 `= max(0, t_ab0 − blank_dur − 0.5)` |
| `blank_end_s` | 同上取窗上界 | s | 气体空白取样窗终点 `= t_ab0 − 0.5` |
| `sample_start_s` | `Tra.t0` | s | 积分窗口起点 `= t_ab0 + trim`（窗口过短时不裁边） |
| `sample_end_s` | `Tra.t1` | s | 积分窗口终点 |

三条口径（写死在产物里，也在 `interval_selection.json` 的 `conventions` 字段重复一遍）：

1. **秒的原点 = 该测点 CSV 时间轴自己的零点**，不是序列的绝对时刻。跨点比较请用
   下标或 τ，不要直接比秒。
2. **所有区间都是闭区间**（`window_mask` 的口径）`[a, b]` 两端都取等号；
   Isoclock 那边是**左闭右开**，所以写进 `Time_setting.npy` 的右端点要落在"末下标 + 1"那一格
   —— 这就是 `iso_sample_end_s` 与 `sample_end_s` 通常不是同一个数的原因。
3. **cycle 下标从 0 开始**（与 `np.flatnonzero` 一致）。τ 的定义与 `depth/windows.py`
   逐字相同：`τ = (t − sample_start_s) / (sample_end_s − sample_start_s)`，
   0 = 坑口、1 = 坑底；**逐 cycle 的 τ 用该 cycle 的时间戳，DRUID 窗口的 τ 用窗口中点**，
   两者定义相同、取值点不同。

Isoclock 那一侧的合同（逐条从 `Isoclock2.0.py` 抄出来，不是猜的）
----------------------------------------------------------------
| 事项 | 出处与事实 |
|---|---|
| 数据格式 | `loaddata()` 的 `ele==0` 支路：`np.loadtxt(..., skiprows=13)`，第 1 行是被跳过的标题、行 14 是表头、行 15 是单位行、行 16 起是数据 ⇒ **与 `io/qtegra.py` 读的 Qtegra 2.8 导出完全同一套偏移**，同一批原始 CSV 双方都能直接读 |
| 键 | **数据文件名（含 `.csv`）**，因为它是 `for name in files` 拿到的 `name` |
| 值 | `num[name] = [Numbers, b0, b1, s0, s1]`（`instructure0` 行 2421 / `instructure1` 行 2727）。`Numbers = len(x)`；其余四个是**秒**：`b0/b1` 背景段起止、`s0/s1` 样品段起止 |
| 秒→下标 | `start_b = int(b0 / Timeinternal)`、`start_s = int(s0 / Timeinternal)` …（`dataprocess` 行 1350-1353）。**截断**，不是四舍五入 |
| 周期 | `Timeinternal = (x[Numbers-1] - x[0]) / Numbers`（行 2358 / 2662）—— 用 `N` 除而不是 `N-1`，且是**全局变量：最后一个载入的文件赢**。本工具按**逐文件**算，见"残余差异" |
| 切数据 | `y[start_b:end_b]` 与 `y[start_s:end_s]`（行 1370-1386）—— **左闭右开**切片，之后各自取 `.mean()` |
| 合法性判据 | `check()`（行 1325-1337）：`s0 - b1 < 0` 或 `s1 - b1 < 0` ⇒ 该点判为错；`s0 + 1 > s1 - 1` ⇒ 区间被就地压塌成零长 |
| 何时写这份文件 | `instructure0/1` 跑完会 `np.save(outputpath + '//' + 'Time_setting', num)`（行 2464）⇒ **加载步会覆盖它**，所以本工具的产物要在"加载"步**之后**放进去 |
| 谁读它 | `dataprocess()` 及之后每一步都从 `Time_setting.npy` 读（`CWD`，即 Isoclock 的 `outputpath`） |

半开切片这一条决定了**秒数该怎么写**：DRUID 的窗口是**闭区间** `[t0, t1]`（`window_mask`
两边都取等号），Isoclock 的切片是**左闭右开**。所以要让两边选中同一批 cycle，写进去的
秒数必须满足

    int(s0_iso / T) = i0        且        int(s1_iso / T) = i1 + 1

（`i0/i1` 是 DRUID 实际选中的首末 cycle 下标）。本工具默认就**反解**这两个式子来生成秒数
（`--no-snap` 可关掉，改成原样写 DRUID 的秒数）。写的是**区间的中点**而不是端点：
`Timeinternal` 在 Isoclock 里是全局量、取值随"最后一个文件"浮动，中点能容忍半个下标的抖动。

**反解不是洁癖，实测有差**：不反解时，示例批次 **83/83 个点的背景段整段早一个 cycle**
（样品段恰好没错，纯属巧合）。原因是 `Timeinternal = span/N` 而不是 `span/(N-1)`，
比真实 cycle 周期小 **0.42%**，而 `int(秒/T)` 的偏差随下标**线性增长** —— 到第 40 个
cycle 就够翻一格。`--no-snap` 的自检会当场拦住这种产物（见 `selfcheck.slice_match`）；
那一次负对照实测：样品段 `[94, 216)` 对、背景段 `[41, 88)` 被算成 `[40, 87)`。

残余差异（都记在 manifest.json 的 `residuals` 里，不藏）
------------------------------------------------------
1. `Timeinternal` 的全局性：Isoclock 拿**最后一个载入文件**的周期去换算**所有**文件的秒数。
   本工具按逐文件周期反解 ⇒ 只有当"最后一个文件"的周期与目标文件不同时下标才会漂。
   实测示例批次 83 个点的周期极差 2.5e-5 s（相对 8e-5），`selfcheck.ti_sensitivity`
   把"换成全批最小/最大周期"之后的最坏下标漂移如实报出来。
2. 本工具**只保证"同一批 cycle"**，不保证"同一套算法"。Isoclock 会做它自己的
   ±2σ 剔离群、普通铅扣除与分馏校正（`samples_method` / `Std*_method`），
   DRUID 也有一套。两边结果不同是**预期的**，正是要对照的东西。
3. `blank_window()` 是**本文件唯一重写的一处逻辑**（复现 `ablation.blank_mean` 取空白用的
   那个时间窗）。它有一条硬校验：用复现出来的窗口重算各通道空白，必须与 `Tra.blank`
   逐通道一致，否则**直接报错退出**（不允许"看起来差不多"）。
   `--from-stage` 那条路没有 `raw` / `blank` 可对账，这条校验跑不了 —— 在 JSON 里
   如实写成 `null`（**没验**，不是"通过"）。
4. `interval_cycles.csv` 的 `seg` 只按**时间**分段，不替使用者做地球化学判断；
   它标的是"这个 cycle 落在你选的哪一段里"，不是"这个 cycle 该不该用"。

用法
----
    python tools/export_isoclock.py --dir examples/EX2022A
    python tools/export_isoclock.py --dir <批次> --out <Isoclock 的 outputpath>
    python tools/export_isoclock.py --dir <批次> --merge-from <已有 Time_setting.npy>
    python tools/export_isoclock.py --from-stage <..._stage>          # 不重跑批次，读已有产物
    python tools/export_isoclock.py --dir <批次> --figure             # 另存一张时间点示意图

产物（`--out` 目录）
--------------------
    Time_setting.npy         → **拷到 Isoclock 的 outputpath**（覆盖它自己生成的那份）
    interval_selection.csv   逐点一行：六个时刻（秒）/ 四个下标 / 四个时间戳 / τ / 诊断
    interval_cycles.csv      逐 cycle 一行：`seg` 段标记 + `tau` + 三个 in_* 掩码
    interval_selection.json  schema、口径、**列名逐列说明**、参数、覆盖、自检、残余差异
    README.txt               上面那张"时间点命名与口径"表 + 三步怎么用（由本文件的
                             `TIME_POINT_LEGEND` / `COLUMN_LEGEND` / `CONVENTIONS` 生成，
                             不另抄一份）
    figure_interval.png      仅 `--figure`（且仅流水线路径）：每个测点一格，把六个时刻画在信号上

⚠ 边界（与仓库其它 `tools/` 脚本同一约定）
------------------------------------------
**只读**：不改 `druid/` 一个字节、不碰任何数值输出，产物只写到 `--out`。
区间值**全部**取自流水线自己的 `Tra`（`load_batch`），本工具一个数都不另算 ——
唯一的例外是上面第 3 条那个空白窗口，它带硬校验。
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from druid import __version__                            # noqa: E402
from druid.console import ensure_utf8_streams            # noqa: E402
from druid.core.constants import CJK_FONTS, ROLE_LABEL_CN  # noqa: E402
from druid.io.qtegra import read_qtegra                  # noqa: E402
from druid.qc import jsonable                            # noqa: E402
from druid.reduction.ablation import blank_mean          # noqa: E402
from druid.reduction.trace import window_mask            # noqa: E402
from druid.workflow import BatchConfig, load_batch       # noqa: E402

#: 本工具产物的契约版本。破坏性改动才升（改文件名 / 列名 / 语义），加列不升。
#: /2：列名从简写（`s0_druid_s` …）改成自解释的六个时刻名（`sample_start_s` …），
#:     并新增 `interval_cycles.csv`、四个 cycle 时间戳、τ 两列。改列名 ⇒ 必须升版。
SCHEMA = "druid.isoclock_bridge/2"

#: 浮点写出用几位有效数字。与 `export_stages.py` 同一约定：这些文件给机器算，
#: 截断到 4 位会引入看不见的偏差。
_FLOAT = "%.12g"

#: `ablation.blank_mean` 在空白段与剥蚀起点之间留的安全间隔 (s)。
#: 它**不是可调参数**，是从那份实现里抄出来的常量（那边写死了 0.5）。
_ISO_BLANK_GAP = 0.5

#: 覆盖检测时要跳过的文件：名单自身也是目录里的一个 `.csv`，它不是测点。
#: 名字由 `BatchConfig.list_file` 的发现规则决定（`<批次名>_LIST.csv` 等）。
_LIST_SUFFIXES = ("_list.csv", "_list.tsv", "_list.xls", "_list.xlsx")

#: `interval_cycles.csv` 的 `seg` 取值，**按时间先后**排列（同时也是自检里核对顺序的依据）。
SEGMENTS = ("before_blank", "blank", "before_ablation", "trimmed_head",
            "integration", "after_integration")

# ─────────────────────────────────────────────────────────────────────────────
# 口径与列名 —— 单一出处：JSON 的 `time_points` / `columns` / `conventions`
# 与 README 的同一张表都由这里生成，不另抄一份（抄一份就是"两份实现"）
# ─────────────────────────────────────────────────────────────────────────────
#: 六个时刻的语义契约。`druid` 列写清它对应 `Tra` 上的哪一个量 ——
#: 读表的人不用回去翻代码就知道这个数是怎么来的。
TIME_POINT_LEGEND = [
    {"name": "ablation_start_s", "druid": "Tra.t_ab0", "unit": "s",
     "meaning": "激光剥蚀起点（find_ablation 自动识别）"},
    {"name": "ablation_end_s", "druid": "Tra.t_ab1", "unit": "s",
     "meaning": "剥蚀终点"},
    {"name": "blank_start_s", "druid": "blank_mean 取窗下界（本工具复现）", "unit": "s",
     "meaning": "气体空白取样窗起点 = max(0, t_ab0 − blank_dur − 0.5)"},
    {"name": "blank_end_s", "druid": "blank_mean 取窗上界（本工具复现）", "unit": "s",
     "meaning": "气体空白取样窗终点 = t_ab0 − 0.5"},
    {"name": "sample_start_s", "druid": "Tra.t0", "unit": "s",
     "meaning": "积分窗口起点 = t_ab0 + trim（窗口过短时不裁边，见 trim_applied_s）"},
    {"name": "sample_end_s", "druid": "Tra.t1", "unit": "s",
     "meaning": "积分窗口终点"},
]

#: 三条口径。写进 JSON 与 README，省得下一个人踩。
CONVENTIONS = [
    "秒的**原点 = 该测点 CSV 时间轴自己的零点**，不是序列的绝对时刻；"
    "跨测点比较请用 cycle 下标或 τ，不要直接比秒。",
    "所有区间都是**闭区间**（window_mask 口径）[a, b] 两端取等号；"
    "Isoclock 那边是左闭右开，所以要写进 Time_setting.npy 的右端点落在『末下标 + 1』那一格 "
    "——这就是 iso_sample_end_s 与 sample_end_s 通常不是同一个数的原因。",
    "cycle 下标**从 0 开始**（np.flatnonzero 口径）；"
    "τ = (t − sample_start_s) / (sample_end_s − sample_start_s)，0 = 坑口、1 = 坑底"
    "（与 druid/depth/windows.py 逐字相同）。",
]

#: `interval_selection.csv` 的列说明（顺序即写出顺序）。
COLUMN_LEGEND = [
    ("idx", "批次内测点序号（从 1 开始，与结果表「序号」列一致）"),
    ("file", "数据文件名（**含 .csv**）—— 同时是 Time_setting.npy 的键"),
    ("sample", "样品名"),
    ("role", "内部角色码（primary_std / secondary_std / unknown / glass）"),
    ("role_cn", "角色的中文标签"),
    ("n_points", "该测点的 cycle 总数（= len(时间轴)，Isoclock 的 Numbers）"),
    ("period_s", "等效周期，Isoclock 的 Timeinternal = (t[-1] − t[0]) / n（**除 n 不是 n−1**）"),
    ("t_exact", "时间轴是不是流水线用的那一条（False = 等间隔重建，见 residuals）"),
    ("t_first_s", "时间轴首点时刻"),
    ("t_last_s", "时间轴末点时刻"),
    ("ablation_start_s", "★ 剥蚀起点 Tra.t_ab0"),
    ("ablation_end_s", "★ 剥蚀终点 Tra.t_ab1"),
    ("trim_applied_s", "实际裁掉的秒数 = sample_start_s − ablation_start_s（窗口过短时 < cfg.trim）"),
    ("blank_start_s", "★ 空白取样窗起点（闭区间）"),
    ("blank_end_s", "★ 空白取样窗终点（闭区间）"),
    ("sample_start_s", "★ 积分窗口起点 Tra.t0（闭区间）"),
    ("sample_end_s", "★ 积分窗口终点 Tra.t1（闭区间）"),
    ("cyc_blank_first", "空白窗内第一个 cycle 的下标（从 0 开始，闭区间）"),
    ("cyc_blank_last", "空白窗内末一个 cycle 的下标"),
    ("cyc_sample_first", "积分窗内第一个 cycle 的下标（从 0 开始，闭区间）"),
    ("cyc_sample_last", "积分窗内末一个 cycle 的下标"),
    ("cyc_blank_first_t_s", "上述首个空白 cycle 的**时间戳**（按秒取数的工具用它复现）"),
    ("cyc_blank_last_t_s", "上述末个空白 cycle 的时间戳"),
    ("cyc_sample_first_t_s", "上述首个积分 cycle 的时间戳"),
    ("cyc_sample_last_t_s", "上述末个积分 cycle 的时间戳"),
    ("tau_sample_first", "首个积分 cycle 的归一化深度 τ（≈0；>0 因为 cycle 未必正落在窗口端点）"),
    ("tau_sample_last", "末个积分 cycle 的归一化深度 τ（≈1）"),
    ("n_cyc_blank", "空白窗内的 cycle 数"),
    ("n_cyc_sample", "积分窗内的 cycle 数（Isoclock 那边是同一批，见 selfcheck.slice_match）"),
    ("iso_blank_start_s", "★ 写进 Time_setting.npy 的 b0（秒；默认取区间中点）"),
    ("iso_blank_end_s", "★ 写进 Time_setting.npy 的 b1"),
    ("iso_sample_start_s", "★ 写进 Time_setting.npy 的 s0"),
    ("iso_sample_end_s", "★ 写进 Time_setting.npy 的 s1（落在『末下标 + 1』那一格）"),
    ("edge_margin_s", "四个窗口端点离最近 cycle 的距离；越小越接近『下标会翻一格』的边界"),
    ("d_i0_raw", "**不吸附**时 Isoclock 会算出的起点下标与 DRUID 的差（默认吸附 ⇒ 0）"),
    ("d_i1_raw", "**不吸附**时右端点下标的差"),
]

#: `interval_cycles.csv` 的列说明。
CYCLE_COLUMN_LEGEND = [
    ("idx", "批次内测点序号（从 1 开始）"),
    ("file", "数据文件名（含 .csv）"),
    ("sample", "样品名"),
    ("role_cn", "角色的中文标签"),
    ("t_exact", "该行的时间轴是不是流水线用的那一条"),
    ("cyc", "cycle 下标（从 0 开始）"),
    ("t_s", "该 cycle 的时间戳（秒，原点是本文件时间轴零点）"),
    ("tau", "归一化深度 τ = (t_s − sample_start_s) / (sample_end_s − sample_start_s)；<1s 的窗口按 1e-9 兜底"),
    ("seg", "该 cycle 落在哪一段：" + " / ".join(SEGMENTS) + "（**互斥**，按时间先后）"),
    ("in_blank", "1 = 落在空白取样窗内（闭区间）"),
    ("in_sample", "1 = 落在积分窗口内（闭区间）—— 与 cyc_sample_first/last 那一对完全一致"),
    ("in_ablation", "1 = 落在剥蚀区间内 [ablation_start_s, ablation_end_s]"),
]

#: `interval_cycles.csv` 的 `seg` 语义（逐个说清，别让下一个人猜 "trimmed_head" 是什么）。
SEG_LEGEND = {
    "before_blank": "空白窗之前（t < blank_start_s）",
    "blank": "气体空白取样窗内（闭区间）",
    "before_ablation": "空白窗之后、剥蚀起点之前（t < ablation_start_s）—— 纯气体段",
    "trimmed_head": "剥蚀已经开始、但被 trim 裁掉的那一段（ablation_start_s ≤ t < sample_start_s）",
    "integration": "积分窗口内（闭区间）—— 真正参与整段比值的那一批 cycle",
    "after_integration": "积分窗口之后（t > sample_end_s）—— 含剥蚀尾部（若有）与剥蚀后气体",
}


# ─────────────────────────────────────────────────────────────────────────────
# 一、Isoclock 的算术（逐字对齐，不"等价改写"）
# ─────────────────────────────────────────────────────────────────────────────
def iso_period(t: np.ndarray) -> float:
    """
    Isoclock 的 `Timeinternal` —— **与它公式逐字相同**：`(x[N-1] - x[0]) / N`。

    ⚠ 是除以 `N` 而不是 `N-1`。看着像笔误，但这是它实际在用的数，
    我们要预测它的下标就必须照抄 —— "修正"它反而会让两边对不上。
    """
    n = len(t)
    return (float(t[-1]) - float(t[0])) / n


def iso_index(seconds: float, period: float) -> int:
    """Isoclock 的秒→下标：`int(b0 / Timeinternal)`（**截断**，不是四舍五入）。"""
    return int(seconds / period)


def iso_center(j: int, period: float) -> float:
    """
    反解：返回一个秒数，使 Isoclock 的 `int(秒 / period)` 恰好等于下标 `j`。

    解集是 `[j·T, (j+1)·T)`，这里取**中点**。为什么不用端点：
    `Timeinternal` 在 Isoclock 里是全局量（最后一个载入的文件赢），
    中点对"半个下标的周期误差"免疫，端点则不然。
    """
    return (j + 0.5) * period


def iso_slice(seconds: float, period: float) -> int:
    """`y[start:end]` 的那个下标 —— 与 `iso_index` 同式，单列出来是为了让调用点读起来是切片语义。"""
    return iso_index(seconds, period)


# ─────────────────────────────────────────────────────────────────────────────
# 二、DRUID 那一侧的事实（全部取自 Tra，不另算）
# ─────────────────────────────────────────────────────────────────────────────
def druid_cycle_span(t: np.ndarray, a: float, b: float):
    """
    DRUID 实际选中的 cycle 下标 —— **闭区间** `[a, b]`（`window_mask` 的口径：两边都取等号）。

    返回 `(i0, i1, n)`；窗口内一个点都没有时返回 `None`。
    用 `flatnonzero` 取两端，而不是假设数据在窗口内连续 —— 万一某天时间轴不单调，
    这里会如实报出首末下标，而不是悄悄给出一个错的跨度。
    """
    m = window_mask(t, a, b)
    if not m.any():
        return None
    i0 = int(np.flatnonzero(m)[0])
    i1 = int(len(m) - 1 - int(np.argmax(m[::-1])))
    return i0, i1, int(m.sum())


def blank_window(t: np.ndarray, t_ab0: float, dur: float):
    """
    复现 `reduction/ablation.blank_mean` 取气体空白用的那个时间窗（闭区间）。

    ⚠ **本文件唯一重写的一处逻辑**。那边只返回一个均值、不返回窗口，
    而这里需要窗口本身（要写进 Isoclock 的背景段）。既然重写了，
    就由 `check_blank_window()` 拿它重算各通道空白、与 `Tra.blank` 逐通道对账，
    对不上直接报错 —— 不允许出现"两份实现、症状为零"的情况。

    返回 `(lo, hi)`，含 `blank_mean` 的退化支路（可用点太少时从 0 取到 ta0-0.5）。
    """
    lo = max(0.0, t_ab0 - dur - _ISO_BLANK_GAP)
    hi = t_ab0 - _ISO_BLANK_GAP
    if int(((t >= lo) & (t <= hi)).sum()) < 5:
        lo = 0.0
    return lo, hi


def check_blank_window(tr, dur: float) -> dict:
    """
    用复现出来的空白窗口重算各通道空白，与 `Tra.blank` 逐通道对账。

    两条断言（任何一条不过 → 由调用方抛错）：
      ① 复现窗口的均值 == `Tra.blank[m]`（这才是在验"窗口复现对了"）；
      ② `blank_mean(...)` 的返回值 == `Tra.blank[m]`（验"调用姿势没写反"）。
    ② 看着多余，但它把"参数顺序写反了"和"窗口算错了"两种错分开报了。
    """
    lo, hi = blank_window(tr.t, tr.t_ab0, dur)
    m = (np.asarray(tr.t, float) >= lo) & (np.asarray(tr.t, float) <= hi)
    worst = 0.0
    worst_m = None
    for mass, arr in tr.raw.items():
        ref = float(tr.blank.get(mass, float("nan")))
        mine = float(np.asarray(arr, float)[m].mean()) if int(m.sum()) >= 2 else 0.0
        call = float(blank_mean(tr.t, arr, tr.t_ab0, dur))
        scale = max(abs(ref), 1.0)
        for got in (mine, call):
            d = abs(got - ref) / scale
            if d > worst:
                worst, worst_m = d, mass
    return {"blank_window_s": [lo, hi], "max_rel_err": worst, "worst_mass": worst_m,
            "ok": bool(worst <= 1e-9)}


def cycles_for(row: dict, t) -> list:
    """
    逐 cycle 一行 —— **把时间点位置标到逐 cycle 分辨率**，这就是"选一次区间、
    多处复算"的那份底表：换个分段方式直接按 `t_s` 重选即可，不必回 DRUID。

    分段判据与 `window_mask` 同一口径（闭区间、从 0 数），按时间先后 if/elif 落到
    唯一一个 `seg`；同时给出三个独立的 `in_*` 掩码（可重叠，用来看区间之间的关系）。

    τ 与 `depth/windows.py` 逐字相同（含 `max(span, 1e-9)` 那个兜底）。
    """
    t = np.asarray(t, float)
    a0, b0 = float(row["blank_start_s"]), float(row["blank_end_s"])
    a_ab, b_ab = float(row["ablation_start_s"]), float(row["ablation_end_s"])
    s0, s1 = float(row["sample_start_s"]), float(row["sample_end_s"])
    span = max(s1 - s0, 1e-9)
    out = []
    for i, ti in enumerate(t):
        ti = float(ti)
        if ti < a0:
            seg = "before_blank"
        elif ti <= b0:
            seg = "blank"
        elif ti < a_ab:
            seg = "before_ablation"
        elif ti < s0:
            seg = "trimmed_head"
        elif ti <= s1:
            seg = "integration"
        else:
            seg = "after_integration"
        out.append({
            "idx": int(row["idx"]), "file": row["file"], "sample": row["sample"],
            "role_cn": row["role_cn"], "t_exact": bool(row["t_exact"]),
            "cyc": int(i), "t_s": ti, "tau": (ti - s0) / span, "seg": seg,
            "in_blank": int(a0 <= ti <= b0),
            "in_sample": int(s0 <= ti <= s1),
            "in_ablation": int(a_ab <= ti <= b_ab),
        })
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 三、逐点算出行（秒数 + 下标 + 反解）
# ─────────────────────────────────────────────────────────────────────────────
def _snap(j: int, period: float, enabled: bool, raw_seconds: float) -> float:
    """把某个下标的秒数写出来：默认反解（吸附），`--no-snap` 时原样写 DRUID 的秒数。"""
    return iso_center(j, period) if enabled else float(raw_seconds)


def _row_from(meta: dict, t, t_exact: bool, t_ab0: float, t_ab1: float,
              t0: float, t1: float, dur: float, snap: bool, period: float) -> dict:
    """
    两个路径（流水线 / from-stage）共用的**唯一**一处行组装 —— 只写一遍，
    保证两条路产出的列名、口径、诊断列完全相同（"两份实现"就是这个模块最怕的东西）。
    """
    t = np.asarray(t, float)
    span = druid_cycle_span(t, t0, t1)
    if span is None:
        raise ValueError(f"{meta['file']}: 积分窗口 [{t0},{t1}] 内一个 cycle 都没有")
    i0, i1, n_druid = span
    lo, hi = blank_window(t, t_ab0, dur)
    bspan = druid_cycle_span(t, lo, hi)
    if bspan is None:
        raise ValueError(f"{meta['file']}: 空白窗口 [{lo},{hi}] 内一个 cycle 都没有")
    bi0, bi1, n_blk = bspan
    tau0, tau1 = float(t0), float(t1)
    tspan = max(tau1 - tau0, 1e-9)
    return {
        **meta,
        "n_points": int(len(t)),
        "period_s": float(period),
        "t_exact": bool(t_exact),
        "t_first_s": float(t[0]),
        "t_last_s": float(t[-1]),
        # ── 六个时刻（名字自解释；对应哪个 Tra 量见 TIME_POINT_LEGEND）──
        "ablation_start_s": float(t_ab0),
        "ablation_end_s": float(t_ab1),
        "trim_applied_s": float(t0 - t_ab0),
        "blank_start_s": float(lo),
        "blank_end_s": float(hi),
        "sample_start_s": float(t0),
        "sample_end_s": float(t1),
        # ── 四个 cycle 下标（从 0 开始，闭区间）──
        "cyc_blank_first": int(bi0), "cyc_blank_last": int(bi1),
        "cyc_sample_first": int(i0), "cyc_sample_last": int(i1),
        # ── 四个 cycle 的**时间戳**：按秒取数的工具靠这四个数复现同一段 ──
        "cyc_blank_first_t_s": float(t[bi0]), "cyc_blank_last_t_s": float(t[bi1]),
        "cyc_sample_first_t_s": float(t[i0]), "cyc_sample_last_t_s": float(t[i1]),
        # ── τ（与 depth/windows.py 同式；逐 cycle 用时间戳，窗口用中点）──
        "tau_sample_first": float((t[i0] - tau0) / tspan),
        "tau_sample_last": float((t[i1] - tau0) / tspan),
        "n_cyc_blank": int(n_blk), "n_cyc_sample": int(n_druid),
        # ── ★ 要写进 Time_setting.npy 的秒数 ──
        #   s1/b1 用 `下标 + 1` 再取中点：Isoclock 那边是**左闭右开**切片，
        #   要让末 cycle 被包含进来，右端点必须落在"末下标 + 1"那一格里。
        "iso_blank_start_s": _snap(bi0, period, snap, lo),
        "iso_blank_end_s": _snap(bi1 + 1, period, snap, hi),
        "iso_sample_start_s": _snap(i0, period, snap, t0),
        "iso_sample_end_s": _snap(i1 + 1, period, snap, t1),
        # ── 诊断：**不吸附**时 Isoclock 会算出什么（差值就是"吸附"修掉的那部分）──
        "d_i0_raw": iso_index(t0, period) - i0,
        "d_i1_raw": iso_index(t1, period) - (i1 + 1),
        # ── 窗口端点离最近 cycle 有多远：默认路径下只是信息量（下标是从 mask 数出来的），
        #    真正起作用的是 --from-stage 且时间轴只能等间隔重建那条路（见 selfcheck.time_axis）──
        "edge_margin_s": float(min(float(np.min(np.abs(t - x)))
                                   for x in (t0, t1, lo, hi))),
    }


def rows_from_spots(spots, cfg: BatchConfig, snap: bool):
    """
    从 `load_batch` 的产物（`Tra`）里取全部区间事实，逐点算好该写什么。

    返回 `(rows, cycles, blank_checks)`。`blank_checks` 是空白窗口的对账结果 ——
    只要有一条 `ok=False`，调用方必须报错退出（那是"两份实现"的信号）。
    """
    rows, cycles, blank_checks = [], [], []
    for tr, _r in spots:
        t = np.asarray(tr.t, float)
        bc = check_blank_window(tr, cfg.blank_dur)
        bc["file"] = Path(tr.path).name
        blank_checks.append(bc)
        meta = {"idx": int(tr.idx) + 1, "file": Path(tr.path).name, "sample": tr.sample,
                "role": tr.role, "role_cn": ROLE_LABEL_CN.get(tr.role, tr.role)}
        row = _row_from(meta, t, True, float(tr.t_ab0), float(tr.t_ab1),
                        float(tr.t0), float(tr.t1), cfg.blank_dur, snap, iso_period(t))
        rows.append(row)
        cycles.extend(cycles_for(row, t))
    return rows, cycles, blank_checks


def _resolve_csv_name(data_dir: Path, given: str) -> str:
    """
    把产物里的 `file` 值还原成**磁盘上的真实文件名**。

    为什么要这一步：`export_stages.py` 的 `file` 列是 `Path(tr.path).stem`
    （**不带 `.csv`**，与序列表的约定一致），而 Isoclock 的键是
    `os.walk` 拿到的**文件名（含 `.csv`）** —— 差一个后缀，`numn[name]` 整条都取不到。
    这里先去数据目录里核对真实文件名，核不上再补 `.csv`（归档场景下数据目录可能不在）。
    """
    name = Path(str(given)).name
    if Path(name).suffix.lower() == ".csv":
        return name
    d = Path(data_dir)
    if d.is_dir():
        for p in d.iterdir():
            if p.is_file() and p.stem == Path(name).stem and p.suffix.lower() == ".csv":
                return p.name
    return name + ".csv"


def _stage_time_axis(data_dir: Path, fname: str, t_first: float, t_last: float, n: int):
    """
    `--from-stage` 要算下标就得有时间轴，但那条路手上只有四个时刻。

    分两级取（**能精确就不近似**）：
      ① 原始 CSV 还在 → 用 `io.qtegra.read_qtegra` 把 `t` 读回来（**精确**，
         与流水线用的是同一个解析器；读时间轴不是"重跑还原"，不动任何数值）；
      ② 原始 CSV 不在了（归档场景）→ 用 `t_first + i·step` 等间隔重建（**近似**）。

    ⚠ ②为什么只是"近似"：实测时间轴**不是**严格等间隔，重建与真实 `t` 的最大偏差
    达 2.0e-3 s（≈0.64% 个 cycle，iCAP 的时间戳本身就有毫秒级抖动）⇒ 窗口端点
    一旦落进"距某个 cycle 边界 2 ms 以内"的危险带，下标就可能差 1。
    返回 `(t, exact)`，由调用方把这件事如实报出去。

    端点校验：读回来的 `t` 若长度或首末时刻与产物对不上，**不采用**，退回重建 ——
    宁可近似并声明，也不拿一份来历不符的时间轴去算下标。
    """
    cands = [Path(data_dir) / fname]
    if not fname.lower().endswith(".csv"):
        cands.append(Path(data_dir) / (fname + ".csv"))
    for cand in cands:
        if not cand.exists():
            continue
        try:
            _title, t, _raw, _cols, _info = read_qtegra(cand)
        except Exception:                                    # noqa: BLE001
            break
        t = np.asarray(t, float)
        if (len(t) == n and abs(float(t[0]) - t_first) < 1e-6
                and abs(float(t[-1]) - t_last) < 1e-6):
            return t, True
        break
    step = (t_last - t_first) / max(n - 1, 1)
    return t_first + np.arange(n) * step, False


def rows_from_stage(csv_path: Path, cfg: BatchConfig, snap: bool, data_dir: Path):
    """
    从 `export_stages.py` 的 `01_blank_and_window.csv` 读区间，**不重跑批次**。

    那条路没有 `Tra`，所以时间轴按 `_stage_time_axis` 分两级取（能读原始 CSV 就精确），
    取到以后**走与流水线路径完全相同的下标逻辑**（同一个 `window_mask` / 同一个
    `blank_window` / 同一个 `_row_from`）。拿不到精确时间轴时，`t_exact=False`，并由
    `selfcheck.time_axis` 报出有多少个点落在"±1 cycle"的危险带里。
    空白窗口无法与 `Tra.blank` 对账（没有 `raw`/`blank`）—— 这一条会写成 `null`，
    **不是"通过"**。
    """
    df = pd.read_csv(csv_path, encoding="utf-8-sig")
    need = {"file", "n_cycles", "t_first_s", "t_last_s", "t_ab0_s", "t0_s", "t1_s"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{csv_path} 缺少列 {sorted(missing)} —— 请用新版 "
                         f"tools/export_stages.py 重新生成")

    rows, cycles = [], []
    for _k, r in df.iterrows():
        n = int(r["n_cycles"])
        t_first, t_last = float(r["t_first_s"]), float(r["t_last_s"])
        period = (t_last - t_first) / n
        # ⚠ 键必须是**磁盘上的文件名（含 .csv）** —— 产物里那一列不带后缀
        fname = _resolve_csv_name(data_dir, str(r["file"]))
        t, t_exact = _stage_time_axis(data_dir, fname, t_first, t_last, n)
        role = str(r.get("role", ""))
        meta = {"idx": int(r["idx"]), "file": fname, "sample": str(r["sample"]),
                "role": role, "role_cn": ROLE_LABEL_CN.get(role, role)}
        row = _row_from(meta, t, t_exact, float(r["t_ab0_s"]),
                        float(r.get("t_ab1_s", r["t_ab0_s"])),
                        float(r["t0_s"]), float(r["t1_s"]), cfg.blank_dur, snap, period)
        rows.append(row)
        cycles.extend(cycles_for(row, t))
    return rows, cycles


# ─────────────────────────────────────────────────────────────────────────────
# 四、组装 Time_setting.npy + 自检
# ─────────────────────────────────────────────────────────────────────────────
def build_time_setting(rows) -> dict:
    """
    按 Isoclock 的 `num[name] = [Numbers, b0, b1, s0, s1]` 组装。

    `Numbers` 填 `len(x)`（实测它只在计算 `Timeinternal` 时用过，之后没人读，
    但既然它写在那儿，我们就填对）。
    """
    out = {}
    for r in rows:
        out[r["file"]] = [int(r["n_points"]),
                          float(r["iso_blank_start_s"]), float(r["iso_blank_end_s"]),
                          float(r["iso_sample_start_s"]), float(r["iso_sample_end_s"])]
    return out


def iso_check(numn: dict) -> list:
    """
    复现 Isoclock `check()` 的两条判据（行 1325-1337），返回不合规的键。

    ⚠ `check()` 会把不合规的条目**就地改写**（把区间压塌成零长）再存回去 ——
    也就是说入口不拦的话，坏值会被静默抹平。所以这里**先拦住**，而不是交给它。
    """
    bad = []
    for name, v in numn.items():
        if v[3] - v[2] < 0 or v[4] - v[2] < 0:
            bad.append((name, "s0-b1<0 或 s1-b1<0"))
        elif v[3] + 1 > v[4] - 1:
            bad.append((name, "s1-s0<=2 会被压塌"))
    return bad


def cycles_match(rows, cycles) -> dict:
    """
    核对逐 cycle 表与逐点表**说的是不是同一段**（这是"标注清楚"的底线）。

    三条：`integration` 段的 cycle 数 == `n_cyc_sample`；其首末 `cyc` ==
    `cyc_sample_first/last`；`blank` 段同理。第二、三条尤其重要 —— 逐点表的四个下标
    是 `window_mask` 数出来的，逐 cycle 表是 if/elif 判出来的，**两条独立路径**，
    对不上就说明分段判据写岔了（而且症状为零）。
    """
    bad = []
    by_file = {}
    for c in cycles:
        by_file.setdefault(c["file"], []).append(c)
    seg_index = {s: i for i, s in enumerate(SEGMENTS)}
    for r in rows:
        cs = by_file.get(r["file"], [])
        if len(cs) != r["n_points"]:
            bad.append({"file": r["file"], "why": f"cycle 行数 {len(cs)} != n_points {r['n_points']}"})
            continue
        # seg 取值合法、且按时间单调（不然 cycles_for 的 if/elif 就没意义）
        unknown = sorted({c["seg"] for c in cs} - set(seg_index))
        if unknown:
            bad.append({"file": r["file"], "why": f"seg 出现未定义取值 {unknown}"})
        elif any(seg_index[cs[i]["seg"]] > seg_index[cs[i + 1]["seg"]]
                 for i in range(len(cs) - 1)):
            bad.append({"file": r["file"], "why": "seg 不是按时间单调的（时间轴可能非单调）"})
        for tag, c0, c1, n_key in (("integration", r["cyc_sample_first"], r["cyc_sample_last"],
                                    "n_cyc_sample"),
                                   ("blank", r["cyc_blank_first"], r["cyc_blank_last"],
                                    "n_cyc_blank")):
            hit = [c["cyc"] for c in cs if c["seg"] == tag]
            if len(hit) != r[n_key] or (hit and (hit[0] != c0 or hit[-1] != c1)):
                bad.append({"file": r["file"], "seg": tag, "n_cycles": len(hit),
                            "want_n": r[n_key], "cyc": [hit[0], hit[-1]] if hit else None,
                            "want_cyc": [c0, c1]})
        hit = [c["cyc"] for c in cs if c["in_sample"]]
        if len(hit) != r["n_cyc_sample"] or (hit and (hit[0] != r["cyc_sample_first"]
                                                      or hit[-1] != r["cyc_sample_last"])):
            bad.append({"file": r["file"], "why": "in_sample 掩码与 cyc_sample_first/last 不一致"})
    return {"ok": not bad, "bad": bad[:10],
            "how": "逐 cycle 表的 seg=='integration' / in_sample 掩码，与逐点表的 "
                   "cyc_sample_first/last + n_cyc_sample 逐点对齐（两条独立路径互证）"}


def selfcheck(rows, cycles, numn: dict, checks: list) -> dict:
    """
    自检：**写进去的秒数，按 Isoclock 的算术解出来的下标，必须与 DRUID 的 cycle 下标一致。**

    这一步看着像同义反复（秒数本来就是反解出来的），但它挡住的是**浮点边界**：
    反解用的是 `(j+0.5)·T`，若某个 `j` 恰好落在浮点误差里，
    `int(把秒数除以 T)` 会掉到隔壁一格 —— 那正是"两边区间差一个 cycle"这种
    零症状错误。另外还量化了 `Timeinternal` 全局性带来的风险（`ti_sensitivity`），
    以及逐 cycle 表与逐点表是否互证（`cycles_match`）。
    """
    bad_slice, bad_order = [], []
    bad_iso = iso_check(numn)
    ti_lo = min(r["period_s"] for r in rows)
    ti_hi = max(r["period_s"] for r in rows)
    worst_ti = {"period_lo": ti_lo, "period_hi": ti_hi, "max_abs_d": 0, "worst_file": None}

    for r in rows:
        v = numn[r["file"]]
        b0, b1, s0, s1 = v[1], v[2], v[3], v[4]
        T = r["period_s"]
        # ① Isoclock 会切出哪一段？（左闭右开）
        sb, eb = iso_slice(b0, T), iso_slice(b1, T)
        ss, es = iso_slice(s0, T), iso_slice(s1, T)
        if (ss != r["cyc_sample_first"]) or (es != r["cyc_sample_last"] + 1) \
                or (sb != r["cyc_blank_first"]) or (eb != r["cyc_blank_last"] + 1):
            bad_slice.append({"file": r["file"],
                              "want_sample": [r["cyc_sample_first"], r["cyc_sample_last"] + 1],
                              "got_sample": [ss, es],
                              "want_blank": [r["cyc_blank_first"], r["cyc_blank_last"] + 1],
                              "got_blank": [sb, eb]})
        # ② 切片得是有序且非空
        if not (es > ss) or not (eb > sb):
            bad_order.append({"file": r["file"], "sample": [ss, es], "blank": [sb, eb]})
        # ③ 周期换成全批最小/最大时的最坏漂移
        for T2 in (ti_lo, ti_hi):
            for got, want in ((iso_slice(s0, T2), r["cyc_sample_first"]),
                              (iso_slice(s1, T2), r["cyc_sample_last"] + 1)):
                d = abs(got - want)
                if d > worst_ti["max_abs_d"]:
                    worst_ti.update(max_abs_d=int(d), worst_file=r["file"],
                                    at=("period_lo" if T2 == ti_lo else "period_hi"))

    bad_blank = [c for c in checks if not c["ok"]]
    # 近似时间轴（只可能出现在 --from-stage 且原始 CSV 不在时）：
    # 重建与真实 t 的偏差实测上限 2.0e-3 s ⇒ 端点落进这个带宽就可能差 1 个 cycle。
    approx = [r for r in rows if not r.get("t_exact", True)]
    band = 2.5e-3
    danger = [r["file"] for r in approx if r["edge_margin_s"] < band]
    return {
        "n_spots": len(rows),
        "time_axis": {
            "n_exact": len(rows) - len(approx),
            "n_reconstructed": len(approx),
            "n_in_danger_band": len(danger),
            "band_s": band,
            "danger_files": danger[:10],
            "note": "重建时间轴的偏差上限 2.0e-3 s；端点距 cycle 边界小于这个带宽时下标可能差 1 "
                    "（即两边区间差一个 cycle）。原始 CSV 在时不会走到这一步。",
        },
        "slice_match": {"ok": not bad_slice, "bad": bad_slice[:10],
                        "rule": "int(s0/T)==i0 且 int(s1/T)==i1+1（Isoclock 左闭右开）"
                                " —— 两侧下标必须完全相等，不容许 ±1"},
        "slice_ordered": {"ok": not bad_order, "bad": bad_order[:10]},
        "cycles_match": (cycles_match(rows, cycles) if cycles
                         else {"ok": None, "note": "--no-cycles：逐 cycle 表没写，这一项**没验**"}),
        "blank_window_reproduced": {
            # ⚠ `--from-stage` 那条路拿不到逐 cycle 数据 ⇒ 这一项**没验**，
            #   所以是 None 而不是 True —— 不能把"没验"写成"通过"。
            "ok": (None if not checks else (not bad_blank)),
            "n_checked": len(checks),
            "max_rel_err": max((c["max_rel_err"] for c in checks), default=None),
            "bad": [{"file": c["file"], "why": f"max_rel_err={c['max_rel_err']:.3g}"}
                    for c in bad_blank[:10]],
            "how": "用复现出来的空白窗口重算各通道均值，与 Tra.blank 逐通道对账（相对误差 <= 1e-9）",
        },
        "iso_check": {"ok": not bad_iso, "bad": bad_iso[:10],
                      "note": "Isoclock check() 会把不合规条目就地压塌再存回去 —— 这里先拦住"},
        "ti_sensitivity": worst_ti,
        "d_i_raw": {"n_nonzero_i0": int(sum(1 for r in rows if r["d_i0_raw"])),
                    "n_nonzero_i1": int(sum(1 for r in rows if r["d_i1_raw"])),
                    "note": "不吸附（--no-snap）时 Isoclock 会算出的下标与 DRUID 差多少；"
                            "默认吸附已把这两列清零"},
    }


# ─────────────────────────────────────────────────────────────────────────────
# 五、覆盖情况：目录里的 CSV 有没有全都被写到
# ─────────────────────────────────────────────────────────────────────────────
def coverage(data_dir: Path, numn: dict) -> dict:
    """
    列出"目录里有、但这份 Time_setting.npy 没有"的 CSV。

    为什么必须报：Isoclock 的 `dataprocess()` 是
    `for name in sampleslist.keys(): b0 = numn[name][1] …` ——
    **样本清单里少一个键就直接 KeyError 停住**。而 DRUID 会跳过玻璃标样与作废点，
    所以"DRUID 载入的测点"天然少于"目录里的 CSV"（示例批次 83 vs 86）。
    这批差值必须让使用者看得见，或者用 `--merge-from` 把 Isoclock 自己
    生成的那几条保底留下来。
    """
    on_disk = sorted(p.name for p in data_dir.glob("*.csv")
                     if not p.name.lower().endswith(_LIST_SUFFIXES))
    have = set(numn)
    return {"csv_on_disk": len(on_disk),
            "covered": len([n for n in on_disk if n in have]),
            "missing": [n for n in on_disk if n not in have],
            "extra_keys": sorted(set(numn) - set(on_disk))}


# ─────────────────────────────────────────────────────────────────────────────
# 六、示意图（可选）：把六个时刻画在信号上，肉眼核对
# ─────────────────────────────────────────────────────────────────────────────
def make_figure(pairs, out_path: Path, panels: int = 6) -> bool:
    """
    一张"时间点落在哪"的示意图 —— **仅流水线路径可用**（只有那条路手上还有逐 cycle raw）。

    `pairs` 是 `[(Tra, row), …]`。蓝底 = 气体空白窗，红底 = 积分窗；
    竖线 = 六个时刻。这张图不是产物的一部分，是给人**核对**用的：
    填进 Isoclock 之前扫一眼，比读 33 列数字快得多。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = CJK_FONTS
    plt.rcParams["axes.unicode_minus"] = False

    n = len(pairs)
    if n == 0:
        return False
    if n <= panels:
        pick = list(range(n))
    else:
        pick = sorted({int(round(k * (n - 1) / (panels - 1))) for k in range(panels)})
    fig, axes = plt.subplots(len(pick), 1, figsize=(11.0, 2.3 * len(pick)), squeeze=False)
    for j, k in enumerate(pick):
        tr, row = pairs[k]
        t = np.asarray(tr.t, float)
        ax = axes[j][0]
        for mass, color in ((238, "#1f6feb"), (206, "#c2410c")):
            if mass in tr.raw:
                y = np.asarray(tr.raw[mass], float)
                ax.plot(t, np.where(y > 0, y, np.nan), lw=0.8, color=color,
                        label=f"{mass} cps")
        ax.set_yscale("log")
        ax.axvspan(row["blank_start_s"], row["blank_end_s"], color="#4c9be8", alpha=0.16)
        ax.axvspan(row["sample_start_s"], row["sample_end_s"], color="#d94a4a", alpha=0.16)
        for x, color, ls in ((row["blank_start_s"], "#4c9be8", ":"),
                             (row["blank_end_s"], "#4c9be8", ":"),
                             (row["ablation_start_s"], "#15803d", "--"),
                             (row["ablation_end_s"], "#15803d", "--"),
                             (row["sample_start_s"], "#d94a4a", "-"),
                             (row["sample_end_s"], "#d94a4a", "-")):
            ax.axvline(x, color=color, lw=1.0, ls=ls)
        ax.set_title(
            f"#{row['idx']} {row['sample']} ({row['role_cn']})   "
            f"cyc_sample={row['cyc_sample_first']}..{row['cyc_sample_last']} "
            f"(n={row['n_cyc_sample']})   cyc_blank={row['cyc_blank_first']}.."
            f"{row['cyc_blank_last']} (n={row['n_cyc_blank']})",
            fontsize=9, loc="left")
        ax.set_xlim(float(t[0]), float(t[-1]))
        ax.set_ylabel("cps", fontsize=8)
        ax.tick_params(labelsize=8)
        if j == 0:
            ax.legend(fontsize=7, loc="upper right", ncol=2)
    axes[-1][0].set_xlabel("时间 (s)     蓝底 = 气体空白窗　红底 = 积分窗", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# 七、README（由上面的口径表生成，不另抄一份）
# ─────────────────────────────────────────────────────────────────────────────
def build_readme() -> str:
    """README 的正文 —— 时刻表与列说明**直接引用**上面的常量，改口径只需改一处。"""
    w = max(len(p["name"]) for p in TIME_POINT_LEGEND)
    lines = [
        "这份目录由 tools/export_isoclock.py 生成 —— 把 DRUID 选好的积分区间交给 Isoclock，",
        f"同时把区间本身写成与软件无关、带完整标注的机读表（schema = {SCHEMA}）。",
        "",
        "【怎么用（三步）】",
        "  1. 在 Isoclock 里跑**加载**步。它会做两件事：把每个 CSV 存成 <名字>_x.npy / _y1..y7.npy，",
        "     并把**它自己**自动识别的区间写进 <outputpath>\\Time_setting.npy。",
        "  2. 把本目录的 Time_setting.npy **拷进 Isoclock 的 outputpath、覆盖**它那份",
        "     （或者本工具直接 --out <Isoclock 的 outputpath>，省掉拷贝）。",
        "  3. 跑后面的处理步。它读的就是这份文件，于是两边的积分区间是同一段。",
        "",
        "【时间点总览（秒；原点 = 该测点自己时间轴的零点）】",
    ]
    lines += [f"  {p['name']:<{w}}  {p['unit']:<2}  ← {p['druid']:<34}  {p['meaning']}"
              for p in TIME_POINT_LEGEND]
    lines += ["", "【三条口径】"]
    lines += [f"  {i}. {c}" for i, c in enumerate(CONVENTIONS, 1)]
    wc = max(len(c) for c, _ in COLUMN_LEGEND)
    lines += ["", "【interval_selection.csv 的列】"]
    lines += [f"  {c:<{wc}}  {d}" for c, d in COLUMN_LEGEND]
    wq = max(len(c) for c, _ in CYCLE_COLUMN_LEGEND)
    lines += ["", "【interval_cycles.csv 的列（逐 cycle；'选一次区间做多种运算'的底表）】"]
    lines += [f"  {c:<{wq}}  {d}" for c, d in CYCLE_COLUMN_LEGEND]
    lines += ["", "【interval_cycles.csv 的 seg 取值】"]
    wseg = max(len(s) for s in SEG_LEGEND)
    lines += [f"  {s:<{wseg}}  {SEG_LEGEND[s]}" for s in SEGMENTS]
    lines += [
        "",
        "【两个坑】",
        "  · 再跑一次**加载**步会把这份覆盖掉（instructure0/1 结尾那句 np.save）。",
        "  · DRUID 会跳过玻璃标样与作废点，所以写进去的条目**天然少于目录里的 CSV**",
        "    （示例批次 83 vs 86）。而 Isoclock 的 dataprocess() 遇到清单里有、",
        "    Time_setting 里没有的键会直接 KeyError。用 --merge-from 把 Isoclock",
        "    自己生成的那几条保底留下来，就不会断。",
        "",
        "【残余差异】",
        "  · Isoclock 的 Timeinternal 是全局量（最后一个载入的文件赢），本工具按逐文件周期反解；",
        "    interval_selection.json 的 selfcheck.ti_sensitivity 给出了最坏漂移。",
        "  · 本工具只保证**同一批 cycle**，不保证同一套算法 —— 两边结果不同是预期的。",
        "  · 逐 cycle 表的 seg 只按**时间**分段，不替你做地球化学判断。",
        "",
        "  详见 interval_selection.json（schema / 口径 / 逐列说明 / 自检 / 残余差异）、",
        "  interval_selection.csv（逐点）与 interval_cycles.csv（逐 cycle）。",
        "",
    ]
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# 八、主入口
# ─────────────────────────────────────────────────────────────────────────────
def _load_existing(path: Path):
    """读一份已有的 Time_setting.npy（Isoclock 自己生成的或上一次本工具生成的）。

    ⚠ 用 `allow_pickle=True` 读 pickle 是会执行代码的 —— 这里只读**你自己**那份
    由 Isoclock / 本工具写的文件，不要指向来路不明的 .npy。
    """
    try:
        d = np.load(path, allow_pickle=True).item()
    except Exception as e:                                   # noqa: BLE001
        raise SystemExit(f"读不了 {path}：{type(e).__name__}: {e}")
    if not isinstance(d, dict):
        raise SystemExit(f"{path} 里不是字典（是 {type(d).__name__}）—— 不像是 Time_setting.npy")
    return d


def main(argv=None) -> int:
    ensure_utf8_streams()
    ap = argparse.ArgumentParser(
        description="把 DRUID 自动选好的积分区间写成 Isoclock 的 Time_setting.npy，"
                    "并另存带完整标注的逐点表 / 逐 cycle 表（只读：不改 druid/ 任何代码、"
                    "不改任何数值输出）。")
    ap.add_argument("--dir", default=str(ROOT / "examples" / "EX2022A"),
                    help="批次目录（含 *_LIST.csv 与逐点 CSV）")
    ap.add_argument("--out", default=None,
                    help="输出目录；默认「批次目录上一级/<批次名>_isoclock」。"
                         "也可以直接指到 Isoclock 的 outputpath")
    ap.add_argument("--from-stage", default=None,
                    help="从 export_stages.py 的 01_blank_and_window.csv 读区间，**不重跑批次**")
    ap.add_argument("--merge-from", default=None,
                    help="已有的 Time_setting.npy（Isoclock 生成的那份）：把它有而我没写的条目留下来")
    ap.add_argument("--blank-dur", type=float, default=None,
                    help="气体空白取样时长 (s)，默认取流水线默认值；--from-stage 时不受 manifest 影响")
    ap.add_argument("--no-snap", action="store_true",
                    help="不反解秒数，原样写 DRUID 的秒数（会让 Isoclock 的切片与 DRUID 差 <=1 cycle）")
    ap.add_argument("--no-cycles", action="store_true",
                    help="不写 interval_cycles.csv（大批次省体积；代价是少一层互证）")
    ap.add_argument("--figure", action="store_true",
                    help="另存 figure_interval.png（把六个时刻画在信号上；仅流水线路径可用）")
    ap.add_argument("--figure-panels", type=int, default=6,
                    help="示意图放几个测点（默认 6，按批次均匀取）")
    args = ap.parse_args(argv)

    cfg = BatchConfig(data_dir=str(args.dir), plot=False, verbose=False)
    if args.blank_dur is not None:
        cfg.blank_dur = float(args.blank_dur)
    snap = not args.no_snap

    data_dir = Path(args.dir)
    out = (Path(args.out) if args.out
           else data_dir.resolve().parent / f"{data_dir.name}_isoclock")
    out.mkdir(parents=True, exist_ok=True)

    # ── 取区间 ──
    spot_list = None
    if args.from_stage:
        stage = Path(args.from_stage)
        csv_path = stage / "01_blank_and_window.csv" if stage.is_dir() else stage
        print(f"[1] 从既有产物读区间：{csv_path}")
        rows, cycles = rows_from_stage(csv_path, cfg, snap, data_dir)
        checks = []
        blank_src = f"assumed(default={cfg.blank_dur} s)"
        n_approx = sum(1 for r in rows if not r["t_exact"])
        print(f"    读入 {len(rows)} 个点；时间轴：精确 {len(rows) - n_approx} 个、"
              f"等间隔重建 {n_approx} 个（原始 CSV {'在' if n_approx == 0 else '不全在'}）")
    else:
        print(f"[1] 装载 {cfg.data_dir}（跑一遍流水线前置步骤，取 Tra 上的区间事实）")
        spot_list, _seq, _skipped = load_batch(cfg)
        if not spot_list:
            raise SystemExit(f"没有任何测点成功还原：{cfg.data_dir}")
        rows, cycles, checks = rows_from_spots(spot_list, cfg, snap)
        bad_blank = [c for c in checks if not c["ok"]]
        if bad_blank:
            worst = max(bad_blank, key=lambda c: c["max_rel_err"])
            raise SystemExit(
                f"空白窗口复现失败：{worst['file']} 上相对误差 {worst['max_rel_err']:.3g}"
                f"（通道 {worst['worst_mass']}）—— ablation.blank_mean 的取窗规则变了，"
                f"本工具必须跟着改，**不写文件**。")
        blank_src = "pipeline(cfg.blank_dur)"
        print(f"    有效测点 {len(rows)} 个；空白窗口已与 Tra.blank 逐通道对账通过"
              f"（最大相对误差 {max(c['max_rel_err'] for c in checks):.2e}）")
    if args.no_cycles:
        cycles = []

    numn_new = build_time_setting(rows)

    # ── 合并已有的条目 ──
    merge_path = Path(args.merge_from) if args.merge_from else out / "Time_setting.npy"
    kept = {}
    if merge_path.exists():
        existing = _load_existing(merge_path)
        kept = {k: v for k, v in existing.items() if k not in numn_new}
        overlap = len(existing) - len(kept)
        print(f"[2] 合并 {merge_path}：保留它独有的 {len(kept)} 条、覆盖同名 {overlap} 条")
        numn = {**kept, **numn_new}
    else:
        numn = dict(numn_new)
        print(f"[2] 没有可合并的 {merge_path}，只写本次的 {len(numn)} 条")

    # ── 自检：不过就不写文件 ──
    chk = selfcheck(rows, cycles, numn, checks)
    for key in ("slice_match", "slice_ordered", "cycles_match",
                "blank_window_reproduced", "iso_check"):
        if chk[key]["ok"] is False:            # None = 这条路没验（stage 路径 / --no-cycles），不是失败
            raise SystemExit(f"自检不通过：{key} -> {json.dumps(chk[key], ensure_ascii=False)[:400]}"
                             f"\n**没有写任何文件**。")
    print(f"[3] 自检通过：{len(rows)} 个点的下标逐点与 DRUID 的 cycle 区间相等；"
          f"周期换用全批极值时最坏漂移 {chk['ti_sensitivity']['max_abs_d']} cycle；"
          f"逐 cycle 表与逐点表互证 "
          f"{'通过' if chk['cycles_match']['ok'] else '未验（--no-cycles）'}")
    if chk["time_axis"]["n_reconstructed"]:
        print(f"    ⚠ 其中 {chk['time_axis']['n_reconstructed']} 个点用的是等间隔重建的时间轴"
              f"（原始 CSV 不在），{chk['time_axis']['n_in_danger_band']} 个落在 ±1 cycle 的危险带里")

    # ── 落盘 ──
    np.save(out / "Time_setting", numn)
    df = pd.DataFrame(rows)
    df.to_csv(out / "interval_selection.csv", index=False,
              encoding="utf-8-sig", float_format=_FLOAT)
    dfc = pd.DataFrame(cycles)
    if not dfc.empty:
        dfc.to_csv(out / "interval_cycles.csv", index=False,
                   encoding="utf-8-sig", float_format=_FLOAT)
    fig_path = None
    if args.figure:
        if args.from_stage or not spot_list:
            print("    ⚠ --figure 只在流水线路径可用（--from-stage 那条路手上没有逐 cycle "
                  "原始信号），已跳过")
        else:
            fig_path = out / "figure_interval.png"
            if make_figure(list(zip([s[0] for s in spot_list], rows)), fig_path,
                           args.figure_panels):
                print(f"    示意图已写出：{fig_path.name}（{args.figure_panels} 个测点）")

    cov = coverage(data_dir, numn)
    payload = {
        "schema": SCHEMA,
        "schema_note": "/2：列名从简写改成自解释的六个时刻名（s0_druid_s → sample_start_s 等），"
                       "并新增 interval_cycles.csv、四个 cycle 时间戳与 tau_sample_first/last。",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "generated_by": "tools/export_isoclock.py",
        "tool": {"name": "druid", "version": __version__},
        "target": {
            "software": "Isoclock (github.com/sndjgm/Isoclock)",
            "file": "Time_setting.npy",
            "entry_format": "[Numbers, b0, b1, s0, s1]  （秒；键 = 数据文件名含 .csv）",
            "place_into": "Isoclock 的 outputpath，**加载步之后**覆盖它自己那份",
        },
        "batch": {"name": data_dir.name, "data_dir": str(data_dir)},
        "config": {"trim": cfg.trim, "blank_dur": cfg.blank_dur, "win": cfg.win,
                   "step": cfg.step, "snap": bool(snap)},
        "counts": {"written": len(numn), "own": len(numn_new), "merged_kept": len(kept),
                   "blank_dur_source": blank_src,
                   "cycle_rows": int(len(dfc)),
                   "intervals_from": ("stage:" + str(args.from_stage)) if args.from_stage
                                     else "pipeline:load_batch"},
        "coverage": cov,
        "time_points": TIME_POINT_LEGEND,
        "conventions": CONVENTIONS,
        "columns": {"interval_selection.csv": [list(x) for x in COLUMN_LEGEND],
                    "interval_cycles.csv": [list(x) for x in CYCLE_COLUMN_LEGEND]},
        "segments": [{"seg": s, "meaning": SEG_LEGEND[s]} for s in SEGMENTS],
        "selfcheck": chk,
        "residuals": [
            "Isoclock 的 Timeinternal 是全局变量（最后一个载入的文件赢），本工具按逐文件"
            "周期反解 ⇒ 只有那个『最后一个文件』的周期与目标文件不同时下标才会漂；"
            "selfcheck.ti_sensitivity 已给出全批极值下的最坏漂移。",
            "秒数写的是**区间中点**（不是端点），以容忍半个下标的周期误差；"
            "--no-snap 可改回原样写 DRUID 的秒数。",
            "本工具只保证**同一批 cycle**，不保证同一套算法：Isoclock 会做它自己的"
            "±2σ 剔离群、普通铅扣除与分馏校正，两边结果不同是预期的。",
            "interval_cycles.csv 的 seg 只按**时间**分段（窗口端点全是秒），"
            "不替使用者做地球化学判断（Th/U、CL 才是做判断的东西）。",
            "逐 cycle 的 τ 用该 cycle 自己的时间戳；DRUID 窗口的 τ 用窗口中点 —— "
            "定义相同（都不是『与 100% 比』），取值点不同。",
            "blank_window() 是唯一的重写逻辑，已用 Tra.blank 逐通道对账（含退化支路）；"
            "intervals_from 是 stage:* 时没有 raw/blank 可对账 —— 该对账为 null"
            "（**没验**，不是通过）。",
            "intervals_from 是 stage:* 且原始 CSV 已不在时，时间轴只能按等间隔重建"
            "（实测重建偏差上限 2.0e-3 s ≈ 0.64% 个 cycle）⇒ 端点落进危险带的点下标"
            "可能差 1 个 cycle；selfcheck.time_axis 报出这一档的个数。原始 CSV 在时走精确路径。",
            "DRUID 跳过的玻璃标样 / 作废点不会出现在本文件里，而 Isoclock 的 "
            "dataprocess() 遇到这种缺键会 KeyError ⇒ 用 --merge-from 保底。",
        ],
    }
    (out / "interval_selection.json").write_text(
        json.dumps(jsonable(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (out / "README.txt").write_text(build_readme(), encoding="utf-8-sig")

    # ── 汇报 ──
    print(f"[4] 写出 {out}/Time_setting.npy（{len(numn)} 条）、"
          f"interval_selection.csv（{len(df)} 行 × {len(df.columns)} 列）、"
          + (f"interval_cycles.csv（{len(dfc)} 行 × {len(dfc.columns)} 列）、" if not dfc.empty
             else "")
          + "interval_selection.json、README.txt")
    if cov["missing"]:
        print(f"    ⚠ 目录里有 {len(cov['missing'])} 个 CSV 没被写到（DRUID 跳过了它们）："
              f"{cov['missing'][:6]}{' …' if len(cov['missing']) > 6 else ''}")
        print("      Isoclock 的 dataprocess() 遇到这些键会 KeyError ⇒ "
              "用 --merge-from <它自己生成的 Time_setting.npy> 保底。")
    if cov["extra_keys"]:
        print(f"    注：文件里有 {len(cov['extra_keys'])} 个键在数据目录里找不到对应 CSV"
              f"（多半是另一次批次的残留）：{cov['extra_keys'][:6]}")
    print(f"    下一步：把 {out.name}/Time_setting.npy 拷到 Isoclock 的 outputpath"
          f"（**加载步之后**）覆盖它自己那份。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
