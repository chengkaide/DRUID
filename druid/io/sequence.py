"""
druid.io.sequence —— 测点序列（LIST.xls / LIST.xlsx）读取与角色判定
=================================================================

序列文件长什么样
----------------
两列、若干行，每行一个测点：

    序号0  20220301CKDB_1    SRM 612
    序号1  20220301CKDB_2    91500        ← Excel 里写成数字 91500.0
    序号2  20220301CKDB_5    Ple
    序号3  20220301CKDB_7    YL-46-1
    …

⚠ 一个陷阱：**91500 在 xls 里存成了浮点数 91500.0**。
   如果直接 str(91500.0) 会得到 "91500.0"，后面的判等就失效了。
   所以下面专门把它还原成 "91500"。
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core.constants import (
    ROLE_GLASS,
    ROLE_LABEL_CN,
    ROLE_PRIMARY,
    ROLE_SECONDARY,
    ROLE_UNKNOWN,
)


# ─────────────────────────────────────────────────────────────────────────────
# 一、样品名 → 角色
# ─────────────────────────────────────────────────────────────────────────────
def sample_role(name: str, primary: str = "91500", secondary: str = "Ple") -> str:
    """
    按样品名判定这个测点在流程里扮演什么角色。

    返回值为 core.constants 里的四种之一：
        primary_std   主标     —— 用来算归一化因子 F（本项目 = 91500）
        secondary_std 监控标样  —— 不参与归一化，只做 QC（本项目 = Plešovice）
        glass         玻璃标样  —— 只测微量元素，U-Pb 完全不参与
        unknown       未知样品  —— 真正要定年的对象

    为什么玻璃标样要单独挑出来
    --------------------------
    NIST SRM 612/610 是硅酸盐玻璃，U 含量人为加标到 ~40 ppm、
    且 Pb 同位素组成是现代普通铅。给它"定年"会得到一个荒谬的结果，
    如果混进 normalization 会直接污染整批数据。必须在入口就拦截。

    primary / secondary 参数
    -------------------------
    允许用户把自己实际使用的标样名告诉程序。例如某批次用 GJ1 当主标、
    用 Temora 当监控标样，只需把 primary="GJ1" / secondary="Temora" 传进来，
    这些测点就会被正确识别，而不是被误判成"未知样品"导致整批失去主标、
    归一化因子曲线为空、最后 np.interp 在空数组上崩溃。
    默认仍是 91500 / Ple，向后兼容。
    """
    s = str(name).strip().lower()
    p = str(primary).strip().lower()
    sec = str(secondary).strip().lower()

    # ① 优先匹配用户显式指定的主标 / 监控标样名（大小写不敏感）
    if p and s == p:
        return ROLE_PRIMARY
    if sec and s == sec:
        return ROLE_SECONDARY
    # ② 默认别名（不指定参数时的兜底）
    if s in ("91500", "91500 zircon", "91500z"):
        return ROLE_PRIMARY
    if s in ("ple", "plesovice", "pl", "plešovice"):
        return ROLE_SECONDARY
    # 玻璃：命中任意一个关键词就算
    if ("srm" in s) or ("nist" in s) or ("612" in s) or ("610" in s):
        return ROLE_GLASS
    return ROLE_UNKNOWN


def normalize_name(value) -> str:
    """
    把序列文件里读到的第二列（样品名）规范化成字符串。

    处理三种 Excel 脏数据情形：
        1) 数字 91500.0        → "91500"（不然字符串判等会失败）
        2) 带首尾空格 " Ple "  → "Ple"
        3) 空值/None           → ""（后续会被当作无效行跳过）
    """
    if value is None:
        return ""
    # 情形 1：浮点整数 → 去掉小数点
    if isinstance(value, float) and abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    # 其余情况统一转字符串并去空白
    return str(value).strip()


# ─────────────────────────────────────────────────────────────────────────────
# 二、读序列文件
# ─────────────────────────────────────────────────────────────────────────────
def _read_rows(path: Path):
    """
    底层读取：按扩展名选择解析器，产出 [(文件名, 样品名), ...]。

    · .xls  → xlrd（老 Excel 二进制格式，pandas 已弃用 xlrd 读 xls 需显式安装）
    · .xlsx → pandas.read_excel（openpyxl 引擎）
    · 其余  → 按 CSV 兜底，容忍度更高
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".xls":
        # xlrd 2.x 只支持 .xls，不支持 .xlsx
        import xlrd
        book = xlrd.open_workbook(str(path))
        sheet = book.sheet_by_index(0)
        out = []
        for r in range(sheet.nrows):
            file_cell = normalize_name(sheet.cell_value(r, 0))
            name_cell = normalize_name(sheet.cell_value(r, 1))
            # 文件名为空的行通常是页眉/空行，直接丢弃
            if file_cell:
                out.append((file_cell, name_cell))
        return out

    if suffix in (".xlsx", ".xlsm"):
        df = pd.read_excel(path, header=None, usecols=[0, 1], engine="openpyxl")
        return [(normalize_name(a), normalize_name(b))
                for a, b in zip(df.iloc[:, 0], df.iloc[:, 1])
                if normalize_name(a)]

    # 兜底：当 CSV 处理（sep=None 让 pandas 自动嗅探分隔符）
    df = pd.read_csv(path, header=None, sep=None, engine="python")
    return [(normalize_name(a), normalize_name(b))
            for a, b in zip(df.iloc[:, 0], df.iloc[:, 1]) if normalize_name(a)]


def read_sequence(path, primary=None, secondary=None) -> pd.DataFrame:
    """
    读 *_LIST.xls(x)，返回标准化后的测点序列表。

    参数
    ----
    primary / secondary : 可选，用户实际使用的主标 / 监控标样名。
        传入后会覆盖 sample_role 里的默认 91500 / Ple 判定（见 sample_role）。

    返回列
    ------
        file   : 不带扩展名的文件名（如 "20220301CKDB_7"，拼 .csv 即得数据文件路径）
        sample : 规范化后的样品名（如 "91500" / "Ple" / "YL-46-1"）
        role   : 角色（见 sample_role）
        order  : 在序列中的原始次序，从 1 开始计数（= Excel 里的序号列）
    """
    p = str(primary) if primary else "91500"
    sec = str(secondary) if secondary else "Ple"
    rows = []
    for i, (fname, sname) in enumerate(_read_rows(Path(path))):
        rows.append(dict(
            # rstrip(".csv")：万一 LIST 里写了带扩展名的文件名，这里统一去掉
            file=str(fname).rsplit(".", 1)[0] if str(fname).lower().endswith(".csv") else str(fname),
            sample=sname,
            role=sample_role(sname, p, sec),
            order=i + 1,
        ))
    return pd.DataFrame(rows, columns=["file", "sample", "role", "order"])


def spot_csv_path(data_dir, fname) -> Path:
    """
    由「批次目录 + 序列登记名」拼出实际 CSV 路径。

    约定（重要）
    ----------
    LIST 里登记的 `file` 列**不带扩展名**（如 "20220301CKDB_7"），
    磁盘上的数据文件是同名的 `.csv`。

    为什么要把这一个字符串拼接单独拎成函数
    ----------------------------------
    这个"补 .csv"的约定如果到处手写，早晚会有人在某个调用点忘记，
    然后在不同地方报出互相矛盾的错误。集中到这里之后：
        · workflow.load_batch 用它；
        · webui.handlers 的完整性检查用它；
    两边永远一致 —— 网页提示"文件齐全"的话，跑起来就一定读得到。
    """
    name = str(fname)
    # 容错：万一传进来的已经带了后缀，不重复拼
    if not name.lower().endswith(".csv"):
        name += ".csv"
    return Path(data_dir) / name


def sequence_summary(seq: pd.DataFrame) -> str:
    """把序列里各角色的点数拼成一行可读文字，用于命令行进度打印。"""
    if seq.empty:
        return "（空序列）"
    parts = [f"{ROLE_LABEL_CN.get(r, r)} {int(n)} 点"
             for r, n in seq["role"].value_counts().items()]
    return " / ".join(parts)
