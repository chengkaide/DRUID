"""
druid.console —— 让中文输出在非 UTF-8 控制台下不崩
==================================================

问题
----
这个包几乎所有的进度与摘要输出都是中文。Windows 上把输出**重定向到文件或管道**时
（`druid-reduce ... > log.txt`），Python 用 locale 编码（英文系统上是 cp1252）
去编码 stdout，于是第一句 `print("...通过")` 就抛：

    UnicodeEncodeError: 'charmap' codec can't encode characters in position 8-9

进程直接死掉，而且是在跑完批处理、准备打印结果的时候 —— 最难受的时机。
交互式控制台不一定会踩到（Python 会为控制台单独用宽字符 API），所以这个 bug
只在"重定向 + 非 UTF-8 locale"这个组合下出现，本地开发很难遇到。

CI 上就是这样挂的：ubuntu 全过，windows 在"不装 pytest 那条自检路"上失败，
因为 pytest 自己的报告器处理了编码，而 `tests/run_all.py` 直接 print。

做法
----
在入口把 stdout/stderr 切成 UTF-8，`errors="replace"` 兜底：
即使某个字符最终编不出来，也只丢那一个字符，不会中断整个流程。

为什么不给每个 print 包 try/except
---------------------------------
那样的噪声远大于收益，而且总会漏掉新加的 print。
输出流是"进程级"的东西，在入口统一处理一次才是它该待的地方。

放在包根目录而不是 `core/`
------------------------
`core/` 按约定是"纯计算、无副作用"，改输出流显然不是纯计算。
"""
from __future__ import annotations

import sys

__all__ = ["ensure_utf8_streams"]


def ensure_utf8_streams(*streams) -> None:
    """
    把给定的文本流重新配置为 UTF-8（默认 sys.stdout / sys.stderr）。

    幂等，可以随便多调几次；对不支持 `reconfigure` 的对象（比如 pytest 的
    捕获对象、`io.StringIO`）静默跳过 —— 那些流本来就按 Unicode 处理，
    不存在编码问题。
    """
    if not streams:
        streams = tuple(s for s in (getattr(sys, "stdout", None),
                                    getattr(sys, "stderr", None)) if s is not None)

    for s in streams:
        reconfigure = getattr(s, "reconfigure", None)
        if reconfigure is None:
            continue                      # StringIO / pytest 捕获对象等
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except (ValueError, OSError):
            # 流已经写过东西、或已被关闭 —— 都不是致命情况，跳过即可
            pass
