"""
druid.cli.serve_ui —— 网页界面启动入口
=====================================

这是给**不用命令行的人**准备的入口：

    python -m druid.cli.serve_ui

回车之后会自动弹出浏览器，剩下的操作全在网页里点。
要退出就关掉那个黑窗口（或 Ctrl+C）。

为什么不动 PYTHONPATH 也能跑
----------------------------
必须以 `python -m druid.cli.serve_ui` 的形式调用，且工作目录在 druid 的**上级**。
包内的相对导入（`from .. import`）会自己找到兄弟模块，不需要 PYTHONPATH。
下面那段 try/except 是为了兼容 `python druid/cli/serve_ui.py` 这种直接跑脚本的写法。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# 与 reduce_batch 同一套判据：__package__ 为空 = 直接跑脚本，需自己补 sys.path。
# 不用 try/except ImportError，理由见 reduce_batch.py 的注释。
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from druid import __version__
from druid.webui.server import serve


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="druid.cli.serve_ui",
        description="启动本地网页界面（默认会自动打开浏览器）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--host", default="127.0.0.1",
                    help="监听地址。**保持 127.0.0.1**，改成 0.0.0.0 会把本机文件暴露到局域网")
    ap.add_argument("--port", type=int, default=0,
                    help="端口号，0 表示自动挑一个空闲端口（8765 起）")
    ap.add_argument("--no-browser", action="store_true",
                    help="不自动打开浏览器（已经手动开着标签页时用）")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    print(f"druid v{__version__}   正在启动网页界面…")
    try:
        # serve() 会一直阻塞到用户 Ctrl+C
        serve(host=args.host, port=args.port, open_browser=not args.no_browser)
    except KeyboardInterrupt:
        print("\n已退出。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
