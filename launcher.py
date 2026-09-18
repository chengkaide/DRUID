"""
打包用的入口脚本（PyInstaller 从它开始分析依赖）。

为什么不能直接用 `druid/cli/serve_ui.py` 当入口
----------------------------------------------
那里面有一句 `try: from .. import / except ImportError: sys.path.insert(...)`，
是给 `python druid/cli/serve_ui.py` 这种直接跑脚本的写法留的兜底。
PyInstaller 冻结后走的是另一套导入器，相对的 `..` 有时解析不出来。

入口脚本就一行直接调用，不做任何路径魔法 —— 打包时最省事。
"""
from druid.cli.serve_ui import main

if __name__ == "__main__":
    raise SystemExit(main())
