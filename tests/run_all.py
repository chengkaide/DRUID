"""
跑 tests/ 下全部自检 —— 不装 pytest 也能用：

    python tests/run_all.py

装过开发依赖的话，等价的命令是 `python -m pytest -q`（CI 用的是后者）。
两条路都必须通：实验室那台机器上只有 numpy/pandas/matplotlib，
而 CI 上想要更好的失败报告。
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import _selftest                                                # noqa: E402


def main() -> int:
    # 本文件与 _selftest 的输出含中文；Windows 上重定向时会踩编码坑，
    # 见 druid/console.py 的说明。CI 的 windows job 就是这么挂的。
    from druid.console import ensure_utf8_streams
    ensure_utf8_streams()

    modules = sorted(p.stem for p in HERE.glob("test_*.py"))
    if not modules:
        print("没找到任何 test_*.py")
        return 1

    rc = 0
    for name in modules:
        print(f"\n=== {name} ===")
        rc |= _selftest.run(vars(importlib.import_module(name)))
    print("\n全部文件跑完" if rc == 0 else "\n有失败项")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
