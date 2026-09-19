"""
检查构建出来的 wheel 里该有的东西都在。

    python -m build --wheel && python tests/check_wheel.py
    python tests/check_wheel.py G:/某个目录        # 指定 wheel 所在目录

为什么要单独一个脚本，而不是写在 CI 的 YAML 里
----------------------------------------------
YAML 里嵌一段 Python 的缩进规则很别扭，出错了也没法在本地复现。
放在这里就能在推之前先跑一遍 —— 而它抓的正是"源码目录里跑得好好的、
装完却坏掉"的那类问题，这种问题最难在本地发现：

  * 静态资源没进 wheel → `pip install` 完打开网页是 404
    （`package-data` 少一行就会这样，而源码目录里永远不会暴露）；
  * 入口点没写进去 → 装完 `druid-reduce` 命令不存在；
  * 模块被 `packages.find` 漏掉 → 开发时能 import，装完 ImportError。
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# 网页界面的三个资源：缺任何一个，界面都起不来
REQUIRED_STATIC = [
    "druid/webui/static/index.html",
    "druid/webui/static/app.js",
    "druid/webui/static/style.css",
]
REQUIRED_SCRIPTS = ["druid-reduce", "druid-gui"]
MIN_MODULES = 25                      # 当前 29 个，留一点余量


def check(wheel: Path) -> list:
    problems = []
    with zipfile.ZipFile(wheel) as z:
        names = z.namelist()

        for need in REQUIRED_STATIC:
            if need not in names:
                problems.append(f"缺静态资源 {need} —— pip 装完网页会 404")

        eps = [n for n in names if n.endswith("entry_points.txt")]
        if len(eps) != 1:
            problems.append(f"entry_points.txt 有 {len(eps)} 个，应为 1")
        else:
            text = z.read(eps[0]).decode("utf-8")
            for cmd in REQUIRED_SCRIPTS:
                if cmd not in text:
                    problems.append(f"缺入口点 {cmd}")

        mods = [n for n in names if n.endswith(".py")]
        if len(mods) < MIN_MODULES:
            problems.append(f"只打进 {len(mods)} 个模块（应 ≥ {MIN_MODULES}）")

        print(f"  {wheel.name}")
        print(f"    {len(names)} 个文件 / {len(mods)} 个模块")
        print(f"    静态资源 {len(REQUIRED_STATIC)}/3"
              f"，入口点 {'有' if len(eps) == 1 else '?'}")

    return problems


def main(argv) -> int:
    if len(argv) > 1:
        target = Path(argv[1])
    else:
        target = ROOT / "dist"

    wheels = sorted(target.glob("*.whl"))
    if not wheels:
        print(f"{target} 下没有 .whl —— 先跑 python -m build --wheel")
        return 1

    print("检查 wheel：")
    problems = []
    for w in wheels:
        problems += check(w)

    if problems:
        print("\n不合格：")
        for p in problems:
            print(f"  ✗ {p}")
        return 1
    print("\nwheel 完整 ✓")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
