"""
打包自检 —— 版本号一致、入口点可调用、模块都能导入、网页资源没丢。

这四条都是"改完不报错、装完才发现"的问题，而且每一条都真的发生过：

  * **版本号写两处**（`pyproject.toml` 与 `druid/__init__.py`），改一处忘一处
    → 打出来的包自称 2.1.0，跑起来打印 2.0.0。网页界面那行 banner 就曾经
    硬编码着 "v2.0"，和 `__version__` 不一致。
  * **`readme = "README.md"` 指向不存在的文件** → `pip install -e .` 直接失败。
    这个仓库的 `pyproject.toml` 就这样躺了很久，直到真的去构建才发现。
  * **静态资源没进 wheel** → `pip` 装完打开网页是 404（就是 `package-data`
    漏声明的后果，而源码目录里跑永远不会暴露）。
  * **某个子模块导入就报错** → 只有真的 import 才知道。

为什么用正则读 pyproject 而不用 tomllib
---------------------------------------
`tomllib` 是 Python 3.11 才进标准库，而这个包支持 3.9（CI 就在跑 3.9）。
为了一个测试去拉 `tomli` 依赖不值当，这里的取值都很规整，正则够用。
"""
from __future__ import annotations

import importlib
import pkgutil
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import _selftest                                                # noqa: E402
import druid                                                    # noqa: E402

PYPROJECT = ROOT / "pyproject.toml"


def _section(name: str) -> str:
    """取 pyproject.toml 里某个 [section] 的正文（到下一个 [ 为止）。"""
    text = PYPROJECT.read_text(encoding="utf-8")
    m = re.search(rf"^\[{re.escape(name)}\]\s*$(.*?)(?=^\[|\Z)", text,
                  re.M | re.S)
    assert m, f"pyproject.toml 里找不到 [{name}]"
    return m.group(1)


# ═════════════════════════════════════════════════════════════════════════════
# 一、版本号
# ═════════════════════════════════════════════════════════════════════════════
def test_version_matches_pyproject():
    m = re.search(r'^version\s*=\s*"([^"]+)"', _section("project"), re.M)
    assert m, "pyproject.toml 的 [project] 里没有 version"
    assert m.group(1) == druid.__version__, (
        f"pyproject 是 {m.group(1)}，druid.__version__ 是 {druid.__version__}")


def test_version_looks_like_semver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", druid.__version__), druid.__version__


def test_installed_metadata_agrees_if_present():
    """装过（pip install -e .）时，元数据里的版本也必须一致。"""
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:                                          # py<3.8 兜底
        return
    try:
        installed = version("druid")
    except PackageNotFoundError:
        return                       # 直接跑源码、没装过 —— 不是失败
    assert installed == druid.__version__, (installed, druid.__version__)


# ═════════════════════════════════════════════════════════════════════════════
# 二、pyproject 自己说的东西必须存在
# ═════════════════════════════════════════════════════════════════════════════
def test_declared_readme_exists():
    m = re.search(r'^readme\s*=\s*"([^"]+)"', _section("project"), re.M)
    if not m:
        return                        # 没声明就不用查
    p = ROOT / m.group(1)
    assert p.is_file(), f"pyproject 声明 readme = {m.group(1)}，但文件不存在"


def test_declared_console_scripts_are_importable():
    """[project.scripts] 里每个 "模块:函数" 都要真的能导入、真的可调用。"""
    body = _section("project.scripts")
    pairs = re.findall(r'^(\S+)\s*=\s*"([\w.]+):(\w+)"', body, re.M)
    assert pairs, "没解析到任何入口点 —— 正则或 pyproject 格式变了"

    for cmd, mod_name, attr in pairs:
        assert cmd.startswith("druid-"), f"入口点命名不符合约定: {cmd}"
        mod = importlib.import_module(mod_name)
        fn = getattr(mod, attr, None)
        assert fn is not None, f"{mod_name} 里没有 {attr}"
        assert callable(fn), f"{mod_name}:{attr} 不可调用"


def test_webui_static_assets_exist_and_are_declared():
    """
    三个静态资源必须同时满足两件事：磁盘上有、且 `package-data` 声明了。
    少第二件时源码目录里跑一切正常，`pip install` 完打开网页就是 404。
    """
    static = ROOT / "druid" / "webui" / "static"
    declared = _section("tool.setuptools.package-data")

    for name in ("index.html", "app.js", "style.css"):
        assert (static / name).is_file(), f"缺静态资源 {name}"
        assert (static / name).stat().st_size > 0, f"{name} 是空文件"
        assert f"webui/static/*{Path(name).suffix}" in declared, (
            f"package-data 没声明 webui/static/*{Path(name).suffix}")


def test_bat_files_are_crlf_on_disk():
    """
    批处理必须是 CRLF。cmd.exe 按行读，LF 结尾的 .bat 会出现"命令明明写在
    文件里、运行时却读不到"这种极难排查的故障 —— 而双击那个 .bat 的人
    （实验室里日常用它启动工具）没有任何办法自己排查。

    ⚠ 为什么这件事值得专门测：
    `.gitattributes` 里的 `*.bat text eol=crlf` **只在 checkout 时生效**。
    如果谁用编辑器或某个脚本按 LF 存了一次，`git status` 依然报"干净"
    —— 因为入库时会规范化，内容被认为没变。于是这个错误可以一直潜伏到
    有人真的去双击它。git 自己看不见，只能直接看磁盘上的字节。
    """
    bats = sorted(ROOT.glob("*.bat")) + sorted(ROOT.glob("*.cmd"))
    assert bats, "一个 .bat/.cmd 都没有 —— 启动脚本不该消失"

    for p in bats:
        data = p.read_bytes()
        assert b"\r\n" in data, f"{p.name} 没有任何 CRLF 行尾"
        lone_lf = data.replace(b"\r\n", b"").count(b"\n")
        assert lone_lf == 0, f"{p.name} 里有 {lone_lf} 行只有 LF —— cmd.exe 可能读不到"


def test_bat_files_stay_ascii():
    """
    cmd.exe 按当前代码页解析批处理，中文注释会变乱码。而解释器路径里
    恰好含中文用户名，所以中文路径必须靠运行时变量（%USERPROFILE%、%~dp0）抵达，
    不能写死在文件里。
    """
    for p in sorted(ROOT.glob("*.bat")) + sorted(ROOT.glob("*.cmd")):
        data = p.read_bytes()
        bad = [(i, b) for i, b in enumerate(data) if b > 127]
        assert not bad, (
            f"{p.name} 第 {bad[0][0]} 字节起出现非 ASCII（{bytes(x for _, x in bad[:8])!r}）"
            " —— 中文注释在 cmd.exe 里会变成乱码")


# ═════════════════════════════════════════════════════════════════════════════
# 三、每个子模块都导得进来
# ═════════════════════════════════════════════════════════════════════════════
def test_every_submodule_imports():
    """
    语法错误、循环导入、写错名字的 from-import —— 只有真的 import 才发现。
    pyflakes 不做导入解析的运行时验证，`python -m druid.cli.serve_ui --help`
    也只会碰到其中一个分支，所以这里遍历全部。
    """
    mods = sorted(m.name for m in
                  pkgutil.walk_packages(druid.__path__, prefix="druid."))
    assert len(mods) >= 25, f"只发现 {len(mods)} 个子模块，是不是少了？{mods}"

    broken = []
    for name in mods:
        try:
            importlib.import_module(name)
        except Exception as e:                                   # noqa: BLE001
            broken.append(f"{name}: {type(e).__name__}: {e}")
    assert not broken, broken


def test_public_api_is_where_it_is_documented():
    """
    `druid/__init__.py` 的文档里点名了这几个地方是"对外用法"，
    所以它们的公开名字不能悄悄搬走 —— 搬走等于破坏使用说明。
    """
    from druid.cli import reduce_batch, serve_ui
    from druid.workflow import BatchConfig, run_batch

    assert callable(run_batch) and callable(BatchConfig)
    assert callable(reduce_batch.main) and callable(serve_ui.main)


# 文本类文件（二进制的不去读）
TEXT_SUFFIX = {".md", ".txt", ".py", ".js", ".toml", ".html", ".yml", ".yaml",
               ".css", ".spec", ".bat", ".csv"}


def test_no_identifying_strings_in_tracked_files():
    """
    仓库里不许再出现真名 —— 地名、样品号、批次号、本机用户名。

    这不是洁癖。这些信息曾经散布在 **13 个文件**里：`pyproject.toml` 的作者名、
    代码注释与文档里的样例路径、以及文档里内嵌那张图的标题（是像素，不是文本）。
    发出去就等于把"这套工具"和"一份未发表的数据"绑在一起，而且收不回来。

    现在示例批次用 `S01`…`S48` / `EX2022A`，文档里的例子也一并改名 ——
    这条测试是防止回流：以后谁顺手贴一条真实路径进来，CI 立刻拦下。

    ⚠ 本文件自己必须把那些模式写出来，所以把自己排除在外。
    """
    import subprocess

    files = subprocess.run(["git", "-c", "core.quotePath=false", "ls-files"],
                           cwd=str(ROOT), capture_output=True, text=True,
                           encoding="utf-8").stdout
    if not files.strip():
        return                    # 不是 git 仓库（例如从源码包解开的），跳过

    me = Path(__file__).name
    patterns = ["云龙", "锡矿", "YL-46", "20220301", "凯凯"]
    hits = []
    for rel in files.splitlines():
        if not rel.strip() or Path(rel).name == me:
            continue
        p = ROOT / rel
        if not p.exists() or p.suffix.lower() not in TEXT_SUFFIX:
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for k in patterns:
            if k in text:
                hits.append(f"{rel}: 含 {k}")
    assert not hits, "仓库里出现了可识别信息，请改成中性示例：" + "; ".join(hits)


if __name__ == "__main__":
    raise SystemExit(_selftest.run(globals()))
