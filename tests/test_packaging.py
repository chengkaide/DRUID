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


def test_webui_server_header_carries_the_live_version():
    """
    HTTP 响应里的 `Server` 头必须来自 `druid.__version__`：BaseHTTPRequestHandler
    把它拼在 `server_version` 后面，所以那一行等于对外的自我声明。

    曾经硬编码成 "druid-webui/2.0" —— 包已经是 2.2.2，浏览器/curl 里还自称 2.0。
    这正是本文件开头说的那个毛病（版本号写两处、改一处忘一处）**第二次复发**：
    第一次是网页界面那行 banner，这次换成了 HTTP 头。所以再钉一遍。
    """
    from druid.webui.server import Handler

    assert Handler.server_version == f"druid-webui/{druid.__version__}", (
        f"Server 头自称 {Handler.server_version}，而包版本是 {druid.__version__}")


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

    ⚠ 扫描范围包含**尚未 git add 的新文件**（`--others --exclude-standard`）。
    只用 `ls-files` 会留一个真实的盲区：新写的文档在提交前是"未跟踪"状态，
    本地跑这条测试会通过、CI 上才失败 —— 这个盲区真的发生过一次：一份新文档
    的说明文字里混进了真实样品号，本地全绿，推送后 CI 立刻拦下。
    """
    import subprocess

    def _git(*args) -> set:
        out = subprocess.run(
            ["git", "-c", "core.quotePath=false", *args],
            cwd=str(ROOT), capture_output=True, text=True,
            encoding="utf-8").stdout
        return {ln for ln in out.splitlines() if ln.strip()}

    # 已入库的 + 尚未 add 但也没被忽略的 —— 后者是"`git add .` 会顺手带走"的东西，
    # 只看 --cached 会留一个真实的盲区（见 docstring）。
    files = _git("ls-files", "--cached", "--others", "--exclude-standard")
    if not files:
        return                    # 不是 git 仓库（例如从源码包解开的），跳过
    tracked = _git("ls-files", "--cached")

    me = Path(__file__).name
    patterns = ["云龙", "锡矿", "YL-46", "20220301", "凯凯"]
    hits = []
    for rel in sorted(files):
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
                where = "已入库" if rel in tracked else "未跟踪"
                hits.append(f"[{where}] {rel}: 含 {k}")
    if not hits:
        return

    n_tracked = sum(1 for h in hits if h.startswith("[已入库]"))
    advice = (
        "这些字会被推上公开仓库，而且收不回来。两种情形，修法不同：\n"
        "  · 这份文件**要进仓库** → 把真名换成中性示例"
        "（地名→`EX2022A`、样品→`S01`…、本机用户名→`<用户名>`）。\n"
        "  · 这份文件**只是本地笔记、本就不该进仓库** → 加进本地排除文件：\n"
        "        echo '<文件名>' >> .git/info/exclude\n"
        "    它与 `.gitignore` 的区别：**.gitignore 入版本库、全队共享；\n"
        "    .git/info/exclude 只对这台机器生效、不进版本库**。"
        "本地笔记要用后者。\n"
        "\n"
        "⚠ 不要改这条测试。它刻意把**未跟踪文件**也算进来，是花过代价的：\n"
        "曾经一份新写的文档在提交前混进了真实样品号，本地跑全绿、推上去 CI 才红。"
        if n_tracked else
        "这些文件还没 `git add`，但 `git add .` 会顺手把它们带走，所以现在就得处理。\n"
        "若它们只是**本地笔记、本就不该进仓库** → 加进本地排除文件：\n"
        "        echo '<文件名>' >> .git/info/exclude\n"
        "（`.gitignore` 会入版本库、全队共享；`.git/info/exclude` 只对这台机器生效。\n"
        " 本地笔记要用后者 —— 换一台机器 clone 下来本就没有这些文件。）\n"
        "\n"
        "⚠ 不要改这条测试。它刻意把未跟踪文件也算进来，是花过代价的：\n"
        "曾经一份新写的文档在提交前混进了真实样品号，本地跑全绿、推上去 CI 才红。"
    )
    assert not hits, ("仓库里出现了可识别信息：\n  " + "\n  ".join(hits)
                      + "\n\n" + advice)


def test_ensure_utf8_streams_survives_a_cp1252_console():
    """
    Windows 上把输出**重定向到文件/管道**时，Python 用 locale 编码
    （英文系统上是 cp1252）编码 stdout，于是第一句中文 print 就
    `UnicodeEncodeError: 'charmap' codec can't encode characters` 把进程带走。

    这不是假想问题：CI 的 windows job 就是这么挂的 —— ubuntu 全过，
    windows 在"不装 pytest 那条自检路"上失败，因为 pytest 自己的报告器
    处理了编码，而 `tests/run_all.py` 直接 print。
    而且它发生在**跑完批处理、准备打印结果**的时候，最难受的时机。

    `druid.console.ensure_utf8_streams()` 在入口把流切成 UTF-8 解决它。
    """
    import io

    from druid.console import ensure_utf8_streams

    # errors 默认是 strict，所以修复前这一句必然抛
    buf = io.TextIOWrapper(io.BytesIO(), encoding="cp1252")
    ensure_utf8_streams(buf)
    buf.write("通过 ✓")                       # 修复前：UnicodeEncodeError
    buf.flush()
    assert buf.encoding.lower().replace("-", "") == "utf8", buf.encoding

    # 幂等
    ensure_utf8_streams(buf)

    # 对没有 reconfigure 的对象要静默跳过（io.StringIO、pytest 的捕获对象）
    ensure_utf8_streams(io.StringIO())
    ensure_utf8_streams()                     # 默认参数 = 真实的 stdout/stderr



# ═════════════════════════════════════════════════════════════════════════════
# 五、文档里写的数字
# ═════════════════════════════════════════════════════════════════════════════
def test_landing_page_states_the_real_selfcheck_count():
    """落地页写着"自检 N 项"，N 必须等于 tests/ 下**真实**的自检项数。

    `_selftest.run()` 收集的是各文件顶层的 `test_*` 函数，所以这里也按 AST 数
    顶层 `def test_`，不 import —— import 一整套测试既慢又可能有副作用。

    这个数在 2.6.0 从 80 涨到 95（其中一条就是本测试），页面上当时还写着 80，
    挂了整整一版才被发现。加测试项的时候顺手把页面上的数字改掉 ——
    否则"自检 95 项"这种话就只是广告词。
    """
    import ast
    n = 0
    for f in sorted((ROOT / "tests").glob("test_*.py")):
        for node in ast.parse(f.read_text(encoding="utf-8")).body:
            if isinstance(node, ast.FunctionDef) and node.name.startswith("test_"):
                n += 1

    doc = (ROOT / "docs" / "index.html").read_text(encoding="utf-8")
    hits = re.findall(r"自检\s*(?:<i>)?(\d+)", doc)
    assert hits, "落地页里找不到「自检 N 项」—— 这条守护失去对象了，请更新它"
    for h in hits:
        assert int(h) == n, f"落地页写「自检 {h} 项」，实际是 {n} 项"

    m = re.search(r'<div class="n">(\d+)</div><div class="d">项自检', doc)
    assert m, "落地页的「基本情况」里没有自检项数，这条守护失去对象了"
    assert int(m.group(1)) == n, f"落地页的「基本情况」写 {m.group(1)} 项，实际是 {n} 项"
if __name__ == "__main__":
    raise SystemExit(_selftest.run(globals()))
