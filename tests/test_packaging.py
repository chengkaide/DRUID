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


def test_webui_ui_params_match_backend_whitelist():
    """
    网页界面的参数入口与后端白名单必须对齐，两个方向各钉一条：

      ① 前端传的键，后端得认。传了白名单外的键会被**静默丢弃**（白名单正是
         为防脏参数而设），症状是"界面上改了、结果一点没变"。
      ② 后端放行的关键参数，界面上得真的有入口。`ref_preset` 就漏过一次：
         `handlers._TUNABLE` 里一直有它、CLI 也有 `--ref-preset`，但
         `static/` 里从来没有控件能把它传出去 —— 网页用户以为只有一套参考值，
         而它其实是五个档。

    只查"有没有这个入口"，不查控件长什么样 —— 后者是排版，会变。
    """
    app_js = (ROOT / "druid" / "webui" / "static" / "app.js").read_text(encoding="utf-8")
    m = re.search(r"const payload = \{(.*?)\n  \};", app_js, re.S)
    assert m, "app.js 里找不到 payload 字面量 —— 抓取正则要跟着改"
    keys = set(re.findall(r"(\w+)\s*:", m.group(1)))

    from druid.webui.handlers import _TUNABLE
    # 这两个键由路由单独处理，不进 BatchConfig，所以本来就不在白名单里
    route_own = {"data_dir", "out_excel"}
    unknown = keys - set(_TUNABLE) - route_own
    assert not unknown, f"前端传了后端不认的参数：{sorted(unknown)}"

    must_have_ui = {"bulk", "primary", "secondary", "ref_preset"}
    missing = must_have_ui - keys
    assert not missing, f"这些参数后端放行了、界面上却没有入口：{sorted(missing)}"


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


def test_shell_launchers_are_lf_and_executable():
    """
    类 Unix 的启动脚本必须是 **LF**，并且在 git 里记着 **100755**。

    `.command` 的第 1 行是 `#!/bin/bash`。若该行以 `\\r` 结尾（按 CRLF 存的），
    内核会去找一个叫 `/bin/bash\\r` 的解释器，报
    `bad interpreter: No such file or directory` —— 报错里看不出是行尾问题，
    双击的人只会说"没反应"。所以只能直接看磁盘字节。

    可执行位同理：Windows 上没有这个概念，新建/拷贝文件时最容易丢。
    丢了之后 Finder 双击报"无法执行"，而 `bash 文件名` 仍然正常 ——
    又是那种"在我机器上没问题"的坑。工作区上（Windows）看不出 mode，
    只能问 git 索引。

    ⚠ 这条与上面 `test_bat_files_are_crlf_on_disk` 互为镜像：那边要 CRLF、
    这边要 LF。两边的 `.gitattributes` 规则**只在 checkout 时生效**，
    谁用编辑器存一次反的，`git status` 依然报"干净"（入库时会规范化），
    错误就一直潜伏到有人真的去双击它。
    """
    import subprocess

    scripts = sorted(ROOT.glob("*.command")) + sorted(ROOT.glob("*.sh"))
    assert scripts, "一个 .command/.sh 都没有 —— macOS/Linux 的启动脚本不该消失"

    for p in scripts:
        data = p.read_bytes()
        crlf = data.count(b"\r\n")
        assert crlf == 0, f"{p.name} 里有 {crlf} 行是 CRLF —— shebang 会变成 bad interpreter"
        assert b"\r" not in data, f"{p.name} 里出现了裸 CR（0x0D）"
        assert data.startswith(b"#!"), f"{p.name} 第 1 行不是 shebang"

    # 可执行位。不是 git 仓库（例如从源码包解开的）就跳过 —— 同
    # test_no_identifying_strings_in_tracked_files 的处理。
    out = subprocess.run(["git", "ls-files", "-s", "--"]
                         + [p.name for p in scripts],
                         cwd=str(ROOT), capture_output=True, text=True,
                         encoding="utf-8").stdout
    if not out.strip():
        return
    for line in out.splitlines():
        mode, _, rest = line.partition(" ")
        name = rest.split("\t", 1)[-1]
        assert mode == "100755", (
            f"{name} 在 git 索引里是 {mode}，不是 100755 —— 丢了可执行位，"
            "Finder 双击会说'无法执行'")


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


def test_a_real_cjk_font_is_available():
    """
    至少有一个候选中文字体**真的能被 matplotlib 找到**，而且真的有汉字字形。

    ⚠ 这条守的是本仓最安静的一个故障：`CJK_FONTS` 里的字体一个都找不到时，
    matplotlib **不报错、不警告**，只是悄悄退回 DejaVu Sans —— 而 DejaVu 没有
    汉字，于是剖面图里的「年龄」「剥蚀时间」全变成「□□□」，批处理照样 exit 0、
    Excel 照样产出。在 Windows 上开发永远看不到（雅黑一定在），换到 macOS 才
    发现，那时已经出了几十张图。CI 的 macos job 就靠这条兜底。

    两个易错点：
      * `DejaVu Sans` **垫在 `CJK_FONTS` 末尾**（只为让 matplotlib 别抛警告），
        必须排除在断言之外 —— 否则这条测试恒真。
      * 只查名字不够：名字在、字形不在，照样画不出字。所以再拿 `ft2font`
        查一个汉字有没有 glyph index（0 = 没有）。
    """
    import platform as _platform

    from matplotlib import font_manager as fm

    from druid.core.constants import CJK_FONTS

    candidates = [f for f in CJK_FONTS if not f.startswith("DejaVu")]
    known = {f.name for f in fm.fontManager.ttflist}
    hit = [f for f in candidates if f in known]

    if _platform.system() == "Linux":
        # 中文字体在 Linux 上是**可选包**（fonts-noto-cjk / fonts-wqy-microhei），
        # 装不装取决于发行版与 CI 镜像 ⇒ 缺了是环境问题、不是代码写错，不在这里
        # 断言。Windows / macOS 的候选字体是**系统自带**，两边都强断言。
        if not hit:
            print("[提示] 本机没有任何候选中文字体；Linux 上可装 fonts-noto-cjk")
        return

    assert hit, (
        f"matplotlib 一个候选中文字体都找不到（候选 {candidates}，"
        f"本机共有 {len(known)} 个字体名）"
        " —— 图里的中文会静默变成「□□□」，而进程照样 exit 0")

    from matplotlib import ft2font
    name = hit[0]
    face = ft2font.FT2Font(fm.findfont(fm.FontProperties(family=name)))
    idx = face.get_char_index(ord("年"))
    assert idx, f"{name} 里没有汉字「年」的字形（glyph index = 0）—— 换一个候选字体"


def test_reveal_picks_the_command_for_each_platform():
    """
    `reveal()` / `list_drives()` 的 Windows / Darwin / Linux 三个分支各走对路。

    ⚠ 动机：开发机只有 Windows，**macOS 的代码路径平时一行都执行不到** ——
    而 `reveal()` 正是 mac 用户在界面上会真点的那个"打开所在文件夹"。
    跑 macos 的 CI job 也盖不住它（没有别的测试调用这个函数），所以只能
    在这里用假 `platform.system()` 把三个分支都走一遍。

    这也是本文件里**唯一**手写还原（而不用 pytest 的 monkeypatch）的地方 ——
    `_selftest.run()` 直接 `fn()` 调，不支持 fixture，见 `tests/_selftest.py`。
    别把它当惯例：能真跑的路径就不要用假平台。

    各平台预期（`reveal` 的语义是"打开**所在文件夹**"）：
      * Windows：目录 → `os.startfile(目录)`；文件 → `explorer /select, 文件`
      * Darwin ：一律 `open <目录>`
      * 其他   ：一律 `xdg-open <目录>`
    """
    import tempfile
    import types

    from druid.webui import handlers

    calls = []
    keep = (handlers.subprocess.Popen, handlers.platform,
            getattr(handlers.os, "startfile", None))
    handlers.subprocess.Popen = lambda cmd, *a, **k: calls.append(list(cmd))
    # os.startfile 只有 Windows 有；补个假的，好让另外两个平台也能走到这条分支
    handlers.os.startfile = lambda p: calls.append(["<startfile>", p])
    try:
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            d = base / "子目录"
            d.mkdir()
            f = base / "a.txt"
            f.write_text("x", encoding="utf-8")

            cases = (
                # 平台, 对目录的期望, 对文件的期望, list_drives 的期望
                ("Darwin", ["open", str(d)], ["open", str(base)], ["/"]),
                ("Linux", ["xdg-open", str(d)], ["xdg-open", str(base)], ["/"]),
                ("Windows", ["<startfile>", str(d)],
                 ["explorer", "/select,", str(f)], None),
            )
            for sysname, want_dir, want_file, want_drives in cases:
                # 默认参数把当前循环值捕获进来，避免闭包晚绑定
                handlers.platform = types.SimpleNamespace(
                    system=lambda s=sysname: s)

                assert handlers.reveal(str(d))["ok"], sysname
                assert calls[-1] == want_dir, (sysname, calls[-1], want_dir)

                assert handlers.reveal(str(f))["ok"], sysname
                assert calls[-1] == want_file, (sysname, calls[-1], want_file)

                got = handlers.list_drives()
                if want_drives is None:
                    # Windows 上返回本机真实存在的盘符；在别的系统上跑这条
                    # 分支时（Path("C:\\") 不存在）会是空表 —— 只查形状。
                    assert isinstance(got, list)
                    assert all(x.endswith(":\\") for x in got), got
                else:
                    assert got == want_drives, (sysname, got)
    finally:
        handlers.subprocess.Popen, handlers.platform = keep[0], keep[1]
        if keep[2] is None:
            del handlers.os.startfile
        else:
            handlers.os.startfile = keep[2]


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
    import re
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
    # 固定子串：本机用户名、项目地名、出现过的样品号与批次号前缀。
    patterns = ["云龙", "锡矿", "YL-17", "YL-46", "YL-48",
                "20220301", "20230501", "20181027", "20181028", "20181030",
                "凯凯"]
    # 光列固定前缀会一直漏 —— 换个批次就换个号。2026-10-02 真的漏过一次：
    # `tools/domain_span_check.py` 的 docstring 里写进了 `20230501ZDD`，
    # 而固定子串表里没有 `20230501`，本地与 CI 全绿。于是再补一条**样式**：
    # 批次号 = 8 位日期 + 2~3 个大写字母。
    # 示例批次 `EX2022A` 不匹配（`20` 后面必须是 6 位数字）。
    batch_like = re.compile(r"20[0-9]{6}[A-Z]{2,3}")
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
        for tok in sorted({m.group(0) for m in batch_like.finditer(text)}):
            where = "已入库" if rel in tracked else "未跟踪"
            hits.append(f"[{where}] {rel}: 含批次号样式 {tok}")
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
    否则"自检 N 项"这种话就只是广告词。（后来又涨到 98：补
    `calibration.primary_rsd` 时带了三条测试；同一天涨到 99：
    补「网页界面参数入口 vs 后端白名单」那条时 +1；同日再涨到 101：
    补「207/206 也要对自己的下限判触界」与「两个下限只能有一个出处」。
    再后来涨到 169：补四个核心链路测试文件 —— 之前 `reduce_interval`
    与整条深度剖面链路是**零测试覆盖**，见 `tests/test_reduction.py`、
    `tests/test_geochronology.py`、`tests/test_depth.py`、
    `tests/test_io_sequence.py`、`tests/test_core_models.py`。
    同轮再涨到 171：`segment` 补一条"下标口径"的**契约回归**（D-1），
    `sample_role` 补一条"作废标记整名相等"的用例。再涨到 172：加 macOS 启动器
    `启动数据处理工具.command` 时补一条「LF 行尾 + 100755 可执行位」的守护，
    与上面 `test_bat_files_are_crlf_on_disk` 互为镜像。174 是 CI 加 macOS
    job 时补「中文字体真的可用」与「reveal/list_drives 三平台分支」两条 ——
    macOS 的代码路径在这台 Windows 开发机上一行都执行不到，只能这么测。
    现为 175：A-22 补「窗口重叠 → 均值的 1σ 按 √K 放大」，只动 se、MSWD 不动。）
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
