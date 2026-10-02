#!/bin/bash
# =====================================================================
#  DRUID -- LA-ICP-MS zircon U-Pb reduction tool   (macOS / Linux launcher)
#
#  这是 Windows 的 `启动数据处理工具.bat` 在 macOS 上的对应物，做同样三件事：
#    ① 找解释器  ② cd 到包根  ③ python -m druid.cli.serve_ui
#
#  用法：
#    在 Finder 里双击本文件（或用 `bash 启动数据处理工具.command`）。
#    终端窗口别关；关掉即停止服务。
#
#  ⚠ 首次双击若被 Gatekeeper 拦下（"来自身份不明的开发者"），二选一：
#      · 右键 → 打开 → 再点一次"打开"；或
#      · 先在终端执行：xattr -d com.apple.quarantine 启动数据处理工具.command
#
#  ⚠ 行尾必须是 LF。CRLF 会让第 1 行的 shebang 带上一个 \r，
#    症状是 macOS 直接报 "bad interpreter"。已在 .gitattributes 里钉死。
# =====================================================================

set -u
cd "$(dirname "$0")" || exit 1

PY=""

# ① WorkBuddy 托管的隔离环境（本机开发用）。目录名有 druid / upb 两种历史命名。
for env in druid upb; do
    cand="$HOME/.workbuddy/binaries/python/envs/$env/bin/python"
    if [ -x "$cand" ]; then PY="$cand"; break; fi
done

# ② 其它账户下的同名环境 —— 覆盖"文件被拷到另一台 Mac、账户名不同"的情形。
if [ -z "$PY" ]; then
    for cand in /Users/*/.workbuddy/binaries/python/envs/druid/bin/python \
                /Users/*/.workbuddy/binaries/python/envs/upb/bin/python; do
        if [ -x "$cand" ]; then PY="$cand"; break; fi
    done
fi

# ③ 系统 python3（Homebrew / python.org 装的都算）。
if [ -z "$PY" ]; then
    PY="$(command -v python3 || true)"
fi

if [ -z "$PY" ]; then
    echo
    echo "  [错误] 找不到 Python 解释器。"
    echo "        装了 Homebrew 的话：brew install python3"
    echo
    read -r -p "  按回车退出…" _
    exit 1
fi

# 依赖自检。macOS 上最可能踩的就是"解释器找到了、包没装"，
# 与其让它跑到一半抛 ImportError，不如在门口说清楚。
if ! "$PY" -c "import numpy, pandas, matplotlib, openpyxl" >/dev/null 2>&1; then
    echo
    echo "  [错误] 这个解释器里缺依赖："
    echo "        $PY"
    echo
    echo "  在包根目录（含 druid/ 的那一层）执行："
    echo "        \"$PY\" -m pip install -e ."
    echo
    read -r -p "  按回车退出…" _
    exit 1
fi

echo "  Python: $PY"
echo
echo "  正在启动本地服务…下面会打印一个 http://127.0.0.1:端口 的地址，"
echo "  用浏览器打开它即可。这个终端窗口保持开着，关掉就停止服务。"
echo

# 与 Windows 版一致地用 --no-browser：把地址打印出来让人自己点，
# 免得 webbrowser 模块挑到意外的浏览器。想自动打开就把这个参数去掉。
exec "$PY" -m druid.cli.serve_ui --port 0 --no-browser
