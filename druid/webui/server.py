"""
druid.webui.server —— 零依赖 HTTP 服务（只用 Python 标准库）
===========================================================

为什么不装 Flask / FastAPI
---------------------------
这是一个**给自己用的本地小工具**，引一个 Web 框架意味着：
    · 多一份依赖要装、要维护、要在别人机器上复现；
    · Python 自带的 http.server 足够跑"一个页面 + 几个接口"；
    · 依赖越少，将来打包成单文件 exe 时体积越小、坑越少。
所以这里只用 http.server + ThreadingHTTPServer。取舍是够用的那一边。

线程模型
--------
    ThreadingHTTPServer：每个请求一个线程        ← 收发 JSON，毫秒级返回
    TaskManager：批处理跑在独立的第 N+1 个线程    ← 分钟级计算
前端通过轮询 /api/status 拿进度，而不是长连接推送 —— 本地场景够用且好调试。

安全边界
--------
只监听 127.0.0.1（本机回环），外部机器连不上。这不是网关服务，
不要 `--host 0.0.0.0` 暴露到局域网。
"""
from __future__ import annotations

import json
import mimetypes
import sys
import threading
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Optional

from .. import __version__
from . import handlers
from .tasks import MANAGER

# 静态资源目录：与 server.py 同级下的 static/
STATIC_DIR = Path(__file__).resolve().parent / "static"

# 常见扩展名 → MIME。mimetypes 在精简版 Python 上偶尔认不全，这里兜一下。
_MIME = {
    ".html": "text/html; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".ico": "image/x-icon",
}


# ═════════════════════════════════════════════════════════════════════════════
# 请求处理器
# ═════════════════════════════════════════════════════════════════════════════
class Handler(BaseHTTPRequestHandler):
    """所有 HTTP 请求的入口。GET 走 do_GET，POST 走 do_POST。"""

    server_version = "druid-webui/2.0"

    # ── 基础工具 ──────────────────────────────────────────────────────────
    def _send(self, code: int, body: bytes, ctype: str, extra: Optional[dict] = None):
        """统一的应答出口，避免每个分支重复写一堆 set_header。"""
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")   # 调试时不想被浏览器缓存坑
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, data: Any, code: int = 200):
        body = json.dumps(data, ensure_ascii=False, default=str).encode("utf-8")
        return self._send(code, body, "application/json; charset=utf-8")

    def _error(self, msg: str, code: int = 400):
        return self._json(dict(ok=False, error=msg), code)

    def _text(self, s: str, code: int = 200, ctype="text/plain; charset=utf-8"):
        return self._send(code, s.encode("utf-8"), ctype)

    def _read_json(self) -> dict:
        """读 POST 的 JSON body。Content-Length 缺失时按 0 处理。"""
        n = int(self.headers.get("Content-Length") or 0)
        if n <= 0:
            return {}
        raw = self.rfile.read(n)
        try:
            return json.loads(raw.decode("utf-8"))
        except Exception:                                # noqa: BLE001
            return {}

    # ── GET ───────────────────────────────────────────────────────────────
    def do_GET(self):                                    # noqa: N802 —— http.server 的命名约定
        u = urllib.parse.urlparse(self.path)
        route, qs = u.path, urllib.parse.parse_qs(u.query)
        one = lambda k, d=None: qs.get(k, [d])[0]        # noqa: E731 —— 取第一个查询值

        try:
            # ── 静态资源 ──
            if route in ("/", "/index.html"):
                return self._file(STATIC_DIR / "index.html")
            if route.startswith("/static/"):
                rel = route[len("/static/"):]
                # 防目录穿越： ../ 一律拒掉
                if ".." in rel or rel.startswith("/"):
                    return self._error("非法路径", 403)
                return self._file(STATIC_DIR / rel)

            # ── API ──
            if route == "/api/drives":
                return self._json(dict(drives=handlers.list_drives()))

            if route == "/api/browse":
                return self._json(handlers.browse(one("path")))

            if route == "/api/inspect":
                d = one("dir")
                if not d:
                    return self._error("缺少 dir 参数")
                return self._json(handlers.inspect_batch(d))

            if route == "/api/tasks":
                return self._json(dict(tasks=MANAGER.list_all()))

            if route == "/api/status":
                tid = one("tid")
                task = MANAGER.get(tid) if tid else None
                if task is None:
                    return self._error("任务不存在", 404)
                since = int(one("since", 0) or 0)
                snap = task.snapshot(since=since)
                snap["running"] = task.is_alive()
                return self._json(snap)

            if route == "/api/download":
                # 只按 tid + 序号取产物，不接受任意路径 —— 避免被当文件服务器用
                task = MANAGER.get(one("tid"))
                if task is None:
                    return self._error("任务不存在", 404)
                try:
                    item = task.outputs[int(one("i", 0) or 0)]
                except Exception:                        # noqa: BLE001
                    return self._error("产物序号不存在", 404)
                return self._file(Path(item["path"]), download=True)

            if route == "/api/ping":
                return self._json(dict(ok=True))

            if route == "/favicon.ico":
                return self._send(204, b"", "image/x-icon")

        except Exception as e:                           # noqa: BLE001
            import traceback
            traceback.print_exc()
            return self._error(f"{type(e).__name__}: {e}", 500)

        return self._error("未找到接口 " + route, 404)

    # ── POST ──────────────────────────────────────────────────────────────
    def do_POST(self):                                   # noqa: N802
        route = urllib.parse.urlparse(self.path).path
        payload = self._read_json()

        try:
            if route == "/api/inspect":
                return self._json(handlers.inspect_batch(str(payload.get("dir", ""))))

            if route == "/api/upload":
                files = payload.get("files") or []
                name = str(payload.get("batch_name") or "上传批次")
                return self._json(dict(ok=True, **handlers.save_uploads(files, name)))

            if route == "/api/run":
                if not payload.get("data_dir"):
                    return self._error("缺少 data_dir")
                cfg = handlers.make_config(payload)

                # 回调签名必须是 (task) —— Task 由 TaskManager 传进来。
                # 不要在这里写 `lambda: run_job(cfg, task)`：线程 start() 之后
                # 可能立刻执行，那时 `task = MANAGER.submit(...)` 还没赋值完，
                # 会抛 UnboundLocalError。
                def job(t):
                    handlers.run_job(cfg, t)

                try:
                    task = MANAGER.submit(job, label=cfg.data_dir.name)
                except RuntimeError as e:
                    # 已有任务在跑 —— 这是用户操作顺序问题，不是故障
                    return self._error(str(e), 409)
                return self._json(dict(ok=True, tid=task.tid, label=task.label,
                                       out_excel=str(cfg.out_excel)))

            if route == "/api/reveal":
                return self._json(handlers.reveal(str(payload.get("path", ""))))

        except Exception as e:                           # noqa: BLE001
            import traceback
            traceback.print_exc()
            return self._error(f"{type(e).__name__}: {e}", 500)

        return self._error("未找到接口 " + route, 404)

    # ── 静态文件输出 ──────────────────────────────────────────────────────
    def _file(self, path: Path, download: bool = False):
        p = Path(path)
        if not p.exists() or not p.is_file():
            return self._error(f"文件不存在：{p.name}", 404)
        ctype = _MIME.get(p.suffix.lower()) or mimetypes.guess_type(str(p))[0] \
            or "application/octet-stream"
        body = p.read_bytes()
        extra = {}
        if download:
            # RFC 5987：中文文件名要用 filename* 这段 UTF-8 编码形式，否则乱码
            fname = urllib.parse.quote(p.name)
            extra["Content-Disposition"] = f"attachment; filename*=UTF-8''{fname}"
        return self._send(200, body, ctype, extra)

    # ── 日志 ──────────────────────────────────────────────────────────────
    def log_message(self, fmt: str, *args):
        """
        默认实现会把每个请求打到 stderr。两个改动的理由：

        1. 轮询请求每 700ms 一次，刷屏会淹掉真正的批处理日志 —— 这里跳过它们。
        2. **必须写 sys.__stdout__ 而不是 sys.stdout**。
           sys.stdout 是进程全局的；任务线程在跑的时候会把它换成 _Tee，
           如果这里用 sys.stdout，HTTP 线程自己的日志就会被吸进任务日志里
           （曾经出现过网页状态和批处理进度串在一起的现象）。
           sys.__stdout__ 是解释器保存的原始 stdout，不会被任何重定向影响。
        """
        msg = fmt % args
        if "/api/status" in msg or "/api/ping" in msg:
            return
        try:
            stream = sys.__stdout__ or sys.stdout
            stream.write(f"  [http] {msg}\n")
            stream.flush()
        except Exception:
            pass


# ═════════════════════════════════════════════════════════════════════════════
# 服务启动
# ═════════════════════════════════════════════════════════════════════════════
def find_free_port(start: int = 8765, tries: int = 40) -> int:
    """
    从 start 开始找第一个能绑定的端口。

    为什么要自动找：上一次没关掉、或者别的软件占了 8765，
    用户双击图标却打不开，是最难自查的一类故障。自动顺延最省事。
    """
    import socket
    for p in range(start, start + tries):
        with socket.socket() as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    raise RuntimeError(f"{start}~{start+tries} 端口全被占用，请检查是否有残留进程")


def serve(host: str = "127.0.0.1", port: int = 0, open_browser: bool = True,
          quiet_open: bool = False) -> None:
    """
    启动 Web 服务并阻塞当前线程。

    参数
    ----
    host : 只建议 127.0.0.1。改成 0.0.0.0 会把本机的文件系统暴露到局域网。
    port : 0 表示自动挑一个空闲端口。
    open_browser : 是否自动拉起默认浏览器
    """
    port = port or find_free_port()
    httpd = ThreadingHTTPServer((host, port), Handler)
    url = f"http://{host}:{port}/"

    print("=" * 68)
    print(f"  锆石 U-Pb 数据处理工具  v{__version__}")
    print(f"  界面地址：{url}")
    print("  关闭这个黑窗口即可退出；也可以在网页里点停止服务")
    print("=" * 68)

    if open_browser and not quiet_open:
        # 延后一点点再开浏览器，等服务真正 listen 起来，避免抢跑变成"无法访问"
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止服务。")
    finally:
        httpd.server_close()
