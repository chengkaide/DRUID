"""
druid.webui.tasks —— 后台任务管理（线程 + 日志捕获）
====================================================

为什么需要这一层
----------------
一次批处理要跑几十秒到几分钟，如果直接在 HTTP 请求的处理线程里跑，
浏览器会一直转圈直到超时。所以必须：

    HTTP 线程（毫秒级响应）  ←→  工作线程（分钟级计算）

两者之间需要一个**线程安全的中间人**，这就是 TaskManager 的职责。

日志是怎么抓到的
----------------
druid/workflow.py 里所有进度都走 `_log()` → `print()` → `sys.stdout`。
这意味着我们**完全不需要改动 workflow 一行代码**，只要在工作线程启动期间
临时把 `sys.stdout` 换成自己的对象，就能截获全部日志。

这个设计很关键：任何 print 都逃不掉，包括 matplotlib、numpy 警告
（stderr 也一并捕获），用户看到的就是真实的运行过程。

进度百分比怎么来的
------------------
workflow.py 的日志里有 `[1]` ~ `[6]` 六个阶段标记，例如：

    [1] 序列：83 个测点  →  ...
    [2] 主标 91500   R68_ref=0.179280 ...
    ...
    [6] QC 二次校正：Ple 实测 ...

用正则抓到最后一个出现的阶段号 n，进度就是 n/6。
这不是精确百分比，但对"还剩多久"的判断足够用，且因为是粗粒度 stage，
在某些阶段卡住很久时进度条会停住 —— 这是诚实的表现，
比用假动画假装在工作要好。
"""
from __future__ import annotations

import io
import re
import sys
import threading
import time
import traceback
import uuid
from collections import deque
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional

# 匹配形如 "[3] 外部重现性..." 的阶段标记，用于推断进度
_STAGE_RE = re.compile(r"^\s*\[([1-6])\]")

# 工作流总共有六个阶段（见 workflow.run_batch）
_TOTAL_STAGES = 6


class TaskState(str, Enum):
    """任务状态。用 str 混合继承，便于直接 JSON 序列化。"""
    PENDING = "pending"      # 已创建，尚未开始
    RUNNING = "running"      # 正在跑
    DONE = "done"            # 正常结束
    ERROR = "error"          # 抛异常结束
    CANCELLED = "cancelled"  # 用户中止


class _Tee(io.TextIOBase):
    """
    把输出同时写到「终端」和「内存队列」的流对象。

    为什么要保留终端输出：万一 Web 界面崩了或用户没开浏览器，
    命令行窗口里仍然能看到完整过程，便于排查。

    为什么要按行切分：HTTP 轮询时按"行号游标"增量取，不用每次传全量日志。
    """

    def __init__(self, sink: Callable[[str], None], original=None):
        self._sink = sink            # 每来一行就回调一次
        self._original = original    # 真正的 stdout
        self._buf = ""               # 跨 write() 调用的半行残留

    def write(self, s: str) -> int:
        # print() 会多次调用 write（内容 + 换行），所以要自己拼行
        self._buf += s
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            self._sink(line)
            if self._original is not None:
                try:
                    self._original.write(line + "\n")
                except Exception:
                    pass  # 终端不可写也不能影响计算
        if self._original is not None:
            try:
                self._original.flush()
            except Exception:
                pass
        return len(s)

    def flush(self):
        # 把不足一行的残留也推给回调，保证最后一句不会丢
        if self._buf.strip():
            self._sink(self._buf)
            self._buf = ""
        if self._original is not None:
            try:
                self._original.flush()
            except Exception:
                pass

    def isatty(self) -> bool:
        """有些库（如 tqdm）据此决定是否输出 ANSI 控制符，这里统一否认。"""
        return False


@dataclass
class Task:
    """一个批处理任务。所有字段都可被 JSON 序列化后直接返回给前端。"""
    tid: str                                  # 任务 ID（uuid4 前 8 位）
    label: str = ""                           # 人类可读的批次名
    state: TaskState = TaskState.PENDING
    logs: deque = field(default_factory=lambda: deque(maxlen=5000))
    t_created: float = field(default_factory=time.time)
    t_started: Optional[float] = None
    t_ended: Optional[float] = None
    stage: int = 0                            # 当前阶段 0~6
    rows: int = 0                             # 已处理行数（预留，暂未用）
    error: Optional[str] = None               # 失败时的 traceback 文本
    outputs: List[Dict[str, str]] = field(default_factory=list)  # 产物清单
    summary: Dict[str, Any] = field(default_factory=dict)        # 结果摘要
    _thread: Optional[threading.Thread] = field(default=None, repr=False)

    # ── 快照：给前端轮询用 ────────────────────────────────────────────────
    def snapshot(self, since: int = 0) -> dict:
        """
        返回任务状态。since 是上次取到的日志行号，只返回之后的新行，
        避免轮询时反复传输全部日志（长任务日志可达数千行）。
        """
        lines = list(self.lines)
        return dict(
            tid=self.tid,
            label=self.label,
            state=self.state.value,
            stage=self.stage,
            # 进度 = 阶段占比；已完成给满
            progress=100.0 if self.state in (TaskState.DONE,) else round(self.stage / _TOTAL_STAGES * 100, 1),
            total_lines=len(lines),
            new_lines=lines[since:] if since < len(lines) else [],
            elapsed=(self.t_ended or time.time()) - (self.t_started or time.time()),
            error=self.error,
            outputs=self.outputs,
            summary=self.summary,
        )

    @property
    def lines(self) -> List[str]:
        # deque 不能直接切片，转 list 后切。日志有上限，代价可接受。
        return list(self.logs)

    def is_alive(self) -> bool:
        return self._thread is not None and self._thread.is_alive()


class TaskManager:
    """
    任务的注册表。全局唯一实例即可 —— 这是个人工具，不是多租户服务器。

    用 threading.Lock 保护字典，因为 HTTP 线程和工作线程会同时访问它。
    """

    def __init__(self, keep: int = 20):
        self._tasks: Dict[str, Task] = {}
        self._lock = threading.Lock()
        self._keep = keep          # 最多保留多少个历史任务，防止内存无限增长

    # ── 创建并立即启动 ────────────────────────────────────────────────────
    def submit(self, func: Callable[[Task], None], label: str = "") -> Task:
        """
        提交一个函数到后台线程执行。

        func 的签名必须是 `def job(task: Task)` —— **Task 由回调参数传进去，
        而不是靠闭包去抓外层变量**。这是个刻意的决定：

            # 危险写法（曾经因此崩过）：
            task = MANAGER.submit(lambda: do(task))   # ← 线程启动可能快过赋值，
                                                      #   运行时 UnboundLocalError
            # 正确写法：
            MANAGER.submit(lambda task: do(task))

        线程一旦 start() 就可能立刻执行，而 submit() 还没返回、外层变量还没
        绑定上。这种竞态不是每次都触发，所以特别难查。把 Task 当参数传进去，
        从根上没有这个问题。

        同一时刻只允许一个任务在跑（见下方 _lock 检查）。批处理是 CPU 密集的，
        并行既没有收益，还会让 stdout 重定向互相打架。
        """
        with self._lock:
            alive = [t for t in self._tasks.values() if t.is_alive()]
            if alive:
                running = alive[0]
                raise RuntimeError(
                    f"已有任务正在运行（{running.label}，#{running.tid}）。"
                    f"批处理是 CPU 密集的，同时跑多个没有收益，等它结束再提交。")

        t = Task(tid=uuid.uuid4().hex[:8], label=label or "未命名任务")
        with self._lock:
            self._tasks[t.tid] = t
            self._trim()

        th = threading.Thread(target=self._wrap(t, func), name=f"task-{t.tid}", daemon=True)
        t._thread = th
        th.start()
        return t

    # ── 线程主函数的包装：负责状态切换、日志重定向、异常兜底 ──────────────
    def _wrap(self, t: Task, func: Callable[[Task], None]) -> Callable[[], None]:
        def runner():
            t.state = TaskState.RUNNING
            t.t_started = time.time()

            def on_line(line: str):
                # 每条日志进来，先更新阶段号，再入队
                m = _STAGE_RE.match(line)
                if m:
                    t.stage = max(t.stage, int(m.group(1)))
                elif line.startswith("    ") and t.stage == 0:
                    # "[1]" 之前还有跳过信息等细节行，也给一点进度感
                    pass
                t.logs.append(line)

            tee_out = _Tee(on_line, sys.stdout)
            tee_err = _Tee(on_line, sys.stderr)

            old_out, old_err = sys.stdout, sys.stderr
            sys.stdout, sys.stderr = tee_out, tee_err
            try:
                func(t)          # ← Task 作为参数传入，避免闭包竞态
                t.state = TaskState.DONE
                t.stage = _TOTAL_STAGES
            except BaseException as e:       # noqa: BLE001 —— 必须兜住所有异常
                t.state = TaskState.ERROR
                t.error = f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}"
                t.logs.append(f"!! 运行失败：{type(e).__name__}: {e}")
            finally:
                tee_out.flush()
                tee_err.flush()
                sys.stdout, sys.stderr = old_out, old_err
                t.t_ended = time.time()

        return runner

    # ── 查询 ──────────────────────────────────────────────────────────────
    def get(self, tid: str) -> Optional[Task]:
        with self._lock:
            return self._tasks.get(tid)

    def list_all(self) -> List[dict]:
        with self._lock:
            tasks = sorted(self._tasks.values(), key=lambda x: x.t_created, reverse=True)
        return [dict(tid=t.tid, label=t.label, state=t.state.value,
                     created=t.t_created) for t in tasks]

    def _trim(self):
        """超出 keep 数量时淘汰最老的任务，防止长时间运行内存膨胀。"""
        if len(self._tasks) <= self._keep:
            return
        oldest = sorted(self._tasks.values(), key=lambda x: x.t_created)[:-self._keep]
        for t in oldest:
            if not t.is_alive():
                self._tasks.pop(t.tid, None)


# 全局单例：整个 Web 服务共用这一个管理器
MANAGER = TaskManager()
