"""
druid.webui —— 网页界面子包
=========================

目标是让**完全不用命令行的人**也能跑这套流水线：
双击一个 bat → 浏览器自动打开 → 点几下选个文件夹 → 出 Excel。

模块分工
--------
    tasks.py      后台任务管理（线程 + stdout 重定向抓日志）
    handlers.py   业务实现（浏览/探测/构造配置/跑任务/上传）—— 不认识 HTTP
    server.py     HTTP 层（路由、静态文件、JSON 应答）—— 不认识 BatchConfig
    static/       前端页面（原生 HTML/CSS/JS，无构建步骤、无 npm）

分层的好处：handler 可以直接被命令行或定时任务复用，
将来想换成 Flask 也只需要重写 server.py。

外部只需 import：
    from druid.webui.handlers import inspect_batch, run_job
    from druid.webui.server import serve
"""

from . import handlers, tasks
from .server import find_free_port, serve
from .tasks import MANAGER, Task, TaskManager, TaskState

__all__ = ["serve", "find_free_port", "handlers", "tasks",
           "MANAGER", "Task", "TaskManager", "TaskState"]
