"""
druid.cli —— 命令行入口
=====================

这里只做两件事：
    1. 把命令行参数翻译成 BatchConfig
    2. 调用 workflow.run_batch 并把结果落到磁盘

**不含任何算法**：要看处理逻辑请去 workflow 与各算法模块。
"""
