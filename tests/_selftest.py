"""
最小自检运行器 —— 让每个测试文件既能被 pytest 收集，也能直接 `python` 跑。

为什么不直接用 pytest
----------------------
这个工具要能在实验室机器上验证自己没坏掉。pytest 是开发依赖（`.[dev]`），
那台"只用来跑数据"的机器上未必有。所以约定：

  * 测试函数只写 `assert`，不用任何 pytest 特性
    （fixture、参数化、tmp_path 一律别用）；
  * 需要临时文件就用标准库的 `tempfile`；
  * 每个测试文件末尾配一行 `raise SystemExit(_selftest.run(globals()))`。

这样 `python tests/test_xxx.py` 和 `pytest` 两种跑法都成立，不需要装任何东西。
"""
from __future__ import annotations


def run(namespace: dict) -> int:
    """执行 namespace 里所有 test_* 函数，返回进程退出码（0 = 全过）。"""
    # 下面的输出里有中文。Windows 上把输出重定向到文件/管道时，stdout 会用
    # locale 编码（cp1252 之类），第一句 print 就 UnicodeEncodeError 把进程带走。
    # 这个坑在 CI 的 windows job 上真实发生过，见 druid/console.py。
    from druid.console import ensure_utf8_streams
    ensure_utf8_streams()

    tests = [(n, f) for n, f in sorted(namespace.items())
             if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  ok    {name}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name}: {e}")
        except Exception as e:                      # noqa: BLE001
            # 意外异常也要报成失败，而不是让它中断整个文件 —— 否则后面
            # 那些本来能过的测试就白跑了，排查面也少了一半。
            failed += 1
            print(f"  ERROR {name}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} 通过")
    return 1 if failed else 0
