# AGENTS.md — 给在这个仓库里干活的人（和 AI 助手）

**读这一页 + `README.md` 就够开始改代码，不必通读全部源码。**
本文件只写"读代码看不出来"的东西：怎么验证、哪些数字不许动、哪里有雷。

---

## 1. 这是什么

LA-ICP-MS 锆石 U-Pb 数据还原工具。把 Qtegra/iCAP 导出的 cps 时间序列，
变成**随深度变化的年龄**（逐窗口，含 1σ），供下游的 R 包 **ADEPT** 判定坪年龄与 MSWD。

一句话分工：**DRUID 测出来，ADEPT 读出来。** 两者关系见 `docs/druid-adept-dataflow.html`。

- 包名 `druid`（历史上叫 `upb`，见到旧叫法都是改名前的引用）
- 版本在 `druid/__init__.py:__version__`
- Python 层：**只有 numpy / pandas / matplotlib / openpyxl / xlrd**。
  **没有 scipy** —— 年龄方程的不动点迭代、二分法、卡方上尾概率都是自己实现的，
  这不是疏漏而是设计：不依赖会变的第三方行为。

## 2. 30 秒上手

```bash
# 必须用这个解释器 —— 系统 Python 没有科学计算包
PY="C:/Users/凯凯/.workbuddy/binaries/python/envs/upb/Scripts/python.exe"

# 跑一个批次
"$PY" -m druid.cli.reduce_batch --dir "<批次目录>" --out "<结果.xlsx>" --plot

# 网页界面（日常用这个）
#   双击 启动数据处理工具.bat

# 自检（不需要 pytest）
"$PY" tests/run_all.py

# 装了开发依赖的话等价于（这是 CI 用的那条路）
pip install -e ".[dev]" && python -m pytest -q
```

虚拟环境目录仍叫 `envs/upb`（历史名），两个启动 bat 同时接受 `envs/druid` 和 `envs/upb`。

## 3. 改完必须做的三件事

```bash
1) "$PY" -m pyflakes druid/            # 必须 0 项
2) "$PY" tests/run_all.py              # 必须全过
3) 跑真实批次，与基线逐列对比（见 §4）
```

**任何一步没过就不要提交。**

### CI 能替你做什么，不能做什么

`.github/workflows/ci.yml`（推上 GitHub 后自动跑，约 1 分钟）：

- ubuntu × py3.9 / 3.12 / 3.13，windows × py3.13；
- 每个平台都跑 pyflakes、`pytest`、`tests/run_all.py`（不含 pytest 的那条路）、
  以及两个命令行入口的 `--help`；
- 另一个 job 构建 wheel，核对静态资源与两个入口点真的在包里。

**它没有真实批次数据，所以 CI 全绿不代表数值没变 —— 第 3 步只能人来做。**
反过来，CI 抓的是"源码目录里跑得好好的、装完却坏掉"那类问题
（网页资源漏进 wheel 就是典型），那种本地很难发现。

写测试时注意约定：**只写 `assert`，不用任何 pytest 特性**（fixture、参数化、
`tmp_path` 都别用），需要临时文件就用 `tempfile`。因为实验室那台机器上
可能只有 numpy/pandas，`tests/run_all.py` 必须能跑。

## 4. 数值不许悄悄动

输出是要写进文章的。除了有意为之的改动，**必须证明数值没变**。

对照基线（批次 `20220301CKDB`，48 个样品测点）：

| 量 | 值 |
|---|---|
| 结果 / 标样QC / 深度剖面域 / 剖面窗口 | 83 / 2 / 91 / 1677 行 |
| 样品 206Pb/238U 中位 | 460.9 Ma |
| 协和度中位 | 101.3%，其中 90–110% 占 94% |
| 91500 偏差 / Ple 偏差 | −0.25% / +1.81% |
| 接 ADEPT 后 | 51 个坪，MSWD 中位 1.25 |

对比做法：改前先 `--out 旧.xlsx` 存一份，改后逐列比。
**比较时用 `fillna('__NA__')` 哨兵法，不要直接 `astype(str)` 比** ——
空单元格在两次读取里可能是 `None` 也可能是 `nan`，会报出假差异（踩过）。

## 5. 分层与依赖方向

```
core/        纯计算，无 I/O，可单独 import 测试   ← 谁都不依赖
io/          读写：Qtegra CSV、序列 LIST、Excel 报表
reduction/   单点还原：装载、剥蚀段识别、同位素比值与 jackknife
depth/       深度剖面：滑窗、F(τ) 分馏校正、BIC 年龄域判别、出图
webui/       网页界面（任务管理 / 业务层 / HTTP 层 / 静态页）
cli/         只做参数解析
workflow.py  编排层：BatchConfig 集中所有可调参数 + 十步 run_batch
```

**依赖只朝一个方向流**：`core ← io ← reduction ← depth ← workflow`。
想改算法去 `core`/`reduction`；想改流程顺序去 `workflow`。

## 6. 改动守则

1. **物理常量只允许在 `core/constants.py` 出现一次。** 同一个衰变常数写两处，
   改了一处忘了另一处不会报错，只会悄悄给出错误的年龄。以前 `workflow.py`
   就把 `L238` 硬编码过一遍。
2. **跨工具对接一律显式传参，不要靠推断。** 曾经做过"有 `Age68_1s` 列就认为数据
   已预处理、自动跳过时间窗口"，结果示例数据从 161 点变成 410 点、年龄从 21.6 变成
   79.8 Ma —— 一个字段的存在说明不了另一个字段的语义。
3. **年龄一律是反比方差加权平均，不是算术平均。** `weighted_mean()` 是唯一实现，
   MSWD 与它的卡方概率都在这个口径上算（与 ADEPT 同口径）。
4. **不要为了"顺手"改数值输出。** 例如 `deadtime_ns` 默认 0 是刻意的：
   设成实测值会让全部历史结果不可比。要改默认值先问人。
5. **改 `druid/` 下的东西不要动 `结果/`。** 那是工具的输出目录，装的是某一次
   真实分析的数据。它现在被 `.gitignore` 挡住了（曾经没有，代价是重写全部历史）。
6. **版本号写在两处**：`pyproject.toml` 的 `version` 与 `druid/__init__.py` 的
   `__version__`。改一处忘一处不会报错，只会让包自称 2.1.0 而跑起来打印 2.0.0
   ——网页界面那行 banner 就这么硬编码过。`tests/test_packaging.py` 会拦下来。
7. **改输出表的列就升版本**，并在 `druid/__init__.py` 的版本说明里写清为什么。
   下游（ADEPT、补充材料里的表）需要靠版本号判断两批结果可不可比。

## 7. 已知地雷

- **204Pb 必须做显著性检验**（`n_sigma_common_pb`，默认 2σ）。204Pb 只有几~几十 cps，
  204Hg 本底几百 cps，不做检验就扣普通铅等于把噪声当铅扣掉，对 207Pb 是毁灭性的。
- **207Pb/206Pb 不套 F(τ)**。它是 Pb-Pb 比值，不随坑深分馏；窗口级 207Pb 只有约 130 个计数，
  套 F(τ) 纯属注入噪声。它只参与**标样归一化**（F76）。
- **QC 二次校正是治标**。主标 91500 约 1e5 cps，样品约 1e6 cps，差一个数量级，
  单一线性因子消不掉非线性残差。正式发表前应缩小信号强度差距或做死时间校正。
- **窗口重叠让 MSWD 系统性偏小**（win=4 s / step=1 s，ρ=0.75，n_eff≈n/2.5）。
  `domains.py` 的判据 `1 + α√(2/(k−1))` 因此偏**宽松**约 0.3（k=30 给 1.53，按 k_eff=12 应 1.85）。
  方向是**漏杀混合窗口，不会错杀**。
- **`druid/deadtime_ns.txt` 里是实测死时间 14.7645 ns，但没有代码读它**。
  要用得显式传 `--deadtime-ns 14.7645`。
- **`.bat` 必须保持纯 ASCII**：cmd.exe 按当前代码页解析，中文注释会变乱码；
  而解释器路径含中文用户名。所有中文路径靠运行时变量（`%USERPROFILE%`、`%~dp0`）抵达。
  `.gitattributes` 已把 `*.bat` 钉成 CRLF。

  ⚠ **而且它真的漂移过。** 属性只在 **checkout** 时生效：一次
  `git filter-branch`（它会重写后重置工作区）就把两个 `.bat` 变成了纯 LF，
  而 `git status` **什么都没报** —— 因为入库时会规范化，git 认为内容没变。
  这种错误会一路潜伏到实验室里有人双击那个 `.bat` 才发现。
  `tests/test_packaging.py` 现在直接看磁盘字节，就是为了接住这类漂移。

- **`git status` 有时会把两个 `.bat` 报成 ` M`，那是 stat 噪音。**
  这台机器 `core.autocrlf=true`（来自 WorkBuddy 自带的 PortableGit，
  **不要去改它**），与 `.gitattributes` 的 `eol=crlf` 叠加后，只要文件被
  git 之外的东西重写过一次，`git status` 就会按 stat 判成"改过"。
  内容其实一致：`git diff` 对它们没有任何输出。

  两件事要记住：

  1. **`git add` 这两个文件是 no-op**，但它会刷新 stat，` M` 就消失了。
     想确认无害就比哈希：

     ```bash
     git rev-parse :打包成exe.bat                       # 索引里现在是什么
     git hash-object --path=打包成exe.bat 打包成exe.bat # git add 之后会变成什么
     ```

     两个哈希相同 → 不改变仓库内容。
  2. **这个 stat 脏状态会挡住 `git filter-branch`**，报
     `Cannot rewrite branches: You have unstaged changes.`。
     重写历史之前先 `git add` 一遍（或 `git update-index --refresh`）把工作区刷干净。
- **`结果/` 曾经被纳管进版本控制**：98 个文件（87 个原始 Qtegra CSV + 11 个结果表），
  而且都在**初始提交**里。等到要发布时，`git rm` 当前提交已经没用，只能重写全部
  历史才剔干净。规律：**工具的仓库不装某一次分析的数据。** 要留档就存
  "原始批次目录 + 运行参数"，那两样足够完整复现。
- **不要用 `git filter-branch` 剔文件后又指望工作区原样**：收尾时它会把工作区
  重置到重写后的 HEAD，于是那些"刚被移出索引"的文件会**从磁盘上一起消失**。
  要保住磁盘文件就先复制一份出来。真的丢了也别急，旧提交对象在 gc 之前都还在，
  `git archive <旧SHA> <目录>` 能完整捞回来（这次就是这么救回来的）。

## 8. 接 ADEPT

`剖面窗口` 表的列名就是按 ADEPT 的 Format 4 约定取的，直接喂过去即可：

```r
adept("<批次>_U-Pb结果.xlsx",
      lower_ablation_time = 0,      # DRUID 已切好自己的剥蚀窗口
      upper_ablation_time = 1e6,    # 上界放到最后一个窗口之后
      smooth = "none",              # 已做 F(τ)，再平滑会抹平真实域边界
      calibration_uncertainty = 0)  # σ 已含外部重现性，再乘 3% 是重复计入
```

四个都不可省：一个 `_1s` 列说明不了剥蚀窗口在哪里结束，ADEPT 刻意不猜。

**`_1s` 是 1σ，不是 2σ。** 结果表里 `s68_1sig` 与 `s68_2sig` 并列，拿错那列 MSWD 差 4 倍。

## 9. 术语对照

| 代码里 | 指什么 |
|---|---|
| `Tra` | 一个测点的全部家当（原始/净信号、时间轴、剥蚀区间、样品名与角色） |
| `τ`(tau) | 归一化剥蚀进度 0→1，等价于坑深 |
| `F(τ)` | 随深度变化的分馏因子（只用于 Pb/U） |
| `f206` | 普通铅在 206 中的占比 |
| `s68` / `s76` | 206Pb/238U 与 207Pb/206Pb 的**相对** 1σ |
| 域 / domain | 深度剖面上一段年龄自洽的区间（DRUID 的用词） |
| 坪 / plateau | 同上，但由 ADEPT 用它的四步过滤判据选出（ADEPT 的用词） |
| `σ_ext` | 外部重现性，由监控标样散度反推 |
