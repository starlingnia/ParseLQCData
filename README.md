### 执行程序使用的源码与核心函数清单

| 组件/文件 | 头文件路径 | 核心函数 / 职责 |
| :--- | :--- | :--- |
| `apps/main.cpp` | `ParseLQCData/Registry.h` | `main()`: `registry::TaskRegistry::instance().all()` 动态匹配并分发命令行子任务 |
| `build/meson_multisrc_task.cpp` | `ParseLQCData/Core/PhysicsSetup.h`, `ParseLQCData/Core/MesonPipeline.h` | 专职多源任务注册：`REGISTER_TASK("meson_multi", ...)` 与 `REGISTER_TASK("meson_multi_all", ...)` |
| `build/meson_singlesrc_task.cpp` | `ParseLQCData/Core/PhysicsSetup.h`, `ParseLQCData/Core/MesonPipeline.h` | 专职单源任务注册：`REGISTER_TASK("meson_single", ...)` 与 `REGISTER_TASK("meson_single_all", ...)` |
| `build/meson_task.cpp` | `ParseLQCData/Core/PhysicsSetup.h`, `ParseLQCData/Core/MesonPipeline.h` | 包含单源与多源的复合任务文件（向后兼容） |
| `build/condensate_task.cpp` | `CondensateAnalysis/CondensateExtractor.h` | 手征凝聚任务注册：`REGISTER_TASK("condensate", ...)` 与 `REGISTER_TASK("condensate_all", ...)` |
| `docs/physics_setup.py` | 全局单一事实源 | 定义 `CHANNEL_CONFIGS` 与 `CONDENSATE_CONFIGS` 并导出为 `docs/physics_setup.json` |
| `src/IOdata/*.cpp` | `IOdata/FileReader.h`, `IOdata/DirectoryScanner.h`, `IOdata/FileWriter.h` | 1. `iodata::scan_natural_sorted()`: 自然排序扫描<br>2. `iodata::read_file_to_string()`: 文件读取<br>3. `iodata::write_matrix_to_csv()`: 矩阵写入 |
| `src/Statistics/*.cpp` | `Statistics/Resampling.h` | 1. `lqcd::stats::compute_block_binning()`: 块平均<br>2. `lqcd::stats::compute_jackknife()`: 留一 Jackknife 矩阵<br>3. `lqcd::stats::compute_folding()`: 时间反演折叠 |
| `src/MesonAnalysis/*.cpp` | `MesonAnalysis/MesonExtractor.h` | 1. `extract_single_file_averaged_block()`: 16 空间源循环平移平均<br>2. `extract_single_file_singlesrc_block()`: 单源抽取 |
| `src/CondensateAnalysis/*.cpp` | `CondensateAnalysis/CondensateExtractor.h` | 1. `extract_pbp_from_xml()`: 标量夸克凝聚 XML 正则抽取<br>2. `process_chiral_condensate()`: 扣除残余质量发散项与重整化 |
| `src/core/MesonPipeline.cpp` | `ParseLQCData/Core/MesonPipeline.h` | `lqcd::meson::run_meson_pipeline()`: 底层多线程并行池调度，全流水线串联归约 |
| `tools/Services/MesonService.cpp` | C ABI 导出 | `extern "C" run_meson_pipeline_c_api`, `run_meson_multisrc_pipeline_c_api`, `run_meson_singlesrc_pipeline_c_api` |
| `tools/Services/CondensateService.cpp` | C ABI 导出 | `extern "C" run_chiral_condensate_c_api`: 导出供 Python ctypes 桥接 |

---

### 由您使用 gc 编译的完整指令

您的环境别名为：`gc` (例如 `g++-mp-15 -std=c++26 -Iinclude`)。

#### 1. 编译原生 CLI 执行程序 (bin/pars)

```bash
mkdir -p bin

# 方案 A: 仅编译多源任务 (极简单组件)
gc apps/main.cpp \
  src/IOdata/*.cpp \
  src/Statistics/*.cpp \
  src/MesonAnalysis/* \
  src/core/* \
  build/meson_multisrc_task.cpp \
  -o bin/pars

# 方案 B: 仅编译单源任务 (极简单组件)
gc apps/main.cpp \
  src/IOdata/*.cpp \
  src/Statistics/*.cpp \
  src/MesonAnalysis/* \
  src/core/* \
  build/meson_singlesrc_task.cpp \
  -o bin/pars

# 方案 C: 全功能 CLI (多源 + 单源 + 手征凝聚全注册)
gc apps/main.cpp \
  src/IOdata/*.cpp \
  src/Statistics/*.cpp \
  src/MesonAnalysis/* \
  src/CondensateAnalysis/* \
  src/core/* \
  build/meson_multisrc_task.cpp \
  build/meson_singlesrc_task.cpp \
  build/condensate_task.cpp \
  -o bin/pars
```

#### 2. 分开编译供 Python 调用的独立共享动态库 (build/liblqcd_meson.dylib 与 build/liblqcd_condensate.dylib)

```bash
mkdir -p build

# 介子关联函数动态库
gc -O3 -Wall -fPIC -shared \
  tools/Services/MesonService.cpp \
  src/core/MesonPipeline.cpp \
  src/MesonAnalysis/MesonExtractor.cpp \
  src/Statistics/Resampling.cpp \
  src/IOdata/*.cpp \
  -o build/liblqcd_meson.dylib

# 手征凝聚动态库
gc -O3 -Wall -fPIC -shared \
  tools/Services/CondensateService.cpp \
  src/CondensateAnalysis/CondensateExtractor.cpp \
  src/Statistics/Resampling.cpp \
  src/IOdata/*.cpp \
  -o build/liblqcd_condensate.dylib
```

#### 3. 编译 C++ 单元测试程序 (仅在需要单独测试底层算子时编译)

```bash
gc -O3 -Wall \
  tests/test_meson_pipeline.cpp \
  src/IOdata/*.cpp \
  src/Statistics/*.cpp \
  src/MesonAnalysis/*.cpp \
  src/CondensateAnalysis/*.cpp \
  src/core/*.cpp \
  -o bin/unit_tests
```

---

### 编译后由您执行的命令

```bash
# 1. 验证命令列表与从 docs/ 加载的信道任务
./bin/pars

# 2. 测量单个信道 (自动从 docs/physics_setup.json 动态读入信道算子结构)
./bin/pars meson_multi 17 AV 4
./bin/pars meson_single 17 S 4

# 3. 批量测量特定 Beta 下所有信道
./bin/pars meson_all 18 multi

# 4. 测量手征凝聚
./bin/pars condensate
```

---

## 编译与执行操作指南 (User Compilation & Execution Guide)

### 1. 编译构建 C++ 独立组件与目标程序

本项目全面模块化，组件彼此完全解耦，提供两种编译方式：

#### 方式 A：使用您的 `gc` 直接编译 (推荐)
参见上方说明。

#### 方式 B：使用 `xmake` 模块化构建各组件

```bash
# 1. 一键编译所有业务组件与产物 (动态库 + CLI 主程序)
xmake

# 2. 单独编译指定原子组件 (如仅编译 I/O 模块或核心管线模块)
xmake build iodata
xmake build statistics
xmake build meson_analysis
xmake build condensate_analysis
xmake build core

# 3. 编译并运行 C++ 单元测试 (tests/ 仅在需要测试时构建)
xmake build unit_tests
xmake run unit_tests
```

---

### 2. 执行数据全量复现与分析流水线 (极简、拆分、无图标)

用户可以通过以下几种方式执行计算：

#### 方式 A：分开执行专职流水线 (强子多源、强子单源与手征完全解耦)
```bash
# 1. 运行多源介子流水线 (抽取 + meff + 拟合)
bash scripts/run_meson_multi.sh
# 或: uv run python scripts/reproduce_meson_multi.py

# 2. 运行单源介子流水线 (抽取 + meff + 拟合)
bash scripts/run_meson_single.sh
# 或: uv run python scripts/reproduce_meson_single.py

# 3. 运行手征凝聚提取流水线 (重整化手征凝聚)
bash scripts/run_condensate.sh
# 或: uv run python scripts/reproduce_condensate.py

# 4. 生成科学图表
uv run python scripts/plot_meson.py
uv run python scripts/plot_condensate.py
```

#### 方式 B：一键端到端总控
```bash
# 运行介子全量流水线 (多源 + 单源 + 绘图)
bash scripts/run_meson.sh

# 运行全项目所有物理流水线
bash scripts/run_all.sh
```

#### 方式 C：独立测试套件 (tests/ 目录下，按需执行)
```bash
# 1. 运行 Binsize 自相关消除与误差饱和度极简测试 (支持指定 beta 与 channel)
uv run python tests/test_binsize_autocorr.py 17 AV

# 2. 运行与 ana 基准 260 项数值精度回归对比
uv run python tests/compare_with_ana.py

# 3. 运行多线程扩展性与极端压力测试
uv run python tests/benchmark_stress_test.py

# 4. 一键执行全部测试套件
bash tests/run_tests.sh
```

---

### 3. C++ 原生 CLI 驱动方式 (任务注册表)

系统支持使用原生可执行文件 `bin/ParseLQCData` 直接派发各项底层计算任务：

```bash
# 查看所有注册的任务及用法
./bin/ParseLQCData

# 1. 运行单信道多源提取任务: [beta] [channel] [binsize]
./bin/ParseLQCData meson_multi 17 AV 4

# 2. 运行单信道单源提取任务: [beta] [channel] [binsize]
./bin/ParseLQCData meson_single 17 S 4

# 3. 批量运行特定 Beta 下全部 6 种物理信道: [beta] [mode: multi|single]
./bin/ParseLQCData meson_all 18 multi

# 4. 运行手征凝聚快速提取任务: [base_dir] [ml] [ms] [mres] [zm] (缺省自动使用本地数据)
./bin/ParseLQCData condensate
```

---

## 数值回归对比与极限承压验证结论

### 1. 全量数值精度回归 (与 ana/ 基准对比)
详细报告见：[`docs/verification_report.md`](file:///Users/junxiongnie/code/build/ParseLQCData/docs/verification_report.md)

| 校验模块 | 测试项总数 | 通过数 | 通过率 | 精度指标 |
| :--- | :--- | :--- | :--- | :--- |
| **多源强子关联函数** | 42 组 (7 Betas x 6 Channels) | 42 | **100.0%** | 最大偏差 $< 1.67 \times 10^{-16}$ (双精度浮点极限) |
| **单源强子关联函数** | 42 组 (7 Betas x 6 Channels) | 42 | **100.0%** | 最大偏差 $< 1.39 \times 10^{-17}$ (双精度浮点极限) |
| **多源有效质量 ($m_{\text{eff}}$)** | 42 组 (7 Betas x 6 Channels) | 42 | **100.0%** | 最大偏差 $< 2.38 \times 10^{-6}$ (非线性求解收敛极限) |
| **单源有效质量 ($m_{\text{eff}}$)** | 42 组 (7 Betas x 6 Channels) | 42 | **100.0%** | 最大偏差 $< 3.70 \times 10^{-6}$ (非线性求解收敛极限) |
| **非线性平台拟合总结** | 84 组拟合参数对比 | 84 | **100.0%** | 最大偏差 $< 1.11 \times 10^{-16}$ (完全一致) |
| **手征凝聚重整化结果** | 8 组温度参数对比 | 8 | **100.0%** | 最大偏差 $< 1.00 \times 10^{-16}$ (完全一致) |
| **总体测试通过率** | **260 项物理测试** | **260 项** | **100.00%** | **全量零回归，100% 严格复现** |

### 2. 高并发与极限负载承压指标
详细报告见：[`docs/stress_test_report.md`](file:///Users/junxiongnie/code/build/ParseLQCData/docs/stress_test_report.md)

* **并行扩展加速比**：单进程多线程并发最高达到 **4.84x 加速**，吞吐量达到 **538.6 cfgs/s**；
* **I/O 与解析吞吐率**：针对 396 个大文本文件（841.01 MB），解析耗时仅 **0.85 秒**，有效 I/O 解析带宽达到 **986.31 MB/秒**；
* **内存稳定性与零泄漏**：连续 12 轮高负荷大规模数据提取与 Jackknife 重采样，物理内存净增长仅 **+0.05 MB**，严格验证了 RAII 内存零泄露。

---

## 项目完整目录树架构

```text
/Users/junxiongnie/code/build/ParseLQCData/
├── apps/
│   └── main.cpp                         # C++ 原生 CLI 驱动入口 (TaskRegistry 统一调度)
├── bin/
│   ├── pars                             # 编译生成的 CLI 执行程序
│   └── unit_tests                       # C++ 单元测试二进制 (按需编译)
├── build/
│   ├── meson_task.cpp                   # 介子测量任务 (从 docs/ 动态载入信道并注册至 Registry)
│   ├── condensate_task.cpp              # 手征凝聚测量任务 (注册至 Registry)
│   ├── liblqcd_meson.dylib              # 介子关联函数专职动态库
│   └── liblqcd_condensate.dylib         # 手征凝聚专职动态库
├── data/                                # 自包含的格点物理原始数据 (6.2 GB)
│   └── readin/
│       ├── 48x16b4.13 ~ 48x16b4.30/     # 7 组温度下的强子关联函数大文本数据
│       └── L32T12beta4.17/              # 手征凝聚 XML 数据
├── dat -> data                          # 目录软链接兼容层
├── docs/                                # 物理信息与元数据中枢 (Single Source of Truth)
│   ├── physics_setup.py                 # 全局物理配置单一事实源 (温度/信道/切片/对称性)
│   ├── physics_setup.json               # 物理配置 JSON 导出版本 (C++ 与 Python 共同读取)
│   ├── physics_inputs.md                # 物理输入构型元数据规范
│   ├── physics_results.md               # 物理输出观测量计算标准
│   ├── verification_report.md           # 260 项数值回归验证报告
│   ├── stress_test_report.md            # 极限负载与高并发性能承压报告
│   └── figures/                         # 出版级科学图像目录
│       ├── channel_comparisons/         # 42 组信道单/多源拟合对比图 (b4.13 ~ b4.30)
│       ├── simulation_results.png       # 有效质量温度演化全景图 (带相变过渡带)
│       ├── symmetry_breaking_ana_setup.png # 手征与轴对称性破缺恢复签名图
│       ├── conden_re_plot.png           # 重整化手征凝聚随温度演化图
│       └── conden_plot.png              # 轻夸克原始手征凝聚图
├── include/                             # 解耦的独立算子头文件与核心管线接口
│   ├── IOdata/                          # 基础文本 I/O 与零拷贝解析头文件
│   │   ├── DirectoryScanner.h
│   │   ├── FastParser.h
│   │   ├── FileReader.h
│   │   ├── FileWriter.h
│   │   └── IOdata.h
│   ├── MesonAnalysis/                   # 基础强子提取算子头文件 (MesonExtractor.h)
│   ├── CondensateAnalysis/              # 基础手征凝聚算子头文件 (CondensateExtractor.h)
│   ├── Statistics/                      # 统计与重采样算子头文件 (Resampling.h)
│   └── ParseLQCData/                    # 项目核心流水线接口与注册表
│       ├── Core/
│       │   ├── PhysicsSetup.h           # 动态读取 docs/physics_setup.json 映射解析器
│       │   ├── LQCDataTypes.h           # 物理强子信道与数据结构定义
│       │   └── MesonPipeline.h          # C++ 高并发调度管线接口
│       └── Registry.h                   # 统一子任务注册驱动器
├── output/                              # 隔离的物理计算数据产物目录
│   ├── pickdata/                        # 多源关联函数 (save_, sym_, dr_, err_)
│   ├── pickdata-singlesrc/              # 单源关联函数 (save_, sym_, dr_, err_)
│   ├── ratio_results/                   # 多源有效质量曲线 (meff_<ch>.csv)
│   ├── ratio_results-singlesrc/         # 单源有效质量曲线 (meff_<ch>.csv)
│   ├── simulateresult/                  # 多源贝叶斯拟合及 all_fits_summary.csv
│   ├── simulateresult-singlesrc/        # 单源贝叶斯拟合及 all_fits_summary.csv
│   └── condensate/                      # 手征凝聚物理结果文本 (results_rm_beta.txt)
├── pyproject.toml                       # Python 工程与依赖管理描述
├── scripts/                             # 极简科学计算流水线 (无图标、完全拆分)
│   ├── reproduce_meson.py               # 介子关联函数全量复现 Python 脚本
│   ├── plot_meson.py                    # 介子科学对比绘图生成 Python 脚本
│   ├── run_meson.sh                     # 介子计算与绘图一键 Shell 脚本
│   ├── reproduce_condensate.py          # 手征凝聚提取与重整化 Python 脚本
│   ├── plot_condensate.py               # 手征凝聚对比绘图生成 Python 脚本
│   ├── run_condensate.sh                # 手征凝聚一键 Shell 脚本
│   └── run_all.sh                       # 极简端到端总控 Shell 脚本
├── src/                                 # C++ 算子实现与 Python 基础管道
│   ├── IOdata/                          # 基础文本解析与文件扫描实现
│   ├── MesonAnalysis/                   # 基础强子关联函数提取实现
│   ├── CondensateAnalysis/              # 基础手征凝聚提取实现
│   ├── Statistics/                      # 重采样与时间反演Folding实现
│   ├── core/                            # C++ 底层统一调度管线实现 (MesonPipeline.cpp)
│   └── parselqcdata/                    # Python 整合基础管道
│       ├── meson_pipeline.py            # 介子端到端分析管道
│       ├── effective_mass.py            # 有效质量超越方程求解器
│       └── plateau_fit.py               # 贝叶斯非线性平台拟合器
├── tests/                               # 纯测试套件 (按需独立执行)
│   ├── test_binsize_autocorr.py         # Binsize 自相关消除与误差饱和度极简测试
│   ├── test_meson_pipeline.cpp          # C++ 核心算子与管线单元测试
│   ├── compare_with_ana.py              # 全量数值精度回归验证测试脚本
│   ├── benchmark_stress_test.py         # 极限承压与高并发性能测试脚本
│   ├── run_all_benchmarks.py            # 架构基准性能压测总控
│   └── run_tests.sh                     # 极简一键回归与基准测试脚本
├── tools/                               # 接口服务层 (完全拆分)
│   ├── meson_orchestrator.py            # 介子关联函数专职调度器 (加载 liblqcd_meson.dylib 等)
│   ├── condensate_orchestrator.py       # 手征凝聚专职调度器 (加载 liblqcd_condensate.dylib 等)
│   ├── lqcd_orchestrator.py             # 统一调用门面封装
│   └── Services/ (Service -> Services)  # C ABI 动态库桥接层 (MesonService.cpp, CondensateService.cpp)
└── xmake.lua                            # 模块化 xmake 现代构建描述文件
```
