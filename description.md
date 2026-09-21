# ParseLQCData 工程可执行程序与最佳实践指南

本项目基于 **Pitchfork Directory Standard (PDS)** 规范，采用 **“底层高性能模块化 C++ (遵循 GSL 规范) + 上层轻量 Python 科学流转流水线”** 的混合架构。

本文档系统性梳理本仓库内**所有可编译生成的执行文件、共享动态库、依赖关系、已注册 Task、功能参数以及各场景下的最佳实践**。

---

## 一、 可编译生成的执行文件全景 (位于 `bin/`)

本工程通过模块化解耦，可独立编译出以下三类可执行二进制程序（统一输出至 [bin/](file:///Users/junxiongnie/code/algo/ParseLQCData/bin/) 目录）：

### 1. `bin/ParseLQCData` —— 原生 CLI 任务调度主程序 (核心驱动)

- **定位**：整个格点量子色动力学（LQCD）数值计算分析的核心命令行驱动入口。
- **入口源码**：[apps/main.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/apps/main.cpp)
- **依赖组件与架构解耦**：
  - **任务组件 (Object Library)**：来自 `build/*_task.cpp`（[build/condensate_task.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/build/condensate_task.cpp)、[build/meson_task.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/build/meson_task.cpp)、[build/meson_singlesrc_task.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/build/meson_singlesrc_task.cpp)、[build/meson_multisrc_task.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/build/meson_multisrc_task.cpp)）。这些组件**单独只编译为 `.o` 目标文件**，通过静态初始化向注册表注入任务。
  - **核心计算流水线 (`core`)**：`src/core/MesonPipeline.cpp`、`src/core/CondensatePipeline.cpp`
  - **业务抽取组件**：`meson_analysis` (`src/MesonAnalysis/MesonExtractor.cpp`)、`condensate_analysis` (`src/CondensateAnalysis/CondensateExtractor.cpp`)
  - **统计重采样 (`statistics`)**：`src/Statistics/Resampling.cpp`
  - **极速 I/O (`iodata`)**：`src/IOdata/DirectoryScanner.cpp`、`src/IOdata/FileReader.cpp`、`src/IOdata/FileWriter.cpp`
  - **底层依赖**：C++26 标准库、系统 `libz` (`-lz`)。
- **构建方式**：
  ```bash
  # 使用 xmake 编译
  xmake build ParseLQCData
  # 或使用 CMake 编译
  cmake --build build_tree --target ParseLQCData
  ```

---

### 2. `bin/unit_tests` —— 介子关联函数与底层算法单元测试程序

- **定位**：底层 I/O、数值解析、Folding 折叠与统计重采样算法的独立验证套件。
- **入口源码**：[tests/test_meson_pipeline.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/tests/test_meson_pipeline.cpp)
- **依赖组件**：`core`, `meson_analysis`, `statistics`, `iodata`
- **验证功能**：
  1. `docs/physics_setup.json` 动态加载与各介子信道空间分量正确映射（如 AV 包含 6 个空间极化分量）；
  2. `FastParser` 零拷贝双精度浮点与整型快速解析精度；
  3. `Jackknife` 统计重采样与误差传播算法的绝对数值自洽性。
- **构建与运行**：
  ```bash
  # xmake 一键编译并运行
  xmake run unit_tests
  ```

---

### 3. `bin/test_condensate` —— 手征凝聚流水线单元测试程序

- **定位**：手征凝聚 XML 抽取、残余质量扣除与 Jackknife 重整化流水线测试程序。
- **入口源码**：[tests/test_condensate_pipeline.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/tests/test_condensate_pipeline.cpp)
- **依赖组件**：`core`, `condensate_analysis`, `statistics`, `iodata`
- **验证功能**：
  1. 验证手征凝聚 8 组温度/耦合常数参数的动态配置加载；
  2. 验证 XML 语法树高速解析器与残余质量扣除管道。
- **构建与运行**：
  ```bash
  # xmake 一键编译并运行
  xmake run test_condensate
  ```

---

## 二、 编译产出的 C-ABI 共享动态库 (位于 `lib/`)

共享动态库专门供 Python 胶水层（`ctypes`）动态加载，**严格输出至 [lib/](file:///Users/junxiongnie/code/algo/ParseLQCData/lib/) 目录，绝不污染 `build/`**：

| 动态库文件名 | 接口桥接文件 | 内部核心依赖 | 导出的 C-ABI 核心接口 | 调用方与使用场景 |
| :--- | :--- | :--- | :--- | :--- |
| **`lib/liblqcd_meson.so`** | `tools/Services/MesonService.cpp` | `core`, `meson_analysis`, `statistics`, `iodata` | `int run_meson_pipeline_c_api(...)` | 被 `tools/meson_orchestrator.py` 载入，执行多源/单源数据极速提取 |
| **`lib/liblqcd_condensate.so`** | `tools/Services/CondensateService.cpp` | `core`, `condensate_analysis`, `statistics`, `iodata` | `int run_condensate_pipeline_c_api(...)` | 被 `tools/condensate_orchestrator.py` 载入，执行 XML 抽取与手征凝聚计算 |
| **`lib/libparselqcdata.so`** | `tools/Services/*.cpp` | 全部计算组件与底层 I/O | 包含介子与手征凝聚全部 C-ABI | 全功能通用动态链接库 |

---

## 三、 `ParseLQCData` 已注册 Task 清单与使用语法

主程序通过 [include/ParseLQCData/Registry.h](file:///Users/junxiongnie/code/algo/ParseLQCData/include/ParseLQCData/Registry.h) 注册表动态挂载以下 7 大核心计算 Task：

```bash
# 查看当前所有可用的 Task 列表及说明
./bin/ParseLQCData
```

### 1. 介子关联函数计算任务

#### ① `meson_multi`：单信道多源关联函数抽取
- **语法**：`./bin/ParseLQCData meson_multi [beta] [channel] [binsize]`
- **示例**：
  ```bash
  ./bin/ParseLQCData meson_multi 17 AV 4
  ```
- **功能**：抽取指定 $\beta=4.17$ 下多源轴矢量介子（AV）的关联函数，进行 Binsize=4 的 Jackknife 重采样并完成 Folding 对称折叠。

#### ② `meson_single`：单信道单源关联函数抽取
- **语法**：`./bin/ParseLQCData meson_single [beta] [channel] [binsize]`
- **示例**：
  ```bash
  ./bin/ParseLQCData meson_single 17 S 4
  ```
- **功能**：抽取指定 $\beta=4.17$ 下单源标量介子（S）的关联函数。

#### ③ `meson_all`：全信道批量测量
- **语法**：`./bin/ParseLQCData meson_all [beta] [multi|single]`
- **示例**：
  ```bash
  # 批量测量 beta=4.17 下全部 6 个信道的多源数据
  ./bin/ParseLQCData meson_all 17 multi
  # 批量测量 beta=4.17 下全部 6 个信道的单源数据
  ./bin/ParseLQCData meson_all 17 single
  ```

#### ④ `meson_multi_all` / `meson_single_all`：独立批量任务
- **语法**：`./bin/ParseLQCData meson_multi_all [beta]` 或 `./bin/ParseLQCData meson_single_all [beta]`

---

### 2. 手征凝聚计算任务

#### ① `condensate`：单数据集手征凝聚测量
- **语法**：`./bin/ParseLQCData condensate [dataset_path_or_beta] [ml] [ms] [mres] [zm]`
- **示例**：
  ```bash
  # 支持自动正则匹配参数 (自动推导 ml, ms, mres, zm)
  ./bin/ParseLQCData condensate L32T12beta4.17
  
  # 支持手动传入物理参数
  ./bin/ParseLQCData condensate test_condensate 0.001 0.0384 0.000339722 0.966247
  ```
- **功能**：扫描数据集中的所有构型 XML 文件，计算夸克单魔环、减去残余质量并计算重整化后的手征凝聚值。

#### ② `condensate_all`：全量批量扫描
- **语法**：`./bin/ParseLQCData condensate_all [readin_dir]`
- **示例**：
  ```bash
  ./bin/ParseLQCData condensate_all data/readin
  ```
- **功能**：自动递归扫描 `data/readin/` 目录下所有包含 `meas.*/PsibarPsi/*.xml` 的数据集，多线程并发全量抽取。

---

## 四、 科学计算最佳实践（告别 `run_all.sh`，按需精准调度）

> ⚠️ **为什么淘汰 `run_all.sh`？**
> `run_all.sh` 盲目地将所有温度点、所有介子信道（多源+单源）以及手征凝聚数十万份 XML 构型强行串行串联，全量运行耗时极长，严重拖慢调试与日常研究效率。
> **强烈推荐以下模块化、阶段化的最佳实践：**

### 最佳实践 1：针对性测量与调试单个物理点（C++ 原生首选）
在排查特定信道或特定温度点时，直接调用 C++ 原生主程序，毫秒至秒级极速响应：
```bash
# 步骤 1: 编译主程序 (增量构建仅需 0.3s)
xmake build ParseLQCData

# 步骤 2: 精准测量目标信道
./bin/ParseLQCData meson_multi 20 PS 4
```

### 最佳实践 2：介子端到端数据流水线（Python 胶水层调度）
需要进行有效质量求解、非线性贝叶斯拟合以及自动生成图谱时，运行专项脚本：
```bash
# 执行多源介子关联函数全量抽取与拟合
bash scripts/run_meson_multi.sh

# 执行单源介子关联函数全量抽取与拟合
bash scripts/run_meson_single.sh
```

### 最佳实践 3：手征凝聚端到端重整化图谱
需要从原始构型 XML 抽取并生成物理相变曲线时：
```bash
# 一键完成手征凝聚 XML 抽取、数据导出与物理图谱绘制
bash scripts/run_condensate.sh
```

### 最佳实践 4：纯出图与可视化呈现 (零重复计算)
当底层数据已经抽取并存在于 `output/` 目录下时，无需重新计算，直接调用画图脚本：
```bash
# 绘制所有介子信道对比图与有效质量谱
uv run python scripts/plot_meson.py

# 绘制手征凝聚与奇异夸克重整化相变曲线
uv run python scripts/plot_condensate.py
```

### 最佳实践 5：高精度数值回归与基准自检
在修改了 C++ 核心算法或重采样逻辑后，执行轻量级回归自检（不跑冗长流水线）：
```bash
# 运行 C++ 底层单元测试 (< 1 秒)
xmake run unit_tests
xmake run test_condensate

# 运行逐点高精回归测试 (与原始 ana 逐点比对，验证 IEEE-754 精度)
uv run python tests/compare_with_ana.py
```

---

## 五、 双构建系统常用命令总结

### 1. 使用 xmake 构建与运行 (推荐)
```bash
# 编译所有默认目标 (CLI、Task 组件 .o 与 lib/ 下的所有动态库)
xmake

# 运行 CLI 查看任务
xmake run ParseLQCData

# 运行单元测试
xmake run unit_tests
xmake run test_condensate
```

### 2. 使用 CMake 构建与运行
```bash
# 1. 导出最新 CMake 配置
xmake project -k cmake .

# 2. 配置与构建
cmake -B build_tree -S .
cmake --build build_tree --target ParseLQCData
cmake --build build_tree --target lqcd_meson lqcd_condensate

# 3. 运行主程序
./bin/ParseLQCData
```

### 3. 原生独立编译脚本 (`scripts/run_build.sh`)
如需在脱离 xmake 的环境下仅重新编译 C-ABI 共享动态库：
```bash
bash scripts/run_build.sh
```
动态库将严格生成并输出至 [lib/](file:///Users/junxiongnie/code/algo/ParseLQCData/lib/) 目录。
