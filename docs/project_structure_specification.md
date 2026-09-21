# ParseLQCData 项目架构与目录结构设计规范 (Pitchfork & Hybrid Pipeline)

本文档系统性阐述 **ParseLQCData** 项目的工程组织架构、目录职责边界、C++ 编码准则以及双构建系统（xmake & CMake）的协同规范。

---

## 一、 总体架构与设计哲学

本项目采用 **“底层模块化 C++ (GSL 标准) + 上层敏捷 Python 流水线”** 的混合高性能科学计算架构：

```
+-------------------------------------------------------------------------+
|                  Python 科学计算胶水层 (Scripts & Orchestrators)          |
|    - 顶层业务流水线调度 (scripts/*.py, tools/*_orchestrator.py)           |
|    - 科学参数组织、数据流编排、贝叶斯/非线性拟合、自动生成图谱 (docs/figures)  |
+-------------------------------------------------------------------------+
                                    │ ctypes C-ABI 动态载入
                                    ▼
+-------------------------------------------------------------------------+
|                C-ABI 动态链接库层 (lib/*.so / lib/*.dylib)               |
|    - 纯 C 接口封装 (tools/Services/*.cpp)                                |
|    - 零拷贝数据传递、跨语言指针与数组映射                                      |
+-------------------------------------------------------------------------+
                                    │ 静态链接 / 内部原子依赖
                                    ▼
+-------------------------------------------------------------------------+
|                高性能原子计算核心 (src/ & include/ & build/)              |
|    - 零拷贝文本解析器 (src/IOdata/)                                       |
|    - 统计重采样与 Folding (src/Statistics/)                             |
|    - 强子特征/XML 抽取 (src/MesonAnalysis/, src/CondensateAnalysis/)     |
|    - 多线程并行流水线 (src/core/)                                        |
|    - 注册表任务组件 (build/*_task.cpp -> .o 目标文件)                     |
+-------------------------------------------------------------------------+
                                    │ 链接
                                    ▼
+-------------------------------------------------------------------------+
|                原生 CLI 命令行驱动入口 (apps/main.cpp -> bin/)           |
|    - 基于注册表 (Registry) 自动挂载各组件的 Task 并在终端执行               |
+-------------------------------------------------------------------------+
```

### 核心设计原则：
1. **职责分离**：
   - **底层的并发度与吞吐交给 C++**：负责大规模磁盘文本的零拷贝极速 I/O、XML 构型批量解析、高负载统计重采样（Jackknife / Bootstrap）以及多线程并行流水线；
   - **高层抽象与流转交给 Python**：负责复杂配置读取、流程调度、非线性贝叶斯回归拟合、数据分析与科学绘图。
2. **现代 C++ 与 GSL (Guideline Support Library) 编程规范**：
   - **现代标准**：采用 C++26 标准，充分利用标准库现代设施；
   - **零裸指针与视图优先**：核心数据传递使用 `std::span`、`std::string_view`，消灭裸指针算术与潜在越界风险，实现跨层零内存拷贝；
   - **控制流防卫式设计**：严格推行**早期返回（Early Exit / Guard Clauses）**，杜绝多层深度的 `if-else` 分支嵌套；消除 `for` 循环体内复杂条件跳转，使异常与边界情况第一时间退回，保持核心逻辑线性清晰；
   - **严格的 RAII 与零内存泄漏**：资源（句柄、文件映射、动态内存）全面由标准容器及智能指针全生命周期托管。
3. **注册表模式 (Registry Pattern) 任务挂载**：
   - CLI 驱动主程序与业务测量任务彻底解耦。各业务任务独立封装为静态自动注册单元，在编译期仅作为 `.o` 目标文件存在，装载时自动注入注册表，主程序按名称进行派发。

---

## 二、 各目录结构的用处与内容规范

项目全面遵循现代 **Pitchfork Directory Standard (PDS)** 规范，各目录界定明确的职责边界：

| 目录名称 | 定位与用处 | 允许放置的内容 | 严禁放置的内容 |
| :--- | :--- | :--- | :--- |
| **`apps/`** | **原生可执行程序入口源码** | 包含 `main()` 函数的 C++ 驱动文件（如 `apps/main.cpp`） | 具体的原子业务计算逻辑、非入口业务 Task |
| **`build/`** | **任务组件源码、中间目标文件与构建目录** | ① 配合 `apps/` 编译的业务 Task 组件源码（`*_task.cpp`，只含注册函数不含 `main`）；<br>② 构建中间产物（`.o` 目标文件、静态库中间缓存、`.build_cache` 等） | **严禁输出动态链接库（.so/.dylib）**！禁止存放最终二进制可执行文件 |
| **`include/`** | **C++ 头文件根目录** | 所有对外与内部公开的 `.h` / `.hpp` 头文件 | 具体的 `.cpp` 实现文件、非源码资源 |
| **`src/`** | **核心 C++ 算法与流水线实现源码** | 核心算法与管道实现（`core/`, `IOdata/`, `Statistics/` 等），以及除 Python 接口以外的管道文件 | 包含 `main()` 的入口文件、纯 C-ABI 包装层代码 |
| **`scripts/`** | **顶层批处理与端到端运行脚本** | 调度整体分析流程的 Bash 与 Python 脚本（`run_all.sh`, `run_meson.sh` 等） | 临时测试数据、编译构建产物 |
| **`tools/`** | **服务中间层与外围辅助工具** | ① `tools/Services/`：C-ABI 动态链接库的接口桥接源码；<br>② `*_orchestrator.py`：Python ctypes 调度控制器 | 核心物理底层算法实现、输入源数据集 |
| **`opt/`** | **可选第三方依赖与优化模块** | 第三方优化配置、可选扩展支持包、外部工具链适配代码 | 项目本身的主线业务源码 |
| **`bin/`** | **编译产出的原生可执行程序** | 构建最终生成的独立可执行二进制（如 `bin/ParseLQCData`） | 动态链接库、中间编译目标文件、源代码 |
| **`lib/`** | **编译产出的 C-ABI 共享动态库** | 构建生成的共享动态链接库（`.so` / `.dylib`） | 最终可执行程序、构建临时缓存文件、源码 |
| **`tests/`** | **测试套件与数值精度回归代码** | C++ 单元测试（`test_*.cpp`）、Python 数值回归对比脚本（`compare_with_ana.py`）、压测脚本 | 生产环境业务脚本、未隔离的实验性改动 |
| **`data/`** | **物理测量原始输入数据** | 构型 XML 测量文件、强子关联函数输入源文件（如 `data/readin/`） | **严禁提交至 Git**（受 `.gitignore` 保护），严禁编译系统在此写入 |
| **`output/`** | **计算与拟合分析落盘产物** | 计算生成的 CSV 表格、有效质量求解结果、拟合汇总文件等 | 静态代码文件、编译库文件 |
| **`docs/`** | **项目工程文档与验证报告** | 架构规范、数学与物理推导文档、全量数值回归报告、压测报告、自动生成的物理图像 | 任何二进制构建产物 |

---

## 三、 重点目录详细解析

### 1. `apps/` 目录：可执行程序唯一入口
- **规范说明**：此目录下是**可以直接编译出最终执行程序的 C++ 源码**。
- **典型文件**：[apps/main.cpp](file:///Users/junxiongnie/code/algo/ParseLQCData/apps/main.cpp)
- **运行逻辑**：
  不硬编码任何具体物理计算逻辑。它在运行时访问 `registry::TaskRegistry::instance()`，列出所有已注册的任务名及描述，并根据命令行参数派发执行对应的任务处理器。

### 2. `build/` 目录：注册表任务组件与目标中间物（防污染红线）
- **规范说明**：
  包含两层含义：
  1. **注册表任务组件（Task Components）**：例如 `condensate_task.cpp`、`meson_task.cpp`。它们配合 `apps/main.cpp` 编译，但**单独只编译出 `.o` 目标文件（或 CMake Object Library）**，绝不包含独立的 `main()` 函数，不单独产生独立的可执行二进制。组件通过 `REGISTER_TASK` 宏在静态初始化时向全局注册表注入可执行函数。
  2. **构建系统中间目录**：存放编译器生成的依赖文件（`.deps`）、对象缓存（`.objs`）、静态链接中间文件。
- **🚨 核心隔离红线**：
  **动态库（`.so` / `.dylib`）必须统一输出到 `lib/` 目录下，严禁输出到 `build/` 下造成目录污染！**

### 3. `include/` 目录：两大核心头文件层级
头文件严格划分为两大板块，各司其职：
- **`include/__project_name__/` (即 `include/ParseLQCData/`)**：
  - **核心注册表引擎**：[include/ParseLQCData/Registry.h](file:///Users/junxiongnie/code/algo/ParseLQCData/include/ParseLQCData/Registry.h)。
    > **规范准则**：**任务注册表头文件必须位于 `include/__project_name__/Registry.h`**，严禁放入 `include/core/`。所有业务 Task 组件与 `apps/main.cpp` 入口统一通过 `#include <__project_name__/Registry.h>`（例如 `#include <ParseLQCData/Registry.h>`）引入注册宏与注册表单例。
  - **项目专属领域头文件**：与本项目强相关的核心计算模型、特征提取器与流水线定义头文件（如 `include/IOdata/`、`include/Statistics/`、`include/MesonAnalysis/`、`include/CondensateAnalysis/` 等）。
- **`include/core/`**：
  通用底层基础功能函数、通用工具类（如字符串/通用切片）、数学工具、并行控制抽象等与具体项目名称解耦的基础设施头文件。

### 4. `src/` 目录：核心 C++ 算法源码
存放核心领域算力实现，与 Python 交互逻辑无关：
- `src/core/`：多线程多信道调度流转（如 `MesonPipeline.cpp`, `CondensatePipeline.cpp`）；
- `src/IOdata/`：系统级零拷贝快速流处理；
- `src/Statistics/`：数值统计重采样算法实现；
- `src/MesonAnalysis/` & `src/CondensateAnalysis/`：物理核心算符抽取实现。

### 5. `tools/` 与 `opt/` 目录：C-ABI 服务桥接与外围工具
- **`tools/Services/`**：
  为 Python 提供纯 C-ABI 函数签名的桥接实现（例如 `MesonService.cpp` 与 `CondensateService.cpp`），内部使用 `extern "C"` 包装，并捕获所有 C++ 异常，转换为整型状态码，保证跨语言调用的绝对稳健。
- **`tools/` 根目录**：
  存放 Python 端的 Orchestrator 封装（`meson_orchestrator.py`, `condensate_orchestrator.py`），利用 `ctypes` 动态加载 `lib/` 下的共享库，向顶层业务脚本提供面向对象的 Python API。
- **`opt/`**：
  预留给后续可选性能调优工具、外部算法库配置。

### 6. `bin/` 与 `lib/` 目录：严密的产物输出边界
- **`bin/`**：仅存放编译生成的可直接启动的终端程序（`ParseLQCData`）；
- **`lib/`**：专职存放动态链接库（`liblqcd_meson.so`, `liblqcd_condensate.so`, `libparselqcdata.so` 等）。

---

## 四、 双构建系统配置与指令规范

本项目原生使用 **xmake** 构建系统，并支持一键导出标准跨平台的 **CMake** 构建配置：

### 1. 使用 xmake 构建与运行

```bash
# 1. 编译全部默认目标 (包含 C++ CLI、任务组件 .o 与 lib/ 下的所有动态库)
xmake

# 2. 仅编译并运行 CLI 主程序
xmake run ParseLQCData

# 3. 运行特定物理任务 (通过注册表)
xmake run ParseLQCData meson_multi 17 AV 4
xmake run ParseLQCData condensate L32T12beta4.17

# 4. 编译并执行单元测试
xmake run unit_tests
```

### 2. 通过 xmake 导出标准 CMake 配置

执行以下命令即可从 `xmake.lua` 生成标准、现代化的 `CMakeLists.txt`：

```bash
xmake project -k cmake .
```

在导出的 [CMakeLists.txt](file:///Users/junxiongnie/code/algo/ParseLQCData/CMakeLists.txt) 中已严格遵循防污染规范：
- 共享动态库配置：
  ```cmake
  set_target_properties(lqcd_meson PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
  set_target_properties(lqcd_condensate PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
  set_target_properties(parselqcdata_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
  ```
- 任务组件作为目标文件编译：
  ```cmake
  add_library(task_components OBJECT "")
  ```
- 原生 CLI 可执行程序输出至 `bin/`：
  ```cmake
  set_target_properties(ParseLQCData PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
  ```

### 3. 使用 CMake 构建

```bash
# 配置构建目录
cmake -B build_tree -S .

# 编译所有动态链接库 (产物自动进入 lib/)
cmake --build build_tree --target lqcd_meson lqcd_condensate parselqcdata_shared

# 编译 CLI 主程序 (产物自动进入 bin/)
cmake --build build_tree --target ParseLQCData

# 运行主程序
./bin/ParseLQCData
```
