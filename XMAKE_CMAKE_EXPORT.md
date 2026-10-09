# xmake 导出 CMakeLists.txt 指南 (XMake to CMake Export Guide)

本文档专门记录与规范本项目通过 **xmake** 自动生成与导出 **`CMakeLists.txt`** 的操作方法、参数说明与维护工作流。

---

## 1. 核心设计原则

- **单一事实源 (Single Source of Truth)**：本项目的原生构建配置为根目录下的 [`xmake.lua`](xmake.lua)。
- **严禁手动编辑 `CMakeLists.txt`**：项目主目录下的 `CMakeLists.txt` 完全由 xmake 自动生成。任何 target、宏定义、编译标志或源文件的变更，均需在 `xmake.lua` 中调整，随后执行导出命令同步覆盖。
- **导出目的**：
  1. 支持 CLion、Visual Studio Code (CMake Tools) 等 IDE 的语义索引、代码补全与调试。
  2. 便于在未安装 xmake 的宿主环境或 CI/CD 机器上使用标准 CMake 进行编译与验证。

---

## 2. 导出命令 (Export Commands)

### 2.1 推荐标准导出命令

在项目根目录下，直接执行以下命令即可一键重新生成最新的 `CMakeLists.txt`：

```bash
xmake project -k cmake
```

> **提示**：`-k cmake` 与 `-k cmakelists` 完全等价。

### 2.2 高级与多配置导出选项

如果需要同时导出 Debug 与 Release 多配置模式，或指定目标架构：

```bash
# 同时导出 Release 与 Debug 模式
xmake project -k cmake -m "release,debug"

# 显式指定架构 (如 macOS Apple Silicon arm64 或 Linux x86_64)
xmake project -k cmake -m "release,debug" -a arm64

# 显式指定输出目录为当前项目根目录
xmake project -k cmake -P . -o .
```

---

## 3. 使用 CMake 编译验证 (CMake Build Workflow)

导出 `CMakeLists.txt` 后，可按标准 CMake 流程进行跨平台编译验证：

```bash
# 1. 创建并配置 CMake 构建目录 (推荐构建目录隔离)
cmake -B build_cmake -DCMAKE_BUILD_TYPE=Release

# 2. 并行编译全部目标
cmake --build build_cmake -j$(sysctl -n hw.logicalcpu || nproc)

# 3. 运行编译生成的二进制验证
./build_cmake/bin/ParseLQCData
```

---

## 4. 本项目 Target 映射关系表

导出的 `CMakeLists.txt` 完整保留了 `xmake.lua` 中定义的依赖拓扑结构：

| Target 名称 | 类型 (Kind) | 源码位置 | 产物输出路径 | 依赖项 (Dependencies) |
| :--- | :--- | :--- | :--- | :--- |
| **`iodata`** | 静态库 (`static`) | `src/IOdata/*.cpp` | 构建缓存 | 无 |
| **`statistics`** | 静态库 (`static`) | `src/Statistics/*.cpp` | 构建缓存 | 无 |
| **`meson_analysis`** | 静态库 (`static`) | `src/MesonAnalysis/*.cpp` | 构建缓存 | `iodata` |
| **`condensate_analysis`** | 静态库 (`static`) | `src/CondensateAnalysis/*.cpp` | 构建缓存 | `iodata` |
| **`core`** | 静态库 (`static`) | `src/core/*.cpp` | 构建缓存 | `iodata`, `statistics`, `meson_analysis` |
| **`task_components`** | 目标对象 (`object`) | `build/*_task.cpp` | 构建缓存 | 全部底层模块 |
| **`lqcd_meson`** | 动态库 (`shared`) | `tools/Services/MesonService.cpp` | `lib/` | `core`, `meson_analysis`, `statistics`, `iodata` |
| **`lqcd_condensate`** | 动态库 (`shared`) | `tools/Services/CondensateService.cpp` | `lib/` | `core`, `condensate_analysis`, `statistics`, `iodata` |
| **`parselqcdata_shared`** | 统一动态库 (`shared`) | `tools/Services/*.cpp` | `lib/` | 全部模块 (供 Python ctypes 调度) |
| **`ParseLQCData`** | CLI 可执行程序 (`binary`) | `apps/*.cpp` | `bin/` | `task_components` + 全部底层模块 |
| **`unit_tests`** | 测试程序 (`binary`) | `tests/test_meson_pipeline.cpp` | `bin/` | 全部核心模块 |
| **`test_condensate`** | 测试程序 (`binary`) | `tests/test_condensate_pipeline.cpp` | `bin/` | `core`, `condensate_analysis`, `statistics`, `iodata` |

---

## 5. 维护与更新工作流 (Maintenance Workflow)

当发生以下情况时，**必须**重新执行 `xmake project -k cmake`：
1. **添加或删除了 C++ 源文件**（例如新增了任务 `build/xxx_task.cpp` 或算法模块）；
2. **修改了头文件包含路径**（`include` 目录变动）；
3. **调整了编译选项或 C++26 语言标准设置**；
4. **新增了外部依赖包或修改了链接库**。

执行导出后，使用 `git diff CMakeLists.txt` 检查差异并一并提交至版本控制库中。
