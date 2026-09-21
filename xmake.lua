add_rules("mode.debug", "mode.release")

-- ==============================================================================
-- 独立原子组件 (Modular Components)
-- ==============================================================================

-- 组件 1: 文本 I/O 与零拷贝解析器
target("iodata")
    set_kind("static")
    set_languages("c++26")
    add_includedirs("include")
    add_files("src/IOdata/*.cpp")

-- 组件 2: 统计重采样与 Folding
target("statistics")
    set_kind("static")
    set_languages("c++26")
    add_includedirs("include")
    add_files("src/Statistics/*.cpp")

-- 组件 3: 介子关联函数特征抽取
target("meson_analysis")
    set_kind("static")
    set_languages("c++26")
    add_includedirs("include")
    add_files("src/MesonAnalysis/*.cpp")
    add_deps("iodata")

-- 组件 4: 手征凝聚 XML 抽取与重整化
target("condensate_analysis")
    set_kind("static")
    set_languages("c++26")
    add_includedirs("include")
    add_files("src/CondensateAnalysis/*.cpp")
    add_deps("iodata")

-- 组件 5: 底层多线程并行调度流水线核心
target("core")
    set_kind("static")
    set_languages("c++26")
    add_includedirs("include")
    add_files("src/core/*.cpp")
    add_deps("iodata", "statistics", "meson_analysis")

-- ==============================================================================
-- 产物目标 (Artifact Targets)
-- ==============================================================================

-- 组件 6: 任务组件 (单独只编译出 .o 目标文件，通过注册表挂载供 CLI 主程序调用)
target("task_components")
    set_kind("object")
    set_languages("c++26")
    add_includedirs("include")
    add_files("build/*_task.cpp")
    add_deps("core", "condensate_analysis", "meson_analysis", "statistics", "iodata")

-- 目标 A1: 介子分析专用动态链接库 (输出至 lib/，严禁污染 build/)
target("lqcd_meson")
    set_kind("shared")
    set_languages("c++26")
    set_basename("lqcd_meson")
    add_includedirs("include")
    add_files("tools/Services/MesonService.cpp")
    add_deps("core", "meson_analysis", "statistics", "iodata")
    set_targetdir("lib")

-- 目标 A2: 手征凝聚专用动态链接库 (输出至 lib/，严禁污染 build/)
target("lqcd_condensate")
    set_kind("shared")
    set_languages("c++26")
    set_basename("lqcd_condensate")
    add_includedirs("include")
    add_files("tools/Services/CondensateService.cpp")
    add_deps("core", "condensate_analysis", "statistics", "iodata")
    set_targetdir("lib")

-- 目标 A3: 全功能统一 C ABI 共享动态库 (供 Python ctypes 统一调度，输出至 lib/)
target("parselqcdata_shared")
    set_kind("shared")
    set_languages("c++26")
    set_basename("parselqcdata")
    add_includedirs("include")
    add_files("tools/Services/*.cpp")
    add_deps("core", "condensate_analysis", "meson_analysis", "statistics", "iodata")
    set_targetdir("lib")

-- 目标 B: 原生 CLI 驱动主程序 (包含 apps/ 下驱动主程序与 build/ 下注册的所有业务 Task)
target("ParseLQCData")
    set_kind("binary")
    set_languages("c++26")
    add_includedirs("include")
    add_files("apps/*.cpp")
    add_deps("task_components", "core", "condensate_analysis", "meson_analysis", "statistics", "iodata")
    set_targetdir("bin")

-- ==============================================================================
-- 测试目标 (Tests Target - 仅在需要测试时单独编译与运行)
-- 默认构建不包含此目标，需使用 `xmake build unit_tests` 或 `xmake run unit_tests`
-- ==============================================================================
target("unit_tests")
    set_default(false)
    set_kind("binary")
    set_languages("c++26")
    add_includedirs("include")
    add_files("tests/test_meson_pipeline.cpp")
    add_deps("core", "condensate_analysis", "meson_analysis", "statistics", "iodata")
    set_targetdir("bin")

target("test_condensate")
    set_default(false)
    set_kind("binary")
    set_languages("c++26")
    add_includedirs("include")
    add_files("tests/test_condensate_pipeline.cpp")
    add_deps("core", "condensate_analysis", "statistics", "iodata")
    set_targetdir("bin")
