add_rules("mode.debug", "mode.release")

-- C++ 共享动态库目标：供 Python ctypes 高性能高并发调用
target("parselqcdata_shared")
    set_kind("shared")
    set_languages("c++26")
    set_basename("parselqcdata")
    add_includedirs("include", "../leetcode/FuncSolv/include")
    add_files(
        "../leetcode/FuncSolv/src/IOdata/*.cpp",
        "src/Statistics/*.cpp",
        "src/MesonAnalysis/*.cpp",
        "tools/Services/*.cpp"
    )
    set_targetdir("bin")

-- 原生可执行文件目标：支持任务注册与 CLI 驱动
target("ParseLQCData")
    set_kind("binary")
    set_languages("c++26")
    add_includedirs("include", "../leetcode/FuncSolv/include")
    add_files(
        "apps/main.cpp",
        "tests/*.cpp",
        "../leetcode/FuncSolv/src/IOdata/*.cpp",
        "src/Statistics/*.cpp",
        "src/MesonAnalysis/*.cpp"
    )
    set_targetdir("bin")
