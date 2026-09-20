{
    depfiles = "build/.objs/ParseLQCData/macosx/x86_64/release/tests/__cpp_demo_task.cpp.cpp:   tests/demo_task.cpp include/ParseLQCData/Registry.h\
",
    files = {
        "tests/demo_task.cpp"
    },
    depfiles_format = "gcc",
    values = {
        "/Library/Developer/CommandLineTools/usr/bin/clang++",
        {
            "-Qunused-arguments",
            "-target",
            "x86_64-apple-macos",
            "-isysroot",
            "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
            "-fvisibility=hidden",
            "-fvisibility-inlines-hidden",
            "-O3",
            "-std=c++26",
            "-Iinclude",
            "-DNDEBUG"
        }
    }
}