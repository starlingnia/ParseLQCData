{
    depfiles = "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/Statistics/__cpp_Resampling.cpp.cpp:   src/Statistics/Resampling.cpp include/Statistics/Resampling.h\
",
    files = {
        "src/Statistics/Resampling.cpp"
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
            "-fPIC",
            "-O3",
            "-std=c++26",
            "-Iinclude",
            "-DNDEBUG"
        }
    }
}