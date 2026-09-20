{
    depfiles = "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/IOdata/__cpp_FileWriter.cpp.cpp:   src/IOdata/FileWriter.cpp include/IOdata/FileWriter.h\
",
    files = {
        "src/IOdata/FileWriter.cpp"
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