{
    depfiles = "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/IOdata/__cpp_FileReader.cpp.cpp:   src/IOdata/FileReader.cpp include/IOdata/FileReader.h\
",
    files = {
        "src/IOdata/FileReader.cpp"
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