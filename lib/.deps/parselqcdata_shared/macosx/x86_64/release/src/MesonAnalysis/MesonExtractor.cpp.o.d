{
    depfiles = "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/MesonAnalysis/__cpp_MesonExtractor.cpp.cpp:   src/MesonAnalysis/MesonExtractor.cpp   include/MesonAnalysis/MesonExtractor.h include/IOdata/FastParser.h   include/IOdata/FileReader.h\
",
    files = {
        "src/MesonAnalysis/MesonExtractor.cpp"
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