{
    depfiles = "build/.objs/ParseLQCData/macosx/x86_64/release/src/CondensateAnalysis/__cpp_CondensateExtractor.cpp.cpp:   src/CondensateAnalysis/CondensateExtractor.cpp   include/CondensateAnalysis/CondensateExtractor.h   include/IOdata/FastParser.h include/IOdata/FileReader.h   include/IOdata/DirectoryScanner.h\
",
    files = {
        "src/CondensateAnalysis/CondensateExtractor.cpp"
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