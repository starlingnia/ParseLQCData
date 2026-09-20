{
    depfiles = "build/.objs/ParseLQCData/macosx/x86_64/release/src/MesonAnalysis/__cpp_MesonPipeline.cpp.cpp:   src/MesonAnalysis/MesonPipeline.cpp   include/ParseLQCData/MesonAnalysis/MesonPipeline.h   include/ParseLQCData/Core/LQCDataTypes.h   include/ParseLQCData/MesonAnalysis/MesonExtractor.h   include/ParseLQCData/Statistics/Resampling.h   ../leetcode/FuncSolv/include/IOdata/DirectoryScanner.h   ../leetcode/FuncSolv/include/IOdata/FileReader.h\
",
    files = {
        "src/MesonAnalysis/MesonPipeline.cpp"
    },
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
            "-I../leetcode/FuncSolv/include",
            "-DNDEBUG"
        }
    },
    depfiles_format = "gcc"
}