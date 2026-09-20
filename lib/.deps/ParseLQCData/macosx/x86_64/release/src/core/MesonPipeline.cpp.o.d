{
    depfiles = "build/.objs/ParseLQCData/macosx/x86_64/release/src/core/__cpp_MesonPipeline.cpp.cpp:   src/core/MesonPipeline.cpp include/ParseLQCData/Core/MesonPipeline.h   include/ParseLQCData/Core/LQCDataTypes.h   include/MesonAnalysis/MesonExtractor.h include/Statistics/Resampling.h   include/IOdata/DirectoryScanner.h include/IOdata/FileReader.h\
",
    files = {
        "src/core/MesonPipeline.cpp"
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