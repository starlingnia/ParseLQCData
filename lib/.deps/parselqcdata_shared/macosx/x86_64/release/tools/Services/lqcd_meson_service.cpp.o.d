{
    depfiles = "build/.objs/parselqcdata_shared/macosx/x86_64/release/tools/Services/__cpp_lqcd_meson_service.cpp.cpp:   tools/Services/lqcd_meson_service.cpp   include/ParseLQCData/Core/MesonPipeline.h   include/ParseLQCData/Core/LQCDataTypes.h   include/CondensateAnalysis/CondensateExtractor.h\
",
    files = {
        "tools/Services/lqcd_meson_service.cpp"
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