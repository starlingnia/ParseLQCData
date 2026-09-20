{
    files = {
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/CondensateAnalysis/CondensateExtractor.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/core/MesonPipeline.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/IOdata/DirectoryScanner.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/IOdata/FileReader.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/IOdata/FileWriter.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/MesonAnalysis/MesonExtractor.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/src/Statistics/Resampling.cpp.o",
        "build/.objs/parselqcdata_shared/macosx/x86_64/release/tools/Services/lqcd_meson_service.cpp.o"
    },
    values = {
        "/Library/Developer/CommandLineTools/usr/bin/clang++",
        {
            "-shared",
            "-target",
            "x86_64-apple-macos",
            "-isysroot",
            "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
            "-lz",
            "-fPIC",
            "-Wl,-x",
            "-Wl,-dead_strip"
        }
    }
}