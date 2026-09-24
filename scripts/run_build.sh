#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

if command -v xmake &>/dev/null; then
    echo "[BUILD] Building project with xmake..."
    xmake -r
else
    echo "[BUILD] xmake not found, compiling manually..."
    CXX=""
    for c in g++-mp-15 clang++ g++; do
        if command -v "$c" &>/dev/null; then
            CXX="$c"
            break
        fi
    done
    if [ -z "$CXX" ]; then
        echo "[ERROR] No C++ compiler found (tried g++-mp-15, clang++, g++)"
        exit 1
    fi
    echo "[BUILD] Using compiler: $CXX"
    mkdir -p lib
    "$CXX" -std=c++26 -Iinclude -O3 -shared -fPIC tools/Services/MesonService.cpp src/IOdata/*.cpp src/Statistics/*.cpp src/MesonAnalysis/MesonExtractor.cpp src/core/MesonPipeline.cpp -o lib/liblqcd_meson.dylib
    "$CXX" -std=c++26 -Iinclude -O3 -shared -fPIC tools/Services/CondensateService.cpp src/IOdata/*.cpp src/Statistics/*.cpp src/CondensateAnalysis/CondensateExtractor.cpp src/core/CondensatePipeline.cpp -o lib/liblqcd_condensate.dylib
    "$CXX" -std=c++26 -Iinclude -O3 -shared -fPIC tools/Services/*.cpp src/IOdata/*.cpp src/Statistics/*.cpp src/MesonAnalysis/MesonExtractor.cpp src/CondensateAnalysis/CondensateExtractor.cpp src/core/*.cpp -o lib/libparselqcdata.dylib
fi
echo "[BUILD] Build completed successfully."
