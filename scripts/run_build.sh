#!/usr/bin/env bash

g++-mp-15 -std=c++26 -Iinclude -O3 -shared -fPIC tools/Services/MesonService.cpp src/IOdata/* src/Statistics/* build/MesonAnalysis/MesonExtractor.cpp src/core/MesonPipeline.cpp -o lib/liblqcd_condensate.so
g++-mp-15 -std=c++26 -Iinclude -O3 -shared -fPIC tools/Services/CondensateService.cpp src/IOdata/* src/Statistics/* build/CondensateAnalysis/CondensateExtractor.cpp -o lib/liblqcd_condensate.so
