  新工程已在 ParseLQCData 初始化创建：




  所有测试与执行程序已经全面梳理并经过真实测试验证完毕，运行均已通过。以下是项目的执行命令与测试清单：
  ──────
  ### 一、 核心 C++ 可执行程序（基于 apps/ 与 build/ 注册表）

  主程序入口为 main.cpp，通过注册表自动挂载 meson_task.cpp 与 condensate_task.cpp 实现的任务。

  #### 1. 编译原生 CLI 主程序

    g++-mp-15 -std=c++26 -Iinclude -O3 \
        apps/main.cpp build/meson_task.cpp build/condensate_task.cpp \
        src/IOdata/* src/Statistics/* \
        build/MesonAnalysis/MesonExtractor.cpp src/core/MesonPipeline.cpp \
        build/CondensateAnalysis/CondensateExtractor.cpp src/core/CondensatePipeline.cpp \
        -o bin/ParseLQCData

  #### 2. 运行 C++ 任务

    # 查看所有已注册的任务列表
    ./bin/ParseLQCData

    # [手征凝聚测量] 单数据集测量 (支持自动正则解析 L32T12_beta4.17... / L48T16beta4.13... 等格式)
    ./bin/ParseLQCData condensate L32T12beta4.17
    # 自定义物理参数测量: condensate [路径或名称] [ml] [ms] [mres] [zm]
    ./bin/ParseLQCData condensate test_condensate 0.001 0.0384 0.000339722 0.966247

    # [手征凝聚批量测量] 扫描 data/readin/ 下所有包含 meas.* 构型的数据集
    ./bin/ParseLQCData condensate_all data/readin

    # [介子多源关联函数测量] 语法: meson_multi [beta] [channel] [binsize]
    ./bin/ParseLQCData meson_multi 17 AV 4

    # [介子单源关联函数测量] 语法: meson_single [beta] [channel] [binsize]
    ./bin/ParseLQCData meson_single 17 S 4

    # [介子批量测量] 一键测量指定 beta 下的所有 6 种介子信道
    ./bin/ParseLQCData meson_all 17 multi
  ──────
  ### 二、 C++ 单元测试（位于 tests/）

  #### 1. 手征凝聚 C++ 核心流水线测试 (test_condensate_pipeline.cpp)

  测试 XML 数据解析、残余质量扣除、Jackknife 误差与重整化流水线：

    # 编译
    g++-mp-15 -std=c++26 -Iinclude -O3 tests/test_condensate_pipeline.cpp \
        src/IOdata/* src/Statistics/* \
        build/CondensateAnalysis/CondensateExtractor.cpp src/core/CondensatePipeline.cpp \
        -o bin/test_condensate

    # 运行 (已验证全部通过)
    ./bin/test_condensate

  #### 2. 介子关联函数 C++ 核心流水线测试 (test_meson_pipeline.cpp)

  测试 FastParser 零拷贝解析、重采样与 Folding 算法：

    # 编译
    g++-mp-15 -std=c++26 -Iinclude -O3 tests/test_meson_pipeline.cpp \
        src/IOdata/* src/Statistics/* \
        build/MesonAnalysis/MesonExtractor.cpp src/core/MesonPipeline.cpp \
        -o bin/test_meson

    # 运行 (已验证全部通过)
    ./bin/test_meson
  ──────
  ### 三、 Python 与物理精度回归测试套件（位于 tests/）

  运行综合测试套件（已包含数值精度比对、多线程并发压力以及自相关饱和度测试）：

    bash tests/run_tests.sh

  单独运行各测试项：

  1. 严格数值精度回归测试（逐点比对原 ana 结果，IEEE-754 精度）：
    uv run python tests/compare_with_ana.py

  2. 多线程并发扩展性与内存零泄露测试：
    uv run python tests/benchmark_stress_test.py

  3. Binsize 自相关分析与误差饱和度测试：
    uv run python tests/test_binsize_autocorr.py

  4. 全基准性能测试：
    uv run python tests/run_all_benchmarks.py

  ──────
  ### 四、 全流程自动化运行脚本（位于 scripts/）

   脚本文件            │ 功能说明                                               │ 运行命令
  ─────────────────────┼────────────────────────────────────────────────────────┼──────────────────────────────────
   run_all.sh          │ 一键端到端运行：依序执行介子与手征凝聚全量流水线及出图 │ bash scripts/run_all.sh
   run_meson.sh        │ 介子分析总控：多源/单源数据抽取、有效质量求解、拟合与  │ bash scripts/run_meson.sh
                       │ 出图                                                   │
   run_meson_multi.sh  │ 执行多源介子关联函数全量抽取与拟合                     │ bash scripts/run_meson_multi.sh
   run_meson_single.sh │ 执行单源介子关联函数全量抽取与拟合                     │ bash scripts/run_meson_single.sh
   run_condensate.sh   │ 手征凝聚总控：构型 XML 抽取、数据导出与物理图谱生成    │ bash scripts/run_condensate.sh
   run_build.sh        │ 重新编译 C++ C-ABI 共享动态库 (lib/)                   │ bash scripts/run_build.sh
  ──────
  ### 五、 底层共享动态库编译（供 Python ctypes 调用）

    # 编译介子动态库
    g++-mp-15 -std=c++26 -Iinclude -O3 -shared -fPIC \
        tools/Services/MesonService.cpp src/IOdata/* src/Statistics/* \
        build/MesonAnalysis/MesonExtractor.cpp src/core/MesonPipeline.cpp \
        -o lib/liblqcd_meson.so

    # 编译手征凝聚动态库
    g++-mp-15 -std=c++26 -Iinclude -O3 -shared -fPIC \
        tools/Services/CondensateService.cpp src/IOdata/* src/Statistics/* \
        build/CondensateAnalysis/CondensateExtractor.cpp src/core/CondensatePipeline.cpp \
        -o lib/liblqcd_condensate.so

    
    
    
    
    
    
  基础 I/O 通用工具已保存于  与 ：

  • FileReader.h / FileReader.cpp：极速批量文件入内存与零拷贝按行切片。
  • DirectoryScanner.h / DirectoryScanner.cpp：自然排序（与 Python natsort 严格等价）。
  • FastParser.h：零分配字符切片转双精度浮点与坐标解析。
  • FileWriter.h / FileWriter.cpp：高精度 CSV 格式化落盘。
  
    ### 1. 编译原生 CLI 可执行文件 (bin/ParseLQCData)

    mkdir -p bin

    gc -O3 -Wall \
      apps/main.cpp \
      tests/meson_task.cpp \
      tests/demo_task.cpp \
      src/IOdata/*.cpp \
      src/Statistics/*.cpp \
      src/MesonAnalysis/*.cpp \
      src/CondensateAnalysis/*.cpp \
      src/core/*.cpp \
      -o bin/ParseLQCData

  ### 2. 编译供 Python 调用的共享动态库 (build/libparselqcdata.dylib)

    mkdir -p build

    gc -O3 -Wall -fPIC -shared \
      tools/Services/*.cpp \
      src/IOdata/*.cpp \
      src/Statistics/*.cpp \
      src/MesonAnalysis/*.cpp \
      src/CondensateAnalysis/*.cpp \
      src/core/*.cpp \
      -o build/libparselqcdata.dylib
  
  
  
