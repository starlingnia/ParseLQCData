  新工程已在 ParseLQCData 初始化创建：

    ParseLQCData/
    ├── include/
    │   └── ParseLQCData/
    │       ├── Core/LQCDataTypes.h           # 基础物理信道配置与结果结构体定义
    │       ├── MesonAnalysis/                # 强子关联函数核心算子头文件
    │       │   ├── MesonExtractor.h          # 零拷贝单文件多源移位平均提取
    │       │   └── MesonPipeline.h           # 多线程并发批处理流水线
    │       └── Statistics/Resampling.h       # 统计重采样算子 (Binning / Jackknife / Folding)
    ├── src/
    │   ├── Statistics/Resampling.cpp         # 纯算术高内聚统计算子实现
    │   ├── MesonAnalysis/
    │   │   ├── MesonExtractor.cpp            # 状态机文本扫描与移位平移
    │   │   └── MesonPipeline.cpp             # 多核多线程并发调度流水线
    │   └── parselqcdata/
    │       ├── __init__.py
    │       └── meson_pipeline.py             # Python 业务层物理信道封装 (AV, S, PS, Vec, Xt, Tt)
    ├── tools/
    │   ├── Services/
    │   │   └── lqcd_meson_service.cpp        # 导出 C ABI (run_meson_pipeline_c_api) 供 ctypes 绑定
    │   └── lqcd_orchestrator.py              # 高层任务编排器 (加载 .dylib/.so，零拷贝内存对接)
    ├── apps/main.cpp                         # Pitchfork 原生注册中心入口
    ├── tests/
    │   ├── meson_task.cpp                    # 注册到 TaskRegistry 的 C++ CLI 任务 (meson_av)
    │   └── verify_against_ground_truth.py    # 端到端精度与性能基准比对测试脚本
    ├── output/b4.17/                         # 运算结果与中间数据存放目录 (隔离原 ana 项目)
    ├── xmake.lua                             # 同时编译 bin/ParseLQCData 与 bin/libparselqcdata.dylib
    └── pyproject.toml                        # uv 管理的 Python 依赖 (numpy, polars)
    
    
  基础 I/O 通用工具已保存于  与 ：

  • FileReader.h / FileReader.cpp：极速批量文件入内存与零拷贝按行切片。
  • DirectoryScanner.h / DirectoryScanner.cpp：自然排序（与 Python natsort 严格等价）。
  • FastParser.h：零分配字符切片转双精度浮点与坐标解析。
  • FileWriter.h / FileWriter.cpp：高精度 CSV 格式化落盘。
