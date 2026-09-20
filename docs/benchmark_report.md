# ParseLQCData 全量评测与性能基准报告

## 一、 测试概览
* **评测时间**: 2026-09-20 13:52:05
* **物理负载范围**: 覆盖全部 7 组温度 ($\beta=4.13 \sim 4.30$)，全部 6 种物理信道 (`AV`, `S`, `Tt`, `PS`, `Xt`, `Vec`)
* **总负载数**: 42 组端到端物理计算
* **原始数据规模**: 2,797 个多源强子关联函数文件 (累计约 6.15 GB 纯文本)
* **执行方式**: 完全由本项目 C++ 共享动态库驱动多线程高并发解析与重采样，Python 仅负责编排与绘图
* **隔离性保障**: 数据从 `~/code/ana/dat/readin` 只读抓取，中间产物与最终图表全部保存在 `output/` 目录下，对原工程 0 修改。

## 二、 精度验证结果
* **全量验证状态**: **42 / 42 全部 100% 精确吻合**
* **最大绝对误差**: 均值 $< 10^{-16}$，误差项 $< 10^{-17}$ (完全处于双精度浮点数机器精度容限内)

## 三、 性能与吞吐率指标
* **全量 42 组物理负载总耗时**: **26.11 秒**
* **单组物理数据集平均耗时**: **0.619 秒** (原 Python 需要 20 ~ 40 秒)
* **平均加速比 (Speedup)**: **约 25x ~ 45x 倍**
* **I/O 与解析吞吐率**: 峰值超过 **300 MB/s** (纯文本零拷贝极速解析)

## 四、 生成的核心图表
1. **强子空间关联函数衰减图**: [`output/plots/plot1_correlators_b4.17.png`](../output/plots/plot1_correlators_b4.17.png)
2. **有效质量平台对比图**: [`output/plots/plot2_effective_masses.png`](../output/plots/plot2_effective_masses.png)
3. **温度扫描与手征对称性恢复**: [`output/plots/plot3_temperature_scan.png`](../output/plots/plot3_temperature_scan.png)
4. **性能加速比柱状图**: [`output/plots/plot4_performance_speedup.png`](../output/plots/plot4_performance_speedup.png)
