#!/usr/bin/env python3
"""
scripts/tasks/meson/04_plot_meson.py
--------------------------------------------------------------------------------
介子测量子任务 4: 出版级强子物理图表绘制与能谱可视化
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 读取由子任务 1~3 导出的强子关联函数 (pickdata)、有效质量 (ratio_results) 与拟合平台 (simulateresult)。
2. 生成符合高能物理顶级期刊 (PRD / JHEP) 标准的可视化图表:
   - 单信道有效质量随时间演化图 (带 cosh 拟合平台带与 Jackknife 误差棒)
   - 多信道空间关联函数衰减对比图
   - 手征对称性破缺与恢复监测对比图 (例如 AV 与 Vec 对比, PS 与 S 对比)
3. 产物输出至 docs/figures/ 目录。

【底层调用的 C++ 功能与关联】:
- 本任务作为分析流水线终端呈现层，直接呈现 C++ 核心库 (target("meson_analysis") + target("core"))
  并发抽取、对称化及平台拟合出的所有物理观测量。
- 自动调用 scripts/plot_meson.py 主绘图引擎。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.plot_meson import main as plot_main


def main() -> None:
    parser = argparse.ArgumentParser(description="介子物理图表绘制小脚本")
    parser.parse_args()

    print("[MESON-TASK-04] 正在调用 plot_meson.py 生成介子全信道图谱...")
    plot_main()
    print("[MESON-TASK-04] 图表生成完毕！产物位于 docs/figures/")


if __name__ == "__main__":
    main()
