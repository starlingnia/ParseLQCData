#!/usr/bin/env python3
"""
scripts/tasks/lcp/04_plot_lcp.py
--------------------------------------------------------------------------------
LCP (常物理线) 子任务 4: 手征凝聚与手征磁化率热力学曲线绘制
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 读取由子任务 1 导出的手征凝聚温度扫描数据 (results_rm_beta.txt, results_beta.txt)。
2. 读取由子任务 2~3 导出的手征磁化率 LCP 规范产物 (results_susceptibility.txt)。
3. 生成顶级期刊标准的高清出版图表:
   - 手征凝聚随温度变化与残余质量散度扣除效果图 (plot_condensate.py)
   - 手征磁化率随温度演化与赝临界相变峰值定位图 (plot_susceptibility.py)
4. 输出至 docs/figures/ 供分析报告与论文直接引用。

【底层调用的 C++ 功能与关联】:
- 本任务作为 LCP 项目可视化呈现层，直接呈现 C++ 核心库提取与重整化的热力学热点数据。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.plot_condensate import main as plot_condensate_main
from scripts.plot_susceptibility import main as plot_susceptibility_main


def main() -> None:
    parser = argparse.ArgumentParser(description="LCP 综合热力学物理曲线绘制小脚本")
    parser.parse_args()

    print("[LCP-TASK-04] 正在生成手征凝聚全景图表 (plot_condensate.py)...")
    plot_condensate_main()

    print("[LCP-TASK-04] 正在生成手征磁化率温度曲线 (plot_susceptibility.py)...")
    plot_susceptibility_main()

    print("[LCP-TASK-04] LCP 综合图表绘制完成！产物位于 docs/figures/")


if __name__ == "__main__":
    main()
