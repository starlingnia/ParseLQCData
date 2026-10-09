#!/usr/bin/env python3
"""
scripts/tasks/ccor/03_compute_delta_m.py
--------------------------------------------------------------------------------
CCOR 对称性破缺子任务 3: 4 组手征/自旋手征对称性信道质量差 ΔM 计算
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 在 Jackknife 样本层级对拟合出的介子质量逐样本作差:
   - ΔM(V - A)   : 手征对称性破缺序参量
   - ΔM(T - Xt)  : 自旋手征对称性信号
   - ΔM(S - PS)  : U(1)_A 轴向反常破缺指示
   - ΔM(Xt - A)  : 高阶自旋对称性
2. 使用 Jackknife 统计给出质量差及其统计误差棒。
3. 对照 ana/dat/ccor 基准数据进行高精度误差检验。

【底层调用的 C++ 功能】:
- C++ 统计核心: target("statistics") (src/Statistics/Resampling.cpp)
  负责两样本关联协方差与逐样本作差统计推断
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.reproduce_ccor_flow import main as ccor_main


def main() -> None:
    parser = argparse.ArgumentParser(description="CCOR 质量差与对称性破缺检验小脚本")
    parser.parse_args()
    print("[CCOR-TASK-03] 正在调度 reproduce_ccor_flow 计算质量差与基准对比...")
    ccor_main()
    print("[CCOR-TASK-03] 对称性信道质量差计算完成！")


if __name__ == "__main__":
    main()
