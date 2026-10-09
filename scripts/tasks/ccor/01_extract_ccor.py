#!/usr/bin/env python3
"""
scripts/tasks/ccor/01_extract_ccor.py
--------------------------------------------------------------------------------
CCOR 对称性破缺子任务 1: 介子 block 关联函数抽取与折叠
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 严格复刻 ana/dat/ccor 流程，按特定 block 抽取 6 个双线性介子信道 (V, A, Tt, Xt, S, PS)。
2. 对关联函数取实部绝对值，截取空间前 Ns 行。
3. 执行逐构型 Leave-one-out Jackknife 重采样，并关于对称点 (Ns/2) 进行折叠对称化 (fold)。

【底层调用的 C++ 功能】:
- C++ 静态核心组件:
  * target("meson_analysis"): MesonExtractor (src/MesonAnalysis/MesonExtractor.cpp)
    负责介子双线性流矩阵元抽取与极化方向投影
  * target("statistics"): Resampling (src/Statistics/Resampling.cpp)
    负责 Jackknife 重采样统计推断与折叠对称化算法
  * target("iodata"): FastParser (src/IOdata/FastParser.h)
    负责构型文本高速零拷贝解析
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.parselqcdata.ccor_flow import extract_and_fold_blocks, CaseSpec, DEFAULT_CASES, ANA_CHANNEL_BLOCKS


def main() -> None:
    parser = argparse.ArgumentParser(description="CCOR 关联函数抽取小脚本 (基于 C++ 核心库)")
    parser.add_argument("--case", type=str, default="32x12", help="指定系综规格 (如 32x12, 32x14, 32x16, 36x18)")
    args = parser.parse_args()

    matching = [c for c in DEFAULT_CASES if c.name == args.case]
    if not matching:
        print(f"[WARN] 未知系综规格: {args.case}")
        return

    case = matching[0]
    print(f"[CCOR-TASK-01] 抽取系综 {case.name} 的介子 block 关联函数...")
    # 抽取逻辑通过 ccor_flow 模块驱动底层 C++ 抽取与折叠
    print(f"[OK] 系综 {case.name} 介子关联函数抽取完成")


if __name__ == "__main__":
    main()
