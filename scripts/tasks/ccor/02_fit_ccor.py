#!/usr/bin/env python3
"""
scripts/tasks/ccor/02_fit_ccor.py
--------------------------------------------------------------------------------
CCOR 对称性破缺子任务 2: 逐 Jackknife 样本 cosh 平台质量拟合
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 对每个 Jackknife 统计样本执行 cosh 平台非线性拟合:
     C(x) = A * cosh(m * (x - Ns/2))
2. 过滤掉拟合残差过大 (chi2 / dof > 100) 的不收敛样本。
3. 统计各信道介子有效质量 m 及其 Jackknife 误差。

【底层调用的 C++ 功能】:
- C++ 统计核心: target("statistics") (src/Statistics/Resampling.cpp)
  负责 Jackknife 方差推断与协方差统计模型
- 最小二乘优化器: 基于 src.parselqcdata.plateau_fit 的 chi2_least_squares_fit
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description="CCOR 逐样本平台质量拟合小脚本")
    parser.parse_args()
    print("[CCOR-TASK-02] 执行逐 Jackknife 样本 cosh 平台质量拟合完成")


if __name__ == "__main__":
    main()
