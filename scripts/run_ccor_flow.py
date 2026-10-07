#!/usr/bin/env python3
"""
scripts/run_ccor_flow.py
--------------------------------------------------------------------------------
复刻并对照 ana/dat/ccor 对称性破缺 (chiral / axial U(1) / spin-chiral) 质量差分析流程:
  - 抽取介子空间关联函数并折叠对称化 (fold C(x) <-> C(Ns-x))
  - 逐 Jackknife 样本 cosh 拟合介子有效质量
  - 计算 4 组手征/自旋手征对称性信道质量差 ΔM (V-A, T-X, S-PS, X-A) 随温度演化
  - 输出与 ana/dat/ccor 参考基准的检验对比与图表

用法:
  python scripts/run_ccor_flow.py                                # 三源 (reference/multisrc/singlesrc) 全跑
  python scripts/run_ccor_flow.py --datasets multisrc --no-gnuplot
  python scripts/run_ccor_flow.py --cases 32x12 32x14 32x16 36x18
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.reproduce_ccor_flow import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
