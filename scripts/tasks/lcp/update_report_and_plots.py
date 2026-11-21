#!/usr/bin/env python3
"""
scripts/tasks/lcp/update_report_and_plots.py
--------------------------------------------------------------------------------
常物理线 (LCP) 图表更新与报告编译调度脚本:
依次调度执行独立的绘图模块与报告编译模块:
1. scripts/tasks/lcp/plot_lcp_figures.py (独立绘制所有出版级图表)
2. scripts/tasks/lcp/compile_report.py (更新报告表格并编译 PDF)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tasks.lcp.plot_lcp_figures import main as run_plotting
from scripts.tasks.lcp.compile_report import main as run_compilation


def main() -> None:
    print("==========================================================================")
    print(">>> 启动 LCP 图表更新与报告编译工作流")
    print("==========================================================================")
    run_plotting()
    print("")
    run_compilation()
    print("\n>>> [ALL COMPLETE] 全套图件与报告更新顺利完成！")


if __name__ == "__main__":
    main()
