#!/usr/bin/env python3
"""
scripts/quick_plot_all_condensates.py
--------------------------------------------------------------------------------
全系综手征凝聚总览绘图入口 (转发至统一的 plot_condensate.py)
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.plot_condensate import main

if __name__ == "__main__":
    main()
