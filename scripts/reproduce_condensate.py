#!/usr/bin/env python3
"""
scripts/reproduce_condensate.py
--------------------------------------------------------------------------------
手征凝聚分析入口 (转发至统一的 run_condensate.py)
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.run_condensate import main

if __name__ == "__main__":
    main()
