#!/usr/bin/env python3
"""
scripts/reproduce_meson.py
--------------------------------------------------------------------------------
介子流水线标准入口 (转发至统一的 run_meson.py)
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.run_meson import run_meson_pipeline, solve_effective_mass_ratio as solve_effective_mass_one_point_exact

if __name__ == "__main__":
    arg = sys.argv[1].lower() if len(sys.argv) > 1 else "all"
    if arg == "single":
        sources = [True]
    elif arg == "multi":
        sources = [False]
    else:
        sources = [False, True]
    run_meson_pipeline(sources=sources)