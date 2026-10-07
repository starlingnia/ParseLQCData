#!/usr/bin/env python3
"""
scripts/reproduce_meson_multi.py
--------------------------------------------------------------------------------
多源介子流水线入口 (转发至统一的 run_meson.py --source multi)
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.run_meson import run_meson_pipeline

if __name__ == "__main__":
    run_meson_pipeline(sources=[False])