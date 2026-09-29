#!/usr/bin/env python3
"""临时端到端驱动：把 OUTPUT_ROOT 重定向到 /tmp，避免污染仓库输出，测各阶段墙钟耗时。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def run() -> None:
    import scripts.reproduce_meson as rm
    out = Path("/tmp/plqc_e2e_output")
    out.mkdir(parents=True, exist_ok=True)
    rm.OUTPUT_ROOT = out
    rm.main()


if __name__ == "__main__":
    run()
