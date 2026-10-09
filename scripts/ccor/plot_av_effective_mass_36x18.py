#!/usr/bin/env python3
"""Run the AV effective-mass sample for the 36x18, beta=4.17 ensemble."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import plot_av_effective_mass_48x18 as sample


sample.ENSEMBLE = "36x18_b4.17_ms0.040m0.0020"
sample.NS = 36
sample.NT = 18
sample.HALF = 18
sample.A_INV_MEV = 2452.96
sample.ML_TAG = "0.0020"
sample.INPUT_DIR = PROJECT_ROOT / "data" / "readin" / sample.ENSEMBLE / "Output"
sample.SAMPLE_OUTPUT_DIR = (
    PROJECT_ROOT
    / "output"
    / "meson_scan"
    / "b4.17"
    / "ccor"
    / "36x18"
)
sample.OUTPUT_PATH = sample.SAMPLE_OUTPUT_DIR / "meff_AV.png"


if __name__ == "__main__":
    raise SystemExit(sample.main())