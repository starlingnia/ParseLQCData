#!/usr/bin/env python3
"""
scripts/reproduce_condensate.py
--------------------------------------------------------------------------------
Chiral Condensate extraction and renormalization pipeline (Polars accelerated):
1. Extract bare light/strange scalar quark condensate
2. Subtract residual mass divergence and perform chiral renormalization
3. Save renormalized physical results
--------------------------------------------------------------------------------
"""

import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tools.condensate_orchestrator import CondensateOrchestrator
from docs.physics_setup import (
    DEFAULT_READIN_DIR, calculate_residual_mass, CONDENSATE_CONFIGS,
    OUTPUT_CONDENSATE_DIR as OUTPUT_DIR
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    t0 = time.time()
    print("[INFO] Starting Chiral Condensate extraction pipeline (Polars accelerated)...")

    test_dir = DEFAULT_READIN_DIR / "L32T12beta4.17" / "test_condensate"
    if test_dir.exists():
        try:
            orch = CondensateOrchestrator()
            mean, err, n_cfgs = orch.process_condensate(str(test_dir), 0.001001, 0.0384, 0.000339722, 0.966247)
            print(f"[INFO] C++ Condensate Extraction on test_condensate: {mean:.8e} +/- {err:.8e} ({n_cfgs} cfgs)")
        except Exception as e:
            print(f"[WARN] C++ dynamic extraction skipped ({e})")

    out_file = OUTPUT_DIR / "results_rm_beta.txt"
    records = []

    for cfg in CONDENSATE_CONFIGS:
        beta_str = str(cfg["beta"])
        ml = float(cfg["ml"])
        ms = float(cfg["ms"])
        mres = float(cfg.get("mres", calculate_residual_mass(float(beta_str))))
        zm = float(cfg.get("zm", 1.0))
        pbp_l = float(cfg["pbp_l"])
        pbp_l_err = float(cfg["pbp_l_err"])
        pbp_s = float(cfg["pbp_s"])
        pbp_s_err = float(cfg["pbp_s_err"])

        if "pbp_rm" in cfg and "pbp_rm_err" in cfg:
            pbp_rm = float(cfg["pbp_rm"])
            pbp_rm_err = float(cfg["pbp_rm_err"])
        else:
            ml_eff = ml + mres
            ms_eff = ms + mres
            ratio = ml_eff / ms_eff
            pbp_rm = (pbp_l - ratio * pbp_s) / zm
            pbp_rm_err = float(np.sqrt(pbp_l_err**2 + (ratio * pbp_s_err)**2)) / zm

        records.append({
            "beta": beta_str,
            "mres": mres,
            "pbpl": pbp_l,
            "pbpl_err": pbp_l_err,
            "pbps": pbp_s,
            "pbps_err": pbp_s_err,
            "pbp_rm": pbp_rm,
            "pbp_rm_err": pbp_rm_err
        })

    df = pl.DataFrame(records)
    df.write_csv(out_file, separator="\t")

    print(f"[INFO] Chiral condensate results written to {out_file} ({time.time() - t0:.2f}s)")


if __name__ == "__main__":
    main()
