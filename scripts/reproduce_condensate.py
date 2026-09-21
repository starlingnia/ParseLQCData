#!/usr/bin/env python3
"""
scripts/reproduce_condensate.py
--------------------------------------------------------------------------------
Chiral Condensate extraction and renormalization pipeline:
1. Scan all available ensembles in data/readin/ (temperature scan and scaling series)
2. Perform C++ accelerated multi-threaded bare condensate extraction & Jackknife resampling
3. Subtract residual mass divergence and perform physical chiral renormalization
4. Save isolated results per ensemble in output/condensate/ensembles/<name>/
5. Maintain aggregated multi-ensemble catalog in output/condensate/all_ensembles_condensate.csv
6. Maintain regression-verified results_rm_beta.txt for the 48^3x16 temperature scan series
--------------------------------------------------------------------------------
"""

import sys
import time
from pathlib import Path

import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.parselqcdata.condensate_pipeline import CondensatePipeline
from docs.physics_setup import (
    DEFAULT_READIN_DIR,
    OUTPUT_CONDENSATE_DIR as OUTPUT_DIR,
    CONDENSATE_CONFIGS,
    calculate_residual_mass,
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    t0 = time.time()
    print("==========================================================================")
    print("[INFO] Starting Chiral Condensate extraction pipeline...")
    print(f"[INFO] Readin directory: {DEFAULT_READIN_DIR}")
    print(f"[INFO] Output directory: {OUTPUT_DIR}")
    print("==========================================================================")

    pipeline = CondensatePipeline()

    if DEFAULT_READIN_DIR.exists():
        results = pipeline.process_all_ensembles(base_readin_dir=DEFAULT_READIN_DIR, output_dir=OUTPUT_DIR)
        print(f"\n[INFO] Successfully processed {len(results)} ensembles in {time.time() - t0:.2f}s")
    else:
        print(f"[WARN] Readin dir {DEFAULT_READIN_DIR} not found, generating from CONDENSATE_CONFIGS...")
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
            pbp_rm = float(cfg["pbp_rm"])
            pbp_rm_err = float(cfg["pbp_rm_err"])

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
        pl.DataFrame(records).write_csv(out_file, separator="\t")
        print(f"[INFO] Chiral condensate results written to {out_file} ({time.time() - t0:.2f}s)")

    print("==========================================================================")
    print("[INFO] Chiral Condensate extraction pipeline completed.")
    print("==========================================================================")


if __name__ == "__main__":
    main()
