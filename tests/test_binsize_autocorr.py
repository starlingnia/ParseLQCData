#!/usr/bin/env python3
"""
tests/test_binsize_autocorr.py
Minimalist test for optimal binsize in LQCD data:
1. Scans binsizes: [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20]
2. Tests auto-correlation elimination via lag-1 autocorrelation rho_1
3. Tests error estimate saturation: sigma(B) / sigma_plateau >= 95%
"""

import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tools.meson_orchestrator import MesonOrchestrator
from docs.physics_setup import CHANNEL_CONFIGS, DEFAULT_READIN_DIR


def compute_lag1_autocorr(series: np.ndarray) -> float:
    """Compute lag-1 autocorrelation coefficient rho_1."""
    n = len(series)
    if n < 3:
        return 0.0
    mean = np.mean(series)
    var = np.sum((series - mean) ** 2)
    if var == 0:
        return 0.0
    cov = np.sum((series[:-1] - mean) * (series[1:] - mean))
    return float(cov / var)


def run_binsize_analysis(beta: str = "17", channel: str = "AV", test_t: int = 10) -> None:
    input_dir = DEFAULT_READIN_DIR / f"48x16b4.{beta}" / "Output"
    if not input_dir.exists():
        print(f"[ERROR] Input directory not found: {input_dir}")
        return

    orch = MesonOrchestrator()
    configs = CHANNEL_CONFIGS[channel]

    test_binsizes = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20]
    results = []

    for b in test_binsizes:
        means, errors, folded_jk, n_bins, n_cfgs = orch.process_channel(
            input_dir=str(input_dir),
            channel_configs=configs,
            binsize=b,
            num_lines=48,
            thread_count=0,
            is_single_source=False
        )

        err_at_t = errors[test_t]

        # Calculate lag-1 autocorrelation from Jackknife samples
        jk_samples = folded_jk[test_t, :]
        rho_1 = compute_lag1_autocorr(jk_samples)
        noise_threshold = 2.0 / np.sqrt(max(n_bins, 1))

        results.append({
            "binsize": b,
            "n_bins": n_bins,
            "err": err_at_t,
            "rho_1": rho_1,
            "noise_threshold": noise_threshold
        })

    # Base error at binsize = 1
    base_err = results[0]["err"]
    # Estimate plateau error from largest valid binsizes (e.g. binsize 8~16)
    plateau_err = np.mean([r["err"] for r in results if 6 <= r["binsize"] <= 16])

    print("=" * 82)
    print(f"Binsize Autocorrelation & Error Saturation Analysis (Beta 4.{beta}, Channel {channel}, t={test_t})")
    print("=" * 82)
    print(f"{'Binsize':<8} {'N_bins':<8} {'Error (sigma)':<15} {'Ratio/B=1':<12} {'Plateau%':<10} {'rho_1':<10} {'Status'}")
    print("-" * 82)

    optimal_binsize = None

    for r in results:
        b = r["binsize"]
        err = r["err"]
        ratio = err / base_err if base_err > 0 else 1.0
        pct_plateau = (err / plateau_err * 100.0) if plateau_err > 0 else 100.0
        rho = r["rho_1"]
        thresh = r["noise_threshold"]

        autocorr_eliminated = abs(rho) < thresh
        error_saturated = pct_plateau >= 95.0

        if error_saturated and autocorr_eliminated:
            status = "Saturated (Optimal)" if optimal_binsize is None else "Saturated"
            if optimal_binsize is None:
                optimal_binsize = b
        elif error_saturated:
            status = "Saturated"
        elif pct_plateau >= 88.0:
            status = "Near Saturation"
        else:
            status = "Underestimated"

        print(f"{b:<8} {r['n_bins']:<8} {err:<15.6e} {ratio:<12.3f} {pct_plateau:<9.1f}% {rho:<+9.3f} {status}")

    print("-" * 82)
    print(f"[SUMMARY] Minimal optimal binsize: {optimal_binsize}")
    print("  - Autocorrelation elimination: |rho_1| < 2/sqrt(N_bins) noise threshold")
    print("  - Error estimate saturation:   sigma(B) >= 95% of plateau error")
    print("=" * 82)


if __name__ == "__main__":
    beta_arg = sys.argv[1] if len(sys.argv) > 1 else "17"
    ch_arg = sys.argv[2] if len(sys.argv) > 2 else "AV"
    run_binsize_analysis(beta=beta_arg, channel=ch_arg)
