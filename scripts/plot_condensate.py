#!/usr/bin/env python3
"""
scripts/plot_condensate.py
Publication-grade plotting for Chiral Condensate observables (Polars accelerated):
1. Renormalized Chiral Condensate vs Temperature (conden_re_plot.png)
2. Bare Light Quark Condensate vs Temperature (conden_plot.png)
"""

import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    TEMP_MAP, TRANSITION_REGION,
    FIGURES_DIR as DOCS_FIGURES_DIR,
    OUTPUT_CONDENSATE_DIR
)
DOCS_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
CONDENSATE_FILE = OUTPUT_CONDENSATE_DIR / "results_rm_beta.txt"


def main() -> None:
    if not CONDENSATE_FILE.exists():
        print(f"[WARN] File not found: {CONDENSATE_FILE}")
        return

    print("[INFO] Generating Chiral Condensate plots (Polars accelerated)...")
    # Read tab- or space-delimited text file
    try:
        df = pl.read_csv(CONDENSATE_FILE, separator="\t")
        if "beta" not in df.columns:
            df = pl.read_csv(CONDENSATE_FILE, separator=" ")
    except Exception:
        df = pl.read_csv(CONDENSATE_FILE, separator=" ")

    df = df.with_columns(
        pl.col("beta").map_elements(
            lambda b: TEMP_MAP.get(str(b), TEMP_MAP.get(f"{float(b):.2f}", np.nan)), return_dtype=pl.Float64
        ).alias("T")
    )

    t_low, t_high = TRANSITION_REGION

    # 1. Renormalized Condensate Plot
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.errorbar(
        df['T'].to_numpy(), df['pbp_rm'].to_numpy(), yerr=df['pbp_rm_err'].to_numpy(),
        fmt='o', color='#2ca02c', ecolor='#2ca02c',
        elinewidth=1.5, capsize=3, capthick=1.5,
        label=r'$\Delta_{\ell, s} = \langle\bar{\psi}\psi\rangle_l - \frac{m_l+m_{\mathrm{res}}}{m_s+m_{\mathrm{res}}}\langle\bar{\psi}\psi\rangle_s$',
        markersize=7, zorder=5
    )
    ax1.axvspan(t_low, t_high, color='gray', alpha=0.15, label='Transition Region')
    ax1.set_xlabel(r'$T$ [MeV]', fontsize=14)
    ax1.set_ylabel(r'$\Delta_{\ell, s}$', fontsize=14)
    ax1.set_title(r'Renormalized Chiral Condensate vs Temperature', fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(loc='best', frameon=True, fontsize=10)

    fig1.tight_layout()
    out_re = DOCS_FIGURES_DIR / "conden_re_plot.png"
    fig1.savefig(out_re, dpi=300)
    plt.close(fig1)
    print(f"[INFO] Renormalized condensate plot saved: {out_re}")

    # 2. Bare Light Condensate Plot
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.errorbar(
        df['T'].to_numpy(), df['pbpl'].to_numpy(), yerr=df['pbpl_err'].to_numpy(),
        fmt='s', color='#d62728', ecolor='#d62728',
        elinewidth=1.5, capsize=3, capthick=1.5,
        label=r'$\langle\bar{\psi}\psi\rangle_l$ (Bare)',
        markersize=7, zorder=5
    )
    ax2.axvspan(t_low, t_high, color='gray', alpha=0.15, label='Transition Region')
    ax2.set_xlabel(r'$T$ [MeV]', fontsize=14)
    ax2.set_ylabel(r'$\langle\bar{\psi}\psi\rangle_l$', fontsize=14)
    ax2.set_title(r'Bare Light Quark Condensate vs Temperature', fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='best', frameon=True, fontsize=10)

    fig2.tight_layout()
    out_bare = DOCS_FIGURES_DIR / "conden_plot.png"
    fig2.savefig(out_bare, dpi=300)
    plt.close(fig2)
    print(f"[INFO] Bare light condensate plot saved: {out_bare}")


if __name__ == "__main__":
    main()
