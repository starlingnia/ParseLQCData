#!/usr/bin/env python3
"""
scripts/plot_condensate.py
--------------------------------------------------------------------------------
Publication-grade plotting for Chiral Condensate observables:
Strictly adheres to ana/ formatting conventions (colors, markers, error bars, annotations):
1. Renormalized Chiral Condensate vs Temperature (conden_re_plot.png, royalblue, 'o')
2. Bare Light & Strange Quark Condensate vs Temperature (conden_plot.png)
3. Finite Volume & Temperature Scaling Series at beta=4.17 from distributed cluster runs
4. Spatial Volume Scaling (1/Ns^3) at T=153.31 MeV (Nt=16)
5. Comprehensive Multi-Ensemble Overview (All 13 Ensembles)
6. Scaled 16*T*Delta_{l,s} Observable (ana/src/tools.py standard)
--------------------------------------------------------------------------------
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    OUTPUT_CONDENSATE_DIR,
    OUTPUT_ROOT
)

DOCS_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
ANA_STYLE_DIR = OUTPUT_ROOT / "plots" / "ana_style"
ANA_STYLE_DIR.mkdir(parents=True, exist_ok=True)

CONDENSATE_RM_FILE = OUTPUT_CONDENSATE_DIR / "results_rm_beta.txt"
CONDENSATE_LIGHT_FILE = OUTPUT_CONDENSATE_DIR / "results_beta.txt"
CONDENSATE_STRANGE_FILE = OUTPUT_CONDENSATE_DIR / "results_strangequark_beta.txt"
ALL_ENSEMBLES_FILE = OUTPUT_CONDENSATE_DIR / "all_ensembles_condensate.csv"
SCALING_FILE = OUTPUT_CONDENSATE_DIR / "scaling_beta4.17.csv"

# ana/ styling palette and marker standards
COLOR_PRIMARY = "royalblue"       # Standard royalblue used in ana/src/tools.py
COLOR_SECONDARY = "#ba68c8"     # Purple used in ana/src/plot_symmetry_breaking.py
COLOR_TERTIARY = "#69b3a2"      # Teal used in ana/src/plot_symmetry_breaking.py
COLOR_ACCENT = "#ffb74d"        # Orange used in ana/src/plot_symmetry_breaking.py
COLOR_BLUE = "#64b5f6"          # Sky blue used in ana/src/plot_symmetry_breaking.py

MARKER_PRIMARY = "o"            # Circle marker standard in ana/src/tools.py
MARKER_SECONDARY = "s"          # Square marker
MARKER_TRIANGLE_UP = "^"        # Up-triangle marker
MARKER_TRIANGLE_DOWN = "v"      # Down-triangle marker
MARKER_DIAMOND = "D"            # Diamond marker


def save_multi_target(fig: plt.Figure, base_filename: str) -> None:
    """Save plot to docs/figures, output/plots/ana_style, and output/condensate in PNG and PDF."""
    targets = [
        DOCS_FIGURES_DIR,
        ANA_STYLE_DIR,
        OUTPUT_CONDENSATE_DIR
    ]
    for target in targets:
        target.mkdir(parents=True, exist_ok=True)
        fig.savefig(target / f"{base_filename}.png", dpi=300, bbox_inches="tight")
        fig.savefig(target / f"{base_filename}.pdf", bbox_inches="tight")
    print(f"[OK] Saved {base_filename}.png/.pdf to docs/figures and output/plots/ana_style")


def plot_renormalized_condensate() -> None:
    """
    Plot Renormalized Chiral Condensate vs Temperature (conden_re_plot.png).
    Strictly follows ana/src/tools.py format:
    - color: royalblue
    - fmt: 'o'
    - capsize=5, capthick=1.5, elinewidth=1.5, markersize=6
    - Point annotations: beta={beta}
    - Transition Region: [155.5, 160.5] MeV gray span
    """
    if not CONDENSATE_RM_FILE.exists():
        print(f"[WARN] File not found: {CONDENSATE_RM_FILE}")
        return

    # Read results_rm_beta.txt
    try:
        df = pl.read_csv(CONDENSATE_RM_FILE, separator="\t")
        if "beta" not in df.columns:
            df = pl.read_csv(CONDENSATE_RM_FILE, separator=" ")
    except Exception:
        df = pl.read_csv(CONDENSATE_RM_FILE, separator=" ")

    df = df.with_columns(
        pl.col("beta").cast(pl.String).alias("beta_str")
    ).with_columns(
        pl.col("beta_str").map_elements(
            lambda b: TEMP_MAP.get(str(b), TEMP_MAP.get(f"{float(b):.2f}", np.nan)),
            return_dtype=pl.Float64
        ).alias("T")
    ).sort("T")

    temps = df["T"].to_numpy()
    pbp_rm = df["pbp_rm"].to_numpy()
    pbp_rm_err = df["pbp_rm_err"].to_numpy()
    betas = df["beta_str"].to_list()

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

    # Shaded Transition Region (ana standard)
    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region")

    # Errorbar plot with strict ana style
    ax.errorbar(
        temps, pbp_rm, yerr=pbp_rm_err,
        fmt=MARKER_PRIMARY,
        color=COLOR_PRIMARY,
        ecolor=COLOR_PRIMARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=6,
        label=r'$\Delta_{\ell, s} = \frac{1}{Z_m} \left( \langle\bar{\psi}\psi\rangle_l - \frac{m_l + m_{\mathrm{res}}}{m_s + m_{\mathrm{res}}} \langle\bar{\psi}\psi\rangle_s \right)$',
        zorder=5
    )

    # Point annotations matching ana/src/tools.py exactly
    for i, b in enumerate(betas):
        ax.annotate(
            f"β={b}",
            (temps[i], pbp_rm[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=8,
            fontweight="medium"
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14)
    ax.set_ylabel(r"Renormalized Condensate $\Delta_{\ell, s}$", fontsize=14)
    ax.set_title(r"Renormalized $\bar{\psi}\psi$ vs Temperature ($48^3 \times 16$ Series)", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", frameon=True, fontsize=10)

    fig.tight_layout()
    save_multi_target(fig, "conden_re_plot")
    plt.close(fig)


def plot_bare_condensate() -> None:
    """
    Plot Bare Light and Strange Quark Condensate vs Temperature (conden_plot.png).
    Strictly follows ana/ formatting:
    - Light quark: royalblue circle
    - Strange quark: #ba68c8 square
    - capsize=5, capthick=1.5, elinewidth=1.5, markersize=6
    """
    if not CONDENSATE_RM_FILE.exists():
        return

    df = pl.read_csv(CONDENSATE_RM_FILE, separator="\t")
    df = df.with_columns(
        pl.col("beta").cast(pl.String).alias("beta_str")
    ).with_columns(
        pl.col("beta_str").map_elements(
            lambda b: TEMP_MAP.get(str(b), TEMP_MAP.get(f"{float(b):.2f}", np.nan)),
            return_dtype=pl.Float64
        ).alias("T")
    ).sort("T")

    temps = df["T"].to_numpy()
    pbpl = df["pbpl"].to_numpy()
    pbpl_err = df["pbpl_err"].to_numpy()
    pbps = df["pbps"].to_numpy()
    pbps_err = df["pbps_err"].to_numpy()
    betas = df["beta_str"].to_list()

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region")

    # Light quark condensate
    ax.errorbar(
        temps, pbpl, yerr=pbpl_err,
        fmt=MARKER_PRIMARY,
        color=COLOR_PRIMARY,
        ecolor=COLOR_PRIMARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=6,
        label=r'$\langle\bar{\psi}\psi\rangle_l$ (Light Quark)',
        zorder=5
    )

    # Strange quark condensate
    ax.errorbar(
        temps, pbps, yerr=pbps_err,
        fmt=MARKER_SECONDARY,
        color=COLOR_SECONDARY,
        ecolor=COLOR_SECONDARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=6,
        label=r'$\langle\bar{\psi}\psi\rangle_s$ (Strange Quark)',
        zorder=5
    )

    for i, b in enumerate(betas):
        ax.annotate(
            f"β={b}",
            (temps[i], pbpl[i]),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=8
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14)
    ax.set_ylabel(r"Bare Condensate $\langle\bar{\psi}\psi\rangle$", fontsize=14)
    ax.set_title(r"Bare Light & Strange Quark Condensate vs Temperature", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="center right", frameon=True, fontsize=10)

    fig.tight_layout()
    save_multi_target(fig, "conden_plot")
    plt.close(fig)


def plot_scaled_16T_observable() -> None:
    """
    Plot exact 16*T*Delta_{l,s} observable conforming to ana/src/tools.py:
    Scaled Value = 16 * T * Delta_{l,s} vs Temperature.
    """
    if not CONDENSATE_RM_FILE.exists():
        return

    df = pl.read_csv(CONDENSATE_RM_FILE, separator="\t")
    df = df.with_columns(
        pl.col("beta").cast(pl.String).alias("beta_str")
    ).with_columns(
        pl.col("beta_str").map_elements(
            lambda b: TEMP_MAP.get(str(b), TEMP_MAP.get(f"{float(b):.2f}", np.nan)),
            return_dtype=pl.Float64
        ).alias("T")
    ).sort("T")

    temps = df["T"].to_numpy()
    pbp_rm = df["pbp_rm"].to_numpy()
    pbp_rm_err = df["pbp_rm_err"].to_numpy()
    betas = df["beta_str"].to_list()

    scale_factor = 16.0 * temps
    scaled_mean = pbp_rm * scale_factor
    scaled_err = pbp_rm_err * scale_factor

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region")

    ax.errorbar(
        temps, scaled_mean, yerr=scaled_err,
        fmt=MARKER_PRIMARY,
        color=COLOR_PRIMARY,
        ecolor=COLOR_PRIMARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=6,
        label=r'$16 \times T \times \Delta_{\ell, s}$',
        zorder=5
    )

    for i, b in enumerate(betas):
        ax.annotate(
            f"β={b}",
            (temps[i], scaled_mean[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=8
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14)
    ax.set_ylabel(r"Scaled Value ($16 \times T \times \Delta_{\ell, s}$) [MeV]", fontsize=14)
    ax.set_title(r"Scaled $\bar{\psi}\psi$ vs Temperature (ana Standard)", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", frameon=True, fontsize=10)

    fig.tight_layout()
    save_multi_target(fig, "pbp_temperature_plot")
    plt.close(fig)


def plot_distributed_scaling_beta417() -> None:
    """
    Extract data calculated by other computers (output/condensate/scaling_beta4.17.csv).
    Plots finite volume and temporal extent evolution at beta=4.17:
    - L32T12 (T=204.41 MeV)
    - L32T14 (T=175.21 MeV)
    - L40T16 (T=153.31 MeV)
    - L48T16 (T=153.31 MeV)
    - L48T18 (T=136.28 MeV, Streams 1 & 2)
    """
    if not SCALING_FILE.exists():
        print(f"[WARN] File not found: {SCALING_FILE}")
        return

    df = pl.read_csv(SCALING_FILE).sort("temperature")

    fig, ax = plt.subplots(figsize=(9, 6), dpi=100)

    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region")

    # Plot each ensemble using ana formatting
    temps = df["temperature"].to_numpy()
    pbp_rm = df["pbp_rm"].to_numpy()
    pbp_rm_err = df["pbp_rm_err"].to_numpy()
    names = df["dataset_name"].to_list()
    ns_list = df["ns"].to_list()
    nt_list = df["nt"].to_list()

    # Differentiate streams or geometries
    ax.errorbar(
        temps, pbp_rm, yerr=pbp_rm_err,
        fmt=MARKER_PRIMARY,
        color=COLOR_PRIMARY,
        ecolor=COLOR_PRIMARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'Scaling Ensembles ($\beta=4.17$)',
        zorder=5
    )

    for i in range(len(names)):
        tag = f"${ns_list[i]}^3 \\times {nt_list[i]}$"
        if "_2" in names[i]:
            tag += " (str2)"
            offset = (28, 6)
        elif nt_list[i] == 18:
            offset = (0, 12)
        elif nt_list[i] == 16:
            offset = (22, 8)
        elif nt_list[i] == 14:
            offset = (0, 12)
        else:
            offset = (0, -18)

        ax.annotate(
            tag,
            (temps[i], pbp_rm[i]),
            textcoords="offset points",
            xytext=offset,
            ha="center",
            fontsize=8,
            fontweight="semibold"
        )

    ax.axhline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14)
    ax.set_ylabel(r"Renormalized Condensate $\Delta_{\ell, s}$", fontsize=14)
    ax.set_title(r"Finite Temporal Scaling Condensate at $\beta=4.17$", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", frameon=True, fontsize=10)

    fig.tight_layout()
    save_multi_target(fig, "conden_scaling_beta4.17")
    plt.close(fig)


def plot_volume_scaling_fixed_temp() -> None:
    """
    Finite spatial volume scaling at fixed beta=4.17, Nt=16 (T=153.31 MeV):
    Extracts Ns=40 and Ns=48 from all_ensembles_condensate.csv to demonstrate
    chiral condensate volume convergence towards the thermodynamic limit.
    """
    if not ALL_ENSEMBLES_FILE.exists():
        return

    df = pl.read_csv(ALL_ENSEMBLES_FILE)
    sub = df.filter(
        (pl.col("beta") == 4.17) & (pl.col("nt") == 16)
    ).sort("ns")

    if sub.height < 2:
        return

    ns = sub["ns"].to_numpy()
    inv_v = 1.0 / (ns ** 3)
    pbp_rm = sub["pbp_rm"].to_numpy()
    pbp_rm_err = sub["pbp_rm_err"].to_numpy()
    pbpl = sub["pbpl"].to_numpy()
    pbpl_err = sub["pbpl_err"].to_numpy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=100)

    # Panel 1: Renormalized Condensate vs 1/Ns^3
    ax1.errorbar(
        inv_v, pbp_rm, yerr=pbp_rm_err,
        fmt=MARKER_PRIMARY,
        color=COLOR_PRIMARY,
        ecolor=COLOR_PRIMARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'$\Delta_{\ell, s}$ vs $1/N_s^3$'
    )
    for i, n in enumerate(ns):
        ax1.annotate(f"$N_s={n}$", (inv_v[i], pbp_rm[i]), textcoords="offset points", xytext=(0, 10), ha="center", fontsize=9)

    ax1.set_xlabel(r"Inverse Spatial Volume $1/N_s^3$", fontsize=13)
    ax1.set_ylabel(r"Renormalized Condensate $\Delta_{\ell, s}$", fontsize=13)
    ax1.set_title(r"Finite Volume Scaling of $\Delta_{\ell, s}$ ($T=153.31$ MeV)", fontsize=13)
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(loc="best", fontsize=10)

    # Panel 2: Bare Light Condensate vs Ns
    ax2.errorbar(
        ns, pbpl, yerr=pbpl_err,
        fmt=MARKER_SECONDARY,
        color=COLOR_SECONDARY,
        ecolor=COLOR_SECONDARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'$\langle\bar{\psi}\psi\rangle_l$ vs $N_s$'
    )
    for i, n in enumerate(ns):
        ax2.annotate(f"$N_s={n}$", (ns[i], pbpl[i]), textcoords="offset points", xytext=(0, 10), ha="center", fontsize=9)

    ax2.set_xlabel(r"Spatial Lattice Extent $N_s$", fontsize=13)
    ax2.set_ylabel(r"Bare Light Condensate $\langle\bar{\psi}\psi\rangle_l$", fontsize=13)
    ax2.set_title(r"Volume Dependence of $\langle\bar{\psi}\psi\rangle_l$ ($T=153.31$ MeV)", fontsize=13)
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(loc="best", fontsize=10)

    fig.tight_layout()
    save_multi_target(fig, "conden_volume_scaling")
    plt.close(fig)


def plot_all_ensembles_overview() -> None:
    """
    Global comprehensive overview comparing all 13 ensembles:
    1. Standard 48^3 x 16 Temperature Scan Series (Physical Quark Masses)
    2. Distributed Scaling Series at beta=4.17 (Quark Mass ml=0.0020, ms=0.040)
    """
    if not ALL_ENSEMBLES_FILE.exists():
        return

    df = pl.read_csv(ALL_ENSEMBLES_FILE)

    # Split into Standard Scan (48x16, physical mass ml <= 0.00105) and Scaling (other computers)
    std_df = df.filter((pl.col("ns") == 48) & (pl.col("nt") == 16)).sort("temperature")
    dist_df = df.filter(~((pl.col("ns") == 48) & (pl.col("nt") == 16))).sort("temperature")

    fig, ax = plt.subplots(figsize=(10, 6.5), dpi=100)

    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region ($155.5-160.5$ MeV)")

    # Standard temperature scan
    ax.errorbar(
        std_df["temperature"].to_numpy(),
        std_df["pbp_rm"].to_numpy(),
        yerr=std_df["pbp_rm_err"].to_numpy(),
        fmt=MARKER_PRIMARY,
        color=COLOR_PRIMARY,
        ecolor=COLOR_PRIMARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'Physical Scan ($48^3 \times 16$, $m_l \approx 0.0010$)',
        zorder=5
    )

    # Distributed cluster runs
    ax.errorbar(
        dist_df["temperature"].to_numpy(),
        dist_df["pbp_rm"].to_numpy(),
        yerr=dist_df["pbp_rm_err"].to_numpy(),
        fmt=MARKER_SECONDARY,
        color=COLOR_SECONDARY,
        ecolor=COLOR_SECONDARY,
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'Distributed Scaling Runs ($\beta=4.17$, $m_l = 0.0020$)',
        zorder=6
    )

    # Annotate standard scan points
    for row in std_df.iter_rows(named=True):
        b = f"{row['beta']:.2f}" if abs(row['beta'] - 4.405) > 1e-4 else "4.405"
        ax.annotate(
            f"β={b}",
            (row["temperature"], row["pbp_rm"]),
            textcoords="offset points",
            xytext=(0, 9),
            ha="center",
            fontsize=8
        )

    # Annotate distributed points
    for row in dist_df.iter_rows(named=True):
        tag = f"${row['ns']}^3\\times{row['nt']}$"
        if "_2" in row["dataset_name"]:
            tag += " (str2)"
            offset = (32, 4)
        elif row["nt"] == 18:
            offset = (-26, 6)
        elif row["nt"] == 16:
            offset = (26, -14)
        elif row["nt"] == 14:
            offset = (0, -18)
        else:
            offset = (0, -18)

        ax.annotate(
            tag,
            (row["temperature"], row["pbp_rm"]),
            textcoords="offset points",
            xytext=offset,
            ha="center",
            fontsize=8,
            color="#6a1b9a",
            fontweight="semibold"
        )

    ax.axhline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14)
    ax.set_ylabel(r"Renormalized Chiral Condensate $\Delta_{\ell, s}$", fontsize=14)
    ax.set_title("Global Chiral Condensate Evolution Across All 13 Ensembles", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", frameon=True, fontsize=10)

    fig.tight_layout()
    save_multi_target(fig, "conden_all_ensembles_overview")
    plt.close(fig)


def main() -> None:
    print("==========================================================================")
    print("[INFO] Starting Chiral Condensate plotting pipeline (Strict ana/ Style)...")
    print(f"[INFO] Source data directory: {OUTPUT_CONDENSATE_DIR}")
    print(f"[INFO] Target figures directory: {DOCS_FIGURES_DIR}")
    print("==========================================================================")

    plot_renormalized_condensate()
    plot_bare_condensate()
    plot_scaled_16T_observable()
    plot_distributed_scaling_beta417()
    plot_volume_scaling_fixed_temp()
    plot_all_ensembles_overview()

    print("==========================================================================")
    print("[INFO] All Chiral Condensate plots generated successfully.")
    print("==========================================================================")


if __name__ == "__main__":
    main()
