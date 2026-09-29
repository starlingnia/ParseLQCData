#!/usr/bin/env python3
"""Plot the AV effective mass for the 48x18, beta=4.17 ensemble."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import CHANNEL_CONFIGS, get_binsize
from src.parselqcdata import (
    compute_effective_mass_matrix_centered,
    detect_plateau_window,
    fit_jackknife_mass_centered,
)
from src.parselqcdata.plateau_fit import LEAST_SQUARES_FITTER
from tools.meson_orchestrator import MesonOrchestrator


ENSEMBLE = "48x18_b4.17_ms0.040m0.0020"
NS = 48
NT = 18
HALF = NS / 2
A_INV_MEV = 2452.96
ML_TAG = "0.0020"
INPUT_DIR = PROJECT_ROOT / "data" / "readin" / ENSEMBLE / "Output"
SAMPLE_OUTPUT_DIR = (
    PROJECT_ROOT
    / "output"
    / "meson_scan"
    / "b4.17"
    / "ccor"
    / "48x18"
)
OUTPUT_PATH = SAMPLE_OUTPUT_DIR / "meff_AV.png"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--a-inv", type=float, default=A_INV_MEV, help="inverse lattice spacing in MeV")
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        parser.error(f"input directory does not exist: {args.input_dir}")

    means, errors, jk_matrix, n_bins, n_cfgs = MesonOrchestrator().process_channel(
        input_dir=str(args.input_dir),
        channel_configs=CHANNEL_CONFIGS["AV"],
        binsize=get_binsize("17", "AV", is_single_source=False),
        num_lines=NS,
        thread_count=args.threads,
        is_single_source=False,
    )
    if jk_matrix.size == 0 or jk_matrix.shape[0] != NS:
        raise RuntimeError("AV correlator extraction returned no 48-site samples")

    meff_jk, meff_mean, meff_err = compute_effective_mass_matrix_centered(jk_matrix, HALF)
    plateau = detect_plateau_window(meff_jk, end=int(HALF), chi2_dof_max=1.5, min_points=4)
    if not plateau["ok"]:
        raise RuntimeError("No AV effective-mass plateau passed the automatic flatness check")

    fit = fit_jackknife_mass_centered(
        np.arange(NS),
        jk_matrix,
        errors,
        int(plateau["start"]),
        int(plateau["end"]),
        HALF,
    )
    if not np.isfinite(fit["mass"]):
        raise RuntimeError("AV cosh mass fit did not produce a finite result")

    n_cfgs = sum(1 for _ in args.input_dir.glob("*_mesons_multi_src"))
    mass_mev = float(fit["mass"]) * args.a_inv
    mass_mev_err = float(fit["mass_err"]) * args.a_inv
    fit_start = int(plateau["start"])
    fit_end = int(plateau["end"])
    output_dir = args.output.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    sym_data_path = output_dir / f"symdatasample{ML_TAG}.csv"
    error_path = output_dir / f"dferr1{ML_TAG}.csv"
    sample_fits_path = output_dir / f"binnedresult{ML_TAG}.csv"
    fit_result_path = output_dir / "mass_AV_fit.json"

    np.savetxt(sym_data_path, jk_matrix, delimiter=",", fmt="%.17e")
    np.savetxt(error_path, errors.reshape(-1, 1), delimiter=",", fmt="%.17e")
    fit_fields = ["chi2", "dof", "massfit_mean", "massfit_err", "fita", "fita_err", "chi2_dof", "tem"]
    with sample_fits_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fit_fields)
        writer.writeheader()
        for sample_fit in fit["per_sample"]:
            if "chi2" not in sample_fit:
                continue
            writer.writerow({
                "chi2": sample_fit["chi2"],
                "dof": sample_fit["dof"],
                "massfit_mean": sample_fit["massfit_mean"],
                "massfit_err": sample_fit["massfit_err"],
                "fita": sample_fit["fita"],
                "fita_err": sample_fit["fita_err"],
                "chi2_dof": sample_fit["chi2_dof"],
                "tem": 1.0 / NT,
            })

    with fit_result_path.open("w", encoding="utf-8") as result_file:
        fit_result = {
            "ensemble": ENSEMBLE,
            "channel": "AV",
            "fit_model": "a*cosh(m*(x-Ns/2))",
            "fit_method": "per_jackknife_correlator_chi2_least_squares",
            "fitter": LEAST_SQUARES_FITTER,
            "fit_window": [fit_start, fit_end],
            "plateau_chi2_dof": float(plateau["chi2_dof"]),
            "mass_parameter_lattice_units": float(fit["mass"]),
            "mass_parameter_error_lattice_units": float(fit["mass_err"]),
            "inverse_lattice_spacing_mev": float(args.a_inv),
            "meson_mass_mev": mass_mev,
            "meson_mass_error_mev": mass_mev_err,
            "n_bins": int(n_bins),
            "n_raw_cfgs": int(n_cfgs),
            "n_fitted_jackknife_samples": int(fit["n_fitted"]),
            "symmetrized_jackknife_csv": str(sym_data_path),
            "jackknife_error_csv": str(error_path),
            "per_sample_fit_csv": str(sample_fits_path),
        }
        json.dump(fit_result, result_file, indent=2)

    x = np.arange(NS)
    valid = (x < HALF) & np.isfinite(meff_mean) & np.isfinite(meff_err)

    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    ax.errorbar(
        x[valid],
        meff_mean[valid],
        yerr=meff_err[valid],
        fmt="o",
        markersize=4,
        capsize=2,
        color="#176b87",
        ecolor="#69a9b8",
        label="AV effective mass",
    )
    fit_x = np.arange(int(plateau["start"]), int(plateau["end"]))
    ax.hlines(
        fit["mass"],
        fit_x[0],
        fit_x[-1],
        color="#c44e52",
        linewidth=2,
        label=f"cosh fit [{fit['x_start']}, {fit['x_end']})",
    )
    ax.axhspan(
        fit["mass"] - fit["mass_err"],
        fit["mass"] + fit["mass_err"],
        color="#c44e52",
        alpha=0.16,
    )
    ax.axvspan(fit_x[0], fit_x[-1], color="#c44e52", alpha=0.06)
    ax.set(
        title=f"AV effective mass: {ENSEMBLE}",
        xlabel="Spatial separation x",
        ylabel=r"Effective mass $m_{\mathrm{eff}}$ (lattice units)",
        xlim=(0, HALF),
    )
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False)
    ax.text(
        0.98,
        0.96,
        f"fit m = {fit['mass']:.5f} +/- {fit['mass_err']:.5f}\n"
        f"bins = {n_bins}, cfgs = {n_cfgs}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.88, "edgecolor": "0.8"},
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(f"AV mass: {fit['mass']:.6f} +/- {fit['mass_err']:.6f} (lattice units)")
    print(f"AV mass: {mass_mev:.2f} +/- {mass_mev_err:.2f} MeV")
    print(f"Fit window: [{fit['x_start']}, {fit['x_end']}), plateau chi2/dof: {plateau['chi2_dof']:.3f}")
    print(f"Saved effective-mass plot: {args.output}")
    print(f"Saved symmetrized Jackknife matrix: {sym_data_path}")
    print(f"Saved Jackknife errors: {error_path}")
    print(f"Saved per-sample fits: {sample_fits_path}")
    print(f"Saved fit result: {fit_result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())