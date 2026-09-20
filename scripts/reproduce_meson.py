#!/usr/bin/env python3
"""
scripts/reproduce_meson.py
--------------------------------------------------------------------------------
Meson correlator data processing pipeline (Polars + NumPy accelerated):
Supports mode argument: 'multi', 'single', or 'all' (default)
1. Extract multi-source & single-source correlators via MesonOrchestrator
2. Solve effective mass (meff) computing Ratio C(t)/C(t+1) strictly matched to ana
3. Perform Bayesian cosh plateau fits via lsqfit (Strict Jackknife alignment)
--------------------------------------------------------------------------------
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys
import time
import warnings
from typing import Dict, List, Tuple

import numpy as np
import polars as pl
from scipy.optimize import fsolve

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    BETAS, CHANNELS, CHANNEL_CONFIGS, MULTI_FIT_SLICES, SINGLE_FIT_SLICES,
    FIT_SLICE_END, get_binsize, DEFAULT_READIN_DIR
)
from src.parselqcdata import fit_single_jackknife_column
from tools.meson_orchestrator import MesonOrchestrator

READIN_DIR = DEFAULT_READIN_DIR
OUTPUT_ROOT = PROJECT_ROOT / "output"


def solve_effective_mass_one_point_exact(y: float, x: int) -> float:
    """完美复刻原 ana 的非线性方程组求根算法"""
    if np.isnan(y):
        return np.nan

    def equation(m):
        a = m * (x - 24.0)
        b = m * (x + 1.0 - 24.0)
        num = np.exp(a) + np.exp(-a)
        den = np.exp(b) + np.exp(-b)
        return num / den - y
    
    try:
        # 完全采用原版参数 x0=0.1, xtol=1e-15, maxfev=5000
        sol = fsolve(equation, x0=0.1, xtol=1e-15, maxfev=5000)
        return float(sol[0])
    except Exception:
        return np.nan


def compute_single_meff(args: Tuple[str, str, str, str]) -> str:
    dr_path_str, out_path_str, beta, ch = args
    dr_path = Path(dr_path_str)
    out_path = Path(out_path_str)

    if not dr_path.exists():
        return f"[WARN] Missing {dr_path}"

    # 此处的 dr 为关联函数 (Correlator)，来自 C++ 输出的 folded_jk_matrix
    dr = pl.read_csv(dr_path, has_header=False).to_numpy()
    nrow, ncol = dr.shape
    J = ncol

    meff = np.zeros((nrow, ncol), dtype=np.float64)

    # 抑制 fsolve 大数值溢出导致的 RuntimeWarning，保持输出整洁
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(nrow):
            x = i
            for j in range(ncol):
                if i < nrow - 1:
                    c_t = dr[i, j]
                    c_next = dr[i+1, j]
                    # 必须计算比值 Ratio！负负得正，无需顾虑负信道
                    if abs(c_next) > 1e-16:
                        y = c_t / c_next
                        meff[i, j] = solve_effective_mass_one_point_exact(y, x)
                    else:
                        meff[i, j] = np.nan
                else:
                    meff[i, j] = np.nan

    # Jackknife 统计：完全复刻原版中的 J 缩放
    meff_mean = np.nanmean(meff, axis=1)
    diffs2_m = (meff - meff_mean[:, None]) ** 2
    meff_err = np.sqrt((J - 1) * np.nansum(diffs2_m, axis=1) / J)

    out_df = pl.DataFrame({"mean": meff_mean, "err": meff_err})
    out_df.write_csv(out_path)
    return f"[OK] meff b4.{beta} {ch}"


def compute_single_fit(args: Tuple[str, str, str, int, int, str, str]) -> dict:
    dr_path_str, err_path_str, out_dir_str, t_start, t_end, beta, ch = args
    dr_path = Path(dr_path_str)
    err_path = Path(err_path_str)
    out_dir = Path(out_dir_str)

    default_result = {
        "beta": beta, "channel": ch,
        "mass_mean": np.nan, "mass_err": np.nan,
        "t_start": t_start, "t_end": t_end
    }

    if not dr_path.exists() or not err_path.exists():
        return default_result

    # 拟合引擎直接读取关联函数进行 cosh 拟合
    dr = pl.read_csv(dr_path, has_header=False).to_numpy()
    err_col = pl.read_csv(err_path, has_header=False).to_series(0).to_numpy()

    x_array = np.arange(t_start, t_end)
    y_slice = dr[t_start:t_end, :]
    err_slice = err_col[t_start:t_end]

    J = y_slice.shape[1]
    fit_records = []
    mass_list = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(J):
            try:
                fit = fit_single_jackknife_column(x_array, y_slice[:, k], err_slice)
                if fit:
                    fit_records.append(fit)
                    m_val = fit.get("massfit_mean", fit.get("mass", fit.get("m", np.nan)))
                    mass_list.append(float(m_val))
                else:
                    mass_list.append(np.nan)
            except Exception:
                mass_list.append(np.nan)

    if fit_records:
        pl.DataFrame(fit_records).write_csv(out_dir / f"fittresult{ch}.csv")

    mass_arr = np.array(mass_list, dtype=np.float64)
    valid_mask = np.isfinite(mass_arr)
    
    if np.any(valid_mask):
        # 原版算法缩放复刻：均值使用 valid 样本，缩放倍率强制除以总数 J
        fit_mass = float(np.mean(mass_arr[valid_mask]))
        diff = mass_arr[valid_mask] - fit_mass
        fit_err = float(np.sqrt((J - 1) * np.sum(diff ** 2) / J))
    else:
        fit_mass, fit_err = np.nan, np.nan

    summary_df = pl.DataFrame({
        "quantity": ["mass"],
        "mean": [fit_mass],
        "jack_err": [fit_err]
    })
    summary_df.write_csv(out_dir / f"summary_fit_{ch}.csv")

    return {
        "beta": beta,
        "channel": ch,
        "mass_mean": fit_mass,
        "mass_err": fit_err,
        "t_start": t_start,
        "t_end": t_end
    }


def main() -> None:
    t_start_total = time.time()
    print("[INFO] Starting Meson correlator reproduction pipeline (Polars + NumPy accelerated)...")

    orch = MesonOrchestrator()

    arg = sys.argv[1].lower() if len(sys.argv) > 1 else "all"
    if arg == "single":
        modes = [True]
    elif arg == "multi":
        modes = [False]
    else:
        modes = [False, True]

    for is_single in modes:
        mode_str = "single" if is_single else "multi"
        print(f"\n[INFO] Mode: {mode_str.upper()}-source extraction")

        for beta in BETAS:
            input_dir = READIN_DIR / f"48x16b4.{beta}" / "Output"
            if not input_dir.exists():
                print(f"[WARN] Directory not found: {input_dir}, skipping")
                continue

            pick_dir = OUTPUT_ROOT / (f"pickdata-singlesrc/b4.{beta}" if is_single else f"pickdata/b4.{beta}")
            pick_dir.mkdir(parents=True, exist_ok=True)

            t0 = time.time()
            for ch in CHANNELS:
                configs = CHANNEL_CONFIGS[ch]
                binsize = get_binsize(beta, ch, is_single)

                # 严格使用 kwargs 调用匹配 C Python 接口参数结构
                means, errors, folded_jk, n_bins, n_cfgs = orch.process_channel(
                    input_dir=str(input_dir),
                    channel_configs=configs,
                    binsize=binsize,
                    num_lines=48,
                    thread_count=0,
                    is_single_source=is_single
                )

                pl.DataFrame({"mean": means, "err": errors}).write_csv(pick_dir / f"save_{ch}.csv")
                pl.DataFrame({"mean": means[:25], "err": errors[:25]}).write_csv(pick_dir / f"sym_{ch}.csv")
                pl.DataFrame(folded_jk).write_csv(pick_dir / f"dr_{ch}.csv", include_header=False)
                pl.DataFrame(errors).write_csv(pick_dir / f"err_{ch}.csv", include_header=False)

            elapsed = time.time() - t0
            print(f"  - Beta 4.{beta}: 6 channels extracted ({elapsed:.2f}s, {n_cfgs} cfgs, {n_bins} bins)")

    print("\n[INFO] Computing effective mass (meff)...")
    meff_tasks = []
    for is_single in modes:
        dir_name = "pickdata-singlesrc" if is_single else "pickdata"
        out_name = "ratio_results-singlesrc" if is_single else "ratio_results"
        for beta in BETAS:
            out_dir = OUTPUT_ROOT / out_name / f"b4.{beta}"
            out_dir.mkdir(parents=True, exist_ok=True)
            for ch in CHANNELS:
                dr_path = OUTPUT_ROOT / dir_name / f"b4.{beta}" / f"dr_{ch}.csv"
                out_path = out_dir / f"meff_{ch}.csv"
                meff_tasks.append((str(dr_path), str(out_path), beta, ch))

    t0_meff = time.time()
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(compute_single_meff, task) for task in meff_tasks]
        for f in as_completed(futures):
            f.result()
    print(f"  - Meff tasks completed in {time.time() - t0_meff:.2f}s")

    print("\n[INFO] Performing Bayesian cosh plateau fits...")
    for is_single in modes:
        src_name = "pickdata-singlesrc" if is_single else "pickdata"
        res_name = "simulateresult-singlesrc" if is_single else "simulateresult"
        slice_map = SINGLE_FIT_SLICES if is_single else MULTI_FIT_SLICES

        fit_tasks = []
        for beta in BETAS:
            out_dir = OUTPUT_ROOT / res_name / f"b4.{beta}"
            out_dir.mkdir(parents=True, exist_ok=True)
            for ch in CHANNELS:
                dr_path = OUTPUT_ROOT / src_name / f"b4.{beta}" / f"dr_{ch}.csv"
                err_path = OUTPUT_ROOT / src_name / f"b4.{beta}" / f"err_{ch}.csv"
                t_start = slice_map[beta][ch]
                t_end = FIT_SLICE_END
                fit_tasks.append((str(dr_path), str(err_path), str(out_dir), t_start, t_end, beta, ch))

        fit_results = []
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(compute_single_fit, task) for task in fit_tasks]
            for f in as_completed(futures):
                fit_results.append(f.result())

        summary_df = pl.DataFrame(fit_results).sort(["beta", "channel"])
        summary_path = OUTPUT_ROOT / res_name / "all_fits_summary.csv"
        summary_df.write_csv(summary_path)
        print(f"  - Fits saved: {summary_path}")

    total_time = time.time() - t_start_total
    print(f"\n[INFO] Meson pipeline completed successfully in {total_time:.2f}s")


if __name__ == "__main__":
    main()