#!/usr/bin/env python3
"""
scripts/run_meson.py
--------------------------------------------------------------------------------
介子强子关联函数标准分析与物理拟合端到端运行脚本:
1. 通过 MesonOrchestrator (C++ 核心库) 批量抽取多源/单源空间关联函数
2. 内存级直接流转求解有效质量方程 (C(t)/C(t+1) = cosh(m*(t-24))/cosh(m*(t-23)))
3. 执行贝叶斯/最小二乘 cosh 平台拟合 (Jackknife 统计推断) 并导出汇总
4. 可选自动生成介子物理图表 (--plot)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from scipy.optimize import fsolve

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    BETAS,
    CHANNELS,
    CHANNEL_CONFIGS,
    DEFAULT_READIN_DIR,
    FIT_SLICE_END,
    MULTI_FIT_SLICES,
    OUTPUT_ROOT,
    SINGLE_FIT_SLICES,
    get_binsize,
)
from src.parselqcdata import fit_single_jackknife_column
from tools.meson_orchestrator import MesonOrchestrator


def solve_effective_mass_ratio(y: float, x: int) -> float:
    """求解单个切片比值的有效质量方程 (与 ana 数值 100% 对齐)"""
    if not np.isfinite(y):
        return np.nan

    def equation(m):
        a = m * (x - 24.0)
        b = m * (x + 1.0 - 24.0)
        num = np.exp(a) + np.exp(-a)
        den = np.exp(b) + np.exp(-b)
        return num / den - y

    try:
        sol = fsolve(equation, x0=0.1, xtol=1e-15, maxfev=5000)
        return float(sol[0])
    except Exception:
        return np.nan


def compute_meff_from_matrix(dr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """从 (num_lines x n_bins) 折叠矩阵直接计算有效质量均值与 Jackknife 误差"""
    nrow, ncol = dr.shape
    meff = np.full((nrow, ncol), np.nan, dtype=np.float64)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(nrow - 1):
            for j in range(ncol):
                c_next = dr[i + 1, j]
                if abs(c_next) > 1e-16:
                    y = dr[i, j] / c_next
                    meff[i, j] = solve_effective_mass_ratio(y, i)

    meff_mean = np.nanmean(meff, axis=1)
    diffs2 = (meff - meff_mean[:, None]) ** 2
    meff_err = np.sqrt((ncol - 1) * np.nansum(diffs2, axis=1) / ncol)
    return meff_mean, meff_err


def fit_plateau_from_matrix(
    dr: np.ndarray,
    errors: np.ndarray,
    t_start: int,
    t_end: int,
    beta: str,
    channel: str,
    out_dir: Optional[Path] = None,
) -> dict:
    """对关联函数切片执行逐 Jackknife 样本 cosh 平台拟合并汇总"""
    x_array = np.arange(t_start, t_end)
    y_slice = dr[t_start:t_end, :]
    err_slice = errors[t_start:t_end]
    n_bins = y_slice.shape[1]

    fit_records = []
    mass_list = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(n_bins):
            try:
                fit = fit_single_jackknife_column(x_array, y_slice[:, k], err_slice)
                if fit:
                    fit_records.append(fit)
                    mass_list.append(float(fit.get("massfit_mean", np.nan)))
                else:
                    mass_list.append(np.nan)
            except Exception:
                mass_list.append(np.nan)

    if out_dir is not None and fit_records:
        pl.DataFrame(fit_records).write_csv(out_dir / f"fittresult{channel}.csv")

    mass_arr = np.array(mass_list, dtype=np.float64)
    valid_mask = np.isfinite(mass_arr)

    if np.any(valid_mask):
        fit_mass = float(np.mean(mass_arr[valid_mask]))
        diff = mass_arr[valid_mask] - fit_mass
        fit_err = float(np.sqrt((n_bins - 1) * np.sum(diff ** 2) / n_bins))
    else:
        fit_mass, fit_err = np.nan, np.nan

    if out_dir is not None:
        summary_df = pl.DataFrame({
            "quantity": ["mass"],
            "mean": [fit_mass],
            "jack_err": [fit_err]
        })
        summary_df.write_csv(out_dir / f"summary_fit_{channel}.csv")

    return {
        "beta": beta,
        "channel": channel,
        "mass_mean": fit_mass,
        "mass_err": fit_err,
        "t_start": t_start,
        "t_end": t_end,
    }


def run_single_channel_workflow(
    orch: MesonOrchestrator,
    input_dir: Path,
    beta: str,
    channel: str,
    is_single: bool,
    slice_start: int,
    pick_dir: Path,
    ratio_dir: Path,
    simulate_dir: Path,
) -> dict:
    """单个信道端到端流式计算：抽取 -> 有效质量 -> 平台拟合 -> 落盘标准产物"""
    configs = CHANNEL_CONFIGS[channel]
    binsize = get_binsize(beta, channel, is_single)

    means, errors, folded_jk, n_bins, n_cfgs = orch.process_channel(
        input_dir=str(input_dir),
        channel_configs=configs,
        binsize=binsize,
        num_lines=48,
        thread_count=0,
        is_single_source=is_single,
    )

    # 1. 保存关联函数标准 CSV
    pl.DataFrame({"mean": means, "err": errors}).write_csv(pick_dir / f"save_{channel}.csv")
    pl.DataFrame({"mean": means[:25], "err": errors[:25]}).write_csv(pick_dir / f"sym_{channel}.csv")
    pl.DataFrame(folded_jk).write_csv(pick_dir / f"dr_{channel}.csv", include_header=False)
    pl.DataFrame(errors).write_csv(pick_dir / f"err_{channel}.csv", include_header=False)

    # 2. 直接在内存中求解有效质量并落盘
    meff_mean, meff_err = compute_meff_from_matrix(folded_jk)
    pl.DataFrame({"mean": meff_mean, "err": meff_err}).write_csv(ratio_dir / f"meff_{channel}.csv")

    # 3. 直接在内存中拟合平台质量并落盘
    fit_record = fit_plateau_from_matrix(
        dr=folded_jk,
        errors=errors,
        t_start=slice_start,
        t_end=FIT_SLICE_END,
        beta=beta,
        channel=channel,
        out_dir=simulate_dir,
    )

    return fit_record


def run_meson_pipeline(
    betas: Sequence[str] = BETAS,
    sources: Sequence[bool] = (False, True),
    readin_dir: Path = DEFAULT_READIN_DIR,
    output_root: Path = OUTPUT_ROOT,
    do_plot: bool = False,
) -> None:
    t0_all = time.time()
    orch = MesonOrchestrator()

    for is_single in sources:
        mode_str = "singlesrc" if is_single else "multisrc"
        dir_pick = "pickdata-singlesrc" if is_single else "pickdata"
        dir_ratio = "ratio_results-singlesrc" if is_single else "ratio_results"
        dir_sim = "simulateresult-singlesrc" if is_single else "simulateresult"
        slice_map = SINGLE_FIT_SLICES if is_single else MULTI_FIT_SLICES

        print(f"\n==========================================================================")
        print(f"[MESON] 启动 {mode_str.upper()} 信道流水线处理...")
        print(f"==========================================================================")

        all_fits = []

        for beta in betas:
            input_dir = readin_dir / f"48x16b4.{beta}" / "Output"
            if not input_dir.exists():
                print(f"  [WARN] 目录不存在，跳过: {input_dir}")
                continue

            pick_dir = output_root / dir_pick / f"b4.{beta}"
            ratio_dir = output_root / dir_ratio / f"b4.{beta}"
            sim_dir = output_root / dir_sim / f"b4.{beta}"

            pick_dir.mkdir(parents=True, exist_ok=True)
            ratio_dir.mkdir(parents=True, exist_ok=True)
            sim_dir.mkdir(parents=True, exist_ok=True)

            t0_beta = time.time()
            for ch in CHANNELS:
                t_start = slice_map[beta][ch]
                fit_rec = run_single_channel_workflow(
                    orch=orch,
                    input_dir=input_dir,
                    beta=beta,
                    channel=ch,
                    is_single=is_single,
                    slice_start=t_start,
                    pick_dir=pick_dir,
                    ratio_dir=ratio_dir,
                    simulate_dir=sim_dir,
                )
                all_fits.append(fit_rec)

            elapsed_beta = time.time() - t0_beta
            print(f"  [OK] Beta 4.{beta}: 6 信道全流程完成 (耗时: {elapsed_beta:.2f}s)")

        if all_fits:
            summary_df = pl.DataFrame(all_fits).sort(["beta", "channel"])
            summary_path = output_root / dir_sim / "all_fits_summary.csv"
            summary_df.write_csv(summary_path)
            print(f"  -> 汇总拟合表已落盘: {summary_path}")

    print(f"\n[OK] 介子全流程分析全部执行完成！总耗时: {time.time() - t0_all:.2f}s")

    if do_plot:
        print("[INFO] 正在调用 plot_meson.py 绘制全信道图谱...")
        import scripts.plot_meson as pm
        pm.main()


def main() -> None:
    parser = argparse.ArgumentParser(description="介子强子关联函数标准分析与物理拟合工具")
    parser.add_argument(
        "--source",
        choices=["all", "multi", "single"],
        default="all",
        help="数据源模式: multi (多源), single (单源), 或 all (全部)",
    )
    parser.add_argument(
        "--beta",
        type=str,
        default="all",
        help="指定分析特定 Beta (如 17 或 4.17)，或 'all'",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="计算完成后自动调用绘图脚本生成出版级物理图表",
    )
    args = parser.parse_args()

    if args.source == "multi":
        sources = [False]
    elif args.source == "single":
        sources = [True]
    else:
        sources = [False, True]

    if args.beta == "all":
        betas = BETAS
    else:
        b_clean = args.beta.replace("b4.", "").replace("4.", "")
        betas = [b_clean]

    run_meson_pipeline(betas=betas, sources=sources, do_plot=args.plot)


if __name__ == "__main__":
    main()
