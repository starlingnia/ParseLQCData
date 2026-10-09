#!/usr/bin/env python3
"""
scripts/lcp/run_lcp_meson.py
--------------------------------------------------------------------------------
【LCP 专区】常物理线 (Line of Constant Physics) 48^3 x 16 介子关联函数分析脚本:
1. 通过 MesonOrchestrator (C++ 核心库) 批量抽取 7 组 Beta 的 6 种信道多源/单源关联函数
2. 调用 effective_mass_solver 内存级求解有效质量方程 (中心对称点 24.0)
3. 调用 cosh_fitter 执行非线性平台拟合与 Jackknife 统计推断
4. 输出结果至 output/pickdata*, output/ratio_results*, output/simulateresult*
5. 可选自动生成 LCP 介子物理图表 (--plot)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys
import time
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
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
from src.parselqcdata.effective_mass_solver import (
    compute_effective_mass_matrix_centered,
    solve_effective_mass_one_point_centered,
)
from src.parselqcdata.cosh_fitter import (
    fit_single_jackknife_column_centered,
)
from tools.meson_orchestrator import MesonOrchestrator


def solve_effective_mass_ratio(y: float, x: int) -> float:
    """求解单个切片比值的有效质量方程 (与 ana 数值 100% 对齐)"""
    return solve_effective_mass_one_point_centered(y, x, half=24.0)


def compute_meff_from_matrix(dr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """从 (num_lines x n_bins) 折叠矩阵直接计算有效质量均值与 Jackknife 误差"""
    _, meff_mean, meff_err = compute_effective_mass_matrix_centered(dr, half=24.0, vectorized=True)
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
                fit = fit_single_jackknife_column_centered(x_array, y_slice[:, k], err_slice, center=24.0)
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


def run_single_task(
    beta: str,
    channel: str,
    is_single: bool,
    readin_dir: Path,
    out_root: Path,
) -> dict:
    """处理单个 (beta, channel, source) 介子计算任务"""
    src_tag = "singlesrc" if is_single else "multisrc"
    suffix = "-singlesrc" if is_single else ""
    slice_map = SINGLE_FIT_SLICES if is_single else MULTI_FIT_SLICES

    input_dir = readin_dir / f"48x16b4.{beta}" / "Output"
    if not input_dir.exists():
        input_dir = readin_dir / f"L48T16beta4.{beta}" / "Output"

    pick_dir = out_root / f"pickdata{suffix}" / f"b4.{beta}"
    ratio_dir = out_root / f"ratio_results{suffix}" / f"b4.{beta}"
    sim_dir = out_root / f"simulateresult{suffix}" / f"b4.{beta}"

    for d in (pick_dir, ratio_dir, sim_dir):
        d.mkdir(parents=True, exist_ok=True)

    orchestrator = MesonOrchestrator()
    binsize = get_binsize(beta, channel, is_single_source=is_single)
    mappings = CHANNEL_CONFIGS[channel]

    means, errors, dr = orchestrator.execute_meson_analysis(
        input_dir=str(input_dir),
        channel_mappings=mappings,
        binsize=binsize,
        num_lines=0,
        return_folded_jk=True,
        is_single_source=is_single,
    )

    # 1. 保存关联函数
    corr_df = pl.DataFrame({
        "t": list(range(len(means))),
        "mean": means,
        "error": errors,
    })
    corr_df.write_csv(pick_dir / f"save_{channel}.csv")

    # 2. 求解有效质量
    meff_mean, meff_err = compute_meff_from_matrix(dr)
    meff_df = pl.DataFrame({
        "t": list(range(len(meff_mean))),
        "meff": meff_mean,
        "error": meff_err,
    })
    meff_df.write_csv(ratio_dir / f"meff_{channel}.csv")

    # 3. 平台拟合
    t_start = slice_map[beta][channel]
    fit_res = fit_plateau_from_matrix(
        dr=dr,
        errors=errors,
        t_start=t_start,
        t_end=FIT_SLICE_END,
        beta=beta,
        channel=channel,
        out_dir=sim_dir,
    )

    fit_res["source"] = src_tag
    return fit_res


def run_lcp_meson_pipeline(
    betas: List[str] = BETAS,
    channels: List[str] = CHANNELS,
    sources: List[bool] = [False, True],
    readin_dir: Path = DEFAULT_READIN_DIR,
    out_root: Path = OUTPUT_ROOT,
    max_workers: int = 4,
    do_plot: bool = False,
) -> List[dict]:
    """常物理线 (48^3 x 16) 介子全量端到端分析"""
    t0 = time.time()
    total_tasks = len(betas) * len(channels) * len(sources)

    print("==========================================================================")
    print(f"[LCP MESON] 启动 48^3x16 介子分析流水线: {total_tasks} 个子任务 (workers={max_workers})")
    print(f"  - Betas: {betas}")
    print(f"  - Channels: {channels}")
    print(f"  - Sources: {['single' if s else 'multi' for s in sources]}")
    print("==========================================================================")

    tasks = []
    for is_single in sources:
        for b in betas:
            for ch in channels:
                tasks.append((b, ch, is_single))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_single_task,
                beta=t[0],
                channel=t[1],
                is_single=t[2],
                readin_dir=readin_dir,
                out_root=out_root,
            ): t
            for t in tasks
        }

        for fut in as_completed(futures):
            t = futures[fut]
            try:
                res = fut.result()
                results.append(res)
                src_label = "single" if t[2] else "multi"
                print(
                    f"  [OK] Beta={t[0]} Ch={t[1]:<3} ({src_label:<6}): "
                    f"mass = {res['mass_mean']:.6f} +/- {res['mass_err']:.6f} (fit [{res['t_start']}:{res['t_end']}])"
                )
            except Exception as e:
                src_label = "single" if t[2] else "multi"
                print(f"  [FAIL] Beta={t[0]} Ch={t[1]} ({src_label}): {e}")

    # 导出汇总
    if results:
        res_df = pl.DataFrame(results).sort(["source", "beta", "channel"])
        summary_path = out_root / "all_fits_summary.csv"
        res_df.write_csv(summary_path)
        print(f"\n[OK] 介子拟合总表已保存至: {summary_path}")

    elapsed = time.time() - t0
    print(f"[OK] LCP 介子全量分析完成，耗时: {elapsed:.2f}s")

    if do_plot:
        print("\n[INFO] 正在生成 LCP 介子物理图表 (plot_lcp_meson.py)...")
        from scripts.lcp.plot_lcp_meson import main as plot_main
        plot_main()

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="【LCP 专区】常物理线 48^3x16 介子分析脚本")
    parser.add_argument("--betas", nargs="+", default=BETAS, help="Beta 编号列表")
    parser.add_argument("--channels", nargs="+", default=CHANNELS, help="物理信道列表")
    parser.add_argument(
        "--source",
        choices=["all", "multi", "single"],
        default="all",
        help="计算源类型 (all: 多源+单源, multi: 仅多源, single: 仅单源)",
    )
    parser.add_argument("--readin-dir", type=Path, default=DEFAULT_READIN_DIR, help="数据目录")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_ROOT, help="输出根目录")
    parser.add_argument("--workers", type=int, default=4, help="并发进程数")
    parser.add_argument("--plot", action="store_true", help="分析完成后自动绘图")
    args = parser.parse_args()

    sources = [False, True] if args.source == "all" else ([True] if args.source == "single" else [False])

    run_lcp_meson_pipeline(
        betas=args.betas,
        channels=args.channels,
        sources=sources,
        readin_dir=args.readin_dir,
        out_root=args.output_dir,
        max_workers=args.workers,
        do_plot=args.plot,
    )


if __name__ == "__main__":
    main()
