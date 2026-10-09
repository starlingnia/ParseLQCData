#!/usr/bin/env python3
"""
scripts/tasks/meson/03_fit_plateau.py
--------------------------------------------------------------------------------
介子测量子任务 3: 基态平台 cosh 拟合与 Jackknife 统计推断
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 读取由子任务 1 (C++ 抽取) 导出的折叠对称化矩阵 dr_{channel}.csv 与误差 err_{channel}.csv。
2. 依据 docs/physics_setup.py 中针对各 Beta 与信道设定的最佳物理饱和时间窗口 [t_start, t_end]，
   截取时间切片区间。
3. 逐 Jackknife 重采样样本执行 cosh 理论模型非线性最小二乘拟合:
     C(t) = A * cosh(m * (t - Nt/2))
4. 汇总 Jackknife 拟合参数，输出:
   - simulateresult[-singlesrc]/b4.{beta}/fittresult{channel}.csv : 逐样本拟合明细 (mass, amp, chi2)
   - simulateresult[-singlesrc]/b4.{beta}/summary_fit_{channel}.csv : 该信道统计平均质量与误差
   - simulateresult[-singlesrc]/all_fits_summary.csv : 全 Beta、全信道汇总大表

【底层调用的 C++ 功能与关联】:
- 数据依赖: C++ 抽取与折叠产物 (target("meson_analysis") + target("core"))。
- 统计误差与协方差: target("statistics") (src/Statistics/Resampling.cpp) 的 Jackknife 方差推断。
- 拟合引擎: 基于 src.parselqcdata.plateau_fit 极速优化器与卡方检验。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
from typing import Sequence
import warnings

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    BETAS,
    CHANNELS,
    FIT_SLICE_END,
    MULTI_FIT_SLICES,
    OUTPUT_ROOT,
    SINGLE_FIT_SLICES,
)
from src.parselqcdata import fit_cosh_plateau


def fit_plateau_from_matrix(
    dr: np.ndarray,
    errors: np.ndarray,
    t_start: int,
    t_end: int,
    beta: str,
    channel: str,
    out_dir: Path | None = None,
) -> dict:
    """调用高层拟合 wrapper 执行逐 Jackknife 样本 cosh 平台拟合并落盘"""
    summary, _ = fit_cosh_plateau(
        dr_matrix=dr,
        errors=errors,
        t_start=t_start,
        t_end=t_end,
        half=24.0,
        channel=channel,
        out_dir=out_dir,
    )
    return {
        "beta": beta,
        "channel": channel,
        "mass_mean": summary["mass_mean"],
        "mass_err": summary["mass_err"],
        "t_start": t_start,
        "t_end": t_end,
    }


def run_fit_task(
    betas: Sequence[str] = BETAS,
    sources: Sequence[bool] = (False, True),
    output_root: Path = OUTPUT_ROOT,
) -> None:
    t0_all = time.time()

    for is_single in sources:
        mode_str = "singlesrc" if is_single else "multisrc"
        dir_pick = "pickdata-singlesrc" if is_single else "pickdata"
        dir_sim = "simulateresult-singlesrc" if is_single else "simulateresult"
        slice_map = SINGLE_FIT_SLICES if is_single else MULTI_FIT_SLICES

        print(f"\n[MESON-TASK-03] 正在对 {mode_str.upper()} 执行平台 cosh 拟合...")
        all_fits = []

        for beta in betas:
            pick_dir = output_root / dir_pick / f"b4.{beta}"
            sim_dir = output_root / dir_sim / f"b4.{beta}"
            sim_dir.mkdir(parents=True, exist_ok=True)

            t0_beta = time.time()
            for ch in CHANNELS:
                dr_parquet = pick_dir / f"dr_{ch}.parquet"
                dr_file = pick_dir / f"dr_{ch}.csv"
                err_file = pick_dir / f"err_{ch}.csv"
                if not err_file.exists():
                    print(f"  [WARN] 找不到拟合所需误差文件: {err_file}")
                    continue

                if dr_parquet.exists():
                    dr_mat = pl.read_parquet(dr_parquet).to_numpy()
                elif dr_file.exists():
                    dr_mat = pl.read_csv(dr_file, has_header=False).to_numpy()
                else:
                    print(f"  [WARN] 找不到拟合所需关联函数: {dr_file}")
                    continue

                err_arr = pl.read_csv(err_file, has_header=False).to_numpy().flatten()
                t_start = slice_map[beta][ch]

                rec = fit_plateau_from_matrix(
                    dr=dr_mat,
                    errors=err_arr,
                    t_start=t_start,
                    t_end=FIT_SLICE_END,
                    beta=beta,
                    channel=ch,
                    out_dir=sim_dir,
                )
                all_fits.append(rec)

            print(f"  [OK] Beta 4.{beta} ({mode_str}): 6 信道平台拟合完成 (耗时: {time.time() - t0_beta:.2f}s)")

        if all_fits:
            new_df = pl.DataFrame(all_fits)
            summary_path = output_root / dir_sim / "all_fits_summary.csv"
            if summary_path.exists() and len(betas) < len(BETAS):
                try:
                    old_df = pl.read_csv(summary_path)
                    filtered_df = old_df.filter(
                        ~pl.col("beta").cast(pl.String).is_in([str(b) for b in betas])
                    )
                    summary_df = pl.concat([filtered_df, new_df]).sort(["beta", "channel"])
                except Exception:
                    summary_df = new_df.sort(["beta", "channel"])
            else:
                summary_df = new_df.sort(["beta", "channel"])
            summary_df.write_csv(summary_path)
            print(f"  -> 全量拟合总表已保存: {summary_path}")

    print(f"[OK] 平台 cosh 拟合子任务执行完成！总耗时: {time.time() - t0_all:.2f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description="介子基态平台 cosh 拟合小脚本")
    parser.add_argument(
        "--source",
        choices=["all", "multi", "single"],
        default="all",
        help="数据源类型: multi, single 或 all",
    )
    parser.add_argument(
        "--beta",
        type=str,
        default="all",
        help="指定 Beta (例如 17 或 4.17)，或 'all'",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="输出根目录",
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

    run_fit_task(
        betas=betas,
        sources=sources,
        output_root=args.output_root,
    )


if __name__ == "__main__":
    main()
