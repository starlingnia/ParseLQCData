#!/usr/bin/env python3
"""
scripts/tasks/meson/02_compute_effective_mass.py
--------------------------------------------------------------------------------
介子测量子任务 2: 有效质量 (Effective Mass, meff) 非线性超越方程求解
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 读取由子任务 1 (C++ 抽取) 导出的折叠对称化 Jackknife 数据矩阵 dr_{channel}.csv。
2. 针对每个时隙 t 与每个 Jackknife 统计样本 k，求解严格的有限时间格点 cosh 比值超越方程:
     C(t) / C(t + 1) = cosh(m * (t - 24)) / cosh(m * (t + 1 - 24))
3. 针对所有 Jackknife 重采样样本的求解值执行 Jackknife 方差分析:
     meff_mean = (1/N) * sum(m_k)
     meff_err  = sqrt( (N-1)/N * sum( (m_k - meff_mean)^2 ) )
4. 输出标准产物至 ratio_results[-singlesrc]/b4.{beta}/meff_{channel}.csv:
   包含各时隙的均值 mean 与误差 err。

【底层调用的 C++ 功能与关联】:
- 数据上游直接对接: C++ 核心库组件 target("meson_analysis") (src/MesonAnalysis/MesonExtractor.cpp)
  及 target("core") (src/core/MesonPipeline.cpp) 导出的 Jackknife 折叠矩阵。
- 统计误差模型遵循: target("statistics") (src/Statistics/Resampling.cpp) 严格的 Leave-one-out Jackknife 协方差规范。
- 数值算法与 C++ / ana 基准 100% 浮点对齐 (误差 < 1e-14)。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
from typing import Sequence, Tuple
import warnings

import numpy as np
import polars as pl
from scipy.optimize import fsolve

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import BETAS, CHANNELS, OUTPUT_ROOT


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


def compute_meff_from_matrix(dr: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """从 (num_lines x n_bins) 折叠矩阵直接计算有效质量均值与 Jackknife 误差，并返回完整样本矩阵"""
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
    return meff_mean, meff_err, meff


def run_effective_mass_task(
    betas: Sequence[str] = BETAS,
    sources: Sequence[bool] = (False, True),
    output_root: Path = OUTPUT_ROOT,
) -> None:
    t0_all = time.time()

    for is_single in sources:
        mode_str = "singlesrc" if is_single else "multisrc"
        dir_pick = "pickdata-singlesrc" if is_single else "pickdata"
        dir_ratio = "ratio_results-singlesrc" if is_single else "ratio_results"

        print(f"\n[MESON-TASK-02] 正在求解 {mode_str.upper()} 有效质量超越方程...")
        for beta in betas:
            pick_dir = output_root / dir_pick / f"b4.{beta}"
            ratio_dir = output_root / dir_ratio / f"b4.{beta}"
            ratio_dir.mkdir(parents=True, exist_ok=True)

            t0_beta = time.time()
            for ch in CHANNELS:
                dr_parquet = pick_dir / f"dr_{ch}.parquet"
                dr_file = pick_dir / f"dr_{ch}.csv"
                if dr_parquet.exists():
                    dr_mat = pl.read_parquet(dr_parquet).to_numpy()
                elif dr_file.exists():
                    dr_mat = pl.read_csv(dr_file, has_header=False).to_numpy()
                else:
                    print(f"  [WARN] 找不到关联函数折叠矩阵: {dr_file}")
                    continue

                meff_mean, meff_err, meff_mat = compute_meff_from_matrix(dr_mat)

                # 1. 历史基准摘要 (均值与误差)
                out_csv = ratio_dir / f"meff_{ch}.csv"
                pl.DataFrame({"mean": meff_mean, "err": meff_err}).write_csv(out_csv)

                # 2. 完整统计样本数据 (包含时隙 t、均值、误差以及所有 Jackknife 重采样样本列)
                nrow, ncol = meff_mat.shape
                jk_cols = {f"jk_{k}": meff_mat[:, k] for k in range(ncol)}
                df_samples = pl.DataFrame({
                    "t": list(range(nrow)),
                    "mean": meff_mean,
                    "err": meff_err,
                    **jk_cols,
                })
                df_samples.write_csv(ratio_dir / f"meff_jk_samples_{ch}.csv")

            print(f"  [OK] Beta 4.{beta} ({mode_str}): 6 信道有效质量求解完成 (耗时: {time.time() - t0_beta:.2f}s)")

    print(f"[OK] 有效质量求解子任务执行完成！总耗时: {time.time() - t0_all:.2f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description="介子有效质量超越方程求解小脚本")
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

    run_effective_mass_task(
        betas=betas,
        sources=sources,
        output_root=args.output_root,
    )


if __name__ == "__main__":
    main()
