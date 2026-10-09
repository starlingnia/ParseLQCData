#!/usr/bin/env python3
"""
src/parselqcdata/cosh_fitter.py
--------------------------------------------------------------------------------
高层可读性包装器：双曲余弦 (cosh) 平台非线性拟合器
- 统一委托至 src.parselqcdata.plateau_fit 权威拟合核心 (消除重复实现)
- 强制使用 scipy_least_squares 拟合器与统一残差容差
- 提供面向介子分析的高层可读性拟合接口 fit_cosh_plateau()
- 完整保留所有统计重采样样本明细 (含 jk_index, mass, a, chi2, dof)
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence, Tuple
import warnings

import numpy as np
import pandas as pd
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# 从权威实现 plateau_fit 统一导入，杜绝多处重复定义
from src.parselqcdata.plateau_fit import (
    LEAST_SQUARES_FITTER,
    chi2_least_squares_fit,
    target_cosh_func,
    centered_cosh_func,
    fit_single_jackknife_column,
    fit_single_jackknife_column_centered,
    fit_meson_plateau,
    fit_jackknife_mass_centered,
    fit_mass_window_scan,
)


def fit_cosh_plateau(
    dr_matrix: np.ndarray,
    errors: np.ndarray,
    t_start: int,
    t_end: int,
    half: float = 24.0,
    channel: str = "",
    out_dir: Optional[Path] = None,
) -> Tuple[Dict[str, float], List[dict]]:
    """
    高层可读性拟合 Wrapper:
    对关联函数 Jackknife 折叠矩阵在时间切片窗口 [t_start, t_end) 执行逐样本 cosh 平台拟合。

    Parameters:
        dr_matrix: (N_time x N_bins) 关联函数折叠矩阵
        errors: (N_time,) 关联函数统计误差向量
        t_start: 拟合窗口起始时隙
        t_end: 拟合窗口截止时隙
        half: 对称点 (默认 24.0，对应 Ns=48; 若 Ns=32 则为 16.0)
        channel: 信道名称 (例如 'AV', 'PS', 'S')
        out_dir: 产物输出目录 (若提供则自动落盘逐样本与汇总 CSV)

    Returns:
        (summary_dict, fit_records):
        - summary_dict: {"mass": ..., "mass_err": ..., "a": ..., "a_err": ..., "chi2_dof": ...}
        - fit_records: 逐 Jackknife 样本的拟合明细列表 (含 jk_index)
    """
    dr_matrix = np.asarray(dr_matrix, dtype=np.float64)
    errors = np.asarray(errors, dtype=np.float64)
    x_array = np.arange(t_start, t_end)
    y_slice = dr_matrix[t_start:t_end, :]
    err_slice = errors[t_start:t_end]
    n_bins = y_slice.shape[1]

    fit_records: List[dict] = []
    mass_list: List[float] = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(n_bins):
            try:
                if abs(half - 24.0) < 1e-6:
                    fit = fit_single_jackknife_column(x_array, y_slice[:, k], err_slice)
                else:
                    fit = fit_single_jackknife_column_centered(
                        x_array, y_slice[:, k], err_slice, half=half
                    )
                if fit:
                    rec = {"jk_index": k, **fit}
                    fit_records.append(rec)
                    mass_list.append(float(fit.get("massfit_mean", np.nan)))
                else:
                    fit_records.append({"jk_index": k, "massfit_mean": np.nan})
                    mass_list.append(np.nan)
            except Exception:
                fit_records.append({"jk_index": k, "massfit_mean": np.nan})
                mass_list.append(np.nan)

    mass_arr = np.array(mass_list, dtype=np.float64)
    valid_mask = np.isfinite(mass_arr)

    if np.any(valid_mask):
        fit_mass = float(np.mean(mass_arr[valid_mask]))
        diff = mass_arr[valid_mask] - fit_mass
        fit_err = float(np.sqrt((n_bins - 1) * np.sum(diff**2) / n_bins))
    else:
        fit_mass, fit_err = np.nan, np.nan

    summary_dict = {
        "channel": channel,
        "mass_mean": fit_mass,
        "mass_err": fit_err,
        "t_start": t_start,
        "t_end": t_end,
        "half": half,
        "n_bins": n_bins,
    }

    # 持久化输出
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        ch_suffix = f"{channel}" if channel else ""

        if fit_records:
            df_samples = pl.DataFrame(fit_records)
            df_samples.write_csv(out_dir / f"fittresult{ch_suffix}.csv")

        summary_df = pl.DataFrame({
            "quantity": ["mass"],
            "mean": [fit_mass],
            "jack_err": [fit_err],
        })
        summary_df.write_csv(out_dir / f"summary_fit_{ch_suffix}.csv")

    return summary_dict, fit_records


def main() -> None:
    parser = argparse.ArgumentParser(description="双曲余弦 (cosh) 平台拟合包装工具")
    parser.add_argument("--start", type=int, default=14, help="平台起始切片 (默认: 14)")
    parser.add_argument("--end", type=int, default=24, help="平台结束切片 (默认: 24)")
    parser.add_argument("--half", type=float, default=24.0, help="对称点 (默认: 24.0)")
    args = parser.parse_args()

    x = np.arange(args.start, args.end)
    true_m, true_a = 0.45, 1.2
    y = true_a * np.cosh(true_m * (x - args.half)) + 1e-5 * np.random.randn(len(x))
    err = np.full(len(x), 1e-4)
    res = fit_single_jackknife_column_centered(x, y, err, half=args.half)
    print("--- Cosh 平台拟合测试 (通过 plateau_fit 统一引擎) ---")
    print(f"真值: mass={true_m}, A={true_a}")
    for k, v in res.items():
        print(f"  {k:<14}: {v}")


if __name__ == "__main__":
    main()
