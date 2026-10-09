#!/usr/bin/env python3
"""
src/parselqcdata/cosh_fitter.py
--------------------------------------------------------------------------------
独立脚本与模块：双曲余弦 (cosh) 平台非线性拟合器
- 拟合形式: f(x) = a * cosh(m * (x - half))
- 统一强制 chi2 最小二乘拟合器 (scipy_least_squares, tol=1e-15, maxit=10000)
- 逐 Jackknife 样本列独立拟合与全量矩阵重采样统计
- 候选窗口全量扫描与 chi2/dof 评估
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence, Tuple
import warnings

import gvar as gv
import lsqfit
import numpy as np
import pandas as pd
import polars as pl

# ==============================================================================
# 拟合器统一入口
# ------------------------------------------------------------------------------
# lsqfit.nonlinear_fit 的 fitter 默认值是 None, 含义是 "按环境自动挑":
# 装了 GSL 时用 'gsl_multifit', 没装时才退到 'scipy_least_squares'。
# 这会让同一份代码在不同机器上走出不同结果 (甚至退化), 因此本仓库**一律显式**
# 指定 chi2 最小二乘拟合器, 所有拟合必须经过 chi2_least_squares_fit()。
# ==============================================================================
from src.parselqcdata.plateau_fit import LEAST_SQUARES_FITTER, chi2_least_squares_fit



def target_cosh_func(x, p):
    """默认对称点 center = 24.0 的 cosh 拟合函数"""
    return p["a"] * np.cosh(p["m"] * (x - 24.0))


def centered_cosh_func(half: float):
    """生成对称点为 half 的 cosh 拟合函数: f(x) = a * cosh(m * (x - half))"""
    def _fcn(x, p):
        return p["a"] * np.cosh(p["m"] * (x - half))
    return _fcn


def _initial_guess(x_data: np.ndarray, y_val: np.ndarray, half: float) -> dict:
    """按首点量级给出初值, 避免 cosh(x - half) 很大时 a 的初值离谱"""
    x0 = float(np.asarray(x_data).ravel()[0])
    y0 = float(np.asarray(y_val).ravel()[0])
    m0 = 0.3
    c0 = np.cosh(m0 * (x0 - half))
    a0 = y0 / c0 if np.isfinite(c0) and c0 != 0.0 else 1.0
    if not np.isfinite(a0) or a0 == 0.0:
        a0 = 1.0
    return {"a": a0, "m": m0}


def fit_single_jackknife_column(
    x_data: np.ndarray, y_val: np.ndarray, y_err: np.ndarray
) -> dict:
    """对单个 Jackknife 样本列执行非线性拟合 (默认 half = 24.0)"""
    y_gv = gv.gvar(y_val, y_err)
    try:
        fit = chi2_least_squares_fit(
            data=(x_data, y_gv),
            fcn=target_cosh_func,
            p0={"a": 1.0, "m": 0.5},
        )
        return {
            "chi2": float(fit.chi2),
            "dof": fit.dof,
            "massfit_mean": abs(float(fit.p["m"].mean)),
            "massfit_err": float(fit.p["m"].sdev),
            "fita": float(fit.p["a"].mean),
            "fita_err": float(fit.p["a"].sdev),
            "chi2_dof": float(fit.chi2 / (fit.dof - 1)) if fit.dof > 1 else float("inf"),
        }
    except Exception:
        return {}


def fit_single_jackknife_column_centered(
    x_data: np.ndarray, y_val: np.ndarray, y_err: np.ndarray, half: float
) -> dict:
    """单个 Jackknife 样本列的非线性 cosh 拟合 (对称点 half = Ns/2)"""
    y_gv = gv.gvar(y_val, y_err)
    try:
        fit = chi2_least_squares_fit(
            data=(x_data, y_gv),
            fcn=centered_cosh_func(half),
            p0=_initial_guess(x_data, y_val, half),
        )
        return {
            "chi2": float(fit.chi2),
            "dof": fit.dof,
            "massfit_mean": abs(float(fit.p["m"].mean)),
            "massfit_err": float(fit.p["m"].sdev),
            "fita": float(fit.p["a"].mean),
            "fita_err": float(fit.p["a"].sdev),
            "chi2_dof": float(fit.chi2 / (fit.dof - 1)) if fit.dof > 1 else float("inf"),
        }
    except Exception:
        return {}


def fit_meson_plateau(
    df_sym: pd.DataFrame,
    err_vals: np.ndarray,
    slice_start: int,
    slice_end: int = 25,
) -> Tuple[pd.DataFrame, Optional[dict]]:
    """对对称折叠矩阵全部 Jackknife 样本进行切片平台拟合 (兼容旧接口)"""
    x_array = np.arange(48)[slice_start:slice_end]
    fits = []

    for col in df_sym.columns:
        y_val = df_sym[col].to_numpy()[slice_start:slice_end]
        y_err = err_vals[slice_start:slice_end]
        f = fit_single_jackknife_column(x_array, y_val, y_err)
        if f and f.get("chi2_dof", float("inf")) <= 1000000.0:
            fits.append(f)

    df_res = pd.DataFrame(fits)
    summary = None

    if not df_res.empty:
        n_samples = len(df_res)
        means = df_res.mean()
        diffs = df_res - means
        sum_sq = (diffs**2).sum()
        jack_err = np.sqrt(((n_samples - 1) / n_samples) * sum_sq)
        summary = {
            "mass_mean": float(means["massfit_mean"]),
            "mass_err": float(jack_err["massfit_mean"]),
            "a_mean": float(means["fita"]),
            "a_err": float(jack_err["fita"]),
            "chi2_dof": float(means["chi2_dof"]),
            "n_samples": int(n_samples),
        }

    return df_res, summary


def fit_mass_window_scan(
    x_data: np.ndarray,
    y_mean: np.ndarray,
    y_err: np.ndarray,
    windows: List[Tuple[int, int]],
    half: float,
) -> List[dict]:
    """
    平台窗口扫描 (对均值关联函数逐窗口拟合), 用于挑选最佳平台窗口。
    返回每个窗口的质量/chi2 记录, 失败窗口以 NaN 占位。
    """
    records: List[dict] = []
    x_data = np.asarray(x_data)
    y_mean = np.asarray(y_mean)
    y_err = np.asarray(y_err)

    for start, end in windows:
        rec = {
            "x_start": int(start), "x_end": int(end), "n_points": int(end - start),
            "mass": np.nan, "mass_err": np.nan, "a": np.nan, "a_err": np.nan,
            "chi2": np.nan, "dof": np.nan, "chi2_dof": np.nan, "ok": False,
        }
        if end - start < 2:
            records.append(rec)
            continue
        fit = fit_single_jackknife_column_centered(
            x_data[start:end], y_mean[start:end], y_err[start:end], half
        )
        if fit:
            rec.update({
                "mass": float(fit["massfit_mean"]),
                "mass_err": float(fit["massfit_err"]),
                "a": float(fit["fita"]),
                "a_err": float(fit["fita_err"]),
                "chi2": float(fit["chi2"]),
                "dof": float(fit["dof"]),
                "chi2_dof": float(fit["chi2_dof"]),
                "ok": True,
            })
        records.append(rec)
    return records


def fit_jackknife_mass_centered(
    x_data: np.ndarray,
    jk_matrix: np.ndarray,
    y_err: np.ndarray,
    start: int,
    end: int,
    half: float,
) -> dict:
    """
    对 Jackknife 折叠矩阵的每一列在窗口 [start, end) 内做 cosh 平台拟合,
    按仓库既有约定汇总 (均值 + Jackknife 误差 sqrt((J-1) * sum(d^2) / J))。
    """
    x_array = np.asarray(x_data)[start:end]
    err_slice = np.asarray(y_err)[start:end]
    jk_slice = np.asarray(jk_matrix)[start:end, :]
    n_samples = jk_slice.shape[1]

    masses: List[float] = []
    amps: List[float] = []
    chi2_dofs: List[float] = []
    per_sample: List[dict] = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(n_samples):
            fit = fit_single_jackknife_column_centered(
                x_array, jk_slice[:, k], err_slice, half
            )
            if fit:
                masses.append(float(fit["massfit_mean"]))
                amps.append(float(fit["fita"]))
                if np.isfinite(fit["chi2_dof"]):
                    chi2_dofs.append(float(fit["chi2_dof"]))
                per_sample.append({"sample": k, **fit})
            else:
                per_sample.append({"sample": k})

    masses_arr = np.asarray(masses, dtype=np.float64)
    amps_arr = np.asarray(amps, dtype=np.float64)

    result = {
        "mass": np.nan, "mass_err": np.nan,
        "a": np.nan, "a_err": np.nan,
        "chi2_dof": np.nan,
        "n_samples": int(n_samples), "n_fitted": int(masses_arr.size),
        "x_start": int(start), "x_end": int(end), "half": float(half),
    }
    if masses_arr.size:
        mass_mean = float(np.mean(masses_arr))
        a_mean = float(np.mean(amps_arr)) if amps_arr.size else np.nan
        result["mass"] = mass_mean
        result["a"] = a_mean
        result["mass_err"] = float(np.sqrt((masses_arr.size - 1) * np.sum((masses_arr - mass_mean) ** 2) / masses_arr.size))
        if amps_arr.size:
            result["a_err"] = float(np.sqrt((amps_arr.size - 1) * np.sum((amps_arr - a_mean) ** 2) / amps_arr.size))
        if chi2_dofs:
            result["chi2_dof"] = float(np.mean(chi2_dofs))
    result["per_sample"] = per_sample
    return result


def select_best_window(
    scan_results: Sequence[dict],
    chi2_max: float = 3.0,
    min_points: int = 3,
) -> Optional[dict]:
    """从窗口扫描结果中挑选 chi2/dof <= chi2_max 且最平坦的窗口。"""
    valid = [
        r for r in scan_results
        if r.get("ok") and np.isfinite(r.get("chi2_dof", np.nan)) and r.get("n_points", 0) >= min_points
    ]
    if not valid:
        return None
    acceptable = [r for r in valid if r["chi2_dof"] <= chi2_max]
    pool = acceptable if acceptable else valid
    return min(pool, key=lambda r: r["chi2_dof"])


def main() -> None:
    parser = argparse.ArgumentParser(description="双曲余弦 (cosh) 平台拟合脚本")
    parser.add_argument("--start", type=int, default=14, help="平台起始切片 (默认: 14)")
    parser.add_argument("--end", type=int, default=24, help="平台结束切片 (默认: 24)")
    parser.add_argument("--half", type=float, default=24.0, help="对称点 (默认: 24.0)")
    args = parser.parse_args()

    x = np.arange(args.start, args.end)
    true_m, true_a = 0.45, 1.2
    y = true_a * np.cosh(true_m * (x - args.half)) + 1e-5 * np.random.randn(len(x))
    err = np.full(len(x), 1e-4)
    res = fit_single_jackknife_column_centered(x, y, err, half=args.half)
    print("--- Cosh 平台拟合测试 ---")
    print(f"真值: mass={true_m}, A={true_a}")
    for k, v in res.items():
        print(f"  {k:<14}: {v}")


if __name__ == "__main__":
    main()
