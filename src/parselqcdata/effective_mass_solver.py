#!/usr/bin/env python3
"""
src/parselqcdata/effective_mass_solver.py
--------------------------------------------------------------------------------
独立脚本与模块：有效质量方程求解器 (Effective Mass Solver)
- 双曲余弦方程:
    C(x) / C(x+1) = cosh(m * (x - half)) / cosh(m * (x + 1 - half))
- 支持对称点 half 参数化配置 (Ns=32->16, 36->18, 40->20, 48->24)
- 结合向量化 Newton-Raphson 极速求解与标度 fsolve 高精度兜底
- 自动计算 Jackknife 均值与统计误差
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Optional, Tuple
import warnings

import numpy as np
import polars as pl
from scipy.optimize import fsolve


def log_cosh(z: np.ndarray) -> np.ndarray:
    """数值稳定的 log(cosh(z))"""
    z = np.asarray(z, dtype=np.float64)
    return np.logaddexp(z, -z)


def cosh_ratio(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """数值稳定的 cosh(a) / cosh(b)"""
    return np.exp(log_cosh(a) - log_cosh(b))


def solve_effective_mass_one_point_centered(y: float, x: int, half: float) -> float:
    """
    求解单个数据点的有效质量方程 (对称点 half = Ns/2):
        C(x) / C(x+1) = cosh(m*(x-half)) / cosh(m*(x+1-half))
    """
    if not np.isfinite(y):
        return np.nan

    if x < half and y <= 1.0:
        return 0.0 if y == 1.0 else np.nan

    a = float(x) - float(half)
    b = float(x) + 1.0 - float(half)

    def equation(m):
        return cosh_ratio(m * a, m * b) - y

    try:
        sol = fsolve(equation, x0=0.1, xtol=1e-12, maxfev=2000)
        root = float(np.ravel(sol)[0])
    except Exception:
        return np.nan

    if not np.isfinite(root):
        return np.nan
    residual = abs(float(np.ravel(cosh_ratio(root * a, root * b))[0]) - y)
    if residual > 1e-8 * max(1.0, abs(y)):
        return np.nan
    return abs(root)


def solve_effective_mass_one_point(y: float, x: int, half: float = 24.0) -> float:
    """求解单个数据点的有效质量方程，默认对称点 half = 24.0 (兼容旧接口)"""
    return solve_effective_mass_one_point_centered(y, x, half=half)


def _newton_column(
    y_col: np.ndarray, x: int, half: float, max_iter: int = 200, tol: float = 1e-14
) -> Tuple[np.ndarray, np.ndarray]:
    """对同一时间切片 x 上的整列 Jackknife 比值做向量化 Newton 迭代。"""
    a = float(x) - float(half)
    b = float(x) + 1.0 - float(half)

    m = np.full(y_col.shape, 0.1, dtype=np.float64)
    converged = np.zeros(y_col.shape, dtype=bool)
    active = np.isfinite(y_col)

    for _ in range(max_iter):
        if not np.any(active):
            break
        idx = np.nonzero(active)[0]
        mm = m[idx]
        am, bm = mm * a, mm * b
        ratio = cosh_ratio(am, bm)
        f = ratio - y_col[idx]
        dlog = a * np.tanh(am) - b * np.tanh(bm)
        deriv = ratio * dlog
        safe = np.abs(deriv) > 1e-300
        step = np.zeros_like(mm)
        step[safe] = f[safe] / deriv[safe]
        step = np.clip(step, -5.0, 5.0)
        mm_new = np.clip(mm - step, -50.0, 50.0)
        m[idx] = np.abs(mm_new)
        resid = np.abs(cosh_ratio(mm_new * a, mm_new * b) - y_col[idx])
        done = (resid <= tol * np.maximum(1.0, np.abs(y_col[idx]))) & (np.abs(step) <= 1e-14)
        converged[idx[done]] = True
        active[idx[done]] = False

    return m, converged


def compute_effective_mass_matrix_centered(
    dr_matrix: np.ndarray, half: float, vectorized: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    输入: dr_matrix (N_lines 行 x J 列 Jackknife 折叠样本), half = Ns/2
    返回: (meff_matrix, meff_mean, meff_err)
    """
    dr_matrix = np.asarray(dr_matrix, dtype=np.float64)
    nrow, ncol = dr_matrix.shape
    meff = np.full((nrow, ncol), np.nan, dtype=np.float64)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(nrow - 1):
            if ncol == 0:
                continue
            ratios = np.full(ncol, np.nan, dtype=np.float64)
            valid = np.abs(dr_matrix[i + 1, :]) > 1e-16
            ratios[valid] = dr_matrix[i, valid] / dr_matrix[i + 1, valid]

            if vectorized and np.any(np.isfinite(ratios)):
                solved, ok = _newton_column(ratios, i, half)
                meff[i, ok] = np.abs(solved[ok])
                fallback = np.nonzero(np.isfinite(ratios) & ~ok)[0]
            else:
                fallback = np.nonzero(np.isfinite(ratios))[0]

            for j in fallback:
                meff[i, j] = solve_effective_mass_one_point_centered(ratios[j], i, half)

    meff = np.abs(meff)
    meff_mean = np.full(nrow, np.nan, dtype=np.float64)
    meff_err = np.full(nrow, np.nan, dtype=np.float64)
    min_k = max(2, int(ncol * 0.8))
    for i in range(nrow):
        finite = meff[i, np.isfinite(meff[i])]
        k = len(finite)
        if k >= min_k:
            m = float(np.mean(finite))
            e = float(np.sqrt((k - 1) * np.sum((finite - m) ** 2) / k))
            meff_mean[i] = m
            meff_err[i] = e

    return meff, meff_mean, meff_err


def compute_effective_mass_matrix(
    dr_matrix: np.ndarray, half: Optional[float] = None, vectorized: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """输入: dr_matrix (N_lines x J 列样本矩阵)，返回: (meff_matrix, meff_mean, meff_err)"""
    if half is None:
        nrow = dr_matrix.shape[0]
        half = 24.0 if nrow in (47, 48) else (nrow + 1) / 2.0
    return compute_effective_mass_matrix_centered(dr_matrix, half=half, vectorized=vectorized)


def main() -> None:
    parser = argparse.ArgumentParser(description="有效质量求解脚本")
    parser.add_argument("--ratio", type=float, default=None, help="单个切片比值 C(x)/C(x+1)")
    parser.add_argument("--x", type=int, default=10, help="切片编号 x (默认: 10)")
    parser.add_argument("--half", type=float, default=24.0, help="对称点 half = Ns/2 (默认: 24.0)")
    parser.add_argument("--matrix-csv", type=str, default=None, help="从 CSV 矩阵文件批量计算")
    args = parser.parse_args()

    if args.matrix_csv:
        path = Path(args.matrix_csv)
        if not path.is_file():
            print(f"[ERROR] 文件不存在: {path}", file=sys.stderr)
            sys.exit(1)
        mat = pl.read_csv(path).to_numpy()
        _, mean, err = compute_effective_mass_matrix_centered(mat, half=args.half)
        print(f"--- 批量有效质量曲线 (half={args.half}) ---")
        for i, (m, e) in enumerate(zip(mean, err)):
            if np.isfinite(m):
                print(f"  x={i:2d}: meff = {m:.8f} +/- {e:.8f}")
            else:
                print(f"  x={i:2d}: meff = NaN")
    elif args.ratio is not None:
        sol = solve_effective_mass_one_point_centered(args.ratio, args.x, args.half)
        print(f"比值 C({args.x})/C({args.x+1}) = {args.ratio}, half = {args.half}")
        print(f"有效质量解: {sol:.8f}")
    else:
        print("[INFO] 请提供 --ratio 或 --matrix-csv 参数进行求解。")


if __name__ == "__main__":
    main()
