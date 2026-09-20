"""
src/parselqcdata/effective_mass.py
--------------------------------------------------------------------------------
Python 整合层：Effective Mass (有效质量) 管道
- 求解双曲余弦方程: C(t)/C(t+1) = cosh(m*(t-24))/cosh(m*(t-23))
- 并发计算各 Jackknife 样本并得出均值与统计误差
--------------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from typing import Tuple


def solve_effective_mass_one_point(y: float, x: int) -> float:
    """求解单个数据点的有效质量方程 (严格匹配 ana/src/effectivemass.py 实现)"""
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


def compute_effective_mass_matrix(dr_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    输入: dr_matrix (47 行 x N_cols 样本)
    返回: (meff_matrix, meff_mean, meff_err)
    """
    nrow, ncol = dr_matrix.shape
    meff = np.zeros((nrow, ncol), dtype=np.float64)

    for i in range(nrow):
        for j in range(ncol):
            meff[i, j] = solve_effective_mass_one_point(dr_matrix[i, j], i)

    meff_mean = np.nanmean(meff, axis=1)
    diffs2 = (meff - meff_mean[:, None]) ** 2
    meff_err = np.sqrt((ncol - 1) * np.nansum(diffs2, axis=1) / ncol)

    return meff, meff_mean, meff_err
