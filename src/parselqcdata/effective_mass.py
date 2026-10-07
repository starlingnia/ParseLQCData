"""
src/parselqcdata/effective_mass.py
--------------------------------------------------------------------------------
Python 整合层：Effective Mass (有效质量) 管道
- 求解双曲余弦方程: C(t)/C(t+1) = cosh(m*(t-24))/cosh(m*(t-23))
- 并发计算各 Jackknife 样本并得出均值与统计误差
--------------------------------------------------------------------------------
"""

import warnings

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from typing import Tuple


def solve_effective_mass_one_point(y: float, x: int, half: float = 24.0) -> float:
    """
    求解单个数据点的有效质量方程:
        C(x) / C(x+1) = cosh(m*(x-half)) / cosh(m*(x+1-half))
    默认对称点 half = 24.0，亦可传入 Ns/2。
    """
    return solve_effective_mass_one_point_centered(y, x, half=half)


def compute_effective_mass_matrix(
    dr_matrix: np.ndarray, half: Optional[float] = None, vectorized: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    输入: dr_matrix (N_lines x J 列样本矩阵)
    返回: (meff_matrix, meff_mean, meff_err)
    """
    if half is None:
        nrow = dr_matrix.shape[0]
        half = 24.0 if nrow in (47, 48) else (nrow + 1) / 2.0
    return compute_effective_mass_matrix_centered(dr_matrix, half=half, vectorized=vectorized)


def detect_plateau_window(
    meff_matrix: np.ndarray,
    end: int,
    chi2_dof_max: float = 2.0,
    min_points: int = 3,
    min_start: int = 1,
) -> dict:
    """
    基于 Jackknife 样本的有效质量平台自动检测。

    判据: 窗口 [s, end) 内所有相邻差分 D(x) = meff(x+1) - meff(x) 在 Jackknife
    误差下与 0 一致, 即 flatness chi2/dof <= chi2_dof_max。因为窗口随 s 减小而变长、
    chi2/dof 单调变差, 第一个通过的 s 就是 "离对称点最近的平坦窗口" (最长平台)。

    返回: {start, end, chi2_dof, n_diffs, n_points, mass_est, mass_err_est, ok}
    """
    meff_matrix = np.asarray(meff_matrix, dtype=np.float64)
    nrow, ncol = meff_matrix.shape
    end = int(min(end, nrow))
    result = {
        "start": None, "end": end, "chi2_dof": np.nan, "n_diffs": 0,
        "n_points": 0, "mass_est": np.nan, "mass_err_est": np.nan, "ok": False,
    }
    if ncol == 0 or end - int(min_start) < int(min_points):
        return result

    # 相邻差分与其 Jackknife 误差
    D = meff_matrix[1:end, :] - meff_matrix[: end - 1, :]
    finite = np.isfinite(D)
    cnt = finite.sum(axis=1)
    dmean = np.where(cnt > 0, np.nansum(D, axis=1) / np.maximum(cnt, 1), np.nan)
    resid = np.where(finite, D - np.where(np.isfinite(dmean), dmean, np.nan)[:, None], 0.0)
    dvar = np.sum(resid ** 2, axis=1)
    derr = np.sqrt(np.maximum(cnt - 1, 0) * dvar / np.maximum(cnt, 1))

    for s in range(int(min_start), end - int(min_points) + 1):
        idx = slice(s, end - 1)          # 差分编号 x = s ... end-2
        dd, de = dmean[idx], derr[idx]
        ok = np.isfinite(dd) & np.isfinite(de) & (de > 0.0)
        n_diff = int(ok.sum())
        if n_diff < 2:
            continue
        chi2 = float(np.sum((dd[ok] / de[ok]) ** 2))
        chi2_dof = chi2 / n_diff
        if chi2_dof <= chi2_dof_max:
            window = meff_matrix[s:end, :]
            m_mean = np.nanmean(window)
            m_err = np.sqrt(
                (ncol - 1) * np.nansum((window - m_mean) ** 2) / ncol
            )
            result.update(
                {
                    "start": s,
                    "end": end,
                    "chi2_dof": chi2_dof,
                    "n_diffs": n_diff,
                    "n_points": end - s,
                    "mass_est": float(m_mean),
                    "mass_err_est": float(m_err),
                    "ok": True,
                }
            )
            return result
    return result


# ==============================================================================
# b4.17 Nt 扫描任务扩展: cosh 对称点可配置 (对称点 = Ns/2, 而非常量 24)
# ------------------------------------------------------------------------------
# 空间关联函数 C(x) = A * cosh(m * (x - Ns/2)), 对 Ns = 32 / 36 / 40 / 48,
# 对称点分别是 16 / 18 / 20 / 24。旧接口写死 24, 只适用于 48^3 ensemble。
# 下面为并存的新实现 (旧接口保持原样, 既有结果可原样复现)。
# ==============================================================================

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
    与 solve_effective_mass_one_point 是同一个方程, 只是把对称点参数化。
    """
    if not np.isfinite(y):
        return np.nan

    if x < half and y <= 1.0:
        return 0.0 if y == 1.0 else np.nan

    a = float(x) - float(half)
    b = float(x) + 1.0 - float(half)

    def equation(m):
        # 注意: fsolve 会以 1 元素数组调用, 此处直接返回数组 (与旧接口写法一致)
        return cosh_ratio(m * a, m * b) - y

    try:
        sol = fsolve(equation, x0=0.1, xtol=1e-12, maxfev=2000)
        root = float(np.ravel(sol)[0])
    except Exception:
        return np.nan

    # 无实根的情形 (如 C(x+1) 与 C(x) 异号) fsolve 会返回伪解, 这里显式判为 NaN
    if not np.isfinite(root):
        return np.nan
    residual = abs(float(np.ravel(cosh_ratio(root * a, root * b))[0]) - y)
    if residual > 1e-8 * max(1.0, abs(y)):
        return np.nan
    return abs(root)


def _newton_column(y_col: np.ndarray, x: int, half: float,
                   max_iter: int = 200, tol: float = 1e-14) -> Tuple[np.ndarray, np.ndarray]:
    """
    对同一时间切片 x 上的整列 Jackknife 比值做向量化 Newton 迭代。
    返回 (解, 收敛掩码); 未收敛元素由调用方回退到标量 fsolve。
    """
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

    先用向量化 Newton 求解 (快), 未收敛/无解的元素回退到标量 fsolve
    (方程与仓库旧实现一致, 仅对称点参数化为 half), 兼顾速度与口径一致性。
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
