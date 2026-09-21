"""
src/parselqcdata/plateau_fit.py
--------------------------------------------------------------------------------
Python 整合层：贝叶斯平台双曲余弦拟合管道
- 拟合形式: f(t) = a * cosh(m * (t - 24))
- 执行 Jackknife 样本批量非线性拟合与参数汇总
--------------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
import lsqfit
import gvar as gv
from typing import Dict, List, Tuple, Optional


def target_cosh_func(x, p):
    return p['a'] * np.cosh(p['m'] * (x - 24))


def fit_single_jackknife_column(x_data: np.ndarray, y_val: np.ndarray, y_err: np.ndarray) -> dict:
    """对单个 Jackknife 样本列执行非线性拟合"""
    y_gv = gv.gvar(y_val, y_err)
    try:
        fit = lsqfit.nonlinear_fit(
            data=(x_data, y_gv),
            fcn=target_cosh_func,
            p0={'a': 1.0, 'm': 0.5},
            fitter='scipy_least_squares',
            debug=False
        )
        return {
            'chi2': float(fit.chi2),
            'dof': fit.dof,
            'massfit_mean': abs(float(fit.p['m'].mean)),
            'massfit_err': float(fit.p['m'].sdev),
            'fita': float(fit.p['a'].mean),
            'fita_err': float(fit.p['a'].sdev),
            'chi2_dof': float(fit.chi2 / (fit.dof - 1)) if fit.dof > 1 else float('inf'),
        }
    except Exception:
        return {}


def fit_meson_plateau(
    df_sym: pd.DataFrame,
    err_vals: np.ndarray,
    slice_start: int,
    slice_end: int = 25
) -> Tuple[pd.DataFrame, Optional[dict]]:
    """
    对对称折叠矩阵全部 Jackknife 样本进行切片平台拟合
    返回: (df_results, summary_dict)
    """
    x_array = np.arange(48)[slice_start:slice_end]
    fits = []

    for col in df_sym.columns:
        y_val = df_sym[col].to_numpy()[slice_start:slice_end]
        y_err = err_vals[slice_start:slice_end]
        f = fit_single_jackknife_column(x_array, y_val, y_err)
        if f and f.get('chi2_dof', float('inf')) <= 1000000.0:
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
            "mass_mean": float(means['massfit_mean']),
            "mass_err": float(jack_err['massfit_mean']),
            "a_mean": float(means['fita']),
            "a_err": float(jack_err['fita']),
            "chi2_dof": float(means['chi2_dof']),
            "N_samples": int(n_samples)
        }

    return df_res, summary
