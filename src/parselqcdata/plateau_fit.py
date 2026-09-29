"""
src/parselqcdata/plateau_fit.py
--------------------------------------------------------------------------------
Python 整合层：贝叶斯平台双曲余弦拟合管道
- 拟合形式: f(t) = a * cosh(m * (t - 24))
- 执行 Jackknife 样本批量非线性拟合与参数汇总
--------------------------------------------------------------------------------
"""

import warnings

import numpy as np
import pandas as pd
import lsqfit
import gvar as gv
from typing import Dict, List, Tuple, Optional

# ==============================================================================
# 拟合器统一入口
# ------------------------------------------------------------------------------
# lsqfit.nonlinear_fit 的 fitter 默认值是 None, 含义是 "按环境自动挑":
# 装了 GSL 时用 'gsl_multifit', 没装时才退到 'scipy_least_squares'。
# 这会让同一份代码在不同机器上走出不同结果 (甚至退化), 因此本仓库**一律显式**
# 指定 chi2 最小二乘拟合器, 所有拟合必须经过 chi2_least_squares_fit()。
# ==============================================================================
LEAST_SQUARES_FITTER: str = "scipy_least_squares"


def chi2_least_squares_fit(data, fcn, p0, **kwargs):
    """
    统一的 chi2 最小二乘拟合入口 (强制 fitter='scipy_least_squares')。

    * 不接受调用方覆盖 fitter (传入也会被忽略), 避免个别调用点退化回环境默认值;
    * 其余关键字参数 (eps / tol / maxit / prior ...) 原样透传。
    """
    kwargs.pop("fitter", None)
    kwargs.setdefault("debug", False)
    return lsqfit.nonlinear_fit(
        data=data, fcn=fcn, p0=p0, fitter=LEAST_SQUARES_FITTER, **kwargs
    )


def target_cosh_func(x, p):
    return p['a'] * np.cosh(p['m'] * (x - 24))


def fit_single_jackknife_column(x_data: np.ndarray, y_val: np.ndarray, y_err: np.ndarray) -> dict:
    """对单个 Jackknife 样本列执行非线性拟合"""
    y_gv = gv.gvar(y_val, y_err)
    try:
        fit = chi2_least_squares_fit(
            data=(x_data, y_gv),
            fcn=target_cosh_func,
            p0={'a': 1.0, 'm': 0.5},
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


# ==============================================================================
# b4.17 Nt 扫描任务扩展: cosh 对称点可配置 (对称点 = Ns/2, 而非常量 24)
# ------------------------------------------------------------------------------
# 空间关联函数满足 C(x) = A * cosh(m * (x - Ns/2)), Ns = 32 / 36 / 40 / 48,
# 因此原 target_cosh_func (写死 t - 24) 不能直接用于这些 ensemble。
# 下面这些函数与上面的旧实现并存, 不修改旧接口以保证既有管道结果可复现。
# ==============================================================================

def centered_cosh_func(half: float):
    """生成对称点为 half 的 cosh 拟合函数: f(x) = a * cosh(m * (x - half))"""
    def _fcn(x, p):
        return p['a'] * np.cosh(p['m'] * (x - half))
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
    return {'a': a0, 'm': m0}


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


def fit_mass_window_scan(
    x_data: np.ndarray,
    y_mean: np.ndarray,
    y_err: np.ndarray,
    windows: List[Tuple[int, int]],
    half: float
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
            'x_start': int(start), 'x_end': int(end), 'n_points': int(end - start),
            'mass': np.nan, 'mass_err': np.nan, 'a': np.nan, 'a_err': np.nan,
            'chi2': np.nan, 'dof': np.nan, 'chi2_dof': np.nan, 'ok': False,
        }
        if end - start < 2:
            records.append(rec)
            continue
        fit = fit_single_jackknife_column_centered(
            x_data[start:end], y_mean[start:end], y_err[start:end], half
        )
        if fit:
            rec.update({
                'mass': float(fit['massfit_mean']),
                'mass_err': float(fit['massfit_err']),
                'a': float(fit['fita']),
                'a_err': float(fit['fita_err']),
                'chi2': float(fit['chi2']),
                'dof': float(fit['dof']),
                'chi2_dof': float(fit['chi2_dof']),
                'ok': True,
            })
        records.append(rec)
    return records


def fit_jackknife_mass_centered(
    x_data: np.ndarray,
    jk_matrix: np.ndarray,
    y_err: np.ndarray,
    start: int,
    end: int,
    half: float
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
                masses.append(float(fit['massfit_mean']))
                amps.append(float(fit['fita']))
                if np.isfinite(fit['chi2_dof']):
                    chi2_dofs.append(float(fit['chi2_dof']))
                per_sample.append({'sample': k, **fit})
            else:
                per_sample.append({'sample': k})

    masses_arr = np.asarray(masses, dtype=np.float64)
    amps_arr = np.asarray(amps, dtype=np.float64)

    result = {
        'mass': np.nan, 'mass_err': np.nan,
        'a': np.nan, 'a_err': np.nan,
        'chi2_dof': np.nan,
        'n_samples': int(n_samples), 'n_fitted': int(masses_arr.size),
        'x_start': int(start), 'x_end': int(end), 'half': float(half),
    }
    if masses_arr.size:
        mass_mean = float(np.mean(masses_arr))
        a_mean = float(np.mean(amps_arr)) if amps_arr.size else np.nan
        result['mass'] = mass_mean
        result['a'] = a_mean
        result['mass_err'] = float(np.sqrt((masses_arr.size - 1) * np.sum((masses_arr - mass_mean) ** 2) / masses_arr.size))
        if amps_arr.size:
            result['a_err'] = float(np.sqrt((amps_arr.size - 1) * np.sum((amps_arr - a_mean) ** 2) / amps_arr.size))
        if chi2_dofs:
            result['chi2_dof'] = float(np.mean(chi2_dofs))
    result['per_sample'] = per_sample
    return result


def select_best_window(
    scan_records: List[dict], min_points: int = 5, chi2_dof_max: float = 3.0
) -> Optional[dict]:
    """
    从窗口扫描结果中挑选最佳平台:
    取 chi2/dof 最小者 (chi2/dof 保留两位小数比较, 同分优先更长窗口),
    并要求 chi2/dof <= chi2_dof_max 且质量为有限正值。
    """
    candidates = [
        r for r in scan_records
        if r.get('ok') and r['n_points'] >= min_points
        and np.isfinite(r['mass']) and r['mass'] > 0.0
        and np.isfinite(r['chi2_dof']) and r['chi2_dof'] <= chi2_dof_max
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda r: (round(float(r['chi2_dof']), 2), -int(r['n_points']), int(r['x_start'])))
