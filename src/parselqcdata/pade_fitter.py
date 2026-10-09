#!/usr/bin/env python3
"""
src/parselqcdata/pade_fitter.py
--------------------------------------------------------------------------------
手征磁化率 Padé 有理式逼近拟合引擎 (基于 lsqfit + gvar)
--------------------------------------------------------------------------------
【物理与算法背景】:
1. 在有限温度格点 QCD 中，手征磁化率 \\chi(T) 在赝临界相变区呈现单峰结构。
2. 采用 Padé 有理式逼近 (Rational Function Approximation) 对 \\chi(T) 建模:
     R(T) = P_m(T) / Q_n(T)
   其中展开中心通常取相变区中心 T_0 (如 157.0 MeV)，有效避免高次多项式在绝对温度下的数值病态。
3. 极点 (Pole) 与峰值 (Peak):
   - 分母 Q_n(T) = 0 的复根给出复温度平面上的奇点 (Lee-Yang / Fisher 极点):
     T_{pole} = Re(T_{pole}) \\pm i Im(T_{pole})
   - 实数轴上的极值点给出赝临界峰值温度 T_{peak} 和峰值磁化率 \\chi_{peak}。
4. Jackknife 统计推断:
   - 每个温度截取相同样本量 N (如 2000)。
   - 对每个 Jackknife 重采样统计切片 k，使用 lsqfit 执行非线性拟合。
   - 提取每个切片的 T_{peak}^{(k)}、Re(T_{pole})^{(k)}、Im(T_{pole})^{(k)} 与 \\chi_{peak}^{(k)}。
   - 最终按 Jackknife 公式计算均值与统计标准误:
       mean = (1/N) * \\sum X^{(k)}
       err  = \\sqrt{N - 1} * std(X^{(k)}, ddof=0)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence, Tuple, Union
import warnings

import numpy as np
import polars as pl
import lsqfit
import gvar as gv
from scipy.optimize import minimize_scalar

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_LCP_DIR, OUTPUT_CONDENSATE_DIR
from src.parselqcdata.plateau_fit import chi2_least_squares_fit


class PadeModel:
    """Padé 有理式逼近模型基类与常用阶数实现"""

    def __init__(self, order: str = "3/2", t0: float = 157.0) -> None:
        self.order = order.strip()
        self.t0 = float(t0)

        if self.order == "0/2":
            self.param_keys = ["p0", "q1", "q2"]
        elif self.order == "1/2":
            self.param_keys = ["p0", "p1", "q1", "q2"]
        elif self.order == "2/2":
            self.param_keys = ["p0", "p1", "p2", "q1", "q2"]
        elif self.order == "3/2":
            self.param_keys = ["p0", "p1", "p2", "p3", "q1", "q2"]
        elif self.order == "0/4":
            self.param_keys = ["p0", "q1", "q2", "q3", "q4"]
        else:
            raise ValueError(f"暂不支持的 Padé 阶数: {order} (支持: '0/2', '1/2', '2/2', '3/2', '0/4')")

    def eval(self, t: Union[float, np.ndarray], p: Dict[str, Union[float, gv.GVar]]) -> Union[float, np.ndarray]:
        """计算 Padé 有理函数值 R(t; p)"""
        dx = t - self.t0
        if self.order == "0/2":
            num = p["p0"]
            den = 1.0 + p["q1"] * dx + p["q2"] * (dx**2)
        elif self.order == "1/2":
            num = p["p0"] + p["p1"] * dx
            den = 1.0 + p["q1"] * dx + p["q2"] * (dx**2)
        elif self.order == "2/2":
            num = p["p0"] + p["p1"] * dx + p["p2"] * (dx**2)
            den = 1.0 + p["q1"] * dx + p["q2"] * (dx**2)
        elif self.order == "3/2":
            num = p["p0"] + p["p1"] * dx + p["p2"] * (dx**2) + p["p3"] * (dx**3)
            den = 1.0 + p["q1"] * dx + p["q2"] * (dx**2)
        elif self.order == "0/4":
            num = p["p0"]
            den = 1.0 + p["q1"] * dx + p["q2"] * (dx**2) + p["q3"] * (dx**3) + p["q4"] * (dx**4)
        else:
            raise ValueError(f"未知阶数 {self.order}")
        return num / den

    def initial_guess(self, t_arr: np.ndarray, y_arr: np.ndarray) -> Dict[str, float]:
        """根据数据特征自动构建稳定的初值"""
        peak_idx = int(np.nanargmax(y_arr))
        p0_est = float(y_arr[peak_idx])
        t_peak_est = float(t_arr[peak_idx])
        dx_peak = t_peak_est - self.t0

        # 分母二阶展开项估算: 宽度通常约 15 ~ 25 MeV
        width_est = 20.0
        q2_est = 1.0 / (width_est**2)
        q1_est = -2.0 * q2_est * dx_peak

        if self.order == "0/2":
            return {"p0": max(p0_est, 0.1), "q1": q1_est, "q2": q2_est}
        elif self.order == "1/2":
            return {"p0": max(p0_est, 0.1), "p1": 0.0, "q1": q1_est, "q2": q2_est}
        elif self.order == "2/2":
            return {"p0": max(p0_est, 0.1), "p1": 0.0, "p2": 0.0, "q1": q1_est, "q2": q2_est}
        elif self.order == "3/2":
            return {"p0": max(p0_est, 0.1), "p1": 0.0, "p2": 0.0, "p3": 0.0, "q1": q1_est, "q2": q2_est}
        elif self.order == "0/4":
            return {"p0": max(p0_est, 0.1), "q1": q1_est, "q2": q2_est, "q3": 0.0, "q4": q2_est**2}
        return {"p0": 1.0}

    def solve_pole(self, p: Dict[str, float]) -> Tuple[float, float]:
        """
        求解分母 Q_n(T) = 0 的最近复根 (Pole):
        Returns:
            (Re(T_pole), Im(T_pole))
        """
        if self.order in ["0/2", "1/2", "2/2", "3/2"]:
            q1 = float(p.get("q1", 0.0))
            q2 = float(p.get("q2", 1e-4))
            if abs(q2) < 1e-14:
                # 线性分母
                if abs(q1) > 1e-14:
                    return self.t0 - 1.0 / q1, 0.0
                return np.nan, np.nan
            disc = q1**2 - 4.0 * q2
            re_part = self.t0 - q1 / (2.0 * q2)
            if disc < 0:
                im_part = float(np.sqrt(-disc) / (2.0 * abs(q2)))
            else:
                im_part = 0.0
            return re_part, im_part

        elif self.order == "0/4":
            # 四次多项式求根
            q1 = float(p.get("q1", 0.0))
            q2 = float(p.get("q2", 0.0))
            q3 = float(p.get("q3", 0.0))
            q4 = float(p.get("q4", 0.0))
            coeffs = [q4, q3, q2, q1, 1.0]
            roots = np.roots(coeffs)
            # 找到虚部绝对值最小且大于 0 的根，或离实轴最近的根
            nearest = min(roots, key=lambda r: abs(r.imag))
            return float(self.t0 + nearest.real), float(abs(nearest.imag))

        return np.nan, np.nan

    def solve_peak(
        self,
        p: Dict[str, float],
        t_bounds: Tuple[float, float] = (140.0, 175.0),
    ) -> Tuple[float, float]:
        """
        求解实数轴上的最大值峰值点 (T_peak, chi_peak)
        """
        if self.order == "0/2":
            # 对于 [0/2], 分母二次多项式最小值点即为最大值点
            q1 = float(p.get("q1", 0.0))
            q2 = float(p.get("q2", 1e-4))
            if q2 > 0 and (4.0 * q2 > q1**2):
                t_pk = self.t0 - q1 / (2.0 * q2)
                chi_pk = float(self.eval(t_pk, p))
                return t_pk, chi_pk

        # 通用优化搜索峰值
        def neg_target(x_val: float) -> float:
            try:
                v = float(self.eval(x_val, p))
                return -v if np.isfinite(v) else 1e9
            except Exception:
                return 1e9

        res = minimize_scalar(neg_target, bounds=t_bounds, method="bounded")
        if res.success:
            return float(res.x), float(-res.fun)

        # 退化保护
        grid_t = np.linspace(t_bounds[0], t_bounds[1], 200)
        grid_y = np.array([float(self.eval(x_val, p)) for x_val in grid_t])
        idx = int(np.nanargmax(grid_y))
        return float(grid_t[idx]), float(grid_y[idx])


def fit_single_sample_pade(
    t_data: np.ndarray,
    y_data: np.ndarray,
    err_data: np.ndarray,
    model: PadeModel,
    p0: Optional[Dict[str, float]] = None,
) -> Optional[dict]:
    """对单组数据样本执行 lsqfit Padé 有理函数非线性拟合"""
    t_data = np.asarray(t_data, dtype=np.float64)
    y_data = np.asarray(y_data, dtype=np.float64)
    err_data = np.asarray(err_data, dtype=np.float64)

    if p0 is None:
        p0 = model.initial_guess(t_data, y_data)

    def fcn(x, p):
        return model.eval(x, p)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = chi2_least_squares_fit(
                data=(t_data, gv.gvar(y_data, err_data)),
                fcn=fcn,
                p0=p0,
            )

        p_fit = {k: float(v.mean) for k, v in fit.p.items()}
        chi2 = float(fit.chi2)
        dof = float(fit.dof)
        chi2_dof = chi2 / dof if dof > 0 else np.nan

        # 求解极点与峰值
        t_pole_re, t_pole_im = model.solve_pole(p_fit)
        t_peak, chi_peak = model.solve_peak(p_fit, t_bounds=(float(np.min(t_data)), float(np.max(t_data))))

        return {
            "ok": True,
            "params": p_fit,
            "chi2": chi2,
            "dof": dof,
            "chi2_dof": chi2_dof,
            "t_peak": t_peak,
            "chi_peak": chi_peak,
            "t_pole_real": t_pole_re,
            "t_pole_imag": t_pole_im,
        }
    except Exception as e:
        return {
            "ok": False,
            "error": str(e),
            "params": p0,
            "chi2": np.nan,
            "dof": np.nan,
            "chi2_dof": np.nan,
            "t_peak": np.nan,
            "chi_peak": np.nan,
            "t_pole_real": np.nan,
            "t_pole_imag": np.nan,
        }


def load_and_align_susceptibility_samples(
    csv_path: Optional[Path] = None,
    max_samples: Optional[int] = None,
    binsize: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[float]]:
    """
    载入并对齐各温度的 Jackknife 样本:
    - 确保每个温度具有完全一致的统计样本量 N = min(N_i) 或 max_samples
    - 支持 binsize 降采样 / 重抽样 (先不用做 binsize 处理时 binsize=1)

    Returns:
        (temps, chi_matrix, errors, betas):
        - temps: 一维温度数组 (K,)
        - chi_matrix: 二维样本矩阵 (K, N)
        - errors: 各温度标准误差 (K,)
        - betas: 对应的 beta 列表
    """
    if csv_path is None or not csv_path.exists():
        cand1 = OUTPUT_LCP_DIR / "susceptibility_jk_samples.csv"
        cand2 = OUTPUT_CONDENSATE_DIR / "results_susceptibility.csv"
        csv_path = cand1 if cand1.exists() else cand2

    if not csv_path.exists():
        raise FileNotFoundError(f"找不到磁化率样本数据文件: {csv_path}")

    df = pl.read_csv(csv_path)

    # 识别列名
    beta_col = "beta" if "beta" in df.columns else "Beta"
    temp_col = "temp" if "temp" in df.columns else "Temp"
    chi_col = "chi_renorm" if "chi_renorm" in df.columns else "Mean_scaled_gev2_renorm"

    # 按 beta / temp 分组提取
    betas_unique = sorted(df[beta_col].unique().to_list())
    temp_list = []
    samples_per_temp = []
    beta_ordered = []

    for b in betas_unique:
        sub = df.filter(pl.col(beta_col) == b)
        if "jk_index" in sub.columns:
            sub = sub.sort("jk_index")
        t_val = float(sub[temp_col][0])
        s_vals = sub[chi_col].to_numpy().astype(np.float64)
        if len(s_vals) > 0:
            temp_list.append(t_val)
            samples_per_temp.append(s_vals)
            beta_ordered.append(b)

    if not samples_per_temp:
        raise ValueError("未能提取到任何有效的磁化率温度样本序列")

    # 1. 确定统一的样本量 N
    min_len = min(len(s) for s in samples_per_temp)
    target_n = min_len if max_samples is None else min(min_len, int(max_samples))

    # 截取前 target_n 个样本以保持各温度样本量严格对齐
    aligned_matrix = np.array([s[:target_n] for s in samples_per_temp], dtype=np.float64)

    # 2. 如果用户指定了 binsize > 1，执行 binning 聚合与重新 Jackknife
    if binsize > 1:
        n_bins = target_n // binsize
        if n_bins < 2:
            raise ValueError(f"binsize={binsize} 过大，导致剩余 bin 数量小于 2 (当前 N={target_n})")
        binned_matrix = np.zeros((len(temp_list), n_bins), dtype=np.float64)
        for i in range(len(temp_list)):
            raw_s = aligned_matrix[i, : n_bins * binsize]
            # 计算各 bin 的局部均值
            bin_means = raw_s.reshape((n_bins, binsize)).mean(axis=1)
            # 对 bin_means 做 Jackknife leave-one-out
            sum_b = np.sum(bin_means)
            binned_matrix[i, :] = (sum_b - bin_means) / (n_bins - 1)
        aligned_matrix = binned_matrix
        target_n = n_bins

    # 3. 计算每个温度的标准误差
    errors = np.zeros(len(temp_list), dtype=np.float64)
    for i in range(len(temp_list)):
        s_row = aligned_matrix[i, :]
        errors[i] = float(np.sqrt(target_n - 1) * np.std(s_row, ddof=0))
        if errors[i] <= 0:
            errors[i] = 1e-4

    return np.array(temp_list), aligned_matrix, errors, beta_ordered


def fit_jackknife_susceptibility_pade(
    temps: np.ndarray,
    chi_matrix: np.ndarray,
    errors: np.ndarray,
    order: str = "3/2",
    t0: float = 157.0,
) -> Tuple[dict, List[dict]]:
    """
    对每个 Jackknife 统计重采样样本执行 Padé 有理函数拟合，并统计极点与峰值的 Jackknife 均值与方差。

    Parameters:
        temps: 一维温度数组 (K,)
        chi_matrix: (K, N) 二维磁化率 Jackknife 样本矩阵
        errors: (K,) 各温度在格点系综上的统计标准误差 (对角不确定度权重)
        order: Padé 逼近阶数 (默认: "3/2")
        t0: 展开参考中心温度 (默认: 157.0 MeV)

    Returns:
        (summary_dict, fit_records):
        - summary_dict: 汇总统计量 (含 mean 与 Jackknife 误差)
        - fit_records: 每个 Jackknife 切片的逐样本拟合明细列表
    """
    temps = np.asarray(temps, dtype=np.float64)
    chi_matrix = np.asarray(chi_matrix, dtype=np.float64)
    errors = np.asarray(errors, dtype=np.float64)

    n_temps, n_samples = chi_matrix.shape
    model = PadeModel(order=order, t0=t0)

    # 1. 首先拟合系综均值，获得全局优良基底参数作为初值
    mean_y = np.mean(chi_matrix, axis=1)
    fit_mean = fit_single_sample_pade(temps, mean_y, errors, model)
    if not fit_mean or not fit_mean["ok"]:
        p_warm = model.initial_guess(temps, mean_y)
    else:
        p_warm = dict(fit_mean["params"])

    # 2. 顺序对每个 Jackknife 样本执行非线性拟合 (热启动加速)
    fit_records: List[dict] = []
    t_peaks: List[float] = []
    chi_peaks: List[float] = []
    t_poles_re: List[float] = []
    t_poles_im: List[float] = []
    chi2_list: List[float] = []

    p_curr = dict(p_warm)

    for k in range(n_samples):
        y_k = chi_matrix[:, k]
        res_k = fit_single_sample_pade(temps, y_k, errors, model, p0=p_curr)

        rec = {
            "jk_index": k,
            "order": order,
            "ok": res_k["ok"],
            "chi2": res_k["chi2"],
            "dof": res_k["dof"],
            "chi2_dof": res_k["chi2_dof"],
            "t_peak": res_k["t_peak"],
            "chi_peak": res_k["chi_peak"],
            "t_pole_real": res_k["t_pole_real"],
            "t_pole_imag": res_k["t_pole_imag"],
        }
        for pk, pv in res_k.get("params", {}).items():
            rec[pk] = pv

        fit_records.append(rec)

        if res_k["ok"] and np.isfinite(res_k["t_peak"]):
            t_peaks.append(res_k["t_peak"])
            chi_peaks.append(res_k["chi_peak"])
            t_poles_re.append(res_k["t_pole_real"])
            t_poles_im.append(res_k["t_pole_imag"])
            chi2_list.append(res_k["chi2_dof"])
            p_curr = dict(res_k["params"])

    # 3. 统计汇总计算
    def compute_jk_stats(arr: Sequence[float]) -> Tuple[float, float]:
        v_arr = np.asarray(arr, dtype=np.float64)
        mask = np.isfinite(v_arr)
        if not np.any(mask):
            return np.nan, np.nan
        valid = v_arr[mask]
        m = float(np.mean(valid))
        n_valid = len(valid)
        if n_valid <= 1:
            return m, np.nan
        err = float(np.sqrt(n_valid - 1) * np.std(valid, ddof=0))
        return m, err

    tpk_mean, tpk_err = compute_jk_stats(t_peaks)
    cpk_mean, cpk_err = compute_jk_stats(chi_peaks)
    pre_mean, pre_err = compute_jk_stats(t_poles_re)
    pim_mean, pim_err = compute_jk_stats(t_poles_im)
    chi2_mean = float(np.nanmean(chi2_list)) if chi2_list else np.nan

    summary_dict = {
        "order": order,
        "t0": t0,
        "n_samples": n_samples,
        "n_fitted": len(t_peaks),
        "chi2_dof_mean": chi2_mean,
        "t_peak_mean": tpk_mean,
        "t_peak_err": tpk_err,
        "chi_peak_mean": cpk_mean,
        "chi_peak_err": cpk_err,
        "t_pole_real_mean": pre_mean,
        "t_pole_real_err": pre_err,
        "t_pole_imag_mean": pim_mean,
        "t_pole_imag_err": pim_err,
    }

    return summary_dict, fit_records


def main() -> None:
    parser = argparse.ArgumentParser(description="手征磁化率 Padé 有理式拟合分析工具")
    parser.add_argument("--order", type=str, default="3/2", choices=["0/2", "1/2", "2/2", "3/2", "0/4"], help="Padé 拟合阶数 (默认: 3/2)")
    parser.add_argument("--t0", type=float, default=157.0, help="展开中心温度 T0 (MeV, 默认: 157.0)")
    parser.add_argument("--binsize", type=int, default=1, help="构型重抽样 binsize (默认: 1，不执行 bin 聚合)")
    parser.add_argument("--max-samples", type=int, default=2000, help="截取对齐的样本总数 (默认: 2000)")
    parser.add_argument("--csv", type=Path, default=None, help="输入的磁化率样本 CSV 路径 (默认自动推断)")
    parser.add_argument("--out-dir", type=Path, default=OUTPUT_LCP_DIR, help="拟合结果导出目录 (默认: output/LCP)")
    args = parser.parse_args()

    print(f"\n==========================================================================")
    print(f"[PADE-FITTER] 启动手征磁化率 Padé [{args.order}] 有理函数拟合...")
    print(f"  - 展开中心 T0: {args.t0} MeV")
    print(f"  - 样本量对齐: {args.max_samples}")
    print(f"  - Binsize 设定: {args.binsize}")
    print(f"==========================================================================")

    temps, chi_matrix, errors, betas = load_and_align_susceptibility_samples(
        csv_path=args.csv,
        max_samples=args.max_samples,
        binsize=args.binsize,
    )

    print(f"[INFO] 成功载入 {len(temps)} 个温度点，样本矩阵维度: {chi_matrix.shape}")
    for t_val, b_val, err_val in zip(temps, betas, errors):
        print(f"  -> Beta {b_val:<6}: T = {t_val:6.2f} MeV, 统计标准误 = {err_val:.6e}")

    summary, records = fit_jackknife_susceptibility_pade(
        temps=temps,
        chi_matrix=chi_matrix,
        errors=errors,
        order=args.order,
        t0=args.t0,
    )

    print("\n--------------------------------------------------------------------------")
    print("【拟合与统计结果 (包含 Jackknife 偏差)】:")
    print(f"  - 拟合通过率:               {summary['n_fitted']} / {summary['n_samples']} ({summary['n_fitted']/summary['n_samples']*100:.1f}%)")
    print(f"  - 平均 \\chi^2/dof:           {summary['chi2_dof_mean']:.3f}")
    print(f"  - 峰值温度 T_peak:           {summary['t_peak_mean']:.4f} +/- {summary['t_peak_err']:.4f} MeV")
    print(f"  - 极点实部 Re(T_pole):       {summary['t_pole_real_mean']:.4f} +/- {summary['t_pole_real_err']:.4f} MeV")
    print(f"  - 极点虚部 Im(T_pole):       {summary['t_pole_imag_mean']:.4f} +/- {summary['t_pole_imag_err']:.4f} MeV")
    print(f"  - 峰值磁化率 \\chi_peak:       {summary['chi_peak_mean']:.4f} +/- {summary['chi_peak_err']:.4f} GeV^2")
    print("--------------------------------------------------------------------------")

    # 持久化输出
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_df = pl.DataFrame([summary])
    summary_path = out_dir / "pade_fit_summary.csv"
    summary_df.write_csv(summary_path)

    samples_df = pl.DataFrame(records)
    samples_path = out_dir / "pade_jk_samples.csv"
    samples_df.write_csv(samples_path)

    print(f"[OK] 拟合结果已成功保存至:\n  - 汇总报告: {summary_path}\n  - 样本明细: {samples_path}\n")


if __name__ == "__main__":
    main()
