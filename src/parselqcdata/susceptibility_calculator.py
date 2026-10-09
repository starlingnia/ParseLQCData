#!/usr/bin/env python3
"""
src/parselqcdata/susceptibility_calculator.py
--------------------------------------------------------------------------------
独立脚本与模块：手征磁化率 Jackknife 统计推断与标度折算计算器
- 算法理论:
    \\chi_{\\text{disc}} = \\langle (\\bar{\\psi}\\psi)^2 \\rangle - \\langle \\bar{\\psi}\\psi \\rangle^2
    a_r = JK(Obar)_r, b_r = JK(O2bar)_r
    \\chi_r = b_r - a_r^2
    mean(\\chi) = (1/N) * \\sum \\chi_r,  err(\\chi) = \\sqrt{N-1} * \\text{std}(\\chi_r, \\text{ddof}=0)
- 执行物理体积标度 (F_vol) 与质量重整化 (除以 Zm^2)
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Optional, Union
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.parselqcdata.scaling_factors import compute_scaling_factors, get_zm_factor


def jackknife_resample(arr: np.ndarray) -> np.ndarray:
    """
    对一维样本序列执行 Jackknife Leave-One-Out 重采样。
    输入大小 N，输出大小 N 的 Jackknife 假样本序列。
    """
    n = len(arr)
    if n <= 1:
        raise ValueError("Jackknife 重采样至少需要 2 个构型样本")
    total_sum = np.sum(arr)
    return (total_sum - arr) / (n - 1)


def compute_jackknife_susceptibility(
    list_obar: List[float],
    list_o2bar: List[float],
    ns: int = 48,
    nt: int = 16,
    temp_mev: float = 157.0,
    beta: Optional[Union[float, str]] = None,
    zm: Optional[float] = None,
) -> Dict[str, float]:
    """
    通过 Jackknife 假样本计算手征磁化率 \\chi = \\langle O^2 \\rangle - \\langle O \\rangle^2
    并根据 Ns, Nt 与温度 T 准确折算标度因子与质量重整化因子 Zm：
    - factor_vol = Ns^3 * Nt
    - factor_scaled = Ns^3 * Nt^3 * T^2
    - chi_disc(lattice unit) = Mean_vol_scaled
    - chi_disc(GeV^2 renormalized) = (Mean_scaled / 1e6) / Zm**2

    Parameters:
        list_obar: 各构型上的算符单体线性均值列表
        list_o2bar: 各构型上的算符两体无偏乘积均值列表
        ns: 空间格点点数
        nt: 时间格点点数
        temp_mev: 物理温度 (MeV)
        beta: 规范耦合常数
        zm: 夸克质量重整化常数 (若未提供则根据 beta 自动查表)

    Returns:
        Dict 包含未标度、格点体积标度、物理温度标度及重整化后的均值与 Jackknife 误差
    """
    n = len(list_obar)
    if n <= 1 or len(list_o2bar) != n:
        raise ValueError(f"样本维度不合法: obar={len(list_obar)}, o2bar={len(list_o2bar)}")

    arr_obar = np.array(list_obar, dtype=np.float64)
    arr_o2bar = np.array(list_o2bar, dtype=np.float64)

    a_jk = jackknife_resample(arr_obar)
    b_jk = jackknife_resample(arr_o2bar)

    # 每一个 Jackknife bin 中的磁化率估计
    chi_jk = b_jk - (a_jk**2)

    mean_unscaled = float(np.mean(chi_jk))
    # 统计标准误: \\sigma = \\sqrt{N-1} * \\text{std}(chi_jk, \\text{ddof}=0)
    err_unscaled = float(np.sqrt(n - 1) * np.std(chi_jk, ddof=0))

    f_vol, f_scaled = compute_scaling_factors(ns, nt, temp_mev)

    mean_vol_scaled = mean_unscaled * f_vol
    error_vol_scaled = err_unscaled * f_vol

    mean_scaled = mean_unscaled * f_scaled
    error_scaled = err_unscaled * f_scaled

    # 物理量纲 [GeV^2] (未重整化)
    mean_scaled_gev2 = mean_scaled / 1e6
    error_scaled_gev2 = error_scaled / 1e6

    # 质量重整化因子 Zm (最新标准口径除以 Zm**2)
    eff_zm = zm if zm is not None else (get_zm_factor(beta) if beta is not None else 1.0)
    inv_zm = (1.0 / eff_zm**2) if abs(eff_zm) > 1e-15 else 1.0

    mean_scaled_renorm = mean_scaled * inv_zm
    error_scaled_renorm = error_scaled * inv_zm
    mean_scaled_gev2_renorm = mean_scaled_gev2 * inv_zm
    error_scaled_gev2_renorm = error_scaled_gev2 * inv_zm

    # 完整统计样本数据序列
    jk_unscaled = chi_jk
    jk_vol = chi_jk * f_vol
    jk_scaled = chi_jk * f_scaled
    jk_gev2 = jk_scaled / 1e6
    jk_renorm = jk_gev2 * inv_zm

    return {
        "mean_unscaled": mean_unscaled,
        "error_unscaled": err_unscaled,
        "factor_vol": f_vol,
        "mean_vol_scaled": mean_vol_scaled,
        "error_vol_scaled": error_vol_scaled,
        "factor_scaled": f_scaled,
        "mean_scaled": mean_scaled,
        "error_scaled": error_scaled,
        "mean_scaled_gev2": mean_scaled_gev2,
        "error_scaled_gev2": error_scaled_gev2,
        "zm": eff_zm,
        "mean_scaled_renorm": mean_scaled_renorm,
        "error_scaled_renorm": error_scaled_renorm,
        "mean_scaled_gev2_renorm": mean_scaled_gev2_renorm,
        "error_scaled_gev2_renorm": error_scaled_gev2_renorm,
        "num_configs": n,
        "ns": ns,
        "nt": nt,
        "temp": temp_mev,
        # 完整统计重采样样本
        "jk_samples_unscaled": jk_unscaled,
        "jk_samples_vol": jk_vol,
        "jk_samples_scaled": jk_scaled,
        "jk_samples_gev2": jk_gev2,
        "jk_samples_renorm": jk_renorm,
    }


def save_susceptibility_jk_samples_csv(
    res: Dict[str, Union[float, int, np.ndarray, list]],
    csv_path: Path,
    ensemble_name: str = "",
) -> None:
    """将单个系综的完整 Jackknife 统计重采样样本序列持久化保存为独立 CSV 文件"""
    import polars as pl

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    jk_renorm = np.asarray(res.get("jk_samples_renorm", []))
    n = len(jk_renorm)
    if n == 0:
        return

    jk_unscaled = np.asarray(res.get("jk_samples_unscaled", np.zeros(n)))
    jk_vol = np.asarray(res.get("jk_samples_vol", np.zeros(n)))
    jk_gev2 = np.asarray(res.get("jk_samples_gev2", np.zeros(n)))

    df_jk = pl.DataFrame({
        "jk_index": list(range(n)),
        "chi_unscaled": jk_unscaled,
        "chi_vol_scaled": jk_vol,
        "chi_scaled_gev2": jk_gev2,
        "chi_renorm": jk_renorm,
    })
    if ensemble_name:
        df_jk = df_jk.with_columns(pl.lit(ensemble_name).alias("ensemble"))
        df_jk = df_jk.select(["ensemble", "jk_index", "chi_unscaled", "chi_vol_scaled", "chi_scaled_gev2", "chi_renorm"])

    df_jk.write_csv(csv_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="手征磁化率 Jackknife 统计计算器")
    parser.add_argument("--ns", type=int, default=48, help="空间格点点数")
    parser.add_argument("--nt", type=int, default=16, help="时间格点点数")
    parser.add_argument("--temp", type=float, default=153.31, help="温度 (MeV)")
    parser.add_argument("--beta", type=str, default="4.17", help="耦合常数 beta")
    parser.add_argument("--zm", type=float, default=None, help="重整化因子 Zm (默认自动查表)")
    args = parser.parse_args()

    # 简单演示或快速测试运行
    print(f"--- 手征磁化率计算器 (Ns={args.ns}, Nt={args.nt}, T={args.temp} MeV, beta={args.beta}) ---")
    dummy_obar = [1.0001, 1.0002, 1.0003, 1.0000, 1.0004]
    dummy_o2bar = [1.000201, 1.000402, 1.000603, 1.000001, 1.000804]
    res = compute_jackknife_susceptibility(
        dummy_obar, dummy_o2bar, ns=args.ns, nt=args.nt, temp_mev=args.temp, beta=args.beta, zm=args.zm
    )
    for k, v in res.items():
        if isinstance(v, float):
            print(f"  {k:<24}: {v:+.8e}")
        else:
            print(f"  {k:<24}: {v}")


if __name__ == "__main__":
    main()
