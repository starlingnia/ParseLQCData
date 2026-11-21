#!/usr/bin/env python3
"""
scripts/tasks/lcp/run_pade_fits.py
--------------------------------------------------------------------------------
常物理线 (LCP) 多系综 Padé 有理函数拟合脚本:
1. 指定输入输出路径
2. 调用 src/parselqcdata/pade_fitter.py 核心拟合函数
3. 导出拟合汇总表与 Jackknife 样本
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_LCP_DIR
from src.parselqcdata.pade_fitter import fit_jackknife_susceptibility_pade

BASE_DIR = PROJECT_ROOT / "22pade"

# 各系综选用的有效拟合阶数
ENSEMBLE_ORDERS = {
    "p2": ["1/2", "2/2", "3/2"],
    "p3": ["0/2", "1/2", "2/2"],
    "q2": ["1/2", "2/2", "3/2", "1/3"],
    "q3": ["0/2", "1/2"],
}


def fit_ensemble(ds_name: str, orders: list[str], t0: float = 157.0) -> None:
    ds_dir = BASE_DIR / ds_name
    susc_df = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
    temps = susc_df["temp"].to_numpy().astype(float)
    errors = susc_df["chi_renorm_err"].to_numpy().astype(float)
    betas = susc_df["beta"].to_numpy().astype(float)

    # 读取 Jackknife 样本构建矩阵 (N_temps, N_samples)
    if ds_name == "q3":
        jk_csv = OUTPUT_LCP_DIR / "susceptibility_jk_samples_L64.csv"
    else:
        jk_csv = OUTPUT_LCP_DIR / ds_name / f"susceptibility_jk_samples_{ds_name}.csv"
    if not jk_csv.exists():
        jk_csv = ds_dir / f"pade22_jk_samples_{ds_name}.csv"

    jk_df = pl.read_csv(jk_csv)
    samples = []
    for b in betas:
        s = jk_df.filter(pl.col("beta") == b).sort("jk_index")["chi_renorm"].to_numpy().astype(float)
        samples.append(s)
    min_len = min(len(s) for s in samples)
    chi_matrix = np.array([s[:min_len] for s in samples], dtype=float)

    # 循环调用脚手架拟合函数
    for ord_val in orders:
        ord_tag = ord_val.replace("/", "")
        summary, records = fit_jackknife_susceptibility_pade(
            temps=temps,
            chi_matrix=chi_matrix,
            errors=errors,
            order=ord_val,
            t0=t0,
        )
        # 写输出地址
        pl.DataFrame([summary]).write_csv(ds_dir / f"pade{ord_tag}_fit_summary_{ds_name}.csv")
        pl.DataFrame(records).write_csv(ds_dir / f"pade{ord_tag}_jk_samples_{ds_name}.csv")
        print(f"  [OK] {ds_name} [{ord_val}]: Tc = {summary['t_peak_mean']:.2f} +/- {summary['t_peak_err']:.2f} MeV, chi2/dof = {summary['chi2_dof_mean']:.3f}")


def main() -> None:
    print("==========================================================================")
    print("[LCP-PADE-FIT] 启动各系综 Padé 有理函数拟合...")
    print("==========================================================================")
    for ds_name, orders in ENSEMBLE_ORDERS.items():
        print(f"\n>>> 拟合系综 {ds_name}...")
        fit_ensemble(ds_name, orders)
    print("\n[DONE] 所有系综拟合完成！")


if __name__ == "__main__":
    main()
