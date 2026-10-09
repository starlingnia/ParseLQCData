#!/usr/bin/env python3
"""
scripts/tasks/lcp/04_fit_pade_susceptibility.py
--------------------------------------------------------------------------------
LCP (常物理线) 子任务 4: 手征磁化率 Padé 有理函数逼近与极点/峰值 Jackknife 推断
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 目的:
   - 拟合手征磁化率随温度变化 \\chi(T)，提取赝临界相变温度 T_{pc} 以及复平面上的有理函数极点 (Pole)。
   - 得到峰值温度与极点位置的完整 Jackknife 统计偏差。
2. 数据准备:
   - 从 output/LCP/susceptibility_jk_samples.csv 载入各温度的重整化磁化率重采样样本。
   - 截取使每个温度具有完全相同的统计样本数 (例如 N=2000)。
   - 不考虑信道/温度间的协方差关联 (独立对角误差权重)。
   - 默认不进行额外 binsize 聚合 (binsize=1)，亦支持 --binsize B 灵活指定。
3. Padé 有理式逼近 (基于 lsqfit + gvar):
   - 模型: R(T) = P_m(T) / Q_n(T) (展开中心 T0 默认 157.0 MeV)。
   - 支持 Padé [1/2], Padé [0/2], Padé [2/2] 等标准阶数。
   - 对每个 Jackknife 重采样统计样本 k 独立求解最优拟合参数。
4. 物理量提取:
   - 峰值位置 T_{peak}^{(k)}: 实温度轴上的极值点。
   - 极点位置 T_{pole}^{(k)}: 分母 Q_n(T) = 0 的最近复根 (Re(T_{pole}) \\pm i Im(T_{pole}))。
   - 峰值高度 \\chi_{peak}^{(k)}: 最大磁化率估计。
5. 产物落盘:
   - output/LCP/pade_fit_summary.csv: 包含均值与 Jackknife 误差的总结表。
   - output/LCP/pade_jk_samples.csv: 包含所有 N 个 Jackknife 样本明细的完整 CSV。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_LCP_DIR
from src.parselqcdata.pade_fitter import (
    load_and_align_susceptibility_samples,
    fit_jackknife_susceptibility_pade,
)


def run_pade_fit_task(
    order: str = "3/2",
    t0: float = 157.0,
    binsize: int = 1,
    max_samples: int = 2000,
    csv_path: Path | None = None,
    output_dir: Path = OUTPUT_LCP_DIR,
) -> tuple[dict, list[dict]]:
    print(f"\n[LCP-TASK-04] 启动手征磁化率 Padé [{order}] 有理函数拟合与极点/峰值抽取...")
    print(f"  - 展开参考中心 T0: {t0} MeV")
    print(f"  - 样本量对齐: 每温度截取前 {max_samples} 个统计样本")
    print(f"  - Binsize 设定: {binsize} (先不用做 binsize 处理，binsize=1)")
    print(f"  - 产物输出目录: {output_dir}")

    # 1. 载入并对齐各温度数据
    temps, chi_matrix, errors, betas = load_and_align_susceptibility_samples(
        csv_path=csv_path,
        max_samples=max_samples,
        binsize=binsize,
    )

    print(f"[INFO] 成功载入 {len(temps)} 个温度点，样本矩阵维度: {chi_matrix.shape}")
    for t_val, b_val, err_val in zip(temps, betas, errors):
        print(f"  -> Beta {b_val:<6}: T = {t_val:6.2f} MeV, 磁化率标准误 = {err_val:.6e}")

    # 2. 执行 Jackknife Padé 拟合
    summary, records = fit_jackknife_susceptibility_pade(
        temps=temps,
        chi_matrix=chi_matrix,
        errors=errors,
        order=order,
        t0=t0,
    )

    print("\n--------------------------------------------------------------------------")
    print("【Padé 拟合核心结果 (Jackknife 统计推断)】:")
    print(f"  - 拟合通过率:               {summary['n_fitted']} / {summary['n_samples']} ({summary['n_fitted']/summary['n_samples']*100:.1f}%)")
    print(f"  - 平均 \\chi^2/dof:           {summary['chi2_dof_mean']:.3f}")
    print(f"  - 峰值温度 T_peak:           {summary['t_peak_mean']:.4f} +/- {summary['t_peak_err']:.4f} MeV")
    print(f"  - 极点实部 Re(T_pole):       {summary['t_pole_real_mean']:.4f} +/- {summary['t_pole_real_err']:.4f} MeV")
    print(f"  - 极点虚部 Im(T_pole):       {summary['t_pole_imag_mean']:.4f} +/- {summary['t_pole_imag_err']:.4f} MeV")
    print(f"  - 峰值磁化率 \\chi_peak:       {summary['chi_peak_mean']:.4f} +/- {summary['chi_peak_err']:.4f} GeV^2")
    print("--------------------------------------------------------------------------")

    # 3. 持久化输出
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_df = pl.DataFrame([summary])
    summary_path = output_dir / "pade_fit_summary.csv"
    summary_df.write_csv(summary_path)

    samples_df = pl.DataFrame(records)
    samples_path = output_dir / "pade_jk_samples.csv"
    samples_df.write_csv(samples_path)

    print(f"[OK] 拟合结果已成功保存至:\n  - 汇总报告: {summary_path}\n  - 样本明细: {samples_path}\n")
    return summary, records


def main() -> None:
    parser = argparse.ArgumentParser(description="LCP 子任务 4: 手征磁化率 Padé 有理函数逼近拟合工具")
    parser.add_argument("--order", type=str, default="3/2", choices=["0/2", "1/2", "2/2", "3/2", "0/4"], help="Padé 阶数 (默认: 3/2)")
    parser.add_argument("--t0", type=float, default=157.0, help="展开中心温度 T0 (MeV, 默认: 157.0)")
    parser.add_argument("--binsize", type=int, default=1, help="构型重抽样 binsize (默认: 1)")
    parser.add_argument("--max-samples", type=int, default=2000, help="截取对齐的样本总数 (默认: 2000)")
    parser.add_argument("--csv", type=Path, default=None, help="输入的磁化率样本 CSV 路径 (可选)")
    parser.add_argument("--out-dir", type=Path, default=OUTPUT_LCP_DIR, help="结果输出目录 (默认: output/LCP)")
    args = parser.parse_args()

    run_pade_fit_task(
        order=args.order,
        t0=args.t0,
        binsize=args.binsize,
        max_samples=args.max_samples,
        csv_path=args.csv,
        output_dir=args.out_dir,
    )


if __name__ == "__main__":
    main()
