#!/usr/bin/env python3
"""
scripts/tasks/lcp/05_plot_pade_susceptibility.py
--------------------------------------------------------------------------------
手征磁化率 Padé 拟合精简绘图脚本:
仅展示:
1. 原数据点及误差棒
2. 拟合的平均曲线
3. 阴影表示拟合范围
4. 平均曲线的约化卡方 chi^2 / dof
5. Tc (峰值 T, T_peak) 的中心值与 Jackknife 偏差估计
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import FIGURES_DIR, OUTPUT_LCP_DIR


def plot_pade_susceptibility(
    lcp_dir: Path = OUTPUT_LCP_DIR,
    out_dir: Path = FIGURES_DIR,
    order: str = "2/2",
    t0: float = 157.0,
) -> Path:
    summary_file = lcp_dir / "pade_fit_summary.csv"
    samples_file = lcp_dir / "pade_jk_samples.csv"
    results_file = lcp_dir / "results_susceptibility.csv"

    if not summary_file.exists() or not samples_file.exists() or not results_file.exists():
        raise FileNotFoundError(
            f"缺少拟合数据产物，请先运行 scripts/tasks/lcp/04_fit_pade_susceptibility.py"
        )

    # 1. 读取拟合统计量与原数据点
    df_sum = pl.read_csv(summary_file)
    t_peak_mean = float(df_sum["t_peak_mean"][0])
    t_peak_err = float(df_sum["t_peak_err"][0])
    chi2_dof = float(df_sum["chi2_dof_mean"][0])
    t0_val = float(df_sum["t0"][0]) if "t0" in df_sum.columns else t0

    df_pts = pl.read_csv(results_file)
    if "Ns" in df_pts.columns:
        df_pts = df_pts.filter(pl.col("Ns") == 48)
    if "Nt" in df_pts.columns:
        df_pts = df_pts.filter(pl.col("Nt") == 16)
    df_pts = df_pts.sort("Temp")
    temps_data = df_pts["Temp"].to_numpy().astype(float)
    chi_data = df_pts["Mean_scaled_gev2_renorm"].to_numpy().astype(float)
    err_data = df_pts["Error_scaled_gev2_renorm"].to_numpy().astype(float)

    # 2. 读取 Jackknife 样本计算平均曲线与拟合范围阴影
    df_samples = pl.read_csv(samples_file).filter(pl.col("ok") == True)
    n_samples = df_samples.height

    t_grid = np.linspace(133.0, 247.0, 500)
    dx_grid = t_grid - t0_val

    num = np.zeros((n_samples, len(t_grid)), dtype=float)
    for p_i in range(5):
        col = f"p{p_i}"
        if col in df_samples.columns:
            pv = df_samples[col].to_numpy().astype(float)
            num += pv[:, None] * (dx_grid[None, :] ** p_i)

    den = np.ones((n_samples, len(t_grid)), dtype=float)
    for q_i in range(1, 5):
        col = f"q{q_i}"
        if col in df_samples.columns:
            qv = df_samples[col].to_numpy().astype(float)
            den += qv[:, None] * (dx_grid[None, :] ** q_i)

    r_matrix = num / den
    curve_mean = np.mean(r_matrix, axis=0)
    curve_err = np.sqrt(n_samples - 1) * np.std(r_matrix, axis=0, ddof=0)

    # 3. 干净绘图
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "mathtext.fontset": "cm",
    })

    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    # ① 阴影表示拟合的范围 (Fit range)
    ax.fill_between(
        t_grid,
        curve_mean - curve_err,
        curve_mean + curve_err,
        color="#1f77b4",
        alpha=0.25,
        label="Fit range",
        zorder=1,
    )

    # ② 拟合的平均曲线 (Fit mean curve)
    ax.plot(
        t_grid,
        curve_mean,
        color="#1f77b4",
        linewidth=2.0,
        label="Fit",
        zorder=2,
    )

    # ③ 原数据点及误差棒 (Data points)
    ax.errorbar(
        temps_data,
        chi_data,
        yerr=err_data,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.4,
        capsize=3.5,
        markersize=5.5,
        label="Data",
        zorder=3,
    )

    # ④ 图例置于右上角
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # ⑤ 数据方框放置在右上角图例下方
    info_text = (
        rf"$T_c = {t_peak_mean:.2f} \pm {t_peak_err:.2f}\ \mathrm{{MeV}}$" + "\n"
        rf"$\chi^2/\mathrm{{dof}} = {chi2_dof:.3f}$"
    )
    ax.text(
        0.96,
        0.72,
        info_text,
        transform=ax.transAxes,
        fontsize=11.5,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        zorder=4,
    )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax.set_xlim(133.0, 247.0)
    ax.set_ylim(-0.02, 1.15)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.tight_layout()

    # 4. 保存产物
    out_dir.mkdir(parents=True, exist_ok=True)
    lcp_dir.mkdir(parents=True, exist_ok=True)

    fig_png = out_dir / "chiral_susceptibility_pade_fit.png"
    fig_pdf = out_dir / "chiral_susceptibility_pade_fit.pdf"
    lcp_png = lcp_dir / "chiral_susceptibility_pade_fit.png"

    plt.savefig(fig_png, dpi=300)
    plt.savefig(fig_pdf)
    plt.savefig(lcp_png, dpi=300)

    # 别名保存
    order_tag = df_sum["order"][0].replace("/", "") if "order" in df_sum.columns else "22"
    plt.savefig(out_dir / f"chiral_susceptibility_pade{order_tag}_fit.png", dpi=300)
    plt.savefig(out_dir / f"chiral_susceptibility_pade{order_tag}_fit.pdf")
    plt.close(fig)

    print(f"[OK] 成功生成拟合图:\n  - {fig_png}\n  - {fig_pdf}\n  - {lcp_png}")
    return fig_png


def main() -> None:
    parser = argparse.ArgumentParser(description="绘制手征磁化率 Padé 拟合图")
    parser.add_argument("--lcp-dir", type=Path, default=OUTPUT_LCP_DIR, help="LCP 数据目录 (默认: output/LCP)")
    parser.add_argument("--out-dir", type=Path, default=FIGURES_DIR, help="图表输出目录 (默认: docs/figures)")
    args = parser.parse_args()

    plot_pade_susceptibility(lcp_dir=args.lcp_dir, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
