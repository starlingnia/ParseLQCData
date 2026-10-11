#!/usr/bin/env python3
"""
scripts/tasks/lcp/05_plot_pade_susceptibility.py
--------------------------------------------------------------------------------
绘制手征磁化率 Padé [3/2] 拟合平均曲线、Jackknife 误差阴影带与赝临界峰值区间
--------------------------------------------------------------------------------
【包含图层】:
1. 8 组格点 QCD 原始测量数据点与误差棒 (包含实际构型统计标准误)。
2. Padé [3/2] 有理函数的均值拟合曲线。
3. 拟合曲线的完整逐样本 Jackknife 误差阴影带 (1-sigma 预测带)。
4. 峰值温度 T_peak 与极点位置 Re(T_pole) 的物理置信区间阴影标记。
5. 局部过渡带高分辨率放大子图 (Inset Zoom)，清晰展示相变峰顶微观结构。
--------------------------------------------------------------------------------
"""

from pathlib import Path
import sys

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_LCP_DIR, FIGURES_DIR


def plot_pade_susceptibility(
    lcp_dir: Path = OUTPUT_LCP_DIR,
    out_dir: Path = FIGURES_DIR,
    order: str = "3/2",
    t0: float = 157.0,
) -> Path:
    summary_file = lcp_dir / "pade_fit_summary.csv"
    samples_file = lcp_dir / "pade_jk_samples.csv"
    results_file = lcp_dir / "results_susceptibility.csv"

    if not summary_file.exists() or not samples_file.exists() or not results_file.exists():
        raise FileNotFoundError(f"缺少拟合数据产物，请先运行 04_fit_pade_susceptibility.py")

    # 1. 载入汇总数据与数据点
    df_sum = pl.read_csv(summary_file)
    t_peak_mean = float(df_sum["t_peak_mean"][0])
    t_peak_err = float(df_sum["t_peak_err"][0])
    pole_re_mean = float(df_sum["t_pole_real_mean"][0])
    pole_re_err = float(df_sum["t_pole_real_err"][0])
    pole_im_mean = float(df_sum["t_pole_imag_mean"][0])
    pole_im_err = float(df_sum["t_pole_imag_err"][0])
    chi_peak_mean = float(df_sum["chi_peak_mean"][0])
    chi_peak_err = float(df_sum["chi_peak_err"][0])
    chi2_dof = float(df_sum["chi2_dof_mean"][0])

    df_pts = pl.read_csv(results_file).sort("Temp")
    temps_data = df_pts["Temp"].to_numpy().astype(float)
    chi_data = df_pts["Mean_scaled_gev2_renorm"].to_numpy().astype(float)
    err_data = df_pts["Error_scaled_gev2_renorm"].to_numpy().astype(float)

    # 2. 载入 2000 个 Jackknife 拟合参数，并在密集温度网格上计算曲线与误差带
    df_samples = pl.read_csv(samples_file).filter(pl.col("ok") == True)
    n_samples = df_samples.height

    t_grid = np.linspace(135.0, 245.0, 600)
    dx_grid = t_grid - t0

    # 提取参数矩阵
    p0 = df_samples["p0"].to_numpy().astype(float)
    p1 = df_samples["p1"].to_numpy().astype(float)
    p2 = df_samples["p2"].to_numpy().astype(float)
    p3 = df_samples["p3"].to_numpy().astype(float)
    q1 = df_samples["q1"].to_numpy().astype(float)
    q2 = df_samples["q2"].to_numpy().astype(float)

    # 计算每个样本在每个温度点上的值: (n_samples, n_grid)
    # R(T) = (p0 + p1*dx + p2*dx^2 + p3*dx^3) / (1 + q1*dx + q2*dx^2)
    num = (
        p0[:, None]
        + p1[:, None] * dx_grid[None, :]
        + p2[:, None] * (dx_grid[None, :] ** 2)
        + p3[:, None] * (dx_grid[None, :] ** 3)
    )
    den = (
        1.0
        + q1[:, None] * dx_grid[None, :]
        + q2[:, None] * (dx_grid[None, :] ** 2)
    )
    r_matrix = num / den

    # Jackknife 均值与方差
    curve_mean = np.mean(r_matrix, axis=0)
    curve_err = np.sqrt(n_samples - 1) * np.std(r_matrix, axis=0, ddof=0)

    # 3. 开始绘图 (双面板结构: 主全景图 + 局部放大图)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10.5,
        "figure.titlesize": 15,
        "mathtext.fontset": "cm",
    })

    fig = plt.figure(figsize=(10.5, 6.5), dpi=300)
    ax = fig.add_subplot(111)

    # 拟合置信带阴影 (1-sigma Jackknife Error Band)
    ax.fill_between(
        t_grid,
        curve_mean - curve_err,
        curve_mean + curve_err,
        color="#1f77b4",
        alpha=0.28,
        label=r"Padé [3/2] 1$\sigma$ Jackknife Band",
    )

    # 拟合均值平滑曲线
    ax.plot(
        t_grid,
        curve_mean,
        color="#0d47a1",
        linewidth=2.2,
        label="Mean Padé [3/2] Fit Curve",
        zorder=3,
    )

    # 峰值温度 T_peak 标记与置信带阴影
    ax.axvspan(
        t_peak_mean - t_peak_err,
        t_peak_mean + t_peak_err,
        color="#d32f2f",
        alpha=0.22,
        label=rf"$T_{{\mathrm{{peak}}}} = {t_peak_mean:.2f} \pm {t_peak_err:.2f}$ MeV",
        zorder=1,
    )
    ax.axvline(
        t_peak_mean,
        color="#b71c1c",
        linestyle="--",
        linewidth=1.8,
        alpha=0.85,
        zorder=2,
    )

    # 极点实部 Re(T_pole) 标记
    ax.axvline(
        pole_re_mean,
        color="#2e7d32",
        linestyle=":",
        linewidth=1.8,
        alpha=0.85,
        label=rf"$\mathrm{{Re}}(T_{{\mathrm{{pole}}}}) = {pole_re_mean:.2f} \pm {pole_re_err:.2f}$ MeV",
        zorder=2,
    )

    # 格点原始数据点
    ax.errorbar(
        temps_data,
        chi_data,
        yerr=err_data,
        fmt="o",
        color="#212121",
        ecolor="#212121",
        elinewidth=1.6,
        capsize=4.5,
        capthick=1.4,
        markersize=6.5,
        label=r"Lattice QCD Data ($N=2000$)",
        zorder=5,
    )

    # 物理参数说明文本框
    info_text = (
        r"$\bf{Pad\acute{e}\ [3/2]\ Model:}$" + "\n"
        r"$R(T) = \frac{p_0 + p_1 \Delta T + p_2 \Delta T^2 + p_3 \Delta T^3}{1 + q_1 \Delta T + q_2 \Delta T^2}$" + "\n"
        rf"$\Delta T = T - 157.0$ MeV" + "\n"
        rf"$\chi^2 / \mathrm{{dof}} = {chi2_dof:.3f}$" + "\n"
        rf"$\mathrm{{Im}}(T_{{\mathrm{{pole}}}}) = {pole_im_mean:.2f} \pm {pole_im_err:.2f}$ MeV" + "\n"
        rf"$\chi_{{\mathrm{{peak}}}} = {chi_peak_mean:.4f} \pm {chi_peak_err:.4f}\ \mathrm{{GeV}}^2$"
    )
    ax.text(
        0.965,
        0.58,
        info_text,
        transform=ax.transAxes,
        fontsize=9.8,
        verticalalignment="center",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#ffffff", edgecolor="#cccccc", alpha=0.92),
        zorder=6,
    )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax.set_title("Chiral Susceptibility Padé [3/2] Rational Fit & Peak/Pole Extraction", pad=12)

    ax.set_xlim(133.0, 247.0)
    ax.set_ylim(-0.03, 1.12)
    ax.grid(True, linestyle="--", alpha=0.45)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="#cccccc")

    # --------------------------------------------------------------------------
    # 局部相变峰顶放大 Inset (145 MeV ~ 172 MeV)
    # --------------------------------------------------------------------------
    ax_inset = ax.inset_axes([0.15, 0.12, 0.40, 0.38])
    sub_mask = (t_grid >= 144.0) & (t_grid <= 172.0)
    ax_inset.fill_between(
        t_grid[sub_mask],
        (curve_mean - curve_err)[sub_mask],
        (curve_mean + curve_err)[sub_mask],
        color="#1f77b4",
        alpha=0.32,
    )
    ax_inset.plot(t_grid[sub_mask], curve_mean[sub_mask], color="#0d47a1", linewidth=2.0)

    # 峰值与极点垂直阴影
    ax_inset.axvspan(
        t_peak_mean - t_peak_err,
        t_peak_mean + t_peak_err,
        color="#d32f2f",
        alpha=0.22,
    )
    ax_inset.axvline(t_peak_mean, color="#b71c1c", linestyle="--", linewidth=1.6)
    ax_inset.axvline(pole_re_mean, color="#2e7d32", linestyle=":", linewidth=1.6)

    # 局部范围数据点
    sub_pts = (temps_data >= 143.0) & (temps_data <= 173.0)
    ax_inset.errorbar(
        temps_data[sub_pts],
        chi_data[sub_pts],
        yerr=err_data[sub_pts],
        fmt="o",
        color="#212121",
        elinewidth=1.4,
        capsize=3.5,
        markersize=5.5,
    )

    ax_inset.set_xlim(144.0, 172.0)
    ax_inset.set_ylim(0.70, 1.05)
    ax_inset.set_title("Transition Region Zoom", fontsize=9.5)
    ax_inset.tick_params(labelsize=8.5)
    ax_inset.grid(True, linestyle=":", alpha=0.5)
    ax.indicate_inset_zoom(ax_inset, edgecolor="#666666", alpha=0.7)

    plt.tight_layout()

    # 保存产物
    out_dir.mkdir(parents=True, exist_ok=True)
    lcp_dir.mkdir(parents=True, exist_ok=True)

    fig_png = out_dir / "chiral_susceptibility_pade32_fit.png"
    fig_pdf = out_dir / "chiral_susceptibility_pade32_fit.pdf"
    lcp_png = lcp_dir / "chiral_susceptibility_pade32_fit.png"

    plt.savefig(fig_png, dpi=300)
    plt.savefig(fig_pdf)
    plt.savefig(lcp_png, dpi=300)
    plt.close()

    print(f"[OK] 成功生成高清拟合效果图:\n  - {fig_png}\n  - {fig_pdf}\n  - {lcp_png}")
    return fig_png


if __name__ == "__main__":
    plot_pade_susceptibility()
