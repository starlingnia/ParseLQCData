#!/usr/bin/env python3
"""
scripts/plot_susceptibility.py
--------------------------------------------------------------------------------
手征磁化率 (Chiral Susceptibility) 高精度学术制图脚本 (纯轻夸克无偏估计，不执行 RM 减除)
严格复现 ana/ttest/ 的画图标准 (颜色、图例、线型、标注、误差线)：

1. sucep_plot.png / .pdf:
   - 提取自 results_susceptibility.csv
   - 纵轴为 4D 体积标度磁化率 F_vol * chi (F_vol = Ns^3 * Nt = 48^3 * 16)
   - 严格对应 ana/ttest/sucep_calc.py:
     - 颜色: 红色系 #ef4444, ecolor: #fca5a5, mfc: #f87171, mec: #ef4444
     - 样式: mew=1.5, ms=7, capsize=4, elinewidth=1.5, fmt='o-'
     - 标注: β={beta}
     - 标题: Disconnected Chiral Susceptibility \\chi vs Temperature

2. pbpchisce.png / .pdf:
   - 温度标度磁化率: F_scaled * chi (F_scaled = Ns^3 * Nt^3 * T^2)
   - 纵轴为 Scaled Susceptibility ($(N_t T)^2 \\times a^2\\chi$) [MeV^2]
   - 严格对应 ana/ttest/plot_chisce.py
   - 渲染相变过渡带 [155, 158] MeV（156.5 ± 1.5 MeV）与 T_pc ≈ 157.0 MeV 极大值峰位
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    FIGURES_DIR,
    OUTPUT_CONDENSATE_DIR,
    OUTPUT_ROOT,
    OUTPUT_LCP_DIR,
    TRANSITION_REGION,
)

ANA_STYLE_DIR = OUTPUT_ROOT / "plots" / "ana_style"
TARGET_DIRS = [FIGURES_DIR, ANA_STYLE_DIR, OUTPUT_CONDENSATE_DIR, OUTPUT_LCP_DIR]


def save_plot_multiformat(fig: plt.Figure, base_name: str) -> None:
    """自动将图表以 300 DPI 高清 PNG 和矢量 PDF 格式同步保存至所有目标目录"""
    for d in TARGET_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        png_path = d / f"{base_name}.png"
        pdf_path = d / f"{base_name}.pdf"
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
    print(f"  [EXPORT] 已同步保存 {base_name}.png/.pdf 至 {len(TARGET_DIRS)} 个目录")


def plot_sucep(df: pl.DataFrame) -> None:
    """严格按照 ana/ttest/sucep_calc.py 格式绘制手征磁化率曲线 (sucep_plot)"""
    print("\n--- 绘制手征磁化率 (sucep_plot) ---")
    temps = df["Temp"].to_list()
    means = df["Mean_vol_scaled"].to_list()
    errs = df["Error_vol_scaled"].to_list()
    betas_lbl = df["Beta"].to_list()

    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    ax.errorbar(
        temps,
        means,
        yerr=errs,
        fmt="o-",
        color="#ef4444",
        ecolor="#fca5a5",
        mfc="#f87171",
        mec="#ef4444",
        mew=1.5,
        ms=7,
        capsize=4,
        elinewidth=1.5,
        label="Chiral Susceptibility",
    )

    for i, beta in enumerate(betas_lbl):
        ax.annotate(
            f"β={beta}",
            (temps[i], means[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=9,
            fontweight="semibold",
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel(r"Disconnected Chiral Susceptibility ($\chi$)", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title(r"Disconnected Chiral Susceptibility $\chi$ vs Temperature", fontsize=13, fontweight="bold", pad=15)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper left")
    plt.tight_layout()

    save_plot_multiformat(fig, "sucep_plot")
    plt.close(fig)


def plot_scaled_pbpchisce(df: pl.DataFrame) -> None:
    """严格按照规范绘制非连通手征磁化率 (pbpchisce) [GeV^2]"""
    print("\n--- 绘制非连通手征磁化率 (pbpchisce) [GeV^2] ---")
    temps = np.array(df["Temp"])
    betas = df["Beta"].to_list()
    # 转换为物理量纲 [GeV^2]
    disc_means = np.array(df["Mean_scaled"]) / 1e6
    disc_errs = np.array(df["Error_scaled"]) / 1e6

    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=150)

    # 1. 渲染相变过渡温带 [155, 158] MeV (156.5 ± 1.5 MeV)
    ax.axvspan(
        TRANSITION_REGION[0],
        TRANSITION_REGION[1],
        color="crimson",
        alpha=0.12,
        label=f"Crossover Band ({TRANSITION_REGION[0]}-{TRANSITION_REGION[1]} MeV)",
        zorder=1,
    )

    # 2. 绘制非连通手征磁化率 (严格使用 ana/ttest/plot_chisce.py 配色与规范)
    color_indigo = "#4f46e5"
    ax.errorbar(
        temps,
        disc_means,
        yerr=disc_errs,
        fmt="o-",
        color=color_indigo,
        ecolor="#a5b4fc",
        mfc="#6366f1",
        mec=color_indigo,
        mew=1.5,
        ms=7,
        capsize=4,
        elinewidth=1.5,
        label=r"Disconnected Susceptibility $\chi_{\mathrm{disc}}$",
        zorder=3,
    )

    for i, beta in enumerate(betas):
        ax.annotate(
            f"β={beta}",
            (temps[i], disc_means[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=9,
            fontweight="semibold",
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel(r"Disconnected Susceptibility $\chi_{\mathrm{disc}}$ [$\mathrm{GeV}^2$]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title(r"Disconnected Chiral Susceptibility $\chi_{\mathrm{disc}}$ vs Temperature", fontsize=13, fontweight="bold", pad=15)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()

    save_plot_multiformat(fig, "pbpchisce")
    plt.close(fig)


def plot_scaling_susceptibility(df: pl.DataFrame) -> None:
    """绘制固定 beta=4.17 的有限时间延展系综非连通磁化率 (pbpchisce_scaling_beta4.17) [GeV^2]"""
    print("\n--- 绘制标度系综 (beta=4.17) 非连通手征磁化率 (pbpchisce_scaling_beta4.17) [GeV^2] ---")
    if "Ensemble" in df.columns:
        df_scaling = (
            df.filter(pl.col("Ensemble").str.contains(r"4\.17|m0\.0020|beta4\.17"))
              .filter(~((pl.col("Ns") == 48) & (pl.col("Nt") == 16)))
              .sort("Temp")
        )
    else:
        df_scaling = (
            df.filter(pl.col("Beta").cast(pl.String).str.contains(r"4\.17"))
              .filter(~((pl.col("Ns") == 48) & (pl.col("Nt") == 16)))
              .sort("Temp")
        )

    if df_scaling.height == 0:
        return

    temps = np.array(df_scaling["Temp"])
    disc_means = np.array(df_scaling["Mean_scaled"]) / 1e6
    disc_errs = np.array(df_scaling["Error_scaled"]) / 1e6
    ensembles = df_scaling["Ensemble"].to_list()
    ns_list = df_scaling["Ns"].to_list()
    nt_list = df_scaling["Nt"].to_list()

    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=150)

    # 1. 渲染相变过渡温带 [155, 158] MeV (156.5 ± 1.5 MeV)
    ax.axvspan(
        TRANSITION_REGION[0],
        TRANSITION_REGION[1],
        color="crimson",
        alpha=0.12,
        label=f"Crossover Band ({TRANSITION_REGION[0]}-{TRANSITION_REGION[1]} MeV)",
        zorder=1,
    )

    # 2. 绘制非连通手征磁化率
    color_indigo = "#4f46e5"
    ax.errorbar(
        temps,
        disc_means,
        yerr=disc_errs,
        fmt="o-",
        color=color_indigo,
        ecolor="#a5b4fc",
        mfc="#6366f1",
        mec=color_indigo,
        mew=1.5,
        ms=7,
        capsize=4,
        elinewidth=1.5,
        label=r"$\chi_{\mathrm{disc}}$ ($\beta=4.17$ series)",
        zorder=3,
    )

    # 标注时空尺寸标签
    for i in range(len(ensembles)):
        tag = f"${ns_list[i]}^3 \\times {nt_list[i]}$"
        if nt_list[i] == 16:
            offset = (0, -18)
        elif nt_list[i] == 14:
            offset = (0, 10)
        else:
            offset = (0, 10)

        ax.annotate(
            tag,
            (temps[i], disc_means[i]),
            textcoords="offset points",
            xytext=offset,
            ha="center",
            fontsize=9,
            fontweight="semibold",
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel(r"Disconnected Susceptibility $\chi_{\mathrm{disc}}$ [$\mathrm{GeV}^2$]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title(r"Disconnected Chiral Susceptibility vs Temperature ($\beta=4.17$)", fontsize=13, fontweight="bold", pad=15)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()

    save_plot_multiformat(fig, "pbpchisce_scaling_beta4.17")
    plt.close(fig)


def plot_renormalized_pbpchisce(df: pl.DataFrame) -> None:
    """绘制经 Zm 质量重整化后的非连通手征磁化率 (pbpchisce_renormalized) [GeV^2]"""
    print("\n--- 绘制重整化非连通手征磁化率 (pbpchisce_renormalized) [GeV^2] ---")
    temps = np.array(df["Temp"])
    betas = df["Beta"].to_list()

    if "Mean_scaled_gev2_renorm" in df.columns:
        disc_means = np.array(df["Mean_scaled_gev2_renorm"])
        disc_errs = np.array(df["Error_scaled_gev2_renorm"])
    else:
        from src.parselqcdata.susceptibility_pipeline import get_zm_factor
        zms = np.array([get_zm_factor(b) for b in betas])
        disc_means = (np.array(df["Mean_scaled"]) / 1e6) / zms
        disc_errs = (np.array(df["Error_scaled"]) / 1e6) / zms

    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=150)

    # 1. 渲染相变过渡温带 [155, 158] MeV
    ax.axvspan(
        TRANSITION_REGION[0],
        TRANSITION_REGION[1],
        color="crimson",
        alpha=0.12,
        label=f"Crossover Band ({TRANSITION_REGION[0]}-{TRANSITION_REGION[1]} MeV)",
        zorder=1,
    )

    color_indigo = "#4f46e5"
    ax.errorbar(
        temps,
        disc_means,
        yerr=disc_errs,
        fmt="o-",
        color=color_indigo,
        ecolor="#a5b4fc",
        mfc="#6366f1",
        mec=color_indigo,
        mew=1.5,
        ms=7,
        capsize=4,
        elinewidth=1.5,
        label=r"Renormalized $\chi_{\mathrm{disc}} / Z_m$",
        zorder=3,
    )

    for i, beta in enumerate(betas):
        ax.annotate(
            f"β={beta}",
            (temps[i], disc_means[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=9,
            fontweight="semibold",
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel(r"Renormalized Susceptibility $\chi_{\mathrm{disc}} / Z_m$ [$\mathrm{GeV}^2$]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title(r"Renormalized Chiral Susceptibility $\chi_{\mathrm{disc}} / Z_m$ vs Temperature", fontsize=13, fontweight="bold", pad=15)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()

    save_plot_multiformat(fig, "pbpchisce_renormalized")
    plt.close(fig)


def main():
    print("=" * 80)
    print("       ParseLQCData 手征磁化率 (Chiral Susceptibility) 绘图脚本       ")
    print("=" * 80)

    csv_path = OUTPUT_CONDENSATE_DIR / "results_susceptibility.csv"

    if not csv_path.exists():
        print("[ERROR] 缺少磁化率输入数据，请先运行 scripts/reproduce_susceptibility.py")
        sys.exit(1)

    df = pl.read_csv(csv_path)
    df_48 = df.filter((pl.col("Ns") == 48) & (pl.col("Nt") == 16)).sort("Temp")

    plot_sucep(df_48)
    plot_scaled_pbpchisce(df_48)
    plot_renormalized_pbpchisce(df_48)

    ensemble_csv_path = OUTPUT_CONDENSATE_DIR / "all_ensembles_susceptibility.csv"
    if ensemble_csv_path.exists():
        fixed_beta_df = pl.read_csv(ensemble_csv_path)
    else:
        fixed_beta_df = df
    plot_scaling_susceptibility(fixed_beta_df)

    print("\n所有磁化率图表已成功生成并同步至:")
    for d in TARGET_DIRS:
        print(f"  - {d}")


if __name__ == "__main__":
    main()
