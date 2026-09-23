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
   - 连续温度标度磁化率: F_scaled * chi (F_scaled = Ns^3 * Nt^3 * T^2)
   - 纵轴为 Scaled Susceptibility ($(16T)^2 \\times a^2\\chi$) [MeV^2]
   - 严格对应 ana/ttest/plot_chisce.py
   - 渲染相变过渡带 [155.5, 160.5] MeV 与 T_pc ≈ 157.0 MeV 极大值峰位
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    FIGURES_DIR,
    OUTPUT_CONDENSATE_DIR,
    OUTPUT_ROOT,
    TRANSITION_REGION,
)

ANA_STYLE_DIR = OUTPUT_ROOT / "plots" / "ana_style"
TARGET_DIRS = [FIGURES_DIR, ANA_STYLE_DIR, OUTPUT_CONDENSATE_DIR]


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
    """严格按照 ana/ttest/plot_chisce.py 格式绘制连续标度化磁化率 (pbpchisce)"""
    print("\n--- 绘制标度化手征磁化率 (pbpchisce) ---")
    temps = np.array(df["Temp"])
    betas = df["Beta"].to_list()
    scaled_means = np.array(df["Mean_scaled"])
    scaled_errs = np.array(df["Error_scaled"])

    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=150)

    # 1. 渲染相变过渡温带 [155.5, 160.5] MeV
    ax.axvspan(
        TRANSITION_REGION[0],
        TRANSITION_REGION[1],
        color="crimson",
        alpha=0.12,
        label=f"Crossover Band ({TRANSITION_REGION[0]}-{TRANSITION_REGION[1]} MeV)",
        zorder=1,
    )

    # 2. 绘制连续标度磁化率 (使用 ana/ttest/plot_chisce.py 经典配色)
    ax.errorbar(
        temps,
        scaled_means,
        yerr=scaled_errs,
        fmt="o-",
        color="#ef4444",
        ecolor="#fca5a5",
        mfc="#f87171",
        mec="#ef4444",
        mew=1.5,
        ms=7,
        capsize=4,
        elinewidth=1.5,
        label="Chiral Susceptibility (Scaled)",
        zorder=3,
    )

    # 标注峰值点与各温度点
    peak_idx = int(np.argmax(scaled_means))
    ax.annotate(
        f"Peak at T={temps[peak_idx]:.1f} MeV\n(β={betas[peak_idx]}, $T_{{pc}}$)",
        (temps[peak_idx], scaled_means[peak_idx]),
        textcoords="offset points",
        xytext=(25, -5),
        fontsize=10,
        fontweight="bold",
        color="#991b1b",
        arrowprops=dict(arrowstyle="->", color="#991b1b", lw=1.5),
    )

    for i, beta in enumerate(betas):
        ax.annotate(
            f"β={beta}",
            (temps[i], scaled_means[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=9,
            fontweight="semibold",
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel(r"Scaled Susceptibility ($(16T)^2 \times a^2\chi$) [$\mathrm{MeV}^2$]", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title("Chiral Susceptibility vs Temperature (Scaled)", fontsize=13, fontweight="bold", pad=15)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper left", framealpha=0.9)
    plt.tight_layout()

    save_plot_multiformat(fig, "pbpchisce")
    plt.close(fig)


def main():
    print("=" * 80)
    print("       ParseLQCData 手征磁化率 (Chiral Susceptibility) 绘图脚本       ")
    print("=" * 80)

    csv_path = OUTPUT_CONDENSATE_DIR / "results_susceptibility.csv"

    if not csv_path.exists():
        print("[ERROR] 缺少磁化率输入数据，请先运行 scripts/reproduce_susceptibility.py")
        sys.exit(1)

    df = pl.read_csv(csv_path).sort("Temp")

    plot_sucep(df)
    plot_scaled_pbpchisce(df)

    print("\n所有磁化率图表已成功生成并同步至:")
    for d in TARGET_DIRS:
        print(f"  - {d}")


if __name__ == "__main__":
    main()
