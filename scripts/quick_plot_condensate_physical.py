#!/usr/bin/env python3
"""
scripts/quick_plot_condensate_physical.py
--------------------------------------------------------------------------------
快速绘制带真实连续物理量纲 [MeV^3] 的手征凝聚随温度演化图：
  Delta_phys = (16 T)^3 * Delta_{l, s}  [MeV^3]
严格遵循 ana/ 风格：
- royalblue 配色
- 与 ana/ 一致的坐标轴格式与标题
- 小方框图例指示栏位于右上角 (loc="upper right")
- 渲染相变过渡温带 [155.5, 160.5] MeV
--------------------------------------------------------------------------------
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 1. 物理参数与输入文件
OUTPUT_CONDENSATE_DIR = PROJECT_ROOT / "output" / "condensate"
DOCS_FIGURES_DIR = PROJECT_ROOT / "docs" / "figures"
DOCS_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

CONDENSATE_FILE = OUTPUT_CONDENSATE_DIR / "results_rm_beta.txt"
SCALING_FILE = OUTPUT_CONDENSATE_DIR / "scaling_beta4.17.csv"

TEMP_MAP = {
    "4.13": 138.23,
    "4.15": 145.75,
    "4.17": 153.31,
    "4.18": 157.03,
    "4.20": 164.55,
    "4.23": 176.43,
    "4.30": 202.52,
    "4.405": 241.60,
}

TRANSITION_REGION = (155.5, 160.5)

def main():
    if not CONDENSATE_FILE.exists():
        print(f"[ERROR] 找不到输入文件: {CONDENSATE_FILE}")
        return

    df = pl.read_csv(CONDENSATE_FILE, separator="\t")
    df = df.with_columns(
        pl.col("beta").cast(pl.String).alias("beta_str")
    ).with_columns(
        pl.col("beta_str").map_elements(
            lambda b: TEMP_MAP.get(str(b), TEMP_MAP.get(f"{float(b):.2f}", np.nan)),
            return_dtype=pl.Float64
        ).alias("T")
    ).sort("T")

    temps = df["T"].to_numpy()
    pbp_rm = df["pbp_rm"].to_numpy()
    pbp_rm_err = df["pbp_rm_err"].to_numpy()
    betas = df["beta_str"].to_list()

    # 乘以正确的物理因子 a^-3 = (16 * T)^3  [MeV^3]
    scale_factor_3 = (16.0 * temps) ** 3
    pbp_phys = pbp_rm * scale_factor_3
    pbp_phys_err = pbp_rm_err * scale_factor_3

    # =========================================================================
    # 图 1: 物理凝聚绝对值 (16T)^3 * Delta_{l, s} [MeV^3]
    # =========================================================================
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

    # 灰色相变过渡温带 (ana/ 标准)
    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region", zorder=1)

    # 主曲线 (采用 ana/ 皇家蓝 royalblue 与经典标记)
    ax.errorbar(
        temps,
        pbp_phys,
        yerr=pbp_phys_err,
        fmt="o-",
        color="royalblue",
        ecolor="royalblue",
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=6,
        label=r'$\langle\bar{\psi}\psi\rangle_{\mathrm{phys}} = (16T)^3 \Delta_{\ell, s}$',
        zorder=5
    )

    # 标注每个温度点的 beta 值
    for i, b in enumerate(betas):
        ax.annotate(
            f"β={b}",
            (temps[i], pbp_phys[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=9,
            fontweight="medium"
        )

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14, labelpad=10)
    ax.set_ylabel(r"Physical Condensate $\langle\bar{\psi}\psi\rangle_{\mathrm{phys}}$ [$\mathrm{MeV}^3$]", fontsize=14, labelpad=10)
    ax.set_title(r"Physical Chiral Condensate vs Temperature ($48^3 \times 16$)", fontsize=14, pad=15)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    
    # 严格将指示栏挪到右上角
    ax.legend(loc="upper right", frameon=True, fontsize=11, framealpha=0.9)
    plt.tight_layout()

    out_mev3 = DOCS_FIGURES_DIR / "conden_physical_mev3.png"
    out_mev3_pdf = DOCS_FIGURES_DIR / "conden_physical_mev3.pdf"
    fig.savefig(out_mev3, dpi=300, bbox_inches="tight")
    fig.savefig(out_mev3_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] 已生成物理单位手性凝聚图: {out_mev3} 及 {out_mev3_pdf}")

    # =========================================================================
    # 图 2: 凝聚开三次方有效能标 [MeV] (如 <pbp>^(1/3) ~ 225.56 MeV)
    # =========================================================================
    fig2, ax2 = plt.subplots(figsize=(8, 6), dpi=150)
    ax2.axvspan(t_low, t_high, color="gray", alpha=0.2, label="Transition Region", zorder=1)

    pbp_cuberoot = np.cbrt(np.maximum(pbp_phys, 0))
    # 传播误差: d(y^(1/3)) = (1/3) * y^(-2/3) * dy
    pbp_cuberoot_err = (1.0 / 3.0) * (np.maximum(pbp_phys, 1e-12) ** (-2.0 / 3.0)) * pbp_phys_err

    ax2.errorbar(
        temps,
        pbp_cuberoot,
        yerr=pbp_cuberoot_err,
        fmt="o-",
        color="royalblue",
        ecolor="royalblue",
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=6,
        label=r'$(\langle\bar{\psi}\psi\rangle_{\mathrm{phys}})^{1/3} = ((16T)^3 \Delta_{\ell, s})^{1/3}$',
        zorder=5
    )

    for i, b in enumerate(betas):
        ax2.annotate(
            f"β={b}\n({pbp_cuberoot[i]:.1f})",
            (temps[i], pbp_cuberoot[i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=8,
            fontweight="medium"
        )

    ax2.set_xlabel("Temperature (T) [MeV]", fontsize=14, labelpad=10)
    ax2.set_ylabel(r"Condensate Scale $(\langle\bar{\psi}\psi\rangle_{\mathrm{phys}})^{1/3}$ [MeV]", fontsize=14, labelpad=10)
    ax2.set_title(r"Physical Chiral Condensate Scale vs Temperature", fontsize=14, pad=15)
    ax2.tick_params(axis="both", which="major", labelsize=12)
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(loc="upper right", frameon=True, fontsize=11, framealpha=0.9)
    plt.tight_layout()

    out_cbrt = DOCS_FIGURES_DIR / "conden_physical_cuberoot.png"
    out_cbrt_pdf = DOCS_FIGURES_DIR / "conden_physical_cuberoot.pdf"
    fig2.savefig(out_cbrt, dpi=300, bbox_inches="tight")
    fig2.savefig(out_cbrt_pdf, bbox_inches="tight")
    plt.close(fig2)
    print(f"[OK] 已生成开三次方手性凝聚能标图: {out_cbrt} 及 {out_cbrt_pdf}")

if __name__ == "__main__":
    main()
