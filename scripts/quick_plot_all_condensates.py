#!/usr/bin/env python3
"""
scripts/quick_plot_all_condensates.py
--------------------------------------------------------------------------------
绘制所有 13 个系综手征凝聚值随温度演化的全局总览图 (All-in-One):
- 8组标准物理质量温扫系综 (48^3 x 16, m_l ≈ 0.0010): royalblue 圆形
- 5组 beta=4.17 有限时间标度系综 (m_l = 0.0020): 紫色正方形
严格遵循 ana/ 风格：
- 右上角小方框指示栏 (loc="upper right")
- 相变过渡带阴影 [155.5, 160.5] MeV
- 坐标轴标签与字体标准
--------------------------------------------------------------------------------
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_CONDENSATE_DIR = PROJECT_ROOT / "output" / "condensate"
DOCS_FIGURES_DIR = PROJECT_ROOT / "docs" / "figures"
DOCS_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

ALL_FILE = OUTPUT_CONDENSATE_DIR / "all_ensembles_condensate.csv"
TRANSITION_REGION = (155.5, 160.5)

def main():
    if not ALL_FILE.exists():
        print(f"[ERROR] 找不到文件: {ALL_FILE}")
        return

    df = pl.read_csv(ALL_FILE)

    # 1. 拆分为两组系综
    # 标准物理质量扫描 (48x16, ml <= 0.00105)
    std_df = df.filter((pl.col("ns") == 48) & (pl.col("nt") == 16)).sort("temperature")
    # beta=4.17 标度系综 (ml = 0.0020)
    scale_df = df.filter(~((pl.col("ns") == 48) & (pl.col("nt") == 16))).sort("temperature")

    fig, ax = plt.subplots(figsize=(9.5, 6.2), dpi=150)

    # 相变过渡带
    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color="gray", alpha=0.18, label="Transition Region (155.5–160.5 MeV)", zorder=1)

    # 绘制标准温扫系综 (royalblue)
    ax.errorbar(
        std_df["temperature"].to_numpy(),
        std_df["pbp_rm"].to_numpy(),
        yerr=std_df["pbp_rm_err"].to_numpy(),
        fmt="o-",
        color="royalblue",
        ecolor="royalblue",
        mfc="#60a5fa",
        mec="royalblue",
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'Physical Quark Mass Scan ($48^3 \times 16$, $m_l \approx 0.0010$)',
        zorder=5
    )

    # 绘制 beta=4.17 标度系综 (紫色 #7c3aed)
    ax.errorbar(
        scale_df["temperature"].to_numpy(),
        scale_df["pbp_rm"].to_numpy(),
        yerr=scale_df["pbp_rm_err"].to_numpy(),
        fmt="s--",
        color="#7c3aed",
        ecolor="#7c3aed",
        mfc="#c084fc",
        mec="#7c3aed",
        capsize=5,
        capthick=1.5,
        elinewidth=1.5,
        markersize=7,
        label=r'Temporal Scaling Series ($\beta=4.17$, $m_l = 0.0020$)',
        zorder=6
    )

    # 标准温扫数据点标注 (beta)
    for row in std_df.iter_rows(named=True):
        b = f"{row['beta']:.2f}" if abs(row['beta'] - 4.405) > 1e-4 else "4.405"
        y_val = row["pbp_rm"]
        offset = (0, 9)
        if abs(row["temperature"] - 153.31) < 0.1:
            offset = (-18, 9)  # 避开与 40x16 碰撞
        ax.annotate(
            f"β={b}",
            (row["temperature"], y_val),
            textcoords="offset points",
            xytext=offset,
            ha="center",
            fontsize=8.5,
            color="#1e3a8a",
            fontweight="normal"
        )

    # 标度系综数据点标注 (时空尺寸)
    for row in scale_df.iter_rows(named=True):
        ns, nt = row["ns"], row["nt"]
        tag = f"${ns}^3\\times{nt}$"
        if "_2" in row["dataset_name"]:
            tag = r"$48^3\times18\text{ (s2)}$"
            offset = (22, 10)
        elif nt == 18:
            tag = r"$48^3\times18\text{ (s1)}$"
            offset = (-24, 10)
        elif nt == 16:
            tag = r"$40^3\times16$"
            offset = (24, 10)
        elif nt == 14:
            offset = (0, 10)
        elif nt == 12:
            offset = (0, -16)
        else:
            offset = (0, 10)

        ax.annotate(
            tag,
            (row["temperature"], row["pbp_rm"]),
            textcoords="offset points",
            xytext=offset,
            ha="center",
            fontsize=8.5,
            color="#581c87",
            fontweight="normal"
        )

    ax.axhline(0, color="black", linestyle=":", linewidth=0.8, alpha=0.5)

    ax.set_xlabel("Temperature (T) [MeV]", fontsize=14, labelpad=10)
    ax.set_ylabel(r"Renormalized Chiral Condensate $\Delta_{\ell, s}$", fontsize=14, labelpad=10)
    ax.set_title(r"Global Chiral Condensate $\Delta_{\ell, s}$ vs Temperature (All 13 Ensembles)", fontsize=14, pad=15)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", alpha=0.55)

    # 严格将小方框图例置于右上角
    ax.legend(loc="upper right", frameon=True, fontsize=10.5, framealpha=0.92)

    plt.tight_layout()

    out_png = DOCS_FIGURES_DIR / "conden_all_ensembles_overview.png"
    out_pdf = DOCS_FIGURES_DIR / "conden_all_ensembles_overview.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] 成功绘制并保存所有系综手征凝聚图: {out_png} 和 {out_pdf}")

if __name__ == "__main__":
    main()
