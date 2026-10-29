#!/usr/bin/env python3
"""
scripts/tasks/lcp/04_plot_lcp.py
--------------------------------------------------------------------------------
LCP (常物理线) 标准极简绘图任务:
1. 手征易感度标准图 (无拟合曲线正常图) -> chiral_susceptibility.png / .pdf
2. 手征凝聚标准图 (同一极简学术风格)   -> chiral_condensate.png / .pdf
--------------------------------------------------------------------------------
【纯净极简统一风格规范】:
- 画布比例: figsize=(8.0, 5.5), dpi=300
- 字体系统: 经典无衬线 sans-serif, LaTeX 数学字体 cm, 标签字号 13, 刻度字号 11
- 数据呈现: 纯数据点加误差棒 (fmt='o', color='black', capsize=3.5)
- 刻度网格: 启用次刻度 AutoMinorLocator(), 浅灰虚线网格 alpha=0.4
- 图例框体: 简洁右上角 frameon=True, edgecolor='#cccccc'
- 绝无多余杂乱元素 (无多余放大框、无冗余公式说明、无每点文字堆叠)
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

from docs.physics_setup import (
    FIGURES_DIR,
    OUTPUT_CONDENSATE_DIR,
    OUTPUT_LCP_DIR,
    TEMP_MAP,
)


def apply_minimalist_style() -> None:
    """应用全局统一极简学术画风配置"""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "mathtext.fontset": "cm",
    })

def _save_figure(fig: plt.Figure, name: str, dirs: list[Path]) -> None:
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / f"{name}.png", dpi=300)
        fig.savefig(d / f"{name}.pdf")
    plt.close(fig)


def plot_chiral_susceptibility(
    lcp_dir: Path = OUTPUT_LCP_DIR,
    out_dir: Path = FIGURES_DIR,
) -> Path:
    """绘制没有拟合曲线的正常手征易感度数据图 (支持 48 尺寸、64 尺寸及对比图)"""
    results_file = lcp_dir / "results_susceptibility.csv"
    if not results_file.exists():
        raise FileNotFoundError(f"缺少易感度数据文件: {results_file}")

    df_pts = pl.read_csv(results_file)
    if "Nt" in df_pts.columns:
        df_pts = df_pts.filter(pl.col("Nt") == 16)

    # 提取 48 尺寸与 64 尺寸子集
    df_48 = df_pts.filter(pl.col("Ns") == 48).sort("Temp") if "Ns" in df_pts.columns else df_pts.sort("Temp")
    df_64 = df_pts.filter(pl.col("Ns") == 64).sort("Temp") if "Ns" in df_pts.columns else pl.DataFrame()

    apply_minimalist_style()
    save_dirs = [out_dir, lcp_dir]

    # 1. 绘制综合对比图 (chiral_susceptibility.png)
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)
    if df_48.height > 0 and df_64.height > 0:
        ax.errorbar(
            df_48["Temp"].to_numpy().astype(float),
            df_48["Mean_scaled_gev2_renorm"].to_numpy().astype(float),
            yerr=df_48["Error_scaled_gev2_renorm"].to_numpy().astype(float),
            fmt="o",
            color="#1e293b",
            mfc="white",
            mec="#1e293b",
            mew=1.4,
            elinewidth=1.4,
            capsize=3.5,
            markersize=6.0,
            label=r"$48^3 \times 16$",
            zorder=3,
        )
        ax.errorbar(
            df_64["Temp"].to_numpy().astype(float),
            df_64["Mean_scaled_gev2_renorm"].to_numpy().astype(float),
            yerr=df_64["Error_scaled_gev2_renorm"].to_numpy().astype(float),
            fmt="s",
            color="#dc2626",
            mfc="#dc2626",
            mec="#dc2626",
            mew=1.4,
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"$64^3 \times 16$",
            zorder=4,
        )
    elif df_48.height > 0:
        ax.errorbar(
            df_48["Temp"].to_numpy().astype(float),
            df_48["Mean_scaled_gev2_renorm"].to_numpy().astype(float),
            yerr=df_48["Error_scaled_gev2_renorm"].to_numpy().astype(float),
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"Data ($\chi_{\mathrm{disc}} / Z_m^2$)",
            zorder=3,
        )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax.set_title("Chiral Susceptibility", pad=12)
    ax.set_xlim(133.0, 247.0)
    ax.set_ylim(-0.02, 1.15)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")
    plt.tight_layout()
    _save_figure(fig, "chiral_susceptibility", save_dirs)

    # 2. 独立绘制 64 尺寸正常图 (chiral_susceptibility_L64.png)
    if df_64.height > 0:
        fig_64, ax_64 = plt.subplots(figsize=(8.0, 5.5), dpi=300)
        ax_64.errorbar(
            df_64["Temp"].to_numpy().astype(float),
            df_64["Mean_scaled_gev2_renorm"].to_numpy().astype(float),
            yerr=df_64["Error_scaled_gev2_renorm"].to_numpy().astype(float),
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"Data ($64^3 \times 16$)",
            zorder=3,
        )
        ax_64.set_xlabel(r"Temperature $T$ [MeV]")
        ax_64.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
        ax_64.set_title(r"Chiral Susceptibility ($64^3 \times 16$)", pad=12)
        ax_64.set_xlim(133.0, 247.0)
        ax_64.set_ylim(-0.02, 1.15)
        ax_64.grid(True, linestyle="--", alpha=0.4)
        ax_64.xaxis.set_minor_locator(AutoMinorLocator())
        ax_64.yaxis.set_minor_locator(AutoMinorLocator())
        ax_64.legend(loc="upper right", frameon=True, edgecolor="#cccccc")
        plt.tight_layout()
        _save_figure(fig_64, "chiral_susceptibility_L64", save_dirs)

    # 3. 独立绘制 48 尺寸正常图 (chiral_susceptibility_L48.png)
    if df_48.height > 0:
        fig_48, ax_48 = plt.subplots(figsize=(8.0, 5.5), dpi=300)
        ax_48.errorbar(
            df_48["Temp"].to_numpy().astype(float),
            df_48["Mean_scaled_gev2_renorm"].to_numpy().astype(float),
            yerr=df_48["Error_scaled_gev2_renorm"].to_numpy().astype(float),
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"Data ($48^3 \times 16$)",
            zorder=3,
        )
        ax_48.set_xlabel(r"Temperature $T$ [MeV]")
        ax_48.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
        ax_48.set_title(r"Chiral Susceptibility ($48^3 \times 16$)", pad=12)
        ax_48.set_xlim(133.0, 247.0)
        ax_48.set_ylim(-0.02, 1.15)
        ax_48.grid(True, linestyle="--", alpha=0.4)
        ax_48.xaxis.set_minor_locator(AutoMinorLocator())
        ax_48.yaxis.set_minor_locator(AutoMinorLocator())
        ax_48.legend(loc="upper right", frameon=True, edgecolor="#cccccc")
        plt.tight_layout()
        _save_figure(fig_48, "chiral_susceptibility_L48", save_dirs)

    print(f"[OK] 成功生成极简手征易感度图表 (含对比图与各尺寸分图) 至: {out_dir}")
    return out_dir / "chiral_susceptibility.png"


def _read_condensate_file(file_path: Path) -> pl.DataFrame:
    try:
        df = pl.read_csv(file_path, separator="\t")
        if "beta" not in df.columns:
            df = pl.read_csv(file_path, separator=" ")
    except Exception:
        df = pl.read_csv(file_path, separator=" ")

    return df.unique(subset=["beta"]).with_columns(
        pl.col("beta").cast(pl.String).alias("beta_str")
    ).with_columns(
        pl.col("beta_str").map_elements(
            lambda b: TEMP_MAP.get(str(b), TEMP_MAP.get(f"{float(b):.2f}", np.nan)),
            return_dtype=pl.Float64,
        ).alias("Temp")
    ).sort("Temp")


def plot_chiral_condensate(
    condensate_dir: Path = OUTPUT_CONDENSATE_DIR,
    out_dir: Path = FIGURES_DIR,
) -> Path:
    """以相同极简风格绘制手征凝聚温度演化图 (支持 48 尺寸、64 尺寸及对比图)"""
    file_48 = condensate_dir / "results_rm_beta_L48.txt"
    file_64 = condensate_dir / "results_rm_beta_L64.txt"
    file_default = condensate_dir / "results_rm_beta.txt"

    df_48 = _read_condensate_file(file_48 if file_48.exists() else file_default)
    df_64 = _read_condensate_file(file_64) if file_64.exists() else pl.DataFrame()

    apply_minimalist_style()
    save_dirs = [out_dir, condensate_dir]

    # 1. 综合对比图 (chiral_condensate.png)
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)
    if df_48.height > 0 and df_64.height > 0:
        ax.errorbar(
            df_48["Temp"].to_numpy().astype(float),
            df_48["pbp_rm"].to_numpy().astype(float),
            yerr=df_48["pbp_rm_err"].to_numpy().astype(float),
            fmt="o",
            color="#1e293b",
            mfc="white",
            mec="#1e293b",
            mew=1.4,
            elinewidth=1.4,
            capsize=3.5,
            markersize=6.0,
            label=r"$48^3 \times 16$",
            zorder=3,
        )
        ax.errorbar(
            df_64["Temp"].to_numpy().astype(float),
            df_64["pbp_rm"].to_numpy().astype(float),
            yerr=df_64["pbp_rm_err"].to_numpy().astype(float),
            fmt="s",
            color="#dc2626",
            mfc="#dc2626",
            mec="#dc2626",
            mew=1.4,
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"$64^3 \times 16$",
            zorder=4,
        )
    elif df_48.height > 0:
        ax.errorbar(
            df_48["Temp"].to_numpy().astype(float),
            df_48["pbp_rm"].to_numpy().astype(float),
            yerr=df_48["pbp_rm_err"].to_numpy().astype(float),
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"Data ($\langle \bar{\psi}\psi \rangle_{\mathrm{sub}}$)",
            zorder=3,
        )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\mathrm{sub}}\ [\mathrm{lattice\ unit}]$")
    ax.set_title("Chiral Condensate", pad=12)
    ax.set_xlim(133.0, 247.0)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")
    plt.tight_layout()
    _save_figure(fig, "chiral_condensate", save_dirs)

    # 2. 独立绘制 64 尺寸手征凝聚图 (chiral_condensate_L64.png)
    if df_64.height > 0:
        fig_64, ax_64 = plt.subplots(figsize=(8.0, 5.5), dpi=300)
        ax_64.errorbar(
            df_64["Temp"].to_numpy().astype(float),
            df_64["pbp_rm"].to_numpy().astype(float),
            yerr=df_64["pbp_rm_err"].to_numpy().astype(float),
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"Data ($64^3 \times 16$)",
            zorder=3,
        )
        ax_64.set_xlabel(r"Temperature $T$ [MeV]")
        ax_64.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\mathrm{sub}}\ [\mathrm{lattice\ unit}]$")
        ax_64.set_title(r"Chiral Condensate ($64^3 \times 16$)", pad=12)
        ax_64.set_xlim(133.0, 247.0)
        ax_64.grid(True, linestyle="--", alpha=0.4)
        ax_64.xaxis.set_minor_locator(AutoMinorLocator())
        ax_64.yaxis.set_minor_locator(AutoMinorLocator())
        ax_64.legend(loc="upper right", frameon=True, edgecolor="#cccccc")
        plt.tight_layout()
        _save_figure(fig_64, "chiral_condensate_L64", save_dirs)

    # 3. 独立绘制 48 尺寸手征凝聚图 (chiral_condensate_L48.png)
    if df_48.height > 0:
        fig_48, ax_48 = plt.subplots(figsize=(8.0, 5.5), dpi=300)
        ax_48.errorbar(
            df_48["Temp"].to_numpy().astype(float),
            df_48["pbp_rm"].to_numpy().astype(float),
            yerr=df_48["pbp_rm_err"].to_numpy().astype(float),
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.4,
            capsize=3.5,
            markersize=5.5,
            label=r"Data ($48^3 \times 16$)",
            zorder=3,
        )
        ax_48.set_xlabel(r"Temperature $T$ [MeV]")
        ax_48.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\mathrm{sub}}\ [\mathrm{lattice\ unit}]$")
        ax_48.set_title(r"Chiral Condensate ($48^3 \times 16$)", pad=12)
        ax_48.set_xlim(133.0, 247.0)
        ax_48.grid(True, linestyle="--", alpha=0.4)
        ax_48.xaxis.set_minor_locator(AutoMinorLocator())
        ax_48.yaxis.set_minor_locator(AutoMinorLocator())
        ax_48.legend(loc="upper right", frameon=True, edgecolor="#cccccc")
        plt.tight_layout()
        _save_figure(fig_48, "chiral_condensate_L48", save_dirs)

    print(f"[OK] 成功生成极简手征凝聚图表 (含对比图与各尺寸分图) 至: {out_dir}")
    return out_dir / "chiral_condensate.png"


def main() -> None:
    parser = argparse.ArgumentParser(description="绘制极简统一风格手征热力学物理曲线")
    parser.add_argument("--lcp-dir", type=Path, default=OUTPUT_LCP_DIR, help="LCP 数据目录 (默认: output/LCP)")
    parser.add_argument("--condensate-dir", type=Path, default=OUTPUT_CONDENSATE_DIR, help="手征凝聚数据目录 (默认: output/condensate)")
    parser.add_argument("--out-dir", type=Path, default=FIGURES_DIR, help="图表输出目录 (默认: docs/figures)")
    args = parser.parse_args()

    plot_chiral_susceptibility(lcp_dir=args.lcp_dir, out_dir=args.out_dir)
    plot_chiral_condensate(condensate_dir=args.condensate_dir, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
