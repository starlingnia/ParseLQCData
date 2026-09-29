#!/usr/bin/env python3
"""
scripts/plot_split_raw.py
--------------------------------------------------------------------------------
把 output/condensate/ 的全部 13 个系综拆成两条独立序列，并各画两个物理量的原始数据图：

  序列 A | 48^3 x 16 温度扫描 : beta = 4.13, 4.15, 4.17, 4.18, 4.20, 4.23, 4.30, 4.405
  序列 B | 固定 beta = 4.17   : L32T12 / L32T14 / L40T16 / L48T18 (stream1&2)

  观测量 1 | 温度 - 手征凝聚   Delta_{l,s} = pbp_rm  (减去残余质量并重整化)
  观测量 2 | 温度 - 手征磁化率 chi_vol = Mean_vol_scaled = chi_raw * Ns^3 * Nt

产出 (全部落在 output/condensate/split/)：
  temperature_scan_48T16__condensate.csv
  temperature_scan_48T16__susceptibility.csv
  fixed_beta4.17__condensate.csv
  fixed_beta4.17__susceptibility.csv
  condensate_vs_T_raw.png / .pdf
  susceptibility_vs_T_raw.png / .pdf

说明：两条序列互不重叠；beta=4.17 的 48^3x16 系综属于序列 A，绘图时在序列 B 面板
      中以空心星号作为"同一温度下的参考点"单独标出（它换了夸克质量 ml=0.001001）。
--------------------------------------------------------------------------------
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
COND_DIR = PROJECT_ROOT / "output" / "condensate"
SPLIT_DIR = COND_DIR / "split"
COND_FILE = COND_DIR / "all_ensembles_condensate.csv"
SUSC_FILE = COND_DIR / "all_ensembles_susceptibility.csv"

TRANSITION_REGION = (155.0, 158.0)

# 配色
C_SCAN = "#1f4e9c"        # 温度扫描序列
C_SCAN_FILL = "#7aa6e8"
C_FIXED = "#7b2d8b"       # 固定 beta=4.17 序列
C_FIXED_FILL = "#c58fd4"
C_REF = "#c1121f"         # 共享参考点

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.linewidth": 1.1,
    "axes.unicode_minus": False,
    "mathtext.fontset": "dejavusans",
})


# ----------------------------------------------------------------------------- 
# 1. 数据装载与序列拆分
# -----------------------------------------------------------------------------
def load_merged() -> pl.DataFrame:
    """把凝聚总表与磁化率总表按系综名合并为一张宽表。"""
    cond = pl.read_csv(COND_FILE)
    susc = pl.read_csv(SUSC_FILE).rename({"Ensemble": "dataset_name"})
    keep = ["dataset_name", "Mean_unscaled", "Error_unscaled",
            "Mean_vol_scaled", "Error_vol_scaled", "Mean_scaled", "Error_scaled"]
    return cond.join(susc.select(keep), on="dataset_name", how="left")


def split_series(df: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """序列 A = 48^3x16 温度扫描; 序列 B = 其余 beta=4.17 系综 (两者不重叠)。"""
    is_scan = (pl.col("ns") == 48) & (pl.col("nt") == 16)
    scan = df.filter(is_scan).sort("temperature")
    fixed = df.filter((pl.col("beta") == 4.17) & ~is_scan).sort("temperature")
    return scan, fixed


def write_condensate_csv(path: Path, df: pl.DataFrame, series: str) -> None:
    pl.DataFrame({
        "series": [series] * df.height,
        "ensemble": df["dataset_name"],
        "beta": df["beta"],
        "temperature_MeV": df["temperature"],
        "ns": df["ns"],
        "nt": df["nt"],
        "ml": df["ml"],
        "ms": df["ms"],
        "num_cfgs": df["num_cfgs"],
        "condensate_renorm": df["pbp_rm"],
        "condensate_renorm_err": df["pbp_rm_err"],
        "condensate_light_bare": df["pbpl"],
        "condensate_light_bare_err": df["pbpl_err"],
        "condensate_strange_bare": df["pbps"],
        "condensate_strange_bare_err": df["pbps_err"],
    }).write_csv(path)


def write_susceptibility_csv(path: Path, df: pl.DataFrame, series: str) -> None:
    pl.DataFrame({
        "series": [series] * df.height,
        "ensemble": df["dataset_name"],
        "beta": df["beta"],
        "temperature_MeV": df["temperature"],
        "ns": df["ns"],
        "nt": df["nt"],
        "num_cfgs": df["num_cfgs"],
        "chi_raw": df["Mean_unscaled"],
        "chi_raw_err": df["Error_unscaled"],
        "chi_vol_scaled": df["Mean_vol_scaled"],
        "chi_vol_scaled_err": df["Error_vol_scaled"],
        "chi_scaled": df["Mean_scaled"],
        "chi_scaled_err": df["Error_scaled"],
    }).write_csv(path)


# -----------------------------------------------------------------------------
# 2. 绘图
# -----------------------------------------------------------------------------
def _panel_style(ax: plt.Axes, show_band: bool, band_label: str | None) -> None:
    t_low, t_high = TRANSITION_REGION
    if show_band:
        ax.axvspan(t_low, t_high, color="#9aa0a6", alpha=0.20,
                   label=band_label, zorder=0)
    ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.65, zorder=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=10.5)


def _plot_scan_panel(ax, df, ycol, yerrcol, ylabel, beta_offsets):
    t = df["temperature"].to_numpy()
    y = df[ycol].to_numpy()
    e = df[yerrcol].to_numpy()

    ax.plot(t, y, "-", color=C_SCAN, linewidth=1.4, alpha=0.55, zorder=3)
    ax.errorbar(t, y, yerr=e, fmt="o", color=C_SCAN, ecolor=C_SCAN,
                mfc=C_SCAN_FILL, mec=C_SCAN, mew=1.3, ms=7.5,
                elinewidth=1.4, capsize=4, capthick=1.3, zorder=4,
                label=r"$48^3\times16$ temperature scan")

    betas = df["beta"].to_list()
    for i, (xi, yi) in enumerate(zip(t, y)):
        off, ha = beta_offsets[i]
        ax.annotate(rf"$\beta={betas[i]:g}$", (xi, yi),
                    textcoords="offset points", xytext=off,
                    ha=ha, va="center", fontsize=8.6, color=C_SCAN,
                    fontweight="bold")

    ax.set_xlabel("Temperature $T$ [MeV]", fontsize=12.5)
    ax.set_ylabel(ylabel, fontsize=12.5)
    _panel_style(ax, show_band=True, band_label="Pseudo-critical region (155–158 MeV)")


def _plot_fixed_panel(ax, df, ref_row, ycol, yerrcol, ylabel, tag_offsets):
    t = df["temperature"].to_numpy()
    y = df[ycol].to_numpy()
    e = df[yerrcol].to_numpy()
    ns = df["ns"].to_list()
    nt = df["nt"].to_list()
    names = df["dataset_name"].to_list()

    # 大体积点更大，便于一眼看出 Ns 趋势
    sizes = [62 if n == 48 else 46 for n in ns]

    # 误差棒单独画，避免 fmt='none' 时图例混乱
    ax.errorbar(t, y, yerr=e, fmt="none", ecolor=C_FIXED,
                elinewidth=1.4, capsize=4, capthick=1.3, zorder=3)
    ax.scatter(t, y, s=sizes, marker="s", facecolor=C_FIXED_FILL,
               edgecolor=C_FIXED, linewidths=1.4, zorder=4,
               label=r"Fixed $\beta=4.17$ ($m_l=0.0020$)")

    # 共享参考点：同温度的 48^3x16 (来自温度扫描, ml 不同)
    if ref_row is not None:
        ax.errorbar([ref_row["temperature"]], [ref_row[ycol]],
                    yerr=[ref_row[yerrcol]], fmt="*", color=C_REF,
                    ecolor=C_REF, ms=15, elinewidth=1.3, capsize=4,
                    capthick=1.3, zorder=5,
                    label=r"$48^3\times16$ ref. ($m_l=0.0010$)")

    for i, name in enumerate(names):
        tag = rf"${ns[i]}^3\times{nt[i]}$"
        off, ha = tag_offsets[i]
        if "_2" in name:
            tag += " (2)"
        ax.annotate(tag, (t[i], y[i]), textcoords="offset points", xytext=off,
                    ha=ha, va="center", fontsize=8.6, color=C_FIXED,
                    fontweight="bold")

    ax.set_xlabel("Temperature $T$ [MeV]", fontsize=12.5)
    ax.set_ylabel(ylabel, fontsize=12.5)
    _panel_style(ax, show_band=True, band_label=None)


def make_figure(scan, fixed, ref_row, ycol, yerrcol, ylabel, suptitle, subtitle,
                baseline: float | None, basename: str,
                ylim: tuple[float, float],
                beta_offsets: list[tuple[tuple[int, int], str]],
                tag_offsets: list[tuple[tuple[int, int], str]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.7), dpi=150)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.795, bottom=0.16, wspace=0.22)

    _plot_scan_panel(axes[0], scan, ycol, yerrcol, ylabel, beta_offsets)
    _plot_fixed_panel(axes[1], fixed, ref_row, ycol, yerrcol, ylabel, tag_offsets)

    for ax in axes:
        ax.set_ylim(*ylim)
        ax.margins(x=0.09)
        if baseline is not None:
            ax.axhline(baseline, color="black", linestyle="--", linewidth=0.9,
                       alpha=0.55, zorder=2)
        ax.legend(loc="upper right", frameon=True, framealpha=0.93,
                  fontsize=9.2, borderpad=0.6)

    axes[0].set_title(r"Series A — $48^3\times16$ temperature scan  ($\beta$ sweep)",
                      fontsize=12, pad=10)
    axes[1].set_title(r"Series B — fixed $\beta=4.17$  ($N_s/N_t$ variation)",
                      fontsize=12, pad=10)

    fig.suptitle(suptitle, fontsize=15.5, fontweight="bold", y=0.965)
    fig.text(0.5, 0.885, subtitle, ha="center", fontsize=10, color="#444444")
    fig.text(0.5, 0.025,
             "Raw ensemble data with Jackknife error bars — no fits applied.  "
             "Source: output/condensate/all_ensembles_{condensate,susceptibility}.csv",
             ha="center", fontsize=8.3, color="#666666")

    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(SPLIT_DIR / f"{basename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(SPLIT_DIR / f"{basename}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  [PLOT] {SPLIT_DIR / (basename + '.png')}")


# -----------------------------------------------------------------------------
def main() -> None:
    if not COND_FILE.exists() or not SUSC_FILE.exists():
        print(f"[ERROR] 缺少输入总表: {COND_FILE} / {SUSC_FILE}")
        sys.exit(1)

    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_merged()
    scan, fixed = split_series(df)

    print("=" * 78)
    print("  系综拆分结果")
    print("=" * 78)
    print(f"  序列 A | 48^3x16 温度扫描 : {scan.height} 个系综")
    for name, t in zip(scan["dataset_name"].to_list(), scan["temperature"].to_list()):
        print(f"      - {name:<42s} T = {t:7.2f} MeV")
    print(f"  序列 B | 固定 beta=4.17   : {fixed.height} 个系综")
    for name, t in zip(fixed["dataset_name"].to_list(), fixed["temperature"].to_list()):
        print(f"      - {name:<42s} T = {t:7.2f} MeV")

    write_condensate_csv(SPLIT_DIR / "temperature_scan_48T16__condensate.csv", scan,
                         "48T16_temperature_scan")
    write_condensate_csv(SPLIT_DIR / "fixed_beta4.17__condensate.csv", fixed,
                         "beta4.17_fixed")
    write_susceptibility_csv(SPLIT_DIR / "temperature_scan_48T16__susceptibility.csv", scan,
                             "48T16_temperature_scan")
    write_susceptibility_csv(SPLIT_DIR / "fixed_beta4.17__susceptibility.csv", fixed,
                             "beta4.17_fixed")
    print(f"\n  已写出 4 张拆分表 -> {SPLIT_DIR}")

    # beta=4.17 的 48^3x16 参考点（属于序列 A）
    ref = scan.filter(pl.col("beta") == 4.17)
    ref_row = None
    if ref.height:
        ref_row = {c: ref[c][0] for c in ref.columns}

    make_figure(
        scan, fixed, ref_row,
        ycol="pbp_rm", yerrcol="pbp_rm_err",
        ylabel=r"Renormalized chiral condensate $\Delta_{\ell,s}$",
        suptitle=r"Chiral Condensate vs Temperature — Raw Ensemble Data",
        subtitle=(r"$\Delta_{\ell,s}=\langle\bar\psi\psi\rangle_{\rm sub}/Z_m$ "
                  r"with residual-mass subtraction;  two series stored separately"),
        baseline=0.0,
        basename="condensate_vs_T_raw",
        ylim=(-0.00015, 0.00250),
        # 序列 A: 单调下降 -> 标签统一放在点右侧
        beta_offsets=[((12, 13), "left"), ((12, 13), "left"), ((-11, -16), "right"),
                      ((12, -14), "left"), ((13, -14), "left"), ((13, -14), "left"),
                      ((12, 13), "left"), ((12, 13), "left")],
        # 序列 B: 两个 L48T18 流重叠 -> 上下错开; 右端 32^3x12 标签改为左上
        tag_offsets=[((-9, 15), "center"), ((16, -13), "left"),
                     ((15, 5), "left"), ((15, 5), "left"), ((-14, 13), "right")],
    )

    make_figure(
        scan, fixed, ref_row,
        ycol="Mean_vol_scaled", yerrcol="Error_vol_scaled",
        ylabel=r"Chiral susceptibility $\chi_{\rm vol}=\chi\cdot N_s^3N_t$",
        suptitle=r"Chiral Susceptibility vs Temperature — Raw Ensemble Data",
        subtitle=(r"$\chi_{\rm vol}$ removes the trivial $1/V$ volume factor, "
                  r"so both series are directly comparable"),
        baseline=None,
        basename="susceptibility_vs_T_raw",
        ylim=(-0.006, 0.172),
        # 序列 A: 峰区三个点挤在一起 -> 上下左右错开
        beta_offsets=[((-12, -14), "right"), ((6, 15), "left"), ((-11, -15), "right"),
                      ((17, 11), "left"), ((15, 6), "left"), ((15, 6), "left"),
                      ((12, 11), "left"), ((12, 11), "left")],
        tag_offsets=[((-9, 16), "center"), ((17, -13), "left"),
                     ((15, 5), "left"), ((15, 5), "left"), ((-14, 13), "right")],
    )

    print("\n完成。")


if __name__ == "__main__":
    main()
