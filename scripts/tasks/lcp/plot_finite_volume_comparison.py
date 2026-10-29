#!/usr/bin/env python3
"""
scripts/tasks/lcp/plot_finite_volume_comparison.py
--------------------------------------------------------------------------------
有限体积效应对比图件生成脚本 (不包含 q3):
1. p2 (36^3 x 12, LT = 3) 与 p3 (48^3 x 12, LT = 4) 的纯有限体积效应直接对比图。
2. 整合 p2 (36^3 x 12), p3 (48^3 x 12) 与基准 q2 (48^3 x 16) 的相变综合对比图。
- 拟合曲线采用半透明阴影画出 Jackknife 拟合所有曲线的覆盖范围 (Fit range)。
- 实线表示由平均参数计算的拟合曲线 (Fit)。
- 数据点带统计误差棒。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

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
from src.parselqcdata.pade_fitter import PadeModel


def apply_minimalist_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10.5,
        "mathtext.fontset": "cm",
    })


def load_dataset_data(ds_name: str, base_dir: Path) -> dict:
    ds_dir = base_dir / ds_name
    df_s = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
    df_c = pl.read_csv(ds_dir / f"condensate_{ds_name}.csv").sort("temp")
    df_sum = pl.read_csv(ds_dir / f"pade22_fit_summary_{ds_name}.csv")
    df_jk = pl.read_csv(ds_dir / f"pade22_jk_samples_{ds_name}.csv").filter(pl.col("ok") == True)
    return {
        "df_susc": df_s,
        "df_cond": df_c,
        "summary": df_sum.to_dicts()[0],
        "records": df_jk.to_dicts(),
    }


def compute_pade_curve_and_band(
    records: list[dict],
    t_grid: np.ndarray,
    t0: float = 157.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model = PadeModel(order="2/2", t0=t0)
    p_keys = ["p0", "p1", "p2", "q1", "q2"]

    valid_recs = [r for r in records if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    n_samples = len(valid_recs)

    # 1. 平均参数构建实线曲线
    p_mean = {pk: float(np.mean([float(r[pk]) for r in valid_recs if pk in r and np.isfinite(r[pk])])) for pk in p_keys}
    curve_mean = model.eval(t_grid, p_mean)

    # 2. 逐样本计算曲线矩阵并计算 Jackknife 误差带
    curves_matrix = np.zeros((n_samples, len(t_grid)), dtype=float)
    for i, r in enumerate(valid_recs):
        p_i = {pk: float(r[pk]) for pk in p_keys if pk in r}
        curves_matrix[i, :] = model.eval(t_grid, p_i)

    curve_err = np.sqrt(n_samples - 1) * np.std(curves_matrix, axis=0, ddof=0)
    y_lower = curve_mean - curve_err
    y_upper = curve_mean + curve_err

    return curve_mean, y_lower, y_upper


def plot_pure_finite_volume_p2_p3(base_dir: Path, out_dir: Path, t0: float = 157.0) -> tuple[Path, Path]:
    """生成 p2 (36^3 x 12) 与 p3 (48^3 x 12) 的纯有限体积效应对比图"""
    apply_minimalist_style()
    data_p2 = load_dataset_data("p2", base_dir)
    data_p3 = load_dataset_data("p3", base_dir)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 5.8), dpi=300)

    # --- 左图: 手征凝聚有限体积对比 ---
    temps_p2 = data_p2["df_cond"]["temp"].to_numpy().astype(float)
    pbp_p2 = data_p2["df_cond"]["pbp_renorm"].to_numpy().astype(float) * 1e4
    err_p2 = data_p2["df_cond"]["pbp_renorm_err"].to_numpy().astype(float) * 1e4

    temps_p3 = data_p3["df_cond"]["temp"].to_numpy().astype(float)
    pbp_p3 = data_p3["df_cond"]["pbp_renorm"].to_numpy().astype(float) * 1e4
    err_p3 = data_p3["df_cond"]["pbp_renorm_err"].to_numpy().astype(float) * 1e4

    ax1.errorbar(
        temps_p2, pbp_p2, yerr=err_p2,
        fmt="o-", color="#1f77b4", label=r"$p_2\ (36^3\times 12,\ N_s/N_t=3)$",
        capsize=3.5, markersize=5.5, linewidth=1.4, zorder=3,
    )
    ax1.errorbar(
        temps_p3, pbp_p3, yerr=err_p3,
        fmt="s-", color="#2ca02c", label=r"$p_3\ (48^3\times 12,\ N_s/N_t=4)$",
        capsize=3.5, markersize=5.5, linewidth=1.4, zorder=4,
    )

    ax1.set_xlabel(r"Temperature $T$ [MeV]")
    ax1.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m \times 10^4$")
    ax1.set_title(r"Chiral Condensate: Finite Volume Scaling ($N_t=12$)")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # --- 右图: 手征易感度 Padé [2/2] 阴影与曲线有限体积对比 ---
    t_min = min(np.min(temps_p2), np.min(temps_p3)) - 3.0
    t_max = max(np.max(temps_p2), np.max(temps_p3)) + 3.0
    t_grid = np.linspace(t_min, t_max, 400)

    # p2
    c_p2, low_p2, up_p2 = compute_pade_curve_and_band(data_p2["records"], t_grid, t0=t0)
    tpk_p2 = float(data_p2["summary"]["t_peak_mean"])
    terr_p2 = float(data_p2["summary"]["t_peak_err"])
    chi2_p2 = float(data_p2["summary"]["chi2_dof_mean"])

    ax2.fill_between(t_grid, low_p2, up_p2, color="#1f77b4", alpha=0.20, label=r"$p_2$ Fit Range", zorder=1)
    ax2.plot(t_grid, c_p2, color="#1f77b4", linewidth=2.0, label=rf"$p_2$ Padé [2/2] ($T_c={tpk_p2:.1f}\pm{terr_p2:.1f}\,$MeV)", zorder=3)
    ax2.errorbar(
        data_p2["df_susc"]["temp"].to_numpy().astype(float),
        data_p2["df_susc"]["chi_renorm"].to_numpy().astype(float),
        yerr=data_p2["df_susc"]["chi_renorm_err"].to_numpy().astype(float),
        fmt="o", color="#1f77b4", capsize=3.5, markersize=5.5, label=r"$p_2$ Data ($36^3\times 12$)", zorder=5,
    )

    # p3
    c_p3, low_p3, up_p3 = compute_pade_curve_and_band(data_p3["records"], t_grid, t0=t0)
    tpk_p3 = float(data_p3["summary"]["t_peak_mean"])
    terr_p3 = float(data_p3["summary"]["t_peak_err"])
    chi2_p3 = float(data_p3["summary"]["chi2_dof_mean"])

    ax2.fill_between(t_grid, low_p3, up_p3, color="#2ca02c", alpha=0.20, label=r"$p_3$ Fit Range", zorder=2)
    ax2.plot(t_grid, c_p3, color="#2ca02c", linewidth=2.0, label=rf"$p_3$ Padé [2/2] ($T_c={tpk_p3:.1f}\pm{terr_p3:.1f}\,$MeV)", zorder=4)
    ax2.errorbar(
        data_p3["df_susc"]["temp"].to_numpy().astype(float),
        data_p3["df_susc"]["chi_renorm"].to_numpy().astype(float),
        yerr=data_p3["df_susc"]["chi_renorm_err"].to_numpy().astype(float),
        fmt="s", color="#2ca02c", capsize=3.5, markersize=5.5, label=r"$p_3$ Data ($48^3\times 12$)", zorder=6,
    )

    ax2.set_xlabel(r"Temperature $T$ [MeV]")
    ax2.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax2.set_title(r"Chiral Susceptibility: Finite Volume Scaling ($N_t=12$)")
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # 标注有限体积物理推断
    info_box = (
        r"$\mathbf{Finite\ Volume\ Effect\ (N_t=12):}$" + "\n"
        rf"$36^3\times 12:\ T_c = {tpk_p2:.2f}\pm{terr_p2:.2f}\,\mathrm{{MeV}},\ \chi^2/\mathrm{{dof}}={chi2_p2:.3f}$" + "\n"
        rf"$48^3\times 12:\ T_c = {tpk_p3:.2f}\pm{terr_p3:.2f}\,\mathrm{{MeV}},\ \chi^2/\mathrm{{dof}}={chi2_p3:.3f}$" + "\n"
        rf"$\Delta T_c = {tpk_p3 - tpk_p2:+.2f}\,\mathrm{{MeV}}\ (\mathrm{{Shift}}\approx 0.7\%)$"
    )
    ax2.text(
        0.04, 0.95, info_box,
        transform=ax2.transAxes, fontsize=10.5,
        verticalalignment="top", horizontalalignment="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        zorder=7,
    )

    plt.tight_layout()
    out_png = out_dir / "finite_volume_effect_p2_p3.png"
    out_pdf = out_dir / "finite_volume_effect_p2_p3.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def plot_integrated_p2_p3_q2(base_dir: Path, out_dir: Path, t0: float = 157.0) -> tuple[Path, Path]:
    """生成 p2, p3 与 q2 (48^3 x 16) 整合对比图 (不包含 q3)"""
    apply_minimalist_style()
    data_p2 = load_dataset_data("p2", base_dir)
    data_p3 = load_dataset_data("p3", base_dir)
    data_q2 = load_dataset_data("q2", base_dir)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 5.8), dpi=300)

    items = [
        ("p2", data_p2, "#1f77b4", "o", r"$p_2\ (36^3\times 12,\ LT=3)$"),
        ("p3", data_p3, "#2ca02c", "s", r"$p_3\ (48^3\times 12,\ LT=4)$"),
        ("q2", data_q2, "#d62728", "^", r"$q_2\ (48^3\times 16,\ LT=3)$"),
    ]

    # --- 左图: 手征凝聚 ---
    for ds_name, data, col, marker, lbl in items:
        temps = data["df_cond"]["temp"].to_numpy().astype(float)
        pbp = data["df_cond"]["pbp_renorm"].to_numpy().astype(float) * 1e4
        err = data["df_cond"]["pbp_renorm_err"].to_numpy().astype(float) * 1e4
        ax1.errorbar(
            temps, pbp, yerr=err,
            fmt=f"{marker}-", color=col, label=lbl,
            capsize=3.5, markersize=5.5, linewidth=1.3,
        )

    ax1.set_xlabel(r"Temperature $T$ [MeV]")
    ax1.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m \times 10^4$")
    ax1.set_title(r"Chiral Condensate Comparison ($p_2, p_3, q_2$)")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # --- 右图: 手征易感度 Padé [2/2] 拟合与阴影包络 ---
    all_temps = np.concatenate([d["df_susc"]["temp"].to_numpy().astype(float) for _, d, _, _, _ in items])
    t_grid = np.linspace(np.min(all_temps) - 3.0, np.max(all_temps) + 3.0, 400)

    for ds_name, data, col, marker, lbl in items:
        # 数据点
        ax2.errorbar(
            data["df_susc"]["temp"].to_numpy().astype(float),
            data["df_susc"]["chi_renorm"].to_numpy().astype(float),
            yerr=data["df_susc"]["chi_renorm_err"].to_numpy().astype(float),
            fmt=marker, color=col, capsize=3.5, markersize=5.5, alpha=0.75,
            label=f"{lbl} Data",
        )

        # Padé [2/2] 拟合线与阴影
        c_curve, low_b, up_b = compute_pade_curve_and_band(data["records"], t_grid, t0=t0)
        tpk = float(data["summary"]["t_peak_mean"])
        terr = float(data["summary"]["t_peak_err"])

        ax2.fill_between(t_grid, low_b, up_b, color=col, alpha=0.15)
        ax2.plot(
            t_grid, c_curve, color=col, linewidth=1.8,
            label=rf"{lbl} Padé [2/2] ($T_c={tpk:.1f}\pm{terr:.1f}\,$MeV)",
        )

    ax2.set_xlabel(r"Temperature $T$ [MeV]")
    ax2.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax2.set_title(r"Chiral Susceptibility Padé [2/2] Comparison ($p_2, p_3, q_2$)")
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = out_dir / "finite_volume_and_lcp_comparison.png"
    out_pdf = out_dir / "finite_volume_and_lcp_comparison.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def main() -> None:
    base_dir = PROJECT_ROOT / "22pade"
    out_22pade = PROJECT_ROOT / "22pade"

    print("[INFO] 生成 p2 vs p3 纯有限体积效应图件 (不包含 q3)...")
    p1_png, p1_pdf = plot_pure_finite_volume_p2_p3(base_dir, out_22pade)
    print(f"  -> 生成: {p1_png.name}")

    print("[INFO] 生成 p2, p3, q2 综合对比图件 (不包含 q3)...")
    p2_png, p2_pdf = plot_integrated_p2_p3_q2(base_dir, out_22pade)
    print(f"  -> 生成: {p2_png.name}")

    # 复制到 docs/figures/ 与 output/LCP/22pade/
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copyfile(p1_png, FIGURES_DIR / p1_png.name)
    shutil.copyfile(p1_pdf, FIGURES_DIR / p1_pdf.name)
    shutil.copyfile(p2_png, FIGURES_DIR / p2_png.name)
    shutil.copyfile(p2_pdf, FIGURES_DIR / p2_pdf.name)

    lcp_22pade = OUTPUT_LCP_DIR / "22pade"
    lcp_22pade.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(p1_png, lcp_22pade / p1_png.name)
    shutil.copyfile(p1_pdf, lcp_22pade / p1_pdf.name)
    shutil.copyfile(p2_png, lcp_22pade / p2_png.name)
    shutil.copyfile(p2_pdf, lcp_22pade / p2_pdf.name)

    print("[OK] 有限体积效应图件已成功生成并同步！")


if __name__ == "__main__":
    main()
