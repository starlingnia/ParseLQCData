#!/usr/bin/env python3
"""
scripts/tasks/lcp/plot_lcp_figures.py
--------------------------------------------------------------------------------
常物理线 (LCP) 出版级图表绘制脚本 (独立专注绘图):
1. Figure 1: 4 系综手征凝聚与易感度综合对比图 (含 Padé 拟合带)
2. Figure 2: 4 系综裸光/奇夸克手征凝聚图
3. Figure 3: 各系综多阶 Padé 拟合对比图 (含 Jackknife 误差阴影带)
4. Figure 4: 物理量纲下的手征凝聚总图 [GeV^3]
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import shutil
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.parselqcdata.pade_fitter import PadeModel

BASE_DIR = PROJECT_ROOT / "22pade"
DOCS_DIR = PROJECT_ROOT / "docs"
FIGURES_DIR = DOCS_DIR / "figures"

# 统一风格配置
STYLES = {
    "p2": {"color": "#1f77b4", "marker": "o", "label": r"$p_2\ (36^3\times 12)$"},
    "p3": {"color": "#2ca02c", "marker": "s", "label": r"$p_3\ (48^3\times 12)$"},
    "q2": {"color": "#d62728", "marker": "^", "label": r"$q_2\ (48^3\times 16)$"},
    "q3": {"color": "#9467bd", "marker": "D", "label": r"$q_3\ (64^3\times 16)$"},
}

NS_NT_LABELS = {
    "p2": r"36^3\times 12",
    "p3": r"48^3\times 12",
    "q2": r"48^3\times 16",
    "q3": r"64^3\times 16",
}


def apply_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 9.5,
        "mathtext.fontset": "cm",
    })


def load_order_records(ds_name: str, order: str) -> tuple[dict, list[dict], np.ndarray, np.ndarray, np.ndarray]:
    """读取指定系综与阶数的拟合结果与数据"""
    ds_dir = BASE_DIR / ds_name
    df_s = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
    temps = df_s["temp"].to_numpy().astype(float)
    errs = df_s["chi_renorm_err"].to_numpy().astype(float)
    chi_vals = df_s["chi_renorm"].to_numpy().astype(float)

    ord_tag = order.replace("/", "")
    sum_df = pl.read_csv(ds_dir / f"pade{ord_tag}_fit_summary_{ds_name}.csv")
    jk_df = pl.read_csv(ds_dir / f"pade{ord_tag}_jk_samples_{ds_name}.csv").filter(pl.col("ok") == True)
    return sum_df.to_dicts()[0], jk_df.to_dicts(), temps, chi_vals, errs


def compute_pade_band(order: str, records: list[dict], t_grid: np.ndarray, t0: float = 157.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算 Padé 拟合均值曲线与 Jackknife 误差带"""
    model = PadeModel(order=order, t0=t0)
    p_keys = model.param_keys
    valid = [r for r in records if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    n_samples = len(valid)

    p_mean = {pk: float(np.mean([float(r[pk]) for r in valid if pk in r and np.isfinite(r[pk])])) for pk in p_keys}
    curve_mean = model.eval(t_grid, p_mean)

    mat = np.zeros((n_samples, len(t_grid)), dtype=float)
    for i, r in enumerate(valid):
        p_i = {pk: float(r[pk]) for pk in p_keys if pk in r}
        mat[i, :] = model.eval(t_grid, p_i)

    err = np.sqrt(n_samples - 1) * np.std(mat, axis=0, ddof=0)
    return curve_mean, curve_mean - err, curve_mean + err


# ==============================================================================
# 1. Figure 1: 4 系综手征凝聚与易感度综合对比图
# ==============================================================================
def plot_figure1_combined() -> tuple[Path, Path]:
    apply_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.0), dpi=300)

    # 左图: 格点单位手征凝聚
    max_c = 0.0
    for ds, st in STYLES.items():
        df_c = pl.read_csv(BASE_DIR / ds / f"condensate_{ds}.csv").sort("temp")
        t = df_c["temp"].to_numpy().astype(float)
        pbp = df_c["pbp_renorm"].to_numpy().astype(float)
        pbp_err = df_c["pbp_renorm_err"].to_numpy().astype(float)
        max_c = max(max_c, np.max(pbp + pbp_err))

        ax1.errorbar(
            t, pbp, yerr=pbp_err,
            fmt=st["marker"], linestyle="none", color=st["color"], label=st["label"],
            capsize=3.5, markersize=5.5, alpha=0.9,
        )

    ax1.set_xlabel(r"Temperature $T$ [MeV]")
    ax1.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\mathrm{sub}} / Z_m$ (lattice units)")
    ax1.set_title(r"Renormalized Chiral Condensate ($p_2, p_3, q_2, q_3$)")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.4f}"))
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.set_ylim(-0.0001, max_c * 1.30)
    ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # 右图: 易感度拟合带与数据点
    orders_map = {"p2": "2/2", "p3": "0/2", "q2": "2/2", "q3": "0/2"}
    all_temps = []
    max_s = 0.0

    for ds, ord_val in orders_map.items():
        st = STYLES[ds]
        _, _, temps, chi, errs = load_order_records(ds, ord_val)
        all_temps.extend(temps)
        ax2.errorbar(
            temps, chi, yerr=errs,
            fmt=st["marker"], color=st["color"], capsize=3.5, markersize=5.5, alpha=0.75,
            label=f"{st['label']} Data",
        )
        max_s = max(max_s, np.max(chi + errs))

    t_grid = np.linspace(min(all_temps) - 3.0, max(all_temps) + 3.0, 500)
    for ds, ord_val in orders_map.items():
        st = STYLES[ds]
        summary, records, _, _, _ = load_order_records(ds, ord_val)
        c_curve, low_b, up_b = compute_pade_band(ord_val, records, t_grid)
        tpk = summary["t_peak_mean"]
        terr = summary["t_peak_err"]

        ax2.fill_between(t_grid, low_b, up_b, color=st["color"], alpha=0.15)
        ax2.plot(
            t_grid, c_curve, color=st["color"], linewidth=1.8,
            label=rf"{st['label']} [{ord_val}] ($T_c={tpk:.1f}\pm{terr:.1f}\,$MeV)",
        )
        max_s = max(max_s, np.max(up_b[np.isfinite(up_b)]))

    ax2.set_xlabel(r"Temperature $T$ [MeV]")
    ax2.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax2.set_title(r"Chiral Susceptibility Comparison ($p_2[2/2], p_3[0/2], q_2[2/2], q_3[0/2]$)")
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f}"))
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.set_ylim(-0.05, max_s * 1.35)
    ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_pdf = FIGURES_DIR / "chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02.pdf"
    out_png = FIGURES_DIR / "chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    shutil.copyfile(out_pdf, BASE_DIR / out_pdf.name)
    shutil.copyfile(out_png, BASE_DIR / out_png.name)
    return out_png, out_pdf


# ==============================================================================
# 2. Figure 2: 4 系综裸 light 与 strange 凝聚数据散点图
# ==============================================================================
def plot_figure2_bare_condensate() -> tuple[Path, Path]:
    apply_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 5.8), dpi=300)
    max_l = 0.0
    max_s = 0.0

    for ds, st in STYLES.items():
        df = pl.read_csv(BASE_DIR / ds / f"condensate_{ds}.csv").sort("temp")
        t = df["temp"].to_numpy().astype(float)
        l_val = df["pbpl"].to_numpy().astype(float)
        l_err = df["pbpl_err"].to_numpy().astype(float)
        s_val = df["pbps"].to_numpy().astype(float)
        s_err = df["pbps_err"].to_numpy().astype(float)

        max_l = max(max_l, np.max(l_val + l_err))
        max_s = max(max_s, np.max(s_val + s_err))

        ax1.errorbar(
            t, l_val, yerr=l_err,
            fmt=st["marker"], linestyle="none", color=st["color"], label=st["label"],
            capsize=3.5, markersize=5.5, alpha=0.9,
        )
        ax2.errorbar(
            t, s_val, yerr=s_err,
            fmt=st["marker"], linestyle="none", color=st["color"], label=st["label"],
            capsize=3.5, markersize=5.5, alpha=0.9,
        )

    ax1.set_xlabel(r"Temperature $T$ [MeV]")
    ax1.set_ylabel(r"$\langle\bar{\psi}\psi\rangle_l$ (lattice units)")
    ax1.set_title(r"Bare Light Quark Condensate $\langle\bar{\psi}\psi\rangle_l$")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.4f}"))
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.set_ylim(0.0, max_l * 1.30)
    ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    ax2.set_xlabel(r"Temperature $T$ [MeV]")
    ax2.set_ylabel(r"$\langle\bar{\psi}\psi\rangle_s$ (lattice units)")
    ax2.set_title(r"Bare Strange Quark Condensate $\langle\bar{\psi}\psi\rangle_s$")
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.3f}"))
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.set_ylim(0.0, max_s * 1.30)
    ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_pdf = FIGURES_DIR / "chiral_condensate_bare_light_strange.pdf"
    out_png = FIGURES_DIR / "chiral_condensate_bare_light_strange.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    shutil.copyfile(out_pdf, BASE_DIR / out_pdf.name)
    shutil.copyfile(out_png, BASE_DIR / out_png.name)
    return out_png, out_pdf


# ==============================================================================
# 3. Figure 3: 各系综多阶 Padé 拟合对比图 (含 Jackknife 误差带)
# ==============================================================================
def plot_figure3_multimodal_fits() -> None:
    apply_plot_style()
    order_specs = {
        "p2": ["1/2", "2/2", "3/2"],
        "p3": ["0/2", "1/2", "2/2"],
        "q2": ["1/2", "2/2", "3/2", "1/3"],
        "q3": ["0/2", "1/2"],
    }
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    linestyles = ["-", "--", "-.", ":"]

    for ds_name, orders in order_specs.items():
        fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=300)
        _, _, temps, chi_vals, chi_errs = load_order_records(ds_name, orders[0])
        t_grid = np.linspace(min(temps) - 4.0, max(temps) + 4.0, 500)

        ax.errorbar(
            temps, chi_vals, yerr=chi_errs,
            fmt="ko", capsize=3.5, markersize=5.5, label="Lattice Data", zorder=5,
        )

        max_y = np.max(chi_vals + chi_errs)
        for i, order in enumerate(orders):
            summary, records, _, _, _ = load_order_records(ds_name, order)
            c_mean, c_low, c_up = compute_pade_band(order, records, t_grid)
            tpk = summary["t_peak_mean"]
            terr = summary["t_peak_err"]
            chi2_dof = summary["chi2_dof_mean"]
            col = colors[i % len(colors)]
            ls = linestyles[i % len(linestyles)]

            ax.fill_between(t_grid, c_low, c_up, color=col, alpha=0.15, zorder=2)
            ax.plot(
                t_grid, c_mean, color=col, linestyle=ls, linewidth=1.8,
                label=rf"Padé [{order}] ($T_c={tpk:.1f}\pm{terr:.1f}\,$MeV, $\chi^2/\mathrm{{dof}}={chi2_dof:.2f}$)",
                zorder=3,
            )
            max_y = max(max_y, np.max(c_up[np.isfinite(c_up)]))

        ax.set_xlabel(r"Temperature $T$ [MeV]")
        ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
        ax.set_title(rf"Chiral Susceptibility Multi-Order Padé Fits: ${NS_NT_LABELS[ds_name]}$")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f}"))
        ax.grid(True, linestyle="--", alpha=0.35)
        ax.set_ylim(-0.05, max_y * 1.35)
        ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

        plt.tight_layout()
        out_pdf = FIGURES_DIR / f"chiral_susceptibility_multimodal_{ds_name}.pdf"
        out_png = FIGURES_DIR / f"chiral_susceptibility_multimodal_{ds_name}.png"
        fig.savefig(out_pdf)
        fig.savefig(out_png, dpi=300)
        plt.close(fig)

        shutil.copyfile(out_pdf, BASE_DIR / ds_name / out_pdf.name)
        shutil.copyfile(out_png, BASE_DIR / ds_name / out_png.name)


# ==============================================================================
# 4. Figure 4: 物理单位重整化手征凝聚总图 [GeV^3]
# ==============================================================================
def plot_figure4_physical_condensate() -> tuple[Path, Path]:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(8.5, 5.8), dpi=300)
    max_c = 0.0

    for ds, st in STYLES.items():
        df = pl.read_csv(BASE_DIR / ds / f"condensate_{ds}.csv").sort("temp")
        t = df["temp"].to_numpy().astype(float)
        c_lat = df["pbp_renorm"].to_numpy().astype(float)
        err_lat = df["pbp_renorm_err"].to_numpy().astype(float)
        a_inv_gev = df["a_inv"].to_numpy().astype(float) / 1000.0
        scale3 = a_inv_gev ** 3
        c_phys = c_lat * scale3
        err_phys = err_lat * scale3
        max_c = max(max_c, np.max(c_phys + err_phys))

        ax.errorbar(
            t, c_phys, yerr=err_phys,
            fmt=st["marker"], linestyle="none", color=st["color"], label=st["label"],
            capsize=3.5, markersize=5.5, alpha=0.9,
        )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\mathrm{sub}} / Z_m\ [\mathrm{GeV}^3]$")
    ax.set_title(r"Renormalized Chiral Condensate in Physical Units ($p_2, p_3, q_2, q_3$)")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.3f}"))
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_ylim(-0.002, max_c * 1.30)
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_pdf = FIGURES_DIR / "chiral_condensate_physical_all.pdf"
    out_png = FIGURES_DIR / "chiral_condensate_physical_all.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    shutil.copyfile(out_pdf, BASE_DIR / out_pdf.name)
    shutil.copyfile(out_png, BASE_DIR / out_png.name)
    return out_png, out_pdf


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    print("==========================================================================")
    print("[LCP-PLOT] 启动常物理线 (LCP) 出版级图表绘制...")
    print("==========================================================================")
    print("  -> 绘制 Figure 1: 综合手征凝聚与易感度对比图...")
    plot_figure1_combined()
    print("  -> 绘制 Figure 2: 裸 light/strange 凝聚图...")
    plot_figure2_bare_condensate()
    print("  -> 绘制 Figure 3: 各系综多阶 Padé 拟合带图...")
    plot_figure3_multimodal_fits()
    print("  -> 绘制 Figure 4: 物理单位手征凝聚总图 [GeV^3]...")
    plot_figure4_physical_condensate()
    print("\n[DONE] 所有图表绘制完成并输出至 docs/figures/ 与 22pade/！")


if __name__ == "__main__":
    main()
