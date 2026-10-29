#!/usr/bin/env python3
"""
scripts/tasks/lcp/update_report_and_plots.py
--------------------------------------------------------------------------------
根据用户最新指令更新图件与 report.tex：
1. 修复/过滤 Padé 拟合中因实轴极点 (Pole) 或零自由度导致的病态拟合：
   - p2: 保留 [1/2], [2/2], [3/2] (剔除 0/2 与 1/3)
   - p3: 保留 [0/2], [1/2], [2/2]
   - q2: 保留 [1/2], [2/2], [3/2], [1/3] (剔除 0/2 与 0/4)
   - q3: 保留 [0/2], [1/2] (剔除 dof=0 且有实根的 2/2 及更高阶)
2. 为各系综的所有有效拟合阶数均绘制 Jackknife 误差阴影带 (Fit error range)。
3. 示意栏 (Legend) 调整位置与画幅 Headroom，坚决杜绝遮挡任何数据点与曲线峰值。
4. 去除独立的 Finite Volume Effect 示意图。
5. 去除介子质量差部分 (原 Section 3)。
6. 报告 (report.tex) 绝对不含任何 title/author/date 以及 section/subsection 标题，纯图表呈现。
7. 表格与坐标轴全面采用物理单位 (GeV^3, GeV^2, MeV) 与科学计数法，杜绝乘常数标头。
8. 编译为 docs/report.pdf。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter, ScalarFormatter
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DOCS_DIR = PROJECT_ROOT / "docs"
FIGURES_DIR = DOCS_DIR / "figures"
BASE_DIR = PROJECT_ROOT / "22pade"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_LCP_DIR
from src.parselqcdata.pade_fitter import (
    PadeModel,
    fit_jackknife_susceptibility_pade,
)


def apply_minimalist_style() -> None:
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


_JK_CACHE: dict[tuple[str, str], tuple[dict, list[dict], np.ndarray, np.ndarray, np.ndarray]] = {}

def compute_jackknife_for_order(
    ds_name: str,
    order: str,
    t0: float = 157.0,
) -> tuple[dict, list[dict], np.ndarray, np.ndarray, np.ndarray]:
    """运行指定系综与指定阶数的 Jackknife 拟合并返回汇总、样本和数据序列 (带内存与磁盘缓存)"""
    key = (ds_name, order)
    if key in _JK_CACHE:
        return _JK_CACHE[key]

    ds_dir = BASE_DIR / ds_name
    df_s = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
    temps = df_s["temp"].to_numpy().astype(float)
    errs = df_s["chi_renorm_err"].to_numpy().astype(float)
    betas = df_s["beta"].to_numpy().astype(float)
    chi_vals = df_s["chi_renorm"].to_numpy().astype(float)

    ord_tag = order.replace("/", "")
    sum_file = ds_dir / f"pade{ord_tag}_fit_summary_{ds_name}.csv"
    jk_file = ds_dir / f"pade{ord_tag}_jk_samples_{ds_name}.csv"

    if sum_file.exists() and jk_file.exists():
        df_sum = pl.read_csv(sum_file)
        df_jk_cached = pl.read_csv(jk_file).filter(pl.col("ok") == True)
        summary = df_sum.to_dicts()[0]
        records = df_jk_cached.to_dicts()
        _JK_CACHE[key] = (summary, records, temps, chi_vals, errs)
        return summary, records, temps, chi_vals, errs

    if ds_name == "q3":
        jk_csv = OUTPUT_LCP_DIR / "susceptibility_jk_samples_L64.csv"
    else:
        jk_csv = OUTPUT_LCP_DIR / ds_name / f"susceptibility_jk_samples_{ds_name}.csv"

    df_jk = pl.read_csv(jk_csv)
    samples = []
    for b in betas:
        s = df_jk.filter(pl.col("beta") == b).sort("jk_index")["chi_renorm"].to_numpy().astype(float)
        samples.append(s)
    min_l = min(len(s) for s in samples)
    chi_mat = np.array([s[:min_l] for s in samples], dtype=float)

    summary, records = fit_jackknife_susceptibility_pade(temps, chi_mat, errs, order=order, t0=t0)
    pl.DataFrame([summary]).write_csv(sum_file)
    pl.DataFrame(records).write_csv(jk_file)

    _JK_CACHE[key] = (summary, records, temps, chi_vals, errs)
    return summary, records, temps, chi_vals, errs


def compute_band_from_records(
    order: str,
    records: list[dict],
    t_grid: np.ndarray,
    t0: float = 157.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算平均拟合曲线与 Jackknife 误差带"""
    model = PadeModel(order=order, t0=t0)
    p_keys = model.param_keys
    valid_recs = [r for r in records if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    n_samples = len(valid_recs)

    p_mean = {pk: float(np.mean([float(r[pk]) for r in valid_recs if pk in r and np.isfinite(r[pk])])) for pk in p_keys}
    curve_mean = model.eval(t_grid, p_mean)

    curves_mat = np.zeros((n_samples, len(t_grid)), dtype=float)
    for i, r in enumerate(valid_recs):
        p_i = {pk: float(r[pk]) for pk in p_keys if pk in r}
        curves_mat[i, :] = model.eval(t_grid, p_i)

    curve_err = np.sqrt(n_samples - 1) * np.std(curves_mat, axis=0, ddof=0)
    return curve_mean, curve_mean - curve_err, curve_mean + curve_err


def generate_multimodal_plots() -> dict[str, list[dict]]:
    """生成每个系综的多阶 Padé 拟合图，包含每个阶数的误差阴影带并排除病态阶数"""
    apply_minimalist_style()
    order_specs = {
        "p2": ["1/2", "2/2", "3/2"],
        "p3": ["0/2", "1/2", "2/2"],
        "q2": ["1/2", "2/2", "3/2", "1/3"],
        "q3": ["0/2", "1/2"],
    }
    ns_nt_labels = {
        "p2": r"p_2\ (36^3 \times 12)",
        "p3": r"p_3\ (48^3 \times 12)",
        "q2": r"q_2\ (48^3 \times 16)",
        "q3": r"q_3\ (64^3 \times 16)",
    }
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    linestyles = ["-", "--", "-.", ":"]

    fit_stats = {}

    for ds_name, orders in order_specs.items():
        print(f"[INFO] 正在为 {ds_name} 生成包含误差阴影的多阶对比图...")
        fig, ax = plt.subplots(figsize=(8.5, 5.8), dpi=300)

        # 加载数据点
        ds_dir = BASE_DIR / ds_name
        df_s = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
        temps = df_s["temp"].to_numpy().astype(float)
        chi_vals = df_s["chi_renorm"].to_numpy().astype(float)
        chi_errs = df_s["chi_renorm_err"].to_numpy().astype(float)

        t_grid = np.linspace(np.min(temps) - 4.0, np.max(temps) + 4.0, 500)

        ax.errorbar(
            temps, chi_vals, yerr=chi_errs,
            fmt="ko", capsize=4.0, markersize=5.5, label="Lattice Data", zorder=5,
        )

        ds_stats = []
        max_y_seen = np.max(chi_vals + chi_errs)

        for i, order in enumerate(orders):
            summary, records, _, _, _ = compute_jackknife_for_order(ds_name, order)
            c_mean, c_low, c_up = compute_band_from_records(order, records, t_grid)

            tpk = summary["t_peak_mean"]
            terr = summary["t_peak_err"]
            chi2_dof = summary["chi2_dof_mean"]
            col = colors[i % len(colors)]
            ls = linestyles[i % len(linestyles)]

            # 绘制阴影带
            ax.fill_between(t_grid, c_low, c_up, color=col, alpha=0.15, zorder=2)
            # 绘制均值拟合曲线
            ax.plot(
                t_grid, c_mean, color=col, linestyle=ls, linewidth=1.8,
                label=rf"Padé [{order}] ($T_c={tpk:.1f}\pm{terr:.1f}\,$MeV, $\chi^2/\mathrm{{dof}}={chi2_dof:.2f}$)",
                zorder=3,
            )

            max_y_seen = max(max_y_seen, np.max(c_up[np.isfinite(c_up)]))

            ds_stats.append({
                "dataset": ds_name,
                "order": order,
                "summary": summary,
            })

        ax.set_xlabel(r"Temperature $T$ [MeV]")
        ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
        ax.set_title(rf"Chiral Susceptibility Multi-Order Padé Fits: ${ns_nt_labels[ds_name]}$")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f}"))
        ax.grid(True, linestyle="--", alpha=0.35)

        # 调整 Headroom 确保示意栏绝不遮挡数据和拟合峰值
        ax.set_ylim(-0.05, max_y_seen * 1.35)
        ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

        plt.tight_layout()
        out_png = FIGURES_DIR / f"chiral_susceptibility_multimodal_{ds_name}.png"
        out_pdf = FIGURES_DIR / f"chiral_susceptibility_multimodal_{ds_name}.pdf"
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_pdf)
        plt.close(fig)

        shutil.copyfile(out_png, BASE_DIR / ds_name / out_png.name)
        shutil.copyfile(out_pdf, BASE_DIR / ds_name / out_pdf.name)

        fit_stats[ds_name] = ds_stats

    return fit_stats


def generate_combined_plots() -> None:
    """生成 4 系综综合对比图 (A: p2 2/2, p3 0/2, q2 2/2, q3 0/2; B: p2 2/2, p3 2/2, q2 2/2, q3 0/2)"""
    apply_minimalist_style()

    styles = {
        "p2": {"color": "#1f77b4", "marker": "o", "label": r"$p_2\ (36^3\times 12)$"},
        "p3": {"color": "#2ca02c", "marker": "s", "label": r"$p_3\ (48^3\times 12)$"},
        "q2": {"color": "#d62728", "marker": "^", "label": r"$q_2\ (48^3\times 16)$"},
        "q3": {"color": "#9467bd", "marker": "D", "label": r"$q_3\ (64^3\times 16)$"},
    }

    cases = [
        (
            {"p2": "2/2", "p3": "0/2", "q2": "2/2", "q3": "0/2"},
            r"($p_2[2/2], p_3[0/2], q_2[2/2], q_3[0/2]$)",
            "chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02",
        ),
        (
            {"p2": "2/2", "p3": "2/2", "q2": "2/2", "q3": "0/2"},
            r"($p_2[2/2], p_3[2/2], q_2[2/2], q_3[0/2]$)",
            "chiral_lcp_comparison_p2_22_p3_22_q2_22_q3_02",
        ),
    ]

    for orders_map, title_suffix, filename_prefix in cases:
        print(f"[INFO] 正在生成 4 系综对比图 {filename_prefix}...")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.0), dpi=300)

        # 1. 左图: 格点单位手征凝聚
        max_cond_seen = 0.0
        for ds_name in ["p2", "p3", "q2", "q3"]:
            st = styles[ds_name]
            df_c = pl.read_csv(BASE_DIR / ds_name / f"condensate_{ds_name}.csv").sort("temp")
            temps = df_c["temp"].to_numpy().astype(float)
            pbp_lat = df_c["pbp_renorm"].to_numpy().astype(float)
            pbp_err_lat = df_c["pbp_renorm_err"].to_numpy().astype(float)
            max_cond_seen = max(max_cond_seen, np.max(pbp_lat + pbp_err_lat))

            ax1.errorbar(
                temps, pbp_lat, yerr=pbp_err_lat,
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
        # Headroom 避免遮挡
        ax1.set_ylim(-0.0001, max_cond_seen * 1.30)
        ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

        # 2. 右图: 易感度拟合带与数据点
        all_temps = []
        max_susc_seen = 0.0
        for ds_name, ord_val in orders_map.items():
            st = styles[ds_name]
            summary, records, temps, chi_vals, errs = compute_jackknife_for_order(ds_name, ord_val)
            all_temps.extend(temps)

            # 数据点
            ax2.errorbar(
                temps, chi_vals, yerr=errs,
                fmt=st["marker"], color=st["color"], capsize=3.5, markersize=5.5, alpha=0.75,
                label=f"{st['label']} Data",
            )
            max_susc_seen = max(max_susc_seen, np.max(chi_vals + errs))

        t_grid = np.linspace(min(all_temps) - 3.0, max(all_temps) + 3.0, 500)

        for ds_name, ord_val in orders_map.items():
            st = styles[ds_name]
            summary, records, _, _, _ = compute_jackknife_for_order(ds_name, ord_val)
            c_curve, low_b, up_b = compute_band_from_records(ord_val, records, t_grid)

            tpk = summary["t_peak_mean"]
            terr = summary["t_peak_err"]

            ax2.fill_between(t_grid, low_b, up_b, color=st["color"], alpha=0.15)
            ax2.plot(
                t_grid, c_curve, color=st["color"], linewidth=1.8,
                label=rf"{st['label']} [{ord_val}] ($T_c={tpk:.1f}\pm{terr:.1f}\,$MeV)",
            )
            max_susc_seen = max(max_susc_seen, np.max(up_b[np.isfinite(up_b)]))

        ax2.set_xlabel(r"Temperature $T$ [MeV]")
        ax2.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
        ax2.set_title(f"Chiral Susceptibility Comparison {title_suffix}")
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f}"))
        ax2.grid(True, linestyle="--", alpha=0.35)
        # Headroom 确保峰值与数据不被 legend 遮挡
        ax2.set_ylim(-0.05, max_susc_seen * 1.35)
        ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

        plt.tight_layout()
        out_png = FIGURES_DIR / f"{filename_prefix}.png"
        out_pdf = FIGURES_DIR / f"{filename_prefix}.pdf"
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_pdf)
        plt.close(fig)

        shutil.copyfile(out_png, BASE_DIR / out_png.name)
        shutil.copyfile(out_pdf, BASE_DIR / out_pdf.name)


def generate_lattice_condensate_plot() -> tuple[Path, Path]:
    """生成 4 系综在格点单位下的重整化手征凝聚对比图 (无任何 LT 字样，科学计数法刻度)"""
    apply_minimalist_style()

    styles = {
        "p2": {"color": "#1f77b4", "marker": "o", "label": r"$p_2\ (36^3\times 12)$"},
        "p3": {"color": "#2ca02c", "marker": "s", "label": r"$p_3\ (48^3\times 12)$"},
        "q2": {"color": "#d62728", "marker": "^", "label": r"$q_2\ (48^3\times 16)$"},
        "q3": {"color": "#9467bd", "marker": "D", "label": r"$q_3\ (64^3\times 16)$"},
    }

    fig, ax = plt.subplots(figsize=(8.5, 5.8), dpi=300)
    max_c_seen = 0.0

    for ds_name, st in styles.items():
        df = pl.read_csv(BASE_DIR / ds_name / f"condensate_{ds_name}.csv").sort("temp")
        t = df["temp"].to_numpy().astype(float)
        c = df["pbp_renorm"].to_numpy().astype(float)
        err = df["pbp_renorm_err"].to_numpy().astype(float)
        max_c_seen = max(max_c_seen, np.max(c + err))

        ax.errorbar(
            t, c, yerr=err,
            fmt=st["marker"], linestyle="none", color=st["color"], label=st["label"],
            capsize=3.5, markersize=5.5, alpha=0.9,
        )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\mathrm{sub}} / Z_m$ (lattice units)")
    ax.set_title(r"Renormalized Chiral Condensate in Lattice Units ($p_2, p_3, q_2, q_3$)")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.4f}"))
    ax.grid(True, linestyle="--", alpha=0.35)

    ax.set_ylim(-0.0001, max_c_seen * 1.30)
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = FIGURES_DIR / "chiral_condensate_lattice_units_all.png"
    out_pdf = FIGURES_DIR / "chiral_condensate_lattice_units_all.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    shutil.copyfile(out_png, BASE_DIR / out_png.name)
    shutil.copyfile(out_pdf, BASE_DIR / out_pdf.name)

    return out_png, out_pdf


def generate_bare_condensate_plot() -> tuple[Path, Path]:
    """生成 4 系综的裸 light 与 strange 凝聚数据散点图 (无连线)"""
    apply_minimalist_style()

    styles = {
        "p2": {"color": "#1f77b4", "marker": "o", "label": r"$p_2\ (36^3\times 12)$"},
        "p3": {"color": "#2ca02c", "marker": "s", "label": r"$p_3\ (48^3\times 12)$"},
        "q2": {"color": "#d62728", "marker": "^", "label": r"$q_2\ (48^3\times 16)$"},
        "q3": {"color": "#9467bd", "marker": "D", "label": r"$q_3\ (64^3\times 16)$"},
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 5.8), dpi=300)
    max_l = 0.0
    max_s = 0.0

    for ds_name, st in styles.items():
        df = pl.read_csv(BASE_DIR / ds_name / f"condensate_{ds_name}.csv").sort("temp")
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
    out_png = FIGURES_DIR / "chiral_condensate_bare_light_strange.png"
    out_pdf = FIGURES_DIR / "chiral_condensate_bare_light_strange.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    shutil.copyfile(out_png, BASE_DIR / out_png.name)
    shutil.copyfile(out_pdf, BASE_DIR / out_pdf.name)
    return out_png, out_pdf


def format_sci(val: float, err: float | None = None, digits: int = 2) -> str:
    """将数值格式化为标准 LaTeX 科学计数法形式 (用 $...$ 包含)"""
    if val is None or not np.isfinite(val):
        return "--"
    if err is not None and np.isfinite(err) and err > 0:
        exp = int(np.floor(np.log10(abs(val)))) if val != 0 else 0
        mant = val / (10**exp)
        mant_err = err / (10**exp)
        if exp == 0:
            return f"${val:.{digits}f} \\pm {err:.{digits}f}$"
        else:
            return f"$({mant:.{digits}f} \\pm {mant_err:.{digits}f}) \\times 10^{{{exp}}}$"
    else:
        exp = int(np.floor(np.log10(abs(val)))) if val != 0 else 0
        mant = val / (10**exp)
        if exp == 0:
            return f"${val:.{digits}f}$"
        else:
            return f"${mant:.{digits}f} \\times 10^{{{exp}}}$"


def generate_pure_report_tex() -> None:
    """生成完全不包含任何 section/title/author、纯图表的 report.tex"""
    tex = []
    tex.append(r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=1.5cm]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{amsmath,amssymb}
\usepackage{subcaption}
\usepackage{hyperref}
\usepackage{microtype}
\usepackage{siunitx}

\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    citecolor=blue,
    urlcolor=blue
}

\graphicspath{{figures/}{../22pade/}}

\pagestyle{plain}

\begin{document}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02.pdf}
    \caption{Renormalized chiral condensate $\langle\bar{\psi}\psi\rangle_{\rm sub}/Z_m$ in lattice units (left) and renormalized chiral susceptibility $\chi_{\mathrm{disc}}/Z_m^2\,[\mathrm{GeV}^2]$ (right) comparison across all four LCP ensembles with Padé configuration: $p_2\,[2/2]$, $p_3\,[0/2]$, $q_2\,[2/2]$, and $q_3\,[0/2]$. Shaded bands indicate Jackknife envelope uncertainty, and solid lines denote mean parameter fits.}
    \label{fig:combined_p3_02_q3_02}
\end{figure}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/chiral_condensate_bare_light_strange.pdf}
    \caption{Bare light quark condensate $\langle\bar{\psi}\psi\rangle_l$ (left) and bare strange quark condensate $\langle\bar{\psi}\psi\rangle_s$ (right) in lattice units across all four LCP ensembles ($p_2, p_3, q_2, q_3$).}
    \label{fig:bare_light_strange}
\end{figure}

\clearpage

\begin{figure}[htbp]
    \centering
    \begin{subfigure}[b]{0.49\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/chiral_susceptibility_multimodal_p2.pdf}
        \caption{$p_2$ ($36^3\times 12$)}
        \label{fig:multi_p2}
    \end{subfigure}
    \hfill
    \begin{subfigure}[b]{0.49\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/chiral_susceptibility_multimodal_p3.pdf}
        \caption{$p_3$ ($48^3\times 12$)}
        \label{fig:multi_p3}
    \end{subfigure}
    \vskip\baselineskip
    \begin{subfigure}[b]{0.49\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/chiral_susceptibility_multimodal_q2.pdf}
        \caption{$q_2$ ($48^3\times 16$)}
        \label{fig:multi_q2}
    \end{subfigure}
    \hfill
    \begin{subfigure}[b]{0.49\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/chiral_susceptibility_multimodal_q3.pdf}
        \caption{$q_3$ ($64^3\times 16$)}
        \label{fig:multi_q3}
    \end{subfigure}
    \caption{Multi-order Padé approximation comparison per ensemble with shaded Jackknife confidence bands. Spurious orders suffering from real-axis poles or zero degrees of freedom have been excluded.}
    \label{fig:multimodal_all}
\end{figure}

\clearpage

\begin{table}[htbp]
    \centering
    \caption{Jackknife fit parameters and thermal transition results for primary Padé configurations.}
    \label{tab:lcp_summary}
    \small
    \begin{tabular}{lccccc}
        \toprule
        \textbf{Ensemble} & \textbf{Dimensions} & \textbf{Order} & $\boldsymbol{T_{pc}}$ [MeV] & $\boldsymbol{\chi_{\mathrm{peak}} / Z_m^2}$ [$\text{GeV}^2$] & $\boldsymbol{\chi^2/\mathrm{dof}}$ \\
        \midrule
""")

    # 填入核心系综数据 (仅保留图1展示的 4 个系综设定)
    core_cases = [
        ("p2", "36^3\\times 12", "2/2"),
        ("p3", "48^3\\times 12", "0/2"),
        ("q2", "48^3\\times 16", "2/2"),
        ("q3", "64^3\\times 16", "0/2"),
    ]
    for ds_name, dims, ord_val in core_cases:
        sm, recs, _, _, _ = compute_jackknife_for_order(ds_name, ord_val)
        tpk = sm["t_peak_mean"]
        terr = sm["t_peak_err"]
        cpk = sm["chi_peak_mean"]
        cerr = sm["chi_peak_err"]
        chi2_dof = sm["chi2_dof_mean"]
        ens_tex = f"${ds_name[0]}_{{{ds_name[1:]}}}$"

        tex.append(f"        {ens_tex} & ${dims}$ & [{ord_val}] & ${tpk:.2f}\\pm {terr:.2f}$ & ${cpk:.3f}\\pm {cerr:.3f}$ & ${chi2_dof:.3f}$ \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{Model comparison across all Padé rational approximations.}
    \label{tab:model_comparison}
    \small
    \begin{tabular}{lcccc}
        \toprule
        \textbf{Ensemble} & \textbf{Order} & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & $\boldsymbol{T_{\mathrm{peak}}}$ [MeV] \\
        \midrule
""")

    # 填入筛选后的各个阶数统计
    selected_fits = [
        ("p2", ["1/2", "2/2", "3/2"]),
        ("p3", ["0/2", "1/2", "2/2"]),
        ("q2", ["1/2", "2/2", "3/2", "1/3"]),
        ("q3", ["0/2", "1/2"]),
    ]
    for ds_name, ords in selected_fits:
        for ord_val in ords:
            sm, _, _, _, _ = compute_jackknife_for_order(ds_name, ord_val)
            model = PadeModel(order=ord_val, t0=157.0)
            k_params = len(model.param_keys)
            chi2_dof = sm["chi2_dof_mean"]
            tpk = sm["t_peak_mean"]
            ens_tex = f"${ds_name[0]}_{{{ds_name[1:]}}}$"
            tex.append(f"        {ens_tex} & [{ord_val}] & {k_params} & ${chi2_dof:.3f}$ & ${tpk:.2f}$ \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\clearpage
""")

    # Table 3: p2 & p3 原始数据
    tex.append(r"""\begin{table}[htbp]
    \centering
    \caption{Lattice measurements for $N_t=12$ ensembles ($p_2$ and $p_3$): Temperature $T$, inverse lattice spacing $a^{-1}$, renormalized chiral condensate $\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m$ in lattice units, and renormalized susceptibility $\chi_{\mathrm{disc}} / Z_m^2$ [$\text{GeV}^2$].}
    \label{tab:raw_p2_p3}
    \small
    \begin{tabular}{lccccc}
        \toprule
        \textbf{Ensemble} & $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{a^{-1}}$ [GeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m}$ & $\boldsymbol{\chi_{\mathrm{disc}} / Z_m^2}$ [$\text{GeV}^2$] \\
        \midrule
""")

    df_p2_c = pl.read_csv(BASE_DIR / "p2/condensate_p2.csv").sort("temp")
    df_p2_s = pl.read_csv(BASE_DIR / "p2/susceptibility_p2.csv").sort("temp")
    p2_c_dict = {float(r["beta"]): r for r in df_p2_c.to_dicts()}

    for r in df_p2_s.to_dicts():
        b = float(r["beta"])
        t_val = float(r["temp"])
        chi_val = float(r["chi_renorm"])
        chi_err = float(r["chi_renorm_err"])
        c_row = p2_c_dict.get(b, {})
        a_inv = float(c_row.get("a_inv", 12 * t_val)) / 1000.0
        pbp_lat = float(c_row.get("pbp_renorm", 0.0))
        pbp_err_lat = float(c_row.get("pbp_renorm_err", 0.0))
        pbp_str = format_sci(pbp_lat, pbp_err_lat, digits=2)
        chi_str = format_sci(chi_val, chi_err, digits=3)
        tex.append(f"        $p_2$ ($36^3\\times 12$) & {b:.3f} & {t_val:.1f} & {a_inv:.3f} & {pbp_str} & {chi_str} \\\\")

    tex.append(r"""        \midrule
""")

    df_p3_c = pl.read_csv(BASE_DIR / "p3/condensate_p3.csv").sort("temp")
    df_p3_s = pl.read_csv(BASE_DIR / "p3/susceptibility_p3.csv").sort("temp")
    p3_c_dict = {float(r["beta"]): r for r in df_p3_c.to_dicts()}

    for r in df_p3_s.to_dicts():
        b = float(r["beta"])
        t_val = float(r["temp"])
        chi_val = float(r["chi_renorm"])
        chi_err = float(r["chi_renorm_err"])
        c_row = p3_c_dict.get(b, {})
        a_inv = float(c_row.get("a_inv", 12 * t_val)) / 1000.0
        pbp_lat = float(c_row.get("pbp_renorm", 0.0))
        pbp_err_lat = float(c_row.get("pbp_renorm_err", 0.0))
        pbp_str = format_sci(pbp_lat, pbp_err_lat, digits=2)
        chi_str = format_sci(chi_val, chi_err, digits=3)
        tex.append(f"        $p_3$ ($48^3\\times 12$) & {b:.3f} & {t_val:.1f} & {a_inv:.3f} & {pbp_str} & {chi_str} \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{Lattice measurements for $N_t=16$ ensembles ($q_2$ and $q_3$): Temperature $T$, inverse lattice spacing $a^{-1}$, renormalized chiral condensate $\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m$ in lattice units, and renormalized susceptibility $\chi_{\mathrm{disc}} / Z_m^2$ [$\text{GeV}^2$].}
    \label{tab:raw_q2_q3}
    \small
    \begin{tabular}{lccccc}
        \toprule
        \textbf{Ensemble} & $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{a^{-1}}$ [GeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m}$ & $\boldsymbol{\chi_{\mathrm{disc}} / Z_m^2}$ [$\text{GeV}^2$] \\
        \midrule
""")

    df_q2_c = pl.read_csv(BASE_DIR / "q2/condensate_q2.csv").sort("temp")
    df_q2_s = pl.read_csv(BASE_DIR / "q2/susceptibility_q2.csv").sort("temp")
    q2_c_dict = {float(r["beta"]): r for r in df_q2_c.to_dicts()}

    for r in df_q2_s.to_dicts():
        b = float(r["beta"])
        t_val = float(r["temp"])
        chi_val = float(r["chi_renorm"])
        chi_err = float(r["chi_renorm_err"])
        c_row = q2_c_dict.get(b, {})
        a_inv = float(c_row.get("a_inv", 16 * t_val)) / 1000.0
        pbp_lat = float(c_row.get("pbp_renorm", 0.0))
        pbp_err_lat = float(c_row.get("pbp_renorm_err", 0.0))
        pbp_str = format_sci(pbp_lat, pbp_err_lat, digits=2)
        chi_str = format_sci(chi_val, chi_err, digits=3)
        tex.append(f"        $q_2$ ($48^3\\times 16$) & {b:.3f} & {t_val:.1f} & {a_inv:.3f} & {pbp_str} & {chi_str} \\\\")

    tex.append(r"""        \midrule
""")

    df_q3_c = pl.read_csv(BASE_DIR / "q3/condensate_q3.csv").sort("temp")
    df_q3_s = pl.read_csv(BASE_DIR / "q3/susceptibility_q3.csv").sort("temp")
    q3_c_dict = {float(r["beta"]): r for r in df_q3_c.to_dicts()}

    for r in df_q3_s.to_dicts():
        b = float(r["beta"])
        t_val = float(r["temp"])
        chi_val = float(r["chi_renorm"])
        chi_err = float(r["chi_renorm_err"])
        c_row = q3_c_dict.get(b, {})
        a_inv = float(c_row.get("a_inv", 16 * t_val)) / 1000.0
        pbp_lat = float(c_row.get("pbp_renorm", 0.0))
        pbp_err_lat = float(c_row.get("pbp_renorm_err", 0.0))
        pbp_str = format_sci(pbp_lat, pbp_err_lat, digits=2)
        chi_str = format_sci(chi_val, chi_err, digits=3)
        tex.append(f"        $q_3$ ($64^3\\times 16$) & {b:.3f} & {t_val:.1f} & {a_inv:.3f} & {pbp_str} & {chi_str} \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\end{document}
""")

    tex_path = DOCS_DIR / "report.tex"
    tex_path.write_text("\n".join(tex), encoding="utf-8")
    print(f"[INFO] 成功写入 {tex_path}")


def compile_pdf() -> None:
    print("[INFO] 正在编译 report.tex 为 PDF...")
    for ext in [".aux", ".out", ".log", ".toc"]:
        f = DOCS_DIR / f"report{ext}"
        if f.exists():
            f.unlink()

    cmd = ["/Library/TeX/texbin/pdflatex", "-interaction=nonstopmode", "report.tex"]
    proc = subprocess.run(cmd, cwd=DOCS_DIR, capture_output=True, text=True)
    if proc.returncode != 0:
        print("[ERROR] pdflatex 编译失败:")
        print(proc.stdout[-2000:])
        sys.exit(1)

    subprocess.run(cmd, cwd=DOCS_DIR, capture_output=True, text=True)
    pdf_path = DOCS_DIR / "report.pdf"
    if pdf_path.exists():
        print(f"[SUCCESS] PDF 编译完成: {pdf_path} (大小: {pdf_path.stat().st_size / 1024:.1f} KB)")
    else:
        print("[ERROR] 未找到生成的 report.pdf")
        sys.exit(1)


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    generate_lattice_condensate_plot()
    generate_bare_condensate_plot()
    generate_multimodal_plots()
    generate_combined_plots()
    generate_pure_report_tex()
    compile_pdf()
    print("[DONE] 所有更新与编译流程全部顺利完成！")


if __name__ == "__main__":
    main()
