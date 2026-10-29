#!/usr/bin/env python3
"""
scripts/tasks/lcp/run_22pade_pipeline.py
--------------------------------------------------------------------------------
1. 手征易感度统一采用 Padé [2/2] 进行 Jackknife 拟合与推断:
   - 拟合曲线用阴影画出拟合出的所有曲线的范围 (Jackknife 拟合包络/误差范围)。
   - 实线表示平均出的参数所计算出的拟合曲线。
   - 原数据点带误差棒，并在图例下方放置信息框展示 Tc 与 chi^2/dof。
2. 覆盖四大核心数据集:
   - p2: 36^3 x 12
   - p3: 48^3 x 12
   - q2: 48^3 x 16
   - q3: 64^3 x 16
3. 所有输出集中归档至 22pade/ (并在 output/LCP/22pade/ 建立统一映射)。
4. 清理中间数据之外的之前做的不需要的历史旧图。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import sys
import warnings

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
from src.parselqcdata.pade_fitter import (
    PadeModel,
    fit_single_sample_pade,
    fit_jackknife_susceptibility_pade,
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


def run_jackknife_pade22_for_dataset(
    dataset_name: str,
    t0: float = 157.0,
    max_samples: int = 2000,
) -> tuple[dict, list[dict], pl.DataFrame, pl.DataFrame]:
    """对单个数据集执行 Padé [2/2] Jackknife 拟合"""
    susc_csv = OUTPUT_LCP_DIR / dataset_name / f"susceptibility_{dataset_name}.csv"
    jk_csv = OUTPUT_LCP_DIR / dataset_name / f"susceptibility_jk_samples_{dataset_name}.csv"
    cond_csv = OUTPUT_LCP_DIR / dataset_name / f"condensate_{dataset_name}.csv"

    if not susc_csv.exists() or not jk_csv.exists() or not cond_csv.exists():
        raise FileNotFoundError(f"找不到数据集 {dataset_name} 的完整数据产物")

    df_susc = pl.read_csv(susc_csv).sort("temp")
    df_jk = pl.read_csv(jk_csv)
    df_cond = pl.read_csv(cond_csv).sort("temp")

    temps = df_susc["temp"].to_numpy().astype(float)
    errors = df_susc["chi_renorm_err"].to_numpy().astype(float)
    betas = df_susc["beta"].to_numpy().astype(float)

    # 提取对齐的样本矩阵
    samples_per_temp = []
    for b in betas:
        s_vals = df_jk.filter(pl.col("beta") == b).sort("jk_index")["chi_renorm"].to_numpy().astype(float)
        samples_per_temp.append(s_vals)

    min_len = min(len(s) for s in samples_per_temp)
    target_n = min(min_len, max_samples)
    chi_matrix = np.array([s[:target_n] for s in samples_per_temp], dtype=float)

    summary, records = fit_jackknife_susceptibility_pade(
        temps=temps,
        chi_matrix=chi_matrix,
        errors=errors,
        order="2/2",
        t0=t0,
    )

    return summary, records, df_susc, df_cond


def plot_pade22_with_shaded_envelope(
    dataset_name: str,
    summary: dict,
    records: list[dict],
    df_susc: pl.DataFrame,
    out_dir: Path,
    t0: float = 157.0,
) -> tuple[Path, Path]:
    """
    绘制易感度 Padé [2/2] 拟合图:
    - 拟合曲线用阴影画出拟合出的所有曲线的范围 (Jackknife 误差包络/曲线极值范围)
    - 实线表示平均出的参数所计算出的拟合曲线
    - 数据点带误差棒
    - 右上角带 Tc 与 chi^2/dof 信息框
    """
    apply_minimalist_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    temps_data = df_susc["temp"].to_numpy().astype(float)
    chi_data = df_susc["chi_renorm"].to_numpy().astype(float)
    err_data = df_susc["chi_renorm_err"].to_numpy().astype(float)

    t_peak_mean = float(summary["t_peak_mean"])
    t_peak_err = float(summary["t_peak_err"])
    chi2_dof = float(summary["chi2_dof_mean"])

    # 1. 过滤有效样本计算所有拟合曲线矩阵
    valid_recs = [r for r in records if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    n_samples = len(valid_recs)

    t_grid = np.linspace(np.min(temps_data) - 4.0, np.max(temps_data) + 4.0, 500)
    dx = t_grid - t0

    model = PadeModel(order="2/2", t0=t0)

    # 计算各参数的平均值 (实线表示平均出的参数表示的曲线)
    p_keys = ["p0", "p1", "p2", "q1", "q2"]
    p_mean_dict: dict[str, float] = {}
    for pk in p_keys:
        p_vals = [float(r[pk]) for r in valid_recs if pk in r and np.isfinite(r[pk])]
        p_mean_dict[pk] = float(np.mean(p_vals)) if p_vals else 0.0

    # 实线: 平均出的参数表示的曲线
    curve_from_mean_params = model.eval(t_grid, p_mean_dict)

    # 逐样本求值构建曲线矩阵
    curves_matrix = np.zeros((n_samples, len(t_grid)), dtype=float)
    for i, r in enumerate(valid_recs):
        p_i = {pk: float(r[pk]) for pk in p_keys if pk in r}
        curves_matrix[i, :] = model.eval(t_grid, p_i)

    # 拟合范围阴影: Jackknife 误差包络范围 (与以前 05_plot_pade_susceptibility 严格对齐)
    jk_curve_err = np.sqrt(n_samples - 1) * np.std(curves_matrix, axis=0, ddof=0)
    y_lower = curve_from_mean_params - jk_curve_err
    y_upper = curve_from_mean_params + jk_curve_err

    # ① 阴影表示拟合的范围 (Fit range)
    ax.fill_between(
        t_grid,
        y_lower,
        y_upper,
        color="#1f77b4",
        alpha=0.25,
        label="Fit range",
        zorder=1,
    )

    # ② 拟合的平均曲线 (Fit mean curve)
    ax.plot(
        t_grid,
        curve_from_mean_params,
        color="#1f77b4",
        linewidth=2.0,
        label="Fit",
        zorder=2,
    )

    # ③ 原数据点及误差棒 (Data points)
    ax.errorbar(
        temps_data,
        chi_data,
        yerr=err_data,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.4,
        capsize=3.5,
        markersize=5.5,
        label="Data",
        zorder=3,
    )

    # ④ 图例置于右上角
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # ⑤ 数据方框放置在右上角图例下方
    chi2_str = f"{chi2_dof:.3f}" if np.isfinite(chi2_dof) else "0.000 (Exact dof=0)"
    info_text = (
        rf"$T_c = {t_peak_mean:.2f} \pm {t_peak_err:.2f}\ \mathrm{{MeV}}$" + "\n"
        rf"$\chi^2/\mathrm{{dof}} = {chi2_str}$"
    )
    ax.text(
        0.96,
        0.72,
        info_text,
        transform=ax.transAxes,
        fontsize=11.5,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        zorder=4,
    )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax.set_title(rf"Chiral Susceptibility Padé [2/2] Fit [{dataset_name.upper()}]")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")

    plt.tight_layout()
    out_png = out_dir / f"chiral_susceptibility_pade22_{dataset_name}.png"
    out_pdf = out_dir / f"chiral_susceptibility_pade22_{dataset_name}.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def plot_condensate_clean(
    dataset_name: str,
    df_cond: pl.DataFrame,
    out_dir: Path,
) -> tuple[Path, Path]:
    """绘制出版级手征凝聚图"""
    apply_minimalist_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    temps = df_cond["temp"].to_numpy().astype(float)
    pbp_ren = df_cond["pbp_renorm"].to_numpy().astype(float)
    pbp_err = df_cond["pbp_renorm_err"].to_numpy().astype(float)

    scale = 1e4
    ax.errorbar(
        temps,
        pbp_ren * scale,
        yerr=pbp_err * scale,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.4,
        capsize=3.5,
        markersize=5.5,
        label=rf"$\langle \bar{{\psi}}\psi \rangle_{{{{\rm sub}}}} / Z_m$ ({dataset_name})",
        zorder=3,
    )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m \times 10^4$")
    ax.set_title(rf"Renormalized Chiral Condensate [{dataset_name.upper()}]")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = out_dir / f"chiral_condensate_{dataset_name}.png"
    out_pdf = out_dir / f"chiral_condensate_{dataset_name}.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def plot_unified_comparison(
    all_results: dict[str, dict],
    out_dir: Path,
    t0: float = 157.0,
) -> tuple[Path, Path]:
    """绘制 4 体系综合对比大图"""
    apply_minimalist_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 5.8), dpi=300)

    colors = {
        "p2": "#1f77b4",
        "p3": "#2ca02c",
        "q2": "#d62728",
        "q3": "#9467bd",
    }
    labels = {
        "p2": r"$p_2\ (36^3\times 12)$",
        "p3": r"$p_3\ (48^3\times 12)$",
        "q2": r"$q_2\ (48^3\times 16)$",
        "q3": r"$q_3\ (64^3\times 16)$",
    }

    # 左图: 手征凝聚
    for ds_name, item in all_results.items():
        df_c = item["df_cond"]
        ax1.errorbar(
            df_c["temp"].to_numpy().astype(float),
            df_c["pbp_renorm"].to_numpy().astype(float) * 1e4,
            yerr=df_c["pbp_renorm_err"].to_numpy().astype(float) * 1e4,
            fmt="o-",
            color=colors[ds_name],
            label=labels[ds_name],
            capsize=3.0,
            markersize=5.0,
            linewidth=1.2,
        )

    ax1.set_xlabel(r"Temperature $T$ [MeV]")
    ax1.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m \times 10^4$")
    ax1.set_title(r"Chiral Condensate Comparison ($p_2, p_3, q_2, q_3$)")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # 右图: 易感度 Padé [2/2] 拟合曲线与数据点
    for ds_name, item in all_results.items():
        df_s = item["df_susc"]
        sum_dict = item["summary"]
        recs = item["records"]

        # 数据点
        ax2.errorbar(
            df_s["temp"].to_numpy().astype(float),
            df_s["chi_renorm"].to_numpy().astype(float),
            yerr=df_s["chi_renorm_err"].to_numpy().astype(float),
            fmt="o",
            color=colors[ds_name],
            capsize=3.0,
            markersize=5.0,
            alpha=0.7,
        )

        # Padé 拟合曲线
        model = PadeModel(order="2/2", t0=t0)
        p_keys = ["p0", "p1", "p2", "q1", "q2"]
        p_mean_dict = {pk: float(np.mean([float(r[pk]) for r in recs if pk in r and np.isfinite(r[pk])])) for pk in p_keys}

        temps_arr = df_s["temp"].to_numpy().astype(float)
        t_grid = np.linspace(np.min(temps_arr) - 3.0, np.max(temps_arr) + 3.0, 300)
        y_curve = model.eval(t_grid, p_mean_dict)

        tpk = float(sum_dict["t_peak_mean"])
        terr = float(sum_dict["t_peak_err"])
        ax2.plot(
            t_grid,
            y_curve,
            color=colors[ds_name],
            linewidth=1.8,
            label=rf"{labels[ds_name]} ($T_{{pc}}={tpk:.1f}\pm{terr:.1f}\,$MeV)",
        )

    ax2.set_xlabel(r"Temperature $T$ [MeV]")
    ax2.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax2.set_title(r"Chiral Susceptibility Padé [2/2] Fit Comparison")
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = out_dir / "chiral_lcp_pade22_all_comparison.png"
    out_pdf = out_dir / "chiral_lcp_pade22_all_comparison.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def cleanup_obsolete_figures() -> None:
    """清理中间数据之外的历史不必要图件"""
    print("\n[CLEANUP] 开始清理历史冗余图件 (保留所有数据表与统计样本)...")

    # docs/figures 冗余图清理
    obsolete_patterns = [
        "chiral_condensate_L48.*",
        "chiral_condensate_L64.*",
        "chiral_susceptibility_L48.*",
        "chiral_susceptibility_L64.*",
        "chiral_susceptibility_pade_fit.*",
        "chiral_susceptibility_pade22_fit.*",
        "chiral_condensate.png",
        "chiral_condensate.pdf",
        "chiral_susceptibility.png",
        "chiral_susceptibility.pdf",
        "chiral_susceptibility_pade_p2.*",
        "chiral_susceptibility_pade_p3.*",
        "chiral_susceptibility_pade_q2.*",
        "chiral_susceptibility_pade_q3.*",
        "conden_all_ensembles_overview.*",
        "conden_plot.*",
        "conden_re_plot.*",
        "conden_scaling_beta4.17.*",
        "conden_volume_scaling.*",
        "pbp_temperature_plot.*",
        "pbpchisce*.*",
        "sucep_plot.*",
    ]

    for pat in obsolete_patterns:
        for f in FIGURES_DIR.glob(pat):
            try:
                f.unlink()
                print(f"  - 已删除 docs/figures/{f.name}")
            except Exception:
                pass

    # output/condensate 冗余图清理
    for pat in obsolete_patterns:
        for f in OUTPUT_CONDENSATE_DIR.glob(pat):
            try:
                f.unlink()
                print(f"  - 已删除 output/condensate/{f.name}")
            except Exception:
                pass

    # output/LCP 根目录旧图清理
    for f in OUTPUT_LCP_DIR.glob("*.png"):
        try:
            f.unlink()
            print(f"  - 已删除 output/LCP/{f.name}")
        except Exception:
            pass
    for f in OUTPUT_LCP_DIR.glob("*.pdf"):
        try:
            f.unlink()
            print(f"  - 已删除 output/LCP/{f.name}")
        except Exception:
            pass

    # output/ 根目录下旧图清理
    for f in (PROJECT_ROOT / "output").glob("*.png"):
        try:
            f.unlink()
            print(f"  - 已删除 output/{f.name}")
        except Exception:
            pass

    print("[CLEANUP] 历史冗余图件清理完毕！")


def main() -> None:
    t0 = 157.0

    # 统一产物目录 22pade/
    out_22pade_lcp = OUTPUT_LCP_DIR / "22pade"
    out_22pade_root = PROJECT_ROOT / "22pade"
    out_22pade_lcp.mkdir(parents=True, exist_ok=True)
    out_22pade_root.mkdir(parents=True, exist_ok=True)

    datasets = ["p2", "p3", "q2", "q3"]
    all_results: dict[str, dict] = {}
    summary_rows = []

    print("\n" + "=" * 80)
    print(">>> [22PADE PIPELINE] 启动手征易感度 Padé [2/2] Jackknife 拟合与产物归档...")
    print("=" * 80)

    for ds in datasets:
        print(f"\n>>> 正在处理 {ds.upper()} 系综 Padé [2/2] 拟合...")
        summary, records, df_susc, df_cond = run_jackknife_pade22_for_dataset(ds, t0=t0)

        # 归档到 22pade/<ds>/
        ds_dir = out_22pade_lcp / ds
        ds_dir.mkdir(parents=True, exist_ok=True)

        # 导出 csv
        summary_df = pl.DataFrame([summary])
        summary_df.write_csv(ds_dir / f"pade22_fit_summary_{ds}.csv")

        samples_df = pl.DataFrame(records)
        samples_df.write_csv(ds_dir / f"pade22_jk_samples_{ds}.csv")

        df_susc.write_csv(ds_dir / f"susceptibility_{ds}.csv")
        df_cond.write_csv(ds_dir / f"condensate_{ds}.csv")

        # 1. 绘制带阴影包络的易感度拟合图
        susc_png, susc_pdf = plot_pade22_with_shaded_envelope(
            dataset_name=ds,
            summary=summary,
            records=records,
            df_susc=df_susc,
            out_dir=ds_dir,
            t0=t0,
        )
        print(f"  [PLOT] 易感度 Padé [2/2] 阴影拟合图已生成: {susc_png.name}")

        # 2. 绘制手征凝聚图
        cond_png, cond_pdf = plot_condensate_clean(
            dataset_name=ds,
            df_cond=df_cond,
            out_dir=ds_dir,
        )
        print(f"  [PLOT] 手征凝聚图已生成: {cond_png.name}")

        all_results[ds] = {
            "summary": summary,
            "records": records,
            "df_susc": df_susc,
            "df_cond": df_cond,
        }

        ns_val = int(df_susc["ns"][0])
        nt_val = int(df_susc["nt"][0])
        summary_rows.append({
            "dataset": ds,
            "ns": ns_val,
            "nt": nt_val,
            "n_pts": df_susc.height,
            "n_samples": summary["n_samples"],
            "n_fitted": summary["n_fitted"],
            "chi2_dof_mean": summary["chi2_dof_mean"],
            "t_peak_mean": summary["t_peak_mean"],
            "t_peak_err": summary["t_peak_err"],
            "chi_peak_mean": summary["chi_peak_mean"],
            "chi_peak_err": summary["chi_peak_err"],
            "t_pole_real_mean": summary["t_pole_real_mean"],
            "t_pole_real_err": summary["t_pole_real_err"],
            "t_pole_imag_mean": summary["t_pole_imag_mean"],
            "t_pole_imag_err": summary["t_pole_imag_err"],
        })

    # 生成全局汇总表
    df_all_sum = pl.DataFrame(summary_rows)
    df_all_sum.write_csv(out_22pade_lcp / "all_datasets_pade22_summary.csv")
    print(f"\n[OK] 全体系 Padé [2/2] 汇总表已写出: all_datasets_pade22_summary.csv")

    # 生成 4 体系综合对比大图
    comp_png, comp_pdf = plot_unified_comparison(all_results, out_22pade_lcp, t0=t0)
    print(f"[PLOT] 四体系手征相变综合对比大图已生成: {comp_png.name}")

    # 将介子与奇异夸克凝聚拟合产物也整合至 22pade/
    src_meson_dir = OUTPUT_LCP_DIR / "meson_strange_fit"
    dst_meson_dir = out_22pade_lcp / "meson_strange_fit"
    if src_meson_dir.exists():
        dst_meson_dir.mkdir(parents=True, exist_ok=True)
        for f in src_meson_dir.iterdir():
            if f.is_file():
                shutil.copyfile(f, dst_meson_dir / f.name)
        print(f"[SYNC] 介子与奇异夸克凝聚拟合产物已同步整合至 22pade/meson_strange_fit/")

    # 同步整个目录结构到根目录 22pade/
    for item in out_22pade_lcp.rglob("*"):
        rel = item.relative_to(out_22pade_lcp)
        target = out_22pade_root / rel
        if item.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        else:
            shutil.copyfile(item, target)
    print(f"[SYNC] 22pade/ 产物目录已完整镜像至根目录: {out_22pade_root}")

    # 清理历史无用图件
    cleanup_obsolete_figures()

    print("\n" + "=" * 80)
    print(">>> [ALL COMPLETE] Padé [2/2] 全面分析与 22pade/ 归档执行完毕！")
    print("=" * 80)


if __name__ == "__main__":
    main()
