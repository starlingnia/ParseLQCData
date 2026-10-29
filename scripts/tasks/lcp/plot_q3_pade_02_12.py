#!/usr/bin/env python3
"""
scripts/tasks/lcp/plot_q3_pade_02_12.py
--------------------------------------------------------------------------------
1. 为 q3 (64^3 x 16) 单独执行 Padé [0/2] 与 Padé [1/2] 的完整 Jackknife 拟合。
2. 绘制符合学术规范的高清易感度拟合图：
   - 阴影区域表示 Jackknife 拟合出的所有曲线的包络范围 (Fit range)
   - 实线表示由平均参数计算出的拟合曲线 (Fit)
   - 原数据点带统计误差棒 (Data)
   - 信息框清晰标注 Tc 与 chi^2/dof
3. 生成单阶数图 (0/2, 1/2) 以及 0/2 与 1/2 的同图对比图。
4. 修复与更新 q2 的手征凝聚数据 (清理混杂的 L32T12 非本征数据点，严格保持 8 个标准点)。
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

from docs.physics_setup import FIGURES_DIR, OUTPUT_LCP_DIR, OUTPUT_CONDENSATE_DIR
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
        "legend.fontsize": 10.5,
        "mathtext.fontset": "cm",
    })


def run_q3_pade_fits(t0: float = 157.0) -> dict[str, dict]:
    """对 q3 执行 0/2 与 1/2 的 Jackknife 拟合"""
    susc_csv = PROJECT_ROOT / "22pade" / "q3" / "susceptibility_q3.csv"
    jk_csv = OUTPUT_LCP_DIR / "susceptibility_jk_samples_L64.csv"

    df_susc = pl.read_csv(susc_csv).sort("temp")
    df_jk = pl.read_csv(jk_csv)

    temps = df_susc["temp"].to_numpy().astype(float)
    errors = df_susc["chi_renorm_err"].to_numpy().astype(float)
    betas = df_susc["beta"].to_numpy().astype(float)

    samples_per_temp = []
    for b in betas:
        s_vals = df_jk.filter(pl.col("beta") == b).sort("jk_index")["chi_renorm"].to_numpy().astype(float)
        samples_per_temp.append(s_vals)

    min_len = min(len(s) for s in samples_per_temp)
    chi_matrix = np.array([s[:min_len] for s in samples_per_temp], dtype=float)

    results = {}
    for order in ["0/2", "1/2"]:
        summary, records = fit_jackknife_susceptibility_pade(
            temps=temps,
            chi_matrix=chi_matrix,
            errors=errors,
            order=order,
            t0=t0,
        )
        results[order] = {
            "summary": summary,
            "records": records,
            "temps": temps,
            "chi": df_susc["chi_renorm"].to_numpy().astype(float),
            "errors": errors,
            "df_susc": df_susc,
        }

    return results


def plot_single_q3_pade(
    order: str,
    data: dict,
    out_dir: Path,
    t0: float = 157.0,
) -> tuple[Path, Path]:
    """绘制单个 Padé 阶数的拟合图 (带阴影与平均实线)"""
    apply_minimalist_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    summary = data["summary"]
    records = data["records"]
    temps_data = data["temps"]
    chi_data = data["chi"]
    err_data = data["errors"]

    t_peak_mean = float(summary["t_peak_mean"])
    t_peak_err = float(summary["t_peak_err"])
    chi2_dof = float(summary["chi2_dof_mean"])

    valid_recs = [r for r in records if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    n_samples = len(valid_recs)

    t_grid = np.linspace(np.min(temps_data) - 4.0, np.max(temps_data) + 4.0, 500)
    model = PadeModel(order=order, t0=t0)
    p_keys = model.param_keys

    # 计算平均参数实线
    p_mean_dict = {pk: float(np.mean([float(r[pk]) for r in valid_recs if pk in r and np.isfinite(r[pk])])) for pk in p_keys}
    curve_mean = model.eval(t_grid, p_mean_dict)

    # 逐样本计算包络阴影
    curves_matrix = np.zeros((n_samples, len(t_grid)), dtype=float)
    for i, r in enumerate(valid_recs):
        p_i = {pk: float(r[pk]) for pk in p_keys if pk in r}
        curves_matrix[i, :] = model.eval(t_grid, p_i)

    jk_curve_err = np.sqrt(n_samples - 1) * np.std(curves_matrix, axis=0, ddof=0)
    y_lower = curve_mean - jk_curve_err
    y_upper = curve_mean + jk_curve_err

    # 阴影包络
    ax.fill_between(
        t_grid,
        y_lower,
        y_upper,
        color="#1f77b4",
        alpha=0.25,
        label="Fit range",
        zorder=1,
    )

    # 实线拟合
    ax.plot(
        t_grid,
        curve_mean,
        color="#1f77b4",
        linewidth=2.0,
        label=f"Fit (Padé [{order}])",
        zorder=2,
    )

    # 数据点
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
        label=r"Data ($q_3:\ 64^3\times 16$)",
        zorder=3,
    )

    # 图例
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # 信息框
    info_text = (
        rf"$T_c = {t_peak_mean:.2f} \pm {t_peak_err:.2f}\ \mathrm{{MeV}}$" + "\n"
        rf"$\chi^2/\mathrm{{dof}} = {chi2_dof:.3f}$"
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
    order_slug = order.replace("/", "")
    ax.set_title(rf"Chiral Susceptibility Padé [{order}] Fit [$q_3:\ 64^3\times 16$]")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")

    plt.tight_layout()
    out_png = out_dir / f"chiral_susceptibility_pade{order_slug}_q3.png"
    out_pdf = out_dir / f"chiral_susceptibility_pade{order_slug}_q3.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def plot_combined_q3_pade(
    results: dict[str, dict],
    out_dir: Path,
    t0: float = 157.0,
) -> tuple[Path, Path]:
    """生成 Padé [0/2] 与 [1/2] 在同一图中的对比图"""
    apply_minimalist_style()
    fig, ax = plt.subplots(figsize=(8.5, 5.8), dpi=300)

    data_02 = results["0/2"]
    data_12 = results["1/2"]

    temps_data = data_02["temps"]
    chi_data = data_02["chi"]
    err_data = data_02["errors"]

    t_grid = np.linspace(np.min(temps_data) - 4.0, np.max(temps_data) + 4.0, 500)

    # 1. Padé [0/2]
    model_02 = PadeModel(order="0/2", t0=t0)
    recs_02 = [r for r in data_02["records"] if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    p_mean_02 = {pk: float(np.mean([float(r[pk]) for r in recs_02 if pk in r and np.isfinite(r[pk])])) for pk in model_02.param_keys}
    curve_02 = model_02.eval(t_grid, p_mean_02)

    mat_02 = np.array([model_02.eval(t_grid, {pk: float(r[pk]) for pk in model_02.param_keys if pk in r}) for r in recs_02])
    err_02 = np.sqrt(len(recs_02) - 1) * np.std(mat_02, axis=0, ddof=0)

    tpk_02 = float(data_02["summary"]["t_peak_mean"])
    terr_02 = float(data_02["summary"]["t_peak_err"])
    chi2_02 = float(data_02["summary"]["chi2_dof_mean"])

    ax.fill_between(t_grid, curve_02 - err_02, curve_02 + err_02, color="#1f77b4", alpha=0.20, label="Padé [0/2] Fit range")
    ax.plot(t_grid, curve_02, color="#1f77b4", linewidth=2.0, label=rf"Padé [0/2] ($T_c={tpk_02:.2f}\pm{terr_02:.2f}\,\mathrm{{MeV}},\ \chi^2/\mathrm{{dof}}={chi2_02:.2f}$)")

    # 2. Padé [1/2]
    model_12 = PadeModel(order="1/2", t0=t0)
    recs_12 = [r for r in data_12["records"] if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    p_mean_12 = {pk: float(np.mean([float(r[pk]) for r in recs_12 if pk in r and np.isfinite(r[pk])])) for pk in model_12.param_keys}
    curve_12 = model_12.eval(t_grid, p_mean_12)

    mat_12 = np.array([model_12.eval(t_grid, {pk: float(r[pk]) for pk in model_12.param_keys if pk in r}) for r in recs_12])
    err_12 = np.sqrt(len(recs_12) - 1) * np.std(mat_12, axis=0, ddof=0)

    tpk_12 = float(data_12["summary"]["t_peak_mean"])
    terr_12 = float(data_12["summary"]["t_peak_err"])
    chi2_12 = float(data_12["summary"]["chi2_dof_mean"])

    ax.fill_between(t_grid, curve_12 - err_12, curve_12 + err_12, color="#d62728", alpha=0.18, label="Padé [1/2] Fit range")
    ax.plot(t_grid, curve_12, color="#d62728", linestyle="--", linewidth=2.0, label=rf"Padé [1/2] ($T_c={tpk_12:.2f}\pm{terr_12:.2f}\,\mathrm{{MeV}},\ \chi^2/\mathrm{{dof}}={chi2_12:.2f}$)")

    # 3. 原数据点
    ax.errorbar(
        temps_data,
        chi_data,
        yerr=err_data,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.4,
        capsize=3.5,
        markersize=6.0,
        label=r"Data ($q_3:\ 64^3\times 16$)",
        zorder=5,
    )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax.set_title(r"Chiral Susceptibility Fit Comparison on $q_3$ ($64^3\times 16$)")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = out_dir / "chiral_susceptibility_pade02_12_q3.png"
    out_pdf = out_dir / "chiral_susceptibility_pade02_12_q3.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def clean_and_replot_q2_condensate() -> tuple[Path, Path]:
    """严格清理 q2 手征凝聚数据 (清理混入的 L32T12 杂质，保留真实 8 点) 并重新绘图"""
    cond_all_csv = OUTPUT_CONDENSATE_DIR / "all_ensembles_condensate.csv"
    df_c_all = pl.read_csv(cond_all_csv)

    # 仅保留真实 L48T16 开头的 8 组基准系综
    df_q2 = (
        df_c_all.filter(
            (pl.col("ns") == 48)
            & (pl.col("nt") == 16)
            & (pl.col("dataset_name").str.starts_with("L48T16"))
        )
        .select([
            pl.lit("q2").alias("dataset"),
            pl.col("ns"),
            pl.col("nt"),
            pl.col("beta"),
            pl.col("temperature").alias("temp"),
            (pl.col("temperature") * 16.0).alias("a_inv"),
            pl.col("ml"),
            pl.col("ms"),
            pl.col("mres"),
            pl.col("zm"),
            pl.col("num_cfgs"),
            pl.col("pbpl"),
            pl.col("pbpl_err"),
            pl.col("pbps"),
            pl.col("pbps_err"),
            pl.col("pbp_rm").alias("pbp_sub"),
            pl.col("pbp_rm_err").alias("pbp_sub_err"),
            (pl.col("pbp_rm") / pl.col("zm")).alias("pbp_renorm"),
            (pl.col("pbp_rm_err") / pl.col("zm")).alias("pbp_renorm_err"),
        ])
        .sort("temp")
    )

    # 写回各目录
    dirs = [
        PROJECT_ROOT / "22pade" / "q2",
        OUTPUT_LCP_DIR / "22pade" / "q2",
        OUTPUT_LCP_DIR / "q2",
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        df_q2.write_csv(d / "condensate_q2.csv")

    apply_minimalist_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    temps = df_q2["temp"].to_numpy().astype(float)
    pbp_ren = df_q2["pbp_renorm"].to_numpy().astype(float) * 1e4
    pbp_err = df_q2["pbp_renorm_err"].to_numpy().astype(float) * 1e4

    ax.errorbar(
        temps,
        pbp_ren,
        yerr=pbp_err,
        fmt="o-",
        color="black",
        ecolor="black",
        elinewidth=1.4,
        capsize=3.5,
        markersize=6.0,
        linewidth=1.3,
        label=r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m$ ($q_2:\ 48^3\times 16$, 8 Ensembles)",
        zorder=3,
    )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m \times 10^4$")
    ax.set_title(r"Renormalized Chiral Condensate [$q_2:\ 48^3\times 16$]")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = PROJECT_ROOT / "22pade" / "q2" / "chiral_condensate_q2.png"
    out_pdf = PROJECT_ROOT / "22pade" / "q2" / "chiral_condensate_q2.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    # 复制到各镜像目录
    for d in [OUTPUT_LCP_DIR / "22pade" / "q2", OUTPUT_LCP_DIR / "q2", FIGURES_DIR]:
        d.mkdir(parents=True, exist_ok=True)
        import shutil
        shutil.copyfile(out_png, d / "chiral_condensate_q2.png")
        shutil.copyfile(out_pdf, d / "chiral_condensate_q2.pdf")

    return out_png, out_pdf


def main() -> None:
    t0 = 157.0
    print("\n" + "=" * 80)
    print(">>> 启动 q3 专属 Padé [0/2] 与 [1/2] Jackknife 拟合及绘图任务...")
    print("=" * 80)

    # 1. 执行 q3 拟合
    results = run_q3_pade_fits(t0=t0)

    q3_dirs = [
        PROJECT_ROOT / "22pade" / "q3",
        OUTPUT_LCP_DIR / "22pade" / "q3",
        OUTPUT_LCP_DIR / "q3",
    ]

    for d in q3_dirs:
        d.mkdir(parents=True, exist_ok=True)
        # 保存汇总与样本
        pl.DataFrame([results["0/2"]["summary"]]).write_csv(d / "pade02_fit_summary_q3.csv")
        pl.DataFrame(results["0/2"]["records"]).write_csv(d / "pade02_jk_samples_q3.csv")
        pl.DataFrame([results["1/2"]["summary"]]).write_csv(d / "pade12_fit_summary_q3.csv")
        pl.DataFrame(results["1/2"]["records"]).write_csv(d / "pade12_jk_samples_q3.csv")

    # 2. 绘制 q3 独立拟合图件
    out_dir = PROJECT_ROOT / "22pade" / "q3"
    png_02, pdf_02 = plot_single_q3_pade("0/2", results["0/2"], out_dir, t0=t0)
    print(f"  [PLOT] q3 Padé [0/2] 阴影拟合图已生成: {png_02.name}")

    png_12, pdf_12 = plot_single_q3_pade("1/2", results["1/2"], out_dir, t0=t0)
    print(f"  [PLOT] q3 Padé [1/2] 阴影拟合图已生成: {png_12.name}")

    png_comp, pdf_comp = plot_combined_q3_pade(results, out_dir, t0=t0)
    print(f"  [PLOT] q3 Padé [0/2] vs [1/2] 对比图已生成: {png_comp.name}")

    # 同步图件到各镜像目录
    for d in [OUTPUT_LCP_DIR / "22pade" / "q3", OUTPUT_LCP_DIR / "q3", FIGURES_DIR]:
        d.mkdir(parents=True, exist_ok=True)
        import shutil
        for p in [png_02, pdf_02, png_12, pdf_12, png_comp, pdf_comp]:
            shutil.copyfile(p, d / p.name)

    # 3. 严格清洗 q2 手征凝聚数据并重新制图
    print("\n>>> 正在清洗并重新生成 q2 手征凝聚数据与图谱 (确保无 L32T12 混杂，严格 8 点)...")
    q2_png, q2_pdf = clean_and_replot_q2_condensate()
    print(f"  [OK] q2 手征凝聚图已更新 (精确 8 点): {q2_png.name}")

    # 4. 重新同步有限体积对比大图
    from scripts.tasks.lcp.plot_finite_volume_comparison import main as plot_fv_main
    print("\n>>> 正在同步更新有限体积整合大图 (已基于清洗后的 q2 8点数据)...")
    plot_fv_main()

    print("\n" + "=" * 80)
    print(">>> 所有 q3 单独拟合与 q2 纯净数据重构已圆满完成！")
    print("=" * 80)


if __name__ == "__main__":
    main()
