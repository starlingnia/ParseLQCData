#!/usr/bin/env python3
"""
scripts/tasks/lcp/generate_custom_pade_comparisons.py
--------------------------------------------------------------------------------
1. 生成定制化 4 系综对比图件：
   (A) p2 [2/2] + p3 [0/2] + q2 [2/2] + q3 [0/2]
   (B) p2 [2/2] + p3 [2/2] + q2 [2/2] + q3 [0/2]
   - 左图: 4 系综手征凝聚随温度变化
   - 右图: 4 系综手征易感度 Padé 拟合实线、Jackknife 阴影包络与数据点
2. 为 p2, p3, q2, q3 各自生成多阶 Padé 拟合对比图 (Multi-order Padé Fits):
   - 包含数据点与各阶 Padé 曲线，标注 Tc 与 chi2/dof
3. 保存所有图件至 22pade/ 及其子目录，并同步至 docs/figures/。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import shutil
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import FIGURES_DIR, OUTPUT_LCP_DIR
from src.parselqcdata.pade_fitter import (
    PadeModel,
    fit_jackknife_susceptibility_pade,
    fit_single_sample_pade,
)


def apply_minimalist_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10.0,
        "mathtext.fontset": "cm",
    })


def ensure_p3_02_jackknife(base_dir: Path, t0: float = 157.0) -> dict:
    """确保 p3 [0/2] 的 Jackknife 拟合结果已生成并缓存"""
    p3_dir = base_dir / "p3"
    sum_path = p3_dir / "pade02_fit_summary_p3.csv"
    jk_path = p3_dir / "pade02_jk_samples_p3.csv"

    if sum_path.exists() and jk_path.exists():
        df_sum = pl.read_csv(sum_path)
        df_jk = pl.read_csv(jk_path)
        return {"summary": df_sum.to_dicts()[0], "records": df_jk.to_dicts()}

    print("[INFO] 正在为 p3 执行 Padé [0/2] Jackknife 拟合...")
    df_s = pl.read_csv(p3_dir / "susceptibility_p3.csv").sort("temp")
    df_jk_raw = pl.read_csv(OUTPUT_LCP_DIR / "p3" / "susceptibility_jk_samples_p3.csv")

    temps = df_s["temp"].to_numpy().astype(float)
    errs = df_s["chi_renorm_err"].to_numpy().astype(float)
    betas = df_s["beta"].to_numpy().astype(float)

    samples = []
    for b in betas:
        s = df_jk_raw.filter(pl.col("beta") == b).sort("jk_index")["chi_renorm"].to_numpy().astype(float)
        samples.append(s)
    min_l = min(len(s) for s in samples)
    chi_mat = np.array([s[:min_l] for s in samples], dtype=float)

    summary, records = fit_jackknife_susceptibility_pade(temps, chi_mat, errs, order="0/2", t0=t0)
    pl.DataFrame([summary]).write_csv(sum_path)
    pl.DataFrame(records).write_csv(jk_path)
    print(f"  -> p3 [0/2] Tc = {summary['t_peak_mean']:.2f} +/- {summary['t_peak_err']:.2f} MeV")

    return {"summary": summary, "records": records}


def load_ensemble_data(
    ds_name: str,
    order: str,
    base_dir: Path,
) -> dict:
    """加载指定系综和指定 Padé 阶数的完整数据"""
    ds_dir = base_dir / ds_name
    df_s = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
    df_c = pl.read_csv(ds_dir / f"condensate_{ds_name}.csv").sort("temp")

    order_tag = order.replace("/", "")
    sum_file = ds_dir / f"pade{order_tag}_fit_summary_{ds_name}.csv"
    jk_file = ds_dir / f"pade{order_tag}_jk_samples_{ds_name}.csv"

    if not sum_file.exists() or not jk_file.exists():
        raise FileNotFoundError(f"Missing fit summary or jk samples for {ds_name} [{order}]: {sum_file}")

    df_sum = pl.read_csv(sum_file)
    df_jk = pl.read_csv(jk_file).filter(pl.col("ok") == True)

    return {
        "df_susc": df_s,
        "df_cond": df_c,
        "summary": df_sum.to_dicts()[0],
        "records": df_jk.to_dicts(),
        "order": order,
    }


def compute_pade_curve_and_band(
    order: str,
    records: list[dict],
    t_grid: np.ndarray,
    t0: float = 157.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """根据 Jackknife 样本计算实线均值拟合与阴影包络误差带"""
    model = PadeModel(order=order, t0=t0)
    p_keys = model.param_keys

    valid_recs = [r for r in records if r.get("ok", False) and np.isfinite(r.get("t_peak", np.nan))]
    n_samples = len(valid_recs)

    p_mean = {pk: float(np.mean([float(r[pk]) for r in valid_recs if pk in r and np.isfinite(r[pk])])) for pk in p_keys}
    curve_mean = model.eval(t_grid, p_mean)

    curves_matrix = np.zeros((n_samples, len(t_grid)), dtype=float)
    for i, r in enumerate(valid_recs):
        p_i = {pk: float(r[pk]) for pk in p_keys if pk in r}
        curves_matrix[i, :] = model.eval(t_grid, p_i)

    curve_err = np.sqrt(n_samples - 1) * np.std(curves_matrix, axis=0, ddof=0)
    y_lower = curve_mean - curve_err
    y_upper = curve_mean + curve_err

    return curve_mean, y_lower, y_upper


def plot_combined_4_ensembles(
    orders_map: dict[str, str],
    title_suffix: str,
    filename_prefix: str,
    base_dir: Path,
    out_dir: Path,
    t0: float = 157.0,
) -> tuple[Path, Path]:
    """
    生成 4 个系综的综合对比图 (Condensate + Susceptibility)
    orders_map: {"p2": "2/2", "p3": "0/2", "q2": "2/2", "q3": "0/2"}
    """
    apply_minimalist_style()

    styles = {
        "p2": {"color": "#1f77b4", "marker": "o", "label": r"$p_2\ (36^3\times 12,\ LT=3)$"},
        "p3": {"color": "#2ca02c", "marker": "s", "label": r"$p_3\ (48^3\times 12,\ LT=4)$"},
        "q2": {"color": "#d62728", "marker": "^", "label": r"$q_2\ (48^3\times 16,\ LT=3)$"},
        "q3": {"color": "#9467bd", "marker": "D", "label": r"$q_3\ (64^3\times 16,\ LT=4)$"},
    }

    datasets_data = {}
    for ds_name, ord_val in orders_map.items():
        datasets_data[ds_name] = load_ensemble_data(ds_name, ord_val, base_dir)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.0), dpi=300)

    # --- 左图: 手征凝聚 ---
    for ds_name in ["p2", "p3", "q2", "q3"]:
        d = datasets_data[ds_name]
        st = styles[ds_name]
        temps = d["df_cond"]["temp"].to_numpy().astype(float)
        pbp = d["df_cond"]["pbp_renorm"].to_numpy().astype(float) * 1e4
        err = d["df_cond"]["pbp_renorm_err"].to_numpy().astype(float) * 1e4

        ax1.errorbar(
            temps, pbp, yerr=err,
            fmt=f"{st['marker']}-", color=st["color"], label=st["label"],
            capsize=3.5, markersize=5.5, linewidth=1.3, alpha=0.9,
        )

    ax1.set_xlabel(r"Temperature $T$ [MeV]")
    ax1.set_ylabel(r"$\langle \bar{\psi}\psi \rangle_{\rm sub} / Z_m \times 10^4$")
    ax1.set_title(r"Chiral Condensate Comparison ($p_2, p_3, q_2, q_3$)")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # --- 右图: 手征易感度 Padé 拟合 ---
    all_temps = np.concatenate([d["df_susc"]["temp"].to_numpy().astype(float) for d in datasets_data.values()])
    t_grid = np.linspace(np.min(all_temps) - 3.0, np.max(all_temps) + 3.0, 500)

    for ds_name in ["p2", "p3", "q2", "q3"]:
        d = datasets_data[ds_name]
        st = styles[ds_name]
        ord_val = orders_map[ds_name]

        # 数据点
        ax2.errorbar(
            d["df_susc"]["temp"].to_numpy().astype(float),
            d["df_susc"]["chi_renorm"].to_numpy().astype(float),
            yerr=d["df_susc"]["chi_renorm_err"].to_numpy().astype(float),
            fmt=st["marker"], color=st["color"], capsize=3.5, markersize=5.5, alpha=0.75,
            label=f"{st['label']} Data",
        )

        # Padé 拟合曲线与阴影
        c_curve, low_b, up_b = compute_pade_curve_and_band(ord_val, d["records"], t_grid, t0=t0)
        tpk = float(d["summary"]["t_peak_mean"])
        terr = float(d["summary"]["t_peak_err"])

        ax2.fill_between(t_grid, low_b, up_b, color=st["color"], alpha=0.15)
        ax2.plot(
            t_grid, c_curve, color=st["color"], linewidth=1.8,
            label=rf"{st['label']} [{ord_val}] ($T_c={tpk:.1f}\pm{terr:.1f}\,$MeV)",
        )

    ax2.set_xlabel(r"Temperature $T$ [MeV]")
    ax2.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax2.set_title(f"Chiral Susceptibility Comparison {title_suffix}")
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = out_dir / f"{filename_prefix}.png"
    out_pdf = out_dir / f"{filename_prefix}.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def plot_multimodal_pade_per_ensemble(
    ds_name: str,
    base_dir: Path,
    out_dir: Path,
    t0: float = 157.0,
) -> tuple[Path, Path]:
    """为单个系综绘制不同阶 Padé 拟合的同图对比 (Multi-order Padé Comparison)"""
    apply_minimalist_style()

    ds_dir = base_dir / ds_name
    df_s = pl.read_csv(ds_dir / f"susceptibility_{ds_name}.csv").sort("temp")
    temps = df_s["temp"].to_numpy().astype(float)
    chi_vals = df_s["chi_renorm"].to_numpy().astype(float)
    chi_errs = df_s["chi_renorm_err"].to_numpy().astype(float)

    # 确定测试的阶数
    if ds_name == "q3":
        orders = ["0/2", "1/2", "2/2"]
    else:
        orders = ["0/2", "1/2", "2/2", "3/2", "1/3"]

    color_cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    linestyle_cycle = ["-", "--", "-.", ":", "-"]

    fig, ax = plt.subplots(figsize=(8.8, 6.0), dpi=300)

    # 绘制数据点
    ax.errorbar(
        temps, chi_vals, yerr=chi_errs,
        fmt="ko", capsize=4.0, markersize=6.0, zorder=5, label="Lattice Data",
    )

    t_grid = np.linspace(np.min(temps) - 4.0, np.max(temps) + 4.0, 500)

    # 逐阶拟合
    for i, order in enumerate(orders):
        model = PadeModel(order=order, t0=t0)
        fit_res = fit_single_sample_pade(temps, chi_vals, chi_errs, model)
        if not fit_res or not fit_res["ok"]:
            continue

        p_fit = fit_res["params"]
        if not p_fit or not all(k in p_fit for k in model.param_keys):
            continue

        y_fit = model.eval(t_grid, p_fit)

        tpk = float(fit_res["t_peak"])
        chi2 = float(fit_res["chi2"])
        dof = float(fit_res["dof"])
        chi2_dof = chi2 / dof if dof > 0 else np.nan

        if not np.isfinite(tpk):
            continue

        lbl = rf"Padé [{order}] ($T_c={tpk:.1f}\,$MeV"
        if np.isfinite(chi2_dof):
            lbl += rf", $\chi^2/\mathrm{{dof}}={chi2_dof:.2f}$)"
        else:
            lbl += r", $\mathrm{dof}=0$)"

        col = color_cycle[i % len(color_cycle)]
        ls = linestyle_cycle[i % len(linestyle_cycle)]

        ax.plot(t_grid, y_fit, color=col, linestyle=ls, linewidth=2.0, label=lbl, zorder=3)

    ns_nt_map = {
        "p2": r"p_2\ (36^3 \times 12,\ N_s/N_t = 3)",
        "p3": r"p_3\ (48^3 \times 12,\ N_s/N_t = 4)",
        "q2": r"q_2\ (48^3 \times 16,\ N_s/N_t = 3)",
        "q3": r"q_3\ (64^3 \times 16,\ N_s/N_t = 4)",
    }

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"$\chi_{\mathrm{disc}} / Z_m^2\ [\mathrm{GeV}^2]$")
    ax.set_title(rf"Chiral Susceptibility Multi-Order Padé Fits: ${ns_nt_map[ds_name]}$")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    # y 轴自适应范围，避免极点发散破坏画幅
    y_max = np.max(chi_vals) * 1.55
    ax.set_ylim(-0.05, y_max)

    plt.tight_layout()
    out_png = out_dir / f"chiral_susceptibility_multimodal_{ds_name}.png"
    out_pdf = out_dir / f"chiral_susceptibility_multimodal_{ds_name}.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return out_png, out_pdf


def main() -> None:
    base_dir = PROJECT_ROOT / "22pade"
    out_22pade = PROJECT_ROOT / "22pade"
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # 1. 确保 p3 [0/2] Jackknife 缓存存在
    ensure_p3_02_jackknife(base_dir)

    # 2. 生成图件 1: 2/2 + q3 0/2 + p3 0/2 (p2 [2/2], p3 [0/2], q2 [2/2], q3 [0/2])
    print("[INFO] 生成图件 1: 2/2 + q3 0/2 + p3 0/2 对比图...")
    f1_orders = {"p2": "2/2", "p3": "0/2", "q2": "2/2", "q3": "0/2"}
    f1_png, f1_pdf = plot_combined_4_ensembles(
        f1_orders,
        title_suffix=r"($p_2[2/2], p_3[0/2], q_2[2/2], q_3[0/2]$)",
        filename_prefix="chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02",
        base_dir=base_dir,
        out_dir=out_22pade,
    )
    shutil.copyfile(f1_png, FIGURES_DIR / f1_png.name)
    shutil.copyfile(f1_pdf, FIGURES_DIR / f1_pdf.name)
    print(f"  -> 生成: {f1_png.name}")

    # 3. 生成图件 2: 其他都是 2/2 就 q3 0/2 (p2 [2/2], p3 [2/2], q2 [2/2], q3 [0/2])
    print("[INFO] 生成图件 2: 其他都是 2/2 就 q3 0/2 对比图...")
    f2_orders = {"p2": "2/2", "p3": "2/2", "q2": "2/2", "q3": "0/2"}
    f2_png, f2_pdf = plot_combined_4_ensembles(
        f2_orders,
        title_suffix=r"($p_2[2/2], p_3[2/2], q_2[2/2], q_3[0/2]$)",
        filename_prefix="chiral_lcp_comparison_p2_22_p3_22_q2_22_q3_02",
        base_dir=base_dir,
        out_dir=out_22pade,
    )
    shutil.copyfile(f2_png, FIGURES_DIR / f2_png.name)
    shutil.copyfile(f2_pdf, FIGURES_DIR / f2_pdf.name)
    print(f"  -> 生成: {f2_png.name}")

    # 4. 生成各个系综的多阶拟合对比图
    print("[INFO] 生成各系综多阶 Padé 拟合对比图...")
    for ds_name in ["p2", "p3", "q2", "q3"]:
        ds_dir = base_dir / ds_name
        p_png, p_pdf = plot_multimodal_pade_per_ensemble(ds_name, base_dir, ds_dir)
        shutil.copyfile(p_png, FIGURES_DIR / p_png.name)
        shutil.copyfile(p_pdf, FIGURES_DIR / p_pdf.name)
        # 也复制一份到 22pade/ 根目录方便查看
        shutil.copyfile(p_png, out_22pade / p_png.name)
        shutil.copyfile(p_pdf, out_22pade / p_pdf.name)
        print(f"  -> 生成 {ds_name}: {p_png.name}")

    print("[SUCCESS] 所有指定图件生成并同步完成！")


if __name__ == "__main__":
    main()
