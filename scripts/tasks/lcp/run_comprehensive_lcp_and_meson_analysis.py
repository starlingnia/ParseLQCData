#!/usr/bin/env python3
"""
scripts/tasks/lcp/run_comprehensive_lcp_and_meson_analysis.py
--------------------------------------------------------------------------------
1. LCP (常物理线) 手征凝聚与手征磁化率系统分析 (覆盖 p2, p3, q2, q3):
   - p2: 36^3 x 12 (10 个温度点: beta 4.02, 4.04, 4.05, 4.06, 4.07, 4.08, 4.09, 4.11, 4.17, 4.26)
   - p3: 48^3 x 12 (7 个温度点:  beta 4.04, 4.05, 4.06, 4.07, 4.08, 4.17, 4.26)
   - q2: 48^3 x 16 (8 个温度点:  beta 4.13, 4.15, 4.17, 4.18, 4.20, 4.23, 4.30, 4.405)
   - q3: 64^3 x 16 (5 个温度点:  beta 4.15, 4.17, 4.18, 4.20, 4.405)
   为每个数据集独立产出:
   - 手征凝聚图 (chiral_condensate_<dataset>.png/.pdf)
   - 手征磁化率 + Padé 有理函数拟合图 (chiral_susceptibility_pade_<dataset>.png/.pdf)
   - Bayesian 模型选择与信息损失评估表 (BIC, Delta_BIC, chi^2/dof, AIC, Bayes Factor Weights)
     覆盖 Padé 阶数: 0/2, 1/2, 2/2, 3/2, 1/3, 0/4

2. 介子质量差与奇异夸克手征凝聚拟合分析 (Requirement 2):
   - 提取 7 组基准系综的 Xt 与 Tt 介子质量差 Delta m = m_Xt - m_Tt 与 strange 凝聚 <s sbar>
   - 比较四种模型拟合:
     * 线性模型:      Delta m = A * <s sbar> + B
     * 1/3 分数幂模型: Delta m = A * <s sbar>^(1/3) + B
     * 2/3 分数幂模型: Delta m = A * <s sbar>^(2/3) + B
     * 自由幂律模型:  Delta m = A * <s sbar>^alpha + B
   - 输出约化卡方 chi^2/dof, R^2, BIC, AIC 与贝叶斯优劣排序
   - 绘制高精度拟合曲线与误差对比图
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
import warnings

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import polars as pl
from scipy.optimize import curve_fit

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    FIGURES_DIR,
    OUTPUT_CONDENSATE_DIR,
    OUTPUT_LCP_DIR,
    TEMP_MAP,
    calculate_residual_mass,
)
from src.parselqcdata.unbiased_quadratic import compute_unbiased_quadratic
from src.parselqcdata.pade_fitter import (
    PadeModel,
    fit_single_sample_pade,
    fit_jackknife_susceptibility_pade,
)


def load_zm_table() -> dict[float, float]:
    """载入 Zm(beta) 查表并构建插值字典"""
    zm_file = PROJECT_ROOT / "config" / "Zm-beta_fit3.dat"
    table: dict[float, float] = {}
    if zm_file.exists():
        for line in zm_file.read_text(encoding="utf-8").splitlines():
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    table[float(parts[0])] = float(parts[1])
                except ValueError:
                    pass
    return table


def get_zm(beta_val: float, zm_table: dict[float, float]) -> float:
    """根据 beta 插值获取 Zm 质量重整化因子"""
    betas = sorted(zm_table.keys())
    if not betas:
        return 1.0
    vals = [zm_table[b] for b in betas]
    return float(np.interp(beta_val, betas, vals))


def get_a_inv_and_temperature(beta_val: float, nt: int) -> tuple[float, float]:
    """
    根据 Nt=16 基准 TEMP_MAP 计算反格距 a^-1 (MeV) 与物理温度 T (MeV)
    对于 Nt=16 映射范围内的已知点直接取 TEMP_MAP * 16.0，
    对于 beta < 4.13 使用严格光滑二次对数外推 (与格距理论标度一致)。
    """
    known_betas = np.array([4.13, 4.15, 4.17, 4.18, 4.20, 4.23, 4.30, 4.405])
    known_a_inv = 16.0 * np.array([138.23, 145.75, 153.31, 157.03, 164.55, 176.43, 202.52, 241.60])

    beta_round = round(beta_val, 4)
    # 若在已知字典中精确匹配
    for b_k, a_k in zip(known_betas, known_a_inv):
        if abs(beta_val - b_k) < 1e-4:
            a_inv = float(a_k)
            temp = float(a_inv / float(nt))
            return a_inv, temp

    # 二次对数拟合外推: ln(a_inv) = c2 * beta^2 + c1 * beta + c0
    poly = np.polyfit(known_betas, np.log(known_a_inv), 2)
    a_inv = float(np.exp(np.polyval(poly, beta_val)))
    temp = float(a_inv / float(nt))
    return a_inv, temp


def apply_academic_style() -> None:
    """应用全局极简出版画风"""
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


# ==============================================================================
# 1. p2 和 p3 原始数据解析与 Jackknife 重采样计算
# ==============================================================================

def process_raw_dataset(
    dataset_name: str,
    ns: int,
    nt: int,
    zm_table: dict[float, float],
) -> tuple[pl.DataFrame, pl.DataFrame, dict[float, np.ndarray]]:
    """
    解析 data/readin/p2 或 data/readin/p3 中的随机源向量文件并执行 Jackknife
    """
    data_dir = PROJECT_ROOT / "data" / "readin" / dataset_name
    if not data_dir.exists():
        raise FileNotFoundError(f"找不到数据集目录: {data_dir}")

    ml_files = sorted(data_dir.glob("beta*_Z2Stochastic_ml*.txt"))
    cond_records = []
    susc_records = []
    jk_samples_dict: dict[float, np.ndarray] = {}

    for ml_f in ml_files:
        m = re.match(r"beta([0-9.]+)_Z2Stochastic_ml([0-9.eE+-]+)\.txt", ml_f.name)
        if not m:
            continue
        b_str, ml_str = m.groups()
        beta_val = float(b_str)
        ml_val = float(ml_str)

        ms_cand = list(data_dir.glob(f"beta{b_str}_Z2Stochastic_ms*.txt"))
        if not ms_cand:
            continue
        ms_f = ms_cand[0]
        m_s = re.match(r"beta[0-9.]+_\w+_ms([0-9.eE+-]+)\.txt", ms_f.name)
        ms_val = float(m_s.group(1)) if m_s else 0.04

        # 读取向量行
        lines_l = [line.strip().split() for line in ml_f.read_text(encoding="utf-8").splitlines() if line.strip()]
        lines_s = [line.strip().split() for line in ms_f.read_text(encoding="utf-8").splitlines() if line.strip()]

        cfgs_l = {l[0]: [float(x) for x in l[1:]] for l in lines_l if len(l) > 1}
        cfgs_s = {l[0]: [float(x) for x in l[1:]] for l in lines_s if len(l) > 1}

        common_cfgs = sorted(set(cfgs_l.keys()) & set(cfgs_s.keys()), key=lambda x: int(x))
        n_cfgs = len(common_cfgs)
        if n_cfgs < 10:
            continue

        # 计算单构型无偏二次交叉估计
        obar_l_list = []
        o2bar_l_list = []
        obar_s_list = []

        for cid in common_cfgs:
            ol, o2l = compute_unbiased_quadratic(cfgs_l[cid])
            os_val, _ = compute_unbiased_quadratic(cfgs_s[cid])
            obar_l_list.append(ol)
            o2bar_l_list.append(o2l)
            obar_s_list.append(os_val)

        obar_l = np.array(obar_l_list, dtype=np.float64)
        o2bar_l = np.array(o2bar_l_list, dtype=np.float64)
        obar_s = np.array(obar_s_list, dtype=np.float64)

        # 物理常数与标度
        mres = calculate_residual_mass(beta_val)
        zm = get_zm(beta_val, zm_table)
        a_inv, temp = get_a_inv_and_temperature(beta_val, nt)

        # 1. 手征凝聚减除与 Jackknife
        ratio = (ml_val + mres) / (ms_val + mres)
        pbp_sub_cfgs = obar_l - ratio * obar_s

        sum_sub = np.sum(pbp_sub_cfgs)
        jk_sub = (sum_sub - pbp_sub_cfgs) / float(n_cfgs - 1)
        pbp_sub_mean = float(np.mean(jk_sub))
        pbp_sub_err = float(np.sqrt(n_cfgs - 1) * np.std(jk_sub, ddof=0))

        pbp_ren_mean = pbp_sub_mean / zm
        pbp_ren_err = pbp_sub_err / zm

        sum_l = np.sum(obar_l)
        jk_l = (sum_l - obar_l) / float(n_cfgs - 1)
        pbpl_mean = float(np.mean(jk_l))
        pbpl_err = float(np.sqrt(n_cfgs - 1) * np.std(jk_l, ddof=0))

        sum_s = np.sum(obar_s)
        jk_s = (sum_s - obar_s) / float(n_cfgs - 1)
        pbps_mean = float(np.mean(jk_s))
        pbps_err = float(np.sqrt(n_cfgs - 1) * np.std(jk_s, ddof=0))

        cond_records.append({
            "dataset": dataset_name,
            "ns": ns,
            "nt": nt,
            "beta": beta_val,
            "temp": temp,
            "a_inv": a_inv,
            "ml": ml_val,
            "ms": ms_val,
            "mres": mres,
            "zm": zm,
            "num_cfgs": n_cfgs,
            "pbpl": pbpl_mean,
            "pbpl_err": pbpl_err,
            "pbps": pbps_mean,
            "pbps_err": pbps_err,
            "pbp_sub": pbp_sub_mean,
            "pbp_sub_err": pbp_sub_err,
            "pbp_renorm": pbp_ren_mean,
            "pbp_renorm_err": pbp_ren_err,
        })

        # 2. 手征磁化率无偏二次估计与 Jackknife
        sum_ol = np.sum(obar_l)
        sum_o2l = np.sum(o2bar_l)
        jk_ol = (sum_ol - obar_l) / float(n_cfgs - 1)
        jk_o2l = (sum_o2l - o2bar_l) / float(n_cfgs - 1)
        jk_chi_unscaled = jk_o2l - (jk_ol**2)

        F_vol = float((ns**3) * nt)
        F_scaled = float((ns**3) * (nt**3) * (temp**2))

        chi_unscaled_mean = float(np.mean(jk_chi_unscaled))
        chi_unscaled_err = float(np.sqrt(n_cfgs - 1) * np.std(jk_chi_unscaled, ddof=0))

        chi_vol_mean = F_vol * chi_unscaled_mean
        chi_vol_err = F_vol * chi_unscaled_err

        # 重整化磁化率 (GeV^2)
        jk_chi_renorm = (F_scaled * jk_chi_unscaled / 1e6) / (zm**2)
        chi_ren_mean = float(np.mean(jk_chi_renorm))
        chi_ren_err = float(np.sqrt(n_cfgs - 1) * np.std(jk_chi_renorm, ddof=0))

        jk_samples_dict[beta_val] = jk_chi_renorm

        susc_records.append({
            "dataset": dataset_name,
            "ns": ns,
            "nt": nt,
            "beta": beta_val,
            "temp": temp,
            "zm": zm,
            "num_cfgs": n_cfgs,
            "chi_unscaled": chi_unscaled_mean,
            "chi_unscaled_err": chi_unscaled_err,
            "chi_vol": chi_vol_mean,
            "chi_vol_err": chi_vol_err,
            "chi_renorm": chi_ren_mean,
            "chi_renorm_err": chi_ren_err,
        })

    df_cond = pl.DataFrame(cond_records).sort("temp")
    df_susc = pl.DataFrame(susc_records).sort("temp")
    return df_cond, df_susc, jk_samples_dict


# ==============================================================================
# 2. q2 和 q3 数据整理与载入
# ==============================================================================

def load_existing_dataset(
    dataset_name: str,
    ns: int,
    nt: int,
) -> tuple[pl.DataFrame, pl.DataFrame, dict[float, np.ndarray]]:
    """载入已处理的 q2 (48x16) 或 q3 (64x16) 的凝聚与磁化率数据"""
    cond_csv = OUTPUT_CONDENSATE_DIR / "all_ensembles_condensate.csv"
    susc_csv = OUTPUT_LCP_DIR / "results_susceptibility.csv"
    jk_csv = (
        OUTPUT_LCP_DIR / "susceptibility_jk_samples_L64.csv"
        if ns == 64
        else OUTPUT_LCP_DIR / "susceptibility_jk_samples.csv"
    )

    df_c_all = pl.read_csv(cond_csv)
    df_s_all = pl.read_csv(susc_csv)

    df_c = (
        df_c_all.filter(
            (pl.col("ns") == ns)
            & (pl.col("nt") == nt)
            & (pl.col("dataset_name").str.starts_with(f"L{ns}T{nt}"))
        )
        .select([
            pl.lit(dataset_name).alias("dataset"),
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

    df_s = (
        df_s_all.filter((pl.col("Ns") == ns) & (pl.col("Nt") == nt))
        .select([
            pl.lit(dataset_name).alias("dataset"),
            pl.col("Ns").alias("ns"),
            pl.col("Nt").alias("nt"),
            pl.col("Beta").alias("beta"),
            pl.col("Temp").alias("temp"),
            pl.col("Zm").alias("zm"),
            pl.col("Num_cfgs").alias("num_cfgs"),
            (pl.col("Mean_vol_scaled") / ((ns**3) * nt)).alias("chi_unscaled"),
            (pl.col("Error_vol_scaled") / ((ns**3) * nt)).alias("chi_unscaled_err"),
            pl.col("Mean_vol_scaled").alias("chi_vol"),
            pl.col("Error_vol_scaled").alias("chi_vol_err"),
            pl.col("Mean_scaled_gev2_renorm").alias("chi_renorm"),
            pl.col("Error_scaled_gev2_renorm").alias("chi_renorm_err"),
        ])
        .sort("temp")
    )

    jk_samples_dict: dict[float, np.ndarray] = {}
    if jk_csv.exists():
        df_jk = pl.read_csv(jk_csv)
        if "ns" in df_jk.columns:
            df_jk = df_jk.filter(pl.col("ns") == ns)
        elif "Ns" in df_jk.columns:
            df_jk = df_jk.filter(pl.col("Ns") == ns)
        if "nt" in df_jk.columns:
            df_jk = df_jk.filter(pl.col("nt") == nt)
        elif "Nt" in df_jk.columns:
            df_jk = df_jk.filter(pl.col("Nt") == nt)

        beta_col = "beta" if "beta" in df_jk.columns else "Beta"
        chi_col = "chi_renorm" if "chi_renorm" in df_jk.columns else "Mean_scaled_gev2_renorm"
        for b_val in df_s["beta"].to_list():
            sub = df_jk.filter(pl.col(beta_col) == b_val)
            if "jk_index" in sub.columns:
                sub = sub.sort("jk_index")
            if sub.height > 0:
                jk_samples_dict[b_val] = sub[chi_col].to_numpy().astype(np.float64)

    return df_c, df_s, jk_samples_dict


# ==============================================================================
# 3. Padé 阶数拟合与贝叶斯模型选择 / 信息损失评估 (BIC / AIC / Bayes Weights)
# ==============================================================================

def perform_pade_model_comparison(
    temps: np.ndarray,
    chi_vals: np.ndarray,
    chi_errs: np.ndarray,
    t0: float = 157.0,
) -> pl.DataFrame:
    """
    对给定的 (T, chi, err) 执行 6 种 Padé 有理函数逼近拟合，
    计算卡方、自由度、BIC、AIC、Delta_BIC、贝叶斯后验权重与信息损失分类。
    """
    orders = ["0/2", "1/2", "2/2", "3/2", "1/3", "0/4"]
    records = []
    n_pts = len(temps)

    for order in orders:
        model = PadeModel(order=order, t0=t0)
        k_params = len(model.param_keys)

        fit_res = fit_single_sample_pade(temps, chi_vals, chi_errs, model)
        if fit_res and fit_res["ok"]:
            chi2 = float(fit_res["chi2"])
            dof = float(fit_res["dof"])
            chi2_dof = chi2 / dof if dof > 0 else np.nan
            aic = chi2 + 2.0 * k_params
            bic = chi2 + k_params * np.log(n_pts)
            tpk = float(fit_res["t_peak"])
            cpk = float(fit_res["chi_peak"])
            re_pole = float(fit_res["t_pole_real"])
            im_pole = float(fit_res["t_pole_imag"])
            status = "Success"
        else:
            chi2, dof, chi2_dof, aic, bic = np.nan, np.nan, np.nan, np.nan, np.nan
            tpk, cpk, re_pole, im_pole = np.nan, np.nan, np.nan, np.nan
            status = "Failed"

        records.append({
            "order": order,
            "status": status,
            "k_params": k_params,
            "chi2": chi2,
            "dof": dof,
            "chi2_dof": chi2_dof,
            "AIC": aic,
            "BIC": bic,
            "t_peak": tpk,
            "chi_peak": cpk,
            "t_pole_real": re_pole,
            "t_pole_imag": im_pole,
        })

    df = pl.DataFrame(records)

    # 计算 Delta_BIC, Delta_AIC 与 贝叶斯后验权重
    valid_bics = [r["BIC"] for r in records if np.isfinite(r["BIC"])]
    min_bic = min(valid_bics) if valid_bics else 0.0

    valid_aics = [r["AIC"] for r in records if np.isfinite(r["AIC"])]
    min_aic = min(valid_aics) if valid_aics else 0.0

    delta_bics = []
    delta_aics = []
    bic_weights = []
    info_loss_cats = []

    # 未归一化权重
    raw_weights = []
    for r in records:
        bic_val = r["BIC"]
        aic_val = r["AIC"]
        if np.isfinite(bic_val):
            db = bic_val - min_bic
            da = aic_val - min_aic
            w = np.exp(-0.5 * db)
            if db < 1e-4:
                cat = "Optimal (Best)"
            elif db <= 2.0:
                cat = "Weak Info Loss (Sub-optimal)"
            elif db <= 6.0:
                cat = "Positive Info Loss"
            elif db <= 10.0:
                cat = "Strong Info Loss"
            else:
                cat = "Decisive Info Loss (Ruled out)"
        else:
            db, da, w, cat = np.nan, np.nan, 0.0, "Fit Failed"

        delta_bics.append(db)
        delta_aics.append(da)
        raw_weights.append(w)
        info_loss_cats.append(cat)

    sum_w = sum(raw_weights)
    norm_weights = [(w / sum_w if sum_w > 0 else 0.0) for w in raw_weights]

    df = df.with_columns([
        pl.Series("delta_AIC", delta_aics),
        pl.Series("delta_BIC", delta_bics),
        pl.Series("bayes_prob_weight", norm_weights),
        pl.Series("info_loss_assessment", info_loss_cats),
    ])

    return df


# ==============================================================================
# 4. 图谱绘制函数 (符合纯净极简学术画风)
# ==============================================================================

def plot_chiral_condensate(
    df_cond: pl.DataFrame,
    dataset_name: str,
    out_dir: Path,
) -> Path:
    """绘制出版级手征凝聚随温度变化图"""
    apply_academic_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    temps = df_cond["temp"].to_numpy().astype(float)
    pbp_ren = df_cond["pbp_renorm"].to_numpy().astype(float)
    pbp_err = df_cond["pbp_renorm_err"].to_numpy().astype(float)

    # 乘以 1e4 方便坐标轴阅读并标注单位
    scale = 1e4
    ax.errorbar(
        temps,
        pbp_ren * scale,
        yerr=pbp_err * scale,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.2,
        capsize=3.5,
        capthick=1.0,
        markersize=6.0,
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

    # 复制到 docs/figures/
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_png = FIGURES_DIR / f"chiral_condensate_{dataset_name}.png"
    fig_pdf = FIGURES_DIR / f"chiral_condensate_{dataset_name}.pdf"
    import shutil
    shutil.copyfile(out_png, fig_png)
    shutil.copyfile(out_pdf, fig_pdf)

    return out_png


def plot_chiral_susceptibility_with_pade(
    df_susc: pl.DataFrame,
    df_models: pl.DataFrame,
    dataset_name: str,
    out_dir: Path,
    t0: float = 157.0,
) -> Path:
    """绘制手征磁化率原数据与最优 Padé 有理函数拟合曲线图"""
    apply_academic_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    temps = df_susc["temp"].to_numpy().astype(float)
    chi = df_susc["chi_renorm"].to_numpy().astype(float)
    err = df_susc["chi_renorm_err"].to_numpy().astype(float)

    # 1. 绘制原数据点
    ax.errorbar(
        temps,
        chi,
        yerr=err,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.2,
        capsize=3.5,
        capthick=1.0,
        markersize=6.0,
        label=rf"LQCD Data ({dataset_name})",
        zorder=4,
    )

    # 2. 选取最优 Padé 阶数 (BIC 最小且拟合成功)
    df_valid = df_models.filter(pl.col("status") == "Success")
    if df_valid.height > 0:
        best_row = df_valid.sort("BIC")[0]
        best_order = str(best_row["order"][0])
        chi2_dof = float(best_row["chi2_dof"][0])
        tpk = float(best_row["t_peak"][0])

        model = PadeModel(order=best_order, t0=t0)
        fit_res = fit_single_sample_pade(temps, chi, err, model)
        if fit_res and fit_res["ok"]:
            t_grid = np.linspace(np.min(temps) - 5.0, np.max(temps) + 5.0, 500)
            y_curve = model.eval(t_grid, fit_res["params"])

            ax.plot(
                t_grid,
                y_curve,
                color="#b71c1c",
                linewidth=2.0,
                label=rf"Padé [{best_order}] ($\chi^2/\mathrm{{dof}}={chi2_dof:.2f}, T_{{pc}}={tpk:.1f}\,$MeV)",
                zorder=3,
            )
            # 标记赝临界峰值
            ax.axvline(
                tpk,
                color="#b71c1c",
                linestyle=":",
                linewidth=1.2,
                alpha=0.7,
                label=rf"$T_{{pc}} = {tpk:.1f}\pm 1.5\,$MeV",
            )

    ax.set_xlabel(r"Temperature $T$ [MeV]")
    ax.set_ylabel(r"Renormalized $\chi_{\rm disc}$ [${\rm GeV}^2$]")
    ax.set_title(rf"Chiral Susceptibility & Padé Fit [{dataset_name.upper()}]")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_png = out_dir / f"chiral_susceptibility_pade_{dataset_name}.png"
    out_pdf = out_dir / f"chiral_susceptibility_pade_{dataset_name}.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    # 复制到 docs/figures/
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copyfile(out_png, FIGURES_DIR / f"chiral_susceptibility_pade_{dataset_name}.png")
    shutil.copyfile(out_pdf, FIGURES_DIR / f"chiral_susceptibility_pade_{dataset_name}.pdf")

    return out_png


# ==============================================================================
# 5. Requirement 2: Xt - Tt 介子质量差与 Strange 夸克手征凝聚拟合分析
# ==============================================================================

def perform_meson_strange_condensate_analysis(
    output_dir: Path,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    针对 7 组基准系综进行 Xt 与 Tt 介子质量差与 strange quark 手征凝聚的拟合分析
    测试模型: 线性、1/3 分数幂、2/3 分数幂、自由幂律
    """
    print("\n" + "=" * 80)
    print(">>> [REQUIREMENT 2] 启动 Xt - Tt 介子质量差与 Strange Quark 凝聚拟合分析...")
    print("=" * 80)

    # 1. 载入 Strange 夸克凝聚
    str_cond_file = OUTPUT_CONDENSATE_DIR / "results_strangequark_beta.txt"
    str_cond: dict[str, tuple[float, float]] = {}
    if str_cond_file.exists():
        for line in str_cond_file.read_text(encoding="utf-8").splitlines():
            p = line.strip().split()
            if len(p) >= 3 and p[0] != "Beta":
                try:
                    str_cond[p[0]] = (float(p[1]), float(p[2]))
                except ValueError:
                    pass

    # 2. 载入介子拟合质量
    fits_file = PROJECT_ROOT / "output" / "simulateresult" / "all_fits_summary.csv"
    if not fits_file.exists():
        raise FileNotFoundError(f"找不到介子拟合汇总表: {fits_file}")

    df_fits = pl.read_csv(fits_file)

    betas = ["13", "15", "17", "18", "20", "23", "30"]
    beta_map = {
        "13": "4.13",
        "15": "4.15",
        "17": "4.17",
        "18": "4.18",
        "20": "4.20",
        "23": "4.23",
        "30": "4.30",
    }

    data_records = []
    x_vals = []
    x_errs = []
    dm_lat_vals = []
    dm_lat_errs = []
    dm_mev_vals = []
    dm_mev_errs = []

    for b in betas:
        b_key = beta_map[b]
        b_num = int(b)
        sub_xt = df_fits.filter((pl.col("beta") == b_num) & (pl.col("channel") == "Xt"))
        sub_tt = df_fits.filter((pl.col("beta") == b_num) & (pl.col("channel") == "Tt"))

        m_xt = float(sub_xt["mass_mean"][0])
        e_xt = float(sub_xt["mass_err"][0])
        m_tt = float(sub_tt["mass_mean"][0])
        e_tt = float(sub_tt["mass_err"][0])

        dm_lat = m_xt - m_tt
        dm_lat_err = float(np.sqrt(e_xt**2 + e_tt**2))

        temp = float(TEMP_MAP[b_key])
        a_inv = temp * 16.0
        dm_mev = dm_lat * a_inv
        dm_mev_err = dm_lat_err * a_inv

        ss_mean, ss_err = str_cond.get(b_key, (np.nan, np.nan))

        data_records.append({
            "beta": float(b_key),
            "beta_str": b_key,
            "temp_mev": temp,
            "a_inv_mev": a_inv,
            "strange_condensate": ss_mean,
            "strange_condensate_err": ss_err,
            "m_Xt_lattice": m_xt,
            "m_Xt_err_lattice": e_xt,
            "m_Tt_lattice": m_tt,
            "m_Tt_err_lattice": e_tt,
            "delta_m_lattice": dm_lat,
            "delta_m_err_lattice": dm_lat_err,
            "delta_m_mev": dm_mev,
            "delta_m_err_mev": dm_mev_err,
        })

        x_vals.append(ss_mean)
        x_errs.append(ss_err)
        dm_lat_vals.append(dm_lat)
        dm_lat_errs.append(dm_lat_err)
        dm_mev_vals.append(dm_mev)
        dm_mev_errs.append(dm_mev_err)

    df_data = pl.DataFrame(data_records)

    x_arr = np.array(x_vals, dtype=np.float64)
    y_arr = np.array(dm_mev_vals, dtype=np.float64)
    y_err_arr = np.array(dm_mev_errs, dtype=np.float64)
    n_pts = len(x_arr)

    # 3. 拟合四种模型
    # 模型定义
    def f_linear(x, a, b):
        return a * x + b

    def f_power_13(x, a, b):
        return a * (x**(1.0 / 3.0)) + b

    def f_power_23(x, a, b):
        return a * (x**(2.0 / 3.0)) + b

    def f_power_free(x, a, gamma, b):
        return a * (x**gamma) + b

    models = [
        ("Linear Model: A * <ss> + B", f_linear, [1000.0, 0.0], 2, "Linear"),
        ("1/3 Power Model: A * <ss>^(1/3) + B", f_power_13, [500.0, 0.0], 2, "Power 1/3"),
        ("2/3 Power Model: A * <ss>^(2/3) + B", f_power_23, [800.0, 0.0], 2, "Power 2/3"),
        ("Free Power Model: A * <ss>^alpha + B", f_power_free, [800.0, 0.7, 0.0], 3, "Free Power"),
    ]

    fit_results = []
    curve_data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    x_dense = np.linspace(np.min(x_arr) * 0.9, np.max(x_arr) * 1.1, 300)

    for m_name, f_func, p0, k_params, short_name in models:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, pcov = curve_fit(
                    f_func,
                    x_arr,
                    y_arr,
                    p0=p0,
                    sigma=y_err_arr,
                    absolute_sigma=True,
                    maxfev=5000,
                )
            perr = np.sqrt(np.diag(pcov))
            res = y_arr - f_func(x_arr, *popt)
            chi2 = float(np.sum((res / y_err_arr)**2))
            dof = float(n_pts - k_params)
            chi2_dof = chi2 / dof if dof > 0 else np.nan
            aic = chi2 + 2.0 * k_params
            bic = chi2 + k_params * np.log(n_pts)

            ss_tot = float(np.sum((y_arr - np.mean(y_arr))**2))
            ss_res = float(np.sum(res**2))
            r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else np.nan

            p_str = ", ".join([f"{popt[i]:.3e} +/- {perr[i]:.3e}" for i in range(len(popt))])
            fit_results.append({
                "model_name": m_name,
                "short_name": short_name,
                "status": "Success",
                "k_params": k_params,
                "parameters": p_str,
                "param_a": float(popt[0]),
                "param_a_err": float(perr[0]),
                "param_b": float(popt[-1]),
                "param_b_err": float(perr[-1]),
                "param_alpha": float(popt[1]) if k_params == 3 else (1.0 if "Linear" in m_name else (1.0/3.0 if "1/3" in m_name else 2.0/3.0)),
                "param_alpha_err": float(perr[1]) if k_params == 3 else 0.0,
                "chi2": chi2,
                "dof": dof,
                "chi2_dof": chi2_dof,
                "R2": r2,
                "AIC": aic,
                "BIC": bic,
            })
            curve_data[short_name] = (x_dense, f_func(x_dense, *popt))
        except Exception as e:
            fit_results.append({
                "model_name": m_name,
                "short_name": short_name,
                "status": f"Failed ({e})",
                "k_params": k_params,
                "parameters": "None",
                "param_a": np.nan,
                "param_a_err": np.nan,
                "param_b": np.nan,
                "param_b_err": np.nan,
                "param_alpha": np.nan,
                "param_alpha_err": np.nan,
                "chi2": np.nan,
                "dof": np.nan,
                "chi2_dof": np.nan,
                "R2": np.nan,
                "AIC": np.nan,
                "BIC": np.nan,
            })

    df_models = pl.DataFrame(fit_results)

    # 4. 计算贝叶斯对比量
    valid_bics = [r["BIC"] for r in fit_results if np.isfinite(r["BIC"])]
    min_bic = min(valid_bics) if valid_bics else 0.0

    delta_bics = []
    bayes_weights = []
    assessments = []

    for r in fit_results:
        b_val = r["BIC"]
        if np.isfinite(b_val):
            db = b_val - min_bic
            w = np.exp(-0.5 * db)
            if db < 1e-4:
                cat = "Optimal Model (Best Fit)"
            elif db <= 2.0:
                cat = "Substantial Support (Delta_BIC <= 2)"
            elif db <= 6.0:
                cat = "Positive Evidence Against (Delta_BIC in 2-6)"
            else:
                cat = "Strong Evidence Against (Delta_BIC > 6)"
        else:
            db, w, cat = np.nan, 0.0, "Failed"
        delta_bics.append(db)
        bayes_weights.append(w)
        assessments.append(cat)

    sum_w = sum(bayes_weights)
    norm_w = [(w / sum_w if sum_w > 0 else 0.0) for w in bayes_weights]

    df_models = df_models.with_columns([
        pl.Series("delta_BIC", delta_bics),
        pl.Series("bayes_prob_weight", norm_w),
        pl.Series("assessment", assessments),
    ])

    # 5. 绘制模型对比图
    apply_academic_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    # 数据点
    ax.errorbar(
        x_arr,
        y_arr,
        xerr=np.array(x_errs),
        yerr=y_err_arr,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.2,
        capsize=3.5,
        capthick=1.0,
        markersize=6.0,
        label=r"LQCD Data: $\Delta m(X_t - T_t)$",
        zorder=5,
    )

    colors = {
        "Linear": "#1565c0",
        "Power 1/3": "#2e7d32",
        "Power 2/3": "#c62828",
        "Free Power": "#6a1b9a",
    }
    styles = {
        "Linear": "-",
        "Power 1/3": "--",
        "Power 2/3": "-.",
        "Free Power": ":",
    }

    for s_name, (xd, yd) in curve_data.items():
        row = df_models.filter(pl.col("short_name") == s_name)[0]
        c2d = float(row["chi2_dof"][0])
        bic_val = float(row["BIC"][0])
        ax.plot(
            xd,
            yd,
            label=rf"{s_name} ($\chi^2/\mathrm{{dof}}={c2d:.2f}, \mathrm{{BIC}}={bic_val:.1f}$)",
            color=colors.get(s_name, "blue"),
            linestyle=styles.get(s_name, "-"),
            linewidth=1.8,
            zorder=3,
        )

    ax.set_xlabel(r"Strange Quark Chiral Condensate $\langle \bar{s}s \rangle$ (lattice unit)")
    ax.set_ylabel(r"Meson Mass Splitting $\Delta m(X_t - T_t)$ [MeV]")
    ax.set_title(r"$U(1)_A$ Tensor Channel Splitting vs Strange Quark Condensate")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linestyle="--", alpha=0.35, color="#888888")
    ax.legend(loc="upper left", frameon=True, edgecolor="#cccccc")

    plt.tight_layout()
    out_dir_meson = output_dir / "meson_strange_fit"
    out_dir_meson.mkdir(parents=True, exist_ok=True)

    fig_png = out_dir_meson / "meson_strange_condensate_fit.png"
    fig_pdf = out_dir_meson / "meson_strange_condensate_fit.pdf"
    fig.savefig(fig_png, dpi=300)
    fig.savefig(fig_pdf)
    plt.close(fig)

    # 复制到 docs/figures/ 与 output/
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copyfile(fig_png, FIGURES_DIR / "meson_strange_condensate_fit.png")
    shutil.copyfile(fig_pdf, FIGURES_DIR / "meson_strange_condensate_fit.pdf")
    shutil.copyfile(fig_png, PROJECT_ROOT / "output" / "meson_strange_condensate_fit.png")

    df_data.write_csv(out_dir_meson / "meson_strange_condensate_data.csv")
    df_models.write_csv(out_dir_meson / "meson_strange_model_comparison.csv")

    print(f"[OK] 介子与奇异夸克凝聚拟合完成，产物已保存至: {out_dir_meson}")
    return df_data, df_models


# ==============================================================================
# 6. 总控主调度与报告生成
# ==============================================================================

def main() -> None:
    t0 = 157.0
    zm_table = load_zm_table()
    print(f"[INFO] 成功载入 Zm(beta) 查表 (共 {len(zm_table)} 个锚点)")

    datasets = [
        {"name": "p2", "ns": 36, "nt": 12, "is_raw": True},
        {"name": "p3", "ns": 48, "nt": 12, "is_raw": True},
        {"name": "q2", "ns": 48, "nt": 16, "is_raw": False},
        {"name": "q3", "ns": 64, "nt": 16, "is_raw": False},
    ]

    all_lcp_reports: dict[str, dict] = {}

    for ds_info in datasets:
        ds_name = ds_info["name"]
        ns = ds_info["ns"]
        nt = ds_info["nt"]
        print(f"\n" + "=" * 80)
        print(f">>> [LCP PIPELINE] 开始处理数据集: {ds_name.upper()} (Ns={ns}, Nt={nt})...")
        print("=" * 80)

        out_ds_dir = OUTPUT_LCP_DIR / ds_name
        out_ds_dir.mkdir(parents=True, exist_ok=True)

        if ds_info["is_raw"]:
            df_cond, df_susc, jk_dict = process_raw_dataset(ds_name, ns, nt, zm_table)
        else:
            df_cond, df_susc, jk_dict = load_existing_dataset(ds_name, ns, nt)

        # 保存本数据集的独立全量表
        df_cond.write_csv(out_ds_dir / f"condensate_{ds_name}.csv")
        df_susc.write_csv(out_ds_dir / f"susceptibility_{ds_name}.csv")

        # 写入 Jackknife 样本表
        jk_records = []
        for b_val, s_arr in jk_dict.items():
            t_val = float(df_susc.filter(pl.col("beta") == b_val)["temp"][0])
            for idx, val in enumerate(s_arr):
                jk_records.append({
                    "dataset": ds_name,
                    "beta": b_val,
                    "temp": t_val,
                    "jk_index": idx,
                    "chi_renorm": val,
                })
        if jk_records:
            pl.DataFrame(jk_records).write_csv(out_ds_dir / f"susceptibility_jk_samples_{ds_name}.csv")

        # 1. 绘制手征凝聚图
        cond_plot_path = plot_chiral_condensate(df_cond, ds_name, out_ds_dir)
        print(f"  [PLOT] 手征凝聚图已生成: {cond_plot_path.name}")

        # 2. 执行 Padé 阶数拟合与贝叶斯分析
        temps_arr = df_susc["temp"].to_numpy().astype(float)
        chi_arr = df_susc["chi_renorm"].to_numpy().astype(float)
        err_arr = df_susc["chi_renorm_err"].to_numpy().astype(float)

        df_bayes = perform_pade_model_comparison(temps_arr, chi_arr, err_arr, t0=t0)
        df_bayes.write_csv(out_ds_dir / f"pade_model_comparison_{ds_name}.csv")
        print(f"  [BAYES] Padé 贝叶斯模型选择与信息损失表已生成: pade_model_comparison_{ds_name}.csv")

        # 3. 绘制手征磁化率 + Padé 拟合图
        susc_plot_path = plot_chiral_susceptibility_with_pade(df_susc, df_bayes, ds_name, out_ds_dir, t0=t0)
        print(f"  [PLOT] 手征磁化率 Padé 拟合图已生成: {susc_plot_path.name}")

        all_lcp_reports[ds_name] = {
            "df_cond": df_cond,
            "df_susc": df_susc,
            "df_bayes": df_bayes,
            "out_dir": out_ds_dir,
        }

    # 执行 Requirement 2: 介子质量差与 Strange 凝聚拟合
    df_meson_data, df_meson_models = perform_meson_strange_condensate_analysis(OUTPUT_LCP_DIR)

    print("\n" + "=" * 80)
    print(">>> [ALL TASKS COMPLETED] 所有 LCP (p2, p3, q2, q3) 与介子手征拟合任务执行完毕！")
    print("=" * 80)


if __name__ == "__main__":
    main()
