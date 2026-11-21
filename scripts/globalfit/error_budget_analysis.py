#!/usr/bin/env python3
"""
scripts/globalfit/error_budget_analysis.py
--------------------------------------------------------------------------------
FLAG 规范手征相变温度与峰值磁化率误差分解表 (FLAG-Style Error Budget Table)
--------------------------------------------------------------------------------
【设计与标准】:
1. 参考基准 (Reference Baseline Anchor):
   严格遵循 FLAG (Flavour Lattice Averaging Group) 标准, 以现有最大体积与最精细格距的
   64^3 x 16 (q3) 系综作为物理参考基准。

2. 关于样本数对齐的物理说明:
   不同系综 (如 p2, p3, q2, q3) 来自相互独立的马尔可夫链蒙特卡洛轨迹, 样本索引之间
   没有点对点的物理配对关联。因此, 每个系综完整保留其自身真实的全部构型样本数
   (p2: 2003, p3: 1923, q2: 2000, q3: 1171), 保证统计量的无损性与数学严谨性。

3. 不确定度来源分解 (Sources of Uncertainty):
   - 统计误差 (Statistical Error, sigma_stat):
     来自基准系综 q3 自身全部 1171 个 Jackknife 样本计算的标准误:
     sigma_stat = sqrt(N-1) * std(X_k)
   - 有限格距/离散化效应 (Discretization / a -> 0, Delta_disc):
     在固定大空间纵横比 L/Nt = 4 下, 比较粗格点 p3 (Nt=12) 与基准细格点 q3 (Nt=16)。
     取样本分布中偏离基准值的最大绝对差值作为保守系统误差上界。
   - 有限体积效应 (Finite Volume / L -> infty, Delta_FV):
     在固定细格距 Nt=16 下, 比较较小体积 q2 (L=48) 与基准大体积 q3 (L=64)。
     取样本分布中偏离基准值的最大绝对差值作为保守系统误差上界。
   - 合成系统误差 (Total Systematic, sigma_sys):
     sigma_sys = sqrt(Delta_disc^2 + Delta_FV^2)
   - 全总误差 (Total Error, sigma_tot):
     sigma_tot = sqrt(sigma_stat^2 + sigma_sys^2)

4. 表格输出:
   严格对齐 FLAG 报告的 Error Budget 标准排版 (绝对误差 + 相对百分比 [%] + FLAG 紧凑括号表述法)。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
from typing import Dict

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = PROJECT_ROOT / "docs"

ENSEMBLES_CONFIG = {
    "p2": {
        "name": "$p_2$",
        "dim": "$36^3 \\times 12$",
        "L": 36,
        "Nt": 12,
        "aspect_ratio": 3.0,
        "order": "[2/2]",
        "summary_csv": "22pade/p2/pade22_fit_summary_p2.csv",
        "samples_csv": "22pade/p2/pade22_jk_samples_p2.csv",
    },
    "p3": {
        "name": "$p_3$",
        "dim": "$48^3 \\times 12$",
        "L": 48,
        "Nt": 12,
        "aspect_ratio": 4.0,
        "order": "[0/2]",
        "summary_csv": "22pade/p3/pade02_fit_summary_p3.csv",
        "samples_csv": "22pade/p3/pade02_jk_samples_p3.csv",
    },
    "q2": {
        "name": "$q_2$",
        "dim": "$48^3 \\times 16$",
        "L": 48,
        "Nt": 16,
        "aspect_ratio": 3.0,
        "order": "[2/2]",
        "summary_csv": "22pade/q2/pade22_fit_summary_q2.csv",
        "samples_csv": "22pade/q2/pade22_jk_samples_q2.csv",
    },
    "q3": {
        "name": "$q_3$",
        "dim": "$64^3 \\times 16$",
        "L": 64,
        "Nt": 16,
        "aspect_ratio": 4.0,
        "order": "[0/2]",
        "summary_csv": "22pade/q3/pade02_fit_summary_q3.csv",
        "samples_csv": "22pade/q3/pade02_jk_samples_q3.csv",
    },
}


def load_full_ensemble_data() -> Dict[str, dict]:
    """读取各系综最优拟合的均值、误差与全部真实有效 Jackknife 样本 (完整保留构型数)"""
    data = {}
    for key, cfg in ENSEMBLES_CONFIG.items():
        sum_p = PROJECT_ROOT / cfg["summary_csv"]
        sam_p = PROJECT_ROOT / cfg["samples_csv"]

        if not sum_p.exists() or not sam_p.exists():
            raise FileNotFoundError(f"找不到系综 {key} 的拟合产物: {sum_p} 或 {sam_p}")

        df_sum = pl.read_csv(sum_p)
        df_sam = pl.read_csv(sam_p)

        valid_mask = (df_sam["ok"] == True) & df_sam["t_peak"].is_not_nan() & df_sam["chi_peak"].is_not_nan()
        df_valid = df_sam.filter(valid_mask)

        t_samples = df_valid["t_peak"].to_numpy().astype(np.float64)
        chi_samples = df_valid["chi_peak"].to_numpy().astype(np.float64)

        t_mean = float(df_sum["t_peak_mean"][0])
        t_err = float(df_sum["t_peak_err"][0])
        chi_mean = float(df_sum["chi_peak_mean"][0])
        chi_err = float(df_sum["chi_peak_err"][0])
        chi2_dof = float(df_sum["chi2_dof_mean"][0]) if "chi2_dof_mean" in df_sum.columns else float("nan")

        data[key] = {
            **cfg,
            "t_mean": t_mean,
            "t_err": t_err,
            "chi_mean": chi_mean,
            "chi_err": chi_err,
            "chi2_dof": chi2_dof,
            "t_samples": t_samples,
            "chi_samples": chi_samples,
            "n_samples": len(t_samples),
        }
    return data


def compute_flag_error_budget(data: Dict[str, dict]) -> dict:
    """严格以 64x16 (q3) 为参考基准，计算 FLAG 规范误差分解"""
    q3 = data["q3"]
    q2 = data["q2"]
    p3 = data["p3"]
    p2 = data["p2"]

    # 1. 基准值与统计误差 (来自 q3)
    T_base = q3["t_mean"]
    T_stat = q3["t_err"]
    chi_base = q3["chi_mean"]
    chi_stat = q3["chi_err"]

    # 2. 有限格距效应 (Discretization / a^2 -> 0):
    # 主考察项: 固定 L/Nt=4 下, p3 (Nt=12) 与基准 q3 (Nt=16) 对比
    # 取样本相对基准中心值的最大绝对差值
    T_disc_max = float(np.max(np.abs(p3["t_samples"] - T_base)))
    chi_disc_max = float(np.max(np.abs(p3["chi_samples"] - chi_base)))
    T_disc_mean = abs(p3["t_mean"] - T_base)
    chi_disc_mean = abs(p3["chi_mean"] - chi_base)

    # 3. 有限体积效应 (Finite Volume / L -> infty):
    # 主考察项: 固定 Nt=16 下, 较小体积 q2 (L=48) 与基准大体积 q3 (L=64) 对比
    T_vol_max = float(np.max(np.abs(q2["t_samples"] - T_base)))
    chi_vol_max = float(np.max(np.abs(q2["chi_samples"] - chi_base)))
    T_vol_mean = abs(q2["t_mean"] - T_base)
    chi_vol_mean = abs(q2["chi_mean"] - chi_base)

    # 4. 合成系统误差与总误差
    T_sys = float(np.sqrt(T_disc_max**2 + T_vol_max**2))
    chi_sys = float(np.sqrt(chi_disc_max**2 + chi_vol_max**2))

    T_tot = float(np.sqrt(T_stat**2 + T_sys**2))
    chi_tot = float(np.sqrt(chi_stat**2 + chi_sys**2))

    # 计算相对百分比不确定度 [%]
    pct = lambda err, base: (err / base) * 100.0

    return {
        "T": {
            "base": T_base,
            "stat": T_stat,
            "stat_pct": pct(T_stat, T_base),
            "disc": T_disc_max,
            "disc_mean": T_disc_mean,
            "disc_pct": pct(T_disc_max, T_base),
            "vol": T_vol_max,
            "vol_mean": T_vol_mean,
            "vol_pct": pct(T_vol_max, T_base),
            "sys": T_sys,
            "sys_pct": pct(T_sys, T_base),
            "tot": T_tot,
            "tot_pct": pct(T_tot, T_base),
        },
        "chi": {
            "base": chi_base,
            "stat": chi_stat,
            "stat_pct": pct(chi_stat, chi_base),
            "disc": chi_disc_max,
            "disc_mean": chi_disc_mean,
            "disc_pct": pct(chi_disc_max, chi_base),
            "vol": chi_vol_max,
            "vol_mean": chi_vol_mean,
            "vol_pct": pct(chi_vol_max, chi_base),
            "sys": chi_sys,
            "sys_pct": pct(chi_sys, chi_base),
            "tot": chi_tot,
            "tot_pct": pct(chi_tot, chi_base),
        },
    }


def export_flag_latex_tables(data: Dict[str, dict], budget: dict, out_tex_path: Path) -> None:
    """输出标准 FLAG 风格的 LaTeX 表格代码"""
    t = budget["T"]
    c = budget["chi"]

    tex_code = []
    tex_code.append("% ==============================================================================")
    tex_code.append("% Auto-generated FLAG-Style Error Budget LaTeX Tables")
    tex_code.append("% Generated by scripts/globalfit/error_budget_analysis.py")
    tex_code.append("% ==============================================================================")
    tex_code.append("")

    # 表 1: LCP 四组核心系综输入与拟合汇总
    tex_code.append(r"\begin{table}[htbp]")
    tex_code.append(r"    \centering")
    tex_code.append(r"    \caption{Summary of lattice geometries, primary Padé fit orders, and full Monte Carlo Jackknife statistics across the four core LCP ensembles.}")
    tex_code.append(r"    \label{tab:flag_ensemble_inputs}")
    tex_code.append(r"    \small")
    tex_code.append(r"    \begin{tabular}{lcccccc}")
    tex_code.append(r"        \toprule")
    tex_code.append(r"        \textbf{Ensemble} & \textbf{Lattice Dimensions} & \textbf{Padé Order} & $\boldsymbol{T_{pc}}$ [MeV] & $\boldsymbol{\chi_{\mathrm{peak}}/Z_m^2}$ [$\text{GeV}^2$] & $\boldsymbol{\chi^2/\mathrm{dof}}$ & \textbf{Jackknife Samples} \\")
    tex_code.append(r"        \midrule")

    for k in ["p2", "p3", "q2", "q3"]:
        d = data[k]
        chi2_str = f"{d['chi2_dof']:.3f}" if np.isfinite(d['chi2_dof']) else "-"
        note = " (Reference Anchor)" if k == "q3" else ""
        tex_code.append(
            f"        {d['name']}{note} & {d['dim']} & {d['order']} & "
            f"${d['t_mean']:.2f} \\pm {d['t_err']:.2f}$ & "
            f"${d['chi_mean']:.3f} \\pm {d['chi_err']:.3f}$ & "
            f"{chi2_str} & {d['n_samples']} \\\\"
        )
    tex_code.append(r"        \bottomrule")
    tex_code.append(r"    \end{tabular}")
    tex_code.append(r"\end{table}")
    tex_code.append("")

    # 表 2: FLAG 标准 Error Budget 表格 (包含绝对误差与相对百分比不确定度，已去除统计误差)
    tex_code.append(r"\begin{table}[htbp]")
    tex_code.append(r"    \centering")
    tex_code.append(r"    \caption{\textbf{FLAG-Style Error Budget Table}: Systematic uncertainty decomposition for the transition temperature $T_{pc}$ and peak susceptibility $\chi_{\mathrm{peak}}/Z_m^2$ anchored at the $64^3 \times 16$ ($q_3$) physical baseline.}")
    tex_code.append(r"    \label{tab:flag_error_budget}")
    tex_code.append(r"    \small")
    tex_code.append(r"    \begin{tabular}{lcccc}")
    tex_code.append(r"        \toprule")
    tex_code.append(r"        \textbf{Uncertainty Source} & \multicolumn{2}{c}{$\boldsymbol{T_{pc}}$} & \multicolumn{2}{c}{$\boldsymbol{\chi_{\mathrm{peak}}/Z_m^2}$} \\")
    tex_code.append(r"        \cmidrule(lr){2-3} \cmidrule(lr){4-5}")
    tex_code.append(r"        & \textbf{Absolute [MeV]} & \textbf{Relative [\%]} & \textbf{Absolute [$\text{GeV}^2$]} & \textbf{Relative [\%]} \\")
    tex_code.append(r"        \midrule")
    tex_code.append(f"        \\textbf{{Central Reference Value ($64^3\\times 16$)}} & \\textbf{{{t['base']:.2f}}} & --- & \\textbf{{{c['base']:.3f}}} & --- \\\\")
    tex_code.append(r"        \midrule")
    tex_code.append(f"        Discretization / $a^2\\to 0$ ($\\Delta_{{\\mathrm{{disc}}}}$) & $\\pm {t['disc']:.2f}$ & {t['disc_pct']:.2f}\\% & $\\pm {c['disc']:.3f}$ & {c['disc_pct']:.2f}\\% \\\\")
    tex_code.append(f"        Finite Volume / $L\\to\\infty$ ($\\Delta_{{\\mathrm{{FV}}}}$) & $\\pm {t['vol']:.2f}$ & {t['vol_pct']:.2f}\\% & $\\pm {c['vol']:.3f}$ & {c['vol_pct']:.2f}\\% \\\\")
    tex_code.append(r"        \midrule")
    tex_code.append(f"        Combined Systematic ($\\sigma_{{\\mathrm{{sys}}}}$) & $\\pm {t['sys']:.2f}$ & {t['sys_pct']:.2f}\\% & $\\pm {c['sys']:.3f}$ & {c['sys_pct']:.2f}\\% \\\\")
    tex_code.append(r"        \midrule")
    tex_code.append(r"        \multicolumn{5}{l}{\textbf{Standard FLAG Result Presentation:}} \\")
    tex_code.append(
        r"        \multicolumn{5}{c}{"
        f"$T_{{pc}} = \\mathbf{{{t['base']:.2f} \\pm {t['disc']:.2f}_{{\\mathrm{{disc}}}} \\pm {t['vol']:.2f}_{{\\mathrm{{FV}}}}\\text{{ MeV}}}} "
        f"= \\mathbf{{{t['base']:.2f} \\pm {t['sys']:.2f}_{{\\mathrm{{sys}}}}}} "
        f"= \\mathbf{{{t['base']:.2f}({int(round(t['sys']*100)):02d})\\text{{ MeV}}}}$"
        r"} \\")
    tex_code.append(
        r"        \multicolumn{5}{c}{"
        f"$\\chi_{{\\mathrm{{peak}}}}/Z_m^2 = \\mathbf{{{c['base']:.3f} \\pm {c['disc']:.3f}_{{\\mathrm{{disc}}}} \\pm {c['vol']:.3f}_{{\\mathrm{{FV}}}}\\text{{ GeV}}^2}} "
        f"= \\mathbf{{{c['base']:.3f} \\pm {c['sys']:.3f}_{{\\mathrm{{sys}}}}}} "
        f"= \\mathbf{{{c['base']:.3f}({int(round(c['sys']*1000)):02d})\\text{{ GeV}}^2}}$"
        r"} \\")
    tex_code.append(r"        \bottomrule")
    tex_code.append(r"    \end{tabular}")
    tex_code.append(r"\end{table}")

    out_tex_path.write_text("\n".join(tex_code), encoding="utf-8")
    print(f"[OK] 成功导出 FLAG 格式 LaTeX 表格: {out_tex_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="FLAG 规范手征相变 Tc 与 chi_peak 误差分解分析")
    parser.add_argument("--compile", action="store_true", default=True, help="是否自动调用 pdflatex 编译 PDF")
    args = parser.parse_args()

    print("==========================================================================")
    print(">>> [FLAG ERROR BUDGET] 启动以 64x16 为基准的 FLAG 标准误差分解分析")
    print("==========================================================================")

    data = load_full_ensemble_data()
    print("\n[INFO] 各系综真实有效样本读取完毕 (无损保留完整构型规模):")
    for k, d in data.items():
        print(f"  - {d['name']} ({d['dim']}): Tc = {d['t_mean']:.2f} +/- {d['t_err']:.2f} MeV, "
              f"chi_peak = {d['chi_mean']:.3f} +/- {d['chi_err']:.3f}, chi2/dof = {d['chi2_dof']:.3f}, "
              f"有效样本数 = {d['n_samples']}")

    budget = compute_flag_error_budget(data)
    t = budget["T"]
    c = budget["chi"]

    print("\n==========================================================================")
    print(">>> [FLAG ERROR BUDGET SUMMARY]")
    print(f"  * T_pc:     {t['base']:.2f} +/- {t['stat']:.2f}(stat, {t['stat_pct']:.2f}%) "
          f"+/- {t['disc']:.2f}(disc, {t['disc_pct']:.2f}%) "
          f"+/- {t['vol']:.2f}(FV, {t['vol_pct']:.2f}%) "
          f"-> 合成sys: +/- {t['sys']:.2f} ({t['sys_pct']:.2f}%), 总tot: +/- {t['tot']:.2f} ({t['tot_pct']:.2f}%)")
    print(f"  * chi_peak: {c['base']:.3f} +/- {c['stat']:.3f}(stat, {c['stat_pct']:.2f}%) "
          f"+/- {c['disc']:.3f}(disc, {c['disc_pct']:.2f}%) "
          f"+/- {c['vol']:.3f}(FV, {c['vol_pct']:.2f}%) "
          f"-> 合成sys: +/- {c['sys']:.3f} ({c['sys_pct']:.2f}%), 总tot: +/- {c['tot']:.3f} ({c['tot_pct']:.2f}%)")
    print("==========================================================================")

    out_tables_tex = DOCS_DIR / "error_budget_tables.tex"
    export_flag_latex_tables(data, budget, out_tables_tex)

    if args.compile:
        main_tex = DOCS_DIR / "error_budget.tex"
        if main_tex.exists():
            print(f"\n[INFO] 正在调用 pdflatex 编译 {main_tex} ...")
            cmd = ["pdflatex", "-interaction=nonstopmode", main_tex.name]
            res = subprocess.run(cmd, cwd=str(DOCS_DIR), capture_output=True, text=True, errors="replace")
            if res.returncode == 0:
                print(f"[SUCCESS] 成功生成 PDF: {DOCS_DIR / 'error_budget.pdf'}")
            else:
                print(f"[WARN] 首次编译退出码 {res.returncode}，再次运行以稳定交叉引用...")
                subprocess.run(cmd, cwd=str(DOCS_DIR), capture_output=True, text=True, errors="replace")
                if (DOCS_DIR / "error_budget.pdf").exists():
                    print(f"[SUCCESS] 最终生成 PDF: {DOCS_DIR / 'error_budget.pdf'}")


if __name__ == "__main__":
    main()
