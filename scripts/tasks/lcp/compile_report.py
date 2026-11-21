#!/usr/bin/env python3
"""
scripts/tasks/lcp/compile_report.py
--------------------------------------------------------------------------------
常物理线 (LCP) LaTeX 报告编译脚本:
1. 读取各系综拟合统计量与测量数据，生成规范的数据表格 (Table 1 - 4)
2. 保持用户在 docs/report.tex 中编写的说明与排版
3. 调用 pdflatex 编译生成 docs/report.pdf
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DOCS_DIR = PROJECT_ROOT / "docs"
BASE_DIR = PROJECT_ROOT / "22pade"


def format_sci(val: float, err: float | None = None, digits: int = 2) -> str:
    """标准 LaTeX 科学计数法格式化"""
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


def build_report_tex() -> str:
    """构建保持纯图表与用户分析文字的完整 report.tex"""
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

\title{Chiral Condensate and Susceptibility along the LCP}
\author{Junxiong Nie}
\date{November 10, 2026}

\begin{document}
\maketitle

\[
\frac{\langle\bar{\psi}\psi\rangle_{\rm sub}}{Z_m}
= \langle\bar{\psi}\psi\rangle_l - \frac{m_l + m_{\rm res}}{m_s + m_{\rm res}}\,\langle\bar{\psi}\psi\rangle_s
\]

\[
a^2\chi_{\rm disc} = V\left[\left\langle\left(\frac{1}{V}\mathrm{Tr}\,D_l^{-1}\right)^2\right\rangle - \left\langle\frac{1}{V}\mathrm{Tr}\,D_l^{-1}\right\rangle^2\right],
\qquad
\chi_{\rm disc}(T) = (N_t T)^2\, a^2\chi_{\rm disc}
\]

\[
\ln\!\left[\frac{a^{-1}(\beta)}{\mathrm{MeV}}\right] = -2.1966\,\beta^2 + 20.7717\,\beta - 40.6167
\]

Shaded bands indicate Jackknife envelope uncertainty, and solid lines denote mean parameter fits.
\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02.pdf}
    \caption{Renormalized chiral condensate $\langle\bar{\psi}\psi\rangle_{\rm sub}/Z_m$ in lattice units (left) and renormalized chiral susceptibility $\chi_{\mathrm{disc}}/Z_m^2\,[\mathrm{GeV}^2]$ (right) comparison across all four LCP ensembles with Padé configuration: $p_2\,[2/2]$, $p_3\,[0/2]$, $q_2\,[2/2]$, and $q_3\,[0/2]$. }
    \label{fig:combined_p3_02_q3_02}
\end{figure}

The value of condensation in high temperature may be negative , light condensate decrease much faster than the strange condensate. 

\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/chiral_condensate_bare_light_strange.pdf}
    \caption{Bare light quark condensate $\langle\bar{\psi}\psi\rangle_l$ (left) and bare strange quark condensate $\langle\bar{\psi}\psi\rangle_s$ (right) in lattice units across all four LCP ensembles ($p_2, p_3, q_2, q_3$).}
    \label{fig:bare_light_strange}
\end{figure}

The height of peak, value of susceptibility, varies with Volume size slowly. The number of data points in q3 ensemble maybe  not so much but it's clearly larger than other ensembles.


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
    \caption{Multi-order Padé approximation comparison per ensemble with shaded Jackknife confidence bands. }
    \label{fig:multimodal_all}
\end{figure}
Spurious orders suffering from real-axis poles or zero degrees of freedom have been excluded.

I picked the order with reduced chi2 close to 1, but the number of parameter, k, is also important.

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

    # Table 1: 核心系综对应图 1 设定
    core_cases = [
        ("p2", "36^3\\times 12", "2/2"),
        ("p3", "48^3\\times 12", "0/2"),
        ("q2", "48^3\\times 16", "2/2"),
        ("q3", "64^3\\times 16", "0/2"),
    ]
    for ds_name, dims, ord_val in core_cases:
        ord_tag = ord_val.replace("/", "")
        sum_df = pl.read_csv(BASE_DIR / ds_name / f"pade{ord_tag}_fit_summary_{ds_name}.csv")
        sm = sum_df.to_dicts()[0]
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
    \begin{tabular}{lccccc}
        \toprule
        \textbf{Ensemble} & \textbf{Order} & $\boldsymbol{N_{\mathrm{data}}}$ & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & $\boldsymbol{T_{\mathrm{peak}}}$ [MeV] \\
        \midrule
""")

    # Table 2: 多阶 Padé 模型比较
    selected_fits = [
        ("p2", ["1/2", "2/2", "3/2"]),
        ("p3", ["0/2", "1/2", "2/2"]),
        ("q2", ["1/2", "2/2", "3/2", "1/3"]),
        ("q3", ["0/2", "1/2"]),
    ]
    n_data_map = {"p2": 10, "p3": 7, "q2": 8, "q3": 5}
    k_params_map = {"0/2": 3, "1/2": 4, "2/2": 5, "3/2": 6, "1/3": 5}

    for ds_name, ords in selected_fits:
        n_pts = n_data_map[ds_name]
        for ord_val in ords:
            ord_tag = ord_val.replace("/", "")
            sum_df = pl.read_csv(BASE_DIR / ds_name / f"pade{ord_tag}_fit_summary_{ds_name}.csv")
            sm = sum_df.to_dicts()[0]
            k_params = k_params_map[ord_val]
            chi2_dof = sm["chi2_dof_mean"]
            tpk = sm["t_peak_mean"]
            ens_tex = f"${ds_name[0]}_{{{ds_name[1:]}}}$"
            tex.append(f"        {ens_tex} & [{ord_val}] & {n_pts} & {k_params} & ${chi2_dof:.3f}$ & ${tpk:.2f}$ \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{\textbf{FLAG-Style Error Budget Table}: Systematic uncertainty decomposition for the transition temperature $T_{pc}$ and peak susceptibility $\chi_{\mathrm{peak}}/Z_m^2$ anchored at the $64^3 \times 16$ ($q_3$) baseline.}
    \label{tab:flag_error_budget}
    \small
    \begin{tabular}{lcccc}
        \toprule
        \textbf{Uncertainty Source} & \multicolumn{2}{c}{$\boldsymbol{T_{pc}}$} & \multicolumn{2}{c}{$\boldsymbol{\chi_{\mathrm{peak}}/Z_m^2}$} \\
        \cmidrule(lr){2-3} \cmidrule(lr){4-5}
        & \textbf{Absolute [MeV]} & \textbf{Relative [\%]} & \textbf{Absolute [$\text{GeV}^2$]} & \textbf{Relative [\%]} \\
        \midrule
        \textbf{Central Reference Value ($64^3\times 16$)} & \textbf{155.12} & --- & \textbf{1.038} & --- \\
        \midrule
        Discretization / $a^2\to 0$ ($\Delta_{\mathrm{disc}}$) & $\pm 1.72$ & 1.11\% & $\pm 0.109$ & 10.47\% \\
        Finite Volume / $L\to\infty$ ($\Delta_{\mathrm{FV}}$) & $\pm 1.16$ & 0.75\% & $\pm 0.097$ & 9.35\% \\
        \midrule
        Combined Systematic ($\sigma_{\mathrm{sys}}$) & $\pm 2.07$ & 1.34\% & $\pm 0.146$ & 14.04\% \\
        \midrule
        \multicolumn{5}{l}{\textbf{Standard FLAG Result Presentation:}} \\
        \multicolumn{5}{c}{$T_{pc} = \mathbf{155.12 \pm 1.72_{\mathrm{disc}} \pm 1.16_{\mathrm{FV}}\text{ MeV}} = \mathbf{155.12 \pm 2.07_{\mathrm{sys}}\text{ MeV}} = \mathbf{155.12(207)\text{ MeV}}$} \\
        \multicolumn{5}{c}{$\chi_{\mathrm{peak}}/Z_m^2 = \mathbf{1.038 \pm 0.109_{\mathrm{disc}} \pm 0.097_{\mathrm{FV}}\text{ GeV}^2} = \mathbf{1.038 \pm 0.146_{\mathrm{sys}}\text{ GeV}^2} = \mathbf{1.038(146)\text{ GeV}^2}$} \\
        \bottomrule
    \end{tabular}
\end{table}

\clearpage

\begin{table}[htbp]
    \centering
    \caption{Lattice measurements for $N_t=12$ ensembles ($p_2$ and $p_3$): Temperature $T$, inverse lattice spacing $a^{-1}$, renormalized chiral condensate $\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m$ in lattice units, and renormalized susceptibility $\chi_{\mathrm{disc}} / Z_m^2$ [$\text{GeV}^2$].}
    \label{tab:raw_p2_p3}
    \small
    \begin{tabular}{lccccc}
        \toprule
        \textbf{Ensemble} & $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{a^{-1}}$ [GeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m}$ & $\boldsymbol{\chi_{\mathrm{disc}} / Z_m^2}$ [$\text{GeV}^2$] \\
        \midrule
""")

    for ds_name, dim_str in [("p2", "36^3\\times 12"), ("p3", "48^3\\times 12")]:
        df_c = pl.read_csv(BASE_DIR / ds_name / f"condensate_{ds_name}.csv").sort("temp")
        df_s = pl.read_csv(BASE_DIR / ds_name / f"susceptibility_{ds_name}.csv").sort("temp")
        c_dict = {float(r["beta"]): r for r in df_c.to_dicts()}

        for r in df_s.to_dicts():
            b = float(r["beta"])
            t_val = float(r["temp"])
            chi_val = float(r["chi_renorm"])
            chi_err = float(r["chi_renorm_err"])
            c_row = c_dict.get(b, {})
            a_inv = float(c_row.get("a_inv", 12 * t_val)) / 1000.0
            pbp_lat = float(c_row.get("pbp_renorm", 0.0))
            pbp_err_lat = float(c_row.get("pbp_renorm_err", 0.0))
            pbp_str = format_sci(pbp_lat, pbp_err_lat, digits=2)
            chi_str = format_sci(chi_val, chi_err, digits=3)
            tex.append(f"        ${ds_name}$ (${dim_str}$) & {b:.3f} & {t_val:.1f} & {a_inv:.3f} & {pbp_str} & {chi_str} \\\\")

        if ds_name == "p2":
            tex.append(r"""        \midrule""")

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

    for ds_name, dim_str in [("q2", "48^3\\times 16"), ("q3", "64^3\\times 16")]:
        df_c = pl.read_csv(BASE_DIR / ds_name / f"condensate_{ds_name}.csv").sort("temp")
        df_s = pl.read_csv(BASE_DIR / ds_name / f"susceptibility_{ds_name}.csv").sort("temp")
        c_dict = {float(r["beta"]): r for r in df_c.to_dicts()}

        for r in df_s.to_dicts():
            b = float(r["beta"])
            t_val = float(r["temp"])
            chi_val = float(r["chi_renorm"])
            chi_err = float(r["chi_renorm_err"])
            c_row = c_dict.get(b, {})
            a_inv = float(c_row.get("a_inv", 16 * t_val)) / 1000.0
            pbp_lat = float(c_row.get("pbp_renorm", 0.0))
            pbp_err_lat = float(c_row.get("pbp_renorm_err", 0.0))
            pbp_str = format_sci(pbp_lat, pbp_err_lat, digits=2)
            chi_str = format_sci(chi_val, chi_err, digits=3)
            tex.append(f"        ${ds_name}$ (${dim_str}$) & {b:.3f} & {t_val:.1f} & {a_inv:.3f} & {pbp_str} & {chi_str} \\\\")

        if ds_name == "q2":
            tex.append(r"""        \midrule""")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\clearpage

\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.88\textwidth]{figures/chiral_condensate_physical_all.pdf}
    \caption{Renormalized chiral condensate $\langle\bar{\psi}\psi\rangle_{\mathrm{sub}}/Z_m$ in physical units [$\text{GeV}^3$] across all four LCP ensembles ($p_2, p_3, q_2, q_3$).}
    \label{fig:condensate_physical_all}
\end{figure}

\end{document}
""")
    return "\n".join(tex)


def compile_pdf() -> Path:
    """写入 report.tex 并调用 pdflatex 编译"""
    tex_content = build_report_tex()
    tex_path = DOCS_DIR / "report.tex"
    tex_path.write_text(tex_content, encoding="utf-8")
    print(f"[INFO] 成功写入 {tex_path}")

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

    # 再次编译以解析交叉引用
    subprocess.run(cmd, cwd=DOCS_DIR, capture_output=True, text=True)

    pdf_path = DOCS_DIR / "report.pdf"
    if not pdf_path.exists():
        print("[ERROR] 未找到生成的 report.pdf")
        sys.exit(1)

    print(f"[SUCCESS] PDF 编译完成: {pdf_path} (大小: {pdf_path.stat().st_size / 1024:.1f} KB)")
    return pdf_path


def main() -> None:
    print("==========================================================================")
    print("[LCP-REPORT] 启动报告表格更新与 PDF 编译...")
    print("==========================================================================")
    compile_pdf()
    print("\n[DONE] 报告编译完成！")


if __name__ == "__main__":
    main()
