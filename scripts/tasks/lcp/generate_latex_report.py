#!/usr/bin/env python3
"""
scripts/tasks/lcp/generate_latex_report.py
--------------------------------------------------------------------------------
自动生成 docs/report.tex 并调用 pdflatex 编译生成 docs/report.pdf。
严格遵循用户要求：“暂时只放图和数据，什么都不加”。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import polars as pl
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DOCS_DIR = PROJECT_ROOT / "docs"
BASE_DIR = PROJECT_ROOT / "22pade"


def format_num(val, fmt="{:.3f}", na_str="--"):
    if val is None or not np.isfinite(val):
        return na_str
    return fmt.format(val)


def build_latex_document() -> str:
    # 1. 读取多模型比较表格数据
    p2_cmp = pl.read_csv(PROJECT_ROOT / "output/LCP/p2/pade_model_comparison_p2.csv")
    p3_cmp = pl.read_csv(PROJECT_ROOT / "output/LCP/p3/pade_model_comparison_p3.csv")
    q2_cmp = pl.read_csv(PROJECT_ROOT / "output/LCP/q2/pade_model_comparison_q2.csv")
    q3_cmp = pl.read_csv(PROJECT_ROOT / "output/LCP/q3/pade_model_comparison_q3.csv")

    # 2. 读取各系综拟合汇总
    p2_sum = pl.read_csv(BASE_DIR / "p2/pade22_fit_summary_p2.csv").to_dicts()[0]
    p3_sum_22 = pl.read_csv(BASE_DIR / "p3/pade22_fit_summary_p3.csv").to_dicts()[0]
    p3_sum_02 = pl.read_csv(BASE_DIR / "p3/pade02_fit_summary_p3.csv").to_dicts()[0]
    q2_sum = pl.read_csv(BASE_DIR / "q2/pade22_fit_summary_q2.csv").to_dicts()[0]
    q3_sum_22 = pl.read_csv(BASE_DIR / "q3/pade22_fit_summary_q3.csv").to_dicts()[0]
    q3_sum_02 = pl.read_csv(BASE_DIR / "q3/pade02_fit_summary_q3.csv").to_dicts()[0]
    q3_sum_12 = pl.read_csv(BASE_DIR / "q3/pade12_fit_summary_q3.csv").to_dicts()[0]

    # 3. 读取介子质量差拟合数据
    meson_cmp = pl.read_csv(BASE_DIR / "meson_strange_fit/meson_strange_model_comparison.csv")
    meson_data = pl.read_csv(BASE_DIR / "meson_strange_fit/meson_strange_condensate_data.csv")

    # 4. 原始数据点
    p2_susc = pl.read_csv(BASE_DIR / "p2/susceptibility_p2.csv").sort("temp")
    p2_cond = pl.read_csv(BASE_DIR / "p2/condensate_p2.csv").sort("temp")
    p3_susc = pl.read_csv(BASE_DIR / "p3/susceptibility_p3.csv").sort("temp")
    p3_cond = pl.read_csv(BASE_DIR / "p3/condensate_p3.csv").sort("temp")
    q2_susc = pl.read_csv(BASE_DIR / "q2/susceptibility_q2.csv").sort("temp")
    q2_cond = pl.read_csv(BASE_DIR / "q2/condensate_q2.csv").sort("temp")
    q3_susc = pl.read_csv(BASE_DIR / "q3/susceptibility_q3.csv").sort("temp")
    q3_cond = pl.read_csv(BASE_DIR / "q3/condensate_q3.csv").sort("temp")

    tex = []
    tex.append(r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=1.8cm]{geometry}
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

\title{\textbf{Lattice QCD Chiral Phase Transition and Meson Analysis}\\[0.3em]
\large Figures and Numerical Data Report}
\author{Antigravity AI Pair Programming Session}
\date{\today}

\begin{document}
\maketitle

\section{Combined Ensemble Comparisons}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/chiral_lcp_comparison_p2_22_p3_02_q2_22_q3_02.pdf}
    \caption{Chiral condensate (left) and renormalized chiral susceptibility (right) comparison across all four LCP ensembles using the Padé configuration: $p_2\,[2/2]$, $p_3\,[0/2]$, $q_2\,[2/2]$, and $q_3\,[0/2]$. Shaded bands represent Jackknife envelope fit ranges, and solid curves denote mean parameter fits.}
    \label{fig:combined_p3_02_q3_02}
\end{figure}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/chiral_lcp_comparison_p2_22_p3_22_q2_22_q3_02.pdf}
    \caption{Chiral condensate (left) and renormalized chiral susceptibility (right) comparison across all four LCP ensembles using the Padé configuration: $p_2\,[2/2]$, $p_3\,[2/2]$, $q_2\,[2/2]$, and $q_3\,[0/2]$.}
    \label{fig:combined_all22_q3_02}
\end{figure}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{figures/finite_volume_effect_p2_p3.pdf}
    \caption{Finite-volume scaling comparison at $N_t = 12$ between $p_2$ ($36^3\times 12, LT=3$) and $p_3$ ($48^3\times 12, LT=4$) using Padé [2/2].}
    \label{fig:finite_volume_p2_p3}
\end{figure}

\clearpage
\section{Multi-Order Padé Fits per Ensemble}

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
    \caption{Comparison of different Padé approximation orders for each individual ensemble.}
    \label{fig:multimodal_all}
\end{figure}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.75\textwidth]{figures/chiral_susceptibility_pade02_12_q3.pdf}
    \caption{Comparison of Padé [0/2] and [1/2] Jackknife fits with shaded uncertainty bands on $q_3$ ($64^3\times 16$).}
    \label{fig:q3_02_12_comp}
\end{figure}

\clearpage
\section{Meson Mass Splitting vs Strange Quark Condensate}

\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.85\textwidth]{figures/meson_strange_condensate_fit.pdf}
    \caption{Fit of the spatial-temporal meson mass splitting $\Delta m = m(X_t) - m(T_t)$ against the strange quark condensate $\langle\bar{s}s\rangle_{\rm sub}/Z_m$, testing Linear, Power $1/3$, Power $2/3$, and Free Power functional forms.}
    \label{fig:meson_strange_fit}
\end{figure}

\begin{table}[htbp]
    \centering
    \caption{Model comparison results for $\Delta m(X_t - T_t)$ vs $\langle\bar{s}s\rangle$.}
    \label{tab:meson_models}
    \small
    \begin{tabular}{lcccccc}
        \toprule
        \textbf{Model} & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & \textbf{AIC} & \textbf{BIC} & $\boldsymbol{\Delta\mathrm{BIC}}$ & $\boldsymbol{w_{\mathrm{Bayes}}}$ \\
        \midrule
""")

    for r in meson_cmp.to_dicts():
        tex.append(f"        {r['short_name']} & {r['k_params']} & {r['chi2_dof']:.3f} & {r['AIC']:.2f} & {r['BIC']:.2f} & {r['delta_BIC']:.2f} & {r['bayes_prob_weight']:.3f} \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{Meson mass splitting and strange condensate data points.}
    \label{tab:meson_data}
    \small
    \begin{tabular}{lccccc}
        \toprule
        $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{\Delta m_{\mathrm{latt}}}$ & $\boldsymbol{\Delta m}$ [MeV] & $\boldsymbol{\langle\bar{s}s\rangle \times 10^4}$ & $\boldsymbol{\sigma(\langle\bar{s}s\rangle) \times 10^4}$ \\
        \midrule
""")

    for r in meson_data.to_dicts():
        tex.append(f"        {r['beta']:.3f} & {r['temp_mev']:.2f} & {r['delta_m_lattice']:.4f}({r['delta_m_err_lattice']:.4f}) & {r['delta_m_mev']:.1f}({r['delta_m_err_mev']:.1f}) & {r['strange_condensate']*1e4:.3f} & {r['strange_condensate_err']*1e4:.3f} \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\clearpage
\section{Padé Fit Summaries and Model Statistics}

\begin{table}[htbp]
    \centering
    \caption{Jackknife fit parameters for key Padé configurations.}
    \label{tab:lcp_summary}
    \small
    \begin{tabular}{lccccccc}
        \toprule
        \textbf{Dataset} & \textbf{Dimensions} & $\boldsymbol{LT}$ & \textbf{Order} & $\boldsymbol{T_{pc}}$ [MeV] & $\boldsymbol{\chi_{\mathrm{peak}} / Z_m^2}$ [$\text{GeV}^2$] & $\boldsymbol{\chi^2/\mathrm{dof}}$ & $\boldsymbol{T_{\mathrm{pole}}}$ [MeV] \\
        \midrule
""")

    # rows: p2 [2/2], p3 [2/2], p3 [0/2], q2 [2/2], q3 [2/2], q3 [0/2], q3 [1/2]
    tex.append(f"        $p_2$ & $36^3\\times 12$ & 3 & [2/2] & ${p2_sum['t_peak_mean']:.2f}\\pm {p2_sum['t_peak_err']:.2f}$ & ${p2_sum['chi_peak_mean']:.3f}\\pm {p2_sum['chi_peak_err']:.3f}$ & ${p2_sum['chi2_dof_mean']:.3f}$ & ${p2_sum['t_pole_real_mean']:.1f} \\pm i {p2_sum['t_pole_imag_mean']:.1f}$ \\\\")
    tex.append(f"        $p_3$ & $48^3\\times 12$ & 4 & [2/2] & ${p3_sum_22['t_peak_mean']:.2f}\\pm {p3_sum_22['t_peak_err']:.2f}$ & ${p3_sum_22['chi_peak_mean']:.3f}\\pm {p3_sum_22['chi_peak_err']:.3f}$ & ${p3_sum_22['chi2_dof_mean']:.3f}$ & ${p3_sum_22['t_pole_real_mean']:.1f} \\pm i {p3_sum_22['t_pole_imag_mean']:.1f}$ \\\\")
    tex.append(f"        $p_3$ & $48^3\\times 12$ & 4 & [0/2] & ${p3_sum_02['t_peak_mean']:.2f}\\pm {p3_sum_02['t_peak_err']:.2f}$ & ${p3_sum_02['chi_peak_mean']:.3f}\\pm {p3_sum_02['chi_peak_err']:.3f}$ & ${p3_sum_02['chi2_dof_mean']:.3f}$ & ${p3_sum_02['t_pole_real_mean']:.1f} \\pm i {p3_sum_02['t_pole_imag_mean']:.1f}$ \\\\")
    tex.append(f"        $q_2$ & $48^3\\times 16$ & 3 & [2/2] & ${q2_sum['t_peak_mean']:.2f}\\pm {q2_sum['t_peak_err']:.2f}$ & ${q2_sum['chi_peak_mean']:.3f}\\pm {q2_sum['chi_peak_err']:.3f}$ & ${q2_sum['chi2_dof_mean']:.3f}$ & ${q2_sum['t_pole_real_mean']:.1f} \\pm i {q2_sum['t_pole_imag_mean']:.1f}$ \\\\")
    tex.append(f"        $q_3$ & $64^3\\times 16$ & 4 & [2/2] & ${q3_sum_22['t_peak_mean']:.2f}\\pm {q3_sum_22['t_peak_err']:.2f}$ & ${q3_sum_22['chi_peak_mean']:.3f}\\pm {q3_sum_22['chi_peak_err']:.3f}$ & -- & ${q3_sum_22['t_pole_real_mean']:.1f} \\pm i {q3_sum_22['t_pole_imag_mean']:.1f}$ \\\\")
    tex.append(f"        $q_3$ & $64^3\\times 16$ & 4 & [0/2] & ${q3_sum_02['t_peak_mean']:.2f}\\pm {q3_sum_02['t_peak_err']:.2f}$ & ${q3_sum_02['chi_peak_mean']:.3f}\\pm {q3_sum_02['chi_peak_err']:.3f}$ & ${q3_sum_02['chi2_dof_mean']:.3f}$ & ${q3_sum_02['t_pole_real_mean']:.1f} \\pm i {q3_sum_02['t_pole_imag_mean']:.1f}$ \\\\")
    tex.append(f"        $q_3$ & $64^3\\times 16$ & 4 & [1/2] & ${q3_sum_12['t_peak_mean']:.2f}\\pm {q3_sum_12['t_peak_err']:.2f}$ & ${q3_sum_12['chi_peak_mean']:.3f}\\pm {q3_sum_12['chi_peak_err']:.3f}$ & ${q3_sum_12['chi2_dof_mean']:.3f}$ & ${q3_sum_12['t_pole_real_mean']:.1f} \\pm i {q3_sum_12['t_pole_imag_mean']:.1f}$ \\\\")

    tex.append(r"""        \bottomrule
    \end{tabular}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{Model selection statistics for $p_2$ ($36^3\times 12$) and $p_3$ ($48^3\times 12$).}
    \label{tab:model_p2_p3}
    \footnotesize
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$p_2$ ($N=10$ points)}
        \begin{tabular}{lcccccc}
            \toprule
            \textbf{Order} & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & \textbf{AIC} & \textbf{BIC} & $\boldsymbol{\Delta\mathrm{BIC}}$ & $\boldsymbol{w_{\mathrm{B}}}$ \\
            \midrule
""")

    for r in p2_cmp.to_dicts():
        tex.append(f"            {r['order']} & {r['k_params']} & {format_num(r['chi2_dof'])} & {format_num(r['AIC'], '{:.1f}')} & {format_num(r['BIC'], '{:.1f}')} & {format_num(r['delta_BIC'], '{:.1f}')} & {format_num(r['bayes_prob_weight'], '{:.3f}')} \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
    \hfill
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$p_3$ ($N=7$ points)}
        \begin{tabular}{lcccccc}
            \toprule
            \textbf{Order} & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & \textbf{AIC} & \textbf{BIC} & $\boldsymbol{\Delta\mathrm{BIC}}$ & $\boldsymbol{w_{\mathrm{B}}}$ \\
            \midrule
""")

    for r in p3_cmp.to_dicts():
        tex.append(f"            {r['order']} & {r['k_params']} & {format_num(r['chi2_dof'])} & {format_num(r['AIC'], '{:.1f}')} & {format_num(r['BIC'], '{:.1f}')} & {format_num(r['delta_BIC'], '{:.1f}')} & {format_num(r['bayes_prob_weight'], '{:.3f}')} \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{Model selection statistics for $q_2$ ($48^3\times 16$) and $q_3$ ($64^3\times 16$).}
    \label{tab:model_q2_q3}
    \footnotesize
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$q_2$ ($N=8$ points)}
        \begin{tabular}{lcccccc}
            \toprule
            \textbf{Order} & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & \textbf{AIC} & \textbf{BIC} & $\boldsymbol{\Delta\mathrm{BIC}}$ & $\boldsymbol{w_{\mathrm{B}}}$ \\
            \midrule
""")

    for r in q2_cmp.to_dicts():
        tex.append(f"            {r['order']} & {r['k_params']} & {format_num(r['chi2_dof'])} & {format_num(r['AIC'], '{:.1f}')} & {format_num(r['BIC'], '{:.1f}')} & {format_num(r['delta_BIC'], '{:.1f}')} & {format_num(r['bayes_prob_weight'], '{:.3f}')} \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
    \hfill
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$q_3$ ($N=5$ points)}
        \begin{tabular}{lcccccc}
            \toprule
            \textbf{Order} & $\boldsymbol{k}$ & $\boldsymbol{\chi^2/\mathrm{dof}}$ & \textbf{AIC} & \textbf{BIC} & $\boldsymbol{\Delta\mathrm{BIC}}$ & $\boldsymbol{w_{\mathrm{B}}}$ \\
            \midrule
""")

    for r in q3_cmp.to_dicts():
        tex.append(f"            {r['order']} & {r['k_params']} & {format_num(r['chi2_dof'])} & {format_num(r['AIC'], '{:.1f}')} & {format_num(r['BIC'], '{:.1f}')} & {format_num(r['delta_BIC'], '{:.1f}')} & {format_num(r['bayes_prob_weight'], '{:.3f}')} \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
\end{table}

\clearpage
\section{Lattice Measurements Data Tables}

\begin{table}[htbp]
    \centering
    \caption{Chiral condensate $\langle\bar{\psi}\psi\rangle_{\rm sub}/Z_m \times 10^4$ and susceptibility $\chi_{\mathrm{disc}} / Z_m^2$ [$\text{GeV}^2$] for $p_2$ and $p_3$ ($N_t=12$).}
    \label{tab:raw_p2_p3}
    \small
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$p_2$ ($36^3\times 12$)}
        \begin{tabular}{ccccc}
            \toprule
            $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle \times 10^4}$ & $\boldsymbol{\chi / Z_m^2}$ \\
            \midrule
""")

    p2_c_dict = {float(r["beta"]): (float(r["pbp_renorm"])*1e4, float(r["pbp_renorm_err"])*1e4) for r in p2_cond.to_dicts()}
    for r in p2_susc.to_dicts():
        b = float(r["beta"])
        c_val, c_err = p2_c_dict.get(b, (np.nan, np.nan))
        tex.append(f"            {b:.3f} & {float(r['temp']):.2f} & {c_val:.2f}({c_err:.2f}) & {float(r['chi_renorm']):.4f}({float(r['chi_renorm_err']):.4f}) \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
    \hfill
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$p_3$ ($48^3\times 12$)}
        \begin{tabular}{ccccc}
            \toprule
            $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle \times 10^4}$ & $\boldsymbol{\chi / Z_m^2}$ \\
            \midrule
""")

    p3_c_dict = {float(r["beta"]): (float(r["pbp_renorm"])*1e4, float(r["pbp_renorm_err"])*1e4) for r in p3_cond.to_dicts()}
    for r in p3_susc.to_dicts():
        b = float(r["beta"])
        c_val, c_err = p3_c_dict.get(b, (np.nan, np.nan))
        tex.append(f"            {b:.3f} & {float(r['temp']):.2f} & {c_val:.2f}({c_err:.2f}) & {float(r['chi_renorm']):.4f}({float(r['chi_renorm_err']):.4f}) \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
\end{table}

\begin{table}[htbp]
    \centering
    \caption{Chiral condensate $\langle\bar{\psi}\psi\rangle_{\rm sub}/Z_m \times 10^4$ and susceptibility $\chi_{\mathrm{disc}} / Z_m^2$ [$\text{GeV}^2$] for $q_2$ and $q_3$ ($N_t=16$).}
    \label{tab:raw_q2_q3}
    \small
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$q_2$ ($48^3\times 16$)}
        \begin{tabular}{ccccc}
            \toprule
            $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle \times 10^4}$ & $\boldsymbol{\chi / Z_m^2}$ \\
            \midrule
""")

    q2_c_dict = {float(r["beta"]): (float(r["pbp_renorm"])*1e4, float(r["pbp_renorm_err"])*1e4) for r in q2_cond.to_dicts()}
    for r in q2_susc.to_dicts():
        b = float(r["beta"])
        c_val, c_err = q2_c_dict.get(b, (np.nan, np.nan))
        tex.append(f"            {b:.3f} & {float(r['temp']):.2f} & {c_val:.2f}({c_err:.2f}) & {float(r['chi_renorm']):.4f}({float(r['chi_renorm_err']):.4f}) \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
    \hfill
    \begin{subtable}[t]{0.49\textwidth}
        \centering
        \caption{$q_3$ ($64^3\times 16$)}
        \begin{tabular}{ccccc}
            \toprule
            $\boldsymbol{\beta}$ & $\boldsymbol{T}$ [MeV] & $\boldsymbol{\langle\bar{\psi}\psi\rangle \times 10^4}$ & $\boldsymbol{\chi / Z_m^2}$ \\
            \midrule
""")

    q3_c_dict = {float(r["beta"]): (float(r["pbp_renorm"])*1e4, float(r["pbp_renorm_err"])*1e4) for r in q3_cond.to_dicts()}
    for r in q3_susc.to_dicts():
        b = float(r["beta"])
        c_val, c_err = q3_c_dict.get(b, (np.nan, np.nan))
        tex.append(f"            {b:.3f} & {float(r['temp']):.2f} & {c_val:.2f}({c_err:.2f}) & {float(r['chi_renorm']):.4f}({float(r['chi_renorm_err']):.4f}) \\\\")

    tex.append(r"""            \bottomrule
        \end{tabular}
    \end{subtable}
\end{table}

\end{document}
""")

    return "\n".join(tex)


def main():
    tex_content = build_latex_document()
    tex_path = DOCS_DIR / "report.tex"
    tex_path.write_text(tex_content, encoding="utf-8")
    print(f"[INFO] 成功写入 {tex_path}")

    print("[INFO] 调用 pdflatex 编译 report.tex...")
    cmd = ["/Library/TeX/texbin/pdflatex", "-interaction=nonstopmode", "report.tex"]
    proc = subprocess.run(cmd, cwd=DOCS_DIR, capture_output=True, text=True)
    if proc.returncode != 0:
        print("[ERROR] pdflatex 编译失败:")
        print(proc.stdout[-1500:])
        sys.exit(1)

    # 运行第二遍以更新引用与编号
    subprocess.run(cmd, cwd=DOCS_DIR, capture_output=True, text=True)

    pdf_path = DOCS_DIR / "report.pdf"
    if pdf_path.exists():
        print(f"[SUCCESS] 编译成功！生成目标文档: {pdf_path} (大小: {pdf_path.stat().st_size / 1024:.1f} KB)")
    else:
        print("[ERROR] 未找到输出的 report.pdf 文件！")
        sys.exit(1)


if __name__ == "__main__":
    main()
