#!/usr/bin/env python3
"""
scripts/plot_meson.py
--------------------------------------------------------------------------------
Publication-grade plotting for meson correlator observables (Polars + NumPy):
1. 42 channel comparison plots (single vs multi with data-coordinate fit bands)
2. Global effective mass thermal evolution across temperature (MeV)
3. Chiral & axial symmetry breaking / restoration signatures
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    BETAS, CHANNELS, TEMP_MAP, MULTI_FIT_SLICES, SINGLE_FIT_SLICES,
    SYMMETRY_PAIRS, TRANSITION_REGION, OUTPUT_ROOT,
    FIGURES_DIR as DOCS_FIGURES_DIR
)
DOCS_FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def normalize_beta_str(b: Any) -> str:
    """统一标准化 beta 键值为不带 '4.' 的短字符串，如 '10', '15'"""
    s = str(b).strip()
    return s[2:] if s.startswith("4.") else s


def get_temperature(beta_val: Any) -> float:
    """根据 beta 获取温度 (MeV)"""
    b_short = normalize_beta_str(beta_val)
    return TEMP_MAP.get(f"4.{b_short}", np.nan)


def load_fit_band(
    result_dir: Path,
    ch: str,
    summary_df: Optional[pl.DataFrame] = None,
    beta_str: str = ""
) -> Tuple[Optional[float], Optional[float]]:
    """
    自适应从 summary_df、summary_fit_{ch}.csv 或 fittresult{ch}.csv 提取拟合参数。
    兼容 mass_mean/mass_err 与 massfit_mean/massfit_err 两种命名规范。
    """
    b_target = normalize_beta_str(beta_str)

    # 1. 优先从预读入的 all_fits_summary.csv 检索
    if summary_df is not None and not summary_df.is_empty():
        # 清理列名首尾空格
        df_clean = summary_df.rename({c: c.strip() for c in summary_df.columns})
        if "beta" in df_clean.columns and "channel" in df_clean.columns:
            sub = df_clean.filter(
                (pl.col("beta").cast(pl.Utf8).map_elements(normalize_beta_str, return_dtype=pl.Utf8) == b_target) &
                (pl.col("channel") == ch)
            )
            if not sub.is_empty():
                m_col = next((c for c in ["mass_mean", "fit_mass", "mean", "massfit_mean"] if c in sub.columns), None)
                e_col = next((c for c in ["mass_err", "fit_err", "err", "massfit_err"] if c in sub.columns), None)
                if m_col and e_col:
                    m = sub[m_col][0]
                    e = sub[e_col][0]
                    if m is not None and np.isfinite(m):
                        return float(m), float(e)

    # 2. 从单独的 summary_fit_{ch}.csv 检索
    sum_file = result_dir / f"summary_fit_{ch}.csv"
    if sum_file.exists():
        try:
            df = pl.read_csv(sum_file)
            df = df.rename({c: c.strip() for c in df.columns})
            # 如果是 key-value 行结构 (quantity, mean, jack_err)
            if "quantity" in df.columns and "mean" in df.columns:
                row = df.filter(pl.col("quantity").is_in(["mass", "massfit_mean", "fit_mass"]))
                if not row.is_empty():
                    err_col = "jack_err" if "jack_err" in row.columns else "err"
                    return float(row["mean"][0]), float(row[err_col][0])
            # 如果是单行宽表结构 (chi2, dof, massfit_mean, massfit_err, ...)
            elif "massfit_mean" in df.columns and "massfit_err" in df.columns:
                return float(df["massfit_mean"][0]), float(df["massfit_err"][0])
            elif "mass_mean" in df.columns and "mass_err" in df.columns:
                return float(df["mass_mean"][0]), float(df["mass_err"][0])
        except Exception:
            pass

    # 3. 从单样本拟合汇总 fittresult{ch}.csv 计算 Jackknife 估计
    fit_file = result_dir / f"fittresult{ch}.csv"
    if fit_file.exists():
        try:
            df = pl.read_csv(fit_file)
            df = df.rename({c: c.strip() for c in df.columns})
            target_col = next((c for c in ["massfit_mean", "mass_mean", "mass"] if c in df.columns), None)
            if target_col and len(df) > 0:
                m_vals = df[target_col].to_numpy().astype(np.float64)
                m_vals = m_vals[np.isfinite(m_vals)]
                n = len(m_vals)
                if n > 1:
                    m_mean = float(np.mean(m_vals))
                    jack_err = float(np.sqrt((n - 1) / n * np.sum((m_vals - m_mean) ** 2)))
                    return m_mean, jack_err
        except Exception:
            pass

    return None, None


def plot_channel_comparisons() -> None:
    """生成 42 组物理信道有效质量对比图 (单源 vs 多源 + 平台拟合带)"""
    comp_dir = DOCS_FIGURES_DIR / "channel_comparisons"
    comp_dir.mkdir(parents=True, exist_ok=True)

    s_summary_path = OUTPUT_ROOT / "simulateresult-singlesrc" / "all_fits_summary.csv"
    m_summary_path = OUTPUT_ROOT / "simulateresult" / "all_fits_summary.csv"
    s_summary = pl.read_csv(s_summary_path) if s_summary_path.exists() else None
    m_summary = pl.read_csv(m_summary_path) if m_summary_path.exists() else None

    for beta in BETAS:
        b_short = normalize_beta_str(beta)
        beta_dir = comp_dir / f"b4.{b_short}"
        beta_dir.mkdir(parents=True, exist_ok=True)

        for ch in CHANNELS:
            s_meff_file = OUTPUT_ROOT / "ratio_results-singlesrc" / f"b4.{b_short}" / f"meff_{ch}.csv"
            m_meff_file = OUTPUT_ROOT / "ratio_results" / f"b4.{b_short}" / f"meff_{ch}.csv"

            if not (s_meff_file.exists() and m_meff_file.exists()):
                continue

            s_df = pl.read_csv(s_meff_file)
            m_df = pl.read_csv(m_meff_file)

            # 只取前 24 个物理有效时间片
            s_mean = s_df['mean'].head(24).to_numpy()
            s_err = s_df['err'].head(24).to_numpy()
            m_mean = m_df['mean'].head(24).to_numpy()
            m_err = m_df['err'].head(24).to_numpy()

            t_axis = np.arange(len(s_mean))
            s_t_start = SINGLE_FIT_SLICES[b_short][ch]
            m_t_start = MULTI_FIT_SLICES[b_short][ch]
            t_fit_end = 23  # 拟合切片右端点数据坐标

            fig, ax = plt.subplots(figsize=(8, 5))

            # 绘制单源与多源散点及误差棒
            ax.errorbar(
                t_axis, s_mean, yerr=s_err,
                fmt='o', color='#d62728', ecolor='#d62728', elinewidth=1.2, capsize=2.5,
                markersize=5, label='Single-Source', zorder=4
            )
            ax.errorbar(
                t_axis, m_mean, yerr=m_err,
                fmt='s', color='#1f77b4', ecolor='#1f77b4', elinewidth=1.2, capsize=2.5,
                markersize=5, label='Multi-Source (16 srcs)', zorder=3
            )

            # 拟合区间数据坐标填充 (彻底替代有坐标系偏差的 axhspan)
            s_val, s_err_val = load_fit_band(
                OUTPUT_ROOT / "simulateresult-singlesrc" / f"b4.{b_short}",
                ch, s_summary, b_short
            )
            if s_val is not None and s_err_val is not None and np.isfinite(s_val):
                t_band = np.array([s_t_start, t_fit_end], dtype=np.float64)
                ax.fill_between(
                    t_band, s_val - s_err_val, s_val + s_err_val,
                    color='#d62728', alpha=0.25, zorder=1,
                    label=f'Single Fit [{s_t_start}, {t_fit_end}]'
                )
                ax.hlines(s_val, s_t_start, t_fit_end, colors='#d62728', linestyles='--', linewidth=1.2, zorder=2)

            m_val, m_err_val = load_fit_band(
                OUTPUT_ROOT / "simulateresult" / f"b4.{b_short}",
                ch, m_summary, b_short
            )
            if m_val is not None and m_err_val is not None and np.isfinite(m_val):
                t_band = np.array([m_t_start, t_fit_end], dtype=np.float64)
                ax.fill_between(
                    t_band, m_val - m_err_val, m_val + m_err_val,
                    color='#1f77b4', alpha=0.25, zorder=1,
                    label=f'Multi Fit [{m_t_start}, {t_fit_end}]'
                )
                ax.hlines(m_val, m_t_start, t_fit_end, colors='#1f77b4', linestyles='--', linewidth=1.2, zorder=2)

            ax.set_ylim([0.0, 2.0])
            ax.set_xlim([-0.5, 24.5])
            ax.set_xlabel(r"$t / a$", fontsize=14)
            ax.set_ylabel(r"$m_{\mathrm{eff}} \cdot a$", fontsize=14)
            t_mev = get_temperature(b_short)
            ax.set_title(f"Channel {ch} — $\\beta=4.{b_short}$ ($T={t_mev:.1f}\\ \\mathrm{{MeV}}$)", fontsize=13)
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.grid(True, linestyle="--", alpha=0.5)
            ax.legend(loc="upper right", frameon=True, fontsize=9)

            fig.tight_layout()
            fig.savefig(beta_dir / f"ratio{ch}.png", dpi=300)
            plt.close(fig)

    print(f"[INFO] 42 组信道对比图已保存至: {comp_dir}")


def plot_simulation_results() -> None:
    """绘制所有强子信道在不同温度下的热演化曲线"""
    fit_summary_file = OUTPUT_ROOT / "simulateresult" / "all_fits_summary.csv"
    if not fit_summary_file.exists():
        print(f"[WARN] 找不到拟合总结文件: {fit_summary_file}")
        return

    df = pl.read_csv(fit_summary_file)
    df = df.rename({c: c.strip() for c in df.columns})

    # 标准化列名映射
    m_col = next((c for c in ["mass_mean", "fit_mass", "massfit_mean"] if c in df.columns), None)
    e_col = next((c for c in ["mass_err", "fit_err", "massfit_err"] if c in df.columns), None)

    if not m_col or not e_col:
        print(f"[ERROR] 无法解析拟合总结字段，当前列为: {df.columns}")
        return

    channel_styles = {
        'AV':  {'color': '#1f77b4', 'marker': 'o', 'label': r'$AV$ ($\bar{\psi}\gamma_5\gamma_\mu\psi$)'},
        'S':   {'color': '#ff7f0e', 'marker': 's', 'label': r'$S$ ($\bar{\psi}\psi$)'},
        'Tt':  {'color': '#2ca02c', 'marker': '^', 'label': r'$T_t$ ($\bar{\psi}\sigma_{0\mu}\psi$)'},
        'PS':  {'color': '#d62728', 'marker': 'D', 'label': r'$PS$ ($\bar{\psi}\gamma_5\psi$)'},
        'Xt':  {'color': '#9467bd', 'marker': 'v', 'label': r'$X_t$ ($\bar{\psi}\gamma_5\sigma_{0\mu}\psi$)'},
        'Vec': {'color': '#8c564b', 'marker': 'p', 'label': r'$Vec$ ($\bar{\psi}\gamma_\mu\psi$)'},
    }

    # 提取并排序温度
    df = df.with_columns([
        pl.col("beta").cast(pl.Utf8).map_elements(get_temperature, return_dtype=pl.Float64).alias("T"),
        pl.col(m_col).cast(pl.Float64).alias("plot_mass"),
        pl.col(e_col).cast(pl.Float64).alias("plot_err"),
    ])

    fig, ax = plt.subplots(figsize=(9, 6))
    for ch in CHANNELS:
        sub = df.filter((pl.col("channel") == ch) & pl.col("T").is_not_nan()).sort("T")
        if sub.is_empty():
            continue

        st = channel_styles.get(ch, {'color': 'black', 'marker': 'o', 'label': ch})
        ax.errorbar(
            sub['T'].to_numpy(), sub['plot_mass'].to_numpy(), yerr=sub['plot_err'].to_numpy(),
            fmt=st['marker'], color=st['color'], ecolor=st['color'],
            elinewidth=1.5, capsize=3, capthick=1.5,
            label=st['label'], markersize=7, zorder=5
        )

    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color='gray', alpha=0.18, label=f'Transition Region ({t_low:.1f}-{t_high:.1f} MeV)')
    ax.set_xlabel(r'$T$ [MeV]', fontsize=14)
    ax.set_ylabel(r'$m_{\mathrm{eff}} \cdot a$', fontsize=14)
    ax.set_title(r'Thermal Evolution of Meson Masses ($48^3 \times 16$ DWF)', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='best', frameon=True, fontsize=10)

    fig.tight_layout()
    out_png = DOCS_FIGURES_DIR / "simulation_results.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"[INFO] 介子质量热演化图已保存: {out_png}")


def plot_symmetry_breaking() -> None:
    """绘制宇称与手征伴侣的质量简并劈裂曲线（表征对称性恢复）"""
    fit_summary_file = OUTPUT_ROOT / "simulateresult" / "all_fits_summary.csv"
    if not fit_summary_file.exists():
        return

    df = pl.read_csv(fit_summary_file)
    df = df.rename({c: c.strip() for c in df.columns})

    m_col = next((c for c in ["mass_mean", "fit_mass", "massfit_mean"] if c in df.columns), None)
    e_col = next((c for c in ["mass_err", "fit_err", "massfit_err"] if c in df.columns), None)

    if not m_col or not e_col:
        return

    df = df.with_columns([
        pl.col("beta").cast(pl.Utf8).map_elements(normalize_beta_str, return_dtype=pl.Utf8).alias("beta_norm"),
        pl.col(m_col).cast(pl.Float64).alias("mass"),
        pl.col(e_col).cast(pl.Float64).alias("err"),
    ])

    fig, ax = plt.subplots(figsize=(8, 5))

    for pair in SYMMETRY_PAIRS:
        ch1, ch2 = pair['ch1'], pair['ch2']
        sub1 = df.filter(pl.col("channel") == ch1).select(["beta_norm", "mass", "err"])
        sub2 = df.filter(pl.col("channel") == ch2).select(["beta_norm", "mass", "err"])

        joined = sub1.join(sub2, on="beta_norm", suffix="_2")
        if joined.is_empty():
            continue

        joined = joined.with_columns(
            pl.col("beta_norm").map_elements(get_temperature, return_dtype=pl.Float64).alias("T")
        ).sort("T")

        t_vals = joined["T"].to_numpy()
        diff_vals = (joined["mass"] - joined["mass_2"]).abs().to_numpy()
        # 误差传播 (无协方差时保守采用方和根)
        diff_errs = np.sqrt(joined["err"].to_numpy()**2 + joined["err_2"].to_numpy()**2)

        ax.errorbar(
            t_vals, diff_vals, yerr=diff_errs,
            fmt=pair['marker'], color=pair['color'], ecolor=pair['color'],
            elinewidth=1.5, capsize=3, capthick=1.5,
            label=pair['label'], markersize=7, zorder=5
        )

    t_low, t_high = TRANSITION_REGION
    ax.axvspan(t_low, t_high, color='gray', alpha=0.15, label='Transition Region')
    ax.axhline(0, color='black', linestyle=':', alpha=0.6)
    ax.set_xlabel(r'$T$ [MeV]', fontsize=14)
    ax.set_ylabel(r'$|\Delta m_{\mathrm{eff}}| \cdot a$', fontsize=14)
    ax.set_title(r'Chiral and Axial Symmetry Restoration Signatures', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='best', frameon=True, fontsize=10)

    fig.tight_layout()
    out_png = DOCS_FIGURES_DIR / "symmetry_breaking_ana_setup.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"[INFO] 手征对称性破缺与恢复图已保存: {out_png}")


def main() -> None:
    print("==========================================================================")
    print("[INFO] 开始生成介子观测物理图谱 (Polars + NumPy Engine)...")
    print("==========================================================================")
    plot_channel_comparisons()
    plot_simulation_results()
    plot_symmetry_breaking()
    print("[INFO] 所有物理图谱生成完毕。")


if __name__ == "__main__":
    main()
