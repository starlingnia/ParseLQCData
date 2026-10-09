#!/usr/bin/env python3
"""
scripts/plot_effective_mass_b417.py
--------------------------------------------------------------------------------
为所有 beta=4.17 的目录绘制 effective mass (有效质量) 图，辅助确定
config/fit_windows.txt 的具体拟合区间 [x_start, x_end)。

功能与产物：
1. 按目录保存 (docs/effective_mass/by_directory/<dir_name>/):
     - 单信道高清图像 (meff_Vec.png, meff_AV.png, meff_Tt.png, meff_Xt.png, meff_S.png, meff_PS.png)
     - 该目录下 6 信道全景组合图 (meff_all_channels.png)
     - 该目录下所有信道的有效质量数据表格 (meff_data.csv)
2. 按频道保存 (docs/effective_mass/by_channel/<channel>/):
     - 各格点规模下 4 种轻夸克质量 (ml=0.0020, 0.0035, 0.0070, 0.0120) 的叠加对比图 (<geom>_compare_ml.png)
3. 格点规模全景概览 (docs/effective_mass/ensemble_overview/):
     - 各格点 6 信道 x 4 质量的大画幅概览卡片 (overview_<geom>.png)
4. 综合索引报告 (docs/effective_mass/README.md):
     - 列出所有目录和信道的文件路径、当前 config/fit_windows.txt 窗口设定、观测平台区间与建议。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.parselqcdata.ccor_flow import FIT_WINDOW_FILE, load_fit_windows

# 物理常数与标度 (153 * 16 = 2448 MeV)
SCALE_UNIT: float = 153.0 * 16.0  # 2448.0 MeV

# 6 个标准信道及其显示名称与配色
CHANNELS: Tuple[str, ...] = ("Vec", "AV", "PS", "Tt", "Xt", "S")

CHANNEL_TITLES: Dict[str, str] = {
    "Vec": "Vector (Vec / V, ρ/ω)",
    "AV": "Axial-Vector (AV / A, a₁)",
    "PS": "Pseudo-Scalar (PS / Ps, π/η)",
    "Tt": "Tensor (Tt / T, b₁)",
    "Xt": "Axial-Tensor (Xt / X, h₁)",
    "S": "Scalar (S, σ/a₀)",
}

CHANNEL_COLORS: Dict[str, str] = {
    "Vec": "#1f77b4",  # blue
    "AV": "#2ca02c",   # green
    "PS": "#ff7f0e",   # orange
    "Tt": "#d62728",   # red
    "Xt": "#9467bd",   # purple
    "S": "#8c564b",    # brown
}

# 4 种轻夸克质量及其对比配色
ML_COLORS: Dict[float, str] = {
    0.0020: "#9400D3",  # DarkViolet
    0.0035: "#1E90FF",  # DodgerBlue
    0.0070: "#2E8B57",  # SeaGreen
    0.0120: "#FF8C00",  # DarkOrange
}

# 信道别名映射 (用于对接 fit_windows.txt)
CHANNEL_ALIASES: Dict[str, Tuple[str, ...]] = {
    "Vec": ("Vec", "V"),
    "AV": ("AV", "A"),
    "PS": ("PS", "Ps", "P"),
    "Tt": ("Tt", "T"),
    "Xt": ("Xt", "X"),
    "S": ("S",),
}


@dataclass(frozen=True)
class DirectorySpec:
    dir_name: str         # 对应 data/readin/ 下的目录名
    case_key: str         # 如 "32x12", "36x18", "48x16"
    ns: int
    nt: int
    ml: float
    ms: float
    temp_mev: float
    data_source_path: Path
    multisrc_case_dir: Optional[Path] = None

    @property
    def half(self) -> int:
        return self.ns // 2


def build_directory_list() -> List[DirectorySpec]:
    """收集全部 25 个 beta=4.17 目录信息"""
    specs: List[DirectorySpec] = []

    # 1. 24 个 Nt 扫描任务目录 (32x12, 32x14, 32x16, 36x18, 40x16, 48x18)
    geoms = [
        ("32x12", 32, 12, "00_32x12"),
        ("32x14", 32, 14, "01_32x14"),
        ("32x16", 32, 16, "02_32x16"),
        ("36x18", 36, 18, "03_36x18"),
        ("40x16", 40, 16, "04_40x16"),
        ("48x18", 48, 18, "05_48x18"),
    ]
    mls = (0.0020, 0.0035, 0.0070, 0.0120)

    for case_key, ns, nt, prefix in geoms:
        temp = SCALE_UNIT / nt  # 2448 / nt MeV
        for ml in mls:
            dir_name = f"{case_key}_b4.17_ms0.040m{ml:.4f}"
            case_folder = f"{prefix}_ms0.040_m{ml:.4f}"
            multisrc_p = PROJECT_ROOT / "output" / "meson_scan" / "b4.17" / "multisrc" / "cases" / case_folder

            specs.append(
                DirectorySpec(
                    dir_name=dir_name,
                    case_key=case_key,
                    ns=ns,
                    nt=nt,
                    ml=ml,
                    ms=0.040,
                    temp_mev=temp,
                    data_source_path=multisrc_p,
                    multisrc_case_dir=multisrc_p,
                )
            )

    # 2. 48x16b4.17 目录 (定标度温度扫描基准点)
    ratio_p = PROJECT_ROOT / "output" / "ratio_results" / "b4.17"
    specs.append(
        DirectorySpec(
            dir_name="48x16b4.17",
            case_key="48x16",
            ns=48,
            nt=16,
            ml=0.001001,
            ms=0.038400,
            temp_mev=153.31,
            data_source_path=ratio_p,
            multisrc_case_dir=None,
        )
    )

    return specs


def get_current_window(windows, case_key: str, channel: str) -> Optional[Tuple[int, int]]:
    """从 fit_windows.txt 查询当前配置的拟合窗口"""
    aliases = CHANNEL_ALIASES.get(channel, (channel,))
    for alias in aliases:
        win = windows.get(case_key, alias)
        if win is not None:
            return win
    return None


def load_channel_meff(spec: DirectorySpec, channel: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    加载指定信道的有效质量数据:
    返回 (x_array, mean_array, err_array)
    """
    p = spec.data_source_path / f"meff_{channel}.csv"
    if not p.exists():
        return np.array([]), np.array([]), np.array([])

    df = pl.read_csv(p)
    if "x" in df.columns:
        x = df["x"].to_numpy().astype(np.float64)
    else:
        x = np.arange(len(df), dtype=np.float64)

    mean = df["mean"].to_numpy().astype(np.float64)
    err = df["err"].to_numpy().astype(np.float64)
    return x, mean, err


def load_summary_fit_dict() -> Dict[Tuple[str, float, str], dict]:
    """从 meson_mass_summary.csv 加载已拟合质量和窗口元数据"""
    summary_path = (
        PROJECT_ROOT
        / "output"
        / "meson_scan"
        / "b4.17"
        / "multisrc"
        / "meson_mass_summary.csv"
    )
    fit_dict: Dict[Tuple[str, float, str], dict] = {}
    if not summary_path.exists():
        return fit_dict

    df = pl.read_csv(summary_path)
    for r in df.iter_rows(named=True):
        key = (str(r["case_key"]), round(float(r["ml"]), 4), str(r["channel"]))
        fit_dict[key] = r
    return fit_dict


def plot_single_channel(
    spec: DirectorySpec,
    channel: str,
    x: np.ndarray,
    mean: np.ndarray,
    err: np.ndarray,
    window: Optional[Tuple[int, int]],
    summary_row: Optional[dict],
    out_path: Path,
) -> None:
    """绘制单个信道的详细有效质量图"""
    half = spec.half
    valid = (x >= 0) & (x <= half) & np.isfinite(mean) & np.isfinite(err)

    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)

    # 1. 绘制有效质量散点与误差棒
    x_val = x[valid]
    y_val = mean[valid]
    y_err = err[valid]

    color = CHANNEL_COLORS.get(channel, "#1f77b4")
    ax.errorbar(
        x_val,
        y_val,
        yerr=y_err,
        fmt="o",
        markersize=5,
        capsize=3,
        color=color,
        ecolor=color,
        elinewidth=1.2,
        label=r"$m_{\mathrm{eff}}(x)$ data",
        zorder=5,
    )

    # 2. 如果当前在 config/fit_windows.txt 中定义了窗口，画出阴影区间
    if window is not None:
        w_start, w_end = window
        ax.axvspan(
            w_start,
            w_end,
            color="#ffa726",
            alpha=0.22,
            label=f"fit_windows.txt: [{w_start}, {w_end})",
            zorder=2,
        )

        # 若存在已拟合的 cosh 质量，画出水平拟合带
        if summary_row and np.isfinite(summary_row.get("mass", np.nan)):
            mass = float(summary_row["mass"])
            mass_err = float(summary_row["mass_err"])
            fit_x = np.arange(w_start, w_end + 1)
            ax.hlines(
                mass,
                w_start,
                w_end,
                color="#e65100",
                linewidth=2.0,
                label=f"Fit $m = {mass:.5f} \\pm {mass_err:.5f}$",
                zorder=6,
            )
            ax.axhspan(
                mass - mass_err,
                mass + mass_err,
                color="#ffa726",
                alpha=0.25,
                zorder=3,
            )

    # 3. 自动检测平台标注 (如有)
    if summary_row and summary_row.get("auto_x_start", -1) >= 0:
        a_start = int(summary_row["auto_x_start"])
        a_end = int(summary_row["auto_x_end"])
        ax.axvspan(
            a_start,
            a_end,
            facecolor="#26a69a",
            edgecolor="#00695c",
            alpha=0.10,
            linestyle="--",
            label=f"Auto detected: [{a_start}, {a_end})",
            zorder=1,
        )

    # 4. 对称点 Ns/2 垂线
    ax.axvline(half, color="#78909c", linestyle=":", linewidth=1.2, label=f"Symmetry point $N_s/2={half}$")

    # 5. 设置坐标范围与刻度网格
    ax.set_xlim(-0.5, half + 0.8)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))

    # 限制 y 轴范围，防止极少奇异点/巨大误差棒撑坏平台视野
    # 选取 x >= 2 的有效数据估计 y 轴范围
    plateau_candidates = y_val[(x_val >= 2) & (y_err < 1.0)]
    if len(plateau_candidates) > 0:
        p_min = np.min(plateau_candidates)
        p_max = np.max(plateau_candidates)
        y_top = min(2.0, max(1.2, p_max * 1.35))
        y_bottom = max(-0.05, p_min - 0.25)
        ax.set_ylim(y_bottom, y_top)
    else:
        ax.set_ylim(-0.05, 2.0)

    ax.grid(True, which="major", color="#e0e0e0", linestyle="-", linewidth=0.8)
    ax.grid(True, which="minor", color="#f5f5f5", linestyle=":", linewidth=0.5)

    # 6. 标签与标题
    ch_title = CHANNEL_TITLES.get(channel, channel)
    ax.set_title(
        f"Effective Mass: {spec.dir_name} | {ch_title}",
        fontsize=12,
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel("Spatial separation $x$ [lattice units]", fontsize=11)
    ax.set_ylabel(r"Effective mass $a m_{\mathrm{eff}}(x)$ [lattice units]", fontsize=11)

    # 7. 信息浮窗卡片
    info_lines = [
        f"Ensemble: {spec.case_key} ($N_s={spec.ns}, N_t={spec.nt}$)",
        f"$m_l = {spec.ml:.4f}$, $m_s = {spec.ms:.3f}$",
        f"Temp $T = {spec.temp_mev:.2f}$ MeV",
        f"Symmetry point $half = {half}$",
    ]
    if window is not None:
        info_lines.append(f"Window: [{window[0]}, {window[1]})")
    else:
        info_lines.append("Window: [Not set in fit_windows.txt]")

    ax.text(
        0.03,
        0.96,
        "\n".join(info_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "alpha": 0.88, "edgecolor": "#cccccc"},
        zorder=10,
    )

    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_all_channels_grid(
    spec: DirectorySpec,
    meff_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    windows,
    summary_dict: Dict[Tuple[str, float, str], dict],
    out_path: Path,
) -> None:
    """绘制该目录 6 信道 2x3 全景对比图"""
    half = spec.half
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.5), constrained_layout=True)
    channel_grid = [
        ["Vec", "AV", "PS"],
        ["Tt", "Xt", "S"],
    ]

    for r in range(2):
        for c in range(3):
            ch = channel_grid[r][c]
            ax = axes[r, c]
            x, mean, err = meff_data.get(ch, (np.array([]), np.array([]), np.array([])))
            window = get_current_window(windows, spec.case_key, ch)
            summary_row = summary_dict.get((spec.case_key, round(spec.ml, 4), ch))

            if len(x) > 0:
                valid = (x >= 0) & (x <= half) & np.isfinite(mean) & np.isfinite(err)
                x_val = x[valid]
                y_val = mean[valid]
                y_err = err[valid]

                color = CHANNEL_COLORS.get(ch, "#1f77b4")
                ax.errorbar(
                    x_val,
                    y_val,
                    yerr=y_err,
                    fmt="o",
                    markersize=4,
                    capsize=2,
                    color=color,
                    ecolor=color,
                    label=r"$m_{\mathrm{eff}}$",
                    zorder=5,
                )

                if window is not None:
                    ax.axvspan(
                        window[0],
                        window[1],
                        color="#ffa726",
                        alpha=0.25,
                        label=f"win: [{window[0]},{window[1]})",
                        zorder=2,
                    )
                    if summary_row and np.isfinite(summary_row.get("mass", np.nan)):
                        mass = float(summary_row["mass"])
                        ax.hlines(mass, window[0], window[1], color="#e65100", linewidth=1.8, zorder=6)

                ax.axvline(half, color="#90a4ae", linestyle=":", linewidth=1.0)
                ax.set_xlim(-0.5, half + 0.8)
                ax.xaxis.set_major_locator(MultipleLocator(2))
                ax.xaxis.set_minor_locator(MultipleLocator(1))

                plateau_candidates = y_val[(x_val >= 2) & (y_err < 1.0)]
                if len(plateau_candidates) > 0:
                    p_min = np.min(plateau_candidates)
                    p_max = np.max(plateau_candidates)
                    y_top = min(2.0, max(1.1, p_max * 1.3))
                    y_bottom = max(-0.05, p_min - 0.2)
                    ax.set_ylim(y_bottom, y_top)
                else:
                    ax.set_ylim(-0.05, 1.8)

            ax.grid(True, which="major", color="#e0e0e0", linestyle="-", linewidth=0.6)
            ax.set_title(CHANNEL_TITLES.get(ch, ch), fontsize=10.5, fontweight="bold")
            ax.set_xlabel("Spatial separation $x$", fontsize=9.5)
            ax.set_ylabel(r"$a m_{\mathrm{eff}}(x)$", fontsize=9.5)
            ax.legend(loc="upper right", fontsize=8.5, framealpha=0.85)

    fig.suptitle(
        f"6-Channel Effective Mass Overview: {spec.dir_name} (T = {spec.temp_mev:.2f} MeV, $m_l = {spec.ml:.4f}$)",
        fontsize=13,
        fontweight="bold",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_by_channel_ml_comparison(
    geom_specs: List[DirectorySpec],
    channel: str,
    windows,
    out_path: Path,
) -> None:
    """按信道比较同一格点几何下 4 种轻夸克质量的有效质量平台 (辅助挑选对所有 ml 都适用的窗口)"""
    fig, ax = plt.subplots(figsize=(8.5, 5.4), constrained_layout=True)

    geom_key = geom_specs[0].case_key
    half = geom_specs[0].half
    window = get_current_window(windows, geom_key, channel)

    # 4 种质量微小横向抖动偏移，防止重叠
    jitter_map = {0.0020: -0.18, 0.0035: -0.06, 0.0070: 0.06, 0.0120: 0.18}

    y_candidates = []
    for spec in sorted(geom_specs, key=lambda s: s.ml):
        x, mean, err = load_channel_meff(spec, channel)
        if len(x) == 0:
            continue
        valid = (x >= 0) & (x <= half) & np.isfinite(mean) & np.isfinite(err)
        x_val = x[valid] + jitter_map.get(spec.ml, 0.0)
        y_val = mean[valid]
        y_err = err[valid]

        color = ML_COLORS.get(spec.ml, "#333333")
        ax.errorbar(
            x_val,
            y_val,
            yerr=y_err,
            fmt="o",
            markersize=4.5,
            capsize=2.5,
            color=color,
            ecolor=color,
            elinewidth=1.0,
            label=f"$m_l = {spec.ml:.4f}$",
            zorder=5,
        )
        cand = mean[(x >= 2) & (x <= half) & (err < 0.8)]
        if len(cand) > 0:
            y_candidates.extend(cand.tolist())

    if window is not None:
        ax.axvspan(
            window[0],
            window[1],
            color="#ffa726",
            alpha=0.22,
            label=f"fit_windows.txt: [{window[0]}, {window[1]})",
            zorder=2,
        )

    ax.axvline(half, color="#90a4ae", linestyle=":", linewidth=1.2, label=f"half={half}")
    ax.set_xlim(-0.5, half + 0.8)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))

    if y_candidates:
        p_min = min(y_candidates)
        p_max = max(y_candidates)
        ax.set_ylim(max(-0.05, p_min - 0.25), min(2.0, max(1.2, p_max * 1.35)))
    else:
        ax.set_ylim(-0.05, 2.0)

    ax.grid(True, which="major", color="#e0e0e0", linestyle="-", linewidth=0.8)
    ax.grid(True, which="minor", color="#f5f5f5", linestyle=":", linewidth=0.5)

    ch_title = CHANNEL_TITLES.get(channel, channel)
    ax.set_title(
        f"Quark Mass Dependence: {geom_key} | {ch_title} (T = {geom_specs[0].temp_mev:.2f} MeV)",
        fontsize=12,
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel("Spatial separation $x$", fontsize=11)
    ax.set_ylabel(r"$a m_{\mathrm{eff}}(x)$", fontsize=11)
    ax.legend(loc="upper right", fontsize=9.5, framealpha=0.9)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_ensemble_overview_card(
    geom_specs: List[DirectorySpec],
    windows,
    out_path: Path,
) -> None:
    """按格点几何绘制 6 信道 x 4 质量的大画幅总览卡片 (overview_<geom>.png)"""
    geom_key = geom_specs[0].case_key
    half = geom_specs[0].half

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 9.0), constrained_layout=True)
    channel_grid = [
        ["Vec", "AV", "PS"],
        ["Tt", "Xt", "S"],
    ]

    jitter_map = {0.0020: -0.18, 0.0035: -0.06, 0.0070: 0.06, 0.0120: 0.18}

    for r in range(2):
        for c in range(3):
            ch = channel_grid[r][c]
            ax = axes[r, c]
            window = get_current_window(windows, geom_key, ch)

            y_cands = []
            for spec in sorted(geom_specs, key=lambda s: s.ml):
                x, mean, err = load_channel_meff(spec, ch)
                if len(x) == 0:
                    continue
                valid = (x >= 0) & (x <= half) & np.isfinite(mean) & np.isfinite(err)
                x_val = x[valid] + jitter_map.get(spec.ml, 0.0)
                y_val = mean[valid]
                y_err = err[valid]

                color = ML_COLORS.get(spec.ml, "#333333")
                ax.errorbar(
                    x_val,
                    y_val,
                    yerr=y_err,
                    fmt="o",
                    markersize=3.5,
                    capsize=2,
                    color=color,
                    ecolor=color,
                    label=f"m={spec.ml}",
                    zorder=5,
                )
                c_valid = mean[(x >= 2) & (x <= half) & (err < 0.8)]
                if len(c_valid) > 0:
                    y_cands.extend(c_valid.tolist())

            if window is not None:
                ax.axvspan(
                    window[0],
                    window[1],
                    color="#ffa726",
                    alpha=0.25,
                    label=f"win: [{window[0]},{window[1]})",
                    zorder=2,
                )

            ax.axvline(half, color="#90a4ae", linestyle=":", linewidth=1.0)
            ax.set_xlim(-0.5, half + 0.8)
            ax.xaxis.set_major_locator(MultipleLocator(2))
            ax.xaxis.set_minor_locator(MultipleLocator(1))

            if y_cands:
                p_min = min(y_cands)
                p_max = max(y_cands)
                ax.set_ylim(max(-0.05, p_min - 0.2), min(2.0, max(1.1, p_max * 1.3)))
            else:
                ax.set_ylim(-0.05, 1.8)

            ax.grid(True, which="major", color="#e0e0e0", linestyle="-", linewidth=0.6)
            ax.set_title(CHANNEL_TITLES.get(ch, ch), fontsize=10.5, fontweight="bold")
            ax.set_xlabel("Spatial separation $x$", fontsize=9.5)
            ax.set_ylabel(r"$a m_{\mathrm{eff}}(x)$", fontsize=9.5)
            ax.legend(loc="upper right", fontsize=8.0, framealpha=0.85)

    fig.suptitle(
        f"Ensemble Overview: {geom_key} ($N_s={geom_specs[0].ns}, N_t={geom_specs[0].nt}$, T = {geom_specs[0].temp_mev:.2f} MeV) Across All 4 Quark Masses",
        fontsize=13,
        fontweight="bold",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def generate_readme_report(
    specs: List[DirectorySpec],
    windows,
    summary_dict: Dict[Tuple[str, float, str], dict],
    base_out_dir: Path,
) -> Path:
    """生成 docs/effective_mass/README.md 综合汇总报告与图片快速导航"""
    readme_path = base_out_dir / "README.md"

    lines = [
        "# Beta 4.17 介子 Effective Mass (有效质量) 全局分析报告与拟合窗口导航",
        "",
        "> **核心目标**：提供所有 beta=4.17 目录下 6 个介子信道（Vec, AV, PS, Tt, Xt, S）的有效质量平台图，作为精确设定与校准 `config/fit_windows.txt` 中各 ensemble 及噪声信道窗口 `[x_start, x_end)` 的判定依据。",
        ">",
        "> **物理标度**：$a^{-1} = 153 \\times 16 = 2448\\text{ MeV}$，$T = 2448 / N_t\\text{ MeV}$。",
        "",
        "## 1. 当前 `config/fit_windows.txt` 配置现状与观察建议",
        "",
        "| Ensemble | Ns x Nt | T (MeV) | 默认窗口 | 噪声道专用窗口 | 观测有效质量平台特征与建议 |",
        "|:---|:---:|:---:|:---:|:---|:---|",
    ]

    # 总结各格点特征
    obs_summary = {
        "32x12": ("`12 17`", "`S: 6 17`", "平台在 x=12~16 表现优良；S 道噪声从 x=7 起增大，故 S 道采用 [6, 17) 可保持稳定。"),
        "32x14": ("`12 17`", "`S: 4 17`", "Vec/AV/PS 在 [12, 17) 稳定；S 道在近相变区信号衰减极快，放宽至 [4, 17) 避免 m→0 退化。"),
        "32x16": ("`12 17`", "无", "处于赝临界点 T=153 MeV 附近，Vec/AV 在 [12, 17) 形成良好平台；S 道涨落显著。"),
        "36x18": ("**待定 (当前注释)**", "**待定**", "Ns=36 对称点 half=18。Vec/AV/PS 在 [12, 19) 或 [14, 19) 展现清晰平坦基态；S/Tt 道需注意晚期涨落。"),
        "40x16": ("`16 21`", "`Tt: 8 21`, `Xt: 8 21`", "Ns=40 对称点 half=20。默认区间 [16, 21) 对 Vec/AV 极佳；Tt/Xt 在小 ml 处早衰，放宽至 [8, 21)。"),
        "48x18": ("`16 25`", "`Tt/Xt: 8 25`, `S: 4 25`", "Ns=48 对称点 half=24。默认 [16, 25) 极平坦；张量与标量道在对称点附近信噪比过低，专用窗口极为关键。"),
        "48x16": ("待设定", "待设定", "定标度物理点基准 (T=153.31 MeV)，对称点 half=24。"),
    }

    seen_geoms = set()
    for s in specs:
        geom = s.case_key
        if geom in seen_geoms:
            continue
        seen_geoms.add(geom)

        def_win = windows.defaults.get(geom, "未设定")
        def_win_str = f"`{def_win[0]} {def_win[1]}`" if isinstance(def_win, tuple) else str(def_win)

        ch_overrides = []
        for ch in ("S", "Tt", "Xt", "V", "A", "PS"):
            win_ch = windows.per_channel.get((geom, ch))
            if win_ch:
                ch_overrides.append(f"`{ch}: {win_ch[0]} {win_ch[1]}`")
        ch_str = ", ".join(ch_overrides) if ch_overrides else "无"

        note = obs_summary.get(geom, ("-", "-", "-"))[2]
        lines.append(f"| **{geom}** | {s.ns}x{s.nt} | {s.temp_mev:.1f} | {def_win_str} | {ch_str} | {note} |")

    lines.extend([
        "",
        "---",
        "",
        "## 2. 格点规模全景概览卡片 (跨 4 种轻夸克质量对比)",
        "",
        "点击链接即可打开大图，用于评估某个格点规模下所有夸克质量是否在同一窗口平坦：",
        "",
    ])

    for geom in sorted(seen_geoms):
        card_p = base_out_dir / "ensemble_overview" / f"overview_{geom}.png"
        lines.append(f"- **{geom} 全景总览卡片**: [{card_p.name}](file://{card_p.resolve()})")

    lines.extend([
        "",
        "---",
        "",
        "## 3. 按目录与频道分类存储的图片索引导航",
        "",
        "### 3.1 按目录划分 (25 个目录)",
        "",
    ])

    for s in specs:
        dir_folder = base_out_dir / "by_directory" / s.dir_name
        all_ch_p = dir_folder / "meff_all_channels.png"
        csv_p = dir_folder / "meff_data.csv"

        lines.append(f"#### 📂 目录：`{s.dir_name}`")
        lines.append(f"- **参数**：格点 ${s.ns}\\times {s.nt}$，半宽 $half={s.half}$，$m_l={s.ml:.4f}$，$T={s.temp_mev:.2f}$ MeV")
        lines.append(f"- **6 信道全景图**：[meff_all_channels.png](file://{all_ch_p.resolve()})")
        lines.append(f"- **数据表格**：[meff_data.csv](file://{csv_p.resolve()})")
        lines.append("- **单信道图**：")
        ch_links = []
        for ch in CHANNELS:
            ch_p = dir_folder / f"meff_{ch}.png"
            ch_links.append(f"[{ch}](file://{ch_p.resolve()})")
        lines.append("  " + " | ".join(ch_links))
        lines.append("")

    lines.extend([
        "---",
        "",
        "### 3.2 按信道划分 (夸克质量对比图)",
        "",
    ])

    for ch in CHANNELS:
        lines.append(f"#### 🏷️ 信道：`{ch}` ({CHANNEL_TITLES.get(ch, ch)})")
        lines.append("在不同格点规模下对比 4 种夸克质量平台：")
        comp_links = []
        for geom in sorted(seen_geoms):
            comp_p = base_out_dir / "by_channel" / ch / f"{geom}_compare_ml.png"
            if comp_p.exists():
                comp_links.append(f"[{geom} 对比图](file://{comp_p.resolve()})")
        lines.append("- " + " | ".join(comp_links))
        lines.append("")

    lines.extend([
        "---",
        "",
        "## 4. 如何更新 `config/fit_windows.txt`",
        "",
        "1. 打开 [`config/fit_windows.txt`](file://" + str(FIT_WINDOW_FILE.resolve()) + ")；",
        "2. 根据上述图表中观测到的平台区间修改或新增行，例如为 `36x18` 指定：",
        "   ```text",
        "   36x18              14       19",
        "   36x18  S            6       19",
        "   ```",
        "3. 修改保存后，直接执行：",
        "   ```bash",
        "   .venv/bin/python scripts/update_ccor_data_and_plots.py",
        "   ```",
        "   系统将自动采用新窗口重拟合介子质量并更新全部图表！",
        "",
    ])

    readme_path.write_text("\n".join(lines), encoding="utf-8")
    return readme_path


def main():
    parser = argparse.ArgumentParser(description="绘制所有 beta=4.17 目录的 effective mass 图像")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "docs" / "effective_mass",
        help="输出基目录 (默认: docs/effective_mass)",
    )
    args = parser.parse_args()
    base_out_dir = args.output_dir.resolve()
    base_out_dir.mkdir(parents=True, exist_ok=True)

    print("================================================================================")
    print("      正在批量绘制所有 beta=4.17 目录的 Effective Mass 图像流水线")
    print(f"      输出根目录: {base_out_dir}")
    print("================================================================================")

    # 1. 加载拟合窗口与拟合总结
    windows = load_fit_windows()
    summary_dict = load_summary_fit_dict()
    specs = build_directory_list()
    print(f"[INIT] 收集到 {len(specs)} 个 beta=4.17 目标目录。")

    # 2. 遍历每个目录绘制单信道图与全景图
    print("\n[STEP 1/3] 正在生成各个目录下的 effective mass 图与数据表格...")
    geom_grouped: Dict[str, List[DirectorySpec]] = {}

    for idx, spec in enumerate(specs, 1):
        geom_grouped.setdefault(spec.case_key, []).append(spec)
        dir_folder = base_out_dir / "by_directory" / spec.dir_name
        dir_folder.mkdir(parents=True, exist_ok=True)

        meff_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        csv_rows = []

        for ch in CHANNELS:
            x, mean, err = load_channel_meff(spec, ch)
            if len(x) == 0:
                continue
            meff_data[ch] = (x, mean, err)

            window = get_current_window(windows, spec.case_key, ch)
            summary_row = summary_dict.get((spec.case_key, round(spec.ml, 4), ch))

            # 绘制单信道图
            ch_plot_path = dir_folder / f"meff_{ch}.png"
            plot_single_channel(spec, ch, x, mean, err, window, summary_row, ch_plot_path)

            # 收集 CSV 数据
            valid = (x >= 0) & (x <= spec.half)
            for xi, mi, ei in zip(x[valid], mean[valid], err[valid]):
                csv_rows.append({
                    "channel": ch,
                    "x": int(xi),
                    "meff_mean": float(mi) if np.isfinite(mi) else None,
                    "meff_err": float(ei) if np.isfinite(ei) else None,
                })

        # 保存目录内 6 信道全景图
        all_channels_path = dir_folder / "meff_all_channels.png"
        plot_all_channels_grid(spec, meff_data, windows, summary_dict, all_channels_path)

        # 保存 CSV
        if csv_rows:
            pl.DataFrame(csv_rows).write_csv(dir_folder / "meff_data.csv")

        print(f"  [{idx:02d}/{len(specs):02d}] 完成目录: {spec.dir_name}")

    # 3. 按信道绘制夸克质量对比图与格点全景卡片
    print("\n[STEP 2/3] 正在生成按信道夸克质量对比图与格点全景卡片...")
    for geom_key, geom_specs in geom_grouped.items():
        # 格点全景卡片
        overview_path = base_out_dir / "ensemble_overview" / f"overview_{geom_key}.png"
        plot_ensemble_overview_card(geom_specs, windows, overview_path)
        print(f"  [OVERVIEW] {overview_path.name}")

        # 各信道质量对比图
        for ch in CHANNELS:
            ch_comp_path = base_out_dir / "by_channel" / ch / f"{geom_key}_compare_ml.png"
            plot_by_channel_ml_comparison(geom_specs, ch, windows, ch_comp_path)

    # 4. 生成 README.md 综合索引与分析导航
    print("\n[STEP 3/3] 正在生成 docs/effective_mass/README.md 综合报告...")
    readme_path = generate_readme_report(specs, windows, summary_dict, base_out_dir)
    print(f"  [README] 报告已保存: {readme_path}")

    print("\n================================================================================")
    print("      全部 Effective Mass 图像与数据已按目录和信道完整生成！")
    print(f"      总览报告: {readme_path}")
    print("================================================================================")


if __name__ == "__main__":
    main()
