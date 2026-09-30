#!/usr/bin/env python3
"""
scripts/measure_meson_b417.py
--------------------------------------------------------------------------------
测量所有以数字开头且包含 beta4.17 的目录中的介子质量，并绘制介子质量随温度演化的图。
物理量纲转换乘法常数: a^-1 = 2453 MeV.
温度: T = 2453 / Nt (MeV).
介子质量: M = (a * M) * 2453 (MeV).

为不同的轻夸克质量 (ml = 0.0020, 0.0035, 0.0070, 0.0120) 分别绘制独立的温度演化图，
并保存到 docs/ 目录下。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

# 确保项目根目录在 sys.path 中
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.meson_scan_setup import (
    BETA_KEY,
    CHANNEL_ORDER,
    MIN_WINDOW_POINTS,
    SCAN_CASES,
    SOURCE_TAGS,
    binsize_for,
    case_output_dir,
    reference_window,
    scan_windows,
)
from docs.physics_setup import CHANNEL_CONFIGS, DOCS_DIR, FIGURES_DIR, TRANSITION_REGION
from src.parselqcdata import (
    compute_effective_mass_matrix_centered,
    detect_plateau_window,
    fit_jackknife_mass_centered,
    fit_mass_window_scan,
    select_best_window,
)
from tools.meson_orchestrator import MesonOrchestrator

# 物理量纲换算常数: a^-1 = 2453 MeV (旧稿的 153*16 = 2448 MeV 不对)
SCALE_UNIT: float = 2453.0

# 介子信道配色与标记样式
CHANNEL_STYLES: Dict[str, dict] = {
    "PS": {
        "color": "#d62728",
        "marker": "o",
        "label": r"$PS$ ($\bar{\psi}\gamma_5\psi$)",
        "linestyle": "-",
        "zorder": 6,
    },
    "Vec": {
        "color": "#1f77b4",
        "marker": "s",
        "label": r"$Vec$ ($\bar{\psi}\gamma_\mu\psi$)",
        "linestyle": "-",
        "zorder": 5,
    },
    "AV": {
        "color": "#2ca02c",
        "marker": "^",
        "label": r"$AV$ ($\bar{\psi}\gamma_5\gamma_\mu\psi$)",
        "linestyle": "--",
        "zorder": 4,
    },
    "S": {
        "color": "#ff7f0e",
        "marker": "D",
        "label": r"$S$ ($\bar{\psi}\psi$)",
        "linestyle": "--",
        "zorder": 3,
    },
    "Tt": {
        "color": "#9467bd",
        "marker": "v",
        "label": r"$T_t$ ($\bar{\psi}\sigma_{0\mu}\psi$)",
        "linestyle": ":",
        "zorder": 2,
    },
    "Xt": {
        "color": "#8c564b",
        "marker": "p",
        "label": r"$X_t$ ($\bar{\psi}\gamma_5\sigma_{0\mu}\psi$)",
        "linestyle": ":",
        "zorder": 2,
    },
}

# 正则匹配数字开头且含 b4.17 / beta4.17 的格点目录
B417_DIR_PATTERN = re.compile(
    r"^(?P<ns>\d+)x(?P<nt>\d+)_b(?:eta)?(?P<beta>4\.17)_ms(?P<ms>[\d.]+)m(?P<ml>[\d.]+)$"
)


def discover_b417_directories(data_dir: Path) -> List[Tuple[str, Path, dict]]:
    """
    扫描数据目录，发现所有数字开头且包含 beta4.17 的目录，并解析物理参数。
    """
    found = []
    if not data_dir.exists():
        print(f"[ERROR] 数据目录不存在: {data_dir}")
        return found

    for child in sorted(data_dir.iterdir()):
        if not child.is_dir():
            continue
        m = B417_DIR_PATTERN.match(child.name)
        if m:
            info = {
                "ns": int(m.group("ns")),
                "nt": int(m.group("nt")),
                "beta": float(m.group("beta")),
                "ms": float(m.group("ms")),
                "ml": float(m.group("ml")),
                "case_key": f"{m.group('ns')}x{m.group('nt')}",
                "label": f"{m.group('ns')}^3x{m.group('nt')}",
            }
            found.append((child.name, child, info))
        elif re.match(r"^\d+.*4\.17.*", child.name):
            print(f"[INFO] 发现数字开头 4.17 相关目录 (非标准命名或对照格点): {child.name}")

    return found


def measure_or_load_meson_mass(
    dir_name: str,
    case_path: Path,
    info: dict,
    orch: Optional[MesonOrchestrator] = None,
    force: bool = False,
    window_policy: str = "auto",
) -> List[dict]:
    """
    对单个格点目录执行介子质量测量（或从已有缓存中加载）。
    """
    ns = info["ns"]
    nt = info["nt"]
    ml = info["ml"]
    half = ns / 2.0
    output_case_name = f"{ns}x{nt}_ms{info['ms']:.3f}_m{ml:.4f}"

    # 优先查找是否有既有落盘
    cached_paths = [
        PROJECT_ROOT / "output" / "meson_scan" / f"b4.{BETA_KEY}" / "multisrc" / "cases" / f"{info['case_key']}_ms0.040_m{ml:.4f}" / "meta.json",
        *list((PROJECT_ROOT / "output" / "meson_scan" / f"b4.{BETA_KEY}" / "multisrc" / "cases").glob(f"*{output_case_name}*/meta.json")),
    ]

    for c_path in cached_paths:
        if c_path.exists() and not force:
            try:
                meta = json.loads(c_path.read_text(encoding="utf-8"))
                rows = meta.get("rows", [])
                if rows and len(rows) == len(CHANNEL_ORDER):
                    return rows
            except Exception:
                pass

    # 若未找到缓存或需要强制重算，调用 MesonOrchestrator 抽取与拟合
    input_dir = case_path / "Output"
    if not input_dir.exists():
        print(f"[WARN] 找不到关联函数输入路径: {input_dir}")
        return []

    if orch is None:
        orch = MesonOrchestrator()

    x_axis = np.arange(ns, dtype=np.float64)
    scan_win_list = scan_windows(
        next(c for c in SCAN_CASES if c.ns == ns and c.nt == nt)
    )

    rows = []
    print(f"[MEASURE] 开始计算格点 {dir_name} (Ns={ns}, Nt={nt}, ml={ml:.4f})...")

    for ch in CHANNEL_ORDER:
        bs = binsize_for(ch, is_single_source=False)
        means, errors, jk, n_bins, n_raw = orch.process_channel(
            input_dir=str(input_dir),
            channel_configs=CHANNEL_CONFIGS[ch],
            binsize=bs,
            num_lines=ns,
            thread_count=4,
            is_single_source=False,
        )

        if jk.size == 0 or jk.shape[0] != ns:
            print(f"[WARN] {dir_name} 通道 {ch} 抽取失败")
            continue

        meff_matrix, meff_mean, meff_err = compute_effective_mass_matrix_centered(jk, half)
        auto = detect_plateau_window(meff_matrix, end=int(half), chi2_dof_max=1.5, min_points=4)
        scan = fit_mass_window_scan(x_axis, means, errors, scan_win_list, half)
        best = select_best_window(scan, min_points=MIN_WINDOW_POINTS, chi2_dof_max=3.0)
        if best is None:
            best = select_best_window(scan, min_points=MIN_WINDOW_POINTS, chi2_dof_max=float("inf"))

        if best is not None:
            scan_fit = fit_jackknife_mass_centered(
                x_axis, jk, errors, int(best["x_start"]), int(best["x_end"]), half
            )
        else:
            scan_fit = {"mass": np.nan, "mass_err": np.nan}

        if auto["ok"]:
            auto_fit = fit_jackknife_mass_centered(
                x_axis, jk, errors, int(auto["start"]), int(auto["end"]), half
            )
        else:
            auto_fit = {"mass": np.nan, "mass_err": np.nan}

        if window_policy == "auto" and auto["ok"]:
            primary = auto_fit
        else:
            primary = scan_fit

        row = {
            "ns": ns,
            "nt": nt,
            "case_key": info["case_key"],
            "label": info["label"],
            "beta": info["beta"],
            "ms": info["ms"],
            "ml": ml,
            "channel": ch,
            "mass": float(primary.get("mass", np.nan)),
            "mass_err": float(primary.get("mass_err", np.nan)),
        }
        rows.append(row)

    return rows


def build_measurement_dataframe(readin_dir: Path, force: bool = False) -> pl.DataFrame:
    """
    收集所有目录的测量结果，并施加物理量纲转换。
    单位比例: a^-1 = 2453 MeV
    """
    dirs = discover_b417_directories(readin_dir)
    print(f"[INFO] 共发现 {len(dirs)} 个符合条件的 beta4.17 构型目录。")

    all_rows = []
    summary_csv = PROJECT_ROOT / "output" / "meson_scan" / f"b4.{BETA_KEY}" / "multisrc" / "meson_mass_summary.csv"

    # 如果有整体汇总表且不强制重新计算，直接加载汇总表
    if summary_csv.exists() and not force:
        print(f"[INFO] 快速加载现有全量拟合汇总表: {summary_csv}")
        df = pl.read_csv(summary_csv)
    else:
        orch = MesonOrchestrator()
        for name, path, info in dirs:
            rows = measure_or_load_meson_mass(name, path, info, orch=orch, force=force)
            all_rows.extend(rows)
        df = pl.DataFrame(all_rows)

    # 实施用户指定的物理量纲转换:
    # 物理量纲需要乘的单位是 a^-1 = 2453 MeV
    df = df.with_columns([
        (pl.col("mass") * SCALE_UNIT).alias("mass_mev_scaled"),
        (pl.col("mass_err") * SCALE_UNIT).alias("mass_mev_err_scaled"),
        (SCALE_UNIT / pl.col("nt")).alias("T_mev_scaled"),
    ])

    return df


def plot_meson_mass_for_quark_mass(
    df_ml: pl.DataFrame,
    ml_val: float,
    output_path: Path,
) -> None:
    """
    针对给定的夸克质量 ml，绘制 6 个介子信道随温度变化的出版级物理曲线。
    """
    fig, ax = plt.subplots(figsize=(9, 6.2), dpi=300)

    # 标记相变过渡区 (Transition region ~ 155 - 158 MeV)
    t_min, t_max = TRANSITION_REGION
    ax.axvspan(
        t_min, t_max,
        color="#ebebeb", alpha=0.7,
        label=f"Transition Region ({t_min:.1f}–{t_max:.1f} MeV)",
        zorder=1,
    )

    # 信道微小水平位移，防止同温度下误差棒互相重叠
    offsets = {
        "PS": -0.6,
        "Vec": -0.2,
        "AV": 0.2,
        "Tt": 0.6,
        "Xt": 1.0,
        "S": -1.0,
    }

    primary_cases = ["32x12", "32x14", "40x16", "48x18"]

    for ch in CHANNEL_ORDER:
        st = CHANNEL_STYLES.get(ch, {"color": "black", "marker": "o", "label": ch, "linestyle": "-"})
        sub = df_ml.filter(pl.col("channel") == ch)
        if sub.is_empty():
            continue

        dx = offsets.get(ch, 0.0)
        main_sub = sub.filter(pl.col("case_key").is_in(primary_cases)).sort("T_mev_scaled")
        other_sub = sub.filter(~pl.col("case_key").is_in(primary_cases)).sort("T_mev_scaled")

        t_main = main_sub["T_mev_scaled"].to_numpy() + dx
        m_main = main_sub["mass_mev_scaled"].to_numpy()
        err_main = main_sub["mass_mev_err_scaled"].to_numpy()

        # 绘制大体积主线与实心数据点
        ax.errorbar(
            t_main, m_main, yerr=err_main,
            fmt=st["marker"],
            color=st["color"],
            ecolor=st["color"],
            elinewidth=1.6,
            capsize=3.5,
            capthick=1.4,
            markersize=7.5,
            linestyle=st.get("linestyle", "-"),
            linewidth=1.4,
            label=st["label"],
            zorder=st.get("zorder", 3),
        )

        # 绘制较小格点体积对比点 (空心标记，小幅度微移)
        if not other_sub.is_empty():
            t_oth = other_sub["T_mev_scaled"].to_numpy() + dx + 0.35
            m_oth = other_sub["mass_mev_scaled"].to_numpy()
            err_oth = other_sub["mass_mev_err_scaled"].to_numpy()
            ax.errorbar(
                t_oth, m_oth, yerr=err_oth,
                fmt=st["marker"],
                markerfacecolor="none",
                markeredgecolor=st["color"],
                markeredgewidth=1.8,
                color=st["color"],
                ecolor=st["color"],
                elinewidth=1.2,
                capsize=3.0,
                capthick=1.2,
                markersize=7.5,
                linestyle="None",
                zorder=st.get("zorder", 3) - 1,
            )

    # 添加体积符号图例说明
    ax.scatter([], [], marker="o", color="gray", label="Large Vol ($40^3, 48^3, 32^3$)", s=45)
    ax.scatter([], [], marker="o", facecolors="none", edgecolors="gray", linewidth=1.5,
               label="Small Vol ($32^3x16, 36^3x18$)", s=45)

    ax.set_xlabel(r"Temperature $T$ [MeV]  ($T = 153 \times 16 / N_t$)", fontsize=13, fontweight="bold")
    ax.set_ylabel(r"Meson Mass $M$ [MeV]  ($M = aM \times 153 \times 16$)", fontsize=13, fontweight="bold")
    ax.set_title(
        rf"Thermal Evolution of Meson Masses ($\beta=4.17$, $m_l={ml_val:.4f}$, $m_s=0.040$)",
        fontsize=14,
        fontweight="bold",
        pad=12,
    )

    ax.tick_params(axis="both", which="major", labelsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(133, 207)
    ax.set_ylim(-50, 2600)
    ax.legend(loc="upper left", frameon=True, fontsize=9.5, framealpha=0.92, ncol=2)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"[SUCCESS] 介子质量温度演化图已生成 (ml={ml_val:.4f}): {output_path}")


def plot_combined_quark_mass_overview(
    df: pl.DataFrame,
    output_path: Path,
) -> None:
    """
    绘制所有 4 种轻夸克质量在 2x2 网格中的全局对比图。
    """
    mls = sorted(df["ml"].unique().to_list())
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=300, sharex=True, sharey=True)
    t_min, t_max = TRANSITION_REGION

    offsets = {
        "PS": -0.6,
        "Vec": -0.2,
        "AV": 0.2,
        "Tt": 0.6,
        "Xt": 1.0,
        "S": -1.0,
    }
    primary_cases = ["32x12", "32x14", "40x16", "48x18"]

    for idx, ml_val in enumerate(mls):
        ax = axes[idx // 2, idx % 2]
        sub_ml = df.filter(pl.col("ml") == ml_val)

        # 标记相变温区
        ax.axvspan(t_min, t_max, color="#e8e8e8", alpha=0.7, zorder=1)

        for ch in CHANNEL_ORDER:
            st = CHANNEL_STYLES.get(ch, {"color": "black", "marker": "o", "label": ch})
            ch_sub = sub_ml.filter(pl.col("channel") == ch)
            if ch_sub.is_empty():
                continue

            dx = offsets.get(ch, 0.0)
            main_sub = ch_sub.filter(pl.col("case_key").is_in(primary_cases)).sort("T_mev_scaled")
            other_sub = ch_sub.filter(~pl.col("case_key").is_in(primary_cases)).sort("T_mev_scaled")

            # 主线 (大体积)
            ax.errorbar(
                main_sub["T_mev_scaled"].to_numpy() + dx,
                main_sub["mass_mev_scaled"].to_numpy(),
                yerr=main_sub["mass_mev_err_scaled"].to_numpy(),
                fmt=st["marker"],
                color=st["color"],
                ecolor=st["color"],
                elinewidth=1.4,
                capsize=3.0,
                markersize=6.5,
                linestyle=st.get("linestyle", "-"),
                linewidth=1.2,
                label=st["label"] if idx == 0 else None,
                zorder=st.get("zorder", 3),
            )

            # 对比点 (小体积)
            if not other_sub.is_empty():
                ax.errorbar(
                    other_sub["T_mev_scaled"].to_numpy() + dx + 0.35,
                    other_sub["mass_mev_scaled"].to_numpy(),
                    yerr=other_sub["mass_mev_err_scaled"].to_numpy(),
                    fmt=st["marker"],
                    markerfacecolor="none",
                    markeredgecolor=st["color"],
                    markeredgewidth=1.6,
                    color=st["color"],
                    ecolor=st["color"],
                    elinewidth=1.0,
                    capsize=2.5,
                    markersize=6.5,
                    linestyle="None",
                    zorder=st.get("zorder", 3) - 1,
                )

        ax.set_title(rf"$m_l = {ml_val:.4f}$ ($m_s = 0.040$)", fontsize=13, fontweight="bold")
        ax.grid(True, linestyle="--", alpha=0.45)
        ax.tick_params(axis="both", which="major", labelsize=11)
        ax.set_ylim(-50, 2600)
        ax.set_xlim(133, 207)
        if idx % 2 == 0:
            ax.set_ylabel(r"Meson Mass $M$ [MeV]", fontsize=12, fontweight="bold")
        if idx // 2 == 1:
            ax.set_xlabel(r"Temperature $T$ [MeV]", fontsize=12, fontweight="bold")

    # 全局图例
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="upper center",
        ncol=len(handles),
        bbox_to_anchor=(0.5, 0.98),
        frameon=True,
        fontsize=10.5,
    )
    fig.suptitle(
        r"Thermal Evolution of Meson Masses across Light Quark Masses ($\beta=4.17$)",
        fontsize=16,
        fontweight="bold",
        y=1.01,
    )

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[SUCCESS] 4 种轻夸克质量全局对比图已保存: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="测量所有数字开头的 beta4.17 目录介子质量并按夸克质量分别绘图"
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "readin",
        help="格点数据根目录 (默认: data/readin)",
    )
    parser.add_argument(
        "--docs-dir",
        type=Path,
        default=DOCS_DIR,
        help="图像输出目录 (默认: docs)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="强制重新进行 C++ 抽取与拟合，不复用已有缓存",
    )
    args = parser.parse_args()

    print("================================================================================")
    print("      LQCD beta=4.17 有限温度介子质量测量与热演化科学绘图流水线")
    print(f"      物理标度常量: 153 * 16 = {SCALE_UNIT:.1f} MeV")
    print("================================================================================")

    # 1. 测量与数据汇总
    df = build_measurement_dataframe(args.data_dir, force=args.force)

    # 2. 针对不同的轻夸克质量分别绘图
    mls = sorted(df["ml"].unique().to_list())
    print(f"\n[INFO] 发现不同轻夸克质量列表: {mls}")

    for ml_val in mls:
        sub_df = df.filter(pl.col("ml") == ml_val)
        out_png = args.docs_dir / f"meson_mass_T_ml{ml_val:.4f}.png"
        plot_meson_mass_for_quark_mass(sub_df, ml_val, out_png)

        # 在 docs/figures/ 也同步保留一份副本
        fig_copy = FIGURES_DIR / f"meson_mass_T_ml{ml_val:.4f}.png"
        plot_meson_mass_for_quark_mass(sub_df, ml_val, fig_copy)

    # 3. 额外生成一张四网格对比大图
    overview_png = args.docs_dir / "meson_mass_T_all_ml.png"
    plot_combined_quark_mass_overview(df, overview_png)
    plot_combined_quark_mass_overview(df, FIGURES_DIR / "meson_mass_T_all_ml.png")

    print("\n================================================================================")
    print("      全部介子质量测量与绘图已圆满完成！图像已成功落盘至 docs/ 目录。")
    print("================================================================================")


if __name__ == "__main__":
    main()
