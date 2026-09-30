#!/usr/bin/env python3
"""
scripts/try_fit_window.py
--------------------------------------------------------------------------------
小工具: 快速试错介子 cosh 平台拟合窗口 (只读, 不改动任何已落盘结果)。

流程 (先便宜后贵):
  1. 直接读 meson_scan 已落盘的数据: corr_<ch>.csv (均值/误差) + jk_<ch>.csv (Jackknife 样本);
  2. 用均值关联函数把所有候选窗口 [x_start, x_end) 扫一遍 -> mass / chi2/dof (很便宜);
  3. 对 1D 网格 (终点固定在对称点+1) 的全部窗口做完整 Jackknife 拟合 -> mass ± err (贵但值得);
  4. 用 "相邻窗口质量差是否在误差内一致" 的判据给出真正的平台起点:
        t_s = (m_s - m_{s+1}) / sqrt(err_s^2 + err_{s+1}^2)
     从长窗口往短窗口找, 取满足所有 |t_k| <= tol 的最长窗口 -> 推荐窗口;
  5. 同时打印 chi2/dof 最优的 top-N 窗口、detect_plateau_window 自动检测结果、
     以及当前生产结果用的窗口, 便于对照。

注意: 强关联数据上 chi2/dof 往往 << 1, 单看 chi2/dof 会一直偏向最短窗口;
      实际挑窗口请以 (4) 的质量稳定性 + (5) 的自动平台检测为准。

用法示例:
  .venv/bin/python scripts/try_fit_window.py --case 48x18 --ml 0.0020
  .venv/bin/python scripts/try_fit_window.py --case 40x16 --ml 0.0020 --channels S Xt Vec
  .venv/bin/python scripts/try_fit_window.py --case 48x18 --ml all --channels AV Vec
  .venv/bin/python scripts/try_fit_window.py --case 48x18 --ml 0.0020 --grid 2d --top 8
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.meson_scan_setup import CHANNEL_ORDER, MIN_WINDOW_POINTS, SCAN_CASES  # noqa: E402
from src.parselqcdata import (  # noqa: E402
    compute_effective_mass_matrix_centered,
    detect_plateau_window,
    fit_jackknife_mass_centered,
    fit_mass_window_scan,
)

SCAN_ROOT: Path = PROJECT_ROOT / "output" / "meson_scan" / "b4.17"
SOURCE_TAGS = ("multisrc", "singlesrc")

#: 生产流程用的阈值 (与 reproduce_meson_b417_nt_scan.py 一致)
CHI2_DOF_MAX: float = 3.0
PLATEAU_CHI2_MAX: float = 1.5


def find_case(case_key: str):
    for case in SCAN_CASES:
        if case.key == case_key:
            return case
    raise SystemExit(f"[ERROR] 未知 ensemble: {case_key} (可选: {[c.key for c in SCAN_CASES]})")


def case_dir(case_key: str, ml: float, source: str) -> Path:
    case = find_case(case_key)
    return SCAN_ROOT / source / "cases" / case.output_name(float(ml))


def load_channel(case_path: Path, channel: str):
    """返回 (x_axis, mean, err, jk_matrix[ns, n_bins]); 缺文件返回 None"""
    corr_p = case_path / f"corr_{channel}.csv"
    jk_p = case_path / f"jk_{channel}.csv"
    if not corr_p.exists() or not jk_p.exists():
        return None
    corr = pl.read_csv(corr_p)
    x_axis = corr["x"].to_numpy().astype(np.float64)
    mean = corr["mean"].to_numpy().astype(np.float64)
    err = corr["err"].to_numpy().astype(np.float64)
    jk = pl.read_csv(jk_p, has_header=False).to_numpy().astype(np.float64)
    if jk.shape[0] != x_axis.size:
        jk = jk.T
    return x_axis, mean, err, jk


def production_window(case_key: str, ml: float, channel: str, source: str) -> str:
    """当前生产结果采用的窗口 (读汇总表)"""
    summary = SCAN_ROOT / source / "meson_mass_summary.csv"
    if not summary.exists():
        return "?"
    df = pl.read_csv(summary).filter(
        (pl.col("case_key") == case_key) & (pl.col("channel") == channel) & (pl.col("ml") == float(ml))
    )
    if df.height == 0:
        return "?"
    r = df.row(0, named=True)
    return f"[{r['x_start']},{r['x_end']}) {r['window_mode']}"


def windows_1d(half: int, min_points: int) -> List[Tuple[int, int]]:
    """终点固定在对称点+1, 扫起点 (生产流程用的网格)"""
    return [(s, half + 1) for s in range(1, max(half - min_points + 1, 1) + 1)]


def windows_2d(half: int, min_points: int) -> List[Tuple[int, int]]:
    """起点 + 终点都扫"""
    out: List[Tuple[int, int]] = []
    for s in range(1, half - min_points + 2):
        for e in range(s + min_points, half + 2):
            out.append((s, e))
    return out


def plateau_from_stability(jk_rows: List[dict], tol: float) -> Optional[int]:
    """
    平台起点判据: 相邻窗口 (s, s+1) 的质量差
        t_s = (m_s - m_{s+1}) / sqrt(err_s^2 + err_{s+1}^2)
    从长窗口往短窗口找, 返回满足所有 |t_k| <= tol 的最长窗口起点。
    """
    ns = [r for r in jk_rows if np.isfinite(r["mass"]) and np.isfinite(r["mass_err"])]
    if len(ns) < 2:
        return None
    starts = [r["x_start"] for r in ns]
    for idx in range(len(ns) - 1):
        ok = True
        for k in range(idx, len(ns) - 1):
            a, b = ns[k], ns[k + 1]
            denom = np.hypot(a["mass_err"], b["mass_err"])
            if denom <= 0:
                ok = False
                break
            if abs(a["mass"] - b["mass"]) / denom > tol:
                ok = False
                break
        if ok:
            return starts[idx]
    return None


def scan_channel(case_key: str, ml: float, channel: str, source: str, grid: str, top: int,
                 min_points: int, plateau_chi2: float, tol: float, full_jk: bool) -> Optional[dict]:
    case = find_case(case_key)
    half = int(case.ns // 2)
    cdir = case_dir(case_key, ml, source)
    loaded = load_channel(cdir, channel)
    if loaded is None:
        print(f"[WARN] {case_key} ml={ml:.4f} {channel}: 缺少 corr_{channel}.csv / jk_{channel}.csv -> 跳过")
        return None
    x_axis, mean, err, jk = loaded

    windows = windows_1d(half, min_points) if grid == "1d" else windows_2d(half, min_points)
    scan = fit_mass_window_scan(x_axis, mean, err, windows, half)

    ok = [
        r for r in scan
        if r.get("ok") and r["n_points"] >= min_points
        and np.isfinite(r["mass"]) and r["mass"] > 0.0 and np.isfinite(r["chi2_dof"])
    ]
    ranked = sorted(ok, key=lambda r: (round(float(r["chi2_dof"]), 2), -int(r["n_points"]), int(r["x_start"])))

    print()
    print("=" * 104)
    print(f" {case_key}  ml={ml:.4f}  {channel:<3}  (Ns={case.ns}, half={half}, n_bins={jk.shape[1]}, "
          f"source={source})")
    print(f" 当前生产窗口: {production_window(case_key, ml, channel, source)}")
    print("=" * 104)

    # ---- 1D: 全窗口 Jackknife + 稳定性判据 ----
    if grid == "1d":
        print(" [1D 扫描: 终点固定 %d; 质量 = 完整 Jackknife 拟合]" % (half + 1))
        print(f" {'start':>6} {'npts':>5} {'mass(JK)':>12} {'err(JK)':>10} {'chi2/dof':>9} "
              f"{'t_next(σ)':>10}   {'mass(mean)':>11}")
        jk_rows: List[dict] = []
        for r in scan:
            s, e = int(r["x_start"]), int(r["x_end"])
            full = fit_jackknife_mass_centered(x_axis, jk, err, s, e, half) if full_jk else None
            row = {
                "x_start": s, "x_end": e, "n_points": r["n_points"],
                "mass": (full or {}).get("mass", np.nan) if full else r["mass"],
                "mass_err": (full or {}).get("mass_err", np.nan) if full else r["mass_err"],
                "chi2_dof": (full or {}).get("chi2_dof", np.nan) if full else r["chi2_dof"],
                "mass_mean": r["mass"] if r.get("ok") else np.nan,
                "n_fitted": (full or {}).get("n_fitted", 0) if full else 0,
                "n_samples": (full or {}).get("n_samples", 0) if full else 0,
            }
            jk_rows.append(row)

        for i, row in enumerate(jk_rows):
            nxt = jk_rows[i + 1] if i + 1 < len(jk_rows) else None
            t_txt = ""
            if nxt is not None and np.isfinite(row["mass"]) and np.isfinite(nxt["mass"]):
                denom = np.hypot(row["mass_err"], nxt["mass_err"])
                t_txt = f"{(row['mass'] - nxt['mass']) / denom:>10.2f}" if denom > 0 else f"{'inf':>10}"
            m_jk = f"{row['mass']:>12.6f}" if np.isfinite(row["mass"]) else f"{'--':>12}"
            e_jk = f"{row['mass_err']:>10.6f}" if np.isfinite(row["mass_err"]) else f"{'--':>10}"
            c2 = f"{row['chi2_dof']:>9.3f}" if np.isfinite(row["chi2_dof"]) else f"{'--':>9}"
            mm = f"{row['mass_mean']:>11.6f}" if np.isfinite(row["mass_mean"]) else f"{'--':>11}"
            print(f" {row['x_start']:>6} {row['n_points']:>5} {m_jk} {e_jk} {c2} {t_txt}   {mm}")

        p_start = plateau_from_stability(jk_rows, tol)
        if p_start is not None:
            rec = next(r for r in jk_rows if r["x_start"] == p_start)
            print(f"\n [稳定性判据 |t|<={tol:g}σ] 最长平台起点 = {p_start}  ->  "
                  f"推荐窗口 [{p_start},{half + 1}) : m = {rec['mass']:.6f} ± {rec['mass_err']:.6f}")
        else:
            print(f"\n [稳定性判据 |t|<={tol:g}σ] 未找到平台 (质量随窗口一直在漂)")

        print(f"\n [chi2/dof 最优 top-{top}]  (注意: 强关联数据 chi2/dof 偏小, 会偏向最短窗口)")
        print(f" {'window':>12} {'npts':>5} {'mass(JK)':>12} {'err(JK)':>10} {'chi2/dof':>9} {'fitted':>9}")
        for r in ranked[:top]:
            full = fit_jackknife_mass_centered(x_axis, jk, err, int(r["x_start"]), int(r["x_end"]), half)
            if not np.isfinite(full["mass"]):
                continue
            print(f" [{full['x_start']:>3},{full['x_end']:>3}) {r['n_points']:>5} {full['mass']:>12.6f} "
                  f"{full['mass_err']:>10.6f} {full['chi2_dof']:>9.3f} "
                  f"{full['n_fitted']:>4}/{full['n_samples']:<4}")

    # ---- 2D: 只列 top ----
    else:
        print(f" [2D 扫描: {len(windows)} 个窗口, 只列 chi2/dof 最优 top-{max(top, 20)}]")
        print(f" {'window':>12} {'npts':>5} {'mass(mean)':>12} {'chi2/dof':>9}")
        for r in ranked[:max(top, 20)]:
            print(f" [{r['x_start']:>3},{r['x_end']:>3}) {r['n_points']:>5} {r['mass']:>12.6f} {r['chi2_dof']:>9.3f}")

    # ---- 自动平台检测 (当前 plateau_auto 的判据) ----
    meff_matrix, _, _ = compute_effective_mass_matrix_centered(jk, float(half))
    auto = detect_plateau_window(meff_matrix, end=half, chi2_dof_max=plateau_chi2, min_points=4)
    auto_txt = (f"[{auto['start']},{auto['end']}) flat chi2/dof={auto['chi2_dof']:.3f}"
                if auto.get("ok") else f"未检出 (flat chi2/dof={auto['chi2_dof']:.3f} > {plateau_chi2})")
    print(f" 自动平台检测 (平坦性 chi2/dof<={plateau_chi2}, min_points=4): {auto_txt}")
    return None


def scan_ana_half_window(case_key: str, ml: float, channels: Sequence[str], source: str,
                         n_min: int, pairs: bool) -> None:
    """
    ccor_flow / ana simulate.py 口径: 先按对称点 fold, 再在 [c-N, c+N+1) 内逐 Jackknife
    样本做 cosh 拟合 (chi2/dof 阈值 100)。这里扫描半宽 N —— 论文 ΔM 图用的就是这个窗口,
    默认 N=7 (Ns=32 -> [9,24), Ns=40 -> [13,28), Ns=48 -> [17,32))。
    """
    from src.parselqcdata.ccor_flow import (  # 局部导入: 只在 ana 模式下需要
        ANA_MIN_MASS,
        ANA_PAIRS,
        channel_masses_from_jk,
        jackknife_mean_err,
        load_fit_windows,
    )

    case = find_case(case_key)
    half = int(case.ns // 2)
    # 生产口径 = config/fit_windows.txt 里该 ensemble 的窗口 (唯一来源)
    prod_wins = {c: load_fit_windows().get(case_key, c) for c in channels}
    cdir = case_dir(case_key, ml, source)

    #: meson_scan 的信道名 -> ana (ccor_flow) 的信道名
    ana_name = {"Vec": "V", "AV": "A", "Tt": "Tt", "Xt": "Xt", "S": "S", "PS": "PS"}
    per_channel: Dict[str, Dict[int, np.ndarray]] = {}
    for ch in channels:
        loaded = load_channel(cdir, ch)
        if loaded is None:
            print(f"[WARN] {case_key} ml={ml:.4f} {ch}: 缺少数据 -> 跳过")
            continue
        jk = loaded[3]
        per_channel[ana_name.get(ch, ch)] = {
            n: channel_masses_from_jk(jk, center=half, half_window=n, min_mass=ANA_MIN_MASS)[0]
            for n in range(n_min, half - 1)
        }

    print()
    print("=" * 104)
    print(f" {case_key}  ml={ml:.4f}  [ccor_flow/ana 口径: 窗口 = [c-N, c+N+1), c={half}]  source={source}")
    print(" 生产窗口 (config/fit_windows.txt): " + ", ".join(
        f"{c}={list(w) if w else '未定义'}" for c, w in prod_wins.items()))
    print("=" * 104)

    for ch, per_n in per_channel.items():
        print(f"\n {ch}:")
        print(f" {'N':>3} {'window':>13} {'mass':>12} {'err':>10} {'n_valid':>8}")
        for n, masses in per_n.items():
            mean, err = jackknife_mean_err(masses)
            tag = "   <- 生产默认" if n == prod_n else ""
            print(f" {n:>3} [{half - n:>3},{half + n + 1:>3}) {mean:>12.6f} {err:>10.6f} "
                  f"{int(np.isfinite(masses).sum()):>8}{tag}")

    if pairs:
        print("\n [对称性破缺质量差 ΔM (逐样本作差 + Jackknife)]")
        for ch1, ch2, title in ANA_PAIRS:
            if ch1 not in per_channel or ch2 not in per_channel:
                continue
            print(f"\n  {title}  ({ch1} - {ch2}):")
            print(f" {'N':>3} {'window':>13} {'dM(lattice)':>13} {'err':>10} {'dM(MeV)':>12}")
            for n, a in per_channel[ch1].items():
                b = per_channel[ch2].get(n)
                if b is None:
                    continue
                k = min(a.size, b.size)
                mean, err = jackknife_mean_err(a[:k] - b[:k])
                tag = "   <- 生产默认" if n == prod_n else ""
                print(f" {n:>3} [{half - n:>3},{half + n + 1:>3}) {mean:>13.6f} {err:>10.6f} "
                      f"{mean * 2453:>12.2f}{tag}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="快速试错 cosh 平台拟合窗口 (只读小工具)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--case", required=True, help="ensemble, 如 48x18 / 40x16 / 32x12")
    p.add_argument("--ml", required=True, help="裸轻夸克质量, 如 0.0020; 或用 all 跑全部 4 个质量")
    p.add_argument("--channels", nargs="*", default=list(CHANNEL_ORDER), help="信道列表")
    p.add_argument("--source", choices=SOURCE_TAGS, default="multisrc", help="多源 / 单源数据")
    p.add_argument("--grid", choices=("1d", "2d"), default="1d",
                   help="1d=终点固定在对称点+1, 只扫起点 (生产口径); 2d=起点终点都扫")
    p.add_argument("--top", type=int, default=5, help="额外列出 chi2/dof 最优的前 N 个窗口")
    p.add_argument("--min-points", type=int, default=MIN_WINDOW_POINTS, help="窗口最小点数")
    p.add_argument("--plateau-chi2", type=float, default=PLATEAU_CHI2_MAX, help="自动平台检测的平坦性阈值")
    p.add_argument("--tol", type=float, default=1.0, help="稳定性判据的 |t| 上限 (单位 σ)")
    p.add_argument("--mode", choices=("scan", "ana"), default="scan",
                   help="scan=meson_scan 拟合窗口扫描; ana=ccor_flow/ana 半宽扫描 (论文 ΔM 口径)")
    p.add_argument("--n-min", type=int, default=3, help="ana 模式下半宽 N 的最小值")
    p.add_argument("--no-pairs", dest="pairs", action="store_false", default=True,
                   help="ana 模式下不打印 4 组 ΔM 随 N 的变化")
    p.add_argument("--no-full-jk", dest="full_jk", action="store_false", default=True,
                   help="1D 网格下跳过全窗口 Jackknife (只用均值拟合, 更快)")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    mls = [0.0020, 0.0035, 0.0070, 0.0120] if args.ml == "all" else [float(args.ml)]
    for ch in args.channels:
        if ch not in CHANNEL_ORDER:
            raise SystemExit(f"[ERROR] 未知信道 {ch} (可选: {list(CHANNEL_ORDER)})")
    if args.mode == "ana":
        for ml in mls:
            scan_ana_half_window(args.case, ml, args.channels, args.source, args.n_min, args.pairs)
        return 0
    for ml in mls:
        for ch in args.channels:
            scan_channel(args.case, ml, ch, args.source, args.grid, args.top,
                         args.min_points, args.plateau_chi2, args.tol, args.full_jk)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
