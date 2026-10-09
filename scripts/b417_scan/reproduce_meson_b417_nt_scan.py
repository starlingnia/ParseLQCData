#!/usr/bin/env python3
"""
scripts/reproduce_meson_b417_nt_scan.py
--------------------------------------------------------------------------------
beta = 4.17 有限温度扫描任务: 关联函数 + 介子质量 (Ns^3 x Nt, Nt = 12/14/16/18)

任务表 (顺序即落盘顺序, 见 docs/meson_scan_setup.py):
    00  32^3 x 12   ms=0.040   ml = 0.0020 / 0.0035 / 0.0070 / 0.0120
    01  32^3 x 14   ms=0.040   ml = 0.0020 / 0.0035 / 0.0070 / 0.0120
    02  32^3 x 16   ms=0.040   ml = 0.0020 / 0.0035 / 0.0070 / 0.0120
    03  36^3 x 18   ms=0.040   ml = 0.0020 / 0.0035 / 0.0070 / 0.0120
    04  40^3 x 16   ms=0.040   ml = 0.0020 / 0.0035 / 0.0070 / 0.0120
    05  48^3 x 18   ms=0.040   ml = 0.0020 / 0.0035 / 0.0070 / 0.0120

每个 ensemble 的处理链:
    1. MesonOrchestrator (C++ 核心库) 抽取 6 个狄拉克双线性信道的空间关联函数
       (num_lines = Ns, Jackknife 折叠矩阵);
    2. 有效质量 meff: C(x)/C(x+1) = cosh(m(x-Ns/2)) / cosh(m(x+1-Ns/2));
    3. 平台窗口扫描 (chi2/dof 选优) + 逐 Jackknife 样本 cosh 拟合 -> 介子质量;
    4. 同时给出 "镜像窗口" 结果: 把 4.17 @ 48^3x16 既有分析的窗口平移到本 ensemble;
    5. 全部结果按任务表顺序写入 output/meson_scan/b4.17/<source>/...

用法:
    .venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --dry-run
    .venv/bin/python scripts/reproduce_meson_b417_nt_scan.py                    # 多源 + 单源, 全部 24 个 ensemble
    .venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --source multi --workers 4
    .venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --cases 48x18 --ml 0.0020 --force
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.meson_scan_setup import (  # noqa: E402
    BETA_KEY,
    BETA_VALUE,
    CHANNEL_ORDER,
    MIN_WINDOW_POINTS,
    MS,
    SCAN_CASES,
    SOURCE_ORDER,
    SOURCE_TAGS,
    ScanCase,
    binsize_for,
    case_output_dir,
    discover_b417_dirs,
    reference_window,
    resolve_tasks,
    scan_root,
    scan_windows,
)
from docs.physics_setup import CHANNEL_CONFIGS  # noqa: E402
from src.parselqcdata import (  # noqa: E402
    LEAST_SQUARES_FITTER,
    compute_effective_mass_matrix_centered,
    detect_plateau_window,
    fit_jackknife_mass_centered,
    fit_mass_window_scan,
    select_best_window,
)
from tools.meson_orchestrator import MesonOrchestrator  # noqa: E402

warnings.filterwarnings("ignore")

#: 平台窗口选优允许的最大 chi2/dof (超限则放宽为 "只要求有限正值")
CHI2_DOF_MAX: float = 3.0

#: 有效质量平坦性检测阈值 (相邻差分与 0 一致的 chi2/dof 上限), 可用 --plateau-chi2 覆盖
PLATEAU_CHI2_MAX: float = 1.5


# ------------------------------------------------------------------------------
# 工具函数
# ------------------------------------------------------------------------------
def find_case(case_key: str) -> ScanCase:
    for case in SCAN_CASES:
        if case.key == case_key or case.label == case_key:
            return case
    raise KeyError(f"未知 ensemble: {case_key} (可选: {[c.key for c in SCAN_CASES]})")


def symmetrize_matrix(matrix: np.ndarray) -> np.ndarray:
    """
    空间关联函数对称化 (关于对称点 Ns/2 折叠, 对任意偶数 Ns 成立):
        C(x) -> [C(x) + C(Ns-x)] / 2,  x = 0 ... Ns-1
    对称点由矩阵行数 Ns 决定 (32->16, 36->18, 40->20, 48->24), 不写死任何尺寸。
    """
    ns = matrix.shape[0]
    idx = (ns - np.arange(ns)) % ns
    return 0.5 * (matrix + matrix[idx, :])


def build_tasks(
    out_root: Path,
    readin_dir: Optional[Path],
    sources: Sequence[bool],
    mls: Optional[Sequence[float]],
    case_keys: Optional[Sequence[str]],
    threads: int,
    window_policy: str = "auto",
) -> List[dict]:
    """生成有序任务列表: (ensemble, source) 笛卡尔积, 顺序严格跟随任务表"""
    ensembles = resolve_tasks(readin_dir=readin_dir, mls=mls, case_keys=case_keys)
    tasks: List[dict] = []
    for case in SCAN_CASES:
        if not any(e["case_key"] == case.key for e in ensembles):
            continue
        windows = [[s, e] for s, e in scan_windows(case)]
        for is_single_source in sources:
            ref_windows = {ch: list(reference_window(case, ch, is_single_source)) for ch in CHANNEL_ORDER}
            binsizes = {ch: binsize_for(ch, is_single_source) for ch in CHANNEL_ORDER}
            for ens in ensembles:
                if ens["case_key"] != case.key:
                    continue
                task = dict(ens)
                task.update(
                    {
                        "is_single_source": is_single_source,
                        "source_tag": SOURCE_TAGS[is_single_source],
                        "windows": windows,
                        "ref_windows": ref_windows,
                        "binsizes": binsizes,
                        "threads": threads,
                        "out_dir": str(case_output_dir(case, ens["ml"], is_single_source, out_root)),
                        "symmetrize": False,
                        "window_policy": window_policy,
                    }
                )
                tasks.append(task)
    return tasks


# ------------------------------------------------------------------------------
# 单个 ensemble 的处理
# ------------------------------------------------------------------------------
def process_ensemble(task: dict) -> dict:
    """抽取关联函数 -> meff -> 平台拟合 (在子进程中执行)"""
    started = time.time()
    out_dir = Path(task["out_dir"])
    base = {
        "order": task["order"],
        "ml_index": task["ml_index"],
        "case_key": task["case_key"],
        "label": task["label"],
        "ns": task["ns"],
        "nt": task["nt"],
        "ms": task["ms"],
        "beta": task["beta"],
        "ml": task["ml"],
        "source": task["source_tag"],
        "is_single_source": task["is_single_source"],
        "out_dir": str(out_dir),
    }

    if task["missing"]:
        return {**base, "status": "missing", "message": f"输入目录不存在: {task['input_dir']}", "rows": []}

    meta_path = out_dir / "meta.json"
    if meta_path.exists() and not task.get("force", False):
        try:
            cached = json.loads(meta_path.read_text(encoding="utf-8"))
            if cached.get("rows"):
                return {**base, "status": "cached", "message": "命中已有结果 (--force 可重算)", "rows": cached["rows"],
                        "n_bins": cached.get("n_bins"), "n_raw_cfgs": cached.get("n_raw_cfgs")}
        except Exception:
            pass

    out_dir.mkdir(parents=True, exist_ok=True)
    orch = MesonOrchestrator()
    ns = int(task["ns"])
    half = float(task["half"])
    x_axis = np.arange(ns, dtype=np.float64)
    is_single = bool(task["is_single_source"])

    corr_rows: List[dict] = []
    meff_rows: List[dict] = []
    scan_rows: List[dict] = []
    mass_rows: List[dict] = []
    n_bins_seen: Optional[int] = None
    n_cfgs_seen: Optional[int] = None

    for channel in CHANNEL_ORDER:
        means, errors, jk, n_bins, n_raw = orch.process_channel(
            input_dir=task["input_dir"],
            channel_configs=CHANNEL_CONFIGS[channel],
            binsize=int(task["binsizes"][channel]),
            num_lines=ns,
            thread_count=int(task["threads"]),
            is_single_source=is_single,
        )
        if jk.size == 0 or jk.shape[0] != ns:
            print(f"[WARN] {task['case_key']} ml={task['ml']:.4f} {channel}: 抽取失败", flush=True)
            continue

        n_bins_seen = int(n_bins)
        n_cfgs_seen = int(n_raw) or task.get("n_configs") or 0

        if task.get("symmetrize"):
            jk = symmetrize_matrix(jk)
            means = jk.mean(axis=1)
            n_jk = jk.shape[1]
            errors = np.sqrt((n_jk - 1) * np.sum((jk - means[:, None]) ** 2, axis=1) / n_jk)

        # ---- 关联函数 (驱动 C(x) 均值 + Jackknife 误差) ----
        pl.DataFrame(
            {"channel": [channel] * ns, "x": x_axis, "mean": means, "err": errors}
        ).write_csv(out_dir / f"corr_{channel}.csv")
        pl.DataFrame(jk).write_csv(out_dir / f"jk_{channel}.csv", include_header=False)
        corr_rows.append((channel, x_axis, means, errors))

        # ---- 有效质量 meff (对称点 = Ns/2) + 平台自动检测 ----
        meff_matrix, meff_mean, meff_err = compute_effective_mass_matrix_centered(jk, half)
        pl.DataFrame({"x": x_axis, "mean": meff_mean, "err": meff_err}).write_csv(
            out_dir / f"meff_{channel}.csv"
        )
        meff_rows.append((channel, meff_mean, meff_err))
        # 平台检测只用到 x < half 的有效质量 (x = half 处 cosh 比值无实根, 不参与平坦性检验)
        auto = detect_plateau_window(
            meff_matrix, end=int(half),
            chi2_dof_max=float(task.get("plateau_chi2", PLATEAU_CHI2_MAX)), min_points=4,
        )

        # ---- 平台窗口扫描 (补充信息, 全窗口落盘) ----
        scan = fit_mass_window_scan(x_axis, means, errors, [tuple(w) for w in task["windows"]], half)
        best = select_best_window(scan, min_points=MIN_WINDOW_POINTS, chi2_dof_max=CHI2_DOF_MAX)
        scan_mode = "scan_best"
        if best is None:
            best = select_best_window(scan, min_points=MIN_WINDOW_POINTS, chi2_dof_max=float("inf"))
            scan_mode = "scan_best_relaxed" if best is not None else "none"

        scan_table = pl.DataFrame(
            [
                {
                    "channel": channel,
                    "x_start": r["x_start"],
                    "x_end": r["x_end"],
                    "n_points": r["n_points"],
                    "mass": r["mass"],
                    "mass_err": r["mass_err"],
                    "chi2": r["chi2"],
                    "dof": r["dof"],
                    "chi2_dof": r["chi2_dof"],
                    "a": r["a"],
                    "a_err": r["a_err"],
                    "ok": r["ok"],
                }
                for r in scan
            ]
        )
        scan_table.write_csv(out_dir / f"mass_scan_{channel}.csv")
        for r in scan:
            scan_rows.append({"channel": channel, **{k: v for k, v in r.items() if k != "ok"}, "ok": bool(r["ok"])})

        empty_fit = {"mass": np.nan, "mass_err": np.nan, "chi2_dof": np.nan, "n_samples": jk.shape[1],
                     "n_fitted": 0, "a": np.nan, "a_err": np.nan}

        # 扫描最优窗口 (逐 Jackknife 样本拟合)
        if best is not None:
            scan_fit = fit_jackknife_mass_centered(
                x_axis, jk, errors, int(best["x_start"]), int(best["x_end"]), half
            )
            scan_start, scan_end = int(best["x_start"]), int(best["x_end"])
        else:
            scan_fit = dict(empty_fit)
            scan_start, scan_end = -1, -1

        # 镜像窗口 = 4.17 @ 48^3x16 既有分析窗口平移到本 ensemble 的对称点
        ref_start, ref_end = (int(v) for v in task["ref_windows"][channel])
        ref_fit = fit_jackknife_mass_centered(x_axis, jk, errors, ref_start, ref_end, half)

        # 自动平台窗口 (有效质量平坦性检测)
        if auto["ok"]:
            auto_fit = fit_jackknife_mass_centered(
                x_axis, jk, errors, int(auto["start"]), int(auto["end"]), half
            )
            auto_start, auto_end = int(auto["start"]), int(auto["end"])
        else:
            auto_fit = dict(empty_fit)
            auto_start, auto_end = -1, -1

        # 主结果按 --window-policy 选择
        policy = task.get("window_policy", "auto")
        if policy == "mirror":
            primary, primary_mode = ref_fit, "mirror_ref"
            primary_start, primary_end = ref_start, ref_end
        elif policy == "scan":
            primary, primary_mode = scan_fit, scan_mode
            primary_start, primary_end = scan_start, scan_end
        elif auto["ok"]:
            primary, primary_mode = auto_fit, "plateau_auto"
            primary_start, primary_end = auto_start, auto_end
        else:
            primary, primary_mode = scan_fit, f"{scan_mode}(no_plateau)"
            primary_start, primary_end = scan_start, scan_end

        def _f(d: dict, key: str) -> float:
            val = float(d.get(key, np.nan))
            return val if np.isfinite(val) else np.nan

        row = {
            **base,
            "channel": channel,
            "channel_index": CHANNEL_ORDER.index(channel),
            "binsize": int(task["binsizes"][channel]),
            "n_bins": n_bins_seen,
            "n_raw_cfgs": n_cfgs_seen,
            "half": half,
            "window_mode": primary_mode,
            "x_start": primary_start,
            "x_end": primary_end,
            "mass": _f(primary, "mass"),
            "mass_err": _f(primary, "mass_err"),
            "chi2_dof": _f(primary, "chi2_dof"),
            "n_fitted": int(primary.get("n_fitted", 0)),
            "mass_auto": _f(auto_fit, "mass"),
            "mass_err_auto": _f(auto_fit, "mass_err"),
            "auto_x_start": auto_start,
            "auto_x_end": auto_end,
            "auto_flat_chi2dof": float(auto["chi2_dof"]) if auto["ok"] else np.nan,
            "mass_ref": _f(ref_fit, "mass"),
            "mass_err_ref": _f(ref_fit, "mass_err"),
            "chi2_dof_ref": _f(ref_fit, "chi2_dof"),
            "ref_x_start": ref_start,
            "ref_x_end": ref_end,
            "mass_scan": _f(scan_fit, "mass"),
            "mass_err_scan": _f(scan_fit, "mass_err"),
            "chi2_dof_scan": _f(scan_fit, "chi2_dof"),
            "scan_x_start": scan_start,
            "scan_x_end": scan_end,
            "scan_mode": scan_mode,
            "temperature_mev": task["temperature_mev"],
            "temp_ref_mev": task["temp_ref_mev"],
            "a_inv_mev": task["a_inv_mev"],
            "mass_mev": _f(primary, "mass") * task["a_inv_mev"],
            "mass_mev_err": _f(primary, "mass_err") * task["a_inv_mev"],
            "mass_auto_mev": _f(auto_fit, "mass") * task["a_inv_mev"],
            "mass_auto_mev_err": _f(auto_fit, "mass_err") * task["a_inv_mev"],
            "mass_ref_mev": _f(ref_fit, "mass") * task["a_inv_mev"],
            "mass_ref_mev_err": _f(ref_fit, "mass_err") * task["a_inv_mev"],
            "mass_scan_mev": _f(scan_fit, "mass") * task["a_inv_mev"],
            "mass_scan_mev_err": _f(scan_fit, "mass_err") * task["a_inv_mev"],
        }
        mass_rows.append(row)

        pl.DataFrame([{k: v for k, v in row.items() if k not in ("out_dir",)}]).write_csv(
            out_dir / f"mass_{channel}.csv"
        )

    # ---- ensemble 级汇总文件 ----
    if corr_rows:
        pl.concat(
            [
                pl.DataFrame({"channel": [ch] * ns, "x": x, "mean": m, "err": e})
                for ch, x, m, e in corr_rows
            ]
        ).write_csv(out_dir / "correlators.csv")
    if meff_rows:
        pl.concat(
            [
                pl.DataFrame({"channel": [ch] * ns, "x": x_axis, "mean": m, "err": e})
                for ch, m, e in meff_rows
            ]
        ).write_csv(out_dir / "meff.csv")
    if scan_rows:
        pl.DataFrame(scan_rows).write_csv(out_dir / "mass_scan.csv")
    if mass_rows:
        pl.DataFrame(mass_rows).write_csv(out_dir / "mass.csv")

    meta = {
        "task": "b4.17 Nt-scan meson correlators & masses",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "case_key": task["case_key"],
        "label": task["label"],
        "order": task["order"],
        "ns": ns,
        "nt": task["nt"],
        "ml": task["ml"],
        "ml_index": task["ml_index"],
        "ms": task["ms"],
        "beta": task["beta"],
        "half": half,
        "source": task["source_tag"],
        "input_dir": task["input_dir"],
        "n_bins": n_bins_seen,
        "n_raw_cfgs": n_cfgs_seen,
        "binsizes": {ch: int(task["binsizes"][ch]) for ch in CHANNEL_ORDER},
        "fitter": LEAST_SQUARES_FITTER,
        "scan_windows": task["windows"],
        "window_policy": task.get("window_policy", "auto"),
        "plateau_chi2_max": float(task.get("plateau_chi2", PLATEAU_CHI2_MAX)),
        "ref_windows": task["ref_windows"],
        "n_channels_done": len(mass_rows),
        "elapsed_sec": round(time.time() - started, 3),
        "rows": mass_rows,
    }
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    return {
        **base,
        "status": "ok",
        "message": f"{len(mass_rows)} 个信道完成",
        "rows": mass_rows,
        "n_bins": n_bins_seen,
        "n_raw_cfgs": n_cfgs_seen,
        "elapsed_sec": meta["elapsed_sec"],
    }


# ------------------------------------------------------------------------------
# 结果汇总 (有序)
# ------------------------------------------------------------------------------
def _build_pivot(df: "pl.DataFrame", value_col: str, err_col: str,
                 win_cols: Tuple[str, str], chi2_col: str) -> List[dict]:
    """把长表转成 '每个 ensemble 一行 x 每个信道一列' 的宽表, 值形如 mass(err×10^5)"""
    pivot_rows: List[dict] = []
    for (order, ml_index), group in df.group_by(["order", "ml_index"], maintain_order=True):
        first = group.row(0, named=True)
        rec = {
            "order": order,
            "ns": first["ns"],
            "nt": first["nt"],
            "ms": first["ms"],
            "ml": first["ml"],
            "temperature_mev": round(first["temperature_mev"], 3),
            "a_inv_mev": round(first["a_inv_mev"], 3),
            "n_raw_cfgs": first["n_raw_cfgs"],
        }
        for ch in CHANNEL_ORDER:
            sub = group.filter(pl.col("channel") == ch)
            if sub.height == 0:
                rec[ch] = ""
                rec[f"{ch}_window"] = ""
                rec[f"{ch}_chi2dof"] = None
                continue
            r = sub.row(0, named=True)
            val, err = r[value_col], r[err_col]
            rec[ch] = f"{val:.5f}({err * 1e5:.0f})" if np.isfinite(val) else "nan"
            rec[f"{ch}_window"] = f"[{r[win_cols[0]]},{r[win_cols[1]]})"
            chi2 = r[chi2_col]
            rec[f"{ch}_chi2dof"] = None if not np.isfinite(chi2) else round(float(chi2), 4)
        pivot_rows.append(rec)
    return pivot_rows


def write_summaries(results: List[dict], source_tag: str, out_root: Path) -> List[Path]:
    """把某个 source 的所有结果按任务表顺序汇总落盘, 返回写出的文件路径"""
    rows: List[dict] = []
    for res in results:
        for row in res.get("rows", []):
            row = dict(row)
            row["source"] = source_tag
            rows.append(row)
    if not rows:
        return []

    df = pl.DataFrame(rows).with_columns(
        pl.col("mass_mev").round(4).alias("mass_mev"),
        pl.col("mass_mev_err").round(4).alias("mass_mev_err"),
        pl.col("mass_auto_mev").round(4).alias("mass_auto_mev"),
        pl.col("mass_auto_mev_err").round(4).alias("mass_auto_mev_err"),
        pl.col("mass_ref_mev").round(4).alias("mass_ref_mev"),
        pl.col("mass_ref_mev_err").round(4).alias("mass_ref_mev_err"),
        pl.col("mass_scan_mev").round(4).alias("mass_scan_mev"),
        pl.col("mass_scan_mev_err").round(4).alias("mass_scan_mev_err"),
    ).sort(["order", "ml_index", "channel_index"])

    src_root = out_root / source_tag
    src_root.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    summary_path = src_root / "meson_mass_summary.csv"
    df.write_csv(summary_path)
    written.append(summary_path)

    # 宽表 1: 主结果 (由 --window-policy 决定)
    pivot_path = src_root / "meson_mass_pivot.csv"
    pl.DataFrame(_build_pivot(df, "mass", "mass_err", ("x_start", "x_end"), "chi2_dof")).write_csv(pivot_path)
    written.append(pivot_path)

    # 宽表 1b: 自动平台窗口结果 (与主结果独立对照)
    pivot_auto_path = src_root / "meson_mass_pivot_auto.csv"
    pl.DataFrame(
        _build_pivot(df, "mass_auto", "mass_err_auto", ("auto_x_start", "auto_x_end"), "auto_flat_chi2dof")
    ).write_csv(pivot_auto_path)
    written.append(pivot_auto_path)

    # 宽表 2: 窗口扫描 chi2/dof 选优结果 (补充对照)
    pivot_scan_path = src_root / "meson_mass_pivot_scan.csv"
    pl.DataFrame(
        _build_pivot(df, "mass_scan", "mass_err_scan", ("scan_x_start", "scan_x_end"), "chi2_dof_scan")
    ).write_csv(pivot_scan_path)
    written.append(pivot_scan_path)

    # 宽表 3: 镜像窗口结果 (与 4.17 @ 48^3x16 既有发布口径对照)
    pivot_ref_path = src_root / "meson_mass_pivot_mirror.csv"
    pl.DataFrame(
        _build_pivot(df, "mass_ref", "mass_err_ref", ("ref_x_start", "ref_x_end"), "chi2_dof_ref")
    ).write_csv(pivot_ref_path)
    written.append(pivot_ref_path)

    return written


def print_ordered_table(results: List[dict]) -> None:
    rows: List[dict] = []
    for res in results:
        rows.extend(res.get("rows", []))
    if not rows:
        return
    df = pl.DataFrame(rows).sort(["order", "ml_index", "channel_index"])
    print("\n================ 介子质量 (lattice units, 平台拟合) ================")
    header = f"{'#':>2} {'ensemble':>10} {'ml':>7} {'T(MeV)':>7} " + " ".join(f"{ch:>16}" for ch in CHANNEL_ORDER)
    print(header)
    print("-" * len(header))
    seen = set()
    for r in df.iter_rows(named=True):
        key = (r["order"], r["ml_index"])
        if key in seen:
            continue
        seen.add(key)
        cells = []
        for ch in CHANNEL_ORDER:
            sub = df.filter((pl.col("order") == r["order"]) & (pl.col("ml_index") == r["ml_index"]) & (pl.col("channel") == ch))
            if sub.height:
                v = sub.row(0, named=True)
                cells.append(f"{v['mass']:>10.5f}({v['mass_err']*1e5:>4.0f})" if np.isfinite(v["mass"]) else f"{'nan':>16}")
            else:
                cells.append(f"{'-':>16}")
        print(f"{r['order']:>2} {r['label']:>10} {r['ml']:>7.4f} {r['temperature_mev']:>7.2f} " + " ".join(cells))
    print("(括号内为 Jackknife 误差 ×10^5)\n")


def write_manifest(all_results: Dict[str, List[dict]], out_root: Path, threads: int, elapsed: float) -> Path:
    used, extras = discover_b417_dirs()
    root = scan_root(out_root)

    def case_entry(case: ScanCase) -> dict:
        return {
            "order": case.order,
            "ns": case.ns,
            "nt": case.nt,
            "ms": case.ms,
            "mls": list(case.mls),
            "temperature_mev": round(case.temperature_mev, 3),
            "temp_ref_mev": case.temp_ref_mev,
            "dirs": [case.dirname(ml) for ml in case.mls],
        }

    manifest = {
        "task": "b4.17_nt_scan_meson",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "beta": BETA_VALUE,
        "beta_key": BETA_KEY,
        "ms": MS,
        "channels": list(CHANNEL_ORDER),
        "source_tags": [SOURCE_TAGS[s] for s in SOURCE_ORDER],
        "cases": [case_entry(c) for c in SCAN_CASES],
        "readin_dirs_used": used,
        "readin_dirs_b417_not_in_task": extras,
        "conventions": {
            "num_lines": "Ns (每信道 spatial:DIR* block 行数)",
            "correlator": "空间关联函数 C(x), x = 0 .. Ns-1 (16 个源位置环形移位平均)",
            "cosh_center": "Ns/2 (32->16, 36->18, 40->20, 48->24)",
            "fitter": f"{LEAST_SQUARES_FITTER} (强制 chi2 最小二乘, 不用 lsqfit 环境默认)",
            "jackknife_error": "sqrt((J-1) * sum((m_j - mean)^2) / J), J = bins 数",
            "reference_window": "由 4.17 @ 48^3x16 的 slice 表平移得到 (保持距对称点距离)",
        },
        "runtime": {"threads_per_task": threads, "elapsed_sec": round(elapsed, 2)},
        "results": {
            tag: [
                {
                    "order": r["order"],
                    "case": r["case_key"],
                    "ml": r["ml"],
                    "status": r["status"],
                    "n_bins": r.get("n_bins"),
                    "n_raw_cfgs": r.get("n_raw_cfgs"),
                    "out_dir": r.get("out_dir"),
                }
                for r in res_list
            ]
            for tag, res_list in all_results.items()
        },
    }
    root.mkdir(parents=True, exist_ok=True)
    path = root / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="b4.17 Nt 扫描: 关联函数 + 介子质量 计算任务",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source", choices=["multi", "single", "both"], default="both",
                        help="多源 (multi) / 单源 (single) / 两者都算")
    parser.add_argument("--cases", nargs="*", default=None,
                        help="只跑指定 ensemble, 如 --cases 32x12 48x18")
    parser.add_argument("--ml", nargs="*", type=float, default=None,
                        help="只跑指定 ml, 如 --ml 0.0020 0.0120")
    parser.add_argument("--workers", type=int, default=4, help="并行进程数")
    parser.add_argument("--threads", type=int, default=2, help="每个抽取任务内 C++ 线程数 (0 = 库默认)")
    parser.add_argument("--out-root", type=Path, default=None, help="结果根目录 (默认 output/meson_scan/b4.17)")
    parser.add_argument("--readin", type=Path, default=None, help="readin 数据根目录 (默认 data/readin)")
    parser.add_argument("--symmetrize", action="store_true", help="先做 C(x) <-> C(Ns-x) 对称化再拟合")
    parser.add_argument("--window-policy", choices=["auto", "mirror", "scan"], default="auto",
                        help="主质量窗口策略: auto=有效质量平坦性自动检测 (默认, 小格子也适用); "
                             "mirror=镜像既有 4.17@48^3x16 窗口; scan=cosh 拟合 chi2/dof 选优")
    parser.add_argument("--plateau-chi2", type=float, default=PLATEAU_CHI2_MAX,
                        help="平台自动检测的平坦性 chi2/dof 阈值 (越小窗口越靠近对称点)")
    parser.add_argument("--force", action="store_true", help="忽略已有结果, 强制重算")
    parser.add_argument("--dry-run", action="store_true", help="只打印任务清单, 不计算")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    t0 = time.time()

    out_root = scan_root(args.out_root)
    sources = {"multi": (False,), "single": (True,), "both": SOURCE_ORDER}[args.source]

    tasks = build_tasks(out_root, args.readin, sources, args.ml, args.cases, args.threads, args.window_policy)
    for task in tasks:
        task["force"] = args.force
        task["symmetrize"] = args.symmetrize
        task["plateau_chi2"] = args.plateau_chi2

    print(f"[TASK] beta=4.17 Nt 扫描 | ensemble {len({t['case_key'] for t in tasks})} 个 x "
          f"source {[SOURCE_TAGS[s] for s in sources]} | 任务数 {len(tasks)}")
    print(f"[TASK] 结果根目录: {out_root}")

    if args.dry_run:
        print("\n序号  ensemble     ml       T(MeV)   source      输入目录")
        for task in tasks:
            flag = "MISSING" if task["missing"] else f"cfgs={task['n_configs']}"
            print(f"{task['order']:>3}  {task['label']:>10}  {task['ml']:.4f}  {task['temperature_mev']:7.2f}  "
                  f"{task['source_tag']:>9}  {task['input_dir']}  [{flag}]")
        used, extras = discover_b417_dirs(args.readin)
        print(f"\n[INFO] 任务表目录数: {len(used)}")
        if extras:
            print(f"[INFO] 数字开头且含 b4.17 但不在任务表的目录 ({len(extras)} 个, 已跳过):")
            for name in extras:
                print(f"       - {name}")
        return 0

    results: List[dict] = []
    done = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [executor.submit(process_ensemble, task) for task in tasks]
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            done += 1
            if res["status"] in ("ok", "cached"):
                print(f"[{done:>2}/{len(tasks)}] {res['source']:<9} {res['label']:>10} ml={res['ml']:.4f} "
                      f"-> {res['status']} ({res['message']})", flush=True)
            else:
                print(f"[{done:>2}/{len(tasks)}] {res['source']:<9} {res['label']:>10} ml={res['ml']:.4f} "
                      f"-> {res['status']}: {res['message']}", flush=True)

    all_results: Dict[str, List[dict]] = {}
    for is_single in sources:
        tag = SOURCE_TAGS[is_single]
        subset = [r for r in results if r["source"] == tag]
        if not subset:
            continue
        all_results[tag] = subset
        print_ordered_table(subset)
        written = write_summaries(subset, tag, out_root)
        if written:
            for path in written:
                print(f"[OK] {tag} 汇总: {path}")

    n_missing = sum(1 for r in results if r["status"] == "missing")
    n_cases = len({t["case_key"] for t in tasks})
    manifest_path = write_manifest(all_results=all_results, out_root=out_root,
                                   threads=args.threads, elapsed=time.time() - t0)
    print(f"[OK] manifest: {manifest_path}")
    if n_missing:
        print(f"[WARN] {n_missing} 个任务因输入目录缺失被跳过")
    print(f"[DONE] {n_cases} 组格点 x {len(tasks) // max(n_cases, 1)} 个任务, 共 {len(tasks)} 个任务, "
          f"用时 {time.time() - t0:.1f}s")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
