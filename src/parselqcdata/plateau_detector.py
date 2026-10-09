#!/usr/bin/env python3
"""
src/parselqcdata/plateau_detector.py
--------------------------------------------------------------------------------
独立脚本与模块：有效质量平台窗口自动检测器 (Plateau Window Detector & Selector)
- 基于 Jackknife 样本的有效质量平坦性检验：
    D(x) = meff(x+1) - meff(x) 在 Jackknife 误差下与 0 一致性检验
    chi2/dof = \\sum (D / sigma_D)^2 / n_diff <= chi2_dof_max
- 自动挑选满足平坦性判据的最长稳定平台窗口
- 支持多种窗口决策策略 (auto, mirror, scan)
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence
import numpy as np
import polars as pl


def detect_plateau_window(
    meff_matrix: np.ndarray,
    end: int,
    chi2_dof_max: float = 2.0,
    min_points: int = 3,
    min_start: int = 1,
) -> dict:
    """
    基于 Jackknife 样本的有效质量平台自动检测。

    判据: 窗口 [s, end) 内所有相邻差分 D(x) = meff(x+1) - meff(x) 在 Jackknife
    误差下与 0 一致, 即 flatness chi2/dof <= chi2_dof_max。因为窗口随 s 减小而变长、
    chi2/dof 单调变差, 第一个通过的 s 就是 "离对称点最近的平坦窗口" (最长平台)。

    Parameters:
        meff_matrix: 有效质量样本矩阵 (nrow x ncol)
        end: 拟合上限点 (通常为 Ns/2)
        chi2_dof_max: 判定平坦性的最大 chi2/dof 阈值
        min_points: 平台最少包含点数
        min_start: 最小起始切片编号

    Returns:
        {start, end, chi2_dof, n_diffs, n_points, mass_est, mass_err_est, ok}
    """
    meff_matrix = np.asarray(meff_matrix, dtype=np.float64)
    nrow, ncol = meff_matrix.shape
    end = int(min(end, nrow))
    result = {
        "start": None,
        "end": end,
        "chi2_dof": np.nan,
        "n_diffs": 0,
        "n_points": 0,
        "mass_est": np.nan,
        "mass_err_est": np.nan,
        "ok": False,
    }
    if ncol == 0 or end - int(min_start) < int(min_points):
        return result

    # 相邻差分与其 Jackknife 误差
    D = meff_matrix[1:end, :] - meff_matrix[: end - 1, :]
    finite = np.isfinite(D)
    cnt = finite.sum(axis=1)
    dmean = np.where(cnt > 0, np.nansum(D, axis=1) / np.maximum(cnt, 1), np.nan)
    resid = np.where(finite, D - np.where(np.isfinite(dmean), dmean, np.nan)[:, None], 0.0)
    dvar = np.sum(resid ** 2, axis=1)
    derr = np.sqrt(np.maximum(cnt - 1, 0) * dvar / np.maximum(cnt, 1))

    for s in range(int(min_start), end - int(min_points) + 1):
        idx = slice(s, end - 1)  # 差分编号 x = s ... end-2
        dd, de = dmean[idx], derr[idx]
        ok = np.isfinite(dd) & np.isfinite(de) & (de > 0.0)
        n_diff = int(ok.sum())
        if n_diff < 2:
            continue
        chi2 = float(np.sum((dd[ok] / de[ok]) ** 2))
        chi2_dof = chi2 / n_diff
        if chi2_dof <= chi2_dof_max:
            window = meff_matrix[s:end, :]
            m_mean = np.nanmean(window)
            m_err = np.sqrt(
                (ncol - 1) * np.nansum((window - m_mean) ** 2) / ncol
            )
            result.update(
                {
                    "start": s,
                    "end": end,
                    "chi2_dof": chi2_dof,
                    "n_diffs": n_diff,
                    "n_points": end - s,
                    "mass_est": float(m_mean),
                    "mass_err_est": float(m_err),
                    "ok": True,
                }
            )
            return result
    return result


def select_best_window(
    scan_results: Sequence[dict],
    chi2_max: float = 3.0,
    min_points: int = 3,
) -> Optional[dict]:
    """
    从窗口扫描结果中选出最佳拟合窗口：
    1. 优先筛选 chi2/dof <= chi2_max 且拟合成功的窗口；
    2. 在满足条件的候选窗口中选择 chi2/dof 最低者；
    3. 若无窗口满足上限，则兜底返回所有有效拟合中 chi2/dof 最低者。
    """
    valid = [
        r for r in scan_results
        if r.get("ok") and np.isfinite(r.get("chi2_dof", np.nan)) and r.get("n_points", 0) >= min_points
    ]
    if not valid:
        return None

    acceptable = [r for r in valid if r["chi2_dof"] <= chi2_max]
    pool = acceptable if acceptable else valid
    return min(pool, key=lambda r: r["chi2_dof"])


def main() -> None:
    parser = argparse.ArgumentParser(description="有效质量平台窗口检测脚本")
    parser.add_argument("--csv", type=str, required=True, help="有效质量矩阵 CSV 文件")
    parser.add_argument("--end", type=int, default=24, help="拟合上限切片位置 (默认: 24)")
    parser.add_argument("--chi2-max", type=float, default=2.0, help="平坦性最大 chi2/dof 阈值 (默认: 2.0)")
    parser.add_argument("--min-points", type=int, default=3, help="平台最少点数 (默认: 3)")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.is_file():
        print(f"[ERROR] 文件不存在: {csv_path}", file=sys.stderr)
        sys.exit(1)

    df = pl.read_csv(csv_path)
    mat = df.to_numpy()
    res = detect_plateau_window(mat, end=args.end, chi2_dof_max=args.chi2_max, min_points=args.min_points)
    print("--- 平台检测结果 ---")
    for k, v in res.items():
        print(f"  {k:<14}: {v}")


if __name__ == "__main__":
    main()
