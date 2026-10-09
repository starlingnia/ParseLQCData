#!/usr/bin/env python3
"""
src/parselqcdata/unbiased_quadratic.py
--------------------------------------------------------------------------------
独立脚本与模块：单构型无偏二次交叉估计器 (Unbiased Quadratic Estimator)
- 数学定义:
    Obar = (1/k) * \\sum_{i=1}^k O_i
    O2bar = [1/(k(k-1))] * \\sum_{i \\neq j} O_i O_j
          = (1/k) * \\sum_{i=1}^k [ 1/(k-1) * O_i * (S - O_i) ]
- 消除有限随机源方差噪声方差污染
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence, Tuple
import numpy as np


def compute_unbiased_quadratic(vals: Sequence[float]) -> Tuple[float, float]:
    """
    计算单个规范构型上的无偏手征凝聚均值与两点方均关联估计：
    - Obar = (1/k) * \\sum_{i=1}^k O_i
    - O2bar = [1/(k(k-1))] * \\sum_{i \\neq j} O_i O_j
            = (1/k) * \\sum_{i=1}^k [ 1/(k-1) * O_i * (S - O_i) ]
    消除了由于有限随机源数目 k 带来的随机噪声方差对算符平方的内积偏差。

    Parameters:
        vals: 单构型上测得的 k 个随机源标量值序列 (k >= 2)

    Returns:
        (obar, obar_sq_unbiased): 算符线性均值与两体无偏方均值
    """
    k = len(vals)
    if k <= 1:
        raise ValueError(f"计算两体无偏关联至少需要 2 个随机源向量，当前仅提供 {k} 个")

    total_s = float(sum(vals))
    obar = total_s / k
    # 向量化无偏方均估计
    obar_sq_unbiased = float(np.mean([(1.0 / (k - 1)) * x * (total_s - x) for x in vals]))
    return obar, obar_sq_unbiased


def main() -> None:
    parser = argparse.ArgumentParser(description="单构型无偏两体乘积计算脚本")
    parser.add_argument("values", nargs="*", type=float, help="随机源测量值列表 (浮点数)")
    parser.add_argument("--file", type=str, default=None, help="从文本文件读取数值 (每行一个浮点数)")
    args = parser.parse_args()

    vals = list(args.values)
    if args.file:
        fpath = Path(args.file)
        if not fpath.is_file():
            print(f"[ERROR] 文件不存在: {fpath}", file=sys.stderr)
            sys.exit(1)
        for line in fpath.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    vals.append(float(line.split()[0]))
                except ValueError:
                    pass

    if len(vals) < 2:
        print("[ERROR] 至少需要提供 2 个测量数值进行无偏两体乘积计算", file=sys.stderr)
        sys.exit(1)

    obar, o2bar = compute_unbiased_quadratic(vals)
    print(f"随机源数量 k: {len(vals)}")
    print(f"Obar (均值)   : {obar:+.10e}")
    print(f"O2bar (无偏两点): {o2bar:+.10e}")
    biased_sq = obar ** 2
    delta = o2bar - biased_sq
    print(f"有偏平方 Obar^2: {biased_sq:+.10e}")
    print(f"方差抵消修正量 : {delta:+.10e}")


if __name__ == "__main__":
    main()
