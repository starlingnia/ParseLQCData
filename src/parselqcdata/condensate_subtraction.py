#!/usr/bin/env python3
"""
src/parselqcdata/condensate_subtraction.py
--------------------------------------------------------------------------------
独立脚本与模块：手征凝聚残余质量相减 (Residual Mass Subtraction) 与 Zm 物理重整化
- 物理公式:
    \\langle \\bar{\\psi}\\psi \\rangle_{\\text{sub}} = Z_m \\cdot \\left[ \\langle \\bar{\\psi}\\psi \\rangle_l - \\frac{m_l + m_{\\text{res}}}{m_s + m_{\\text{res}}} \\langle \\bar{\\psi}\\psi \\rangle_s \\right]
- 消除 Domain Wall Fermion 有限第五维引入的加法手征破缺残余质量发散
- 施加有限重整化常数 Zm
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, Optional, Tuple, Union
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def compute_condensate_subtraction(
    pbp_l: float,
    pbp_s: float,
    ml: float,
    ms: float,
    mres: float,
    zm: float = 1.0,
    pbp_l_err: Optional[float] = None,
    pbp_s_err: Optional[float] = None,
) -> Tuple[float, Optional[float]]:
    """
    计算残余质量扣除并施加 Zm 重整化后的手征凝聚：
    pbp_sub = Zm * (pbp_l - [(ml + mres) / (ms + mres)] * pbp_s)

    Parameters:
        pbp_l: 裸光夸克手征凝聚
        pbp_s: 裸奇夸克手征凝聚
        ml: 光夸克质量
        ms: 奇夸克质量
        mres: 残余夸克质量
        zm: 质量重整化因子 (默认: 1.0)
        pbp_l_err: 光夸克手征凝聚统计误差 (可选)
        pbp_s_err: 奇夸克手征凝聚统计误差 (可选)

    Returns:
        (pbp_sub, pbp_sub_err): 重整化手征凝聚均值与误差 (若未提供分量误差则返回 None)
    """
    ratio = (ml + mres) / (ms + mres)
    raw_sub = pbp_l - ratio * pbp_s
    val_sub = zm * raw_sub

    err_sub = None
    if pbp_l_err is not None and pbp_s_err is not None:
        raw_var = (pbp_l_err**2) + ((ratio * pbp_s_err)**2)
        err_sub = zm * float(np.sqrt(raw_var))

    return val_sub, err_sub


def main() -> None:
    parser = argparse.ArgumentParser(description="手征凝聚残余质量扣除与 Zm 重整化工具")
    parser.add_argument("--pbpl", type=float, required=True, help="裸光夸克凝聚 pbp_l")
    parser.add_argument("--pbps", type=float, required=True, help="裸奇夸克凝聚 pbp_s")
    parser.add_argument("--ml", type=float, required=True, help="光夸克质量 ml")
    parser.add_argument("--ms", type=float, required=True, help="奇夸克质量 ms")
    parser.add_argument("--mres", type=float, required=True, help="残余夸克质量 mres")
    parser.add_argument("--zm", type=float, default=1.0, help="重整化因子 Zm (默认: 1.0)")
    parser.add_argument("--pbpl-err", type=float, default=None, help="pbp_l 误差")
    parser.add_argument("--pbps-err", type=float, default=None, help="pbp_s 误差")
    args = parser.parse_args()

    val, err = compute_condensate_subtraction(
        args.pbpl, args.pbps, args.ml, args.ms, args.mres, args.zm, args.pbpl_err, args.pbps_err
    )

    print("--- 手征凝聚重整化计算结果 ---")
    print(f"pbp_l: {args.pbpl:.8e}, pbp_s: {args.pbps:.8e}")
    print(f"ml: {args.ml:.6e}, ms: {args.ms:.6e}, mres: {args.mres:.6e}, Zm: {args.zm:.6f}")
    if err is not None:
        print(f"<pbp_sub>: {val:+.8e} +/- {err:.8e}")
    else:
        print(f"<pbp_sub>: {val:+.8e}")


if __name__ == "__main__":
    main()
