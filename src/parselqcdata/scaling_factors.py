#!/usr/bin/env python3
"""
src/parselqcdata/scaling_factors.py
--------------------------------------------------------------------------------
独立脚本与模块：格点物理标度因子、质量重整化常数 Zm 与目录元数据解析器
- 计算四维格点体积因子 F_vol = Ns^3 * Nt
- 计算无量纲标度因子 F_scaled = Ns^3 * Nt^3 * T^2 = F_vol * (Nt * T)^2
- 耦合常数 Beta -> 夸克质量重整化常数 Zm 查询
- 从系综目录名称中解析格点几何尺寸与参数
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Dict, Optional, Tuple, Union

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def compute_scaling_factors(ns: int, nt: int, temp_mev: float) -> Tuple[float, float]:
    """
    计算手征磁化率所需的体积因子与标度因子：
    - factor_vol = Ns^3 * Nt
    - factor_scaled = Ns^3 * Nt^3 * T^2 = factor_vol * (Nt * T)^2 (单位: MeV^2)

    Parameters:
        ns: 空间格点点数
        nt: 时间格点点数
        temp_mev: 物理温度 (MeV)

    Returns:
        (factor_vol, factor_scaled): 四维格点体积与无量纲标度因子
    """
    f_vol = float((ns**3) * nt)
    f_scaled = float((ns**3) * (nt**3) * (temp_mev**2))
    return f_vol, f_scaled


def get_zm_factor(beta_val: Union[float, str]) -> float:
    """
    根据规范耦合常数 Beta 查询获取质量重整化因子 Zm(beta)。
    默认查阅 docs.physics_setup.MRES_TABLE 标准表，缺失时平滑兜底为 1.0。
    """
    try:
        from docs.physics_setup import MRES_TABLE
        b_float = float(beta_val)
        b_str = f"{b_float:.2f}" if abs(b_float - 4.405) > 1e-4 else "4.405"
        if str(beta_val) in MRES_TABLE:
            return float(MRES_TABLE[str(beta_val)]["zm"])
        if b_str in MRES_TABLE:
            return float(MRES_TABLE[b_str]["zm"])
    except Exception:
        pass
    return 1.0


def parse_ensemble_meta_from_dir(dir_name: str) -> Dict[str, Union[int, float, str]]:
    """
    从目录名称中解析格点几何尺寸 (Ns, Nt) 与耦合常数 Beta。
    支持格式:
      - L48T16beta4.18ms0.037265m0.001022
      - 48x16b4.18
      - L32T12beta4.17
      - L40T16_beta4.17
    """
    meta: Dict[str, Union[int, float, str]] = {
        "ns": 48,
        "nt": 16,
        "beta": "4.17",
    }

    # 1. 匹配时空几何 Ns, Nt
    geom_match = re.search(r"L(\d+)T(\d+)", dir_name, re.IGNORECASE)
    if geom_match:
        meta["ns"] = int(geom_match.group(1))
        meta["nt"] = int(geom_match.group(2))
    else:
        geom_x = re.search(r"(\d+)x(\d+)", dir_name)
        if geom_x:
            meta["ns"] = int(geom_x.group(1))
            meta["nt"] = int(geom_x.group(2))

    # 2. 匹配 Beta
    beta_match = re.search(r"b(?:eta)?([0-9.]+)", dir_name, re.IGNORECASE)
    if beta_match:
        meta["beta"] = beta_match.group(1)

    return meta


def is_valid_condensate_dir(dir_path: Path | str) -> bool:
    """
    检查目录是否为有效的手征凝聚测量目录。
    自动识别并跳过强子/介子关联函数目录（例如包含 Output/test1_lhadrons_* 且无 PsibarPsi 的目录）。
    """
    path = Path(dir_path)
    if not path.is_dir():
        return False
    # 介子关联函数输出目录通常包含 Output/ 且没有 meas.*
    if (path / "Output").exists():
        has_meas = any((path / m).is_dir() for m in ["test_condensate"] if (path / m).exists()) or list(path.glob("meas.*"))
        if not has_meas:
            return False

    # 检查根目录下是否有 meas.*/PsibarPsi
    for m in path.glob("meas.*"):
        if (m / "PsibarPsi").exists():
            return True

    # 检查子目录下是否有 test_condensate/meas.*/PsibarPsi
    sub_cand = path / "test_condensate"
    if sub_cand.exists():
        for m in sub_cand.glob("meas.*"):
            if (m / "PsibarPsi").exists():
                return True

    return False


def main() -> None:
    parser = argparse.ArgumentParser(description="格点物理标度与参数查询工具")
    parser.add_argument("--ns", type=int, default=48, help="空间格点 Ns (默认: 48)")
    parser.add_argument("--nt", type=int, default=16, help="时间格点 Nt (默认: 16)")
    parser.add_argument("--temp", type=float, default=153.31, help="温度 T (MeV, 默认: 153.31)")
    parser.add_argument("--beta", type=str, default="4.17", help="规范耦合常数 beta (默认: 4.17)")
    parser.add_argument("--dir-name", type=str, default=None, help="解析指定系综目录名称")
    args = parser.parse_args()

    if args.dir_name:
        meta = parse_ensemble_meta_from_dir(args.dir_name)
        print(f"解析目录 '{args.dir_name}':")
        for k, v in meta.items():
            print(f"  {k}: {v}")
        args.ns = int(meta["ns"])
        args.nt = int(meta["nt"])
        args.beta = str(meta["beta"])

    f_vol, f_scaled = compute_scaling_factors(args.ns, args.nt, args.temp)
    zm = get_zm_factor(args.beta)

    print("\n--- 物理标度参数 ---")
    print(f"Ns: {args.ns}, Nt: {args.nt}, T: {args.temp:.2f} MeV, Beta: {args.beta}")
    print(f"体积因子 F_vol   : {f_vol:.6e} (lattice units)")
    print(f"标度因子 F_scaled: {f_scaled:.6e} (MeV^2)")
    print(f"重整化因子 Zm    : {zm:.6f}")
    print(f"重整化除数 Zm^2  : {zm**2:.6f}")


if __name__ == "__main__":
    main()
