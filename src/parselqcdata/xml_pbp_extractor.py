#!/usr/bin/env python3
"""
src/parselqcdata/xml_pbp_extractor.py
--------------------------------------------------------------------------------
独立脚本与模块：XML 手征凝聚随机源测量值抽取器
- 从单个 XML 测量文件中使用正则抽取 <pbp>(real, imag)</pbp> 的实部数值
- 批量扫描构型目录并提取全部随机源向量
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import List, Optional, Tuple


def extract_pbp_from_xml(file_path: Path | str) -> Optional[float]:
    """
    使用高效正则从 XML 文件中提取 <pbp>(real, imag)</pbp> 的实部数值。

    Parameters:
        file_path: XML 文件路径

    Returns:
        float: 提取的实部浮点数；若提取失败则返回 None
    """
    path = Path(file_path)
    if not path.is_file():
        return None
    try:
        content = path.read_text(encoding="utf-8", errors="ignore")
        match = re.search(r"<pbp>\(([^,]+),", content)
        if match:
            return float(match.group(1))
    except Exception:
        pass
    return None


def extract_pbp_from_directory(
    meas_dir: Path | str,
    target_quark: str = "light",
) -> Tuple[List[float], List[str]]:
    """
    扫描测量目录 (如 meas.XXXXXX/) 下的 XML 文件，提取特定夸克类型的全部随机源数值。

    Parameters:
        meas_dir: 包含测量 XML 文件的目录
        target_quark: 目标夸克类别，"light" 对应 pbp_l，"strange" 对应 pbp_s

    Returns:
        (values, file_names): 提取的数值列表与对应的文件名列表
    """
    dir_path = Path(meas_dir)
    if not dir_path.is_dir():
        return [], []

    pbp_sub = dir_path / "PsibarPsi"
    scan_dir = pbp_sub if pbp_sub.is_dir() else dir_path

    xml_files = sorted(scan_dir.glob("*.xml"))
    values: List[float] = []
    file_names: List[str] = []

    for f in xml_files:
        name_lower = f.name.lower()
        if target_quark == "light" and "strange" in name_lower:
            continue
        if target_quark == "strange" and "strange" not in name_lower:
            continue

        val = extract_pbp_from_xml(f)
        if val is not None:
            values.append(val)
            file_names.append(f.name)

    return values, file_names


def main() -> None:
    parser = argparse.ArgumentParser(description="XML 手征凝聚测量值提取脚本")
    parser.add_argument("path", type=str, help="XML 文件或测量目录路径")
    parser.add_argument("--quark", type=str, default="light", choices=["light", "strange"], help="夸克类型 (默认: light)")
    args = parser.parse_args()

    target = Path(args.path)
    if target.is_file():
        val = extract_pbp_from_xml(target)
        if val is not None:
            print(f"{target.name}: {val:+.8e}")
        else:
            print(f"[ERROR] 无法从 {target} 中提取 <pbp> 数值", file=sys.stderr)
            sys.exit(1)
    elif target.is_dir():
        vals, names = extract_pbp_from_directory(target, target_quark=args.quark)
        print(f"从目录 {target} 提取到 {len(vals)} 个 {args.quark} 夸克测量值:")
        for name, v in zip(names, vals):
            print(f"  {name}: {v:+.8e}")
    else:
        print(f"[ERROR] 路径不存在: {target}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
