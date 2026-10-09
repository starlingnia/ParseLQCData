#!/usr/bin/env python3
"""
scripts/tasks/lcp/03_export_lcp.py
--------------------------------------------------------------------------------
LCP (常物理线) 子任务 3: 质量因子 Zm^2 重整化与 LCP 标准格式导出
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 读取由子任务 2 计算出的格点手征磁化率数据 (output/condensate/results_susceptibility.csv)。
2. 根据规范耦合常数 Beta 查询质量重整化因子 Zm(beta)，应用物理重整化:
     chi_ren = (chi_scaled / 1e6) / (Zm^2)   [单位: GeV^2]
3. 严格按照高能物理实验组与 ana 规范导出标准格式文本文件:
   - results_susceptibility.txt
   - results_susceptibility_lcp.txt
   - results_susceptibility_all.txt
   - results_susceptibility_scaling.txt
   标头为:
   # beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error
4. 输出至标准 output/LCP/，保证目录规范统一。

【底层调用的 C++ 功能与关联】:
- 输入数据直接源自 C++ 磁化率计算组件 target("condensate_analysis") 与 target("statistics") 的产物。
- 物理常数与质量重整化参数与 include/core/PhysicsSetup.h 严格对应。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    OUTPUT_CONDENSATE_DIR,
    OUTPUT_LCP_DIR,
)
from src.parselqcdata.susceptibility_pipeline import SusceptibilityPipeline


def run_export_lcp_task(
    input_csv: Path | None = None,
    output_lcp_dir: Path = OUTPUT_LCP_DIR,
) -> None:
    print(f"\n[LCP-TASK-03] 启动 LCP 标准格式文本产物导出与 Zm^2 重整化...")
    
    pipeline = SusceptibilityPipeline(output_dir=OUTPUT_CONDENSATE_DIR)
    
    if input_csv is not None and input_csv.exists():
        df = pl.read_csv(input_csv)
    else:
        df = pipeline.get_full_scan_results()

    # 规范导出至 output/LCP/
    target_dirs = [output_lcp_dir]
    pipeline.export_lcp_outputs(df, target_dirs=target_dirs)
    print(f"[OK] LCP 规范产物已成功导出至: {[str(d) for d in target_dirs]}")


def main() -> None:
    parser = argparse.ArgumentParser(description="LCP 标准格式导出小脚本")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=OUTPUT_CONDENSATE_DIR / "results_susceptibility.csv",
        help="输入磁化率 CSV 文件",
    )
    parser.add_argument(
        "--output-lcp-dir",
        type=Path,
        default=OUTPUT_LCP_DIR,
        help="LCP 产物输出目录",
    )
    args = parser.parse_args()

    in_csv = args.input_csv if args.input_csv.exists() else None
    run_export_lcp_task(
        input_csv=in_csv,
        output_lcp_dir=args.output_lcp_dir,
    )


if __name__ == "__main__":
    main()
