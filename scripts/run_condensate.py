#!/usr/bin/env python3
"""
scripts/run_condensate.py
--------------------------------------------------------------------------------
手征凝聚 (Chiral Condensate) 端到端全量分析与物理标度运行脚本:
1. 扫描 data/readin/ 中所有可用格点系综 (涵盖温度扫描与有限体积标度系列)
2. 调度 C++ 核心库多线程并发抽取裸光/奇夸克手征凝聚与 Jackknife 重采样
3. 扣除残余质量散度并执行 Zm 物理重整化
4. 独立隔离保存单系综结果，并维护全局总表 all_ensembles_condensate.csv 与 results_rm_beta.txt
5. 可选自动生成物理曲线图表 (--plot)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    DEFAULT_READIN_DIR,
    OUTPUT_CONDENSATE_DIR,
    CONDENSATE_CONFIGS,
    calculate_residual_mass,
)
from src.parselqcdata.condensate_pipeline import CondensatePipeline


def run_condensate_analysis(
    readin_dir: Path = DEFAULT_READIN_DIR,
    output_dir: Path = OUTPUT_CONDENSATE_DIR,
    target_ensemble: str | None = None,
    do_plot: bool = False,
) -> None:
    t0 = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("==========================================================================")
    print("[CONDENSATE] 启动手征凝聚数据抽取与重整化流水线...")
    print(f"  - 输入目录: {readin_dir}")
    print(f"  - 输出目录: {output_dir}")
    print("==========================================================================")

    pipeline = CondensatePipeline()

    if target_ensemble:
        print(f"[INFO] 正在单独处理指定系综: {target_ensemble}")
        rec = pipeline.process_ensemble(target_ensemble, base_readin_dir=readin_dir, output_dir=output_dir)
        print(f"  [OK] {rec['dataset_name']}: <pbp_sub>={rec['pbp_rm']:+.6e} +/- {rec['pbp_rm_err']:.6e} ({rec['num_cfgs']} cfgs)")
    elif readin_dir.exists():
        results = pipeline.process_all_ensembles(base_readin_dir=readin_dir, output_dir=output_dir)
        print(f"\n[OK] 成功处理 {len(results)} 组格点系综，耗时: {time.time() - t0:.2f}s")
    else:
        print(f"[WARN] 数据目录 {readin_dir} 不存在，从预置 CONDENSATE_CONFIGS 导出标准基准表...")
        out_file = output_dir / "results_rm_beta.txt"
        records = []
        for cfg in CONDENSATE_CONFIGS:
            beta_str = str(cfg["beta"])
            mres = float(cfg.get("mres", calculate_residual_mass(float(beta_str))))
            records.append({
                "beta": beta_str,
                "mres": mres,
                "pbpl": float(cfg["pbp_l"]),
                "pbpl_err": float(cfg["pbp_l_err"]),
                "pbps": float(cfg["pbp_s"]),
                "pbps_err": float(cfg["pbp_s_err"]),
                "pbp_rm": float(cfg["pbp_rm"]),
                "pbp_rm_err": float(cfg["pbp_rm_err"]),
            })
        import polars as pl
        pl.DataFrame(records).write_csv(out_file, separator="\t")
        print(f"  -> 已写出基准表: {out_file}")

    if do_plot:
        print("\n[INFO] 正在生成手征凝聚全景图表 (plot_condensate.py)...")
        from scripts.plot_condensate import main as plot_main
        plot_main()


def main() -> None:
    parser = argparse.ArgumentParser(description="手征凝聚抽取与重整化分析工具")
    parser.add_argument(
        "--readin-dir",
        type=Path,
        default=DEFAULT_READIN_DIR,
        help="数据根目录 (默认 data/readin)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_CONDENSATE_DIR,
        help="结果输出根目录 (默认 output/condensate)",
    )
    parser.add_argument(
        "--ensemble",
        type=str,
        default=None,
        help="指定仅分析特定系综目录名",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="计算完成后自动生成图表",
    )
    args = parser.parse_args()

    run_condensate_analysis(
        readin_dir=args.readin_dir,
        output_dir=args.output_dir,
        target_ensemble=args.ensemble,
        do_plot=args.plot,
    )


if __name__ == "__main__":
    main()
