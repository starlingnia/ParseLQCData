#!/usr/bin/env python3
"""
scripts/tasks/lcp/01_extract_condensate.py
--------------------------------------------------------------------------------
LCP (常物理线) 子任务 1: 手征凝聚数据抽取与残余质量减除重整化
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 扫描格点系综目录中测量的 XML 随机源，提取光夸克与奇夸克手征凝聚:
     pbp_l = <psibar psi>_l,  pbp_s = <psibar psi>_s
2. 执行残余质量散度减除 (Residual Mass Subtraction):
     <pbp_sub> = pbp_l - [(ml + mres) / (ms + mres)] * pbp_s
3. 施加 Zm 乘法质量重整化，通过 Jackknife 评估各系综统计不确定度。
4. 独立隔离保存单系综结果至 output/condensate/ensembles/<name>/，
   并维护全局汇总总表 all_ensembles_condensate.csv、scaling_beta4.17.csv
   与严格对齐的历史基准表 results_rm_beta.txt。

【底层调用的 C++ 功能】:
- C++ 动态链接库: lib/libparselqcdata.dylib / lib/liblqcd_condensate.dylib
- C++ 导出 C ABI 函数: run_chiral_condensate_c_api(...) (定义于 tools/Services/CondensateService.cpp)
- C++ 底层静态核心组件:
  * target("condensate_analysis"): CondensateExtractor (src/CondensateAnalysis/CondensateExtractor.cpp)
    负责 PsibarPsi 随机源 XML 抽取与多随机源均值计算
  * target("core"): CondensatePipeline (src/core/CondensatePipeline.cpp)
    负责多线程并发系综调度与参数推导
  * target("statistics"): Resampling (src/Statistics/Resampling.cpp)
    负责构型级与系综级 Jackknife 重采样误差传递
  * target("iodata"): FastParser (src/IOdata/FastParser.h)
    负责极速无拷贝 XML 浮点数解析
- 对应 C++ 原生 CLI 任务:
  * ./bin/ParseLQCData condensate [dataset_path_or_beta] [ml] [ms] [mres] [zm]
  * ./bin/ParseLQCData condensate_all [readin_dir]
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    CONDENSATE_CONFIGS,
    DEFAULT_READIN_DIR,
    OUTPUT_CONDENSATE_DIR,
    calculate_residual_mass,
)
from src.parselqcdata.condensate_pipeline import CondensatePipeline


def run_extract_condensate_task(
    readin_dir: Path = DEFAULT_READIN_DIR,
    output_dir: Path = OUTPUT_CONDENSATE_DIR,
    target_ensemble: str | None = None,
) -> None:
    t0 = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[LCP-TASK-01] 启动手征凝聚数据抽取与重整化...")
    print(f"  - 输入目录: {readin_dir}")
    print(f"  - 输出目录: {output_dir}")

    pipeline = CondensatePipeline()

    if target_ensemble:
        print(f"[INFO] 正在单独处理指定系综: {target_ensemble}")
        rec = pipeline.process_ensemble(target_ensemble, base_readin_dir=readin_dir, output_dir=output_dir)
        print(f"  [OK] {rec['dataset_name']}: <pbp_sub>={rec['pbp_rm']:+.6e} +/- {rec['pbp_rm_err']:.6e} ({rec['num_cfgs']} cfgs)")
    elif readin_dir.exists():
        results = pipeline.process_all_ensembles(base_readin_dir=readin_dir, output_dir=output_dir)
        print(f"[OK] 成功处理 {len(results)} 组格点系综，耗时: {time.time() - t0:.2f}s")
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

    print(f"[OK] 手征凝聚抽取子任务执行完成！耗时: {time.time() - t0:.2f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description="LCP 手征凝聚抽取与残余质量减除小脚本 (基于 C++ 核心库)")
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
    args = parser.parse_args()

    run_extract_condensate_task(
        readin_dir=args.readin_dir,
        output_dir=args.output_dir,
        target_ensemble=args.ensemble,
    )


if __name__ == "__main__":
    main()
