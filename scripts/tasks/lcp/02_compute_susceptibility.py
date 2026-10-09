#!/usr/bin/env python3
"""
scripts/tasks/lcp/02_compute_susceptibility.py
--------------------------------------------------------------------------------
LCP (常物理线) 子任务 2: 手征磁化率 (Chiral Susceptibility) 无偏二次估计与标度
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 物理理论定义:
     chi_disc = (V_3 / T) * [ <(psibar psi)^2> - <psibar psi>^2 ]
   其中时空四维体积 V_4 = V_3 / T = Ns^3 * Nt * a^4。
2. 单构型无偏二次交叉估计器 (Unbiased Quadratic Estimator):
     O_bar = (1/k) * sum_{i=1}^k O_i
     O2_bar = [1/(k*(k-1))] * sum_{i != j} O_i * O_j
   消除有限随机源数目 k 引入的对角虚假噪声方差抬升。
3. 构型级 Jackknife 误差分析:
     a_r = JK(O_bar)_r,  b_r = JK(O2_bar)_r
     chi_r = b_r - a_r^2
4. 物理标度因子折算:
     - 四维格点体积标度: F_vol = Ns^3 * Nt,  chi_vol = F_vol * chi_unscaled
     - 物理温度标度:     F_scaled = Ns^3 * Nt^3 * T^2,  chi_scaled = F_scaled * chi_unscaled (MeV^2)
5. 导出全量汇总表至 output/condensate/results_susceptibility.csv 与 all_ensembles_susceptibility.csv。

【底层调用的 C++ 功能】:
- C++ 静态核心组件:
  * target("condensate_analysis"): CondensateExtractor (src/CondensateAnalysis/CondensateExtractor.cpp)
    其成员函数 compute_unbiased_pbp 与 extract_sources_from_xml 负责随机源交叉积与对角消噪
  * target("statistics"): Resampling (src/Statistics/Resampling.cpp)
    负责构型级 Jackknife 重采样方差分析与无偏统计量折算
  * target("iodata"): FastParser (src/IOdata/FastParser.h)
    负责微观随机源向量的快速解析
- 对应 C++ 原生 CLI 任务:
  * ./bin/ParseLQCData susceptibility [readin_dir] [output_dir]
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

from docs.physics_setup import DEFAULT_READIN_DIR, OUTPUT_CONDENSATE_DIR
from src.parselqcdata.susceptibility_pipeline import SusceptibilityPipeline


def run_compute_susceptibility_task(
    readin_dir: Path | None = None,
    output_dir: Path = OUTPUT_CONDENSATE_DIR,
) -> None:
    t0 = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[LCP-TASK-02] 启动手征磁化率无偏交叉估计与标度计算...")
    print(f"  - 输入目录: {readin_dir if readin_dir else '自动定位/基准验证模式'}")
    print(f"  - 输出目录: {output_dir}")

    pipeline = SusceptibilityPipeline(output_dir=output_dir)

    if readin_dir and readin_dir.exists():
        print(f"[INFO] 正在扫描格点系综构型目录: {readin_dir}")
        df = pipeline.run_cluster_scan(readin_dir)
    else:
        df = pipeline.get_full_scan_results(force_recompute=True)

    print(f"[OK] 磁化率计算完成 (共 {df.height} 组系综)，耗时: {time.time() - t0:.2f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description="LCP 手征磁化率无偏估计与标度计算小脚本 (基于 C++ 核心库)")
    parser.add_argument(
        "--readin-dir",
        type=Path,
        default=None,
        help="数据根目录 (默认自动探测)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_CONDENSATE_DIR,
        help="结果输出根目录 (默认 output/condensate)",
    )
    args = parser.parse_args()

    readin = args.readin_dir if (args.readin_dir and args.readin_dir.exists()) else (
        DEFAULT_READIN_DIR if DEFAULT_READIN_DIR.exists() else None
    )

    run_compute_susceptibility_task(
        readin_dir=readin,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
