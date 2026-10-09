#!/usr/bin/env python3
"""
scripts/lcp/run_lcp_susceptibility.py
--------------------------------------------------------------------------------
【LCP 专区】常物理线 (Line of Constant Physics) 手征磁化率分析脚本:
1. 驱动 C++ 核心库提取纯轻夸克微观随机源向量
2. 计算单构型无偏两体乘积并执行构型级 Jackknife 误差分析
3. 准确折算四维时空体积因子 F_vol 与温度标度因子 F_scaled，施加 Zm^2 质量重整化
4. 导出 LCP 规范结果至 output/LCP/
5. 可选自动生成 LCP 手征磁化率物理图表 (--plot)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import DEFAULT_READIN_DIR, OUTPUT_CONDENSATE_DIR, OUTPUT_LCP_DIR
from src.parselqcdata.susceptibility_pipeline import SusceptibilityPipeline


def run_lcp_susceptibility(
    readin_dir: Path | None = None,
    output_dir: Path = OUTPUT_CONDENSATE_DIR,
    do_plot: bool = False,
) -> None:
    t0 = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("==========================================================================")
    print("[LCP SUSCEPTIBILITY] 启动常物理线手征磁化率数据抽取与物理标度流水线...")
    print(f"  - 输入目录: {readin_dir if readin_dir else '基准/本地验证模式'}")
    print(f"  - 输出目录: {output_dir}")
    print("==========================================================================")

    pipeline = SusceptibilityPipeline(output_dir=output_dir)

    if readin_dir and readin_dir.exists():
        print(f"[INFO] 正在从真实构型目录执行扫描: {readin_dir}")
        df = pipeline.run_cluster_scan(readin_dir)
    else:
        df = pipeline.get_full_scan_results(force_recompute=True)

    print(f"[OK] LCP 磁化率计算完成 (共 {df.height} 组数据)，耗时: {time.time() - t0:.2f}s")

    if do_plot:
        print("\n[INFO] 正在生成 LCP 手征磁化率图表 (plot_lcp_susceptibility.py)...")
        from scripts.lcp.plot_lcp_susceptibility import main as plot_main
        plot_main()


def main() -> None:
    parser = argparse.ArgumentParser(description="【LCP 专区】常物理线手征磁化率分析工具")
    parser.add_argument(
        "--readin-dir",
        type=Path,
        default=None,
        help="数据根目录 (默认自动定位或使用基准数据)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_CONDENSATE_DIR,
        help="结果输出根目录 (默认 output/condensate)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="计算完成后自动生成图表",
    )
    args = parser.parse_args()

    readin = args.readin_dir if (args.readin_dir and args.readin_dir.exists()) else (
        DEFAULT_READIN_DIR if DEFAULT_READIN_DIR.exists() else None
    )

    run_lcp_susceptibility(
        readin_dir=readin,
        output_dir=args.output_dir,
        do_plot=args.plot,
    )


if __name__ == "__main__":
    main()
