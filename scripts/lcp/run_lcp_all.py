#!/usr/bin/env python3
"""
scripts/lcp/run_lcp_all.py
--------------------------------------------------------------------------------
【LCP 专区】常物理线 (Line of Constant Physics) 全流程一键端到端总入口:
一键依次调度 LCP 的三大核心物理分析与可视化流程：
1. 手征凝聚 (Chiral Condensate): 8 个 LCP 系综抽取、残余质量相减、Zm 重整化
2. 手征磁化率 (Chiral Susceptibility): 无偏二次量、Jackknife 统计、除 Zm^2 导出规范 output/LCP/
3. 介子温度演化 (Meson Mass & Correlators): 48^3x16 各温度介子有效质量与平台拟合
4. 综合物理可视化出图 (--plot 默认启用)
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

from scripts.lcp.run_lcp_condensate import run_lcp_condensate
from scripts.lcp.run_lcp_susceptibility import run_lcp_susceptibility
from scripts.lcp.run_lcp_meson import run_lcp_meson_pipeline
from scripts.lcp.plot_lcp_condensate import main as plot_condensate_main
from scripts.lcp.plot_lcp_susceptibility import main as plot_susceptibility_main
from scripts.lcp.plot_lcp_meson import main as plot_meson_main
from docs.physics_setup import DEFAULT_READIN_DIR, OUTPUT_ROOT, OUTPUT_CONDENSATE_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="【LCP 专区】常物理线全流程端到端运行总脚本")
    parser.add_argument("--readin-dir", type=Path, default=DEFAULT_READIN_DIR, help="数据输入根目录")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_ROOT, help="输出根目录")
    parser.add_argument("--workers", type=int, default=4, help="并发工作进程数")
    parser.add_argument("--skip-condensate", action="store_true", help="跳过手征凝聚分析")
    parser.add_argument("--skip-susceptibility", action="store_true", help="跳过手征磁化率分析")
    parser.add_argument("--skip-meson", action="store_true", help="跳过介子分析")
    parser.add_argument("--no-plot", action="store_true", help="计算完成后不自动绘图")
    args = parser.parse_args()

    t_start = time.time()
    print("==========================================================================")
    print("🚀 启动常物理线 (Line of Constant Physics, LCP) 全流程分析")
    print(f"  - 输入目录: {args.readin_dir}")
    print(f"  - 输出目录: {args.output_dir}")
    print("==========================================================================")

    # 1. 手征凝聚
    if not args.skip_condensate:
        print("\n>>> [1/3] 执行 LCP 手征凝聚分析 (run_lcp_condensate)...")
        cond_out = args.output_dir / "condensate"
        run_lcp_condensate(readin_dir=args.readin_dir, output_dir=cond_out, do_plot=False)

    # 2. 手征磁化率
    if not args.skip_susceptibility:
        print("\n>>> [2/3] 执行 LCP 手征磁化率分析 (run_lcp_susceptibility)...")
        readin_sucep = args.readin_dir if (args.readin_dir and args.readin_dir.exists()) else None
        run_lcp_susceptibility(readin_dir=readin_sucep, output_dir=args.output_dir / "condensate", do_plot=False)

    # 3. 介子关联函数与质量
    if not args.skip_meson:
        print("\n>>> [3/3] 执行 LCP 介子热演化分析 (run_lcp_meson)...")
        run_lcp_meson_pipeline(
            readin_dir=args.readin_dir,
            out_root=args.output_dir,
            max_workers=args.workers,
            do_plot=False,
        )

    # 4. 可选自动绘图
    if not args.no_plot:
        print("\n==========================================================================")
        print("🎨 正在生成全套 LCP 物理曲线与报告图表...")
        print("==========================================================================")
        try:
            print("  -> 生成 LCP 手征凝聚图表...")
            plot_condensate_main()
        except Exception as e:
            print(f"  [WARN] 手征凝聚绘图跳过: {e}")

        try:
            print("  -> 生成 LCP 手征磁化率标度图表...")
            plot_susceptibility_main()
        except Exception as e:
            print(f"  [WARN] 手征磁化率绘图跳过: {e}")

        try:
            print("  -> 生成 LCP 介子热演化对比图表...")
            plot_meson_main()
        except Exception as e:
            print(f"  [WARN] 介子绘图跳过: {e}")

    elapsed = time.time() - t_start
    print("\n==========================================================================")
    print(f"✅ LCP 全流程分析完成! 总耗时: {elapsed:.2f} 秒")
    print(f"  - LCP 规范结果输出目录: {args.output_dir}/LCP/")
    print(f"  - 手征凝聚输出目录    : {args.output_dir}/condensate/")
    print(f"  - 介子关联函数与质量表: {args.output_dir}/simulateresult/")
    print(f"  - 物理图表输出目录    : docs/figures/")
    print("==========================================================================")


if __name__ == "__main__":
    main()
