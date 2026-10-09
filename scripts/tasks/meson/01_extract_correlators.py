#!/usr/bin/env python3
"""
scripts/tasks/meson/01_extract_correlators.py
--------------------------------------------------------------------------------
介子测量子任务 1: 强子空间关联函数抽取与对称化折叠
--------------------------------------------------------------------------------
【物理与算法说明】:
1. 从格点 QCD 原始测量文本 (Output/test1_lhadrons_*) 扫描并解析介子关联函数。
2. 投影到 6 大物理信道 (AV, S, Tt, PS, Xt, Vec)，支持多源 (multisrc) 与单源 (singlesrc)。
3. 在空间/时间方向执行折叠对称化 (fold): C(t) 与 C(Nt - t) 对称求均值。
4. 运用 Jackknife 重采样算法构建统计协方差与误差矩阵。
5. 输出标准产物:
   - save_{channel}.csv : 完整全时隙均值与误差
   - sym_{channel}.csv  : 对称折叠前半段均值与误差 (0 ~ 24)
   - dr_{channel}.csv   : (num_lines x n_bins) Jackknife 子样本折叠矩阵
   - err_{channel}.csv  : 各时隙 Jackknife 误差向量

【底层调用的 C++ 功能】:
- C++ 动态链接库: lib/libparselqcdata.dylib / lib/liblqcd_meson.dylib
- C++ 导出 C ABI 函数: run_meson_pipeline_c_api(...) (定义于 tools/Services/MesonService.cpp)
- C++ 底层静态核心组件:
  * target("meson_analysis"): MesonExtractor (src/MesonAnalysis/MesonExtractor.cpp) - 负责狄拉克矩阵投影与极化态张量求和
  * target("core"): MesonPipeline (src/core/MesonPipeline.cpp) - 多线程并行调度、流水线并发执行
  * target("iodata"): DirectoryScanner / FastParser (src/IOdata/DirectoryScanner.cpp) - 极速文本扫描与零拷贝解析
  * target("statistics"): Resampling (src/Statistics/Resampling.cpp) - Jackknife 折叠与方差计算
- 对应 C++ 原生 CLI 任务:
  * ./bin/ParseLQCData meson_multi [beta] [channel] [binsize]
  * ./bin/ParseLQCData meson_single [beta] [channel] [binsize]
  * ./bin/ParseLQCData meson_all [beta] [multi|single]
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
from typing import Sequence

import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    BETAS,
    CHANNELS,
    CHANNEL_CONFIGS,
    DEFAULT_READIN_DIR,
    OUTPUT_ROOT,
    get_binsize,
)
from tools.meson_orchestrator import MesonOrchestrator


def extract_correlators_for_channel(
    orch: MesonOrchestrator,
    input_dir: Path,
    beta: str,
    channel: str,
    is_single: bool,
    pick_dir: Path,
) -> None:
    """调用 C++ 引擎抽取单个信道的强子关联函数并落盘标准 CSV"""
    configs = CHANNEL_CONFIGS[channel]
    binsize = get_binsize(beta, channel, is_single)

    # 调用底层 C++ 核心库通过 C ABI 完成高并发数据流抽取与折叠
    means, errors, folded_jk, n_bins, n_cfgs = orch.process_channel(
        input_dir=str(input_dir),
        channel_configs=configs,
        binsize=binsize,
        num_lines=48,
        thread_count=0,
        is_single_source=is_single,
    )

    # 规范产物落盘
    pick_dir.mkdir(parents=True, exist_ok=True)
    pl.DataFrame({"mean": means, "err": errors}).write_csv(pick_dir / f"save_{channel}.csv")
    pl.DataFrame({"mean": means[:25], "err": errors[:25]}).write_csv(pick_dir / f"sym_{channel}.csv")
    pl.DataFrame(folded_jk).write_csv(pick_dir / f"dr_{channel}.csv", include_header=False)
    pl.DataFrame(folded_jk).write_parquet(pick_dir / f"dr_{channel}.parquet")
    pl.DataFrame(errors).write_csv(pick_dir / f"err_{channel}.csv", include_header=False)


def run_extract_task(
    betas: Sequence[str] = BETAS,
    sources: Sequence[bool] = (False, True),
    readin_dir: Path = DEFAULT_READIN_DIR,
    output_root: Path = OUTPUT_ROOT,
) -> None:
    t0_all = time.time()
    orch = MesonOrchestrator()

    for is_single in sources:
        mode_str = "singlesrc" if is_single else "multisrc"
        dir_pick = "pickdata-singlesrc" if is_single else "pickdata"

        print(f"\n[MESON-TASK-01] 正在抽取 {mode_str.upper()} 关联函数...")
        for beta in betas:
            input_dir = readin_dir / f"48x16b4.{beta}" / "Output"
            if not input_dir.exists():
                print(f"  [WARN] 目录不存在，跳过: {input_dir}")
                continue

            pick_dir = output_root / dir_pick / f"b4.{beta}"
            t0_beta = time.time()
            for ch in CHANNELS:
                extract_correlators_for_channel(
                    orch=orch,
                    input_dir=input_dir,
                    beta=beta,
                    channel=ch,
                    is_single=is_single,
                    pick_dir=pick_dir,
                )
            print(f"  [OK] Beta 4.{beta} ({mode_str}): 6 信道关联函数抽取完毕 (耗时: {time.time() - t0_beta:.2f}s)")

    print(f"[OK] 关联函数抽取子任务执行完成！总耗时: {time.time() - t0_all:.2f}s")


def main() -> None:
    parser = argparse.ArgumentParser(description="介子关联函数抽取与折叠小脚本 (基于 C++ 核心引擎)")
    parser.add_argument(
        "--source",
        choices=["all", "multi", "single"],
        default="all",
        help="数据源类型: multi, single 或 all",
    )
    parser.add_argument(
        "--beta",
        type=str,
        default="all",
        help="指定 Beta (例如 17 或 4.17)，或 'all'",
    )
    parser.add_argument(
        "--readin-dir",
        type=Path,
        default=DEFAULT_READIN_DIR,
        help="数据源根目录",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="输出根目录",
    )
    args = parser.parse_args()

    if args.source == "multi":
        sources = [False]
    elif args.source == "single":
        sources = [True]
    else:
        sources = [False, True]

    if args.beta == "all":
        betas = BETAS
    else:
        b_clean = args.beta.replace("b4.", "").replace("4.", "")
        betas = [b_clean]

    run_extract_task(
        betas=betas,
        sources=sources,
        readin_dir=args.readin_dir,
        output_root=args.output_root,
    )


if __name__ == "__main__":
    main()
