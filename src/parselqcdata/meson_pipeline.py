"""
src/parselqcdata/meson_pipeline.py
--------------------------------------------------------------------------------
Python 整合层：Meson 强子关联函数端到端分析基础管道
- 导入 docs/physics_setup 中的全局配置
- 调度 tools/lqcd_orchestrator.py 驱动 C++ 核心库
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import CHANNEL_CONFIGS, get_binsize, BETAS, CHANNELS, DEFAULT_READIN_DIR
from tools.meson_orchestrator import MesonOrchestrator


class MesonPipeline:
    """Meson 数据处理与统计分析管道"""

    def __init__(self, orchestrator: Optional[MesonOrchestrator] = None):
        self.orchestrator = orchestrator or MesonOrchestrator()

    def process_channel(
        self,
        beta: str,
        channel: str,
        base_readin_dir: Optional[str] = None,
        is_single_source: bool = False,
        binsize: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        处理单个 (beta, channel) 组合，返回: (means, errors, folded_jk_matrix)
        """
        if base_readin_dir is None:
            base_readin_dir = str(DEFAULT_READIN_DIR)
        if channel not in CHANNEL_CONFIGS:
            raise ValueError(f"未知信道: {channel}，可选信道: {list(CHANNEL_CONFIGS.keys())}")

        if binsize is None:
            binsize = get_binsize(beta, channel, is_single_source)

        input_dir = Path(base_readin_dir) / f"48x16b4.{beta}" / "Output"
        if not input_dir.exists():
            raise FileNotFoundError(f"输入构型目录不存在: {input_dir}")

        mappings = CHANNEL_CONFIGS[channel]
        means, errors, folded_jk = self.orchestrator.execute_meson_analysis(
            input_dir=str(input_dir),
            channel_mappings=mappings,
            binsize=binsize,
            num_lines=48,
            return_folded_jk=True,
            is_single_source=is_single_source
        )

        return means, errors, folded_jk

    def run_ensemble(
        self,
        beta: str,
        base_readin_dir: Optional[str] = None,
        is_single_source: bool = False
    ) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """批量处理特定 Beta 下的所有 6 种物理信道"""
        if base_readin_dir is None:
            base_readin_dir = str(DEFAULT_READIN_DIR)
        results = {}
        for ch in CHANNELS:
            results[ch] = self.process_channel(beta, ch, base_readin_dir, is_single_source)
        return results


PRESET_CHANNELS = CHANNEL_CONFIGS


def process_channel(
    channel: str,
    beta: str,
    binsize: Optional[int] = None,
    is_single_source: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """快捷单信道分析接口，返回 (means, errors)"""
    pipeline = MesonPipeline()
    means, errors, _ = pipeline.process_channel(beta, channel, binsize=binsize, is_single_source=is_single_source)
    return means, errors

