import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

from tools.lqcd_orchestrator import LQCDOrchestrator

# 预设的标准物理信道多方向映射表 (与原项目 scripts/run.sh 严格对齐)
PRESET_CHANNELS: Dict[str, List[dict]] = {
    'AV': [
        {'type': 'AVector2', 'dir': 'DIRX'},
        {'type': 'AVector3', 'dir': 'DIRX'},
        {'type': 'AVector3', 'dir': 'DIRY'},
        {'type': 'AVector1', 'dir': 'DIRY'},
        {'type': 'AVector1', 'dir': 'DIRZ'},
        {'type': 'AVector2', 'dir': 'DIRZ'},
    ],
    'S': [
        {'type': 'TVector4', 'dir': 'DIRZ'},
        {'type': 'TVector4', 'dir': 'DIRX'},
        {'type': 'TVector4', 'dir': 'DIRY'},
    ],
    'Tt': [
        {'type': 'TVector1', 'dir': 'DIRX'},
        {'type': 'TVector3', 'dir': 'DIRZ'},
        {'type': 'TVector2', 'dir': 'DIRY'},
    ],
    'PS': [
        {'type': 'TAVector4', 'dir': 'DIRZ'},
        {'type': 'TAVector4', 'dir': 'DIRX'},
        {'type': 'TAVector4', 'dir': 'DIRY'},
    ],
    'Xt': [
        {'type': 'TAVector1', 'dir': 'DIRX'},
        {'type': 'TAVector3', 'dir': 'DIRZ'},
        {'type': 'TAVector2', 'dir': 'DIRY'},
    ],
    'Vec': [
        {'type': 'Vector2', 'dir': 'DIRX'},
        {'type': 'Vector3', 'dir': 'DIRX'},
        {'type': 'Vector3', 'dir': 'DIRY'},
        {'type': 'Vector1', 'dir': 'DIRY'},
        {'type': 'Vector1', 'dir': 'DIRZ'},
        {'type': 'Vector2', 'dir': 'DIRZ'},
    ],
}


def process_channel(
    channel_name: str,
    beta_str: str,
    binsize: int = 4,
    base_readin_dir: str = "/Users/junxiongnie/code/ana/dat/readin",
    orchestrator: Optional[LQCDOrchestrator] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    高层核心接口：驱动 C++ 底层高并发解析与统计计算
    """
    if channel_name not in PRESET_CHANNELS:
        raise ValueError(f"未知信道: {channel_name}，可选信道: {list(PRESET_CHANNELS.keys())}")

    mappings = PRESET_CHANNELS[channel_name]
    input_dir = os.path.join(base_readin_dir, f"48x16b4.{beta_str}", "Output")

    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"输入数据目录不存在: {input_dir}")

    if orchestrator is None:
        orchestrator = LQCDOrchestrator()

    means, errors, _ = orchestrator.execute_meson_analysis(
        input_dir=input_dir,
        channel_mappings=mappings,
        binsize=binsize,
        num_lines=48,
        thread_count=0
    )

    return means, errors
