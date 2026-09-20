#!/usr/bin/env python3
"""
tools/lqcd_orchestrator.py
--------------------------------------------------------------------------------
Unified bridge wrapping MesonOrchestrator and CondensateOrchestrator
--------------------------------------------------------------------------------
"""

from typing import List, Optional, Tuple
import numpy as np

from tools.meson_orchestrator import MesonOrchestrator
from tools.condensate_orchestrator import CondensateOrchestrator


class LQCDOrchestrator:
    """Unified Orchestrator delegating to specialized Meson and Condensate bridges."""

    def __init__(self, library_path: Optional[str] = None) -> None:
        self.meson = MesonOrchestrator(library_path)
        self.condensate = CondensateOrchestrator(library_path)

    def process_channel(
        self,
        input_dir: str,
        channel_configs: List[dict],
        binsize: int = 4,
        num_lines: int = 48,
        thread_count: int = 0,
        is_single_source: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
        return self.meson.process_channel(
            input_dir=input_dir,
            channel_configs=channel_configs,
            binsize=binsize,
            num_lines=num_lines,
            thread_count=thread_count,
            is_single_source=is_single_source
        )

    def process_chiral_condensate(
        self,
        base_dir: str,
        m_light: float = 0.001,
        m_strange: float = 0.050,
        m_residual: float = 0.0003,
        zm_factor: float = 1.0
    ) -> Tuple[float, float, int]:
        return self.condensate.process_condensate(
            base_dir=base_dir,
            m_light=m_light,
            m_strange=m_strange,
            m_residual=m_residual,
            zm_factor=zm_factor
        )

    def execute_meson_analysis(
        self,
        input_dir: str,
        channel_mappings: Optional[List[dict]] = None,
        channel_configs: Optional[List[dict]] = None,
        binsize: int = 4,
        num_lines: int = 48,
        thread_count: int = 0,
        is_single_source: bool = False,
        return_folded_jk: bool = True
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.meson.execute_meson_analysis(
            input_dir=input_dir,
            channel_mappings=channel_mappings,
            channel_configs=channel_configs,
            binsize=binsize,
            num_lines=num_lines,
            thread_count=thread_count,
            is_single_source=is_single_source,
            return_folded_jk=return_folded_jk
        )
