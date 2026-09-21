#!/usr/bin/env python3
"""
tools/condensate_orchestrator.py

Chiral Condensate Pipeline C ABI Bridge
Loads liblqcd_condensate.dylib / libparselqcdata.dylib
"""

import ctypes
from pathlib import Path
from typing import Optional, Tuple

class CondensateOrchestrator:
    """Python bridge for Chiral Condensate C++ pipeline."""
    def __init__(self, library_path: Optional[str] = None) -> None:
        if library_path is None:
            library_path = self._find_library()
        self.library_path = str(library_path)
        self._lib = ctypes.CDLL(self.library_path)
        self._setup_bindings()
    
    def _find_library(self) -> Path:
        current_dir = Path(__file__).resolve().parent
        project_root = current_dir.parent
        candidates = [
            project_root / "build" / "libparselqcdata.dylib",
            project_root / "lib" / "libparselqcdata.dylib",
            project_root / "lib" / "liblqcd_condensate.so",
        ]
        for c in candidates:
            if c.exists():
                return c
        raise FileNotFoundError(
            "Condensate shared library not found. Checked:\n"
            + "\n".join(f"  - {c}" for c in candidates)
        )
    
    def _setup_bindings(self) -> None:
        # int run_chiral_condensate_c_api(...)
        self._lib.run_chiral_condensate_c_api.argtypes = [
            ctypes.c_char_p,                 # base_dir
            ctypes.c_double,                 # m_light
            ctypes.c_double,                 # m_strange
            ctypes.c_double,                 # m_residual
            ctypes.c_double,                 # zm_factor
            ctypes.POINTER(ctypes.c_double), # out_mean
            ctypes.POINTER(ctypes.c_double), # out_error
            ctypes.POINTER(ctypes.c_int),    # out_num_cfgs
        ]
        self._lib.run_chiral_condensate_c_api.restype = ctypes.c_int

        # int run_chiral_condensate_full_c_api(...)
        if hasattr(self._lib, "run_chiral_condensate_full_c_api"):
            self._lib.run_chiral_condensate_full_c_api.argtypes = [
                ctypes.c_char_p,                 # base_dir
                ctypes.c_double,                 # m_light
                ctypes.c_double,                 # m_strange
                ctypes.c_double,                 # m_residual
                ctypes.c_double,                 # zm_factor
                ctypes.POINTER(ctypes.c_double), # out_pbp_sub_mean
                ctypes.POINTER(ctypes.c_double), # out_pbp_sub_error
                ctypes.POINTER(ctypes.c_double), # out_pbp_l_mean
                ctypes.POINTER(ctypes.c_double), # out_pbp_l_error
                ctypes.POINTER(ctypes.c_double), # out_pbp_s_mean
                ctypes.POINTER(ctypes.c_double), # out_pbp_s_error
                ctypes.POINTER(ctypes.c_int),    # out_num_cfgs
            ]
            self._lib.run_chiral_condensate_full_c_api.restype = ctypes.c_int
    
    def process_condensate(
        self,
        base_dir: str,
        m_light: float = 0.001,
        m_strange: float = 0.050,
        m_residual: float = 0.0003,
        zm_factor: float = 1.0
    ) -> Tuple[float, float, int]:
        out_mean = ctypes.c_double(0.0)
        out_error = ctypes.c_double(0.0)
        out_num_cfgs = ctypes.c_int(0)
    
        base_dir_bytes = str(base_dir).encode('utf-8')
    
        ret = self._lib.run_chiral_condensate_c_api(
            base_dir_bytes,
            ctypes.c_double(m_light),
            ctypes.c_double(m_strange),
            ctypes.c_double(m_residual),
            ctypes.c_double(zm_factor),
            ctypes.byref(out_mean),
            ctypes.byref(out_error),
            ctypes.byref(out_num_cfgs)
        )
    
        if ret != 0:
            raise RuntimeError(f"run_chiral_condensate_c_api failed with code {ret}")
    
        return float(out_mean.value), float(out_error.value), int(out_num_cfgs.value)

    def process_condensate_full(
        self,
        base_dir: str,
        m_light: float = 0.001,
        m_strange: float = 0.050,
        m_residual: float = 0.0003,
        zm_factor: float = 1.0
    ) -> dict:
        if not hasattr(self._lib, "run_chiral_condensate_full_c_api"):
            mean, err, n = self.process_condensate(base_dir, m_light, m_strange, m_residual, zm_factor)
            return {
                "pbp_sub_mean": mean,
                "pbp_sub_error": err,
                "pbp_l_mean": 0.0,
                "pbp_l_error": 0.0,
                "pbp_s_mean": 0.0,
                "pbp_s_error": 0.0,
                "num_cfgs": n,
            }

        out_sub_mean = ctypes.c_double(0.0)
        out_sub_error = ctypes.c_double(0.0)
        out_l_mean = ctypes.c_double(0.0)
        out_l_error = ctypes.c_double(0.0)
        out_s_mean = ctypes.c_double(0.0)
        out_s_error = ctypes.c_double(0.0)
        out_num_cfgs = ctypes.c_int(0)

        base_dir_bytes = str(base_dir).encode('utf-8')

        ret = self._lib.run_chiral_condensate_full_c_api(
            base_dir_bytes,
            ctypes.c_double(m_light),
            ctypes.c_double(m_strange),
            ctypes.c_double(m_residual),
            ctypes.c_double(zm_factor),
            ctypes.byref(out_sub_mean),
            ctypes.byref(out_sub_error),
            ctypes.byref(out_l_mean),
            ctypes.byref(out_l_error),
            ctypes.byref(out_s_mean),
            ctypes.byref(out_s_error),
            ctypes.byref(out_num_cfgs)
        )

        if ret != 0:
            raise RuntimeError(f"run_chiral_condensate_full_c_api failed with code {ret}")

        return {
            "pbp_sub_mean": float(out_sub_mean.value),
            "pbp_sub_error": float(out_sub_error.value),
            "pbp_l_mean": float(out_l_mean.value),
            "pbp_l_error": float(out_l_error.value),
            "pbp_s_mean": float(out_s_mean.value),
            "pbp_s_error": float(out_s_error.value),
            "num_cfgs": int(out_num_cfgs.value),
        }
