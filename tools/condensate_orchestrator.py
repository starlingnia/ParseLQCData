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
            project_root / "lib" / "libparselqcdata.dylib",
            project_root / "lib" / "liblqcd_condensate.dylib",
            project_root / "build" / "libparselqcdata.dylib",
            project_root / "lib" / "libparselqcdata.so",
            project_root / "lib" / "liblqcd_condensate.so",
            project_root / "build" / "libparselqcdata.so",
        ]
        errors = []
        for c in candidates:
            if c.exists():
                try:
                    test_lib = ctypes.CDLL(str(c))
                    if hasattr(test_lib, "run_chiral_condensate_c_api"):
                        return c
                except OSError as e:
                    errors.append(f"{c.name}: {e}")
        # Fallback to first existing if test loading did not identify one
        for c in candidates:
            if c.exists():
                return c
        raise FileNotFoundError(
            "Condensate shared library not found. Checked:\n"
            + "\n".join(f"  - {c}" for c in candidates)
            + (f"\nLoad errors:\n" + "\n".join(f"  - {e}" for e in errors) if errors else "")
        )

    def is_available(self) -> bool:
        return self._lib is not None

    def has_susceptibility(self) -> bool:
        return hasattr(self._lib, "run_chiral_susceptibility_c_api")
    
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

        if hasattr(self._lib, "run_chiral_susceptibility_c_api"):
            self._lib.run_chiral_susceptibility_c_api.argtypes = [
                ctypes.c_char_p,                 # base_dir
                ctypes.c_int,                    # ns
                ctypes.c_int,                    # nt
                ctypes.c_double,                 # temp_mev
                ctypes.POINTER(ctypes.c_double), # out_mean_unscaled
                ctypes.POINTER(ctypes.c_double), # out_error_unscaled
                ctypes.POINTER(ctypes.c_double), # out_mean_vol_scaled
                ctypes.POINTER(ctypes.c_double), # out_error_vol_scaled
                ctypes.POINTER(ctypes.c_double), # out_mean_scaled
                ctypes.POINTER(ctypes.c_double), # out_error_scaled
                ctypes.POINTER(ctypes.c_int),    # out_num_cfgs
            ]
            self._lib.run_chiral_susceptibility_c_api.restype = ctypes.c_int
    
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

    def process_susceptibility(
        self,
        base_dir: str,
        ns: int = 48,
        nt: int = 16,
        temp_mev: float = 157.0,
    ) -> dict:
        """
        通过 C++ 高性能多线程引擎提取纯轻夸克手征磁化率并进行物理标度计算。
        """
        if not hasattr(self._lib, "run_chiral_susceptibility_c_api"):
            raise NotImplementedError("Dynamic library does not export run_chiral_susceptibility_c_api")

        out_mean_unscaled = ctypes.c_double(0.0)
        out_error_unscaled = ctypes.c_double(0.0)
        out_mean_vol_scaled = ctypes.c_double(0.0)
        out_error_vol_scaled = ctypes.c_double(0.0)
        out_mean_scaled = ctypes.c_double(0.0)
        out_error_scaled = ctypes.c_double(0.0)
        out_num_cfgs = ctypes.c_int(0)

        base_dir_bytes = str(base_dir).encode('utf-8')

        ret = self._lib.run_chiral_susceptibility_c_api(
            base_dir_bytes,
            ctypes.c_int(ns),
            ctypes.c_int(nt),
            ctypes.c_double(temp_mev),
            ctypes.byref(out_mean_unscaled),
            ctypes.byref(out_error_unscaled),
            ctypes.byref(out_mean_vol_scaled),
            ctypes.byref(out_error_vol_scaled),
            ctypes.byref(out_mean_scaled),
            ctypes.byref(out_error_scaled),
            ctypes.byref(out_num_cfgs)
        )

        if ret != 0:
            raise RuntimeError(f"run_chiral_susceptibility_c_api failed with code {ret}")

        f_vol = float((ns**3) * nt)
        f_scaled = float((ns**3) * (nt**3) * (temp_mev**2))

        return {
            "mean_unscaled": float(out_mean_unscaled.value),
            "error_unscaled": float(out_error_unscaled.value),
            "factor_vol": f_vol,
            "mean_vol_scaled": float(out_mean_vol_scaled.value),
            "error_vol_scaled": float(out_error_vol_scaled.value),
            "factor_scaled": f_scaled,
            "mean_scaled": float(out_mean_scaled.value),
            "error_scaled": float(out_error_scaled.value),
            "num_configs": int(out_num_cfgs.value),
            "ns": ns,
            "nt": nt,
            "temp": temp_mev,
        }

