import ctypes
from pathlib import Path
from typing import List, Optional, Tuple
import numpy as np

class MesonOrchestrator:
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
            project_root / "lib" / "liblqcd_meson.so",
        ]
        for c in candidates:
            if c.exists():
                return c
        raise FileNotFoundError(
            "Meson shared library not found. Checked:\n"
            + "\n".join(f"  - {c}" for c in candidates)
        )
    
    def _setup_bindings(self) -> None:
        # int run_meson_pipeline_c_api(...)
        self._lib.run_meson_pipeline_c_api.argtypes = [
            ctypes.c_char_p,                                    # input_dir
            ctypes.c_char_p,                                    # channel_types_csv
            ctypes.c_char_p,                                    # channel_dirs_csv
            ctypes.c_int,                                       # binsize
            ctypes.c_int,                                       # num_lines
            ctypes.c_int,                                       # thread_count
            ctypes.c_int,                                       # is_single_source
            ctypes.POINTER(ctypes.c_double),                    # out_means
            ctypes.POINTER(ctypes.c_double),                    # out_errors
            ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),    # out_folded_jk
            ctypes.POINTER(ctypes.c_int),                       # out_n_bins
            ctypes.POINTER(ctypes.c_int),                       # out_n_raw_cfgs
        ]
        self._lib.run_meson_pipeline_c_api.restype = ctypes.c_int
    
        # void free_lqcd_buffer(void* ptr)
        self._lib.free_lqcd_buffer.argtypes = [ctypes.c_void_p]
        self._lib.free_lqcd_buffer.restype = None
    
    def process_channel(
        self,
        input_dir: str,
        channel_configs: List[dict],
        binsize: int = 4,
        num_lines: int = 48,
        thread_count: int = 0,
        is_single_source: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
        types = [c['type'] for c in channel_configs]
        dirs = [c['dir'] for c in channel_configs]
    
        types_csv = ",".join(types).encode('utf-8')
        dirs_csv = ",".join(dirs).encode('utf-8')
        input_dir_bytes = str(input_dir).encode('utf-8')
    
        out_means = (ctypes.c_double * num_lines)()
        out_errors = (ctypes.c_double * num_lines)()
        out_folded_jk_ptr = ctypes.POINTER(ctypes.c_double)()
        out_n_bins = ctypes.c_int(0)
        out_n_raw_cfgs = ctypes.c_int(0)
    
        ret = self._lib.run_meson_pipeline_c_api(
            input_dir_bytes,
            types_csv,
            dirs_csv,
            ctypes.c_int(binsize),
            ctypes.c_int(num_lines),
            ctypes.c_int(thread_count),
            ctypes.c_int(1 if is_single_source else 0),
            out_means,
            out_errors,
            ctypes.byref(out_folded_jk_ptr),
            ctypes.byref(out_n_bins),
            ctypes.byref(out_n_raw_cfgs)
        )
    
        if ret != 0:
            raise RuntimeError(f"run_meson_pipeline_c_api failed with code {ret}")
    
        means_arr = np.ctypeslib.as_array(out_means, shape=(num_lines,)).copy()
        errors_arr = np.ctypeslib.as_array(out_errors, shape=(num_lines,)).copy()
    
        n_bins = out_n_bins.value
        n_raw_cfgs = out_n_raw_cfgs.value
    
        if bool(out_folded_jk_ptr) and n_bins > 0:
            folded_jk_arr = np.ctypeslib.as_array(out_folded_jk_ptr, shape=(num_lines, n_bins)).copy()
            self._lib.free_lqcd_buffer(ctypes.cast(out_folded_jk_ptr, ctypes.c_void_p))
        else:
            folded_jk_arr = np.zeros((num_lines, 0), dtype=np.float64)
    
        return means_arr, errors_arr, folded_jk_arr, n_bins, n_raw_cfgs
    
    def process_channel_multi(
        self,
        input_dir: str,
        channel_configs: List[dict],
        binsize: int = 4,
        num_lines: int = 48,
        thread_count: int = 0
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
        """专门用于多源（16点源空间循环平移平均）信道抽取"""
        return self.process_channel(
            input_dir=input_dir,
            channel_configs=channel_configs,
            binsize=binsize,
            num_lines=num_lines,
            thread_count=thread_count,
            is_single_source=False
        )
    
    def process_channel_single(
        self,
        input_dir: str,
        channel_configs: List[dict],
        binsize: int = 4,
        num_lines: int = 48,
        thread_count: int = 0
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
        """专门用于单源（原点 0/0/0/0 单点源）信道抽取"""
        return self.process_channel(
            input_dir=input_dir,
            channel_configs=channel_configs,
            binsize=binsize,
            num_lines=num_lines,
            thread_count=thread_count,
            is_single_source=True
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
        """与早期分析调用保持兼容的接口"""
        configs = channel_mappings if channel_mappings is not None else (channel_configs or [])
        means, errors, folded_jk, _, _ = self.process_channel(
            input_dir=input_dir,
            channel_configs=configs,
            binsize=binsize,
            num_lines=num_lines,
            thread_count=thread_count,
            is_single_source=is_single_source
        )
        return means, errors, folded_jk