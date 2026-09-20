import ctypes
import os
from pathlib import Path
from typing import List, Optional, Tuple
import numpy as np

class LQCDOrchestrator:
    """
    Python 高层任务编排器：负责参数组织与驱动底层 C++ 多线程并发数据分析引擎
    """
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
            project_root / "bin" / "libparselqcdata.dylib",
            project_root / "bin" / "libparselqcdata.so",
            project_root / "build" / "libparselqcdata.dylib",
            project_root / "build" / "libparselqcdata.so",
        ]
        for c in candidates:
            if c.exists():
                return c
        raise FileNotFoundError(f"未找到共享库文件！尝试查找路径:\n" + "\n".join(f"  - {c}" for c in candidates))

    def _setup_bindings(self) -> None:
        # int run_meson_pipeline_c_api(...)
        self._lib.run_meson_pipeline_c_api.argtypes = [
            ctypes.c_char_p,                                    # input_dir
            ctypes.c_char_p,                                    # channel_types_csv
            ctypes.c_char_p,                                    # channel_dirs_csv
            ctypes.c_int,                                       # binsize
            ctypes.c_int,                                       # num_lines
            ctypes.c_int,                                       # thread_count
            ctypes.POINTER(ctypes.c_double),                    # out_means
            ctypes.POINTER(ctypes.c_double),                    # out_errors
            ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),    # out_folded_jk
            ctypes.POINTER(ctypes.c_int),                       # out_n_bins
            ctypes.POINTER(ctypes.c_int),                       # out_n_raw_cfgs
        ]
        self._lib.run_meson_pipeline_c_api.restype = ctypes.c_int

        self._lib.free_lqcd_buffer.argtypes = [ctypes.c_void_p]
        self._lib.free_lqcd_buffer.restype = None

    def execute_meson_analysis(
        self,
        input_dir: str,
        channel_mappings: List[dict],
        binsize: int = 4,
        num_lines: int = 48,
        thread_count: int = 0,
        return_folded_jk: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """
        驱动 C++ 引擎执行端到端 Meson 数据分析与统计计算
        返回: (means, errors, [folded_jk_matrix])
        """
        types = [m['type'] for m in channel_mappings]
        dirs = [m['dir'] for m in channel_mappings]
        types_str = ",".join(types).encode('utf-8')
        dirs_str = ",".join(dirs).encode('utf-8')
        input_dir_bytes = str(input_dir).encode('utf-8')

        out_means = (ctypes.c_double * num_lines)()
        out_errors = (ctypes.c_double * num_lines)()
        out_n_bins = ctypes.c_int(0)
        out_n_raw_cfgs = ctypes.c_int(0)

        c_folded_ptr = ctypes.POINTER(ctypes.c_double)() if return_folded_jk else None
        p_c_folded_ptr = ctypes.byref(c_folded_ptr) if return_folded_jk else None

        ret = self._lib.run_meson_pipeline_c_api(
            input_dir_bytes,
            types_str,
            dirs_str,
            ctypes.c_int(binsize),
            ctypes.c_int(num_lines),
            ctypes.c_int(thread_count),
            out_means,
            out_errors,
            p_c_folded_ptr,
            ctypes.byref(out_n_bins),
            ctypes.byref(out_n_raw_cfgs)
        )

        if ret != 0:
            raise RuntimeError(f"C++ 底层计算失败，错误码: {ret}")

        means = np.array([out_means[i] for i in range(num_lines)], dtype=np.float64)
        errors = np.array([out_errors[i] for i in range(num_lines)], dtype=np.float64)

        folded_jk = None
        if return_folded_jk and c_folded_ptr:
            n_bins = out_n_bins.value
            total_elements = num_lines * n_bins
            # 零拷贝读取连续内存
            raw_buf = np.ctypeslib.as_array(c_folded_ptr, shape=(num_lines, n_bins))
            folded_jk = raw_buf.copy()
            self._lib.free_lqcd_buffer(c_folded_ptr)

        return means, errors, folded_jk
