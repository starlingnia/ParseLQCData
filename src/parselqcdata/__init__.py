"""
src/parselqcdata/__init__.py
--------------------------------------------------------------------------------
ParseLQCData Python 科学计算管道导出包
--------------------------------------------------------------------------------
"""

# 基础管道类
from .meson_pipeline import MesonPipeline
from .condensate_pipeline import CondensatePipeline
from .susceptibility_pipeline import SusceptibilityPipeline, SUSCEPTIBILITY_TEMP_MAP

# 独立原子计算脚本与模块 (按职责拆分)
from .xml_pbp_extractor import extract_pbp_from_xml, extract_pbp_from_directory
from .unbiased_quadratic import compute_unbiased_quadratic
from .scaling_factors import (
    compute_scaling_factors,
    get_zm_factor,
    parse_ensemble_meta_from_dir,
    is_valid_condensate_dir,
)
from .susceptibility_calculator import (
    jackknife_resample,
    compute_jackknife_susceptibility,
)
from .lcp_exporter import export_lcp_outputs
from .condensate_subtraction import compute_condensate_subtraction

# 有效质量与拟合工具
from .effective_mass_solver import (
    solve_effective_mass_one_point,
    compute_effective_mass_matrix,
    solve_effective_mass_one_point_centered,
    compute_effective_mass_matrix_centered,
    cosh_ratio,
    log_cosh,
)
from .plateau_detector import detect_plateau_window, select_best_window
from .cosh_fitter import (
    LEAST_SQUARES_FITTER,
    chi2_least_squares_fit,
    fit_meson_plateau,
    fit_single_jackknife_column,
    fit_single_jackknife_column_centered,
    fit_mass_window_scan,
    fit_jackknife_mass_centered,
)

__all__ = [
    # 核心管道
    "MesonPipeline",
    "CondensatePipeline",
    "SusceptibilityPipeline",
    "SUSCEPTIBILITY_TEMP_MAP",
    # 原子工具函数
    "extract_pbp_from_xml",
    "extract_pbp_from_directory",
    "compute_unbiased_quadratic",
    "compute_scaling_factors",
    "get_zm_factor",
    "parse_ensemble_meta_from_dir",
    "is_valid_condensate_dir",
    "jackknife_resample",
    "compute_jackknife_susceptibility",
    "export_lcp_outputs",
    "compute_condensate_subtraction",
    # 有效质量与平台拟合
    "solve_effective_mass_one_point",
    "compute_effective_mass_matrix",
    "detect_plateau_window",
    "solve_effective_mass_one_point_centered",
    "compute_effective_mass_matrix_centered",
    "cosh_ratio",
    "log_cosh",
    "select_best_window",
    "LEAST_SQUARES_FITTER",
    "chi2_least_squares_fit",
    "fit_meson_plateau",
    "fit_single_jackknife_column",
    "fit_single_jackknife_column_centered",
    "fit_mass_window_scan",
    "fit_jackknife_mass_centered",
]
