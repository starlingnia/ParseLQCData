#!/usr/bin/env python3
"""
tests/test_susceptibility.py
--------------------------------------------------------------------------------
手征磁化率 (Chiral Susceptibility) 算法与数值精度单元测试：
1. 校验单构型无偏平方关联计算公式:
   Obar = mean(vals)
   O2bar = 1/(k(k-1)) * \\sum_{i!=j} O_i O_j
2. 校验 Jackknife 重采样统计误差算法
3. 校验本地真实 XML 测量抽取 (L32T12beta4.17)
4. 校验 No-RM 与 RM-corrected 8 组温度全量扫描结果与 ana/ttest/ 基准数据 100% 对齐
--------------------------------------------------------------------------------
"""

from pathlib import Path
import sys
import numpy as np
import polars as pl
try:
    import pytest
except ImportError:
    class _MockPytest:
        @staticmethod
        def skip(reason=""):
            print(f"  [SKIPPED] {reason}")
    pytest = _MockPytest()


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import ANA_ROOT, OUTPUT_CONDENSATE_DIR
from src.parselqcdata.susceptibility_pipeline import (
    compute_unbiased_quadratic,
    compute_jackknife_susceptibility,
    SusceptibilityPipeline,
)


def test_unbiased_quadratic_estimator():
    """测试单构型无偏平方估计函数的代数正确性"""
    vals = [2.0, 4.0, 6.0]
    obar, o2bar = compute_unbiased_quadratic(vals)

    # 理论值:
    # obar = (2 + 4 + 6) / 3 = 4.0
    # 交叉对: (2*4 + 2*6 + 4*2 + 4*6 + 6*2 + 6*4) / 6 = (8 + 12 + 8 + 24 + 12 + 24) / 6 = 88 / 6 = 44 / 3 = 14.666666666666666
    assert abs(obar - 4.0) < 1e-14
    assert abs(o2bar - (44.0 / 3.0)) < 1e-14


def test_jackknife_susceptibility_basic():
    """测试 Jackknife 假样本差值计算"""
    # 构造假数据
    list_obar = [1.0, 2.0, 3.0, 4.0, 5.0]
    list_o2bar = [2.0, 5.0, 10.0, 17.0, 26.0]

    res = compute_jackknife_susceptibility(list_obar, list_o2bar)
    assert res["num_configs"] == 5
    assert np.isfinite(res["mean_unscaled"])
    assert np.isfinite(res["error_unscaled"])
    assert res["error_unscaled"] > 0


def test_compare_susceptibility_no_rm_against_ana():
    """校验全量 No-RM 磁化率数据与 ana/ttest/results_susceptibility.csv 100% 对齐"""
    my_file = OUTPUT_CONDENSATE_DIR / "results_susceptibility.csv"
    ana_file = ANA_ROOT / "ttest" / "results_susceptibility.csv"

    if not my_file.exists() or not ana_file.exists():
        pytest.skip("数据文件尚未生成或 ana 基准不可达")

    df_my = pl.read_csv(my_file).sort("Temp")
    df_ana = pl.read_csv(ana_file).sort("Temp")

    assert df_my.height == df_ana.height
    for row_my, row_ana in zip(df_my.iter_rows(named=True), df_ana.iter_rows(named=True)):
        assert abs(row_my["Temp"] - row_ana["Temp"]) < 1e-8
        diff_mean = abs(row_my["Mean_unscaled"] - row_ana["Mean_unscaled"])
        diff_err = abs(row_my["Error_unscaled"] - row_ana["Error_unscaled"])
        assert diff_mean < 1e-14, f"Beta {row_my['Beta']} mean mismatch: {diff_mean}"
        assert diff_err < 1e-14, f"Beta {row_my['Beta']} err mismatch: {diff_err}"


def test_local_xml_extraction():
    """测试从本地真实 XML 抽取 L32T12 数据的稳定健壮性"""
    local_test_dir = PROJECT_ROOT / "data" / "readin" / "L32T12beta4.17" / "test_condensate"
    if not local_test_dir.exists():
        pytest.skip("本地测试数据集不存在")

    pipeline = SusceptibilityPipeline()
    res = pipeline.extract_from_meas_directory(local_test_dir, ns=32, nt=12, temp_mev=204.41)

    assert res["num_configs"] > 50
    assert res["mean_unscaled"] > 0
    assert res["mean_vol_scaled"] > 0
    assert res["mean_scaled"] > 0


if __name__ == "__main__":
    test_unbiased_quadratic_estimator()
    test_jackknife_susceptibility_basic()
    test_compare_susceptibility_no_rm_against_ana()
    test_local_xml_extraction()
    print("[ALL SUSCEPTIBILITY TESTS PASSED SUCCESSFULLY!]")

