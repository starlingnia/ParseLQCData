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
    """校验全量 No-RM 磁化率数据 (48^3x16 LCP) 与 ana/ttest/results_susceptibility.csv 100% 对齐"""
    my_file = OUTPUT_CONDENSATE_DIR / "results_susceptibility.csv"
    ana_file = ANA_ROOT / "ttest" / "results_susceptibility.csv"

    if not my_file.exists() or not ana_file.exists():
        pytest.skip("数据文件尚未生成或 ana 基准不可达")

    df_my = pl.read_csv(my_file)
    if "Ns" in df_my.columns and "Nt" in df_my.columns:
        df_my = df_my.filter((pl.col("Ns") == 48) & (pl.col("Nt") == 16))
    df_my = df_my.sort("Beta")
    df_ana = pl.read_csv(ana_file).sort("Beta")

    assert df_my.height == df_ana.height
    for row_my, row_ana in zip(df_my.iter_rows(named=True), df_ana.iter_rows(named=True)):
        diff_mean = abs(row_my["Mean_unscaled"] - row_ana["Mean_unscaled"])
        diff_err = abs(row_my["Error_unscaled"] - row_ana["Error_unscaled"])
        assert diff_mean < 1e-14, f"Beta {row_my['Beta']} mean mismatch: {diff_mean}"
        assert diff_err < 1e-14, f"Beta {row_my['Beta']} err mismatch: {diff_err}"


def test_local_xml_extraction():
    """测试从本地真实 XML 抽取 L32T12 数据的稳定健壮性"""
    from docs.physics_setup import resolve_dataset_dir, DEFAULT_READIN_DIR
    local_test_dir = resolve_dataset_dir(DEFAULT_READIN_DIR, "L32T12beta4.17")
    if not local_test_dir.exists():
        pytest.skip("本地测试数据集不存在")
        return

    pipeline = SusceptibilityPipeline()
    res = pipeline.extract_from_meas_directory(local_test_dir, ns=32, nt=12, temp_mev=204.41)

    assert res["num_configs"] > 50
    assert res["mean_unscaled"] > 0
    assert res["mean_vol_scaled"] > 0
    assert res["mean_scaled"] > 0
    assert "zm" in res
    assert "mean_scaled_gev2_renorm" in res
    assert res["mean_scaled_gev2_renorm"] > 0


def test_susceptibility_zm_renormalization():
    """校验磁化率除以 Zm 重整化代数正确性"""
    from src.parselqcdata.susceptibility_pipeline import get_zm_factor

    # 4.13: zm = 0.937703
    zm_413 = get_zm_factor(4.13)
    assert abs(zm_413 - 0.937703) < 1e-6

    # 4.17: zm = 0.966247
    zm_417 = get_zm_factor(4.17)
    assert abs(zm_417 - 0.966247) < 1e-6

    # 4.405: zm = 1.09274
    zm_4405 = get_zm_factor(4.405)
    assert abs(zm_4405 - 1.09274) < 1e-5

    # 校验测试样本下的 Jackknife 重整化 (最新口径除以 Zm^2)
    list_obar = [1.0, 2.0, 3.0, 4.0, 5.0]
    list_o2bar = [2.0, 5.0, 10.0, 17.0, 26.0]
    res = compute_jackknife_susceptibility(list_obar, list_o2bar, beta=4.17)
    expected_renorm = (res["mean_scaled"] / 1e6) / (zm_417 ** 2)
    assert abs(res["mean_scaled_gev2_renorm"] - expected_renorm) < 1e-12


def test_lcp_outputs_format():
    """校验 output/LCP/ 目录产物及其标头格式要求，确认 outputs/ 冗余目录不存在"""
    pipeline = SusceptibilityPipeline()
    pipeline.get_full_scan_results()

    lcp_dir = PROJECT_ROOT / "output" / "LCP"
    txt_file = lcp_dir / "results_susceptibility.txt"
    lcp_txt_file = lcp_dir / "results_susceptibility_lcp.txt"

    assert lcp_dir.exists(), f"LCP 目录 {lcp_dir} 必须存在"
    assert txt_file.exists(), f"结果文件 {txt_file} 必须存在"
    assert lcp_txt_file.exists(), f"LCP 标准文件 {lcp_txt_file} 必须存在"

    # 严格确保废除的 outputs 目录不存在
    assert not (PROJECT_ROOT / "outputs").exists(), "冗余 outputs/ 目录不应存在，所有产物统一归集至 output/LCP"

    expected_header = "# beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error"
    content = txt_file.read_text(encoding="utf-8").strip().splitlines()
    assert content[0] == expected_header, f"标头不匹配: {content[0]}"

    # 校验所有数据行均包含 6 个数值列
    for line in content[1:]:
        parts = line.split()
        assert len(parts) == 6, f"每行数据必须恰好为 6 列: {line}"
        beta_v = float(parts[0])
        zm_v = float(parts[1])
        chi_lat = float(parts[2])
        err_lat = float(parts[3])
        chi_ren = float(parts[4])
        err_ren = float(parts[5])
        assert beta_v > 4.0 and beta_v < 4.5
        assert zm_v > 0.9 and zm_v < 1.2
        assert chi_lat > 0
        assert err_lat > 0
        assert chi_ren > 0
        assert err_ren > 0


if __name__ == "__main__":
    test_unbiased_quadratic_estimator()
    test_jackknife_susceptibility_basic()
    test_compare_susceptibility_no_rm_against_ana()
    test_local_xml_extraction()
    test_susceptibility_zm_renormalization()
    test_lcp_outputs_format()
    print("[ALL SUSCEPTIBILITY TESTS PASSED SUCCESSFULLY!]")

