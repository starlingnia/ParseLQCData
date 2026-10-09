#!/usr/bin/env python3
"""
tests/test_pade_fitter.py
--------------------------------------------------------------------------------
手征磁化率 Padé 有理函数拟合与极点/峰值 Jackknife 统计推断单元测试套件
--------------------------------------------------------------------------------
"""

from pathlib import Path
import sys

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_LCP_DIR
from src.parselqcdata.pade_fitter import (
    PadeModel,
    fit_single_sample_pade,
    load_and_align_susceptibility_samples,
    fit_jackknife_susceptibility_pade,
)
import importlib

lcp_task_04 = importlib.import_module("scripts.tasks.lcp.04_fit_pade_susceptibility")
run_pade_fit_task = lcp_task_04.run_pade_fit_task


def test_pade_model_eval_and_poles():
    """验证 Padé [0/2] 与 [1/2] 模型的数值求值、极点求解与实轴峰值解析"""
    t0 = 157.0
    model02 = PadeModel(order="0/2", t0=t0)
    
    # 构造参数: p0=1.0, q1=0.0, q2=0.01 -> 峰值正好在 t=157.0, 分母 1 + 0.01*(t-157)^2
    p = {"p0": 1.0, "q1": 0.0, "q2": 0.01}
    assert abs(model02.eval(157.0, p) - 1.0) < 1e-12
    assert abs(model02.eval(167.0, p) - 0.5) < 1e-12  # 1 / (1 + 0.01*100) = 0.5
    
    # 极点检验: 1 + 0.01*dx^2 = 0 -> dx = +/- 10 i -> T_pole = 157.0 +/- 10.0 i
    pole_re, pole_im = model02.solve_pole(p)
    assert abs(pole_re - 157.0) < 1e-10
    assert abs(pole_im - 10.0) < 1e-10
    
    # 峰值检验: 在实轴上最大值在 t=157.0, 高度=1.0
    t_pk, chi_pk = model02.solve_peak(p)
    assert abs(t_pk - 157.0) < 1e-6
    assert abs(chi_pk - 1.0) < 1e-6


def test_align_susceptibility_samples_same_count():
    """确保 load_and_align_susceptibility_samples 能够使所有温度拥有严格一致的样本量"""
    csv_cand = OUTPUT_LCP_DIR / "susceptibility_jk_samples.csv"
    if not csv_cand.exists():
        print(f"  [SKIP] 样本文件 {csv_cand} 尚未生成")
        return

    # 1. 默认对齐 (N=2000)
    temps, mat, errs, betas = load_and_align_susceptibility_samples(csv_cand, max_samples=2000, binsize=1)
    assert len(temps) == 8
    assert mat.shape == (8, 2000), "所有 8 组温度必须截取相同且完整的 2000 组样本"
    assert len(errs) == 8
    assert np.all(errs > 0)
    
    # 2. binsize 聚合测试: N=2000, binsize=5 -> 400 bins
    temps_b, mat_b, errs_b, _ = load_and_align_susceptibility_samples(csv_cand, max_samples=2000, binsize=5)
    assert mat_b.shape == (8, 400), "binsize=5 时应聚合为 400 个 Jackknife 重抽样切片"


def test_fit_single_sample_with_lsqfit():
    """验证 fit_single_sample_pade 能正确调用 lsqfit 拟合合成数据并提取正确极点"""
    t0 = 157.0
    model = PadeModel(order="1/2", t0=t0)
    
    # 合成数据
    t_grid = np.array([138.0, 145.0, 153.0, 157.0, 164.0, 176.0, 202.0])
    true_p = {"p0": 0.95, "p1": -0.01, "q1": 0.002, "q2": 0.002}
    y_true = np.array([model.eval(t, true_p) for t in t_grid])
    y_noise = y_true + 1e-4 * np.ones_like(y_true)
    errs = np.full_like(t_grid, 0.01)
    
    res = fit_single_sample_pade(t_grid, y_noise, errs, model)
    assert res is not None
    assert res["ok"] is True
    assert res["chi2"] < 1.0
    assert 150.0 < res["t_peak"] < 160.0
    assert 150.0 < res["t_pole_real"] < 165.0
    assert res["t_pole_imag"] > 0


def test_full_pade_jackknife_workflow():
    """验证完整的 Jackknife 样本拟合与 CSV 产物落盘"""
    csv_cand = OUTPUT_LCP_DIR / "susceptibility_jk_samples.csv"
    if not csv_cand.exists():
        print(f"  [SKIP] 样本文件 {csv_cand} 尚未生成")
        return

    # 对前 50 个样本进行快速集成测试 (验证 Padé [3/2] 阶数)
    summary, records = run_pade_fit_task(
        order="3/2",
        t0=157.0,
        binsize=1,
        max_samples=50,
        csv_path=csv_cand,
        output_dir=OUTPUT_LCP_DIR,
    )

    assert summary["n_samples"] == 50
    assert summary["n_fitted"] == 50
    assert np.isfinite(summary["t_peak_mean"])
    assert summary["t_peak_err"] > 0, "Jackknife 峰值温度偏差必须大于 0"
    assert np.isfinite(summary["t_pole_real_mean"])
    assert summary["t_pole_real_err"] > 0, "Jackknife 极点实部偏差必须大于 0"
    assert summary["t_pole_imag_mean"] > 0, "极点虚部必须为正"

    # 验证输出 CSV 文件存在
    summary_file = OUTPUT_LCP_DIR / "pade_fit_summary.csv"
    samples_file = OUTPUT_LCP_DIR / "pade_jk_samples.csv"
    assert summary_file.exists(), f"汇总表 {summary_file} 必须存在"
    assert samples_file.exists(), f"明细表 {samples_file} 必须存在"

    df_samples = pl.read_csv(samples_file)
    assert df_samples.height >= 50
    assert "jk_index" in df_samples.columns
    assert "t_peak" in df_samples.columns
    assert "t_pole_real" in df_samples.columns
    assert "t_pole_imag" in df_samples.columns
    assert "chi_peak" in df_samples.columns
    assert "chi2_dof" in df_samples.columns


def main() -> None:
    print("[TEST] 运行 Padé 有理函数逼近与极点推断测试套件...")
    test_pade_model_eval_and_poles()
    print("  [PASS] test_pade_model_eval_and_poles")
    test_align_susceptibility_samples_same_count()
    print("  [PASS] test_align_susceptibility_samples_same_count")
    test_fit_single_sample_with_lsqfit()
    print("  [PASS] test_fit_single_sample_with_lsqfit")
    test_full_pade_jackknife_workflow()
    print("  [PASS] test_full_pade_jackknife_workflow")
    print("[ALL PADE FITTER TESTS PASSED SUCCESSFULLY!]")


if __name__ == "__main__":
    main()
