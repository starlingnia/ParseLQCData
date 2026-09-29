#!/usr/bin/env python3
"""
tests/test_meson_b417_scan.py
--------------------------------------------------------------------------------
b4.17 Nt 扫描计算任务 (关联函数 + 介子质量) 的设定与数值回归测试:
1. 任务表有序性: 6 组格点 x 4 组 ml 的顺序 / 目录名 / 温度与用户给定表格一致
2. 镜像窗口: Ns=48 时必须原样复现 4.17 @ 48^3x16 既有 slice 表
3. 对称点参数化: 对 48^3x16 的真实数据, 新接口必须与仓库既有发布结果逐位一致
   (meff_PS.csv 与 summary_fit_PS.csv)
4. 平台自动检测: 平坦序列通过, 含污染漂移的序列给出更靠后的窗口
5. 结果顺序: resolve_tasks / case 输出名严格按任务表排序
--------------------------------------------------------------------------------
"""

from pathlib import Path
import sys

import numpy as np
import polars as pl

try:
    import pytest
except ImportError:  # pragma: no cover
    class _MockPytest:
        @staticmethod
        def skip(reason=""):
            print(f"  [SKIPPED] {reason}")

    pytest = _MockPytest()


def close(a, b, rel: float = 1e-9, abs_tol: float = 0.0) -> bool:
    """不依赖 pytest 的近似比较 (仓库依赖里没有 pytest)"""
    return abs(float(a) - float(b)) <= max(abs_tol, rel * max(abs(float(a)), abs(float(b))))

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.meson_scan_setup import (  # noqa: E402
    BETA_KEY,
    CHANNEL_ORDER,
    REFERENCE_HALF,
    SCAN_CASES,
    SINGLE_FIT_SLICES,
    MULTI_FIT_SLICES,
    binsize_for,
    case_output_dir,
    discover_b417_dirs,
    readin_root,
    reference_window,
    resolve_tasks,
    scan_windows,
)
from docs.physics_setup import CHANNEL_CONFIGS, OUTPUT_ROOT  # noqa: E402
from src.parselqcdata import (  # noqa: E402
    compute_effective_mass_matrix_centered,
    detect_plateau_window,
    fit_jackknife_mass_centered,
    solve_effective_mass_one_point_centered,
)
from src.parselqcdata.effective_mass import solve_effective_mass_one_point  # noqa: E402
from tools.meson_orchestrator import MesonOrchestrator  # noqa: E402

REFERENCE_ENSEMBLE = readin_root() / "48x16b4.17" / "Output"
REFERENCE_MEFF = OUTPUT_ROOT / "ratio_results" / "b4.17" / "meff_PS.csv"
REFERENCE_FIT = OUTPUT_ROOT / "simulateresult" / "b4.17" / "summary_fit_PS.csv"


def test_scan_table_matches_spec():
    """任务表必须与用户给定表格逐行一致, 且顺序即结果顺序"""
    spec = [
        (32, 12, 204.00),
        (32, 14, 174.86),
        (32, 16, 153.00),
        (36, 18, 136.00),
        (40, 16, 153.00),
        (48, 18, 136.00),
    ]
    assert len(SCAN_CASES) == len(spec)
    assert [c.order for c in SCAN_CASES] == list(range(len(spec)))
    for case, (ns, nt, temp_ref) in zip(SCAN_CASES, spec):
        assert (case.ns, case.nt) == (ns, nt)
        assert close(case.temp_ref_mev, temp_ref, 0.0, 1e-6)
        assert close(case.ms, 0.040, 0.0, 1e-12)
        assert case.mls == (0.0020, 0.0035, 0.0070, 0.0120)
        assert case.half == ns // 2
        assert case.dirname(0.0020) == f"{ns}x{nt}_b4.17_ms0.040m0.0020"


def test_output_names_are_ordered():
    """输出目录名前缀携带任务序号, 文件系统排序 == 任务表排序"""
    names = []
    for case in SCAN_CASES:
        for ml in case.mls:
            names.append(case_output_dir(case, ml, False).name)
    assert names == sorted(names)
    assert names[0].startswith("00_32x12_")
    assert names[-1].startswith("05_48x18_")


def test_resolve_tasks_order_and_count():
    tasks = resolve_tasks()
    assert len(tasks) == 24
    keys = [(t["order"], t["ml_index"]) for t in tasks]
    assert keys == sorted(keys)
    assert all(t["input_dir"].endswith("Output") for t in tasks)
    # 每组格点 4 个 ml, 顺序与任务表一致
    for case in SCAN_CASES:
        sub = [t for t in tasks if t["case_key"] == case.key]
        assert [t["ml"] for t in sub] == list(case.mls)


def test_directory_discovery_counts():
    used, extras = discover_b417_dirs()
    assert len(used) == 24
    # 数字开头 + 含 b4.17 但不在任务表内的目录应被显式标记跳过
    assert "48x16b4.17" in extras
    assert not any(name in used for name in extras)


def test_reference_window_reproduces_published_slices():
    """Ns=48 (half=24) 时镜像窗口必须与既有 slice 表完全相同"""
    case48 = next(c for c in SCAN_CASES if c.ns == 48)
    for ch in CHANNEL_ORDER:
        for is_single, table in ((False, MULTI_FIT_SLICES), (True, SINGLE_FIT_SLICES)):
            start, end = reference_window(case48, ch, is_single)
            assert (start, end) == (table[BETA_KEY][ch], 25)
    assert REFERENCE_HALF == 24


def test_scan_windows_are_well_formed():
    for case in SCAN_CASES:
        windows = scan_windows(case)
        assert len(windows) >= 5
        for start, end in windows:
            assert 1 <= start < case.half
            assert end == case.half + 1
            assert end - start >= 5


def test_centered_solver_matches_legacy_at_half_24():
    """
    对称点取 24 时, 新求解器必须与旧接口给出同一根;
    x < half 时 cosh 比值 >= 1, 因此 y < 1 无实根 -> 新接口显式返回 NaN (旧接口返回伪解)。
    """
    for y, x in ((4.0, 1), (2.0, 8), (1.2, 20), (1.05, 22)):
        new = solve_effective_mass_one_point_centered(y, x, 24.0)
        old = solve_effective_mass_one_point(y, x)
        assert close(new, old, 1e-10, 1e-12), f"x={x} y={y}: {new} != {old}"

    # 无实根 -> NaN (旧接口此处返回 1.379e-05 之类的伪解)
    assert not np.isfinite(solve_effective_mass_one_point_centered(0.5, 14, 24.0))


def _synthetic_jackknife_meff(profile: np.ndarray, n_bins: int = 60, seed: int = 7) -> np.ndarray:
    """
    构造接近真实 Jackknife 结构的有效质量矩阵: 各 bin 共享同一条 "物理" 曲线,
    bin 间只有很小的乘性抖动 + 噪声 (真实数据里 meff 的 bin 间差分远小于点间漂移)。
    """
    rng = np.random.default_rng(seed)
    jitter = 1.0 + 0.002 * rng.standard_normal((1, n_bins))
    noise = 5e-4 * rng.standard_normal((profile.size, n_bins))
    return 0.25 + profile[:, None] * jitter + noise


def test_plateau_detection_on_synthetic_data():
    """平坦有效质量序列 -> 检出长窗口; 含未消失污染漂移 -> 窗口起点后移"""
    half = 16
    x = np.arange(half, dtype=float)

    flat = _synthetic_jackknife_meff(np.zeros_like(x))
    det_flat = detect_plateau_window(flat, end=half, chi2_dof_max=1.5, min_points=4)
    assert det_flat["ok"], "平坦序列必须检出平台"
    assert det_flat["start"] <= 3
    assert close(det_flat["mass_est"], 0.25, 0.0, 0.01)

    drift_profile = np.where(x < 8, 0.03 * (8 - x), 0.0)
    drifted = _synthetic_jackknife_meff(drift_profile)
    det_drift = detect_plateau_window(drifted, end=half, chi2_dof_max=1.5, min_points=4)
    assert det_drift["ok"]
    assert det_drift["start"] >= 8, f"污染段必须被排除, 实际 start={det_drift['start']}"
    assert det_drift["start"] > det_flat["start"]


def test_binsize_mapping():
    for ch in CHANNEL_ORDER:
        assert binsize_for(ch, False) == 4
        assert binsize_for(ch, True) == 4


def test_against_published_b417_results():
    """真实数据回归: 48^3x16 的 meff 与 PS 平台质量必须复现既有发布结果"""
    if not REFERENCE_ENSEMBLE.is_dir() or not REFERENCE_MEFF.exists():
        pytest.skip(f"缺少 48^3x16 参考数据: {REFERENCE_ENSEMBLE}")

    orch = MesonOrchestrator()
    means, errors, jk, n_bins, _ = orch.process_channel(
        input_dir=str(REFERENCE_ENSEMBLE),
        channel_configs=CHANNEL_CONFIGS["PS"],
        binsize=4,
        num_lines=48,
        thread_count=0,
        is_single_source=False,
    )
    assert jk.shape[0] == 48

    _, meff_mean, meff_err = compute_effective_mass_matrix_centered(jk, 24.0)
    ref = pl.read_csv(REFERENCE_MEFF).to_numpy()
    assert np.nanmax(np.abs(meff_mean - ref[:, 0])) < 1e-10
    assert np.nanmax(np.abs(meff_err - ref[:, 1])) < 1e-10

    fit = fit_jackknife_mass_centered(np.arange(48, dtype=float), jk, errors, 8, 25, 24.0)
    ref_summary = pl.read_csv(REFERENCE_FIT)
    ref_mass = float(ref_summary["mean"][0])
    ref_err = float(ref_summary["jack_err"][0])
    # 与既有发布值一致到优化器容差 (同一 cosh 模型, 仅初值策略不同)
    assert close(fit["mass"], ref_mass, 1e-6), f"{fit['mass']} != {ref_mass}"
    assert close(fit["mass_err"], ref_err, 1e-5), f"{fit['mass_err']} != {ref_err}"


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"  [PASS] {name}")
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"  [FAIL] {name}: {type(exc).__name__}: {exc}")
    print(f"[TEST] b4.17 Nt 扫描测试完成, 失败 {failures} 项")
    raise SystemExit(1 if failures else 0)
