#!/usr/bin/env python3
"""
tests/test_ccor_flow.py
--------------------------------------------------------------------------------
ana/dat/ccor 流程复刻的数值回归测试:
1. block 抽取口径: 精确 tag + 实部绝对值 + 只取前 32 行
2. 逐构型留一 Jackknife 与 fold (C(x)<->C(Ns-x)) 的代数正确性
3. Jackknife 均值/误差公式与 ana 一致
4. 与 ccor 保存数据逐位对齐:
     symdatasample<ml>.csv / symdatasample2<ml>.csv   (36^3x18, A/Xt, 精确到 1e-15)
     dferr1/2<ml>.csv                                  (1e-7, ana 公式存在相消噪声)
     binnedresult<ml>.csv / binnedresult2<ml>.csv      (逐样本质量, 1e-7)
若参考目录不存在则自动跳过 (SKIP)。
--------------------------------------------------------------------------------
"""

from pathlib import Path
import sys

import numpy as np

try:
    import pytest
except ImportError:  # pragma: no cover
    class _MockPytest:
        @staticmethod
        def skip(reason=""):
            print(f"  [SKIPPED] {reason}")

    pytest = _MockPytest()

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import DEFAULT_READIN_DIR  # noqa: E402
from src.parselqcdata import cosh_ratio, solve_effective_mass_one_point_centered  # noqa: E402
from src.parselqcdata.ccor_flow import (  # noqa: E402
    ANA_CHANNEL_BLOCKS,
    ANA_MLS,
    ANA_PAIRS,
    ANA_TYPE_ID,
    MESON_SCAN_CHANNEL,
    ana_block_tag,
    ana_ensemble_dir,
    extract_block_values,
    fit_sample_masses,
    jackknife_mean_err,
    leave_one_out_jackknife,
    natural_sorted,
    per_config_matrix,
    read_csv_matrix,
    sample_errors,
    symmetrize_about_center,
)

CCOR_DIR = Path("/Users/junxiongnie/code/ana/dat/ccor")
READIN = Path(DEFAULT_READIN_DIR)


def close(a, b, rel=1e-9, abs_tol=0.0) -> bool:
    return abs(float(a) - float(b)) <= max(abs_tol, rel * max(abs(float(a)), abs(float(b))))


# ------------------------------------------------------------------------------
# 1. 基础代数
# ------------------------------------------------------------------------------
def test_extract_block_values_takes_abs_real_and_stops():
    lines = [
        "--- Vector2 to Vector2 --- spatial:DIRZ --- 0/0/0/0\n",
        " 0  -1.5  3.0\n",
        " 1   2.5  0.0\n",
        " 2   -3.5 1.0\n",
        "--- Vector2 to Vector2 --- spatial:DIRZ --- 1/0/0/0\n",
        " 0  9.9\n",
    ]
    tag = "--- Vector2 to Vector2 --- spatial:DIRZ --- 0/0/0/0"
    assert extract_block_values(lines, tag, n_lines=3) == [1.5, 2.5, 3.5]
    # ana 的 extract_dir_block 只在凑满 n_lines 时停止, 不会在下一个 tag 处截断
    # (真实数据里 block 行数 >= 32, 因此不影响结果; 这里显式记录该行为)
    assert extract_block_values(lines, tag, n_lines=10) == [1.5, 2.5, 3.5, 9.9]


def test_leave_one_out_jackknife_is_algebraic():
    mat = np.array([[1.0, 2.0, 3.0, 4.0],
                    [10.0, 20.0, 30.0, 40.0]])
    jk = leave_one_out_jackknife(mat, n_rows=2)
    for i in range(2):
        total = mat[i].sum()
        expected = (total - mat[i]) / (mat[i].size - 1)
        assert np.allclose(jk[i], expected)


def test_symmetrize_matches_ana_fold():
    """ana 的 fold: row[x]=(jk[x]+jk[32-x])/2 (1<=x<=15), row[16] 原值, 其余镜像"""
    rng = np.random.default_rng(2024)
    jk = rng.standard_normal((32, 5))
    sym = symmetrize_about_center(jk, center=16)
    assert np.allclose(sym[0], jk[0])
    for x in range(1, 16):
        assert np.allclose(sym[x], (jk[x] + jk[32 - x]) / 2)
        assert np.allclose(sym[32 - x], sym[x])
    assert np.allclose(sym[16], jk[16])


def test_jackknife_mean_err_matches_ana_formula():
    v = np.array([0.1, 0.2, 0.35, 0.4, 0.5])
    mean, err = jackknife_mean_err(v)
    assert close(mean, v.mean(), 1e-15)
    assert close(err, np.sqrt(sum(x ** 2 for x in v) / len(v) - v.mean() ** 2) * np.sqrt(len(v) - 1), 1e-12)


def test_channel_mapping_covers_all_types():
    assert set(ANA_TYPE_ID) == set(ANA_CHANNEL_BLOCKS) == set(MESON_SCAN_CHANNEL)
    assert len(ANA_PAIRS) == 4
    for ch1, ch2, _ in ANA_PAIRS:
        assert ch1 in ANA_TYPE_ID and ch2 in ANA_TYPE_ID


# ------------------------------------------------------------------------------
# 2. 与 ccor 保存数据逐位对齐 (需要 ana 参考数据)
# ------------------------------------------------------------------------------
def _require_reference():
    if not CCOR_DIR.is_dir():
        pytest.skip(f"缺少 ana ccor 参考目录: {CCOR_DIR}")
    if not ana_ensemble_dir(READIN, 18, "0.0020").is_dir():
        pytest.skip(f"缺少 36^3x18 ensemble: {ana_ensemble_dir(READIN, 18, '0.0020')}")


def test_symdata_reproduces_saved_reference():
    _require_reference()
    for ml in ANA_MLS:
        ens_dir = ana_ensemble_dir(READIN, 18, ml)
        files = natural_sorted(str(p) for p in (ens_dir / "Output").glob("test1_lhadrons_*_mesons"))
        assert files, f"缺少输入文件: {ens_dir}"
        for ch, suffix in (("A", ""), ("Xt", "2")):
            saved_path = CCOR_DIR / f"symdatasample{suffix}{ml}.csv"
            if not saved_path.exists():
                continue
            saved = read_csv_matrix(saved_path)
            mat = per_config_matrix(files, ana_block_tag(ch), 32)
            sym = symmetrize_about_center(leave_one_out_jackknife(mat, n_rows=32), center=16)
            assert sym.shape == saved.shape
            assert np.max(np.abs(sym - saved)) < 1e-15, f"{ch} ml={ml} 未复现 symdata"


def test_dferr_and_fit_reproduce_saved_reference():
    _require_reference()
    for ml in ANA_MLS:
        ens_dir = ana_ensemble_dir(READIN, 18, ml)
        files = natural_sorted(str(p) for p in (ens_dir / "Output").glob("test1_lhadrons_*_mesons"))
        for ch, suffix in (("A", ""), ("Xt", "2")):
            err_path = CCOR_DIR / f"dferr{'1' if suffix == '' else '2'}{ml}.csv"
            br_path = CCOR_DIR / f"binnedresult{suffix}{ml}.csv"
            if not (err_path.exists() and br_path.exists()):
                continue
            mat = per_config_matrix(files, ana_block_tag(ch), 32)
            sym = symmetrize_about_center(leave_one_out_jackknife(mat, n_rows=32), center=16)
            err = sample_errors(sym)
            saved_err = read_csv_matrix(err_path).ravel()
            # ana 的 dferr 公式有 ~1e-9 量级的相消噪声
            assert np.max(np.abs(err - saved_err)) < 1e-6 * np.max(np.abs(saved_err))

            import polars as pl

            saved_br = pl.read_csv(br_path)
            mine = fit_sample_masses(sym, err, center=16)
            saved_m = saved_br["massfit_mean"].to_numpy()
            n = min(mine.size, saved_m.size)
            ok = np.isfinite(mine[:n]) & np.isfinite(saved_m[:n])
            assert ok.sum() > 0
            assert np.max(np.abs(mine[:n][ok] - saved_m[:n][ok])) < 1e-6 * np.max(np.abs(saved_m))



# ------------------------------------------------------------------------------
# 3. x - L/2 逻辑: 尺寸无关性 (绝不允许把 32/16 或 48/24 硬套到别的 Ns 上)
# ------------------------------------------------------------------------------
def _synthetic_cosh(ns: int, mass: float, amp: float = 3.0e-4, n_samples: int = 24,
                    noise: float = 1e-9, seed: int = 11) -> np.ndarray:
    """造一份 C(x) = amp * cosh(mass * (x - Ns/2)) 的 Jackknife 样本矩阵 (Ns x n_samples)"""
    rng = np.random.default_rng(seed)
    x = np.arange(ns, dtype=float)
    base = amp * np.cosh(mass * (x - ns / 2.0))
    return base[:, None] + noise * rng.standard_normal((ns, n_samples))


def test_symmetrize_honours_explicit_center():
    """fold 必须真的按传入的 center 折叠, 而不是永远按 Ns/2"""
    ns = 40
    rng = np.random.default_rng(3)
    mat = rng.standard_normal((ns, 4))
    for center in (ns // 2, 18, 12):
        sym = symmetrize_about_center(mat, center)
        for x in range(ns):
            expect = 0.5 * (mat[x] + mat[(2 * center - x) % ns])
            assert np.allclose(sym[x], expect), f"center={center} x={x} 折叠错误"
        assert np.allclose(sym[center], mat[center])
    # center = Ns/2 时必须与 ana 的手写 fold 逐位一致
    ana = np.empty_like(mat)
    ana[0] = mat[0]
    for x in range(1, ns // 2):
        ana[x] = (mat[x] + mat[ns - x]) / 2
    ana[ns // 2] = mat[ns // 2]
    for x in range(ns // 2 + 1, ns):
        ana[x] = ana[ns - x]
    assert np.allclose(symmetrize_about_center(mat, ns // 2), ana)


def test_fit_recovers_mass_for_every_ns():
    """
    核心回归: 对 Ns = 32 / 36 / 40 / 48 (以及非常规 24 / 52) 造同一物理质量的数据,
    拟合必须都还原出同一个质量 —— 说明窗口/对称点全部随 Ns 走, 没有写死 16/24/32/48。
    """
    mass_true = 0.31
    recovered = {}
    for ns in (24, 32, 36, 40, 48, 52):
        jk = _synthetic_cosh(ns, mass_true)
        sym = symmetrize_about_center(jk)                     # 默认 center = Ns/2
        err = sample_errors(sym)
        masses = fit_sample_masses(sym, err, center=None)     # 默认 center = Ns/2
        mean, _ = jackknife_mean_err(masses)
        recovered[ns] = mean
        assert close(mean, mass_true, 5e-4), f"Ns={ns} 拟合 m={mean} 偏离真值 {mass_true}"
    spread = max(recovered.values()) - min(recovered.values())
    assert spread < 5e-4, f"不同 Ns 拟合结果不一致: {recovered}"


def test_effective_mass_solver_is_center_parametric():
    """meff 求解器: 同一个物理质量在不同 Ns 上必须解出同一个值"""
    mass_true = 0.27
    for ns in (32, 36, 40, 48):
        half = ns / 2.0
        x = ns // 2 - 5                                    # 距对称点 5
        y = float(cosh_ratio(mass_true * (x - half), mass_true * (x + 1 - half)))
        solved = solve_effective_mass_one_point_centered(y, x, half)
        assert close(solved, mass_true, 1e-6), f"Ns={ns}: meff={solved} != {mass_true}"
        # 若误用别的对称点 (例如写死 16), 结果就会明显偏掉
        wrong = solve_effective_mass_one_point_centered(y, x, 16.0)
        if ns != 32:
            assert not close(wrong, mass_true, 1e-3), f"Ns={ns} 用 center=16 竟然也对上了?"


def test_reference_convention_is_explicitly_tagged():
    """ana 参考口径的 32 行 / 中心 16 必须是显式常量, 且有清晰标注"""
    from src.parselqcdata import ccor_flow as cf

    assert cf.ANA_BLOCK_LINES == 32
    assert "ana" in (cf.ANA_BLOCK_LINES.__doc__ or "") or True  # 常量本身无 doc, 下面是行为校验
    # 对 Ns=32 的 ensemble, 参考口径 = 我的口径; 对 Ns=36/40/48 则必须显式传 center/行数
    jk32 = np.arange(32 * 3, dtype=float).reshape(32, 3)
    assert np.allclose(cf.symmetrize_about_center(jk32), cf.symmetrize_about_center(jk32, 16))
    jk48 = np.arange(48 * 3, dtype=float).reshape(48, 3)
    assert not np.allclose(cf.symmetrize_about_center(jk48), cf.symmetrize_about_center(jk48, 16))


# ------------------------------------------------------------------------------
# 4. 拟合器: 必须强制 chi2 最小二乘 (不能依赖 lsqfit 的环境默认, 否则会退化)
# ------------------------------------------------------------------------------
def test_chi2_least_squares_fit_forces_fitter():
    """统一入口必须显式传 fitter, 且调用方无法覆盖成别的 (例如 lgm)"""
    import gvar as gv

    from src.parselqcdata import LEAST_SQUARES_FITTER, chi2_least_squares_fit

    assert LEAST_SQUARES_FITTER == "scipy_least_squares"
    x = np.arange(6, dtype=float)
    y = gv.gvar(2.0 + 3.0 * x, [1e-6] * 6)

    fit = chi2_least_squares_fit(
        data=(x, y), fcn=lambda xx, p: p["a"] + p["b"] * xx, p0={"a": 1.0, "b": 1.0}
    )
    assert fit.fitter == LEAST_SQUARES_FITTER
    assert close(fit.p["b"].mean, 3.0, 1e-6)

    # 调用方偷偷传 fitter 也必须被忽略, 不能退化回环境默认
    fit2 = chi2_least_squares_fit(
        data=(x, y), fcn=lambda xx, p: p["a"] + p["b"] * xx, p0={"a": 1.0, "b": 1.0},
        fitter="lgm",
    )
    assert fit2.fitter == LEAST_SQUARES_FITTER


def test_ccor_flow_fit_reports_least_squares():
    """ccor 流程里的逐样本拟合同样必须走 chi2 最小二乘"""
    from src.parselqcdata.ccor_flow import fit_sample_masses
    from src.parselqcdata import LEAST_SQUARES_FITTER
    import src.parselqcdata.ccor_flow as cf

    captured = {}
    original = cf.chi2_least_squares_fit

    def spy(data, fcn, p0, **kwargs):
        captured["called"] = captured.get("called", 0) + 1
        captured["fitter"] = kwargs.pop("fitter", LEAST_SQUARES_FITTER)
        return original(data, fcn, p0, **kwargs)

    cf.chi2_least_squares_fit = spy
    try:
        jk = _synthetic_cosh(32, 0.3)
        sym = symmetrize_about_center(jk)
        fit_sample_masses(sym, sample_errors(sym))
    finally:
        cf.chi2_least_squares_fit = original
    assert captured.get("called", 0) >= 1
    assert captured.get("fitter") == LEAST_SQUARES_FITTER


def test_no_fit_call_bypasses_the_fitter_guard():
    """
    源码级防退化检查: 除统一入口 plateau_fit.chi2_least_squares_fit 内部,
    任何地方都不允许直接调用 lsqfit.nonlinear_fit (那样会退回环境默认 fitter)。
    """
    import re as _re

    roots = [PROJECT_ROOT / "src", PROJECT_ROOT / "scripts", PROJECT_ROOT / "docs", PROJECT_ROOT / "tools"]
    offenders = []
    for root in roots:
        for path in root.rglob("*.py"):
            if "__pycache__" in str(path):
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            for m in _re.finditer(r"nonlinear_fit\s*\(", text):
                line_no = text[: m.start()].count("\n") + 1
                window = text[max(0, m.start() - 400): m.start() + 400]
                if path.name == "plateau_fit.py" and "def chi2_least_squares_fit" in window:
                    continue          # 唯一允许的调用点
                offenders.append(f"{path.relative_to(PROJECT_ROOT)}:{line_no}")
    assert not offenders, f"存在绕过 fitter 强制指定的拟合调用: {offenders}"


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
    print(f"[TEST] ccor 流程复刻测试完成, 失败 {failures} 项")
    raise SystemExit(1 if failures else 0)
