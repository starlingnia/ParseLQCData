"""
src/parselqcdata/ccor_flow.py
--------------------------------------------------------------------------------
ana/dat/ccor 流程的严格复刻 (Python 侧)

ana 仓库里 `dat/ccor/` 保存的是一套 "两信道质量差 -> 对称性破缺信号" 的分析:
    1. 32x10_Beta4.70_m0.01_mesons.py:
       从 <readin>/32x12_b4.17_ms0.040m0.0020/Output/test1_lhadrons_*_mesons
       中按精确 tag 取 block, **取实部绝对值**, 只取前 32 行, 逐构型留一 Jackknife,
       再关于对称点 fold (C(x) <-> C(32-x));
    2. simulate.py: 对每个 Jackknife 样本做 cosh 平台拟合 a*cosh(m*(x-16)),
       窗口 x in [9, 24), 逐样本质量 m, 丢掉 chi2/dof > 100 的样本;
    3. dfiltered.py: 两个信道的逐样本质量作差, 用 Jackknife 统计给出质量差与误差;
    4. pplot.py / newdata.py / *.gp: 质量差 vs ml, 质量 vs 温度 (MeV = 质量 * 2640)。

本模块把这套口径参数化 (对称点 = Ns/2), 用于:
  * 复现 ccor 里保存的数据 (数值校验, 见 tests/test_ccor_flow.py);
  * 用 ParseLQCData 自己算出来的关联函数跑同一套流程 (同流程出图);
  * 两者对照。

注意: ana 脚本里 `for mapping in data_blocks_to_list1:` 只保留了列表最后一个 mapping,
      因此每个信道实际用的是固定 block (见 ANA_CHANNEL_BLOCKS)。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import glob
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gvar as gv
import numpy as np

from src.parselqcdata.plateau_fit import LEAST_SQUARES_FITTER, chi2_least_squares_fit

# ------------------------------------------------------------------------------
# 1. ana ccor 常量 (与 ana/dat/ccor/*.py 一一对应)
# ------------------------------------------------------------------------------
#: 6 个狄拉克双线性信道 -> ana 实际使用的 (vector_type, direction) block
ANA_CHANNEL_BLOCKS: Dict[str, Tuple[str, str]] = {
    "V": ("Vector2", "DIRZ"),      # trial.sh 第 1 块 list1 的最后一个 mapping
    "A": ("AVector2", "DIRZ"),     # 第 1 块 list2 / 第 4 块 list1
    "Tt": ("TVector2", "DIRY"),    # 第 2 块 list1
    "Xt": ("TAVector2", "DIRY"),   # 第 2 块 list2 / 第 4 块 list2
    "S": ("TVector4", "DIRY"),     # 第 3 块 list1
    "PS": ("TAVector4", "DIRY"),   # 第 3 块 list2
}

#: 信道 -> newdata.py 里的 type 编号
ANA_TYPE_ID: Dict[str, int] = {"V": 1, "A": 2, "Tt": 3, "Xt": 4, "S": 5, "PS": 6}

#: ana 信道名 -> ParseLQCData meson_scan 的信道名
MESON_SCAN_CHANNEL: Dict[str, str] = {
    "V": "Vec", "A": "AV", "Tt": "Tt", "Xt": "Xt", "S": "S", "PS": "PS",
}

#: 任务表里的 4 组 ml (字符串形式, 与 ana 文件名一致)
ANA_MLS: Tuple[str, ...] = ("0.0020", "0.0035", "0.0070", "0.0120")

#: 4 组对称性破缺信道对: (信道1, 信道2, 图标题/文件名)
ANA_PAIRS: Tuple[Tuple[str, str, str], ...] = (
    ("V", "A", "SU(2)XSU(2) V-A"),
    ("Tt", "Xt", "U(1)A T-X"),
    ("S", "PS", "U(1) S-PS"),
    ("A", "Xt", "SU(2)spinxchiral X-A"),
)

#: pplot.py / newdata.py 用的能量标度 a^-1 = 2640 MeV (T = 2640 / Nt)
ANA_SCALE_MEV: float = 2640.0

#: simulate.py 的拟合窗口与选样阈值 (窗口中心 = 对称点, 半宽 = 7 -> [c-7, c+8))
#: 说明: N=7 是**复刻旧 48^3x16 (ana) 口径**用的固定值, reference 数据源仍然照用它;
#:       本仓库自己的 multisrc/singlesrc 数据源改用 half_window_for(Ns) (见下)。
ANA_HALF_WINDOW: int = 7
ANA_CHI2_DOF_MAX: float = 100.0

#: 逐样本质量下限 (格点单位): 低于它的解基本是 m->0 的退化解 (cosh 退化成常数,
#: 例如 48x18 的 Xt 在 N=7 下给出 m=0.000000), 只在非 reference 数据源上启用。
ANA_MIN_MASS: float = 1.0e-3


def half_window_for(ns: int, min_half_window: int = ANA_HALF_WINDOW) -> int:
    """
    拟合半宽随格子大小缩放: 窗口 = [center-N, center+N+1), N = round(Ns/3)。

    用 scripts/try_fit_window.py --mode ana 扫过 N 的依赖: ΔM 在 N ≈ Ns/3 处进入平台
    (Ns=32 -> 11, 40 -> 13, 48 -> 16), 而固定 N=7 在 Ns=40/48 上还没收敛
    (48x18 的 Xt 甚至会退化成 m=0), 于是低温柔给出虚高的 ΔM 和巨大的误差。

    注意: 这是"名义"半宽, 小格子上可能宽到拟合崩掉 (32^3 上 N=11 -> [5,28) 时
    V/A 全部样本 chi2/dof 爆掉), 实际使用的半宽请用 choose_half_window()。
    """
    return max(int(min_half_window), int(round(ns / 3.0)))


def choose_half_window(
    jk_by_channel: Dict[str, np.ndarray],
    ns: int,
    min_valid_samples: int = 2,
    min_half_window: int = ANA_HALF_WINDOW,
) -> int:
    """
    选该 ensemble 实际可用的半宽: 从名义值 half_window_for(ns) 往下退,
    直到**每个信道**都至少有 min_valid_samples 个有效样本
    (有效 = chi2/dof <= 100 且 m > ANA_MIN_MASS)。

    为什么要"每个信道都有效": 6 个信道两两组成 4 组 ΔM (V-A / Tt-Xt / S-PS / A-Xt),
    只要有一个信道全灭, 对应的 ΔM 就是 NaN。实测:
      32^3 上名义 N=11 (窗口 [5,28)) 会让 V/A 全灭 -> 退到 N=10;
      48^3x18 ml=0.0020 的 S 道本身报废, 但 N=16 时仍有 2 个样本 -> 保持 N=16
      (若为救 S 而退到 N=10, T-X / X-A 的平台反而被破坏)。
    """
    keys = [ch for ch, jk in jk_by_channel.items() if jk is not None and np.asarray(jk).size]
    if not keys:
        return int(min_half_window)
    for n in range(half_window_for(ns, min_half_window), int(min_half_window) - 1, -1):
        ok = True
        for ch in keys:
            masses, _ = channel_masses_from_jk(
                np.asarray(jk_by_channel[ch], dtype=np.float64),
                center=ns // 2, half_window=n, min_mass=ANA_MIN_MASS,
            )
            if int(np.isfinite(masses).sum()) < int(min_valid_samples):
                ok = False
                break
        if ok:
            return n
    return int(min_half_window)

#: simulate.py 的粗略初值 (拟合失败时回退到自适应初值)
ANA_P0: Dict[str, float] = {"a": 2.5e-5, "m": 0.43058515595986924}

#: ana 从 block 里固定只读 32 行 (对 Ns > 32 的 ensemble 也是 32 行)
ANA_BLOCK_LINES: int = 32


# ------------------------------------------------------------------------------
# 2. 基础工具
# ------------------------------------------------------------------------------
def natural_key(path: str) -> List[int]:
    """自然排序 key (替代 ana 依赖的 natsort)"""
    return [int(t) for t in re.findall(r"\d+", Path(path).name)]


def natural_sorted(paths: Iterable[str]) -> List[str]:
    return sorted(paths, key=natural_key)


def extract_block_values(
    lines: Sequence[str], tag: str, n_lines: int = ANA_BLOCK_LINES, apply_abs: bool = True
) -> List[float]:
    """
    完全复刻 ana 的 extract_dir_block: 精确匹配 tag 后取紧随其后 n_lines 行的第二列 (实部)。
    apply_abs=True 时取绝对值 (ana 的写法); False 时保留符号 (ParseLQCData 管道写法)。
    """
    out: List[float] = []
    in_block = False
    for raw in lines:
        s = raw.strip()
        if s == tag:
            in_block = True
            continue
        if in_block:
            if len(out) >= n_lines:
                break
            parts = re.split(r"\s+", s)
            if len(parts) >= 2:
                try:
                    value = float(parts[1])
                except (ValueError, IndexError):
                    continue
                out.append(abs(value) if apply_abs else value)
    return out


def ana_block_tag(channel: str, source: str = "0/0/0/0") -> str:
    vtype, direction = ANA_CHANNEL_BLOCKS[channel]
    return f"--- {vtype} to {vtype} --- spatial:{direction} --- {source}"


def ana_dir_pattern(intem: int) -> str:
    """ana 使用的目录 glob: 3*x{intem}_b4.17_ms0.040m{ml}"""
    return f"3*x{intem}_b4.17_ms0.040m{{ml}}"


def ana_ensemble_dir(readin_dir: Path, intem: int, ml: str) -> Path:
    """
    ana 的 glob `3*x{intem}` 在本仓库会同时匹配 36x18 与 48x18, 这里显式解析:
    取 config 数最多的候选之外的默认顺序 —— 优先 3{intem}x{intem} (32/36), 否则报错。
    """
    readin_dir = Path(readin_dir)
    if intem in (12, 14, 16):
        return readin_dir / f"32x{intem}_b4.17_ms0.040m{ml}"
    if intem == 18:
        return readin_dir / f"36x18_b4.17_ms0.040m{ml}"
    raise ValueError(f"未知 intem={intem} (ana ccor 只覆盖 12/14/16/18)")


# ------------------------------------------------------------------------------
# 3. ana 三步流程 (提取 -> Jackknife/fold -> 逐样本拟合)
# ------------------------------------------------------------------------------
def per_config_matrix(files: Sequence[str], tag: str,
                      n_lines: int = ANA_BLOCK_LINES, apply_abs: bool = True) -> np.ndarray:
    """(n_lines, n_cfg) 逐构型关联函数矩阵 (ana 口径; apply_abs=False 则保留符号)"""
    columns = []
    for f in files:
        with open(f, "r", encoding="utf-8", errors="ignore") as fh:
            columns.append(extract_block_values(fh.readlines(), tag, n_lines, apply_abs))
    if not columns:
        return np.zeros((0, 0))
    width = min(len(c) for c in columns)
    if width < n_lines:
        columns = [c for c in columns if len(c) == n_lines]
    return np.array(columns, dtype=np.float64).T


def leave_one_out_jackknife(matrix: np.ndarray, n_rows: Optional[int] = None) -> np.ndarray:
    """ana: 对每一行 (时间/空间切片) 做逐构型留一 Jackknife"""
    matrix = np.asarray(matrix, dtype=np.float64)
    n_rows = matrix.shape[0] if n_rows is None else n_rows
    out = np.empty((n_rows, matrix.shape[1]), dtype=np.float64)
    for i in range(n_rows):
        row = matrix[i]
        total = row.sum()
        out[i] = (total - row) / (row.size - 1)
    return out


def symmetrize_about_center(matrix: np.ndarray, center: Optional[int] = None) -> np.ndarray:
    """
    fold: row[x] = (row[x] + row[(2*center - x) mod Ns]) / 2, 即关于 **center** 对称折叠。

    * center = Ns/2 (默认) 时与 ana 的 fold 完全等价, 也与 ParseLQCData 里
      symmetrize_matrix 等价 —— 适用于 Ns = 32 / 36 / 40 / 48 等任意偶数尺寸;
    * center 可以是别的值 (例如复刻 ana 对 36^3 数据仍按 32 行/中心 16 的旧口径),
      此时按 (2*center - x) mod Ns 折叠, row[center] 保持不变。
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError(f"symmetrize_about_center 需要 2D 矩阵 (Ns x N_samples), 得到 shape={matrix.shape}")
    ns = matrix.shape[0]
    center = ns // 2 if center is None else int(center)
    if not (0 <= center < ns):
        raise ValueError(f"对称点 center={center} 超出矩阵行数 Ns={ns}")
    idx = (2 * center - np.arange(ns)) % ns
    # center = Ns/2 时该式自动退化为 ana 的写法: sym[0]=M[0], sym[c]=M[c], 其余取镜像平均
    return 0.5 * (matrix + matrix[idx, :])


def jackknife_mean_err(values: np.ndarray) -> Tuple[float, float]:
    """ana 的 Jackknife 统计: mean 与 sqrt((n-1) * mean((x-mean)^2))"""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    mean = float(v.mean())
    err = float(np.sqrt((v.size - 1) * np.mean((v - mean) ** 2)))
    return mean, err


def _adaptive_p0(x: np.ndarray, y: np.ndarray, center: int) -> Dict[str, float]:
    y0 = float(y[0]) if y.size else 1e-5
    c0 = float(np.cosh(0.3 * (x[0] - center)))
    a0 = y0 / c0 if np.isfinite(c0) and c0 != 0.0 else 1.0
    if not np.isfinite(a0) or a0 == 0.0:
        a0 = 1.0
    return {"a": a0, "m": 0.3}


def fit_sample_masses(
    sym: np.ndarray,
    err: np.ndarray,
    center: Optional[int] = None,
    half_window: int = ANA_HALF_WINDOW,
    chi2_dof_max: float = ANA_CHI2_DOF_MAX,
    abs_values: bool = True,
    p0: Optional[Dict[str, float]] = None,
    min_mass: Optional[float] = None,
) -> np.ndarray:
    """
    对 sym 的每一列 (Jackknife 样本) 做 cosh 拟合, 返回逐样本质量数组 (无效为 NaN)。

    拟合形式与窗口与 ana simulate.py 完全一致: f = a * cosh(m * (x - center)),
    窗口 [center - half_window, center + half_window + 1)。
    ana 用固定初值 ANA_P0; 这里先用它, 失败/发散时回退到自适应初值,
    并同样用 chi2/dof <= chi2_dof_max 过滤。
    min_mass 不为 None 时, 额外丢弃 m <= min_mass 的退化解 (m->0 时 cosh 退化成常数)。
    """
    sym = np.asarray(sym, dtype=np.float64)
    ns, n_samples = sym.shape
    center = ns // 2 if center is None else int(center)
    lo, hi = center - half_window, center + half_window + 1
    if lo < 0 or hi > ns:
        raise ValueError(f"窗口 [{lo},{hi}) 超出 Ns={ns}")

    x = np.arange(ns, dtype=np.float64)[lo:hi]
    prior = dict(ANA_P0 if p0 is None else p0)
    masses = np.full(n_samples, np.nan, dtype=np.float64)

    for j in range(n_samples):
        y = sym[lo:hi, j].copy()
        if abs_values:
            y = np.abs(y)
        if not np.all(np.isfinite(y)):
            continue
        for guess in (prior, _adaptive_p0(x, y, center)):
            try:
                # 强制 chi2 最小二乘 (fitter 由 chi2_least_squares_fit 统一指定, 不可覆盖)
                fit = chi2_least_squares_fit(
                    data=(x, gv.gvar(y, err[lo:hi])),
                    fcn=lambda xx, p: p["a"] * np.cosh(p["m"] * (xx - center)),
                    p0=guess,
                )
                chi2_dof = float(fit.chi2 / fit.dof) if fit.dof > 0 else np.inf
                m = abs(float(fit.p["m"].mean))
                if (np.isfinite(m) and chi2_dof <= chi2_dof_max
                        and (min_mass is None or m > float(min_mass))):
                    masses[j] = m
                break
            except Exception:
                continue
    return masses


def sample_errors(sym: np.ndarray) -> np.ndarray:
    """
    ana 的 dferr: 逐行 sqrt(mean(v^2) - mean(v)^2) * sqrt(n-1)
    (与 jackknife_mean_err 代数等价, 这里刻意保留 ana 的写法以便逐位复现其保存结果)
    """
    out = np.empty(sym.shape[0], dtype=np.float64)
    for i in range(sym.shape[0]):
        v = sym[i]
        v = v[np.isfinite(v)]
        if v.size == 0:
            out[i] = np.nan
            continue
        mean = v.mean()
        # ana 原式 mean(v^2) - mean(v)^2 在相消时可能略小于 0, 夹到 0 再开方
        var = max(float(np.mean(v ** 2) - mean ** 2), 0.0)
        out[i] = np.sqrt(var) * np.sqrt(v.size - 1)
    return out


def channel_masses_from_jk(
    jk_matrix: np.ndarray,
    center: Optional[int] = None,
    half_window: int = ANA_HALF_WINDOW,
    chi2_dof_max: float = ANA_CHI2_DOF_MAX,
    abs_values: bool = True,
    min_mass: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    输入 Jackknife 样本矩阵 (Ns x n_samples), 先 fold 再逐样本 cosh 拟合。
    返回 (逐样本质量, fold 后的矩阵)。
    """
    jk_matrix = np.asarray(jk_matrix, dtype=np.float64)
    sym = symmetrize_about_center(jk_matrix, center)
    err = sample_errors(sym)
    masses = fit_sample_masses(sym, err, center, half_window, chi2_dof_max, abs_values,
                               min_mass=min_mass)
    return masses, sym


# ------------------------------------------------------------------------------
# 4. 参考数据 (ana/dat/ccor) 读取
# ------------------------------------------------------------------------------
def reference_symdata(ccor_dir: Path, ml: str, channel_index: int) -> Optional[np.ndarray]:
    """读取 ana 保存的 symdatasample<ml>.csv (channel_index=0) / symdatasample2<ml>.csv (=1)"""
    suffix = "" if channel_index == 0 else "2"
    path = Path(ccor_dir) / f"symdatasample{suffix}{ml}.csv"
    if not path.exists():
        return None
    return read_csv_matrix(path)


def read_csv_matrix(path: Path) -> np.ndarray:
    import polars as pl

    return pl.read_csv(path, has_header=False, infer_schema_length=0).cast(
        pl.Float64, strict=False
    ).to_numpy()


def load_reference_channel(
    readin_dir: Path, ccor_dir: Path, intem: int, ml: str, channel: str,
    n_lines: int = ANA_BLOCK_LINES,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    复刻 ana 的提取: 从 <ensemble>/Output/test1_lhadrons_*_mesons 取 block,
    abs + 前 n_lines 行 + 逐构型 Jackknife + fold。
    返回 (fold 后矩阵, 逐行误差)。
    """
    ens_dir = ana_ensemble_dir(Path(readin_dir), intem, ml)
    files = natural_sorted(glob.glob(str(ens_dir / "Output" / "test1_lhadrons_*_mesons")))
    if not files:
        raise FileNotFoundError(f"找不到 ana 输入文件: {ens_dir}/Output/test1_lhadrons_*_mesons")
    mat = per_config_matrix(files, ana_block_tag(channel), n_lines)
    jk = leave_one_out_jackknife(mat, n_rows=n_lines)
    sym = symmetrize_about_center(jk, center=n_lines // 2)
    return sym, sample_errors(sym)


# ------------------------------------------------------------------------------
# 5. ParseLQCData meson_scan 结果读取
# ------------------------------------------------------------------------------
def ana_channel_components(ana_channel: str) -> List[Tuple[str, str]]:
    """ana 脚本里给该信道的完整 (type, dir) 列表 (与仓库 CHANNEL_CONFIGS 一致)"""
    from docs.physics_setup import CHANNEL_CONFIGS

    name = MESON_SCAN_CHANNEL.get(ana_channel, ana_channel)
    return [(c["type"], c["dir"]) for c in CHANNEL_CONFIGS[name]]


def per_config_matrix_channels(
    files: Sequence[str],
    components: Sequence[Tuple[str, str]],
    n_lines: int = ANA_BLOCK_LINES,
    apply_abs: bool = True,
    average: bool = True,
) -> np.ndarray:
    """
    多分量逐构型矩阵: 对 components 里每个 (type, dir) 取 block;
    average=True 时对所有分量求平均 (ParseLQCData 信道定义), 否则只用最后一个分量
    (ana 脚本实际用到的那个)。
    """
    mats = [
        per_config_matrix(files, f"--- {t} to {t} --- spatial:{d} --- 0/0/0/0", n_lines, apply_abs)
        for t, d in components
    ]
    if not mats:
        return np.zeros((0, 0))
    if average:
        return sum(mats) / len(mats)
    return mats[-1]


def load_meson_scan_channel(case_dir: Path, channel: str) -> Optional[np.ndarray]:
    """
    读取 output/meson_scan 里保存的 Jackknife 折叠矩阵 (Ns x n_bins)。
    列 = 逐 bin 留一 Jackknife 样本, 与 ana 的逐构型 Jackknife 语义一致。
    channel 既可用 ana 短名 (V/A/Tt/Xt/S/PS), 也可用 meson_scan 名 (Vec/AV/...)。
    """
    import polars as pl

    name = MESON_SCAN_CHANNEL.get(channel, channel)
    path = Path(case_dir) / f"jk_{name}.csv"
    if not path.exists():
        return None
    return pl.read_csv(path, has_header=False, infer_schema_length=0).cast(
        pl.Float64, strict=False
    ).to_numpy()
