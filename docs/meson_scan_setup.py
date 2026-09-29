"""
docs/meson_scan_setup.py
--------------------------------------------------------------------------------
b4.17 有限温度扫描任务 (Nt = 12 / 14 / 16 / 18) 的物理设定与任务清单
(Single Source of Truth for the beta=4.17 Nt-scan task)

任务范围: ms = 0.040, ml = 0.0020 / 0.0035 / 0.0070 / 0.0120

    序号  Ns x Nt     ms       ml                                     T_ref(MeV)
    00   32 x 12     0.040    0.0020, 0.0035, 0.0070, 0.0120          204.00
    01   32 x 14     0.040    0.0020, 0.0035, 0.0070, 0.0120          174.86
    02   32 x 16     0.040    0.0020, 0.0035, 0.0070, 0.0120          153.00
    03   36 x 18     0.040    0.0020, 0.0035, 0.0070, 0.0120          136.00
    04   40 x 16     0.040    0.0020, 0.0035, 0.0070, 0.0120          153.00
    05   48 x 18     0.040    0.0020, 0.0035, 0.0070, 0.0120          136.00

表格顺序即结果落盘顺序 (order 前缀保证文件系统与汇总表同序)。

约定说明 (与仓库既有 4.17 @ 48^3x16 分析严格一致, 已数值复核):
  * C++ 核心库扫描 "<case>/Output/test1_lhadrons_*_mesons_multi_src";
  * 每个 multi_src 文件含 16 个源位置 (4x4 空间网格), 每个信道 4 个 block:
        temporal (Nt 行) / spatial:DIRX,DIRY,DIRZ (各 Ns 行);
  * 抽取信道只使用 spatial:DIR* block, 因此 num_lines == Ns,
    得到的关联函数是沿极化方向的空间关联函数 C(x), x = 0 ... Ns-1;
  * 环形移位平均后 C(x) 满足周期边界 cosh 形式
        C(x) = A * cosh(m * (x - Ns/2)),
    故有效质量与平台拟合的对称点必须取 Ns/2 (32->16, 36->18, 40->20, 48->24),
    而不是 48^3 分析里写死的 24。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from docs.physics_setup import (
    CHANNELS,
    DEFAULT_READIN_DIR,
    ENSEMBLE_CONFIGS,
    FIT_SLICE_END,
    MRES_TABLE,
    MULTI_FIT_SLICES,
    OUTPUT_ROOT,
    SINGLE_FIT_SLICES,
    TEMP_MAP,
    calculate_temperature,
    get_binsize,
)

# ------------------------------------------------------------------------------
# 1. 任务常量
# ------------------------------------------------------------------------------
BETA_KEY: str = "17"           # 与 MULTI_FIT_SLICES / SINGLE_FIT_SLICES 的键一致
BETA_VALUE: float = 4.17
MS: float = 0.040

#: 48^3x16 (beta=4.17) 基准分析的 cosh 对称点, 现有 slice 表均以它为参考
REFERENCE_HALF: int = 24

#: 结果根目录: output/meson_scan/b4.17
SCAN_ROOT: Path = OUTPUT_ROOT / "meson_scan" / f"b4.{BETA_KEY}"

SOURCE_TAGS: Dict[bool, str] = {False: "multisrc", True: "singlesrc"}
SOURCE_ORDER: Tuple[bool, ...] = (False, True)   # 多源优先 (已与既有结果数值复核)


@dataclass(frozen=True)
class ScanCase:
    """一组格点规模 (Ns x Nt) 对应的任务条目"""

    order: int
    ns: int
    nt: int
    mls: Tuple[float, ...]
    temp_ref_mev: float
    ms: float = MS
    beta: float = BETA_VALUE

    @property
    def key(self) -> str:
        return f"{self.ns}x{self.nt}"

    @property
    def label(self) -> str:
        return f"{self.ns}^3x{self.nt}"

    @property
    def half(self) -> int:
        return self.ns // 2

    @property
    def temperature_mev(self) -> float:
        """仓库统一口径的物理温度 T = 1/(a(beta) * Nt)"""
        return calculate_temperature(self.beta, self.nt)

    @property
    def a_inv_mev(self) -> float:
        """由 TEMP_MAP 反推的格距倒数 a^-1 = T(Nt=16) * 16 (MeV)"""
        return float(TEMP_MAP[f"{self.beta:.2f}"]) * 16.0

    def dirname(self, ml: float) -> str:
        return f"{self.ns}x{self.nt}_b{self.beta:.2f}_ms{self.ms:.3f}m{ml:.4f}"

    def output_name(self, ml: float) -> str:
        return f"{self.order:02d}_{self.ns}x{self.nt}_ms{self.ms:.3f}_m{ml:.4f}"

    def sub_cases(self) -> Iterator[Tuple[int, float]]:
        """按任务表给定的 ml 顺序产出 (ml_index, ml)"""
        for ml_index, ml in enumerate(self.mls):
            yield ml_index, ml


#: 任务表 (顺序即最终结果顺序, 与用户给定表格逐行对应)
SCAN_CASES: Tuple[ScanCase, ...] = (
    ScanCase(order=0, ns=32, nt=12, mls=(0.0020, 0.0035, 0.0070, 0.0120), temp_ref_mev=204.00),
    ScanCase(order=1, ns=32, nt=14, mls=(0.0020, 0.0035, 0.0070, 0.0120), temp_ref_mev=174.86),
    ScanCase(order=2, ns=32, nt=16, mls=(0.0020, 0.0035, 0.0070, 0.0120), temp_ref_mev=153.00),
    ScanCase(order=3, ns=36, nt=18, mls=(0.0020, 0.0035, 0.0070, 0.0120), temp_ref_mev=136.00),
    ScanCase(order=4, ns=40, nt=16, mls=(0.0020, 0.0035, 0.0070, 0.0120), temp_ref_mev=153.00),
    ScanCase(order=5, ns=48, nt=18, mls=(0.0020, 0.0035, 0.0070, 0.0120), temp_ref_mev=136.00),
)

#: 本任务使用的 6 个狄拉克双线性信道 (顺序 = 汇总表内信道顺序)
CHANNEL_ORDER: Tuple[str, ...] = tuple(CHANNELS)


# ------------------------------------------------------------------------------
# 2. 目录解析 (数字开头 + b4.17 参数)
# ------------------------------------------------------------------------------
#: 任务表目录名: 32x12_b4.17_ms0.040m0.0020
SCAN_DIR_RE = re.compile(
    r"^(?P<ns>\d+)x(?P<nt>\d+)_b(?P<beta>[\d.]+)_ms(?P<ms>[\d.]+)m(?P<ml>[\d.]+)$"
)
#: 任意含 b4.17 的目录名, 用于 manifest 里登记 "发现但未纳入任务表" 的目录
ANY_B417_RE = re.compile(r"^[0-9].*4\.17.*$")


def readin_root(readin_dir: Optional[Path] = None) -> Path:
    return Path(readin_dir) if readin_dir is not None else Path(DEFAULT_READIN_DIR)


def case_input_dir(case: ScanCase, ml: float, readin_dir: Optional[Path] = None) -> Path:
    """ensemble 目录: <readin>/32x12_b4.17_ms0.040m0.0020"""
    return readin_root(readin_dir) / case.dirname(ml)


def case_output_dir(
    case: ScanCase, ml: float, is_single_source: bool, out_root: Optional[Path] = None
) -> Path:
    """结果目录: <scan_root>/<source>/cases/00_32x12_ms0.040_m0.0020"""
    root = Path(out_root) if out_root is not None else SCAN_ROOT
    return root / SOURCE_TAGS[is_single_source] / "cases" / case.output_name(ml)


def scan_root(out_root: Optional[Path] = None) -> Path:
    return Path(out_root) if out_root is not None else SCAN_ROOT


def discover_b417_dirs(readin_dir: Optional[Path] = None) -> Tuple[List[str], List[str]]:
    """返回 (已纳入任务表的目录名, 数字开头且含 b4.17 但未纳入任务表的目录名)"""
    root = readin_root(readin_dir)
    used: List[str] = []
    for case in SCAN_CASES:
        for _, ml in case.sub_cases():
            used.append(case.dirname(ml))
    used_set = set(used)

    extras: List[str] = []
    if root.exists():
        for child in sorted(root.iterdir()):
            if not child.is_dir():
                continue
            name = child.name
            if name in used_set:
                continue
            if ANY_B417_RE.match(name):
                extras.append(name)
    return used, extras


def resolve_tasks(
    readin_dir: Optional[Path] = None,
    mls: Optional[Sequence[float]] = None,
    case_keys: Optional[Sequence[str]] = None,
) -> List[dict]:
    """
    生成有序任务列表 (每个元素 = 一个 ensemble), 缺失目录会被标记 missing 而不是静默跳过。
    """
    root = readin_root(readin_dir)
    want_mls = {round(float(m), 4) for m in mls} if mls else None
    want_cases = set(case_keys) if case_keys else None

    tasks: List[dict] = []
    for case in SCAN_CASES:
        if want_cases is not None and case.key not in want_cases:
            continue
        for ml_index, ml in case.sub_cases():
            if want_mls is not None and round(ml, 4) not in want_mls:
                continue
            in_dir = case_input_dir(case, ml, root)
            out_dir = in_dir / "Output"
            tasks.append(
                {
                    "order": case.order,
                    "ml_index": ml_index,
                    "ns": case.ns,
                    "nt": case.nt,
                    "half": case.half,
                    "ms": case.ms,
                    "beta": case.beta,
                    "ml": ml,
                    "case_key": case.key,
                    "label": case.label,
                    "output_name": case.output_name(ml),
                    "dir_name": case.dirname(ml),
                    "input_dir": str(out_dir),
                    "temperature_mev": case.temperature_mev,
                    "temp_ref_mev": case.temp_ref_mev,
                    "a_inv_mev": case.a_inv_mev,
                    "missing": not out_dir.is_dir(),
                    "n_configs": _count_configs(out_dir),
                }
            )
    return tasks


def _count_configs(output_dir: Path) -> int:
    """统计 multi_src 构型文件数 (C++ 库读取的文件集合)"""
    if not output_dir.is_dir():
        return 0
    return sum(1 for _ in output_dir.glob("*_mesons_multi_src"))


# ------------------------------------------------------------------------------
# 3. 平台拟合窗口策略
# ------------------------------------------------------------------------------
#: 窗口扫描的最小点数 (含两端)
MIN_WINDOW_POINTS: int = 5


def fit_slice_table(is_single_source: bool) -> Dict[str, int]:
    return SINGLE_FIT_SLICES[BETA_KEY] if is_single_source else MULTI_FIT_SLICES[BETA_KEY]


def reference_window(case: ScanCase, channel: str, is_single_source: bool) -> Tuple[int, int]:
    """
    镜像窗口: 把 4.17 @ 48^3x16 既有分析的 [start, 25) 平移到本 ensemble 的对称点,
    即保持 "距 cosh 对称点的距离" 不变, 起点夹紧到 >= 1。
    """
    start_ref = fit_slice_table(is_single_source)[channel]
    half = case.half
    start = max(1, half - (REFERENCE_HALF - start_ref))
    end = half + (FIT_SLICE_END - REFERENCE_HALF)     # = half + 1
    end = min(max(end, start + 1), half + 1)
    return start, end


def scan_windows(case: ScanCase, min_points: int = MIN_WINDOW_POINTS) -> List[Tuple[int, int]]:
    """窗口扫描网格: 固定终点 = 对称点 + 1, 起点从 1 扫到 Ns/2 - min_points + 1"""
    half = case.half
    end = half + 1
    last_start = half - min_points + 1
    return [(s, end) for s in range(1, max(last_start, 1) + 1)]


def binsize_for(channel: str, is_single_source: bool) -> int:
    return get_binsize(BETA_KEY, channel, is_single_source)


def residual_mass() -> float:
    """DWF 残余质量 (beta = 4.17)"""
    return float(MRES_TABLE[f"{BETA_VALUE:.2f}"]["m_residual"])


def zm_factor() -> float:
    return float(MRES_TABLE[f"{BETA_VALUE:.2f}"]["zm"])


def known_ensemble_dirs() -> List[str]:
    """已登记在 ENSEMBLE_CONFIGS 里的 b4.17 相关目录名 (供 manifest 参考)"""
    return sorted(name for name in ENSEMBLE_CONFIGS if "4.17" in name)
