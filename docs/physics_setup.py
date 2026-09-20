"""
docs/physics_setup.py
--------------------------------------------------------------------------------
格点 QCD 全局物理参数、构型字典与对称性映射总表 (Single Source of Truth)
所有 scripts/, tests/, tools/, src/ 统一自此处导入物理设定，杜绝硬编码与参数碎片化。
--------------------------------------------------------------------------------
"""

from typing import Dict, List, Tuple
import json
from pathlib import Path

# ==============================================================================
# 0. 全局目录与路径统一数据配置 (Single Source of Truth for Directories)
# ==============================================================================
PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
DOCS_DIR: Path = PROJECT_ROOT / "docs"
FIGURES_DIR: Path = DOCS_DIR / "figures"
REPORT_PATH: Path = DOCS_DIR / "verification_report.md"

ANA_ROOT: Path = Path("/Users/junxiongnie/code/ana")
DEFAULT_DATA_DIR: Path = PROJECT_ROOT / "data"
DEFAULT_READIN_DIR: Path = DEFAULT_DATA_DIR / "readin" if (DEFAULT_DATA_DIR / "readin").exists() else Path("/Users/junxiongnie/code/ana/dat/readin")

OUTPUT_ROOT: Path = PROJECT_ROOT / "output"
OUTPUT_PICKDATA_DIR: Path = OUTPUT_ROOT / "pickdata"
OUTPUT_PICKDATA_SINGLE_DIR: Path = OUTPUT_ROOT / "pickdata-singlesrc"
OUTPUT_RATIO_DIR: Path = OUTPUT_ROOT / "ratio_results"
OUTPUT_RATIO_SINGLE_DIR: Path = OUTPUT_ROOT / "ratio_results-singlesrc"
OUTPUT_SIMULATE_DIR: Path = OUTPUT_ROOT / "simulateresult"
OUTPUT_SIMULATE_SINGLE_DIR: Path = OUTPUT_ROOT / "simulateresult-singlesrc"
OUTPUT_CONDENSATE_DIR: Path = OUTPUT_ROOT / "condensate"

# 1. 耦合常数 (Beta) 与物理温度 (MeV) 对应字典
TEMP_MAP: Dict[str, float] = {
    "4.13": 138.23,
    "4.15": 145.75,
    "4.17": 153.31,
    "4.18": 157.03,  # 赝临界相变区附近
    "4.20": 164.55,
    "4.23": 176.43,
    "4.30": 202.52,
    "4.405": 241.60,
}

# 核心分析覆盖的 7 组 Beta 编号 (字符串形式)
BETAS: List[str] = ["13", "15", "17", "18", "20", "23", "30"]

# 覆盖的 6 种狄拉克介子流信道
CHANNELS: List[str] = ["AV", "S", "Tt", "PS", "Xt", "Vec"]

# 相变过渡带温度范围 (Transition Region, MeV)
TRANSITION_REGION: Tuple[float, float] = (155.5, 160.5)

# 2. 介子流算子及其空间极化投影映射字典 (Channel Configurations)
CHANNEL_CONFIGS: Dict[str, List[dict]] = {
    "AV": [
        {'type': 'AVector2', 'dir': 'DIRX'},
        {'type': 'AVector3', 'dir': 'DIRX'},
        {'type': 'AVector3', 'dir': 'DIRY'},
        {'type': 'AVector1', 'dir': 'DIRY'},
        {'type': 'AVector1', 'dir': 'DIRZ'},
        {'type': 'AVector2', 'dir': 'DIRZ'}
    ],
    "S": [
        {'type': 'TVector4', 'dir': 'DIRZ'},
        {'type': 'TVector4', 'dir': 'DIRX'},
        {'type': 'TVector4', 'dir': 'DIRY'}
    ],
    "Tt": [
        {'type': 'TVector1', 'dir': 'DIRX'},
        {'type': 'TVector3', 'dir': 'DIRZ'},
        {'type': 'TVector2', 'dir': 'DIRY'}
    ],
    "PS": [
        {'type': 'TAVector4', 'dir': 'DIRZ'},
        {'type': 'TAVector4', 'dir': 'DIRX'},
        {'type': 'TAVector4', 'dir': 'DIRY'}
    ],
    "Xt": [
        {'type': 'TAVector1', 'dir': 'DIRX'},
        {'type': 'TAVector3', 'dir': 'DIRZ'},
        {'type': 'TAVector2', 'dir': 'DIRY'}
    ],
    "Vec": [
        {'type': 'Vector2', 'dir': 'DIRX'},
        {'type': 'Vector3', 'dir': 'DIRX'},
        {'type': 'Vector3', 'dir': 'DIRY'},
        {'type': 'Vector1', 'dir': 'DIRY'},
        {'type': 'Vector1', 'dir': 'DIRZ'},
        {'type': 'Vector2', 'dir': 'DIRZ'}
    ]
}

# 3. 对称性破缺与恢复物理监测信道对字典 (Symmetry Breaking Signatures)
SYMMETRY_PAIRS: List[dict] = [
    {
        "name": "Chiral SU(2)_L x SU(2)_R",
        "ch1": "Vec",
        "ch2": "AV",
        "label": r"$SU(2)_L \times SU(2)_R$ (Vec - AV)",
        "color": "#69b3a2",
        "marker": "o"
    },
    {
        "name": "Axial Anomaly U(1)_A (Scalar-Pseudoscalar)",
        "ch1": "S",
        "ch2": "PS",
        "label": r"$U(1)_A$ (S - PS)",
        "color": "#ba68c8",
        "marker": "s"
    },
    {
        "name": "Axial Anomaly U(1)_A (Tensor Channels)",
        "ch1": "Xt",
        "ch2": "Tt",
        "label": r"$U(1)_A$ (Xt - Tt)",
        "color": "#64b5f6",
        "marker": "^"
    },
    {
        "name": "Chiral-Spin SU(2)_CS",
        "ch1": "AV",
        "ch2": "Tt",
        "label": r"$SU(2)_{CS}$ (AV - Tt)",
        "color": "#ffb74d",
        "marker": "v"
    }
]

# 4. 平台拟合时间切片区间字典 (Plateau Fit Intervals: [start, 25])
FIT_SLICE_END: int = 25

MULTI_FIT_SLICES: Dict[str, Dict[str, int]] = {
    '13': {'AV': 16, 'S': 14, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
    '15': {'AV': 13, 'S': 10, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 13},
    '17': {'AV': 14, 'S': 8,  'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 14},
    '18': {'AV': 16, 'S': 18, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
    '20': {'AV': 15, 'S': 19, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 15},
    '23': {'AV': 16, 'S': 8,  'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
    '30': {'AV': 16, 'S': 8,  'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
}

SINGLE_FIT_SLICES: Dict[str, Dict[str, int]] = {
    '13': {'AV': 16, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
    '15': {'AV': 13, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 13},
    '17': {'AV': 14, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 14},
    '18': {'AV': 16, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
    '20': {'AV': 15, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 15},
    '23': {'AV': 16, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
    '30': {'AV': 16, 'S': 8, 'Tt': 10, 'PS': 8, 'Xt': 10, 'Vec': 16},
}

# 5. 格点 Domain Wall Fermion 质量参数与残余质量查表 (Single Source of Truth)
MRES_TABLE: Dict[str, dict] = {
    "4.13": {
        "m_light": 0.000805,
        "m_strange": 0.043547,
        "m_residual": 0.000731464,
        "zm": 0.937703,
    },
    "4.15": {
        "m_light": 0.00093,
        "m_strange": 0.040843,
        "m_residual": 0.000498493,
        "zm": 0.95239,
    },
    "4.17": {
        "m_light": 0.001001,
        "m_strange": 0.0384,
        "m_residual": 0.000339722,
        "zm": 0.966247,
    },
    "4.18": {
        "m_light": 0.001022,
        "m_strange": 0.037265,
        "m_residual": 0.000280451,
        "zm": 0.972899,
    },
    "4.20": {
        "m_light": 0.001041,
        "m_strange": 0.03515,
        "m_residual": 0.000191127,
        "zm": 0.98571,
    },
    "4.23": {
        "m_light": 0.001033,
        "m_strange": 0.032315,
        "m_residual": 0.000107528,
        "zm": 1.00385,
    },
    "4.30": {
        "m_light": 0.000939,
        "m_strange": 0.02693,
        "m_residual": 0.0000280963,
        "zm": 1.04224,
    },
    "4.405": {
        "m_light": 0.00076,
        "m_strange": 0.021032,
        "m_residual": 0.00000375269,
        "zm": 1.09274,
    },
}

# 5.1 格点数据集目录命名结构字典 (Ensemble Configurations & Directory Patterns)
ENSEMBLE_CONFIGS: Dict[str, dict] = {
    "L32T12_beta4.17ms0.040m0.0020": {
        "ns": 32, "nt": 12, "beta": 4.17, "ms": 0.040, "ml": 0.0020,
        "mres": 0.000339722, "zm": 0.966247, "aliases": ["L32T12beta4.17"]
    },
    "L32T14_beta4.17ms0.040m0.0020": {
        "ns": 32, "nt": 14, "beta": 4.17, "ms": 0.040, "ml": 0.0020,
        "mres": 0.000339722, "zm": 0.966247, "aliases": []
    },
    "L40T16_beta4.17ms0.040m0.0020": {
        "ns": 40, "nt": 16, "beta": 4.17, "ms": 0.040, "ml": 0.0020,
        "mres": 0.000339722, "zm": 0.966247, "aliases": []
    },
    "L48T18_beta4.17ms0.040m0.0020": {
        "ns": 48, "nt": 18, "beta": 4.17, "ms": 0.040, "ml": 0.0020,
        "mres": 0.000339722, "zm": 0.966247, "aliases": []
    },
    "L48T18_beta4.17ms0.040m0.0020_2": {
        "ns": 48, "nt": 18, "beta": 4.17, "ms": 0.040, "ml": 0.0020,
        "mres": 0.000339722, "zm": 0.966247, "aliases": []
    },
    "L48T16beta4.13ms0.043547m0.000805": {
        "ns": 48, "nt": 16, "beta": 4.13, "ms": 0.043547, "ml": 0.000805,
        "mres": 0.000731464, "zm": 0.937703, "aliases": ["48x16b4.13", "13", "4.13"]
    },
    "L48T16beta4.15ms0.040843m0.000930": {
        "ns": 48, "nt": 16, "beta": 4.15, "ms": 0.040843, "ml": 0.000930,
        "mres": 0.000498493, "zm": 0.95239, "aliases": ["48x16b4.15", "15", "4.15"]
    },
    "L48T16beta4.17ms0.038400m0.001001": {
        "ns": 48, "nt": 16, "beta": 4.17, "ms": 0.038400, "ml": 0.001001,
        "mres": 0.000339722, "zm": 0.966247, "aliases": ["48x16b4.17", "17", "4.17"]
    },
    "L48T16beta4.18ms0.037265m0.001022": {
        "ns": 48, "nt": 16, "beta": 4.18, "ms": 0.037265, "ml": 0.001022,
        "mres": 0.000280451, "zm": 0.972899, "aliases": ["48x16b4.18", "18", "4.18"]
    },
    "L48T16beta4.20ms0.035150m0.001041": {
        "ns": 48, "nt": 16, "beta": 4.20, "ms": 0.035150, "ml": 0.001041,
        "mres": 0.000191127, "zm": 0.98571, "aliases": ["48x16b4.20", "20", "4.20"]
    },
    "L48T16beta4.23ms0.032315m0.001033": {
        "ns": 48, "nt": 16, "beta": 4.23, "ms": 0.032315, "ml": 0.001033,
        "mres": 0.000107528, "zm": 1.00385, "aliases": ["48x16b4.23", "23", "4.23"]
    },
    "L48T16beta4.30ms0.026930m0.000939": {
        "ns": 48, "nt": 16, "beta": 4.30, "ms": 0.026930, "ml": 0.000939,
        "mres": 0.0000280963, "zm": 1.04224, "aliases": ["48x16b4.30", "30", "4.30"]
    },
    "L48T16beta4.405ms0.021032m0.000760": {
        "ns": 48, "nt": 16, "beta": 4.405, "ms": 0.021032, "ml": 0.000760,
        "mres": 0.00000375269, "zm": 1.09274, "aliases": ["48x16b4.405", "405", "4.405"]
    }
}

# 手征凝聚物理配置列表 (Chiral Condensate Configurations)
CONDENSATE_CONFIGS: List[dict] = [
    {
        "beta": "4.13",
        "ml": 0.000805,
        "ms": 0.043547,
        "mres": 0.000731464,
        "zm": 0.937703,
        "pbp_l": 0.00166880,
        "pbp_l_err": 0.00000634,
        "pbp_s": 0.01942524,
        "pbp_s_err": 0.00000135,
        "pbp_rm": 0.00106082,
        "pbp_rm_err": 0.00000674,
    },
    {
        "beta": "4.15",
        "ml": 0.00093,
        "ms": 0.040843,
        "mres": 0.000498493,
        "zm": 0.952390,
        "pbp_l": 0.00135913,
        "pbp_l_err": 0.00000644,
        "pbp_s": 0.01808069,
        "pbp_s_err": 0.00000128,
        "pbp_rm": 0.00077109,
        "pbp_rm_err": 0.00000675,
    },
    {
        "beta": "4.17",
        "ml": 0.001001,
        "ms": 0.0384,
        "mres": 0.000339722,
        "zm": 0.966247,
        "pbp_l": 0.00107396,
        "pbp_l_err": 0.00000642,
        "pbp_s": 0.01686886,
        "pbp_s_err": 0.00000127,
        "pbp_rm": 0.00050728,
        "pbp_rm_err": 0.00000662,
    },
    {
        "beta": "4.18",
        "ml": 0.001022,
        "ms": 0.037265,
        "mres": 0.000280451,
        "zm": 0.972899,
        "pbp_l": 0.00097061,
        "pbp_l_err": 0.00000642,
        "pbp_s": 0.01629887,
        "pbp_s_err": 0.00000134,
        "pbp_rm": 0.00041649,
        "pbp_rm_err": 0.00000659,
    },
    {
        "beta": "4.20",
        "ml": 0.001041,
        "ms": 0.03515,
        "mres": 0.000191127,
        "zm": 0.985710,
        "pbp_l": 0.00076708,
        "pbp_l_err": 0.00000554,
        "pbp_s": 0.01524584,
        "pbp_s_err": 0.00000127,
        "pbp_rm": 0.00023896,
        "pbp_rm_err": 0.00000561,
    },
    {
        "beta": "4.23",
        "ml": 0.001033,
        "ms": 0.032315,
        "mres": 0.000107528,
        "zm": 1.003850,
        "pbp_l": 0.00058018,
        "pbp_l_err": 0.00000383,
        "pbp_s": 0.01384239,
        "pbp_s_err": 0.00000120,
        "pbp_rm": 0.00009289,
        "pbp_rm_err": 0.00000381,
    },
    {
        "beta": "4.30",
        "ml": 0.000939,
        "ms": 0.02693,
        "mres": 0.0000280963,
        "zm": 1.042240,
        "pbp_l": 0.00041210,
        "pbp_l_err": 0.00000139,
        "pbp_s": 0.01124911,
        "pbp_s_err": 0.00000090,
        "pbp_rm": 0.00000820,
        "pbp_rm_err": 0.00000133,
    },
    {
        "beta": "4.405",
        "ml": 0.00076,
        "ms": 0.021032,
        "mres": 0.00000375269,
        "zm": 1.092740,
        "pbp_l": 0.00031129,
        "pbp_l_err": 0.00000003,
        "pbp_s": 0.00856161,
        "pbp_s_err": 0.00000027,
        "pbp_rm": 0.00000040,
        "pbp_rm_err": 0.00000002,
    },
]

# 6. 重采样 Binsize 配置函数
def get_binsize(beta: str, channel: str, is_single_source: bool) -> int:
    """根据模式、温度和信道严格返回对应的 binsize"""
    if is_single_source:
        return 5 if beta in ['18', '20', '23', '30'] else 4
    else:
        return 5 if (beta == '18' and channel == 'S') else 4

# 7. 残余手征质量函数 (Domain Wall Fermion Residual Mass)
def calculate_residual_mass(beta_val) -> float:
    """计算 DWF 残余质量: 优先查表，缺失时使用拟合公式 m_res = 2.547e28 * exp(-17.559 * beta)"""
    b_float = float(beta_val)
    beta_str = f"{b_float:.2f}" if abs(b_float - 4.405) > 1e-4 else "4.405"
    if str(beta_val) in MRES_TABLE:
        return float(MRES_TABLE[str(beta_val)]["m_residual"])
    if beta_str in MRES_TABLE:
        return float(MRES_TABLE[beta_str]["m_residual"])
    import numpy as np
    return float(2.547e28 * np.exp(-17.559 * b_float))

def calculate_temperature(beta_val: float, nt: int = 16) -> float:
    """根据有限温度公式 T = 1 / (a(beta) * nt) 计算物理温度 (MeV)"""
    b_val = float(beta_val)
    b_key = f"{b_val:.2f}" if abs(b_val - 4.405) > 1e-4 else "4.405"
    t_16 = TEMP_MAP.get(b_key, TEMP_MAP.get(str(beta_val), 153.31))
    return float(t_16 * 16.0 / float(nt)) if nt > 0 else float(t_16)

# 8. 格点目录结构正则解析与智能定位函数 (Directory Pattern Matcher & Resolver)
import re

ENSEMBLE_REGEX = re.compile(
    r"^L(?P<L>\d+)T(?P<T>\d+)(?:_)?beta(?P<beta>[\d\.]+)ms(?P<ms>[\d\.]+)m(?:0)?(?P<ml>[\d\.]+)(?:_(?P<subid>\w+))?$"
)
LEGACY_REGEX = re.compile(
    r"^(?P<L>\d+)x(?P<T>\d+)b(?P<beta>[\d\.]+)$"
)

def parse_ensemble_dirname(dirname: str) -> dict:
    """从目录名称解析格点规模与夸克质量参数及对应温度"""
    # 1. 优先查表精确匹配
    if dirname in ENSEMBLE_CONFIGS:
        info = dict(ENSEMBLE_CONFIGS[dirname])
        info["dirname"] = dirname
        info["temperature"] = calculate_temperature(info["beta"], info.get("nt", 16))
        return info

    # 2. 正则解析新格式: L32T12_beta4.17ms0.040m0.0020 或 L48T16beta4.13ms0.043547m0.000805
    m = ENSEMBLE_REGEX.match(dirname)
    if m:
        gd = m.groupdict()
        b_val = float(gd["beta"])
        t_val = int(gd["T"])
        return {
            "dirname": dirname,
            "ns": int(gd["L"]),
            "nt": t_val,
            "beta": b_val,
            "ms": float(gd["ms"]),
            "ml": float(gd["ml"]),
            "mres": calculate_residual_mass(b_val),
            "zm": float(MRES_TABLE.get(f"{b_val:.2f}", {}).get("zm", 1.0)),
            "temperature": calculate_temperature(b_val, t_val),
        }

    # 3. 解析遗留格式: 48x16b4.17 或 L32T12beta4.17
    m_leg = LEGACY_REGEX.match(dirname)
    if m_leg:
        gd = m_leg.groupdict()
        b_val = float(gd["beta"])
        b_key = f"{b_val:.2f}"
        mres_info = MRES_TABLE.get(b_key, {})
        return {
            "dirname": dirname,
            "ns": int(gd["L"]),
            "nt": int(gd["T"]),
            "beta": b_val,
            "ms": float(mres_info.get("m_strange", 0.04)),
            "ml": float(mres_info.get("m_light", 0.001)),
            "mres": float(mres_info.get("m_residual", calculate_residual_mass(b_val))),
            "zm": float(mres_info.get("zm", 1.0)),
            "temperature": calculate_temperature(b_val, int(gd["T"])),
        }

    return {}

def resolve_dataset_dir(base_dir: Path, target: str) -> Path:
    """在 base_dir (如 data/readin/) 检索目标数据集，支持别名、beta 和完整目录名"""
    # 直接存在
    direct = base_dir / target
    if direct.exists():
        return direct

    # 别名/名称在 ENSEMBLE_CONFIGS 中的匹配
    for name, cfg in ENSEMBLE_CONFIGS.items():
        if target == name or target in cfg.get("aliases", []):
            candidate = base_dir / name
            if candidate.exists():
                return candidate
            for alias in cfg.get("aliases", []):
                if (base_dir / alias).exists():
                    return base_dir / alias

    # 模糊匹配 beta
    target_clean = target.replace("beta", "").replace("b", "")
    for child in base_dir.iterdir():
        if not child.is_directory():
            continue
        cname = child.name
        if target in cname:
            return child
        if target_clean in cname:
            return child

    return base_dir / target

# 导出为标准 JSON 便于 C++ 及其他语言读取
def export_to_json(out_path: Path) -> None:
    data = {
        "DIRECTORIES": {
            "PROJECT_ROOT": str(PROJECT_ROOT),
            "DOCS_DIR": str(DOCS_DIR),
            "FIGURES_DIR": str(FIGURES_DIR),
            "REPORT_PATH": str(REPORT_PATH),
            "ANA_ROOT": str(ANA_ROOT),
            "DEFAULT_DATA_DIR": str(DEFAULT_DATA_DIR),
            "DEFAULT_READIN_DIR": str(DEFAULT_READIN_DIR),
            "OUTPUT_ROOT": str(OUTPUT_ROOT),
            "OUTPUT_PICKDATA_DIR": str(OUTPUT_PICKDATA_DIR),
            "OUTPUT_PICKDATA_SINGLE_DIR": str(OUTPUT_PICKDATA_SINGLE_DIR),
            "OUTPUT_RATIO_DIR": str(OUTPUT_RATIO_DIR),
            "OUTPUT_RATIO_SINGLE_DIR": str(OUTPUT_RATIO_SINGLE_DIR),
            "OUTPUT_SIMULATE_DIR": str(OUTPUT_SIMULATE_DIR),
            "OUTPUT_SIMULATE_SINGLE_DIR": str(OUTPUT_SIMULATE_SINGLE_DIR),
            "OUTPUT_CONDENSATE_DIR": str(OUTPUT_CONDENSATE_DIR),
        },
        "TEMP_MAP": TEMP_MAP,
        "BETAS": BETAS,
        "CHANNELS": CHANNELS,
        "TRANSITION_REGION": list(TRANSITION_REGION),
        "CHANNEL_CONFIGS": CHANNEL_CONFIGS,
        "CONDENSATE_CONFIGS": CONDENSATE_CONFIGS,
        "MRES_TABLE": MRES_TABLE,
        "ENSEMBLE_CONFIGS": ENSEMBLE_CONFIGS,
        "SYMMETRY_PAIRS": SYMMETRY_PAIRS,
        "FIT_SLICE_END": FIT_SLICE_END,
        "MULTI_FIT_SLICES": MULTI_FIT_SLICES,
        "SINGLE_FIT_SLICES": SINGLE_FIT_SLICES
    }
    out_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")

if __name__ == "__main__":
    json_target = Path(__file__).parent / "physics_setup.json"
    export_to_json(json_target)
    print(f"[OK] Physics setup exported to JSON: {json_target}")

