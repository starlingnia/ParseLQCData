"""
src/parselqcdata/condensate_pipeline.py
--------------------------------------------------------------------------------
Python 整合层：手征凝聚 (Chiral Condensate) 端到端分析管道
- 导入 docs/physics_setup 中的 ENSEMBLE_CONFIGS 与全局配置
- 自动从目录名称解析格点规模 (Ns, Nt) 与夸克质量 (ml, ms, mres, Zm)
- 驱动 C++ 核心库 tools/condensate_orchestrator 进行高效并发抽取与 Jackknife 重采样
--------------------------------------------------------------------------------
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Union

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    ENSEMBLE_CONFIGS,
    parse_ensemble_dirname,
    resolve_dataset_dir,
    DEFAULT_READIN_DIR,
    calculate_residual_mass,
    MRES_TABLE,
)
from tools.condensate_orchestrator import CondensateOrchestrator


class CondensatePipeline:
    """手征凝聚数据处理与统计分析管道"""

    def __init__(self, orchestrator: Optional[CondensateOrchestrator] = None):
        self.orchestrator = orchestrator or CondensateOrchestrator()

    def process_ensemble(
        self,
        target_or_dir: Union[str, Path],
        base_readin_dir: Optional[Union[str, Path]] = None,
        ml: Optional[float] = None,
        ms: Optional[float] = None,
        mres: Optional[float] = None,
        zm: Optional[float] = None,
    ) -> Dict[str, Union[float, int, str, Path]]:
        """
        处理单个手征凝聚数据集，支持自动解析目录名中编码的格点与物理参数：
        例如: L32T12_beta4.17ms0.040m0.0020 或 L48T16beta4.13ms0.043547m0.000805
        """
        if base_readin_dir is None:
            base_readin_dir = DEFAULT_READIN_DIR
        base_readin = Path(base_readin_dir)

        target_path = Path(target_or_dir)
        if not target_path.exists():
            target_path = resolve_dataset_dir(base_readin, str(target_or_dir))

        if not target_path.exists():
            raise FileNotFoundError(f"找不到手征凝聚数据集目录: {target_or_dir} (基准目录: {base_readin})")

        # 智能匹配 test_condensate 子目录
        run_path = target_path / "test_condensate" if (target_path / "test_condensate").exists() else target_path

        # 解析物理参数 (优先使用用户显式传入参数，其次自 docs/ 字典或目录名解析)
        meta = parse_ensemble_dirname(target_path.name)

        eff_ml = ml if ml is not None else float(meta.get("ml", 0.001001))
        eff_ms = ms if ms is not None else float(meta.get("ms", 0.0384))
        eff_beta = float(meta.get("beta", 4.17))
        eff_mres = mres if mres is not None else float(meta.get("mres", calculate_residual_mass(eff_beta)))
        eff_zm = zm if zm is not None else float(meta.get("zm", 0.966247))

        # 调用 C++ 底层流水线
        mean, error, n_cfgs = self.orchestrator.process_condensate(
            base_dir=str(run_path),
            m_light=eff_ml,
            m_strange=eff_ms,
            m_residual=eff_mres,
            zm_factor=eff_zm,
        )

        return {
            "dataset_name": target_path.name,
            "dataset_path": str(run_path),
            "beta": eff_beta,
            "ns": meta.get("ns", 0),
            "nt": meta.get("nt", 0),
            "temperature": float(meta.get("temperature", 153.31)),
            "ml": eff_ml,
            "ms": eff_ms,
            "mres": eff_mres,
            "zm": eff_zm,
            "pbp_sub_mean": mean,
            "pbp_sub_error": error,
            "num_cfgs": n_cfgs,
        }

    def process_all_ensembles(
        self,
        base_readin_dir: Optional[Union[str, Path]] = None,
    ) -> List[Dict[str, Union[float, int, str, Path]]]:
        """
        批量扫描并运行基准目录下所有可识别的手征凝聚数据集
        """
        if base_readin_dir is None:
            base_readin_dir = DEFAULT_READIN_DIR
        base_readin = Path(base_readin_dir)

        results = []
        if not base_readin.exists():
            return results

        for child in sorted(base_readin.iterdir()):
            if not child.is_directory():
                continue
            # 判断是否具有 meas 目录或 test_condensate 目录
            is_condensate_dir = (child / "test_condensate").exists() or any(
                p.is_directory() and p.name.starts_with("meas.") for p in child.iterdir()
            )
            if is_condensate_dir:
                try:
                    res = self.process_ensemble(child, base_readin_dir=base_readin)
                    results.append(res)
                except Exception as e:
                    print(f"[WARN] 数据集 {child.name} 处理跳过: {e}")

        return results
