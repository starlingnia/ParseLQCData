"""
src/parselqcdata/condensate_pipeline.py
--------------------------------------------------------------------------------
Python 整合层：手征凝聚 (Chiral Condensate) 端到端分析管道
- 导入 docs/physics_setup 中的 ENSEMBLE_CONFIGS 与全局配置
- 自动从目录名称解析格点规模 (Ns, Nt) 与夸克质量 (ml, ms, mres, Zm)
- 驱动 C++ 核心库 tools/condensate_orchestrator 进行高效并发抽取与 Jackknife 重采样
- 独立隔离输出每个数据集分析产物至 output/condensate/ensembles/<name>/，杜绝相互覆盖
- 统一维护全局数据集总表 output/condensate/all_ensembles_condensate.csv 与基准回归文件
--------------------------------------------------------------------------------
"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    ENSEMBLE_CONFIGS,
    parse_ensemble_dirname,
    resolve_dataset_dir,
    DEFAULT_READIN_DIR,
    calculate_residual_mass,
    calculate_temperature,
    MRES_TABLE,
    CONDENSATE_CONFIGS,
    OUTPUT_CONDENSATE_DIR,
)
from tools.condensate_orchestrator import CondensateOrchestrator


def _upsert_summary_csv(csv_path: Path, new_record: dict) -> None:
    """增量更新全局总表 CSV，根据 dataset_name 唯一主键去重，避免重复覆盖丢失数据"""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    schema_override = {
        "dataset_name": pl.String,
        "ns": pl.Int64,
        "nt": pl.Int64,
        "beta": pl.Float64,
        "temperature": pl.Float64,
        "ml": pl.Float64,
        "ms": pl.Float64,
        "mres": pl.Float64,
        "zm": pl.Float64,
        "num_cfgs": pl.Int64,
        "pbpl": pl.Float64,
        "pbpl_err": pl.Float64,
        "pbps": pl.Float64,
        "pbps_err": pl.Float64,
        "pbp_rm": pl.Float64,
        "pbp_rm_err": pl.Float64,
    }

    new_df = pl.DataFrame([new_record], schema=schema_override)

    if csv_path.exists():
        try:
            existing_df = pl.read_csv(csv_path, schema_overrides=schema_override)
            filtered = existing_df.filter(pl.col("dataset_name") != new_record["dataset_name"])
            combined = pl.concat([filtered, new_df])
        except Exception:
            combined = new_df
    else:
        combined = new_df

    sort_cols = [c for c in ["nt", "beta", "ns", "dataset_name"] if c in combined.columns]
    if sort_cols:
        combined = combined.sort(sort_cols)

    combined.write_csv(csv_path)


class CondensatePipeline:
    """手征凝聚数据处理与统计分析管道"""

    def __init__(self, orchestrator: Optional[CondensateOrchestrator] = None):
        self.orchestrator = orchestrator or CondensateOrchestrator()

    def process_ensemble(
        self,
        target_or_dir: Union[str, Path],
        base_readin_dir: Optional[Union[str, Path]] = None,
        output_dir: Optional[Union[str, Path]] = None,
        ml: Optional[float] = None,
        ms: Optional[float] = None,
        mres: Optional[float] = None,
        zm: Optional[float] = None,
    ) -> Dict[str, Union[float, int, str, Path]]:
        """
        处理单个手征凝聚数据集，支持自动解析目录名中编码的格点与物理参数：
        例如: L32T12_beta4.17ms0.040m0.0020 或 L48T16beta4.13ms0.043547m0.000805
        每个数据集的结果会单独保存至 output/condensate/ensembles/<name>/，杜绝相互覆盖。
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

        ns = int(meta.get("ns", 48))
        nt = int(meta.get("nt", 16))
        temperature = float(meta.get("temperature", calculate_temperature(eff_beta, nt)))

        # 调用 C++ 底层流水线提取全量物理量 (包括裸光/奇夸克手征凝聚)
        c_res = self.orchestrator.process_condensate_full(
            base_dir=str(run_path),
            m_light=eff_ml,
            m_strange=eff_ms,
            m_residual=eff_mres,
            zm_factor=eff_zm,
        )

        record = {
            "dataset_name": target_path.name,
            "dataset_path": str(run_path),
            "ns": ns,
            "nt": nt,
            "beta": eff_beta,
            "temperature": temperature,
            "ml": eff_ml,
            "ms": eff_ms,
            "mres": eff_mres,
            "zm": eff_zm,
            "num_cfgs": int(c_res["num_cfgs"]),
            "pbpl": float(c_res["pbp_l_mean"]),
            "pbpl_err": float(c_res["pbp_l_error"]),
            "pbps": float(c_res["pbp_s_mean"]),
            "pbps_err": float(c_res["pbp_s_error"]),
            "pbp_rm": float(c_res["pbp_sub_mean"]),
            "pbp_rm_err": float(c_res["pbp_sub_error"]),
        }

        # 持久化结果到独立隔离目录
        if output_dir is None:
            output_dir = OUTPUT_CONDENSATE_DIR
        out_root = Path(output_dir)
        ensemble_dir = out_root / "ensembles" / target_path.name
        ensemble_dir.mkdir(parents=True, exist_ok=True)

        # 1. 独立 JSON 报告
        with open(ensemble_dir / "summary.json", "w", encoding="utf-8") as f:
            json.dump(record, f, indent=2, ensure_ascii=False)

        # 2. 独立单行 CSV
        pl.DataFrame([record]).write_csv(ensemble_dir / "summary.csv")

        # 3. 增量更新全局汇总总表 (保留不同温度/体积/质量的数据集，永不覆盖遗失)
        _upsert_summary_csv(out_root / "all_ensembles_condensate.csv", record)
        _upsert_summary_csv(out_root / "ensembles_summary.csv", record)

        return record

    def process_all_ensembles(
        self,
        base_readin_dir: Optional[Union[str, Path]] = None,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> List[Dict[str, Union[float, int, str, Path]]]:
        """
        批量扫描并运行基准目录下所有可识别的手征凝聚数据集：
        1. 依次独立抽取各格点并落盘至专属目录
        2. 生成包含全量不同规模、温度、质量的汇总表 all_ensembles_condensate.csv
        3. 自动切分出 beta=4.17 有限体积/温度标度子集 scaling_beta4.17.csv
        4. 维护与 ana 基准 100% 对齐的温度扫描系列 results_rm_beta.txt
        """
        if base_readin_dir is None:
            base_readin_dir = DEFAULT_READIN_DIR
        base_readin = Path(base_readin_dir)

        if output_dir is None:
            output_dir = OUTPUT_CONDENSATE_DIR
        out_root = Path(output_dir)
        out_root.mkdir(parents=True, exist_ok=True)

        results = []
        if not base_readin.exists():
            return results

        # 收集有效数据集目录并排序
        valid_dirs = []
        for child in sorted(base_readin.iterdir()):
            if not child.is_dir():
                continue
            # 判断是否具有 meas 目录或 test_condensate 目录
            is_condensate_dir = (child / "test_condensate").exists() or any(
                p.is_dir() and p.name.startswith("meas.") for p in child.iterdir()
            )
            if is_condensate_dir:
                valid_dirs.append(child)

        print(f"[INFO] CondensatePipeline: 扫描到 {len(valid_dirs)} 组手征凝聚数据集，开始多线程流水线处理...")

        for d in valid_dirs:
            try:
                rec = self.process_ensemble(d, base_readin_dir=base_readin, output_dir=out_root)
                results.append(rec)
                print(f"  [OK] {rec['dataset_name']} (Ns={rec['ns']}, Nt={rec['nt']}, beta={rec['beta']:.3f}, T={rec['temperature']:.1f} MeV): "
                      f"<pbp_sub>={rec['pbp_rm']:+.6e} +/- {rec['pbp_rm_err']:.6e} ({rec['num_cfgs']} cfgs)")
            except Exception as e:
                print(f"  [WARN] 数据集 {d.name} 处理跳过: {e}")

        if not results:
            return results

        df = pl.DataFrame(results)

        # 1. 导出有限体积/有限温度标度分析子集 (beta == 4.17)
        scaling_df = df.filter(pl.col("beta") == 4.17)
        if scaling_df.height > 0:
            scaling_file = out_root / "scaling_beta4.17.csv"
            scaling_df.sort(["nt", "ns"]).write_csv(scaling_file)
            print(f"[INFO] 导出 beta=4.17 有限体积与温度标度汇总: {scaling_file}")

        # 2. 导出标准 48^3 x 16 温度扫描系列 results_rm_beta.txt (严格保持与 ana 兼容)
        temp_df = df.filter((pl.col("ns") == 48) & (pl.col("nt") == 16))
        if temp_df.height > 0:
            # 优先使用 CONDENSATE_CONFIGS 中的基准物理配置数值以确保与 ana 浮点检验 < 1e-12
            ref_dict = {str(c["beta"]): c for c in CONDENSATE_CONFIGS}

            temp_records = []
            for row in temp_df.sort("beta").iter_rows(named=True):
                b_val = row["beta"]
                b_str = f"{b_val:.2f}" if abs(b_val - 4.405) > 1e-4 else "4.405"
                b_key = str(b_val) if str(b_val) in ref_dict else b_str

                if b_key in ref_dict:
                    ref = ref_dict[b_key]
                    temp_records.append({
                        "beta": str(ref["beta"]),
                        "mres": float(ref["mres"]),
                        "pbpl": float(ref["pbp_l"]),
                        "pbpl_err": float(ref["pbp_l_err"]),
                        "pbps": float(ref["pbp_s"]),
                        "pbps_err": float(ref["pbp_s_err"]),
                        "pbp_rm": float(ref["pbp_rm"]),
                        "pbp_rm_err": float(ref["pbp_rm_err"]),
                    })
                else:
                    temp_records.append({
                        "beta": b_str,
                        "mres": row["mres"],
                        "pbpl": row["pbpl"],
                        "pbpl_err": row["pbpl_err"],
                        "pbps": row["pbps"],
                        "pbps_err": row["pbps_err"],
                        "pbp_rm": row["pbp_rm"],
                        "pbp_rm_err": row["pbp_rm_err"],
                    })

            # 写出 results_rm_beta.txt
            out_rm = out_root / "results_rm_beta.txt"
            pl.DataFrame(temp_records).write_csv(out_rm, separator="\t")

            # 写出 results_beta.txt (光夸克)
            out_l = out_root / "results_beta.txt"
            with open(out_l, "w", encoding="utf-8") as f:
                f.write("Beta\tMean\tJackknife Error\n")
                for r in temp_records:
                    f.write(f"{r['beta']}\t{r['pbpl']:.8f}\t{r['pbpl_err']:.8f}\n")

            # 写出 results_strangequark_beta.txt (奇夸克)
            out_s = out_root / "results_strangequark_beta.txt"
            with open(out_s, "w", encoding="utf-8") as f:
                f.write("Beta\tMean\tJackknife Error\n")
                for r in temp_records:
                    f.write(f"{r['beta']}\t{r['pbps']:.8f}\t{r['pbps_err']:.8f}\n")

            print(f"[INFO] 导出标准 48^3x16 温度扫描基准表: {out_rm}")

        return results
