"""
src/parselqcdata/susceptibility_pipeline.py
--------------------------------------------------------------------------------
Python 整合层：手征磁化率 (Chiral Susceptibility, \\chi) 端到端分析管道
严格复现 ana/ttest/sucep_calc.py 与 plot_chisce.py 的算法标准：

1. 物理理论定义 (Continuous Definition):
   \\chi_{\\text{disc}} = (V_3 / T) * [ \\langle (\\bar{\\psi}\\psi)^2 \\rangle - \\langle \\bar{\\psi}\\psi \\rangle^2 ]
   其中 V_3 = (Ns * a)^3, T = 1 / (Nt * a), 故 V_3 / T = Ns^3 * Nt * a^4 = V_4。

2. 单构型无偏二次交叉估计器 (Unbiased quadratic estimator on single configuration):
   Obar = (1/k) * \\sum_{i=1}^k O_i
   O2bar = [1/(k(k-1))] * \\sum_{i \\neq j} O_i O_j
         = (1/k) * \\sum_{i=1}^k [ 1/(k-1) * O_i * (\\sum_j O_j - O_i) ]
   数学证明：E[O_i O_j] = [Tr(D^{-1})]^2 + \\sigma_{noise}^2 * \\delta_{ij}。
   对角项消除彻底去除了有限随机源数目 k 引入的虚假噪声方差抬升。

3. 构型级 Jackknife 重采样统计推断:
   a_r = JK(Obar)_r, b_r = JK(O2bar)_r
   \\chi_r = b_r - a_r^2
   mean(\\chi) = (1/N) * \\sum \\chi_r,  err(\\chi) = \\sqrt{N-1} * \\text{std}(\\chi_r, \\text{ddof}=0)

4. 物理标度因子 (Scaling Factors with Ns, Nt, Temperature T):
   - 四维格点体积因子: F_vol = Ns^3 * Nt
     \\chi_{vol} = F_vol * \\chi_{unscaled}
   - 连续温度标度因子: F_scaled = Ns^3 * Nt^3 * T^2 = F_vol * (Nt * T)^2
     \\chi_{scaled} = F_scaled * \\chi_{unscaled}  (单位: MeV^2)
     对应 plot_chisce.py 中的 (16T)^2 * a^2 * \\chi 标准。

5. 跨机器移植性:
   既支持在包含全量原始构型向量的集群计算机上全自动遍历提取，
   又支持在本地轻量环境中调用基准回归数据进行快速验证与画图。
--------------------------------------------------------------------------------
"""

import argparse
import os
from pathlib import Path
import re
import sys
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import (
    ANA_ROOT,
    DEFAULT_READIN_DIR,
    OUTPUT_CONDENSATE_DIR,
    TEMP_MAP,
)

# 8组温度与 Beta 标准映射表 (Nt=16基准, 单位 MeV)
SUSCEPTIBILITY_TEMP_MAP: Dict[str, float] = {
    "4.13": 138.0,
    "4.15": 145.0,
    "4.17": 153.0,
    "4.18": 157.0,
    "4.20": 164.5,
    "4.23": 176.0,
    "4.30": 202.0,
    "4.405": 241.6,
}


def parse_ensemble_meta_from_dir(dir_name: str) -> Dict[str, Union[int, float, str]]:
    """
    从目录名称中解析格点几何尺寸 (Ns, Nt) 与耦合常数 Beta。
    支持格式:
      - L48T16beta4.18ms0.037265m0.001022
      - 48x16b4.18
      - L32T12beta4.17
      - L40T16_beta4.17
    """
    meta: Dict[str, Union[int, float, str]] = {
        "ns": 48,
        "nt": 16,
        "beta": "4.17",
    }

    # 1. 匹配时空几何 Ns, Nt
    geom_match = re.search(r"L(\d+)T(\d+)", dir_name, re.IGNORECASE)
    if geom_match:
        meta["ns"] = int(geom_match.group(1))
        meta["nt"] = int(geom_match.group(2))
    else:
        geom_x = re.search(r"(\d+)x(\d+)", dir_name)
        if geom_x:
            meta["ns"] = int(geom_x.group(1))
            meta["nt"] = int(geom_x.group(2))

    # 2. 匹配 Beta
    beta_match = re.search(r"b(?:eta)?([0-9.]+)", dir_name, re.IGNORECASE)
    if beta_match:
        meta["beta"] = beta_match.group(1)

    return meta


def compute_scaling_factors(ns: int, nt: int, temp_mev: float) -> Tuple[float, float]:
    """
    计算磁化率所需的体积因子与连续标度因子：
    - factor_vol = Ns^3 * Nt
    - factor_scaled = Ns^3 * Nt^3 * T^2 = factor_vol * (Nt * T)^2
    """
    f_vol = float((ns**3) * nt)
    f_scaled = float((ns**3) * (nt**3) * (temp_mev**2))
    return f_vol, f_scaled


def extract_pbp_from_xml(file_path: Path) -> Optional[float]:
    """
    使用高效正则从 XML 中提取 <pbp>(real, imag)</pbp> 的实部数值。
    """
    try:
        content = file_path.read_text(encoding="utf-8", errors="ignore")
        match = re.search(r"<pbp>\(([^,]+),", content)
        if match:
            return float(match.group(1))
    except Exception:
        pass
    return None


def compute_unbiased_quadratic(vals: List[float]) -> Tuple[float, float]:
    """
    计算单个规范构型上的无偏手征凝聚均值与两点方均关联估计：
    - Obar = (1/k) * \\sum_{i=1}^k O_i
    - O2bar = [1/(k(k-1))] * \\sum_{i \\neq j} O_i O_j
            = (1/k) * \\sum_{i=1}^k [ 1/(k-1) * O_i * (S - O_i) ]
    消除了由于有限随机源数目 k 带来的随机噪声方差对算符平方的内积偏差。
    """
    k = len(vals)
    if k <= 1:
        raise ValueError(f"计算两体无偏关联至少需要 2 个随机源向量，当前仅提供 {k} 个")

    total_s = float(sum(vals))
    obar = total_s / k
    # 向量化无偏方均估计
    obar_sq_unbiased = float(np.mean([(1.0 / (k - 1)) * x * (total_s - x) for x in vals]))
    return obar, obar_sq_unbiased


def jackknife_resample(arr: np.ndarray) -> np.ndarray:
    """
    对一维样本序列执行 Jackknife Leave-One-Out 重采样。
    输入大小 N，输出大小 N 的 Jackknife 假样本序列。
    """
    n = len(arr)
    if n <= 1:
        raise ValueError("Jackknife 重采样至少需要 2 个构型样本")
    total_sum = np.sum(arr)
    return (total_sum - arr) / (n - 1)


def compute_jackknife_susceptibility(
    list_obar: List[float],
    list_o2bar: List[float],
    ns: int = 48,
    nt: int = 16,
    temp_mev: float = 157.0,
) -> Dict[str, float]:
    """
    通过 Jackknife 假样本计算手征磁化率 \\chi = \\langle O^2 \\rangle - \\langle O \\rangle^2
    并根据 Ns, Nt 与温度 T 准确折算标度因子：
    - factor_vol = Ns^3 * Nt
    - factor_scaled = Ns^3 * Nt^3 * T^2
    """
    n = len(list_obar)
    if n <= 1 or len(list_o2bar) != n:
        raise ValueError(f"样本维度不合法: obar={len(list_obar)}, o2bar={len(list_o2bar)}")

    arr_obar = np.array(list_obar, dtype=np.float64)
    arr_o2bar = np.array(list_o2bar, dtype=np.float64)

    a_jk = jackknife_resample(arr_obar)
    b_jk = jackknife_resample(arr_o2bar)

    # 每一个 Jackknife bin 中的磁化率估计
    chi_jk = b_jk - (a_jk**2)

    mean_unscaled = float(np.mean(chi_jk))
    # 统计标准误: \\sigma = \\sqrt{N-1} * \\text{std}(chi_jk, \\text{ddof}=0)
    err_unscaled = float(np.sqrt(n - 1) * np.std(chi_jk, ddof=0))

    f_vol, f_scaled = compute_scaling_factors(ns, nt, temp_mev)

    return {
        "mean_unscaled": mean_unscaled,
        "error_unscaled": err_unscaled,
        "factor_vol": f_vol,
        "mean_vol_scaled": mean_unscaled * f_vol,
        "error_vol_scaled": err_unscaled * f_vol,
        "factor_scaled": f_scaled,
        "mean_scaled": mean_unscaled * f_scaled,
        "error_scaled": err_unscaled * f_scaled,
        "num_configs": n,
        "ns": ns,
        "nt": nt,
        "temp": temp_mev,
    }


class SusceptibilityPipeline:
    """手征磁化率端到端处理与统计提取管道（直接使用光夸克随机源向量，不进行向量级 RM 减除）"""

    def __init__(self, ana_root: Optional[Path] = None, output_dir: Optional[Path] = None):
        self.ana_root = Path(ana_root) if ana_root else ANA_ROOT
        self.output_dir = Path(output_dir) if output_dir else OUTPUT_CONDENSATE_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def extract_from_meas_directory(
        self,
        ensemble_dir: Path,
        ns: Optional[int] = None,
        nt: Optional[int] = None,
        temp_mev: Optional[float] = None,
    ) -> Dict[str, float]:
        """
        从包含 meas.* 文件夹的实际目录中遍历解析轻夸克随机源 XML 文件，
        提取单构型量与全系综 Jackknife 磁化率结果。
        自动结合 Ns, Nt, Temperature 计算相应标度因子。
        """
        meas_dirs = sorted(ensemble_dir.glob("meas.*"))
        if not meas_dirs:
            sub_cand = ensemble_dir / "test_condensate"
            if sub_cand.exists():
                meas_dirs = sorted(sub_cand.glob("meas.*"))

        if not meas_dirs:
            raise FileNotFoundError(f"在目录 {ensemble_dir} 中未找到任何 meas.* 测量子目录")

        meta = parse_ensemble_meta_from_dir(ensemble_dir.name)
        eff_ns = int(ns if ns is not None else meta["ns"])
        eff_nt = int(nt if nt is not None else meta["nt"])
        eff_beta = str(meta["beta"])

        if temp_mev is not None:
            eff_temp = float(temp_mev)
        else:
            base_temp = SUSCEPTIBILITY_TEMP_MAP.get(eff_beta, 157.0)
            eff_temp = base_temp * (16.0 / eff_nt)

        list_obar: List[float] = []
        list_o2bar: List[float] = []

        for m_dir in meas_dirs:
            # 读取光夸克随机源向量
            vals_light = []
            for f in sorted(m_dir.glob("PsibarPsi/Z2Stochastic_pbp_l_vec*.xml")):
                val = extract_pbp_from_xml(f)
                if val is not None:
                    vals_light.append(val)

            if len(vals_light) > 1:
                obar, o2bar = compute_unbiased_quadratic(vals_light)
                list_obar.append(obar)
                list_o2bar.append(o2bar)

        if not list_obar:
            raise ValueError(f"在目录 {ensemble_dir} 中未能成功提取到有效随机源向量数据")

        return compute_jackknife_susceptibility(
            list_obar, list_o2bar, ns=eff_ns, nt=eff_nt, temp_mev=eff_temp
        )

    def run_cluster_scan(self, readin_dir: Path) -> pl.DataFrame:
        """
        在拥有全量构型测量数据的计算机/集群上执行全温度扫描计算。
        自动遍历 readin_dir 下的所有 L48T16beta* 或 48x16b* 目录，
        计算手征磁化率，导出至 output/condensate/results_susceptibility.csv。
        """
        betas = sorted(list(SUSCEPTIBILITY_TEMP_MAP.keys()), key=float)
        results = []

        print(f"[SusceptibilityPipeline] 正在从数据目录 {readin_dir} 执行集群全量扫描...")

        for beta in betas:
            patterns = [
                f"L48T16beta{beta}*",
                f"48x16b{beta}*",
                f"L*beta{beta}*",
            ]
            matched = []
            for pat in patterns:
                matched.extend(list(readin_dir.glob(pat)))
            matched = [d for d in set(matched) if d.is_dir()]

            if not matched:
                print(f"  [WARN] 未找到 beta={beta} 的数据目录")
                continue

            target_dir = matched[0]
            temp = SUSCEPTIBILITY_TEMP_MAP[beta]
            print(f"  -> 处理 beta={beta} (T={temp} MeV) 目录: {target_dir.name}")

            res = self.extract_from_meas_directory(target_dir, temp_mev=temp)
            results.append({
                "Beta": float(beta),
                "Temp": temp,
                "Ns": res["ns"],
                "Nt": res["nt"],
                "Mean_unscaled": res["mean_unscaled"],
                "Error_unscaled": res["error_unscaled"],
                "Mean_vol_scaled": res["mean_vol_scaled"],
                "Error_vol_scaled": res["error_vol_scaled"],
                "Mean_scaled": res["mean_scaled"],
                "Error_scaled": res["error_scaled"],
            })

        df_res = pl.DataFrame(results).sort("Temp")

        out_csv = self.output_dir / "results_susceptibility.csv"
        out_parquet = self.output_dir / "results_susceptibility.parquet"

        df_res.write_csv(out_csv)
        df_res.write_parquet(out_parquet)

        print(f"[OK] 集群计算完成，结果已保存至 {out_csv}")
        return df_res

    def get_full_scan_results(
        self,
        force_recompute: bool = False,
    ) -> pl.DataFrame:
        """
        获取全温度扫描 (8 组温度) 的手征磁化率完整数据集。
        若集群基准数据存在，保证与 ana/ttest 完全一致；
        同时精确计算并补充 Ns, Nt, T, Mean_vol_scaled 与 Mean_scaled 列。
        """
        out_csv = self.output_dir / "results_susceptibility.csv"
        out_parquet = self.output_dir / "results_susceptibility.parquet"
        ana_csv = self.ana_root / "ttest" / "results_susceptibility.csv"

        if not force_recompute and out_csv.exists():
            return pl.read_csv(out_csv)

        # 优先读取基准数据
        if ana_csv.exists():
            df_norm = pl.read_csv(ana_csv)
        else:
            raise FileNotFoundError(f"未找到基准数据: {ana_csv}")

        # 标准化添加 Ns=48, Nt=16, 体积因子与连续标度因子
        rows = []
        for r in df_norm.iter_rows(named=True):
            b = r["Beta"]
            t = float(r["Temp"])
            mu = float(r["Mean_unscaled"])
            eu = float(r["Error_unscaled"])
            f_vol, f_scaled = compute_scaling_factors(48, 16, t)
            rows.append({
                "Beta": b,
                "Temp": t,
                "Ns": 48,
                "Nt": 16,
                "Mean_unscaled": mu,
                "Error_unscaled": eu,
                "Mean_vol_scaled": mu * f_vol,
                "Error_vol_scaled": eu * f_vol,
                "Mean_scaled": mu * f_scaled,
                "Error_scaled": eu * f_scaled,
            })

        df_enriched = pl.DataFrame(rows).sort("Temp")
        df_enriched.write_csv(out_csv)
        df_enriched.write_parquet(out_parquet)

        print(f"[SusceptibilityPipeline] 成功同步全量磁化率数据集至 {self.output_dir}")
        return df_enriched


def main():
    parser = argparse.ArgumentParser(description="手征磁化率数据抽取与统计分析管道")
    parser.add_argument("--readin-dir", type=str, default=None, help="包含原始格点构型 meas.* 的数据根目录")
    parser.add_argument("--output-dir", type=str, default=str(OUTPUT_CONDENSATE_DIR), help="输出结果目录")
    args = parser.parse_args()

    pipeline = SusceptibilityPipeline(output_dir=Path(args.output_dir))

    if args.readin_dir and Path(args.readin_dir).exists():
        print(f"正在集群模式下扫描: {args.readin_dir}")
        pipeline.run_cluster_scan(Path(args.readin_dir))
    else:
        df_res = pipeline.get_full_scan_results(force_recompute=True)
        print("=== Chiral Susceptibility (8 温度扫描) ===")
        print(df_res)


if __name__ == "__main__":
    main()
