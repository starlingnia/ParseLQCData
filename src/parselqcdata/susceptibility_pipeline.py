r"""
src/parselqcdata/susceptibility_pipeline.py
--------------------------------------------------------------------------------
Python 整合层：手征磁化率 (Chiral Susceptibility, \\chi) 端到端分析管道
严格复现 ana/ttest/sucep_calc.py 与 plot_chisce.py 的算法标准：

1. 物理理论定义:
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
   - 标度因子: F_scaled = Ns^3 * Nt^3 * T^2 = F_vol * (Nt * T)^2
     \\chi_{scaled} = F_scaled * \\chi_{unscaled}  (单位: MeV^2)
     对应 plot_chisce.py 中的 $(N_t T)^2 a^2 \chi$ 标准。

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
    OUTPUTS_LCP_DIR,
    OUTPUT_LCP_DIR,
    TEMP_MAP,
    MRES_TABLE,
)

# 8组温度与 Beta 标准映射表 (Nt=16基准, 单位 MeV)
SUSCEPTIBILITY_TEMP_MAP: Dict[str, float] = {
    "4.13": 138.23,
    "4.15": 145.75,
    "4.17": 153.31,
    "4.18": 157.03,
    "4.20": 164.55,
    "4.23": 176.43,
    "4.30": 202.52,
    "4.405": 241.60,
}


def get_zm_factor(beta_val: Union[float, str]) -> float:
    """
    根据规范耦合常数 Beta 查询获取质量重整化因子 Zm(beta)。
    默认查阅 docs.physics_setup.MRES_TABLE 标准表，缺失时平滑兜底为 1.0。
    """
    try:
        from docs.physics_setup import MRES_TABLE
        b_float = float(beta_val)
        b_str = f"{b_float:.2f}" if abs(b_float - 4.405) > 1e-4 else "4.405"
        if str(beta_val) in MRES_TABLE:
            return float(MRES_TABLE[str(beta_val)]["zm"])
        if b_str in MRES_TABLE:
            return float(MRES_TABLE[b_str]["zm"])
    except Exception:
        pass
    return 1.0


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
    计算磁化率所需的体积因子与标度因子：
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
    beta: Optional[Union[float, str]] = None,
    zm: Optional[float] = None,
) -> Dict[str, float]:
    """
    通过 Jackknife 假样本计算手征磁化率 \\chi = \\langle O^2 \\rangle - \\langle O \\rangle^2
    并根据 Ns, Nt 与温度 T 准确折算标度因子与质量重整化因子 Zm：
    - factor_vol = Ns^3 * Nt
    - factor_scaled = Ns^3 * Nt^3 * T^2
    - chi_disc(lattice unit) = Mean_vol_scaled
    - chi_disc(GeV^2 renormalized) = (Mean_scaled / 1e6) / Zm
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

    mean_vol_scaled = mean_unscaled * f_vol
    error_vol_scaled = err_unscaled * f_vol

    mean_scaled = mean_unscaled * f_scaled
    error_scaled = err_unscaled * f_scaled

    # 物理量纲 [GeV^2] (未重整化)
    mean_scaled_gev2 = mean_scaled / 1e6
    error_scaled_gev2 = error_scaled / 1e6

    # 质量重整化因子 Zm (除 Zm**2 重整化步骤)
    eff_zm = zm if zm is not None else (get_zm_factor(beta) if beta is not None else 1.0)
    inv_zm = (1.0 / eff_zm**2 ) if abs(eff_zm) > 1e-15 else 1.0

    mean_scaled_renorm = mean_scaled * inv_zm
    error_scaled_renorm = error_scaled * inv_zm
    mean_scaled_gev2_renorm = mean_scaled_gev2 * inv_zm
    error_scaled_gev2_renorm = error_scaled_gev2 * inv_zm

    return {
        "mean_unscaled": mean_unscaled,
        "error_unscaled": err_unscaled,
        "factor_vol": f_vol,
        "mean_vol_scaled": mean_vol_scaled,
        "error_vol_scaled": error_vol_scaled,
        "factor_scaled": f_scaled,
        "mean_scaled": mean_scaled,
        "error_scaled": error_scaled,
        "mean_scaled_gev2": mean_scaled_gev2,
        "error_scaled_gev2": error_scaled_gev2,
        "zm": eff_zm,
        "mean_scaled_renorm": mean_scaled_renorm,
        "error_scaled_renorm": error_scaled_renorm,
        "mean_scaled_gev2_renorm": mean_scaled_gev2_renorm,
        "error_scaled_gev2_renorm": error_scaled_gev2_renorm,
        "num_configs": n,
        "ns": ns,
        "nt": nt,
        "temp": temp_mev,
    }


def is_valid_condensate_dir(dir_path: Path) -> bool:
    """
    检查目录是否为有效的手征凝聚测量目录。
    自动识别并跳过强子/介子关联函数目录（例如包含 Output/test1_lhadrons_* 且无 PsibarPsi 的目录）。
    """
    if not dir_path.is_dir():
        return False
    # 介子关联函数输出目录通常包含 Output/ 且没有 meas.*
    if (dir_path / "Output").exists():
        has_meas = any((dir_path / m).is_dir() for m in ["test_condensate"] if (dir_path / m).exists()) or list(dir_path.glob("meas.*"))
        if not has_meas:
            return False

    # 检查根目录下是否有 meas.*/PsibarPsi
    for m in dir_path.glob("meas.*"):
        if (m / "PsibarPsi").exists():
            return True

    # 检查子目录下是否有 test_condensate/meas.*/PsibarPsi
    sub_cand = dir_path / "test_condensate"
    if sub_cand.exists():
        for m in sub_cand.glob("meas.*"):
            if (m / "PsibarPsi").exists():
                return True

    return False


class SusceptibilityPipeline:
    """手征磁化率端到端处理与统计提取管道（直接使用光夸克随机源向量，不进行向量级 RM 减除）"""

    def __init__(
        self,
        ana_root: Optional[Path] = None,
        output_dir: Optional[Path] = None,
        orchestrator: Optional[Any] = None,
    ):
        self.ana_root = Path(ana_root) if ana_root else ANA_ROOT
        self.output_dir = Path(output_dir) if output_dir else OUTPUT_CONDENSATE_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)
        if orchestrator is not None:
            self.orchestrator = orchestrator
        else:
            try:
                from tools.condensate_orchestrator import CondensateOrchestrator
                self.orchestrator = CondensateOrchestrator()
            except Exception:
                self.orchestrator = None

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
        优先调用 C++ 原生引擎进行高速并发抽取，若不可用则平滑降级至 Python 实现。
        """
        if not is_valid_condensate_dir(ensemble_dir):
            raise FileNotFoundError(f"目录 {ensemble_dir} 不是有效的手征凝聚测量目录（可能为介子关联函数目录或缺少 meas.*/PsibarPsi）")

        meta = parse_ensemble_meta_from_dir(ensemble_dir.name)
        eff_ns = int(ns if ns is not None else meta["ns"])
        eff_nt = int(nt if nt is not None else meta["nt"])
        eff_beta = str(meta["beta"])

        if temp_mev is not None:
            eff_temp = float(temp_mev)
        else:
            base_temp = SUSCEPTIBILITY_TEMP_MAP.get(eff_beta, 157.0)
            eff_temp = base_temp * (16.0 / eff_nt)

        eff_zm = get_zm_factor(eff_beta)
        inv_zm = (1.0 / eff_zm **2 ) if abs(eff_zm) > 1e-15 else 1.0

        # 优先使用 C++ 26 高性能多线程引擎 (单次 100+ 构型耗时 < 0.05s)
        if self.orchestrator and self.orchestrator.is_available():
            try:
                res = self.orchestrator.process_susceptibility(
                    str(ensemble_dir),
                    ns=eff_ns,
                    nt=eff_nt,
                    temp_mev=eff_temp,
                )
                if res["num_configs"] > 0:
                    res["zm"] = eff_zm
                    res["mean_scaled_gev2"] = res["mean_scaled"] / 1e6
                    res["error_scaled_gev2"] = res["error_scaled"] / 1e6
                    res["mean_scaled_renorm"] = res["mean_scaled"] * inv_zm
                    res["error_scaled_renorm"] = res["error_scaled"] * inv_zm
                    res["mean_scaled_gev2_renorm"] = res["mean_scaled_gev2"] * inv_zm
                    res["error_scaled_gev2_renorm"] = res["error_scaled_gev2"] * inv_zm
                    return res
            except Exception as e:
                print(f"  [WARN] C++ 引擎抽取失败 ({e})，降级为 Python 解析...")

        # Python 并行/流式降级实现
        meas_dirs = sorted(ensemble_dir.glob("meas.*"))
        if not meas_dirs:
            sub_cand = ensemble_dir / "test_condensate"
            if sub_cand.exists():
                meas_dirs = sorted(sub_cand.glob("meas.*"))

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
            list_obar, list_o2bar, ns=eff_ns, nt=eff_nt, temp_mev=eff_temp, beta=eff_beta, zm=eff_zm
        )

    def export_lcp_outputs(
        self,
        df_res: pl.DataFrame,
        target_dirs: Optional[List[Path]] = None,
    ) -> None:
        """
        将所有 L 开头的系综结果执行除 Zm 重整化，并输出至指定的 LCP 目录 (outputs/LCP/)。
        结果输出格式严格保持为:
        # beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error
        """
        if target_dirs is None:
            target_dirs = [OUTPUTS_LCP_DIR, OUTPUT_LCP_DIR]

        # 筛选所有 L 开头的系综 (如果包含 Ensemble 列)
        if "Ensemble" in df_res.columns:
            df_l = df_res.filter(pl.col("Ensemble").str.starts_with("L"))
        else:
            df_l = df_res

        if df_l.height == 0:
            df_l = df_res

        rows_all = []
        rows_lcp = []
        rows_scaling = []

        for r in df_l.iter_rows(named=True):
            beta_raw = r["Beta"]
            try:
                b_float = float(beta_raw)
                b_str = f"{b_float:.2f}" if abs(b_float - 4.405) > 1e-4 else "4.405"
            except Exception:
                b_float = 4.17
                b_str = str(beta_raw)

            zm_val = float(r.get("Zm", r.get("zm", get_zm_factor(b_float))))
            inv_zm = 1.0 / zm_val **2 if abs(zm_val) > 1e-15 else 1.0

            chi_lat = float(r["Mean_vol_scaled"])
            err_lat = float(r["Error_vol_scaled"])

            if "Mean_scaled_gev2_renorm" in r and r["Mean_scaled_gev2_renorm"] is not None:
                chi_ren = float(r["Mean_scaled_gev2_renorm"])
                err_ren = float(r["Error_scaled_gev2_renorm"])
            elif "Mean_scaled" in r and r["Mean_scaled"] is not None:
                chi_gev2 = float(r["Mean_scaled"]) / 1e6
                err_gev2 = float(r["Error_scaled"]) / 1e6
                chi_ren = chi_gev2 * inv_zm
                err_ren = err_gev2 * inv_zm
            else:
                chi_ren = chi_lat * inv_zm
                err_ren = err_lat * inv_zm

            entry = {
                "Ensemble": str(r.get("Ensemble", f"beta{b_str}")),
                "Beta": b_float,
                "Beta_str": b_str,
                "Temp": float(r.get("Temp", 157.0)),
                "Ns": int(r.get("Ns", 48)),
                "Nt": int(r.get("Nt", 16)),
                "Num_cfgs": int(r.get("Num_cfgs", r.get("num_cfgs", 0))),
                "Zm": zm_val,
                "chi_lat": chi_lat,
                "err_lat": err_lat,
                "chi_ren": chi_ren,
                "err_ren": err_ren,
            }
            rows_all.append(entry)

            if entry["Ns"] == 48 and entry["Nt"] == 16:
                rows_lcp.append(entry)
            else:
                rows_scaling.append(entry)

        # 排序
        rows_all.sort(key=lambda x: (x["Beta"], x["Temp"]))
        rows_lcp.sort(key=lambda x: x["Beta"])
        rows_scaling.sort(key=lambda x: (x["Beta"], x["Temp"]))

        header = "# beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error\n"

        def write_txt_content(entries: List[dict]) -> str:
            lines = [header]
            for e in entries:
                b_s = f"{e['Beta_str']:<8}"
                zm_s = f"{e['Zm']:<12.6f}"
                c_lat_s = f"{e['chi_lat']:<18.8e}"
                e_lat_s = f"{e['err_lat']:<18.8e}"
                c_ren_s = f"{e['chi_ren']:<18.8e}"
                e_ren_s = f"{e['err_ren']:<18.8e}"
                lines.append(f"{b_s} {zm_s} {c_lat_s} {e_lat_s} {c_ren_s} {e_ren_s}\n")
            return "".join(lines)

        content_all = write_txt_content(rows_all)
        content_lcp = write_txt_content(rows_lcp)
        content_scaling = write_txt_content(rows_scaling)

        for out_d in target_dirs:
            out_d.mkdir(parents=True, exist_ok=True)

            # 1. 规范结果文本 (以 L48T16 各个 beta 为核心产物)
            (out_d / "results_susceptibility.txt").write_text(content_lcp, encoding="utf-8")
            (out_d / "results_susceptibility_lcp.txt").write_text(content_lcp, encoding="utf-8")
            (out_d / "results_susceptibility_all.txt").write_text(content_all, encoding="utf-8")
            if rows_scaling:
                (out_d / "results_susceptibility_scaling.txt").write_text(content_scaling, encoding="utf-8")
            (out_d / "chi_disc.txt").write_text(content_lcp, encoding="utf-8")
            (out_d / "results_beta.txt").write_text(content_lcp, encoding="utf-8")

            # 2. 导出完整 CSV 与 Parquet 格式
            df_export_lcp = pl.DataFrame(rows_lcp).rename({
                "chi_lat": "Mean_vol_scaled",
                "err_lat": "Error_vol_scaled",
                "chi_ren": "Mean_scaled_gev2_renorm",
                "err_ren": "Error_scaled_gev2_renorm",
            })
            df_export_lcp.write_csv(out_d / "results_susceptibility.csv")
            df_export_lcp.write_parquet(out_d / "results_susceptibility.parquet")

            df_export_all = pl.DataFrame(rows_all).rename({
                "chi_lat": "Mean_vol_scaled",
                "err_lat": "Error_vol_scaled",
                "chi_ren": "Mean_scaled_gev2_renorm",
                "err_ren": "Error_scaled_gev2_renorm",
            })
            df_export_all.write_csv(out_d / "all_ensembles_susceptibility.csv")

        print(f"[OK] 成功导出 L48T16 各个 beta 重整化手征磁化率数据产物至: {[str(d) for d in target_dirs]}")

    def run_cluster_scan(self, readin_dir: Path) -> pl.DataFrame:
        """
        在拥有全量构型测量数据的计算机/集群上执行全温度与全系综扫描计算。
        自动遍历 readin_dir 下的所有格点目录，智能过滤介子关联函数目录，
        计算手征磁化率，导出至 output/condensate/results_susceptibility.csv 及 outputs/LCP/。
        """
        all_dirs = sorted([d for d in readin_dir.iterdir() if d.is_dir()])
        results = []

        print(f"[SusceptibilityPipeline] 正在从数据目录 {readin_dir} 执行集群全量扫描，共发现 {len(all_dirs)} 个子目录...")

        for d in all_dirs:
            if not is_valid_condensate_dir(d):
                if (d / "Output").exists():
                    print(f"  [跳过介子目录] {d.name} (包含 Output/ 强子关联函数，无 PsibarPsi 随机源)")
                continue

            meta = parse_ensemble_meta_from_dir(d.name)
            try:
                from docs.physics_setup import parse_ensemble_dirname as p_ed, calculate_temperature as c_temp
                rich_meta = p_ed(d.name)
                if rich_meta:
                    meta.update(rich_meta)
            except Exception:
                pass

            eff_ns = int(meta.get("ns", 48))
            eff_nt = int(meta.get("nt", 16))
            eff_beta_str = str(meta.get("beta", "4.17"))
            try:
                eff_beta_val = float(eff_beta_str)
            except ValueError:
                eff_beta_val = 4.17

            # 温度计算：优先使用标准几何比例 T = T(Nt=16) * 16 / Nt
            try:
                from docs.physics_setup import calculate_temperature
                eff_temp = float(calculate_temperature(eff_beta_val, eff_nt))
            except Exception:
                base_temp = SUSCEPTIBILITY_TEMP_MAP.get(eff_beta_str, 153.31)
                eff_temp = float(base_temp * (16.0 / eff_nt))

            eff_zm = get_zm_factor(eff_beta_str)
            inv_zm = (1.0 / eff_zm ** 2) if abs(eff_zm) > 1e-15 else 1.0

            print(f"  -> 正在处理 [{d.name}] (Ns={eff_ns}, Nt={eff_nt}, beta={eff_beta_str}, T={eff_temp:.2f} MeV, Zm={eff_zm:.6f})...")
            try:
                res = self.extract_from_meas_directory(d, ns=eff_ns, nt=eff_nt, temp_mev=eff_temp)
                m_gev2 = res["mean_scaled"] / 1e6
                e_gev2 = res["error_scaled"] / 1e6
                results.append({
                    "Ensemble": d.name,
                    "Beta": eff_beta_val,
                    "Temp": eff_temp,
                    "Ns": eff_ns,
                    "Nt": eff_nt,
                    "Num_cfgs": res.get("num_configs", 0),
                    "Zm": eff_zm,
                    "Mean_unscaled": res["mean_unscaled"],
                    "Error_unscaled": res["error_unscaled"],
                    "Mean_vol_scaled": res["mean_vol_scaled"],
                    "Error_vol_scaled": res["error_vol_scaled"],
                    "Mean_scaled": res["mean_scaled"],
                    "Error_scaled": res["error_scaled"],
                    "Mean_scaled_gev2": m_gev2,
                    "Error_scaled_gev2": e_gev2,
                    "Mean_scaled_gev2_renorm": m_gev2 * inv_zm,
                    "Error_scaled_gev2_renorm": e_gev2 * inv_zm,
                })
            except Exception as e:
                print(f"  [WARN] 处理 {d.name} 失败: {e}")

        if not results:
            print("[WARN] 未能在给定数据目录中提取到任何有效数据，载入已有全量基准数据...")
            return self.get_full_scan_results()

        df_res = pl.DataFrame(results).sort(["Beta", "Temp"])

        out_csv = self.output_dir / "results_susceptibility.csv"
        out_parquet = self.output_dir / "results_susceptibility.parquet"
        out_all_csv = self.output_dir / "all_ensembles_susceptibility.csv"

        df_res.write_csv(out_csv)
        df_res.write_parquet(out_parquet)
        df_res.write_csv(out_all_csv)

        # 导出至 outputs/LCP/ 与 output/LCP/
        self.export_lcp_outputs(df_res)

        print(f"[OK] 集群计算完成，结果已保存至 {out_csv} 及 {out_all_csv}")
        return df_res


    def get_full_scan_results(
        self,
        force_recompute: bool = False,
    ) -> pl.DataFrame:
        """
        获取手征磁化率完整数据集 (含 Zm 重整化与物理标度)。
        优先载入包含全量真实测量的 all_ensembles_susceptibility.csv 并注入 Zm 重整化；
        若不存在则自 ana 基准数据对齐同步，并自动导出至 outputs/LCP/。
        """
        out_csv = self.output_dir / "results_susceptibility.csv"
        out_parquet = self.output_dir / "results_susceptibility.parquet"
        out_all_csv = self.output_dir / "all_ensembles_susceptibility.csv"
        ana_csv = self.ana_root / "ttest" / "results_susceptibility.csv"

        # 1. 优先从已有全量系综汇总表载入并计算重整化列
        if out_all_csv.exists():
            df_exist = pl.read_csv(out_all_csv)
            if "Ensemble" in df_exist.columns and df_exist.height > 0:
                rows = []
                for r in df_exist.iter_rows(named=True):
                    b = float(r["Beta"])
                    zm = float(r.get("Zm", r.get("zm", get_zm_factor(b))))
                    inv_zm = (1.0 / zm **2 ) if abs(zm) > 1e-15 else 1.0
                    m_scaled = float(r["Mean_scaled"])
                    e_scaled = float(r["Error_scaled"])
                    m_gev2 = m_scaled / 1e6
                    e_gev2 = e_scaled / 1e6
                    row_dict = dict(r)
                    row_dict["Beta"] = b
                    row_dict["Zm"] = zm
                    row_dict["Mean_scaled_gev2"] = m_gev2
                    row_dict["Error_scaled_gev2"] = e_gev2
                    row_dict["Mean_scaled_gev2_renorm"] = m_gev2 * inv_zm
                    row_dict["Error_scaled_gev2_renorm"] = e_gev2 * inv_zm
                    rows.append(row_dict)
                df_enriched = pl.DataFrame(rows).sort(["Beta", "Temp"])
                df_enriched.write_csv(out_csv)
                df_enriched.write_parquet(out_parquet)
                df_enriched.write_csv(out_all_csv)
                self.export_lcp_outputs(df_enriched)
                print(f"[SusceptibilityPipeline] 成功同步并重整化全量 {df_enriched.height} 组格点系综磁化率数据")
                return df_enriched

        if not force_recompute and out_csv.exists():
            return pl.read_csv(out_csv)

        # 2. 兜底读取 ana/ttest 基准数据
        if ana_csv.exists():
            df_norm = pl.read_csv(ana_csv)
        else:
            raise FileNotFoundError(f"未找到基准数据: {ana_csv}")

        # 标准化添加 Ns=48, Nt=16, 体积因子, 标度因子与 Zm 重整化
        rows = []
        for r in df_norm.iter_rows(named=True):
            b = r["Beta"]
            t = float(r["Temp"])
            mu = float(r["Mean_unscaled"])
            eu = float(r["Error_unscaled"])
            f_vol, f_scaled = compute_scaling_factors(48, 16, t)
            zm = get_zm_factor(b)
            inv_zm = (1.0 / zm**2) if abs(zm) > 1e-15 else 1.0
            m_scaled = mu * f_scaled
            e_scaled = eu * f_scaled
            m_gev2 = m_scaled / 1e6
            e_gev2 = e_scaled / 1e6

            b_str = f"{float(b):.2f}" if abs(float(b) - 4.405) > 1e-4 else "4.405"
            rows.append({
                "Ensemble": f"L48T16beta{b_str}",
                "Beta": float(b),
                "Temp": t,
                "Ns": 48,
                "Nt": 16,
                "Num_cfgs": int(r.get("Num_cfgs", 2000)),
                "Zm": zm,
                "Mean_unscaled": mu,
                "Error_unscaled": eu,
                "Mean_vol_scaled": mu * f_vol,
                "Error_vol_scaled": eu * f_vol,
                "Mean_scaled": m_scaled,
                "Error_scaled": e_scaled,
                "Mean_scaled_gev2": m_gev2,
                "Error_scaled_gev2": e_gev2,
                "Mean_scaled_gev2_renorm": m_gev2 * inv_zm,
                "Error_scaled_gev2_renorm": e_gev2 * inv_zm,
            })

        df_enriched = pl.DataFrame(rows).sort("Temp")
        df_enriched.write_csv(out_csv)
        df_enriched.write_parquet(out_parquet)
        df_enriched.write_csv(out_all_csv)
        self.export_lcp_outputs(df_enriched)

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
