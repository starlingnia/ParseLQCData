#!/usr/bin/env python3
"""
scripts/reproduce_susceptibility.py
--------------------------------------------------------------------------------
手征磁化率 (Chiral Susceptibility) 算法复现与数据抽取执行脚本 (纯轻夸克无偏估计，不执行 RM 减除)
1. 具备完整的微观随机源抽取逻辑、单构型无偏二次交叉乘积估计器与 Jackknife 误差分析
2. 精确计算格点几何与温度因子:
   - 4D 时空体积因子: F_vol = Ns^3 * Nt
   - 标度因子: F_scaled = Ns^3 * Nt^3 * T^2 = F_vol * (Nt * T)^2
3. 支持两种工作模式:
   - 跨机器/集群模式: 指定 --readin-dir 遍历拥有真实构型向量的集群目录，全自动提取全量结果
   - 本地验证模式: 对本地现有小规模数据集进行端到端全真 XML 抽取，并与基准数据对齐
4. 导出规范化 CSV/Parquet 产物至 output/condensate/
--------------------------------------------------------------------------------
"""

import argparse
import os
from pathlib import Path
import sys

import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.physics_setup import OUTPUT_CONDENSATE_DIR
from src.parselqcdata.susceptibility_pipeline import SusceptibilityPipeline


def main():
    parser = argparse.ArgumentParser(description="ParseLQCData 手征磁化率抽取与复现工具")
    parser.add_argument(
        "--readin-dir",
        type=str,
        default=None,
        help="包含各温度或格点构型 meas.* 的原始数据根目录（用于有全量数据的集群机器）",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(OUTPUT_CONDENSATE_DIR),
        help="结果 CSV 与 Parquet 文件的保存目录",
    )
    parser.add_argument(
        "--lcp-dir",
        type=str,
        default=str(PROJECT_ROOT / "output" / "LCP"),
        help="重整化 LCP 规范结果输出目录 (默认 output/LCP/)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="计算完成后自动调用 scripts/plot_susceptibility.py 生成高清图表",
    )
    args = parser.parse_args()

    print("=" * 85)
    print("      ParseLQCData 手征磁化率 (Chiral Susceptibility) 计算与物理标度工具      ")
    print("=" * 85)

    pipeline = SusceptibilityPipeline(output_dir=Path(args.output_dir))

    # 1. 如果用户显式传入了集群数据根目录，则在全量真实数据上跑遍全谱
    if args.readin_dir and Path(args.readin_dir).exists():
        readin_path = Path(args.readin_dir)
        print(f"\n[运行模式] 检测到集群原始数据目录: {readin_path}")
        print("正在启动全温度多构型端到端流式解析与 Jackknife 重采样...")
        df_res = pipeline.run_cluster_scan(readin_path)
    else:
        # 本地模式：对本地已有 XML 进行全真逻辑提取，并同步全量数据
        from docs.physics_setup import resolve_dataset_dir, DEFAULT_READIN_DIR
        local_test_dir = resolve_dataset_dir(DEFAULT_READIN_DIR, "L32T12beta4.17")
        if local_test_dir.exists():
            print(f"\n[本地验证] 正在从本地真实 XML 向量中测试抽取 L32T12beta4.17 数据...")
            r = pipeline.extract_from_meas_directory(
                ensemble_dir=local_test_dir,
                ns=32,
                nt=12,
                temp_mev=204.41,
            )
            print(f"  -> 构型数: {r['num_configs']}, Ns={r['ns']}, Nt={r['nt']}, T={r['temp']:.1f} MeV, Zm={r.get('zm', 0.966247):.6f}")
            print(f"     原始晶格量:  χ_unscaled = {r['mean_unscaled']:.8e} ± {r['error_unscaled']:.8e}")
            print(f"     4D体积标度:  χ_vol      = {r['mean_vol_scaled']:.6f} ± {r['error_vol_scaled']:.6f} (F_vol = {r['factor_vol']:.0f})")
            print(f"     标度量:      χ_scaled   = {r['mean_scaled']:.2e} ± {r['error_scaled']:.2e} MeV²")
            print(f"     重整化标度:  χ_renorm   = {r.get('mean_scaled_gev2_renorm', 0.0):.6f} ± {r.get('error_scaled_gev2_renorm', 0.0):.6f} GeV²")

        print(f"\n[数据同步] 正在载入并计算全量格点系综标准标度因子与质量重整化因子 Zm...")
        df_res = pipeline.get_full_scan_results(force_recompute=True)

    target_lcp_dirs = [Path(args.lcp_dir), PROJECT_ROOT / "outputs" / "LCP", PROJECT_ROOT / "output" / "LCP"]
    target_lcp_dirs = list(dict.fromkeys(target_lcp_dirs))
    pipeline.export_lcp_outputs(df_res, target_dirs=target_lcp_dirs)

    # 打印格式化物理总表 (聚焦 L48T16 各个 Beta 的 LCP 序列)
    df_l48 = df_res.filter((pl.col("Ns") == 48) & (pl.col("Nt") == 16)) if "Ns" in df_res.columns else df_res
    print("\n" + "=" * 115)
    print("                 L48T16 手征磁化率各个 Beta 规范计算结果 (含 Zm 重整化)")
    print("=" * 115)
    print(f"{'Ensemble':<35} {'Beta':<6} {'T[MeV]':<8} {'Size':<8} {'Zm':<8} {'χ_vol (lat unit)':<19} {'χ_renorm [GeV²]':<22}")
    print("-" * 115)
    for r in df_l48.iter_rows(named=True):
        ens_name = r.get("Ensemble", f"beta{r['Beta']}")
        size_str = f"{r['Ns']}x{r['Nt']}"
        zm_val = float(r.get("Zm", 1.0))
        chi_vol = float(r["Mean_vol_scaled"])
        err_vol = float(r["Error_vol_scaled"])
        chi_ren = float(r.get("Mean_scaled_gev2_renorm", (float(r["Mean_scaled"]) / 1e6) / zm_val))
        err_ren = float(r.get("Error_scaled_gev2_renorm", (float(r["Error_scaled"]) / 1e6) / zm_val))
        print(
            f"{ens_name:<35} {r['Beta']:<6} {r['Temp']:<8.1f} {size_str:<8} {zm_val:<8.4f} "
            f"{chi_vol:<9.4e}±{err_vol:<7.1e} "
            f"{chi_ren:<10.6e}±{err_ren:<9.2e}"
        )
    print("=" * 115)

    print(f"\n[输出路径] 数据产物已保存至:")
    print(f"  - output/LCP/ 规范结果:")
    print(f"    * 规范文本格式:    {PROJECT_ROOT / 'output' / 'LCP' / 'results_susceptibility.txt'}")
    print(f"    * 标头约定:        # beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error")
    print(f"    * 完整 CSV:        {PROJECT_ROOT / 'output' / 'LCP' / 'results_susceptibility.csv'}")
    print(f"    * 镜像目录:        {PROJECT_ROOT / 'outputs' / 'LCP' / 'results_susceptibility.txt'}")

    if args.plot:
        from scripts.plot_susceptibility import main as plot_main
        plot_main()


if __name__ == "__main__":
    main()
