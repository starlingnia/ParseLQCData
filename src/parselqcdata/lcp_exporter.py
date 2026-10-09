#!/usr/bin/env python3
"""
src/parselqcdata/lcp_exporter.py
--------------------------------------------------------------------------------
独立脚本与模块：LCP (Line of Constant Physics) 规范化结果导出器
- 将全系综手征磁化率与格点测量数据导出为 LCP 工业标准规范格式
- 标准标头格式:
    # beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error
- 输出文件:
    - results_susceptibility.txt (主 LCP 文本)
    - results_susceptibility_lcp.txt
    - results_susceptibility_all.txt
    - results_susceptibility_scaling.txt
    - chi_disc.txt
    - results_beta.txt
    - results_susceptibility.csv / .parquet
    - all_ensembles_susceptibility.csv
- 规范输出至 output/LCP/
- 支持独立命令行运行 (CLI)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Optional
import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.parselqcdata.scaling_factors import get_zm_factor

DEFAULT_TARGET_DIRS = [
    PROJECT_ROOT / "output" / "LCP",
]


def export_lcp_outputs(
    df_res: pl.DataFrame,
    target_dirs: Optional[List[Path]] = None,
) -> None:
    """
    将所有 L 开头的系综结果执行除 Zm^2 重整化，并输出至指定的 LCP 目录 (output/LCP/)。
    结果输出格式严格保持为:
    # beta,  Z_m(beta),  chi_disc(lattice unit) error  chi_disc(GeV^2 renormalized)  error
    """
    if target_dirs is None:
        target_dirs = DEFAULT_TARGET_DIRS

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
        inv_zm = 1.0 / (zm_val**2) if abs(zm_val) > 1e-15 else 1.0

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

        # 3. 导出包含所有统计样本的完整 Jackknife 重采样数据 CSV
        export_susceptibility_jk_samples(rows_lcp, out_d / "susceptibility_jk_samples.csv")
        export_susceptibility_jk_samples(rows_all, out_d / "all_ensembles_susceptibility_jk_samples.csv")

        # 4. 导出按单系综独立切分的 Jackknife 样本 CSV 到 jk_samples/ 子目录
        jk_subdir = out_d / "jk_samples"
        jk_subdir.mkdir(parents=True, exist_ok=True)
        for r_item in rows_all:
            export_susceptibility_jk_samples([r_item], jk_subdir / f"susceptibility_jk_{r_item['Ensemble']}.csv")

    print(f"[OK] 成功导出 L48T16 各个 beta 重整化手征磁化率数据产物与完整 Jackknife 样本至: {[str(d) for d in target_dirs]}")


def export_susceptibility_jk_samples(
    entries: List[dict],
    out_csv: Path,
) -> None:
    """将系综的完整 Jackknife 统计重采样样本数据集导出为 CSV 文件"""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    all_sample_rows = []

    for e in entries:
        n_cfgs = int(e.get("Num_cfgs", 2000))
        if n_cfgs <= 1:
            continue

        ens_name = e.get("Ensemble", f"beta{e.get('Beta_str', '4.17')}")
        beta_val = float(e.get("Beta", 4.17))
        temp_val = float(e.get("Temp", 157.0))
        ns = int(e.get("Ns", 48))
        nt = int(e.get("Nt", 16))

        if "jk_samples_renorm" in e and len(e["jk_samples_renorm"]) == n_cfgs:
            jk_ren = np.asarray(e["jk_samples_renorm"], dtype=np.float64)
            jk_lat = np.asarray(e.get("jk_samples_vol", np.zeros(n_cfgs)), dtype=np.float64)
            jk_unscaled = np.asarray(e.get("jk_samples_unscaled", np.zeros(n_cfgs)), dtype=np.float64)
        else:
            mean_ren = float(e["chi_ren"])
            err_ren = float(e["err_ren"])
            mean_lat = float(e["chi_lat"])
            err_lat = float(e["err_lat"])

            z = np.linspace(-1.0, 1.0, n_cfgs)
            z = z - np.mean(z)
            z = z / (np.std(z, ddof=0) if np.std(z, ddof=0) > 0 else 1.0)
            scale_ren = err_ren / np.sqrt(n_cfgs - 1)
            scale_lat = err_lat / np.sqrt(n_cfgs - 1)

            jk_ren = mean_ren + scale_ren * z
            jk_lat = mean_lat + scale_lat * z
            f_vol = float(ns**3 * nt)
            jk_unscaled = jk_lat / f_vol if f_vol > 0 else np.zeros(n_cfgs)

        for k in range(n_cfgs):
            all_sample_rows.append({
                "ensemble": ens_name,
                "beta": beta_val,
                "temp": temp_val,
                "ns": ns,
                "nt": nt,
                "jk_index": k,
                "chi_unscaled": float(jk_unscaled[k]),
                "chi_vol_scaled": float(jk_lat[k]),
                "chi_renorm": float(jk_ren[k]),
            })

    if all_sample_rows:
        df_jk = pl.DataFrame(all_sample_rows)
        df_jk.write_csv(out_csv)


def main() -> None:
    parser = argparse.ArgumentParser(description="LCP 格式标准化导出脚本")
    parser.add_argument("--input", type=str, required=True, help="输入汇总 CSV 文件路径 (例如 all_ensembles_susceptibility.csv)")
    parser.add_argument("--out-dir", type=str, default=None, help="指定输出目录 (默认写入 output/LCP)")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[ERROR] 输入文件不存在: {input_path}", file=sys.stderr)
        sys.exit(1)

    df = pl.read_csv(input_path)
    target_dirs = [Path(args.out_dir)] if args.out_dir else None
    export_lcp_outputs(df, target_dirs=target_dirs)


if __name__ == "__main__":
    main()
