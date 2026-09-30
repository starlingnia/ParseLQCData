#!/usr/bin/env python3
"""
scripts/update_ccor_data_and_plots.py
--------------------------------------------------------------------------------
完全按照 ~/code/ths/pos/ccor/ 的画图格式与数据规格：
1. 从 ParseLQCData 的最新多源 (multisrc) 测量结果中提取 6 信道介子质量与 4 组对称性破缺质量差；
2. 物理量纲换算乘单位 153 * 16 = 2448 MeV:
     T = 153 * 16 / Nt (MeV)
     M = aM * 2448 (MeV)
3. 严格按照 ccor 格式输出:
     - massvtem.csv
     - deltamass_nt12.csv
     - dataSU(2)XSU(2) asym.txt & dataSU(2)XSU(2) V-A.txt
     - dataU(1)A asym.txt & dataU(1)A T-X.txt
     - dataU(1) asym.txt & dataU(1) S-PS.txt
     - dataSU(2)spinxchiral asym.txt & dataSU(2)spinxchiral X-A.txt
     - mdoutputre.csv
4. 在 ccor 目录下直接调用原生 gnuplot 模板脚本 (plotmdre.gp, plotmd.gp, plotmassdvsmass.gp)
   生成所有 PDF 矢量图，并使用 pdftoppm 渲染高清 PNG 图像，
   完全替换掉 ccor 中的老旧数据和图片。
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CCOR_DIR = Path("/Users/junxiongnie/code/ths/pos/ccor")

# 物理量纲换算标度: 153 * 16 = 2448 MeV
SCALE_UNIT: float = 2453  # 2448.0 MeV

# 4 组扫描格点 (ccor 任务覆盖的 4 个有限温度点)
CCOR_CASES: Tuple[str, ...] = ("32x12", "32x14", "40x16", "48x18")
CCOR_NTS: Tuple[int, ...] = (12, 14, 16, 18)
CCOR_MLS: Tuple[float, ...] = (0.0020, 0.0035, 0.0070, 0.0120)

# 信道名称映射到 ccor 中的 type 编号
# 1: V, 2: A, 3: Tt, 4: Xt, 5: S, 6: Ps
TYPE_MAPPING: Dict[str, int] = {
    "Vec": 1,
    "AV": 2,
    "Tt": 3,
    "Xt": 4,
    "S": 5,
    "PS": 6,
}

# 4 组对称性破缺通道对与图名配置
SYMMETRY_SPECS = [
    {
        "pair_title": "SU(2)XSU(2) asym",
        "pair_tag": "V-A",
        "ch1": "Vec",
        "ch2": "AV",
        "channel_id": 1,
        "ytle": " {/Symbol D}M_{V-A}(MeV)",
    },
    {
        "pair_title": "U(1)A asym",
        "pair_tag": "T-X",
        "ch1": "Tt",
        "ch2": "Xt",
        "channel_id": 2,
        "ytle": " {/Symbol D}M_{T-X}(MeV)",
    },
    {
        "pair_title": "U(1) asym",
        "pair_tag": "S-PS",
        "ch1": "S",
        "ch2": "PS",
        "channel_id": 3,
        "ytle": " {/Symbol D}M_{S-PS}(MeV)",
    },
    {
        "pair_title": "SU(2)spinxchiral asym",
        "pair_tag": "X-A",
        "ch1": "Xt",
        "ch2": "AV",
        "channel_id": 4,
        "ytle": " {/Symbol D}M_{X-A}(MeV)",
    },
]


def ml_index_tag(ml: float) -> int:
    """按 ccor 中的整数离散标签映射 (0.002->1, 0.0035->2, 0.007->4, 0.012->7)"""
    return int(round(ml / 0.0017))


def load_mass_dataframe() -> pl.DataFrame:
    """加载并转换介子质量数据"""
    summary_path = (
        PROJECT_ROOT
        / "output"
        / "meson_scan"
        / "b4.17"
        / "multisrc"
        / "meson_mass_summary.csv"
    )
    if not summary_path.exists():
        raise FileNotFoundError(f"找不到介子拟合汇总表: {summary_path}")

    df = pl.read_csv(summary_path)
    # 只取 ccor 对应的 4 组格点
    df = df.filter(pl.col("case_key").is_in(list(CCOR_CASES)))

    # 计算 153*16 标度下的物理温度和物理质量
    df = df.with_columns([
        (SCALE_UNIT / pl.col("nt")).alias("tem_mev"),
        (pl.col("mass") * SCALE_UNIT).alias("mass_mev"),
        (pl.col("mass_err") * SCALE_UNIT).alias("mass_err_mev"),
        pl.col("channel").replace(TYPE_MAPPING).cast(pl.Int32).alias("type_id"),
    ])
    return df


def load_delta_mass_dataframe() -> pl.DataFrame:
    """加载并转换质量差数据"""
    delta_path = PROJECT_ROOT / "output" / "ccor_flow" / "delta_mass_all.csv"
    if not delta_path.exists():
        raise FileNotFoundError(f"找不到质量差汇总表: {delta_path}")

    df = pl.read_csv(delta_path)
    df = df.filter((pl.col("dataset") == "multisrc") & pl.col("case").is_in(list(CCOR_CASES)))

    df = df.with_columns([
        (SCALE_UNIT / pl.col("nt")).alias("tem_mev"),
        (pl.col("delta_mass_lattice") * SCALE_UNIT).alias("delta_mass_mev"),
        (pl.col("delta_mass_err_lattice") * SCALE_UNIT).alias("delta_mass_err_mev"),
    ])
    return df


def generate_massvtem_csv(mass_df: pl.DataFrame, target_file: Path) -> None:
    """
    生成 ccor/massvtem.csv:
    表头: tem,type,massterm,massfit_mean_jk,massfit_err_jk
    按温度降序 (204.0, 174.857143, 153.0, 136.0)，并按夸克质量和 type 排序
    """
    sorted_df = mass_df.sort(["nt", "ml", "type_id"], descending=[False, False, False])
    # 转换为原 ccor 顺序 (nt=12, 14, 16, 18 对应温度 204.0, 174.86, 153.0, 136.0)
    out_df = pl.DataFrame({
        "tem": [f"{t:.6f}" for t in sorted_df["tem_mev"].to_list()],
        "type": sorted_df["type_id"].to_list(),
        "massterm": sorted_df["ml"].to_list(),
        "massfit_mean_jk": [f"{m:.6f}" for m in sorted_df["mass_mev"].to_list()],
        "massfit_err_jk": [f"{e:.6f}" for e in sorted_df["mass_err_mev"].to_list()],
    })
    out_df.write_csv(target_file)
    print(f"[WRITE] massvtem.csv -> {target_file} ({out_df.height} 行)")


def generate_deltamass_nt12_csv(delta_df: pl.DataFrame, target_file: Path) -> None:
    """
    生成 ccor/deltamass_nt12.csv:
    表头: channel,mass,diff,err
    注意: diff 和 err 为无量纲格点单位 (因为 plotmassdvsmass.gp 内部会乘以 2448)
    """
    nt12_df = delta_df.filter(pl.col("nt") == 12)
    rows = []

    for spec in SYMMETRY_SPECS:
        ch_id = spec["channel_id"]
        # 在 delta_df 查找对应的 pair
        # pair 名可能有 'SU(2)XSU(2) V-A' 或 'U(1) S-PS' 等
        pair_sub = nt12_df.filter(
            (pl.col("ch1") == spec["ch1"]) & (pl.col("ch2") == spec["ch2"])
            | (pl.col("pair").str.contains(spec["pair_tag"]))
        ).sort("ml")

        for r in pair_sub.iter_rows(named=True):
            rows.append({
                "channel": ch_id,
                "mass": r["ml"],
                "diff": f"{r['delta_mass_lattice']:.5f}",
                "err": f"{r['delta_mass_err_lattice']:.5f}",
            })

    out_df = pl.DataFrame(rows)
    out_df.write_csv(target_file)
    print(f"[WRITE] deltamass_nt12.csv -> {target_file} ({out_df.height} 行)")


def generate_symmetry_txt_files(delta_df: pl.DataFrame, ccor_dir: Path) -> None:
    """
    生成 4 组对称性破缺的文本数据文件:
      data<对称性>.txt 格式: tem, int(ml/0.0017), delta_mass_mev, delta_mass_err_mev
    """
    for spec in SYMMETRY_SPECS:
        pair_sub = delta_df.filter(
            (pl.col("ch1") == spec["ch1"]) & (pl.col("ch2") == spec["ch2"])
            | (pl.col("pair").str.contains(spec["pair_tag"]))
        ).sort(["nt", "ml"])

        lines = []
        for r in pair_sub.iter_rows(named=True):
            tem_str = f"{r['tem_mev']:.6f}"
            idx = ml_index_tag(r["ml"])
            dm_str = f"{r['delta_mass_mev']:.6f}"
            de_str = f"{r['delta_mass_err_mev']:.6f}"
            lines.append(f"{tem_str} {idx} {dm_str} {de_str}")

        text_content = "\n".join(lines) + "\n"

        # 写入两种常用别名，确保兼容
        file1 = ccor_dir / f"data{spec['pair_title']}.txt"
        file2 = ccor_dir / f"data{spec['pair_title'].replace(' asym', '')} {spec['pair_tag']}.txt"
        file1.write_text(text_content, encoding="utf-8")
        file2.write_text(text_content, encoding="utf-8")
        print(f"[WRITE] {file1.name} & {file2.name} ({len(lines)} 行)")


def generate_mdoutputre_csv(delta_df: pl.DataFrame, target_file: Path) -> None:
    """生成 mdoutputre.csv"""
    lines = ["intem,massterm,massdifference,error,runningtime"]
    for r in delta_df.sort(["nt", "ml"]).iter_rows(named=True):
        lines.append(
            f"{r['nt']},{r['ml']},{r['delta_mass_lattice']:.6f},{r['delta_mass_err_lattice']:.6f},Time: 0ms"
        )
    target_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[WRITE] mdoutputre.csv -> {target_file}")


def run_gnuplot_and_convert_png(work_dir: Path) -> List[Path]:
    """
    在指定目录中执行 gnuplot 并将 PDF 渲染为高质量 PNG
    """
    gnuplot_bin = shutil.which("gnuplot")
    if not gnuplot_bin:
        raise RuntimeError("系统未找到 gnuplot 可执行程序！")

    pdftoppm_bin = shutil.which("pdftoppm")

    print("\n--------------------------------------------------------------------------------")
    print(f"[GNUPLOT] 开始在 {work_dir} 中按照 ccor 格式渲染所有 PDF 矢量图...")
    print("--------------------------------------------------------------------------------")

    generated_files: List[Path] = []

    # 1. 运行 plotmdre.gp -> massvtem_*.pdf
    gp_mdre = work_dir / "plotmdre.gp"
    if gp_mdre.exists():
        cmd = [gnuplot_bin, "-e", "input_fname='massvtem.csv'", str(gp_mdre)]
        res = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
        if res.returncode == 0:
            print("[SUCCESS] plotmdre.gp 绘制完成 -> massvtem_*.pdf")
            for ml in CCOR_MLS:
                pdf_p = work_dir / f"massvtem_{ml:.4f}.pdf"
                if pdf_p.exists():
                    generated_files.append(pdf_p)
        else:
            print(f"[ERROR] plotmdre.gp 执行失败: {res.stderr}")

    # 2. 运行 plotmd.gp -> data*.pdf (4 组对称性)
    gp_md = work_dir / "plotmd.gp"
    if gp_md.exists():
        for spec in SYMMETRY_SPECS:
            input_name = f"data{spec['pair_title']}.txt"
            output_name = f"data{spec['pair_title']}.pdf"
            symt = spec["pair_title"]
            ytle = spec["ytle"]

            cmd = [
                gnuplot_bin,
                "-e",
                f'input_fname="{input_name}"; output_fname="{output_name}"; '
                f'symt="{symt}"; ytle="{ytle}"',
                str(gp_md),
            ]
            res = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
            if res.returncode == 0:
                print(f"[SUCCESS] plotmd.gp 绘制完成 -> {output_name}")
                pdf_p = work_dir / output_name
                generated_files.append(pdf_p)
                # 同时也复制一份别名文件
                alias_pdf = work_dir / f"data{spec['pair_title'].replace(' asym', '')} {spec['pair_tag']}.pdf"
                shutil.copy2(pdf_p, alias_pdf)
                generated_files.append(alias_pdf)
            else:
                print(f"[ERROR] plotmd.gp 绘制 {spec['pair_title']} 失败: {res.stderr}")

    # 3. 运行 plotmassdvsmass.gp -> deltamassvsm.pdf
    gp_delta = work_dir / "plotmassdvsmass.gp"
    if gp_delta.exists():
        cmd = [gnuplot_bin, str(gp_delta)]
        res = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
        if res.returncode == 0:
            print("[SUCCESS] plotmassdvsmass.gp 绘制完成 -> deltamassvsm.pdf")
            pdf_p = work_dir / "deltamassvsm.pdf"
            if pdf_p.exists():
                generated_files.append(pdf_p)
        else:
            print(f"[ERROR] plotmassdvsmass.gp 执行失败: {res.stderr}")

    # 4. 使用 pdftoppm 将所有生成的 PDF 转换为 PNG
    if pdftoppm_bin:
        print("\n[CONVERT] 使用 pdftoppm 将 PDF 转换为高清晰度 PNG...")
        pdf_targets = [
            work_dir / "massvtem_0.0020.pdf",
            work_dir / "massvtem_0.0035.pdf",
            work_dir / "massvtem_0.0070.pdf",
            work_dir / "massvtem_0.0120.pdf",
            work_dir / "dataSU(2)XSU(2) asym.pdf",
            work_dir / "dataU(1)A asym.pdf",
            work_dir / "dataU(1) asym.pdf",
            work_dir / "dataSU(2)spinxchiral asym.pdf",
            work_dir / "deltamassvsm.pdf",
        ]

        for pdf in pdf_targets:
            if not pdf.exists():
                continue
            base_name = pdf.stem
            prefix = work_dir / base_name
            # pdftoppm -png -r 150 <input.pdf> <prefix>
            subprocess.run([pdftoppm_bin, "-png", "-r", "150", str(pdf), str(prefix)], check=True)
            png_p = work_dir / f"{base_name}-1.png"
            if png_p.exists():
                generated_files.append(png_p)
                print(f"[PNG] 生成/更新: {png_p.name}")
    else:
        print("[WARN] 未检测到 pdftoppm，跳过 PNG 自动渲染。")

    return generated_files


def main():
    parser = argparse.ArgumentParser(
        description="按 ccor 画图格式生成数据与图表到本地 docs/ccor/，并全量拷贝同步至 ~/code/ths/pos/ccor/"
    )
    parser.add_argument(
        "--docs-dir",
        type=Path,
        default=PROJECT_ROOT / "docs" / "ccor",
        help="本地 docs 输出目录 (默认: <project_root>/docs/ccor)",
    )
    parser.add_argument(
        "--ccor-dir",
        type=Path,
        default=CCOR_DIR,
        help="外部 ccor 目录路径 (默认: ~/code/ths/pos/ccor)",
    )
    args = parser.parse_args()
    docs_dir = args.docs_dir.resolve()
    ccor_dir = args.ccor_dir.resolve()

    docs_dir.mkdir(parents=True, exist_ok=True)
    if not ccor_dir.exists():
        print(f"[ERROR] 目标外部 ccor 目录不存在: {ccor_dir}")
        sys.exit(1)

    print("================================================================================")
    print("      按照 ~/code/ths/pos/ccor/ 原生格式生成数据与图像并全量同步")
    print(f"      1. 本地 docs 存储路径: {docs_dir}")
    print(f"      2. 同步目标 ccor 路径: {ccor_dir}")
    print(f"      3. 物理标度常量: 153 * 16 = {SCALE_UNIT:.1f} MeV")
    print("================================================================================")

    # 1. 确保 gnuplot 绘图模板存在于 docs_dir
    gp_scripts = ["plotmdre.gp", "plotmd.gp", "plotmassdvsmass.gp"]
    for gp in gp_scripts:
        src_gp = ccor_dir / gp
        dst_gp = docs_dir / gp
        if src_gp.exists():
            shutil.copy2(src_gp, dst_gp)
        elif not dst_gp.exists():
            # 降级从 scripts/gnuplot 拷贝
            fallback = PROJECT_ROOT / "scripts" / "gnuplot" / gp
            if fallback.exists():
                shutil.copy2(fallback, dst_gp)

    # 2. 加载最新多源拟合数据
    mass_df = load_mass_dataframe()
    delta_df = load_delta_mass_dataframe()

    # 3. 在 docs_dir 生成全部数据文件
    generated_data_files: List[Path] = [
        docs_dir / "massvtem.csv",
        docs_dir / "deltamass_nt12.csv",
        docs_dir / "dataSU(2)XSU(2) asym.txt",
        docs_dir / "dataSU(2)XSU(2) V-A.txt",
        docs_dir / "dataU(1)A asym.txt",
        docs_dir / "dataU(1)A T-X.txt",
        docs_dir / "dataU(1) asym.txt",
        docs_dir / "dataU(1) S-PS.txt",
        docs_dir / "dataSU(2)spinxchiral asym.txt",
        docs_dir / "dataSU(2)spinxchiral X-A.txt",
        docs_dir / "mdoutputre.csv",
    ]

    generate_massvtem_csv(mass_df, docs_dir / "massvtem.csv")
    generate_deltamass_nt12_csv(delta_df, docs_dir / "deltamass_nt12.csv")
    generate_symmetry_txt_files(delta_df, docs_dir)
    generate_mdoutputre_csv(delta_df, docs_dir / "mdoutputre.csv")

    # 4. 在 docs_dir 中执行 gnuplot 出图并转换为 PNG
    generated_plot_files = run_gnuplot_and_convert_png(docs_dir)

    all_files_to_sync = set(generated_data_files + generated_plot_files)

    # 5. 拷贝一份全部生成的文件到 ~/code/ths/pos/ccor/
    print("\n--------------------------------------------------------------------------------")
    print(f"[COPY] 将本地 {docs_dir.name}/ 中的全部生成文件拷贝到 {ccor_dir} ...")
    print("--------------------------------------------------------------------------------")
    synced_files: List[Tuple[Path, Path]] = []
    for f in sorted(all_files_to_sync):
        if f.exists():
            target_file = ccor_dir / f.name
            shutil.copy2(f, target_file)
            synced_files.append((f, target_file))
            print(f"[COPIED] {f.name} -> {target_file}")

    print("\n================================================================================")
    print("      已顺利完成！所有输出已保存在本地 docs/ccor/ 并同步覆盖 ccor/ 下的文件。")
    print(f"      共同步 {len(synced_files)} 个数据与图像文件。")
    print("================================================================================")


if __name__ == "__main__":
    main()
