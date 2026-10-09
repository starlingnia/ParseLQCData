#!/usr/bin/env python3
"""
scripts/update_ccor_data_and_plots.py
--------------------------------------------------------------------------------
完全按照 ~/code/ths/pos/ccor/ 的画图格式与数据规格：
1. 从 ParseLQCData 的最新多源 (multisrc) 测量结果中提取 6 信道介子质量与 4 组对称性破缺质量差；
2. 物理量纲换算统一用 a^-1 = 2453 MeV 
     T = 2453 / Nt (MeV)
     M = aM * 2453 (MeV)
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
拟合窗口 (唯一来源):
  config/fit_windows.txt —— 每个 ensemble 一行 "<ens> <x_start> <x_end>",
  介子质量与对称性 ΔM 都用它, 代码里不再自动检测窗口、也没有窗口相关的开关。
  改了窗口文件后直接重跑本脚本: 检测到文件比 delta 表新就会自动重跑 ccor_flow。
--------------------------------------------------------------------------------
图例位置 (按图定死, 见 LEGEND_POSITIONS):
  对称性图 (plotmd.gp / plotmassdvsmass.gp) -> 右上角 (set key right top)
  介子质量图 (plotmdre.gp)                  -> 右下角 (set key right bottom)
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:      # 让 docs.* / src.* 可导入
    sys.path.insert(0, str(PROJECT_ROOT))
CCOR_DIR = Path("/Users/junxiongnie/code/ths/pos/ccor")

# 物理量纲换算标度: a^-1 = 2453 MeV (旧稿的 153*16 = 2448 MeV 不对, 已废弃)
SCALE_UNIT: float = 2453.0

# 4 组扫描格点 (ccor 任务覆盖的 4 个有限温度点)
CCOR_CASES: Tuple[str, ...] = ("32x12", "32x14", "40x16", "48x18")
CCOR_NTS: Tuple[int, ...] = (12, 14, 16, 18)
CCOR_MLS: Tuple[float, ...] = (0.0020, 0.0035, 0.0070, 0.0120)

# 上游数据来源
#: 介子质量从哪里现算 (与 ΔM 同一套 jk 关联函数)
MASS_SOURCE: str = "multisrc"
MESON_CASE_ROOT: Path = PROJECT_ROOT / "output" / "meson_scan" / "b4.17" / MASS_SOURCE / "cases"
DELTA_PATH: Path = PROJECT_ROOT / "output" / "ccor_flow" / "delta_mass_all.csv"
#: 拟合窗口的唯一来源 (独立文本文件)
FIT_WINDOW_FILE: Path = PROJECT_ROOT / "config" / "fit_windows.txt"
FLOW_SCRIPT: Path = PROJECT_ROOT / "scripts" / "reproduce_ccor_flow.py"
FLOW_TMP_ROOT: Path = PROJECT_ROOT / "output" / "ccor_flow_regen"

# 4 组对称性破缺通道对与图名配置
# pair: delta_mass_all.csv 里 pair 列的完整名字 (精确匹配, 不再依赖子串)
# ch1/ch2: ccor_flow 的通道记号 (V/A/Tt/Xt/S/PS), 仅作记录
SYMMETRY_SPECS = [
    {
        "pair": "SU(2)XSU(2) V-A",
        "pair_title": "SU(2)XSU(2) asym",
        "pair_tag": "V-A",
        "ch1": "V",
        "ch2": "A",
        "channel_id": 1,
        "ytle": " {/Symbol D}M_{V-A}(MeV)",
        "key_pos": "right bottom",
        "extra_cmd": "set yrange [*:40];",
    },
    {
        "pair": "U(1)A T-X",
        "pair_title": "U(1)A asym",
        "pair_tag": "T-X",
        "ch1": "Tt",
        "ch2": "Xt",
        "channel_id": 2,
        "ytle": " {/Symbol D}M_{T-X}(MeV)",
        "key_pos": "right bottom",
        "extra_cmd": "",
    },
    {
        "pair": "U(1) S-PS",
        "pair_title": "U(1) asym",
        "pair_tag": "S-PS",
        "ch1": "S",
        "ch2": "PS",
        "channel_id": 3,
        "ytle": " {/Symbol D}M_{S-PS}(MeV)",
        "key_pos": "right top",
        "extra_cmd": "",
    },
    {
        "pair": "SU(2)spinxchiral X-A",
        "pair_title": "SU(2)spinxchiral asym",
        "pair_tag": "X-A",
        "ch1": "A",
        "ch2": "Xt",
        "channel_id": 4,
        "ytle": " {/Symbol D}M_{X-A}(MeV)",
        "key_pos": "right top",
        "extra_cmd": "",
    },
]


def ml_index_tag(ml: float) -> int:
    """按 ccor 中的整数离散标签映射 (0.002->1, 0.0035->2, 0.007->4, 0.012->7)"""
    return int(round(ml / 0.0017))


# ------------------------------------------------------------------------------
# 上游数据检查 / 补齐
# ------------------------------------------------------------------------------
def _column_values(path: Path, column: str) -> set:
    """安全地读取某个 CSV 的某一列, 返回去重集合 (文件不存在/读不动时返回空集)"""
    if not path.exists():
        return set()
    try:
        df = pl.read_csv(path, columns=[column])
    except Exception:
        return set()
    return set(df[column].to_list())


def repair_delta_mass_table(cases: Sequence[str]) -> None:
    """
    保证 output/ccor_flow/delta_mass_all.csv 覆盖 cases (对称性部分的数据源)。

    做法: 重跑 ccor_flow 里"本仓库自己的"数据源 (multisrc / singlesrc, 即除 ana 复刻用的
    reference* 之外的全部), 输出到临时目录, 然后
      * 这些数据源的行    -> 用最新结果整体替换 (拟合窗口口径变化时会被刷新)
      * reference* 的旧行 -> 原样保留 (它复刻的是 ana 的固定 N=7 口径)
    这样既补齐了新 case, 又不会破坏 reference 的既有内容。
    """
    old = pl.read_csv(DELTA_PATH) if DELTA_PATH.exists() else None
    have = set(old.filter(pl.col("dataset") == "multisrc")["case"].to_list()) if old is not None else set()
    need = sorted(have | set(cases))

    regen_datasets = ["multisrc"]
    if old is not None:
        regen_datasets += sorted(
            {str(d) for d in old["dataset"].to_list()
             if not str(d).startswith("reference") and str(d) != "multisrc"}
        )

    if FLOW_TMP_ROOT.exists():
        shutil.rmtree(FLOW_TMP_ROOT)

    cmd = [
        sys.executable,
        str(FLOW_SCRIPT),
        "--datasets",
        *regen_datasets,
        "--cases",
        *need,
        "--out-root",
        str(FLOW_TMP_ROOT),
        "--no-verify",
        "--no-figures",
        "--no-gnuplot",
    ]
    print(f"[REGEN] delta_mass_all.csv 缺 case 或需刷新窗口口径, 自动重跑 ccor_flow: {regen_datasets}")
    print("        " + " ".join(cmd))
    res = subprocess.run(cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            "ccor_flow 流程执行失败 (exit {}):\n{}\n{}".format(
                res.returncode, res.stdout[-1500:], res.stderr[-1500:]
            )
        )

    fresh = pl.read_csv(FLOW_TMP_ROOT / "delta_mass_all.csv")
    columns = list(old.columns) if old is not None else list(fresh.columns)
    parts = [fresh.select(columns)]
    if old is not None:
        keep = old.filter(~pl.col("dataset").is_in(regen_datasets))
        if keep.height:
            parts.append(keep.select(columns))
    merged = pl.concat(parts, how="vertical_relaxed")
    merged.write_csv(DELTA_PATH)
    shutil.rmtree(FLOW_TMP_ROOT, ignore_errors=True)
    print(
        "[REGEN] 完成: {} (case={}, dataset={})".format(
            DELTA_PATH,
            sorted(set(merged["case"].to_list())),
            sorted(set(merged["dataset"].to_list())),
        )
    )


def ensure_delta_cases(cases: Sequence[str]) -> None:
    """
    检查 delta 表是否需要刷新: 缺 case, 或者 config/fit_windows.txt 比它新
    (改了窗口文件就直接重跑, 不需要任何额外开关)。
    """
    have = _column_values(DELTA_PATH, "case")
    missing = [c for c in cases if c not in have]
    stale = (DELTA_PATH.exists() and FIT_WINDOW_FILE.exists()
             and FIT_WINDOW_FILE.stat().st_mtime > DELTA_PATH.stat().st_mtime)

    if missing or stale:
        why = "缺 case {}".format(missing) if missing else "拟合窗口文件已更新"
        print(f"[REGEN] {why} -> 重跑 ccor_flow 刷新 ΔM 数据")
        repair_delta_mass_table(cases)
        have = _column_values(DELTA_PATH, "case")
        missing = [c for c in cases if c not in have]

    if missing:
        raise RuntimeError(
            "对称性数据源 {} 缺少 case: {}\n        当前包含: {}\n"
            "        请先运行: .venv/bin/python scripts/reproduce_ccor_flow.py "
            "--datasets multisrc --cases {}".format(
                DELTA_PATH, missing, sorted(have), " ".join(cases))
        )


def load_mass_dataframe() -> pl.DataFrame:
    """
    介子质量表: 用 config/fit_windows.txt 里的统一窗口, 从落盘的 Jackknife 关联函数现算
    —— 与 ΔM 用同一套窗口、同一套拟合代码, 不再读 meson_scan 的自动窗口结果。
    """
    from docs.meson_scan_setup import SCAN_CASES
    from src.parselqcdata.ccor_flow import (
        ANA_MIN_MASS,
        ANA_TYPE_ID,
        channel_masses_from_jk,
        jackknife_mean_err,
        load_fit_windows,
        load_meson_scan_channel,
    )

    windows = load_fit_windows()
    rows: List[dict] = []
    for key in CCOR_CASES:
        case = next((c for c in SCAN_CASES if c.key == key), None)
        if case is None:
            raise RuntimeError(f"未知 ensemble: {key}")
        for ml in CCOR_MLS:
            case_dir = MESON_CASE_ROOT / case.output_name(float(ml))
            for ch in ANA_TYPE_ID:
                win = windows.get(key, ch)
                if win is None:
                    raise RuntimeError(f"{key}/{ch} 未在 {FIT_WINDOW_FILE} 中定义拟合窗口")
                jk = load_meson_scan_channel(case_dir, ch)
                if jk is None or jk.size == 0:
                    continue
                masses, _ = channel_masses_from_jk(
                    jk, center=case.ns // 2, window=win, min_mass=ANA_MIN_MASS
                )
                mean, err = jackknife_mean_err(masses)
                rows.append({
                    "case_key": key, "ns": case.ns, "nt": case.nt, "ml": float(ml),
                    "channel": ch, "type_id": int(ANA_TYPE_ID[ch]),
                    "tem_mev": SCALE_UNIT / case.nt,
                    "mass_mev": mean * SCALE_UNIT,
                    "mass_err_mev": err * SCALE_UNIT,
                    "x_start": int(win[0]), "x_end": int(win[1]),
                })
    df = pl.DataFrame(rows)
    print("[DATA] 介子质量: {} 行; 窗口来自 {}: {}".format(
        df.height, FIT_WINDOW_FILE.name,
        ", ".join(f"{k}=[{v[0]},{v[1]})" for k, v in windows.defaults.items()
                  if k in set(CCOR_CASES))
        + (("; 信道覆盖: " + ", ".join(f"{k[0]}/{k[1]}=[{v[0]},{v[1]})"
                                       for k, v in windows.per_channel.items()))
           if windows.per_channel else "")))
    return df


def load_delta_mass_dataframe() -> pl.DataFrame:
    """加载并转换质量差数据 (必要时自动重跑 ccor_flow)"""
    ensure_delta_cases(CCOR_CASES)

    df = pl.read_csv(DELTA_PATH)
    df = df.filter((pl.col("dataset") == "multisrc") & pl.col("case").is_in(list(CCOR_CASES)))

    df = df.with_columns([
        (SCALE_UNIT / pl.col("nt")).alias("tem_mev"),
        (pl.col("delta_mass_lattice") * SCALE_UNIT).alias("delta_mass_mev"),
        (pl.col("delta_mass_err_lattice") * SCALE_UNIT).alias("delta_mass_err_mev"),
    ])

    got = set(df["case"].to_list())
    missing = [c for c in CCOR_CASES if c not in got]
    if missing:
        raise RuntimeError(
            "质量差数据缺少 case: {} (dataset=multisrc); 现有 {}".format(missing, sorted(got))
        )
    print(
        "[DATA] 对称性数据源: {} -> 4 个温度点 {}".format(
            DELTA_PATH.name, sorted(df["tem_mev"].unique().to_list(), reverse=True)
        )
    )
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


def _symmetry_subset(delta_df: pl.DataFrame, spec: dict) -> pl.DataFrame:
    """按 pair 精确名取出某一组对称性破缺的 (temperature, ml, dM, err)"""
    return delta_df.filter(pl.col("pair") == spec["pair"]).sort(["nt", "ml"])


def generate_deltamass_nt12_csv(delta_df: pl.DataFrame, target_file: Path) -> None:
    """
    生成 ccor/deltamass_nt12.csv:
    表头: channel,mass,diff,err
    注意: diff 和 err 为无量纲格点单位 (plotmassdvsmass.gp 内部会乘以 SCALE_UNIT=2453 换成 MeV)
    """
    nt12_df = delta_df.filter(pl.col("nt") == 12)
    rows = []

    for spec in SYMMETRY_SPECS:
        pair_sub = _symmetry_subset(nt12_df, spec)
        for r in pair_sub.iter_rows(named=True):
            rows.append({
                "channel": spec["channel_id"],
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
        pair_sub = _symmetry_subset(delta_df, spec)

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


#: 论文里 4 张固定 beta 表格 (T, am_l, aΔM, ΔM[MeV]) 的数据文件: 由本脚本生成并同步到 pos/ccor
#: 对应 section/data_modi_analysis.tex 里 tab:su2_l_su2_r_asym_data / u1_a_tx / u1_s_ps / su2_cs_xa
THESIS_TABLES: Dict[str, str] = {
    "SU(2)XSU(2) V-A": "thesis_V-A.csv",
    "U(1)A T-X": "thesis_T-X.csv",
    "U(1) S-PS": "thesis_S-PS.csv",
    "SU(2)spinxchiral X-A": "thesis_X-A.csv",
}


def generate_thesis_table_csvs(delta_df: pl.DataFrame, out_dir: Path) -> None:
    r"""
    写论文 4 张表格用的 CSV: 每行 = T_ref, am_l, aΔM ± err, ΔM ± err (MeV), 行尾 hline 标记
    (标记列给每个温度组的最后一行加 \hline, 最后一组不加, 由 .tex 的 late after last line 收尾)
    """
    from docs.meson_scan_setup import SCAN_CASES

    tref = {c.key: float(c.temp_ref_mev) for c in SCAN_CASES}
    for pair, fname in THESIS_TABLES.items():
        sub = delta_df.filter(pl.col("pair") == pair).sort(["nt", "ml"])
        rows = list(sub.iter_rows(named=True))
        lines = []
        for i, r in enumerate(rows):
            last_of_group = (i + 1 == len(rows)) or (rows[i + 1]["nt"] != r["nt"])
            hline = "\\hline" if (last_of_group and i + 1 != len(rows)) else ""
            lines.append("{:.2f},{:.4f},${:.5f} \\pm {:.5f}$,${:.2f} \\pm {:.2f}$,{}".format(
                tref.get(r["case"], r["tem_mev"]), r["ml"],
                r["delta_mass_lattice"], r["delta_mass_err_lattice"],
                r["delta_mass_lattice"] * SCALE_UNIT, r["delta_mass_err_lattice"] * SCALE_UNIT,
                hline))
        (out_dir / fname).write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"[WRITE] {fname} ({len(lines)} 行)")


#: 各图的图例位置: 对称性图 (ΔM) -> 右上角; 介子质量图 -> 右下角
LEGEND_POSITIONS: Dict[str, str] = {
    "plotmd.gp": "right top",             # ΔM vs T (对称性)
    "plotmassdvsmass.gp": "right top",    # ΔM vs quark mass (对称性)
    "plotmdre.gp": "right bottom",        # 介子质量 vs T
}


def apply_legend_position(gp_path: Path, position: str = "right bottom") -> None:
    """
    把 gnuplot 模板的图例统一放到指定位置 (默认右下角, 避免挡住曲线)。

    模板里原有的 `set key ...` 行会被删掉, 然后在第一条绘图命令
    (plot / splot / do for) 之前插入 `set key <position>`,
    因此无论模板原来怎么写, 最终生效的都是这里指定的位置。
    """
    if not gp_path.exists():
        return
    text = gp_path.read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if not ln.strip().startswith("set key")]

    out: List[str] = []
    inserted = False
    for ln in lines:
        stripped = ln.strip()
        if not inserted and (
            stripped.startswith("plot")
            or stripped.startswith("splot")
            or stripped.startswith("do for")
        ):
            out.append(f"set key {position}")
            inserted = True
        out.append(ln)
    if not inserted:
        out.append(f"set key {position}")

    gp_path.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"[LEGEND] {gp_path.name}: set key {position}")


#: plotmassdvsmass.gp 中把格点单位换成 MeV 的写法, 例如 ($3*2448):($4*2448)
_UNIT_SCALE_PATTERN = re.compile(r"(\$\d+\s*\*\s*)(\d+(?:\.\d+)?)")


def apply_unit_scale(gp_path: Path, scale: float = SCALE_UNIT) -> None:
    """
    把 gnuplot 模板里硬编码的 MeV 换算因子统一成 SCALE_UNIT。

    plotmassdvsmass.gp 读入的是格点单位的 dM, 需要乘 a^-1 才变成 MeV;
    不同来源的模板里可能写死 2448 / 2640 / 2452.96 等旧标度。
    这里统一改成 SCALE_UNIT (2453), 保证与 massvtem.csv / 对称性数据同一套物理量纲。
    """
    if not gp_path.exists():
        return
    text = gp_path.read_text(encoding="utf-8")
    new_text, n = _UNIT_SCALE_PATTERN.subn(lambda m: f"{m.group(1)}{scale:g}", text)
    if n:
        gp_path.write_text(new_text, encoding="utf-8")
        print(f"[SCALE] {gp_path.name}: {n} 处 MeV 换算因子 -> {scale:g}")


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
            extra_cmd = spec.get("extra_cmd", "")
            key_pos = spec.get("key_pos", "right top")

            cmd = [
                gnuplot_bin,
                "-e",
                f'input_fname="{input_name}"; output_fname="{output_name}"; '
                f'symt="{symt}"; ytle="{ytle}"; set key {key_pos}; {extra_cmd}',
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
    print(f"      3. 物理标度常量: a^-1 = {SCALE_UNIT:.0f} MeV (T = {SCALE_UNIT:.0f}/Nt, M = aM * {SCALE_UNIT:.0f})")
    print(f"      4. 图例位置: {LEGEND_POSITIONS}")
    print(f"      5. 拟合窗口文件: {FIT_WINDOW_FILE}")
    print(f"      6. 对称性数据 case: {list(CCOR_CASES)}")
    print("================================================================================")

    # 1. 确保 gnuplot 绘图模板存在于 docs_dir, 并把图例挪到指定位置
    gp_scripts = ["plotmdre.gp", "plotmd.gp", "plotmassdvsmass.gp"]
    for gp in gp_scripts:
        src_gp = PROJECT_ROOT / "scripts" / "gnuplot" / gp
        dst_gp = docs_dir / gp
        if src_gp.exists():
            shutil.copy2(src_gp, dst_gp)
        elif (ccor_dir / gp).exists():
            shutil.copy2(ccor_dir / gp, dst_gp)
        if gp != "plotmd.gp":
            apply_legend_position(dst_gp, LEGEND_POSITIONS.get(gp, "right top"))
        if gp == "plotmassdvsmass.gp":
            apply_unit_scale(dst_gp, SCALE_UNIT)

    # 2. 加载最新多源拟合数据 (对称性部分会自动补齐缺的 case)
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
        docs_dir / "thesis_V-A.csv",
        docs_dir / "thesis_T-X.csv",
        docs_dir / "thesis_S-PS.csv",
        docs_dir / "thesis_X-A.csv",
    ]

    generate_massvtem_csv(mass_df, docs_dir / "massvtem.csv")
    generate_deltamass_nt12_csv(delta_df, docs_dir / "deltamass_nt12.csv")
    generate_symmetry_txt_files(delta_df, docs_dir)
    generate_mdoutputre_csv(delta_df, docs_dir / "mdoutputre.csv")
    generate_thesis_table_csvs(delta_df, docs_dir)

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
