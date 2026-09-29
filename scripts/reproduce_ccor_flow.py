#!/usr/bin/env python3
"""
scripts/reproduce_ccor_flow.py
--------------------------------------------------------------------------------
对照 `ana/dat/ccor` 保存的数据, 并用**同一套流程**出图。

三个数据源 (--datasets):
  reference : 严格复刻 ana 的提取 (test1_lhadrons_*_mesons -> block -> abs -> 前 32 行
              -> 逐构型 Jackknife -> fold), 用来和 ccor 里保存的中间/最终结果对数;
  multisrc  : ParseLQCData 新算的多源关联函数 (16 个源位置平均, 保留符号, Ns 行);
  singlesrc : ParseLQCData 单源关联函数 (multi_src 文件里 0/0/0/0 源块)。

流程 (与 ana/dat/ccor/{32x10_*.py, simulate.py, dfiltered.py, pplot.py, newdata.py} 一致):
  fold C(x)<->C(Ns-x) -> 逐 Jackknife 样本 cosh 拟合 a*cosh(m(x-Ns/2)),
  窗口 [Ns/2-7, Ns/2+8), 丢掉 chi2/dof > 100 的样本
  -> 两信道逐样本质量作差 -> Jackknife 均值/误差
  -> data<对称性>.txt + massvtem.csv + gnuplot (plotmd.gp / plotmdre.gp / plotmassdvsmass.gp)

输出: output/ccor_flow/
  verification.csv                 # 与 ccor 保存数据的逐项对照 (max|diff|)
  comparison_mass.csv              # 三个数据源的介子质量对照
  comparison_delta.csv             # 三个数据源的质量差对照
  compare_delta_mass.png/.pdf      # 对照图: dM vs T (4 组对称性)
  compare_correlators.png/.pdf     # 对照图: 关联函数 (A / Xt)
  <dataset>/data*.txt|.pdf, massvtem.csv, massvtem_*.pdf, deltamassvsm.pdf

用法:
  .venv/bin/python scripts/reproduce_ccor_flow.py                 # 三源全跑 + 出图 + 对数
  .venv/bin/python scripts/reproduce_ccor_flow.py --datasets multisrc --no-gnuplot
  .venv/bin/python scripts/reproduce_ccor_flow.py --cases 32x12 32x14 32x16 36x18 40x16 48x18
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from docs.meson_scan_setup import SCAN_CASES, SOURCE_TAGS  # noqa: E402
from docs.physics_setup import DEFAULT_READIN_DIR, OUTPUT_ROOT  # noqa: E402
from src.parselqcdata import LEAST_SQUARES_FITTER  # noqa: E402
from src.parselqcdata.ccor_flow import (  # noqa: E402
    ANA_CHANNEL_BLOCKS,
    ana_channel_components,
    per_config_matrix_channels,
    ANA_CHI2_DOF_MAX,
    ANA_HALF_WINDOW,
    ANA_MLS,
    ANA_PAIRS,
    ANA_SCALE_MEV,
    ANA_TYPE_ID,
    ana_block_tag,
    ana_ensemble_dir,
    channel_masses_from_jk,
    fit_sample_masses,
    jackknife_mean_err,
    load_meson_scan_channel,
    natural_sorted,
    per_config_matrix,
    leave_one_out_jackknife,
    symmetrize_about_center,
    sample_errors,
    read_csv_matrix,
)

#: ccor 里保存的数据集 (intem = 12/14/16/18) 对应的本仓库 ensemble
CCOR_CASES: Tuple[str, ...] = ("32x12", "32x14", "32x16", "36x18")
#: 额外可加入的 ensemble (没有 ccor 对照)
EXTRA_CASES: Tuple[str, ...] = ("40x16", "48x18")

MESON_SCAN_ROOT: Path = OUTPUT_ROOT / "meson_scan" / "b4.17"
CCOR_DIR_DEFAULT: Path = Path("/Users/junxiongnie/code/ana/dat/ccor")
OUT_ROOT_DEFAULT: Path = OUTPUT_ROOT / "ccor_flow"
GNUPLOT_DIR: Path = PROJECT_ROOT / "scripts" / "gnuplot"

#: 数据源: reference = ana 逐位口径; ref_signed / ref_allcomp = 用于归因的对照口径
DATASETS: Tuple[str, ...] = ("reference", "multisrc", "singlesrc")
EXTRA_DATASETS: Tuple[str, ...] = ("ref_signed", "ref_allcomp")


# ------------------------------------------------------------------------------
# 任务构造
# ------------------------------------------------------------------------------
@dataclass(frozen=True)
class Ensemble:
    key: str
    order: int
    ns: int
    nt: int

    @property
    def mls(self) -> Tuple[str, ...]:
        return ANA_MLS

    def meson_case_dir(self, dataset: str, ml: str) -> Path:
        case = next(c for c in SCAN_CASES if c.key == self.key)
        return MESON_SCAN_ROOT / dataset / "cases" / case.output_name(float(ml))


def build_ensembles(case_keys: Sequence[str]) -> List[Ensemble]:
    by_key = {c.key: c for c in SCAN_CASES}
    out: List[Ensemble] = []
    for key in case_keys:
        case = by_key[key]
        out.append(Ensemble(key=key, order=case.order, ns=case.ns, nt=case.nt))
    return out


def ml_index(ml: str) -> int:
    """ana pplot.py: int(massterm / 0.0017)"""
    return int(float(ml) / 0.0017)


def temperature_mev(nt: int, scale: float = ANA_SCALE_MEV) -> float:
    return scale / float(nt)


# ------------------------------------------------------------------------------
# 单个 (dataset, ensemble, ml) 的 6 信道逐样本质量
# ------------------------------------------------------------------------------
def channel_masses_for_task(task: dict) -> dict:
    dataset = task["dataset"]
    ens = Ensemble(**task["ensemble"])
    ml = task["ml"]
    out: Dict[str, List[float]] = {}

    if dataset.startswith("reference") or dataset.startswith("ref_"):
        # ana 口径: 固定读 32 行 + 固定对称点 16 —— 只在 Ns=32 的 ensemble 上物理成立。
        # 对 Ns > 32 的 ensemble 这里仍然照做 (为了与 ccor 保存的数据逐位对齐),
        # 但必须显式告警, 避免把 32/16 当成通用逻辑套到别的尺寸上。
        n_lines = 32
        if ens.ns != n_lines:
            print(f"[WARN] {dataset}: ensemble {ens.key} 的 Ns={ens.ns} != 32, "
                  f"ana 参考口径 (前 32 行 + 对称点 16) 仅为复刻 ccor 而保留, 物理上不可信; "
                  f"请用 multisrc / singlesrc 数据源 (按 Ns/2 自动对称化)", flush=True)
        apply_abs = dataset != "ref_signed"
        average_all = dataset == "ref_allcomp"
        ens_dir = ana_ensemble_dir(Path(task["readin"]), ens.nt, ml)
        files = natural_sorted(str(p) for p in (ens_dir / "Output").glob("test1_lhadrons_*_mesons"))
        if not files:
            return {"dataset": dataset, "case": ens.key, "ml": ml, "masses": {}, "error": "no files"}
        for ch in ANA_TYPE_ID:
            mat = per_config_matrix_channels(
                files, ana_channel_components(ch), n_lines, apply_abs, average=average_all
            )
            if mat.size == 0:
                out[ch] = []
                continue
            jk = leave_one_out_jackknife(mat, n_rows=n_lines)
            sym = symmetrize_about_center(jk, center=n_lines // 2)
            err = sample_errors(sym)
            masses = fit_sample_masses(sym, err, center=n_lines // 2)
            out[ch] = [float(v) for v in masses]
    else:
        case_dir = ens.meson_case_dir(dataset, ml)
        for ch in ANA_TYPE_ID:
            jk = load_meson_scan_channel(case_dir, ch)
            if jk is None or jk.size == 0:
                out[ch] = []
                continue
            masses, _ = channel_masses_from_jk(jk, center=ens.ns // 2)
            out[ch] = [float(v) for v in masses]

    return {"dataset": dataset, "case": ens.key, "ml": ml, "masses": out, "error": None}


# ------------------------------------------------------------------------------
# 汇总: 质量表 / 质量差表 / ana 格式输出
# ------------------------------------------------------------------------------
def build_tables(records: List[dict], ensembles: Sequence[Ensemble]) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """返回 (质量表[长表], 质量差表[长表])"""
    mass_rows: List[dict] = []
    delta_rows: List[dict] = []

    by_key: Dict[Tuple[str, str, str], dict] = {
        (r["dataset"], r["case"], r["ml"]): r for r in records
    }

    for ens in ensembles:
        for ml in ens.mls:
            for dataset in sorted({r["dataset"] for r in records}):
                rec = by_key.get((dataset, ens.key, ml))
                if rec is None:
                    continue
                masses = {ch: np.array(v, dtype=float) for ch, v in rec["masses"].items()}
                for ch, arr in masses.items():
                    mean, err = jackknife_mean_err(arr)
                    mass_rows.append(
                        {
                            "dataset": dataset,
                            "case": ens.key,
                            "ns": ens.ns,
                            "nt": ens.nt,
                            "ml": float(ml),
                            "ml_str": ml,
                            "ml_index": ml_index(ml),
                            "tem_mev_ana": round(temperature_mev(ens.nt), 4),
                            "tem_mev_repo": round(2452.96 / ens.nt, 4),
                            "channel": ch,
                            "type_id": ANA_TYPE_ID[ch],
                            "n_samples": int(np.isfinite(arr).sum()),
                            "mass_lattice": mean,
                            "mass_err_lattice": err,
                            "mass_mev": mean * ANA_SCALE_MEV,
                            "mass_err_mev": err * ANA_SCALE_MEV,
                        }
                    )
                for ch1, ch2, title in ANA_PAIRS:
                    a1, a2 = masses.get(ch1), masses.get(ch2)
                    if a1 is None or a2 is None or a1.size == 0 or a2.size == 0:
                        continue
                    n = min(a1.size, a2.size)
                    diff = a1[:n] - a2[:n]
                    mean, err = jackknife_mean_err(diff)
                    delta_rows.append(
                        {
                            "dataset": dataset,
                            "case": ens.key,
                            "ns": ens.ns,
                            "nt": ens.nt,
                            "ml": float(ml),
                            "ml_str": ml,
                            "ml_index": ml_index(ml),
                            "tem_mev_ana": round(temperature_mev(ens.nt), 4),
                            "tem_mev_repo": round(2452.96 / ens.nt, 4),
                            "pair": title,
                            "ch1": ch1,
                            "ch2": ch2,
                            "n_samples": int(np.isfinite(diff).sum()),
                            "delta_mass_lattice": mean,
                            "delta_mass_err_lattice": err,
                            "delta_mass_mev": mean * ANA_SCALE_MEV,
                            "delta_mass_err_mev": err * ANA_SCALE_MEV,
                        }
                    )
    mass_df = pl.DataFrame(mass_rows).sort(["dataset", "nt", "ml", "type_id"]) if mass_rows else pl.DataFrame()
    delta_df = pl.DataFrame(delta_rows).sort(["dataset", "nt", "ml", "pair"]) if delta_rows else pl.DataFrame()
    return mass_df, delta_df


def write_ana_format(mass_df: pl.DataFrame, delta_df: pl.DataFrame, dataset: str, out_dir: Path) -> None:
    """写出 ana 同格式文件: massvtem.csv / data<对称性>.txt / mdoutputre.csv"""
    out_dir.mkdir(parents=True, exist_ok=True)
    sub = mass_df.filter(pl.col("dataset") == dataset)

    # --- massvtem.csv (newdata.py 格式: tem,type,massterm,massfit_mean_jk,massfit_err_jk) ---
    if sub.height:
        pl.DataFrame(
            {
                "tem": sub["tem_mev_ana"].to_list(),
                "type": sub["type_id"].to_list(),
                "massterm": sub["ml"].to_list(),
                "massfit_mean_jk": sub["mass_mev"].to_list(),
                "massfit_err_jk": sub["mass_err_mev"].to_list(),
            }
        ).write_csv(out_dir / "massvtem.csv")

    # --- massvtem_lattice.csv (同一流程, 格点单位 + 两种温度口径) ---
    if sub.height:
        sub.select(
            ["case", "ns", "nt", "ml", "ml_index", "tem_mev_ana", "tem_mev_repo",
             "channel", "type_id", "n_samples", "mass_lattice", "mass_err_lattice",
             "mass_mev", "mass_err_mev"]
        ).write_csv(out_dir / "massvtem_lattice.csv")

    # --- data<对称性>.txt (pplot.py 现版本格式: T, int(ml/0.0017), dM*2640, err*2640) ---
    dsub = delta_df.filter(pl.col("dataset") == dataset)
    for _, _, title in ANA_PAIRS:
        pair = dsub.filter(pl.col("pair") == title).sort(["tem_mev_ana", "ml"])
        if pair.height == 0:
            continue
        lines = [
            "{} {} {} {}".format(
                r["tem_mev_ana"], r["ml_index"], r["delta_mass_mev"], r["delta_mass_err_mev"]
            )
            for r in pair.iter_rows(named=True)
        ]
        (out_dir / f"data{title}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # --- mdoutputre.csv (plotmassdvsmass.gp 用的旧格式: pair_type, ml, dM, err, runningtime) ---
    pair_ids = {title: i + 1 for i, (_, _, title) in enumerate(ANA_PAIRS)}
    if dsub.height:
        rows = [
            f"{pair_ids[r['pair']]},{r['ml']},{r['delta_mass_lattice']},{r['delta_mass_err_lattice']},"
            f"Time taken: 0ms"
            for r in dsub.sort(["pair", "tem_mev_ana", "ml"]).iter_rows(named=True)
        ]
        (out_dir / "mdoutputre.csv").write_text(
            "intem,massterm,massdifference,error,runningtime\n" + "\n".join(rows) + "\n",
            encoding="utf-8",
        )


# ------------------------------------------------------------------------------
# gnuplot (与 ana 相同的 .gp 模板)
# ------------------------------------------------------------------------------
def run_gnuplot(dataset: str, out_dir: Path, pairs_titles: Sequence[str],
                ylabels: Dict[str, str]) -> List[str]:
    exe = shutil.which("gnuplot")
    if exe is None:
        return ["[WARN] 未找到 gnuplot, 跳过 PDF 出图 (apt/brew install gnuplot)"]

    logs: List[str] = []
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) data<对称性>.txt -> data<对称性>.pdf
    for title in pairs_titles:
        src = out_dir / f"data{title}.txt"
        if not src.exists():
            continue
        ytle = ylabels.get(title, " {/Symbol D}M(MeV)")
        gnuplot_cmd = [
            exe, "-e",
            f"input_fname=\"{src}\"; output_fname=\"{out_dir / f'data{title}.pdf'}\"; "
            f"symt=\"{title.replace('_', ' ')}\"; ytle=\"{ytle}\"",
            str(GNUPLOT_DIR / "plotmd.gp"),
        ]
        res = subprocess.run(gnuplot_cmd, capture_output=True, text=True)
        logs.append(f"[GNUPLOT] {dataset} {title}: rc={res.returncode} {res.stderr.strip()[:120]}")

    # 2) massvtem.csv -> massvtem_<ml>.pdf
    if (out_dir / "massvtem.csv").exists():
        res = subprocess.run(
            [exe, "-e", "input_fname='massvtem.csv'", str(GNUPLOT_DIR / "plotmdre.gp")],
            capture_output=True, text=True, cwd=str(out_dir),
        )
        logs.append(f"[GNUPLOT] {dataset} massvtem: rc={res.returncode} {res.stderr.strip()[:120]}")

    # 3) mdoutputre.csv -> deltamassvsm.pdf
    if (out_dir / "mdoutputre.csv").exists():
        res = subprocess.run(
            [exe, "-e", "input_fname='mdoutputre.csv'", str(GNUPLOT_DIR / "plotmassdvsmass.gp")],
            capture_output=True, text=True, cwd=str(out_dir),
        )
        logs.append(f"[GNUPLOT] {dataset} deltamassvsm: rc={res.returncode} {res.stderr.strip()[:120]}")

    return logs


# ------------------------------------------------------------------------------
# 与 ccor 保存数据对数
# ------------------------------------------------------------------------------
def verify_against_ccor(records: List[dict], readin: Path, ccor_dir: Path,
                        ensembles: Sequence[Ensemble]) -> List[dict]:
    """把 reference 数据源的结果与 ccor 目录里保存的中间/最终产物逐项对比"""
    rows: List[dict] = []
    ref = {(r["case"], r["ml"]): r for r in records if r["dataset"] == "reference" and not r["error"]}

    #: 逐项容差:
    #:   symdata  -> 1e-12 (逐位复刻)
    #:   dferr    -> 1e-7  (ana 的 sqrt(mean(v^2)-mean(v)^2) 存在 ~1e-9 量级的相消噪声, 无法逐位对齐)
    #:   fit_mass -> 1e-7  (lsqfit 优化器容差)
    #:   其余聚合量 -> 1e-6
    def tolerance(item: str) -> float:
        if item.startswith("symdata"):
            return 1e-12
        if item.startswith("dferr"):
            return 1e-7
        if item.startswith("fit_mass"):
            return 1e-7
        return 1e-6

    def add(item: str, detail: str, mine: float, saved: float, scale: Optional[float] = None,
            note: str = "") -> None:
        if not (np.isfinite(mine) and np.isfinite(saved)):
            rows.append({"item": item, "detail": detail, "mine": mine, "saved": saved,
                         "abs_diff": np.nan, "rel_diff": np.nan, "tolerance": tolerance(item),
                         "match": False, "note": note or "无有限值"})
            return
        diff = abs(mine - saved)
        ref_scale = abs(saved) if scale is None else max(abs(scale), 1e-300)
        rel = diff / ref_scale
        rows.append(
            {
                "item": item,
                "detail": detail,
                "mine": mine,
                "saved": saved,
                "abs_diff": diff,
                "rel_diff": rel,
                "tolerance": tolerance(item),
                "match": bool(rel <= tolerance(item)),
                "note": note,
            }
        )

    # 1) intem=18 (36x18) 的 symdata / dferr / 逐样本质量
    ens18 = next((e for e in ensembles if e.nt == 18 and e.key in CCOR_CASES), None)
    if ens18 is not None:
        files = natural_sorted(str(p) for p in (ana_ensemble_dir(readin, 18, "0.0020") / "Output").glob(
            "test1_lhadrons_*_mesons"))
        for ml in ANA_MLS:
            ml_files = natural_sorted(str(p) for p in (ana_ensemble_dir(readin, 18, ml) / "Output").glob(
                "test1_lhadrons_*_mesons"))
            if not ml_files:
                continue
            for ch, suffix in (("A", ""), ("Xt", "2")):
                sym_path = ccor_dir / f"symdatasample{suffix}{ml}.csv"
                if not sym_path.exists():
                    continue
                saved = read_csv_matrix(sym_path)
                mat = per_config_matrix(ml_files, ana_block_tag(ch), 32)
                jk = leave_one_out_jackknife(mat, n_rows=32)
                sym = symmetrize_about_center(jk, center=16)
                scale = float(np.max(np.abs(saved)))
                add(f"symdata[{ch}]", f"ml={ml} 矩阵 {sym.shape} (max|saved|={scale:.3e})",
                    float(np.max(np.abs(sym - saved))), 0.0, scale=scale)
                err_path = ccor_dir / f"dferr{1 if suffix == '' else 2}{ml}.csv"
                if err_path.exists():
                    saved_err = read_csv_matrix(err_path).ravel()
                    add(f"dferr[{ch}]", f"ml={ml}",
                        float(np.max(np.abs(sample_errors(sym) - saved_err))), 0.0,
                        scale=float(np.max(np.abs(saved_err))))

                rec = ref.get((ens18.key, ml))
                br_path = ccor_dir / f"binnedresult{suffix}{ml}.csv"
                if rec is not None and br_path.exists() and rec["masses"].get(ch):
                    saved_br = pl.read_csv(br_path)
                    mine = np.array(rec["masses"][ch], dtype=float)
                    saved_m = saved_br["massfit_mean"].to_numpy()
                    n = min(mine.size, saved_m.size)
                    finite = np.isfinite(mine[:n]) & np.isfinite(saved_m[:n])
                    add(
                        f"fit_mass[{ch}]",
                        f"ml={ml} 共 {int(finite.sum())}/{n} 个样本 (max|saved|={np.max(np.abs(saved_m)):.3e})",
                        float(np.max(np.abs(mine[:n][finite] - saved_m[:n][finite]))) if finite.any() else np.nan,
                        0.0,
                        scale=float(np.max(np.abs(saved_m))) if saved_m.size else 1.0,
                    )

    # 2) massvtem.csv (trial.sh 产物, 全部 16 个 ensemble x 6 信道)
    mv_path = ccor_dir / "massvtem.csv"
    if mv_path.exists():
        saved_mv = pl.read_csv(mv_path)
        # 参考侧自查: 同一 (tem, massterm) 下是否有两个 type 完全同值 (ana 多次运行残留的典型症状)
        dup_note: Dict[Tuple[float, float, int], str] = {}
        groups: Dict[Tuple[float, float], List[dict]] = {}
        for r in saved_mv.iter_rows(named=True):
            groups.setdefault((float(r["tem"]), float(r["massterm"])), []).append(r)
        for (tem, ml_v), items in groups.items():
            for i, a in enumerate(items):
                for b in items[i + 1:]:
                    va, vb = float(a["massfit_mean_jk"]), float(b["massfit_mean_jk"])
                    if abs(va - vb) <= 1e-9 * max(abs(va), 1.0):
                        msg = f"参考 massvtem 中 type{a['type']} 与 type{b['type']} 数值完全相同 (ana 多次运行残留), 该行不可信"
                        dup_note[(tem, ml_v, int(a["type"]))] = msg
                        dup_note[(tem, ml_v, int(b["type"]))] = msg

        for r in saved_mv.iter_rows(named=True):
            nt = int(round(ANA_SCALE_MEV / float(r["tem"])))
            case = {12: "32x12", 14: "32x14", 16: "32x16", 18: "36x18"}.get(nt)
            ch = next((k for k, v in ANA_TYPE_ID.items() if v == r["type"]), None)
            if case is None or ch is None:
                continue
            rec = ref.get((case, f"{r['massterm']:.4f}"))
            if rec is None or not rec["masses"].get(ch):
                continue
            mean, err = jackknife_mean_err(np.array(rec["masses"][ch], dtype=float))
            note = dup_note.get((float(r["tem"]), float(r["massterm"]), int(r["type"])), "")
            saved_val = float(r["massfit_mean_jk"])
            if not note and abs(mean * ANA_SCALE_MEV - saved_val) > 1e-6 * max(abs(saved_val), 1.0):
                # 参考值是否等于**另一个** type 的重算值 (ana 多次运行残留导致的标签错位)
                for ch_other, tid_other in ANA_TYPE_ID.items():
                    if ch_other == ch or not rec["masses"].get(ch_other):
                        continue
                    m_other, _ = jackknife_mean_err(np.array(rec["masses"][ch_other], dtype=float))
                    if abs(m_other * ANA_SCALE_MEV - saved_val) <= 1e-6 * max(abs(saved_val), 1.0):
                        note = (f"参考该行数值等于重算的 type{tid_other}({ch_other}) 结果 (标签错位, "
                                f"ana 多次运行残留), 该行不可信")
                        break
            add(f"massvtem[{ch}]", f"case={case} ml={r['massterm']:.4f} tem={r['tem']:g}", mean * ANA_SCALE_MEV,
                saved_val, note=note)
            add(f"massvtem_err[{ch}]", f"case={case} ml={r['massterm']:.4f} tem={r['tem']:g}", err * ANA_SCALE_MEV,
                float(r["massfit_err_jk"]), note=note)

    # 3) data<对称性>.txt (旧版 pplot 产物: 第 3 列是格点单位 dM)
    for _, _, title in ANA_PAIRS:
        path = ccor_dir / f"data{title}.txt"
        if not path.exists():
            continue
        raw_rows = [l.split() for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) < 4:
                continue
            tem, idx, dm, _err = float(parts[0]), int(parts[1]), float(parts[2]), float(parts[3])
            nt = int(round(ANA_SCALE_MEV / tem))
            case = {12: "32x12", 14: "32x14", 16: "32x16", 18: "36x18"}.get(nt)
            ml = [m for m in ANA_MLS if ml_index(m) == idx]
            if case is None or not ml:
                continue
            rec = ref.get((case, ml[0]))
            if rec is None:
                continue
            ch1, ch2, _ = next(p for p in ANA_PAIRS if p[2] == title)
            a1, a2 = rec["masses"].get(ch1), rec["masses"].get(ch2)
            if not a1 or not a2:
                continue
            n = min(len(a1), len(a2))
            mean, err = jackknife_mean_err(np.array(a1[:n]) - np.array(a2[:n]))
            # 质量差本身可能很小, 相对误差参照两个信道的质量量级
            scale = max(float(np.nanmean(np.abs(a1))), float(np.nanmean(np.abs(a2))))
            note = ""
            for other in raw_rows:
                if other is line.split():
                    continue
                if len(other) >= 4 and abs(float(other[2]) - dm) <= 1e-9 * max(abs(dm), 1e-12) \
                        and float(other[0]) != tem:
                    note = (f"参考 data{title}.txt 该行与 T={float(other[0]):g} 行数值完全相同 "
                            f"(ana 多次运行残留), 该行不可信")
                    break
            add(f"delta_mass[{title}]", f"case={case} ml={ml[0]} tem={tem:g} (格点单位)", mean, dm,
                scale=scale, note=note)

    return rows


# ------------------------------------------------------------------------------
# 关联函数数值对比 + 口径归因
# ------------------------------------------------------------------------------
def compare_correlators(mass_df: pl.DataFrame, ensembles: Sequence[Ensemble],
                        readin: Path, datasets: Sequence[str]) -> pl.DataFrame:
    """
    直接把两边的关联函数放在一起比: 参考 (ana 口径) vs 我的 meson_scan 数据。
    逐 (case, ml, channel) 给出折叠后曲线在拟合窗口内的相对差与 pull。
    """
    rows: List[dict] = []
    mine_ds = [d for d in ("multisrc", "singlesrc") if d in datasets]
    if "reference" not in datasets or not mine_ds:
        return pl.DataFrame()
    for ens in ensembles:
        if ens.key not in CCOR_CASES:
            continue
        ens_dir = ana_ensemble_dir(readin, ens.nt, "0.0020")
        files = natural_sorted(str(p) for p in (ens_dir / "Output").glob("test1_lhadrons_*_mesons"))
        for ml in ens.mls:
            ml_suffix = ana_ensemble_dir(readin, ens.nt, ml).name.split("m")[-1]
            files = natural_sorted(str(p) for p in (ana_ensemble_dir(readin, ens.nt, ml) / "Output").glob(
                "test1_lhadrons_*_mesons"))
            if not files:
                continue
            for ch in ANA_TYPE_ID:
                mat = per_config_matrix(files, ana_block_tag(ch), 32)
                ref_sym = symmetrize_about_center(leave_one_out_jackknife(mat, n_rows=32), center=16)
                ref_mean = np.array([jackknife_mean_err(ref_sym[i])[0] for i in range(32)])
                ref_err = sample_errors(ref_sym)
                xr = np.arange(32)
                for ds in mine_ds:
                    jk = load_meson_scan_channel(ens.meson_case_dir(ds, ml), ch)
                    if jk is None:
                        continue
                    sym = symmetrize_about_center(jk, center=ens.ns // 2)
                    mean = np.array([jackknife_mean_err(sym[i])[0] for i in range(sym.shape[0])])
                    err = sample_errors(sym)
                    # 两边统一到 "距对称点的偏移" 再比较
                    lo, hi = -ANA_HALF_WINDOW, ANA_HALF_WINDOW
                    off = np.arange(lo, hi + 1)
                    a = ref_mean[16 + off]
                    b = mean[ens.ns // 2 + off]
                    ae = ref_err[16 + off]
                    be = err[ens.ns // 2 + off]
                    finite = np.isfinite(a) & np.isfinite(b) & (np.abs(a) > 0)
                    # 只在 "信号点" 上做相对差统计 (窗口靠对称点的一侧信噪比可能 < 1)
                    signal = finite & (np.abs(a) > 2.0 * ae)
                    if not finite.any():
                        continue
                    rel_all = np.abs(np.abs(b[finite]) - np.abs(a[finite])) / np.abs(a[finite])
                    pull = np.abs(np.abs(b[finite]) - np.abs(a[finite])) / np.sqrt(ae[finite] ** 2 + be[finite] ** 2)
                    if signal.any():
                        rel = np.abs(np.abs(b[signal]) - np.abs(a[signal])) / np.abs(a[signal])
                        pull_sig = np.abs(np.abs(b[signal]) - np.abs(a[signal])) / np.sqrt(ae[signal] ** 2 + be[signal] ** 2)
                    else:
                        rel = rel_all
                        pull_sig = pull
                    sub = mass_df.filter((pl.col("case") == ens.key) & (pl.col("ml") == float(ml)) & (pl.col("channel") == ch))
                    m_ref = sub.filter(pl.col("dataset") == "reference")["mass_mev"]
                    m_mine = sub.filter(pl.col("dataset") == ds)["mass_mev"]
                    rows.append(
                        {
                            "case": ens.key,
                            "nt": ens.nt,
                            "ml": float(ml),
                            "channel": ch,
                            "mine_dataset": ds,
                            "n_points": int(finite.sum()),
                            "n_signal": int(signal.sum()),
                            "max_rel_diff_signal": float(rel.max()),
                            "median_rel_diff_signal": float(np.median(rel)),
                            "mean_pull_signal": float(pull_sig.mean()),
                            "max_pull": float(pull.max()),
                            "mass_ref_mev": float(m_ref[0]) if m_ref.len() else np.nan,
                            "mass_mine_mev": float(m_mine[0]) if m_mine.len() else np.nan,
                            "mass_sign_ref": float(np.sign(ref_mean[16])),
                            "mass_sign_mine": float(np.sign(mean[ens.ns // 2])),
                        }
                    )
    return pl.DataFrame(rows) if rows else pl.DataFrame()


def write_attribution(delta_df: pl.DataFrame, mass_df: pl.DataFrame, out_root: Path) -> Optional[Path]:
    """当 ref_signed / ref_allcomp 都在时, 输出口径归因表 (差值分解)"""
    present = set(mass_df["dataset"].to_list())
    if not {"reference", "ref_signed", "ref_allcomp"} <= present:
        return None
    piv = mass_df.pivot(on="dataset", index=["case", "ns", "nt", "ml", "ml_index", "channel"],
                        values="mass_mev", aggregate_function="first")
    piv = piv.with_columns(
        (pl.col("ref_signed") - pl.col("reference")).alias("d_signed_minus_abs"),
        (pl.col("ref_allcomp") - pl.col("reference")).alias("d_allcomp_minus_last"),
    )
    if "singlesrc" in present:
        piv = piv.with_columns(
            ((pl.col("singlesrc") - pl.col("ref_signed")) / pl.col("ref_signed")).alias("rel_mine_vs_refsigned"),
            ((pl.col("multisrc") - pl.col("singlesrc")) / pl.col("singlesrc")).alias("rel_multisrc_vs_singlesrc"),
        )
    piv = piv.with_columns(
        ((pl.col("ref_signed") - pl.col("reference")) / pl.col("reference")).alias("rel_abs_effect"),
        ((pl.col("ref_allcomp") - pl.col("reference")) / pl.col("reference")).alias("rel_component_effect"),
    ).sort(["nt", "ml", "channel"])
    path = out_root / "attribution.csv"
    piv.write_csv(path)
    return path


# ------------------------------------------------------------------------------
# matplotlib 对照图
# ------------------------------------------------------------------------------
def make_comparison_figures(mass_df: pl.DataFrame, delta_df: pl.DataFrame, out_dir: Path) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    datasets = [d for d in DATASETS if d in set(mass_df["dataset"].to_list())]
    colors = {"reference": "#444444", "multisrc": "#1f77b4", "singlesrc": "#d62728",
              "ref_signed": "#2ca02c", "ref_allcomp": "#9467bd"}
    markers = {"reference": "o", "multisrc": "s", "singlesrc": "^",
               "ref_signed": "v", "ref_allcomp": "D"}
    ml_colors = {"0.0020": "#69b3a2", "0.0035": "#ba68c8", "0.0070": "#64b5f6", "0.0120": "#ffb74d"}

    # ---- 图 1: dM vs T, 4 组对称性 ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    for ax, (_, _, title) in zip(axes.ravel(), ANA_PAIRS):
        sub = delta_df.filter(pl.col("pair") == title)
        for dataset in datasets:
            ds = sub.filter(pl.col("dataset") == dataset)
            for ml_str, color in ml_colors.items():
                rows = ds.filter(pl.col("ml_str") == ml_str).sort("tem_mev_ana")
                if rows.height == 0:
                    continue
                ax.errorbar(
                    rows["tem_mev_ana"].to_numpy(),
                    rows["delta_mass_mev"].to_numpy(),
                    yerr=rows["delta_mass_err_mev"].to_numpy(),
                    marker=markers[dataset],
                    color=color,
                    mfc="none" if dataset == "reference" else color,
                    ls="none",
                    ms=5,
                    capsize=2,
                    alpha=0.9 if dataset == "multisrc" else 0.6,
                )
        ax.axvspan(155, 158, color="gray", alpha=0.3)
        ax.axhline(0.0, color="black", lw=0.8)
        ax.set_title(title)
        ax.set_xlabel("temperature (MeV)")
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel(r"$\Delta M$ (MeV)")
    axes[1, 0].set_ylabel(r"$\Delta M$ (MeV)")
    handles = [
        plt.Line2D([], [], marker=markers[d], ls="none", color=colors[d], label=d) for d in datasets
    ] + [
        plt.Line2D([], [], marker="o", ls="none", color=c, label=f"ml={m}") for m, c in ml_colors.items()
    ]
    axes[0, 1].legend(handles=handles, fontsize=8, ncol=2)
    fig.suptitle("ccor flow: mass difference vs temperature (open = ana/ccor reference)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = out_dir / f"compare_delta_mass.{ext}"
        fig.savefig(p, dpi=140)
        written.append(p)
    plt.close(fig)

    # ---- 图 2: 关联函数对照 (每个 case 的 A / Xt) ----
    cases = sorted(set(delta_df["case"].to_list()))
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 2.2 * len(cases)), squeeze=False)
    for i, case in enumerate(cases):
        for j, ch in enumerate(("A", "Xt")):
            ax = axes[i][j]
            for dataset in datasets:
                curves = _load_curves_for_plot(case, ch, dataset)
                if curves is None:
                    continue
                x, y, e = curves
                ax.errorbar(x, y, yerr=e, marker=".", ls="none", ms=3, color=colors[dataset],
                            label=dataset, alpha=0.8)
            ax.set_yscale("log")
            ax.set_title(f"{case}  {ch} (ml=0.0020)")
            ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=8)
    fig.suptitle("correlator C(x) comparison (multi-source vs ana single-source reference)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = out_dir / f"compare_correlators.{ext}"
        fig.savefig(p, dpi=140)
        written.append(p)
    plt.close(fig)
    return written


def _load_curves_for_plot(case_key: str, channel: str, dataset: str):
    """取某个 ensemble/channel 的关联函数曲线 (mean, err) 用于对照图"""
    ens = next((e for e in build_ensembles(CCOR_CASES + EXTRA_CASES) if e.key == case_key), None)
    if ens is None:
        return None
    if dataset == "reference":
        try:
            ens_dir = ana_ensemble_dir(Path(DEFAULT_READIN_DIR), ens.nt, "0.0020")
        except ValueError:
            return None
        files = natural_sorted(str(p) for p in (ens_dir / "Output").glob("test1_lhadrons_*_mesons"))
        if not files:
            return None
        mat = per_config_matrix(files, ana_block_tag(channel), 32)
        if mat.size == 0:
            return None
        sym = symmetrize_about_center(leave_one_out_jackknife(mat, n_rows=32), center=16)
        x = np.arange(32)
    else:
        jk = load_meson_scan_channel(ens.meson_case_dir(dataset, "0.0020"), channel)
        if jk is None:
            return None
        sym = symmetrize_about_center(jk, center=ens.ns // 2)
        x = np.arange(ens.ns)
    mean = np.array([jackknife_mean_err(sym[i])[0] for i in range(sym.shape[0])])
    err = sample_errors(sym)
    return x, mean, err


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="复刻 ana/dat/ccor 流程: 对照数据并同流程出图",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--datasets", nargs="*", default=list(DATASETS),
                   choices=list(DATASETS) + list(EXTRA_DATASETS),
                   help="要跑的数据源: reference / multisrc / singlesrc (+ 归因用 ref_signed / ref_allcomp)")
    p.add_argument("--cases", nargs="*", default=list(CCOR_CASES),
                   help=f"ensemble 列表 (默认 {list(CCOR_CASES)}; 可加 {list(EXTRA_CASES)})")
    p.add_argument("--readin", type=Path, default=Path(DEFAULT_READIN_DIR), help="readin 数据根目录")
    p.add_argument("--ccor-dir", type=Path, default=CCOR_DIR_DEFAULT, help="ana ccor 参考数据目录")
    p.add_argument("--out-root", type=Path, default=OUT_ROOT_DEFAULT, help="输出根目录")
    p.add_argument("--workers", type=int, default=6, help="并行进程数")
    p.add_argument("--gnuplot", dest="gnuplot", action="store_true", default=True, help="调用 gnuplot 出 PDF")
    p.add_argument("--no-gnuplot", dest="gnuplot", action="store_false", help="跳过 gnuplot")
    p.add_argument("--no-verify", dest="verify", action="store_false", default=True, help="跳过与 ccor 的对数")
    p.add_argument("--no-figures", dest="figures", action="store_false", default=True, help="跳过 matplotlib 对照图")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    t0 = time.time()
    ensembles = build_ensembles(args.cases)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"[CCOR] datasets={args.datasets} ensembles={[e.key for e in ensembles]}")
    print(f"[CCOR] 参考数据目录: {args.ccor_dir}")

    tasks = [
        {"dataset": ds, "ensemble": {"key": e.key, "order": e.order, "ns": e.ns, "nt": e.nt},
         "ml": ml, "readin": str(args.readin)}
        for ds in args.datasets
        for e in ensembles
        for ml in e.mls
    ]
    print(f"[CCOR] 共 {len(tasks)} 个任务 (dataset x ensemble x ml)")

    records: List[dict] = []
    done = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = [ex.submit(channel_masses_for_task, t) for t in tasks]
        for fut in as_completed(futs):
            rec = fut.result()
            records.append(rec)
            done += 1
            status = rec["error"] or f"{len([v for v in rec['masses'].values() if v])}/6 信道"
            print(f"[{done:>3}/{len(tasks)}] {rec['dataset']:<10} {rec['case']:<7} ml={rec['ml']} -> {status}",
                  flush=True)

    mass_df, delta_df = build_tables(records, ensembles)
    print(f"[CCOR] 质量表 {mass_df.height} 行, 质量差表 {delta_df.height} 行")

    # 口径留档 (可复现性): 拟合器 / 窗口 / 过滤阈值 / 标度
    import json as _json

    (out_root / "flow_config.json").write_text(
        _json.dumps(
            {
                "datasets": list(args.datasets),
                "cases": [e.key for e in ensembles],
                "fitter": LEAST_SQUARES_FITTER,
                "fitter_note": "强制 chi2 最小二乘, 不使用 lsqfit 的环境默认 fitter",
                "cosh_center": "Ns/2 (我的数据源) / 16 与 32 行 (ana reference 数据源)",
                "window": f"[center-{ANA_HALF_WINDOW}, center+{ANA_HALF_WINDOW}]",
                "chi2_dof_max": ANA_CHI2_DOF_MAX,
                "scale_mev": ANA_SCALE_MEV,
                "jackknife_error": "sqrt((n-1)*mean((x-mean)^2))",
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    # 逐数据集落盘 + gnuplot
    ylabels = {
        "SU(2)XSU(2) V-A": " {/Symbol D}M_{V-A}(MeV)",
        "U(1)A T-X": " {/Symbol D}M_{T-X}(MeV)",
        "U(1) S-PS": " {/Symbol D}M_{S-PS}(MeV)",
        "SU(2)spinxchiral X-A": " {/Symbol D}M_{X-A}(MeV)",
    }
    for ds in args.datasets:
        ds_dir = out_root / ds
        write_ana_format(mass_df, delta_df, ds, ds_dir)
        if args.gnuplot:
            for line in run_gnuplot(ds, ds_dir, [t for _, _, t in ANA_PAIRS], ylabels):
                print(line)
    mass_df.write_csv(out_root / "meson_mass_all.csv")
    delta_df.write_csv(out_root / "delta_mass_all.csv")

    # 关联函数直接对比 + 口径归因
    if "reference" in args.datasets and any(d in args.datasets for d in ("multisrc", "singlesrc")):
        corr_df = compare_correlators(mass_df, ensembles, Path(args.readin), args.datasets)
        if corr_df.height:
            corr_df.write_csv(out_root / "comparison_correlator.csv")
            print(f"[CCOR] 关联函数对比: {out_root / 'comparison_correlator.csv'} "
                  f"(信号点中位相对差 {corr_df['median_rel_diff_signal'].median():.2%}, "
                  f"信号点平均 pull {corr_df['mean_pull_signal'].mean():.2f})")
    attr = write_attribution(delta_df, mass_df, out_root)
    if attr is not None:
        print(f"[CCOR] 口径归因表: {attr}")

    # 多数据源对照
    if len(args.datasets) > 1:
        pivot = delta_df.pivot(
            on="dataset", index=["case", "nt", "ml", "pair"], values="delta_mass_mev", aggregate_function="first"
        )
        pivot.write_csv(out_root / "comparison_delta.csv")
        pivot_m = mass_df.pivot(
            on="dataset", index=["case", "nt", "ml", "channel"], values="mass_mev", aggregate_function="first"
        )
        pivot_m.write_csv(out_root / "comparison_mass.csv")
        print(f"[CCOR] 对照表: {out_root / 'comparison_delta.csv'}, {out_root / 'comparison_mass.csv'}")

    # 与 ccor 保存数据对数
    if args.verify and "reference" in args.datasets:
        rows = verify_against_ccor(records, Path(args.readin), Path(args.ccor_dir), ensembles)
        if rows:
            vdf = pl.DataFrame(rows)
            vdf.write_csv(out_root / "verification.csv")
            n_ok = int(vdf["match"].sum())
            bad = vdf.filter(~pl.col("match"))
            n_ref_bug = int(bad.filter(pl.col("note") != "").height)
            print(f"[CCOR] 与 ccor 保存数据对数: {n_ok}/{vdf.height} 项一致, "
                  f"{n_ref_bug} 项为参考数据自身缺陷 (verification.csv)")
            for r in bad.sort("abs_diff", descending=True).head(6).iter_rows(named=True):
                print(f"        {r['item']:<22} {r['detail']:<40} |diff|={r['abs_diff']:.3e} "
                      f"{('<- ' + r['note']) if r['note'] else '<- 需检查'}")

    if args.figures and mass_df.height and delta_df.height:
        figs = make_comparison_figures(mass_df, delta_df, out_root)
        for p in figs:
            print(f"[CCOR] 对照图: {p}")

    print(f"[DONE] ccor 流程复刻完成, 用时 {time.time() - t0:.1f}s -> {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
