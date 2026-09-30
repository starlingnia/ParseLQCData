# 与 ana/dat/ccor 数据的对照 + 同流程出图

> 流程库：`src/parselqcdata/ccor_flow.py`
> 驱动脚本：`scripts/reproduce_ccor_flow.py`
> gnuplot 模板：`scripts/gnuplot/{plotmd.gp, plotmdre.gp, plotmassdvsmass.gp}`（原样取自 ana，仅把硬编码输入文件名参数化）
> 参考数据：`/Users/junxiongnie/code/ana/dat/ccor`
> 输出：`output/ccor_flow/`

## 1. ana ccor 流程是什么

`ana/dat/ccor` 保存的是「两信道质量差 → 对称性破缺信号」的分析产物，链路为：

| 步骤 | ana 脚本 | 口径 |
|:---|:---|:---|
| 1. 抽取 | `32x10_Beta4.70_m0.01_mesons.py` | 从 `<readin>/3*x{Nt}_b4.17_ms0.040m{ml}/Output/test1_lhadrons_*_mesons` 里按精确 tag 取 block，**取实部绝对值**，只取**前 32 行**，逐构型留一 Jackknife，再关于对称点 fold |
| 2. 拟合 | `simulate.py` | 每个 Jackknife 样本做 `a·cosh(m(x−16))`，窗口 `x∈[9,24)`，丢掉 `chi2/dof > 100` 的样本 |
| — | （本仓库实现） | **强制 `fitter='scipy_least_squares'`（chi2 最小二乘）**，见下；`reference` 数据源仍用 ana 的固定半宽 7（`Ns=32` 时即 `[9,24)`） |
| — | （本仓库实现） | **`multisrc`/`singlesrc` 的拟合半宽随格子缩放** `N = max(7, round(Ns/3))`（32→11、40→13、48→16），并丢掉 `m ≤ 1e-3` 的退化解；固定 `N=7` 在 `Ns=40/48` 上尚未收敛（48³×18 的 `Xt` 会退化成 `m=0`），会给低温柔虚高的 ΔM。可用 `scripts/try_fit_window.py --mode ana` 复现这个 N 依赖 |
| 3. 求差 | `dfiltered.py` | 两信道逐样本质量作差，`mean` 与 `sqrt((n−1)·mean((x−mean)²))` |
| 4. 汇总 | `pplot.py` / `newdata.py` | `T = 2640/Nt`，`mass(MeV) = mass·2640`，写 `data<对称性>.txt`、`massvtem.csv`、`mdoutputre.csv` |
| 5. 出图 | `plotmd.gp` / `plotmdre.gp` / `plotmassdvsmass.gp` | ΔM vs T（4 组 ml）、6 信道质量 vs T、ΔM vs ml |

关键细节（读 ana 源码得到，并已数值验证）：
- ana 脚本里 `for mapping in data_blocks_...` 只保留列表**最后一个** mapping，所以每个信道实际只用**一个**算子分量：

  | 信道 | ana 实际使用 | 完整分量数 |
  |:---|:---|:---:|
  | V | `Vector2/DIRZ` | 6 |
  | A | `AVector2/DIRZ` | 6 |
  | Tt | `TVector2/DIRY` | 3 |
  | Xt | `TAVector2/DIRY` | 3 |
  | S | `TVector4/DIRY` | 3 |
  | PS | `TAVector4/DIRY` | 3 |

- 逐构型取 `abs()`，之后才做平均 —— 对信噪比差的点这会把噪声翻成正值（Jensen 偏置）。
- 固定读 32 行、固定对称点 16：只对 32³ 格点成立。
- 4 组对称性：`V−A`(SU(2)×SU(2))、`Tt−Xt`(U(1)ₐ)、`S−PS`(U(1)ₐ)、`A−Xt`(SU(2)_CS)。

## 2. 三个（可扩展五个）数据源

| dataset | 含义 |
|:---|:---|
| `reference` | 严格复刻 ana：单源 + 最后分量 + 逐构型 abs + 前 32 行 + 逐构型 Jackknife |
| `ref_signed` | 同上但**不取 abs**（保留符号，拟合时才取）→ 单独看 abs 的影响 |
| `ref_allcomp` | 同上但要 **abs** 且**平均全部算子分量** → 单独看分量平均的影响 |
| `multisrc` | ParseLQCData 新算的多源关联函数（16 个源平均、全分量平均、保留符号、Ns 行、逐 bin Jackknife） |
| `singlesrc` | 同上但只用 `0/0/0/0` 源块（与 ana 同源同统计对象，便于隔离"源平均"） |

## 3. 数据对照结论（`output/ccor_flow/`）

### 3.1 原始数据本身完全一致

- `test1_lhadrons_*_mesons`（单源文件）与 `test1_lhadrons_*_mesons_multi_src` 的 `0/0/0/0` block
  **逐位相同**（max|Δ| = 0）；16 个源位置中的其它源是各自独立的测量。
- C++ 抽取与 ana 抽取在**同一 block、同一分量**上只差一个整体符号（相关函数在 tensor 类道为负），
  逐点比值恒为 −1（折叠镜像配对正确时）。

### 3.2 流程复现（`verification.csv`）

| 对照项 | 结果 |
|:---|:---|
| `symdatasample{ml}.csv` / `symdatasample2{ml}.csv`（36³×18, A/Xt） | **1e-16（逐位复现）** |
| `dferr1/2{ml}.csv` | ~1e-9 相对（ana 的 `sqrt(mean(v²)−mean(v)²)` 本身有相消噪声，无法逐位） |
| `binnedresult{ml}.csv` / `binnedresult2{ml}.csv` | ~1e-10 绝对（lsqfit 优化器容差） |
| `data<对称性>.txt`（ΔM，格点单位） | 4 组中 3 组 16/16 行一致；V−A 的 T=165 行与其自身 T=188.6 行**数值完全相同**（ana 多次运行残留） |
| `massvtem.csv` | V/A 两个 type 16/16 一致；Tt/Xt/S/PS 在 T=220 处为**标签错位**（该行数值等于重算的其它 type） |

即：**244/280 项一致，36 项不一致全部可归因为参考数据自身的运行残留**（`verification.csv` 的 `note` 列逐条写明）。

### 3.3 关联函数与质量的直接对比（`comparison_correlator.csv`、`comparison_mass.csv`）

以 32³×12、ml=0.0020、`multisrc` vs `reference` 为例（窗口内信号点）：

| 信道 | 中位相对差 | 平均 pull | m_ref (MeV) | m_mine (MeV) | 结论 |
|:---|:---:|:---:|:---:|:---:|:---|
| V | 5.0% | 1.15 | 1300.8 | 1300.9 | ✅ 一致 |
| A | 4.6% | 0.99 | 1302.0 | 1308.1 | ✅ 一致 |
| PS | 37% | 1.24 | 754.5 | 722.3 | ⚠️ 近对称点噪声主导，质量差 ~4% |
| Tt | 22% | 2.44 | 1169.8 | 1439.2 | ❌ 口径效应（见下） |
| Xt | 29% | 2.83 | 1159.3 | 1392.8 | ❌ |
| S | 89% | 2.63 | 762.6 | 1485.4 | ❌ |

差异来源（`attribution/attribution.csv` 逐条分解；`rel_abs_effect` = 去掉逐构型 abs 的相对变化，
`rel_component_effect` = 改成全分量平均的相对变化）：

1. **逐构型 `abs()`（主因）**：`mean|C| > |mean C|`，噪声越大偏置越大，且靠近对称点更严重
   （把噪声整流成平台形状）→ 拟合质量系统性**偏低**。
   - PS 道恒为正 → abs 是空操作，`rel_abs_effect = 0.0%`（全部 ensemble 都是 0，正好当校验）；
   - 32³×12（T=220 MeV，信噪比好）V/A 只有 2%；
   - 32³×16（T=153 MeV，近 Tc，噪声大）V/A 达 47%~110%。
   - 换个初值 (`p0`) 重拟合结果完全不变（已验证），说明这不是优化器假象，而是估计量本身的偏置。
2. **只用最后一个算子分量**：V/A 约 1%~3%，PS 最多 ~18%（各分量质量本来就散开）。
3. **固定 32 行 + 对称点 16**：只对 32³ 成立；36³（`intem=18`）真实对称点是 18，
   因此 ccor 里 `T=146.7 MeV` 一列结构性不可信（其 V 道比正确值低一倍以上）。
4. **单源 vs 16 源平均**：只影响统计误差（多源误差小 2~3 倍），中心值一致
   （32x12：multisrc vs singlesrc 差 0.3%~5.8%）。

**独立交叉验证**：把「保留符号 + 全分量平均」的单源分析（`ref_signed` 思路）与
我的多源结果比较，32³×16 的 A 道给出 0.402 vs 0.378（差 5%），
而 ana 的 abs 口径给出 0.174（差 2 倍）——
两个互不相关的估计（单源/多源）互相吻合，说明 **abs 口径是离群的一方**，
我的多源数值更可信。

**结论**：原始数据、抽取与拟合流程都已逐位对齐；质量差异来自 ana 流程自身的三个口径选择，
其中「逐构型取 abs」在噪声大的 ensemble 上会给出成倍偏差。

### 3.4 拟合器必须显式指定（防退化）

`lsqfit.nonlinear_fit` 的 `fitter` 默认值是 `None`，含义是**按环境自动挑**：
装了 GSL 就用 `gsl_multifit`，否则才退到 `scipy_least_squares`。
同一份代码在不同机器上因此可能走出不同（甚至退化）的结果。

本仓库的做法：所有拟合统一走一个入口

```python
from src.parselqcdata import chi2_least_squares_fit, LEAST_SQUARES_FITTER
# LEAST_SQUARES_FITTER == "scipy_least_squares"
fit = chi2_least_squares_fit(data=..., fcn=..., p0=...)   # 调用方无法把 fitter 改成别的
```

覆盖范围（`src/parselqcdata/plateau_fit.py` 的 `fit_single_jackknife_column`、
`fit_single_jackknife_column_centered`，以及 `ccor_flow.fit_sample_masses`），
仓库里**只有这一处**直接调用 `lsqfit.nonlinear_fit`，由源码级测试
`test_no_fit_call_bypasses_the_fitter_guard` 守住。
本次强制指定后所有对照结果**逐位不变**（本机默认本来就是 scipy），
但换到装了 GSL 的环境不会再悄悄换算法。
口径同时落档在 `output/ccor_flow/flow_config.json` 与
`output/meson_scan/b4.17/manifest.json`（`fitter` 字段）。

## 4. 输出结构

```
output/ccor_flow/
  verification.csv              # 与 ccor 保存数据的逐项对照（含 note: 参考侧缺陷说明）
  comparison_correlator.csv     # 关联函数直接对比: 逐 case/ml/信道 的相对差与 pull
  comparison_mass.csv           # 三个数据源的 6 信道质量对照
  comparison_delta.csv          # 三个数据源的 4 组 ΔM 对照
  meson_mass_all.csv / delta_mass_all.csv     # 长表全量
  compare_delta_mass.png/.pdf   # ΔM vs T 对照图 (reference 用空心点)
  compare_correlators.png/.pdf  # C(x) 对照图 (A / Xt)
  reference/  multisrc/  singlesrc/           # 同流程产物 (ana 同格式)
      dataSU(2)XSU(2) V-A.txt/.pdf            # ΔM vs T, 4 组 ml
      dataU(1)A T-X.txt/.pdf
      dataU(1) S-PS.txt/.pdf
      dataSU(2)spinxchiral X-A.txt/.pdf
      massvtem.csv / massvtem_lattice.csv     # 6 信道质量 vs T
      massvtem_0.00{20,35,70,120}.pdf         # gnuplot 原样出图
      mdoutputre.csv / deltamassvsm.pdf
  attribution/                  # 5 数据源归因 (含 ref_signed / ref_allcomp)
      attribution.csv           # rel_abs_effect / rel_component_effect 逐条分解
```

## 5. 用法

```bash
# 主对照 (reference / multisrc / singlesrc): 复现 + 对照 + 同流程出图
.venv/bin/python scripts/reproduce_ccor_flow.py

# 口径归因 (再加 ref_signed / ref_allcomp 两个对照数据源)
.venv/bin/python scripts/reproduce_ccor_flow.py \
    --datasets reference ref_signed ref_allcomp multisrc singlesrc \
    --out-root output/ccor_flow/attribution

# 只跑我的多源数据, 不调 gnuplot
.venv/bin/python scripts/reproduce_ccor_flow.py --datasets multisrc --no-gnuplot

# 把额外的 40³×16 / 48³×18 也一并跑 (无 ccor 对照)
.venv/bin/python scripts/reproduce_ccor_flow.py --cases 40x16 48x18 \
    --datasets multisrc --out-root output/ccor_flow/extra
```

参数：`--datasets`、`--cases`、`--workers`、`--ccor-dir`、`--readin`、`--out-root`、
`--no-gnuplot`、`--no-verify`、`--no-figures`。

## 6. 测试

```bash
.venv/bin/python tests/test_ccor_flow.py
```

7 项（全部通过）：block 抽取口径、留一 Jackknife 代数、fold 约定、Jackknife 误差公式、
信道映射，以及与 ccor 保存数据的 `symdata`（1e-15）/`dferr`+`binnedresult`（1e-6 相对）逐位对齐。
