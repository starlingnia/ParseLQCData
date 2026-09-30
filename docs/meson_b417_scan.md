# b4.17 有限温度扫描任务：关联函数 + 介子质量

> 任务脚本：`scripts/reproduce_meson_b417_nt_scan.py`
> 一键运行：`bash scripts/run_meson_b417_nt_scan.sh`（或 `--dry-run` 先看清单）
> 设定单一事实来源：`docs/meson_scan_setup.py`
> 结果目录：`output/meson_scan/b4.17/`

## 1. 任务范围

读取 `data/readin/` 下**数字开头、含 b4.17 参数**的目录，共 6 组格点 × 4 组
light quark mass = 24 个 ensemble（`ms = 0.040`）：

| 序号 | ensemble | 格点目录 | ml | T_ref (MeV) |
|:---:|:---:|:---|:---|:---:|
| 00 | 32³×12 | `32x12_b4.17_ms0.040m0.00XX` | 0.0020 / 0.0035 / 0.0070 / 0.0120 | 204.00 |
| 01 | 32³×14 | `32x14_b4.17_ms0.040m0.00XX` | 同上 | 174.86 |
| 02 | 32³×16 | `32x16_b4.17_ms0.040m0.00XX` | 同上 | 153.00 |
| 03 | 36³×18 | `36x18_b4.17_ms0.040m0.00XX` | 同上 | 136.00 |
| 04 | 40³×16 | `40x16_b4.17_ms0.040m0.00XX` | 同上 | 153.00 |
| 05 | 48³×18 | `48x18_b4.17_ms0.040m0.00XX` | 同上 | 136.00 |

- `T_ref` 为用户给定表格中的参考温度（等价于 `a⁻¹ = 2453 MeV`）。
- 汇总表里另给仓库统一口径的温度 `temperature_mev = T(β=4.17, Nt=16) × 16 / Nt`
  （`TEMP_MAP['4.17'] = 153.31 MeV` ⇒ `a⁻¹ = 2452.96 MeV`），论文出图统一取 `a⁻¹ = 2453 MeV`。
- 数据目录里另外两个数字开头且含 b4.17 的目录
  （`16-z-32x12_b4.17_ms0.040m0.0020`、`48x16b4.17`）**不在任务表内**，
  会在 `manifest.json` 的 `readin_dirs_b417_not_in_task` 中登记并跳过
  （`48x16b4.17` 已由既有 `reproduce_meson_multi.py` 管道覆盖）。

## 2. 计算链

1. **关联函数抽取**：`tools/meson_orchestrator.py` → C++ 核心库
   `run_meson_pipeline_c_api`，扫描 `<ensemble>/Output/test1_lhadrons_*_mesons_multi_src`。
   每个文件含 16 个源位置（4×4 空间网格）× 每个信道 4 个 block：
   `temporal`(Nt 行) / `spatial:DIRX, DIRY, DIRZ`(各 Ns 行)。
   6 个狄拉克双线性信道 `AV, S, Tt, PS, Xt, Vec` 只使用 `spatial:DIR*`，
   因此 **`num_lines = Ns`**，得到的是沿极化方向的**空间关联函数 `C(x)`, x = 0…Ns-1**
   （16 个源位置环形移位平均后按 bin 折叠成 Jackknife 样本矩阵）。
2. **有效质量**：解
   `C(x)/C(x+1) = cosh(m(x−Ns/2)) / cosh(m(x+1−Ns/2))`，
   对称点取 **Ns/2**（32→16, 36→18, 40→20, 48→24），
   与既有 4.17@48³×16 分析（写死 24）在 Ns=48 时完全等价。
   先用向量化 Newton 求解，未收敛/无实根的点回退到标量 `fsolve`，
   并对残差做校验（无实根 → NaN，不再输出伪解）。
3. **介子质量**：单参数 cosh 平台拟合 `f(x) = a·cosh(m(x−Ns/2))`，
   **拟合器强制 `fitter='scipy_least_squares'`（chi2 最小二乘，统一入口
   `chi2_least_squares_fit`，不依赖 lsqfit 的环境默认，避免在装了 GSL 的机器上退化）**，
   对每个 Jackknife 样本列独立拟合，按仓库既有约定汇总：
   `mean = ⟨m_j⟩`，`err = sqrt((J−1)·Σ(m_j−mean)²/J)`。
4. **落盘顺序**：目录名前缀 = 任务表序号；汇总表按
   `(order, ml_index, channel_index)` 排序，与用户表格逐行对应。

## 3. 窗口策略（三种结果同时给出，`--window-policy` 决定主结果）

| 策略 | 说明 | 适用性 |
|:---|:---|:---|
| `auto`（默认） | 用有效质量平坦性自动检测平台：窗口 `[s, Ns/2)` 内所有相邻差分 `meff(x+1)−meff(x)` 在 Jackknife 误差下与 0 一致的 `chi2/dof ≤ --plateau-chi2`（默认 1.5），取满足条件的最长窗口，再在 `[s, Ns/2+1)` 做 cosh 拟合 | 各格点规模都可用，推荐 |
| `mirror` | 把 4.17@48³×16 既有 slice 表按"距对称点距离"平移到本 ensemble（起点夹紧到 ≥1） | 仅 Ns=48 时等价于既有发布口径；32³ 上会落到接触项区域，**不可用** |
| `scan` | 对均值关联函数做 cosh 拟合，取 `chi2/dof` 最小的窗口 | 相关性强的数据上 chi2/dof 区分度差，仅作对照 |

`auto` 的平坦性检验用的是**逐 bin 差分的 Jackknife 误差**（而不是把各点误差当独立），
因此对点间强相关不过度乐观；阈值可用 `--plateau-chi2` 调紧/调松，
`mass_scan_<ch>.csv` 与 `meson_mass_pivot_scan.csv` 保留了全部窗口结果，
方便事后改判而不必重算。

## 4. 输出结构（`output/meson_scan/b4.17/`）

```
manifest.json                      # 任务表、目录解析、约定、每个 ensemble 的状态
multisrc/                          # 多源（默认主结果，已与既有 48³×16 发布值逐位复核）
  meson_mass_summary.csv           # 长表: 有序 + 全部列（主结果 / auto / mirror / scan）
  meson_mass_pivot.csv             # 宽表: 每 ensemble 一行, 每信道一列 "mass(err×10⁵)"
  meson_mass_pivot_auto.csv        # 同上, 自动平台窗口
  meson_mass_pivot_mirror.csv      # 同上, 镜像窗口（与既有发布口径对照）
  meson_mass_pivot_scan.csv        # 同上, chi2/dof 选优窗口
  cases/00_32x12_ms0.040_m0.0020/
    correlators.csv                # 6 信道关联函数 C(x): channel, x, mean, err
    corr_<ch>.csv                  # 单信道关联函数
    jk_<ch>.csv                    # Jackknife 折叠矩阵 (Ns × n_bins, 无表头)
    meff.csv / meff_<ch>.csv       # 有效质量曲线 (无实根的点为 NaN)
    mass.csv / mass_<ch>.csv       # 该 ensemble 的介子质量（主结果 + 三个窗口口径）
    mass_scan.csv / mass_scan_<ch>.csv  # 全窗口扫描明细 (含 cosh 拟合 chi2/dof)
    meta.json                      # 输入目录 / bins / 窗口 / 耗时
singlesrc/                         # 单源（同结构；单源取 multi_src 文件里 0/0/0/0 源块，与既有 single 管道一致）
```

> 宽表里 `<ch>_chi2dof` 的含义：`meson_mass_pivot.csv` / `_mirror` / `_scan` 给的是该窗口
> **cosh 平台拟合**的 chi2/dof；`meson_mass_pivot_auto.csv` 给的是平台**平坦性检验**的
> chi2/dof（相邻有效质量差分与 0 的一致性）。

## 5. 用法

```bash
# 只列清单（不计算）: 检查 24 个目录是否齐全
.venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --dry-run

# 全量（多源 + 单源，24 ensemble × 2 = 48 个任务，4 进程）
bash scripts/run_meson_b417_nt_scan.sh

# 只跑多源、指定 ensemble / ml
.venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --source multi \
    --cases 48x18 32x12 --ml 0.0020 0.0120 --workers 4

# 收紧平台判据（窗口更靠近对称点、统计误差更大）
.venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --plateau-chi2 0.8

# 主结果改用镜像窗口口径
.venv/bin/python scripts/reproduce_meson_b417_nt_scan.py --window-policy mirror
```

- `--symmetrize`：先做 `C(x) ↔ C(Ns−x)` 对称化再拟合（默认关闭，保持与既有管道一致）。
- 结果可续算：已完成的 ensemble 命中 `meta.json` 直接复用，`--force` 强制重算。

## 6. 数值校验

`tests/test_meson_b417_scan.py`（10 项，全部通过）：

- 任务表、目录名、温度、输出目录顺序与用户表格逐行一致；
- Ns=48 时镜像窗口必须原样复现既有 slice 表 `[14,25)/[8,25)/[10,25)`；
- 对真实 48³×16 数据，新接口必须复现既有发布结果：
  `output/ratio_results/b4.17/meff_PS.csv`（逐点 ≤ 1e-15）
  与 `output/simulateresult/b4.17/summary_fit_PS.csv`
  （`m_PS = 0.0835528531 ± 0.0013779298`）；
- 平台自动检测：平坦序列检出长窗口，含污染漂移的序列把窗口起点后移。

```bash
.venv/bin/python tests/test_meson_b417_scan.py
```

## 7. 注意事项

- 32³ 系列（Nt = 12/14/16）在高温下空间关联函数的平台区很窄，
  `auto` 检出的窗口可能仍带少量激发态污染；请对照 `meff.csv` 与 `mass_scan.csv`
  判断，必要时用 `--plateau-chi2` 收紧窗口。
- `S`/`Xt` 道在部分 ensemble 上信噪比差（关联函数在对称点附近过零），
  Jackknife 误差可能很大甚至窗口检不出平台（此时主结果自动回退到 `scan`，
  并在 `window_mode` 中标注 `no_plateau`）。
- `mass_mev*` 列由 `a⁻¹ = 2453 MeV` 换算（仓库 `TEMP_MAP['4.17'] = 153.31 MeV`
  ⇒ `a⁻¹ = 2452.96 MeV`，取整为 2453 MeV）。
