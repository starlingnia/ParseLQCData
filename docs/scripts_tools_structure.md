# ParseLQCData 脚本与工具调用结构梳理

> 本文档梳理当前仓库中 `scripts/`、`tools/`、`src/parselqcdata/` 以及原生 CLI `bin/ParseLQCData` 的职责、调用关系和使用方式。
> 目标：说明“哪个脚本调用了哪个工具、基础功能是什么、应该怎么运行”。

---

## 1. 总体调用链

```text
scripts/*.sh
  └─ scripts/run_*.py、reproduce_*.py
       ├─ docs/physics_setup.py / docs/meson_scan_setup.py    # 全局物理参数、任务表
       ├─ src/parselqcdata/*.py                               # Python 管道、拟合、统计
       ├─ tools/meson_orchestrator.py                         # Meson ctypes 桥
       │    └─ lib/liblqcd_meson.* / libparselqcdata.*
       ├─ tools/condensate_orchestrator.py                    # Condensate / Susceptibility ctypes 桥
       │    └─ lib/liblqcd_condensate.* / libparselqcdata.*
       └─ scripts/plot_*.py                                   # 读取 output 落盘数据并出图

原生 C++ 路径：
bin/ParseLQCData
  └─ apps/main.cpp
       └─ build/*_task.cpp  (REGISTER_TASK 注册)
            └─ src/core、src/IOdata、src/Statistics、
               src/MesonAnalysis、src/CondensateAnalysis
```

一句话总结：

- `scripts/` 负责批处理调度、参数组织、拟合汇总和绘图；
- `tools/` 负责 Python 与 C++ 动态库之间的 ctypes / C ABI 桥接；
- `src/parselqcdata/` 负责 Python 侧可复用的物理分析管道；
- C++ 动态库负责高性能文本/XML I/O、多线程抽取、binning 和 Jackknife 等底层计算；
- `bin/ParseLQCData` 是另一条原生 CLI 路径，直接通过注册表任务调用 C++ 核心。

---

## 2. `tools/`：Python 调 C++ 的桥接层

| 工具 | 基础功能 | 主要接口 |
|---|---|---|
| `tools/meson_orchestrator.py` | 加载 `libparselqcdata.dylib`、`liblqcd_meson.dylib` 等动态库；调用介子关联函数 C ABI | `process_channel()`、`process_channel_multi()`、`process_channel_single()`、`execute_meson_analysis()` |
| `tools/condensate_orchestrator.py` | 加载手征凝聚/磁化率动态库；调用相关 C ABI | `process_condensate()`、`process_condensate_full()`、`process_susceptibility()` |
| `tools/lqcd_orchestrator.py` | 统一包装 Meson 与 Condensate 两个桥接器 | `process_channel()`、`process_chiral_condensate()`、`execute_meson_analysis()` |
| `tools/Services/MesonService.cpp` | 导出介子分析纯 C ABI | `run_meson_pipeline_c_api`、`run_meson_multisrc_pipeline_c_api`、`run_meson_singlesrc_pipeline_c_api`、`free_lqcd_buffer` |
| `tools/Services/CondensateService.cpp` | 导出手征凝聚/磁化率纯 C ABI | `run_chiral_condensate_c_api`、`run_chiral_condensate_full_c_api`、`run_chiral_susceptibility_c_api` |

关键点：

1. Python 侧通过 `ctypes.CDLL` 查找 `lib/` 下的动态库。
2. C ABI 返回值统一使用 `int` 状态码；结果数组/矩阵通过指针回传。
3. `free_lqcd_buffer` 负责释放 C++ 分配给 Python 的缓冲区。
4. 业务脚本通常直接使用 `MesonOrchestrator` 或 `CondensateOrchestrator`；`LQCDOrchestrator` 是统一兼容封装。

---

## 3. `src/parselqcdata/`：Python 独立原子模块与基础管道

项目推行**“单脚本单职责 (Single Responsibility Principle)”**规范，将原本集中在管道内的各类函数剥离为独立可执行的原子脚本，并由顶层 Pipeline 进行组装：

### 3.1 独立原子计算脚本与模块 (可直接 CLI 运行)

| 脚本文件 | 基础功能 | CLI 使用示例 |
|---|---|---|
| `xml_pbp_extractor.py` | 正则抽取 XML 中 `<pbp>` 随机源实部测量值，支持单文件与目录批量扫描 | `uv run python src/parselqcdata/xml_pbp_extractor.py <meas_dir> --quark light` |
| `unbiased_quadratic.py` | 单构型无偏二次交叉估计器 $O_{2\text{bar}} = \frac{1}{k(k-1)}\sum_{i\neq j} O_i O_j$，消除有限随机源方差偏差 | `uv run python src/parselqcdata/unbiased_quadratic.py 1.0 2.0 3.0` |
| `scaling_factors.py` | 计算四维体积因子 $F_{\text{vol}} = N_s^3 N_t$、标度因子 $F_{\text{scaled}} = (N_t T)^2 F_{\text{vol}}$、查询 $Z_m(\beta)$ | `uv run python src/parselqcdata/scaling_factors.py --ns 48 --nt 16 --temp 153.31 --beta 4.17` |
| `susceptibility_calculator.py` | 构型级 Jackknife 统计重采样求解 $\chi = \langle O_2 \rangle - \langle O \rangle^2$ 并折算 $Z_m^2$ 重整化 | `uv run python src/parselqcdata/susceptibility_calculator.py --ns 48 --nt 16 --temp 153.31 --beta 4.17` |
| `lcp_exporter.py` | 规范化导出 LCP 标准产物 (`results_susceptibility.txt`, `.csv`, `.parquet`) 至 `output/LCP/` | `uv run python src/parselqcdata/lcp_exporter.py --input <results.csv>` |
| `condensate_subtraction.py` | 手征凝聚残余质量相减 $\langle\bar{\psi}\psi\rangle_{\text{sub}} = Z_m (\langle\bar{\psi}\psi\rangle_l - \frac{m_l+m_{\text{res}}}{m_s+m_{\text{res}}}\langle\bar{\psi}\psi\rangle_s)$ | `uv run python src/parselqcdata/condensate_subtraction.py --pbpl ... --pbps ... --ml ... --ms ... --mres ...` |
| `effective_mass_solver.py` | 求解双曲余弦比值有效质量方程，支持参数化对称点 $N_s/2$，向量化 Newton-Raphson 极速求解 | `uv run python src/parselqcdata/effective_mass_solver.py --ratio 1.25 --x 10 --half 24` |
| `plateau_detector.py` | 基于 Jackknife 差分的平坦性 $\chi^2/\text{dof}$ 自动平台窗口检测与选优 | `uv run python src/parselqcdata/plateau_detector.py --csv <meff_matrix.csv> --end 24` |
| `cosh_fitter.py` | 统一 chi2 最小二乘非线性平台拟合器，逐 Jackknife 样本拟合与统计推断 | `uv run python src/parselqcdata/cosh_fitter.py --start 14 --end 24 --half 24` |

### 3.2 顶层业务管道

| 模块 | 基础功能 |
|---|---|
| `meson_pipeline.py` | `MesonPipeline.process_channel()` / `run_ensemble()`，按 `beta + channel` 调度 MesonOrchestrator |
| `condensate_pipeline.py` | `CondensatePipeline.process_ensemble()` / `process_all_ensembles()`，独立落盘、增量总表、导出 scaling 子表 |
| `susceptibility_pipeline.py` | `SusceptibilityPipeline`：组装 XML 抽取、无偏估计、Jackknife、LCP 输出与集群扫描 |
| `ccor_flow.py` | 严格复刻 `ana/dat/ccor` 的“两信道质量差 → 对称性破缺”流程 |

---

## 4. `scripts/`：子项目划分与脚本导航

为了彻底解决“LCP 处理脚本与其他项目混杂”的问题，`scripts/` 已按**物理对象与子项目独立分目录管理**：

```text
scripts/
├── lcp/             # 【LCP 专区】常物理线 (Line of Constant Physics, L48T16 8温度点系列)
├── b417_scan/       # 【b4.17 扫描专区】固定 beta=4.17 的 24 个系综有限温度 (Nt) 与体积 (Ns) 扫描
├── ccor/            # 【ccor 破缺流专区】复刻 ana/dat/ccor 介子两信道质量差与对称性破缺流
└── run_build.sh     # C++ 底层动态库构建入口
```

### 4.1 【LCP 专区】常物理线核心分析与绘图 (`scripts/lcp/`)

> 详见 [scripts/lcp/README.md](file:///Users/junxiongnie/code/build/ParseLQCData/scripts/lcp/README.md)

| 脚本 | 功能与说明 | 核心输出 |
|---|---|---|
| `run_lcp_all.sh` / `run_lcp_all.py` | **一键运行 LCP 全流程**（凝聚 + 磁化率 + 介子 + 出图） | `output/LCP/`, `output/condensate/`, `docs/figures/` |
| `run_lcp_condensate.py` | 扫描 8 组 LCP 系综，扣除残余质量散度，施加 $Z_m$ 物理重整化 | `output/condensate/results_rm_beta.txt` |
| `run_lcp_susceptibility.py` | 抽取轻夸克随机源向量，计算无偏二次量与 Jackknife，折算 $Z_m^2$ 重整化 | `output/LCP/results_susceptibility.txt`, `output/LCP/results_susceptibility.csv` |
| `run_lcp_meson.py` | 48³×16 上 7 组 $\beta$ 的 6 种信道介子有效质量求解与 cosh 平台拟合 | `output/pickdata*`, `output/ratio_results*`, `output/simulateresult*` |
| `plot_lcp_condensate.py` | 绘制重整化手征凝聚、裸凝聚、物理三次根、有限体积标度图 | `docs/figures/conden_re_plot.png` 等 |
| `plot_lcp_susceptibility.py` | 绘制手征磁化率随温度演化曲线与 $T_{pc} \approx 157.0\text{ MeV}$ 极大值峰位 | `docs/figures/pbpchisce_renormalized.png` 等 |
| `plot_lcp_meson.py` | 绘制 42 组信道对比、介子质量热演化、矢量-轴矢量手征恢复与对称性破缺 | `docs/figures/channel_comparisons/` 等 |

### 4.2 【b4.17 扫描专区】固定 $\beta=4.17$ 有限温度与体积扫描 (`scripts/b417_scan/`)

> 详见 [scripts/b417_scan/README.md](file:///Users/junxiongnie/code/build/ParseLQCData/scripts/b417_scan/README.md)

| 脚本 | 功能与说明 | 核心输出 |
|---|---|---|
| `run_meson_b417_nt_scan.sh` | 一键多进程批处理 Shell 脚本 | `output/meson_scan/b4.17/` |
| `reproduce_meson_b417_nt_scan.py` | 核心扫描脚本：24 系综 $\times$ 6 信道空间关联函数抽取、对称点 $N_s/2$ 有效质量求解、平台拟合 | `output/meson_scan/b4.17/multisrc/`, `singlesrc/` |
| `run_meson_b417_scan.py` | 顶层调度入口，支持 `--source`、`--workers`、`--force`、`--plot` | 同上 |
| `measure_meson_b417.py` | 介子质量测量与温度演化图 | `docs/meson_mass_T_*.png` |
| `plot_effective_mass_b417.py` | 批量绘制有效质量曲线与拟合窗口可视化诊断 | `docs/effective_mass/` |
| `try_fit_window.py` | 针对特定系综的拟合窗口试算与平坦性诊断 | 控制台诊断打印 |
| `plot_split_raw.py` | 拆分“48³×16 温度扫描”和“固定 beta=4.17 序列”并绘制原始数据图 | `output/condensate/split/` |

### 4.3 【ccor 破缺流专区】对称性破缺恢复与论文复刻 (`scripts/ccor/`)

> 详见 [scripts/ccor/README.md](file:///Users/junxiongnie/code/build/ParseLQCData/scripts/ccor/README.md)

| 脚本 | 功能与说明 | 核心输出 |
|---|---|---|
| `run_ccor_flow.py` | 端到端对称性破缺流调度入口 | `output/ccor_flow/` |
| `reproduce_ccor_flow.py` | 对照 `ana/dat/ccor`，输出数据表、PDF、PNG、verification CSV | `output/ccor_flow/data*.txt`, `verification_summary.csv` |
| `update_ccor_data_and_plots.py` | 生成 ccor 格式并同步至外部论文/展示目录 | `massvtem.csv`, `mdoutputre.csv`, thesis CSV |
| `plot_av_effective_mass_48x18.py` | 48x18 轴矢量 (AV) 信道现算有效质量并落盘 | `output/meson_scan/b4.17/ccor/48x18/` |
| `plot_av_effective_mass_36x18.py` | 36x18 轴矢量 (AV) 信道现算有效质量并落盘 | `output/meson_scan/b4.17/ccor/36x18/` |
| `gnuplot/` | 经典 gnuplot 绘图模版与脚本 | 矢量图 PDF / PNG |

### 4.4 根目录兼容性转发入口 (Backwards Compatibility)

根目录既有脚本全部保留并作为轻量转发入口，确保 CI 与已有调用习惯完全正常运作：
- `scripts/run_condensate.py` / `.sh` $\to$ `scripts/lcp/run_lcp_condensate.py`
- `scripts/run_susceptibility.py` $\to$ `scripts/lcp/run_lcp_susceptibility.py`
- `scripts/run_meson.py` / `.sh` $\to$ `scripts/lcp/run_lcp_meson.py`
- `scripts/plot_condensate.py` $\to$ `scripts/lcp/plot_lcp_condensate.py`
- `scripts/plot_susceptibility.py` $\to$ `scripts/lcp/plot_lcp_susceptibility.py`
- `scripts/plot_meson.py` $\to$ `scripts/lcp/plot_lcp_meson.py`
- `scripts/run_meson_b417_nt_scan.sh` $\to$ `scripts/b417_scan/run_meson_b417_nt_scan.sh`
- `scripts/run_meson_b417_scan.py` $\to$ `scripts/b417_scan/run_meson_b417_scan.py`
- `scripts/run_ccor_flow.py` $\to$ `scripts/ccor/run_ccor_flow.py`
- `scripts/reproduce_ccor_flow.py` $\to$ `scripts/ccor/reproduce_ccor_flow.py`

这些脚本通常只是导入另一个脚本的 `main` 或 `run_*_pipeline`：

- `reproduce_condensate.py` → `run_condensate.main()`
- `reproduce_susceptibility.py` → `run_susceptibility.main()`
- `reproduce_meson_multi.py` → `run_meson.run_meson_pipeline(sources=[False])`
- `reproduce_meson_single.py` → `run_meson.run_meson_pipeline(sources=[True])`
- `reproduce_meson.py` → 根据第一个位置参数 `single/multi/all` 调 `run_meson_pipeline`
- `run_ccor_flow.py` → `reproduce_ccor_flow.main()`
- `quick_plot_all_condensates.py`、`quick_plot_condensate_physical.py` → `plot_condensate.main()`

### 4.4 绘图与诊断脚本

| 脚本 | 输入 | 输出 | 说明 |
|---|---|---|---|
| `plot_condensate.py` | `output/condensate/*` | `docs/figures/` 等 | 手征凝聚重整化、裸凝聚、标度、总览图 |
| `plot_meson.py` | `output/simulateresult*`、`ratio_results*` | `docs/figures/` | 42 组信道对比、有效质量热演化、对称性破缺 |
| `plot_susceptibility.py` | `results_susceptibility.csv`、`all_ensembles_susceptibility.csv` | `docs/figures/`、`output/plots/` 等 | 手征磁化率多种标度图 |
| `plot_effective_mass_b417.py` | `output/meson_scan/...`、`config/fit_windows.txt` | `docs/effective_mass/` | 批量画有效质量、拟合窗口诊断、README 报告 |
| `plot_split_raw.py` | `all_ensembles_condensate.csv`、`all_ensembles_susceptibility.csv` | `output/condensate/split/` | 拆分“48^3×16 温度扫描”和“固定 beta=4.17 序列”并画原始数据图 |
| `plot_av_effective_mass_48x18.py` | 真实 AV 关联函数 | `output/meson_scan/b4.17/ccor/48x18/` | AV 信道示例，现算有效质量并落盘 |
| `plot_av_effective_mass_36x18.py` | 同上 | 36x18 输出 | 通过改模块常量复用 48x18 脚本 |
| `try_fit_window.py` | `output/meson_scan` 已落盘数据 | 仅打印 | 只读试算拟合窗口，支持 `--mode scan` / `--mode ana` |
| `measure_meson_b417.py` | 真实目录或 meson_scan 缓存 | `docs/meson_mass_T_*.png` | 测量 + 按 ml 出温度演化图 |

---

## 5. C++ 原生 CLI：`bin/ParseLQCData`

先看任务列表：

```bash
./bin/ParseLQCData
```

已注册任务：

| 任务 | 参数 | 说明 |
|---|---|---|
| `condensate` | `[dataset_path_or_beta] [ml] [ms] [mres] [zm]` | 单数据集手征凝聚 |
| `condensate_all` | `[readin_dir]` | 批量扫描手征凝聚数据集 |
| `meson_multi` | `[beta] [channel] [binsize]` | 多源介子关联函数 |
| `meson_single` | `[beta] [channel] [binsize]` | 单源介子关联函数 |
| `meson_all` | `[beta] [multi|single]` | 批量所有信道 |
| `meson_multi_all` | `[beta]` | 批量多源全信道 |
| `meson_single_all` | `[beta]` | 批量单源全信道 |
| `susceptibility` | `[readin_dir] [output_dir]` | 手征磁化率 C++ 抽取 |

示例：

```bash
./bin/ParseLQCData condensate L32T12_beta4.17ms0.040m0.0020
./bin/ParseLQCData condensate_all data/readin
./bin/ParseLQCData meson_multi 17 AV 4
./bin/ParseLQCData meson_single 17 S 4
./bin/ParseLQCData meson_all 18 multi
./bin/ParseLQCData susceptibility data/readin output/condensate
```

说明：`bin/pars` 也存在，但只列出 5 个任务，缺少 `susceptibility` 和部分批量任务，可视为旧/精简 CLI；正式使用应优先 `bin/ParseLQCData`。

---

## 6. 常用使用方式

```bash
# 1. 构建 C++ 动态库
bash scripts/run_build.sh

# 2. 手征凝聚全量分析 + 画图
bash scripts/run_condensate.sh

# 只跑某个系综
.venv/bin/python scripts/run_condensate.py \
    --ensemble L32T12_beta4.17ms0.040m0.0020 --plot

# 3. 介子旧 beta 扫描（48^3×16，beta=4.13~4.30）
bash scripts/run_meson.sh
.venv/bin/python scripts/run_meson.py --source multi --beta 17 --plot

# 4. b4.17 Nt 扫描
bash scripts/run_meson_b417_nt_scan.sh --dry-run
bash scripts/run_meson_b417_nt_scan.sh --source multi --cases 48x18 --ml 0.0020 0.0120
.venv/bin/python scripts/run_meson_b417_scan.py --source both --plot --workers 4

# 5. 手征磁化率
.venv/bin/python scripts/run_susceptibility.py --readin-dir data/readin --plot

# 6. 仅出图
.venv/bin/python scripts/plot_condensate.py
.venv/bin/python scripts/plot_meson.py
.venv/bin/python scripts/plot_susceptibility.py

# 7. ccor 对照
.venv/bin/python scripts/reproduce_ccor_flow.py --datasets multisrc --no-gnuplot

# 8. 拟合窗口试算
.venv/bin/python scripts/try_fit_window.py --case 48x18 --ml 0.0020 --channels PS
```

---

## 7. 注意事项

1. **多数脚本要从项目根目录运行**，因为代码大量使用 `PROJECT_ROOT = Path(__file__).resolve().parent.parent` 和相对路径。
2. Python 环境推荐 `.venv/bin/python` 或 `uv run`；Shell 脚本两种都有使用。
3. 很多脚本有硬编码 fallback 路径，例如 `~/code/ana/dat/readin`、`~/code/ths/pos/ccor`；本地 `data/readin` 不存在时会回退这些路径。
4. `update_ccor_data_and_plots.py` 默认会覆盖外部 ccor 目录中的同名文件，**先改 `--ccor-dir` 到安全目录**。
5. 物理标度常量不完全统一：部分脚本用 2453 MeV，部分用 `TEMP_MAP` 推导的 2452.96 MeV，`plot_effective_mass_b417.py` 里又有 `153×16 = 2448`。比较结果时要注意口径。
6. `run_build.sh` 的手工 fallback 输出固定 `.dylib` 名称；在 Linux 上如需 `.so`，要自行调整。
7. `scripts/` 本身主要是调度和绘图，真正高性能的 I/O、XML/文本抽取、binning、Jackknife 在 C++ 动态库中；Python 层负责参数组织、拟合、汇总和出图。
