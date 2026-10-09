# ParseLQCData 脚本目录总览 (Scripts Navigation Map)

本项目包含三大科学研究与分析子项目，按职责与物理对象清晰划分至不同子目录：

```text
scripts/
├── lcp/             # 【LCP 专区】常物理线 (Line of Constant Physics, L48T16 8温度点系列)
├── b417_scan/       # 【b4.17 扫描专区】固定 beta=4.17 的 24 个系综有限温度 (Nt) 与有限体积 (Ns) 介子扫描
├── ccor/            # 【ccor 破缺流专区】复刻 ana/dat/ccor 介子两信道质量差与手征对称性破缺流
└── run_build.sh     # C++ 底层动态库构建入口
```

---

## 一、 子项目划分与物理定位

### 1. `scripts/lcp/` —— 常物理线 (Line of Constant Physics) 分析专区
- **物理目标**：固定夸克质量在物理点，在 $48^3 \times 16$ 上扫描 8 个温度点 ($\beta = 4.13 \sim 4.405$)，观测手征相变与解禁闭。
- **涵盖系综**：`L48T16beta4.13...` 至 `L48T16beta4.405...`，以及强子关联函数 `48x16b4.13` 至 `48x16b4.30`。
- **核心脚本**：
  - `run_lcp_all.sh` / `run_lcp_all.py`：**一键运行全套 LCP 分析并出图**
  - `run_lcp_condensate.py`：手征凝聚抽取、扣残余质量、Zm 物理重整化
  - `run_lcp_susceptibility.py`：手征磁化率无偏二次量、Jackknife 与重整化，输出至 `output/LCP/`
  - `run_lcp_meson.py`：48x16 介子关联函数全量拟合与有效质量求解
  - `plot_lcp_condensate.py`、`plot_lcp_susceptibility.py`、`plot_lcp_meson.py`：全套物理图表绘制
- **详见**：[scripts/lcp/README.md](file:///Users/junxiongnie/code/build/ParseLQCData/scripts/lcp/README.md)

---

### 2. `scripts/b417_scan/` —— 固定 $\beta=4.17$ 有限温度与体积扫描
- **物理目标**：固定晶格间距 $a(\beta=4.17)$，通过改变 $N_t \in \{12, 14, 16, 18\}$ 改变温度，通过改变 $N_s \in \{32, 36, 40, 48\}$ 研究有限体积标度。
- **涵盖系综**：24 个系综 ($6\text{ 组格点} \times 4\text{ 组 } m_l \in \{0.0020, 0.0035, 0.0070, 0.0120\}$)。
- **核心脚本**：
  - `reproduce_meson_b417_nt_scan.py`：核心扫描计算脚本
  - `run_meson_b417_scan.py`：顶层调度包装入口
  - `run_meson_b417_nt_scan.sh`：一键批处理 Shell 脚本
  - `measure_meson_b417.py`：介子质量测量与温度演化图
  - `plot_effective_mass_b417.py`：有效质量与拟合窗口可视化
  - `try_fit_window.py`：拟合窗口平坦性试算与诊断
- **详见**：[scripts/b417_scan/README.md](file:///Users/junxiongnie/code/build/ParseLQCData/scripts/b417_scan/README.md)

---

### 3. `scripts/ccor/` —— 对称性破缺流与 Thesis 基准复刻
- **物理目标**：严格复刻 `ana/dat/ccor` 流程，通过两信道介子质量差观测 $SU(2)_L \times SU(2)_R$, $U(1)_A$, $SU(2)_{CS}$ 对称性恢复行为。
- **核心脚本**：
  - `run_ccor_flow.py`：端到端对称性破缺流调度
  - `reproduce_ccor_flow.py`：质量差与 Jackknife 统计计算
  - `update_ccor_data_and_plots.py`：格式转换与外部同步
  - `gnuplot/`：经典论文矢量图绘制模版
- **详见**：[scripts/ccor/README.md](file:///Users/junxiongnie/code/build/ParseLQCData/scripts/ccor/README.md)

---

## 二、 根目录兼容性转发入口 (Backwards-Compatibility Wrappers)

为了保证既有命令行习惯、自动化测试与 CI 流水线不发生断裂，根目录中的脚本一律作为轻量转发入口，自动指向对应的子项目脚本：

| 根目录旧入口 | 实际转发指向 | 所属项目 |
|---|---|---|
| `scripts/run_condensate.py` | `scripts/lcp/run_lcp_condensate.py` | LCP |
| `scripts/run_condensate.sh` | `uv run python scripts/lcp/run_lcp_condensate.py --plot` | LCP |
| `scripts/run_susceptibility.py` | `scripts/lcp/run_lcp_susceptibility.py` | LCP |
| `scripts/run_meson.py` | `scripts/lcp/run_lcp_meson.py` | LCP |
| `scripts/run_meson.sh` | `uv run python scripts/lcp/run_lcp_meson.py --source all --plot` | LCP |
| `scripts/plot_condensate.py` | `scripts/lcp/plot_lcp_condensate.py` | LCP |
| `scripts/plot_susceptibility.py`| `scripts/lcp/plot_lcp_susceptibility.py` | LCP |
| `scripts/plot_meson.py` | `scripts/lcp/plot_lcp_meson.py` | LCP |
| `scripts/run_meson_b417_nt_scan.sh` | `scripts/b417_scan/run_meson_b417_nt_scan.sh` | b417_scan |
| `scripts/run_meson_b417_scan.py` | `scripts/b417_scan/run_meson_b417_scan.py` | b417_scan |
| `scripts/run_ccor_flow.py` | `scripts/ccor/run_ccor_flow.py` | ccor |
| `scripts/reproduce_ccor_flow.py` | `scripts/ccor/reproduce_ccor_flow.py` | ccor |
