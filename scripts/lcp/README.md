# LCP (Line of Constant Physics / 常物理线) 分析专区

> **定位说明**：本文档与本目录 (`scripts/lcp/`) 专门归集**格点 QCD 常物理线 (Line of Constant Physics)** 的处理、拟合与绘图脚本。
> 在格点 QCD 计算中，常物理线保持夸克质量调谐在物理点 ($m_\pi \approx 135\text{ MeV}$, $m_K \approx 495\text{ MeV}$)，通过改变规范耦合常数 $\beta$ 改变格距 $a(\beta)$，从而在固定的有限时间切片点数 $N_t=16$ 下实现**温度热演化扫描** ($T = 1/(N_t a)$)。

---

## 1. LCP 覆盖的系综与数据集

| 系综目录名 | 格点尺寸 | $\beta$ | $m_l$ | $m_s$ | $m_{\text{res}}$ | $Z_m(\beta)$ | $T$ (MeV) | 物理区段 |
|---|---|---|---|---|---|---|---|---|
| `L48T16beta4.13ms0.043547m0.000805` | $48^3 \times 16$ | 4.13 | 0.000805 | 0.043547 | 0.000731464 | 0.937703 | 138.23 | 手征破缺相 (低温区) |
| `L48T16beta4.15ms0.040843m0.000930` | $48^3 \times 16$ | 4.15 | 0.000930 | 0.040843 | 0.000498493 | 0.952390 | 145.75 | 手征破缺相 |
| `L48T16beta4.17ms0.038400m0.001001` | $48^3 \times 16$ | 4.17 | 0.001001 | 0.038400 | 0.000339722 | 0.966247 | 153.31 | 相变过渡带前沿 |
| `L48T16beta4.18ms0.037265m0.001022` | $48^3 \times 16$ | 4.18 | 0.001022 | 0.037265 | 0.000280451 | 0.972899 | 157.03 | **赝临界相变区 ($T_{pc}$)** |
| `L48T16beta4.20ms0.035150m0.001041` | $48^3 \times 16$ | 4.20 | 0.001041 | 0.035150 | 0.000191127 | 0.985710 | 164.55 | 手征恢复相 (高温区) |
| `L48T16beta4.23ms0.032315m0.001033` | $48^3 \times 16$ | 4.23 | 0.001033 | 0.032315 | 0.000107528 | 1.003850 | 176.43 | 手征恢复相 |
| `L48T16beta4.30ms0.026930m0.000939` | $48^3 \times 16$ | 4.30 | 0.000939 | 0.026930 | 0.000028096 | 1.042240 | 202.52 | QGP 夸克胶子等离子体 |
| `L48T16beta4.405ms0.021032m0.000760`| $48^3 \times 16$ | 4.405| 0.000760 | 0.021032 | 0.000003753 | 1.092740 | 241.60 | 极高温区 |

> 强子/介子关联函数在 `data/readin/48x16b4.{beta}` 目录下。

---

## 2. 脚本清单与职责分工

| 脚本文件 | 基础功能 | 对应产物路径 |
|---|---|---|
| `run_lcp_all.py` / `run_lcp_all.sh` | **LCP 全流程总入口**：依次运行凝聚、磁化率与介子分析并生成全套图表 | `output/LCP/`, `output/condensate/`, `docs/figures/` |
| `run_lcp_condensate.py` | 抽取 8 个 LCP 系综的裸光/奇手征凝聚，扣除残余质量散度并乘 $Z_m$ 重整化 | `output/condensate/results_rm_beta.txt`, `output/condensate/all_ensembles_condensate.csv` |
| `run_lcp_susceptibility.py` | 抽取轻夸克微观随机源向量，计算单构型无偏两体乘积与 Jackknife，施加 $Z_m^2$ 重整化 | `output/LCP/results_susceptibility.txt`, `output/LCP/results_susceptibility.csv` |
| `run_lcp_meson.py` | 抽取 7 组 $\beta$ 的 6 种狄拉克介子信道空间关联函数，求解有效质量并做 cosh 平台拟合 | `output/pickdata*`, `output/ratio_results*`, `output/simulateresult*`, `output/all_fits_summary.csv` |
| `plot_lcp_condensate.py` | 绘制重整化手征凝聚、裸凝聚、物理三次根、有限体积标度图 | `docs/figures/conden_re_plot.png`, `conden_physical_cuberoot.png` 等 |
| `plot_lcp_susceptibility.py` | 绘制手征磁化率随温度演化曲线、过渡带渲染与 $T_{pc} \approx 157.0\text{ MeV}$ 极大值峰位 | `docs/figures/pbpchisce_renormalized.png`, `pbp_temperature_plot.png` 等 |
| `plot_lcp_meson.py` | 绘制 42 组信道对比、介子质量热演化、矢量-轴矢量手征恢复与对称性破缺信号 | `docs/figures/channel_comparisons/`, `simulation_results.png` 等 |

---

## 3. 快速运行指南

```bash
# 1. 一键运行 LCP 全流程 (手征凝聚 + 手征磁化率 + 介子分析 + 出图)
bash scripts/lcp/run_lcp_all.sh

# 2. 单独运行 LCP 手征凝聚并出图
uv run python scripts/lcp/run_lcp_condensate.py --plot

# 3. 单独运行 LCP 手征磁化率并出图
uv run python scripts/lcp/run_lcp_susceptibility.py --plot

# 4. 单独运行 LCP 介子温度扫描 (多进程并发)
uv run python scripts/lcp/run_lcp_meson.py --workers 4 --plot

# 5. 仅出图 (不重复计算数据)
uv run python scripts/lcp/plot_lcp_condensate.py
uv run python scripts/lcp/plot_lcp_susceptibility.py
uv run python scripts/lcp/plot_lcp_meson.py
```
