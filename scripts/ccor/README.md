# ccor 对称性破缺流项目 (ccor_flow)

> **定位说明**：本目录 (`scripts/ccor/`) 专门归集**严格复刻 `ana/dat/ccor`** 的手征与手征-自旋对称性破缺流分析。
> 该流程通过两信道介子有效质量差与逐样本 Jackknife 统计，观测 QCD 相变过程中各对称性的恢复行为。

---

## 1. 监测的 4 组物理信道对

1. **手征对称性 $SU(2)_L \times SU(2)_R$**:
   - 对比信道: 矢量与轴矢量信道质量差 ($V - A$)
   - 文件名/标头: `SU(2)XSU(2) V-A`
2. **轴矢反常 $U(1)_A$ (标量-赝标量)**:
   - 对比信道: 标量与赝标量质量差 ($S - PS$)
   - 文件名/标头: `U(1) S-PS`
3. **轴矢反常 $U(1)_A$ (张量信道)**:
   - 对比信道: 张量与赝张量质量差 ($T - X$)
   - 文件名/标头: `U(1)A T-X`
4. **手征-自旋对称性 $SU(2)_{CS}$**:
   - 对比信道: 轴矢量与赝张量质量差 ($X - A$)
   - 文件名/标头: `SU(2)spinxchiral X-A`

---

## 2. 脚本清单与职责分工

| 脚本文件 | 基础功能 | 对应产物路径 |
|---|---|---|
| `run_ccor_flow.py` | 统一包装入口，调用 `reproduce_ccor_flow.py` 执行端到端分析 | `output/ccor_flow/` |
| `reproduce_ccor_flow.py` | 抽取关联函数、拟合质量、计算两信道质量差、Jackknife 统计推断、生成 verification CSV | `output/ccor_flow/data*.txt`, `verification_summary.csv` |
| `update_ccor_data_and_plots.py` | 生成 ccor 格式并可配置同步至外部论文/展示目录 | `massvtem.csv`, `mdoutputre.csv`, thesis CSV |
| `plot_av_effective_mass_48x18.py` | 48x18 轴矢量 (AV) 信道现算有效质量并落盘出图 | `output/meson_scan/b4.17/ccor/48x18/` |
| `plot_av_effective_mass_36x18.py` | 36x18 轴矢量 (AV) 信道现算有效质量并落盘出图 | `output/meson_scan/b4.17/ccor/36x18/` |
| `gnuplot/` | 经典 gnuplot 绘图模版 (`plotmassdvsmass.gp`, `plotmd.gp`, `plotmdre.gp`) | 矢量图 PDF / PNG |

---

## 3. 运行示例

```bash
# 1. 运行 ccor 流程 (无需 gnuplot 依赖模式)
uv run python scripts/ccor/reproduce_ccor_flow.py --datasets multisrc --no-gnuplot

# 2. 运行完整 ccor 流程并生成全部 gnuplot 图表
uv run python scripts/ccor/run_ccor_flow.py
```
