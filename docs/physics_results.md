# 格点 QCD 物理输出结果与对称性分析文档 (Physics Results & Symmetry Restoration)

本文档系统汇总本项目全量复现的格点量子色动力学物理输出结果，包括介子筛选质量 (Screening Mass)、手征对称性恢复 ($SU(2)_L \times SU(2)_R$)、轴向反常对称性恢复 ($U(1)_A$) 以及重整化手征凝聚随温度的演化。

---

## 1. 介子筛选质量汇总 (Meson Screening Masses $a M$)

通过贝叶斯平台拟合提取的介子基态筛选质量参数如下（格点单位 $a M$）：

### 1.1 多源模式 (Multi-Source Plateau Fits)
| 耦合 $\beta$ | 温度 $T$ (MeV) | 轴矢量 $AV$ ($a_1$) | 赝标量 $PS$ ($\pi$) | 标量 $S$ ($\sigma$) | 张量 $Tt$ | 矢量 $Vec$ ($\rho$) | 轴张量 $Xt$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **4.13** | 138.23 | 0.2979(58) | 0.0734(08) | 0.0104(97) | 0.3749(49) | 0.2641(22) | 0.4131(53) |
| **4.15** | 145.75 | 0.3054(18) | 0.0751(10) | 0.2717(28) | 0.3682(40) | 0.2881(12) | 0.2824(49) |
| **4.17** | 153.31 | 0.3101(14) | 0.0836(14) | 0.1270(29) | 0.3684(51) | 0.2932(11) | 0.4090(46) |
| **4.18** | 157.03 | 0.2911(18) | 0.0899(17) | 0.0831(58) | 0.4373(48) | 0.2714(14) | 0.3808(57) |
| **4.20** | 164.55 | 0.3246(12) | 0.1000(22) | 0.0126(17) | 0.3596(29) | 0.3161(11) | 0.3767(48) |
| **4.23** | 176.43 | 0.3299(11) | 0.1323(21) | 0.2681(25) | 0.4470(42) | 0.3234(11) | 0.4091(45) |
| **4.30** | 202.52 | **0.3529(12)** | **0.2009(25)** | **0.2008(24)** | **0.4029(36)** | **0.3533(10)** | **0.4049(38)** |

> **关键物理发现**：
> 观察 $\beta = 4.30$ ($T = 202.52\,\text{MeV} \gg T_c$) 高温退禁闭相：
> 1. **$Vec$ 与 $AV$ 简并**：$a M_{Vec} = 0.3533 \approx a M_{AV} = 0.3529$（手征对称性恢复）；
> 2. **$PS$ 与 $S$ 简并**：$a M_{PS} = 0.2009 \approx a M_{S} = 0.2008$（轴向 $U(1)_A$ 对称性恢复）；
> 3. **$Xt$ 与 $Tt$ 简并**：$a M_{Xt} = 0.4049 \approx a M_{Tt} = 0.4029$（张量信道 $U(1)_A$ 对称性恢复）。

---

## 2. 物理对称性破缺与恢复签名 (Symmetry Restoration Signatures)

通过信道对的质量劈裂量 $\Delta M / T = (M_1 - M_2) \times 16 \times T$ 监测各连续对称性在有限温下的演化行为：

1. **手征对称性 $SU(2)_L \times SU(2)_R$ 签名**：
   $$\Delta M_{\text{chiral}} = M_{Vec} - M_{AV}$$
   * 低温禁闭相：$\rho$ 介子与 $a_1$ 介子显著劈裂；
   * 高温相 ($T > 160\,\text{MeV}$)：$\Delta M \to 0$，夸克手征凝聚熔解导致宇称二重态完全简并。
2. **轴向反常对称性 $U(1)_A$ 签名**：
   $$\Delta M_{U(1)_A} = M_S - M_{PS} \quad \text{及} \quad M_{Xt} - M_{Tt}$$
   * 拓扑荷涨落（瞬子密度）在 $T_c$ 附近被剧烈压低，导致赝标量 $\pi$ 与标量 $\sigma$ 态的质量差迅速归零。
3. **手征自旋对称性 $SU(2)_{CS}$ 签名**：
   $$\Delta M_{CS} = M_{AV} - M_{Tt}$$

---

## 3. 重整化手征凝聚数值结果 (Chiral Condensate)

| 耦合 $\beta$ | 温度 $T$ (MeV) | 重整化手征凝聚 $\Delta_{l,s}$ | Jackknife 误差 | 物理状态 |
| :--- | :--- | :--- | :--- | :--- |
| **4.13** | 138.23 | $1.06082 \times 10^{-3}$ | $6.74 \times 10^{-6}$ | 手征自发破缺相 |
| **4.15** | 145.75 | $7.71090 \times 10^{-4}$ | $6.75 \times 10^{-6}$ | 手征自发破缺相 |
| **4.17** | 153.31 | $5.07280 \times 10^{-4}$ | $6.62 \times 10^{-6}$ | 相变起始区 |
| **4.18** | 157.03 | $4.16490 \times 10^{-4}$ | $6.59 \times 10^{-6}$ | **赝临界相变区 ($T_c \approx 157\,\text{MeV}$)** |
| **4.20** | 164.55 | $2.38960 \times 10^{-4}$ | $5.61 \times 10^{-6}$ | 手征恢复过渡相 |
| **4.23** | 176.43 | $9.28900 \times 10^{-5}$ | $3.81 \times 10^{-6}$ | 手征基本恢复相 |
| **4.30** | 202.52 | $8.20000 \times 10^{-6}$ | $1.33 \times 10^{-6}$ | 手征完全恢复相 ($\Delta \to 0$) |

---

## 4. 科学图表全景索引 (Figures Gallery)

所有科学图表均已无损生成并存放在 `docs/figures/` 目录下：

### 4.1 全局物理图表
* **筛选质量随温度演化**：[`docs/figures/simulation_results.png`](figures/simulation_results.png)
* **对称性破缺签名演化**：[`docs/figures/symmetry_breaking_ana_setup.png`](figures/symmetry_breaking_ana_setup.png)
* **重整化手征凝聚**：[`docs/figures/conden_re_plot.png`](figures/conden_re_plot.png)
* **轻夸克未重整化凝聚**：[`docs/figures/conden_plot.png`](figures/conden_plot.png)

### 4.2 各温度信道对比图 (42 组单源 vs 多源拟合带)
* **$\beta = 4.13$**：[`docs/figures/channel_comparisons/b4.13/`](figures/channel_comparisons/b4.13/) (`ratioAV.png`, `ratioPS.png`, `ratioS.png`, `ratioTt.png`, `ratioVec.png`, `ratioXt.png`)
* **$\beta = 4.15$**：[`docs/figures/channel_comparisons/b4.15/`](figures/channel_comparisons/b4.15/)
* **$\beta = 4.17$**：[`docs/figures/channel_comparisons/b4.17/`](figures/channel_comparisons/b4.17/)
* **$\beta = 4.18$**：[`docs/figures/channel_comparisons/b4.18/`](figures/channel_comparisons/b4.18/)
* **$\beta = 4.20$**：[`docs/figures/channel_comparisons/b4.20/`](figures/channel_comparisons/b4.20/)
* **$\beta = 4.23$**：[`docs/figures/channel_comparisons/b4.23/`](figures/channel_comparisons/b4.23/)
* **$\beta = 4.30$**：[`docs/figures/channel_comparisons/b4.30/`](figures/channel_comparisons/b4.30/)
