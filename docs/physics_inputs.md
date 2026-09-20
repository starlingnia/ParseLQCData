# 格点 QCD 物理输入与构型元数据文档 (Physics Inputs & Configuration Metadata)

本文档系统整理本项目分析的所有格点量子色动力学 (LQCD) 规范场构型参数、强子算子定义、热力学温度刻度以及统计重采样设置。

---

## 1. 规范场构型与格点几何结构 (Ensemble Geometry)

* **费米子作用量**：手征对称性保持优异的畴壁费米子 (Domain Wall Fermions, DWF)
* **格点空间几何**：$N_s^3 \times N_t = 48^3 \times 16$（具有高空间体积纵横比 $N_s / N_t = 3$，有效抑制有限体积效应）
* **时间反演对称性周期**：$N_t = 48$（时间片 $t \in [0, 47]$，反转对称中心位于 $t = 24$）

### 1.1 耦合常数 $\beta$ 与温度刻度 $T$ (MeV)

物理温度通过反向格距 $a^{-1}(\beta)$ 与时间格点数 $N_t$ 由关系式 $T = 1 / (N_t a)$ 确定：

| 耦合常数 $\beta$ | 物理温度 $T$ (MeV) | 相对转变区域位置 | 原始数据路径 | 原始构型数 |
| :--- | :--- | :--- | :--- | :--- |
| **4.13** | 138.23 | 低温禁闭相 ($T < T_c$) | `data/readin/48x16b4.13/Output/` | 392 cfgs |
| **4.15** | 145.75 | 低温禁闭相 ($T < T_c$) | `data/readin/48x16b4.15/Output/` | 396 cfgs |
| **4.17** | 153.31 | 相变前夕区 ($T \lesssim T_c$) | `data/readin/48x16b4.17/Output/` | 396 cfgs |
| **4.18** | 157.03 | **赝临界相变区 ($T \approx T_c$)** | `data/readin/48x16b4.18/Output/` | 400 cfgs |
| **4.20** | 164.55 | 早期退禁闭相 ($T \gtrsim T_c$) | `data/readin/48x16b4.20/Output/` | 404 cfgs |
| **4.23** | 176.43 | 高温退禁闭相 ($T > T_c$) | `data/readin/48x16b4.23/Output/` | 392 cfgs |
| **4.30** | 202.52 | 高温退禁闭相 ($T \gg T_c$) | `data/readin/48x16b4.30/Output/` | 404 cfgs |

> **相变过渡带 (Transition Region)** 标定为：$T \in [155.5, 160.5]\,\text{MeV}$（以 $\beta \approx 4.18$ 为对称中心）。

---

## 2. 介子流算子与空间投影信道映射 (Meson Channels)

系统解析以下 6 类标准狄拉克介子流（双夸克束缚态或等效能量标度算子）：

| 信道代号 | 物理态对应 | 算子矩阵结构 | 空间方向组合 (Spatial Directions & Types) |
| :--- | :--- | :--- | :--- |
| **AV** | 轴矢量介子 ($a_1$) | $\bar{\psi} \gamma_5 \gamma_\mu \psi$ | AVector2 (DIRX), AVector3 (DIRX), AVector3 (DIRY), AVector1 (DIRY), AVector1 (DIRZ), AVector2 (DIRZ) |
| **S** | 标量介子 ($\sigma / f_0$) | $\bar{\psi} \psi$ | TVector4 (DIRZ), TVector4 (DIRX), TVector4 (DIRY) |
| **Tt** | 张量介子 ($t$-投影) | $\bar{\psi} \sigma_{0\mu} \psi$ | TVector1 (DIRX), TVector3 (DIRZ), TVector2 (DIRY) |
| **PS** | 赝标量介子 ($\pi$) | $\bar{\psi} \gamma_5 \psi$ | TAVector4 (DIRZ), TAVector4 (DIRX), TAVector4 (DIRY) |
| **Xt** | 轴张量介子 ($t$-投影) | $\bar{\psi} \gamma_5 \sigma_{0\mu} \psi$ | TAVector1 (DIRX), TAVector3 (DIRZ), TAVector2 (DIRY) |
| **Vec** | 矢量介子 ($\rho$) | $\bar{\psi} \gamma_\mu \psi$ | Vector2 (DIRX), Vector3 (DIRX), Vector3 (DIRY), Vector1 (DIRY), Vector1 (DIRZ), Vector2 (DIRZ) |

### 2.1 源抽取模式 (Source Modes)
1. **多源叠加平均模式 (Multi-Source)**：
   * 文本特征匹配：`--- {type} to {type} --- spatial:{dir} ---`
   * 取前 16 个源块，根据各源空间坐标偏移量 $\vec{x}_{\text{src}}$ 进行循环移位对齐（Roll），求 16 源代数平均以压低规范噪声；
2. **单源点源模式 (Single-Source)**：
   * 文本特征匹配：`--- {type} to {type} --- spatial:{dir} --- 0/0/0/0`
   * 直接提取位于格点原点的点源强子关联函数。

---

## 3. 夸克质量与手征重整化输入参数 (Chiral Parameters)

在畴壁费米子形式下，由于第五维有限长度 $L_s$，存在非零残余质量 $m_{\text{res}}$ 破坏精确手征对称性：

* **残余质量公式**：
  $$m_{\text{res}}(\beta) = 2.547 \times 10^{28} \exp(-17.559 \beta)$$
* **裸夸克质量**：
  * 轻夸克裸质量：$m_l$
  * 奇异夸克裸质量：$m_s$
* **重整化减除手征凝聚定义**：
  $$\Delta_{l, s} = \frac{1}{Z_m} \left( \langle \bar{\psi}\psi \rangle_l - \frac{m_l + m_{\text{res}}}{m_s + m_{\text{res}}} \langle \bar{\psi}\psi \rangle_s \right)$$
  该组合消除了二次发散项与残余手征破坏项，直接表征自发手征对称破缺在有限温下的恢复。

---

## 4. 平台拟合时间切片区间设置 (Fit Plateau Windows)

针对各信道在各温度下的信噪比特征，平台拟合区间统一设定为 $[t_{\min}, 25]$：

| Beta | AV ($t_{\min}$) | S ($t_{\min}$, multi/single) | Tt ($t_{\min}$) | PS ($t_{\min}$) | Xt ($t_{\min}$) | Vec ($t_{\min}$) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **13** | 16 | 14 / 8 | 10 | 8 | 10 | 16 |
| **15** | 13 | 10 / 8 | 10 | 8 | 10 | 13 |
| **17** | 14 | 8 / 8 | 10 | 8 | 10 | 14 |
| **18** | 16 | 18 / 8 | 10 | 8 | 10 | 16 |
| **20** | 15 | 19 / 8 | 10 | 8 | 10 | 15 |
| **23** | 16 | 8 / 8 | 10 | 8 | 10 | 16 |
| **30** | 16 | 8 / 8 | 10 | 8 | 10 | 16 |

* **统计 Binsize 设置**：
  * **多源 (Multi-Source)**：所有信道与温度默认 `binsize = 4`（唯 $\beta=4.18$ 信道 S 使用 `binsize = 5`）；
  * **单源 (Single-Source)**：$\beta \in \{13, 15, 17\}$ 使用 `binsize = 4`；$\beta \in \{18, 20, 23, 30\}$ 使用 `binsize = 5`。
