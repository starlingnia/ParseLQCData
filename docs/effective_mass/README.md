# Beta 4.17 介子 Effective Mass (有效质量) 全局分析报告与拟合窗口导航

> **核心目标**：提供所有 beta=4.17 目录下 6 个介子信道（Vec, AV, PS, Tt, Xt, S）的有效质量平台图，作为精确设定与校准 `config/fit_windows.txt` 中各 ensemble 及噪声信道窗口 `[x_start, x_end)` 的判定依据。
>
> **物理标度**：$a^{-1} = 153 \times 16 = 2448\text{ MeV}$，$T = 2448 / N_t\text{ MeV}$。

## 1. 当前 `config/fit_windows.txt` 配置现状与观察建议

| Ensemble | Ns x Nt | T (MeV) | 默认窗口 | 噪声道专用窗口 | 观测有效质量平台特征与建议 |
|:---|:---:|:---:|:---:|:---|:---|
| **32x12** | 32x12 | 204.0 | `8 17` | 无 | 平台在 x=12~16 表现优良；S 道噪声从 x=7 起增大，故 S 道采用 [6, 17) 可保持稳定。 |
| **32x14** | 32x14 | 174.9 | `6 14` | 无 | Vec/AV/PS 在 [12, 17) 稳定；S 道在近相变区信号衰减极快，放宽至 [4, 17) 避免 m→0 退化。 |
| **32x16** | 32x16 | 153.0 | `8 17` | 无 | 处于赝临界点 T=153 MeV 附近，Vec/AV 在 [12, 17) 形成良好平台；S 道涨落显著。 |
| **36x18** | 36x18 | 136.0 | 未设定 | 无 | Ns=36 对称点 half=18。Vec/AV/PS 在 [12, 19) 或 [14, 19) 展现清晰平坦基态；S/Tt 道需注意晚期涨落。 |
| **40x16** | 40x16 | 153.0 | `8 21` | 无 | Ns=40 对称点 half=20。默认区间 [16, 21) 对 Vec/AV 极佳；Tt/Xt 在小 ml 处早衰，放宽至 [8, 21)。 |
| **48x18** | 48x18 | 136.0 | `8 25` | 无 | Ns=48 对称点 half=24。默认 [16, 25) 极平坦；张量与标量道在对称点附近信噪比过低，专用窗口极为关键。 |
| **48x16** | 48x16 | 153.3 | 未设定 | 无 | 定标度物理点基准 (T=153.31 MeV)，对称点 half=24。 |

---

## 2. 格点规模全景概览卡片 (跨 4 种轻夸克质量对比)

点击链接即可打开大图，用于评估某个格点规模下所有夸克质量是否在同一窗口平坦：

- **32x12 全景总览卡片**: [overview_32x12.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_32x12.png)
- **32x14 全景总览卡片**: [overview_32x14.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_32x14.png)
- **32x16 全景总览卡片**: [overview_32x16.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_32x16.png)
- **36x18 全景总览卡片**: [overview_36x18.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_36x18.png)
- **40x16 全景总览卡片**: [overview_40x16.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_40x16.png)
- **48x16 全景总览卡片**: [overview_48x16.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_48x16.png)
- **48x18 全景总览卡片**: [overview_48x18.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/ensemble_overview/overview_48x18.png)

---

## 3. 按目录与频道分类存储的图片索引导航

### 3.1 按目录划分 (25 个目录)

#### 📂 目录：`32x12_b4.17_ms0.040m0.0020`
- **参数**：格点 $32\times 12$，半宽 $half=16$，$m_l=0.0020$，$T=204.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0020/meff_S.png)

#### 📂 目录：`32x12_b4.17_ms0.040m0.0035`
- **参数**：格点 $32\times 12$，半宽 $half=16$，$m_l=0.0035$，$T=204.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0035/meff_S.png)

#### 📂 目录：`32x12_b4.17_ms0.040m0.0070`
- **参数**：格点 $32\times 12$，半宽 $half=16$，$m_l=0.0070$，$T=204.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0070/meff_S.png)

#### 📂 目录：`32x12_b4.17_ms0.040m0.0120`
- **参数**：格点 $32\times 12$，半宽 $half=16$，$m_l=0.0120$，$T=204.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x12_b4.17_ms0.040m0.0120/meff_S.png)

#### 📂 目录：`32x14_b4.17_ms0.040m0.0020`
- **参数**：格点 $32\times 14$，半宽 $half=16$，$m_l=0.0020$，$T=174.86$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0020/meff_S.png)

#### 📂 目录：`32x14_b4.17_ms0.040m0.0035`
- **参数**：格点 $32\times 14$，半宽 $half=16$，$m_l=0.0035$，$T=174.86$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0035/meff_S.png)

#### 📂 目录：`32x14_b4.17_ms0.040m0.0070`
- **参数**：格点 $32\times 14$，半宽 $half=16$，$m_l=0.0070$，$T=174.86$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0070/meff_S.png)

#### 📂 目录：`32x14_b4.17_ms0.040m0.0120`
- **参数**：格点 $32\times 14$，半宽 $half=16$，$m_l=0.0120$，$T=174.86$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x14_b4.17_ms0.040m0.0120/meff_S.png)

#### 📂 目录：`32x16_b4.17_ms0.040m0.0020`
- **参数**：格点 $32\times 16$，半宽 $half=16$，$m_l=0.0020$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0020/meff_S.png)

#### 📂 目录：`32x16_b4.17_ms0.040m0.0035`
- **参数**：格点 $32\times 16$，半宽 $half=16$，$m_l=0.0035$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0035/meff_S.png)

#### 📂 目录：`32x16_b4.17_ms0.040m0.0070`
- **参数**：格点 $32\times 16$，半宽 $half=16$，$m_l=0.0070$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0070/meff_S.png)

#### 📂 目录：`32x16_b4.17_ms0.040m0.0120`
- **参数**：格点 $32\times 16$，半宽 $half=16$，$m_l=0.0120$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/32x16_b4.17_ms0.040m0.0120/meff_S.png)

#### 📂 目录：`36x18_b4.17_ms0.040m0.0020`
- **参数**：格点 $36\times 18$，半宽 $half=18$，$m_l=0.0020$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0020/meff_S.png)

#### 📂 目录：`36x18_b4.17_ms0.040m0.0035`
- **参数**：格点 $36\times 18$，半宽 $half=18$，$m_l=0.0035$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0035/meff_S.png)

#### 📂 目录：`36x18_b4.17_ms0.040m0.0070`
- **参数**：格点 $36\times 18$，半宽 $half=18$，$m_l=0.0070$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0070/meff_S.png)

#### 📂 目录：`36x18_b4.17_ms0.040m0.0120`
- **参数**：格点 $36\times 18$，半宽 $half=18$，$m_l=0.0120$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/36x18_b4.17_ms0.040m0.0120/meff_S.png)

#### 📂 目录：`40x16_b4.17_ms0.040m0.0020`
- **参数**：格点 $40\times 16$，半宽 $half=20$，$m_l=0.0020$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0020/meff_S.png)

#### 📂 目录：`40x16_b4.17_ms0.040m0.0035`
- **参数**：格点 $40\times 16$，半宽 $half=20$，$m_l=0.0035$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0035/meff_S.png)

#### 📂 目录：`40x16_b4.17_ms0.040m0.0070`
- **参数**：格点 $40\times 16$，半宽 $half=20$，$m_l=0.0070$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0070/meff_S.png)

#### 📂 目录：`40x16_b4.17_ms0.040m0.0120`
- **参数**：格点 $40\times 16$，半宽 $half=20$，$m_l=0.0120$，$T=153.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/40x16_b4.17_ms0.040m0.0120/meff_S.png)

#### 📂 目录：`48x18_b4.17_ms0.040m0.0020`
- **参数**：格点 $48\times 18$，半宽 $half=24$，$m_l=0.0020$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0020/meff_S.png)

#### 📂 目录：`48x18_b4.17_ms0.040m0.0035`
- **参数**：格点 $48\times 18$，半宽 $half=24$，$m_l=0.0035$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0035/meff_S.png)

#### 📂 目录：`48x18_b4.17_ms0.040m0.0070`
- **参数**：格点 $48\times 18$，半宽 $half=24$，$m_l=0.0070$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0070/meff_S.png)

#### 📂 目录：`48x18_b4.17_ms0.040m0.0120`
- **参数**：格点 $48\times 18$，半宽 $half=24$，$m_l=0.0120$，$T=136.00$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x18_b4.17_ms0.040m0.0120/meff_S.png)

#### 📂 目录：`48x16b4.17`
- **参数**：格点 $48\times 16$，半宽 $half=24$，$m_l=0.0010$，$T=153.31$ MeV
- **6 信道全景图**：[meff_all_channels.png](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_all_channels.png)
- **数据表格**：[meff_data.csv](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_data.csv)
- **单信道图**：
  [Vec](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_Vec.png) | [AV](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_AV.png) | [PS](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_PS.png) | [Tt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_Tt.png) | [Xt](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_Xt.png) | [S](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_directory/48x16b4.17/meff_S.png)

---

### 3.2 按信道划分 (夸克质量对比图)

#### 🏷️ 信道：`Vec` (Vector (Vec / V, ρ/ω))
在不同格点规模下对比 4 种夸克质量平台：
- [32x12 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/32x12_compare_ml.png) | [32x14 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/32x14_compare_ml.png) | [32x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/32x16_compare_ml.png) | [36x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/36x18_compare_ml.png) | [40x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/40x16_compare_ml.png) | [48x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/48x16_compare_ml.png) | [48x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Vec/48x18_compare_ml.png)

#### 🏷️ 信道：`AV` (Axial-Vector (AV / A, a₁))
在不同格点规模下对比 4 种夸克质量平台：
- [32x12 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/32x12_compare_ml.png) | [32x14 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/32x14_compare_ml.png) | [32x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/32x16_compare_ml.png) | [36x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/36x18_compare_ml.png) | [40x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/40x16_compare_ml.png) | [48x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/48x16_compare_ml.png) | [48x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/AV/48x18_compare_ml.png)

#### 🏷️ 信道：`PS` (Pseudo-Scalar (PS / Ps, π/η))
在不同格点规模下对比 4 种夸克质量平台：
- [32x12 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/32x12_compare_ml.png) | [32x14 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/32x14_compare_ml.png) | [32x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/32x16_compare_ml.png) | [36x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/36x18_compare_ml.png) | [40x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/40x16_compare_ml.png) | [48x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/48x16_compare_ml.png) | [48x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/PS/48x18_compare_ml.png)

#### 🏷️ 信道：`Tt` (Tensor (Tt / T, b₁))
在不同格点规模下对比 4 种夸克质量平台：
- [32x12 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/32x12_compare_ml.png) | [32x14 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/32x14_compare_ml.png) | [32x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/32x16_compare_ml.png) | [36x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/36x18_compare_ml.png) | [40x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/40x16_compare_ml.png) | [48x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/48x16_compare_ml.png) | [48x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Tt/48x18_compare_ml.png)

#### 🏷️ 信道：`Xt` (Axial-Tensor (Xt / X, h₁))
在不同格点规模下对比 4 种夸克质量平台：
- [32x12 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/32x12_compare_ml.png) | [32x14 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/32x14_compare_ml.png) | [32x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/32x16_compare_ml.png) | [36x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/36x18_compare_ml.png) | [40x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/40x16_compare_ml.png) | [48x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/48x16_compare_ml.png) | [48x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/Xt/48x18_compare_ml.png)

#### 🏷️ 信道：`S` (Scalar (S, σ/a₀))
在不同格点规模下对比 4 种夸克质量平台：
- [32x12 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/32x12_compare_ml.png) | [32x14 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/32x14_compare_ml.png) | [32x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/32x16_compare_ml.png) | [36x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/36x18_compare_ml.png) | [40x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/40x16_compare_ml.png) | [48x16 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/48x16_compare_ml.png) | [48x18 对比图](file:///Users/junxiongnie/code/algo/ParseLQCData/docs/effective_mass/by_channel/S/48x18_compare_ml.png)

---

## 4. 如何更新 `config/fit_windows.txt`

1. 打开 [`config/fit_windows.txt`](file:///Users/junxiongnie/code/algo/ParseLQCData/config/fit_windows.txt)；
2. 根据上述图表中观测到的平台区间修改或新增行，例如为 `36x18` 指定：
   ```text
   36x18              14       19
   36x18  S            6       19
   ```
3. 修改保存后，直接执行：
   ```bash
   .venv/bin/python scripts/update_ccor_data_and_plots.py
   ```
   系统将自动采用新窗口重拟合介子质量并更新全部图表！
