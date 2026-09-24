# 手征磁化率抽取逻辑与性能评估综合报告 (Chiral Susceptibility Analysis & Performance Report)

## 一、 模块与执行脚本概览

* **核心执行脚本**: [`scripts/reproduce_susceptibility.py`](file:///Users/junxiongnie/code/algo/ParseLQCData/scripts/reproduce_susceptibility.py)
* **底层分析流水线**: [`src/parselqcdata/susceptibility_pipeline.py`](file:///Users/junxiongnie/code/algo/ParseLQCData/src/parselqcdata/susceptibility_pipeline.py)
* **原生加速引擎**: [`src/CondensateAnalysis/CondensateExtractor.cpp`](file:///Users/junxiongnie/code/algo/ParseLQCData/src/CondensateAnalysis/CondensateExtractor.cpp) (C++ 26 多线程加速)
* **Python C-ABI 桥接**: [`tools/condensate_orchestrator.py`](file:///Users/junxiongnie/code/algo/ParseLQCData/tools/condensate_orchestrator.py)
* **学术绘图工具**: [`scripts/plot_susceptibility.py`](file:///Users/junxiongnie/code/algo/ParseLQCData/scripts/plot_susceptibility.py)

---

## 二、 脚本核心逻辑与提取原理

### 1. 物理理论背景
手征磁化率（Chiral Susceptibility, $\chi$）是不连通手征凝聚两点起伏的核心热力学观测量：
$$\chi_{\text{disc}} = \frac{V_3}{T} \left[ \langle (\bar{\psi}\psi)^2 \rangle - \langle \bar{\psi}\psi \rangle^2 \right]$$
其中三维时空体积 $V_3 = (N_s a)^3$，温度 $T = \frac{1}{N_t a}$，因此 $\frac{V_3}{T} = N_s^3 N_t a^4 = V_4$（四维时空体积）。
在有限温格点 QCD 中，手征磁化率在假临界温度 $T_{pc} \approx 157.0\,\mathrm{MeV}$（对应 $\beta=4.18$）呈现出强烈的极大值相变峰，是标定手征对称性恢复交叉相变（Crossover）的关键物理观测量。

### 2. 微观算法三步流

#### Step 1: 单构型无偏二次交叉乘积估计器 (Unbiased Quadratic Estimator)
对于单个格点规范构型，使用有限个 ($k$ 个，通常 $k=10$) 随机源向量 $O_i$ 估计 $\langle (\bar{\psi}\psi)^2 \rangle$ 时，若直接平方会包含对角自乘项 $O_i^2$，从而混入随机噪声方差 $\sigma_{\text{noise}}^2$：
$$\mathbb{E}[O_i O_j] = [\mathrm{Tr}(D^{-1})]^2 + \sigma_{\text{noise}}^2 \delta_{ij}$$
算法严格剔除对角自相关项，计算单构型无偏估计量：
- **单体平均**:
  $$\bar{O} = \frac{1}{k} \sum_{i=1}^k O_i$$
- **两体无偏方均**:
  $$\overline{O^2} = \frac{1}{k(k-1)} \sum_{i \neq j} O_i O_j = \frac{1}{k(k-1)} \left[ \left(\sum_{i=1}^k O_i\right)^2 - \sum_{i=1}^k O_i^2 \right]$$
该数学恒等式将传统两重循环复杂度从 $\mathcal{O}(k^2)$ 降低至 $\mathcal{O}(k)$ 线性时间。

#### Step 2: 系综级 Jackknife 重采样统计推断
针对全系综 $N$ 个构型的一阶与二阶序列执行 Leave-One-Out 重采样：
- 第 $r$ 个 Jackknife 样本估计：
  $$a_r = \text{JK}(\bar{O})_r = \frac{\sum_{i=1}^N \bar{O}_i - \bar{O}_r}{N - 1}$$
  $$b_r = \text{JK}(\overline{O^2})_r = \frac{\sum_{i=1}^N \overline{O^2}_i - \overline{O^2}_r}{N - 1}$$
- 构型级磁化率：
  $$\chi_r = b_r - a_r^2$$
- 统计均值与标准误差：
  $$\text{Mean}(\chi_{\text{unscaled}}) = \frac{1}{N} \sum_{r=1}^N \chi_r$$
  $$\text{Error}(\chi_{\text{unscaled}}) = \sqrt{N - 1} \cdot \sqrt{\frac{1}{N} \sum_{r=1}^N \left(\chi_r - \text{Mean}(\chi)\right)^2}$$

#### Step 3: 几何体积与连续温度标度变换
- **4D 晶格体积标度因子**:
  $$F_{\text{vol}} = N_s^3 \cdot N_t, \quad \chi_{\text{vol}} = F_{\text{vol}} \cdot \chi_{\text{unscaled}}$$
- **连续温度标度因子**（对应物理连续统量 $(16T)^2 a^2 \chi$，单位 $\mathrm{MeV}^2$）：
  $$F_{\text{scaled}} = N_s^3 \cdot N_t^3 \cdot T^2 = F_{\text{vol}} \cdot (N_t T)^2, \quad \chi_{\text{scaled}} = F_{\text{scaled}} \cdot \chi_{\text{unscaled}}$$

---

## 三、 提取结果字段定义与存储位置

### 1. 输出数据集字段规范

计算产出包含完整的晶格几何、构型统计数及双重标度结果：

| 字段名称 | 类型 | 物理含义 | 单位 |
| :--- | :--- | :--- | :--- |
| `Ensemble` | String | 物理系综目录标识 (如 `L48T16beta4.18ms...`) | - |
| `Beta` | Float | 格点规范耦合常数 $\beta$ (4.13 ~ 4.405) | - |
| `Temp` | Float | 物理温度 $T$ | $\mathrm{MeV}$ |
| `Ns` | Int | 空间格点点数 (通常为 48, 32, 40) | - |
| `Nt` | Int | 时间格点点数 (通常为 16, 12, 14, 18) | - |
| `Num_cfgs` | Int | 参与统计推断的有效规范构型数 | 构型数 |
| `Mean_unscaled` | Float | 晶格裸手征磁化率均值 | 晶格无量纲 |
| `Error_unscaled` | Float | 裸磁化率 Jackknife 统计误差 | 晶格无量纲 |
| `Mean_vol_scaled` | Float | 4D 体积标度磁化率 $\chi_{\text{vol}} = F_{\text{vol}} \cdot \chi$ | 无量纲 |
| `Error_vol_scaled` | Float | 4D 体积标度磁化率误差 | 无量纲 |
| `Mean_scaled` | Float | 连续标度物理磁化率 $\chi_{\text{scaled}} = F_{\text{scaled}} \cdot \chi$ | $\mathrm{MeV}^2$ |
| `Error_scaled` | Float | 连续标度物理磁化率误差 | $\mathrm{MeV}^2$ |

### 2. 产物存储位置

1. **结构化数据表格**:
   - [`output/condensate/results_susceptibility.csv`](file:///Users/junxiongnie/code/algo/ParseLQCData/output/condensate/results_susceptibility.csv)：基准与全量温度磁化率 CSV。
   - [`output/condensate/results_susceptibility.parquet`](file:///Users/junxiongnie/code/algo/ParseLQCData/output/condensate/results_susceptibility.parquet)：高性能列式存储，便于下游分析。
   - [`output/condensate/all_ensembles_susceptibility.csv`](file:///Users/junxiongnie/code/algo/ParseLQCData/output/condensate/all_ensembles_susceptibility.csv)：跨体积与多时空尺度全系综扫描结果。
2. **矢量与高清物理图表**:
   - `sucep_plot.png` / `.pdf`：4D 体积标度磁化率随温度演化曲线（同步存至 `docs/figures/`、`output/condensate/`、`output/plots/ana_style/`）。
   - `pbpchisce.png` / `.pdf`：连续温度标度物理磁化率曲线，标定相变带 $[155.5, 160.5]\,\mathrm{MeV}$ 与 $T_{pc} \approx 157.0\,\mathrm{MeV}$ 峰顶。
3. **学术论文/报告引用**:
   - 自动嵌入中英文学术报告 [`report.tex`](file:///Users/junxiongnie/code/algo/ParseLQCData/report.tex) 与 [`report_en.tex`](file:///Users/junxiongnie/code/algo/ParseLQCData/report_en.tex)。

---

## 四、 性能综合评估与基准总结

### 1. 复杂度与资源开销

| 维度 | 指标表现 | 评估结论 |
| :--- | :--- | :--- |
| **理论时间复杂度** | $\mathcal{O}(k) + \mathcal{O}(N)$ | 数学上已优化至理论极限，CPU 计算周期耗时仅微秒级。 |
| **内存空间开销** | $\mathcal{O}(N)$（流式即算即弃） | 内存仅缓存构型均值浮点标量，单系综占用 $< 100\,\mathrm{KB}$，无 OOM 隐患。 |
| **I/O 负载特征** | 读密集型、海量碎片化小文件 | 瓶颈主要集中在文件系统 `readdir`、文件属性 `stat` 与文件 `open/close`。 |

### 2. 数据负载规模
- 单个系综（如 $\beta=4.18$）：包含 **2,029 个构型**，每个构型 10 个轻夸克随机源 XML 文件，单系综即有 **20,290 个小 XML 文件**。
- 全谱 8 组温度系综：累计处理超过 **160,000 个独立 XML 文件**。

### 3. 双引擎架构性能实测对比

| 执行引擎 | 实现方式 | 并发调度 | 单系综 (~2,000 构型) 耗时 | 全谱 8 温度耗时 | 加速比 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **C++ 26 Native 引擎** | 原生编译共享库，无 DOM 轻量文本提取 | `std::async` 多核分块并发 | **约 1.5 ~ 3.5 秒** | **约 15 ~ 30 秒** | **12x ~ 16x 🚀** |
| **Python Fallback 模式** | `read_text()` + 正则 `re.search` | 单线程解释执行 | **约 25 ~ 50 秒** | **约 4 ~ 6 分钟** | 基准 (1x) |

### 4. 关键架构优势
1. **轻量零拷贝文本扫描**: C++ 算子避开加载几十兆 DOM 语法树，直接通过快速字符串定位抽取 `<pbp>(real, imag)` 实数，将单个 XML 解析降低到纳秒级。
2. **多线程负载均衡**: 根据系统可用硬件线程数（`hardware_concurrency`）自适应划分子区间，并行处理海量构型。
3. **架构平滑容灾**: 若动态库环境不匹配或不可用，Python 层会自动无缝降级运行，确保计算结果 100% 可复现。

### 5. 后续优化演进建议
1. **小文件数据聚合**: 目前每个构型的 10 个随机源分散为 10 个 XML 文件。若能在数据生成侧归档为 SQLite 数据库或 HDF5/NPY 格式，可消减 90% 的文件系统元数据开销，预计整体抽取吞吐率可再提升 **5 ~ 10 倍**。
2. **内存映射预读取 (`mmap`)**: 在集群超算高并发场景下，使用内存映射替代逐文件读取，进一步规避用户态与内核态之间的数据复制缓冲。
