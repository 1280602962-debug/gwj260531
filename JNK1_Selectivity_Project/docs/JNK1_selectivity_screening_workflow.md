# JNK1 选择性抑制剂筛选：CADD/AIDD 完整流程思路

> **版本**: 1.0  
> **适用场景**: 基于 ChEMBL 等公开数据，通过计算与机器学习方法发现 JNK1 亚型选择性小分子抑制剂  
> **核心策略**: 先比较 JNK1/2/3 数据集相似性 → 逐步整合 5 类 SAR 方法 → 构建选择性预测模型 → 百万分子库虚拟筛选 → SHAP 可解释性验证

---

## 目录

1. [项目背景与设计原则](#1-项目背景与设计原则)
2. [Phase 0：数据获取与标准化](#2-phase-0数据获取与标准化)
3. [Phase 1：三靶点数据集相似性比较](#3-phase-1三靶点数据集相似性比较)
4. [Phase 2：五类 SAR 方法逐步整合](#4-phase-2五类-sar-方法逐步整合)
5. [Phase 3：选择性模型构建与评估](#5-phase-3选择性模型构建与评估)
6. [Phase 4：百万分子库虚拟筛选](#6-phase-4百万分子库虚拟筛选)
7. [Phase 5：可解释性分析（SHAP 等）](#7-phase-5可解释性分析shap-等)
8. [Phase 6：实验验证建议](#8-phase-6实验验证建议)
9. [决策树与里程碑检查点](#9-决策树与里程碑检查点)
10. [参考文献](#10-参考文献)

---

## 1. 项目背景与设计原则

### 1.1 生物学与成药背景

JNK（c-Jun N-terminal kinase）家族包含三个高度同源的丝/苏氨酸激酶亚型：

| 亚型 | 基因 | 组织分布 | 成药考量 |
|------|------|----------|----------|
| JNK1 | MAPK8 | 广泛表达 | 代谢、炎症、纤维化主靶点 |
| JNK2 | MAPK9 | 广泛表达 | 免疫/CNS 相关，选择性 off-target |
| JNK3 | MAPK10 | 主要 CNS | 神经保护相关，外周选择性通常期望 JNK3 无活性 |

**JNK1 选择性抑制的逻辑**：JNK1 在胰岛素抵抗、IPF、NASH 等适应症中有明确证据，而 pan-JNK 或 JNK2 偏向抑制可能带来不必要的系统毒性。已知 JNK1 偏向化合物 **CC-90001** 已进入临床（IPF），说明该方向具有成药可行性 [Bennett et al., 2021]。

### 1.2 计算策略核心原则

| 原则 | 说明 |
|------|------|
| **配对优先** | 同一分子在多个亚型上的测定数据直接编码选择性信息 |
| **迁移学习** | JNK1 数据少时，从 JNK2/3 共享表示中迁移 |
| **Scaffold split** | 按骨架划分训练/测试集，避免化学相似性泄漏 [Wu & Rajpal, 2020] |
| **选择性显式建模** | 不仅预测绝对活性，更预测 ΔpActivity |
| **Benchmark 回测** | CC-90001、SP600125 等已知化合物必须被正确排序 |
| **可解释性闭环** | SHAP 结果需与 SAR/MMP 及结构差异残基一致 |

---

## 2. Phase 0：数据获取与标准化

### 2.1 ChEMBL 靶点定义

| 亚型 | ChEMBL ID | 基因 | UniProt | 靶点名称 |
|------|-----------|------|---------|----------|
| JNK1 | **CHEMBL2276** | MAPK8 | P45983 | Mitogen-activated protein kinase 8 |
| JNK2 | **CHEMBL4179** | MAPK9 | P45984 | Mitogen-activated protein kinase 9 |
| JNK3 | **CHEMBL2637** | MAPK10 | P53779 | Mitogen-activated protein kinase 10 |

> **注意**：搜索 "JNK" 会返回 pan-JNK 或未区分亚型的测定；必须按 `target_chembl_id` 精确过滤。

### 2.2 数据下载流程

```
ChEMBL API / chembl_webresource_client
    ↓
按 CHEMBL2276 / 4179 / 2637 分别拉取 activity 表
    ↓
合并 assay、compound、target 元数据
    ↓
输出: jnk1_raw.csv, jnk2_raw.csv, jnk3_raw.csv
```

**推荐保留字段**：
- `molecule_chembl_id`, `canonical_smiles`
- `standard_type`, `standard_value`, `standard_units`
- `pchembl_value`, `assay_type`, `assay_description`
- `assay_confidence_score`, `bao_label`

**补充数据源**（强烈建议）：
- **BindingDB** [Gilson et al., 2016]
- **PubChem BioAssay** [Kim et al., 2023]
- **Davis kinase panel** [Davis et al., 2011]（含多激酶选择性矩阵）
- 文献/专利手工整理（CC-90001、JNK-IN-8 等）

### 2.3 数据清洗规则

#### 2.3.1 活性类型统一

保留：`IC50`, `Ki`, `Kd`, `EC50`（biochemical binding/functional）

统一单位 → **nM**，计算 **pActivity = −log10(M)**

#### 2.3.2 质量过滤

| 规则 | 阈值 |
|------|------|
| Assay confidence | ≥ 6（ChEMBL 单蛋白靶点） |
| 精确关系 | 优先 `=`；`>` 标记为 censored 或剔除 |
| 重复测定 | 同分子-靶点-条件取 **几何平均** |
| 无效 SMILES | RDKit 无法解析则剔除 [Landrum, 2024] |
| PAINS / 反应性 | 剔除 [Baell & Holloway, 2010] |

#### 2.3.3 二分类标签（用于选择性分类）

```text
Active:     pActivity ≥ 6.5  (~316 nM)
Inactive:     pActivity < 5.5  (~3.16 µM)
Intermediate: 5.5–6.5（可单独处理或剔除）

JNK1-selective (Class 1):
    pActivity_JNK1 ≥ 6.5
    AND pActivity_JNK2 < 5.5
    AND pActivity_JNK3 < 5.5
    OR Δ12 = pAct_JNK1 − pAct_JNK2 ≥ 1.0
    OR Δ13 = pAct_JNK1 − pAct_JNK3 ≥ 1.0
```

#### 2.3.4 构建配对数据集（Paired Set）

对至少在 **2 个 JNK 亚型**上有测定数据的分子：

```text
paired_set.csv:
  molecule_chembl_id | smiles | pAct_JNK1 | pAct_JNK2 | pAct_JNK3 | delta_12 | delta_13 | sel_class
```

这是后续 **选择性模型** 和 **MMP 分析** 的核心输入。

---

## 3. Phase 1：三靶点数据集相似性比较

> **目的**：在建模前理解 JNK1/2/3 数据集的化学空间重叠程度、活性分布差异和数据可迁移性，为后续 MTL 和选择性建模提供依据。

### 3.1 分析维度总览

```
                ┌─────────────────────────────────────┐
                │     JNK1 / JNK2 / JNK3 数据集        │
                └─────────────────────────────────────┘
                    │           │           │
         ┌──────────┘           │           └──────────┐
         ▼                      ▼                      ▼
  3.2 化学空间重叠        3.3 骨架/系列重叠      3.4 活性分布比较
         │                      │                      │
         └──────────────────────┼──────────────────────┘
                                ▼
                    3.5 配对化合物选择性分布
                                ▼
                    3.6 数据可迁移性结论 → 指导建模策略
```

### 3.2 化学空间相似性

#### 方法 A：分子指纹 + Tanimoto 相似性

- 指纹：**Morgan (ECFP4, radius=2, 2048 bits)** [Rogers & Hahn, 2010]
- 对每个数据集随机采样 ≥5000 分子（若不足则全取）
- 计算 **数据集间平均 Tanimoto 相似性矩阵**：

```text
Sim(JNK1, JNK2) = mean_{i∈JNK1, j∈JNK2} Tanimoto(fp_i, fp_j)
```

- 可视化：三数据集两两相似性热图

#### 方法 B：降维可视化（UMAP / PCA）

- 输入：Morgan FP 或 **MAP4** 指纹
- 将 JNK1/2/3 全部分子投影到 2D
- 着色：靶点来源 + 活性（active/inactive）
- **解读**：
  - 高度重叠 → 适合 MTL / 迁移学习
  - JNK1 独占区域 → 需警惕外推风险（applicability domain）

#### 方法 C：Fréchet ChemNet Distance（可选进阶）

比较两个数据集在深度学习化学空间中的分布距离 [Preuer et al., 2018]。

### 3.3 骨架与系列重叠

#### Bemis-Murcko 骨架分析 [Bemis & Murcko, 1996]

```text
对每个分子提取 Murcko scaffold
统计:
  - 各数据集 unique scaffold 数
  - 三数据集共享 scaffold 数及比例
  - JNK1 独有 scaffold 比例
```

#### 系列（Series）聚类

- 按 scaffold 或 Butina 聚类（Tanimoto > 0.7）
- 比较三靶点在各系列中的 **活性富集差异**
- 识别 **JNK1-selective series**（某系列在 JNK1 高活性、JNK2/3 低活性）

### 3.4 活性分布比较

| 统计量 | JNK1 | JNK2 | JNK3 | 意义 |
|--------|------|------|------|------|
| 化合物数 | n₁ | n₂ | n₃ | 数据量差异 |
| pActivity 中位数 | | | | 整体活性水平 |
| Active 比例 (≥6.5) | | | | 命中率 |
| Assay type 组成 | B vs F | | | 测定类型偏差 |

**Kolmogorov-Smirnov 检验**：比较 JNK1 vs JNK2、JNK1 vs JNK3 的 pActivity 分布是否显著不同。

### 3.5 配对化合物分析（最关键）

对 `paired_set` 中的分子：

```text
Δ12 分布直方图
Δ13 分布直方图
JNK1-selective / pan-JNK / JNK2-biased 比例
```

**活性 cliffs（活性悬崖）**：
- 定义：Tanimoto > 0.8 但 Δselectivity > 1 log unit 的分子对
- 这些 cliffs 是 **选择性模型的关键训练信号**

### 3.6 相似性分析结论模板

完成 Phase 1 后，填写以下决策表：

| 观察 | 建模策略调整 |
|------|--------------|
| JNK1/2/3 化学空间高度重叠 (>0.4 mean Tanimoto) | 启用 MTL，共享编码器 |
| JNK1 数据量 << JNK2/3 | JNK1 任务权重上调；使用 transfer learning |
| 配对化合物 > 200 | 可直接训练选择性分类器 |
| 配对化合物 < 100 | 以回归 Δ 为主，分类为辅；加大 SBDD 权重 |
| 共享 scaffold 比例高 | Scaffold split 必须严格执行 |
| JNK1 存在独有 chemotype | 单独标注，SHAP 重点分析 |

**对应脚本**: `scripts/02_dataset_similarity.py`

---

## 4. Phase 2：五类 SAR 方法逐步整合

> **逻辑**：由浅入深、由解释到预测，每一类的输出作为下一类的输入或约束。

```
① 经典 SAR/MMP ──→ 选择性药效团假设、关键取代基
        ↓
② 单靶点 ML ──→ 各亚型 baseline 预测能力
        ↓
③ 多任务 MTL ──→ 共享表示 + 迁移 learning
        ↓
④ 选择性专用模型 ──→ Δ12/Δ13 回归 + 分类
        ↓
⑤ SBDD + ML 联合 ──→ 对接特征 + 2D/3D 描述符融合
```

---

### 4.1 方法一：经典 SAR / MMP 分析（解释层）

**目的**：从已知数据中提炼 "什么结构变化导致 JNK1 选择性" 的可解释规则。

#### 4.1.1 R-group 分解

- 识别系列共同母核（Murcko scaffold / 最大公共子结构）
- 分解 R1、R2 等取代基
- 绘制 **R-group × pActivity** 热图（分 JNK1/2/3 三张）

#### 4.1.2 Matched Molecular Pair (MMP) 分析

- 使用 **mmpdb** 或 RDKit 实现 [Dalke et al., 2002; Dossetter & Griffen, 2013]
- 在 paired_set 中找 MMP：`CC(=O)N → CC(=O)NH` 等
- 计算 **ΔpActivity** 和 **Δselectivity** 的 transform 规则

```text
示例输出:
  Transform: [*:1]F → [*:1]Cl  |  median Δ12 = +0.8  |  n = 15
  Transform: [*:1]Me → [*:1]Et  |  median Δ13 = −0.6  |  n = 8
```

#### 4.1.3 活性悬崖 (Activity Cliff) 识别

- 定义：Tanimoto ≥ 0.85 且 |ΔpAct| ≥ 1.0 的分子对
- 分别找 JNK1 cliffs 和 **selectivity cliffs**
- 3D 构象叠合分析 cliff 对的结合模式差异

#### 4.1.4 输出

- `sar_mmp_rules.csv`：选择性相关 transform 规则
- `sar_pharmacophore_hypothesis.json`：JNK1 选择性药效团
- 为 Phase 5 SHAP 分析提供 **先验假设**

**对应脚本**: `scripts/03_sar_analysis.py`

---

### 4.2 方法二：单靶点机器学习（Baseline 层）

**目的**：评估各亚型独立模型的预测上限，作为 MTL 的对比 baseline。

#### 4.2.1 特征

| 特征类型 | 具体方法 |
|----------|----------|
| 2D 指纹 | Morgan (2048), MACCS, RDKit FP |
| 描述符 | RDKit 2D descriptors (MW, LogP, TPSA, HBD/HBA...) |
| 3D 描述符（可选） | ROCS shape, USR |

#### 4.2.2 模型

| 模型 | 用途 |
|------|------|
| **Random Forest** | 小样本 baseline |
| **XGBoost / LightGBM** | 主力 baseline，支持 SHAP |
| **SVM** | 对照 |

#### 4.2.3 验证策略

- **Scaffold split**（必须）：按 Bemis-Murcko 骨架 80/10/10 划分 [Wu & Rajpal, 2020]
- **Random split**（仅作对照，不作为主结论）
- **Nested CV** 用于超参数优化

#### 4.2.4 评估指标

| 任务 | 指标 |
|------|------|
| 回归 (pActivity) | RMSE, MAE, R², Spearman ρ |
| 分类 (active/inactive) | AUC-ROC, AUC-PR, F1, MCC |
| 一致性 | 与实验值误差 < 0.5 log unit 的比例 |

#### 4.2.5 预期结果

- JNK2/JNK3 模型 R² 通常高于 JNK1（数据量差异）
- 若 JNK1 R² < 0.3 → 确认必须引入 MTL

---

### 4.3 方法三：多任务学习 MTL（核心预测层）

**目的**：利用 JNK 家族 SAR 高度相关的特点，共享底层分子表示，缓解 JNK1 数据不足。

#### 4.3.1 架构

```text
Input: SMILES
    ↓
Shared Encoder (GNN or FP-MLP)
    ├── Task Head: JNK1 pActivity (regression)
    ├── Task Head: JNK2 pActivity (regression)
    └── Task Head: JNK3 pActivity (regression)
```

#### 4.3.2 推荐实现

| 框架 | 模型 | 参考 |
|------|------|------|
| **Chemprop** | D-MPNN multitask | [Yang et al., 2019; Heid et al., 2024] |
| **DeepChem** | GCN multitask | [Ramsundar et al., 2019] |
| **XGBoost** | 3-output multi-output regressor | [Chen & Guestrin, 2016] |

#### 4.3.3 缺失值处理

- JNK1/2/3 数据天然不完整（同一分子不一定三个亚型都测过）
- Chemprop 原生支持 **missing target values**
- 训练时每个 task 仅用该 task 有标签的样本计算 loss

#### 4.3.4 任务权重

```text
loss = w1 * L_JNK1 + w2 * L_JNK2 + w3 * L_JNK3

建议: w1 = 2.0, w2 = 1.0, w3 = 1.0  (上调 JNK1 权重)
```

#### 4.3.5 评估

- 分 task 报告 scaffold-split 指标
- **Benchmark 回测**：CC-90001、SP600125、JNK-IN-8 的预测 vs 实验

---

### 4.4 方法四：选择性专用模型（目标层）

**目的**：直接预测 JNK1 选择性，而非间接相减。

#### 4.4.1 选择性回归

```text
Target:
  delta_12 = pAct_JNK1 − pAct_JNK2
  delta_13 = pAct_JNK1 − pAct_JNK3
  delta_min = pAct_JNK1 − max(pAct_JNK2, pAct_JNK3)

Input: SMILES (+ 可选 MTL 预测值作为 feature)
Model: XGBoost / LightGBM / GNN
Data: 仅 paired_set 中有 ≥2 个亚型数据的分子
```

#### 4.4.2 选择性分类

```text
Class 1 (JNK1-selective): delta_min ≥ 1.0 AND pAct_JNK1 ≥ 6.5
Class 0 (non-selective):  其余

Model: XGBoost classifier / Random Forest
Metric: AUC-ROC, AUC-PR, Precision@Top100
```

#### 4.4.3 MMP-based 选择性模型（补充）

- 输入：MMP transform + 母核描述符
- 输出：Δselectivity
- 优势：可解释性极强，适合与 SHAP 交叉验证

#### 4.4.4 集成策略（推荐最终模型）

```text
Final_Score = α * MTL_pAct_JNK1
            + β * SelModel_delta_min
            + γ * DockScore_JNK1
            − δ * max(MTL_pAct_JNK2, MTL_pAct_JNK3)
```

权重 α, β, γ, δ 在 validation set 上 grid search 优化。

---

### 4.5 方法五：结构导向 + ML 联合（精筛层）

**目的**：利用 JNK1/2/3 蛋白结构差异，为选择性提供 3D 层面解释。

#### 4.5.1 结构准备

| 亚型 | 推荐 PDB | 共晶配体 | 备注 |
|------|----------|----------|------|
| JNK1 | 3ELJ, 4L7F | GS7; AX13587 (1V5) | DFG-in；3ELJ (Q=9.3) + 4L7F (Q=9.8, 1.95 Å, 单链) |
| JNK2 | 3E7O | indazole (35F) | DFG-in；**仅 chain A**（sole receptor） |
| JNK3 | 3TTI, 4WHZ | CC-930 (KBI); pyrazole (3NL) | DFG-in；4WHZ 占 selectivity pocket |

- 序列比对找 **差异残基**（尤其 gatekeeper、hinge、back pocket）
- CC-90001 为 JNK1 偏选择性临床候选 [Bennett et al., 2021]，但**目前无公开 JNK1 共晶 PDB**；对接验证可借用同系列 CC-930/3TTI 结合模式

#### 4.5.2 Ensemble Docking

```text
对每个候选分子:
  Dock → JNK1 (top pose, score_1)
  Dock → JNK2 (top pose, score_2)
  Dock → JNK3 (top pose, score_3)

Selectivity docking score:
  ΔDock = score_JNK2/3 − score_JNK1  (越正越 JNK1-selective)
```

推荐工具：**GNINA** [McNutt et al., 2021]（CNN rescoring）或 **Glide SP/XP** [Friesner et al., 2004]

#### 4.5.3 Interaction Fingerprint (IFP)

- 提取每 pose 的 IFP（氢键、疏水、盐桥）
- 拼接 `[IFP_JNK1 | IFP_JNK2 | IFP_JNK3 | Morgan FP]` 作为 ML 输入
- 训练 **Activity Cliff-aware** 模型 [Czarnecki et al., 2018]

#### 4.5.4 选择性结构假设验证

- 对 SHAP 识别的重要子结构，检查其在 JNK1 vs JNK2/3 口袋中的 **空间冲突/互补**
- 形成 "结构–SAR–ML" 三角验证闭环

---

## 5. Phase 3：选择性模型构建与评估

### 5.1 推荐建模路线（综合 Phase 1–2 结果）

```text
Stage A: XGBoost baseline (单靶点 × 3)           → baseline 指标
Stage B: Chemprop MTL (JNK1/2/3 回归)            → 共享表示
Stage C: XGBoost selective (Δ_min 回归 + 分类)  → 选择性
Stage D: Stacking ensemble                        → 最终模型
Stage E: SBDD rescoring (Top 5000)               → 精筛
```

### 5.2 数据划分（严格）

```text
Split strategy: Scaffold-based (Bemis-Murcko)
  Train:  80%  (scaffold groups)
  Val:    10%  (hyperparameter tuning)
  Test:   10%  (final report, touch once)

注意: 配对化合物中同一 scaffold 的所有分子必须在同一 fold
```

### 5.3 评估指标矩阵

| 模型 | 回归指标 | 分类指标 | 选择性指标 |
|------|----------|----------|------------|
| JNK1 MTL | RMSE, R² | AUC (active) | — |
| Sel-Δ regression | RMSE on Δ12, Δ13 | — | Spearman ρ |
| Sel-classifier | — | AUC-ROC, F1 | Precision@K |
| Ensemble | 综合 | 综合 | **Benchmark rank** |

### 5.4 Benchmark 化合物回测（必须通过）

| 化合物 | 预期 JNK1 | 预期 JNK2/3 | 预期分类 |
|--------|-----------|-------------|----------|
| CC-90001 | 高 | 低 | JNK1-selective ✓ |
| SP600125 | 高 | 高 | pan-JNK ✓ |
| JNK-IN-8 | 高 | 高 (JNK2) | pan/covalent ✓ |
| BIRB-796 | 中 | 高 (JNK2) | JNK2-biased ✓ |

```text
Benchmark Pass Criteria:
  - CC-90001 预测 delta_min ≥ 1.0
  - SP600125 预测 delta_min < 0.5
  - Benchmark 排序 Spearman ρ ≥ 0.8
```

### 5.5 Applicability Domain (AD)

对百万库筛选后的 hits，评估 **是否在训练化学空间内** [Sheridan, 2013]：

```text
AD score = mean Tanimoto(hits, training_set)
Flag: AD score < 0.3 → 高外推风险，降低优先级
```

**对应脚本**: `scripts/04_train_selectivity_model.py`

---

## 6. Phase 4：百万分子库虚拟筛选

### 6.1 化合物库选择

| 库 | 规模 | 说明 |
|----|------|------|
| **Enamine REAL** | ~1–3.5 billion (pre-made subset ~1M–10M) | 可购买，推荐 [Enamine, 2023] |
| **ZINC22** |  billions | 免费下载 [Irwin et al., 2020] |
| **Mcule** |  millions | 按需下载 |

推荐：先用 **Enamine REAL Space 预计算构象子集 (~1M)** 做 pilot，再扩展。

### 6.2 筛选漏斗

```text
~1,000,000 compounds
    │
    ├─ [F1] 结构预处理
    │     RDKit sanitize, 去盐, 中性化
    │     PAINS/Brenk 过滤 [Baell & Holloway, 2010]
    │     → ~900,000
    │
    ├─ [F2] 类药过滤 (Lipinski/Ro5)
    │     MW 200–600, LogP −1–5, HBD ≤5, HBA ≤10
    │     → ~700,000
    │
    ├─ [F3] MTL 活性预测
    │     pred_pAct_JNK1 ≥ 7.0
    │     → ~50,000
    │
    ├─ [F4] 选择性预测
    │     pred_delta_min ≥ 1.0
    │     AND pred_pAct_JNK2 < 6.0
    │     AND pred_pAct_JNK3 < 6.0
    │     → ~5,000
    │
    ├─ [F5] 合成可行性
    │     SA score ≤ 4 [Ertl & Schuffenhauer, 2009]
    │     QED ≥ 0.4 [Bickerton et al., 2012]
    │     → ~2,000
    │
    ├─ [F6] ADMET 过滤 (in silico)
    │     hERG, CYP, 溶解度预测
    │     → ~500
    │
    ├─ [F7] Ensemble docking (可选)
    │     JNK1/2/3 三靶点对接
    │     ΔDock ≥ threshold
    │     → ~200
    │
    └─ [F8] 多样性选择 + AD 过滤
          Butina 聚类 (Tanimoto < 0.7), 每簇 Top 1
          Applicability domain check
          → ~50–100 最终候选
```

### 6.3 综合评分函数

```python
FinalScore = (
    0.35 * norm(pred_pAct_JNK1)
  + 0.30 * norm(pred_delta_min)
  + 0.15 * norm(dock_score_JNK1)
  + 0.10 * norm(QED)
  - 0.10 * norm(SA_score)
)
```

权重在 validation set 上用 **已知活性分子 rank recovery** 优化。

### 6.4 输出

```text
results/screening/
  ├── all_predictions.csv          # 全库预测 (可压缩)
  ├── top500.csv                   # 漏斗后 Top 500
  ├── top100_diverse.csv           # 多样性选择后 Top 100
  └── screening_report.md          # 筛选统计报告
```

**对应脚本**: `scripts/06_virtual_screening.py`

---

## 7. Phase 5：可解释性分析（SHAP 等）

> **原则**：可解释性不是可选项；JNK1 选择性模型的每一个关键预测都需有化学意义支撑。

### 7.1 SHAP 分析（树模型主选）

#### 7.1.1 适用模型

- **XGBoost / LightGBM / Random Forest** → SHAP TreeExplainer [Lundberg & Lee, 2017]
- 特征：**Morgan FP bits** 或 **MACCS keys**

#### 7.1.2 分析内容

| 分析 | 输出 | 意义 |
|------|------|------|
| **Global SHAP summary** | beeswarm plot | 哪些子结构 globally 促/抑 JNK1 选择性 |
| **SHAP dependence plot** | 特定 bit vs SHAP | 非线性 SAR 关系 |
| **Per-compound SHAP** | waterfall plot | 单个 hit 的选择性归因 |
| **Class comparison** | JNK1-sel vs pan-JNK SHAP 对比 | 选择性差异特征 |
| **Scaffold-level SHAP** | 按 Murcko 分组聚合 | 系列水平的选择性驱动 |

#### 7.1.3 指纹 bit → 子结构映射

```text
SHAP 高贡献 bit → RDKit bitInfo 提取对应 substructure
→ 与 MMP transform rules 交叉验证
→ 与 JNK1/2/3 口袋残基差异对照
```

### 7.2 GNN 可解释性（Chemprop 补充）

| 方法 | 说明 | 参考 |
|------|------|------|
| **Grad-CAM / Saliency** | 原子/键级别贡献 | [McCloskey et al., 2019] |
| **Integrated Gradients** | 路径积分归因 | [Sundararajan et al., 2017] |
| **Attention weights** | 若模型含 attention 层 | 直接读取 |

### 7.3 MMP-SHAP 交叉验证（推荐独特分析）

```text
For each top MMP transform rule (Phase 2):
  1. 找含该 transform 的化合物子集
  2. 计算该子集 SHAP 值是否显著高于 background
  3. 若一致 → "MMP-SHAP validated rule"
  4. 若不一致 → 需进一步调查 (assay noise / 多机制)
```

### 7.4 可解释性报告模板

```markdown
## SHAP Analysis Report

### Global drivers of JNK1 selectivity (Top 10 substructures)
1. [substructure SMILES] — mean |SHAP| = 0.42 — promotes JNK1 selectivity
2. ...

### Validated MMP-SHAP rules
- Rule: Ar-F → Ar-Cl | MMP Δ12 = +0.9 | SHAP bit 1234 | ✓ Consistent

### Benchmark compound attribution
- CC-90001: Top SHAP features = [hinge binder, back pocket substituent]
- SP600125: Top SHAP features = [pan-JNK pharmacophore elements]

### Structural mapping
- SHAP feature X maps to JNK1 Leu155 pocket (not present in JNK2 Ile)
```

**对应脚本**: `scripts/05_model_interpretation.py`

---

## 8. Phase 6：实验验证建议

### 8.1 最低验证集

| 实验 | 方法 | 通过标准 |
|------|------|----------|
| JNK1/2/3 酶活 IC50 | ADP-Glo (Promega V4070/V4080/V4090) | JNK1 IC50 < 1 µM |
| 选择性指数 | SI = IC50(JNK2/3) / IC50(JNK1) | SI > 10 |
| ATP 竞争性 | 不同 ATP 浓度 IC50 shift | 确认结合模式 |
| 细胞验证 | HTRF p-c-Jun (S63) | JNK1 功能抑制 |

### 8.2 推荐验证化合物

- Top 10 计算 hits（多样性选择后）
- 2–3 个 benchmark 回测分子（阳性对照）
- 1–2 个 SHAP 高度归因的 designed analog（可选）

---

## 9. 决策树与里程碑检查点

```
M0: 数据下载完成
    └─ Check: JNK1 ≥ 200, JNK2 ≥ 500, JNK3 ≥ 300 条 curated records
    └─ Check: paired_set ≥ 100 molecules

M1: 相似性分析完成
    └─ Check: 化学空间重叠报告
    └─ Decision: MTL vs single-task

M2: SAR/MMP 分析完成
    └─ Check: ≥ 5 条 validated MMP selectivity rules

M3: MTL 模型训练完成
    └─ Check: JNK1 scaffold-test R² ≥ 0.4
    └─ Check: Benchmark rank ρ ≥ 0.7

M4: 选择性模型完成
    └─ Check: Sel-classifier AUC ≥ 0.75
    └─ Check: CC-90001 被正确分类

M5: SHAP 分析完成
    └─ Check: ≥ 3 MMP-SHAP validated rules
    └─ Check: 结构映射完成

M6: 百万库筛选完成
    └─ Check: Top 100 diverse hits 输出
    └─ Check: AD filter applied

M7: 实验验证 (optional)
    └─ Check: ≥ 1 compound with SI > 10
```

---

## 10. 参考文献

完整格式化文献见 [REFERENCES.md](REFERENCES.md)。核心引用：

1. Zdrazil B, et al. The ChEMBL Database in 2023. *Nucleic Acids Res.* 2024;52(D1). [PMID: 37933841]
2. Bennett BL, et al. Discovery of CC-90001, a JNK1-Selective Inhibitor. *J Med Chem.* 2021;64(3). [PMID: 33404223]
3. Yang K, et al. Analyzing Learned Molecular Representations for Property Prediction. *J Chem Inf Model.* 2019;59(8). [PMID: 31361430]
4. Heid E, et al. Chemprop 2. *J Chem Inf Model.* 2024. [PMID: 38421620]
5. Lundberg SM, Lee SI. A Unified Approach to Interpreting Model Predictions. *NeurIPS.* 2017.
6. Wu Z, Rajpal DK. Machine Learning in Drug Discovery. *J Chem Inf Model.* 2020;60(12). [PMID: 33089927]
7. Bemis GW, Murcko MA. The Properties of Known Drugs. *J Med Chem.* 1996;39(15). [PMID: 8709122]
8. Dossetter AG, Griffen EJ. Matched Molecular Pair Analysis. *MedChemComm.* 2013;4. [DOI: 10.1039/C2MD20081B]
9. Davis MI, et al. Comprehensive Analysis of Kinase Inhibitor Selectivity. *Nat Biotechnol.* 2011;29(11). [PMID: 22037378]
10. Manning BD, Davis RJ. Targeting JNK for Therapeutic Benefit. *Nat Rev Drug Discov.* 2003;2(7). [PMID: 12838265]
11. Rogers D, Hahn M. Extended-Connectivity Fingerprints. *J Chem Inf Model.* 2010;50(5). [PMID: 20426451]
12. Sheridan RP. Applicability Domain. *J Chem Inf Model.* 2013;53(4). [PMID: 23560694]
13. Baell JB, Holloway GA. PAINS. *J Med Chem.* 2010;53(7). [PMID: 20131845]
14. Ertl P, Schuffenhauer A. SA Score. *J Cheminform.* 2009;1(1). [PMID: 20298526]
15. Bickerton GR, et al. Quantifying Druglikeness (QED). *Nat Chem.* 2012;4(2). [PMID: 22270643]
16. Friesner RA, et al. Glide. *J Med Chem.* 2004;47(7). [PMID: 15027865]
17. McNutt AT, et al. GNINA 1.0. *J Cheminform.* 2021;13(1). [PMID: 34719433]
18. Sundararajan M, et al. Axiomatic Attribution for Deep Networks. *ICML.* 2017.
19. Gilson MK, et al. BindingDB in 2015. *Nucleic Acids Res.* 2016;44(D1). [PMID: 26496856]
20. Irwin JJ, et al. ZINC20. *J Chem Inf Model.* 2020;60(12). [PMID: 33107064]
21. Landrum G. RDKit: Open-Source Cheminformatics. https://www.rdkit.org (2024).
22. Chen T, Guestrin C. XGBoost. *KDD.* 2016.
23. Preuer K, et al. Fréchet ChemNet Distance. *J Chem Inf Model.* 2018;58(9). [PMID: 30118593]
24. Czarnecki WM, et al. Activity Cliff-aware QSAR. *J Chem Inf Model.* 2018;58(6). [PMID: 29792674]
25. McCloskey K, et al. Using Attribution to Probe DL Models in Chemistry. *Chem Sci.* 2019;10. [PMID: 30842837]

---

## 附录 A：文件与脚本对应关系

| Phase | 脚本 | 输入 | 输出 |
|-------|------|------|------|
| 0 | `01_download_chembl_data.py` | ChEMBL API | `data/raw/jnk{1,2,3}_raw.csv` |
| 1 | `02_dataset_similarity.py` | raw/processed | `results/similarity/` |
| 2 | `03_sar_analysis.py` | processed + paired | `results/sar/` |
| 3 | `04_train_selectivity_model.py` | processed + paired | `models/` |
| 5 | `05_model_interpretation.py` | model + data | `results/shap/` |
| 4 | `06_virtual_screening.py` | model + library | `results/screening/` |

## 附录 B：预期时间线

| 阶段 | 预计时间 |
|------|----------|
| Phase 0–1 (数据 + 相似性) | 1–2 周 |
| Phase 2 (SAR + baseline ML) | 2–3 周 |
| Phase 3 (MTL + 选择性模型) | 2–3 周 |
| Phase 5 (SHAP) | 1 周 |
| Phase 4 (百万库筛选) | 1–2 天 (计算) |
| Phase 6 (实验) | 4–8 周 |
| **总计** | **约 3–4 个月** |

---

*本文档由 CADD/AIDD 专家流程设计，对应代码见 `scripts/` 目录。*
