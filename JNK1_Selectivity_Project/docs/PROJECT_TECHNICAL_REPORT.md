# JNK1/2/3 选择性抑制剂 CADD 项目技术报告

> **版本**: 1.0  
> **项目路径**: `JNK1_Selectivity_Project/`  
> **更新日期**: 2026-06

---

## 一、项目目标与总体策略

本项目针对 **JNK（c-Jun N-terminal kinase）三个亚型 JNK1、JNK2、JNK3**，基于 ChEMBL 公开生物活性数据，构建 **计算机辅助药物设计（CADD）/ AI 药物发现（AIDD）** 流程，核心目标为：

1. 清洗并比较 JNK1/2/3 三个数据集的结构与活性差异  
2. 建立 **单靶点活性预测模型**（QSAR），用于百万分子库 **粗筛**  
3. 用文献 benchmark 化合物校准筛选阈值  
4. 通过虚拟筛选漏斗获得候选分子，**亚型选择性** 留待结构对接（F3）验证  

**设计原则**：配对数据弱、isoform 选择性难以靠 ML 直接回归 → **活性用 ML 粗筛，选择性用对接**。

---

## 二、`results/` 文件夹内容说明

```
results/
├── calibration/              # F1 阈值校准（9 个文献 benchmark）
├── model_comparison/         # XGBoost vs Chemprop 对比
├── similarity/               # JNK1/2/3 数据集相似性分析
├── training/                 # 多任务/选择性模型训练曲线与报告
├── shap/                     # 可解释性分析（SHAP）
├── screening/                # 旧版筛选结果（JNK1≥7 + ML 选择性，已弃用）
├── screening_v2/             # 当前推荐：F1 p_family≥6.0 筛选（demo 库）
└── screening_v2_smoke/       # 小规模测试跑
```

### 各子目录详解

| 目录 | 主要文件 | 含义 |
|------|----------|------|
| **calibration/** | `threshold_scan.csv`, `benchmark_predictions.csv`, `threshold_recommendation.json` | 在 9 个文献抑制剂上扫描 F1 阈值；**推荐 p_family ≥ 6.0**（9/9 召回） |
| **model_comparison/** | `comparison.json`, `MODEL_COMPARISON_REPORT.md` | 三亚型 XGBoost 5-fold 骨架 CV + holdout 性能；**XGBoost 优于 Chemprop** |
| **similarity/** | `cross_similarity_matrix.csv`, `scaffold_overlap.json`, 多张 PNG | 三数据集化学空间、骨架重叠、活性分布比较 |
| **training/** | `training_report.json`, 多张 parity 图 | MTL 多任务模型 + 选择性 Δ 回归（**样本极少，仅供参考**） |
| **shap/** | `top_shap_features.csv`, beeswarm/waterfall 图 | 选择性分类模型的特征重要性（Bit_182、MolWt、NumHeteroatoms 等） |
| **screening_v2/** | `all_hits.csv`, `top500.csv`, `top100_diverse.csv`, `screening_report.json`, 漏斗图 | **当前筛选流程**在 demo 库（1835 分子）上的完整输出 |
| **screening/** | 旧 funnel 结果 | 历史版本，勿与 v2 混用 |

> **说明**：百万 Taosu 库筛选结果输出在 `results/screening_taosu_1M/`（本地运行后生成，默认不入 Git）。

---

## 三、训练数据：来源与清洗

### 3.1 数据来源

| 亚型 | ChEMBL ID | 清洗后化合物数 |
|------|-----------|----------------|
| JNK1 | CHEMBL2276 | **444** |
| JNK2 | CHEMBL4179 | **610** |
| JNK3 | CHEMBL2637 | **1147** |

原始导出文件：`docs/JNK1.csv`, `JNK2.csv`, `JNK3.csv`。

### 3.2 清洗规则（`scripts/utils_ml.py` → `curate_isoform_raw`）

为保证 QSAR 质量，仅保留：

- **测定类型**：IC50，`Standard Relation = '='`（精确值）  
- **Assay Type = B**（生化测定，非细胞）  
- **pActivity（pChEMBL）** 范围：**4.0–10.0**（10 nM – 0.1 M）  
- **冲突去除**：同一 SMILES 多条记录时，若 std > 0.5 log 或 range > 1.0 log，丢弃  
- **Assay harmonization**（按亚型调整）：
  - JNK1：assay 内 ≥ **50** 个化合物  
  - JNK2：≥ **8**（数据异质性大，放宽以保留样本）  
  - JNK3：≥ **20**  
- 多值取 **中位数** 作为该分子 pActivity  

### 3.3 数据集相似性（`results/similarity/`）

- 交叉 Tanimoto 中位相似度：JNK1–JNK2 ≈ **0.17**，JNK1–JNK3 ≈ **0.15**（化学空间部分重叠但可区分）  
- 三数据集共享 Murcko 骨架：**38** 个  
- JNK1/JNK2/JNK3 独特骨架：194 / 327 / 610  
- **结论**：三靶点数据相关但不冗余 → 适合 **分亚型单靶点模型**，不宜强行用一个 pan-JNK 模型替代  

### 3.4 数据划分

- **Murcko 骨架分组** 80/10/10 → train / val / test  
- 5-fold **骨架 GroupKFold** 做交叉验证  
- **目的**：避免随机划分导致化学相似分子同时出现在训练/测试集而造成 R² 虚高  

---

## 四、分子描述符（活性模型输入特征）

每个分子被编码为 **2060 维** 特征向量（`featurize_smiles`，`morgan_bits=2048`）：

### 4.1 Morgan 指纹（2048 维）

| 参数 | 值 |
|------|-----|
| 类型 | **Morgan 指纹（ECFP4 等价）** |
| 半径 | **2** |
| 位数 | **2048**（二进制 bit） |
| 工具 | RDKit `GetMorganFingerprintAsBitVect` |

作用：编码 **子结构/药效团** 信息，是激酶 QSAR 最常用的 2D 拓扑描述符。

### 4.2 RDKit 理化描述符（12 维）

| 序号 | 描述符 | 含义 |
|------|--------|------|
| 1 | MolWt | 分子量 |
| 2 | MolLogP | 脂水分配系数 |
| 3 | TPSA | 拓扑极性表面积 |
| 4 | NumHDonors | 氢键供体数 |
| 5 | NumHAcceptors | 氢键受体数 |
| 6 | NumRotatableBonds | 可旋转键数 |
| 7 | RingCount | 环数 |
| 8 | NumAromaticRings | 芳香环数 |
| 9 | FractionCSP3 | sp³ 碳比例 |
| 10 | BertzCT | 分子复杂度 |
| 11 | NumHeteroatoms | 杂原子数 |
| 12 | QED | 药物相似性定量估计 |

作用：补充 **整体理化性质**，改善极性、大小、复杂度与活性的关系；筛选阶段也单独用 QED/SA 过滤。

### 4.3 为何选「Morgan + 理化」而非仅用图神经网络？

- 数据量（JNK1 仅 444 条）对 **Chemprop/GNN** 偏少  
- Morgan+XGBoost 在激酶 benchmark 上稳健、可解释、训练快  
- 项目内对比：**XGBoost 平均 holdout R² = 0.68**，Chemprop 约 **0.53**（同数据）  

---

## 五、活性模型：怎么建、为什么这样建

### 5.1 模型架构

**三个独立的 XGBoost 回归器**，分别预测：

- `pAct_JNK1`、`pAct_JNK2`、`pAct_JNK3`  

**标签**：pActivity = −log₁₀(IC50 [M]) = pChEMBL Value  

**不用**单一多任务模型做筛选的原因：

- 配对分子仅 **322** 个（≥2 个亚型有数据），MTL 测试集极小  
- MTL 在 holdout 上 JNK2/JNK3 R² 仅 **~0.25**  
- 单靶点模型 JNK1/JNK3 holdout R² **> 0.69**  

### 5.2 XGBoost 超参数（`config/targets.yaml`）

| 参数 | 值 |
|------|-----|
| n_estimators | 2500 |
| max_depth | 7 |
| learning_rate | 0.015 |
| subsample | 0.85 |
| colsample_bytree | 0.65 |
| min_child_weight | 3 |
| reg_alpha | 0.8 |
| reg_lambda | 2.5 |
| early_stopping_rounds | 100 |

训练脚本：`scripts/07_compare_models.py`  
输出模型：`models/xgboost/xgboost_jnk{1,2,3}.joblib`（本地生成，不入 Git）

### 5.3 辅助模型（`results/training/`，非筛选用）

| 阶段 | 内容 | 用途 |
|------|------|------|
| Stage B MTL | Chemprop/多任务联合表 | 探索性；性能差 |
| Stage C Δ 回归 | 预测 JNK1 − max(JNK2,JNK3) | 配对样本 n=4（测试），**不可用于筛选** |
| Stage C 分类 | JNK1 选择性二分类 | 正样本训练仅 8 个，测试 0 个正例 |

---

## 六、模型性能

### 6.1 五折骨架交叉验证（XGBoost）

| 亚型 | Mean R² | Std | Mean Spearman |
|------|---------|-----|---------------|
| JNK1 | **0.662** | 0.086 | 0.772 |
| JNK2 | **0.443** | 0.074 | 0.678 |
| JNK3 | **0.633** | 0.089 | 0.793 |
| **平均** | **0.579** | — | **0.748** |

### 6.2 骨架 Holdout 测试集（80/10/10）

| 亚型 | n_test | R² | RMSE | Spearman |
|------|--------|-----|------|----------|
| JNK1 | 31 | **0.697** | 0.626 | **0.858** |
| JNK2 | 67 | **0.574** | 0.639 | **0.780** |
| JNK3 | 98 | **0.774** | 0.711 | **0.869** |
| **平均** | — | **0.682** | — | **0.836** |

### 6.3 性能解读

**优点**：

- JNK1、JNK3 达到 **R² ≈ 0.70–0.77**，Spearman **> 0.85** → 适合 **活性排序** 和 **粗筛**  
- RMSE ≈ **0.6 log**（约 4 倍浓度误差）对预筛可接受  

**局限**：

- JNK2 CV R² 仅 **0.44**（测定异质性、assay 来源杂）  
- Holdout 测试集较小（JNK1 仅 31 个）  
- 模型 **不能可靠预测 isoform 方向**（见 benchmark 验证）  

### 6.4 Benchmark 验证（9 个文献化合物）

| 分子 | 预期 profile | F1 通过 (≥6.0) | ML 亚型方向 |
|------|--------------|----------------|-------------|
| SP600125 | pan-JNK | ✅ | 大致合理 |
| CC-930 | JNK2/3 偏向 | ✅ | ✅ |
| E1 | JNK1 偏向 | ✅ | ❌（预测 JNK2 最高） |
| TCS JNK 6O | JNK1 偏向 | ✅ | ❌ |

→ **F1 用于“有没有 JNK 活性”可靠；用于“JNK1 选择性”不可靠**。

---

## 七、为什么这样筛选（虚拟筛选漏斗 v2）

### 7.1 漏斗设计

```
输入 SMILES
    ↓ F0 预处理（RDKit 标准化、去重）
    ↓ F2 Lipinski 类药（MW 200–600, logP −1~5, HBD≤5, HBA≤10）
    ↓ F1 ML 预筛：p_family = max(pred_JNK1, pred_JNK2, pred_JNK3) ≥ 6.0
    ↓ F5 成药性：SA ≤ 6.0, QED ≥ 0.35
    ↓ 综合评分排序 → Top N → Butina 多样性挑选
```

### 7.2 为何用 `p_family = max(JNK1,JNK2,JNK3)`？

- 百万库第一步只需去掉 **明显无 JNK 活性** 的分子  
- pan-JNK 工具化合物（SP600125、CC-90001）也应保留  
- 9 个 benchmark 在阈值 **6.0** 时 **100% 召回**（`results/calibration/threshold_recommendation.json`）  

### 7.3 为何不用 ML 选择性过滤？

旧版 `screening/` 使用 JNK1 ≥ 7.0 + ΔpActivity 选择性，已被证明：

- E1、TCS JNK 6O 的 **亚型方向预测错误**  
- ChEMBL 配对选择性标签仅 8 个正例，分类模型 F1=0  

### 7.4 综合评分公式

```
final_score = 0.55 × (p_family/10)
            + 0.15 × (pred_JNK1/10)
            + 0.20 × QED
            + 0.10 × (10 − SA)/10
```

略向 JNK1 倾斜，但 **不硬筛选择性**。

### 7.5 一次筛选、三个模型都算

每个分子 **只遍历一遍** 库，同时得到 3 个预测值，取 max 过 F1 — **不是三个模型各筛一遍**。

---

## 八、已得到的筛选结果

### 8.1 Demo 库（`results/screening_v2/`，1835 个 ChEMBL 分子）

| 阶段 | 数量 | 通过率 |
|------|------|--------|
| 输入 | 1,835 | — |
| 预处理 | 1,835 | 100% |
| Lipinski 类药 | 1,541 | 84% |
| **F1 p_family ≥ 6.0** | **1,292** | 84% |
| SA/QED 通过 | **1,211** | 94% |

输出：

- `all_hits.csv`：**1,211** 个 hit  
- `top500.csv` / `top100_diverse.csv`  
- Top 分子多为 **氯代嘧啶-苯甲酰胺** 类 ChEMBL 已知 JNK 抑制剂骨架  

> Demo 库来自 JNK 活性数据本身，通过率高 **不代表** 百万随机库也有同样比例；商业库预期 F1 通过率约 **5–15%**。

### 8.2 百万 Taosu 库（本地 WSL 运行）

```bash
cd JNK1_Selectivity_Project
source .venv/bin/activate
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop
python3 scripts/06_virtual_screening.py \
  --library "/path/to/taosu_20210823_100w_asteroid_murcko_protonized.csv" \
  --output results/screening_taosu_1M \
  --batch-size 50000 \
  --top-n 5000 \
  --diverse-n 500
```

- **没有 `--max-rows`** → 读取完整 CSV（约 100 万行）  
- 预期耗时：**1–3 小时**（视 CPU）  
- 预期输出：`all_hits.csv` + `top5000.csv` + `top500_diverse.csv`  

---

## 九、SHAP 可解释性（`results/shap/`）

基于 **322 个配对分子** 的选择性分类模型，Top 特征：

| 排名 | 特征 | 方向 |
|------|------|------|
| 1 | Bit_182（Morgan 子结构） | 抑制选择性标签 |
| 2 | NumHeteroatoms | 促进 |
| 3 | MolWt | 促进 |
| 4–5 | Bit_656, Bit_231 | 抑制 |

用途：理解 **哪些结构模式与“标注的选择性”相关**，但因正样本极少，结论需谨慎，不能替代对接验证。

---

## 十、方法学总结（论文 Methods 要点）

| 环节 | 方法 | 理由 |
|------|------|------|
| 数据 | ChEMBL 生化 IC50，按亚型清洗 | 公开可重复 |
| 特征 | Morgan-2/2048 + 12 RDKit 2D 描述符 | 激酶 QSAR 标准组合 |
| 模型 | 三独立 XGBoost 回归（pActivity） | 数据量适配、性能优于 GNN |
| 验证 | Murcko 骨架 CV + holdout | 避免相似性泄漏 |
| 阈值 | 9 个文献抑制剂校准 F1=6.0 | 数据驱动、可解释 |
| 筛选 | p_family + Lipinski + SA/QED | 粗筛活性与成药性 |
| 选择性 | **不用 ML**；规划 F3 对接 | 配对数据与 benchmark 证明 ML 方向不可靠 |

---

## 十一、局限与下一步

| 局限 | 下一步 |
|------|--------|
| JNK2 模型 CV 偏弱 | 收紧 assay 或引入外部数据 |
| ML 无法判 isoform 方向 | JNK1/2/3 三结构对接 + ΔIFP |
| 百万库结果未做对接 | Top 5000 → 对接 → Top 100 |
| 无实验验证 | 同批次 JNK1/2/3 IC50 |

---

## 十二、相关脚本与配置文件

| 文件 | 功能 |
|------|------|
| `scripts/00_prepare_user_data.py` | ChEMBL CSV → 清洗数据 |
| `scripts/02_dataset_similarity.py` | 三靶点相似性分析 |
| `scripts/07_compare_models.py` | 训练三亚型 XGBoost |
| `scripts/calibrate_threshold.py` | Benchmark F1 阈值校准 |
| `scripts/06_virtual_screening.py` | 虚拟筛选漏斗 v2 |
| `config/targets.yaml` | 靶点、清洗、训练、筛选参数 |
| `data/benchmarks/literature_benchmarks.csv` | 9 个文献参考抑制剂 |

---

## 十三、一句话结论

本项目用 **2048-bit Morgan + 12 个 RDKit 描述符** 训练 **三个 XGBoost 单靶点模型**（JNK1/JNK3 holdout R² ≈ 0.70/0.77），以 **p_family ≥ 6.0** 对百万库做 **JNK 家族活性粗筛**，benchmark 召回 100%；**亚型选择性不在此步判断**。`results/screening_v2` 与 `results/calibration` 为当前有效结果；`results/screening` 为已弃用旧版。

---

## 参考文献

详见 [REFERENCES.md](REFERENCES.md) 与 [JNK1_selectivity_screening_workflow.md](JNK1_selectivity_screening_workflow.md)。
