# JNK1/2/3 亚型抑制剂计算筛选项目报告

> **版本**: 2.0  
> **日期**: 2026-07-03  
> **原则**: 本报告所有数值均来自仓库内可复现文件、对接工作区归档 CSV 或 MD QC 结果；未在数据中出现的结论一律不作断言。

---

## 摘要

本项目以 **JNK（c-Jun N-terminal kinase）三个亚型 JNK1、JNK2、JNK3** 为对象，构建 **ML 活性粗筛 → Glide XP 5000 化合物 VSW → MM-GBSA/ADMET 短名单 → Desmond MD pose QC** 四级计算漏斗，最终推荐 **10 个分子** 进入同批次 JNK1/2/3 酶学 IC50 湿实验。

**端到端漏斗（有数据支撑的各阶段）**：

| 阶段 | 数量 | 数据来源 |
|------|------|----------|
| ML 初筛后对接库（F0） | **4983** | `md_shortlist_report_23c8.md` |
| Glide XP VSW 有效记录 | **4979** | `JNK1_SELECTIVITY_FINAL_REPORT_41d9.md` |
| VSW pass_selectivity（Δsel_dock + MM-GBSA） | **233** | 同上 |
| VSW Tier 1 | **57** | 同上 |
| MD shortlist（ADMET 后） | **25** | `md_shortlist_report_23c8.md` |
| MD pose QC 输入 | **16** | `MD_QC_report_cf26.md` |
| 最终采购推荐 | **10** | `data/purchase/purchase_after_md.csv` |

**核心结论**：

1. **受体准备可信**：5/5 共晶再对接 RMSD < 2 Å（`redocking_summary_7725.csv`）。
2. **Glide XP 对接选择性方向不可靠**：Spearman(Δsel_dock, −ΔpIC50_sel) = **0.786**（p = 0.021），但 isoform 方向准确率仅 **29%**（2/7）；TCS JNK 6O、SP600125 方向预测失败（`direction_confusion_27c3.csv`）。
3. **ML 同样不能预测亚型方向**（E1、TCS JNK 6O 预测错误）；ChEMBL **JNK1-selective 标注仅 8 个**，选择性分类器测试 F1 = **0**。
4. **Gly87（KLIFS b.l.37）占据策略回顾性失败**：5/5 benchmark `occ_JNK1=True`，但配体距 Gly87 **0.59–1.18 Å**，无选择性判别力（`gly87_selfcheck_16be.csv`）。
5. MD QC：**G1 3/4、G2 0/6** 通过 `pass_md_overall`；**没有任何分子的 JNK1 选择性可被计算确认**。
6. 采购 10 个分子的理由：**G3 实验校准 + G1/G2 chemotype 假说检验 + MD pose 可信度分层**，而非“已算出选择性 hit”。

---

## 1. 项目背景与目标

### 1.1 生物学与成药背景

| 亚型 | 基因 | ChEMBL Target ID |
|------|------|------------------|
| JNK1 | MAPK8 | CHEMBL2276 |
| JNK2 | MAPK9 | CHEMBL4179 |
| JNK3 | MAPK10 | CHEMBL2637 |

JNK1 在 IPF、NASH 等疾病中有明确证据；**CC-90001** 为 JNK1 功能偏向临床候选（Bennett et al., *J. Med. Chem.* 2021, PMID: 33404223）。**SP600125** 为经典 pan-JNK 工具药（Bennett et al., *PNAS* 2001, PMID: 11717429）。**E1** 为 Pan et al. 报道的 JNK1 强抑制剂（IC50 = 2.7 nM，*J. Med. Chem.* 2024, doi:10.1021/acs.jmedchem.4c01764）。

### 1.2 策略演进

| 阶段 | 目标 | 数据结论 |
|------|------|----------|
| 初期 | 计算筛选 JNK1 选择性 hit | ML + 对接 benchmark 均否定 |
| 中期 | Glide Δsel + MM-GBSA 选择性分级 | 233 个严格通过，但方向准确率 29% |
| 后期 | MD pose QC + 湿实验 | 以 **pan-JNK 家族结合剂** 为主假说，选择性靠酶学 |

---

## 2. 训练数据与 ML 模型

### 2.1 数据清洗结果

| 亚型 | 化合物数 | Holdout R² | Holdout Spearman ρ |
|------|----------|------------|-------------------|
| JNK1 | 444 | **0.697** | **0.858** |
| JNK2 | 610 | **0.574** | **0.780** |
| JNK3 | 1147 | **0.774** | **0.869** |

数据来源：`data/processed/data_summary.json`，`results/model_comparison/comparison.json`。

### 2.2 选择性标签稀缺性

- 配对分子（≥2 亚型）：**322** 个
- JNK1-selective 标注：**8** 个（`sel_class_counts.csv`）
- 选择性分类器：训练正例 8，测试正例 0，**F1 = 0**（`training_report.json`）

### 2.3 ML 虚拟筛选（F1）

9 个文献 benchmark 在 **p_family ≥ 6.0** 时 **9/9 全部通过**（`threshold_recommendation.json`）。

Demo 库（1835 分子）漏斗：`screening_v2/screening_report.json`

| 阶段 | 数量 |
|------|------|
| Lipinski 通过 | 1541 |
| F1 通过 | 1292 |
| SA/QED 通过 | 1211 |

**ML 用途**：去除无 JNK 家族活性潜力的分子；**不用于 isoform 方向判断**。

### 2.4 ML vs 对接：benchmark 方向对比

| 化合物 | 实验 profile | ML 预测最高亚型 | 对接 Δsel 预测方向 | 实验方向（IC50） |
|--------|--------------|-----------------|-------------------|------------------|
| E1 | JNK1 偏好 | **JNK2**（7.56） | **JNK1**（Δsel +3.05） | JNK1 |
| TCS JNK 6O | JNK1 偏好 | **JNK3**（6.97） | **JNK23**（Δsel −1.18） | JNK1 |
| CC-930 | JNK2/3 偏好 | JNK2（7.47） | JNK23（Δsel −4.90） | JNK23 |
| SP600125 | pan-JNK | JNK1（6.13） | JNK23（Δsel −2.57） | pan |

→ ML 与对接在关键对照上**均不能一致预测 isoform 方向**；对接对 E1 方向正确、ML 错误；对 TCS JNK 6O 两者均错误。

---

## 3. Glide XP 结构对接与 5000 化合物 VSW

### 3.1 结构 Ensemble（`config/docking_ensemble.yaml`）

| Isoform | PDB | 聚合 | VSW 主结构 |
|---------|-----|------|------------|
| JNK1 | 3ELJ, 4L7F | 均值 | 3ELJ |
| JNK2 | 3E7O | 单结构 | 3E7O |
| JNK3 | 3TTI, 4WHZ | 均值 | 3TTI |

对接得分：**Glide XP `r_i_glide_gscore`**（非 SP 中间分）。

### 3.2 共晶再对接验证

**5/5 通过**（RMSD 阈值 2.0 Å，`results/docking_validation/redocking_summary_7725.csv`）：

| PDB | 靶标 | 配体 | Glide XP 得分 | RMSD (Å) |
|-----|------|------|---------------|----------|
| 3ELJ | JNK1 | GS7 | −12.79 | **0.66** |
| 4L7F | JNK1 | AX13587 | −12.21 | **0.92** |
| 3E7O | JNK2 | 35F | −11.27 | **0.26** |
| 3TTI | JNK3 | CC-930 | −12.86 | **1.50** |
| 4WHZ | JNK3 | 3NL | −10.09 | **1.88** |

**结论**：prepared 受体网格可信，可用于 pose 比较与 VSW。

### 3.3 选择性指标定义

```
Δsel_dock = min(score_JNK2, score_JNK3) − score_JNK1
```

- 得分越负 = 结合越强  
- **Δsel_dock > 0** → 计算上偏向 JNK1

严格选择性门槛：`pass_selectivity` = Δsel_dock > 0 **且** Δsel_MMGBSA ≥ 2.0 kcal/mol。  
**注意**：9 个 benchmark **未运行 Prime MM-GBSA**，MM-GBSA 选择性门槛**未经 benchmark 标定**（`validation_report.md` §3.4）。

### 3.4 5000 化合物 VSW 漏斗

数据来源：`JNK1_SELECTIVITY_FINAL_REPORT_41d9.md`，`candidates_ranked_befe.csv`（4979 条）

| 阶段 | 数量 |
|------|------|
| 输入库 | **4979** |
| pass_pose | 3234 |
| pass_potency（score_JNK1 ≤ 中位数 **−6.65** kcal/mol） | 1681 |
| **pass_selectivity**（Δsel + MM-GBSA 双通过） | **233** |
| has_selectivity_contact（IFP 代理） | 63 |

**Tier 分布**：

| Tier | 数量 | 含义 |
|------|------|------|
| **Tier 1** | **57** | pose + potency + selectivity + consistency + contact 代理 |
| Tier 2 | 92 | 能量选择性通过，机制/一致性待确认 |
| Tier 3 | 1191 | JNK1 有活性，选择性在噪声区 |
| Tier 0 | 3639 | 未达 Tier3 |

**泛 JNK + 计算 JNK1 偏好子集**（`panJNK_JNK1bias_ba7c.csv`）：

| 子集 | 数量 |
|------|------|
| pass_potency ∧ Δsel_dock > 0 | **679** |
| 上述 + 三 isoform score 均 ≤ −6 | **431**（更可能 pan-JNK） |
| pass_selectivity（严格双门槛） | **233** |

> “计算 JNK1 偏好”仅表示 Δsel_dock > 0，**不等于实验 JNK1 选择性**。

### 3.5 Top 选择性候选（pass_selectivity，Δsel 降序前 5）

| compound_id | Δsel_dock | score_JNK1 | Tier |
|-------------|-----------|------------|------|
| 4931 | 3.67 | −7.97 | 1 |
| 2627 | 3.65 | −6.83 | 1 |
| 1941 | 3.37 | −7.01 | 1 |
| 2749 | 3.09 | −10.21 | 1 |
| 2760 | 2.68 | −10.19 | 1 |

`top_selective_f4a0.csv`：50 个 Butina 聚类代表（Tanimoto 0.5），其中 Tier1 = **15**。

---

## 4. 文献 Benchmark 对接选择性验证

### 4.1 定量结果（9 化合物，`benchmark_deltas_51c1.csv`）

| 指标 | 数值 | 阈值 | 达标 |
|------|------|------|------|
| Spearman(Δsel_dock, −ΔpIC50_sel) | **0.786** (p = 0.021, n = 8) | \|ρ\| ≥ 0.35 | ✓ |
| 方向准确率（全部有 IC50） | **22%** (2/9) | ≥ 55% | ✗ |
| 方向准确率（JNK1/JNK23/PAN 标签） | **29%** (2/7) | ≥ 55% | ✗ |

**Spearman 与方向准确率的“分裂”**：连续变量秩相关尚可，但离散 isoform 方向标签预测失败率高 → **Glide Δscore 不宜作为 isoform 分型决策依据**。

### 4.2 四个关键对照（`direction_confusion_27c3.csv`）

| 化合物 | 实验 profile | Δsel_dock | ΔpIC50_sel | 实验方向 | 预测方向 | 匹配 |
|--------|--------------|-----------|------------|----------|----------|------|
| **E1** | JNK1 偏好 | **+3.05** | −0.85 | JNK1 | JNK1 | **✓** |
| **CC-930** | JNK2/3 偏好 | **−4.90** | +0.94 | JNK23 | JNK23 | **✓** |
| TCS JNK 6O | JNK1 偏好 | −1.18 | −0.55 | JNK1 | JNK23 | ✗ |
| SP600125 | pan-JNK | −2.57 | −0.35 | pan | JNK23 | ✗ |

### 4.3 各 isoform 活性排序（对接分 vs pIC50）

| Isoform | Spearman ρ | n |
|---------|------------|---|
| JNK1 | **−0.429** | 8 |
| JNK2 | +0.190 | 8 |
| JNK3 | +0.371 | 6 |

JNK1 排序呈**负相关**，说明跨 PDB 绝对得分比较存在**系统偏差**（`isoform_rank_correlations_299a.csv`）。

### 4.4 验证结论

**对接选择性方向不能可靠用于 isoform 分型。** 建议将 VSW 命中视为 **JNK 家族活性/构象假设**，选择性须由 **同批次 JNK1/2/3 IC50** 确认。

---

## 5. 选择性策略尝试与失败记录

### 5.1 尝试一览

| # | 策略 | 结果 | 数据依据 |
|---|------|------|----------|
| 1 | ML 三模型 ΔpActivity 选择性过滤 | **失败** | E1/TCS JNK 6O 方向错误 |
| 2 | 选择性二分类模型 | **失败** | 正例 n=8，F1=0 |
| 3 | Glide Δsel_dock 排序 | **不可靠** | 方向准确率 29% |
| 4 | Δsel_dock + MM-GBSA ≥ 2 kcal/mol | **未标定** | Benchmark 无 MM-GBSA |
| 5 | KLIFS 非保守位点 / Gly87 IFP | **放弃** | 见 §5.2 |
| 6 | MD pose QC（RMSD + hinge HB） | **部分可用** | 验证 pose，非选择性 |
| 7 | Tier 1 + FEP+ 推荐 | **待做** | 15 个 Tier1 候选含 690 |

### 5.2 Gly87（b.l.37）占据策略 — 回顾性自检失败

`results/docking_validation/gly87_selfcheck_16be.csv`（5 个 benchmark）：

| 配体 | d_Gly87 (Å) | occ_JNK1 | 预测 JNK1 选择性 | 实验 profile | 匹配 |
|------|-------------|----------|------------------|--------------|------|
| E1 | 0.744 | True | False | JNK1 偏好（7.0×） | False |
| TCS JNK 6O | 0.899 | True | False | JNK1 偏好（3.6×） | False |
| CC-930 | 1.180 | True | False | JNK2/3 偏好 | True* |
| SP600125 | 1.009 | True | False | pan | True* |
| CC-90001 | 0.590 | True | False | 近 pan（2.8×） | True* |

\* CC-930/SP600125/CC-90001 的 `match=True` 来自反向/近 pan 标签，**不是** JNK1 选择性验证。

**失败原因**：所有 benchmark 的配体距 JNK1 Gly87（JNK2 对应 Ser87）仅 **0.59–1.18 Å**，铰链区高度保守，**b.l.37 位点无法区分 isoform**。该策略在 MD shortlist 阶段被明确排除（`md_shortlist_report_23c8.md`）。

### 5.3 方法局限（来自 `JNK1_SELECTIVITY_FINAL_REPORT_41d9.md`）

1. ATP 口袋高度保守，Glide 得分差 1–3 kcal/mol 常处噪声水平  
2. 跨 PDB 蛋白准备差异引入 spurious Δsel  
3. MM-GBSA 选择性门槛未 benchmark 标定  
4. `has_selectivity_contact` 仅为铰链 H-bond + Δsel 启发式，非完整 IFP  
5. 5000 库未在备用结构 4L7F/4WHZ 上重跑 VSW；`pass_consistency` 为占位  
6. JNK-IN-8（共价）等不符合简单 Δsel 逻辑  

---

## 6. MD 短名单与 Pose QC

### 6.1 MD 短名单漏斗（`md_shortlist_report_23c8.md`）

**声明**：未使用 Δsel_dock/MMGBSA 选择性方向、Gly87 IFP 作为硬筛；短名单为 **JNK 家族结合剂** pose QC，非最终采购单。

| 阶段 | 数量 |
|------|------|
| 输入（F0 后） | **4983** |
| F1 pose QC 通过 | 3125 |
| F2 活性 + 配体效率通过 | 182 |
| **F1 ∧ F2 通过** | **157** |
| F7 QikProp ADMET 剔除 | 9 |
| ADMET backfill | 9 |
| **ADMET 后 shortlist** | **25** |

**门槛**：
- score_JNK1 ≤ **−7.43**（活性 benchmark 中位数）
- MMGBSA_JNK1 ≤ **−51.60**

**F2（基础成药性）**：MW、cLogP、HBD/HBA、QED、SA、PAINS  
**F7（QikProp @3ELJ）**：hERG、口服吸收、Caco-2、溶解度、#stars ≤ 0  
**G3 对照**：ADMET 豁免保留

**分组（25 个 shortlist）**：

| 组 | 数量 | 进入 MD（16） |
|----|------|---------------|
| G1 骨架模仿 | 9 | **4**（690, 2232, 2157, 2389） |
| G2 新骨架 | 10 | **6** |
| G3 对照 | 4 | **4**（全部） |
| G4 阴性锚点 | 2 | **2**（全部） |

chemotype_sim（ECFP4 Tanimoto vs E1/Q63/TCS JNK 6O）：
- G1 mean = **0.217**（median 0.213）
- G2 mean = **0.120**（median 0.114）

### 6.2 MD QC 方法与结果（`MD_QC_report_cf26.md`）

- 48 个 Desmond 任务（16 × 3 PDB：3ELJ / 3E7O / 3TTI）
- pass_md_JNK1：RMSD ≤ 3 Å + hinge HB ≥ 30%
- pass_md_overall：JNK1 pass + (JNK2 **或** JNK3 pass)

| 阶段 | 数量 |
|------|------|
| MD 输入 | 16 |
| pass_md_JNK1 | 6 |
| pass_md_overall | 5 |
| 采购推荐 | 10 |

| 组 | n | pass_overall | pose grade |
|----|---|--------------|------------|
| G1 | 4 | **3/4** | A×3, F×1 |
| G2 | 6 | **0/6** | C×1, F×5 |
| G3 | 4 | 2/4 | — |
| G4 | 2 | **0/2** | F×2（阴性验证 ✓） |

### 6.3 采购分子对接背景（VSW 数据）

| ID | 组 | Tier | Δsel_dock | score_JNK1 | pass_selectivity | MD pass_overall |
|----|-----|------|-----------|------------|------------------|-----------------|
| **690** | G1 | **1** | **+1.08** | −7.76 | **Yes** | **Yes**（三 isoform 全 pass） |
| 2232 | G1 | 1 | +1.32 | −8.13 | Yes | Yes |
| 2157 | G1 | 3 | −1.05 | −8.46 | No | Yes |
| 2231 | G2 | 3 | +3.37 | −11.22 | No | No（JNK1-only） |
| 4795 | G2 | 1 | +1.57 | −8.38 | Yes | No |
| 1280 | G2 | 3 | +0.89 | −7.85 | No | No |

690 同时出现在：**Tier 1**、**top_selective 聚类代表**、**FEP+ 推荐 15 清单**、**panJNK_JNK1bias 子集**（`candidates_ranked_befe.csv`）。

> 注意：2157 的 Δsel_dock 为负（计算偏向 JNK2/3），但 MD 仍 pass overall——再次说明 **Δsel 与 MD pose 可不一致**，不宜混用。

---

## 7. 采购清单与花钱理由

完整表格：`data/purchase/purchase_after_md.csv`（10 分子，SMILES 经 RDKit 验证）

### 7.1 采购结构

| 类别 | n | 化合物 | 理由 |
|------|---|--------|------|
| G3 对照 | 4 | SP600125, CC-90001, CC-930, E1 | **酶学校准尺**（无论 MD 是否通过） |
| G1 主力 | 3 | 690, 2232, 2157 | MD pass_overall + 对接 Tier1/活性 |
| G2 探索 | 3 | 2231, 1280, 4795 | G2 最优 + off-target pose 假说 |

### 7.2 花钱的逻辑链（可向合作者说明）

```
4979 化合物 VSW
  → 233 严格选择性通过（但方向 benchmark 仅 29% 准确）
  → 157 活性+pose 通过
  → 25 ADMET shortlist
  → 16 MD QC
  → 10 采购（含 4 个已知活性对照）
```

**花钱买的不是“选择性 hit”**，而是：

1. **验证计算管线**：G3 对照建立 IC50 vs MD 的相关性  
2. **检验 chemotype 假说**：G1（Tc~0.22）是否比 G2（Tc~0.12）更易出 JNK 活性  
3. **捕捉最有信息量的候选**：690（Tier1 + MD 三 isoform pass + FEP+ 推荐）  
4. **探索性 backup**：2231（G2 中 JNK1 MD 最好）、1280/4795（JNK2/3 pose 稳、JNK1 不稳）

---

## 8. 当前选择性状况（诚实评估）

### 8.1 计算层面

| 问题 | 答案 |
|------|------|
| 能否计算确认 JNK1 选择性？ | **不能** |
| 最强计算证据是什么？ | 690：Tier1 + Δsel>0 + MD 三 isoform pass → 更支持 **pan-JNK 结合** |
| 对接“233 个选择性通过”有意义吗？ | 仅作 **家族内优先级**，不能作 isoform 标签 |

### 8.2 各分子选择性先验（待实验，非结论）

| 分子 | 先验假说 | 依据 |
|------|----------|------|
| 690 | pan-JNK 或弱 JNK1 偏好 | Tier1, Δsel+1.08, MD 三 isoform pass |
| 2232 | 可能 pan-JNK | MD JNK1/JNK2 极好 |
| 2157 | 未知 | Δsel 为负但 MD pass |
| 2231 | 未知 | 仅 JNK1 MD A 级 |
| 1280/4795 | 可能 JNK2/3 ≥ JNK1 | JNK1 MD fail |

---

## 9. 湿实验预测

### 9.1 必做

**同批次 JNK1 + JNK2 + JNK3 重组酶 IC50**（10 分子 + DMSO 空白）

### 9.2 保守预测

| 场景 | 可能性 | 若成立的意义 |
|------|--------|--------------|
| G3 有活性但 MD-fail（CC-930, SP600125） | **较可能** | hinge HB 门槛偏严；MD 不能代替活性 |
| G1 ≥1 个 IC50 < 1 µM | **中等** | 三级漏斗后合理命中率 |
| G1 活性 > G2 | **不确定** | G2 MD 0/6 overall pass |
| ≥10× JNK1 选择性 | **低** | 所有计算选择性方法均失败 |
| 690 为 pan-JNK | **需实验区分** | 与 MD/对接数据一致 |
| G4 无活性 | **预期** | 阴性锚点验证 |

### 9.3 后续（若有 hit）

kinome 面板 + 细胞 p-c-Jun；对 top 1–2 考虑 **FEP+**（690 已在推荐清单）。

---

## 10. 数据文件索引

### 10.1 仓库内（Git）

| 路径 | 内容 |
|------|------|
| `docs/JNK1_PROJECT_REPORT.md` | **本报告** |
| `data/benchmarks/literature_benchmarks.csv` | 9 个文献 benchmark |
| `results/calibration/` | ML F1 阈值校准 |
| `results/screening_v2/` | ML 虚拟筛选 demo |
| `results/model_comparison/` | XGBoost 性能 |
| `results/docking_validation/` | 再对接、benchmark Δ、Gly87 自检 |
| `config/docking_ensemble.yaml` | Ensemble 与门槛配置 |
| `data/purchase/purchase_after_md.csv` | 10 分子采购表 |

### 10.2 对接工作区归档（本地 + 已摘要入本报告）

| 文件 | 内容 |
|------|------|
| `JNK1_SELECTIVITY_FINAL_REPORT_41d9.md` | VSW 综合报告 |
| `candidates_ranked_befe.csv` | 4979 化合物全量排名 |
| `top_selective_f4a0.csv` | 50 聚类代表 |
| `panJNK_JNK1bias_ba7c.csv` | 679 pan-JNK + 计算 JNK1 偏好 |
| `md_shortlist_report_23c8.md` | MD 短名单漏斗 |
| `MD_QC_report_cf26.md` | MD pose QC |
| `md_pose_qc_summary_5ffb.csv` | 16 化合物 MD 汇总 |

---

## 11. 参考文献

1. Zdrazil B, et al. ChEMBL 2023. *Nucleic Acids Res.* 2024;52(D1):D1180-D1192. doi:10.1093/nar/gkad1004  
2. Bennett BL, et al. CC-90001. *J. Med. Chem.* 2021;64(3):1776-1795. doi:10.1021/acs.jmedchem.0c01843  
3. Bennett BL, et al. SP600125. *PNAS* 2001;98(24):13681-13686. doi:10.1073/pnas.251194298  
4. Pan X, et al. JNK1 inhibitors for IPF (compound E1). *J. Med. Chem.* 2024. doi:10.1021/acs.jmedchem.4c01764  
5. Szczepankiewicz BG, et al. TCS JNK 6o. *J. Med. Chem.* 2006;49(14):3563-3566. doi:10.1021/jm060150w  
6. Plantevin-Krenitsky V, et al. 3TTI/CC-930 co-crystal. *Bioorg. Med. Chem. Lett.* 2012;22(3):1433-1438  
7. Friesner RA, et al. Glide. *J. Med. Chem.* 2004;47(7):1739-1749. doi:10.1021/jm0306430  
8. Manning BD, Davis RJ. Targeting JNK. *Nat. Rev. Drug Discov.* 2003;2(7):554-565. doi:10.1038/nrd1132  
9. Chen T, Guestrin C. XGBoost. *Proc. 22nd ACM SIGKDD* 2016. doi:10.1145/2939672.2939785  

---

## 附录 A：端到端筛选流程图

```mermaid
flowchart TD
    A[ChEMBL 444/610/1147] --> B[XGBoost 三靶点模型]
    B --> C[ML F1: p_family >= 6.0]
    C --> D[4983 化合物 F0]
    D --> E[Glide XP VSW 4979]
    E --> F{pass_selectivity?}
    F -->|233| G[Tier1=57]
    F -->|679 panJNK+bias| H[计算JNK1偏好 未验证]
    G --> I[F1+F2: 157]
    I --> J[QikProp ADMET: 25]
    J --> K[MD QC: 16]
    K --> L[采购: 10]
    L --> M[湿实验 JNK1/2/3 IC50]
    
    B -.->|方向错误| X1[ML选择性 放弃]
    E -.->|方向29%| X2[Δsel 不作硬筛]
    J -.->|Gly87失败| X3[IFP策略 放弃]
```

## 附录 B：一句话答辩版

> 我们从 ChEMBL 训练 JNK1/2/3 活性模型（holdout R² 0.70/0.57/0.77），经 ML 粗筛后对 **4979** 个化合物做 Glide XP 对接，再对接 benchmark 证明 **isoform 方向准确率仅 29%**；因此将 hit 定位为 **pan-JNK 家族结合剂**，经 ADMET 缩至 25 个、MD 验证 16 个，最终采购 10 个（含 4 个文献对照）做同批次三 isoform IC50——**选择性只能由实验回答，不能由计算采购**。
