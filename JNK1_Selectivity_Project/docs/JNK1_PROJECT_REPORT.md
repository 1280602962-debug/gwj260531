# JNK1/2/3 亚型抑制剂计算筛选项目报告

> **版本**: 1.0  
> **日期**: 2026-07-03  
> **原则**: 本报告所有数值均来自仓库内可复现文件或已归档的 MD QC 结果；未在数据中出现的结论一律不作断言。

---

## 摘要

本项目以 **JNK（c-Jun N-terminal kinase）三个亚型 JNK1、JNK2、JNK3** 为对象，基于 ChEMBL 公开生化 IC50 数据构建机器学习（ML）活性预测模型，对商业化合物库进行 **JNK 家族活性粗筛**；随后在本地 Schrödinger 环境中完成结构对接、MM-GBSA 与 ADMET 过滤，并对 **16 个候选分子** 开展 Desmond 分子动力学（MD）pose 质量控制（QC）。最终推荐 **10 个分子** 进入同批次 JNK1/2/3 酶学 IC50 湿实验。

**核心结论**：

1. **ML 模型适合判断“是否具有 JNK 家族结合潜力”**（9/9 文献 benchmark 在 p_family ≥ 6.0 时全部通过 F1 预筛），但 **不能可靠预测亚型选择性方向**（例如 E1、TCS JNK 6O 的预测最高活性亚型与实验不符）。
2. ChEMBL 配对数据中 **JNK1-selective 标注仅 8 个化合物**，选择性分类模型在测试集上 **F1 = 0**，不具备筛选选择性先导物的统计基础。
3. MD QC 表明：**G1（与文献 JNK1 抑制剂 chemotype 相近的类似物）3/4 通过整体 pose 质控，G2（新骨架）0/6 通过**；**没有任何候选分子的选择性可被计算确认**。
4. 花钱购买 10 个分子的理由不是“已算出 JNK1 选择性 hit”，而是：**用 G3 对照校准实验、用 G1 验证 chemotype 假说、用 G2 探索新骨架并结合 MD 信息价值最大化**。

---

## 1. 项目背景与目标

### 1.1 生物学与成药背景

JNK 家族包含三个高度同源的丝/苏氨酸激酶：

| 亚型 | 基因 | ChEMBL Target ID |
|------|------|------------------|
| JNK1 | MAPK8 | CHEMBL2276 |
| JNK2 | MAPK9 | CHEMBL4179 |
| JNK3 | MAPK10 | CHEMBL2637 |

JNK1 在代谢、炎症与纤维化等疾病中具有明确病理角色；临床上 **CC-90001** 为 JNK1 功能偏向的临床候选化合物（特发性肺纤维化 IPF），见 Bennett et al., *J. Med. Chem.* 2021（PMID: 33404223）。经典 pan-JNK 工具药 **SP600125** 见 Bennett et al., *Proc. Natl. Acad. Sci. USA* 2001（PMID: 11717429）。

### 1.2 项目初始目标与策略调整

**初始目标**：发现 **JNK1 亚型选择性** 小分子抑制剂。

**策略调整（基于数据）**：当 ML 与对接均无法可靠排序 isoform 选择性后，项目目标转为：

- 发现 **结构多样、计算上结合模式可信的 JNK 家族结合剂**；
- 通过 **G1（已知 chemotype 类似物）vs G2（新骨架）** 分组设计湿实验；
- 用 **G3 文献对照** 校准“计算 pose 质量 vs 真实活性”的相关性。

---

## 2. 训练数据与 ML 模型

### 2.1 数据来源与清洗

数据来自 ChEMBL 导出（`docs/JNK1.csv`, `JNK2.csv`, `JNK3.csv`），经 `scripts/00_prepare_user_data.py` 清洗。清洗后化合物数（`data/processed/data_summary.json`）：

| 亚型 | 化合物数 | 训练/验证/测试 |
|------|----------|----------------|
| JNK1 | 444 | 384 / 29 / 31 |
| JNK2 | 610 | 477 / 66 / 67 |
| JNK3 | 1147 | 966 / 83 / 98 |

清洗规则（详见 `docs/PROJECT_TECHNICAL_REPORT.md`）：

- 仅保留生化测定（Assay Type = B）、精确 IC50（Relation = '='）
- pActivity 范围 4.0–10.0
- 同 SMILES 多测定冲突：std > 0.5 或 range > 1.0 log 时丢弃
- 按亚型 assay harmonization：JNK1 assay ≥ 50 化合物、JNK2 ≥ 8、JNK3 ≥ 20

### 2.2 化学空间与选择性标签稀缺性

`results/similarity/` 分析显示三数据集化学空间部分重叠但可区分：

- 交叉 Tanimoto 中位相似度：JNK1–JNK2 **0.172**，JNK1–JNK3 **0.149**
- 三数据集共享 Murcko 骨架：**38** 个（`scaffold_overlap.json`）
- 配对分子（≥2 个亚型有数据）：**322** 个；其中 **JNK1-selective 仅 8 个**（`sel_class_counts.csv`）

这说明：**公开数据中“JNK1 选择性”标签极度稀缺**，是后续选择性 ML 失败的根本原因之一。

### 2.3 模型架构与性能

采用 **三个独立 XGBoost 回归器**，输入为 Morgan-2/2048 bit 指纹 + 12 个 RDKit 2D 描述符（共 2060 维），预测各亚型 pActivity。性能（`results/model_comparison/comparison.json`，骨架 holdout）：

| 亚型 | Holdout R² | Holdout Spearman ρ | n_test |
|------|------------|-------------------|--------|
| JNK1 | **0.697** | **0.858** | 31 |
| JNK2 | **0.574** | **0.780** | 67 |
| JNK3 | **0.774** | **0.869** | 98 |
| 平均 | **0.682** | **0.836** | — |

5-fold 骨架 CV 平均 R²：JNK1 **0.662**，JNK2 **0.443**，JNK3 **0.633**（`MODEL_COMPARISON_REPORT.md`）。

**解读**：模型对 **绝对活性排序** 在 JNK1/JNK3 上较可靠，适合 **粗筛**；JNK2 因 assay 异质性 CV 偏弱。XGBoost 平均 holdout R² **0.682**，优于 Chemprop（**0.532**，同报告）。

### 2.4 多任务/选择性模型（未用于筛选）

`results/training/training_report.json`：

| 任务 | 结果 | 结论 |
|------|------|------|
| MTL（Chemprop）JNK2/JNK3 holdout R² | ~0.25 | 不可用 |
| Δ 选择性回归（测试集） | n = **4**，R² = 0.64 | 样本过少，不具泛化性 |
| JNK1 选择性二分类 | 训练正例 **8**，测试正例 **0**，F1 = **0** | 不可用 |

---

## 3. 虚拟筛选漏斗（ML 阶段）

### 3.1 F1 阈值校准

在 **9 个文献 benchmark**（`data/benchmarks/literature_benchmarks.csv`）上扫描 p_family = max(pred_JNK1, pred_JNK2, pred_JNK3) 阈值（`results/calibration/threshold_scan.csv`）：

| 阈值 | Benchmark 通过率 |
|------|------------------|
| 6.0 | **9/9（100%）** |
| 6.5 | 7/9 |
| 7.0 | 5/9 |

推荐阈值：**p_family ≥ 6.0**（`threshold_recommendation.json`）。

### 3.2 漏斗设计（v2）

```
输入 SMILES
  → F0 预处理（RDKit 标准化、去重）
  → F2 Lipinski 类药（MW 200–600, logP −1~5, HBD≤5, HBA≤10）
  → F1 ML：p_family ≥ 6.0
  → F5 成药性：SA ≤ 6.0, QED ≥ 0.35
  → 综合评分 + Butina 多样性挑选
```

综合评分（`screening_v2/screening_report.json`）：

```
final_score = 0.55×(p_family/10) + 0.15×(pred_JNK1/10) + 0.20×QED + 0.10×(10−SA)/10
```

**设计理由**：用 p_family 保留 pan-JNK 工具化合物；**不用 ML 选择性硬筛**（见第 4 节）。

### 3.3 Demo 库筛选结果（仓库内基准）

在 ChEMBL 衍生的 demo 库（1835 分子）上（`screening_v2/screening_report.json`）：

| 阶段 | 数量 |
|------|------|
| 输入 | 1835 |
| Lipinski 通过 | 1541 |
| F1（p_family ≥ 6.0）通过 | 1292 |
| SA/QED 通过（最终 hit） | **1211** |

> 注：百万级 Enamine 库 ML 初筛在本地 WSL 完成（`docs/PROJECT_TECHNICAL_REPORT.md` 第 8.2 节），完整漏斗统计未纳入本仓库 Git；后续对接与 MD 以本地工作流为准。

---

## 4. 亚型选择性：问题、尝试与失败原因

### 4.1 问题一：ML 无法预测 isoform 方向

9 个 benchmark 的 ML 预测与实验 profile 对比（`benchmark_predictions_detailed.csv`）：

| 化合物 | 实验 profile | 实验 IC50 关系 | ML 预测最高亚型 | 方向是否正确 |
|--------|--------------|----------------|-----------------|--------------|
| E1 | JNK1-preferring | JNK1 2.7 nM；JNK2 19 nM（**7.0×**） | **JNK2**（pred 7.56） | **否** |
| TCS JNK 6O | JNK1-preferring | JNK2/JNK1 = 3.6× | **JNK3**（pred 6.97） | **否** |
| CC-930 | JNK2/3-biased | JNK1/JNK2 ~8.7× | JNK2（pred 7.47） | 是 |
| SP600125 | pan-JNK | 各亚型相近 | JNK1（pred 6.13） | 大致合理 |
| CC-90001 | pan-JNK（酶学） | JNK2/JNK1 = 2.8×（弱） | JNK2（pred 7.22） | 部分合理 |

E1 的实验数据来自 Pan et al., *J. Med. Chem.* **2024**, doi:[10.1021/acs.jmedchem.4c01764](https://doi.org/10.1021/acs.jmedchem.4c01764)（化合物 E1，JNK1 IC50 = 2.7 nM）。TCS JNK 6O 的 JNK1/JNK2 IC50（45/160 nM）见项目 benchmark 表，与文献/vendor 数据一致；原始报道见 Szczepankiewicz et al., *J. Med. Chem.* 2006, 49, 3563（PMID: 16759099）。

**结论**：F1 预筛可靠；**亚型方向预测不可靠**，因此 v2 漏斗明确放弃 ML 选择性过滤。

### 4.2 问题二：选择性训练标签极度不平衡

- 配对数据 322 个分子中，JNK1-selective 仅 **8** 个
- 选择性分类器：训练正例 8 个，测试正例 **0** 个 → F1 = 0
- SHAP 分析（`results/shap/top_shap_features.csv`）基于极少正例，仅作探索性参考

### 4.3 问题三：结构对接选择性评分（规划 vs 执行）

仓库中规划了 ensemble docking 方案（`config/docking_ensemble.yaml`）：

- JNK1：3ELJ + 4L7F 平均分
- JNK2：3E7O
- JNK3：3TTI + 4WHZ 平均分
- 选择性指标：Δsel = Score_JNK1 − max(Score_JNK2, Score_JNK3)

**MD QC 阶段明确不使用 Δsel、MD-ΔΔG 或 Gly87 IFP 进行选择性排序**（`MD_QC_report_cf26.md`），原因是前期 benchmark 验证表明对接选择性指标对 isoform 方向判别不足（该验证在本地对接报告中完成，具体数值未纳入本 Git 仓库）。

实际 MD 使用单结构：**3ELJ（JNK1）、3E7O（JNK2）、3TTI（JNK3）**，共 48 个 Desmond 任务（16 化合物 × 3 PDB）。

### 4.4 问题四：Gly87 / KLIFS 非保守位点策略未落地

曾规划利用 JNK1 特有 Gly87（KLIFS b.l.37）设计选择性接触，但 occupancy 自检未通过（该分析在本地完成，未纳入仓库）。MD QC 报告亦将 Gly87 IFP 排除在决策指标之外。

---

## 5. 结构筛选与 MD 短名单

### 5.1 进入 MD 的 16 个化合物分组

在完成 **ML 初筛 + 本地 Glide XP 对接 + MM-GBSA + QikProp ADMET** 后，**16 个化合物** 进入 MD pose QC（`MD_QC_report_cf26.md`）：

| 组别 | 数量 | 定义 |
|------|------|------|
| G1 | 4 | 与 E1/Q63/TCS JNK 6O chemotype 相似度较高的类似物（butina Tanimoto ~0.22–0.26） |
| G2 | 6 | 新骨架（butina 聚类不同，chemotype_sim 低） |
| G3_control | 4 | 文献对照：SP600125, CC-90001, CC-930, E1 |
| G4_anchor | 2 | 对接打分差的阴性锚点：3237, 3411 |

G1 的 chemotype_sim 来自 `md_pose_qc_summary_5ffb.csv`（以 E1 为参照时 E1 = 1.0；G1 约 **0.23**）。

### 5.2 MD QC 方法

- 软件：Desmond MD；分析 Amber cpptraj/MMPBSA
- 蛋白对齐后配体 RMSD：丢弃前 20% 平衡，分析剩余后 50%
- **通过标准**：
  - pass_md_JNK1：RMSD_JNK1 ≤ **3.0 Å** 且 hinge H-bond occupancy ≥ **30%**
  - pass_md_JNK2/3：RMSD ≤ **4.0 Å** 且 hinge HB ≥ 30%
  - pass_md_overall：pass_md_JNK1 **且**（pass_md_JNK2 **或** pass_md_JNK3）

---

## 6. MD Pose QC 结果

### 6.1 总漏斗

| 阶段 | 数量 |
|------|------|
| MD 输入化合物 | **16** |
| pass_md_JNK1 | **6** |
| pass_md_overall | **5** |
| 最终采购推荐 | **10** |

### 6.2 分组表现

| 组别 | n | pass_md_JNK1 | pass_md_overall | pose grade 分布 |
|------|---|--------------|-----------------|-----------------|
| G1 | 4 | 3 | **3** | A×3, F×1 |
| G2 | 6 | 1 | **0** | C×1, F×5 |
| G3_control | 4 | 2 | 2 | A×1, B×1, F×2 |
| G4_anchor | 2 | 0 | **0** | F×2 |

- G1 中位 RMSD_JNK1 = **0.64 Å**（3/4 通过 overall）
- G2 中位 RMSD_JNK1 = **0.65 Å**，但 **0/6 通过 overall**（pose 稳定性不足以支持三靶点结合模式一致）
- G4 两组均 fail → **阴性验证成立**（对接差 + MD 差）

### 6.3 G3 对照：MD 与活性的校准意义

| 对照 | 已知活性 | pass_md_overall | 主要失败原因 |
|------|----------|-----------------|--------------|
| CC-90001 | JNK 家族活性（酶学 pan-JNK，细胞 JNK1 偏向） | **通过** | — |
| E1 | JNK1 IC50 = 2.7 nM（强 JNK1 偏好） | **通过** | JNK3 hinge HB 偏低 |
| CC-930 | JNK2/3 强（7/6 nM） | 未通过 | JNK1 hinge HB ≈ 0% |
| SP600125 | pan-JNK（40/40/90 nM） | 未通过 | JNK1/JNK2 hinge HB < 30% |

**解读**：hinge HB ≥ 30% 门槛对 **经典 hinge binder** 合理，但对 **SP600125（吡唑并喹啉酮，非典型 hinge 氢键）** 和 **CC-930** 偏严。因此 **G3 必须全部购买**——它们的作用是 **酶学活性标尺**，而非 MD-QC 标尺。

### 6.4 G1 候选（MD 支持度最高）

| ID | pass_md_overall | RMSD_JNK1 (Å) | hinge_occ_JNK1 | JNK2/JNK3 pass | 备注 |
|----|-----------------|---------------|----------------|----------------|------|
| **690** | **是** | 0.72 | **100%** | **两者均通过** | 唯一三 isoform 全 pass 的 hit |
| 2232 | 是 | 0.57 | 100% | JNK2 pass | JNK3 hinge HB 低 |
| 2157 | 是 | 0.49 | 84.9% | JNK2 pass | JNK3 hinge HB 低 |
| 2389 | 否 | 1.00 | 28.2% | 均未 pass | JNK1 hinge 差 1.8% 至阈值 |

### 6.5 G2 候选

| ID | pass_md_overall | 说明 |
|----|-----------------|------|
| 2231 | 否（grade C） | 唯一 JNK1 grade A，但 JNK2/3 hinge 未通过 |
| 1280, 4795 | 否 | JNK1 fail，JNK2/3 pose 相对稳定（“off-target backup”假说） |
| 其余 3 个 | 否 | 三靶点 pose 均不可靠 |

---

## 7. 采购清单与花钱理由

完整表格见 `data/purchase/purchase_after_md.csv`（10 个分子，SMILES 经 RDKit 验证）。

### 7.1 采购结构

| 类别 | 数量 | 化合物 | 花钱理由 |
|------|------|--------|----------|
| **G3 对照** | 4 | SP600125, CC-90001, CC-930, E1 | **实验校准必需**：建立“酶学 IC50 vs MD pose”参考系；无论 MD 是否通过均购买 |
| **G1 主力** | 3 | 690, 2232, 2157 | **MD pass_md_overall = True，grade A**；验证“已知 chemotype 类似物是否更易出 JNK 活性” |
| **G2 探索** | 3 | 2231, 1280, 4795 | 2231 为 G2 中最优（JNK1-only pass）；1280/4795 用于检验“JNK2/3 pose 稳、JNK1 pose 不稳”是否仍对应 pan-JNK 活性 |

### 7.2 为什么不是“买选择性 hit”

必须向合作者/审稿人明确：

1. **没有任何一个采购分子的 JNK1 选择性被计算确认**；MD pass 仅表示 **结合位姿在动力学模拟中可信**。
2. 采购的核心信息目标是：
   - G3：校准实验体系
   - G1 vs G2：检验 **chemotype 假说**（G1 chemotype_sim ~0.23，与 E1 仍属低-中等相似，并非 E1 类似物）
   - 690：优先验证是否 **pan-JNK** 或 **JNK1 偏好**（三 isoform MD 全 pass）

### 7.3 预算优化建议（可选）

若需压缩至 8 个：可暂缓 **1280** 和 **4795**（MD 假说性强、JNK1 pose 已 fail）。

---

## 8. 当前“选择性”状况评估

### 8.1 计算层面：选择性不可判定

| 方法 | 对选择性的支持 | 依据 |
|------|----------------|------|
| ML 三模型 ΔpActivity | **不支持** | E1、TCS JNK 6O 方向错误 |
| 选择性分类器 | **不支持** | 正例 n=8，测试 F1=0 |
| 对接 Δsel | **未采用** | 本地 benchmark 方向判别差；MD 阶段明确排除 |
| MD pose QC | **不支持选择性** | 仅验证 pose；G3 中 2/4 MD-fail 仍有已知活性 |
| Gly87 IFP | **未采用** | 自检未通过 |

**诚实表述**：截至目前，**没有任何计算管线能对“JNK1 选择性”给出可采购级别的排序**；项目产出的是 **JNK 家族结合剂候选 + 实验设计**。

### 8.2 各候选分子的选择性“先验”

以下仅为 **待实验检验的假说**，不是计算结论：

| 分子 | 选择性先验 | 依据 |
|------|------------|------|
| 690 | 可能 pan-JNK 或弱 JNK1 偏好 | 三 isoform MD 均 pass；无实验 IC50 |
| 2232, 2157 | 未知；可能 pan-JNK | JNK1/JNK2 MD 好，JNK3 hinge 弱 |
| 2231 | 未知 | 仅 JNK1 MD pass |
| 1280, 4795 | 可能 JNK2/3 ≥ JNK1 | JNK1 MD fail、JNK2/3 相对稳定 |

---

## 9. 湿实验方案与结果预测

### 9.1 必做实验

**同批次 JNK1、JNK2、JNK3 重组酶生化 IC50**（10 个分子 + 建议复测 1 个 DMSO 空白），检测格式与 G3 文献条件尽量一致。

### 9.2 基于数据的预测（保守）

| 场景 | 预测 | 若成立的意义 |
|------|------|--------------|
| G3 有活性但 MD-fail（CC-930, SP600125） | **较可能** | 证明 hinge HB 门槛不能代替活性判断；计算 QC 需降级为辅助 |
| G1（690/2232/2157）至少 1 个 IC50 < 1 µM | **中等可能** | ML F1 + 对接 + MD 三级过滤后，具合理命中率 |
| G1 活性优于 G2 | **不确定** | G2 的 MD 支持度明显弱于 G1（0/6 vs 3/4 overall pass） |
| 出现 JNK1 选择性 ≥ 10× | **低可能** | 所有计算选择性方法均失败；无先验支持 |
| 690 为 pan-JNK | **需实验区分** | 三 isoform MD 全 pass 更支持“多亚型结合”而非“选择性” |
| G4 阴性锚点若被测 | 应无活性或 IC50 > 10 µM | G4 MD 0/2 pass，作阴性对照 |

### 9.3 后续（若有 hit）

对 top 1–2 分子考虑：**kinome 面板** + **细胞 p-c-Jun**（与 JNK 药理相关，见 Manning & Davis, *Nat. Rev. Drug Discov.* 2003, PMID: 12838265）。

---

## 10. 方法学贡献与局限

### 10.1 可写入论文的正面贡献

1. **激酶 isoform 选择性可行性评估框架**：从数据稀缺性（8 个 selective 标签）到 ML、对接、MD 的分层否定，系统记录 **何种计算证据不足以支持选择性采购**。
2. **benchmark 驱动的漏斗设计**：p_family ≥ 6.0 保证 9/9 文献化合物 F1 召回。
3. **MD pose QC + G3 校准实验设计**：将“pose 可信”与“选择性”解耦，避免过度解读。

### 10.2 主要局限

| 局限 | 说明 |
|------|------|
| JNK1 测试集小 | holdout n=31，R² 置信区间宽 |
| JNK2 模型偏弱 | CV R² = 0.443 |
| 本地对接漏斗未入 Git | 4983→…→16 等中间统计需从本地工作区补充归档 |
| MD 仅用单 PDB/isoform | 未执行 ensemble 平均 |
| 无湿实验数据 | 本报告预测均为假说 |

---

## 11. 数据与文件索引

| 内容 | 路径 |
|------|------|
| 文献 benchmark（9 个） | `data/benchmarks/literature_benchmarks.csv` |
| ML 阈值校准 | `results/calibration/` |
| 虚拟筛选 v2 | `results/screening_v2/` |
| 模型性能 | `results/model_comparison/` |
| 对接 ensemble 配置 | `config/docking_ensemble.yaml` |
| MD QC 汇总（本地归档） | `md_pose_qc_summary_5ffb.csv` |
| MD QC 报告（本地归档） | `MD_QC_report_cf26.md` |
| 采购清单 | `data/purchase/purchase_after_md.csv` |
| 参考文献列表 | `docs/REFERENCES.md` |

---

## 12. 参考文献

1. Zdrazil B, et al. The ChEMBL Database in 2023. *Nucleic Acids Res.* 2024;52(D1):D1180-D1192. doi:10.1093/nar/gkad1004  
2. Bennett BL, et al. Discovery of CC-90001, a JNK1-Selective Inhibitor for IPF. *J. Med. Chem.* 2021;64(3):1776-1795. doi:10.1021/acs.jmedchem.0c01843  
3. Bennett BL, et al. SP600125, an anthrapyrazolone inhibitor of JNK. *Proc. Natl. Acad. Sci. USA* 2001;98(24):13681-13686. doi:10.1073/pnas.251194298  
4. Pan X, et al. Structure Optimization of c-Jun N-terminal Kinase 1 Inhibitors for Treating Idiopathic Pulmonary Fibrosis. *J. Med. Chem.* 2024. doi:10.1021/acs.jmedchem.4c01764（化合物 **E1**）  
5. Szczepankiewicz BG, et al. Aminopyridine-based c-Jun N-terminal kinase inhibitors with cellular activity. *J. Med. Chem.* 2006;49(14):3563-3566. doi:10.1021/jm060150w（**TCS JNK 6o**）  
6. Plantevin-Krenitsky V, et al. CC-930/Tanzisertib JNK3 co-crystal (3TTI). *Bioorg. Med. Chem. Lett.* 2012;22(3):1433-1438  
7. Friesner RA, et al. Glide docking. *J. Med. Chem.* 2004;47(7):1739-1749. doi:10.1021/jm0306430  
8. Manning BD, Davis RJ. Targeting JNK for therapeutic benefit. *Nat. Rev. Drug Discov.* 2003;2(7):554-565. doi:10.1038/nrd1132  
9. Chen T, Guestrin C. XGBoost. *Proc. 22nd ACM SIGKDD* 2016. doi:10.1145/2939672.2939785  
10. Rogers D, Hahn M. Extended-Connectivity Fingerprints. *J. Chem. Inf. Model.* 2010;50(5):742-754. doi:10.1021/ci100050t  

---

## 附录：一句话给合作者/答辩用

> 我们用 ChEMBL 训练了 JNK1/2/3 活性 ML 模型（holdout R² 约 0.70/0.57/0.77），以文献 benchmark 校准的 p_family ≥ 6.0 对化合物库做 JNK 家族粗筛，再经本地 Glide 对接与 Desmond MD pose QC，从 16 个分子中选出 10 个采购。**计算无法确认 JNK1 选择性**——采购是为了用 G3 对照校准实验，并检验 G1 chemotype 类似物是否比 G2 新骨架更易获得 JNK 活性；**选择性只能由同批次 JNK1/2/3 IC50 回答**。
