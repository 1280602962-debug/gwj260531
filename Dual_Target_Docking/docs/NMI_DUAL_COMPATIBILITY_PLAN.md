# 双靶兼容性打分冲击 NMI：完整思考、分析与执行方案

> **状态（重要）：** 本文档保留「分数融合 / dual-compatibility」历史主线。  
> **现行 NMI 投稿规划**已切换为 **双药效团 / passenger 污染 + moiety-resolved 评测**：  
> → [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md)  
> 本文中的 fusion / softmin 等内容降为 **基线与组件**，不再作为封面级创新。

本文档归纳本课题早期关于 **Nature Machine Intelligence（NMI）导向** 的思考：问题定位、与协同/DTI/对接算法的区别、数据与对接协议、难点、细胞/PK 数据用法等。

相关已有文档：

- **现行投稿规划：** [`NMI_SUBMISSION_PLAN_MOIETY.md`](NMI_SUBMISSION_PLAN_MOIETY.md)
- 问答总结：[`DUAL_TARGET_DOCKING_QA_SUMMARY.md`](DUAL_TARGET_DOCKING_QA_SUMMARY.md)
- 实现蓝图：[`DUAL_TARGET_SCORING_IMPLEMENTATION.md`](DUAL_TARGET_SCORING_IMPLEMENTATION.md)
- 文献与分子：[`REFERENCES_AND_MOLECULES.md`](REFERENCES_AND_MOLECULES.md)
- **NMI/高分对标流程：** [`NMI_REFERENCE_PAPER_PLAYBOOK.md`](NMI_REFERENCE_PAPER_PLAYBOOK.md)
- 共晶说明：[`DUAL_TARGET_COCRYSTAL_CATALOG_NOTES.md`](DUAL_TARGET_COCRYSTAL_CATALOG_NOTES.md)

---

## 0. 一句话定调（全文最重要）

**现行定调见 moiety 规划。** 下文为历史表述：

**计算主线不是新对接算法，也不是药物协同预测，而是：**

> **双靶兼容性打分与排序（dual-compatibility scoring）**  
> + **配套开放评测基准（Dual-VSDS）**  
> + **以 NLRP3/JNK1 fused & linked 细胞数据作时间盲测外部验证**

普通双靶分子 = 同一配体分别形成两个独立二元复合物，再融合两侧信号。  
**不是** PROTAC 三元桥接；**不是**两药联用协同（DeepSynergy 那条线）。

---

## 1. 文章在算什么、不算什么

| 概念 | 含义 | 本课题 |
|------|------|--------|
| 对接采样器 | 姿态搜索（Vina/DiffDock 等） | **复用现成引擎**，不自研 |
| DTI / 亲和力预测 | 单靶结合强弱 | 作**组件/校准输入** |
| 药物协同 | 两药联用 1+1>2 | **不做** |
| **双靶兼容性** | 同一分子两端是否同时够好，并压过 A-only/B-only | **正文核心** |

文章结构建议：**约 80% 方法+基准，20% 应用验证（公开靶点对 + 私有盲测）**。

---

## 2. 为何对标 NMI：高分刊买什么

近三年 NMI（含药学相关）常见创新模式：

| 模式 | 含义 | 本课题对应 |
|------|------|------------|
| A 纠偏/打假 | 证明主流做法系统性误导 | 两端分数相加/平均 → 假双靶、假阳性放大 |
| B 任务重定义 | 新目标函数贴近真实决策 | dual-compatibility + dual-vs-single |
| C 组合式方法 | 零件已知，组合解决新约束 | 校准 + softmin + 硬负样本 + 类型条件 |
| D 资源+闭环 | 开放基准 / 工具 / 外测 | Dual-VSDS + NLRP3/JNK1 holdout |

**零件可以全部已知**（GNINA、RTMScore、PoseBusters、Z-score、softmin）。  
**创新 = 第一次把双靶评成「校准后的短板约束任务」，并证明朴素融合失败。**

可写进 cover letter 的主张：

> We do not introduce a new docking sampler. We identify that dual-target assessment is an ill-posed score-aggregation problem, introduce a dual-compatibility learning objective with architecture-aware evaluation, and release a leakage-controlled benchmark showing that standard fusion systematically promotes single-target–biased molecules.

---

## 3. 课题真正难点（专家视角）

1. **标签语义**：细胞 ≠ 双靶直接结合；缺亲和力时必须分层标注。  
2. **配对数据稀疏**：两端都有可靠活性的分子少；fused/linked 标注更少。  
3. **硬负样本**：A-only/B-only 才是双靶难例；随机 decoy 会虚高。  
4. **跨靶分数不可比**：口袋不同，原始分不能加。  
5. **活性不均衡**：药效不对称 + 标签不均衡 + 难度不对称，导致强端霸权。  
6. **fused ≠ linked**：物理机制不同，不能一套头硬套；linked ≠ PROTAC 三元。  
7. **评测易泄漏**：随机划分、同系列泄漏会制造假 SOTA。  
8. **细胞好、PK 差**：必须解释暴露脱钩，不能假装预测成药性。

过不了「标签 / 硬负样本+严格划分 / 类型分治」三关，文章会退化成分数加权案例。

---

## 4. 跨靶分数不可比 & 活性不均衡：协议级对策

### 4.1 分数不可比

- **最低合格**：每靶 rank / percentile / Z-score。  
- **更好（主协议）**：每靶用已知活/非活校准到 \(\hat p_t(\mathrm{active})\)（isotonic / Platt / 小模型）。  
- **再融合**：短板敏感，例如  
  \(m_t=\hat p_t-\theta_t\)，\(S_{\mathrm{dual}}=\mathrm{softmin}(m_A,m_B)\)。  
- 消融必须含：raw / z / rank / \(p(\mathrm{active})\)。

### 4.2 活性不均衡（三种问题分开）

| 名称 | 含义 | 后果 | 对策 |
|------|------|------|------|
| Potency asymmetry | 一端 nM、一端 μM | 均值被强端绑架 | softmin、端特异 \(\theta_t\) |
| Label imbalance | dual 少、A-only 多 | 指标虚高 | 重加权、PR-AUC、hard negatives |
| Difficulty asymmetry | 一端对接更噪 | 放大噪声端 | 不确定性门控、分靶报告校准 |

---

## 5. 算法贡献：什么算新

### 5.1 前人已做（较初级）

- Perez-Castillo：多打分融合双靶 VS  
- Zhou 2013：双激酶对接假阳性高  
- PARP1–BRD4 2024：单靶点对 merged 案例  
- DualDiff/FuseDiff：双靶**生成**，非兼容性判别  

**缺口：** 校准 + 短板目标 + dual-vs-single 硬负样本 + fused/linked 条件 + 开放配对基准 + 严格泛化，尚未被系统做成统一框架。

### 5.2 主张的创新（选 1–2 条主打）

1. 双靶兼容性任务形式化与损失（短板 + dual>A-only/B-only）  
2. 结构类型条件专家（fused/linked；你有独家数据）  
3. Dual-VSDS 开放基准与泄漏控制协议（资源型）  
4. （次要）多层活性解耦：靶点层 vs 细胞层；PK 仅案例分析  

**不当主创新：** 再训通用 pose scorer / 新 sampler。

---

## 6. 数据准备方案（要什么、从哪来）

### 6.1 四主表 + 一盲测

| ID | 名称 | 作用 |
|----|------|------|
| **D1** | 姿态/共晶集 | RMSD、PB-valid、PLIF |
| **D2** | 配对活性集 | 主训练/主评测（dual / A-only / B-only / inactive） |
| **D3** | 设计类型元数据 | merged/fused/linked |
| **D4** | 负例协议 | TrueNegative / Hard / RandomDecoy |
| **D5** | NLRP3/JNK1 私有 holdout | 时间盲测；永不训练 |

### 6.2 D2 标签规则（锁死）

- `dual`：两端 biochemical/binding，`pAct ≥ θ`（如 6）  
- `A_only` / `B_only`：一端强、另一端有实测且弱  
- `inactive`：两端都有实测且都弱（**禁止把未测当 inactive**）  
- cell-only 文献数据：弱标签/附录，不与 binding 混训  

### 6.3 数据来源

| 来源 | 用途 |
|------|------|
| ChEMBL / BindingDB / PubChem | 配对活性骨干（程序化） |
| PDB / PDBbind | D1 共晶 |
| **已发表双靶论文** | **必须**：标 design_type、确认真双靶设计、补系列/PDB |
| 实验室 NLRP3/JNK1 | D5 盲测（细胞分层；缺亲和力可接受） |

**文献要不要抠？要。**  
库解决“两端数字”；文献解决“是否真双靶设计 + fused/linked 标注”。  
推荐：**D2-public（大）+ D2-curated 文献金集（小而真）**。

### 6.4 私有细胞数据分层

| 层 | 内容 | 角色 |
|----|------|------|
| L1 | 直接结合 IC50/Kd | 最好；暂缺则用公开 L1 替代训练 |
| L2 | 通路读出 | 辅助/盲测 |
| L3 | 细胞表型 | 盲测，不当双靶金标准 |
| L4 | ADME/PK | 脱钩案例分析，不训 PK 预测 |

细胞好、PK 差 → 写成：兼容性/细胞转化可预测，暴露失败导致体内转化失败；**不声称候选药物或 PK 预测模型**。

### 6.5 划分协议（防泄漏）

scaffold；leave-series-out；leave-target-pair-out；fused↔linked 外推；D5 time split。

### 6.6 主指标

dual-vs-single pairwise accuracy；EF@0.5/1/5%；PR-AUC；校准（Brier/ECE）；姿态 RMSD+PB-valid；**两套 decoy 设定都报**。

---

## 7. 对接方法选择

### 7.1 主协议（Methods 默认）

- 采样：**GNINA**（Vina/smina 搜索）或 AutoDock Vina  
- 重打分对照：**RTMScore**（或 GenScore）  
- 门控：**PoseBusters**  
- 每靶 top-K 姿态；指定口袋 box；两端**独立**对接  
- 输出特征：vina / cnnscore / cnnaffinity / rtmscore / pb_valid → 再校准融合  

### 7.2 对照引擎

DiffDock / CarsiDock / KarmaDock 等作消融，**不作唯一主引擎**（VSDS-VD：精度≠物理合理≠TrueDecoy 富集）。

### 7.3 受体

每靶 1–3 个 holo 构象；统一质子化 SOP；NLRP3/JNK1 用现有 PDB。

---

## 8. 最小模型结构（算法落地）

```text
分子 → 对接A/打分A + 对接B/打分B + design_type
     → 每靶校准 p_A, p_B
     → softmin / 阈值边距 +（可选）类型专家
     → S_dual
损失：L_A + L_B + λ1 L_dual-vs-single + λ2 L_cell(辅助) + λ3 L_calibration
```

基线：raw mean / min / z-mean / rank fusion / RTMScore 双侧简单融合。

---

## 9. 外部实验现实约束下的写法

- 无力补亲和力/动物药效：**可以投方法文**，主张限定在兼容性与细胞转化。  
- D5 = 时间盲测；公开 D2 binding = 主结合证据。  
- 已有细胞数据尽量挖：母核 vs 联用 vs fused vs linked 对照；毒性并行；细胞–PK 排序脱钩。  
- 若只能补很少实验：CETSA/占有 > 溶渗稳 > 一轮前瞻细胞验证 ≫ 动物药效。

---

## 10. 仓库内课题优先级（纯计算冲高分）

| 课题 | 判断 |
|------|------|
| **Dual_Target_Docking / Dual-VSDS** | **主线**：叙事最贴 NMI |
| URAT1–NLRP3 DualTarget | 可作同一框架的第二案例，指标对齐 |
| JNK1 Selectivity | 工程最全；适合先发稳妥 JCIM，或作双靶框架的激酶案例 |
| 共价表/综述 | 非本方法文主线 |

策略：双靶兼容性为主；JNK/激酶对作公开案例；NLRP3/JNK1 私有盲测；避免五线并行。

---

## 11. 后续执行方案（建议顺序）

### Phase 0｜立项锁定（文档）

- [x] 问题定调、创新模式、数据/对接原则（本文档）  
- [ ] 冻结标签阈值 θ、主指标、主引擎版本  

### Phase 1｜数据

1. ~~选定 2–3 个公开靶点对~~ → **已冻结** PIK3CA/mTOR、EGFR/HER2、Mcl-1/Bcl-xL（见 [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)）；下一步 ChEMBL 配对抽取 → D2-public  
2. 文献 curated 30–80 分子 → D2-curated + design_type 
3. 扩展 D1 共晶 catalog  
4. 入库 D5（私有细胞）+ 与训练 InChIKey 去重  
5. 实现划分与负例脚本  

### Phase 2｜对接协议

1. YAML 固定受体/box/随机种子  
2. 跑通 GNINA（± RTMScore）+ PoseBusters  
3. 每靶校准曲线；导出 \(p_t\)  

### Phase 3｜打假表（创新种子）

在严格 split 上报告朴素融合：

| 方法 | dual 召回 | A-only 混入 Top1% | fused/linked 差 |
|------|-----------|-------------------|-----------------|
| mean | | | |
| min | | | |
| z-mean | | | |

**A-only 混入高 → 纠偏贡献成立。**

### Phase 4｜融合头 + 消融

校准、softmin、hard negative、类型条件逐步加；只认严格 split 上的显著提升。

### Phase 5｜盲测与叙事

D5 前瞻；细胞–PK 脱钩图；开源 Zenodo/GitHub；先瞄准可投稿包（亦可以 JCIM 为台阶）。

### Phase 6｜投稿包图表（NMI 向）

1. 任务定义与泄漏/标签问题  
2. 朴素融合系统性失败  
3. 方法总图  
4. Dual-VSDS 多 split 结果  
5. 类型消融与可解释  
6. NLRP3/JNK1 时间盲测  
7. 细胞–PK 脱钩 / Pareto  

---

## 12. 数据表头模板（落地用）

### D2（配对活性）建议列

```text
compound_id, inchikey, smiles,
target_A, uniprot_A, chembl_target_A,
target_B, uniprot_B, chembl_target_B,
pAct_A, assay_type_A, assay_id_A, doc_id_A,
pAct_B, assay_type_B, assay_id_B, doc_id_B,
label_quad, design_type, series_id, source, year
```

### D5（私有盲测）建议列

```text
compound_id, inchikey, smiles, design_type, synth_date, batch,
endpoint_name, value, unit, conc, n_rep, cytotoxicity,
has_warhead_control, has_combo_control, split=holdout
```

---

## 13. 红线与可声称边界

| 可以声称 | 不要声称 |
|----------|----------|
| 双靶兼容性排序优于朴素融合 | 发现可开发临床候选 |
| 降低 A-only 假双靶混入 | 模型可靠预测动物 PK |
| 细胞层盲测外推 | 细胞活性 = 两端直接结合（除非有 L1） |
| linked PK 风险的计算解释 | linked = PROTAC 三元问题 |

---

## 14. 文档维护

- 本文档随 Phase 推进更新勾选状态与关键数字（混入率、EF、校准误差）。  
- 与 `DUAL_TARGET_SCORING_IMPLEMENTATION.md` 冲突时，以**本文档的任务定义与评测原则**为准；实现细节以 implementation 文档补充。  
- 文献条目继续追加至 `REFERENCES_AND_MOLECULES.md`。
