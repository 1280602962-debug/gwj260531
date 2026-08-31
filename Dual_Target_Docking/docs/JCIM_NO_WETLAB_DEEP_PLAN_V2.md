# 无湿实验条件下冲刺 JCIM：独立验证与非激酶/PPI 扩展计划 v2

> 版本日期：2026-08-26  
> 适用稿件：*A Four-Pair Formulation Audit of Docking-Based Dual-Target Recognition*  
> 原则：本文是方法学评价稿，不声称发现新药或预测新活性；新计算必须检验预先写明的问题，不能按结果重抽样或调参。

## 1. 当前证据与剩余主风险

当前分支已经完成：四对方向性 pocket-matched AUROC、wrong-pocket 对照、配体性质与 ECFP 基线、scaffold-aware CV、ChEMBL max/median 与高置信记录敏感性、同批次 unused-pool holdout、两个 PIK3CA 相关靶对的受体实现敏感性、两对独立 GNINA 姿态生成、部分姿态/PLIF 审计及失败配体审计。

这些结果足以支持“评测构造和化学混淆会改变对双靶对接的解释”，但仍有四个明显审稿风险：

1. 所有实验状态主标签来自 ChEMBL；现有 BindingDB/PubChem 仅做供给计数，不是外部验证。
2. 主结果的明显 formulation gap 主要来自 EGFR/HER2，四对均不足以建立目标普遍规律。
3. 四对中三对涉及 kinase ATP 口袋，且两对共享 PIK3CA；AChE/BChE 已是非激酶对，因此新增 MCL1/Bcl-xL 的准确价值是增加 **PPI/BH3 groove** 域，而不是“首次加入非激酶”。
4. MCL1 与 Bcl-xL 同属 BCL-2 家族、具有同源折叠；不得称为 disparate-fold pair。它们的价值是非 ATP、浅而柔性的蛋白–蛋白相互作用槽，以及同一配体的双端晶体结构。

无湿实验版本的最强可行证据链应为：

> ChEMBL 开发/审计面板 → BindingDB-native、文献与结构去重的真正独立切片 → PPI/BH3 槽域外推 → 物理姿态与受体实现检查 → 全流程冻结与可复现发布。

## 2. 工作包 A：BindingDB 真正独立切片（最高优先级）

### A0. 为什么现有 BindingDB 计数不能直接使用

BindingDB 会导入 ChEMBL 数据，ChEMBL 也包含部分 BindingDB 策展数据；PubChem又接收二者沉积。因此“数据库名不同”不代表来源独立。当前 `jcim_supply_crossdb_v0` 的 REST 结果按 `monomerid` 聚合，只保留最大 p-activity，未保留来源 DOI/PMID、策展来源、assay、物种/突变、结构跨库映射，故只能支持供给审计。

独立切片必须从 BindingDB 的版本化 TSV/归档数据重新建造，不能从当前 pmax JSON 反推。

### A1. 预先冻结的独立性定义

每条候选记录必须同时满足：

- `Curation/DataSource` 为 BindingDB 原生策展的 article 或 patent；排除 ChEMBL、PubChem、PDSP、CSAR 导入行。
- 人源、目标序列可映射到指定 UniProt；主分析仅纳入 wild-type。突变体、融合构建和无法确定构建者单列敏感性，不混入主标签。
- 仅 IC50、Ki、Kd，且关系为 `=`；EC50、百分比抑制、模糊范围和不等式不进入主分析。
- 两个靶点均有实测值；缺测不作 inactive。
- 标准化结构按去盐最大有机片段生成 Standard InChIKey、canonical SMILES 和 Murcko scaffold。
- 来源 DOI/PMID/专利号不得出现在 ChEMBL 开发面板或其标签来源文献中。
- Standard InChIKey 不得出现在主 K=4、unused-pool holdout、PM110 或任何用于方法选择的分子集合中。
- 主外部集要求 max ECFP4 Tanimoto 对全部开发分子 `< 0.70`；`0.50–0.70` 与 `<0.50` 分层报告。若严格过滤导致样本不足，放宽结果只能作为 sensitivity，不得仍称 chemical-independent primary slice。

独立性应逐层报告，而非只给最终 n：原始 BindingDB-native → 去共享文献 → 去共享结构 → 去近邻化学系列 → 可成四状态/三状态面板。

### A2. 活性标签规则

避免沿用“跨 endpoint 最大值”作为唯一规则。建议：

1. 每个 ligand–target–endpoint 内取 exact measurements 的中位数。
2. 若同一 endpoint 的记录跨越强/弱两侧，标记 `discordant` 并从主分析排除。
3. 主标签沿用论文统一阈值 `theta = 6.0`，保持与主表直接可比。
4. strict 6.5/5.5 作为高置信敏感性；灰区剔除。
5. Ki/Kd-only、IC50-only、same-document/same-assay-family 分层作为 assay 稳健性，不在看到 docking 结果后选择“最好”的层作为主结果。

所有筛选、聚合和排除原因必须写入逐记录 provenance 表。

### A3. 靶对优先顺序与最低门槛

建议按下列次序做供给审计，供给门槛必须在任何新 docking 前冻结：

| 优先级 | 靶对 | 目的 | 现有 equal-only 供给提示 | 预定角色 |
|---|---|---|---:|---|
| 1 | EGFR/HER2 | 直接复核全文最强 formulation-gap 案例 | strict A/B 34/31（去重前） | 若去重后每臂 ≥15，作为 thin external replication；≥20 更理想 |
| 2 | AChE/BChE | 不共享 PIK3CA、非激酶、供给充足 | strict A/B 181/92（去重前） | primary external pair |
| 3 | PIK3CA/mTOR | 复核当前相对较好的方向性信号 | strict A/B 76/96（去重前） | EGFR 不足时的第二 primary pair |
| 4 | PIK3CA/PIK3CB | 同源激酶敏感性 | strict A/B 58/75（去重前） | 备选，不优先增加 PIK3CA 权重 |

Primary external pair 的最低门槛：dual、A-only、B-only 各 `n >= 20`，每类至少 3 个独立来源文献/专利，且最大单一文献的 ligand fraction 不高于 50%。Neither `n >= 10` 才做 Dual-vs-neither formulation contrast；否则只做方向性验证。EGFR/HER2 若每臂 15–19，可作为预先标记的 thin replication，但不进入跨对总结性主估计。

最小合格组合：AChE/BChE + PIK3CA/mTOR 两个 primary external pairs，并尽力加入 EGFR/HER2 thin replication。理想组合：EGFR/HER2 + AChE/BChE 均达 primary 门槛。

### A4. 冻结、盲化与 docking

在对接前完成并提交：

- `protocol/external_slice_contract.yaml`：数据版本、全部过滤、阈值、配额、随机种子、主/次指标、停止规则。
- `tables/external_candidate_flow.csv`：每层过滤后的 n。
- `tables/external_panel_*.csv`：含标签和 provenance 的锁定面板。
- `tables/external_panel_*.sha256`：冻结校验。
- `analysis/EXTERNAL_SLICE_FREEZE.md`：只讨论供给与可行性，不看 docking 分数。

对接复用冻结受体、盒、Vina 版本和 ligand preparation；主协议统一 exhaustiveness = 8。不得因某个外部对表现差而更换受体或调盒。受体替换仅作为预先规定的 sensitivity。

### A5. 统计分析

主指标仍为两个方向性 pocket-matched AUROC 及 `summary_min`，但补充：

- 每个 arm 同时报 AUROC、AUPRC、n、失败率和方向。
- 以 source document/patent 为 cluster 的 bootstrap 95% CI；另报 scaffold-cluster bootstrap，避免把同一 SAR 系列当独立重复。
- Dual-vs-neither 与 directional 指标使用不同负类，只报预先定义的 descriptive delta，不称为配对显著性检验。
- 报 ligand-only ECFP/physchem 基线及 docking 对 ECFP 的增量；外部集不重新训练或选择 docking score。
- 所有 docking failure 纳入 rank-extreme bounds，并报告其按类别、大小、柔性和来源的分布。
- 跨靶对只作 forest plot 和异质性描述；K 很小时不把随机效应均值包装成普遍规律。

预先写明判读，不把“显著”设为唯一成功条件：

- 若外部方向和效应大小与主面板相近：支持可转移的 failure mode。
- 若外部结果接近 0.5：支持“结论依赖面板/来源”，仍与论文主题一致，但须降低外推。
- 若外部结果反向：这是有价值的边界条件；检查来源、assay 和化学空间后如无错误，必须保留并据实改写。

## 3. 工作包 B：MCL1/Bcl-xL PPI/BH3 groove 扩展

### B0. 科学定位

MCL1/Bcl-xL 不是用来“增加一个好看的 AUROC”，而是检验四状态评测在非 ATP、PPI 槽、较大且柔性配体上是否仍暴露相同的 formulation/chemistry/receptor 依赖。AChE/BChE 已经是非激酶，因此该扩展的卖点必须写成 **binding-site class/domain shift**。

3WIY 与 3WIZ 均含同一 compound 10（PDB ligand LC6），分别结合 MCL1 与 Bcl-xL，是难得的双端 pose-gold。文献中的 compound 11 才是进一步优化后的更强双抑制剂，不要把 10/11 混写。

### B1. 面板供给与功效边界

现有 ChEMBL 审计在 theta = 6.0 下约为 dual 82、A-only 77、B-only 24、neither 122；strict 6.5/5.5 的最小 hard-negative 仅 12。因此：

- 主构建可采用 theta = 6.0，目标 24/24/24/24；B-only 将接近穷举，必须明示无独立同库 holdout。
- strict 分析仅作为薄样本 sensitivity，不得把其宽 CI 解读为机制差异。
- 优先先审计 BindingDB-native 是否能补充独立 B-only；若新增分子与 ChEMBL/LC6 系列高度重叠，不得以数量替代独立性。
- 对样本量做可检测效应模拟，正文明确该面板只能排除多大的效应。

### B2. 受体冻结与 pose-gold gate

Primary receptors：MCL1 3WIY、Bcl-xL 3WIZ。两者在挑选链时记录分辨率、缺失残基、配体完整度、晶体接触和 pocket residues；不能只任取 chain A。LC6 较大且柔性，RCSB ligand-validation 质量也应披露。

在任何面板 docking 前，必须通过：

1. LC6 双端 self-redocking，symmetry-aware heavy-atom RMSD 报 Top-1、Top-3、best-of-all。
2. PoseBusters 或等价物理有效性检查：键长/角、内部碰撞、蛋白–配体碰撞、口袋内位置。
3. 晶体 PLIF recovery：氢键/盐桥/芳香/疏水接触，报告 top pose 与 best-RMSD pose，避免只证明采样到了而排序失败。
4. 至少两个确定性随机种子；若结论对 seed 不稳定，扩大重复而不是挑最佳 seed。

Gate：两端均需 best-of-top3 RMSD < 2.0 Å 且通过物理有效性。若任一端失败，不把该对当成标准 screening-performance 证据；可作为“协议超出适用域”的预先定义 stress-test，但必须降级表述。不得为了过 gate 在看到面板标签/分数后反复调盒。

### B3. 受体实现与 PPI 槽特有控制

- 每个靶预先选一个替代 holo structure；选择依据只用结构质量、wild-type/construct、口袋完整性和配体位置，不用面板 AUROC。
- 报 primary/alternate 的 pocket C-alpha RMSD、体积、关键侧链和水分子差异。
- 主分析仍用刚性受体以保持跨对一致；关键水、侧链柔性和 protonation 的结果作为预先规定的敏感性，而不是为某个靶对定制最佳协议。
- 对 LC6/大分子单独报告 preparation 成功率、rotatable bonds、formal charge、tautomer/protomer 数和 docking 超时。

### B4. 分析必须与 K=4 完全同构

MCL1/Bcl-xL 需要生成与四对相同的全部核心输出：directional pocket-matched、wrong-pocket、Dual-vs-neither、mean/worst/min aggregation、ligand efficiency、descriptor/ECFP baseline、chemotype proximity、scaffold/document bootstrap、mixed-library enrichment、failure bounds 和 receptor sensitivity。不要只加入一张“第五对 AUROC”表，否则它不能解决泛化质疑。

## 4. 无湿实验还能显著提高质量的工作

按收益排序：

1. **真正独立的 BindingDB-native 验证**：唯一能直接改变“无外部验证”判断的工作。
2. **来源/assay 可审计标签**：逐记录 DOI/PMID、entry、assay、endpoint、relation、construct、聚合前后值和排除原因；同文献成簇统计。
3. **PPI/BH3 域外推**：在通过 pose-gold gate 后加入 MCL1/Bcl-xL，降低 ATP-pocket 主导性。
4. **全姿态物理有效性与 PLIF recovery**：不仅检查 cognate RMSD，还统计每个面板 top pose 的 PoseBusters pass rate、clash 和 pocket occupancy。近期 docking benchmark 已将物理有效性与 RMSD并列。
5. **预注册式冻结**：在新 docking 前提交 machine-readable contract、面板哈希、primary endpoint 和 failure handling。即使不是正式注册，也能证明没有按结果调数据。
6. **发布级复现**：Zenodo DOI、release tag、完整输入/输出/姿态、容器或锁定环境、checksum manifest、CI 运行零对接分析和一小套 docking smoke test。
7. **功效与不确定性**：document/scaffold cluster CI、可检测效应模拟、AUPRC、失败极端界；不再依赖点估计故事。
8. **外部文献对照表**：逐项比较 Zhou 2013、DUD-E、LIT-PCBA、CASF、DOCKSTRING 与本文在 negative definition、paired measurements、chemical controls、receptor sensitivity 和 externality 上的差异。

低收益或可能稀释文章的工作：继续堆多个重打分模型、把短 MD 当活性验证、只增加相似 kinase 对、只报告 MM/GBSA 点值、根据结果优化阈值/盒/受体、把 PubChem 当第二个独立数据库重复 BindingDB 数据。

## 5. 推荐执行顺序和停止规则

### Phase 1：只做供给与独立性审计，不 docking

1. 下载并锁定 BindingDB archived TSV、assay mapping 与 checksum。
2. 重建四对 + MCL1/Bcl-xL 的逐记录 provenance。
3. 完成 source/DOI/structure/scaffold/ECFP 去重流图。
4. 冻结两个 external primary pairs；若 EGFR 达 thin 门槛一并冻结。
5. 冻结 MCL1/Bcl-xL panel 和两个 primary receptor chains。

停止规则：若去重后不足两个 primary external pairs，先做 document-blocked ChEMBL 分析并把稿件定位维持为 formulation audit；不得把不足门槛的 BindingDB 集包装成 external validation。

### Phase 2：先跑小规模结构 gate

1. MCL1/Bcl-xL LC6 双端 redock、PoseBusters、PLIF、seed sensitivity。
2. 外部两对各抽不带标签的 5–10 个分子做 preparation/docking smoke test，只查失败和协议兼容性，不计算 AUROC。

停止规则：LC6 gate 失败则 MCL1/Bcl-xL 降为 applicability stress-test；不要先跑完整面板再寻找能过的受体。

### Phase 3：冻结后完整对接

1. 至少两个 BindingDB independent pairs。
2. MCL1/Bcl-xL primary receptor panel。
3. 预定 alternate receptor sensitivity。
4. 一次性执行冻结分析脚本，生成全部表和失败清单。

### Phase 4：稿件重构

若 external 证据合格，标题可升级为：

> *Benchmark Formulation, Chemical Confounding, and External Transfer in Docking-Based Dual-Target Recognition*

若 external 结果高度异质，保留当前 formulation-audit 标题，把异质性作为边界结论。任何情况下都不要写“validated biological recognition”“prospective utility”或“general reliability”。

## 6. 文献依据（用于方案，不替代最终逐条核对）

- JCIM scope：cutting-edge methodologies and applications in chemical informatics and molecular modeling。<https://pubs.acs.org/jcisd8/pages/info-for-authors>
- Zhou, Li, Hou 2013，双激酶 docking 基准：<https://doi.org/10.1021/ci400065e>
- DUD-E，property-matched decoys：<https://doi.org/10.1021/jm300687e>
- LIT-PCBA，实验 assay 标签与偏倚控制：<https://doi.org/10.1021/acs.jcim.0c00155>
- CASF-2016，多维 docking/scoring benchmark：<https://doi.org/10.1021/acs.jcim.8b00545>
- DOCKSTRING，标准化流程、规模和 target diversity：<https://doi.org/10.1021/acs.jcim.1c01334>
- DataSAIL，生物数据 leakage-reduced splitting：<https://doi.org/10.1038/s41467-025-58606-8>
- BindingDB 2024，数据来源、文献/专利 entry、与 ChEMBL 双向关系：<https://doi.org/10.1093/nar/gkae1075>
- Tanaka et al.，MCL1/Bcl-xL compound 10 双端晶体与 compound 11 活性：<https://doi.org/10.1021/jm401170c>
- Schaller et al.，kinase cross-docking 与受体实现：<https://doi.org/10.1021/acs.jcim.4c00905>
- 近期 JCIM pose sampling/physical-validity 对照：<https://doi.org/10.1021/acs.jcim.5c00380>

## 7. 完成定义

只有同时满足以下条件，才可在摘要中写 “external evaluation”：

- 至少两对通过 A1 独立性过滤及 A3 primary 门槛；
- 面板在 docking 前冻结并有 Git commit/哈希；
- 主结果含 document/scaffold cluster uncertainty 和完整失败处理；
- 数据、标签来源、代码、受体、盒、姿态和环境以不可变 DOI 发布；
- 结果无论支持、接近随机或反向均完整报告。

MCL1/Bcl-xL 只有在 LC6 双端 pose-gold gate 通过后，才可作为 target-domain extension；否则只能作为明确标注的 PPI groove applicability failure/stress-test。

