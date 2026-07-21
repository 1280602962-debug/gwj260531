# 中文初稿（痛风 URAT1–NLRP3 双节点临床药物重定位）

> **⛔ 已过时 · 请勿再按本稿写作或投稿**（2026-07-21）  
> 本稿全文仍以 **Schrödinger Glide XP** 为主叙事，与当前 **V2 协议筛选（Vina/gnina/RTMScore → Π\*）** 路线冲突。  
> **请改用**：  
> - 总规划：[`MANUSCRIPT_OUTLINE_V2.md`](MANUSCRIPT_OUTLINE_V2.md)  
> - 中文大纲：[`MANUSCRIPT_OUTLINE_V2_CN_DRAFT.md`](MANUSCRIPT_OUTLINE_V2_CN_DRAFT.md)  
> - 引言：[`INTRO_DRAFT_CN.md`](INTRO_DRAFT_CN.md)  
> - Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)  
> 本文件仅作历史参考，保留 Glide 开发跑数字备查。

---

> ~~目标期刊：*Journal of Molecular Modeling*~~（已取消作为首投）  
> ~~对接引擎：Schrödinger Glide XP~~（已降级；非当前默认）  
> **结构**：URAT1 **9DKB** + NLRP3 **7ALV**（结构设定仍有效）  
> **MD**：代表体系仍可参考；lead 叙事以 V2（审计后假说）为准  
> **定位**：不对称双证据漏斗 + Pareto 短名单 + 代表药 MD；**不声称**发现首个双靶抑制剂

---

## 文章思路总览（一页纸）

### 科学问题

痛风同时涉及 **高尿酸（URAT1）** 与 **NLRP3 炎症轴**。能否在 **全临床药物库** 上，用标准化、可复现的计算漏斗，筛出 **双节点上均不被单轴支配** 的重定位候选，并用 MD 给出可讨论的构象假设？

### 核心叙事（五幕）

| 幕 | 内容 | 信息图 |
|----|------|--------|
| 1 | 数据不对称：URAT1 ML 不可靠作主筛，NLRP3 ML 宜缩库 | Fig 1 流程 |
| 2 | NLRP3 ML：8319→1588；痛风药/colchicine 对照 | Fig 2 |
| 3 | 8973 独立轨：证明 URAT1 **应对接** | Fig 3 |
| 4 | Glide XP 双靶对接 + Pareto → 6 分子；文献裁决 → **EGCG** | Fig 4 + 表 1 |
| 5 | 五组 MD 验证 benchmark + EGCG 双口袋 pose 稳定性 | Fig 5–6 |

### 三套数据（禁止混用）

| 数据集 | n | 用途 |
|--------|---|------|
| ChEMBL 临床 manifest | 8319 | NLRP3 ML → Glide XP@9DKB+7ALV → Pareto |
| 8973 distill | 8973 | **仅** URAT1 回顾富集 |
| Benchmark + MD | 5 体系 | 方法锚定 + EGCG 机制讨论 |

### 一句话结论（预定调）

我们建立了 **NLRP3 ML 缩库 + Glide XP 双靶对接百分位 Pareto** 的可复现重定位漏斗，并叠加 **结构警报/ADMET 过滤、模型 y-scrambling 与适用域验证、Pareto 稳健性** 等非对接计算证据；在短名单中，**EGCG** 因命中 PAINS 且落在 ML 适用域外被主动降级,漏斗转而提示 **canagliflozin** 等结构更干净的双节点候选;五组 MD 为 URAT1 基准药与 NLRP3 工具药提供构象参照，**不构成实验活性证明**。

---

## 题名（中文）

痛风相关 URAT1 与 NLRP3 靶点的临床药物重定位：NLRP3 机器学习预筛、Glide XP 双靶对接与代表抑制剂分子动力学研究

**英文题名（投稿用）**  
*Clinical drug repurposing for gout-related URAT1 and NLRP3 targets: NLRP3 machine-learning prescreening, Glide XP dual-target docking, and molecular dynamics of benchmark and lead compounds*

**Running title**：URAT1–NLRP3 双节点重定位漏斗

---

## 摘要

**背景** 痛风由尿酸沉积与 NLRP3 炎性小体激活共同驱动，临床降尿酸药与抗炎药分属不同机制节点。在 ChEMBL 临床化合物库中系统筛选兼具 URAT1 与 NLRP3 相关证据的重定位候选，尚缺乏可复现的一体化计算流程。

**方法** 对 8,319 个临床阶段小分子先进行 NLRP3 集成机器学习预筛（P(active)≥0.5，n=1,588），再对 inward-open URAT1（PDB 9DKB）与 NLRP3 NACHT（PDB 7ALV）并行 **Glide XP** 分子对接；以 URAT1 对接百分位（S_U）与 max(NLRP3 ML 百分位, NLRP3 对接百分位)（S_N）构建 Pareto 非支配前沿。URAT1 机器学习回归模型因 benchmark 回收不足（2/4）**不用于**主库排序；8,973 化合物蒸馏集单独用于 URAT1 对接富集回顾。对 Pareto 六分子进行文献裁决后，选定表没儿儿儿茶素没食子酸酯（EGCG）为双节点代表；并对 benzbromarone、dotinurad（URAT1 基准）、MCC950（NLRP3 基准）及 EGCG 在 9DKB/7ALV 上共 **五组** 分子动力学（MD）模拟（【待填：时长、力场、软件版本】）。

**结果** NLRP3 模型五折骨架交叉验证 AUROC 为 0.89。双靶 Glide XP 合并后 n=【待填：Glide XP 双靶合并数，开发参考 1451】，Pareto 前沿 **6** 个化合物；已知 URAT1 药 lesinurad、verinurad 的 S_U 较高但未落入 Pareto 前沿，colchicine 呈高 ML、低双轴特征。8973 回顾轨 ROC-AUC【待填：Glide XP 重算值；开发参考 0.705】，支持 URAT1 对接主导排序。六分子中 **EGCG** 为唯一 Phase 3 且具痛风/MSU–NLRP3 上游文献支撑的候选。MD 初步显示【待填：RMSD、关键残基接触、MM-PBSA 定性结论】。

**结论** 所构建的 **不对称双证据、标准化可复现** 漏斗可压缩临床库并产生可检验的双节点假说；结合化学过滤与模型验证,漏斗能主动识别并降级 EGCG 这类泛干扰命中,并提示更可跟进的候选。结果适用于计算药理学与重定位假设生成，仍需体外 URAT1 转运与 NLRP3 炎症实验验证。

**关键词** URAT1；NLRP3；痛风；药物重定位；Glide XP；机器学习；Pareto 优化；分子动力学

---

## 1 引言

痛风是全球最常见的炎症性关节炎之一，病理核心包括血尿酸升高、尿酸盐晶体沉积及 NLRP3 炎性小体活化所介导的 IL-1β 释放。临床治疗长期分为 **降尿酸**（如黄嘌呤氧化酶抑制剂、URAT1 促尿酸排泄药）与 **抗炎**（如秋水仙碱）两条路径，尚缺乏经充分验证的、同时作用于代谢与炎症双节点的单一药物。

URAT1（SLC22A12）是肾近端小管尿酸重吸收的关键转运体。近年来 inward-open 构象的 cryo-EM 结构（如 PDB **9DKB**）为结构导向的 URAT1 抑制剂研究提供了模板。NLRP3 方面，药物发现则面临检测体系异质、构象多样及作用机制不一（直接 NACHT 结合剂 vs 微管/上游调节剂）等挑战；MCC950 等磺酰脲类变构抑制剂与秋水仙碱等间接调节剂在机制上并不等价。

若对全临床药物库同时构建 URAT1 与 NLRP3 的机器学习排序模型，两类靶点的数据深度与标签噪声高度不对称：URAT1 回归模型对已知尿酸药的 benchmark 回收不完整，而 NLRP3 分类模型虽可用于缩库，却可能将 colchicine 等高 ML 分数的间接调节剂排在前列。因此，本研究采用 **不对称双证据漏斗**：以 NLRP3 ML **缩小化学空间**，以 **Glide XP 结构对接** 提供 URAT1 与 NLRP3 的并行结构证据，再以 **Pareto 优化** 筛选在双轴上均不被支配的分子，并通过独立的 8973 URAT1 回顾轨验证「URAT1 应对接、不对 ML 单独排序」的设计选择。

本工作 **不旨在** 发现首个双口袋 URAT1–NLRP3 抑制剂，而是提供一条 **标准化、可复现** 的临床库重定位流程（对接用 Schrödinger Glide XP,并配有开源脚本供无许可复现），并以 Pareto 与化学/模型过滤后的候选为讨论对象，辅以五组 MD 对基准药与 lead 的结合模式进行定性讨论。

---

## 2 材料与方法

### 2.1 临床重定位化合物库

自 ChEMBL 提取 max phase≥1 和/或具有 ATC 分类的临床化合物，去重并规范化为 canonical SMILES，得到 **8,319** 个小分子（`repurposing_manifest.csv`）。文献基准药（lesinurad、benzbromarone、verinurad、dotinurad、colchicine、allopurinol、MCC950 等）仅用于回溯评价与 MD，不参与重定位模型训练。

### 2.2 机器学习模型及其角色

**NLRP3 预筛分类器**：基于五个 NLRP3 相关生物测定数据，采用 Morgan 指纹与 RDKit 描述符的 XGBoost 集成模型，骨架分组五折交叉验证 AUROC **0.89**（AUPRC 0.91）。以 P(active)≥**0.5** 为操作阈值，得到对接池 **n=1,588**。该模型仅用于 **库缩减**，不作为 Pareto 层 NLRP3 证据的唯一来源。

**URAT1 回归模型（回顾用）**：五折 OOF R²=0.51，Spearman ρ=0.73；在严格阈值下仅回收 **2/4** 个 scaffold-novel 尿酸药 benchmark（lesinurad、dotinurad 未通过）。故 **不对 8,319 临床库做 URAT1 ML 排序**；该模型仅用于 8973 蒸馏集上的 ML–对接对比及方法学“负面对照”。

### 2.3 Glide XP 分子对接

**软件与参数** 采用 **Schrödinger Glide**（Extra Precision, XP）。受体经 Protein Preparation Wizard 预处理（加氢、指派键级、pH 7.0±2.0 质子化、限制性最小化）；配体经 LigPrep（Epik 质子化态、立体异构体枚举）准备。对接为柔性配体–刚性受体,每分子取 **最优 XP GScore** 作为对接分（kcal/mol，越低越好）。

**受体准备** 自 PDB **9DKB**（URAT1，inward-open）与 **7ALV**（NLRP3 NACHT，MCC950 类似物共晶口袋）提取蛋白链、去除水与非必要杂原子。

**网格盒** 以各自共晶配体质心定义 receptor grid（9DKB：lesinurad 中心；7ALV：共晶类似物中心），见 `config/docking_ensemble.yaml`。

**打分与排序** 每分子最优 pose 的 XP GScore 作为 `glide_score_xp`，仅在 **同一引擎、同一靶点、同一配体池内** 计算百分位，**不跨软件混比**。

**批量执行** Maestro/Glide 导出 → `normalize_canvas_docking_export.py` → `merge_docking_pareto.py`。

> **可复现性说明** 本文对接使用商业软件 Schrödinger Glide；*Journal of Molecular Modeling* 接受商业对接工具。为便于无许可复现,仓库同时提供开源 Vina/gnina 脚本(`run_vina_batch.py`、`run_gnina_batch.py`),但**本文所有报告数值均来自 Glide XP**,不与开源分数混用。共晶配体 redock 验证见 §2.9。

### 2.4 双靶合并与 Pareto 分析

对同时具备 URAT1 与 NLRP3 Glide XP 分数的分子（n=【待填】）：

- **S_U**：9DKB 对接分数的池内百分位（越高越好）  
- **S_N**：max( NLRP3 ML 百分位, 7ALV Glide XP 百分位 )（`--sn-mode both`）

在 (S_U, S_N) 平面上取 **非支配解** 作为 Pareto 前沿。附加阈值 `--min-su`/`--min-sn` 默认为 0，即短名单等于前沿全体。

### 2.5 URAT1 8973 蒸馏集回顾（独立轨）

8,973 个化合物（`distill_manifest.csv`）采用 **与主流程相同的 Glide XP@9DKB 协议**【若尚未重算则标注待办】，在活性子集 A（n=822）与诱饵子集 D（n=【待填】）上计算 ROC-AUC、EF@5%/10%，并报告 lesinurad、benzbromarone、verinurad、dotinurad 的对接与 ML 百分位。**该轨不参与 NLRP3 ML 或 Pareto 合并。**

### 2.6 分子动力学模拟

**体系（五组，进行中）**

| 编号 | 体系名称 | 蛋白 | 配体 | 角色 |
|------|----------|------|------|------|
| 1 | `benzbromarone_9DKB` | 9DKB | benzbromarone | URAT1 阳性基准 |
| 2 | `dotinurad_9DKB` | 9DKB | dotinurad | URAT1 阳性基准 |
| 3 | `9DKB_EGCG` | 9DKB | EGCG | Pareto 重定位 lead（URAT1） |
| 4 | `7ALV_MCC950` | 7ALV | MCC950 | NLRP3 直接抑制剂基准 |
| 5 | `7ALV_EGCG` | 7ALV | EGCG | Pareto lead（NLRP3 探索性 pose） |

**初始结构** 配体坐标取自 Glide XP 最优 pose；9DKB/7ALV 蛋白取自【待填：PDB 来源及预处理步骤】。  
**软件与力场**【待填：GROMACS / AMBER 版本；CHARMM36、AMBER ff14SB、OPLS 等】。  
**模拟参数**【待填：时长 ns；盒子类型；溶剂；离子浓度；温度/压强；约束】。  
**分析** 蛋白骨架 RMSD、配体 RMSD、关键距离（URAT1：Phe 笼、Arg477 等【待填残基编号】；7ALV：Walker B / 磺酰脲口袋【待填】）、氢键占有率、可选 MM-PBSA/MM-GBSA 相对结合自由能（定性比较，不作绝对亲和力）。

**说明** EGCG 在 7ALV 上的 MD 为 **探索性**：文献表明 EGCG 可能不以 MCC950 式直接占据 NACHT 磺酰脲口袋，而经上游通路调节 NLRP3；该体系用于讨论对接 pose 的 **动力学稳定性假说**，不单独作为直接抑制证据。

### 2.7 图与统计

- Fig 1：双节点 + 漏斗流程（8319→1588→双靶 Glide XP→Pareto→MD）  
- Fig 2：NLRP3 ML 分布与漏斗  
- Fig 3：8973 URAT1 富集  
- Fig 4：Pareto 散点（标注 lesinurad、colchicine、EGCG）  
- Fig 5：URAT1 MD（benz、dotinurad、EGCG）  
- Fig 6：NLRP3 MD（MCC950、EGCG）  
- Fig 7：非对接证据（化学空间 PCA + 适用域 + 结构警报/ADMET + Pareto 稳健性 + 候选提名）

统计：Pareto 为确定性非支配排序；ML 为五折 OOF；富集为 ROC-AUC / EF；ML 与对接 Spearman 相关用于方法辩护。

### 2.8 非对接计算过滤与模型验证

在不改动对接分与 Pareto 归属的前提下,对候选叠加以下 **下游注释与验证**（脚本见 `scripts/09–14`，汇总见 `docs/NON_DOCKING_COMPUTATION.md`）：

- **结构警报过滤**（模块 B）：RDKit FilterCatalog 计算 PAINS(A/B/C)、Brenk、NIH 命中,并用 logP/芳香环/TPSA 启发式标记胶体聚集风险。
- **ADMET/类药性**（模块 C）：MW、cLogP、TPSA、HBD/HBA、QED 及 Lipinski/Veber/Egan/Ghose 规则；作为 SwissADME/ADMETlab 的可复现替代。
- **化学空间与新颖性**（模块 D）：ECFP4 + PCA 投影;候选对 URAT1/NLRP3 已知活性的最近邻 Tanimoto。
- **模型严谨性**（模块 A）：URAT1 回归与 NLRP3 分类的 **y-scrambling 置换检验**（经验 p 值）、**适用域（AD）** 阈值(训练集内 NN Tanimoto 5% 分位)与 NLRP3 概率校准（Brier）。
- **Pareto 稳健性**（模块 E）：top-k% 双百分位交集、阈值门(τ=85/90/95)、bootstrap(500)前沿成员频率。
- **候选提名**（模块 F）：以 τ 双门放宽薄前沿,叠加"无 PAINS/Brenk + Lipinski/Veber 通过 + NLRP3 结构轴支撑"的透明过滤并排序,标注已知参照药以区分真正的新候选。

### 2.9 共晶配体 redock 验证【待补数据】

对 9DKB(lesinurad)与 7ALV(共晶类似物)做 **自对接**,报告最优 pose 相对共晶的 RMSD(目标 ≤2 Å),用于确认 Glide XP pose 的几何可信度;若偏大则在讨论中说明。【待填:各体系 redock RMSD】

---

## 3 结果

### 3.1 数据不对称支持漏斗设计

URAT1 ML 在四个尿酸药 benchmark 中仅稳定回收 2/4，而 NLRP3 分类 AUROC 达 0.89。ChEMBL 临床库与 8973 蒸馏集 **SMILES 零重叠**，故 URAT1 证据必须依赖独立对接回顾轨，而不能假设单一 merged 训练表。上述结果支持：**NLRP3 ML 缩库 + URAT1 Glide XP 过滤** 的不对称设计。

### 3.2 NLRP3 ML 预筛（对应 Fig 2）

8,319 个化合物中 **1,588**（19.1%）满足 P(active)≥0.5。痛风常用药 allopurinol、febuxostat、benzbromarone、dotinurad 的 NLRP3 概率接近 0，符合其非炎性小体机制；colchicine、verinurad 概率升高（约 0.92），提示 **表型/适应症混淆**。高 ML 分分子中 Phase I/II 占比偏高，需在讨论中作为方法局限而非生物学发现。

经双靶 Glide XP 后，**【待填】** 个分子同时具备 9DKB 与 7ALV 分数（开发参考 **1451**；缺失主要源于配体准备或对接失败）。

### 3.3 URAT1 8973 回顾（对应 Fig 3）

在子集 A vs D 上，Glide XP@9DKB 的 ROC-AUC 为 **【待填】**（开发参考 **0.705**），EF@5% 为 **【待填】**（开发参考 **4.23**）。dotinurad 呈现 **对接百分位高、ML 百分位极低**，与「URAT1 必须对接主导」一致。benzbromarone 在 8973 上对接回收良好，但因 NLRP3 ML≈0 **未进入** 主漏斗 P≥0.5 池，故通过单独 MD 与回顾轨表征。

### 3.4 Glide XP 双靶 Pareto（对应 Fig 4、表 1）

**表 1. Pareto 前沿六分子（Glide XP 百分位；数值【待填 Glide XP 终版】，下表为开发阶段参考）**

| 化合物 | S_U (9DKB %) | S_N (%) | 最高临床阶段 | 文献适配（主文判断） |
|--------|--------------|---------|--------------|----------------------|
| SLV-334 | 99.9 | 92.1 | 2 | 弱/无关痛风双节点【待补文献句】 |
| LANPROSTON | 99.9 | 96.8 | 2 | 弱/无关【待补】 |
| LASALOCID | 99.7 | 98.3 | 2 | **不推荐**：离子载体，文献提示可 **激活** NLRP3【待补 DOI】 |
| **EGCG** | 99.2 | 99.7 | **3** | **推荐 nominee**：MSU/痛风 NLRP3 上游 + URAT1 表达文献；多酚脱靶与 PK 局限 |
| FOSIGOTIFATOR | 98.7 | 99.8 | 2 | 磺酰脲样骨架；临床与痛风无关【待补】 |
| FOSRAVUCONAZOLE | 96.9 | 99.9 | 2 | 抗真菌；与双节点无关【待补】 |

**池内对照药（非 Pareto）**

| 药物 | S_U | S_N | 解读 |
|------|-----|-----|------|
| lesinurad | 91.6 | 95.0 | URAT1 对接高；S_N 被 ML 抬高 |
| verinurad | 77.7 | 97.9 | 同上 |
| colchicine | 30.7 | 50.1 | ML 高、URAT1 差 → 漏斗按设计剔除 |

Spearman ρ( P(active), 7ALV Glide XP 分 ) = **【待填】**（开发参考 **−0.04**），支持 S_N 取 max(ML, 对接)。

**EGCG 的 Pareto 定位**  
EGCG 在 1,588 池中 NLRP3 ML 百分位仅约 **12%**，S_N 几乎完全由 **7ALV Glide XP 对接** 驱动（S_N^dock≈99.7%），属「结构轴抬升、ML 轴一般」的 **对接主导型** Pareto 点，而非 ML 假阳性。

### 3.5 基准药与漏斗行为是否合理？

1. **URAT1**：8973 富集与 lesinurad 高 S_U 支持 Glide XP@9DKB；强 URAT1 药因 NLRP3 ML 未进池，由回顾轨 + MD 覆盖，而非指责其「Pareto 失败」。  
2. **NLRP3**：colchicine 高 ML 低双轴，证明漏斗能 **暴露 ML 混淆**。  
3. **六分子输出**：为 **计算非支配解**，非临床验证的双靶药；需文献剔除 lasalocid 等机制相反者。  
4. **与 canagliflozin 等**：非 Pareto 但文献双节点强的分子，可放讨论作「Pareto≠临床最优」对照【待补】。

### 3.6 分子动力学（对应 Fig 5–6）【进行中，待填数据】

**URAT1 @ 9DKB**

| 体系 | 骨架 RMSD (Å) | 配体 RMSD (Å) | 关键接触 | MM-PBSA 趋势 |
|------|---------------|---------------|----------|--------------|
| benzbromarone_9DKB | 【待填】 | 【待填】 | 【待填】 | 【待填】 |
| dotinurad_9DKB | 【待填】 | 【待填】 | 【待填】 | 【待填】 |
| 9DKB_EGCG | 【待填】 | 【待填】 | 【待填】 | 【待填】 |

**NLRP3 @ 7ALV**

| 体系 | 骨架 RMSD (Å) | 配体 RMSD (Å) | 关键接触 | MM-PBSA 趋势 |
|------|---------------|---------------|----------|--------------|
| 7ALV_MCC950 | 【待填】 | 【待填】 | 【待填】 | 【待填】 |
| 7ALV_EGCG | 【待填】 | 【待填】 | 【待填】 | 【待填】 |

**初步解读（占位，待模拟完成后改写）**  
- 若 benchmark 药 RMSD 稳定而 EGCG 漂移大 → 强调 EGCG 为 **假设性 pose**，需实验验证。  
- 若 9DKB_EGCG 稳定、7ALV_EGCG 不稳定 → 支持「URAT1 结构证据较强、NLRP3 以上游调节为主」的 **双节点调节剂** 叙事，与 Lee 2019 等一致【待补引用】。  
- MCC950@7ALV 作为直接抑制剂构象参照。

### 3.7 非对接计算证据与候选提名（对应 Fig 7、表 2）

**模型可信度（模块 A）** y-scrambling 显示两模型均学到真实信号:URAT1 真实 Spearman **0.732** vs 置换最大 **0.065**(经验 p≈0.048);NLRP3 真实 AUROC **0.891** vs 置换最大 **0.559**(p≈0.048);NLRP3 概率校准 Brier **0.128**。URAT1 训练集适用域阈值(NN Tanimoto 5% 分位)为 **0.578**,而 Pareto 六分子对训练集最近邻 Tanimoto 仅 **0.15–0.26**,即 **全部落在 URAT1 ML 适用域之外**——从数据上支持"URAT1 应对接主导、不对这些新颖分子做 ML 排序"。

**化学过滤（模块 B/C）** 短名单中 **仅 EGCG 命中 PAINS**(PAINS_B) 且同时命中 Brenk;其 ADMET 谱(TPSA 197、HBD 8、HBA 11,Lipinski/Veber/口服吸收均不过)印证多酚类"氢键过多、口服利用度低"的转化短板。对接池整体 41% 至少命中一类结构警报,提示过滤对短名单质量的必要性。

**Pareto 稳健性（模块 E）** 原始 6 分子前沿 bootstrap 频率 ~0.60–0.67(中等稳定);top-5%/10% 双百分位交集可将可讨论短名单扩至 **10/46**。

**候选提名（模块 F，表 2）** 在 top-10% 双门(46)上叠加"无 PAINS/Brenk + Lipinski/Veber 通过"过滤后得 **16** 个干净候选,排除已知参照药后 **15** 个,其中 **6** 个同时具备 **NLRP3 结构轴(7ALV 对接)支撑**（而非仅 ML 抬升）。据临床阶段与双结构轴平衡排序,靠前的干净新候选包括:

**表 2. 过滤后的替代双节点候选（top clean candidates；Glide XP 开发跑百分位）**

| 化合物 | 最高阶段 | S_U (%) | S_N 对接 (%) | P(active)_NLRP3 | QED | NLRP3 结构轴支撑 | 备注 |
|--------|:---:|:---:|:---:|:---:|:---:|:---:|------|
| **canagliflozin** | 4(上市) | 98.1 | 94.0 | 0.92 | 0.49 | 是 | SGLT2i,已有降尿酸/抗炎文献,双轴结构均高 |
| caficrestat | 3 | 91.5 | 95.0 | 0.57 | 0.54 | 是 | 醛糖还原酶抑制剂 |
| nelutroctiv | 2 | 92.4 | 94.4 | 0.92 | 0.56 | 是 | — |
| SLV-334 | 2 | 99.9 | 92.1 | 0.92 | 0.32 | 是 | 原 Pareto 前沿成员 |
| fulimetibant | 2 | 90.3 | 90.1 | 1.00 | 0.43 | 是 | — |
| PF-06952229 | 1 | 97.1 | 94.6 | 0.92 | 0.41 | 是 | — |

**参照药定位**:lesinurad(上市 URAT1 药)在同一双门内(S_U 91.6、S_N 95.0),作为阳性对照落在候选区,佐证提名逻辑合理。

> 提名为 **计算假设清单**,非临床推荐;canagliflozin 等仍需 URAT1 摄取与 NLRP3/IL-1β 实验验证,且 SGLT2 主适应症与痛风的关系需在讨论中审慎表述。

---

## 4 讨论

### 4.1 不对称双证据漏斗的方法学意义

本研究的核心贡献是 **在数据深度不均的条件下，仍能以可复现方式整合 URAT1 与 NLRP3 证据**。NLRP3 ML 高效缩库；URAT1 由 Glide XP 百分位主导；Pareto 层避免单轴支配。该设计与 colchicine、dotinurad 等对照药的表现相一致，回应了「为何不用一个统一 ML 模型筛完全库」的质疑。

### 4.2 Glide XP 对接的定位

采用 Schrödinger Glide XP 完成双靶虚筛,是转运体/变构口袋结构筛选的标准做法,*Journal of Molecular Modeling* 接受商业对接工具。对接分数仅用于 **池内相对排序**；本文 **不将** Glide XP GScore 等同于实验 K_i。为便于无许可复现,仓库提供 Vina/gnina 开源脚本,但报告数值均来自 Glide XP,不跨软件混比。共晶配体 redock（§2.9）用于确认 pose 可信度。

### 4.3 Pareto 六分子的生物学可读性

六分子均为 **计算前沿**，不等于临床推荐。经文献初筛：  
- **LASALOCID** 因可能 **激活** NLRP3，与抗炎叙事相反，应剔除。  
- **EGCG** 在六者中临床阶段最高，且有痛风/MSU 炎症与 URAT1 相关报道，适合作为 **双节点调节剂（dual-node modulator）** 讨论对象，而非 **双口袋直接抑制剂**。  
- 其余四分子与痛风双节点关联弱，主文列表报告、讨论中说明 **不宜作 lead**。

### 4.4 EGCG：盲筛入围、文献锚定、MD 边界

EGCG 由 **盲筛 Pareto（1/6）** 进入短名单，文献用于 **六选一可读性**，而非事后用文献反向驱动对接排名（post-hoc 辩护）。  
7ALV 对接针对 **NACHT 磺酰脲口袋**；文献中 EGCG 可不结合 PYD/ASC，而经 mtDNA/ROS 等 **上游** 调节 NLRP3——**不矛盾**。故 `7ALV_EGCG` MD 为 **pose 稳定性探索**；`9DKB_EGCG` MD 对 URAT1 重定位叙事更关键。  
EGCG 已知多酚脱靶、口服生物利用度低，讨论中须强调 **PK/安全性** 限制。

### 4.5 MD 五体系的设计逻辑

| 体系 | 目的 |
|------|------|
| benz + dot @ 9DKB | 锚定 URAT1 对接/MD 管线合理性 |
| EGCG @ 9DKB | 主文 lead 的 URAT1 结合假说 |
| MCC950 @ 7ALV | NLRP3 直接抑制剂阳性对照 |
| EGCG @ 7ALV | 探索性：对接 pose 能否在磺酰脲口袋稳定 |

未对 Pareto 其余五分子做 MD，符合「代表药 + 假说生成」篇幅定位。

### 4.6 EGCG 之外如何寻找候选（提名策略）

原始 Pareto 前沿仅 6 分子且成分与痛风双节点关联弱,单靠"前沿 + 文献裁决"易过度依赖 EGCG。我们改用一条 **可复现的多层提名策略**(模块 F)在不重对接的前提下扩展并净化候选:

1. **放宽前沿为双百分位门**:用 top-k%(如 90/95 双门)代替严格非支配集,把可讨论候选从 6 扩到 10–46,规避薄前沿的偶然性。
2. **结构与类药性净化**:剔除 PAINS/Brenk 命中、要求 Lipinski + Veber 通过,直接过滤掉 EGCG 这类泛干扰/低成药性分子。
3. **区分证据类型**:优先保留 **NLRP3 结构轴(7ALV 对接)也高** 的分子,而非仅靠 ML 百分位抬升者,降低适应症混淆风险。
4. **标注参照药**:已知尿酸药(如 lesinurad)落在候选区可作阳性对照,验证提名逻辑;真正的新候选单列。
5. **临床阶段 + 双结构轴平衡排序**:优先上市/后期、且两轴结构证据均衡的分子。

该策略提示 **canagliflozin**(上市 SGLT2 抑制剂,已有降尿酸与抗炎报道)、caficrestat、nelutroctiv 等结构更干净、双轴结构证据更均衡的替代候选。它们同样是 **计算假设**,须经 URAT1 摄取与 NLRP3/IL-1β 实验验证;canagliflozin 的 SGLT2 主适应症与痛风的机制关系亦需审慎讨论,避免过度外推。该提名流程完全脚本化、可复现,是对"仅报告 EGCG"的方法学补强。

### 4.7 局限性与未来工作

1. URAT1 ML benchmark 不完整；NLRP3 ML 对 colchicine 过敏感。  
2. Glide XP 为计算 pose 假设，无体外 URAT1 转运或 NLRP3 抑制实验。  
3. 8973 与临床库无重叠；8973 富集需与 Glide XP 终版分数同步更新。  
4. 膜蛋白 + 变构口袋的对接打分泛化性未知;Glide XP 对转运体柔性建模有限。  
5. EGCG 不能外推为「首个双靶痛风新药」。  

**未来**：URAT1 摄取实验、MSU 诱导巨噬细胞 IL-1β、候选(canagliflozin 等)体外验证与 PK 评估；对 9DKB/7ALV 的 redock 与网格敏感性分析、共识对接。

---

## 5 结论

本研究建立了面向痛风 **URAT1–NLRP3 双节点** 的 **NLRP3 ML 预筛 + Glide XP 双靶对接 + Pareto 优化** 可复现重定位漏斗，将 8,319 个临床化合物压缩为非支配前沿,并叠加 **结构警报/ADMET 过滤、模型 y-scrambling 与适用域验证、Pareto 稳健性** 等非对接计算证据；漏斗能主动识别并降级 EGCG 这类泛干扰命中,提示 canagliflozin 等更干净的候选;五组 MD 为 URAT1 基准药与 NLRP3 工具药提供构象层面的补充讨论。本工作为 *Journal of Molecular Modeling* 风格的计算重定位研究，强调 **可复现流程与可检验假说**，而非经实验验证的双靶药物发现。

---

## 支持信息（计划）

- NLRP3 OOF ROC/PR 曲线  
- URAT1 ML–对接不对称示意图  
- 全量 `pareto_merged_scores.csv`（Glide XP 终版）  
- Glide XP 参数、网格定义、配体/受体准备脚本  
- 五组 MD 拓扑、参数、轨迹分析补充图  
- Phase≥3 子集敏感性（可选）  

---

## 数据与代码可用性

GitHub：`1280602962-debug/gwj260531`，目录 `URAT1_NLRP3_DualTarget_AIDD_Project/`。Glide XP 对接输出、Pareto 表与作图脚本见 `data/repurposing/pareto/` 与 `scripts/`。

---

## 作者贡献 / 利益冲突

【待填】

---

# 附录 A：待你逐项补充的清单

## A. Glide XP 对接（优先）

- [ ] Glide/Maestro 版本号、XP 精度、网格中心与尺寸、redock RMSD  
- [ ] 双靶合并成功数 n（替代 1451）  
- [ ] Pareto 六分子 **Glide XP 终版** S_U、S_N 表  
- [ ] lesinurad / colchicine / EGCG 的 Glide XP 原始分与百分位  
- [ ] Spearman( ML, 7ALV Glide XP )  
- [ ] 8973 轨 Glide XP@9DKB：AUC、EF@5%、四药百分位（若已重算）  
- [ ] Fig 4 用 Glide XP 数据重绘确认  

## B. 五组 MD（你正在跑）

- [ ] 软件：GROMACS/AMBER 版本  
- [ ] 力场、水模型、离子浓度、模拟时长 (ns)  
- [ ] 各体系：蛋白/配体 RMSD 图数据  
- [ ] URAT1：Phe 笼 / Arg477 等关键距离–时间曲线  
- [ ] 7ALV：Walker B / 磺酰脲口袋关键氢键占有率  
- [ ] MM-PBSA 相对排序（benz vs dot vs EGCG；MCC950 vs EGCG）  
- [ ] 代表性 snapshot 图（各 1–2 张）  

## C. 图表

- [ ] **Fig 1** 流程图（8319→1588→Glide XP→Pareto→MD 五体系）  
- [ ] Fig 2–4：确认 Glide XP 版数字  
- [ ] **Fig 5** URAT1：benz、dot、EGCG 三体系 RMSD + 接触  
- [ ] **Fig 6** NLRP3：MCC950、EGCG 两体系  

## D. 文献（Discussion / 表 1）

- [ ] EGCG：痛风 MSU、NLRP3 上游（Lee 2019 等）  
- [ ] EGCG 与 URAT1 表达/尿酸（如有）  
- [ ] LASALOCID 激活 NLRP3（PNAS 2024 等）  
- [ ] 六分子其余 4 个：适应症/脱靶一句判词  
- [ ] canagliflozin：非 Pareto 但双节点文献（讨论用，可选）  
- [ ] 9DKB、7ALV 结构原始文献 DOI  

## E. 投稿元数据

- [ ] 作者单位、通讯作者  
- [ ] 英文 Abstract 定稿（由本中文稿翻译）  
- [ ] 致谢、基金号  

---

*本初稿与英文底稿 `MANUSCRIPT_DRAFT_CURRENT.md` 并行维护；对接数字以 Glide XP 终版为准后替换文中【待填】项。*
