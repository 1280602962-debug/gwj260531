# 参考文献链接与双靶分子信息汇总

本文件汇总 `Dual_Target_Docking/` 课题中引用的**全部文献链接**（含内容简介），以及已搜集的**双靶分子/共晶种子目录**全文信息。  
机器可读表：[`../data/dual_target_structures/dual_target_cocrystal_catalog.csv`](../data/dual_target_structures/dual_target_cocrystal_catalog.csv)。

---

## 一、方法学与综述文献（按主题）

### 1.1 多靶配体设计分类与成药性

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R01 | Morphy R, Rankovic Z. *The physicochemical challenges of designing multiple ligands.* **J. Med. Chem.** 2006. | [doi:10.1021/jm0603015](https://doi.org/10.1021/jm0603015) | 经典综述：提出/巩固 **linked / fused / merged** 设计分类，系统讨论多靶配体的理化性质、分子量与成药性挑战。本课题设计类型标注的基础。 |
| R02 | Proschak E et al. *Polypharmacology: A Systematic Investigation of Dual-Target-Directed Ligands.* **J. Med. Chem.** 2024. | [doi:10.1021/acs.jmedchem.4c00838](https://doi.org/10.1021/acs.jmedchem.4c00838) | 近年对双靶导向配体（DTDL）的系统梳理：靶点组合、化学策略与数据层面观察，适合作为课题背景与靶点对选择参考。 |

### 1.2 双靶对接 / 虚拟筛选方法学

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R03 | Zhou S et al. *Feasibility of Using Molecular Docking-Based Virtual Screening for Searching Dual Target Kinase Inhibitors.* **JCIM** 2013. | [doi:10.1021/ci400065e](https://doi.org/10.1021/ci400065e) | **关键对照文献**：在 CDK2–GSK3B、EGFR–Src、Lck–Src、Lck–VEGFR2 等激酶对上评估对接 VS 找双靶抑制剂。结论：对接能识别单靶抑制剂，但对双靶抑制剂**假阳性高、富集有限**，需联用其他方法。支撑本课题“假阳性放大 / 需任务级融合”论点。 |
| R04 | Jaiteh M et al. *Docking Screens for Dual Inhibitors of Disparate Drug Targets for Parkinson’s Disease.* **J. Med. Chem.** 2018/2019. | [PMC6716773](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6716773/) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/30525601/) | 针对差异较大靶点对（A2A 腺苷受体与 MAO-B）做双靶对接筛选的案例研究，说明异构靶对上“分别对接再筛选”的可行路径与局限。 |
| R05 | Perez-Castillo Y et al. *Toward a Better Scoring Function: Fusing Docking Scoring Functions Improves Virtual Screening for Dual Target Ligands of Parkinson’s Disease.* 2017. | [PMC5725543](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5725543/) | 将多种对接打分函数**融合**用于帕金森相关双靶配体 VS，显示共识/融合打分相对单函数更稳。直接启发本课题“跨靶融合头 + 校准”设计。 |
| R06 | Sivakumar KC et al. *Multitarget approach for natural products targeting inflammation by molecular docking and molecular dynamics.* **Drug Dev. Res.** 2020. | [doi:10.1002/ddr.21673](https://doi.org/10.1002/ddr.21673) | 对接 + 短 MD 联用做多靶天然产物筛选的工作流示例，说明动力学可作姿态/结合稳定性后处理。 |
| R07 | Ferreira LG et al. *Molecular Modeling Techniques Applied to Multitarget Drug Design.* **Curr. Top. Med. Chem.** 2022. | [doi:10.2174/1568026621666211129140958](https://doi.org/10.2174/1568026621666211129140958) | 多靶药物设计中的分子建模技术综述（药效团、对接、QSAR 等），适合写 Related Work 总览。 |

### 1.3 Linked / Bivalent / de novo 专用工具

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R08 | Bai L et al. *TwistDock: Twist-and-Dock for Bivalent Ligand Binding.* **Drug Des. Devel. Ther.** 2019（XIAP Smac mimetics 等）. | [doi:10.2147/DDDT.S194276](https://doi.org/10.2147/DDDT.S194276) · [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6503218/) | **TwistDock**：先固定两端 warhead 对接，再对 linker 单键扭转采样构象系综。适用于**同蛋白双结构域 bivalent**，启发 linked 路线的 `bridge_ok` / 应变指标，但**不是**普通异构双靶通解。 |
| R09 | Yuan Y, Pei J, Lai L. *LigBuilder V3: Multi-Target de novo Drug Design.* **Front. Chem.** 2020. | [doi:10.3389/fchem.2020.00142](https://doi.org/10.3389/fchem.2020.00142) · [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7066513/) | **LigBuilder V3**：多靶 de novo；支持两端片段独立生长后的 **ensemble linking** 与 framework combination。偏设计生成，而非对接打分评价。 |

### 1.4 Merged 设计 + 对接案例

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R10 | *Combining Data-Driven and Structure-Based Approaches in Designing Dual PARP1–BRD4 Inhibitors.* **JCIM** 2024. | [doi:10.1021/acs.jcim.4c01421](https://doi.org/10.1021/acs.jcim.4c01421) | PARP1–BRD4 双抑制剂：数据驱动 + 结构方法；明确主张优先 **merged 公共药效团** 而非简单 linking（linking 易抬高 MW/logP），并用对接优先筛选，有实验双靶活性分子。 |
| R11 | CDK4/6–芳香化酶双抑制剂虚拟筛选案例（*Molecules* 等，2023 前后）. | 见调研文 Related Work；检索关键词 `CDK4/6 aromatase dual inhibitor docking` | 合并药效团/对接的双靶 VS 案例，说明 merged 路线在激酶–酶异构对上的常见写法。 |

### 1.5 评价基准、物理合理性与共识打分（含近年更新）

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R12 | Su M et al. *Comparative Assessment of Scoring Functions: The CASF-2016 Update.* **JCIM** 2019. | [doi:10.1021/acs.jcim.8b00545](https://doi.org/10.1021/acs.jcim.8b00545) · [CASF 主页](http://www.pdbbind.org.cn/casf.php) | **CASF-2016**：scoring / ranking / docking / screening power 四维评测；285 高质量复合物。双靶研究可迁移组件指标，但**尚无**官方“双靶 CASF”。 |
| R12b | Li Y et al. *Assessing protein–ligand interaction scoring functions with the CASF-2013 benchmark.* **Nat. Protoc.** 2018. | [doi:10.1038/nprot.2017.114](https://doi.org/10.1038/nprot.2017.114) | CASF-2013 操作协议：如何按统一脚本评测打分函数；写方法学 Methods 时的标准引用。 |
| R13 | Ericksen SS et al. *Machine Learning Consensus Scoring Improves Performance Across Targets.* **JCIM** 2017. | [doi:10.1021/acs.jcim.7b00153](https://doi.org/10.1021/acs.jcim.7b00153) | 机器学习共识打分提升结构基 VS 跨靶稳健性；可嵌入双靶流水线的单靶侧打分层。 |
| R14 | Consensus docking 综述 / 调查（多程序降低靶点依赖方差）. | [PMC9821981](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9821981/) | 共识对接降低单程序偏差；双靶场景可对每个靶做多引擎共识后再跨靶融合。 |
| R21 | Buttenschoen M et al. *PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences.* **Chem. Sci.** 2024. | [doi:10.1039/D3SC04185A](https://doi.org/10.1039/D3SC04185A) · [arXiv:2308.05777](https://arxiv.org/abs/2308.05777) · [文档](https://posebusters.readthedocs.io/) | **关键现代对接评测**：在 Astex Diverse + 2021 后 PoseBusters Benchmark（~308）上对比 DiffDock/EquiBind/TankBind/Uni-Mol 与 Vina/GOLD。结论：仅看 RMSD 不够；加入**物理合理性（PB-valid）**后，经典对接仍更稳，多数 DL 方法在未见复合物上有效预测少。本课题姿态门控应参考 PB 检查。 |
| R22 | Harris C et al. *Benchmarking Generated Poses: How Rational is Structure-based Drug Design with Generative Models?*（PoseCheck）. | [arXiv:2308.07413](https://arxiv.org/abs/2308.07413) · [GitHub](https://github.com/cch1999/posecheck) | **PoseCheck**：评测生成式 SBDD 姿态的应变能、蛋白–配体 clash、相互作用；显示许多生成模型姿态“看起来对”但物理不合理，常需最小化/重对接。与 PoseBusters 互补。 |
| R23 | Morehead A et al. *Assessing the potential of deep learning for protein-ligand docking*（**PoseBench**）. | [arXiv:2405.14108](https://arxiv.org/abs/2405.14108) · [GitHub](https://github.com/BioinfoMachineLearning/PoseBench) | 综合 benchmark：apo→holo、盲对接、多配体；对比 Vina、DiffDock、DynamicBind、NeuralPLexer、AF3/Chai/Boltz 等。结论：co-folding 整体更强，但对新颖姿态仍困难；强调化学特异性与结构精度的权衡。含 DockGen-E 扩展。 |
| R24 | Corso G et al. **DockGen**（随 DiffDock-L / BindingMOAD 推广）：未见口袋泛化测试集. | [Zenodo DockGen](https://zenodo.org/records/10656052) · 见 DiffDock 仓库说明 | 强调训练/测试口袋不重叠的泛化评测；近年 PocketVina、PoseBench、Bento 等均引用。双靶论文写“泛化”时应避免仅用同源口袋。 |
| R25 | *Bento: Benchmarking Classical and AI Docking on Drug Design–Relevant Data.* bioRxiv 2025/2026. | [bioRxiv 10.64898/2025.12.30.696741](https://www.biorxiv.org/content/10.64898/2025.12.30.696741v1) | 在药化更相关子集上对比经典与 AI 对接；强调药物设计相关过滤后，物理感知方法（如 Gnina）在 PB-valid 上仍有优势。 |
| R26 | *Assessing interaction recovery of predicted protein-ligand poses.* **J. Cheminform.** 2025. | [doi:10.1186/s13321-025-01011-6](https://doi.org/10.1186/s13321-025-01011-6) | 在 PoseBusters 上增加 **PLIF（相互作用指纹）恢复率**：RMSD≤2 Å 且 PB-valid 仍可能丢关键相互作用；GOLD 等经典方法在 PLIF 上仍领先多数 ML。启发本课题用相互作用恢复作双端门控。 |

### 1.5b VS 视角对接评测（**VSDS-VD 标准**及同系文献）

> **漏检说明（为何上次未收录 Gu et al. 2025）：**  
> 上一轮检索关键词偏重 `PoseBusters / CASF / DiffDock / docking RMSD`，未覆盖 `virtual screening enrichment`、`TrueDecoy`、`VSDS-VD`、`Nature Machine Intelligence docking`。该文主贡献是 **VS 富集协议 + 物理合理性联评**，而非再做一个 RMSD leaderboard，故被漏掉。以下按该文标准补全：  
> **(i)** 再对接 RMSD；**(ii)** PoseBusters 物理合理性；**(iii)** TrueDecoy（实验活性/非活）上的 EF；**(iv)** RandomDecoy / MassiveDecoy（更接近真实库）上的 EF；**(v)** AI 对接 vs 物理对接 vs AI 重打分；**(vi)** 层级 VS 策略。

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| **R34** | **Gu S, Shen C, Zhang X et al.** *Benchmarking AI-powered docking methods from the perspective of virtual screening.* **Nat. Mach. Intell.** 2025, 7, 509–520. | [doi:10.1038/s42256-025-00993-0](https://doi.org/10.1038/s42256-025-00993-0) · [GitHub VSDS-VD](https://github.com/shukai1997/VSDS-VD) · [Zenodo 数据](https://doi.org/10.5281/zenodo.13684010) · [Zenodo 代码](https://doi.org/10.5281/zenodo.14649209) | **本系列标杆文**。构建 **VSDS-VD**：TrueDecoy / RandomDecoy / MassiveDecoy。评测 4 个 AI 对接（CarsiDock、KarmaDock、DiffDock、FlexPose）、4 个物理对接（Glide、Surflex、rDock、LeDock）、2 个 AI 重打分（RTMScore、EquiScore）。要点：① KarmaDock/CarsiDock 再对接精度常高于物理法；② 物理法 **PoseBusters 结构合理性**明显更好（CarsiDock 短板多在分子间有效性）；③ TrueDecoy 上 Glide 系 EF 最高，RTMScore 重打分有效；④ RandomDecoy（更像真实 VS）上 AI 工具明显优于 Glide；⑤ 配体后处理对构象/VS 帮助弱甚至负；⑥ 提出**层级 VS** 平衡通量与精度。双靶论文的单靶后端评测协议应直接对齐此文。 |
| R35 | Cai H et al. **CarsiDock**: large-scale pre-training docking + screening. **Chem. Sci.** 2024. | [doi:10.1039/D3SC05552C](https://doi.org/10.1039/D3SC05552C) · [GitHub](https://github.com/carbonsilicon-ai/CarsiDock) | VSDS-VD 中精度领先的 AI 对接之一：~9M 复合物预训练 → 距离矩阵 → 几何优化构象；DEKOIS2.0 上与 RTMScore 联用早期识别强于 Glide SP。 |
| R36 | Cao D et al. **EquiScore**: physical prior + data augmentation scoring. **Nat. Mach. Intell.** 2024. | [doi:10.1038/s42256-024-00849-z](https://doi.org/10.1038/s42256-024-00849-z) · [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.06.18.545464) · [GitHub](https://github.com/CAODH/EquiScore) | VSDS-VD 重打分对照；等变异构图 + PDBscreen；在 DEKOIS2.0/DUD-E 未见蛋白上优于多种打分；对不同对接姿态的重打分稳健。 |
| R37 | Dong T et al. **FlexPose**: equivariant flexible protein–ligand pose modeling. **J. Chem. Theory Comput.** 2023. | [doi:10.1021/acs.jctc.3c00273](https://doi.org/10.1021/acs.jctc.3c00273) · [GitHub](https://github.com/tiejundong/FlexPose) | VSDS-VD 纳入的 AI 对接之一；几何深度学习直接柔性建模结合姿态。 |
| R38 | Zhang X et al. *Advancing ligand docking through deep learning: challenges and prospects in virtual screening.* **Acc. Chem. Res.** 2024. | [doi:10.1021/acs.accounts.4c00093](https://doi.org/10.1021/acs.accounts.4c00093) | 侯廷军组观点文：DL 对接在 VS 中的挑战（泛化、物理合理性、打分与排序脱节）与前景；与 VSDS-VD 叙事一致。 |
| R39 | *SCORCH2*: heterogeneous consensus for high-enrichment interaction-based VS. **Adv. Sci.** 2025. | [doi:10.1002/advs.202508318](https://doi.org/10.1002/advs.202508318) · [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.31.646332) | 在 **VSDS-vd TrueDecoy** 等上评测的共识重打分；强调相互作用特征与早期富集；与 RTMScore 并列提升 EF。 |
| R40 | *PoseX*: self-docking + **cross-docking** open benchmark（23 方法）. | [arXiv:2505.01700](https://arxiv.org/abs/2505.01700) · [GitHub](https://github.com/CataAI/PoseX) | 718 self + 1312 cross；物理 / AI docking / co-folding；强调松弛后处理可大幅消除 AI clash；口袋指定显著提升。补强 VSDS-VD 未充分展开的**交叉对接**维度。 |
| R41 | *UniDock-Pro*: GPU 高通量 SBVS/LBVS/hybrid VS. **JCIM** 2025/2026. | [doi:10.1021/acs.jcim.5c02587](https://doi.org/10.1021/acs.jcim.5c02587) | 在 DUDE-Z 与 **VSDS-vd TrueDecoy** 上报告早期富集；代表工程化高通量 VS 平台对照。 |
| R42 | Mysinger et al. **DUD-E**；Bauer et al. **DEKOIS 2.0**；Tran-Nguyen et al. **LIT-PCBA**；Stein et al. property-unmatched decoys. | [DUD-E 10.1021/jm300687e](https://doi.org/10.1021/jm300687e) · [DEKOIS2.0 10.1021/ci400115b](https://doi.org/10.1021/ci400115b) · [LIT-PCBA 10.1021/acs.jcim.0c00155](https://doi.org/10.1021/acs.jcim.0c00155) · [property-unmatched 10.1021/acs.jcim.0c00598](https://doi.org/10.1021/acs.jcim.0c00598) | VSDS-VD 之前的经典 VS 基准；注意 **DUD-E 隐藏偏差**（Chen et al. *PLoS ONE* 2019, [10.1371/journal.pone.0220113](https://doi.org/10.1371/journal.pone.0220113)）与 Sieg et al. bias control（[10.1021/acs.jcim.8b00712](https://doi.org/10.1021/acs.jcim.8b00712)）。TrueDecoy/RandomDecoy 正是为更贴近真实 VS、降低误导性 EF 而设计。 |

**按 VSDS-VD 标准写本课题评测时的最小清单：**

| 维度 | 指标 | 说明 |
|------|------|------|
| 姿态精度 | RMSD ≤ 2 Å success | 双靶则两端分别报告 |
| 物理合理性 | PoseBusters PB-valid | AI 对接必报；不可只报 RMSD |
| 严格 VS | TrueDecoy 类：实验活性 vs 实验非活 | EF@0.5/1/5%、AUROC/BEDROC |
| 真实库 VS | RandomDecoy / 商用库随机负例 | 与 TrueDecoy 结论可能**相反**（VSDS-VD 已证明） |
| 协议消融 | 对接原生分 vs RTMScore/EquiScore 重打分 | 采样与打分解耦 |
| 策略 | 层级筛选（快 AI → 物理/重打分精筛） | 大规模库实用路径 |

### 1.5c 侯廷军组系谱 + 同档高质量期刊对照

> 检索策略：以 Hou/Kang/Shen/Zhang（浙大药学院 + CarbonSilicon）的 **AI 对接 / 打分 / VS 评测 / 层级筛选** 工作为中心，再补同档期刊（*Nat. Mach. Intell.*、*Nat. Methods*、*Nat. Comput. Sci.*、*Chem. Sci.*、*JCIM*、*JCTC*、*J. Cheminform.*、*Acc. Chem. Res.*、*TiPS*、*APSB*）上主题相近的竞品与评测文。

#### A. 侯组 / CarbonSilicon 核心产出（按主题）

| # | 文献 | 期刊 / 链接 | 内容简介 |
|---|------|-------------|----------|
| H01 | Gu et al. **VSDS-VD**（见 R34） | **Nat. Mach. Intell.** 2025 · [10.1038/s42256-025-00993-0](https://doi.org/10.1038/s42256-025-00993-0) | VS 视角系统评测 AI vs 物理对接；TrueDecoy/RandomDecoy；层级 VS。 |
| H02 | Zhang et al. **KarmaDock**（见 R28） | **Nat. Comput. Sci.** 2023 · [10.1038/s43588-023-00511-5](https://doi.org/10.1038/s43588-023-00511-5) | 超大库 AI 对接：EGNN 姿态 + MDN 打分；实验验证 LTK 抑制剂。 |
| H03 | Cai et al. **CarsiDock**（见 R35） | **Chem. Sci.** 2024 · [10.1039/D3SC05552C](https://doi.org/10.1039/D3SC05552C) | ~9M 预训练距离矩阵对接；与 RTMScore 联用 VS 早期识别强。 |
| H04 | Shen et al. **RTMScore**（见 R20） | **J. Med. Chem.** 2022 · [10.1021/acs.jmedchem.2c00991](https://doi.org/10.1021/acs.jmedchem.2c00991) | 残基–原子距离似然 + Graph Transformer；CASF docking/screening 强。 |
| H05 | Shen et al. **GenScore**（见 R20b） | **Chem. Sci.** 2023 · [10.1039/D3SC02044D](https://doi.org/10.1039/D3SC02044D) | 可调亲和项使 scoring/ranking/docking/screening 更均衡。 |
| H06 | Zhang et al. DL docking 观点文（见 R38） | **Acc. Chem. Res.** 2024 · [10.1021/acs.accounts.4c00093](https://doi.org/10.1021/acs.accounts.4c00093) | 评指标、应用场景、物理合理性、生成 vs 回归；VS 实用导向。 |
| H07 | Zhang, Shen, Hsieh, Hou. *Harnessing deep learning for enhanced ligand docking.* | **Trends Pharmacol. Sci.** 2024 · [10.1016/j.tips.2023.12.004](https://doi.org/10.1016/j.tips.2023.12.004) | 短评：DLLD 相对传统 search–score 范式的机会与局限。 |
| H08 | Shen et al. **DrugFlow**：一站式 AI 药物发现平台. | **JCIM** 2024 · [10.1021/acs.jcim.4c00621](https://doi.org/10.1021/acs.jcim.4c00621) · [drugflow.com](https://www.drugflow.com) | 对接 / QSAR / 生成 / ADMET / VS 流水线工程化；体现「层级筛选」落地。 |
| H09 | Zhang et al. **TB-IECS**：理论能量项 + XGBoost VS 打分. | **J. Cheminform.** 2023 · [10.1186/s13321-023-00731-x](https://doi.org/10.1186/s13321-023-00731-x) | Smina/NNScore2 能量组分组合；DUD-E / LIT-PCBA / ChemDiv 靶点集上优于/对标 Glide SP。 |
| H10 | Cai, Shen, Hou et al. **CarsiDock-Cov**：共价对接与筛选. | **Acta Pharm. Sin. B** 2025 · [10.1016/j.apsb.2025.07.043](https://doi.org/10.1016/j.apsb.2025.07.043) | 将 CarsiDock 范式扩展到共价配体；公开共价对接/筛选基准上有竞争力。 |

#### B. 同档高质量对照（竞品 / 评测 / 高通量）

| # | 文献 | 期刊 / 链接 | 与侯组关系 |
|---|------|-------------|------------|
| H11 | Cao, Zheng et al. **EquiScore**（见 R36） | **Nat. Mach. Intell.** 2024 · [10.1038/s42256-024-00849-z](https://doi.org/10.1038/s42256-024-00849-z) | 中科院上海药物所；VSDS-VD 重打分对照；物理先验 + 数据增强。 |
| H12 | Cao, Zheng et al. **SurfDock**（见 R29） | **Nat. Methods** 2025 · [10.1038/s41592-024-02516-y](https://doi.org/10.1038/s41592-024-02516-y) | 表面引导扩散对接；强调物理约束与真实 VS 发现案例（ALDH1B1）。 |
| H13 | Li, Cao et al. *Decoding the limits of deep learning in molecular docking…* | **Chem. Sci.** 2025 · [10.1039/D5SC05395A](https://doi.org/10.1039/D5SC05395A) | 五维系统对比：姿态精度、物理合理性、相互作用恢复、VS、泛化；扩散/回归/混合优缺点与失败机制。与 VSDS-VD / PoseBusters 互补的「诊断型」综述评测。 |
| H14 | Corso et al. **DiffDock-L / DockGen** 泛化 | **ICLR / arXiv** · [2402.18396](https://arxiv.org/abs/2402.18396) · DiffDock [2210.01776](https://arxiv.org/abs/2210.01776) | 未见口袋泛化；VSDS-VD 与 PoseBench 常作 AI 对接基线。 |
| H15 | Yu et al. **Uni-Dock**：GPU 超大库物理对接. | **JCTC** 2023 · [10.1021/acs.jctc.2c01145](https://doi.org/10.1021/acs.jctc.2c01145) · [GitHub](https://github.com/dptech-corp/Uni-Dock) | 千倍加速 Vina 系；与 KarmaDock「AI 加速」并列的高通量物理路径；UniDock-Pro（R41）为其平台扩展。 |
| H16 | Moon et al. **PIGNet2**：物理知情 GNN 打分 + VS. | **Digit. Discov.** 2024 · [10.1039/D3DD00149K](https://doi.org/10.1039/D3DD00149K) | 与 GenScore/RTMScore 同台比较；强调多任务均衡与数据增强。 |
| H17 | Xia / Hou 系相关：**DiffDock-NMDN** 端到端盲对接 + VS. | **JCIM** 2024 · [10.1021/acs.jcim.4c01014](https://doi.org/10.1021/acs.jcim.4c01014) | DiffDock 采样 + 归一化距离似然打分（NMDN）；LIT-PCBA 上报告 EF；与「采样–打分解耦」路线一致。 |
| H18 | Méndez-Lucio et al. **DeepDock**（距离似然早期代表）. | **Nat. Mach. Intell.** 2021 · [10.1038/s42256-021-00409-9](https://doi.org/10.1038/s42256-021-00409-9) | RTMScore/GenScore 方法论前身之一；几何深度学习预测结合构象。 |
| H19 | Buttenschoen et al. **PoseBusters**（见 R21） | **Chem. Sci.** 2024 | 物理有效性评测标杆；VSDS-VD 再对接合理性分析直接依赖。 |
| H20 | McNutt / Koes **GNINA**（见 R19/R19b） | **J. Cheminform.** 2021/2025 | 「经典采样 + DL 重打分」工程范式；与 RTMScore 重打分哲学相近、实现不同。 |

#### C. 读这些文章时的对照表（写 Related Work 用）

| 问题 | 优先读 |
|------|--------|
| AI 对接准不准、物理是否合理、VS 富集如何 | **H01 VSDS-VD** + **H13 Decoding limits** + **H19 PoseBusters** |
| 要可插拔重打分 | **H04 RTMScore** / **H05 GenScore** / **H11 EquiScore** / **H20 GNINA** |
| 要超大库通量 | **H02 KarmaDock**（AI）或 **H15 Uni-Dock**（GPU 物理）→ **H08 DrugFlow** 流水线 |
| 要精度优先姿态 | **H03 CarsiDock** / **H12 SurfDock** / DiffDock-L |
| 方法叙事与缺口 | **H06 Acc. Chem. Res.** + **H07 TiPS** |
| 共价扩展 | **H10 CarsiDock-Cov** |

### 1.6 生成式双靶（前沿算法）

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R15 | DualDiff / CompDiff. *Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design.* **NeurIPS** 2024. | [arXiv:2410.20688](https://arxiv.org/abs/2410.20688) | 将单靶扩散模型重编程到双靶：口袋对齐策略（center / RMSD-anchor / score-anchor）；评估用 P1/P2 Vina、**Max Vina**、Dual High Affinity、双姿态 RMSD 等。偏生成，指标可借鉴。 |
| R16 | FuseDiff. *Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design.* | [arXiv:2603.05567](https://arxiv.org/abs/2603.05567) | 对称保持的联合扩散双靶 SBDD；与 DualDiff 同属生成式前沿，非本课题首选落地路径。 |

### 1.7 新对接算法与单靶打分（引擎候选 / 对照，非双靶专用）

| # | 文献 / 工具 | 链接 | 内容简介 |
|---|-------------|------|----------|
| R17 | Zheng L et al. *OnionNet-2: … Residue-Atom Contacting Shells.* **Front. Chem.** 2021. | [doi:10.3389/fchem.2021.753002](https://doi.org/10.3389/fchem.2021.753002) | 残基–原子多层壳接触 CNN 预测亲和力；CASF-2016 上 Pearson r≈0.86。可作单靶侧亲和 scorer，**不是**双靶融合方法。 |
| R18 | Wang Z et al. *DeepRMSD+Vina*：可微姿态优化框架. **Brief. Bioinform.** 2023. | [doi:10.1093/bib/bbac520](https://doi.org/10.1093/bib/bbac520) · [arXiv:2206.13345](https://arxiv.org/abs/2206.13345) | MLP 预测姿态 RMSD 并与 Vina 混合；CASF-2016 docking power 很高。适合作**姿态质量门控/重打分**，非主创新点。 |
| R19 | McNutt A et al. *GNINA 1.0: molecular docking with deep learning.* **J. Cheminform.** 2021. | [doi:10.1186/s13321-021-00522-2](https://doi.org/10.1186/s13321-021-00522-2) · [GitHub](https://github.com/gnina/gnina) | Vina/Smina 采样 + CNN 姿态/亲和重打分；再对接与交叉对接均优于纯 Vina。**推荐可插拔引擎。** |
| R19b | McNutt A et al. *GNINA 1.3: the next increment…* **J. Cheminform.** 2025. | [doi:10.1186/s13321-025-00973-x](https://doi.org/10.1186/s13321-025-00973-x) · [PMC11874439](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11874439/) | GNINA 后续版本与默认 CNN 集成更新；CACHE 挑战等应用见配套文。 |
| R20 | Shen C et al. **RTMScore**：residue–atom 距离似然 + Graph Transformer. **J. Med. Chem.** 2022. | [doi:10.1021/acs.jmedchem.2c00991](https://doi.org/10.1021/acs.jmedchem.2c00991) · [GitHub](https://github.com/sc8668/RTMScore) | CASF-2016 上 docking/screening power 突出；可作交叉对接与大规模 VS 重打分。本课题推荐单靶侧对照 scorer。 |
| R20b | Shen C et al. **GenScore**：广义蛋白–配体打分框架（RTMScore 扩展）. **Chem. Sci.** 2023. | [doi:10.1039/D3SC02044D](https://doi.org/10.1039/D3SC02044D) · [GitHub](https://github.com/sc8668/GenScore) | 在 scoring/ranking/docking/screening 间更均衡；适合与 RTMScore 并列作消融。 |
| R27 | Corso G et al. *DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking.* **ICLR** 2023. | [arXiv:2210.01776](https://arxiv.org/abs/2210.01776) · [GitHub](https://github.com/gcorso/DiffDock) | 将对接建模为配体姿态流形上的扩散生成；PDBBind 上 top-1 RMSD\<2 Å 显著高于当时回归式 DL 与部分经典法；提供置信度。后续 PoseBusters 显示物理有效性仍不足，宜与经典/重打分配套。 |
| R28 | Zhang X et al. *KarmaDock: … ultra-large library docking.* **Nat. Comput. Sci.** 2023. | [doi:10.1038/s43588-023-00511-5](https://doi.org/10.1038/s43588-023-00511-5) · [GitHub](https://github.com/schrojunzhang/KarmaDock) | EGNN 姿态生成/校正 + MDN 结合强度估计；面向超大库 VS，有实验验证案例。可作吞吐对照。 |
| R29 | Cao D et al. *SurfDock is a surface-informed diffusion generative model…* **Nat. Methods** 2025. | [doi:10.1038/s41592-024-02516-y](https://doi.org/10.1038/s41592-024-02516-y) | 表面信息引导的扩散对接/复合物预测；代表最新生成式对接方向之一。 |
| R30 | Wang Z et al. **IGModel**：几何 GNN 同时预测姿态 RMSD 与 pKd. **Brief. Bioinform.** 2024. | [doi:10.1093/bib/bbae145](https://doi.org/10.1093/bib/bbae145) · [bioRxiv](https://doi.org/10.1101/2023.11.01.565115) · [GitHub](https://github.com/zchwang/IGModel) | CASF-2016 docking power 极高；交叉对接与 AF2 结构上表现稳健。适合姿态排序/门控。 |
| R31 | *Integrating ML-Based Pose Sampling with Established Scoring Functions for VS.* **JCIM** 2025. | [doi:10.1021/acs.jcim.5c00380](https://doi.org/10.1021/acs.jcim.5c00380) · [PMC12117556](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12117556/) | **DiffDock-L 采样 + Vina/GNINA/RTMScore 打分** 在 DUDE-Z 上的 VS 评估：采样与打分解耦；打分函数选择强烈影响 VS 成败。与本课题“引擎可插拔 + 融合头”高度同构（单靶版）。 |
| R32 | *CompassDock / Compass*：PoseCheck + AA-Score 综合评估与微调. | [arXiv:2406.06841](https://arxiv.org/abs/2406.06841) | 将物理/化学检查与亲和经验分统一进 DiffDock 推理与微调；强调仅 RMSD 不够。 |
| R33 | *PocketVina*：多口袋条件 + 搜索式物理有效对接. | [arXiv:2506.20043](https://arxiv.org/abs/2506.20043) | 口袋预测 + 系统多口袋搜索；在 PoseBusters/DockGen/Astex 等上强调 **PB-valid** 成功率；无需任务特异训练，适合大规模筛选对照。 |

### 1.8 调研流程工具（非学术论文）

| # | 资源 | 链接 | 说明 |
|---|------|------|------|
| T01 | PaperSpine | [github.com/WUBING2023/PaperSpine](https://github.com/WUBING2023/PaperSpine) | 研究问题结构化、SOTA gap、证据分层 |
| T02 | ACS Writer v2 | [github.com/Caosmart1979/acs-writer-v2](https://github.com/Caosmart1979/acs-writer-v2) | 文献探索与方法学评价工作流 |

### 1.9 对本课题的直接启示（对接算法 / benchmark 更新后）

1. **单靶引擎选型**：优先 Vina / **GNINA** / **RTMScore（重打分）**；DiffDock-L / CarsiDock / KarmaDock 可作采样对照，但必须过 PoseBusters / PLIF 门控。  
2. **评价协议（对齐 VSDS-VD）**：RMSD + **PB-valid** + TrueDecoy 式 EF + RandomDecoy 式 EF；双靶再叠加 dual-vs-single 与短板敏感融合。**不要只在一种 decoy 设定上宣称 VS 优势。**  
3. **不要**把再训通用 pose scorer 当主创新；CASF / PoseBusters / VSDS-VD / DUDE-Z 上单靶打分已高度内卷。  
4. **可发表空白仍在任务层**：跨靶校准、阈值-边距 softmin、硬负样本、类型条件融合——可用上述新引擎作可插拔后端。  
5. **层级 VS**：大规模筛选可走「快 AI 对接 → RTMScore/EquiScore 重打分 → 物理法精筛」；双靶版应对两端分别过阈值后再融合。

---

## 二、双靶分子 / 共晶相关文献（按分子体系）

下列文献同时支撑 **第三节分子目录** 中的条目。

| # | 体系 | 文献 | 链接 | 内容简介 |
|---|------|------|------|----------|
| M01 | Mcl-1 / Bcl-xL 杂交抑制剂 | Tanaka Y et al. *Discovery of Potent Mcl-1/Bcl-xL Dual Inhibitors…* **J. Med. Chem.** 2013. | [doi:10.1021/jm401170c](https://doi.org/10.1021/jm401170c) | 通过杂交/拴系得到 Mcl-1 与 Bcl-xL 双抑制剂；**compound 10** 在两端均有共晶（3WIY / 3WIZ，配体 LC6）。优化类似物 11：Mcl-1 IC50 ≈ 0.088 µM，Bcl-xL IC50 ≈ 0.0037 µM。**Tier A + linked 金标准。** |
| M02 | LpxA / LpxD 双结合配体 | *Scientific Reports* 2019（铜绿假单胞菌 LpxA/LpxD 与 Q5M 等）. | [doi:10.1038/s41598-019-51844-z](https://doi.org/10.1038/s41598-019-51844-z) | 底物模拟型酰基口袋配体；同一化学组分 **Q5M** 分别与 LpxD（6UEC）、LpxA（6UEE/6UEG）成复合物。**Tier A + merged。** |
| M03 | PknA / PknB 双抑制剂 | PDB 主引文：6B2P / 6B2Q 系列（结核分枝杆菌丝氨酸/苏氨酸激酶）. | [6B2P](https://www.rcsb.org/structure/6B2P) · [6B2Q](https://www.rcsb.org/structure/6B2Q) | 同一配体 **CJJ** 分别与 PknB、PknA 结晶；系列 Ki 可达约 nM 级双抑制。**Tier A + merged。** |
| M04 | EGFR / HER2 — TAK-285 | Ishikawa T et al. / JBC 结构与药化工作（TAK-285）. | [doi:10.1074/jbc.M110.209817](https://doi.org/10.1074/jbc.M110.209817) | 临床阶段双 EGFR/HER2 TKI **TAK-285**：EGFR [3POZ](https://www.rcsb.org/structure/3POZ)、HER2 [3RCD](https://www.rcsb.org/structure/3RCD)；相关突变体见 3W2O。**Tier A + merged。** |
| M05 | BET / HDAC — 4-酰基吡咯异羟肟酸杂合体 | *J. Med. Chem.* 2021（LSH 系列 hybrids）. | [doi:10.1021/acs.jmedchem.1c00733](https://doi.org/10.1021/acs.jmedchem.1c00733) | BET–HDAC 双靶杂合体；**7AXR** 为 BRD4 BD1 共晶（配体 S7T），HDAC 端未见对等沉积。**Tier B + linked。** |
| M06 | BET / HDAC — NB 系列 merged | *ACS Chem. Biol.* 2023（NanoBRET 等细胞结合）. | [doi:10.1021/acschembio.3c00427](https://doi.org/10.1021/acschembio.3c00427) | Merged 型 BET / I 类 HDAC 双靶系列；PDB **8P9F–8P9L** 均为 BRD4 端；优化类似物细胞结合可达约 100 nM 级。**Tier B + merged。** |
| M07 | ER / 碳酸酐酶（CA）双调节剂 | *Eur. J. Med. Chem.* 2022. | [doi:10.1016/j.ejmech.2022.115011](https://doi.org/10.1016/j.ejmech.2022.115011) | 首批双 CA/ER 调节剂报道；ERα 复合物 [8EV1](https://www.rcsb.org/structure/8EV1)、[8EV2](https://www.rcsb.org/structure/8EV2)。hCA 端 PDB 待核对，暂 **Tier B**。 |
| M08 | MurD / MurE 双抑制剂 | PDB [2Y1O](https://www.rcsb.org/structure/2Y1O) 主引文（噻唑烷类 T26）. | [2Y1O](https://www.rcsb.org/structure/2Y1O) | 设计为 Mur 连接酶双抑制剂；种子集中仅有 MurD 复合物。**Tier B。** |
| M09 | PD-L1 / VISTA — P17 | PDB [9INU](https://www.rcsb.org/structure/9INU) 主引文. | [9INU](https://www.rcsb.org/structure/9INU) | 双通路免疫小分子 P17；与 PD-L1 成晶，VISTA 端未在种子集。**Tier B。** |
| M10 | GyrB / ParE（Topo IV）抗菌系列 | PDB [4GEE](https://www.rcsb.org/structure/4GEE)、[4HY1](https://www.rcsb.org/structure/4HY1) 同系列主引文. | [4GEE](https://www.rcsb.org/structure/4GEE) · [4HY1](https://www.rcsb.org/structure/4HY1) | 双 GyrB/TopoIV 抗菌项目；两端配体**不完全相同**，仅作系列级交叉对接压力测试。**Tier C。** |

---

## 三、已搜集双靶分子信息一览

分类约定：

- **Tier A**：同一双靶配体两端均有共晶 → 双端姿态 RMSD 金标准  
- **Tier B**：已发表双靶分子，但目前仅一端有 PDB  
- **Tier C**：同系列不同配体分别结晶在两靶  
- **Tier D**：谱系说明用，不作双靶姿态金标准  

设计类型：`linked` / `fused` / `merged`（Morphy）。

### 3.1 Tier A（两端共晶）

#### DT-A-001 — Mcl-1 / Bcl-xL · compound 10（LC6）

| 字段 | 内容 |
|------|------|
| 配体 | compound_10_Tanaka2014；chem_comp **LC6** |
| 靶点 A/B | Mcl-1 (UniProt Q07820) / Bcl-xL (Q07817) |
| PDB | [3WIY](https://www.rcsb.org/structure/3WIY) (2.15 Å) / [3WIZ](https://www.rcsb.org/structure/3WIZ) (2.45 Å) |
| 设计类型 | **linked**（杂交/拴系） |
| 活性 | 优化类似物 11：Mcl-1 IC50 0.088 µM；Bcl-xL IC50 0.0037 µM |
| 文献 | [10.1021/jm401170c](https://doi.org/10.1021/jm401170c) |
| 评估用途 | `pose_both_ends` 金标准 |

#### DT-A-002 — LpxD / LpxA · Q5M

| 字段 | 内容 |
|------|------|
| 配体 | compound_1_Lpx_dual；**Q5M** = 4-(naphthalen-1-yl)-4-oxobutanoic acid |
| 靶点 | LpxD_Pa (Q9HXY6) / LpxA_Pa (Q9HVF4) |
| PDB | [6UEC](https://www.rcsb.org/structure/6UEC) (2.60 Å) / [6UEE](https://www.rcsb.org/structure/6UEE) (2.10 Å) |
| 设计类型 | **merged** |
| 活性/备注 | SPR 双结合；底物模拟酰基口袋 |
| 文献 | [10.1038/s41598-019-51844-z](https://doi.org/10.1038/s41598-019-51844-z) |

#### DT-A-003 — PknB / PknA · CJJ

| 字段 | 内容 |
|------|------|
| 配体 | **CJJ**（氯代吡唑氨基嘧啶–噻吩磺酰胺类） |
| 靶点 | PknB_Mtb (P0A5Z1) / PknA_Mtb (P0A5Z0) |
| PDB | [6B2P](https://www.rcsb.org/structure/6B2P) (3.01 Å) / [6B2Q](https://www.rcsb.org/structure/6B2Q) (2.88 Å) |
| 设计类型 | **merged** |
| 活性 | 系列约 Ki ≈ 5 nM 级双 PknA/PknB |
| 文献 | PDB 主引文 6B2P/6B2Q |

#### DT-A-004 — EGFR / HER2 · TAK-285

| 字段 | 内容 |
|------|------|
| 配体 | **TAK-285**（chem_comp 03Q） |
| 靶点 | EGFR (P00533) / HER2 (P04626) |
| PDB | [3POZ](https://www.rcsb.org/structure/3POZ) (1.50 Å) / [3RCD](https://www.rcsb.org/structure/3RCD) (3.21 Å)；相关 [3W2O](https://www.rcsb.org/structure/3W2O) |
| 设计类型 | **merged** |
| 活性 | 临床阶段双 EGFR/HER2 TKI |
| 文献 | [10.1074/jbc.M110.209817](https://doi.org/10.1074/jbc.M110.209817) |

#### DT-A-005 — LpxA / LpxD · Q5M（重复结构对）

| 字段 | 内容 |
|------|------|
| 说明 | 与 DT-A-002 同化学型；LpxA 用 [6UEG](https://www.rcsb.org/structure/6UEG) (2.00 Å) 配对 LpxD [6UEC](https://www.rcsb.org/structure/6UEC) |
| 用途 | **结构重复对**，统计独立体系时勿与 A-002 重复计数 |

### 3.2 Tier B（仅一端共晶）

#### BET–HDAC linked：DT-B-001

| case_id | 配体 | PDB（仅 BRD4） | 第二靶 | 类型 | 文献 |
|---------|------|----------------|--------|------|------|
| DT-B-001 | LSH24 / S7T | [7AXR](https://www.rcsb.org/structure/7AXR) 1.50 Å | HDAC1–3/6（无 PDB） | linked | [10.1021/acs.jmedchem.1c00733](https://doi.org/10.1021/acs.jmedchem.1c00733) |

亚微摩尔级 BET 与 HDAC 抑制（文中 hybrids 49/61 等相关）。

#### BET–HDAC merged 系列：DT-B-002–008

| case_id | 配体名 | PDB | chem_comp / 备注 |
|---------|--------|-----|------------------|
| DT-B-002 | NB161 | [8P9F](https://www.rcsb.org/structure/8P9F) | 系列起点 |
| DT-B-003 | NB390 | [8P9G](https://www.rcsb.org/structure/8P9G) | 同系列 |
| DT-B-004 | NB437 | [8P9H](https://www.rcsb.org/structure/8P9H) | 同系列 |
| DT-B-005 | NB462 | [8P9I](https://www.rcsb.org/structure/8P9I) | 同系列 |
| DT-B-006 | NB500 | [8P9J](https://www.rcsb.org/structure/8P9J) | 同系列 |
| DT-B-007 | NB503 | [8P9K](https://www.rcsb.org/structure/8P9K) 1.25 Å；**X9U** | lead；高分辨 |
| DT-B-008 | NB512 | [8P9L](https://www.rcsb.org/structure/8P9L) | 同系列 |

- 第二靶：HDAC1/2（Q13547 等），**无双配体 HDAC 共晶**  
- 文献：[10.1021/acschembio.3c00427](https://doi.org/10.1021/acschembio.3c00427)  
- 活性：优化类似物细胞 NanoBRET ~100 nM 级  

#### ER / CA：DT-B-009–010

| case_id | 配体 | ER PDB | CA 端 | 文献 |
|---------|------|--------|-------|------|
| DT-B-009 | WVE / LYQ 等 | [8EV2](https://www.rcsb.org/structure/8EV2) 2.01 Å | 待确认 | [10.1016/j.ejmech.2022.115011](https://doi.org/10.1016/j.ejmech.2022.115011) |
| DT-B-010 | 同研究相关配体 | [8EV1](https://www.rcsb.org/structure/8EV1) 1.83 Å | 待确认 | 同上 |

确认 CA PDB 后方可升为 Tier A。

#### 其他 Tier B

| case_id | 名称 | 靶点对 | PDB | 类型 | 备注 |
|---------|------|--------|-----|------|------|
| DT-B-011 | T26 | MurD / MurE | [2Y1O](https://www.rcsb.org/structure/2Y1O) 1.49 Å（仅 MurD） | merged | 双 Mur 连接酶设计 |
| DT-B-012 | P17 | PD-L1 / VISTA | [9INU](https://www.rcsb.org/structure/9INU) 2.70 Å（仅 PD-L1） | merged | 双通路免疫；体内优于单靶（文内） |

### 3.3 Tier C（系列相关，非同一配体）

| case_id | 体系 | PDB A / B | 说明 |
|---------|------|-----------|------|
| DT-C-001 | GyrB / ParE 吡咯并嘧啶抗菌系列 | [4GEE](https://www.rcsb.org/structure/4GEE) / [4HY1](https://www.rcsb.org/structure/4HY1) | 两端配体不同 → 仅 `crossdock_stress`，不可做同一配体双端 RMSD |

### 3.4 Tier D（谱系注释）

| case_id | 名称 | PDB | 说明 |
|---------|------|-----|------|
| DT-NOTE-001 | compound 4（Mcl-1 前体） | [3WIX](https://www.rcsb.org/structure/3WIX) | 杂交设计前的 Mcl-1 选择性前体；**不作双靶姿态金标准**；文献同 M01 |

### 3.5 目录统计（v0.1）

| Tier | 条目数 | 独立化学体系（约） |
|------|--------|-------------------|
| A | 5（含 1 条 Lpx 重复对） | ~4（Mcl-1/Bcl-xL；Lpx；Pkn；EGFR/HER2） |
| B | 12 | ~5（BET–HDAC linked；BET–HDAC merged；ER/CA；Mur；PD-L1/VISTA） |
| C | 1 | 1（GyrB/ParE） |
| D | 1 | 谱系用 |
| **合计** | **19** | — |

> 注意：Lpx 多 PDB、BET–HDAC 多类似物会**膨胀行数**；报告成功率时请按**独立靶点对/化学型**去重。

---

## 四、快速链接索引（仅 URL）

### 方法学 DOI / PMC / arXiv

- https://doi.org/10.1021/jm0603015  
- https://doi.org/10.1021/acs.jmedchem.4c00838  
- https://doi.org/10.1021/ci400065e  
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6716773/  
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5725543/  
- https://doi.org/10.1002/ddr.21673  
- https://doi.org/10.2174/1568026621666211129140958  
- https://doi.org/10.2147/DDDT.S194276  
- https://doi.org/10.3389/fchem.2020.00142  
- https://doi.org/10.1021/acs.jcim.4c01421  
- https://doi.org/10.1021/acs.jcim.8b00545  
- https://doi.org/10.1038/nprot.2017.114  
- https://doi.org/10.1021/acs.jcim.7b00153  
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9821981/  
- https://doi.org/10.1039/D3SC04185A  
- https://arxiv.org/abs/2308.05777  
- https://arxiv.org/abs/2308.07413  
- https://arxiv.org/abs/2405.14108  
- https://zenodo.org/records/10656052  
- https://www.biorxiv.org/content/10.64898/2025.12.30.696741v1  
- https://doi.org/10.1186/s13321-025-01011-6  
- https://doi.org/10.1038/s42256-025-00993-0  
- https://github.com/shukai1997/VSDS-VD  
- https://doi.org/10.5281/zenodo.13684010  
- https://doi.org/10.5281/zenodo.14649209  
- https://doi.org/10.1039/D3SC05552C  
- https://doi.org/10.1038/s42256-024-00849-z  
- https://doi.org/10.1021/acs.jctc.3c00273  
- https://doi.org/10.1021/acs.accounts.4c00093  
- https://doi.org/10.1002/advs.202508318  
- https://arxiv.org/abs/2505.01700  
- https://doi.org/10.1021/acs.jcim.5c02587  
- https://doi.org/10.1021/jm300687e  
- https://doi.org/10.1021/ci400115b  
- https://doi.org/10.1021/acs.jcim.0c00155  
- https://doi.org/10.1371/journal.pone.0220113  
- https://arxiv.org/abs/2410.20688  
- https://arxiv.org/abs/2603.05567  

### 对接算法 / 打分（近年补充）

- https://doi.org/10.3389/fchem.2021.753002  
- https://doi.org/10.1093/bib/bbac520  
- https://arxiv.org/abs/2206.13345  
- https://doi.org/10.1186/s13321-021-00522-2  
- https://doi.org/10.1186/s13321-025-00973-x  
- https://doi.org/10.1021/acs.jmedchem.2c00991  
- https://doi.org/10.1039/D3SC02044D  
- https://arxiv.org/abs/2210.01776  
- https://doi.org/10.1038/s43588-023-00511-5  
- https://doi.org/10.1038/s41592-024-02516-y  
- https://doi.org/10.1093/bib/bbae145  
- https://doi.org/10.1101/2023.11.01.565115  
- https://doi.org/10.1021/acs.jcim.5c00380  
- https://arxiv.org/abs/2406.06841  
- https://arxiv.org/abs/2506.20043  
- https://github.com/gnina/gnina  
- https://github.com/gcorso/DiffDock  
- https://github.com/sc8668/RTMScore  
- https://github.com/sc8668/GenScore  
- https://github.com/schrojunzhang/KarmaDock  
- https://github.com/carbonsilicon-ai/CarsiDock  
- https://github.com/CAODH/EquiScore  
- https://github.com/tiejundong/FlexPose  
- https://github.com/CataAI/PoseX  
- https://github.com/BioinfoMachineLearning/PoseBench  
- https://posebusters.readthedocs.io/  

### 分子/结构文献 DOI

- https://doi.org/10.1021/jm401170c  
- https://doi.org/10.1038/s41598-019-51844-z  
- https://doi.org/10.1074/jbc.M110.209817  
- https://doi.org/10.1021/acs.jmedchem.1c00733  
- https://doi.org/10.1021/acschembio.3c00427  
- https://doi.org/10.1016/j.ejmech.2022.115011  

### 种子集 PDB（RCSB）

3WIY, 3WIZ, 3WIX, 6UEC, 6UEE, 6UEG, 6B2P, 6B2Q, 3POZ, 3RCD, 3W2O, 7AXR, 8P9F–8P9L, 8EV1, 8EV2, 2Y1O, 9INU, 4GEE, 4HY1  

统一入口：`https://www.rcsb.org/structure/<PDBID>`

---

## 五、维护说明

- 新增文献：在第一节按主题追加编号（Rxx / Mxx），并同步调研文参考文献列表。  
- 新增分子：先写入 `dual_target_cocrystal_catalog.csv`，再在本文件第三节补条目与链接。  
- 本文件为**人类可读总表**；CSV 为评估流水线输入源。
