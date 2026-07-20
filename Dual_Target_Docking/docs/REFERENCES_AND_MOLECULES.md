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

### 1.5 评价基准与共识打分

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R12 | Li Y et al. **CASF** 系列（如 CASF-2016）：scoring / ranking / docking / screening power. | [CASF-2016 doi:10.1021/acs.jcim.7b00650](https://doi.org/10.1021/acs.jcim.7b00650) · [CASF 主页](http://www.pdbbind.org.cn/casf.php) | 单靶打分函数标准评测协议。双靶研究可迁移其组件指标（姿态成功率、排序、富集），但**尚无**官方“双靶 CASF”。 |
| R13 | Ericksen SS et al. *Machine Learning Consensus Scoring Improves Performance Across Targets.* **JCIM** 2017. | [doi:10.1021/acs.jcim.7b00153](https://doi.org/10.1021/acs.jcim.7b00153) | 机器学习共识打分提升结构基 VS 跨靶稳健性；可嵌入双靶流水线的单靶侧打分层。 |
| R14 | Consensus docking 综述 / 调查（多程序降低靶点依赖方差）. | [PMC9821981](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9821981/) | 共识对接降低单程序偏差；双靶场景可对每个靶做多引擎共识后再跨靶融合。 |

### 1.6 生成式双靶（前沿算法）

| # | 文献 | 链接 | 内容简介 |
|---|------|------|----------|
| R15 | DualDiff / CompDiff. *Reprogramming Pretrained Target-Specific Diffusion Models for Dual-Target Drug Design.* **NeurIPS** 2024. | [arXiv:2410.20688](https://arxiv.org/abs/2410.20688) | 将单靶扩散模型重编程到双靶：口袋对齐策略（center / RMSD-anchor / score-anchor）；评估用 P1/P2 Vina、**Max Vina**、Dual High Affinity、双姿态 RMSD 等。偏生成，指标可借鉴。 |
| R16 | FuseDiff. *Symmetry-Preserving Joint Diffusion for Dual-Target Structure-Based Drug Design.* | [arXiv:2603.05567](https://arxiv.org/abs/2603.05567) | 对称保持的联合扩散双靶 SBDD；与 DualDiff 同属生成式前沿，非本课题首选落地路径。 |

### 1.7 讨论中提及的单靶深度打分（对比用，非双靶专用）

| # | 文献 / 工具 | 链接 | 内容简介 |
|---|-------------|------|----------|
| R17 | OnionNet / OnionNet-2（蛋白–配体亲和力 CNN） | 检索 [OnionNet protein-ligand](https://pubmed.ncbi.nlm.nih.gov/?term=OnionNet+protein-ligand) | 基于旋转对称壳层特征的亲和力预测；可作为单靶侧 scorer 候选，但**不是**双靶融合方法本身。 |
| R18 | DeepRMSD / SFCT 等姿态质量相关方法 | 检索对应关键词 | 偏姿态质量/打分校正；本课题建议作门控特征，而非主创新。 |
| R19 | GNINA（CNN 辅助对接） | [github.com/gnina/gnina](https://github.com/gnina/gnina) | Vina 系采样 + CNN 亲和/姿态分；推荐作为可插拔引擎。 |
| R20 | RTMScore 等现代结构基打分 | 检索 `RTMScore docking` | 可用于单靶侧替换/对照；创新仍在跨靶融合。 |

### 1.8 调研流程工具（非学术论文）

| # | 资源 | 链接 | 说明 |
|---|------|------|------|
| T01 | PaperSpine | [github.com/WUBING2023/PaperSpine](https://github.com/WUBING2023/PaperSpine) | 研究问题结构化、SOTA gap、证据分层 |
| T02 | ACS Writer v2 | [github.com/Caosmart1979/acs-writer-v2](https://github.com/Caosmart1979/acs-writer-v2) | 文献探索与方法学评价工作流 |

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
- https://doi.org/10.1021/acs.jcim.7b00650  
- https://doi.org/10.1021/acs.jcim.7b00153  
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9821981/  
- https://arxiv.org/abs/2410.20688  
- https://arxiv.org/abs/2603.05567  

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
