# 课题确立分析报告：双靶筛选中的「假双靶」问题与 Dual-VSDS

> **文档性质：** 说明「双靶存在什么问题」这一判断是怎么来的、为何据此立项、与最初「按类型把对接做准」有何取舍。  
> **日期：** 2026-07-23（本版补全可核对文献）  
> **状态：** 问题诊断与立项依据已齐；**主结论（朴素拼分有害）仍属待证假说**，须由本项目打假表证实或证伪。  
> **引用原则：** 下文凡写「文献表明 / 已报道」者，均给出可点击 DOI 或 PMC；**实验室内部动机**单独标注，不伪装成已发表结论。未找到可核对出处的论断一律降级为「本课题假说」或删除。  
> **关联：** [`PROJECT_MASTER_PLAN.md`](PROJECT_MASTER_PLAN.md) · [`REFERENCES_AND_MOLECULES.md`](REFERENCES_AND_MOLECULES.md) · [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)

---

## 1. 报告要回答的三个问题

1. 「双靶计算筛选有问题」这一说法，**依据什么推出**，哪些已坐实、哪些还只是假说？  
2. **为什么要开展本课题**（而不是继续调对接、或只做实验室个案）？  
3. 相对「让不同类型双靶分子对接更准」的最初想法，**为何主线改到现在这样**？

---

## 2. 「双靶存在的问题」——结论是如何得到的

### 2.1 先分清：问题有两层，不能混谈

| 层级 | 问题表述 | 主要伤害 | 证据成熟度（立项时） | 主要文献支撑 |
|------|----------|----------|----------------------|--------------|
| **层甲：对接与类型** | 融合 / 连接 / 拼合等设计类型不同，对接难度与成药性约束不同；连接型构象自由度更高 | 姿态难、类型难分治 | **部分已知** | Morphy & Rankovic 设计分类与理化挑战[^1]；TwistDock 等说明「双端+连接」需专用策略且场景有限[^2]；本仓库共晶编目显示双端同配体共晶极少（见 §2.3、[^3]） |
| **层乙：分数怎么合** | 两端各对接一次后，用排名/分数聚合挑「双靶」候选，可能混入假双靶 | 名单污染 | **文献有近亲与警示；本项目打假表尚未产出** | 双靶对接 VS 假阳性高[^4]；两端对接再取交集/融合为实务路径[^5][^6]；单靶 VS 中「负例设定」可翻转结论[^7] |

本课题**主攻层乙**；层甲降为子问题与加分项。

### 2.2 推导链（观察 → 机制 → 可证伪命题）

#### 步骤一：实务默认做法（有文献记录的操作）

已发表工作表明，结构基双靶虚拟筛选常采用：

1. **同一分子对两个靶分别对接**（异构靶对上的明确案例：A2A / MAO-B）[^5]；  
2. 再要求分子在**两端排序都靠前**，或对两端（或多套）打分做**融合 / 共识**后再排序[^5][^6]；  
3. 取靠前分子做实验验证。

Pérez-Castillo 等明确写出：要提高双靶配体富集，**必须组合两端对接信息**，且往往还要融合多套打分函数的排名；单打分不够用[^6]。  
Jaiteh 等则从 540 万分子库两端对接后，取「两端都高排」的 24 个分子做实验，得到 4 个双活分子——说明「两端对接 → 取交集式高排」是可发表的标准路径，但**命中率仍有限**（24 测 4 双活）[^5]。

> **本课题据此将「两端独立对接 + 朴素聚合（均值 / 求和 / 排名平均等）」列为必须对比的基线。**  
> 这不是猜测实验室习惯，而是对上述已发表流程的抽象。

#### 步骤二：为何可能系统性出错（机制 + 文献警示）

三条机制中，**①③有直接文献；②为药化常识并与设计综述一致，但「平均分抬高假双靶」本身仍待本项目检验。**

1. **两边分数/排名尺度未必可直接算术合并**  
   共识/融合打分文献强调：单函数不稳定，需融合多排名才改善富集——侧面说明「随便用一个原始分相加」并不可靠[^6][^8]。跨靶时尺度问题更严重（本课题假说的核心之一）。

2. **双靶分子两端活性不必对称；设计还受分子量等约束**  
   Morphy & Rankovic 系统讨论多靶配体的理化与设计约束（linked / fused / merged）[^1]；Proschak 等近年梳理双靶导向配体的靶组合与化学策略[^9]。这些支持「双靶 ≠ 两端分数一样高」，但**不直接证明**平均分会抬高假双靶。

3. **评测设定会改变「方法看起来好不好」**  
   Zhou 等在多组激酶对上评估对接 VS 搜寻双靶抑制剂，结论包括：对接可识别单靶抑制剂，但对双靶抑制剂**假阳性高、富集有限**，需联用其他方法[^4]。  
   Gu 等（VSDS-VD）在单靶虚拟筛选视角证明：**TrueDecoy 与 RandomDecoy 上的富集结论可以相反**；只报一套负例会误导[^7]。  
   → 双靶任务若只用「随机诱饵」、不设「只强一端」硬负例，同样可能虚高——这是立项采用四类标签与双负例协议的直接文献依据[^7]。

由此得到立项用的**核心假说（本项目数据尚未证实）**：

> 在「两端都有活性标签」的分子集合上，朴素拼分会相对抬高「只强甲 / 只强乙」分子的名次；先按靶校准再强调短板的合并，应能降低这种污染。

#### 步骤三：为何「值得立项」而不是空想（证据表）

| 证据类型 | 可核对内容 | 支撑什么 | 不支撑什么 |
|----------|------------|----------|------------|
| **已发表双靶对接流程** | Jaiteh 2018[^5]；Pérez-Castillo 2017[^6] | 「两端对接 + 聚合/取高排」是真实方法路线 | 不证明你们实验室名单已被污染 |
| **双靶对接 VS 警示** | Zhou 2013[^4] | 双靶对接 VS **假阳性高、富集有限** | 不自动等于「平均分有害」的定量证明 |
| **评测设定可翻转结论** | VSDS-VD 2025[^7]；CleanSplit / 泄漏讨论 2025[^10]；PoseBusters 2024[^11] | 必须严肃设计负例、泄漏与物理合理性；不能只看 RMSD | 这些是单靶/打分文，不是双靶打假表 |
| **生成多、判别尺子相对弱** | POLYGON 等可生成并合成验证多靶分子[^12]；本课题定位是排序/评测 | 「判别尺子」有独立价值 | 不表示生成工作无意义 |
| **设计类型与共晶稀缺** | Morphy[^1]；Tanaka Mcl-1/Bcl-xL 双端共晶个案[^3]；本仓库 catalog | 类型问题真实，但双端姿态金标准极少 → 不宜单独立主 | 不否定类型化对接作为子问题 |
| **数据可测性（本仓库）** | [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)（ChEMBL 四类统计，2026-07-23） | 假说有可操作靶点对 | ChEMBL 配对 ≠ 与所选晶体条件完全一致 |
| **实验室语境** | **本课题组内部需求**（NLRP3/JNK1 等双靶分子与细胞读数；资源有限） | 说明应用出口 | **非公开文献，不得写成领域共识** |

#### 步骤四：现在仍不能写进「结果」的句子

- 「我们已经证明简单平均有害」——**否**（待打假表）；  
- 「校准短板一定优于朴素融合」——**否**；  
- 「对接排序等于真实双结合强弱」——**否**（公开活性与晶体条件错配是已知局限，见 §5.3）。

立项阶段允许的结论：

- **问题结构有文献锚点、可检验、有可测数据与应用出口，因此值得开展。**

### 2.3 层甲（类型与对接准度）的文献位置

- **设计分类：** linked / fused / merged 及成药性挑战，见 Morphy & Rankovic[^1]。  
- **连接/双价特殊对接：** TwistDock 针对同蛋白双结构域 bivalent，先固定两端再扭 linker，**不是**普通异构双靶通解[^2]。  
- **双端同配体共晶金标准稀缺：** 典型正面例子为 Tanaka 等 Mcl-1/Bcl-xL 杂交抑制剂 compound 10（PDB 3WIY / 3WIZ）[^3]；本仓库编目仅数个 Tier A 体系（见 [`REFERENCES_AND_MOLECULES.md`](REFERENCES_AND_MOLECULES.md)）。  
→ 「按类型全面把对接做准并大规模证明」缺姿态标签，**不宜作为主命题**；宜作分层与质检。

---

## 3. 为什么要开展这个课题

### 3.1 科学理由（对应文献）

1. **默认双靶对接 VS 已有「假阳性高」的公开警示**[^4]，但缺少以「只强一端」为硬负例、跨多靶点对、可复现的**配对评测基准**——相对 Pérez-Castillo 的单靶点对融合案例[^6]与 VSDS-VD 的单靶 VS 审计[^7]，仍有任务级缺口。  
2. **评测纪律已被证明能改变方法排名**（负例类型[^7]、训练泄漏[^10]、物理有效性[^11]）→ 双靶拼分若无纪律，结论不可信。  
3. **生成类工作增多**[^12]，更需要可公开复现的**判别/排序尺子**，否则无法比较「谁真的更双靶」。

### 3.2 应用理由

- **本课题组内部：** 需要计算侧门槛，减少把单靶偏倚分子推入合成；细胞数据宜作盲测排序参考，不宜在缺亲和力时宣称双结合（内部约束，非文献）。  
- **与资源匹配：** 方法 + 开放基准不依赖立刻补齐动物药效。

### 3.3 可行性理由（已完成、可核对）

| 动作 | 依据 |
|------|------|
| 冻结三组公开靶点对 | [`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)；ChEMBL / RCSB 审计 |
| 文献候选池 | [`DUAL_TARGET_PAPERS_CATALOG.md`](DUAL_TARGET_PAPERS_CATALOG.md) |
| 引擎可复现原则 | GNINA 为已发表开源对接程序（Vina/smina 采样 + CNN 打分）[^13]；PoseBusters 作物理门控[^11]；不自研采样器 |

### 3.4 发表与风险（诚实边界）

- **叙事可对标：** 诊断型评测 + 开放基准（同构于 VSDS-VD 的「设定可翻转」贡献[^7]，但任务改为配对双靶）。  
- **主要风险：** 配对标签稀缺；若所有融合头一起失败 → 按方案视为可证伪结果（对接信息不足），不是隐瞒。  
- **禁止过声称：** 不得把细胞活性写成双结合金标准；不得无检索限定写「全球首个」。

---

## 4. 与最初想法的对比（文献如何支持取舍）

| 维度 | 最初：类型化对接更准 | 现在：防假双靶拼分诊断 |
|------|----------------------|-------------------------|
| 文献锚点 | Morphy 分类[^1]；TwistDock 特例[^2]；双端共晶个案[^3] | Zhou 假阳性警示[^4]；Jaiteh/Pérez-Castillo 流程[^5][^6]；VSDS-VD 评测纪律[^7] |
| 金标准密度 | 双端共晶极少[^3] | 公开库配对活性可规模化（有错配局限） |
| 主线选择 | 证据形态弱，易成调参工程 | **更值得作主线**；类型问题作子问题 |

---

## 5. 课题边界

### 5.1 承诺回答

- 在冻结靶点对与统一对接协议下，朴素拼分是否更易抬高「只强一端」？  
- 校准短板合并能否改善双强相对只强一端的排序？  
- 换靶点对、换 TrueDecoy / RandomDecoy 类设定后结论是否大体同向？（设定敏感性的文献先例见[^7]）

### 5.2 不承诺回答

- 生理条件下一定双结合；对接分 = 亲和力绝对值；公开活性一定对应所选晶体；已解决全部类型姿态精度；动物药效预测。

### 5.3 已知局限（文献与方法学共识层面）

- 对接 VS 本身受打分与设定影响大[^4][^7][^11]；  
- 公开活性—结构条件错配是结构基筛选的普遍风险（故主张必须写窄）；  
- 公开配对集多数无设计类型标签（类型靠文献 curated 补）。

---

## 6. 立项判决

| 项 | 内容 |
|----|------|
| **问题** | 双靶筛选常用「两端对接 + 分数/排名聚合」[^5][^6]；在双靶任务上已有假阳性高的报道[^4]；若负例与评测不严，方法优劣会被误判[^7]。 |
| **本课题假说** | 朴素聚合会相对抬高「只强一端」；校准短板可纠偏——**待证**。 |
| **为何开展** | 有文献锚点的流程与警示 + 可测的公开配对数据 + 实验室应用出口；比「全面类型化对接精度」更具备可完成证据形态。 |
| **成功标准** | 打假表成立；或诚实报告融合头集体失败。 |

**最终一句话：**  
立项依据来自**已发表的双靶对接流程与假阳性/评测设定警示**，而不是来自「本项目已经证明平均分有害」；开展本课题是为了把「假双靶」风险变成可证伪实验，而不是把假说提前写成结论。

---

## 7. 论断—文献对照总表（防胡编）

| 编号 | 报告中的论断 | 文献 / 来源 | 证据强度 |
|------|--------------|-------------|----------|
| A1 | 双靶 VS 常两端分别对接再取高排/融合 | Jaiteh et al., *J. Med. Chem.* 2018[^5]；Pérez-Castillo et al., *Curr. Neuropharmacol.* 2017[^6] | 强（方法描述） |
| A2 | 双靶对接 VS 假阳性高、富集有限 | Zhou et al., *JCIM* 2013[^4] | 强（直接结论） |
| A3 | 负例/评测设定可翻转方法排名 | Gu et al., *Nat. Mach. Intell.* 2025（VSDS-VD）[^7] | 强（单靶 VS；迁移为协议依据） |
| A4 | 只看 RMSD 不够，需物理合理性 | Buttenschoen et al., *Chem. Sci.* 2024（PoseBusters）[^11] | 强 |
| A5 | 训练—测试泄漏可虚高打分表现 | Graber et al., *Nat. Mach. Intell.* 2025（CleanSplit）[^10] | 强（亲和预测；作纪律类比） |
| A6 | 多靶配体有 linked/fused/merged 等类型与理化挑战 | Morphy & Rankovic, *J. Med. Chem.* 2006[^1]；Proschak et al., *J. Med. Chem.* 2024[^9] | 强 |
| A7 | 双端同配体共晶存在但稀少（有明确个案） | Tanaka et al., *J. Med. Chem.* 2013[^3]；本仓库 catalog | 中—强 |
| A8 | 连接/双价可有专用对接策略，非通解 | Bai et al., TwistDock, 2019[^2] | 中 |
| A9 | 多靶分子生成+实验验证已有高分刊工作 | Munson et al., *Nat. Commun.* 2024（POLYGON）[^12] | 强（生成侧） |
| A10 | GNINA 为可引用的开源对接实现 | McNutt et al., *J. Cheminform.* 2021[^13] | 强（工具） |
| A11 | 「朴素平均会系统性抬高只强一端」 | **本课题假说** | 待证（**无**现成定量论文可直接等同） |
| A12 | 本实验室名单已被简单拼分坑过 | **内部经验** | 非公开文献 |

---

## 8. 参考文献（仅列正文引用者）

[^1]: Morphy R, Rankovic Z. The physicochemical challenges of designing multiple ligands. *J. Med. Chem.* 2006. [doi:10.1021/jm0603015](https://doi.org/10.1021/jm0603015)

[^2]: Bai L, et al. TwistDock: Twist-and-Dock for Bivalent Ligand Binding. *Drug Des. Devel. Ther.* 2019. [doi:10.2147/DDDT.S194276](https://doi.org/10.2147/DDDT.S194276)

[^3]: Tanaka Y, et al. Discovery of Potent Mcl-1/Bcl-xL Dual Inhibitors… *J. Med. Chem.* 2013. [doi:10.1021/jm401170c](https://doi.org/10.1021/jm401170c)（PDB：3WIY / 3WIZ）

[^4]: Zhou S, et al. Feasibility of Using Molecular Docking-Based Virtual Screening for Searching Dual Target Kinase Inhibitors. *J. Chem. Inf. Model.* 2013. [doi:10.1021/ci400065e](https://doi.org/10.1021/ci400065e)

[^5]: Jaiteh M, et al. Docking Screens for Dual Inhibitors of Disparate Drug Targets for Parkinson’s Disease. *J. Med. Chem.* 2018, 61, 5269–5278. [doi:10.1021/acs.jmedchem.8b00204](https://doi.org/10.1021/acs.jmedchem.8b00204) · [PMC6716773](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6716773/)

[^6]: Pérez-Castillo Y, et al. Fusing Docking Scoring Functions Improves the Virtual Screening Performance for Discovering Parkinson’s Disease Dual Target Ligands. *Curr. Neuropharmacol.* 2017, 15, 1107–1116. [doi:10.2174/1570159X15666170109143757](https://doi.org/10.2174/1570159X15666170109143757) · [PMC5725543](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5725543/)

[^7]: Gu S, Shen C, Zhang X, et al. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat. Mach. Intell.* 2025, 7, 509–520. [doi:10.1038/s42256-025-00993-0](https://doi.org/10.1038/s42256-025-00993-0)

[^8]: Ericksen SS, et al. Machine Learning Consensus Scoring Improves Performance Across Targets. *J. Chem. Inf. Model.* 2017. [doi:10.1021/acs.jcim.7b00153](https://doi.org/10.1021/acs.jcim.7b00153)

[^9]: Proschak E, et al. Polypharmacology: A Systematic Investigation of Dual-Target-Directed Ligands. *J. Med. Chem.* 2024. [doi:10.1021/acs.jmedchem.4c00838](https://doi.org/10.1021/acs.jmedchem.4c00838)

[^10]: Graber D, et al. Resolving data bias improves generalization in binding affinity prediction（PDBbind CleanSplit）. *Nat. Mach. Intell.* 2025. [doi:10.1038/s42256-025-01124-5](https://doi.org/10.1038/s42256-025-01124-5)

[^11]: Buttenschoen M, et al. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chem. Sci.* 2024. [doi:10.1039/D3SC04185A](https://doi.org/10.1039/D3SC04185A)

[^12]: Munson BP, et al. De novo generation of multi-target compounds using deep generative chemistry（POLYGON）. *Nat. Commun.* 2024, 15, 3636. [doi:10.1038/s41467-024-47120-y](https://doi.org/10.1038/s41467-024-47120-y)

[^13]: McNutt AT, et al. GNINA 1.0: molecular docking with deep learning. *J. Cheminform.* 2021. [doi:10.1186/s13321-021-00522-2](https://doi.org/10.1186/s13321-021-00522-2)

完整文献库与共晶编目见 [`REFERENCES_AND_MOLECULES.md`](REFERENCES_AND_MOLECULES.md)。

---

## 9. 附件索引

- 靶点对硬门槛审计：[`PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)  
- 执行总览：[`PROJECT_MASTER_PLAN.md`](PROJECT_MASTER_PLAN.md)  
- Idea / Scoop 审计：[`researchstudio_audit/RESEARCHSTUDIO_AUDIT.md`](researchstudio_audit/RESEARCHSTUDIO_AUDIT.md)  
- RQ / 审稿攻击面：[`ars_audit/ARS_TOPIC_ANALYSIS.md`](ars_audit/ARS_TOPIC_ANALYSIS.md)
