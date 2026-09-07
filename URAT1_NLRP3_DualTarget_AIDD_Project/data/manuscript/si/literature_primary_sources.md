# 原始文献核对（主张收缩）

核对日期：2026-09-07。只收缩已有主张，不新增靶点或暴露结论。12/21 名单不改。

## 1. PF-03882845 是 MR 拮抗剂（Meyers 等，2010）

- Meyers MJ 等，*J. Med. Chem.* 2010, 53, 5979–6002. DOI [10.1021/jm100505n](https://doi.org/10.1021/jm100505n)。
- 标题化合物即 `(3S,3aR)-2-(3-chloro-4-cyanophenyl)-3-cyclopentyl-3,3a,4,5-tetrahydro-2H-benzo[g]indazole-7-carboxylic acid`，文中编号 `3S,3aR-27d` / PF-3882845。与临床库 SMILES / CHEMBL1215331 一致。
- 文中：MR 结合 IC50 2.7 nM；PR 结合 IC50 310 nM（约 115 倍选择）；Dahl 盐敏感大鼠上降压、降尿白蛋白并有肾脏保护，因而进入临床。
- 该文**没有** URAT1、NLRP3 或痛风实验。

## 2. 肾损伤模型中的炎症基因表达（Orena 等，2013）

- Orena S 等，*Front. Pharmacol.* 2013, 4:115. PMC [3796291](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3796291/)。
- 单肾切除 SD 大鼠 + 高盐 + 醛固酮泵。PF-03882845 抑制肾脏 *Il-6*、*Tgf-β1*、*Icam-1*、osteopontin 等基因表达，并降低 UACR。
- 这支持**上游炎症通路混杂**（醛固酮–MR），不能当作 NACHT 结合或 IL-1β 下降即直接 NLRP3 抑制。

## 3. 已有 URAT1/NLRP3 双节点化学型，堵住“首次”主张

**HNW005.** Sun 等，*Eur. J. Med. Chem.* 2025. DOI [10.1016/j.ejmech.2025.117644](https://doi.org/10.1016/j.ejmech.2025.117644)，PMID 40286449。标题与摘要将其定义为 NLRP3 inflammasome 与 URAT1 双抑制剂，用于痛风性关节炎。厂商页引用该文给出 NLRP3 激活 IC50 1.7 µM、KD 204.6 nM、URAT1 转运 IC50 6.4 µM。本次**未取到论文 SI 结构式**；Tc 表里的 HNW005 SMILES 仍标为厂商来源，不得写成 SI。

**Compound 32.** Zhang 等，*Nat. Commun.* 2025, 16:7430. DOI [10.1038/s41467-025-62645-6](https://doi.org/10.1038/s41467-025-62645-6)。正文：HEK293-URAT1 对 14C-尿酸 IC50 3.81 ± 0.99 µM；NLRP3 SPR KD 27.8 µM；OAT4 IC50 6.08 ± 2.58 µM；另报 GLUT9 抑制与 BMDM IL-1β。SI（MOESM1）化合物 **(32)** 为 `1-((4-bromonaphthalen-1-yl)methyl)-N-(thiophen-2-ylsulfonyl)-1H-indole-2-carboxamide`，C24H17BrN2O3S2。由此重建的 SMILES 与原先厂商串一致，**以后引用 SI IUPAC 名，不引厂商页**。

这两篇已经存在实验双节点/多靶抗痛风分子。本工作的贡献是临床库 + 几何门框架，不是“首次双靶策略”。

## 明确删除或禁止的句子

- 不得由糖尿病肾病适应症或肾脏保护动物模型推断 URAT1 **顶膜游离暴露**。
- 不得把 PubMed 对 `PF-03882845 AND (URAT1 OR NLRP3 OR gout)` 零命中写成**已确认新靶点**。零命中只说明公开检索未把该药连到这些词，不是专利或实验阴性证明。
- 不得写 dual inhibitor / 首次发现。

临床终止：NCT01488877 等因策略/入组停止，不进入效力或靶点论证。
