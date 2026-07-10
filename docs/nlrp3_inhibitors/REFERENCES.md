# NLRP3 / BAL 位点抑制剂发现——参考文献（含链接）

> 按主题分类；优先列出本课题直接引用的文献。  
> 最后更新：2026-07-10

---

## A. BAL 系列与 Glu-switch 结合位点（核心）

| # | 文献 | 链接 | 备注 |
|---|------|------|------|
| A1 | Hartman G et al. The discovery of novel and potent indazole NLRP3 inhibitors enabled by DNA-encoded library screening. *Bioorg Med Chem Lett* 2024;102:129675 | https://doi.org/10.1016/j.bmcl.2024.129675 | BAL-0028 发现；DEL 筛选；**事后 CADD** 预测结合位点 |
| A2 | Wilhelmsen K et al. Discovery of potent and selective inhibitors of human NLRP3 with a novel mechanism of action. *J Exp Med* 2025 | https://doi.org/10.1084/jem.20242403 | BAL-0028 机制；nanoDSF；与 MCC950 不同位点 |
| A3 | Torp J et al. Inhibition of NLRP3 by a CNS-penetrating indazole scaffold. bioRxiv 2025 | https://doi.org/10.1101/2025.07.01.662566 | **BAL-1516** 结构；PDB 9IHN/9Q8V；7PZC 为精修起点 |
| A4 | Wilhelmsen K et al. Discovery of a Potent and Selective Inhibitor of Human NLRP3 with a Novel Binding Modality. bioRxiv 2024 | https://doi.org/10.1101/2024.12.21.629867 | BAL-0028 机制预印本 |
| A5 | WO2025207644 — Indazole NLRP3 inhibitors (BAL 专利) | https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2025207644 | 专利数据来源之一 |

---

## B. NLRP3 结构生物学

| # | 文献 | 链接 | PDB |
|---|------|------|-----|
| B1 | Hochheiser IV et al. Structure of the NLRP3 decamer bound to CRID3. *Nature* 2022 | https://doi.org/10.1038/s41586-022-04467-w | [7PZC](https://www.rcsb.org/structure/7PZC) |
| B2 | Dekker C et al. Crystal Structure of NLRP3 NACHT with inhibitor. *J Mol Biol* 2021 | https://doi.org/10.1016/j.jmb.2021.167309 | [7ALV](https://www.rcsb.org/structure/7ALV) |
| B3 | Ohto U et al. Structural basis for oligomerization-mediated regulation of NLRP3. *PNAS* 2022 | https://doi.org/10.1073/pnas.2121353119 | [7VTP](https://www.rcsb.org/structure/7VTP), [7VTQ](https://www.rcsb.org/structure/7VTQ) |
| B4 | Murray JM, Johnson MC. CryoEM Structure of NLRP3 NACHT with G2394. *J Med Chem* 2022 | https://doi.org/10.1021/acs.jmedchem.2c01250 | [8ETR](https://www.rcsb.org/structure/8ETR) |
| B5 | BAL-1516 cryo-EM (Torp 2025) | https://doi.org/10.1101/2025.07.01.662566 | [9IHN](https://www.rcsb.org/structure/9IHN), [9Q8V](https://www.rcsb.org/structure/9Q8V) (**HPUB**) |

---

## C. 其他 NLRP3 结合位点（对照）

| # | 文献 | 链接 | 位点 |
|---|------|------|------|
| C1 | Coll RC et al. MCC950 directly targets NLRP3 ATP-hydrolysis motif. *Nat Chem Biol* 2019 | https://doi.org/10.1038/s41589-019-0277-7 | MCC950 口袋 |
| C2 | Jiang H et al. Identification of selective NLRP3 inhibitor CY-09. *J Exp Med* 2017 | https://doi.org/10.1084/jem.20171419 | Walker A |
| C3 | He H et al. Oridonin is a covalent NLRP3 inhibitor. *Nat Commun* 2018 | https://doi.org/10.1038/s41467-018-04947-6 | Cys279 |
| C4 | Huang Y et al. Tranilast directly targets NLRP3. *EMBO Mol Med* 2018 | https://doi.org/10.15252/emmm.201708689 | 寡聚化界面 |
| C5 | Hooftman A et al. Itaconate modifies NLRP3. *Cell Metab* 2020 | https://doi.org/10.1016/j.cmet.2020.07.016 | Cys548 (mouse) |
| C6 | Chen Y et al. RRx-001 covalent NLRP3 inhibitor. *Cell Mol Immunol* 2021 | https://doi.org/10.1038/s41423-021-00683-y | Cys409 |
| C7 | Stanton C et al. Covalent targeting NLRP3 inflammasome. *ACS Chem Biol* 2024 | https://doi.org/10.1021/acschembio.3c00330 | VLX1570 交联 |

---

## D. AI 共折叠 / 蛋白-配体预测基准

| # | 文献 | 链接 | 比较内容 |
|---|------|------|----------|
| D1 | Xu Y et al. Benchmarking all-atom biomolecular structure prediction with FoldBench. *Nat Commun* 2025 | https://doi.org/10.1038/s41467-025-67127-3 | **AF3 vs Boltz vs Chai vs Protenix**；558 蛋白-配体 |
| D2 | FoldBench 代码与数据 | https://github.com/BEAM-Labs/FoldBench | 基准复现 |
| D3 | Abramson J et al. Accurate structure prediction with AlphaFold 3. *Nature* 2024 | https://doi.org/10.1038/s41586-024-07487-w | AF3 原文；PoseBusters 对比 |
| D4 | Buttenschoen M et al. PoseBusters: AI docking methods fail physical validity. *Chem Sci* 2024 | https://doi.org/10.1039/D3SC04185A | DL 对接 vs Vina/GOLD |
| D5 | PoseBusters 工具文档 | https://posebusters.readthedocs.io/en/latest/ | PB-valid 检验 |
| D6 | SiteAF3: Accurate site-specific folding via conditional diffusion. 2025 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12595467/ | Pocket 约束提升对接 |
| D7 | Decoding the Allosteric Paradox (AI cofolding on allosteric sites). bioRxiv 2026 | https://doi.org/10.64898/2026.02.24.707829 | **变构位点预测短板** |
| D8 | Inductive Bio: Approaching AF3 docking accuracy (strong baseline) | https://www.inductive.bio/blog/strong-baseline-for-alphafold-3-docking | Vina+口袋 ≈ AF3 盲对接 |

---

## E. 无共晶时建模与对接（方法学先例）

| # | 文献 | 链接 | 相关性 |
|---|------|------|--------|
| E1 | Biscetti F et al. In Silico Insights towards NLRP3 Druggable Hot Spots. *Int J Mol Sci* 2019 | https://doi.org/10.3390/ijms20204974 | NLRC4 同源模建 + Glide 对接 MCC950 |
| E2 | Assessment of AI-Based Structure Prediction for NLRP3. *Molecules* 2022 | https://doi.org/10.3390/molecules27185797 | 7PZC 发布前 AF 模型 + MD 精修 |
| E3 | Ligand-guided homology modeling for H3 receptor. *PLoS ONE* 2019 | https://doi.org/10.1371/journal.pone.0218820 | 已知配体引导建模 + 虚拟筛选 |
| E4 | Ligand-directed modeling for GPCR virtual screening. *PLoS Comput Biol* 2017 | https://doi.org/10.1371/journal.pcbi.1005819 | LDM 方法 |
| E5 | Juliana C et al. Bay11-7082 direct inflammasome inhibitor. *J Biol Chem* 2010 | https://doi.org/10.1074/jbc.M109.082305 | NLRP3 ATP 酶抑制 |

---

## F. 活性机器学习与小数据筛选（方法参考）

| # | 文献 | 链接 | 备注 |
|---|------|------|------|
| F1 | 用户项目 XGBoost vs Chemprop（JNK1, n≈444） | 见 `JNK1_Selectivity_Project/` | **不可外推**至 NLRP3；小数据活性学习参考 |
| F2 | Cocco M et al. INF58 NLRP3 inhibitor. *J Med Chem* 2016 | https://doi.org/10.1021/acs.jmedchem.6b00452 | 共价抑制剂 SAR |
| F3 | Cocco M et al. INF39 NLRP3 inhibitor. *J Med Chem* 2017 | https://doi.org/10.1021/acs.jmedchem.7b00120 | 丙烯酸酯 SAR |

---

## G. 专利数据源

| 专利号 | 说明 | 原始数据 |
|--------|------|----------|
| WO2025207644 | BAL indazole 系列（主） | 用户 uploads CSV |
| WO2022204227 | +/++/+++ 活性分类 | 用户 uploads CSV |
| WO2024064655 | 部分活性标签 | 用户 uploads CSV |
| WO2023147468 | 需 N-氧化物 SMILES 修复 | 用户 uploads CSV |
| WO2026054623 | 9 条 nM IC50 | 用户 uploads CSV |

---

## H. 工具与数据库链接

| 资源 | 链接 |
|------|------|
| RCSB PDB | https://www.rcsb.org |
| AlphaFold Server | https://alphafoldserver.com |
| Boltz | https://github.com/jwohlwend/boltz |
| Chai-1 | https://github.com/chaidiscovery/chai-lab |
| GNINA | https://github.com/gnina/gnina |
| HADDOCK | https://wenmr.science.uu.nl/haddock2.4/ |
| Enamine REAL | https://enamine.net/compound-collections/real-compounds |
| SwissADME | http://www.swissadme.ch |
| PubChem | https://pubchem.ncbi.nlm.nih.gov |

---

## 推荐阅读顺序（本课题）

1. **A1 → A2 → A3**（BAL 发现与结构完整故事）
2. **B1**（7PZC 结构，建模起点）
3. **D1, D7**（AI 方法能力边界，尤其变构位点）
4. **D6**（如何加 pocket 约束）
5. **E1, E3**（无共晶建模先例）
6. 本仓库 **BAL_SITE_INHIBITOR_DISCOVERY_PROJECT.md**（实施方案）
