# 双靶分子文献清单（JMC100 + 高分刊 / EJMC 补充）

> 更新：2026-07-23  
> **核心池：** 原 JMC 优先 100 篇（不变）  
> **本次补充：** JACS / Angew / Nat Commun / Cell Chem Biol / Chem Sci 等高分刊 **22** 篇 + EJMC 较新补充 **16** 篇  
> **合并机读表：** [`../data/literature/dual_target_papers_merged.csv`](../data/literature/dual_target_papers_merged.csv)（共 **138** 篇）  
> 原 100 篇细表仍见 [`DUAL_TARGET_PAPERS_JMC100.md`](DUAL_TARGET_PAPERS_JMC100.md)

---

## 0. 期刊现实：JACS / Nature 子刊会有双靶药化文吗？

**会有，但数量远少于 JMC，且体裁不同。**

| 期刊层级 | 双靶药化文常见吗？ | 典型形态 |
|----------|--------------------|----------|
| **J. Med. Chem. / Eur. J. Med. Chem.** | **主战场** | 设计–合成–两端–两端 SAR 双靶抑制剂 |
| **JACS / Angew** | **有，但少** | 概念创新强：双机制、结构导向多药理、少见靶组合；通常不是长 SAR 表 |
| **Nat. Commun.** | **有** | 结构生物学双通路抑制剂、AI 设计+湿实验验证、转化向双靶先导 |
| **Nat. Chem. Biol. / Cell Chem. Biol.** | **有（偏化学生物）** | 探针、共价双抑制、多药理学机制解析多于「药化系列优化」 |
| **Nature 正刊** | **极少** | 几乎不发传统双靶 SAR；偶见机制/靶点生物学，不宜当药化金标准池 |
| **Chem. Sci. / ACS Cent. Sci.** | **偶见** | 方法/化学生物/特殊骨架；需甄别是否真有合成+两端活性 |

**对本课题的用法：**  
- 抄「文献双靶金标准」分子 → **仍以 JMC/EJMC 为主**（结构+两端 IC50 最齐）。  
- JACS/Angew/Nat Commun → 作 **高概念案例 / 设计类型旁证 / AI+实验对标**，不要指望篇篇都有完整两端 SAR 表。

---

## 1. 合并库统计

| 来源分桶 | 篇数 |
|----------|------|
| jmc100_core（原清单） | 100 |
| high_impact（高分刊补充） | 22 |
| ejmc_supplement（EJMC 新补） | 16 |
| **合计** | **138** |

### 期刊分布（合并后）

| 期刊 | 篇数 |
|------|------|
| J Med Chem | 96 |
| Eur J Med Chem | 20 |
| Cell Chem Biol | 7 |
| Nat Commun | 5 |
| Chem Sci | 4 |
| JACS | 3 |
| Angew Chem Int Ed | 3 |

### 年份分布（合并后）

| 年份 | 篇数 |
|------|------|
| 2026 | 43 |
| 2025 | 29 |
| 2024 | 23 |
| 2023 | 8 |
| 2022 | 10 |
| 2021 | 14 |
| 2020 | 8 |
| 2019 | 3 |

---

## 2. 高分刊补充（22 篇，优先阅读）

机读：[`../data/literature/dual_target_papers_high_impact.csv`](../data/literature/dual_target_papers_high_impact.csv)

| # | 年 | 期刊 | 题目 | DOI链接 | 备注 |
|---|----|------|------|---------|------|
| 101 | 2025 | Chem Sci | Multi-target macrocycles: pyrogallol derivatives to control multiple pathological factors associated with Alzheimer's disease. | [10.1039/d4sc06417h](https://doi.org/10.1039/d4sc06417h) | 多靶大环/化学生物；Chem Sci：焦棓酚衍生物多病理因子 |
| 102 | 2025 | JACS | A Dual-Target and Dual-Mechanism Design Strategy by Combining Inhibition and Degradation Together. | [10.1021/jacs.4c11930](https://doi.org/10.1021/jacs.4c11930) | 降解剂/PROTAC/双机制；JACS：mTOR抑制+GSPT1降解 |
| 103 | 2025 | JACS | Dual Inhibitors of SARS-CoV-2 3CL Protease and Human Cathepsin L Containing Glutamine Isosteres Are Anti-CoV-2 Agents. | [10.1021/jacs.4c11620](https://doi.org/10.1021/jacs.4c11620) | 原创双靶小分子；JACS：3CLpro + cathepsin L |
| 104 | 2025 | Nat Commun | Discovery of multi-target anti-gout agents from Eurycoma longifolia Jack through phenotypic screening and structural optimization. | [10.1038/s41467-025-62645-6](https://doi.org/10.1038/s41467-025-62645-6) | 天然产物多靶；Nat Commun：东革阿里多靶抗痛风，需核小分子结构 |
| 105 | 2025 | Nat Commun | Targeting vascular adhesion protein-1 and myeloperoxidase with a dual inhibitor SNT-8370 in preclinical models of inflammatory disease. | [10.1038/s41467-025-58454-6](https://doi.org/10.1038/s41467-025-58454-6) | 原创双靶小分子；Nat Commun：VAP-1/MPO 双抑制 SNT-8370 |
| 106 | 2024 | Cell Chem Biol | Differential network analysis of ROS1 inhibitors reveals lorlatinib polypharmacology through co-targeting PYK2. | [10.1016/j.chembiol.2023.09.011](https://doi.org/10.1016/j.chembiol.2023.09.011) | 多药理学谱系；Cell Chem Biol：lorlatinib 多药理网络 |
| 107 | 2024 | Chem Sci | Structure-aware dual-target drug design through collaborative learning of pharmacophore combination and molecular simulation. | [10.1039/d4sc00094c](https://doi.org/10.1039/d4sc00094c) | 计算方法(少湿实验)；Chem Sci：结构感知双靶设计方法，作方法参考 |
| 108 | 2024 | Nat Commun | Automated design of multi-target ligands by generative deep learning. | [10.1038/s41467-024-52060-8](https://doi.org/10.1038/s41467-024-52060-8) | AI设计+合成验证；Nat Commun：CLM多靶配体，FXR/sEH等 |
| 109 | 2024 | Nat Commun | De novo generation of multi-target compounds using deep generative chemistry. | [10.1038/s41467-024-47120-y](https://doi.org/10.1038/s41467-024-47120-y) | AI设计+合成验证；Nat Commun：POLYGON，MEK1/mTOR等 |
| 110 | 2023 | Nat Commun | Structure-based discovery of dual pathway inhibitors for SARS-CoV-2 entry. | [10.1038/s41467-023-42527-5](https://doi.org/10.1038/s41467-023-42527-5) | 原创双靶小分子；Nat Commun：TMPRSS2+CTSL/CTSB |
| 111 | 2022 | Cell Chem Biol | The non-canonical target PARP16 contributes to polypharmacology of the PARP inhibitor talazoparib and its synergy with WEE1 inhibitors. | [10.1016/j.chembiol.2021.07.008](https://doi.org/10.1016/j.chembiol.2021.07.008) | 多药理学谱系；Cell Chem Biol：talazoparib 非经典靶 PARP16 |
| 112 | 2022 | Chem Sci | Repurposing of intestinal defensins as multi-target, dual-function amyloid inhibitors via cross-seeding. | [10.1039/d2sc01447e](https://doi.org/10.1039/d2sc01447e) | 肽/防御素多靶；Chem Sci：肠防御素多靶淀粉样+抗菌 |
| 113 | 2021 | Angew Chem Int Ed | Structure-Guided Design of G-Protein-Coupled Receptor Polypharmacology. | [10.1002/anie.202101478](https://doi.org/10.1002/anie.202101478) | 多药理学设计；Angew：GPCR多药理学结构导向 |
| 114 | 2021 | Angew Chem Int Ed | The Design of a GLP-1/PYY Dual Acting Agonist. | [10.1002/anie.202016464](https://doi.org/10.1002/anie.202016464) | 双靶激动/调节；Angew：GLP-1/PYY双激动 |
| 115 | 2021 | Cell Chem Biol | Evolution of kinase polypharmacology across HSP90 drug discovery. | [10.1016/j.chembiol.2021.05.004](https://doi.org/10.1016/j.chembiol.2021.05.004) | 多药理学谱系；Cell Chem Biol：HSP90研发中的激酶多药理 |
| 116 | 2021 | Cell Chem Biol | Exploiting polypharmacology to dissect host kinases and kinase inhibitors that modulate endothelial barrier integrity. | [10.1016/j.chembiol.2021.06.004](https://doi.org/10.1016/j.chembiol.2021.06.004) | 多药理学谱系；Cell Chem Biol：宿主激酶多药理 |
| 117 | 2021 | Chem Sci | Antimicrobial α-defensins as multi-target inhibitors against amyloid formation and microbial infection. | [10.1039/d1sc01133b](https://doi.org/10.1039/d1sc01133b) | 肽/防御素多靶；Chem Sci：α-防御素多靶 |
| 118 | 2020 | Angew Chem Int Ed | Potent Dual BET/HDAC Inhibitors for Efficient Treatment of Pancreatic Cancer. | [10.1002/anie.201915896](https://doi.org/10.1002/anie.201915896) | 原创双靶小分子；Angew：双BET/HDAC |
| 119 | 2020 | Cell Chem Biol | Discovery of Covalent MKK4/7 Dual Inhibitor. | [10.1016/j.chembiol.2020.08.014](https://doi.org/10.1016/j.chembiol.2020.08.014) | 原创双靶小分子；Cell Chem Biol：共价 MKK4/7 |
| 120 | 2019 | Cell Chem Biol | Divergent Polypharmacology-Driven Cellular Activity of Structurally Similar Multi-Kinase Inhibitors through Cumulative Effects on Individual Targets. | [10.1016/j.chembiol.2019.06.003](https://doi.org/10.1016/j.chembiol.2019.06.003) | 多药理学谱系；Cell Chem Biol：结构相近多激酶抑制谱差异 |
| 121 | 2019 | Cell Chem Biol | Multiomics Profiling Establishes the Polypharmacology of FDA-Approved CDK4/6 Inhibitors and the Potential for Differential Clinical Activity. | [10.1016/j.chembiol.2019.05.005](https://doi.org/10.1016/j.chembiol.2019.05.005) | 多药理学谱系；Cell Chem Biol：CDK4/6 抑制剂多组学多药理 |
| 122 | 2019 | JACS | Data-Driven Construction of Antitumor Agents with Controlled Polypharmacology. | [10.1021/jacs.9b08660](https://doi.org/10.1021/jacs.9b08660) | 多药理学设计；JACS：数据驱动多靶抗肿瘤分子 |

### 高分刊阅读优先级（合成向）

1. **明确双靶小分子 + 活性：** JACS 3CLpro/CTSL；Angew BET/HDAC；Nat Commun TMPRSS2/CTSL、VAP-1/MPO；Cell Chem Biol MKK4/7  
2. **AI + 已合成验证：** Nat Commun POLYGON / CLM（有湿实验，可作方法叙事对标，分子是否进金标准另核）  
3. **双机制/降解：** JACS mTOR 抑制 + GSPT1 降解 → **另表**，勿与普通双抑制混用  
4. **多药理学谱系文（Cell Chem Biol 多篇）：** 适合理解「单靶药其实多靶」，**不一定**是按双靶设计的 fused/linked 分子  

---

## 3. EJMC 较新补充（16 篇）

机读：[`../data/literature/dual_target_papers_ejmc_supplement.csv`](../data/literature/dual_target_papers_ejmc_supplement.csv)

| # | 年 | 期刊 | 题目 | DOI链接 | 备注 |
|---|----|------|------|---------|------|
| 123 | 2026 | Eur J Med Chem | An adamantylureido-benzylamide aniline as FLAP/sEH dual inhibitor: Rational design, in vitro and in vivo lipidomic profiling. | [10.1016/j.ejmech.2025.118338](https://doi.org/10.1016/j.ejmech.2025.118338) | 高分刊双靶小分子 |
| 124 | 2026 | Eur J Med Chem | Design and synthesis of c-Met/PARP dual-target inhibitors for the treatment of BRCA wild-type TNBC. | [10.1016/j.ejmech.2026.118722](https://doi.org/10.1016/j.ejmech.2026.118722) | 高分刊双靶小分子 |
| 125 | 2026 | Eur J Med Chem | Design, synthesis and bioactive evaluation of novel quinoline-linked sulfonamide-pyridine derivatives as PI3K/HDAC dual-target inhibitors. | [10.1016/j.ejmech.2025.118170](https://doi.org/10.1016/j.ejmech.2025.118170) | 高分刊双靶小分子 |
| 126 | 2026 | Eur J Med Chem | Design, synthesis and bioevaluation of novel combretastatin A-4 based derivatives as potent tubulin/HDAC6 dual-target inhibitors for cancer therapy. | [10.1016/j.ejmech.2025.118456](https://doi.org/10.1016/j.ejmech.2025.118456) | 高分刊双靶小分子 |
| 127 | 2026 | Eur J Med Chem | Design, synthesis and biological evaluation of novel guanidine-containing matrine derivatives as Topo I/II dual target inhibitors. | [10.1016/j.ejmech.2026.118884](https://doi.org/10.1016/j.ejmech.2026.118884) | 高分刊双靶小分子 |
| 128 | 2026 | Eur J Med Chem | Design, synthesis, and evaluation of dual-target inhibitors of acetylcholinesterase (AChE) and soluble epoxide hydrolase (sEH) for the treatment of Alzheimer's disease. | [10.1016/j.ejmech.2026.118844](https://doi.org/10.1016/j.ejmech.2026.118844) | 高分刊双靶小分子 |
| 129 | 2026 | Eur J Med Chem | Discovery of 3-indolealkylamines as novel dual-target σ(1)R/H(3)R ligands with potent analgesia. | [10.1016/j.ejmech.2026.118616](https://doi.org/10.1016/j.ejmech.2026.118616) | 双靶激动/调节为主 |
| 130 | 2026 | Eur J Med Chem | Discovery of novel bis-aryl urea-linked triazine derivatives as dual PI3K/mTOR inhibitors via scaffold hopping strategy and biological activity evaluations. | [10.1016/j.ejmech.2026.118856](https://doi.org/10.1016/j.ejmech.2026.118856) | 高分刊双靶小分子 |
| 131 | 2026 | Eur J Med Chem | Discovery of the first dual PD-L1/JAK inhibitor with enhanced in vivo antitumor immunity. | [10.1016/j.ejmech.2026.118605](https://doi.org/10.1016/j.ejmech.2026.118605) | 高分刊双靶小分子 |
| 132 | 2026 | Eur J Med Chem | Dual BRD4/AKT inhibition overcomes c-MYC-driven resistance in metastatic castration-resistant prostate cancer. | [10.1016/j.ejmech.2026.118824](https://doi.org/10.1016/j.ejmech.2026.118824) | 高分刊双靶小分子 |
| 133 | 2026 | Eur J Med Chem | Exploring dual inhibitors Carbonic Anhydrases and Phosphodiesterase 5 as potential agents for treatment Alzheimer's disease. | [10.1016/j.ejmech.2025.118404](https://doi.org/10.1016/j.ejmech.2025.118404) | 高分刊双靶小分子 |
| 134 | 2026 | Eur J Med Chem | Multi-target pyrazolopyrimidine-coumarin derivatives as potent CA IX/XII and tubulin polymerization inhibitors: Design, synthesis, and biological evaluation. | [10.1016/j.ejmech.2026.118789](https://doi.org/10.1016/j.ejmech.2026.118789) | 高分刊双靶小分子 |
| 135 | 2026 | Eur J Med Chem | Quinazoline-based dual-target inhibitors disrupt influenza virus RNP complex: Rational design, synthesis and mechanistic validation of potent anti-influenza agents. | [10.1016/j.ejmech.2025.118185](https://doi.org/10.1016/j.ejmech.2025.118185) | 高分刊双靶小分子 |
| 136 | 2026 | Eur J Med Chem | Rational design of dual PB2/JAK2 inhibitors achieving balanced antiviral and host-directed immunomodulatory effects. | [10.1016/j.ejmech.2026.118642](https://doi.org/10.1016/j.ejmech.2026.118642) | 高分刊双靶小分子 |
| 137 | 2026 | Eur J Med Chem | Synthesis and biological evaluation of phenanthridine derivatives as dual-target inhibitors of DNA topoisomerase IB (TOP1) and tyrosyl-DNA phosphodiesterase 1 (TDP1), and potential antitumor agents. | [10.1016/j.ejmech.2025.118541](https://doi.org/10.1016/j.ejmech.2025.118541) | 高分刊双靶小分子 |
| 138 | 2025 | Eur J Med Chem | Achieving dual-target fluorescent probes for tracing and inhibiting BRD4/PLK1 in tumor cells and tissues synchronously. | [10.1016/j.ejmech.2025.117886](https://doi.org/10.1016/j.ejmech.2025.117886) | 高分刊双靶小分子 |

---

## 4. 与原 100 篇如何一起用

1. 金标准抄录顺序建议：  
   - 原 JMC100 中 **B/C 类 + 有 IC50 数字**（见摘要分类）  
   - 本节 **高分刊「原创双靶小分子」**  
   - EJMC 补充中与冻结公开对相关者（如 PI3K/mTOR、AChE 等）  
2. 目标仍是 **30–80 个分子**；合并库是候选池扩容，不是要求全抄。  
3. PROTAC / 双机制 / 肽与防御素 / 纯方法文：分表或降优先级。

## 5. 本次高分刊检索可复现查询（PubMed）

```
("dual inhibitor"[Title] OR "dual-target"[Title] OR "dual-acting"[Title] OR polypharmacology[Title])
AND ("Journal of the American Chemical Society"[Journal]
  OR "Angewandte Chemie (International ed. in English)"[Journal]
  OR "Nature Communications"[Journal]
  OR "Nature Chemical Biology"[Journal]
  OR "Nature Chemistry"[Journal]
  OR "Cell Chemical Biology"[Journal]
  OR "Chemical Science"[Journal]
  OR "ACS Central Science"[Journal])
AND 2018:2026[DP]
```

另对 EJMC：`dual inhibitor/dual-target` + 2024:2026，去重后取较新且题录明确者。

检索日：2026-07-23。
