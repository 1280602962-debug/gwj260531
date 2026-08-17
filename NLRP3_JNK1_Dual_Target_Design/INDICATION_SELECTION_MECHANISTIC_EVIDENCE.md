# NLRP3 / JNK1 双靶点：适应症选择（机制文献核对版）

本文只讨论**两个靶点均已有机制实验**（遗传敲除/敲入、磷酸化位点突变、或选择性药理抑制）的适应症。  
引用条目以文末「本文参考文献」及 `REFERENCES.md` 交叉编号为准；PMID 均经 PubMed/DOI 核对。

> **使用说明**：本文替代主文档第 11 节与文献综述第 7.2 节中较简略、未区分证据等级的适应症排序。那两处已改为指向本文。

---

## 1. 入选标准

一个适应症只有在**同时满足**下列两条时，才能写成「JNK1 + NLRP3 双靶适应症」：

| 条件 | 必须具备的证据 |
|------|----------------|
| **NLRP3 在该病中必要** | 该疾病模型中：`Nlrp3` 敲除、炎症小体组分缺失、或 NLRP3 选择性抑制剂改变表型 |
| **JNK1 在该病中必要** | 该疾病模型中：`Jnk1`（*Mapk8*）敲除、细胞类型特异性缺失、骨髓嵌合体、或可归因于 JNK1 的药理抑制改变表型 |

更强一档（**偶联已在该模型直接验证**）：同一模型中证明 **JNK1 → NLRP3-Ser194（小鼠）/ Ser198（人）磷酸化** 决定炎症表型。公开文献中，这一档**仅出现在少数模型**（见第 3 节）。

**不作为入选依据的材料**：综述性「可能相关」、仅表达量升高、仅 pan-JNK 抑制剂（如 SP600125，选择性不足）、无该疾病模型的遗传数据。

---

## 2. 核心偶联机制（所有适应症的共同上游）

**Song 等，*Mol Cell* 2017** [8] 是目前把两个靶点连在同一分子事件上的关键论文。经核对，第一作者为 **Nan Song**（不是 Song H）。实验事实包括：

1. **JNK1（不是 JNK2）** 在 priming 阶段直接磷酸化小鼠 NLRP3 **Ser194**（人 **Ser198**）；
2. 该磷酸化促进 NLRP3 去泛素化与自寡聚，是炎症小体活化的必要 priming 事件；
3. `Jnk1` 缺失巨噬细胞中 NLRP3 炎症小体活化被抑制；JNK 抑制剂抑制 NLRP3，**不抑制 AIM2 / NLRC4**；
4. **NLRP3-S194A 敲入小鼠**在体内被保护于：
   - **MSU 诱导的腹膜炎**；
   - **LPS 诱导的内毒素血症**。

因此：凡主张「双靶优于单靶」，生物学上最硬的表述是——**打断 JNK1 对 NLRP3-S194 的 priming 许可，并直接抑制 NLRP3 组装**。能直接沿用这句话的疾病，只有 Song 做过的两类模型；其余疾病只能写「两靶各自被证明参与该病；偶联机制来自巨噬细胞 priming，尚未在该病中用 S194A 复现」。

---

## 3. 第一档：偶联已在该模型被直接验证（优先）

### 3.1 痛风 / MSU 晶体性炎症

这是目前唯一同时满足：**NLRP3 致病、JNK1 致病、JNK1–NLRP3-S194 偶联均已在同一类模型中验证**。

| 靶点 / 环节 | 机制实验 | 文献 |
|-------------|----------|------|
| **NLRP3** | `Nlrp3` / `Asc` / `Casp1` 缺失巨噬细胞对 MSU 无 IL-1β；缺失小鼠 MSU 腹膜炎中性粒浸润下降 | Martinon 等，*Nature* 2006 [54] |
| **NLRP3 药理** | MCC950 抑制 NLRP3 依赖的晶体炎症 | Coll 等，*Nat Med* 2015 [1] |
| **JNK1 → NLRP3** | S194A 敲入小鼠抵抗 **MSU 腹膜炎**；磷酸化由 JNK1 介导 | Song 等，*Mol Cell* 2017 [8] |
| **JNK（溶酶体破裂 / MSU）** | 溶酶体破裂经 Ca²⁺–CaMKII–TAK1–JNK 调控 ASC 寡聚与 NLRP3 活化；JNK 抑制降低 LLME/MSU 诱导的炎症小体活化 | Okada 等，*J Biol Chem* 2014 [55] |

**可写入开题的句子：**  
MSU 通过 NLRP3 炎症小体驱动痛风样炎症 [54]；JNK1 对 NLRP3-S194 的磷酸化是该活化的 priming 许可，S194A 小鼠在 MSU 腹膜炎中被保护 [8]。因此痛风 / 晶体性腹膜炎是**机制上最闭合**的双靶适应症。

**必须写明的限制：**

- Song 用的是 **MSU 腹腔注射**，不是人关节腔痛风石模型。
- Okada 2014 证明的是 **TAK1–JNK → ASC 寡聚**，不是 S194 位点本身；与 Song 2017 互补，但不是同一实验。
- JNK1 在「临床痛风关节」中的必要性，尚未达到 NLRP3 那样的遗传证据强度。

体内 POC 应优先复现 **MSU 腹膜炎 / 气囊 / 关节腔注射**，并设 NLRP3 抑制剂与 JNK1 偏向抑制剂对照。

### 3.2 LPS 诱导的全身炎症 / 内毒素血症

| 靶点 | 机制实验 | 文献 |
|------|----------|------|
| **JNK1 → NLRP3** | S194A 敲入小鼠抵抗 **LPS 内毒素血症** | Song 等，*Mol Cell* 2017 [8] |
| **NLRP3 药理** | MCC950 改善 NLRP3 依赖的炎症模型 | Coll 等，*Nat Med* 2015 [1] |

**定位：** 这是验证「双靶分子是否打到 priming 许可位点」的标准体内读出，适合**机制 POC**。不宜作为主要临床适应症（脓毒症病因混杂非经典 caspase-11 通路，转化风险高）。

---

## 4. 第二档：两靶在同一疾病中各自有遗传 / 药理证据，但 S194 偶联尚未在该病复现

开题必须写成：

> 两靶分别被证明驱动该病；JNK1–NLRP3-S194 偶联来自巨噬细胞 priming 的一般机制 [8]，**尚未**在该疾病模型中用 S194A 敲入复现。

### 4.1 饮食诱导的肥胖–胰岛素抵抗 / 代谢性炎症

**NLRP3 端证据扎实。JNK1 端「全身敲除」扎实，但「髓系 / 巨噬细胞 JNK1」存在文献分歧，不得写成定论。**

| 靶点 | 机制实验 | 文献 | 证据等级 |
|------|----------|------|----------|
| **JNK1（全身）** | `Jnk1` 缺失减轻肥胖与胰岛素抵抗 | Hirosumi 等，*Nature* 2002 [56] | 高 |
| **NLRP3** | `Nlrp3` 缺失阻止肥胖诱导的脂肪 / 肝脏炎症小体活化，改善胰岛素信号 | Vandanmagsar 等，*Nat Med* 2011 [47] | 高 |
| **JNK1（造血细胞，BMT）** | 造血细胞缺 JNK1 不减肥胖，但减轻 HFD 炎症与胰岛素抵抗 | Solinas 等，*Cell Metab* 2007 [16] | **有争议** |
| **JNK1（造血细胞，BMT，反证）** | 骨髓源 JNK1 缺失**不足以**改变巨噬细胞浸润或全身胰岛素敏感性；实质细胞缺 JNK1 才有效 | Vallerie 等，*PLoS One* 2008 [58] | 与 [16] 直接矛盾 |
| **巨噬细胞 JNK1+JNK2** | 髓系同时去除 JNK1 与 JNK2 后，肥胖相似但胰岛素敏感 | Han 等，*Science* 2013 [57] | 高（但是**双亚型**，不能单归因于 JNK1） |

**正确表述：**

- JNK1 驱动肥胖相关代谢炎症与胰岛素抵抗：**总体上成立** [56]。
- 「造血细胞 / 巨噬细胞内的 **JNK1 单独**」是否足够：Solinas 2007 支持 [16]，Vallerie 2008 不支持 [58]；Han 2013 表明需 **JNK1+JNK2 联合缺失** 才在巨噬细胞水平看到稳定保护 [57]。
- 因此**不能**把「两靶在同一巨噬细胞内耦合」写成该适应症的已证机制；只能写：NLRP3 在肥胖炎症中必要 [47]；JNK1 在全身代谢炎症中必要 [56]；细胞类型是否完全重合仍有争议。

JNK1 另有独立机制（如 IRS-1 丝氨酸磷酸化、拮抗胰岛素信号）[56]。这反而支持双靶覆盖**两条不完全重叠的轴**，但不能把获益全部算成 NLRP3 priming。

### 4.2 NASH / 脂肪性肝炎与肝纤维化

两靶都在饮食性肝炎模型中被验证。NLRP3 端已有 Cre-lox 髓系特异性证据；JNK1 造血细胞证据主要来自骨髓移植，证据等级不对等。

| 靶点 | 机制实验 | 文献 | 方法学备注 |
|------|----------|------|------------|
| **JNK1（不是 JNK2）** | MCD 模型：`jnk1` 而非 `jnk2` 缺失减轻脂肪性肝炎 | Schattenberg 等，*Hepatology* 2006 [59] | 全身敲除 |
| **JNK1 vs JNK2（HFD）** | HFD：JNK1 促进脂肪变与肝炎；JNK2 抑制肝细胞死亡 | Singh 等，*Hepatology* 2009 [67] | 全身敲除 + ASO |
| **JNK1（造血 / Kupffer）** | CDAA：造血细胞缺 JNK1 减轻肝炎与纤维化，脂肪变相似 | Kodama 等，*Gastroenterology* 2009 [60] | **BMT**；尚缺乏同模型 Cre-lox 独立复现 |
| **肝细胞 JNK1（脂毒性）** | 游离胆固醇脂毒性经 JNK1 介导的线粒体损伤；`Jnk1−/−` 肝细胞抵抗 | Gan 等，*J Hepatol* 2014 [61] | 体外肝细胞 + 遗传；**不经过 NLRP3** |
| **NLRP3 遗传** | `Nlrp3` 缺失保护 CDAA 肝损伤、巨噬细胞浸润与纤维化 | Wree 等，*J Mol Med* 2014 [62] | 全身敲除 |
| **NLRP3 细胞类型** | 髓系而非肝细胞 / HSC 特异性敲除 NLRP3，减轻 CDAHFD / Western diet 炎症与纤维化 | Kaufmann 等，*Cell Mol Gastroenterol Hepatol* 2022 [64] | **Cre-lox**，证据较强 |
| **NLRP3 药理** | MCC950 减轻肥胖糖尿病 NASH 及 MCD 小鼠的肝炎与纤维化；胆固醇结晶激活 Kupffer NLRP3 | Mridha 等，*J Hepatol* 2017 [63] | 药理 |

**可写的句子：**  
在 MCD / CDAA 脂肪性肝炎中，JNK1 [59,60] 与髓系 NLRP3 [62–64] 都被证明驱动炎症向纤维化进展。双靶的合理细胞落点是 **Kupffer / 巨噬细胞**，但 JNK1 的造血细胞定位目前主要依赖 BMT [60]。

**必须写明的缺口：**

- 没有 NASH 模型的 NLRP3-S194A 数据。
- 肝细胞 JNK1 有**不经过 NLRP3** 的线粒体损伤通路 [61]。NASH 中抑制 JNK1 的获益不能全部算成 NLRP3 priming。
- 另有文献表明炎症小体缺陷可通过菌群失调**加重** NAFLD [66]。这不否定髓系 NLRP3 促纤维化 [64]，但说明适应症叙事必须限定在**髓系执行端**，不能写成「抑制 NLRP3 在所有 NAFLD 模型中一律有益」。

---

## 5. 第三档：不能作为双靶适应症主张

### 5.1 IBD / DSS 结肠炎 — 不推荐

| 问题 | 证据 |
|------|------|
| **JNK1 不是致病必需；上皮 JNK2 保护屏障** | 非造血 **JNK2 而非 JNK1** 保护 DSS 结肠炎（减少上皮凋亡与屏障破坏） | Mandić 等，*Mucosal Immunol* 2017 [65] |
| **NLRP3 在代谢 / 肠稳态中可呈保护性** | NLRP3 / NLRP6 炎症小体与 IL-18 负向调节 NAFLD/NASH 进展（菌群依赖） | Henao-Mejia 等，*Nature* 2012 [66] |

Song 2017 **没有**做 DSS / IBD。把 IBD 写成双靶适应症，与现有 JNK 亚型遗传学冲突：非选择性抑制 JNK 可能伤及保护性 JNK2 [65]。

### 5.2 神经退行 / 中枢神经炎症 — 现阶段不能主张双靶

- NLRP3 在部分神经炎症模型中有单靶证据，但这不是 JNK1–NLRP3-S194 闭环。
- 神经元损伤更多涉及 **JNK3**，不是巨噬细胞 priming 轴。
- **没有** S194A 或 `Jnk1` 敲除与 `Nlrp3` 敲除在 AD/PD 中的平行验证。
- 另需血脑屏障暴露，与当前大分子双靶先导的 PK 不匹配（成药问题，不是机制证据）。

### 5.3 类风湿关节炎 — 证据层级不够

滑膜中 MAPK/JNK 与 NLRP3 均有表达或药理报道，但缺少与痛风、HFD、NASH 同等级的：`Jnk1` KO 与 `Nlrp3` KO 在同一关节炎模型中的平行遗传证据，也没有 S194A 数据。暂不列入机制闭合适应症。

---

## 6. 按机制闭合程度排序（供课题使用）

```
闭合程度
高  │  痛风 / MSU 晶体炎症     Song 2017 S194A + Martinon 2006 + Okada 2014
    │  LPS 内毒素血症          Song 2017 S194A（机制 POC，非临床主适应症）
中  │  NASH / 肝纤维化         JNK1：全身 KO 扎实；造血细胞为 BMT
    │                          NLRP3：全身 KO + 髓系 Cre-lox + MCC950
    │  肥胖–胰岛素抵抗         NLRP3 扎实；JNK1 全身扎实；
    │                          髓系 JNK1 文献互相矛盾，不得写成定论
低  │  IBD、CNS、RA            亚型冲突、单靶证据或无 S194 数据 → 不主张
```

**课题表述建议：**

1. **机制验证适应症（体内 POC）：** MSU 晶体性腹膜炎 / 关节炎 —— 唯一把 JNK1–NLRP3 偶联做到疾病模型的。
2. **转化适应症（两靶均有疾病遗传证据，需在正文声明缺口）：** NASH（优先于「巨噬细胞 JNK1」尚未闭合的单纯 HFD 胰岛素抵抗）。
3. **不写进双靶适应症：** IBD、神经退行、RA（除非补做对应遗传实验）。

组织暴露、口服 PK、分子量是**成药性**问题，不能用来把未验证的疾病「升级」为机制已闭合适应症。

---

## 7. 文献核对记录（相对早期草稿）

以下错误已在本文与 `REFERENCES.md` 中更正，请勿沿用旧 PMID / 作者 / 期刊：

| 错误 | 更正 |
|------|------|
| Song 2017 写作 Song H | 第一作者为 **Nan Song** [8] |
| Solinas 2007 写作 *J Clin Invest*，PMID 17549254 | 正确为 ***Cell Metab.* 2007;6(5):386-397**，PMID **17983584** [16] |
| Okada 2014 PMID 25271163 | 正确 PMID **25288801** [55] |
| Schattenberg 2006 PMID 16496325 | 正确 PMID **16374858** [59] |
| Kodama 2009 PMID 19549524 | 正确 PMID **19549522** [60] |
| Gan 2014 PMID 25060693 | 正确 PMID **25064435** [61] |
| Kaufmann 2022 写作 *J Hepatol* | 正确期刊 ***Cell Mol Gastroenterol Hepatol.*** [64] |
| IBD 机制文写作「Chorus 等」 | 正确第一作者 **Mandić AD** [65] |
| 将 Solinas 2007 的造血细胞 JNK1 写成无争议定论 | 必须并列 Vallerie 2008 [58] 与 Han 2013 [57] |
| 将 IBD / RA / 神经炎症列入双靶优先适应症 | 机制闭合不足，降至不主张 |

---

## 本文参考文献

编号 [1][8][47] 见 `REFERENCES.md`。下列为适应症专文新增或更正后的条目，亦已写入 `REFERENCES.md`。

**[16]** Solinas G, Vilcu C, Neels JG, et al. JNK1 in hematopoietically derived cells contributes to diet-induced inflammation and insulin resistance without affecting obesity. *Cell Metab.* 2007;6(5):386-397. doi:[10.1016/j.cmet.2007.09.011](https://doi.org/10.1016/j.cmet.2007.09.011). PMID: [17983584](https://pubmed.ncbi.nlm.nih.gov/17983584/).

**[54]** Martinon F, Pétrilli V, Mayor A, Tardivel A, Tschopp J. Gout-associated uric acid crystals activate the NALP3 inflammasome. *Nature.* 2006;440(7081):237-241. doi:[10.1038/nature04516](https://doi.org/10.1038/nature04516). PMID: [16407889](https://pubmed.ncbi.nlm.nih.gov/16407889/).

**[55]** Okada M, Matsuzawa A, Yoshimura A, Ichijo H. The lysosome rupture-activated TAK1-JNK pathway regulates NLRP3 inflammasome activation. *J Biol Chem.* 2014;289(47):32926-32936. doi:[10.1074/jbc.M114.579961](https://doi.org/10.1074/jbc.M114.579961). PMID: [25288801](https://pubmed.ncbi.nlm.nih.gov/25288801/). PMCID: [PMC4239639](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239639/).

**[56]** Hirosumi J, Tuncman G, Chang L, et al. A central role for JNK in obesity and insulin resistance. *Nature.* 2002;420(6913):333-336. doi:[10.1038/nature01137](https://doi.org/10.1038/nature01137). PMID: [12447443](https://pubmed.ncbi.nlm.nih.gov/12447443/).

**[57]** Han MS, Jung DY, Morel C, et al. JNK expression by macrophages promotes obesity-induced insulin resistance and inflammation. *Science.* 2013;339(6116):218-222. doi:[10.1126/science.1227568](https://doi.org/10.1126/science.1227568). PMID: [23223452](https://pubmed.ncbi.nlm.nih.gov/23223452/). PMCID: [PMC3835653](https://pmc.ncbi.nlm.nih.gov/articles/PMC3835653/).

**[58]** Vallerie SN, Furuhashi M, Fucho R, Hotamisligil GS. A predominant role for parenchymal c-Jun amino terminal kinase (JNK) in the regulation of systemic insulin sensitivity. *PLoS One.* 2008;3(9):e3151. doi:[10.1371/journal.pone.0003151](https://doi.org/10.1371/journal.pone.0003151). PMID: [18773088](https://pubmed.ncbi.nlm.nih.gov/18773088/).

**[59]** Schattenberg JM, Singh R, Wang Y, et al. JNK1 but not JNK2 promotes the development of steatohepatitis in mice. *Hepatology.* 2006;43(1):163-172. doi:[10.1002/hep.20999](https://doi.org/10.1002/hep.20999). PMID: [16374858](https://pubmed.ncbi.nlm.nih.gov/16374858/).

**[60]** Kodama Y, Kisseleva T, Iwaisako K, et al. c-Jun N-terminal kinase-1 from hematopoietic cells mediates progression from hepatic steatosis to steatohepatitis and fibrosis in mice. *Gastroenterology.* 2009;137(4):1467-1477.e5. doi:[10.1053/j.gastro.2009.06.045](https://doi.org/10.1053/j.gastro.2009.06.045). PMID: [19549522](https://pubmed.ncbi.nlm.nih.gov/19549522/). PMCID: [PMC2757473](https://pmc.ncbi.nlm.nih.gov/articles/PMC2757473/).

**[61]** Gan LT, Van Rooyen DM, Koina ME, et al. Hepatocyte free cholesterol lipotoxicity results from JNK1-mediated mitochondrial injury and is HMGB1 and TLR4-dependent. *J Hepatol.* 2014;61(6):1376-1384. doi:[10.1016/j.jhep.2014.07.024](https://doi.org/10.1016/j.jhep.2014.07.024). PMID: [25064435](https://pubmed.ncbi.nlm.nih.gov/25064435/).

**[62]** Wree A, McGeough MD, Peña CA, et al. NLRP3 inflammasome activation is required for fibrosis development in NAFLD. *J Mol Med (Berl).* 2014;92(10):1069-1082. doi:[10.1007/s00109-014-1170-1](https://doi.org/10.1007/s00109-014-1170-1). PMID: [24804980](https://pubmed.ncbi.nlm.nih.gov/24804980/).

**[63]** Mridha AR, Wree A, Robertson AAB, et al. NLRP3 inflammasome blockade reduces liver inflammation and fibrosis in experimental NASH in mice. *J Hepatol.* 2017;66(5):1037-1046. doi:[10.1016/j.jhep.2017.01.022](https://doi.org/10.1016/j.jhep.2017.01.022). PMID: [28167322](https://pubmed.ncbi.nlm.nih.gov/28167322/). PMCID: [PMC6536116](https://pmc.ncbi.nlm.nih.gov/articles/PMC6536116/).

**[64]** Kaufmann B, Kui L, Reca A, et al. Cell-specific deletion of NLRP3 inflammasome identifies myeloid cells as key drivers of liver inflammation and fibrosis in murine steatohepatitis. *Cell Mol Gastroenterol Hepatol.* 2022;14(4):751-767. doi:[10.1016/j.jcmgh.2022.06.007](https://doi.org/10.1016/j.jcmgh.2022.06.007). PMID: [35787975](https://pubmed.ncbi.nlm.nih.gov/35787975/).

**[65]** Mandić AD, Bennek E, Verdier J, et al. c-Jun N-terminal kinase 2 promotes enterocyte survival and goblet cell differentiation in the inflamed intestine. *Mucosal Immunol.* 2017;10(5):1211-1223. doi:[10.1038/mi.2016.125](https://doi.org/10.1038/mi.2016.125). PMID: [28098247](https://pubmed.ncbi.nlm.nih.gov/28098247/).

**[66]** Henao-Mejia J, Elinav E, Jin C, et al. Inflammasome-mediated dysbiosis regulates progression of NAFLD and obesity. *Nature.* 2012;482(7384):179-185. doi:[10.1038/nature10809](https://doi.org/10.1038/nature10809). PMID: [22297845](https://pubmed.ncbi.nlm.nih.gov/22297845/).

**[67]** Singh R, Wang Y, Xiang Y, Tanaka KE, Gaarde WA, Czaja MJ. Differential effects of JNK1 and JNK2 inhibition on murine steatohepatitis and insulin resistance. *Hepatology.* 2009;49(1):87-96. doi:[10.1002/hep.22578](https://doi.org/10.1002/hep.22578). PMID: [19053047](https://pubmed.ncbi.nlm.nih.gov/19053047/).
