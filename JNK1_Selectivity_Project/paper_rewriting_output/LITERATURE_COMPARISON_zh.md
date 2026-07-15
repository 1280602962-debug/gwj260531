# 相近文献对照（高质量期刊优先）

> **选刊原则**：主对照只用 **非 OA 或质量足够的期刊**（ACS *JCIM* / *J. Med. Chem.* / *ACS Med. Chem. Lett.*，Elsevier *Eur. J. Med. Chem.*）。  
> **不作为主对照**：低分 MDPI OA（*Molecules*、*Pharmaceuticals* 等）、*ChemistrySelect* 等应用型弱刊——可作“弱 hit 也有人发”的背景，**不要当作对标档位**。  
> 湿实验约束：仅 **JNK1/2/3 酶促 IC50**；阳性 E1、CC-90001。

---

## 0. 结论（按你的选刊标准）

| 用途 | 推荐主对照 | 期刊质量 |
|------|------------|----------|
| **方法/选择性叙事（最该学）** | Kinase-Bench (*JCIM* 2024)；JAK IFP-VS (*JCIM* 2016)；FEP 选择性 (*JCIM* 2020) | ACS 旗舰计算化学，hybrid/OA 但质量够 |
| **流程+商业库+酶活闭环** | RET ML+对接 VS (*JCIM* 2020)；CCR2/5 (*JCIM* 2025) | 同上 |
| **同靶点 JNK VS→酶活→选择性（天花板，勿对标结果）** | Dou 2019 / 吲唑 2022 (*J. Med. Chem.*) | ACS 药物化学旗舰，传统订阅/hybrid |
| **领域阳性与课题动机** | CC-90001 (*JMC* 2021)；E1 (*JMC* 2024) | 同上 |
| **降级：仅知“µM 也能发 OA”** | *Molecules* 2022 Tricin 等 | **不引用为对标** |

你课题若坚持高质量对照，**投稿也不宜以 *Molecules*/*Pharmaceuticals* 为心理锚点**；更合理目标是 *JCIM*（方法+负结果基准）、*ChemMedChem* / *ACS Med. Chem. Lett.*（若出活性）、或 *Eur. J. Med. Chem.*（若后续有化学优化）。

---

## 1. 主对照文献（推荐全部精读）

### A. 选择性 / 方法学（与你的“预测器失败”最贴）

#### [A1] *JCIM* 2024 — Kinase-Bench：激酶选择性 VS 基准 + 湿实验验证

- **链接**：https://doi.org/10.1021/acs.jcim.4c01830  
- **期刊**：*J. Chem. Inf. Model.*（ACS；计算化学一流）  
- **内容摘要**：针对 ATP 口袋高度保守、对接难做选择性，构建 75 激酶选择性配体/decoy 基准；Glide HTVS/SP + 多打分 + **定制蛋白–配体相互作用过滤**；前瞻筛 JAK1 vs TYK2。  
- **计算**：基准集构建；Glide；交互指纹/残基过滤；回顾+前瞻 VS。  
- **生物**：命中化合物酶活——如 Compound 2：JAK1 IC50 ≈ **980 nM**，TYK2 ≈ **4.5 µM**；另报 Capivasertib 的 JAK1/TYK2 差。  
- **与你**：  
  - **科学问题最近**：对接/交互规则能否买到 isoform/family 选择性。  
  - 他们做成了**正结果基准工具**；你已有 **JNK1/2/3 上 Δsel/Gly87/ML 的量化负结果**——可对标为“同一问题的另一面”。  
  - 生物上他们做双激酶 IC50；你可做 **三亚型**，但缺他们那种大基准集叙事时要写清 scope。

---

#### [A2] *JCIM* 2016 — 亚型选择性 JAK 抑制剂的结构 VS

- **链接**：https://doi.org/10.1021/acs.jcim.5b00634  
- **期刊**：*JCIM*  
- **内容摘要**：蛋白 ensemble 对接 + **Interaction Fingerprint (IFP)** 估亲和与选择性；回顾验证后前瞻筛商业库；6 个化合物体外确认，命中率约 **11%**；吲唑类显示 JAK2 vs JAK1 偏好。  
- **计算**：ensemble docking；IFP 选择性打分；回顾/前瞻。  
- **生物**：少量购买化合物的 JAK 亚型酶活确认。  
- **与你**：你也用了对接差/IFP 类思路但 **失败**——引用此文时强调：IFP 在 JAK 上可工作，在 **JNK 近同源口袋上未通过你的校准**，这本身是可发表信息。

---

#### [A3] *JCIM* 2020 — SBDD 是否准备好做选择性优化？（FEP）

- **链接**：https://doi.org/10.1021/acs.jcim.0c00815  
- **预印本**：https://doi.org/10.1101/2020.07.02.185132  
- **期刊**：*JCIM*（无新化合物湿实验）  
- **计算**：CDK2/CDK9、CDK2/ERK2 上 alchemical FEP；Bayesian 分解系统/统计误差；讨论选择性预测是否比亲和更准。  
- **生物**：无新测。  
- **与你**：Discussion 级背书——即便 FEP 都要专门检验选择性；你的对接/ML 负结果是更“筛选层”的补充，勿宣称已做 FEP。

---

### B. 管线形态（ML/对接 → 购买 → 酶活）

#### [B1] *JCIM* 2020 — 统计模型 + 结构模型接力找化学新颖激酶抑制剂（RET 等）

- **链接**：https://doi.org/10.1021/acs.jcim.9b01204  
- **PMC**：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7525794/  
- **期刊**：*JCIM*  
- **内容摘要**：QSAR/PCM 与对接/SPLIF 等结构方法系统基准后接力 VS；购买 **46** 个测 RET；单点 10 / 0.1 µM，再对阳性做 IC50；命中偏结构/共识打分而非纯统计。  
- **计算**：多激酶统计模型；对接；共识打分；部分 Binding Pose Metadynamics。  
- **生物**：放射性激酶酶活；5 个在 10 µM 有抑制等。  
- **与你**：  
  - **流程最近的高质量对照**：ML→结构→买→酶活，且诚实比较哪种打分带来 hit。  
  - 他们 n_buy=46；你 n_new≈2 ——功率弱，须把贡献放在 **方法校准/失败**，不是 hit 率统计。  
  - 他们测单靶；你测 **JNK1/2/3** 是差异化。

---

#### [B2] *JCIM* 2025 — AlphaFold + 聚焦库：CCR2/CCR5（µM、选择性有限）

- **链接**：https://doi.org/10.1021/acs.jcim.5c01596  
- **期刊**：*JCIM*（ACS；质量足够，可采用）  
- **内容摘要**：组合库对接 AF 精修口袋；合成 18 个；CCR2 Ki **1.3–6 µM**（3 hits），CCR5 IC50 **10.8 µM**（1 hit）；明确写 hit 率约 5–15%、选择性有限。  
- **计算**：库生成；AF 模型精修；对接；结构解释。  
- **生物**：结合 + 功能双靶。  
- **与你**：**措辞与期望管理**最好的高质量范文——µM、低命中、不硬吹选择性。你缺合成/双 assay 类型，但酶活三亚型可替代“双靶对照”角色。

---

### C. 同靶点 JNK（高质量；结果天花板，用于定位而非对标活性）

#### [C1] *J. Med. Chem.* 2019 — 多阶段筛选：吲哚啉酮类 JNK3 异构体选择性抑制剂

- **链接**：https://doi.org/10.1021/acs.jmedchem.9b00537  
- **期刊**：*J. Med. Chem.*（药物化学旗舰）  
- **内容摘要**：比较性结构 VS → 命中 → SAR；**J30-8** JNK3 IC50 ≈ **40 nM**，对 JNK1/2 **>2500×**；kinome 选择性；细胞神经保护 + 体内 AD 相关药效。  
- **计算**：多阶段/比较对接 VS；后续结构优化。  
- **生物**：JNK1/2/3 酶活；激酶 panel；体外/体内。  
- **与你**：同“VS→三亚型酶活”起点，但后续有 **合成 SAR + kinome + 体内**。你只能停在酶活矩阵——文中定位为 **early enrichment / method audit**，不是 tool compound 发现。

---

#### [C2] *J. Med. Chem.* 2022 — 吲唑类 JNK3 异构体选择性（PD）

- **链接**：https://doi.org/10.1021/acs.jmedchem.2c01410  
- **期刊**：*J. Med. Chem.*  
- **内容摘要**：对接 VS + SAR → **25c** JNK3 IC50 ≈ **85 nM**，对 JNK1/2 **>100×**；376 激酶 panel；细胞/体内 PD 模型；BBB。  
- **与你**：领域天花板；只引“异构体选择性极难、成功案例长什么样”，**不对标实验深度**。

---

#### [C3] *ACS Med. Chem. Lett.* 2020 — 噻吩–吡唑脲 JNK3 选择性 + 共晶

- **链接**：https://doi.org/10.1021/acsmedchemlett.0c00533  
- **期刊**：*ACS Med. Chem. Lett.*（ACS，质量好）  
- **内容摘要**：SAR + 酶/细胞 + DMPK；化合物 17：JNK3 IC50 ≈ **35 nM**，374 激酶几乎只打 JNK3；共晶解释 hinge/疏水口袋。  
- **与你**：说明真正选择性常靠 **化学优化+结构**，不是一轮商业库对接。

---

#### [C4] *Eur. J. Med. Chem.* 2020 — 喹喔啉酮 JNK3：从 VS hit 到选择性优化

- **链接**：https://doi.org/10.1016/j.ejmech.2020.112416（文题：Rational modification… dihydroquinoxalin-2-one… JNK3）  
- **期刊**：*European Journal of Medicinal Chemistry*（Elsevier 主流药化，非低分 OA）  
- **内容摘要**：商业库 VS 得到 J46（JNK3 IC50 ≈ 0.25 µM，已有亚型差，但 DDR1/EGFR 差、物化差）→ 合成优化 → **J46-37** 强效且对 JNK1/2 **>500×**，改善脱靶与物化。  
- **与你**：清晰画出你现在所处阶段 = **他们的 J46 之前/之初**；没有合成就不要写成他们的终态。

---

### D. 领域动机与阳性对照出处（必引，非 VS 流程对标）

#### [D1] *J. Med. Chem.* 2021 — CC-90001 发现

- **链接**：https://doi.org/10.1021/acs.jmedchem.1c01716  
- **内容**：JNK1 bias 临床候选；nM 酶活、kinome、细胞 p-c-Jun、PK/体内。  
- **你的用法**：阳性对照与生物学动机；**不是**你的 VS 流程模板。

#### [D2] *J. Med. Chem.* 2024 — E1（JNK1，IPF）

- **链接**：https://doi.org/10.1021/acs.jmedchem.4c01764  
- **内容**：SBDD+SAR；E1 JNK1 IC50 ≈ **2.7 nM**（文献亦报 JNK2/3 nM 级）；MD；PK；纤维化模型。  
- **你的用法**：同批阳性；承认你的新分子预期远弱于 E1。

---

## 2. 总表（仅高质量主对照）

| ID | 论文 | 期刊 | 计算核心 | 生物核心 | 与你距离 |
|----|------|------|----------|----------|----------|
| A1 | Kinase-Bench 2024 | **JCIM** | 选择性基准+Glide+交互过滤 | JAK1/TYK2 IC50（µM–nM） | **方法问题★★★★★** |
| A2 | JAK IFP-VS 2016 | **JCIM** | ensemble+IFP 选择性 | 6 化合物亚型酶活 | 选择性打分★★★★★ |
| A3 | FEP selectivity 2020 | **JCIM** | FEP 误差分析 | 无新测 | Discussion★★★★ |
| B1 | RET ML+SBVS 2020 | **JCIM** | QSAR+对接共识 | 46 买→酶活 | **管线形态★★★★★** |
| B2 | CCR2/5 2025 | **JCIM** | AF+对接 | µM 双靶，诚实选择性 | 期望管理★★★★★ |
| C1 | Dou JNK3 2019 | **JMC** | 多阶段 VS | nM + >2500× + 体内 | 同靶天花板★★★ |
| C2 | 吲唑 JNK3 2022 | **JMC** | VS+SAR | nM + >100× + panel | 天花板★★ |
| C3 | 噻吩脲 JNK3 2020 | **ACS MCL** | SAR+共晶 | nM + kinome | 结构生物学★★ |
| C4 | 喹喔啉 EJMC 2020 | **EJMC** | VS hit→合成 | 0.25 µM→高选择性 | 阶段定位★★★★ |
| D1/D2 | CC-90001 / E1 | **JMC** | SAR/SBDD | 临床/先导全套 | 阳性/动机★★★★★ |

---

## 3. 明确降级、不作为对标的文献

| 论文 | 为何降级 | 若提及怎么用 |
|------|----------|--------------|
| *Molecules* 2022 Tricin（MDPI OA） | 低分 OA；仅 JNK1、µM | 可一句带过“弱酶活闭环先例”，**不作档位锚点** |
| *Pharmaceuticals* 2023 JNK3 DL-VS（MDPI OA） | 同上 | 可知 DL 重打分管线存在；不对标 |
| *ChemistrySelect* PIM1 / Haspin ChemBridge | Wiley 应用型弱刊 | 证明商业库买 10 个出 µM——档位不够 |
| *Molecules* 2025 BCR-ABL | MDPI OA | 同上 |

---

## 4. 与高质量文相比：你缺什么、能补什么

### 4.1 生物（只做 JNK1/2/3 酶活时）

| 高质量文常见 | 你 | 处理 |
|--------------|----|------|
| 目标±近缘亚型 IC50 | **能做三亚型** | 主实验；预注册 SI |
| 买几十个算 hit 率（B1） | n≈2 | **不要报 hit 率**；报校准/enrichment 个案 |
| Kinome panel（C1/C2/C3） | 不能 | Limitation + 计算脱靶推断 |
| 细胞 / 体内（C/D） | 不能 | 不对标；引用 D1/D2 作领域背景 |
| 双靶功能 assay（B2） | 不能 | 用 JNK1/2/3 矩阵代替“选择性读出” |

### 4.2 计算（相对 A1/A2/B1 仍建议补）

| ID | 缺口（审稿人对照 JCIM 会问） | 补法 | 优先级 |
|----|------------------------------|------|--------|
| **C1** | chemotype 相对已知 JNK 是否新颖 | ChEMBL maxTc / Murcko vs SP600125/E1/CC-90001/CC-930 | P0 |
| **C2** | pose 可重复性 | ≥3 seed 重对接 | P0 |
| **C3** | 单轨迹 MD | 购买集 × 3 亚型 × ≥2 seed | P0 |
| **C4** | 终点预注册 | 主：任亚型 IC50；次：SI≥3 vs JNK2 **且** JNK3 | P0 |
| **C5** | 负结果是否主文化（对标 A1/A2） | 主文表：Δsel 方向准确率、Gly87、ML F1 | **P0（差异化）** |
| C6 | 有无类似 Kinase-Bench 的回顾基准 | 用你已有文献 JNK 集整理回顾表（不必新建 75 激酶） | P1 |
| C7 | 脱靶 | SEA/SwissTargetPrediction（假设） | P1 |
| C8 | FEP（A3） | 不强制；引用即可 | P2 |

### 4.3 你相对这些高质量文仍可主张的独特点

1. **JNK1/2/3 上多种“便宜选择性过滤器”的联合负校准**（对接 Δ、Gly87、ML）——A1/A2 偏正结果/工具，你补 **失败案例**。  
2. **MD QC 与酶选择性解耦**（E1 / SP600125 反例）——多数 JMC 成功文不会强调。  
3. **购买决策与选择性打分解耦**的诚实管线——对齐 B2 的语气，内容对齐 A1 的问题。

---

## 5. 写作引用策略（高质量档）

1. **Introduction 动机**：D1（CC-90001）、D2（E1）；选择性难 → C1/C2 一句。  
2. **方法学问题**：A1、A2、A3。  
3. **管线合理性**：B1、B2。  
4. **结果期望**：明确写——本研究停留在 B1/B2 的 **early biochemical confirmation** 层，不达到 C1–C3 的 tool-compound 层。  
5. **禁止**：用 *Molecules*/*Pharmaceuticals* 当“我们达到同类水平”的依据。

---

## 6. 半天阅读顺序（只含高质量）

1. **A1** Kinase-Bench（全文）  
2. **A2** JAK IFP-VS（方法+命中表）  
3. **B1** RET ML+SBVS（购买与酶活段）  
4. **B2** CCR2/5（Conclusion 措辞）  
5. **C1** Dou 2019（VS→三亚型表，知天花板）  
6. **D1/D2** 扫摘要即可  
