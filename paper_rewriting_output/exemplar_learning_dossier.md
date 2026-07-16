# Exemplar Learning Dossier

**Tier:** pro（目标场景样例 + 领域 SOTA 样例）  
**Purpose:** 学结构与修辞，不抄结果。

## Exemplar Inventory

| ID | Paper | Year | Role | What to Learn | What NOT to Copy Blindly |
|----|-------|------|------|---------------|--------------------------|
| E1 | Lu et al., JMC — YL5084/YL2056 | 2023 | 目标场景强样例 | Intro：口袋同源→共价+动力学选择性；Results：kinact/KI、共晶、C116；Discussion：诚实写 off-target / JNK2-independent phenotype | 不要复制 MM.1S 成药叙事 |
| E2 | Wydra et al., JMC — 56d ligand-first | 2025 | 目标场景强样例 | “先可逆 Leu106 fit，再共价”的论证链；linker 几何（meta vs para） | 不要假装已有 56d 级 kinome 数据 |
| E3 | Zhang et al., Chem Biol — JNK-IN-8 | 2012 | 领域奠基 | Cys116 可药化、pan-JNK 共价探针范式 | 不是亚型选择性故事 |
| E4 | Shamir et al., COValid/AF3 | 2025–26 | 方法学样例 | decoy 设计（性质匹配+拓扑不相似）；mPAE 排序；校准 vs 前瞻分层 | 不要把 BTK 前瞻成功率外推到 JNK2 未校准大库 |
| E5 | London/Shoichet DOCKovalent | 2014+ | 库筛样例 | 分 warhead 商业库；高 hit rate 来自库匹配+位点清晰 | 大库需要服务器/集群 |
| E6 | Biology 14, 2025 — stressed epithelium / isoform specificity | 2025 | 叙事样例 | 明确写出：**DSS 等模型缺 isoform-specific 工具，故无法判断 JNK1 vs JNK2** | 不要据此宣称 JNK2 已是 IBD 主靶 |
| E7 | CC-90001 clinical/IPF design papers | 2021– | 对照叙事 | JNK1-bias 临床化合物；纤维化叙事属于 **JNK1 轴** | 勿挪用为 JNK2 适应症证据 |
| E8 | Enamine Cys/acrylamide library notes | 2023–24 | 库工程样例 | 固定 warhead、去过反应、Ro3/Ro5、可采购优先 | 商业库≠靶点定制邻域库 |

## Structural Patterns

1. **Problem → Progress → Gap → Tool/Method → Evidence → Boundary**（Lu/Wydra 共有）
2. **化学贡献在前，疾病展望在后**（成功 JMC 共价探针文）
3. **方法学文：benchmark → enrichment → prospective case**（COValid）
4. **失败生物学结果被写成边界，而不是删掉**（Lu MM.1S）

## Rhetorical Patterns

- 用 **“probe / tool compound / kinetic selectivity”** 而不是 **“therapeutic for X disease”** 作主标题语气。
- 用 **kinact/KI + C116S + washout** 作为选择性金标准句式。
- 疾病段用 **“enables isoform-resolved questions in …”** 而不是 **“treats …”**。

## Language Patterns（中文报告可用）

- 「首个/唯一」仅在有文献支持时使用（YL5084：首个充分表征的 JNK2>JNK1 共价动力学选择性）。
- 「虚筛发现」改为「校准通过后的候选富集」。
- 「优势适应症」改为「尚缺亚型拆分工具的争议生物学场景」。
