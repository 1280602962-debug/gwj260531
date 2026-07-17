# ARS Academic Paper Reviewer + PaperSpine Logic Diagnosis

**Mode:** full（condensed；仅 Results & Discussion 可用）  
**Manuscript:** `/tmp/manuscript_polished_en.md`（Polished EN R&D extract）  
**Prior scientific flags:** `/workspace/JNK1_Selectivity_Project/1_Results_and_Discussion_translated_checked.md` Part 1  
**PaperSpine mode:** draft-only → **logic diagnosis + objection register**（无修订稿，无法做 full rewrite transfer）  
**Date:** 2026-07-17  

**Caveat（全局）:** 本稿仅有 Results & Discussion；无 Abstract / Introduction / Methods / References / Supplementary 全文。评分与决定仅针对已呈现内容；缺失章节记为 **submission blockers**。不虚构 citation、实验或数值。

---

## 1. Field Analyst（学科与审稿人配置）

| 维度 | 判定 |
|------|------|
| **Primary discipline** | Computational medicinal chemistry / CADD（ligand-based + structure-based virtual screening） |
| **Secondary disciplines** | Cheminformatics & ML-QSAR；kinase–inflammasome multi-target discovery；MD / MM-GBSA binding analysis |
| **Research paradigm** | In silico discovery pipeline（data curation → QSAR → dual-target funnel → ADMET/MPO → MD） |
| **Methodology type** | Computational / predictive modeling + hierarchical docking + all-atom MD（无 wet-lab 结果于本稿） |
| **Paper maturity** | Mid-draft polished R&D；科学笔误部分已修，但仍有结构性 overclaim 与内部不一致；全文未齐 |
| **Suggested venue band** | Mid-specialty CADD / cheminformatics（如 *JCIM*、*J. Chem. Inf. Model.* 类；*Molecules* / *Comput. Biol. Chem.* 等）。无实验验证时不宜冲击高影响力综合期刊 |
| **Maturity vs claim** | 技术流程完整度尚可；“dual-target inhibitor / micromolar potency” 类措辞超前于证据 |

### Suggested reviewer personas（5）

| # | Role | Persona |
|---|------|---------|
| 1 | **EIC** | Editor, *Journal of Chemical Information and Modeling*（或同等 CADD 专刊）；关注 workflow novelty、claim–evidence 匹配、读者可复现性 |
| 2 | **Methodology** | QSAR / ML validation specialist（scaffold split、AD、Y-scrambling、UQ calibration、enrichment statistics） |
| 3 | **Domain** | Structural pharmacology of PLK1 hinge / NLRP3 NACHT（PDB 编号、hinge Cys、口袋几何、reference ligand 公平比较） |
| 4 | **Perspective** | Translational multi-target oncology–inflammation discovery（MPO、hERG、合成可及性、从 in silico 到 assay 的路径） |
| 5 | **Devil's Advocate** | Overclaim / logic-gap stress tester（“inhibitor” 语言、enrichment 定义、ESM-2 是否真正约束筛选、对接分数≠活性） |

---

## 2. Rubric Scores（仅评已有 R&D）+ Decision Mapping

**Calibration note:** 无校准集；分数作序数比较用，不作基数“录用保证”。缺全文会系统性压低 Evidence / Coherence。

| Dimension | Weight | Score (0–100) | Descriptor | Rationale（简） |
|-----------|--------|---------------|------------|----------------|
| Originality | 20% | **68** | Adequate→Strong 边缘 | PLK1–NLRP3 dual-target + UQ 过滤 + multi-anchor NLRP3 策略有一定新意；整体仍是标准 VS funnel 的增量组合 |
| Methodological Rigor | 25% | **52** | Weak→Adequate | Scaffold split / Y-scramble / UQ 方向正确；但 Williams \(p=270\) 与 50-PC 模型矛盾、enrichment/intersection 统计 framing 脆弱、MPO 权重和不归一、ESM-2 与打分链路脱节 |
| Evidence Sufficiency | 25% | **48** | Weak | 无 wet-lab；无 Intro/Methods/文献网；多处从 docking/MD 跃迁至 inhibitor / potency 结论；ESM UMAP 未证明筛选增益 |
| Argument Coherence | 15% | **55** | Weak→Adequate | 漏斗叙事可读，但多处 metric-dump；Cys133 vs Cys97 自相矛盾；147-fold enrichment、MCC950 “outperformed” 等 claim-leap |
| Writing Quality | 15% | **72** | Adequate→Strong | 英文化程度良好；仍有 “and and”、过度营销措辞、个别统计句子含混 |

**Weighted estimate:**

\[
0.20\times68 + 0.25\times52 + 0.25\times48 + 0.15\times55 + 0.15\times72 = \mathbf{57.65}
\]

| Mapping | Decision |
|---------|----------|
| ≥80 Accept / 65–79 Minor / 50–64 Major / <50 Reject | **Major Revision（条件性）** |

**Editorial decision mapping（带 caveat）:**

- 对**已呈现的 R&D**：落在 **Major Revision** 区间（≈58）。
- 对**投稿完备性**：缺 Abstract / Introduction / Methods / References / 完整 SI → **Incomplete submission / not review-ready**。
- Devil's Advocate 识别多项 **CRITICAL** overclaim → 按 ARS 规则 **不可 Accept**。
- 综合裁定：**Major Revision — Incomplete Manuscript**（先补齐全文并消化 Must-Fix，再进入正式同行评议）。

### Optional dimensions（叙事用，不入加权）

| Optional | Score | Note |
|----------|-------|------|
| Literature Integration | N/A（≈40 占位） | R&D 中几乎无文献定位；Intro 缺失 |
| Significance & Impact | ≈58 | 双靶概念有转化吸引力；无实验则 impact 上限明显 |

### Prior scientific flags → Polished EN status

| Prior flag | Polished EN status | Review action |
|------------|-------------------|---------------|
| ESM-2 480 vs 640 dims | **已修正为 640**（`esm2_t30_150M_UR50D`） | Nice-to-verify in Methods/logs；非 Must-Fix |
| Hindi `मैच` artifact | **未见残留** | Closed |
| MM/GBSA 46–50 ns vs 50–100 ns | **已写 50–100 ns** | Closed（Methods 需与脚本一致） |
| NLRP3 “浅表口袋” | **已改为 buried / NACHT cavity** | Closed |
| Cys97 vs Cys133 hinge | **仍矛盾**：§3 写 Cys133；§9 写 Cys97 为 2RKU hinge | **Must-Fix** |

---

## 3. Five Short Reviewer Reports

### 3.1 EIC Review（JCIM-class）

**Verdict on available material:** Interesting dual-target computational story, but **not submission-ready**; claims outrun evidence; full paper absent.

**Major**
1. 缺少 Abstract / Introduction / Methods — 无法判断 novelty framing、可复现协议与 journal fit。
2. 标题级承诺（dual-target **inhibitors**）与证据层级（docking + MD + predicted ADMET）不匹配；需全面降调为 *candidate / predicted binder*。
3. Enrichment（≈147-fold）与 “outperforming MCC950” 易被读者/审稿人视为 marketing；须定义 baseline 或删除。
4. 内部结构生物学不一致（Cys133 vs Cys97）损害信任。

**Minor**
1. MPO 选 Mol_997197（总分第 3）需更清晰的 decision rule，避免 post-hoc 印象。
2. 写作大体流畅；去掉夸饰形容词即可。

---

### 3.2 Methodology Reviewer

**Major**
1. **Williams AD：** \(h^*=3(p+1)/n\approx0.713\) 用 \(p=270\)，但冠军模型输入为 **50 PCA features** — \(p\) 定义与模型不一致，AD 结论不可信。
2. **Intersection “random expectation = 25%”：** \(500/1000\times500/1000\) 假设两靶 top-1000 集合相同且独立均匀抽样，与前文 “structurally distinct receptors” 叙事冲突；统计学 framing 错误或至少未证明。
3. **Enrichment ≈147-fold：** 未给出 null model / decoy rate / 何为 “success”；\(10000/68\approx147\) 仅是压缩比，**不是** enrichment factor。
4. **ESM-2：** 生成 640-d embeddings + UMAP 聚类，但筛选实际依赖 **hybrid ECFP4 + physchem similarity**；ESM 是否进入评分函数未证伪也未证实 → 方法贡献不清晰。
5. **UQ：** \(\sigma_i\) 与误差相关 \(r=0.33\) 偏弱却写 “highly significant / self-confidence awareness”；需校准图（reliability）与效应量表述。
6. **MPO 权重和 = 1.08：** 归一化失效；18 参数未完整列出。
7. **Paired t-test on 5 CV folds：** \(n=5\) 功效极低；应用嵌套 CV 或置换检验，并报告多重比较语境。

**Minor**
1. “and and PC1” 笔误。
2. Baseline 12 模型 “default hyperparameters” 对树模型不公平；需说明是否仅作粗筛。
3. H-bond occupancy 0.089 / 0.045 极低 — 需报告 cutoff、是否含水桥、与“稳定结合”叙事如何共存。

---

### 3.3 Domain Reviewer

**Major**
1. **Hinge 残基：** 领域共识 PLK1 hinge 多为 **Cys133**；本稿 §3 正确写 Cys133，§9 却将 **Cys97** 标为 2RKU hinge — 必须核对 PDB numbering / chain mapping，统一术语（hinge vs 邻近接触残基）。
2. **Docking score 比较 MCC950（−5.70）vs Mol_997197（−8.06）：** 不同 chemotype、同一 scoring function 下分数差 **不能** 直接解读为更优 NLRP3 抑制剂；MCC950 为已知活性对照，该比较易误导。
3. **MM/GBSA → “micromolar or sub-micromolar inhibitory potency”：** 绝对 MM/GBSA 值与实验 IC50 无可靠绝对换算；属 **claim-leap**。
4. **NLRP3 “anchor compatibility”：** 多锚点相似合理作 prefilter，但不能替代活性学习或结构药效团验证；需承认 chemotype bias。
5. **Residue numbering in decomposition（ALA45, ILE228 等）：** 需确认相对 7ALV 全长/截短编号，避免审稿人质疑 “wrong pocket”。

**Minor**
1. Imidazo[1,2-a]pyridine ↔ Cys133 叙事与 Mol_997197 实际支架是否一致需在姿态分析中明示。
2. hERG 高淘汰率解释合理，但“large conjugated aromatic scaffolds required” 略绝对化。

---

### 3.4 Perspective Reviewer

**Major**
1. 转化路径断层：从 68 → 41 → lead 全为计算；无 biochemical / cellular inflammasome 或 PLK1 assay 计划陈述（因 Intro/Discussion 闭环缺失而更明显）。
2. Dual-target 的生物学动机（肿瘤–炎症交叉）在 R&D 中几乎未回收 — Discussion 更像结果复述。
3. 安全性（hERG）已识别，但未讨论 polypharmacology / kinase selectivity panel — 对 PLK1 候选是关键 gap。

**Minor**
1. SA / Lipinski 高通过率是正面工程信号，可保留为 translational strength。
2. 开源工具链（Vina/Gnina/Amber）有利于可重复与社区采用 — 应在 Methods 强化。

---

### 3.5 Devil's Advocate

**Strongest counter-argument**  
本稿核心主张是“发现了有前景的 PLK1–NLRP3 **双靶抑制剂** lead（Mol_997197）”。然而：（i）NLRP3 臂几乎无监督活性模型，仅靠相似性 + docking；（ii）“优于 MCC950”主要来自 docking/MM-GBSA，而 MCC950 是已验证分子；（iii）PLK1 QSAR 虽较扎实，但不能外推双靶抑制功能。最强反对意见是：这是一次 **精心设计的计算优先排序练习**，不是抑制剂发现的证据闭环——在缺乏正交实验或至少严格 decoy/enrichment 基准时，主结论不成立。

**CRITICAL**
1. “Dual-target inhibitor / lock-and-key / inhibitory potency” 语言超证据。
2. 147-fold “enrichment” 定义错误或未定义。
3. ESM-2 被叙述为 pocket compatibility 的物理约束，但未显示进入决策函数 → 可能为装饰性方法。

**MAJOR**
1. Cys133 / Cys97 矛盾。
2. Williams \(p\) 与特征维度不一致。
3. Intersection 25% null 不成立。
4. MCC950 docking 劣势可能反映评分函数偏差而非真实亲和。
5. 极低 H-bond + 仍宣称高稳定性 — 需机制一致性解释。

**MINOR**
1. 选第 3 名 MPO 分子的叙事风险（cherry-picking 观感）。
2. \(r=0.33\) 被写成强 UQ 效用。

**Ignored alternatives**
- 单靶优化 + 后续 crosstalk 药理学，可能比强行 dual-pocket ligand 更可行。
- NLRP3 用 docking-only / pharmacophore / free-energy 绝对结合，而非 similarity pooling。
- 对 ChemDiv/Taosu 输出做 prospective assay 子集，而非仅 MD 两体系。

**“So what?” test**  
若删除 ESM-2 段与 147-fold 句，故事仍基本成立为“ML-QSAR + docking funnel + MD lead triage”——说明部分卖点未通过必要性测试。

---

## 4. Logic Diagnosis（PaperSpine adapted：draft-only）

### 4.1 Executive Verdict（logic diagnosis）

| Item | Content |
|------|---------|
| **Logic transfer verdict** | **N/A — no revised manuscript**；本报告产出 **original logic diagnosis + objection register** |
| **原稿试图做什么** | 用 PLK1 监督 QSAR + NLRP3 低数据相容性约束 + 多级对接/ADMET/MPO + 100 ns MD，论证 Mol_997197 为双靶 lead |
| **最强论证段** | §3.1–3.2 scaffold split + SVR 泛化与残差诊断（相对自洽） |
| **最弱论证段** | §3.5（ESM 角色）、§3.6（enrichment/intersection）、§3.8–3.9（从分数到 inhibitor / potency） |
| **缺陷标签汇总** | `metric-dump`, `claim-leap`, `method-recipe`（ESM）, `discussion-repeat`（收束段）, 局部 `gap-vague`（双靶生物学动机未在 R&D 回收） |

### 4.2 Whole-paper argument lineage（Results subsections）

| Subsection | Promise / question tested | Evidence | Interpretation | Implication / transition | Verdict |
|------------|---------------------------|----------|----------------|--------------------------|---------|
| **3.1 Dataset & chemical space** | 数据是否足以支撑可泛化 QSAR？ | 1426 mol；pIC50 分布；314 scaffolds；Tanimoto≈0.17；scaffold 8:1:1 | 化学空间分散、属性分布对齐 | → 可做严格外推基准 | **Pass**（叙事清晰） |
| **3.2 Baseline → champion SVR** | 非线性模型是否必要？SVR 是否显著更优？ | 12 模型；Optuna；CV MSE；配对 t；test \(R^2=0.74\) | 存在非线性 SAR；SVR 冠军 | → 可用于大规模 VS | **Pass with caution**（默认超参；t-test n=5） |
| **3.3 Interpretability PCA** | 模型是否学到与 PLK1 结构生物学一致的特征？ | Permutation importance；Morgan bits → imidazo[1,2-a]pyridine / Cys133 | 特征与 hinge/疏水口袋一致 | → 生物学合理性 | **Partial**（载荷→药效团映射需防过度解读） |
| **3.4 Robustness / AD / UQ** | 是否非偶然相关？UQ 能否降假阳性？ | Y-scramble；Williams；bootstrap \(r=0.33\)；排除 15% → MSE↓25.5%；ChemDiv 过滤 | 模型稳健；UQ 有用 | → 输出高置信 PLK1 候选 | **Fail elements**（\(p=270\)；UQ 措辞过强） |
| **3.5 NLRP3 ESM + anchors** | 无统一 QSAR 时如何约束 NLRP3？ | ESM 640-d UMAP；五锚点 hybrid similarity；0.5:0.5 融合至 10k | “anchor compatibility” 可预筛 | → 对接漏斗 | **Fail**（ESM 未证入打分；`method-recipe` + `claim-leap`） |
| **3.6 Docking funnel** | 双靶共优先能否富集？ | Vina→Gnina→MM/GBSA；68 交集；13.6% vs 25%；“147-fold” | 双靶难；漏斗高效 | → ADMET | **Fail**（统计与 enrichment 定义） |
| **3.7 ADMET cascade** | 68 个是否药学可用？ | MM/GBSA≤−50；Lipinski/SA/logS/HIA；hERG 淘汰 42% | hERG 主障碍 | → MPO | **Pass**（工程叙事好；阈值来源需 Methods） |
| **3.8 MPO → Mol_997197** | 谁是 lead？ | 加权 18 参数；雷达图；对接分 vs BI2536/MCC950 | 第 3 名但更平衡；“double lock-and-key” | → MD | **Fail elements**（权重和 1.08；MCC950 比较；选第 3 名需规则） |
| **3.9 MD / MM-GBSA** | 动态是否稳定？是否支持双靶？ | 100 ns RMSD/Rg/SASA/Hbond；ΔG；残基分解 | 稳定；可与对照比较；优化路径（Glu398） | → 合成与生测 | **Partial→Fail on claims**（稳定性证据尚可；potency/inhibitor 语言越界；Cys97） |

### 4.3 Shallow-patch / metric-dump / overclaim checklist

| Warning sign | Observed? | Evidence | Required fix |
|--------------|-----------|----------|--------------|
| Same paragraph order with sentence-level polish only | **Partial** | 相对中文稿多为 polish + 已知科学纠错，论证骨架未重构 | 对 §3.5–3.6、§3.8–3.9 做 **logic rewrite**，非仅润色 |
| Mostly append-only revision | **Likely** | 红字 Table/Figure 标签增补感强 | 将 SI 引用纳入完整论证，而非挂件 |
| New claims without evidence bank support | **Yes** | micromolar potency；enrichment 147×；ESM “physical constraints” for screening | 删除或改为假设；补 null model / ablation |
| Results still metric-dump style | **Yes（多段）** | §3.7 通过率列表；§3.9 指标清单式 | 每节用 promise→evidence→interpretation 重写首尾句 |
| Discussion repeats Results without resolving motivation | **Yes** | 收束段复述 MD/MM-GBSA，未回收双靶疾病动机 | 需真正 Discussion：回答“为何 PLK1+NLRP3”、局限、下一步 assay |
| LaTeX/formatting displaced logic | **Low** | 主要为 EN polish | — |
| Overclaim: computational score = biological superiority | **Yes** | vs MCC950 docking/MD | 改为 *predicted score under scoring function X*；禁止 “outperforming inhibitor” |
| Internal contradiction | **Yes** | Cys133 vs Cys97；50-PC vs \(p=270\) | Must-Fix 统一 |
| Decorative method component | **Suspected** | ESM-2 UMAP | Ablation：无 ESM 时漏斗是否变化；或降级为 exploratory figure |

### 4.4 Objection Register（供后续 rewrite 对账）

| ID | Objection | Location | Severity | Disposition needed |
|----|-----------|----------|----------|-------------------|
| O1 | 用语 “inhibitor / potency / lock-and-key” 超证据 | §3.8–3.9 | CRITICAL | 全局降调 + 局限段 |
| O2 | Enrichment factor 未定义 / 实为压缩比 | §3.6 | CRITICAL | 重算或删除 |
| O3 | ESM-2 未证明参与筛选决策 | §3.5 | CRITICAL | Ablation 或降级叙述 |
| O4 | Cys97 vs Cys133 | §3.3 vs §3.9 | CRITICAL（信任） | 核对 PDB 后统一 |
| O5 | Williams \(p=270\) vs 50 PCs | §3.4 | MAJOR | 重算 AD 或改正 \(p\) |
| O6 | Intersection null=25% 不合理 | §3.6 | MAJOR | 重做统计或改为描述性 |
| O7 | MCC950 docking 比较不公平 | §3.8 | MAJOR | 删除优劣措辞或加控制实验说明 |
| O8 | MM/GBSA→µM 活性 | §3.9 | MAJOR | 删除绝对换算 |
| O9 | MPO Σw=1.08；选第 3 名 | §3.8 | MAJOR | 归一化 + 预设规则 |
| O10 | UQ \(r=0.33\) 过度解读 | §3.4 | MINOR→MAJOR | 改措辞 + 校准图 |
| O11 | H-bond≪1 与“紧密锚定”张力 | §3.9 | MAJOR | 解释非键主导或检查分析 |
| O12 | 缺 Intro/Methods/Abstract | 全局 | BLOCKER | 补齐后方可送审 |

---

## 5. Consensus Must-Fix vs Nice-to-Fix Revision Roadmap

### Submission blockers（非科学细节，但是硬门槛）

1. 补齐 **Abstract, Introduction, Methods, References** 与完整 **SI**（表格/图在正文已引用）。
2. Methods 中固化：特征维度、Williams \(p\)、MM/GBSA 窗口、力场、对接网格、随机种子、软件版本。

### Must-Fix（共识；含仍相关的 prior flags）

| Priority | Item | Action |
|----------|------|--------|
| M1 | **Cys97 vs Cys133** | 对照 2RKU 编号与标准 hinge（Cys133）统一全文；若 Cys97 仅为局部接触，禁止称 “hinge” |
| M2 | **Claim 降调** | 全局替换/限定：*predicted dual-target binder / computational lead*；删除 micromolar 断语与 “outperforming MCC950” 抑制剂表述 |
| M3 | **Enrichment 147-fold** | 给出正式 enrichment 定义与 baseline，或改为 library reduction factor 并避免 enrichment 一词 |
| M4 | **Williams \(p\)** | 与真实模型描述符数一致后重算 \(h^*\) 与 AD 统计 |
| M5 | **ESM-2 角色** | 证明其进入相容性分数（公式/伪代码），或移出主线改为补充可视化 + 诚实局限 |
| M6 | **Intersection 统计学** | 修正 null model；勿在“受体结构迥异”前提下假设集合等同 |
| M7 | **MPO** | 权重归一化为 1.00；列出全部 18 参数；写明为何选 rank-3 Mol_997197 的预先规则 |
| M8 | **MM/GBSA 解释边界** | 仅作相对排序；禁止 IC50 推断；保留 50–100 ns 与脚本一致 |

*Prior flags 状态提醒：* ESM dims / Hindi / MM-GBSA window / pocket depth 在 polished EN 中已处理；**仅 Cys hinge 仍为 Must-Fix。**

### Nice-to-Fix

| Item | Action |
|------|--------|
| N1 | 修正 “and and PC1” |
| N2 | UQ：可靠性图、效应量谦逊表述；报告 B=30 敏感性 |
| N3 | H-bond 分析参数与疏水主导叙事对齐 |
| N4 | 残基编号对照表（7ALV / 2RKU） |
| N5 | Kinase selectivity / off-target 讨论展望 |
| N6 | 嵌套 CV 或置换检验强化 SVR vs XGBoost |
| N7 | 对 Gnina CNNscore 与 Vina 一致性做简短相关分析 |

### Suggested revision sequence

1. 科学一致性格（M1, M4, M7）→ 2. 主张与统计重写（M2, M3, M6, M8）→ 3. ESM 去留决策（M5）→ 4. 补齐全文与 SI → 5. Discussion 动机回收（非结果复述）→ 6. Nice-to-Fix 抛光。

---

## 6. Editorial Decision Letter（约 1 页）

**Manuscript:** PLK1–NLRP3 dual-target virtual screening / MD（Results & Discussion extract only）  
**Decision:** **Major Revision — Incomplete Manuscript**  
**Weighted score (R&D only):** **≈57.7 / 100**  

Dear Authors,

We thank you for submitting this polished Results & Discussion draft describing a dual-target computational pipeline against PLK1 and NLRP3, culminating in MD evaluation of `Mol_997197`. The panel recognizes a coherent screening narrative and several methodological strengths on the PLK1 arm (scaffold splitting, Y-scrambling, residual diagnostics, and an attempt at uncertainty-aware filtering).

However, **the package is not yet suitable for formal peer review or acceptance consideration**. Critical sections (Abstract, Introduction, Methods, References, and a complete Supplement) are absent. Independently, the available Results text contains **critical claim–evidence mismatches** and internal inconsistencies that would block acceptance even in a complete submission.

**Consensus required revisions**

1. **Unify PLK1 hinge annotation (Cys133 vs Cys97)** after PDB verification; eliminate contradictory structural claims.  
2. **Retone all biological efficacy language** to computational prediction; remove absolute potency inferences from MM/GBSA and “outperformance” of validated inhibitors based solely on docking/MD scores.  
3. **Correct or remove the ~147-fold enrichment claim** and the 25% intersection null model; replace with transparent, well-defined statistics.  
4. **Reconcile applicability-domain leverage parameter \(p\)** with the actual 50-PC SVR representation.  
5. **Clarify or demote ESM-2**: either demonstrate its quantitative role in NLRP3 compatibility scoring (with ablation) or present it as exploratory context only.  
6. **Normalize and fully specify the MPO function**, including the pre-specified rule for selecting a non–top-1 compound.  
7. **Supply the missing manuscript sections** with reproducible Methods aligned to the reported numbers (including MM/GBSA 50–100 ns).

Several prior scientific errors (ESM-2 dimensionality, trajectory window wording, NLRP3 pocket depth, and the Hindi artifact) appear corrected in the English extract; please maintain those corrections in the full paper and Methods logs.

We encourage resubmission of a **complete manuscript** after addressing the Must-Fix list. Upon return, the paper will be evaluated as a computational methods / virtual screening study whose contribution is a prioritization workflow—not a demonstrated dual-target inhibitor—unless experimental validation is added.

Sincerely,  
**Editorial Synthesizer**（ARS panel）  
Decision: **Major Revision — Incomplete Manuscript**

---

## Appendix A — Score snapshot

| Dimension | Score |
|-----------|------:|
| Originality | 68 |
| Methodological Rigor | 52 |
| Evidence Sufficiency | 48 |
| Argument Coherence | 55 |
| Writing Quality | 72 |
| **Weighted** | **57.65** |
| **Decision** | **Major Revision (Incomplete)** |

## Appendix B — Materials reviewed

- `/tmp/manuscript_polished_en.md`
- `/workspace/JNK1_Selectivity_Project/1_Results_and_Discussion_translated_checked.md`（Part 1 flags + section headings）
- ARS: `WORKFLOW.md`, `quality_rubrics.md`, agent role cues（EIC / Methodology / Domain / Devil's Advocate）
- PaperSpine: `logic-transfer-audit.md`（adapted to diagnosis + objection register）
