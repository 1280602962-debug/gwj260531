# Language Optimization Report / 语言优化报告

**Manuscript unit reviewed:** Results & Discussion only (polished English extract)  
**Source:** `/tmp/manuscript_polished_en.md` (~3,500 words of R&D prose)  
**Target venue (assumed):** *Journal of Medicinal Chemistry* or similar CADD / computational medicinal chemistry journal  
**Playbooks applied:** PaperSpine `humanize` (general platform profile, medium tier advisory) + ARS `academic-paper` Writing Quality Check  
**Companion science flags:** Part 1 of `/workspace/JNK1_Selectivity_Project/1_Results_and_Discussion_translated_checked.md`  
**Date:** 2026-07-17  

---

## 1. Scope Note / 范围说明

| EN | 中文 |
|----|------|
| This review covers **only the Results & Discussion (R&D) excerpt**. Introduction, Methods, Abstract, Conclusions, and full IMRaD integration were **not** available. Judgments on journal readiness therefore apply to **R&D language quality**, not to the complete manuscript package (figures, SI consistency, Methods reproducibility, citation compliance). | 本审查**仅覆盖结果与讨论（R&D）摘录**。引言、方法、摘要、结论及完整 IMRaD 整合**不可用**。期刊就绪度判断因此仅针对 **R&D 语言质量**，不代表全文投稿包（图、SI 一致性、方法可复现性、引用合规）已就绪。 |
| No AIGC detector was run. **No fabricated detector % scores** are reported. Findings map to PaperSpine risk dimensions D1–D5 and ARS writing-quality categories A–E as qualitative/count-based signals only. | **未运行任何 AIGC 检测器**，**不虚构检测百分比**。发现仅作为 PaperSpine D1–D5 与 ARS A–E 类别的定性/计数信号。 |
| Part 1 scientific issues (ESM-2 dimensionality, Hindi artifact, MM/GBSA window, pocket depth, Cys97 vs Cys133) remain **science blockers** for submission even when English is polished; several are already corrected in the polished extract, but hinge-residue numbering is still internally inconsistent. | Part 1 科学问题仍是投稿的**科学阻断项**；抛光英译本已修正若干条，但铰链残基编号 internally 仍不一致。 |

---

## 2. Executive Verdict / 总体判定

**Verdict: Language-near-ready for R&D draft circulation; not yet submission-ready without a focused humanize + claim-tone pass.**

| Dimension | EN | 中文 |
|-----------|----|------|
| Overall | Prose is fluent, technically dense, and largely journal-register. Data reporting is concrete. Remaining issues are **template cadence**, **inflated claim diction**, a few **non-academic metaphors**, and **residual typos**—not broken English. | 行文流畅、技术密度高、大体符合期刊语域。数据陈述具体。主要问题是**模板化节奏**、**措辞夸大**、少量**非学术隐喻**与**残留笔误**，而非英语不通。 |
| Readiness | Suitable as an internal R&D draft after Top-10 fixes. For J. Med. Chem.–class submission, expect **1 dedicated language revision** (humanize matrix + ARS A/C/D/E) plus resolution of Part 1 science flags. | 完成 Top-10 修订后可作为内部 R&D 草稿。若投 J. Med. Chem. 级别，还需 **1 轮专项语言修订**，并消化 Part 1 科学问题。 |
| Detector disclaimer | We do **not** claim the text will pass any AIGC checker. Several high-visibility AI-risk patterns remain (purpose-opener stacks, mirrored MD subsections, “Crucially/This…” cadence). | **不承诺**可通过任何 AIGC 检测。若干高可见度风险模式仍在（目的状语开篇堆叠、MD 镜像小节、“Crucially/This…”节奏）。 |

**Blocking before submission (language + carry-over science):**
1. Fix typo `and and PC1`; replace `outclasses`, `champion model`, `"self-confidence awareness"`, overused `lock-and-key`.
2. Diversify section/paragraph openers (16× sentence-initial `To …`).
3. Resolve Part 1: Cys97 vs Cys133 hinge annotation; confirm MM/GBSA window and ESM-2 dim already match methods/scripts.

---

## 3. humanize_matrix / 去 AI 痕迹矩阵

Severity: **High** = claim tone / template stamp likely to draw reviewer or detector attention; **Med** = clear pattern, fix in revision; **Low** = polish.

| Row ID | Manuscript Unit | AI/Language Pattern | Severity | Suggested rewrite | Rationale |
|--------|-----------------|---------------------|----------|-------------------|-----------|
| H01 | §2 model selection | Inflated verb: “statistically confirming that SVR's predictive accuracy significantly **outclasses** that of the XGBoost model” | High | “…confirming that SVR’s test-set and CV MSE were lower than XGBoost under a paired *t*-test (*p* = 0.00611).” | “Outclasses” is sports/marketing register; JMC prefers quantitative understatement. |
| H02 | §2 / §4 | Anthropomorphic metaphor: model exhibits **“self-confidence awareness”** | High | “higher σᵢ when prediction error is large, consistent with useful uncertainty ranking” | Quotation marks do not salvage non-scientific anthropomorphism; invites reviewer pushback. |
| H03 | §2, §4 | Sports metaphor: **“champion SVR model”** / “champion model's hyperparameters” | Med | “final SVR model” / “selected SVR hyperparameters” | “Champion” is competition jargon; standard: selected / final / best-performing. |
| H04 | §8, §10 | Cliché metaphor ×2: **“double lock-and-key matching”** / “dual lock-and-key matching mode” | Med | “complementary binding poses in both pockets” / “stable dual-pocket occupancy” | Lock-and-key is overused and imprecise for flexible NLRP3 NACHT site. |
| H05 | §2–§9 openings | **Purpose-infinitive stack** (16× sentence-initial `To …`): e.g. “To select…”, “To overcome…”, “To rule out…”, “To evaluate…” | High | Vary: “We next compared…”, “Permutation importance then…”, “Y-scrambling…”, “100 ns MD of Mol_997197…” | Classic AI/template cadence (PaperSpine D2/D4); 5/10 sections open with `To`. |
| H06 | §6, §7, §9 | Discourse tic: **“Crucially,”** ×3 | Med | Delete adverb; lead with the fact (“The intersection fraction was 13.6%…”) | ARS flagged-family *crucial*; thrice in one R&D section is throat-clearing. |
| H07 | Multiple paras | Template closer: **“This [noun] indicates/demonstrates/confirms/supports…”** (~10 instances) | Med | Fold interpretation into the data sentence, or use specific subject (“Scaffold independence did not…”). | Reduces D2 paragraph-template similarity; keeps claim attached to evidence. |
| H08 | §1 | Flagged term (non-exempt): “To establish a **robust** QSAR model” | Low | “To fit a QSAR regression model” / “a predictive QSAR model” | Not statistical “robust estimator”; ARS A flags vague quality claim. |
| H09 | §10 summary | Flagged term: “achieves a **robust** dynamic balance” | Med | “maintains comparable MM/GBSA binding free energies for PLK1 and NLRP3 (−32.68 vs −32.15 kcal/mol)” | Replace vague “robust balance” with numbers already in text. |
| H10 | §5 | Flagged verb: “This pre-filter **streamlined** the library to 10,000…” | Low | “reduced the library to 10,000…” | ARS A: *streamline* = non-academic business verb. |
| H11 | §9 vs §10 | **Mirror structure**: identical 1–5 headers (RMSD / *R*g–SASA / H-bonds / MM/GBSA / per-residue) | Med | Keep parallel metrics, but vary prose order/emphasis (e.g., lead NLRP3 with ligand RMSD advantage; shorten PLK1 H-bond item). | ARS D Mirror Structure + PaperSpine D2; dual-target comparison can stay parallel without cloning outline. |
| H12 | §7 | **Rule-of-three** ADMET list (Lipinski / SA / solubility–HIA) then “Crucially, hERG…” | Low–Med | Keep items if all three are needed; break list into prose or add hERG as equal peer, not dramatic coda. | Trio + “Crucially” coda reads template-stamped. |
| H13 | §3 | Typo: “shuffling PC3, PC2, **and and** PC1” | High | “shuffling PC3, PC2, and PC1” | Hard typo; must fix before any submission. |
| H14 | §2–§10 | Verb monoculture: *exhibit\** ×10, *indicat\** ×10, *confirm\** ×7 | Med | Rotate: *showed*, *was*, *gave*, *was consistent with*, *supported*. | D3 lexical diversity; academic clarity prefers verb variation without synonym cycling of technical nouns. |
| H15 | §3, §4 | Connector density: **Concurrently** ×2, **Additionally**, **Consequently** ×2, **Specifically**, **Notably** | Med | Cut half; join clauses with content logic (“…Cys133. The sulfonamide arene…”) | PaperSpine D4 connector overuse. |
| H16 | §4 UQ | Quotation + slogan tone around uncertainty | Med | State Pearson *r* and filtration ΔMSE only; drop “self-confidence” framing. | Keeps UQ scientific; removes AI-flavored flourish. |
| H17 | §9.1 | Mid-paragraph **“Crucially,”** before ligand RMSD | Low | “Ligand RMSD of Mol_997197 converged to…” | Emphasis adverb unnecessary when numbers follow. |
| H18 | §8 docking | Emphatic bold + comparative overclaim risk: docking score “**outperforming** MCC950” | Med | “more favorable docking score than MCC950 (−8.06 vs −5.70 kcal/mol); docking scores are not binding free energies” | Softens claim; docking ≠ potency (science+tone). |
| H19 | §1–§4 | Uniform “methods→numbers→This interprets” paragraph recipe | Med | Alternate: lead with result number; end some paras without meta-sentence. | D2 adjacent-paragraph similarity. |
| H20 | Part 1 carry-over | Hinge residue **Cys133** (§3) vs **CYS97** (§9.1 note) | High (science) | Align numbering to PDB `2RKU` / UniProt; one residue ID throughout. | Not “AI style,” but Part 1 flag still visible in polished EN. |

**Exempt (do not “humanize” away):** *leverage* in Williams-plot / hat-matrix sense (`h`, `h*`) — standard QSAR statistics terminology (ARS Exception Rule).

**PaperSpine D1–D5 (qualitative only; no detector %):**

| Dim | Observation (R&D extract) | Status (advisory) |
|-----|---------------------------|-------------------|
| D1 Burstiness | ~130 sentences; mean ~30 words; CV ≈ 0.45; several 5-sentence mid-length runs; few ≤12-word punches in Discussion-like stretches | Mixed — improve in MD summary & UQ claims |
| D2 Paragraph template | Many `To…` opens; mirrored §9/§10; recurring `This…` closers | Needs revision |
| D3 Lexical diversity | Technical vocabulary good; verb set narrow (*exhibit/indicate/confirm*) | Needs light pass |
| D4 Connectors | Crucially/Concurrently/Consequently/Additionally stack | Needs revision |
| D5 Generic context | Mostly anchored by numbers; weak spots are metaphors (champion, lock-and-key, self-confidence) | Targeted fixes |

---

## 4. ARS Writing Quality Check / 写作质量检查（分类）

### A. High-Frequency Term Warnings / 高频措辞预警

| Term | Count in R&D extract | Notes / 说明 |
|------|----------------------|--------------|
| crucial / Crucially | **3** | All sentence-level emphasis; replace with bare claims. |
| robust | **2** | QSAR setup + MD summary; neither is “robust estimator.” |
| streamline / streamlined | **1** | “streamlined the library…” → *reduced* / *narrowed*. |
| leverage | **3** | **Exempt** — statistical leverage in Williams plot. |
| delve, tapestry, landscape, pivotal, foster, showcase, testament, navigate, realm, embark, underscore, multifaceted, nuanced, comprehensive, intricate, cornerstone, paradigm, synergy, holistic, cutting-edge, groundbreaking | **0** | Clean. |

**Related inflation (not on ARS table but same spirit):** `outclasses` (1), `champion` (2), `self-confidence awareness` (1), `lock-and-key` (2).

### B. Punctuation Pattern Control / 标点模式

| Rule | Finding | Assessment |
|------|----------|------------|
| Em dash (—) ≤3 / paper | **0** em dashes in extract | Clean |
| En dash (–) | **5** (mostly ranges / Figure S8–S10) | Acceptable for scientific ranges |
| Semicolons ≤2 / 1000 words | **22** semicolons ≈ **6.3 / 1000 words** | **Over limit** if applied to extract alone; majority are `Figure X; Table Y` citation pairs — prefer commas or “and” in running text (`Figure 3 and Table 1`). True clause-chaining semis are fewer (e.g., docking score sentence with “meanwhile”). |
| Colon-list sequences | One intentional numbered ADMET block; MD uses numbered metric headers | Acceptable; avoid adding more colon-led enumerations in adjacent paragraphs |

### C. Throat-Clearing Openers / 清嗓式开场

| Checklist phrase | Count |
|------------------|------:|
| In the realm of… / It’s important to note… / It is worth mentioning… / In today’s rapidly evolving… / testament… / goes without saying… / In order to… / It should be noted… / As a matter of fact… / When it comes to… / At the end of the day… / With that being said… / This section will discuss… / We now turn… | **0** |

**Functional equivalent still present:** purpose stacks (`To overcome…`, `To rule out…`, `To evaluate…`) and emphasis adverbs (`Crucially`, `Notably`). Treat as soft throat-clearing under C + D4.

### D. Structure Pattern Warnings / 结构模式

| Pattern | Finding |
|---------|---------|
| Rule of Three | ADMET pass-rate trio (Lipinski / SA / solubility–HIA) then dramatic hERG turn; MD metrics forced into five mirrored bullets for both targets |
| Uniform paragraph length | Long workflow paragraphs dominate early sections; MD uses clipped numbered items — variation exists but early R&D is formulaic |
| Synonym cycling | Mild (`candidates` / `molecules` / `compounds`); technical nouns mostly stable — **good**. Prefer keeping *compound* or *candidate* per subsection |
| Binary contrast | `rather than` ×2, `instead of` ×1 — within ≤2 “Not X. Y.” spirit; OK |
| Mirror structure | **Strong** between PLK1 MD (§9.1) and NLRP3 MD (§10): identical five headings |

### E. Burstiness (Sentence Rhythm) / 句长起伏

| Metric | Approx. value |
|--------|----------------|
| Sentence count | ~130 |
| Mean length | ~30 words |
| SD / CV | ~13.7 / ~0.45 |
| Short sentences (≤12 words) | ~6 (sparse for Discussion-level emphasis) |
| Flagged mid-length runs (≥5 consecutive ~18–34 words) | ≥4 stretches (scaffold-split discussion; SVR stats wrap-up; ADMET list vicinity; MD setup) |

**Guidance:** Insert 1–2 short factual sentences per major subsection (e.g., after enrichment factor; after hERG attrition; after each MM/GBSA headline). Discussion-level R&D should show higher burstiness than Methods.

### Additional language notes / 其他语言要点

- **Register:** Mostly appropriate for computational chemistry R&D.
- **Hedging:** Occasional over-certainty (“statistically confirming… outclasses”; “significantly outperforming MCC950” for docking/MM/GBSA without experimental IC₅₀). Soften computational claims.
- **Typo:** `and and PC1` must be fixed.
- **Part 1 science (still relevant to polished EN):** Cys133 vs Cys97 inconsistency remains; ESM-2 **640**-dim and **50–100 ns** MM/GBSA and **buried** NLRP3 pocket appear corrected vs original Chinese issues.

---

## 5. Priority Fix List (Top 10) / 优先修订清单（Top 10）

Exact original snippet → suggested rewrite.

### 1. Typo — PC importance sentence

**Original:**
> shuffling PC3, PC2, and and PC1 caused the test $R^2$ to drop by 0.22, 0.18, and 0.15, respectively.

**Suggested:**
> shuffling PC3, PC2, and PC1 lowered test $R^2$ by 0.22, 0.18, and 0.15, respectively.

---

### 2. “Outclasses” claim tone

**Original:**
> statistically confirming that SVR's predictive accuracy significantly outclasses that of the XGBoost model.

**Suggested:**
> supporting a lower CV MSE for SVR than for XGBoost under a paired $t$-test ($p = 0.00611$).

---

### 3. Anthropomorphic UQ slogan

**Original:**
> This correlation shows that the model exhibits "self-confidence awareness," automatically generating higher uncertainty estimates when evaluating structurally novel chemical spaces.

**Suggested:**
> Thus $\sigma_i$ tracks absolute error (Pearson $r = 0.3266$, $p = 7.29 \times 10^{-5}$), so high-uncertainty compounds can be deprioritized in screening.

---

### 4. “Champion model”

**Original:**
> The finalized champion SVR model achieved the following generalization metrics on the independent test set: $R^2 = 0.74$, $MSE = 0.42$, and $MAE = 0.44$

**Suggested:**
> The selected SVR model reached $R^2 = 0.74$, MSE $= 0.42$, and MAE $= 0.44$ on the independent test set

---

### 5. Purpose-opener diversification (section lead)

**Original:**
> To overcome the limitations of the "black-box" nature of machine learning models in rational drug design, we performed permutation importance analysis (10 permutations) on the 50 PCA features of the SVR model (**Figure 6**).

**Suggested:**
> Permutation importance (10 shuffles) on the 50 PCA inputs of the SVR model ranked features that drive PLK1 $\mathrm{pIC}_{50}$ prediction (**Figure 6**).

---

### 6. “Crucially” + dual-target intersection

**Original:**
> Crucially, the proportion of dual-target candidates appearing in the intersection of both pools was 13.6%.

**Suggested:**
> The dual-target intersection was 13.6% of each 500-compound pool.

---

### 7. “Streamlined”

**Original:**
> This pre-filter streamlined the library to 10,000 high-confidence candidates for downstream docking, reducing computational cost.

**Suggested:**
> The pre-filter reduced the library to 10,000 candidates for docking and cut the compute load accordingly.

---

### 8. Lock-and-key metaphor (docking)

**Original:**
> In its docked conformation, `Mol_997197` exhibited a balanced **double lock-and-key matching** profile

**Suggested:**
> In the docked poses, `Mol_997197` occupied both pockets with comparable docking scores

---

### 9. Lock-and-key + “robust” summary closer

**Original:**
> In summary, 100 ns molecular dynamics simulations and MM/GBSA calculations demonstrate that `Mol_997197` achieves a robust dynamic balance between PLK1 and NLRP3. This compound binds both targets in a dual lock-and-key matching mode, making it a promising lead candidate for chemical synthesis and biological evaluation.

**Suggested:**
> Over 100 ns, MM/GBSA averages were $-32.68 \pm 3.06$ kcal/mol (PLK1) and $-32.15 \pm 2.10$ kcal/mol (NLRP3), with ligand RMSDs below 1 Å in both pockets. These results support advancing `Mol_997197` to synthesis and assay.

---

### 10. “Crucially” hERG + “robust” QSAR opener (pair fix)

**Original (a):**
> To establish a robust quantitative structure-activity relationship (QSAR) model, we retrieved and curated PLK1 inhibitory activity data from the ChEMBL database.

**Suggested (a):**
> PLK1 $\mathrm{pIC}_{50}$ records were retrieved from ChEMBL and curated for QSAR modeling.

**Original (b):**
> Crucially, hERG K+ channel inhibition was identified as the primary barrier to safety, eliminating 42% of the candidate molecules

**Suggested (b):**
> Predicted hERG liability removed 42% of candidates (probability $> 0.5$) and was the largest single ADMET attrition filter

---

## 6. Part 1 Science Flags Still Touching Language / 仍影响表述的科学问题

| Part 1 issue | Status in polished EN extract | Language action |
|--------------|-------------------------------|-----------------|
| ESM-2 480 vs 640 | States **640-dimensional** for `esm2_t30_150M_UR50D` | Keep; verify against code |
| Hindi `मैच` artifact | Not present in EN | OK |
| MM/GBSA 46–50 ns vs 50–100 ns | States **50–100 ns** | Keep; verify scripts |
| “Shallow” NLRP3 pocket | “buried, highly flexible NACHT binding pocket” | OK |
| Cys97 vs Cys133 | §3 uses **Cys133**; §9.1 still discusses **CYS97** as hinge contact in `2RKU` | **Must reconcile** before submission |

---

## 7. Recommended Revision Order / 建议修订顺序

1. Mechanical: `and and`; Cys97/Cys133; figure/table semicolon style.  
2. Claim tone: outclasses / champion / self-confidence / lock-and-key / robust balance.  
3. Cadence: cut ~half of `To…` openers and `Crucially/Concurrently/Additionally`.  
4. Rhythm: add short emphasis sentences in UQ, docking enrichment, ADMET, MD energy headlines.  
5. Re-read §9–§10 for mirrored outline without losing dual-target comparability.  

**Explicit non-claims:** No AIGC detector percentage is reported or implied. Passing any commercial detector is not guaranteed by these edits.

---

*End of report / 报告结束*
