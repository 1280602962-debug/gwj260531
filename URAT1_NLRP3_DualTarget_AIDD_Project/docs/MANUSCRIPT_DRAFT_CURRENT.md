# Manuscript draft (JMM · 9DKB + 7ALV workflow)

**Title:** Clinical drug repurposing for gout-related URAT1 and NLRP3 targets: NLRP3 machine-learning prescreening, dual-target docking at 9DKB and 7ALV, and molecular dynamics of benchmark inhibitors

**Running title:** URAT1–NLRP3 dual-node repurposing funnel

**Target journal:** *Journal of Molecular Modeling*

---

## Abstract

Gout involves hyperuricemia (urate transporter URAT1, SLC22A12) and inflammasome-driven inflammation (NLRP3). We asked whether an asymmetric computational funnel—fast NLRP3 machine learning (ML) on a clinical drug library, followed by structure-based filtering at inward-open URAT1 (PDB 9DKB) and NLRP3 NACHT (PDB 7ALV)—could prioritize repurposing candidates without claiming de novo dual-target discovery. From 8,319 ChEMBL clinical compounds, NLRP3 ensemble classification (scaffold cross-validation AUROC 0.89) yielded 1,588 compounds with P(active) ≥ 0.5. Glide XP docking produced dual scores for 1,451 compounds; Pareto analysis on URAT1 percentile (S_U) versus max(NLRP3 ML percentile, NLRP3 docking percentile) identified six non-dominated leads. Known URAT1 inhibitors lesinurad and verinurad ranked high on S_U but not on the Pareto front; colchicine showed high ML probability yet poor dual-axis balance, consistent with indirect NLRP3 modulation. URAT1 evidence was independently validated on a 8,973-compound distill set (ROC-AUC 0.705; enrichment factor at 5% = 4.23), where dotinurad recovered ~89th docking percentile despite ML failure. Spearman correlation between NLRP3 ML scores and 7ALV XP scores within the docking pool was negligible (ρ ≈ −0.04), supporting a dual-evidence scoring axis. We recommend four 50–100 ns molecular dynamics simulations (benzbromarone and dotinurad @ 9DKB; MCC950 and one Pareto lead @ 7ALV) to complete mechanistic interpretation. The workflow is computationally tractable and benchmark-consistent; Pareto hits are hypothesis-generating early-phase agents requiring literature and experimental follow-up.

**Keywords:** URAT1; NLRP3; drug repurposing; molecular docking; machine learning; gout; Pareto optimization

---

## 1. Introduction

Gout is the most common inflammatory arthritis in adults, driven by monosodium urate crystal deposition and NLRP3 inflammasome activation with interleukin-1β release. Clinically, urate-lowering therapy (xanthine oxidase inhibitors, uricosurics targeting URAT1) and anti-inflammatory agents (e.g., colchicine) address distinct nodes of the disease. A single approved agent that simultaneously lowers urate and suppresses NLRP3-driven inflammation remains elusive; computational repurposing of existing clinical molecules offers a lower-risk route to generate testable hypotheses.

Structure-guided campaigns on URAT1 have accelerated following inward-open cryo-electron microscopy structures (e.g., PDB 9DKB). NLRP3 drug discovery, by contrast, relies on heterogeneous cell-based assays and multiple conformational states; direct NACHT binders such as MCC950 differ mechanistically from microtubule-directed inflammasome modulators. These asymmetries imply that a unified ML classifier cannot be expected to rank URAT1 uricosurics and NLRP3 tool compounds with equal reliability.

We therefore designed a **dual-evidence funnel**: (i) NLRP3 ML on the full clinical library to reduce chemical space; (ii) parallel Glide XP at URAT1 9DKB and NLRP3 7ALV on the ML-enriched pool; (iii) Pareto ranking on URAT1 docking percentile versus a composite NLRP3 score; (iv) independent URAT1 retrospective enrichment on a large distill docking set (PDB 8973-related panel, n = 8,973); and (v) targeted molecular dynamics (MD) on benchmark inhibitors plus one repurposing lead. This manuscript reports the completed ML, retrospective URAT1, and dual-docking stages; MD figures are specified for local completion.

---

## 2. Methods

### 2.1 Clinical repurposing library

ChEMBL-derived clinical compounds (max phase ≥ 1 and/or ATC-assigned) were merged into a manifest of 8,319 unique small molecules (`repurposing_manifest.csv`). Compounds were standardized to canonical SMILES. Literature benchmarks (lesinurad, benzbromarone, verinurad, dotinurad, colchicine, allopurinol, MCC950) were annotated but not used to train repurposing scores.

### 2.2 NLRP3 machine learning prescreen

An ensemble classifier (Chemprop-style graph models on five NLRP3-related bioassays; scaffold-split five-fold cross-validation) predicted P(active) for each library compound. The operating threshold P(active) ≥ 0.5 defined the **docking pool** (n = 1,588). Model quality on held-out folds: AUROC 0.89, supporting use as a **library reduction** step rather than a sole ranking metric.

### 2.3 Molecular docking

**URAT1:** Glide XP against PDB **9DKB** (inward-open URAT1).  
**NLRP3:** Glide XP against PDB **7ALV** (NACHT domain with MCC950-class sulfonylurea analog in the binding site).

Ligands from the docking pool were prepared in Schrödinger Maestro (LigPrep, OPLS4). XP docking used default precision; the best pose per compound (lowest XP score) was retained. Canvas batch exports were normalized with `scripts/normalize_canvas_docking_export.py` (join on `repurposing_id`).

### 2.4 Dual-target merge and Pareto analysis

For compounds with both URAT1 and NLRP3 scores (n = 1,451), we computed:

- **S_U:** percentile rank of 9DKB XP score (more negative = better).  
- **S_N:** max( NLRP3 ML percentile, 7ALV XP percentile ) when `--sn-mode both`.

Non-dominated points maximizing both axes defined the **Pareto front** (`scripts/merge_docking_pareto.py`). Six compounds formed the shortlist.

### 2.5 URAT1 retrospective track (8973 distill)

Separately, 8,973 compounds with URAT1-related docking on the 8973 structural ensemble were used **only** for URAT1 enrichment (active vs decoy labels from curated pActivity). Metrics: ROC-AUC, enrichment factor at 5% (EF@5%), and benchmark percentiles for four uricosurics. This track did not enter NLRP3 ML or Pareto merging.

### 2.6 Planned molecular dynamics

Four systems (50–100 ns, Desmond or GROMACS):  
**URAT1 @ 9DKB:** benzbromarone, dotinurad (poses from standalone redocking; not in P≥0.5 pool).  
**NLRP3 @ 7ALV:** MCC950 (redock to analog template), epigallocatechin gallate (EGCG) or fosigotifator (Pareto representative).  
Analysis: backbone RMSD, key residue distances (URAT1 Phe cage; NLRP3 Walker B / sulfonylurea pocket), qualitative MM-GBSA.

---

## 3. Results

### 3.1 Asymmetric data justify the funnel design

URAT1 regression ML on curated ChEMBL data recovered only two of four benchmark uricosurics under strict activity thresholds, whereas NLRP3 classification achieved AUROC 0.89 on scaffold splits. Zero compound overlap exists between the 8973 distill set and the 8,319 clinical manifest, mandating parallel URAT1 retrospective validation rather than a single merged training table.

### 3.2 NLRP3 ML prescreen (Fig. 2)

Of 8,319 compounds, 1,588 (19.1%) exceeded P(active) ≥ 0.5. Gout co-medications allopurinol, febuxostat, benzbromarone, and dotinurad showed low NLRP3 probabilities (P ≈ 0), as expected for non-inflammasome mechanisms. Conversely, colchicine and verinurad scored high (P ≈ 0.92), illustrating phenotypic and off-target confounding. Phase composition skewed toward Phase I/II among high-scoring compounds (Discussion). After dual docking, **1,451** compounds retained both XP scores (137 pool members lacked 9DKB and/or 7ALV poses).

### 3.3 URAT1 8973 retrospective (Fig. 3)

On subset A (n = 822 actives) versus subset D decoys (n = 7,957 docked), 9DKB XP achieved ROC-AUC **0.705** and **EF@5% = 4.23**. Dotinurad ranked ~89th percentile by docking but ~5th by ML pActivity, demonstrating that **URAT1 ranking in this project must be docking-led**. Lesinurad and verinurad showed moderate-to-high docking recovery; benzbromarone was well recovered on 8973 but excluded from the NLRP3 funnel by design.

### 3.4 Dual docking Pareto @ 9DKB + 7ALV (Fig. 4)

| Compound | S_U (%) | S_N (%) | max_phase |
|----------|---------|---------|-----------|
| SLV-334 | 99.9 | 92.1 | 2 |
| LANPROSTON | 99.9 | 96.8 | 2 |
| LASALOCID | 99.7 | 98.3 | 2 |
| Epigallocatechin gallate | 99.2 | 99.7 | 3 |
| Fosigotifator | 98.7 | 99.8 | 2 |
| Fosravuconazole | 96.9 | 99.9 | 2 |

**Benchmarks within the dual merge (not Pareto):**

| Drug | S_U | S_N | Interpretation |
|------|-----|-----|----------------|
| Lesinurad | 91.6 | 95.0 | Strong URAT1 docking among P≥0.5 set; NLRP3 axis inflated by ML |
| Verinurad | 77.7 | 97.9 | Similar pattern |
| Colchicine | 30.7 | 50.1 | High ML but poor URAT1 and moderate NLRP3 docking—consistent with indirect mechanism |

Spearman ρ(P(active), 7ALV XP) = **−0.036** (p = 0.17) across the merged set, indicating weak linear coupling between ML prescreen and structural NLRP3 scores.

### 3.5 Comparison to benchmarks—is the outcome reasonable?

**Yes, for a methods-and-repurposing-hypothesis paper**, with explicit caveats:

1. **URAT1 axis:** 8973 enrichment and in-pool lesinurad/verinurad percentiles align with known transporter pharmacology; potent uricosurics excluded from P≥0.5 are handled via the retrospective track and planned MD—not as Pareto failures.  
2. **NLRP3 axis:** Weak ML–docking correlation justifies S_N = max(ML, dock). Colchicine fails Pareto despite high ML, as required for a deconfounded funnel.  
3. **Repurposing output:** Six Pareto compounds are **early-phase** (mostly Phase 2; EGCG Phase 3). They are dual-high **computational** hits, not validated dual-target drugs.  
4. **Versus “better” benchmarks:** The funnel does not outperform dedicated URAT1 or NLRP3 single-target campaigns on approved drugs; it **combines** evidence under complementary constraints. That is the intended contribution.

### 3.6 Molecular dynamics (to be completed locally)

Table 2 (planned): RMSD stability, key contacts, and relative MM-GBSA for benzbromarone, dotinurad, MCC950, and EGCG (or fosigotifator). These simulations will anchor Discussion of binding modes without extending MD to all six Pareto molecules.

---

## 4. Discussion

We present a **gout-relevant dual-node repurposing workflow** that respects unequal data depth for URAT1 and NLRP3. NLRP3 ML efficiently shrinks the library; URAT1 filtering via 9DKB removes chemotypes with poor transporter complementarity; Pareto analysis surfaces compounds that are not dominated on either axis.

**7ALV versus 8ETR:** 7ALV provides an MCC950-class sulfonylurea pocket template appropriate for direct NACHT binders. Alternative NLRP3 structures (e.g., 8ETR) may be compared in supplementary material but are not required for the main narrative.

**Colchicine and phase bias:** High NLRP3 ML scores for colchicine and verinurad reflect assay and indication confounding, not reliable direct NACHT inhibition. Enrichment of Phase I/II agents among top ML scores further cautions against interpreting clinical phase as biological validation.

**Pareto shortlist:** EGCG offers the most mature clinical profile among the six leads but carries known pharmacokinetic limitations; fosigotifator represents a modern sulfonylurea-like scaffold with extreme docking ranks. All six require independent target engagement and safety review before any experimental claim.

**Limitations:** (i) 204 pool compounds missing dual scores; (ii) docking scores are not binding affinities; (iii) no experimental validation in this study; (iv) benchmark uricosurics outside P≥0.5 cannot appear on the main Pareto plot; (v) MD section pending.

We do **not** claim discovery of a first-in-class dual URAT1–NLRP3 inhibitor. We claim a **reproducible, benchmark-aware funnel** that yields a small hypothesis list and a clear MD follow-up.

---

## 5. Conclusions

An asymmetric NLRP3 ML prescreen plus 9DKB/7ALV dual docking reduced 8,319 clinical compounds to six Pareto-front repurposing candidates while preserving interpretable behavior for lesinurad, verinurad, and colchicine controls. URAT1 retrospective analysis on 8,973 compounds supports docking-led uricosuric ranking. The results are **sufficient and reasonable to proceed** to four benchmark-and-lead MD simulations and final figure assembly for *Journal of Molecular Modeling* submission.

---

## Supporting Information (planned)

- NLRP3 OOF ROC/PR curves  
- URAT1 OOF parity and data-asymmetry schematic  
- Phase ≥ 3 sensitivity subset (optional re-docking)  
- Full `pareto_merged_scores.csv` (1,451 rows)

---

## Data availability

Scripts and processed tables: GitHub repository `1280602962-debug/gwj260531`, directory `URAT1_NLRP3_DualTarget_AIDD_Project/`. Pareto outputs: `data/repurposing/pareto/`.

---

## Author contributions

*[To be completed]*

## Conflicts of interest

*[To be completed]*

---

*Numerical source of truth: `docs/RESULTS_DOCKING_9DKB_7ALV.md`, `results/repurposing/pareto_benchmark_report.json`.*
