# PROJECT SCIENTIFIC MAP

Date: 2026-09-19  
Tree: PR #39 `cursor/scientific-freeze-c7cc` / Dual_Target_Docking  
Rule: this map describes the scientific chain. Folder names are not the chain.

## Chain

**Research question**  
Can rigid-receptor docking distinguish dual-active ligands from single-target-active ligands in dual-target virtual screening, and how does the negative-class definition change that judgment?

This is a retrospective four-state discrimination question, not a prospective hit-finding campaign and not a claim that docking “works” in general.

↓

**Data sources**  
ChEMBL activity (three kinds, not one): `chembl_assay_adjudicated` (EGFR/HER2, PIK3CA/mTOR, most AChE/BChE); `chembl37_dump_panel` (Track-B five pairs); `panel_pchembl_no_audit_rows` (five AChE ligands; only AB_056 in primary).  
Structures: 14 deposited receptor PDBQT slots.  
Ligands: panel SMILES; exact production PDBQT confirmed only for current EGFR/HER2 uniform prep.

↓

**Target-pair screening**  
Universe census → strict 6.5/5.5 supply / human holo / site interpretability / noncovalent Vina feasibility.  
Runtime roster is the hard-coded `PRIMARY_PAIRS` tuple.  
EGFR/HER2 is a documented supply-limited exception.  
PIK3CA/PIK3CB fails human-holo gate (human PIK3CB has zero PDB xrefs).  
Eligibility CSV contains no AUROC / GNINA / ECFP tokens.

↓

**Panel construction**  
Six pairs built under strict 6.5/5.5 quotas; EGFR/HER2 and PIK3CA/mTOR under θ=6.0 supply.  
Construction quotas (Table 1) are not Table 2 n_scored.  
EH120_059 and AB_087 are `unresolved_missing_arm` and are excluded from primary.

↓

**Activity labels**  
Primary estimands use θ=6.0 four-state only: dual / A-only / B-only / neither.  
Strict 6.5/5.5 is construction / supply audit, not a second Table 2 definition.  
Independent replay: 0 primary θ=6.0 class mismatches versus `pA`/`pB`.

↓

**Ligand input**  
EGFR/HER2 current: RDKit ETKDGv3 + MMFF + Meeko; 110/110 PDBQT deposited.  
PIK3CA / AChE / Track-B: recovered scripts exist; exact production PDBQT **not in tree**.  
Cannot claim eight-pair uniform ligand preparation as a current fact.

↓

**Receptor input**  
14 deposited PDBQT files are the production receptors.  
Protonation, pH, loop rebuild, and a common Meeko/PPW campaign are **NOT_RECOVERABLE**.  
Shared files: JAK1 6N7A; PPARA 6LXA.

↓

**Docking box**  
Current rule: cognate heavy-atom AABB + 5 Å, min edge 20 Å.  
EGFR uses `*_box_corrected.json`. Legacy `3POZ_box.json` / `3RCD_box.json` are superseded historical snapshots in `boxes/archive/`.  
PIK3CA boxes lack a construction field (`METADATA_ONLY`).

↓

**Docking**  
PRIMARY: AutoDock Vina 1.2.7; seed 20260727; E=16 only for PIK3CA/mTOR; E=8 otherwise.  
SENSITIVITY: PIK3CA E=8, PM110, 4JPS, 5DXT, 4JSX, five-seed, label/aggregation.  
INDEPENDENT GNINA: EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2 only.  
SAME-POSE: RTM / GNINA CNN of Vina poses — not independent docking.

↓

**Score extraction**  
`score_S = −energy` (higher-better).  
40/40 random source→master traces MATCH (5 ligands × 8 pairs).  
Timeouts skipped at 600 s (EGFR EH40_31 all seeds); not filled with 0 or extremes.

↓

**Primary estimand**  
D vs A-only → pocket B (`score_B`)  
D vs B-only → pocket A (`score_A`)  
`summary_min` = weaker of the two directions  
Two-pocket mean is used only for D vs neither and ranking (Table 3), not directional Table 2.

↓

**Sensitivity analyses**  
Five-seed (available-case + fixed-membership), max vs median, label thresholds, cluster/non-stratified bootstrap, receptor swap, unused-pool holdout, protocol E=8/PM110.

↓

**Ligand-only baseline**  
ECFP4 GroupKFold (radius 2, 2048 bits) vs docking-only vs ECFP4+docking.  
Descriptors: full-panel best = descriptive univariate screen; nested scaffold-CV = predictive chemistry control.

↓

**Statistics**  
Class-stratified shared-dual percentile bootstrap, B=2000, seed 20260729.  
Point estimate from the original sample, not the bootstrap mean.  
Independent audit-only replay: 192/192 MATCH (tolerance 1e-4).

↓

**Canonical tables**  
`results/canonical/*.csv` written by `compute_canonical_results.py` from `current_score_master.csv`.  
`data/processed/current_score_master.csv` was a superseded identical duplicate; current pointer is `data/processed/CURRENT_SCORE_MASTER.md`.

↓

**Figures**  
`figures/jcim_article/plotted_values_postfix.json` is the plotted-value lock (339 records).  
Figure 5C still stores a fuzzy `five_seed_comparable=1` for AChE while `same_membership_as_primary=0`.

↓

**Submission pack**  
`submission_pack/` mirrors manuscript, SI, figures, and canonical tables.  
It currently copies the SI Table S6 prose error.

↓

**Manuscript claims**  
Allowed: negative-class dependence; pair-specific directional AUROCs; chemistry already separates classes; docking adds little incremental ECFP information; matched-pocket advantage is pair-specific (AChE and EGFR exclude 0 on main panels); no external docking validation (0/8 eligibility).  
Not allowed: eight independent replicates; uniform receptor/ligand prep; “only AChE excludes 0”; observed power; external validation performed.

## Experiment classification

| Experiment | Class |
|------------|--------|
| Directional AUROC + summary_min (Table 2) | ESSENTIAL |
| Fixed-score Δ (neither vs selective) | ESSENTIAL |
| ECFP4 incremental | ESSENTIAL |
| Matched/mismatched pocket | ESSENTIAL |
| Two-pocket mean ranking / Top10% / AND (Table 3) | SUPPORTING |
| Nested descriptor CV | SUPPORTING |
| Full-panel best descriptor | DIAGNOSTIC (descriptive) |
| Five-seed / fixed-membership | DIAGNOSTIC |
| Receptor substitution / E=8 / PM110 | DIAGNOSTIC |
| Independent GNINA (3 pairs) | DIAGNOSTIC |
| Same-pose RTM/CNN | DIAGNOSTIC |
| Cognate RMSD | DIAGNOSTIC (search coverage) |
| External eligibility 0/8 | SUPPORTING (eligibility only) |
| Detectable-effect simulation | SUPPORTING (not observed power) |
| Historical EGFR ablation / LigPrep / PIK3CA/PIK3CB figures | HISTORICAL |
| Old 0.3237 / 0.0112 / 9.505 Å as current | UNNECESSARY / SUPERSEDED |
