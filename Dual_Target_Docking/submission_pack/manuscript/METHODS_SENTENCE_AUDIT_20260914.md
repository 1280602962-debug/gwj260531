# Methods sentence audit (pasted 2.1–2.7)

Date: 2026-09-14  
Branch checked: `cursor/jcim-language-polish-0b1a` at `dbe7ac1d`  
Text checked: the English Methods 2.1–2.7 supplied in this request, not the manuscript draft in the branch.  
Evidence: locked CSVs, docking/prep scripts, `ENV_PIN.md`, `STATISTICAL_LOCK_V1.md`, `external_slice_contract.yaml`, and the current typeset SI headings (Tables S1–S9).

## Revision (same day): no REST / dump split in the article

`parity_v1.csv` shows 0 missing ends and 0 mismatches for all eight scored panels against local ChEMBL 37 (tolerance 0.015). Table 2 max pChEMBL therefore matches ChEMBL 37. The article now treats ChEMBL 37 as the sole activity source. The eight-pair max-versus-median relabeling on those same records is the only typeset aggregation sensitivity (Table S3). The 2026-08-26 API snapshot and the three-pair high-confidence field screen are archived and are not typeset.

Do **not** restore “2026-07-23 ChEMBL REST harvest” versus “ChEMBL 37 dump” wording in Methods, Results, Discussion, or SI. Historical harvest paths remain in repository CSVs only.

Superseded items below: §1 (split-source harvest) and §3 (three-pair high-confidence in Methods). The implemented 2.2 text is: activities from ChEMBL 37; max pChEMBL; IC50/Ki/Kd/EC50/Potency/IC50app/Ki app; max-versus-median on the same records for all eight pairs.

Only sentences that need a factual or protocol change are listed. Unchanged sentences are accurate against the implementation.

SI numbering used below is the current typeset set:

| Table | Role |
|---|---|
| S1 | software, seeds, Vina knobs, bootstrap, ligand-prep recipe |
| S2 | receptors, boxes, CalcRMS cognate QC |
| S3 | θ / max–median label sensitivity on ChEMBL 37 |
| S4 | fixed-score dual-versus-neither (pocket held fixed) |
| S5 | ligand-chemistry baselines |
| S6 | matched vs mismatched pocket + unused-pool holdout |
| S7 | alt receptor, independent GNINA, PPARG rescore |
| S8 | BindingDB eligibility after independence filters |
| S9 | target-pair eligibility audit |

---

## 2.2 Bioactivity data

### 1. Table S3 is not the retrieval table (superseded: article now uses ChEMBL 37 only)

**Original**

> Experimental activity data were obtained from ChEMBL. Database versions, retrieval information, and activity-processing checks are provided in Table S3.

**Problem**

Table S3 is label/aggregation sensitivity, not the harvest record. Production labels are split-source. Table 2 does not come from a single ChEMBL download.

**Change to**

> Experimental activity data were obtained from ChEMBL. For EGFR/HER2, AChE/BChE, and PIK3CA/mTOR, Table 2 labels used a 2026-07-23 ChEMBL REST harvest. For F2/F10, JAK1/TYK2, JAK1/JAK2, PPARG/PPARA, and PPARA/PPARD, labels were taken from the ChEMBL 37 SQLite dump (release 2026-05-01). The pair-universe census used the same dump. A later dump join found that all eight production maximum pChEMBL values agree with ChEMBL 37 (tolerance 0.015); that check does not replace Table 2. Activity-label and aggregation sensitivity analyses are in Table S3.

**Evidence:** `parity_v1.csv` (`production_label_source`); `analyze_eight_pair_dump_gated_v1.py`; current Table S3 heading.

### 2. Activity-type list is incomplete

**Original**

> Quantitative records with IC50, Ki, Kd, EC50, or Potency measurements were retained when a pChEMBL value was available.

**Problem**

The locked harvest also keeps `IC50app` and `Ki app`.

**Change to**

> Quantitative records with IC50, Ki, Kd, EC50, Potency, IC50app, or Ki app measurements were retained when a pChEMBL value was available.

**Evidence:** `run_checks_v1.json` `standard_ok`; `chembl_exhaustive_pair_census_v1.py` `STANDARD_OK`.

### 3. High-confidence subset is not eight-pair (superseded: dropped from typeset Methods/SI)

**Original**

> Activity-processing sensitivity analyses compared maximum and median aggregation and examined a high-confidence human single-protein subset.

**Problem**

The high-confidence rebuild covers only EGFR/HER2, AChE/BChE, and PIK3CA/mTOR. It is an automatic field screen (human SINGLE PROTEIN, confidence ≥ 8), not a paper-by-paper reread, and it does not cover the five dump pairs.

**Change to**

> Activity-processing sensitivity analyses compared maximum and median aggregation. A high-confidence human SINGLE PROTEIN field screen (confidence score ≥ 8) was applied only to EGFR/HER2, AChE/BChE, and PIK3CA/mTOR (Table S3).

**Evidence:** `high_confidence_summary_v1.csv` (three rows); `high_confidence_run_meta_v1.json`.

---

## 2.3 Target-pair screening and panel construction

### 4. EGFR-only θ = 6.0 wording contradicts the next paragraph

**Original**

> Seven target pairs had sufficient single-target-active supply under the strict 6.5/5.5 rule, while EGFR/HER2 was evaluated from the θ=6.0 pool because it provided adequate directional sample sizes and compatible structural data.

**Problem**

EGFR/HER2 failed the 50/50 strict gate (dump-level B-only = 7). PIK3CA/mTOR met that gate (80/81) but its docked panel was still drawn from the θ = 6.0 pool. The next paragraph already says both EGFR/HER2 and PIK3CA/mTOR came from θ = 6.0 pools.

**Change to**

> EGFR/HER2 did not meet the 50/50 strict selective-supply gate and was retained as a supply-limited exception because it had compatible human holo structures and adequate dual, A-only, and B-only counts at θ = 6.0. PIK3CA/mTOR met dump-level strict supply, but its evaluation panel was still sampled from the θ = 6.0 pool. The other six pairs were sampled from the strict 6.5/5.5 pools.

**Evidence:** `parity_v1.csv` `dump_pool_strict_D_A_B` (EGFR 951/39/7; PIK3CA 1552/80/81); `build_egfr_her2_panel120.py` `THETA = 6.0`; `build_pik3ca_mtor_panel48.py` `THETA = 6.0`.

### 5. Scaffold caps for the other six pairs are unstated

**Original**

> Within each class, a Bemis–Murcko scaffold was represented at most five times for EGFR/HER2 and twice for PIK3CA/mTOR.

**Problem**

True for those two panels. The other six primary panels had no additional Murcko cap. Unused-pool holdouts used a Murcko cap of 3.

**Change to**

> Within each class, a Bemis–Murcko scaffold was represented at most five times for EGFR/HER2 and twice for PIK3CA/mTOR. The other six primary panels had no additional scaffold cap. Unused-pool holdout draws used a Murcko cap of 3.

**Evidence:** `build_egfr_her2_panel120.py` `MAX_PER_SCAFFOLD = 5`; `build_pik3ca_mtor_panel48.py` max 2; holdout builder Murcko cap 3.

---

## 2.4 Receptor and ligand preparation, docking, and scoring

### 6. Shared receptors and box geometry are missing

**Original**

> Human experimental structures and cognate ligands used for docking are listed in Tables 1 and S2. Binding sites were defined from the bound cognate ligands.

**Problem**

JAK1 6N7A is used for both JAK1 pairs; PPARA 6LXA is used for both PPARA pairs. The production box is cognate AABB + 5 Å per axis, then any edge < 20 Å is extended to 20 Å. Alternate locations without an altLoc flag, or marked A, were retained.

**Change to**

> Human experimental structures and cognate ligands used for docking are listed in Tables 1 and S2. JAK1 6N7A was used for both JAK1/TYK2 and JAK1/JAK2. PPARA 6LXA was used for both PPARG/PPARA and PPARA/PPARD. Binding sites were defined from the bound cognate ligands. The docking box was the Cartesian range of cognate heavy atoms, expanded by 5 Å on each side along x, y, and z; any edge that remained shorter than 20 Å was extended to 20 Å. Alternate locations without an altLoc flag, or marked A, were retained.

**Evidence:** Table S1 “Box: Cognate AABB + 5 Å; minimum edge 20 Å”; Table S2 (6N7A, 6LXA shared).

### 7. Ligand-prep recipe is incomplete

**Original**

> Processed ChEMBL SMILES served as ligand input. RDKit added hydrogens and generated one three-dimensional conformer with ETKDGv3. Ligands were geometry-optimized with MMFF and converted to PDBQT with Meeko.

**Problem**

The frozen recipe is: largest organic fragment → AddHs → ETKDGv3 seed 20260727 → MMFF ≤ 200 (failures skipped) → Meeko PDBQT. Protonation states and tautomers were not enumerated. Unspecified stereochemistry followed RDKit defaults.

**Change to**

> Processed ChEMBL SMILES served as ligand input. The largest organic fragment was retained, RDKit added hydrogens, and ETKDGv3 generated one three-dimensional conformer with seed 20260727. Ligands were geometry-optimized with MMFF for at most 200 iterations and converted to PDBQT with Meeko. Protonation states and tautomers were not enumerated systematically; unspecified stereochemistry followed RDKit’s default handling of the given SMILES.

**Evidence:** `prep_track_b_ligands_v1.py`; Table S1 ligand-prep row; `ENV_PIN.md`.

### 8. `energy_range` is omitted; two receptors saved eight poses

**Original**

> Up to nine poses were retained for each ligand–receptor pair. Exhaustiveness was 16 for PIK3CA/mTOR and 8 for the other target pairs.

**Problem**

Production docking used `num_modes = 9` and `energy_range = 3` kcal mol⁻¹. AChE 4EY7 and TYK2 3LXP saved 8 poses, so their lowest RMSD is not a uniform best-of-nine. PIK3CA/mTOR used exhaustiveness 16 because mTOR 4JT6 failed the < 2.0 Å search-coverage gate at exhaustiveness 8.

**Change to**

> At most nine poses were requested for each ligand–receptor pair, with `energy_range` set to 3 kcal mol⁻¹. Exhaustiveness was 16 for PIK3CA/mTOR, because mTOR 4JT6 failed the near-native search-coverage gate at exhaustiveness 8, and 8 for the other target pairs. AChE 4EY7 and TYK2 3LXP retained eight saved poses.

**Evidence:** `ENV_PIN.md`; `all14_cognate_rmsd_calcrrms_v1.csv` (`n_modes` 8 for 4EY7 and 3LXP); Table S2 note on 4JT6.

### 9. The 2.0 Å cutoff is search coverage, not top-1 ranking

**Original**

> The RMSD of the top-ranked pose and the lowest RMSD among saved poses were recorded, with 2.0 Å used as the near-native recovery cutoff.

**Problem**

The receptor gate is lowest saved-pose RMSD < 2.0 Å (search coverage), not top-1 RMSD < 2.0 Å. Several primary receptors fail top-1 and still pass the gate.

**Change to**

> The RMSD of the top-ranked pose and the lowest RMSD among saved poses were recorded. The near-native recovery cutoff was 2.0 Å and was applied to the lowest RMSD among saved poses (search coverage), not to the top-ranked pose.

**Evidence:** Table S1 “Cognate gate: lowest heavy-atom RMSD among all saved poses < 2.0 Å”; `all14_cognate_rmsd_calcrrms_v1.csv` (`pass_best_lt2` vs `pass_top1_lt2`).

### 10. Independent GNINA is not a complete-case copy of the Vina panels

**Original**

> Independent GNINA docking was also performed for EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2 using the same receptors, ligand sets, and docking boxes as the primary analysis. Detailed scoring settings, docking parameters, and complete-case sample sizes are provided in Table S7.

**Problem**

Table S7 is the correct SI table for this analysis. Independent GNINA did not return both-end scores for every ligand, so complete-case n differs from Vina. Exhaustiveness followed the Vina settings (8 / 16 / 8). Primary Vina knobs belong in Table S1, not S7. Independent GNINA used rank-1 `minimizedAffinity` (sign-flipped for statistics).

**Change to**

> Independent GNINA 1.3.2 docking was also performed for EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2 using the same receptors, ligand inputs, docking boxes, and exhaustiveness settings as the primary analysis (8, 16, and 8). Rank-1 minimizedAffinity was used as the ligand-level score. Independent GNINA did not return both-end scores for every ligand: complete-case counts were 28/38/32/11 (EGFR/HER2), 18/13/12/4 (PIK3CA/mTOR), and 30/32/29/14 (JAK1/TYK2). Detailed settings and counts are in Table S7. Primary Vina parameters are in Table S1.

**Evidence:** `independent_dock_formulation_v1.csv`; Table S7a; `analyze_independent_dock_v1.py`.

---

## 2.5 Evaluation metrics and statistical analysis

### 11. Dual-versus-neither mixes two different estimands

**Original**

> To isolate the effect of changing the comparison population from changes in score aggregation or target-score assignment, the same target-specific docking score was used for dual-versus-single-target-active and dual-versus-neither comparisons. Pocket A scores were used for dual-versus-B-only and dual-versus-neither, whereas pocket B scores were used for dual-versus-A-only and dual-versus-neither.
>
> A two-pocket mean score, Smean=SA+SB2， was reported only as a descriptive aggregate ranking for dual-versus-neither and was not used as the primary directional endpoint.

**Problem**

These are two analyses:

1. Fixed-score dual-versus-neither, pocket held fixed — Table S4. This is the isolation analysis.
2. Two-pocket mean \(S_{\mathrm{mean}}=(S_A+S_B)/2\) — Table 3. This is the conventional dual-versus-neither ranking. It is not a directional endpoint, but it is the locked Table 3 estimand, not an unnamed extra.

The formula text is also garbled (`Smean=SA+SB2` and the Chinese comma).

**Change to**

> To isolate the effect of changing the comparison population from changes in score aggregation or target-score assignment, the same target-specific docking score was used for dual-versus-single-target-active and dual-versus-neither comparisons (Table S4). Pocket A scores were used for dual-versus-B-only and for the paired dual-versus-neither contrast; pocket B scores were used for dual-versus-A-only and for the paired dual-versus-neither contrast.
>
> A two-pocket mean score, \(S_{\mathrm{mean}}=(S_A+S_B)/2\), was the conventional dual-versus-neither ranking reported in Table 3. It was not used as a primary directional endpoint.

**Evidence:** `STATISTICAL_LOCK_V1.md` Table 3 = `vina_mean`; `formulation_equal_score_negative_v1.csv`; `analyze_five_pair_stack_v1.py` `boot_single` on `vina_mean`.

### 12. `summary_min` formula is unreadable

**Original**

> The lower of the two directional AUROCs was reported as a descriptive weaker-arm summary,summarymin=minAUCD/A(B),AUCD/B(A).This value was used only to summarize the weaker of the two directional results.

**Change to**

> The lower of the two directional AUROCs was reported as a descriptive weaker-arm summary, \(\mathrm{summary}_{\min}=\min(\mathrm{AUROC}_{D/A}(B),\mathrm{AUROC}_{D/B}(A))\). This value was used only to summarize the weaker of the two directional results and is not a separate scoring function.

**Evidence:** `STATISTICAL_LOCK_V1.md`; `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json`.

### 13. Dual-versus-neither bootstrap needs to name Table 3

**Original**

> Dual-versus-neither analyses used class-stratified bootstrap resampling.

**Problem**

True for Table 3 (`boot_single`: resample dual and neither within class) and for Table S4 Δ. Readers may think it also applies to Table 2; Table 2 is pooled non-stratified. The class-stratified `summary_min` check is archived and does not replace Table 2.

**Change to**

> Table 3 dual-versus-neither confidence intervals used class-stratified resampling of dual and neither ligands. The class-stratified `summary_min` analysis is a sensitivity check and does not replace the Table 2 intervals.

**Evidence:** `analyze_five_pair_stack_v1.py` `boot_pm_ci` vs `boot_single`; Table S1 bootstrap row.

---

## 2.6 Diagnostic and sensitivity analyses

### 14. Logistic-regression settings are missing

**Original**

> ECFP4, docking-only, and ECFP4 plus docking models were evaluated with logistic regression using the same Bemis-Murcko scaffold-grouped cross-validation splits. Out-of-fold predictions were used to calculate AUROC.

**Problem**

Primary models are unscaled `LogisticRegression(C=1.0, max_iter=4000)` with default non-shuffled GroupKFold. The number of splits is min(5, n_scaffolds, the two class counts). StandardScaler on each training fold is sensitivity only (Table S5).

**Change to**

> ECFP4, docking-only, and ECFP4 plus docking models were evaluated with unscaled logistic regression (C = 1.0, maximum 4000 iterations) using the same Bemis–Murcko scaffold-grouped GroupKFold splits. The number of splits was the minimum of five, the number of scaffolds, and the two class counts. Out-of-fold predictions were used to calculate AUROC. A sensitivity analysis fitted StandardScaler on each training fold of the same splits (Table S5).

**Evidence:** Table S1 GroupKFold/logistic row; `ecfp4_docking_scaler_sensitivity_v1.py`.

### 15. Holdout, five-seed, and exhaustiveness scope are overstated

**Original**

> Robustness analyses included alternative receptor structures for PIK3CA/mTOR, internal holdout ligands, five fixed Vina random seeds, and alternative exhaustiveness settings. Internal holdout sets were constructed from unused candidates after main-panel ligands were removed.

**Problem**

- EGFR/HER2 has no holdout.
- Holdout panels are directional only (dual / A-only / B-only); they do not contain a neither class.
- Alternate crystals are PIK3CA/mTOR only (4JPS, 5DXT, 4JSX), one pocket at a time (Table S7c).
- Exhaustiveness 8 vs 16 is the PIK3CA/mTOR PM48 check.
- Five-seed Vina (20260727, 20260811–20260814) is repository QC, not a typeset SI table.

**Change to**

> Robustness analyses included alternative receptor structures for PIK3CA/mTOR (Table S7), unused-pool holdout ligands (Table S6), five fixed Vina random seeds (repository archive), and exhaustiveness 8 versus 16 on PIK3CA/mTOR. Internal holdout sets were constructed from unused candidates after main-panel ligands were removed for the seven pairs with remaining supply; EGFR/HER2 has no holdout. Holdout panels contain dual, A-only, and B-only ligands only.

**Evidence:** Table S6 (EGFR no holdout row; holdout n columns have no neither); Table S7c; `STATISTICAL_LOCK_V1.md` five-seed archive note.

---

## 2.7 External evaluation data

### 16. BindingDB processing is not the ChEMBL max-pChEMBL recipe

**Original**

> BindingDB was used to construct candidate external evaluation sets for all eight target pairs. PubChem was searched separately as an additional check on paired-data availability. Experimental states were assigned using θ=6.0. The development set included the primary panels, the expanded PIK3CA/mTOR panel, and the internal holdout sets. Molecules sharing literature sources with the development set or duplicating development-set structures were removed. Molecules with an ECFP4 Tanimoto similarity ≥ 0.70 to any development molecule were also excluded. A target pair was eligible for external evaluation only if the dual, A-only, and B-only classes each contained at least 20 ligands from at least three independent sources.

**Problem**

The independence-filtered remainder is Table S8. BindingDB used a different activity contract than ChEMBL Table 2: human wild-type, exact IC50/Ki/Kd, median aggregation, missing not inactive. Structure filtering also dropped molecules sharing ChEMBL IDs or InChIKeys with the pair’s θ = 6.0 ChEMBL maps, not only development-set structures. An additional gate required that no single document contribute more than 50% of a class. PubChem was a paired-data availability check and is not the Table S8 remainder. No pair passed, so no external docking was performed.

**Change to**

> BindingDB was used to construct candidate external evaluation sets for all eight target pairs (Table S8). PubChem was searched separately as an additional check on paired-data availability and was not used as the independence-filtered remainder. BindingDB experimental states were assigned at θ = 6.0 from exact human wild-type IC50, Ki, or Kd values, using median aggregation within a ligand–target pair; missing measurements were not treated as low activity. The development set included the primary panels, the expanded PIK3CA/mTOR PM110 panel, and the internal holdout sets. Molecules sharing literature sources with the development set, duplicating development-set or same-pair ChEMBL-map structures, or having ECFP4 Tanimoto similarity ≥ 0.70 to any development molecule were removed. A target pair was eligible for external evaluation only if the dual, A-only, and B-only classes each contained at least 20 ligands from at least three independent sources. No pair met those gates, and no external docking was performed.

**Evidence:** `external_slice_contract.yaml`; `external_slice_summary_v1.csv` (`packaged_as_external_evaluation = 0`); `bindingdb_native_slice_eight_pairs_v1.py`.

---

## Optional wording (not a data error)

These sentences are not false, but they hide protocol details that the SI already typesets. Add them only if the journal wants a self-contained Methods:

- 2.3 screening paragraph does not list later census gates (qHTS hubs, CYP ADME panels, zinc enzymes, ≥5 human holo structures per end, molecular-size filters). Those are in Table S9.
- 2.4.2 “RTMScore and GNINA 1.3.2 CNN models were used to rescore Vina-generated poses” is true; PPARG reporting used RTMScore over all saved poses and GNINA CNN affinity (Table S7b). Independent GNINA in the next sentence is a new search, not that rescore.
- 2.6.2 citing Table S6 for pocket correspondence is correct; S6 also contains the holdout table.

---

## Sentences that can stay

2.1 four-state assignment, pocket-matched directional arms, Figure 1A,B, dual as the positive class, and the role of ligand-only / pocket-correspondence / sensitivity analyses.  
2.2 four-state inequalities, θ = 6.0 primary labels, max pChEMBL as the primary aggregate, salt fragmentation, both-end measurements required, missing ≠ inactive.  
2.3 n ≥ 10 directional gate, 6.5/5.5 with 50/50 selectives, eight-pair list, post-construction relabel at θ = 6.0, both-end scores required for directional AUROC.  
2.4 PG08-NL retained in 9V8H; Meeko receptor prep; Vina 1.2.7 default scoring function; mode-1 affinity as the primary ligand score; CalcRMS as the formal RMSD.  
2.5 sign flip of Vina affinity; B = 2000; pooled non-stratified Table 2 bootstrap with shared dual resamples; omit empty-class replicates; 2.5th/97.5th percentiles; paired bootstrap for method or pocket swaps; cluster bootstrap of the two largest fixed-score differences; pointwise unadjusted CIs.  
2.6 MW, heavy atoms, cLogP, TPSA, ECFP4 radius 2 / 2048 bits; best-descriptor-on-same-panel treated as descriptive; wrong-pocket swap without redocking; Δ from `summary_min`; 2018 earliest-document year subset reported only when both directional arms met sample requirements.  
2.7 closing pointer to Table S1.
