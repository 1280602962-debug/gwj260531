# FINAL FREEZE REPORT

Branch: `cursor/methods-sentence-audit-c7cc`
Directory: `Dual_Target_Docking/`
Date: 2026-09-17
Verdict: **READY_TO_FREEZE**

Scientific numbers below come only from current scripts, current inputs, independent recomputation, and re-checked `results/canonical` tables. Manuscript/README prose was not used as a source of results.

---

## A. CURRENT SCIENTIFIC CHAIN

```
panel membership + pChEMBL
    EGFR: data/egfr_her2_panel120_v0/tables/panel_v0_120.csv
          + ablation_ligand_scores.csv (vina_3POZ_hb, vina_3RCD_hb; 3POZ/3RCD heavy-atom box)
    AChE: data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv
          (vina_ACHE_hb, vina_BCHE_hb; no-ChEMBL-ID-prefix panel; includes AB_056)
    six Track B / PM48 pairs: data/jcim_novelty_v0/tables/review_scored_membership_v1.csv
        → results/canonical/current_score_master.csv  (n=799, analysis_set=main only)
        → class_counts.csv
        → primary_directional_auroc.csv
        → primary_summary_min.csv
        → fixed_score_negative_class_delta.csv
        → two_pocket_mean_ranking.csv
        → top10_operating_points.csv
        → matched_minus_mismatched.csv
        → figure-source CSVs (unified_threshold theta_6.0; five_pair table2)
        → figures/jcim_article (plot_jcim_article_figures_v3.py)
```

Zero-dock recomputation this freeze: labels already in panels; AUROC; class-stratified bootstrap B=2000 seed=20260729; ranking; max-vs-median point/CI overlay; scaffold/document cluster (current scores); equal-score cluster (EGFR rebuilt; JAK scaffold rebuilt); detectable-effect 8-pair simulation; time-split AUROC overlay; figure regen.

Docking was not repeated. Production Vina outputs were reused.

Score convention: stored `score_A`/`score_B` are already favorable-higher (`*_hb` or negated raw affinity). Dual-vs-A_only always uses pocket B; dual-vs-B_only always uses pocket A; `summary_min` is min of those two.

---

## B. CURRENT TARGET PAIRS

Confirmed from `current_score_master` (not hard-coded). Production classes:

| pair | n_dual | n_A_only | n_B_only | n_neither | n_complete_case |
|---|---:|---:|---:|---:|---:|
| EGFR/HER2 | 28 | 38 | 32 | 12 | 110 |
| JAK1/JAK2 | 32 | 32 | 32 | 14 | 110 |
| JAK1/TYK2 | 31 | 32 | 32 | 14 | 109 |
| PIK3CA/mTOR | 18 | 14 | 12 | 4 | 48 |
| AChE/BChE | 27 | 26 | 28 | 15 | 96 |
| F2/F10 | 31 | 32 | 32 | 12 | 107 |
| PPARG/PPARA | 32 | 31 | 32 | 14 | 109 |
| PPARA/PPARD | 32 | 32 | 32 | 14 | 110 |

AChE production labels: strict 6.5/5.5. All other primary pairs: θ=6.0. `class` vs `class_from_pchembl`: 0 mismatches (AChE gray not counted as mismatch). Duplicate `(pair, ligand_id)`: 0. Missing scores: 0. Holdout is not in this master (`analysis_set=main`).

Withdrawn PIK3CA/PIK3CB is not a current pair.

---

## C. PRIMARY RESULTS VERIFIED

Independent SciPy Mann–Whitney U/(n_pos n_neg) on `current_score_master` matches freeze pairwise AUROC to 4 decimals for all eight pairs. Two implementations inside `recompute_current_primary.py` (pairwise vs midrank) also agree. Seed 20260729 bootstrap percentile CI for EGFR `summary_min` is identical on repeat: [0.1919, 0.4632].

| pair | D vs A (pocket B) | D vs B (pocket A) | summary_min [CI] | Δ pocket A (neither−selective) | matched−mismatched [CI excludes 0?] | EF_dual,10% |
|---|---:|---:|---|---:|---|---:|
| EGFR/HER2 | 0.6607 | 0.3237 | 0.3237 [0.1919, 0.4632] | 0.462 | 0.0558 [no] | 0.357 |
| JAK1/JAK2 | 0.5884 | 0.7275 | 0.5884 [0.4482, 0.7158] | −0.004 | −0.019 [no] | 1.875 |
| JAK1/TYK2 | 0.5751 | 0.3649 | 0.3649 [0.2329, 0.5051] | 0.444 | −0.065 [no] | 0.320 |
| PIK3CA/mTOR | 0.7143 | 0.6921 | 0.6921 [0.4801, 0.8016] | −0.220 | 0.090 [no] | 2.133 |
| AChE/BChE | 0.6524 | 0.6058 | 0.6058 [0.4392, 0.7354] | −0.016 | 0.177 [yes] | 1.778 |
| F2/F10 | 0.4133 | 0.3448 | 0.3448 [0.2157, 0.4819] | 0.163 | −0.031 [no] | 1.255 |
| PPARG/PPARA | 0.6492 | 0.7061 | 0.6492 [0.5111, 0.7460] | 0.053 | 0.030 [no] | 2.168 |
| PPARA/PPARD | 0.6465 | 0.4463 | 0.4463 [0.3008, 0.5898] | 0.038 | 0.012 [no] | 1.562 |

Do not claim an EGFR matched-pocket advantage: the 0.0558 CI includes 0.

Table 2 figure-source CSVs match these class-stratified CIs (`scripts/primary/bootstrap_primary.py` PASS). Prior non-stratified pooled CIs are superseded and are no longer the publication interval.

---

## D. SUPPORTING ANALYSES VERIFIED

| analysis | verdict | evidence |
|---|---|---|
| Label θ grid (EGFR, PIK3CA/mTOR, five-pair) | PASS | Existing grid retained; θ=6.0 rows patched to current scores/CIs |
| Label θ grid (AChE) | PASS with limitation | All AChE unified rows are production 27/26/28. Not a real θ reassignment. SI must not call this AChE θ-sensitivity |
| max vs median | PASS | See § known issue 1. AChE max membership now 27/26/28 |
| Scaffold cluster | PASS | Rebuilt on current scores. EGFR D/B point 0.3237; AChE D/A 0.6524 after AB_056 |
| Document cluster (EGFR, AChE, PM) | PASS with provenance note | Grouping rebuilt from committed `high_confidence_activity_audit_v1.csv`, not live ChEMBL 37 sqlite. Scores overlaid from master. AB_056 is a singleton group (no harvested documents) |
| Equal-score cluster EGFR | PASS | Δ=0.4621, document and scaffold CIs exclude 0 |
| Equal-score cluster JAK1/TYK2 scaffold | PASS | Rebuilt from current Track B scores/SMILES; Δ=0.4438 |
| Equal-score cluster JAK1/TYK2 document | PASS with limitation | Point Δ uses current JAK scores. CI kept from committed grouping; sqlite absent so grouping was not rebuilt |
| Descriptor controls | PASS as current SI | PhysChem from current ligand structures in existing tables; AChE fig3C n=27/26/28 |
| ECFP4 / incremental | PASS | Radius 2, 2048 bits, LR C=1. Ligand-only EGFR D/B remains frozen 0.8895. Dock ranks patched to current Vina. GroupKFold: 0 scaffold overlap and 0 ligand overlap on EGFR and AChE D/B. Sensitivity scaler is Pipeline-on-train-fold only |
| Incremental docking information | PASS | Same GroupKFold / OOF framework; compares ligand-only vs ligand-only+dock |
| Receptor substitution | PASS as current SI | Alternate PDB replaces one pocket; other pocket from primary. Fig4B 4JPS/5DXT/4JSX |
| GNINA independent pose generation | PASS | Distinct from CNN rescoring of Vina poses. Authoritative EGFR `summary_min`=0.2645 from `independent_dock_formulation_v1.csv` (MASTER 0.2199 was stale and was synced) |
| GNINA/RTM rescoring of Vina poses | PASS as SI | Separate analysis; must not be written as the same experiment as independent dock |
| Exhaustiveness / multiseed | PASS as SI | Same ligand set/receptor/labels in existing five-seed and E=8 tables; not re-docked |
| Cognate RMSD | PASS as QC | `all14_cognate_rmsd_calcrrms_v1.csv` is the current table. Not re-run this freeze (requires docking/poses) |
| Holdout | PASS | Holdout ligands are not in `current_score_master`. Not used for θ, receptor, or ECFP hyperparameter selection |
| Detectable-effect | PASS | 8-pair simulation from `class_counts.csv`; AChE 27/26/28; no PIK3CB; not observed power |
| Time-split | PASS as SI feasibility | 2018 cutoff: 0 evaluable pairs. Not external validation. Cutoff not moved |
| BindingDB | PASS as negative feasibility | 0 eligible pairs; not external validation |

---

## E. STALE ANALYSES

Obsolete / superseded (do not delete in this freeze; do not cite as current):

- PIK3CA/PIK3CB panel scripts and any leftover PIK3CB figure ticks (`figures/jcim_article/scripts/plot_jcim_si_composites_v1.py`)
- Pre-correction EGFR scores 10.605/9.745 (legacy ablation / old `canonical_ligand_level` before overlay)
- Old AChE complete-case 27/25/28 and docking census listing AB_056 as a failure
- Old three-pair detectable-effect + meta text claiming non-stratified Table 2
- `plot_jcim_article_figures_v1.py` / `v2.py`
- `build_checksum_manifest_v1.py` and other release-engineering wrappers (not part of the scientific chain)
- MCL1/BCL-xL and other non-eight-pair docking campaigns
- `submission_pack/` copies not rebuilt this freeze (auxiliary)

---

## F. DUPLICATE SCRIPTS

| item | keep | why |
|---|---|---|
| Primary AUROC + bootstrap | `scripts/freeze/recompute_current_primary.py` → `results/canonical/*` | Single current generator; independent of sklearn |
| Table 2 lock | `scripts/primary/bootstrap_primary.py` | Now verifies figure-source CSVs against canonical; no longer hard-codes 0.4297 |
| Article figures | `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py` | v1/v2 SUPERSEDED |
| Detectable-effect | `detectable_effect_simulation_v1.py` | Reads `class_counts.csv`; 8 pairs |
| Cluster | `scaffold_cluster_bootstrap_v1.py`, `document_blocked_cv_v1.py`, `equal_score_cluster_bootstrap_v1.py` | KEEP_CURRENT; JAK document = KEEP with provenance limitation |
| max vs median | `assay_aggregation_max_vs_median_v1.py` | KEEP_CURRENT (hb columns); freeze overlay repaired AB_056 without API |
| Time-split | `time_split_validation_v1.py` | KEEP as SI feasibility, not current eight-pair supporting |
| Checksum/pack scripts | SUPERSEDED / UNUSED for freeze | Do not expand |

---

## G. DATA VERSION CONFLICTS

Resolved this freeze:

- Ligand-level EGFR scores: membership, `canonical_ligand_level_v1.csv`, `high_confidence_labels_v1.csv` now match current 3POZ/3RCD hb (EH40_01 = 10.438/9.917).
- AChE AB_056 in master, membership, max-vs-median, cluster labels, class counts 27/26/28.
- Table 2 CIs unified to class-stratified shared-dual.
- Detectable-effect class sizes = current eight-pair counts; PIK3CB removed.
- MASTER GNINA EGFR rows synced to `independent_dock_formulation_v1.csv` (0.2645, not 0.2199).
- EGFR API-max `summary_min` on current scores is 0.3136, not the old 0.417 (old value mixed stale scores with API labels).

Still present, not publication-authoritative:

- `data/jcim_novelty_v0/tables/docking_failure_census_v1.csv` still lists AB_056 as fail.
- `submission_pack/` not regenerated.
- `jcim_figure_style.PAIR_ORDER` still the old three names; detectable SI plot now uses `PAIR_ORDER_DETECT` (eight pairs).

---

## H. UNRESOLVED ISSUES

Not BLOCKER. Concrete:

1. JAK1/TYK2 **document-cluster** grouping cannot be rebuilt here: `/tmp/chembl/chembl_37/.../chembl_37.db` is absent. WHAT CAN BE VERIFIED: point Δ=0.4438 from current Track B Vina scores; scaffold-cluster CI rebuilt. WHAT CANNOT: live sqlite harvest of JAK document connected components. WHETHER IT BLOCKS: no; ligand-level Δ is the official equal-score interval; document CI must be labeled as committed grouping.

2. AChE `unified_threshold_sensitivity_v2.csv` θ=5.5/6.0/6.5/strict rows are the production 27/26/28 panel, not re-labelled. Do not write AChE θ-sensitivity from this file.

3. ChEMBL activity API was not re-fetched for max-vs-median. AB_056 uses panel pChEMBL as max/median (single-value). Other ligands use the committed ligand-level API table joined to current scores.

4. `submission_pack/` and `docking_failure_census_v1.csv` remain stale relative to AB_056 success.

5. Cognate RMSD, holdout docking, GNINA independent search, receptor-swap docking, exhaustiveness, and multiseed were not re-executed (zero-dock freeze). Provenance is the existing pose/score tables.

6. ECFP4 ligand-only OOF for EGFR D/B is intentionally frozen at 0.8895 (pre-fix chemistry model). Vina rank AUROC is post-fix 0.3237.

---

## I. MANUSCRIPT READINESS

**READY_TO_FREEZE**

No unresolved BLOCKER against the 18 freeze conditions.

Suggested placement (not manuscript text):

- PRIMARY (Methods): eight pairs; θ=6.0 except AChE strict 6.5/5.5; complete-case; directional AUROC and `summary_min`; class-stratified bootstrap B=2000 seed=20260729; Vina 1.2.7 current receptors; holdout unused for selection.
- SUPPORTING (one sentence + SI): fixed-score Δ; matched vs mismatched (EGFR CI includes 0); Top-10%; max vs median; cluster bootstrap; ECFP4 increment; receptor swap; independent GNINA vs rescoring distinction; detectable-effect simulation.
- SI-ONLY: θ grid; exhaustiveness; multiseed; cognate RMSD; scaler sensitivity; document-blocked CV.
- REPOSITORY-ONLY / do not call validation: time-split 2018 (0 evaluable pairs); BindingDB feasibility (0 eligible pairs).
- OBSOLETE: PIK3CA/PIK3CB; pre-correction EGFR 0.4297; old AChE 27/25/28.

---

## J. FINAL AUTHORITATIVE FILE MAP

Current master data:

- `results/canonical/current_score_master.csv`

Primary result tables:

- `results/canonical/class_counts.csv`
- `results/canonical/primary_directional_auroc.csv`
- `results/canonical/primary_summary_min.csv`
- `results/canonical/fixed_score_negative_class_delta.csv`
- `results/canonical/two_pocket_mean_ranking.csv`
- `results/canonical/top10_operating_points.csv`
- `results/canonical/matched_minus_mismatched.csv`

Publication figure-source CSVs (must match canonical; v3 checksum PASS):

- `data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv` (θ=6.0)
- `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv`
- `data/jcim_novelty_v0/tables/formulation_*`
- `data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv`
- `data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv`

Supporting:

- `data/jcim_novelty_v0/tables/assay_max_vs_median_*.csv`
- `data/jcim_novelty_v0/tables/scaffold_cluster_bootstrap_v1.csv`
- `data/jcim_novelty_v0/tables/document_cluster_bootstrap_v1.csv`
- `data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv`
- `data/jcim_novelty_v0/tables/detectable_effect_simulation_v1.csv`
- `data/jcim_novelty_v0/tables/time_split_*.csv`
- `data/jcim_novelty_v0/analysis/BINDINGDB_EXTERNAL_FEASIBILITY_SI_V1.md`

Current scripts:

- `scripts/freeze/recompute_current_primary.py`
- `scripts/freeze/phase2_independent_audit.py`
- `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`
- `data/jcim_novelty_v0/scripts/detectable_effect_simulation_v1.py`
- cluster / max-vs-median / time-split scripts listed above

Inventory: `remediation_outputs/freeze_audit/phase1_analysis_inventory.csv`

---

## Fixes this freeze (objective only)

| problem | cause | files | numerical result changed? | interpretation changed? |
|---|---|---|---|---|
| No current_score_master; EGFR ligand scores still pre-fix 0.4297 | version desync | `scripts/freeze/recompute_current_primary.py`; canonical tables; membership; unified; formulation; MASTER; ligand-level | EGFR D/B 0.4297→0.3237; CIs now class-stratified | EGFR weak arm still chance; matched-pocket CI still includes 0 |
| AChE 26 vs 25 A-only in max-vs-median | AB_056 (CHEMBL3960861) missing from stale n=95 table | ligand table + auroc overlay; hb columns in aggregation script | AChE max membership 27/25/28→27/26/28; D/A 0.6504→0.6524 | Aggregation biology did not cause the missing A-only |
| Detectable-effect old 3 pairs + PIK3CB + 27/25/28 | hardcoded stale sizes | `detectable_effect_simulation_v1.py` + regenerated csv/meta | full 8-pair grid recomputed | Still simulation, not observed power |
| Cluster EGFR D/B 0.4297 | high_confidence vina columns stale | overlay + rebuild | EGFR cluster point 0.3237 | Weak arm remains near chance under cluster resample |
| Table 2 lock 0.4297 | hardcoded LOCKED | `bootstrap_primary.py` | n/a (now reads canonical) | n/a |
| Figure checksum CIs | old non-stratified | `plot_jcim_article_figures_v3.py` | plotted CIs = class-stratified | Interval slightly different; point AUROCs unchanged for six pairs |
| MASTER GNINA 0.2199 vs formulation 0.2645 | stale aggregator | MASTER_RESULTS_TABLE.csv | MASTER now 0.2645 | Independent-dock gap remains vs current Vina 0.3237 |

---

## Freeze-condition checklist

1. Eight pairs and membership determined — PASS  
2. Master traceable to current ablation / Track B — PASS  
3. Primary AUROC/delta/ranking recompute — PASS (independent MWU)  
4. Bootstrap definition unified for publication Table 2 — PASS (class-stratified, shared dual, B=2000, seed=20260729)  
5. Corrected EGFR and AChE throughout current analyses — PASS  
6. max-vs-median population mismatch explained and repaired — PASS  
7. detectable-effect synced to current eight-pair sizes — PASS  
8. time-split classified SI feasibility, not external validation — PASS  
9. ECFP4/descriptor no GroupKFold leakage — PASS  
10. Cluster bootstrap on current scores — PASS (JAK document grouping limitation in H.1)  
11. Receptor/GNINA/RTM/exhaustiveness/multiseed provenance — PASS as existing SI tables  
12. Holdout excluded from selection — PASS  
13. External analysis correctly bounded — PASS  
14. Publication figures read current figure-source CSVs — PASS (`v3 verification OK`)  
15. One current primary statistic generator — PASS  
16. Stochastic seed repeat — PASS  
17. No current figure depends on superseded 0.4297 / 27/25/28 / PIK3CB — PASS  
18. No unresolved BLOCKER — PASS  
