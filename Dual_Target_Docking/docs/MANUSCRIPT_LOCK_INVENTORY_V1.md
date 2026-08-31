# Manuscript lock inventory (no prose rewrite)

Date: 2026-08-31  
HEAD checked: `8d5f2eeb` (`cursor/pik3ca-mtor-structure-freeze-0b1a`)  
Scope: DualFourClass JCIM article assets. **Manuscript text was not edited.**

This note freezes what the article currently uses, what is coeval but conflicting, and what is legacy. Companion CSVs:

- `data/manuscript_lock/table2_ci_source_conflict_v1.csv`
- `data/manuscript_lock/table3_vs_multiseed_estimand_v1.csv`
- `data/manuscript_lock/ARTICLE_ASSET_INDEX_v1.csv`

---

## Verdict

Two reproducibility defects are real. Frozen docking scores and four-state labels are intact. Point AUROCs for Table 2 are identical across current executable tables. The problems are (1) **bootstrap CIs exist in more than one current executable stream**, and (2) **multi-seed Dual-versus-neither is a different estimand from Table 3**.

One correction to the previous diagnosis: `pocket_matched_directional_v1.csv` is **not newer** than the Table 2 lock file. Both CSVs were last committed **2026-08-26 14:47** (`531f55a1`). They differ because each script hashes a different `stable_offset` key into the bootstrap seed:

| Stream | Script | Offset key | Role today |
|---|---|---|---|
| Table 2 / Figure 3 / validator / MASTER | `build_t0_strengthen_v1.py` | `stable_offset(pair, "theta_6.0")` | **current manuscript lock** |
| Pocket-matched variant table (S6, Fig 6A) | `build_pocket_matched_diagnostics_v1.py` | `stable_offset(pair, "pocket_matched_vina")` | coeval executable; **do not overwrite Table 2 CIs yet** |
| `PRIMARY_METRIC_V2.md` | prose snapshot | `seed=20260729` with no offset | **deprecated** (2026-07-29) |

Both bootstrap implementations resample ligands **without class stratification** (`rng.choice` over the dual+A-only+B-only pool). Cluster bootstrap remains sensitivity only (`document_cluster_bootstrap_v1.csv`, `scaffold_cluster_bootstrap_v1.csv`).

Do **not** hand-replace the four manuscript CIs with the pocket-matched column. Figure 3, `plotted_values.json`, `MASTER_RESULTS_TABLE.csv`, and `validate_revision_v1.py` are already pinned to `unified_threshold_sensitivity_v2.csv`. Swapping one file would desynchronize figures and CI.

---

## 1. Table 2: three CI sets, one locked point estimate

Point estimates (identical in unified_threshold θ=6.0 and pocket_matched_vina):

| Pair | n (D/A/B) | D vs A (pocket B) | D vs B (pocket A) | summary_min |
|---|---:|---:|---:|---:|
| EGFR/HER2 | 28/38/32 | 0.6664 | 0.4297 | 0.4297 |
| AChE/BChE | 27/25/28 | 0.6504 | 0.6058 | 0.6058 |
| PIK3CA/PIK3CB | 28/27/28 | 0.6905 | 0.5000 | 0.5000 |
| PIK3CA/mTOR | 18/14/12 | 0.7143 | 0.6921 | 0.6921 |

CIs (B=2000 ligand bootstrap; all still include 0.5):

| Pair | PRIMARY_METRIC_V2.md (2026-07-29) | Manuscript / unified_threshold (lock) | pocket_matched_directional_v1.csv |
|---|---|---|---|
| EGFR/HER2 | [0.281, 0.576] | **[0.2818, 0.5775]** → text 0.430 [0.282, 0.578] | [0.2840, 0.5774] |
| AChE/BChE | [0.442, 0.737] | **[0.4370, 0.7303]** → text 0.606 [0.437, 0.730] | [0.4447, 0.7431] |
| PIK3CA/PIK3CB | [0.340, 0.648] | **[0.3502, 0.6495]** → text 0.500 [0.350, 0.650] | [0.3444, 0.6522] |
| PIK3CA/mTOR | [0.457, 0.813] | **[0.4702, 0.8133]** → text 0.692 [0.470, 0.813] | [0.4722, 0.8000] |

Additional CI drift on the “original receptor” rows of Table S30 (`receptor_realization_two_pair_v1.csv`): PIK3CA/mTOR 4L23/4JT6 is [0.4638, 0.8015]; PIK3CA/PIK3CB 4L23/2WXF is [0.3468, 0.648]. Those rows should be used for **swap deltas**, not as a second Table 2.

Legacy pooled `vina_mean` directional (not primary): `bootstrap_directional_ci_v1.csv` EGFR summary_min **0.2824** [0.1577, 0.4219]. Same numbers live in `forest_summary_min_ci_v1.csv` vina_mean rows. Keep as a control; do not treat the filename as Table 2.

Descriptor columns currently in Table 2 (from `descriptor_all_four_directional_v1.csv`, rounded to 3 d.p.):

| Pair | heavy | MW | cLogP | TPSA | best single |
|---|---:|---:|---:|---:|---|
| EGFR/HER2 | 0.369 | 0.416 | 0.482 | 0.427 | cLogP 0.4821 |
| AChE/BChE | 0.582 | 0.579 | 0.467 | 0.733 | TPSA 0.7333 |
| PIK3CA/PIK3CB | 0.622 | 0.620 | 0.595 | 0.418 | heavy 0.6217 |
| PIK3CA/mTOR | 0.463 | 0.448 | 0.310 | 0.260 | heavy 0.463 |

---

## 2. Table 3 vs multi-seed: different Dual-versus-neither estimands

Table 3 lock (`formulation_conventional_vs_directional_v1.csv`): per-ligand `vina_mean = (S_A+S_B)/2`, then one AUROC. Recomputed from the four frozen ablation score tables: **exact match**.

| Pair | n_neither | Table 3 AUC(vina_mean) | 95% CI | Dual vs all non-duals |
|---|---:|---:|---|---:|
| EGFR/HER2 | 12 | 0.7560 | [0.5625, 0.9197] | 0.5514 [0.4429, 0.6664] |
| AChE/BChE | 15 | 0.6494 | [0.4840, 0.8123] | 0.5792 [0.4417, 0.7157] |
| PIK3CA/PIK3CB | 16 | 0.5592 | [0.3728, 0.7456] | 0.5558 [0.4369, 0.6720] |
| PIK3CA/mTOR | 4 (underpowered) | 0.5139 | [0.2222, 0.8056] | 0.6741 [0.5148, 0.8167] |

`analyze_multiseed_vina_v1.py` instead computes `mean(AUC_A, AUC_B)` on Dual vs neither. Recomputed on the same primary-seed scores:

| Pair | AUC(vina_mean) = Table 3 | mean(AUC_A, AUC_B) = multi-seed column | Δ |
|---|---:|---:|---:|
| EGFR/HER2 | 0.7560 | 0.7641 | +0.0081 |
| AChE/BChE | 0.6494 | 0.6494 | 0 (coincidence) |
| PIK3CA/PIK3CB | 0.5592 | 0.5798 | +0.0206 |
| PIK3CA/mTOR | 0.5139 | 0.5278 | +0.0139 |

Independent GNINA Table S32 already uses `score_mean` (ligand-level mean, then one AUC): EGFR Dual vs neither 0.7825 [0.6104, 0.9222]; directional summary_min 0.2199; PIK3CA/mTOR summary_min 0.6325. That matches the Table 3 definition, not the multi-seed definition.

Multi-seed **directional** columns (D vs A, D vs B, summary_min) are pocket-matched and valid. Only the Dual-vs-neither / gap columns are the wrong estimand. Scores in `multiseed_scores_long_v1.csv` do not need re-docking; the analyzer does.

---

## 3. Frozen inputs (do not refresh because they are old)

| Pair | Receptor A/B | Exhaustiveness | n_panel | n_scored D/A/B | Score table |
|---|---|---:|---:|---|---|
| PIK3CA/mTOR | 4L23 / 4JT6 | 16 | 48 | 18/14/12 | `data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv` |
| AChE/BChE | 4EY7 / 4BDS | 8 | 100 | 27/25/28 | `data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv` |
| PIK3CA/PIK3CB | 4L23 / 2WXF | 8 | 100 | 28/27/28 | `data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv` |
| EGFR/HER2 | 3POZ / 3RCD | 8 | 110 | 28/38/32 | `data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv` |

Primary labels: unified θ = 6.0, max pChEMBL, freeze date 2026-07-23. 2026-08-26 API refetch is sensitivity (`assay_max_vs_median_auroc_v1.csv`; EGFR 0.417 vs 0.424). High-confidence view: 352/352 class match (`high_confidence_summary_v1.csv`).

---

## 4. Main-text numbers that already have a single CSV source

These are consistent with MASTER / validator / current Results (3 d.p. rounding in prose):

- Formulation gap, fixed-score negative class: EGFR pocket A Δ = 0.3783 [0.2050, 0.5469] (`formulation_equal_score_negative_v1.csv`).
- Mixed-library EGFR `vina_mean` Top-10: 1 dual / 5 A-only / 4 B-only; hard-negative fraction 0.90 (`mixed_library_enrichment_v1.csv`).
- Complete-case fractions: 0.145119 / 0.340172 / 0.233349 / 0.265252 (`complete_case_usable_pchembl_overlap_v1.csv`) → text 14.5%–34.0%.
- Document-blocked: EGFR D/B still 0.4297; document-cluster CI [0.3214, 0.6171]; PIK3CA/mTOR D/B `cannot_stably_estimate`.
- Receptor points: PIK3CA/mTOR 0.6921 → 0.4861 / 0.5046 (4JPS/5DXT); 4JSX 0.6389. PIK3CA/PIK3CB 0.5000 → 0.6905 / 0.6849.
- AND filter EGFR Dual-median `vina_worst`: precision 0.2979, hard-negative fraction 0.7021.
- Full-map ECFP4 EGFR: Dual vs neither 0.9214; Dual vs B-only 0.8636.
- θ=6.0 census: 49 pairs, 17 with directional n≥10, 4 docked.
- Assay-context: 179 include / 7 uncertain / 0 exclude. Ligand-level source-reading ledger is **not** in the article (`ASSAY_CONTEXT_SOURCE_READING_V1.md`).
- BindingDB native slice: zero pairs through the pre-frozen gate; not docked.

---

## 5. Figures (submission files)

All live in `figures/jcim_article/`. Captions: `figures/jcim_article/CAPTIONS.md`. Numeric lock: `figures/jcim_article/plotted_values.json` (Vina CIs = unified_threshold).

| Figure | File stem | Reads |
|---|---|---|
| 1 | Fig1_task_schematic | schematic |
| 2 | Fig2_hardneg_supply | J0 supply + Table S12 |
| 3 | Fig3_formulation_comparison | unified_threshold + formulation CSV |
| 4 | Fig4_confounds | unified_threshold + AChE assembled scores |
| 5 | Fig5_receptor_realization | receptor_realization_two_pair + alt CSVs |
| 6 | Fig6_wrong_pocket_paradox | pocket_matched_directional + holdout |
| 7 | Fig7_confound_anatomy | ML / matched-subset tables |
| 8 | Fig8_diagnostic_workflow | schematic |
| TOC | TOC_graphic | no AUROCs |
| S1 | FigS1_protocol_sensitivity | unified_threshold grid + GNINA rescore + PM110 |
| S2 | FigS2_equal_relation_and_sampling | crossdb + holdout shift |
| S3 | FigS3_paired_delta_bootstrap | paired-delta CSVs |
| S4 | FigS_pocket_matched_forest | unified Vina CIs + forest descriptor CIs |
| S5 | FigS_unused_pool_holdout | holdout_pocket_matched |
| S6 | FigS_detectable_effect | detectable_effect_simulation |
| S7 | FigS_formulation_upgrades_v1 | census / AND / full-map (**png only**) |
| S8 | FigS_bindingdb_native_slice_v1 | external_candidate_flow |

Regenerate command (does not change scores):  
`python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py`

---

## 6. SI table map (S1–S53)

The machine-readable map is the “溯源” table in `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` plus `ARTICLE_ASSET_INDEX_v1.csv`. English compressed outline: `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`.

MCL1/Bcl-xL (S50–S53) is SI-only after formal demotion: Vina Dual vs neither 0.628 [0.462, 0.786]; Dual vs A 0.793; Dual vs B 0.609; summary_min 0.609. Not Table 2.

---

## 7. Files that look official but are not Table 2

| File | Why it misleads | Keep? |
|---|---|---|
| `PRIMARY_METRIC_V2.md` | 2026-07-29 CIs; claims pocket_matched CSV as source but was not rebuilt | Keep as dated snapshot; mark deprecated |
| `bootstrap_directional_ci_v1.csv` | pooled vina_mean; EGFR 0.2824 | Keep as legacy control |
| `forest_summary_min_ci_v1.csv` vina_mean rows | same pooled 0.2824 | Keep for descriptor forest CIs |
| `directional_with_baselines_v1.csv` (if present) | pre-pocket-matched | Archive |
| `multiseed_auroc_by_seed_v1.csv` Dual-vs-neither column | mean of pocket AUCs | Keep scores; reanalyze metrics |
| `MANUSCRIPT_JCIM_*.md` | assembler output | Never edit by hand |
| `RESULTS_DRAFT_ZH_JCIM_V1.md` | latest prose, numbers from unified_threshold | Edit only after a statistical freeze |

---

## 8. Integration-branch status (2026-08-31)

Done on `cursor/jcim-final-integration-0b1a` without changing Table 2 CIs:

1. `docs/STATISTICAL_LOCK_V1.md` now defines the estimands on one page.
2. Table 2 CIs remain `unified_threshold_sensitivity_v2.csv`. A read-only loader is `scripts/primary/bootstrap_primary.py`. Rebuilding a shared hash-offset bootstrap (item 2 below) is still future work.
3. `analyze_multiseed_vina_v2.py` writes `AUC(vina_mean)` Dual-versus-neither and keeps `mean_marginal_pocket_auroc_D_vs_neither` as a named diagnostic. Primary seed recovered Table 3. Table S54 cites v2 only.
4. MASTER, checksum, and assembled manuscripts are regenerated on this branch. Figures were not redrawn because Table 2 CIs did not change.
5. Table 2 CIs in Results are unchanged.

Still later, and **not** this integration:

- Rebuild `unified_threshold_sensitivity_v2.csv` and `pocket_matched_directional_v1.csv` from one seed function so those two executable streams cannot fork.

Until that later freeze, the article continues to cite **unified_threshold θ=6.0** for Table 2 CIs and **formulation_conventional_vs_directional_v1.csv** for Table 3.

