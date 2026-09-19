# PR #39 scientific data audit

**VERDICT: READY_FOR_SCIENTIFIC_FREEZE**

Scope: scientific design, zero-dock reproducibility, pinned analysis environment, canonical computation, Table 1–3 / S1–S10 numbers, Figure 1–5 / S1–S5 sources. Manuscript prose was not reviewed. Two SI Table S7 interval cells (and the matching results-sentence digits) were aligned to the predetermined class-stratified bootstrap after same-pose RTM/CNN scores were recomputed; no wording was rewritten.

Base: `17435413410290960da3437de640509705f00b51`  
Working branch: `cursor/scientific-freeze-c7cc`

---

## Environment

| Item | Value |
|------|--------|
| Declared pins (`requirements-analysis.txt`) | numpy 2.5.2; pandas 3.0.5; scipy 1.18.1; scikit-learn 1.9.0; rdkit 2026.3.5; matplotlib 3.11.1; meeko 0.7.1; gemmi 0.7.5 |
| Actual freeze environment | `/tmp/dual_target_pr39_env` Python 3.12.3 |
| Final freeze packages | numpy 2.5.2; pandas 3.0.5; scipy 1.18.1; scikit-learn 1.9.0; rdkit 2026.03.5 (`importlib.metadata` 2026.3.5) |
| Exact-match status | **PASS** (`scripts/check_analysis_env.py` Result: `PIN_MATCH`, exit 0) |
| Unpinned system Python | **FAIL** (numpy 2.4.4; scikit-learn 1.9.1; rdkit 2026.03.6). Mismatch is no longer reported as COMPATIBLE. |
| matplotlib / Pillow | Figure-render only; not used to accept a scientific mismatch |
| Meeko / Gemmi | Optional; zero-dock chain does not call them |

`scripts/check_analysis_env.py` now fails (exit 1) if numpy, pandas, scipy, scikit-learn, or rdkit differ from the exact pin. No SHA/hash check was added.

Adjudication is step `00_adjudicate` in `scripts/analysis/rebuild_freeze.py`. The committed adjudication tables remain on disk; a full rebuild regenerates them rather than treating them as a hidden intermediate. Default rebuild output is `/tmp/dual_target_freeze_rebuild`.

---

## ECFP reproducibility

| Item | Value |
|------|--------|
| Base `17435413` max \|Δ\| | **0.0112** at JAK1/JAK2 D vs A |
| Pinned-environment rebuild max \|Δ\| | **0.0234** at PPARA/PPARD D vs B |
| Final canonical max \|Δ\| | **0.0234** at PPARA/PPARD D vs B |
| Fold rows compared vs base | 930 / 930 |
| Scaffold-string mismatch | **3** |
| Scaffold same, `fold_id` changed | **491** |
| Identical scaffold + fold | **436** |
| Same-env repeated run | `/tmp/pr39_ecfp_run_a` and `/tmp/pr39_ecfp_run_b` byte-identical for `model_fold_assignments.csv`, `ecfp4_oof_predictions.csv`, `ecfp4_incremental_information.csv` |

Decision rule: the repository-declared pin is truth. The pinned environment reproduced **0.0234**, not 0.0112 and not a third value. Canonical keeps 0.0234. `requirements-analysis.txt` was not changed to preserve either number.

The three scaffold-string changes are RDKit canonicalization of the same Bemis–Murcko graph (JAK1/JAK2 `J1J2_034`; JAK1/TYK2 `J1TYK2_065` on both arms). The 491 fold-only changes are GroupKFold / sklearn partition differences relative to the base commit. Independent OOF AUROC replay, Δ replay, unique-ligand, unsplit-scaffold, shared ligand rows, and shared fold IDs across ECFP4 / docking / ECFP4+docking all **PASS**.

JAK1/TYK2 best-descriptor `summary_min` under the pin is **0.5796** (typeset 0.580).

---

## Primary data

Eight-pair complete-case activity-eligible n (θ = 6.0; Vina mode-1; `score_S = −energy`; higher better):

| Pair | n_scored D / A / B | n_neither | weaker arm | summary_min [95% CI] |
|------|-------------------:|----------:|------------|----------------------|
| EGFR/HER2 | 28 / 37 / 32 | 12 | D vs B pocket A | 0.3237 [0.1953, 0.4744] |
| JAK1/JAK2 | 32 / 32 / 32 | 14 | D vs A pocket B | 0.5884 [0.4482, 0.7158] |
| JAK1/TYK2 | 31 / 32 / 32 | 14 | D vs B pocket A | 0.3649 [0.2329, 0.5051] |
| PIK3CA/mTOR | 18 / 14 / 12 | 4 | D vs B pocket A | 0.6921 [0.4801, 0.8016] |
| AChE/BChE | 27 / 26 / 28 | 14 | D vs B pocket A | 0.6058 [0.4392, 0.7354] |
| F2/F10 | 31 / 32 / 32 | 12 | D vs B pocket A | 0.3448 [0.2157, 0.4819] |
| PPARG/PPARA | 32 / 31 / 32 | 14 | D vs A pocket B | 0.6492 [0.5111, 0.7460] |
| PPARA/PPARD | 32 / 32 / 32 | 14 | D vs B pocket A | 0.4463 [0.3008, 0.5898] |

Independent replay from `current_score_master`: directional AUROC, shared-dual class-stratified CIs (B = 2000, seed 20260729), and `summary_min` all match canonical. EGFR quota 28 / 38 / 32 / 12 is not mixed with n_scored 28 / 37 / 32. AChE primary n_scored 27 / 26 / 28 is not mixed with the max-vs-median qualified intersection n = 94.

`verify_freeze_rebuild.py --outdir /tmp/pr39_pinned_rebuild` → **result: PASS**.

---

## Robustness

| Check | Result |
|-------|--------|
| Fixed-score 8 × 2 directions from master | PASS. EGFR pocket A: D vs B-only 0.3237, D vs neither 0.7857, Δ 0.4621. JAK1/TYK2 pocket A Δ 0.4438 |
| Ranking `score_mean`, tie by ligand ID, k = ceil(0.10 n) | PASS. Table 3 / Figure 2C / 2D / Table S10 share `two_pocket_mean_ranking.csv` and `top10_operating_points.csv` |
| AND filter | Canonical `and_filter_operating_points.csv`; neither excluded; threshold = median dual `score_worst` |
| Matched vs mismatched | Only AChE/BChE main-panel interval excludes 0. EGFR includes 0. Seven holdouts include 0. EGFR has no holdout |
| Independent GNINA | EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2 from deposited pose scores + current labels. Distinct from CNN rescoring of Vina poses |
| Receptor substitution | 4JPS / 5DXT / 4JSX. `primary_summary_min` copied from rebuilt Table 2 (0.6921), not hard-coded |
| Five-seed Vina | Canonical `five_seed_summary_min.csv`. EGFR/HER2 from corrected-box per-seed scores + current activity-eligible labels (n = 28/37/32/12). Production seed matches Table 2 (0.3237). Other pairs unchanged: deposited per-seed scores + current labels |
| Same-pose RTM / CNN | Recomputed from `scores_rtm_best9_v1.csv` and `scores_gnina_cnn_best9_v1.csv` into `computational_robustness.csv`. PPARG/PPARA RTM 0.3691 [0.2369, 0.4776]; CNN 0.5000 [0.3467, 0.6321]. Notes mark “not independent pose generation” |
| RMSD | 14 primary slots from unified CalcRMS. EGFR 3POZ top-1 1.019 Å; HER2 3RCD top-1 1.947 Å. 0.760 Å / 9.505 Å absent. JAK2, mTOR, BChE, PPARG, PPARA: top-1 ≥ 2 Å and lowest saved pose < 2 Å (search coverage, not top-1 ranking) |

---

## Target-pair selection

`data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv` has no AUROC, `summary_min`, matched/mismatched, GNINA, or ranking columns. Inclusion uses activity supply, target identity, structural availability, protocol compatibility, and site interpretation.

EGFR/HER2 is recorded as a **supply-limited exception**: it missed the strict 6.5/5.5 selective-supply gate (minimum selective-class count = 7) but had suitable human holo structures, a cognate-defined site, and enough θ = 6.0 dual / A-only / B-only ligands for directional evaluation. That is a design-time supply/structure decision, not an outcome-driven add. **No MAJOR DESIGN ISSUE.**

---

## Figure matrix

Official writer: `figures/jcim_article/scripts/update_figures_pr32.py`. Helper library `plot_jcim_article_figures_v3.py` no longer writes artwork or `plotted_values.json`.

| Figure | Panel | Authoritative source | Status |
|--------|-------|----------------------|--------|
| 1 | C | frozen `universe_census_summary_v1.csv` | PASS |
| 2 | A | `results/canonical/fixed_score_negative_class_delta.csv` | PASS |
| 2 | B | `results/canonical/primary_directional_auroc.csv` | PASS |
| 2 | C | `primary_summary_min.csv` + `two_pocket_mean_ranking.csv` | PASS |
| 2 | D | `results/canonical/top10_operating_points.csv` | PASS |
| 3 | A/B | `results/canonical/ecfp4_incremental_information.csv` | PASS |
| 4 | A | `results/canonical/matched_mismatched_pocket.csv` | PASS |
| 4 | B | `results/canonical/holdout_metrics.csv` | PASS |
| 5 | A | `computational_robustness.csv` (independent GNINA) + `primary_summary_min.csv` / `two_pocket_mean_ranking.csv` | PASS |
| 5 | B | `results/canonical/receptor_substitution.csv` | PASS |
| 5 | C | `results/canonical/five_seed_summary_min.csv` | PASS |
| S1 | A | `descriptor_baselines.csv` + `primary_summary_min.csv` | PASS |
| S1 | B | TPSA from `current_score_master` SMILES + RDKit | PASS |
| S2 | A/B | `results/canonical/protocol_sensitivity.csv` (PM110 and E8 rebuilt from frozen scores) | PASS |
| S3 | A | `results/canonical/cognate_rmsd.csv` | PASS |
| S4 | A | `results/canonical/label_aggregation_sensitivity.csv` | PASS |
| S4 | B | `results/canonical/cluster_bootstrap_sensitivity.csv` | PASS |
| S5 | A/B | `results/canonical/external_eligibility.csv` (copy of frozen eligibility input) | PASS |

`plotted_values_postfix.json`: 323 records; no `inputs_sha256` / `artwork_git_head` authority fields; every analysis-derived raw value is in the declared canonical or frozen census file. Sole plotted-value authority. `plotted_values.json` is not written.

Allowed sources only: `results/canonical/*` or an explicit immutable frozen experimental/eligibility input.

---

## Table matrix

Displayed values use three-decimal half-up rounding unless the cell is an integer count or an Å RMSD copied at the stored precision.

| Table | Authoritative source | Status |
|-------|----------------------|--------|
| 1 | class counts + construction quotas / PDB / exhaustiveness from the panel design record | PASS. n_scored ≠ quota (EGFR 28/37/32 vs 28/38/32/12) |
| 2 | `primary_directional_auroc.csv`, `primary_summary_min.csv`, `class_counts.csv` | PASS |
| 3 | `two_pocket_mean_ranking.csv`, `top10_operating_points.csv`, `class_counts.csv` | PASS |
| S1 | `requirements-analysis.txt` + docking protocol constants | PASS. Analysis freeze is Python 3.12.3; SI Table S1 now splits analysis vs docking environments. |
| S2 | frozen box JSON + `cognate_rmsd.csv` | PASS. Coverage vs top-1 not mixed |
| S3 | `label_aggregation_sensitivity.csv`; max/median `max_vs_median_sensitivity.csv` (AChE n = 94) | PASS |
| S4 | `fixed_score_negative_class_delta.csv`; cluster `cluster_bootstrap_sensitivity.csv`; detectable-effect `detectable_effect_simulation.csv` | PASS. JAK document cluster remains `unresolved_mapping_unavailable` |
| S5 | `ecfp4_incremental_information.csv`; `descriptor_baselines.csv`; scaler `ecfp4_scaler_sensitivity.csv` | PASS |
| S6 | `matched_mismatched_pocket.csv`; `holdout_metrics.csv` | PASS |
| S7 | `computational_robustness.csv`; `receptor_substitution.csv`; `five_seed_summary_min.csv` | PASS after scheme-B RTM/CNN intervals |
| S8 | `external_eligibility.csv` / frozen `external_slice_summary_v1.csv` | PASS. 0/8 is eligibility, not external validation |
| S9 | frozen `pair_eligibility_audit_s14_v1.csv` | PASS |
| S10 | `top10_operating_points.csv`; `and_filter_operating_points.csv` | PASS |

No typeset numeric cell in Table 1–3 / S2–S10 is a hand entry without a current source file.

---

## Code simplification

- Removed hashlib / `DATA_SNAPSHOT_COMMIT` / `snapshot_path` / `inputs_sha256` / git-head gates from the figure writer
- Removed `baseline_sha_expected` as a science gate (`git_head` remains a log line only)
- Removed manuscript token blacklists and magic `fig3B_max_abs > 0.05` “implausibly large” checks
- Figure builders renamed; `main()` order is Figure 1–5, S1–S5, TOC
- `promote_freeze_to_canonical.py` requires the expected CSV filename set, replaces canonical, and asserts set equality. No hashes
- `close_publication_from_canonical.py` is retired and does not rewrite historical publication CSVs

---

## Historical cleanup

Deleted from the current tree (git history retains them; historical / unavailable; not current sources):

- `data/jcim_novelty_v0/tables/incremental_information_v1.csv` (historical / unavailable)
- `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv` (historical / unavailable)
- `data/jcim_strengthen_t0t1_v0/tables/ligand_ml_baseline_scaffold_cv_v1.csv` (historical / unavailable)
- `data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv` (historical / unavailable)
- `data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv` (historical / unavailable)
- `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv` (historical / unavailable)
- `data/jcim_strengthen_t0t1_v0/tables/pocket_matched_vs_best_descriptor_delta_v1.csv` (historical / unavailable)
- `data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/descriptor_paired_delta_s19_v1.csv` (historical / unavailable)
- `data/jcim_novelty_v0/tables/ecfp4_docking_scaler_sensitivity_v1.csv` (historical / unavailable)
- `docs/HISTORICAL_LEFTOVER_FILES.md` (historical / unavailable)
- `docs/FINAL_PRE_SUBMISSION_AUDIT.md` (historical / unavailable)
- `results/freeze_rebuild/` (historical / unavailable; duplicate authoritative CSV tree)
- `figures/jcim_article/plotted_values.json` (historical / unavailable; must be absent)

Frozen experimental inputs that remain: ablation / GNINA pose scores, receptor-substitution scores, five-seed score tables, RTM/CNN pose scores, cognate RMSD, universe census, pair eligibility, external eligibility.

---

## Submission pack

- `submission_pack/tables/canonical/`: 29 CSVs, byte-identical to `results/canonical/*.csv`
- Artwork aliases only: `Figure1–5`, `FigureS1–5`, `TOC` (pdf/png/tif as produced). Source stems remain under `figures/jcim_article/`
- README points to `docs/LOCAL_RUN.md` for the full zero-dock chain (including `fit_ecfp4_models.py`)
- `scripts/qa/check_current_chain.py`: PASS

---

## Outstanding scientific limitations

These are genuine unresolved boundaries, not open calculation bugs:

1. JAK1/TYK2 document-cluster bootstrap remains `unresolved_mapping_unavailable`. No old JAK document CI is used as current.
2. External eligibility 0/8 is not external validation. No external docking was added.
3. The eight pair observations are not eight independent biological replicates: JAK1 and PPARA are reused.
4. Pointwise CIs are not multiplicity-adjusted. No eight-pair global inference was added.
5. PIK3CA/PIK3CB is not in the primary eight-pair panel.
6. EGFR/HER2 five-seed scores are the corrected-box realization restored from commit `40edc431` and recomputed with current activity-eligible labels. They are comparable to Table 2.
7. Table S1 analysis freeze is Python 3.12.3. Docking versions are written only with direct evidence.
8. Some SI “Source:” filename strings still name deleted historical CSVs. Typeset numbers come from `results/canonical`. Prose was not rewritten.

No new docking, model, pair, descriptor, threshold, or bootstrap scheme is required to close the current manuscript’s scientific numbers.
