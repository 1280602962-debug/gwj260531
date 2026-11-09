# FIGURE / TABLE LOCK POSTFIX V4

**File:** `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`

**Status:** CURRENT publication-facing figure/table numbering and source provenance.  
**Supersedes (historical / unavailable):** `docs/FIGURE_PANEL_LOCK_V3.md`, `docs/FIGURE_TABLE_PLACEMENT_REPORT_ZH.md`. Captions live in `figures/jcim_article/CAPTIONS.md` and `figures/jcim_article/MANUSCRIPT_FIGURE_CAPTIONS.md`.  
**Plotter:** `figures/jcim_article/scripts/update_figures_pr32.py` (single official pipeline).  
**Value lock:** `figures/jcim_article/plotted_values_postfix.json`.  
**Canonical results:** `results/canonical/`.  
Historical `*_v1.csv` names under `data/jcim_*` are **superseded aliases**, not current Unique sources.

All plotted AUROC, CI, n, EF, Top-10 counts, and RMSD values are read from canonical post-fix CSV/JSON. Hard-coded values are limited to axis limits, fonts, panel geometry, and reference lines at 0, 0.5, and 2 Å.

Main text: **5 Figures + 3 Tables**. SI: **Figures S1–S5** and **Tables S1–S10**. There is no main-text Figure 6.

Narrative order: directional ranking → candidate-ranking consequence → ligand chemistry → pocket correspondence → robustness.

## Main figures

| Figure | Panel | Content | Unique source |
|---|---|---|---|
| 1 | A | Four experimental states: dual / A-only / B-only / neither | schematic |
| 1 | B | Directional evaluation: D vs A-only → score B; D vs B-only → score A | schematic |
| 1 | C | Paired-data supply census 2,164,618 / 63,790 / 5,253 / 86; separate box for Primary evaluation: 8 target pairs. Not a direct 86→8 funnel | `universe_census_summary_v1.csv` |
| 2 | A | Fixed-score ΔAUROC = AUROC(D vs neither) − AUROC(D vs matched single-target-active); 8 pairs × 2 directions; 95% CI. EGFR pocket A Δ = 0.462 | `results/canonical/fixed_score_negative_class_delta.csv` |
| 2 | B | Two directional AUROCs: D vs A-only using B; D vs B-only using A; 8 pairs | `results/canonical/primary_directional_auroc.csv` |
| 2 | C | Directional summary_min vs two-pocket mean D-vs-neither. Not a fixed-score comparison. PIK3CA/mTOR neither n=4 annotation | `results/canonical/primary_summary_min.csv`; `results/canonical/two_pocket_mean_ranking.csv`; `results/canonical/class_counts.csv` |
| 2 | D | Eight-pair top-10% 100% stacked bars (dual / A-only / B-only / neither), vina_mean rank, k=ceil(0.10 n) at right. No EF on the figure | `results/canonical/top10_operating_points.csv` |
| 3 | A | Vina directional rank AUROC vs ECFP4 scaffold-GroupKFold OOF AUROC; 8 pairs / 16 directions. Not an algorithm leaderboard | `results/canonical/ecfp4_incremental_information.csv` |
| 3 | B | ΔAUROC = ECFP4+docking − ECFP4; same folds, same ligands. Freeze max \|Δ\| = 0.0234 (PPARA/PPARD D vs B), written as 0.023 | `results/canonical/ecfp4_incremental_information.csv` |
| 4 | A | Main-panel matched-minus-mismatched Δsummary_min; 8 pairs; 95% CI. Color = main, not significance | `results/canonical/matched_mismatched_pocket.csv` |
| 4 | B | Unused-pool internal holdout matched-minus-mismatched Δsummary_min; 7 pairs; 95% CI. EGFR/HER2 has no holdout. Color = holdout | `results/canonical/holdout_metrics.csv` |
| 5 | A | Independent GNINA pose generation: EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2; summary_min and D-vs-neither. Corrected EGFR box: D-vs-B-only ≈ 0.265; D-vs-neither ≈ 0.737 | `results/canonical/computational_robustness.csv` |
| 5 | B | PIK3CA/mTOR receptor substitution: primary 4L23 / 4JT6, 4JPS, 5DXT, 4JSX | alt-receptor CSVs + Table 2 |
| 5 | C | Five-seed summary_min range; 8 pairs; EGFR/HER2 uses corrected-box scores. AChE available-case membership changes; fixed-membership n=88 is a sensitivity | `results/canonical/five_seed_summary_min.csv`; `results/canonical/five_seed_fixed_membership_sensitivity.csv` |

AChE/BChE TPSA jitter is **Figure S1B**, not a main Figure 3 panel.  
Main-versus-holdout summary_min is **Table S6**, not a third Figure 4 panel.  
RTMScore / same-pose CNN rescoring is **Table S7**, not a main panel.  
Threshold, cluster bootstrap, BindingDB eligibility, and PIK3CA protocol-size/exhaustiveness are SI only.

Figure 4 caption requirement: AChE/BChE and EGFR/HER2 main-panel intervals exclude zero; all seven holdout intervals include zero.

## Main tables

| Table | Content | Source |
|---|---|---|
| 1 | Pair, candidate pool, quota D/A/B/N, PDB A/B, n_scored D/A/B, Vina exhaustiveness. AChE/BChE n_scored = 27 / 26 / 28 | `results/canonical/class_counts.csv`; panel CSVs; `scripts/analysis/analysis_config.py` `RECEPTORS` |
| 2 | Pair, n_scored, D vs A-only pocket B [95% CI], D vs B-only pocket A [95% CI], summary_min [95% CI] | `results/canonical/primary_summary_min.csv`; `results/canonical/primary_directional_auroc.csv`; `results/canonical/class_counts.csv` |
| 3 | Pair, two-pocket mean D-vs-neither AUROC [95% CI], n_neither, Top 10% D/A/B/N, EFdual,10%. No summary_min column. PIK3CA/mTOR n_neither=4 footnote | `results/canonical/two_pocket_mean_ranking.csv`; `results/canonical/top10_operating_points.csv` |

Figure 2C is jointly readable from Table 2 + Table 3. Mean-score D-vs-pooled(A+B) is Table S10 secondary diagnostic.

## SI figures

| Figure | Panel | Content | Source |
|---|---|---|---|
| S1 | A | Best single descriptor vs Vina summary_min forest | `results/canonical/descriptor_baselines.csv` + Table 2 |
| S1 | B | AChE/BChE TPSA distribution | `results/canonical/current_score_master.csv` |
| S2 | A | PM48 vs PM110 | `results/canonical/protocol_sensitivity.csv` |
| S2 | B | E=16 vs E=8 | `results/canonical/protocol_sensitivity.csv` |
| S3 | — | Cognate redocking RMSD, 14 primary slots; top-1 and lowest saved pose; 2 Å line. Post-fix EGFR 3POZ top-1 = 1.019 Å; HER2 3RCD top-1 = 1.947 Å | `results/canonical/cognate_rmsd.csv` |
| S4 | A | Theta threshold sensitivity | `results/canonical/label_aggregation_sensitivity.csv` |
| S4 | B | Cluster-resampling intervals, EGFR/HER2 ligand-level official + JAK1/TYK2 cluster. EGFR cluster recomputed with corrected-box scores and frozen groupings | `results/canonical/cluster_bootstrap_sensitivity.csv` |
| S5 | A | BindingDB independent-filter compound counts | `results/canonical/external_eligibility.csv` |
| S5 | B | Independent-source counts; 0/8 passed full external docking gate. Eligibility screen, not external validation | `results/canonical/external_eligibility.csv` |

Not typeset as SI figures: JAK1/TYK2 AND-filter bar chart; duplicate holdout summary; duplicate matched/mismatched bars; duplicate BindingDB matrices; old PIK3CA/PIK3CB figures; historical J0 plots.

## SI tables

| Table | Content |
|---|---|
| S1 | Computational settings and statistical definitions |
| S2 | Receptors, canonical docking boxes from `phase1_boxes/*_box_corrected.json`, cognate-redocking QC |
| S3 | Activity-threshold and max-vs-median aggregation sensitivity |
| S4 | Fixed-score control-class comparisons + cluster bootstrap |
| S5 | Ligand chemistry: single descriptors, ECFP4, ECFP4+docking, incremental AUROC, feature-scaling sensitivity |
| S6 | Matched vs mismatched pocket, main + unused-pool holdout |
| S7 | Computational realization: independent GNINA, receptor substitution, RTMScore / CNN rescoring, five-seed summary |
| S8 | External-data eligibility after independence filters |
| S9 | Final target-pair eligibility audit |
| S10 | Eight-pair candidate-ranking operating points: n_ranked D/A/B/N, k, Top10% D/A/B/N, dual fraction, panel dual base rate, EFdual,10%, AND filter; optional mean-score D-vs-pooled A+B as secondary diagnostic |

Do not create S11/S12 unless a genuinely new submission-required analysis appears.

## Canonical post-fix flagship values

EGFR/HER2: D vs A-only pocket B = 0.663; D vs B-only pocket A = 0.334; summary_min = 0.334 [0.197, 0.471].  
Same pocket-A score: D vs B-only = 0.334; D vs neither = 0.780; fixed-score ΔAUROC = 0.446 [0.259, 0.632].  
Two-pocket mean D-vs-neither = 0.759 [0.551, 0.926].  
Matched-minus-mismatched = 0.107 [0.006, 0.220]; CI excludes zero.  
AChE/BChE: n_scored = 27 / 26 / 28; matched-minus-mismatched = 0.177 [0.053, 0.291]; CI excludes 0. EGFR/HER2 and AChE/BChE are the two main-panel intervals that exclude zero.  
EGFR/HER2 Top 10% = 1 / 4 / 6 / 0; EFdual,10% = 0.351.  
AChE/BChE Top 10% = 5 / 3 / 1 / 1; EFdual,10% = 1.759.  
GNINA EGFR: D-vs-B-only = 0.227 [0.104, 0.373]; D-vs-neither = 0.705; n=89.  
3POZ box center = 18.816, 31.837, 11.725; size = 20.931, 20.000, 21.342 (from JSON).  
3RCD box center = 12.552, 2.982, 28.152; size = 21.385, 22.378, 20.000 (from JSON).

Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are not current publication results.  
Current canonical ECFP max \|Δ\| is taken from `results/canonical/ecfp4_incremental_information.csv` (0.0234 at PPARA/PPARD D vs B). The old 0.0112 value is a superseded historical-environment result; background is the superseded snapshot `docs/archive/PR39_SCIENTIFIC_DATA_AUDIT.md` (not current truth). Current authority: `results/canonical/ecfp4_incremental_information.csv`.

## Results mapping

| Section | Figures / tables |
|---|---|
| 3.1 Paired experimental data and four-state evaluation panels | Figure 1 |
| 3.2 Directional docking performance and candidate-ranking consequence | Figure 2; Tables 2 and 3 |
| 3.3 Ligand-chemistry baselines and incremental docking discrimination | Figure 3 |
| 3.4 Pocket correspondence | Figure 4 |
| 3.5 Computational robustness | Figure 5 |
| 3.6 Sensitivity and external-data scope | SI only (S2–S5, S3, S4, S8) |
