# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate main figures and S4/S5/S7/S8 with `python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`; regenerate S1–S3/S9/S10 with `python3 data/jcim_bench_v0/scripts/plot_jcim_si_composites_v1.py`; regenerate S6 with `python3 data/jcim_novelty_v0/scripts/plot_detectable_effect_and_workflow_v1.py`. The same scripts are copied under `scripts/` for figure-folder organization.

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four activity states defined from experimental measurements at targets A and B. A-only and B-only serve as single-target selective controls. (B) Pocket-matched directional evaluation: Dual versus A-only is scored at target B, whereas Dual versus B-only is scored at target A. The lower of the two directional AUROCs is reported as the descriptive summary$_{\mathrm{min}}$. (C) Assembly of the eight primary target-pair rows. The initial ChEMBL audit contained 49 candidate pairs; four met the hard-negative supply gate, HDAC1/HDAC6 was excluded because metal-dependent docking was outside the protocol, and supply-limited EGFR/HER2 was retained. PIK3CA/PIK3CB was subsequently withdrawn after receptor-identity review, and five pairs from the later census were added. Bars show the smaller of the A-only and B-only pools under the strict activity rule. The dashed line marks the supply gate of 50 compounds; the upper and lower bar groups derive from the initial and later audits, respectively.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Results are shown for the eight primary target pairs at θ = 6.0. Horizontal gray rules separate the three initially evaluated pairs from the five pairs added after the census. (A) AUROC for Dual versus A-only scored at target B (circles) and Dual versus B-only scored at target A (squares). (B) Directional summary$_{\mathrm{min}}$ (circles) and Dual versus neither using the mean Vina score across the two targets (squares). Error bars are ligand-level bootstrap 95% confidence intervals. The diamond marks PIK3CA/mTOR Dual versus neither, for which the neither class contained four ligands. (C) Change in AUROC after replacing B-only with neither while holding the target-A score fixed. Error bars are 95% confidence intervals; diamonds indicate comparisons with an underpowered neither class. The dashed vertical line indicates no change.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Pocket-matched Vina rank AUROC and ECFP4 logistic-regression AUROC under Bemis–Murcko scaffold-grouped cross-validation for both directional contrasts. (B) Change in scaffold-grouped cross-validated AUROC after adding the corresponding Vina score to ECFP4. Positive values favor ECFP4 plus Vina; circles and squares denote Dual versus A-only and Dual versus B-only, respectively. (C) Topological polar surface area (TPSA) of the AChE/BChE ligands by activity class. Points represent individual ligands; horizontal and vertical black lines show the median and interquartile range, respectively.

## Figure 4. Computational realization.

(A) Comparison of Vina and independent GNINA 1.3.2 pose generation for the three evaluated target pairs. Circles and squares denote Vina and GNINA, respectively; filled blue symbols show directional summary$_{\mathrm{min}}$, and open orange symbols show Dual versus neither. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ after substituting the PIK3CA structure (4JPS or 5DXT) or the mTOR structure (4JSX) for the primary structures (PIK3CA 4L23 and mTOR 4JT6). Error bars are ligand-level bootstrap 95% confidence intervals. (C) Directional summary$_{\mathrm{min}}$ across five Vina random seeds. Horizontal segments show the range, circles the median, and diamonds the primary seed; the dashed line marks AUROC = 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Matched-pocket scoring uses target B for Dual versus A-only and target A for Dual versus B-only; the mismatched control exchanges these score channels without redocking. (A) Difference in summary$_{\mathrm{min}}$ between matched and mismatched scoring in the primary panels. (B) The same paired difference in unused-pool holdouts. Points and horizontal lines show the estimate and ligand-level bootstrap 95% confidence interval; blue intervals exclude zero and gray intervals include zero. EGFR/HER2 had no holdout. (C) Primary-panel (circles) and holdout (squares) summary$_{\mathrm{min}}$ estimates with 95% confidence intervals. Dashed vertical lines indicate zero in panels A and B and AUROC = 0.5 in panel C.

## Figure 6. Robustness checks and evidence boundary.

(A) Directional summary$_{\mathrm{min}}$ under three single activity thresholds and a strict two-threshold rule. The strict rule is a separate categorical definition rather than a continuation of θ. Daggers mark settings with fewer than 10 ligands in at least one directional class. (B) PIK3CA/mTOR summary$_{\mathrm{min}}$ in the PM48 and PM110 panels. (C) PIK3CA/mTOR PM48 summary$_{\mathrm{min}}$ at Vina exhaustiveness 16 and 8. In panels B and C, connected symbols are descriptive point-estimate comparisons and the dashed line marks AUROC = 0.5. (D) Outcome of the preregistered BindingDB external-set gate for the four historical contract pairs. None provided sufficient Dual, A-only, and B-only compounds after source, structure-identity, and chemical-similarity filtering; no external docking evaluation was performed.

## Figure S1. Protocol and panel sensitivities.

Original-set protocol grid, GNINA CNN rescoring of Vina poses, PM48 versus PM110, and exhaustiveness. This SI record still includes the later-withdrawn PIK3CA/PIK3CB row where that is what the source CSVs contain. Independent GNINA pose generation is Figure 4A.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged original-scrape sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`. Five-pair ChEMBL 37 dump supply is Figure 1C (right group).

## Figure S3. Additional paired bootstrap differences.

Paired bootstrap comparisons on the original docked set: matched versus mismatched scoring in the primary and holdout panels, Vina versus the strongest single-descriptor baseline, and ECFP4 estimates under scaffold-grouped and random folds. The later-withdrawn PIK3CA/PIK3CB row is marked with a dagger.

## Figure S4. Pocket-matched summary_min forest.

Vina CIs and the best single-descriptor reference on the eight primary rows. PIK3CA/PIK3CB is omitted.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` on the seven pairs that have a holdout. EGFR/HER2 has no holdout. Mismatched-pocket Δ CIs are Figure 5B.

## Figure S6. Detectable-effect simulation.

Probability that a ligand-level 95% confidence interval for `summary_min` excludes 0.5 under the simulated true AUROC values and the observed class sizes. This is a power-oriented diagnostic, not an estimate of observed docking performance.

## Figure S7. Post-hoc formulation and screening diagnostics.

θ = 6.0 J0 pair census (`docked_in_this_paper` = 4 includes the later-withdrawn PIK3CA/PIK3CB row), current primary n = 8, AND-like dual filter on the original three, and ligand-only full-map ECFP4 on the original three. Not docking upgrades and not a replacement for Table 2.

## Figure S8. BindingDB-native slice.

Filter cascade and remaining four-state counts after literature, structure, and ECFP4 < 0.70 on the original four-pair contract (`external_slice_summary_v1.csv`). Zero of four pairs meet the pre-frozen external gate; nothing was docked. Five census pairs were not re-opened as BindingDB external.

## Figure S9. Additional ligand-structure controls.

Additional ligand-structure controls on the original docked set: ECFP4 versus Vina directional AUROCs, Vina and prespecified single-descriptor baselines, covariate-adjusted Dual versus B-only logistic models, and potency- or size-matched subsets. The Vina-only logistic AUROC in panel C is distinct from the rank-based AUROC in the primary analysis.

## Figure S10. Matched versus mismatched point estimates.

Matched and mismatched scoring point estimates for the original primary panels, unused-pool holdouts, potency- and size-matched holdout subsets, and contact-count controls. Paired confidence intervals for the primary and holdout differences are shown in Figure 5.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs and no decorative arrows.
