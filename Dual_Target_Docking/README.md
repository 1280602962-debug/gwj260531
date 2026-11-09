# Dual-target docking evaluation (eight primary pairs)

This directory is the current article, analysis, and reproduction tree for an evaluation of dual-target docking. It is not a new scoring function, a general virtual-screening leaderboard, or a wet-lab study.

**Question.** When ligands have experimental measurements at both targets, does a favorable dual-versus-neither docking readout also separate dual-active ligands from single-target selectives?

**Primary pairs.** EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, PPARA/PPARD.

**Primary labels and tasks.** θ = 6.0 four-state labels (dual / A-only / B-only / neither). Dual vs A-only uses pocket B; dual vs B-only uses pocket A. `summary_min` is the weaker of those two AUROCs (a descriptive summary, not a new ligand score). Two-pocket mean scores are used only for dual-versus-neither and four-state ranking. EGFR/HER2 uses corrected cognate-heavy-atom boxes; AChE/BChE uses the corrected panel.

## Where the current numbers live

| What | Path |
|------|------|
| Per-ligand scores and classes | `results/canonical/current_score_master.csv` (pointer: `data/processed/CURRENT_SCORE_MASTER.md`) |
| Canonical result tables | `results/canonical/` |
| Detectable-effect simulation | `results/canonical/detectable_effect_simulation.csv` |
| English / Chinese manuscripts | `docs/MANUSCRIPT_JCIM_EN.md`, `docs/MANUSCRIPT_JCIM_ZH.md` |
| SI | `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`, `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` |
| Figures 1–5 and S1–S5 | `figures/jcim_article/` |
| Submission slice | `submission_pack/` |
| Common statistics | `scripts/analysis/bootstrap_metrics.py` |

Primary uncertainty is a **class-stratified nonparametric percentile bootstrap** (B = 2000, seed 20260729). Dual ligands are drawn once per replicate and reused across both directional AUROCs, `summary_min`, and fixed-score Δ. Non-stratified ligand bootstrap and scaffold/document cluster bootstrap are sensitivity analyses only.

## Zero-dock recomputation (no Vina/GNINA)

```bash
cd Dual_Target_Docking
python3 scripts/check_analysis_env.py
python3 scripts/analysis/adjudicate_activity_records.py
python3 scripts/analysis/build_current_score_master.py
python3 scripts/analysis/compute_canonical_results.py
python3 scripts/analysis/compute_leave_one_document.py
python3 scripts/analysis/compute_class_chemistry.py
python3 scripts/analysis/compute_detectable_effect.py
python3 scripts/analysis/fit_ecfp4_models.py
python3 scripts/analysis/compute_descriptor_baselines.py
python3 scripts/analysis/close_publication_from_canonical.py
python3 scripts/analysis/patch_publication_text.py
python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root .
python3 scripts/analysis/rebuild_submission_pack.py
python3 scripts/qa/check_current_chain.py
```

`compute_canonical_results.py` rebuilds Table 2 / Table 3 / Figure 2 sources (`primary_directional_auroc.csv`, `primary_summary_min.csv`, `fixed_score_negative_class_delta.csv`, `two_pocket_mean_ranking.csv`, `top10_operating_points.csv`). Figure 3 uses `ecfp4_oof_predictions.csv` / `ecfp4_incremental_information.csv`. Figure 4 uses `matched_mismatched_pocket.csv`. Figure 5 uses `computational_robustness.csv`, receptor-substitution CSVs, `five_seed_summary_min.csv`, and `five_seed_fixed_membership_sensitivity.csv`.

Python dependencies: `requirements-analysis.txt`.

## What requires docking

Docking is required only to (re)generate poses or scores: production Vina, independent GNINA pose generation, PIK3CA/mTOR receptor substitution, five-seed Vina, PPARG same-pose RTMScore/CNN rescoring, and cognate redocking RMSD. Those deposited scores already exist. Do not redock unless a raw score file is missing.

GNINA is another pose-generation realization, not a claim that GNINA is better or worse than Vina.

## Formal sensitivity analyses in the article

- Label rules: θ = 5.5 / 6.0 / 6.5 and strict 6.5/5.5; max vs median pChEMBL (ChEMBL 37)
- Unused-pool internal holdout
- Independent GNINA on EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2
- PIK3CA/mTOR receptor substitution (4JPS, 5DXT, 4JSX)
- Five-seed Vina
- PPARG/PPARA same-pose RTMScore / GNINA CNN
- Cognate redocking RMSD (search-coverage gate: lowest saved RMSD < 2 Å)
- Non-stratified bootstrap and scaffold/document cluster bootstrap
- Detectable-effect simulation (binormal; current eight-pair class sizes; same class-stratified shared-dual bootstrap as Table 2)

JNK1_Selectivity_Project is a separate tree and is not part of this submission.
