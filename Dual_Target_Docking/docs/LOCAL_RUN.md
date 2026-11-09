# Local run

This tree supports two usage levels:

1. **Zero-dock reproducibility**: rebuild canonical tables, figures, and the submission pack from committed scores.
2. **Docking / rescoring**: regenerate poses. Required only if a deposited score file is missing. Do not redock for the current paper.

## Zero-dock chain

From `Dual_Target_Docking/`:

```bash
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

Python packages: `requirements-analysis.txt`.

Authoritative per-ligand table: `results/canonical/current_score_master.csv`.  
Authoritative statistics: `results/canonical/`.  
Authoritative manuscripts: `docs/MANUSCRIPT_JCIM_EN.md`, `docs/MANUSCRIPT_JCIM_ZH.md`.  
Authoritative SI: `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`.

## What this run does not do

It does not call AutoDock Vina, GNINA, RTMScore, or Open Babel, and it does not download ChEMBL.

Docking binaries (optional): `scripts/check_docking_env.py`.
