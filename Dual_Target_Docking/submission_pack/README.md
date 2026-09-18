# JCIM submission pack

Current English/Chinese manuscripts, SI, figures, and canonical result tables.
Authoritative per-ligand scores: `tables/canonical/current_score_master.csv`.
Figure 3 / Table S5 ECFP incremental numbers come from `tables/canonical/ecfp4_incremental_information.csv` (freeze max |Δ|=0.0234 at PPARA/PPARD D vs B; three-decimal claim 0.023).
Regenerate from Dual_Target_Docking with `python3 scripts/analysis/compute_canonical_results.py` then `python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root Dual_Target_Docking`.
