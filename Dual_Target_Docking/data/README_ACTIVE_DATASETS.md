# Current datasets

Authoritative per-ligand table: `results/canonical/current_score_master.csv`.

## Primary pairs

| Pair | Production panel / scores |
|------|---------------------------|
| EGFR/HER2 | `data/egfr_her2_panel120_v0/` (corrected-box Vina in `tables/ablation_ligand_scores.csv`; boxes in `boxes/*_box_corrected.json`) |
| AChE/BChE | `data/ache_bche_panel_v0/` (corrected panel) |
| PIK3CA/mTOR | `data/pik3ca_mtor_panel48_rdkit_v0/` (PM48 primary) |
| JAK1/JAK2, JAK1/TYK2, F2/F10, PPARG/PPARA, PPARA/PPARD | `data/jcim_chembl_universe_v0/` Track B panels and scores |

## Sensitivity inputs (no new pairs)

- Holdout: `data/jcim_holdout_v0/`
- Independent GNINA: `data/jcim_independent_dock_v0/`
- Five-seed Vina: `data/jcim_multiseed_v0/`
- Receptor substitution: `data/jcim_structure_robust_v0/`
- PM110 protocol size: `data/pik3ca_mtor_panel110_rdkit_v0/`

Canonical statistics are regenerated into `results/canonical/`. Do not mix those tables with historical MASTER_RESULTS_TABLE or pre-fix membership scores.
