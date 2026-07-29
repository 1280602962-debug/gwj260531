# jcim_bench_v0 — DualFourClass-Bench pack

JCIM evaluation/benchmark aggregation after docking-phase scoring.

Authorization: [`../protocols/PAIR_ROLES_APPROVED_JCIM.yaml`](../protocols/PAIR_ROLES_APPROVED_JCIM.yaml)  
Claim ceiling: [`CLAIM_CEILING.md`](CLAIM_CEILING.md)  
GNINA: [`analysis/GNINA_STATUS.md`](analysis/GNINA_STATUS.md) (**DONE**, mode_01 CNN rescore)  
Directional forest: [`tables/directional_forest_v0.csv`](tables/directional_forest_v0.csv) · **preferred** [`tables/directional_with_baselines_v1.csv`](tables/directional_with_baselines_v1.csv) · [`analysis/DIRECTIONAL_FOREST.md`](analysis/DIRECTIONAL_FOREST.md) · [`analysis/POST_DOCKING_VERDICT.md`](analysis/POST_DOCKING_VERDICT.md)

| Pair | Pack | Vina | RTM | GNINA |
|------|------|------|-----|-------|
| EGFR/HER2 | `../egfr_her2_panel120_v0/` | existing | existing | mode_01 done |
| PIK3CA/mTOR | `../pik3ca_mtor_panel48_rdkit_v0/` | 96/96 | done + LigPrep Δ | mode_01 done |
| AChE/BChE | `../ache_bche_panel_v0/` | 191/200 | done | mode_01 done |
| PIK3CA/PIK3CB | `../pik3ca_pik3cb_panel_v0/` | 199/200 | done | mode_01 done |

Primary large pose workspaces remain under  
`/mnt/d/CADD paper exercise/dual target docking/results/` (not all poses committed).

## Reproduce main tables
```bash
# directional forest already in tables/directional_forest_v0.csv
# per-pack scores: ablation_ligand_scores.csv + scores_gnina_best.csv
```

## Status
Docking-phase **scoring complete** (Vina + RTM + GNINA). Next: manuscript / trivial baselines polish if needed for submission.
