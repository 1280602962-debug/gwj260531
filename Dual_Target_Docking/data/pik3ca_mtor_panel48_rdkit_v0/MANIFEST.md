# MANIFEST — pik3ca_mtor_panel48_rdkit_v0

## Role
JCIM K=4 development pair. Same 48 ligands as LigPrep panel48; **RDKit ETKDG + meeko** reprep for primary tables.

## Dock
- Receptors: 4L23 (PIK3CA), 4JT6 (mTOR) — frozen from prior pack
- Vina 1.2.7, **E=16**, n_modes=9, seed=20260727
- Jobs: **96/96 success**
- RTM best-of-9 complete

## Key tables (synced from results workspace)
- `tables/ablation_ligand_scores.csv` — primary RDKit scores
- `tables/prep_delta_vs_ligprep.csv` — LigPrep Δ
- `tables/directional_by_prep.csv` — directional AUROC by prep
- `analysis/PREP_DELTA.md`

## Directional snapshot (summary_min)
| prep | arm | D/A | D/B | min |
|------|-----|-----|-----|-----|
| ligprep_old | vina_mean | 0.698 | 0.597 | 0.597 |
| ligprep_old | rtm_min_z | 0.611 | 0.792 | 0.611 |
| rdkit_new | vina_mean | 0.722 | 0.671 | **0.671** |
| rdkit_new | rtm_min_z | 0.520 | 0.671 | **0.520** |

LigPrep poses remain sensitivity-only; not mixed into primary RDKit tables.
