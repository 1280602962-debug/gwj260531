# Recovered AChE/BChE scripts

Restored from git commit `c213c485` as method provenance. They are not a license to re-dock the production panel or re-prepare 4BDS.

| file | role |
|------|------|
| `dock_panel.py` | RDKit ETKDGv3 + Meeko + Vina recipe used for the AChE panel. Exact production ligand PDBQT are **not** in git. |
| `fetch_panel_smiles.py` | SMILES fetch helper. |
| `freeze_receptors_cognate_qc_v2.py` | Historical receptor/cognate QC. Lists BChE as **6ZWI**, not production **4BDS**. Do not run it against current receptors. |

Current scores: `tables/ablation_ligand_scores.csv`. Current receptors: `ACHE_receptor.pdbqt` / `BCHE_receptor.pdbqt`.
