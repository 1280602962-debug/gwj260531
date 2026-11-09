# Ligand-input reproducibility

Date: 2026-09-18  
Rule: recover scripts / PDBQT / status tables from git when they exist. Do not re-dock production panels to fill missing ligand PDBQT. Do not claim that the repository can fully rerun all original docking from SMILES unless that is true.

## Status by pack

| Pack | Prep script | Status table | Exact ligand PDBQT | Production scores | What is reproducible now |
|------|-------------|--------------|--------------------|-------------------|--------------------------|
| EGFR/HER2 historical deposited ablation | original prep not recoverable; LigPrep job log is not proof of current PDBQT | `data/jcim_multiseed_v0/tables/multiseed_ligand_prep_v1.csv` is not the production EH40 chain | **not_recoverable** (deleted; never in git as production set) | `data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv` **historical / unavailable** as Table 2 | analysis chain from deposited scores only; not current |
| EGFR/HER2 uniform RDKit/Meeko (current) | `data/egfr_her2_uniform_rdkit_v1/scripts/prep_uniform_ligands.py` **confirmed** | `data/egfr_her2_uniform_rdkit_v1/tables/ligand_prep_status.csv` **confirmed** (110/110) | `data/egfr_her2_uniform_rdkit_v1/ligands_pdbqt/` **confirmed** | `tables/scores_vina_mode1_fiveseed.csv` **confirmed** (1090 ok / 10 timeout_skipped on EH40_31); independent GNINA `tables/gnina_dock_scores_EGFR_HER2.csv` **confirmed** (185 ok / 33 timeout_skipped / 2 fail on EH120_109) | this is the only current EGFR ligand-prep and production-score truth |
| PIK3CA/mTOR PM48 | recovered `data/pik3ca_mtor_panel48_rdkit_v0/scripts/pm48_rdkit_prep_dock.py` (git `c213c485`) | not deposited as a separate status CSV | **not_recoverable** (`ligands_pdbqt/` absent from git history) | `tables/ablation_ligand_scores.csv` **confirmed** | analysis chain from deposited scores; script is a method snapshot, not a bit-identical rerun ticket |
| AChE/BChE | recovered `data/ache_bche_panel_v0/scripts/dock_panel.py` and `fetch_panel_smiles.py` (git `c213c485`) | not deposited | **not_recoverable** | `tables/ablation_ligand_scores.csv` **confirmed** | analysis chain from deposited scores |
| Track-B five pairs | recovered `data/jcim_chembl_universe_v0/scripts/prep_track_b_ligands_v1.py` | recovered `data/jcim_chembl_universe_v0/tables/track_b_ligand_prep_status_v1.csv` | **not_recoverable** (`ligands_pdbqt/` never in git) | `local_track_b_v0/tables/scores_vina_mode1_v1.csv` **confirmed**; production **poses** deposited under `local_track_b_v0/poses/` | analysis chain from deposited scores; poses exist; input ligand PDBQT do not |

## What the recovered scripts are

- RDKit ETKDGv3 seed `20260727`, largest organic fragment, MMFF ≤200, Meeko `mk_prepare_ligand.py`.
- PIK3CA script hard-codes a historical ROOT (`/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0`). That path is **historical / unavailable** as a current instruction.
- AChE `freeze_receptors_cognate_qc_v2.py` is a receptor/cognate QC snapshot that lists BChE as **6ZWI**, not production **4BDS**. It is not a current receptor-rebuild instruction.

## Claims that are false

- “The repository can fully rerun all original docking from SMILES.”
- “Every original ligand PDBQT is deposited.”
- “EGFR deposited scores were produced by the recovered Track-B RDKit/Meeko script.” (they were not; CASE 2 rebuild is a new uniform prep.)

Zero-dock statistics remain reproducible from the deposited score tables and `results/canonical/current_score_master.csv`.
