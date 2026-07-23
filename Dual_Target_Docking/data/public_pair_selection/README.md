# Public target-pair selection artifacts

Audit date: **2026-07-23**.

| File | Role |
|------|------|
| `FROZEN_PUBLIC_PAIRS.yaml` | Locked three public Dual-VSDS pairs |
| `chembl_pair_fourclass.csv` | Hard-gate four-class counts |
| `egfr_her2_fourclass_chembl_ids.csv` | EGFR/HER2 molecule-level dual / A_only / B_only / neither IDs |
| `pdb_holo_counts.csv` | RCSB holo structure rough counts |
| `chembl_target_ids.json` | UniProt → ChEMBL SINGLE PROTEIN IDs |
| `mols_*.json` | Cached molecule→best pChEMBL maps (regenerable) |

Human-readable report: [`../../docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md`](../../docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md)

Regenerate ChEMBL table:

```bash
python Dual_Target_Docking/scripts/audit_public_target_pairs.py
# or force refresh:
python Dual_Target_Docking/scripts/audit_public_target_pairs.py --refetch
```

Export EGFR/HER2 labeled molecule list (for diagnostic panel):

```bash
python Dual_Target_Docking/scripts/export_egfr_her2_fourclass.py
# → data/public_pair_selection/egfr_her2_fourclass_chembl_ids.csv
# class=dual means both ends measured and pChEMBL ≥ 6 (operational true dual)
```