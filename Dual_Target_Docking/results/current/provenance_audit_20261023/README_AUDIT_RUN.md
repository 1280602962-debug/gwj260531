# Provenance audit run

Read-only. Outputs stay in this directory.

```bash
python3 Dual_Target_Docking/results/current/provenance_audit_20261023/audit_q110_chembl37.py
python3 Dual_Target_Docking/results/current/provenance_audit_20261023/audit_seven_pair_selection.py
```

`audit_q110_chembl37.py` opens the local ChEMBL 37 SQLite in read-only mode and does not rewrite `activity_label_source_audit_20261008/`.

`audit_seven_pair_selection.py` reads historical panels from git commit `a0ef2726` and replays only the AChE/BChE strict draw from `mols_ACHE.json` and `mols_BCHE.json`. It does not call the ChEMBL API and does not rerun the Track B extractor.
