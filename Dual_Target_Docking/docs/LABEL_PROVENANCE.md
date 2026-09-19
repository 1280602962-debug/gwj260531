# Label provenance

Primary analysis labels are the **θ = 6.0 four-state definition**. That is the only Table 2 class rule.

The strict **6.5 / 5.5** split is a **panel-construction / supply-audit** rule. It is not a second primary definition and is not a second Table 2.

## Three activity-source kinds

Recorded in `results/canonical/current_score_master.csv` as `activity_source_kind`:

| kind | Pairs / ligands | Meaning |
|------|-----------------|--------|
| `chembl_assay_adjudicated` | EGFR/HER2, PIK3CA/mTOR, most AChE/BChE | Track-A assay-level adjudication (`data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv`) |
| `chembl37_dump_panel` | Track-B five pairs (F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, PPARA/PPARD) | ChEMBL 37 dump pChEMBL on the deposited panel CSV; **not** assay-level adjudicated |
| `panel_pchembl_no_audit_rows` | five AChE ligands (`AB_001`, `AB_053`, `AB_054`, `AB_056`, `AB_097`); only **AB_056** is in the primary set | panel pChEMBL without assay-audit rows |

Do **not** claim that all eight pairs are assay-level adjudicated.

`EH120_059` and `AB_087` are `unresolved_missing_arm` and are excluded from primary.
