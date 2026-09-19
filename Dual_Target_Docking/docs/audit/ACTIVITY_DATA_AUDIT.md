# ACTIVITY DATA AUDIT

Date: 2026-09-19  
Machine table: `data/provenance/activity_provenance_matrix.csv`  
Independent check: recompute θ=6.0 from `pA`/`pB` on `current_score_master.csv`.

## Three activity kinds (must not be collapsed)

| kind | Where | What it is |
|------|--------|------------|
| `chembl_assay_adjudicated` | EGFR/HER2, PIK3CA/mTOR, most AChE/BChE | Track-A assay-level aggregate |
| `chembl37_dump_panel` | F2/F10, JAK1/JAK2, JAK1/TYK2, PPARG/PPARA, PPARA/PPARD | ChEMBL 37 dump pChEMBL on deposited panel. **Not** assay-level adjudicated |
| `panel_pchembl_no_audit_rows` | AB_001, AB_053, AB_054, AB_056, AB_097 | Panel pChEMBL without assay-audit rows. Only **AB_056** is in primary |

Track-B must not be described as assay-level adjudicated. Current manuscript/SI already distinguish the three kinds.

## θ=6.0 vs strict 6.5/5.5

Independent replay: **0** primary rows where `primary_class_theta6` ≠ four-state assignment from `pA`/`pB` at 6.0.

Table 2 n_scored matches θ=6.0 complete-case counts, not construction quotas.

Strict 6.5/5.5 appears in panel files (`label_rule=strict_6.5_5.5` on AChE) and Table 1 quotas. It is not the Table 2 class rule.

## Missing arms

| ligand | status | in primary? |
|--------|--------|-------------|
| EH120_059 | `unresolved_missing_arm` | no |
| AB_087 | `unresolved_missing_arm` | no |

Both remain in the master as ineligible rows. They do not enter directional AUROC.

## AChE AB_056

Reason it is in primary: it is the only `panel_pchembl_no_audit_rows` ligand that is activity-eligible and complete-case under θ=6.0.

`verify_freeze_rebuild.py` records an exclude-AB_056 sensitivity: directional 0.6504 / 0.6058; `summary_min` unchanged at 0.6058.  
This audit did not re-dock. The exclusion test is a label/membership sensitivity on deposited scores.

## Duplicate identity

Primary packs: no duplicate `ligand_id`, SMILES, or ChEMBL ID within a pair.

Cross-pair shared chemistry is expected (JAK1 appears in two pairs) and is a **shared-target** limitation, not a within-comparison duplicate.

Salt/stereo/canonicalization collisions across different ChEMBL IDs were not exhaustively InChI-keyed in this pass. No intra-pair collision was found on deposited IDs/SMILES.

## Max vs median

Canonical `max_vs_median_sensitivity.csv`:

- EGFR: 5 class flips; summary_min 0.3341 → 0.3393
- AChE: 1 flip (CHEMBL659); n=94 vs primary 95 because AB_056 has no adjudicated max/median pair
- PPARA/PPARD: 1 flip; summary_min unchanged 0.4463
- Other pairs: 0 flips

Max aggregation can systematically prefer the more potent reported value and therefore can inflate dual or selective calls. The sensitivity exists and does not replace Table 2.

This audit did not re-aggregate raw ChEMBL assay rows from the dump. Track-B max/median is dump-gated. That boundary is **METADATA_ONLY** for raw assay multiplicity.

## Historical vs adjudicated

Master stores `historical_pA`/`historical_pB` beside current `pA`/`pB`. Primary classes use current values.

## Verdict

Activity labels used in Table 2 are internally consistent with θ=6.0 on the master.  
Provenance is **not uniform**. That is a documented limitation, not a silent mix-up of Table 2 labels.
