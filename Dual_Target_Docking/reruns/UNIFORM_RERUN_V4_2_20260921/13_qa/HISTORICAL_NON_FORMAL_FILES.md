# HISTORICAL_NON_FORMAL_FILES

These files remain as provenance. They are not formal V4.2 analysis inputs.
See also `HISTORICAL_PATH_MAPPING.csv`. Historical audit documents still cite the original names; that is intentional.

## Old M1

- Original: `09_secondary_scoring/GNINA_VINA_POSE_RESCORE_MASTER.csv`
- Now: `09_secondary_scoring/HISTORICAL_NON_FORMAL_GNINA_VINA_POSE_RESCORE_MASTER.csv`
- Why not formal: CG0 / input representation failures; 61 original input-format failures; not the verified M1 table.
- Formal replacement: `GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv`

- Original: `09_secondary_scoring/GNINA_VINA_POSE_RESCORE_POSE_LONG.csv`
- Now: `HISTORICAL_NON_FORMAL_GNINA_VINA_POSE_RESCORE_POSE_LONG.csv`
- Formal replacement: `GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv`

- Directory `09_secondary_scoring/gnina_rescore/` was not moved. It is old M1 raw output. See its `README_HISTORICAL_NON_FORMAL.md`. Formal raw dir: `gnina_rescore_representation_v2/`.

## Old / intermediate M3

- Original: `10_gnina_independent_docking/GNINA_SEED42_PRODUCTION_MASTER.csv`
- Now: `HISTORICAL_NON_FORMAL_GNINA_SEED42_PRODUCTION_MASTER.csv`
- Why not formal: historical CG0-to-C / drop-G0 conversion did not restore frozen SDF ring closure for 29 ligands / 61 physical jobs. Topology repair was frozen before redocking. Old master is not official M3.
- Formal replacement: `GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv`

- Original: `GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv`
- Now: `HISTORICAL_NON_FORMAL_GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv`
- Why not formal: intermediate reparse used before topology-verified official merge.

- Original: `GNINA_SEED42_PRODUCTION_LEDGER.csv`
- Now: `HISTORICAL_NON_FORMAL_GNINA_SEED42_PRODUCTION_LEDGER.csv`
- Why not formal: ledger / intermediate production table, not official_m3.

## Old M0 subset

- Original: `13_qa/official_primary_seed42_score_master_formal4.csv`
- Now: `HISTORICAL_NON_FORMAL_official_primary_seed42_score_master_formal4.csv`
- Why not formal: 4-pair subset. Formal M0 is the 8-pair table.

## Pre-repair M3 audit

- Original: `13_qa/GNINA_INDEPENDENT_COMPLETENESS_AUDIT.md`
- Now: `GNINA_INDEPENDENT_COMPLETENESS_AUDIT_PRE_TOPOLOGY_REPAIR.md`
- Reports 1590 SUCCESS / 2 TIMEOUT. That is pre-repair. Final authority is `GNINA_INDEPENDENT_COMPLETENESS_AUDIT_FINAL.md` generated from the topology-verified master: 1586 / 6.

Formal statistics must not read any HISTORICAL_NON_FORMAL file.
