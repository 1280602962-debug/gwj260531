# Active versus legacy datasets

This file is the allowlist for **current** DualFourClass / JCIM eight-pair analyses.
Historical directories remain in the tree for provenance. They must not be pulled
into primary tables by globbing.

## ACTIVE PRIMARY DATASETS (eight pairs)

| Pair | Panel / score source |
|------|----------------------|
| EGFR/HER2 | `data/egfr_her2_panel40_v0/` + `data/egfr_her2_panel120_v0/` (complete-case n=110; θ=6.0) |
| JAK1/JAK2 | `data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_JAK2_v1.csv` |
| JAK1/TYK2 | `data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_TYK2_v1.csv` |
| PIK3CA/mTOR | `data/pik3ca_mtor_panel48_v0/` (PM48 primary; PM110 is protocol sensitivity only) |
| AChE/BChE | `data/ache_bche_panel_v0/` |
| F2/F10 | `data/jcim_chembl_universe_v0/tables/track_b_panels/panel_F2_F10_v1.csv` |
| PPARG/PPARA | `data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARG_PPARA_v1.csv` |
| PPARA/PPARD | `data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARA_PPARD_v1.csv` |

Canonical ligand-level table used by the pre-submission audit:
`audit_outputs/canonical_ligand_table.csv`

Canonical numeric lock:
`data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv`

## ACTIVE SENSITIVITY DATASETS

- Unused-pool internal holdouts: `data/jcim_holdout_v0/` and Track B `holdout_panel_HO*_v1.csv`
- Five-seed Vina: `data/jcim_multiseed_v0/` and Track B `five_pair_local_channels_v1/`
- Independent GNINA: `data/jcim_independent_dock_v0/` and JAK1/TYK2 GNINA scores under Track B
- Receptor substitution (PIK3CA/mTOR): `data/jcim_structure_robust_v0/`
- PM110 protocol-size panel: `data/pik3ca_mtor_panel110_rdkit_v0/` (not a holdout)

## LEGACY / EXCLUDED DATASETS

### `data/pik3ca_pik3cb_panel_v0/`

**LEGACY — excluded from current eight-pair primary analysis.**

Contains historical receptor/target mismatch materials. Must not be included by
current analysis globbing. Do not treat PIK3CA/PIK3CB as a ninth primary pair.

### Backup CSVs (`*_backup.csv`)

Examples:

- `data/egfr_her2_panel120_v0/tables/scores_gnina_best_mode01_backup.csv`
- `data/egfr_her2_panel120_v0/tables/scores_gnina_long_mode01_backup.csv`
- analogous `*_backup.csv` under AChE and PIK3CA/mTOR packs

These are historical archives of a previous GNINA mode-01 extraction. Current
primary Vina analysis must not read `*_backup.csv`. Leave the files in place.

### Other non-primary packs

- `data/mcl1_bclxl_panel_v0/` — demoted / not in the eight-pair primary set
- `data/egfr_her2_panel40_v0/` exhaustiveness-sensitivity poses — subset QC, not a second primary panel

## Globbing rule

Do **not** use `data/*panel*/` as an input glob for primary metrics.
Current figure and Table 2 generators take **explicit pair paths**.
Historical scripts under `data/jcim_bench_v0/` that still name PIK3CA/PIK3CB or
`*_backup.csv` are archival comparators, not the JCIM primary pipeline.
