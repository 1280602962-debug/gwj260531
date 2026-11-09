# Backward provenance audit

Date: 2026-09-18  
Direction: canonical result → score source → frozen docking table → ligand input → receptor → box → protocol.

Machine companion: `docs/SCIENTIFIC_PROVENANCE_MATRIX.csv` and `results/qa/SCIENTIFIC_PROVENANCE_MATRIX.csv`.

## Chain that holds (zero-dock statistics)

For every primary pair:

`results/canonical/primary_summary_min.csv`
→ `results/canonical/current_score_master.csv`
→ pair score file in `analysis_config.py` (`EGFR_SCORE_SOURCE`, `ACHE_SCORE_SOURCE`, `PIK3CA_SCORE_SOURCE`, `TRACK_B_SCORE_SOURCE`)
→ receptor PDBQT on disk
→ box JSON on disk
→ `scripts/analysis/compute_canonical_results.py` + `bootstrap_metrics.py`

Activity labels:

Track-A (EGFR, PIK3CA, AChE) → `data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv` (`chembl_assay_adjudicated`, plus AChE `panel_pchembl_no_audit_rows` for five ligands).  
Track-B → panel CSV `chembl37_dump_panel`.  
`unresolved_missing_arm` (EH120_059, AB_087) are not in the primary gate.

## Break points

| Pair | Broken link | Evidence | Status |
|------|-------------|----------|--------|
| EGFR/HER2 | ligand 3D input for deposited Vina/GNINA | No `ligands_pdbqt/` under panel120; LigPrep maegz never in git; protocol is not proof | **not_recoverable** (P0) |
| EGFR/HER2 | ligand-prep method for deposited scores | CASE 2 audit; mixed LigPrep/RDKit yaml was metadata | **not_recoverable** (P0) |
| PIK3CA/mTOR | ligand PDBQT | only receptor PDBQT in pack | scores confirmed; cannot redock from current tree (P1) |
| AChE/BChE | ligand PDBQT | only receptor/cognate PDBQT | same (P1) |
| Track-B | input `ligands_pdbqt/` | absent; production `poses/` present; yaml + historical prep script | scores confirmed (P1 for missing input PDBQT) |
| All 14 receptors | protonation / loops / PPW | PDBQT exist; no prep script | file confirmed; method **not_recoverable** (P1) |
| Table S9 footnotes | `TIER1_DOCKING_ROSTER_V1.md`, `FEASIBLE_PAIR_LADDER_V1.md` | cited, missing | historical citation (P1) |
| SI scaler filename | `ecfp4_docking_scaler_sensitivity_v1.csv` | missing; canonical `ecfp4_scaler_sensitivity.csv` exists | stale name (fixed in lock this round) |

## Software-name claims vs files

| Claim | Direct evidence in current tree? |
|-------|----------------------------------|
| RDKit ETKDGv3 + Meeko for Track-B | yaml `no_ligprep: true`; historical `prep_track_b_ligands_v1.py` in git; uniform EGFR rebuild script now present | metadata + git script; input PDBQT absent |
| RDKit for PIK3CA/AChE primary | protocol.yaml `rdkit_etkdg_meeko`; primary score tables | metadata_only for 3D inputs |
| LigPrep as current EGFR or PIK3CA primary | **No.** Local job logs prove jobs ran; they do not prove current PDBQT | must not be stated as fact |
| Schrödinger Protein Preparation Wizard | **No** production script/log | not_recoverable |
| OpenBabel | SI notes CCD fallback for 2JKH/BI7 cognate QC | not a ligand-prep pipeline for panels |
| GNINA independent | deposited independent-search CSVs for EGFR, PIK3CA, JAK1/TYK2 | files exist; EGFR ligand inputs not_recoverable |
| RTMScore | `scores_rtm_best9_v1.csv` | Vina-pose rescore, not independent docking |

## Canonical figure/table cells

See `results/qa/figure_table_source_audit.csv`. Official plotter `update_figures_pr32.py` reads `results/canonical/`. Lock Unique sources were remapped to those files this round.
