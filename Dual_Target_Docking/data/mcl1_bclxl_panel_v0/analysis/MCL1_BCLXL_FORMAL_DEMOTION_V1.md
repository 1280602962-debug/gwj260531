# MCL1_BCLXL_FORMAL_DEMOTION_V1

Updated: 2026-08-27

## Decision

**Option B — formal demotion.** MCL1/Bcl-xL is **not** a pose-gold-validated fifth pair, **not** a non-kinase domain extension claim, and **not** external validation.

## Why Option B (not Option A)

1. **Hungarian element-matched RMSD** (prior session) reported 3WIY top1/best3 = 1.689 Å and 3WIZ best-of-top3 = 2.011 Å. That metric does **not** preserve molecular graph isomorphism and cannot certify a topology-aware pose-gold gate.
2. **RDKit `GetBestRMS` topology-aware recompute** on the same cognate outputs failed to produce usable RMSDs: Open Babel PDBQT→SDF conversion of docked LC6 did not yield a chemical graph isomorphic to the crystal LC6 reference (`natoms_mismatch` / no successful graph match for either end). Archive: `tables/cognate_qc_lc6_topology_rmsd_v1.json`.
3. Prespecified secondary gate elements (PoseBusters / physical validity, key-interaction recovery, second random seed) were **not** completed. Without topology-aware RMSD, they would not salvage a formal gate.
4. Stop rule: do **not** retune box, receptor, or exhaustiveness to rescue AUROC after seeing labels.

## Allowed uses

- Repository / SI archive entry documenting an **exploratory applicability stress-test** on a homologous BCL-2 BH3-groove pair.
- Honest statement that LC6 formal pose-gold was **not** established.

## Forbidden uses (manuscript)

- Do **not** report Hungarian or failed topology RMSD numbers as gate-pass evidence in Results.
- Do **not** report panel AUROCs (`summary_min` 0.609, Dual-vs-neither 0.628, etc.) in the main Results as confirmatory.
- Do **not** call MCL1/Bcl-xL a fifth main pair, first non-kinase pair, or PPI domain-extension success.
- Do **not** use MCL1 to support target-general claims.

## Freeze metadata

`mcl1_bclxl_panel_freeze_v1.csv`:
- `docked=1` (Vina panel scores exist in `data/mcl1_bclxl_panel_v0/`)
- `pose_gold_gate=formal_demotion_option_B; hungarian_was_preliminary_only; topology_GetBestRMS_not_computable`
- Analysis hierarchy: **post-hoc exploratory**

## Files retained (archive only)

- `data/mcl1_bclxl_panel_v0/` scores, poses, scripts
- `analysis/MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1.md` (superseded for formal claims by this demotion)
- `tables/formulation_auroc_MBX_v1.csv` (exploratory only)
