# M4 — Unified prep (M4-min: EH40 RDKit re-dock)

## STATUS: **Go** (M4-min completed)

- Protocol: 3POZ/3RCD, E=8, seed=20260727, n_modes=9, RTM best-of-9
- Prep: **RDKit ETKDG + meeko** (same as panel120 new70)
- Scope: EH40 only (`from_panel40=yes`); M4-full not run
- Local results: `results/egfr_her2_panel40_reprep_rdkit_v0/`; repo tables under `data/egfr_her2_panel40_reprep_rdkit_v0/`

## Directional AUROC by prep (same 40 ligands)

| prep | arm | D/A | D/B | pooled | summary_min |
|------|-----|-----|-----|--------|-------------|
| ligprep_old | vina_mean | 0.628 | 0.453 | 0.548 | 0.453 |
| ligprep_old | rtm_min_z | 0.800 | 0.607 | 0.712 | 0.607 |
| rdkit_new | vina_mean | 0.600 | 0.347 | 0.485 | 0.347 |
| rdkit_new | rtm_min_z | 0.656 | 0.467 | 0.570 | 0.467 |

## Per-ligand |Δscore| (RDKit − LigPrep)

- |Δvina_mean| median ≈ 0.10
- |Δrtm_min_z| median ≈ 0.21

## Interpretation

- Under **both** preps, EH40 `rtm_min_z` still beats `vina_mean` on `summary_min`, so the LigPrep→RDKit switch does **not** fully erase an RTM edge on this small old40 set.
- Absolute RTM is **strongly prep-sensitive**: LigPrep D/A 0.80 → RDKit 0.66; D/B 0.61 → 0.47. Mixed panel120 old40/new70 RTM splits therefore **must not** be written as confirmed method conclusions.
- Correct Track A / S1 wording: treat prep confound as a first-class caveat; prefer unified-prep numbers (this pack) when discussing RTM on EH40.

Tables: `stage_m_v0/tables/m4_old40_prep_delta.csv`, `m4_directional_by_prep.csv`.
