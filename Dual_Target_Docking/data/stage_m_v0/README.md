# Stage M v0 — Measurement audit pack

Zero-docking audit (M1–M3, M5) plus **M4-min** (EH40 unified RDKit re-dock) after M2=Weak.

## Relation to `plan_v2_redteam_v0`

- Reuses / rechecks redteam directional numbers (±0.005).
- Expands into formal tables, margin/noise/threshold tests, baseline gate, cleaned arm list, and a single `STAGE_M_VERDICT.md`.

## Layout

| path | content |
|------|---------|
| `analysis/STAGE_M_VERDICT.md` | **Total gate** (Track B Weak) |
| `analysis/M1_*.md` … `M5_*.md` | Per-gate writeups |
| `tables/` | All numeric CSVs |
| `scripts/run_m1_*.py` … | Reproducible zero-dock scripts |
| `scripts/m4_*.py` | M4-min prep/dock/RTM (local paths) |

## M4

- Executed **M4-min** (not skipped): old 40 ligands, RDKit ETKDG+meeko, 3POZ/3RCD E=8.
- Full pose dumps live under local `results/egfr_her2_panel40_reprep_rdkit_v0/`; score tables synced to `data/egfr_her2_panel40_reprep_rdkit_v0/`.

## Track B

See verdict: **Weak** — document planning only; do not launch mass docking.
