# Acid track operations (Amendment A2)

## Status (2026-08-27)

| Step | Status |
|------|--------|
| A1 freeze (24 exploratory) | DONE — `acid_dual_a1_frozen/` |
| Amendment A2 preregistration | DONE — `00_preregistration/AMENDMENT_A2_*.yaml` |
| A2 reference validation | DONE — 5/6 carboxylate refs pass |
| Retrospective acid-gate benchmark | DONE — A1 OR≈3.2; A2 OR≈1.0 (not discriminative) |
| A2 clinical seed 42 rescore | DONE — 59 dual keep |
| Seeds 43/44 docking | IN PROGRESS — `acid_dual_a2_multiseed.log` |
| Competition shortlist | DONE (seed42) — `acid_shortlist_a2_competition.csv` |
| MD (L7) | CLOSED — `md_authorized=false` |

## Pose selection

- **A1 (retired for nomination):** CNNscore Top-1 → geometry check
- **A2 (active):** geometry-compatible poses → CNNscore tie-break

## Do not

- Overwrite `acid_dual_a1_frozen/`
- Use A2 gate as activity retrieval (benchmark shows OR≈1)
- Open MD before multiseed stability + explicit authorization
