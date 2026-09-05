# C5 docking handoff status

Updated 2026-09-05 after shortlist freeze (`scripts/freeze_c5_shortlist.py`).
Names: `data/campaigns/c5/04_shortlist_frozen/SHORTLIST_FROZEN.md`.

## Done

1. W1 benzbromarone@9DKA gate: **fail** (Top-1 pose_rmsd ≈ 3.59 Å; search_ok_selection_fail).
2. W4 146 jobs: **131 ok / 15 fail**; structural gate **pass** vs 40 decoys.
3. W2 IFP on 228 vs 64: **gate_pass=true** (OR=3.15) but IFP is a stricter subset of A1.
4. **Shortlist frozen**: 12 primary + 21 backup. No name shuffle by IFP.
   GSK-3008348 = control. Three cephalosporins dropped from reportable backup.

## Next (no required docking)

1. Rewrite manuscript drafts (`stale_docs` in campaign yaml) using the frozen tables.
2. Optional SI: retry 15 failed W4 decoy jobs, or remaining 29 W1 cells.
3. W5 MD only after explicit authorize (`md_authorized` still false). Slots: 2–3
   from the 12 primary names + crystal controls.

## Do not

- Reopen URAT1 docking-score ranking
- Promote/demote frozen names by W2 IFP or CNNscore
- Declare the 2.0 Å W1 gate passed via GetBestRMS
- Grid-search W2 or W4 thresholds
- Treat remaining 29 W1 jobs as required
- Start MD before flipping `md_authorized`
