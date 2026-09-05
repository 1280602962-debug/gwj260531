# C5 docking handoff status

Updated 2026-09-05 after upload `c37f99ad` and zero-dock W4 scoring.
Analysis: `data/campaigns/c5/04_next/W4_W2_UPLOAD_ANALYSIS.md`.

## Done

1. W1 benzbromarone@9DKA gate: **fail** (Top-1 pose_rmsd ≈ 3.59 Å; search_ok_selection_fail).
2. W4 146 jobs launched: **131 ok / 15 fail** (5 decoys × 3 seeds, empty SDF).
3. W4 structural gate scored (frozen C1 thresholds): **pass** vs 40 decoys
   (seed42 structural 9/9 vs 11/40, p=8.2e-5; loose 9/9 vs 15/40).
   Clinical-acid background: loose==structural still (M2 not fixed on that set).
4. W2 IFP on 228 vs 64: **gate_pass=true** (OR=3.15, CI 1.72–7.08) but IFP is a
   stricter subset of A1 (0 extra actives; 7 A1-passers fail IFP). Do not retune.

## Next (no required docking)

1. Freeze shortlist: annotate existing clinical 9DKB/7ALV SDFs with W2 IFP;
   drop `beta_lactam_flag` from tier-2; do not invent a new ranker.
2. Rewrite manuscript drafts (stale_docs in campaign yaml).
3. Optional SI: retry 15 failed W4 decoy jobs, or remaining 29 W1 cells.
4. W5 MD only after shortlist freeze (`md_authorized` still false).

## Do not

- Reopen URAT1 docking-score ranking
- Declare the 2.0 Å W1 gate passed via GetBestRMS
- Grid-search W2 or W4 thresholds
- Treat remaining 29 W1 jobs as required
