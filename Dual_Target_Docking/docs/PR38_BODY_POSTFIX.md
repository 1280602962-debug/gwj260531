# PR #38 body (post-fix consolidation)

Independent audit of the dual-target docking record is complete. Canonical EGFR/HER2 docking boxes are now the cognate heavy-atom AABB + 5 Å construction with a 20 Å minimum edge (3POZ / 3RCD). EGFR/HER2 ligands were redocked under those boxes. The AChE/BChE panel was rebuilt after removing the erroneous ChEMBL-ID-prefix Murcko cap (n_scored = 27 / 26 / 28).

Core conclusions are retained: changing the experimental-state control at a fixed score channel changes apparent dual-target docking discrimination in a pair- and direction-dependent way, and dual-versus-neither ranking does not establish discrimination against single-target selectives.

The EGFR/HER2 matched-pocket claim is withdrawn. The post-fix matched-minus-mismatched interval is 0.056 [−0.044, 0.160] and includes zero. Only the AChE/BChE main-panel interval excludes zero.

Canonical 20-class invariant audit: **98 PASS / 0 WARNING / 0 FAIL** (`remediation_outputs/POST_FIX_AUDIT_REPORT.md`).

Publication-facing files now follow post-fix canonical data:
- EGFR D vs B-only (pocket A) = 0.324; D vs neither = 0.786; fixed-score ΔAUROC = 0.462
- 3POZ top-1 RMSD = 1.019 Å; 3RCD top-1 = 1.947 Å
- SI Table S2 boxes are generated from `phase1_boxes/*_box_corrected.json`
- Main text is 5 figures + 3 tables (no main Figure 6)
- Unique figure/table lock: `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`

The submission pack was regenerated (`submission_pack/`, built from `submission_pack_postfix/`). Five-round post-fix audit: `docs/SUBMISSION_AUDIT_FIVE_ROUNDS_V2_POSTFIX.md`. Final consolidation verdict: READY FOR LANGUAGE/EDITORIAL POLISH (`docs/FINAL_POSTFIX_SUBMISSION_AUDIT.md`).
