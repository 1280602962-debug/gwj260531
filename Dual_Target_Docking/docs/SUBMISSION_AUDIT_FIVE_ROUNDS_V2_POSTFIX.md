# Five-round JCIM submission audit (V2 post-fix)

Date: 2026-09-16
Branch: `cursor/methods-sentence-audit-c7cc`
Commit: `34a511a537e8af3be98843b2276455435f27e4f1`
This audit verifies post-fix canonical values. It does not replace `docs/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md` (pre-remediation).

Summary: **101 PASS**, **0 FAIL**.

| Round | Status | Finding |
|---|---|---|
| 20-class | PASS | canonical_audit: remediation_outputs/POST_FIX_AUDIT_REPORT.md 98/0/0 |
| 20-class | PASS | root_pointer: root report is a pointer |
| provenance | PASS | Figure 2A:EGFR/HER2:fixed_delta_pocketA: 0.4621 (expect 0.462) |
| provenance | PASS | Figure 2B:EGFR/HER2:AUROC_D_vs_B_pocketA: 0.3237 (expect 0.324) |
| provenance | PASS | Figure 2C:EGFR/HER2:summary_min: 0.3237 (expect 0.324) |
| provenance | PASS | Figure 2D:EGFR/HER2:top_dual: 1.0 (expect 1) |
| provenance | PASS | Figure 2D:AChE/BChE:top_dual: 5.0 (expect 5) |
| provenance | PASS | Figure S3A:EGFR 3POZ:top1_rmsd_A: 1.019 (expect 1.019) |
| provenance | PASS | Figure 5A:EGFR/HER2:gnina_summary_min: 0.2645 (expect 0.265) |
| provenance | PASS | Figure 5A:EGFR/HER2:gnina_D_vs_neither: 0.737 (expect 0.737) |
| manuscript-number | PASS | EN EGFR 0.324:  |
| manuscript-number | PASS | EN EGFR 0.786:  |
| manuscript-number | PASS | EN EGFR 0.462:  |
| manuscript-number | PASS | EN EGFR matched 0.056:  |
| manuscript-number | PASS | EN AChE 27 / 26 / 28:  |
| manuscript-number | PASS | EN AChE 0.177:  |
| manuscript-number | PASS | EN 3POZ 1.019:  |
| manuscript-number | PASS | EN 3RCD 1.947:  |
| manuscript-number | PASS | EN range 0.324 to 0.728:  |
| manuscript-number | PASS | ZH range 0.324:  |
| manuscript-number | PASS | EN no Figure 6:  |
| manuscript-number | PASS | ZH no Figure 6:  |
| manuscript-number | PASS | EN Figure 2D eight-pair:  |
| manuscript-number | PASS | lock V4 unique:  |
| manuscript-number | PASS | S2 18.816:  |
| manuscript-number | PASS | S2 no 18.680:  |
| manuscript-number | PASS | S2b 1.019:  |
| manuscript-number | PASS | S2b no current 9.505 unlabeled:  |
| EN/ZH | PASS | 0.324: EN=True ZH=True |
| EN/ZH | PASS | 0.462: EN=True ZH=True |
| EN/ZH | PASS | 0.056: EN=True ZH=True |
| EN/ZH | PASS | 0.177: EN=True ZH=True |
| EN/ZH | PASS | 1.019: EN=True ZH=True |
| EN/ZH | PASS | 1.947: EN=True ZH=True |
| EN/ZH | PASS | 27 / 26 / 28: EN=True ZH=True |
| EN/ZH | PASS | 0.357: EN=True ZH=True |
| EN/ZH | PASS | 1.778: EN=True ZH=True |
| xref | PASS | EN Figure 1:  |
| xref | PASS | EN Figure 2:  |
| xref | PASS | EN Figure 3:  |
| xref | PASS | EN Figure 4:  |
| xref | PASS | EN Figure 5:  |
| xref | PASS | EN no Fig6:  |
| xref | PASS | EN Figure S1:  |
| xref | PASS | EN Figure S2:  |
| xref | PASS | EN Figure S3:  |
| xref | PASS | EN Figure S4:  |
| xref | PASS | EN Figure S5:  |
| xref | PASS | ZH Figure 1:  |
| xref | PASS | ZH Figure 2:  |
| xref | PASS | ZH Figure 3:  |
| xref | PASS | ZH Figure 4:  |
| xref | PASS | ZH Figure 5:  |
| xref | PASS | ZH no Fig6:  |
| xref | PASS | ZH Figure S1:  |
| xref | PASS | ZH Figure S2:  |
| xref | PASS | ZH Figure S3:  |
| xref | PASS | ZH Figure S4:  |
| xref | PASS | ZH Figure S5:  |
| caption | PASS | CAPTIONS.md superseded or V4:  |
| caption | PASS | MANUSCRIPT_FIGURE_CAPTIONS.md superseded or V4:  |
| caption | PASS | FIGURE_PANEL_LOCK_V3.md superseded or V4:  |
| caption | PASS | no JAK1-only 2D as current def:  |
| caption | PASS | no 0.378 as current:  |
| caption | PASS | no 9.505 as current unlabeled:  |
| pack | PASS | EN manuscript checksum:  |
| pack | PASS | EN SI checksum:  |
| pack | PASS | Figure1.pdf:  |
| pack | PASS | Figure1.png:  |
| pack | PASS | Figure1.tif:  |
| pack | PASS | Figure2.pdf:  |
| pack | PASS | Figure2.png:  |
| pack | PASS | Figure2.tif:  |
| pack | PASS | Figure3.pdf:  |
| pack | PASS | Figure3.png:  |
| pack | PASS | Figure3.tif:  |
| pack | PASS | Figure4.pdf:  |
| pack | PASS | Figure4.png:  |
| pack | PASS | Figure4.tif:  |
| pack | PASS | Figure5.pdf:  |
| pack | PASS | Figure5.png:  |
| pack | PASS | Figure5.tif:  |
| pack | PASS | FigureS1.pdf:  |
| pack | PASS | FigureS1.png:  |
| pack | PASS | FigureS1.tif:  |
| pack | PASS | FigureS2.pdf:  |
| pack | PASS | FigureS2.png:  |
| pack | PASS | FigureS2.tif:  |
| pack | PASS | FigureS3.pdf:  |
| pack | PASS | FigureS3.png:  |
| pack | PASS | FigureS3.tif:  |
| pack | PASS | FigureS4.pdf:  |
| pack | PASS | FigureS4.png:  |
| pack | PASS | FigureS4.tif:  |
| pack | PASS | FigureS5.pdf:  |
| pack | PASS | FigureS5.png:  |
| pack | PASS | FigureS5.tif:  |
| pack | PASS | TOC.png:  |
| pack | PASS | TOC.tif:  |
| stale-token | PASS | blocking: [] |
| stale-token | PASS | labeled_or_unrelated: n=25 |

Post-fix tokens verified:

- EGFR D-vs-B-only 0.324; D-vs-neither 0.786; Δ 0.462
- EGFR matched-minus-mismatched 0.056; CI includes zero
- AChE/BChE n_scored 27/26/28; matched 0.177
- Top10 EGFR 1/5/5/0 EF=0.357; AChE 5/3/1/1 EF≈1.778
- Figure 2D eight-pair stacked composition
- cognate 3POZ top-1 1.019 Å; 3RCD 1.947 Å
- boxes 18.816 / 12.552 from canonical JSON

Target: 0 FAIL.

