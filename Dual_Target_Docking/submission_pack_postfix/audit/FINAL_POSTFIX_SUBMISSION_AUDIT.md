# FINAL POSTFIX SUBMISSION AUDIT

Commit: `34a511a537e8af3be98843b2276455435f27e4f1`
Branch: `cursor/methods-sentence-audit-c7cc`
Canonical audit: `remediation_outputs/POST_FIX_AUDIT_REPORT.md` (98 PASS / 0 WARNING / 0 FAIL)
Figure lock: `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`

## Audits

1. 20-class invariant audit — canonical report 98/0/0; root file is a pointer
2. figure-value provenance audit — plotted_values_postfix.json vs post-fix CSVs
3. manuscript-number trace — EN/ZH/SI vs canonical post-fix values
4. EN/ZH consistency audit — flagship tokens present in both languages
5. SI/main cross-reference audit — Figures 1–5 and S1–S5; no main Figure 6
6. caption/panel-number audit — V4 lock is the unique current lock
7. submission-pack checksum audit — pack copies match docs/figures

## Result table

| Class | Name | Status | Detail |
|---|---|---|---|
| 20-class | canonical_audit | PASS | remediation_outputs/POST_FIX_AUDIT_REPORT.md 98/0/0 |
| 20-class | root_pointer | PASS | root report is a pointer |
| provenance | Figure 2A:EGFR/HER2:fixed_delta_pocketA | PASS | 0.4621 (expect 0.462) |
| provenance | Figure 2B:EGFR/HER2:AUROC_D_vs_B_pocketA | PASS | 0.3237 (expect 0.324) |
| provenance | Figure 2C:EGFR/HER2:summary_min | PASS | 0.3237 (expect 0.324) |
| provenance | Figure 2D:EGFR/HER2:top_dual | PASS | 1.0 (expect 1) |
| provenance | Figure 2D:AChE/BChE:top_dual | PASS | 5.0 (expect 5) |
| provenance | Figure S3A:EGFR 3POZ:top1_rmsd_A | PASS | 1.019 (expect 1.019) |
| provenance | Figure 5A:EGFR/HER2:gnina_summary_min | PASS | 0.2645 (expect 0.265) |
| provenance | Figure 5A:EGFR/HER2:gnina_D_vs_neither | PASS | 0.737 (expect 0.737) |
| manuscript-number | EN EGFR 0.324 | PASS |  |
| manuscript-number | EN EGFR 0.786 | PASS |  |
| manuscript-number | EN EGFR 0.462 | PASS |  |
| manuscript-number | EN EGFR matched 0.056 | PASS |  |
| manuscript-number | EN AChE 27 / 26 / 28 | PASS |  |
| manuscript-number | EN AChE 0.177 | PASS |  |
| manuscript-number | EN 3POZ 1.019 | PASS |  |
| manuscript-number | EN 3RCD 1.947 | PASS |  |
| manuscript-number | EN range 0.324 to 0.728 | PASS |  |
| manuscript-number | ZH range 0.324 | PASS |  |
| manuscript-number | EN no Figure 6 | PASS |  |
| manuscript-number | ZH no Figure 6 | PASS |  |
| manuscript-number | EN Figure 2D eight-pair | PASS |  |
| manuscript-number | lock V4 unique | PASS |  |
| manuscript-number | S2 18.816 | PASS |  |
| manuscript-number | S2 no 18.680 | PASS |  |
| manuscript-number | S2b 1.019 | PASS |  |
| manuscript-number | S2b no current 9.505 unlabeled | PASS |  |
| EN/ZH | 0.324 | PASS | EN=True ZH=True |
| EN/ZH | 0.462 | PASS | EN=True ZH=True |
| EN/ZH | 0.056 | PASS | EN=True ZH=True |
| EN/ZH | 0.177 | PASS | EN=True ZH=True |
| EN/ZH | 1.019 | PASS | EN=True ZH=True |
| EN/ZH | 1.947 | PASS | EN=True ZH=True |
| EN/ZH | 27 / 26 / 28 | PASS | EN=True ZH=True |
| EN/ZH | 0.357 | PASS | EN=True ZH=True |
| EN/ZH | 1.778 | PASS | EN=True ZH=True |
| xref | EN Figure 1 | PASS |  |
| xref | EN Figure 2 | PASS |  |
| xref | EN Figure 3 | PASS |  |
| xref | EN Figure 4 | PASS |  |
| xref | EN Figure 5 | PASS |  |
| xref | EN no Fig6 | PASS |  |
| xref | EN Figure S1 | PASS |  |
| xref | EN Figure S2 | PASS |  |
| xref | EN Figure S3 | PASS |  |
| xref | EN Figure S4 | PASS |  |
| xref | EN Figure S5 | PASS |  |
| xref | ZH Figure 1 | PASS |  |
| xref | ZH Figure 2 | PASS |  |
| xref | ZH Figure 3 | PASS |  |
| xref | ZH Figure 4 | PASS |  |
| xref | ZH Figure 5 | PASS |  |
| xref | ZH no Fig6 | PASS |  |
| xref | ZH Figure S1 | PASS |  |
| xref | ZH Figure S2 | PASS |  |
| xref | ZH Figure S3 | PASS |  |
| xref | ZH Figure S4 | PASS |  |
| xref | ZH Figure S5 | PASS |  |
| caption | CAPTIONS.md superseded or V4 | PASS |  |
| caption | MANUSCRIPT_FIGURE_CAPTIONS.md superseded or V4 | PASS |  |
| caption | FIGURE_PANEL_LOCK_V3.md superseded or V4 | PASS |  |
| caption | no JAK1-only 2D as current def | PASS |  |
| caption | no 0.378 as current | PASS |  |
| caption | no 9.505 as current unlabeled | PASS |  |
| pack | EN manuscript checksum | PASS |  |
| pack | EN SI checksum | PASS |  |
| pack | Figure1.pdf | PASS |  |
| pack | Figure1.png | PASS |  |
| pack | Figure1.tif | PASS |  |
| pack | Figure2.pdf | PASS |  |
| pack | Figure2.png | PASS |  |
| pack | Figure2.tif | PASS |  |
| pack | Figure3.pdf | PASS |  |
| pack | Figure3.png | PASS |  |
| pack | Figure3.tif | PASS |  |
| pack | Figure4.pdf | PASS |  |
| pack | Figure4.png | PASS |  |
| pack | Figure4.tif | PASS |  |
| pack | Figure5.pdf | PASS |  |
| pack | Figure5.png | PASS |  |
| pack | Figure5.tif | PASS |  |
| pack | FigureS1.pdf | PASS |  |
| pack | FigureS1.png | PASS |  |
| pack | FigureS1.tif | PASS |  |
| pack | FigureS2.pdf | PASS |  |
| pack | FigureS2.png | PASS |  |
| pack | FigureS2.tif | PASS |  |
| pack | FigureS3.pdf | PASS |  |
| pack | FigureS3.png | PASS |  |
| pack | FigureS3.tif | PASS |  |
| pack | FigureS4.pdf | PASS |  |
| pack | FigureS4.png | PASS |  |
| pack | FigureS4.tif | PASS |  |
| pack | FigureS5.pdf | PASS |  |
| pack | FigureS5.png | PASS |  |
| pack | FigureS5.tif | PASS |  |
| pack | TOC.png | PASS |  |
| pack | TOC.tif | PASS |  |
| stale-token | blocking | PASS | [] |
| stale-token | labeled_or_unrelated | PASS | n=25 |

## Non-blocking stale-token notes (current chemistry / labeled historical)

- `0.808` in `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: BChE | D vs A | 0.895 | 0.893 | −0.002 | 0.652 | | AChE/BChE | D vs B | 0.821 | 0.808 | −0.013 | 0.606 | | PIK3CA/mTOR | D vs A | 0.762 | 0.742 | −0.020 | 0.714
- `0.170` in `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: 0.060, 0.295] | no | | JAK1/JAK2 | pocket A | 0.728 | 0.723 | −0.004 | [−0.177, 0.170] | no | | JAK1/JAK2 | pocket B | 0.588 | 0.730 | 0.142 | [−0.058, 0.334] |
- `0.170` in `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: OR | heavy | 0.463 | 0.229 | [−0.011, 0.435] | no | | F2/F10 | cLogP | 0.515 | −0.170 | [−0.324, −0.012] | yes | | JAK1/TYK2 | cLogP | 0.580 | −0.215 | [−0.374,
- `9.505` in `docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: lots use the unified chemically mapped CalcRMS table (Table S2). The historical 9.505 Å value is not a current result.  ### Figure S4. Label and source robustne
- `0.808` in `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`: BChE | D vs A | 0.895 | 0.893 | −0.002 | 0.652 | | AChE/BChE | D vs B | 0.821 | 0.808 | −0.013 | 0.606 | | PIK3CA/mTOR | D vs A | 0.762 | 0.742 | −0.020 | 0.714
- `0.170` in `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`:  | [−0.060, 0.295] | 否 | | JAK1/JAK2 | 口袋 A | 0.728 | 0.723 | −0.004 | [−0.177, 0.170] | 否 | | JAK1/JAK2 | 口袋 B | 0.588 | 0.730 | 0.142 | [−0.058, 0.334] | 否 | 
- `0.170` in `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`: mTOR | 重原子数 | 0.463 | 0.229 | [−0.011, 0.435] | 否 | | F2/F10 | cLogP | 0.515 | −0.170 | [−0.324, −0.012] | 是 | | JAK1/TYK2 | cLogP | 0.580 | −0.215 | [−0.374, −
- `9.505` in `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`:  Å，HER2 3RCD 的 top-1 为 1.947 Å，均低于 2 Å。14 个受体统一使用化学对应 CalcRMS 表（Table S2）。不得将历史 9.505 Å 作为当前结果。  ### Figure S4. 标签与来源稳健性  ![Figure S4](../figures/jcim_article/F
- `0.430` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: .552, 2.982, 28.152; size = 21.385, 22.378, 20.000 (from JSON).  Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are n
- `0.808` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: 982, 28.152; size = 21.385, 22.378, 20.000 (from JSON).  Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are not curre
- `0.378` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: 152; size = 21.385, 22.378, 20.000 (from JSON).  Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are not current publi
- `0.170` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: e = 21.385, 22.378, 20.000 (from JSON).  Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are not current publication r
- `0.161` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: 85, 22.378, 20.000 (from JSON).  Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are not current publication results. 
- `9.505` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`: 78, 20.000 (from JSON).  Pre-fix tokens 0.430 / 0.808 / 0.378 / 0.170 / 0.161 / 9.505 and old EGFR box coordinates are not current publication results.  ## Resu
- `Figure 6` in `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`:  3 Tables**. SI: **Figures S1–S5** and **Tables S1–S10**. There is no main-text Figure 6.  Narrative order: directional ranking → candidate-ranking consequence 
- `0.378` in `figures/jcim_article/CAPTIONS.md`: LE_LOCK_POSTFIX_V4.md`  The captions below match V4. Historical pre-fix tokens (0.378, 9.505 Å, JAK1-only Figure 2D, main-text Figure 6) are not current results
- `9.505` in `figures/jcim_article/CAPTIONS.md`: _POSTFIX_V4.md`  The captions below match V4. Historical pre-fix tokens (0.378, 9.505 Å, JAK1-only Figure 2D, main-text Figure 6) are not current results.  Rege
- `9.505` in `figures/jcim_article/CAPTIONS.md`: . Post-fix EGFR 3POZ top-1 = 1.019 Å; HER2 3RCD top-1 = 1.947 Å. The historical 9.505 Å value is not current.  ## Figure S4. Label and source robustness.  (A) T
- `Figure 6` in `figures/jcim_article/CAPTIONS.md`: h V4. Historical pre-fix tokens (0.378, 9.505 Å, JAK1-only Figure 2D, main-text Figure 6) are not current results.  Regenerate: `python3 figures/jcim_article/sc
- `9.505` in `figures/jcim_article/MANUSCRIPT_FIGURE_CAPTIONS.md`: atom RMSD. EGFR 3POZ top-1 = 1.019 Å; HER2 3RCD top-1 = 1.947 Å. The historical 9.505 Å value is not current.  ## Figure S4. Label and source robustness.  (A) A
- `Figure 6` in `figures/jcim_article/MANUSCRIPT_FIGURE_CAPTIONS.md`: otein-system order throughout Figures 2–5 and Tables 1–3. There is no main-text Figure 6.  ## Figure 1. Four-state dual-target evaluation and data supply.  (A) 
- `0.808` in `submission_pack_postfix/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: BChE | D vs A | 0.895 | 0.893 | −0.002 | 0.652 | | AChE/BChE | D vs B | 0.821 | 0.808 | −0.013 | 0.606 | | PIK3CA/mTOR | D vs A | 0.762 | 0.742 | −0.020 | 0.714
- `0.170` in `submission_pack_postfix/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: 0.060, 0.295] | no | | JAK1/JAK2 | pocket A | 0.728 | 0.723 | −0.004 | [−0.177, 0.170] | no | | JAK1/JAK2 | pocket B | 0.588 | 0.730 | 0.142 | [−0.058, 0.334] |
- `0.170` in `submission_pack_postfix/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: OR | heavy | 0.463 | 0.229 | [−0.011, 0.435] | no | | F2/F10 | cLogP | 0.515 | −0.170 | [−0.324, −0.012] | yes | | JAK1/TYK2 | cLogP | 0.580 | −0.215 | [−0.374,
- `9.505` in `submission_pack_postfix/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md`: lots use the unified chemically mapped CalcRMS table (Table S2). The historical 9.505 Å value is not a current result.  ### Figure S4. Label and source robustne

## Verdict

**READY FOR LANGUAGE/EDITORIAL POLISH**

