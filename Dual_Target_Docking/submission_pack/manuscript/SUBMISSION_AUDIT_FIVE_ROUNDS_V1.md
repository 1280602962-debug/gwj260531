# Five-round JCIM submission audit

Date: 2026-09-09
Branch: `cursor/jcim-submission-pack-0b1a`
Script: `scripts/audit/audit_submission_five_rounds_v1.py`

This audit compares assembled manuscripts and SI tables to frozen CSVs.
It does not re-bootstrap and does not mint a DOI.

Summary: **111 PASS**, **0 FAIL**, **4 NOTE**.

| Round | Status | Finding |
|---|---|---|
| R1 | PASS | EN Table 2 EGFR/HER2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 AChE/BChE matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 PIK3CA/mTOR matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 F2/F10 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 JAK1/TYK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 JAK1/JAK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 PPARG/PPARA matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 2 PPARA/PPARD matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 EGFR/HER2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 AChE/BChE matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 PIK3CA/mTOR matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 F2/F10 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 JAK1/TYK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 JAK1/JAK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 PPARG/PPARA matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 2 PPARA/PPARD matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 EGFR/HER2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 AChE/BChE matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 PIK3CA/mTOR matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 F2/F10 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 JAK1/TYK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 JAK1/JAK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 PPARG/PPARA matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | EN Table 3 PPARA/PPARD matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 EGFR/HER2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 AChE/BChE matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 PIK3CA/mTOR matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 F2/F10 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 JAK1/TYK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 JAK1/JAK2 matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 PPARG/PPARA matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | ZH Table 3 PPARA/PPARD matches locked CSV (3-dp half-up or half-even) |
| R1 | PASS | Table 1 PIK3CA/mTOR construction/n_scored |
| R1 | PASS | Table 1 AChE/BChE construction/n_scored |
| R1 | PASS | Table 1 EGFR/HER2 construction/n_scored |
| R1 | PASS | Table 1 F2/F10 construction/n_scored |
| R1 | PASS | Table 1 JAK1/TYK2 construction/n_scored |
| R1 | PASS | Table 1 JAK1/JAK2 construction/n_scored |
| R1 | PASS | Table 1 PPARG/PPARA construction/n_scored |
| R1 | PASS | Table 1 PPARA/PPARD construction/n_scored |
| R1 | PASS | STATISTICAL_LOCK EGFR/HER2 |
| R1 | PASS | STATISTICAL_LOCK AChE/BChE |
| R1 | PASS | STATISTICAL_LOCK PIK3CA/mTOR |
| R1 | PASS | STATISTICAL_LOCK F2/F10 |
| R1 | PASS | STATISTICAL_LOCK JAK1/TYK2 |
| R1 | PASS | STATISTICAL_LOCK JAK1/JAK2 |
| R1 | PASS | STATISTICAL_LOCK PPARG/PPARA |
| R1 | PASS | STATISTICAL_LOCK PPARA/PPARD |
| R2 | PASS | Table S4 EGFR/HER2 pocket A points 0.430 / 0.808 / 0.378 |
| R2 | PASS | Table S4 EGFR/HER2 pocket A CI [0.205, 0.547] |
| R2 | PASS | Table S4 JAK1/TYK2 pocket A points 0.365 / 0.809 / 0.444 |
| R2 | PASS | Table S4 JAK1/TYK2 pocket A CI [0.263, 0.620] |
| R2 | PASS | EGFR document-cluster CI present |
| R2 | PASS | EGFR scaffold-cluster CI present |
| R2 | PASS | JAK1 document-cluster CI present |
| R2 | PASS | JAK1 scaffold-cluster CI present |
| R2 | PASS | Table S6 EGFR/HER2 Δ 0.170 |
| R2 | PASS | Table S6 AChE/BChE Δ 0.161 |
| R2 | PASS | Table S6 PIK3CA/mTOR Δ 0.090 |
| R2 | PASS | Table S6 F2/F10 Δ -0.031 |
| R2 | PASS | Table S6 JAK1/TYK2 Δ -0.065 |
| R2 | PASS | Table S6 JAK1/JAK2 Δ -0.019 |
| R2 | PASS | Table S6 PPARG/PPARA Δ 0.030 |
| R2 | PASS | Table S6 PPARA/PPARD Δ 0.012 |
| R2 | PASS | S11b EGFR/HER2 180 / 10 / 20 |
| R2 | PASS | S11b AChE/BChE 4 / 8 / 14 |
| R2 | PASS | S11b PIK3CA/mTOR 91 / 4 / 1 |
| R2 | PASS | S11b F2/F10 46 / 15 / 16 |
| R2 | PASS | S11b JAK1/TYK2 323 / 7 / 103 |
| R2 | PASS | S11b JAK1/JAK2 928 / 40 / 14 |
| R2 | PASS | S11b PPARG/PPARA 0 / 0 / 1 |
| R2 | PASS | S11b PPARA/PPARD 0 / 1 / 0 |
| R2 | PASS | ECFP increment max \|Δ\| = 0.0234 (manuscript 0.023) |
| R2 | PASS | EGFR independent GNINA 0.220 / 0.783 (0.7825) |
| R2 | NOTE | EGFR GNINA Dual-vs-neither n_neg=11 (EH120_109 failed); Vina Table 3 uses n=12 |
| R2 | PASS | JAK1/TYK2 independent GNINA 0.317 / 0.705 |
| R2 | PASS | JAK1 GNINA CI [0.183, 0.463] |
| R2 | PASS | Table S9 EGFR/HER2 median 0.373 |
| R2 | PASS | Table S9 F2/F10 median 0.366 |
| R2 | PASS | Table S9 PPARG/PPARA median 0.651 |
| R2 | PASS | AChE holdout summary_min 0.618 |
| R2 | PASS | PPARG holdout 0.535 [0.350, 0.717] from table2_comparable_by_channel |
| R3 | PASS | EGFR/HER2 EN/ZH table-number sets match |
| R3 | PASS | AChE/BChE EN/ZH table-number sets match |
| R3 | PASS | PIK3CA/mTOR EN/ZH table-number sets match |
| R3 | PASS | F2/F10 EN/ZH table-number sets match |
| R3 | PASS | JAK1/TYK2 EN/ZH table-number sets match |
| R3 | PASS | JAK1/JAK2 EN/ZH table-number sets match |
| R3 | PASS | PPARG/PPARA EN/ZH table-number sets match |
| R3 | PASS | PPARA/PPARD EN/ZH table-number sets match |
| R3 | PASS | flagship token 0.378 in both manuscripts |
| R3 | PASS | flagship token 0.444 in both manuscripts |
| R3 | PASS | flagship token 0.205 in both manuscripts |
| R3 | PASS | flagship token 0.547 in both manuscripts |
| R3 | PASS | flagship token 0.263 in both manuscripts |
| R3 | PASS | flagship token 0.620 in both manuscripts |
| R3 | PASS | flagship token 0.083 in both manuscripts |
| R3 | PASS | flagship token 0.529 in both manuscripts |
| R3 | PASS | flagship token 0.168 in both manuscripts |
| R3 | PASS | flagship token 0.562 in both manuscripts |
| R3 | PASS | flagship token 0.234 in both manuscripts |
| R3 | PASS | flagship token 0.633 in both manuscripts |
| R3 | PASS | ZH JAK1 document-cluster lower bound present |
| R4 | PASS | fig1C labels no longer include a PIK3CB tick without a bar |
| R4 | PASS | fig1C complete-case range uses the three original maps only |
| R4 | NOTE | plot_jcim_si_composites_v1.py still ticks withdrawn PIK3CB (S1–S3/S9/S10 archive figures; do not submit as eight-row primaries) |
| R4 | NOTE | v1 figure script retained as archive; submission figures must come from v3 only |
| R5 | PASS | EN manuscript has no PIK3CA/PIK3CB primary-set mention |
| R5 | PASS | ZH manuscript has no PIK3CA/PIK3CB primary-set mention |
| R5 | PASS | bootstrap_primary OK |
| R5 | PASS | validate_revision_v1 OK |
| R5 | PASS | checksum --check OK |
| R5 | PASS | MASTER Table 2 has exactly the eight primary pairs |
| R5 | PASS | n_scored D+A+B matches Table 2 class counts |
| R5 | NOTE | J0 θ=6.0 census docked_in_this_paper != 4 (historical J0-era count may include withdrawn PIK3CB) |

## Blocking rule

Publication-facing numbers must match the locked CSVs after three-decimal rounding.
Fig1C must not plot a withdrawn PIK3CA/PIK3CB bar or read that pair from a CSV that no longer contains it.
Assembled EN/ZH manuscripts must not list PIK3CA/PIK3CB as a primary pair.

## Canonical sources

- Table 2 original three: `unified_threshold_sensitivity_v2.csv` `theta_6.0`
- Table 2 / Table 3 five pairs: `five_pair_stack_v1/table2_comparable_theta6_v1.csv`
- Table 3 original three: `formulation_conventional_vs_directional_v1.csv`
- Table S4: `formulation_equal_score_negative_v1.csv` + `equal_score_negative_s34_v1.csv`
- Table S6: `wrong_pocket_paired_delta_bootstrap_v1.csv` + `wrong_pocket_by_channel_v1.csv`
- Table S11b: `external_slice_summary_v1.csv`
- JAK1/TYK2 independent GNINA: `table2_comparable_by_channel_v1.csv` `gnina_independent_jak1_tyk2`
- EGFR/PIK3CA independent GNINA: `independent_dock_formulation_v1.csv`

## Remediation on this branch

- Table 2 EGFR/HER2 TPSA displayed 0.427 for CSV 0.4275, which is not valid under half-up or half-even; corrected to **0.428** in EN/ZH Table 2.
- Fig1C had five x-tick labels (including withdrawn PIK3CA/PIK3CB) but only four J0 bars, and read complete-case overlap for a pair no longer in that CSV. Bars/labels now match `j0_strict_label_supply.csv`; overlap uses the three original maps (14.5%–34.0%).
- Figure verify locks now match the current J0 scrape (48 pairs, 3 thick) and θ=6.0 census (`directional_n10` = 16, `docked_in_this_paper` = 3). The formal Figure 1 funnel shows J0 scrape → min selective ≥50 → eight-pair evaluation set and no longer records historically docked 4 → PIK3CB withdrawal.
- `plot_jcim_si_composites_v1.py` S1–S3/S9/S10 remain original-set archive figures (may still tick PIK3CB) and were moved to `figures/jcim_article/archive_original_set/`.
- Submission slice: `submission_pack/`.

