# Descriptor result impact

Date: 2026-09-20

## Old vs corrected nested `summary_min`

| Pair | Old invalid pooled-raw OOF | Corrected P(dual) OOF | Absolute change | Corrected 95% CI |
|------|---------------------------:|----------------------:|----------------:|------------------|
| EGFR/HER2 | 0.4171 | 0.6325 | +0.2154 | [0.4666, 0.6926] |
| JAK1/JAK2 | 0.5996 | 0.4380 | −0.1616 | [0.2939, 0.5913] |
| JAK1/TYK2 | 0.5630 | 0.4587 | −0.1043 | [0.3019, 0.5393] |
| PIK3CA/mTOR | 0.3542 | 0.3889 | +0.0347 | [0.1943, 0.5927] |
| AChE/BChE | 0.7422 | 0.6496 | −0.0926 | [0.4972, 0.7964] |
| F2/F10 | 0.4768 | 0.3538 | −0.1230 | [0.2278, 0.4925] |
| PPARG/PPARA | 0.4595 | 0.3579 | −0.1016 | [0.2266, 0.4990] |
| PPARA/PPARD | 0.4756 | 0.4438 | −0.0318 | [0.3008, 0.5410] |

Old numbers pooled heterogeneous raw descriptors. They are not a target to retain.

## Comparison with other rankings (corrected)

| Pair | Vina `summary_min` | Descriptive best-descriptor `summary_min` | Nested descriptor OOF `summary_min` | ECFP4 weaker-arm OOF |
|------|-------------------:|------------------------------------------:|------------------------------------:|---------------------:|
| EGFR/HER2 | 0.334 | 0.474 (cLogP) | 0.632 | 0.822 / 0.879 |
| JAK1/JAK2 | 0.588 | 0.578 (heavy) | 0.438 | 0.915 / 0.968 |
| JAK1/TYK2 | 0.365 | 0.580 (cLogP) | 0.459 | 0.846 / 0.902 |
| PIK3CA/mTOR | 0.692 | 0.463 (heavy) | 0.389 | 0.762 / 0.889 |
| AChE/BChE | 0.606 | 0.742 (TPSA) | 0.650 | 0.895 / 0.821 |
| F2/F10 | 0.345 | 0.509 (cLogP) | 0.354 | 0.943 / 0.693 |
| PPARG/PPARA | 0.649 | 0.627 (TPSA) | 0.358 | 0.833 / 0.668 |
| PPARA/PPARD | 0.446 | 0.564 (cLogP) | 0.444 | 0.932 / 0.858 |

ECFP max |Δ AUROC| remains 0.0234 (ECFP4 pipeline not shared).

## Claim

Former implied reading: “simple ligand properties explain a substantial portion of the apparent discrimination.”

**Qualitative conclusion changes: yes, it is weakened.**

- Full-panel TPSA on AChE/BChE (0.742 / 0.801) remains a valid *descriptive* observation.
- The selection-adjusted nested OOF is lower on several pairs. JAK1/JAK2 falls to 0.438 [0.294, 0.591], below Vina 0.588. PPARG/PPARA nested `summary_min` is 0.358 [0.227, 0.499]. Those docking summaries are not explained by the nested single-descriptor baseline.
- AChE/BChE nested OOF 0.650 [0.497, 0.796] is no longer interchangeable with the full-panel TPSA screen, and its interval includes 0.5.
- ECFP4 remains the stronger ligand-based comparator.

## Affected manuscript sentences

| Location | Action |
|----------|--------|
| Methods 2.6.1 nested-CV description | WORDING_UPDATE — inner CV + P(dual) scale |
| Results 3.3 descriptive TPSA / heavy / Δ CIs | UNCHANGED (descriptive screen) |
| Results 3.3 after descriptive paragraph | NUMERIC_UPDATE — added nested OOF and weakened claim |
| Conclusion “simple physicochemical descriptors already separated some states on AChE/BChE” | UNCHANGED as a descriptive statement |
| SI Table S5 nested block | NUMERIC_UPDATE — new table and method prose |
| Table S7a GNINA weaker-arm CI | NUMERIC_UPDATE — directional CI, not copied `summary_min` CI |

Core docking conclusion (negative-class / ranking paper) is unchanged. Table 2 and Table 3 are unchanged.
