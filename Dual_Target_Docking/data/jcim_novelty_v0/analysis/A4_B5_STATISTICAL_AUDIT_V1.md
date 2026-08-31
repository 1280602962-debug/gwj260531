# A4 / B5 statistical audit (2026-08-24)

Internal verification. Manuscript numbers must match these deposited tables.
Master index: `tables/MASTER_RESULTS_TABLE.csv`.

## A4 — max vs median pChEMBL (full frozen panels)

Sources: `assay_max_vs_median_{summary,auroc,flips,ligand,agreement}_v1.csv`.
Script: `scripts/assay_aggregation_max_vs_median_v1.py` (API re-fetch; frozen Vina scores).

**Do not mix frozen Table 2 with A4 API-max on EGFR/HER2.**

| Pair | n scored | class flips (θ=6.0) | label agreement | numeric max≠median | frozen `summary_min` | API-max min | API-median min |
|---|---:|---:|---:|---:|---:|---:|---:|
| EGFR/HER2 | 110 | 7 | 103/110 = 93.6% | 40/110 | **0.430** (28/38/32) | **0.417** (29/37/32) | **0.424** (26/35/33) |
| AChE/BChE | 95 | 1 | 94/95 = 98.9% | 13/95 | 0.606 | 0.606 | 0.629 |
| PIK3CA/PIK3CB | 99 | 1 | 98/99 = 99.0% | 25/99 | 0.500 | 0.500 | 0.500 |
| PIK3CA/mTOR | 48 | 0 | 48/48 = 100% | 27/48 | 0.692 | 0.692 | 0.692 |

Denominators are **scored** n (AChE 95/100; PIK3CB 99/100), not construction n.

EGFR cache/API mismatch (n = 1): `EH120_060` / CHEMBL24828. Frozen class A_only (cached pB = 5.58); API max pB = 6.85 → API-max class dual. `n_cache_matches_both_max` = 109/110 on EGFR; 95/95, 99/99, 48/48 on the other pairs.

Allowed claim: primary pair-level conclusions were insensitive to replacing maximum pChEMBL with the median among repeated measurements. Assay-level heterogeneity remains a limitation because pChEMBL values are not assay-equivalent (numeric max≠median is common; class flips at θ = 6.0 are not).

Forbidden: promoting the 27-ligand diagnostic; writing “7 flips” without 7/110; writing Table 2 0.430 → 0.424 as if that were the A4 native contrast (native is API-max 0.417 → median 0.424; frozen vs median is a second, still small, contrast).

## B5 — PIK3CA/PIK3CB receptor realization (2WXF frozen)

Sources: `pocket_matched_PAB_alt{4JPS,5DXT}_v1.csv`; combined `receptor_realization_two_pair_v1.csv`.
Use **deposited CSV CIs**, not ad-hoc recomputes.

| Pair | PIK3CA | kept B | n attempted / success / fail | D/A | D/B | `summary_min` [95% CI] | Δ vs original |
|---|---|---|---|---:|---:|---|---:|
| PIK3CA/mTOR | 4L23 | 4JT6 | 48/48/0 | 0.714 | 0.692 | **0.692 [0.464, 0.802]** | — |
| PIK3CA/mTOR | 4JPS | 4JT6 | 48/48/0 | 0.714 | 0.486 | **0.486 [0.259, 0.692]** | −0.206 |
| PIK3CA/mTOR | 5DXT | 4JT6 | 48/48/0 | 0.714 | 0.505 | **0.505 [0.292, 0.696]** | −0.187 |
| PIK3CA/mTOR | 4L23 | 4JSX | 48/48/0 | 0.639 | 0.692 | **0.639 [0.418, 0.776]** | −0.053 |
| PIK3CA/PIK3CB | 4L23 | 2WXF | **100/99/1** | 0.691 | 0.500 | **0.500 [0.347, 0.648]** | — |
| PIK3CA/PIK3CB | 4JPS | 2WXF | **100/99/1** | 0.691 | 0.707 | **0.691 [0.516, 0.779]** | +0.191 |
| PIK3CA/PIK3CB | 5DXT | 2WXF | **100/99/1** | 0.691 | 0.685 | **0.685 [0.506, 0.768]** | +0.185 |

Weak-arm switch on PIK3CA/PIK3CB: original bottleneck is D/B (4L23) = 0.500; 4JPS raises D/B to 0.707 so the bottleneck becomes frozen D/A (2WXF) = 0.691; 5DXT D/B = 0.685 ≈ D/A.

## PAB_034

- Ligand `PAB_034` / CHEMBL5089694 / **A_only**.
- Original 4L23: `timeout_900s_torsdof=23_E=4` (2WXF success). Table 2 already uses 99/100.
- 4JPS: `timeout_600s`, 668.6 s. 5DXT: `timeout_600s`, 665.0 s.
- Failure is a **docking timeout**, not an experimental label filter, and is **not unique to the alternate crystals**.
- No 100-ligand directional AUROC exists for this pair under any PIK3CA crystal in this protocol. 99-vs-100 sensitivity is not computable from existing scores.

## Aggregation (`summary_min` vs arithmetic, geometric, and harmonic means)

Table S26 now reports all four aggregators. Pair ranking is unchanged (PM > AChE > PIK3CB > EGFR). This is a deterministic sensitivity calculation on the deposited directional AUROCs, not a new docking experiment.

## Claim ceiling after this audit

- Receptor replacement is a **realization-effect / sensitivity analysis**, not robustness, and not a unidirectional collapse.
- Two opposite-direction examples share PIK3CA; this is not a universal law over K = 4.
- max pChEMBL is a **controlled limitation**, not an unresolved validity threat.
