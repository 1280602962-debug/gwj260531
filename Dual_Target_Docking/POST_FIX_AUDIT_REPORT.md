# POST_FIX_AUDIT_REPORT

Summary: **95 PASS**, **3 WARNING**, **0 FAIL**.

Canonical protocol: EGFR/HER2 3POZ/3RCD cognate heavy-atom AABB+5 Å (min 20 Å); AChE/BChE no-ChEMBL-ID-prefix-cap panel.
Pre-fix boxes/results remain in `data/_legacy_archive/` and are not used as official results.

## Invariant tests

- **PASS** TEST 1 EGFR/HER2 D/A uses B score — AUROC=0.660714
- **PASS** TEST 2 EGFR/HER2 D/B uses A score — AUROC=0.323661
- **PASS** TEST 3 EGFR/HER2 dual is positive class
- **PASS** TEST 4 EGFR/HER2 AUROC(y,-score)≈1-AUROC — sum=1
- **PASS** TEST 5 EGFR/HER2 summary_min=min(two directional points) — 0.323661
- **PASS** TEST 6 EGFR/HER2 pocket A Δ uses identical score channel — delta=0.462054
- **PASS** TEST 7 EGFR/HER2 Smean=(SA+SB)/2
- **PASS** TEST 8 EGFR/HER2 Top10% D+A+B+N=k — 1/5/5/0 k=11
- **PASS** TEST 9 EGFR/HER2 k=ceil(0.10 n) — n=110 k=11
- **PASS** TEST 10 EGFR/HER2 EF formula — 0.357143
- **PASS** TEST 18 EGFR/HER2 matched/mismatched only swaps assignment — matched=0.3237 mismatched=0.2679
- **PASS** TEST 1 AChE/BChE D/A uses B score — AUROC=0.652422
- **PASS** TEST 2 AChE/BChE D/B uses A score — AUROC=0.605820
- **PASS** TEST 3 AChE/BChE dual is positive class
- **PASS** TEST 4 AChE/BChE AUROC(y,-score)≈1-AUROC — sum=1
- **PASS** TEST 5 AChE/BChE summary_min=min(two directional points) — 0.605820
- **PASS** TEST 6 AChE/BChE pocket A Δ uses identical score channel — delta=-0.015697
- **PASS** TEST 7 AChE/BChE Smean=(SA+SB)/2
- **PASS** TEST 8 AChE/BChE Top10% D+A+B+N=k — 5/3/1/1 k=10
- **PASS** TEST 9 AChE/BChE k=ceil(0.10 n) — n=96 k=10
- **PASS** TEST 10 AChE/BChE EF formula — 1.777778
- **PASS** TEST 18 AChE/BChE matched/mismatched only swaps assignment — matched=0.6058 mismatched=0.4288
- **PASS** TEST 1 JAK1/JAK2 D/A uses B score — AUROC=0.588379
- **PASS** TEST 2 JAK1/JAK2 D/B uses A score — AUROC=0.727539
- **PASS** TEST 3 JAK1/JAK2 dual is positive class
- **PASS** TEST 4 JAK1/JAK2 AUROC(y,-score)≈1-AUROC — sum=1.000000000000
- **PASS** TEST 5 JAK1/JAK2 summary_min=min(two directional points) — 0.588379
- **PASS** TEST 7 JAK1/JAK2 Smean=(SA+SB)/2
- **PASS** TEST 6 JAK1/JAK2 pocket A Δ uses identical score channel
- **PASS** TEST 8 JAK1/JAK2 Top10% D+A+B+N=k — 6/5/0/0 k=11
- **PASS** TEST 9 JAK1/JAK2 k=ceil(0.10 n) — n=110 k=11
- **PASS** TEST 10 JAK1/JAK2 EF formula — 1.875000
- **PASS** TEST 18 JAK1/JAK2 matched/mismatched only swaps assignment — Δsummary_min=-0.0190 CI=[-0.088,0.054]
- **PASS** TEST 1 JAK1/TYK2 D/A uses B score — AUROC=0.575101
- **PASS** TEST 2 JAK1/TYK2 D/B uses A score — AUROC=0.364919
- **PASS** TEST 3 JAK1/TYK2 dual is positive class
- **PASS** TEST 4 JAK1/TYK2 AUROC(y,-score)≈1-AUROC — sum=1.000000000000
- **PASS** TEST 5 JAK1/TYK2 summary_min=min(two directional points) — 0.364919
- **PASS** TEST 7 JAK1/TYK2 Smean=(SA+SB)/2
- **PASS** TEST 6 JAK1/TYK2 pocket A Δ uses identical score channel
- **PASS** TEST 8 JAK1/TYK2 Top10% D+A+B+N=k — 1/3/7/0 k=11
- **PASS** TEST 9 JAK1/TYK2 k=ceil(0.10 n) — n=109 k=11
- **PASS** TEST 10 JAK1/TYK2 EF formula — 0.319648
- **PASS** TEST 18 JAK1/TYK2 matched/mismatched only swaps assignment — Δsummary_min=-0.0645 CI=[-0.151,0.033]
- **PASS** TEST 1 PIK3CA/mTOR D/A uses B score — AUROC=0.714286
- **PASS** TEST 2 PIK3CA/mTOR D/B uses A score — AUROC=0.692130
- **PASS** TEST 3 PIK3CA/mTOR dual is positive class
- **PASS** TEST 4 PIK3CA/mTOR AUROC(y,-score)≈1-AUROC — sum=1.000000000000
- **PASS** TEST 5 PIK3CA/mTOR summary_min=min(two directional points) — 0.692130
- **PASS** TEST 7 PIK3CA/mTOR Smean=(SA+SB)/2
- **PASS** TEST 6 PIK3CA/mTOR pocket A Δ uses identical score channel
- **PASS** TEST 8 PIK3CA/mTOR Top10% D+A+B+N=k — 4/1/0/0 k=5
- **PASS** TEST 9 PIK3CA/mTOR k=ceil(0.10 n) — n=48 k=5
- **PASS** TEST 10 PIK3CA/mTOR EF formula — 2.133333
- **PASS** TEST 18 PIK3CA/mTOR matched/mismatched only swaps assignment — Δsummary_min=0.0903 CI=[-0.117,0.286]
- **PASS** TEST 1 F2/F10 D/A uses B score — AUROC=0.413306
- **PASS** TEST 2 F2/F10 D/B uses A score — AUROC=0.344758
- **PASS** TEST 3 F2/F10 dual is positive class
- **PASS** TEST 4 F2/F10 AUROC(y,-score)≈1-AUROC — sum=1.000000000000
- **PASS** TEST 5 F2/F10 summary_min=min(two directional points) — 0.344758
- **PASS** TEST 7 F2/F10 Smean=(SA+SB)/2
- **PASS** TEST 6 F2/F10 pocket A Δ uses identical score channel
- **PASS** TEST 8 F2/F10 Top10% D+A+B+N=k — 4/1/6/0 k=11
- **PASS** TEST 9 F2/F10 k=ceil(0.10 n) — n=107 k=11
- **PASS** TEST 10 F2/F10 EF formula — 1.255132
- **PASS** TEST 18 F2/F10 matched/mismatched only swaps assignment — Δsummary_min=-0.0307 CI=[-0.115,0.046]
- **PASS** TEST 1 PPARG/PPARA D/A uses B score — AUROC=0.649194
- **PASS** TEST 2 PPARG/PPARA D/B uses A score — AUROC=0.706055
- **PASS** TEST 3 PPARG/PPARA dual is positive class
- **PASS** TEST 4 PPARG/PPARA AUROC(y,-score)≈1-AUROC — sum=1.000000000000
- **PASS** TEST 5 PPARG/PPARA summary_min=min(two directional points) — 0.649194
- **PASS** TEST 7 PPARG/PPARA Smean=(SA+SB)/2
- **PASS** TEST 6 PPARG/PPARA pocket A Δ uses identical score channel
- **PASS** TEST 8 PPARG/PPARA Top10% D+A+B+N=k — 7/3/0/1 k=11
- **PASS** TEST 9 PPARG/PPARA k=ceil(0.10 n) — n=109 k=11
- **PASS** TEST 10 PPARG/PPARA EF formula — 2.167614
- **PASS** TEST 18 PPARG/PPARA matched/mismatched only swaps assignment — Δsummary_min=0.0296 CI=[-0.079,0.162]
- **PASS** TEST 1 PPARA/PPARD D/A uses B score — AUROC=0.646484
- **PASS** TEST 2 PPARA/PPARD D/B uses A score — AUROC=0.446289
- **PASS** TEST 3 PPARA/PPARD dual is positive class
- **PASS** TEST 4 PPARA/PPARD AUROC(y,-score)≈1-AUROC — sum=1.000000000000
- **PASS** TEST 5 PPARA/PPARD summary_min=min(two directional points) — 0.446289
- **PASS** TEST 7 PPARA/PPARD Smean=(SA+SB)/2
- **PASS** TEST 6 PPARA/PPARD pocket A Δ uses identical score channel
- **PASS** TEST 8 PPARA/PPARD Top10% D+A+B+N=k — 5/2/3/1 k=11
- **PASS** TEST 9 PPARA/PPARD k=ceil(0.10 n) — n=110 k=11
- **PASS** TEST 10 PPARA/PPARD EF formula — 1.562500
- **PASS** TEST 18 PPARA/PPARD matched/mismatched only swaps assignment — Δsummary_min=0.0117 CI=[-0.088,0.143]
- **PASS** TEST 15 all primary analyses use θ=6.0 labels — EGFR θ=6.0; AChE sampled strict 6.5/5.5 coincides on deposited panel
- **PASS** TEST 16 strict 6.5/5.5 not used as primary labels
- **PASS** TEST 17 primary Vina score = mode 1 — canonical corrected-box mode-1 affinities
- **PASS** TEST 12 main/holdout ligand overlap = 0 — n_overlap=0
- **PASS** TEST 11 ECFP4 train/test scaffold overlap = 0 — leaked=0
- **PASS** TEST boxes canonical heavy-atom AABB+5Å min20 — 3POZ=canonical_post_remediation 3RCD=canonical_post_remediation
- **WARNING** TEST 19 receptor substitution and independent GNINA are not mixed — independent GNINA still running
- **PASS** TEST 20 reported CI procedures match Methods for summary_min — B=2000 pooled dual+A+B min-inside-replicate
- **WARNING** TEST 13 Figure/Table n values correspond to real ligand IDs — manuscript/figures not yet overwritten
- **WARNING** TEST 14 manuscript numeric claims have result-table provenance — manuscript not yet overwritten

## Core qualitative conclusions

a. directional docking performance varies across pairs/directions — **HOLDS** (EGFR D/B now 0.324 with CI excluding 0.5 below; other pairs unchanged).
b. dual-vs-neither does not reliably represent exclusion of single-target-active ligands — **HOLDS** (EGFR fixed-score ΔAUROC 0.462, CI excludes 0).
c. ligand-only chemistry carries substantial class information — **HOLDS** (AChE ECFP4 OOF 0.897 / 0.843; EGFR ligand-only ECFP4 not rerun because membership unchanged).
d. docking adds limited incremental discrimination to ECFP4 — **HOLDS** (AChE increment ≈ 0; EGFR increment pending staging incremental table).
e. matched-pocket advantage is not consistently reproduced — **HOLDS, pair-level EGFR inference updated**: corrected EGFR 95% CI includes 0, so EGFR/HER2 is **not** reported as having a matched-pocket advantage. AChE main-panel CI still excludes 0; AChE holdout CI includes 0.
f. top-ranking dual enrichment is target-pair dependent — **HOLDS** (EGFR EF_dual,10%=0.357; AChE 1.778).

## CORE STATUS

**SUPPORTED BUT NUMERIC/PAIR-LEVEL INTERPRETATION UPDATED**

This is not a core-project failure. Six core conclusions remain. The EGFR/HER2 matched-minus-mismatched interval no longer excludes 0; official text must follow the corrected interval.

