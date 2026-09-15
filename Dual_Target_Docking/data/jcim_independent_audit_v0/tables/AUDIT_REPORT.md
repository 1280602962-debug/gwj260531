# Independent Dual-Target Docking Statistical Audit

This report recomputes primary metrics from ligand-level pChEMBL and raw Vina affinities.
It does not import production analysis modules.

- Canonical ligand table: `data/jcim_independent_audit_v0/tables/canonical_ligand_level_v1.csv`
- Machine-readable metrics: `data/jcim_independent_audit_v0/tables/audit_master_metrics.csv`
- Finding log: `data/jcim_independent_audit_v0/tables/audit_findings_v1.csv`

Summary: **142 PASS**, **15 WARNING**, **0 FAIL**.

## Method

- Primary labels: dual if both pChEMBL ≥ 6.0, A-only / B-only / neither otherwise (θ = 6.0).
- Strict 6.5/5.5 is used only as the sampling rule for six pairs and for unused-pool holdouts.
- Scores: S = −E_Vina (higher is better). D vs A-only uses pocket B; D vs B-only uses pocket A.
- summary_min is min of the two directional point AUROCs. Bootstrap (B=2000) resamples the dual+A+B pool once per replicate and recomputes both arms then the min; dual ligands are shared inside a replicate.
- Fixed-score ΔAUROC = AUROC(D vs neither) − AUROC(D vs single-target-active) at one pocket, with a shared dual draw.
- Two-pocket mean uses only complete-case ligands. Top-10% uses k=⌈0.10 n⌉, descending S_mean, ties by ligand_id.
- EF_dual,10% = (D_top/k) / (D_total/n) on the full four-state panel.
- ECFP4: Morgan radius 2, 2048 bits; GroupKFold on Bemis–Murcko scaffolds; logistic C=1.0; AUROC from out-of-fold scores. ECFP4 and ECFP4+docking share the same splits.

## 1. Directional AUROC (independent recompute)

| Pair | n (D/A/B) | D vs A-only (pocket B) [95% CI] | D vs B-only (pocket A) [95% CI] | summary_min [95% CI] |
|------|-----------|--------------------------------|--------------------------------|----------------------|
| EGFR/HER2 | 28 / 38 / 32 | 0.666 [0.524, 0.793] | 0.430 [0.282, 0.579] | 0.430 [0.282, 0.578] |
| JAK1/JAK2 | 32 / 32 / 32 | 0.588 [0.444, 0.729] | 0.728 [0.595, 0.848] | 0.588 [0.444, 0.725] |
| JAK1/TYK2 | 31 / 32 / 32 | 0.575 [0.434, 0.725] | 0.365 [0.231, 0.505] | 0.365 [0.231, 0.503] |
| PIK3CA/mTOR | 18 / 14 / 12 | 0.714 [0.506, 0.899] | 0.692 [0.495, 0.874] | 0.692 [0.470, 0.813] |
| AChE/BChE | 27 / 25 / 28 | 0.650 [0.483, 0.801] | 0.606 [0.442, 0.751] | 0.606 [0.437, 0.730] |
| F2/F10 | 31 / 32 / 32 | 0.413 [0.259, 0.562] | 0.345 [0.214, 0.486] | 0.345 [0.211, 0.477] |
| PPARG/PPARA | 32 / 31 / 32 | 0.649 [0.507, 0.778] | 0.706 [0.569, 0.833] | 0.649 [0.504, 0.751] |
| PPARA/PPARD | 32 / 32 / 32 | 0.646 [0.504, 0.776] | 0.446 [0.296, 0.586] | 0.446 [0.296, 0.584] |

## 2–3. Fixed-score ΔAUROC (independent recompute)

| Pair | pocket | D vs selective | D vs neither | Δ |
|------|--------|---------------:|-------------:|--:|
| EGFR/HER2 | B | 0.666 | 0.720 | +0.054 |
| EGFR/HER2 | A | 0.430 | 0.808 | +0.378 |
| JAK1/JAK2 | B | 0.588 | 0.730 | +0.142 |
| JAK1/JAK2 | A | 0.728 | 0.723 | -0.004 |
| JAK1/TYK2 | B | 0.575 | 0.705 | +0.130 |
| JAK1/TYK2 | A | 0.365 | 0.809 | +0.444 |
| PIK3CA/mTOR | B | 0.714 | 0.583 | -0.131 |
| PIK3CA/mTOR | A | 0.692 | 0.472 | -0.220 |
| AChE/BChE | B | 0.650 | 0.709 | +0.058 |
| AChE/BChE | A | 0.606 | 0.590 | -0.016 |
| F2/F10 | B | 0.413 | 0.528 | +0.115 |
| F2/F10 | A | 0.345 | 0.508 | +0.163 |
| PPARG/PPARA | B | 0.649 | 0.644 | -0.005 |
| PPARG/PPARA | A | 0.706 | 0.759 | +0.053 |
| PPARA/PPARD | B | 0.646 | 0.665 | +0.019 |
| PPARA/PPARD | A | 0.446 | 0.484 | +0.038 |

## 4. Two-pocket mean dual-vs-neither

| Pair | n_dual / n_neither | D vs neither (S_mean) |
|------|--------------------|------------------------|
| EGFR/HER2 | 28 / 12 | 0.756 |
| JAK1/JAK2 | 32 / 14 | 0.730 |
| JAK1/TYK2 | 31 / 14 | 0.770 |
| PIK3CA/mTOR | 18 / 4 | 0.514 |
| AChE/BChE | 27 / 15 | 0.649 |
| F2/F10 | 31 / 12 | 0.519 |
| PPARG/PPARA | 32 / 14 | 0.685 |
| PPARA/PPARD | 32 / 14 | 0.565 |

## 6. Incremental ECFP4 (independently retrained OOF)

| Pair | contrast | n | ECFP4 OOF | ECFP4+docking OOF | Δ | deposited Δ |
|------|----------|--:|----------:|------------------:|--:|------------:|
| EGFR/HER2 | D_vs_A | 66 | 0.745 | 0.751 | +0.0056 | +0.0056 |
| EGFR/HER2 | D_vs_B | 60 | 0.890 | 0.887 | -0.0022 | -0.0022 |
| JAK1/JAK2 | D_vs_A | 64 | 0.915 | 0.916 | +0.0015 | +0.0014 |
| JAK1/JAK2 | D_vs_B | 64 | 0.968 | 0.969 | +0.0010 | +0.0010 |
| JAK1/TYK2 | D_vs_A | 63 | 0.846 | 0.840 | -0.0060 | -0.0061 |
| JAK1/TYK2 | D_vs_B | 63 | 0.902 | 0.903 | +0.0010 | +0.0010 |
| PIK3CA/mTOR | D_vs_A | 32 | 0.762 | 0.742 | -0.0198 | -0.0198 |
| PIK3CA/mTOR | D_vs_B | 30 | 0.889 | 0.898 | +0.0093 | +0.0092 |
| AChE/BChE | D_vs_A | 52 | 0.895 | 0.893 | -0.0015 | -0.0015 |
| AChE/BChE | D_vs_B | 55 | 0.821 | 0.808 | -0.0132 | -0.0132 |
| F2/F10 | D_vs_A | 63 | 0.936 | 0.929 | -0.0071 | -0.0071 |
| F2/F10 | D_vs_B | 63 | 0.665 | 0.686 | +0.0212 | +0.0212 |
| PPARG/PPARA | D_vs_A | 63 | 0.833 | 0.820 | -0.0131 | -0.0131 |
| PPARG/PPARA | D_vs_B | 64 | 0.668 | 0.676 | +0.0078 | +0.0078 |
| PPARA/PPARD | D_vs_A | 64 | 0.932 | 0.928 | -0.0039 | -0.0039 |
| PPARA/PPARD | D_vs_B | 64 | 0.858 | 0.835 | -0.0234 | -0.0234 |

## 7. Matched-minus-mismatched

| Pair | set | matched summary_min | mismatched | Δ |
|------|-----|--------------------:|-----------:|--:|
| EGFR/HER2 | main | 0.430 | 0.260 | +0.170 |
| JAK1/JAK2 | main | 0.588 | 0.607 | -0.019 |
| JAK1/TYK2 | main | 0.365 | 0.429 | -0.065 |
| PIK3CA/mTOR | main | 0.692 | 0.602 | +0.090 |
| AChE/BChE | main | 0.606 | 0.444 | +0.161 |
| F2/F10 | main | 0.345 | 0.376 | -0.031 |
| PPARG/PPARA | main | 0.649 | 0.620 | +0.030 |
| PPARA/PPARD | main | 0.446 | 0.435 | +0.012 |
| AChE/BChE | holdout | 0.618 | 0.642 | -0.025 |
| PIK3CA/mTOR | holdout | 0.765 | 0.787 | -0.022 |
| JAK1/JAK2 | holdout | 0.619 | 0.611 | +0.008 |
| JAK1/TYK2 | holdout | 0.475 | 0.450 | +0.025 |
| F2/F10 | holdout | 0.392 | 0.471 | -0.079 |
| PPARG/PPARA | holdout | 0.535 | 0.529 | +0.006 |
| PPARA/PPARD | holdout | 0.445 | 0.295 | +0.150 |

## 8–9. Top-10% composition and EF_dual,10%

| Pair | n | k=ceil(0.10 n) | D / A / B / N in top-10% | EF_dual,10% |
|------|--:|---------------:|--------------------------|------------:|
| EGFR/HER2 | 110 | 11 | 1 / 5 / 5 / 0 | 0.357 |
| JAK1/JAK2 | 110 | 11 | 6 / 5 / 0 / 0 | 1.875 |
| JAK1/TYK2 | 109 | 11 | 1 / 3 / 7 / 0 | 0.320 |
| PIK3CA/mTOR | 48 | 5 | 4 / 1 / 0 / 0 | 2.133 |
| AChE/BChE | 95 | 10 | 5 / 3 / 1 / 1 | 1.759 |
| F2/F10 | 107 | 11 | 4 / 1 / 6 / 0 | 1.255 |
| PPARG/PPARA | 109 | 11 | 7 / 3 / 0 / 1 | 2.168 |
| PPARA/PPARD | 110 | 11 | 5 / 2 / 3 / 1 | 1.562 |

## PASS

- **primary experimental states equal theta=6.0 from pChEMBL** (`EGFR/HER2`) independent=110
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`EGFR/HER2`)
- **Bemis-Murcko scaffold present for every scored ligand** (`EGFR/HER2`) independent=110
- **primary experimental states equal theta=6.0 from pChEMBL** (`JAK1/JAK2`) independent=110
- **strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel** (`JAK1/JAK2`)
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`JAK1/JAK2`)
- **Bemis-Murcko scaffold present for every scored ligand** (`JAK1/JAK2`) independent=110
- **primary experimental states equal theta=6.0 from pChEMBL** (`JAK1/TYK2`) independent=109
- **strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel** (`JAK1/TYK2`)
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`JAK1/TYK2`)
- **Bemis-Murcko scaffold present for every scored ligand** (`JAK1/TYK2`) independent=109
- **primary experimental states equal theta=6.0 from pChEMBL** (`PIK3CA/mTOR`) independent=48
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`PIK3CA/mTOR`)
- **Bemis-Murcko scaffold present for every scored ligand** (`PIK3CA/mTOR`) independent=48
- **primary experimental states equal theta=6.0 from pChEMBL** (`AChE/BChE`) independent=95
- **strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel** (`AChE/BChE`)
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`AChE/BChE`)
- **Bemis-Murcko scaffold present for every scored ligand** (`AChE/BChE`) independent=95
- **primary experimental states equal theta=6.0 from pChEMBL** (`F2/F10`) independent=107
- **strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel** (`F2/F10`)
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`F2/F10`)
- **Bemis-Murcko scaffold present for every scored ligand** (`F2/F10`) independent=107
- **primary experimental states equal theta=6.0 from pChEMBL** (`PPARG/PPARA`) independent=109
- **strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel** (`PPARG/PPARA`)
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`PPARG/PPARA`)
- **primary experimental states equal theta=6.0 from pChEMBL** (`PPARA/PPARD`) independent=110
- **strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel** (`PPARA/PPARD`)
- **raw Vina affinity sign-corrected S=-E matches membership scores** (`PPARA/PPARD`)
- **summary_min equals min(D/A, D/B) point estimates** (`EGFR/HER2`) independent=0.4296875
- **sign-flip AUROC satisfies 1-AUC** (`EGFR/HER2`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`EGFR/HER2`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`EGFR/HER2`) independent=smin_hi=0.577505 weaker_arm_hi=0.579012
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`EGFR/HER2`) independent=[0.2818,0.5775]
- **EGFR/HER2 pocket-A dual-vs-B-only 0.430 vs neither 0.808, Δ=0.378** (`EGFR/HER2`) independent=0.3783482142857143
- **matched assignment is D/A pocket B and D/B pocket A** (`EGFR/HER2`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`EGFR/HER2`) independent=1/5/5/0 of k=11
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`EGFR/HER2`) independent=0.3571428571428572
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`EGFR/HER2`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`JAK1/JAK2`) independent=0.58837890625
- **sign-flip AUROC satisfies 1-AUC** (`JAK1/JAK2`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`JAK1/JAK2`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`JAK1/JAK2`) independent=smin_hi=0.724579 weaker_arm_hi=0.729207
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`JAK1/JAK2`) independent=[0.4444,0.7246]
- **matched assignment is D/A pocket B and D/B pocket A** (`JAK1/JAK2`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`JAK1/JAK2`) independent=6/5/0/0 of k=11
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`JAK1/JAK2`) independent=1.875
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`JAK1/JAK2`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`JAK1/TYK2`) independent=0.3649193548387097
- **sign-flip AUROC satisfies 1-AUC** (`JAK1/TYK2`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`JAK1/TYK2`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`JAK1/TYK2`) independent=smin_hi=0.503032 weaker_arm_hi=0.505059
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`JAK1/TYK2`) independent=[0.2306,0.5030]
- **JAK1/TYK2 pocket-A fixed-score ΔAUROC = 0.444** (`JAK1/TYK2`) independent=0.4438364055299539
- **matched assignment is D/A pocket B and D/B pocket A** (`JAK1/TYK2`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`JAK1/TYK2`) independent=1/3/7/0 of k=11
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`JAK1/TYK2`) independent=0.3196480938416422
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`JAK1/TYK2`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`PIK3CA/mTOR`) independent=0.6921296296296297
- **sign-flip AUROC satisfies 1-AUC** (`PIK3CA/mTOR`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`PIK3CA/mTOR`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`PIK3CA/mTOR`) independent=smin_hi=0.813337 weaker_arm_hi=0.873808
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`PIK3CA/mTOR`) independent=[0.4702,0.8133]
- **matched assignment is D/A pocket B and D/B pocket A** (`PIK3CA/mTOR`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`PIK3CA/mTOR`) independent=4/1/0/0 of k=5
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`PIK3CA/mTOR`) independent=2.1333333333333333
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`PIK3CA/mTOR`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`AChE/BChE`) independent=0.6058201058201058
- **sign-flip AUROC satisfies 1-AUC** (`AChE/BChE`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`AChE/BChE`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`AChE/BChE`) independent=smin_hi=0.730321 weaker_arm_hi=0.750994
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`AChE/BChE`) independent=[0.4370,0.7303]
- **matched assignment is D/A pocket B and D/B pocket A** (`AChE/BChE`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`AChE/BChE`) independent=5/3/1/1 of k=10
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`AChE/BChE`) independent=1.7592592592592593
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`AChE/BChE`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`F2/F10`) independent=0.34475806451612906
- **sign-flip AUROC satisfies 1-AUC** (`F2/F10`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`F2/F10`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`F2/F10`) independent=smin_hi=0.477305 weaker_arm_hi=0.485571
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`F2/F10`) independent=[0.2109,0.4773]
- **matched assignment is D/A pocket B and D/B pocket A** (`F2/F10`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`F2/F10`) independent=4/1/6/0 of k=11
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`F2/F10`) independent=1.2551319648093844
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`F2/F10`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`PPARG/PPARA`) independent=0.6491935483870968
- **sign-flip AUROC satisfies 1-AUC** (`PPARG/PPARA`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`PPARG/PPARA`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`PPARG/PPARA`) independent=smin_hi=0.750838 weaker_arm_hi=0.778483
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`PPARG/PPARA`) independent=[0.5045,0.7508]
- **matched assignment is D/A pocket B and D/B pocket A** (`PPARG/PPARA`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`PPARG/PPARA`) independent=7/3/0/1 of k=11
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`PPARG/PPARA`) independent=2.1676136363636362
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`PPARG/PPARA`) independent=2000
- **summary_min equals min(D/A, D/B) point estimates** (`PPARA/PPARD`) independent=0.4462890625
- **sign-flip AUROC satisfies 1-AUC** (`PPARA/PPARD`)
- **Figure 2 plotted directional AUROCs match independent recompute** (`PPARA/PPARD`)
- **summary_min CI is not identical to either directional CI (min taken inside each replicate)** (`PPARA/PPARD`) independent=smin_hi=0.584130 weaker_arm_hi=0.585857
- **independent non-stratified bootstrap recovers locked summary_min CI within 0.005** (`PPARA/PPARD`) independent=[0.2958,0.5841]
- **matched assignment is D/A pocket B and D/B pocket A** (`PPARA/PPARD`)
- **Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id** (`PPARA/PPARD`) independent=5/2/3/1 of k=11
- **EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel** (`PPARA/PPARD`) independent=1.5625
- **Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw** (`PPARA/PPARD`) independent=2000
- **Figure 2D JAK1/TYK2 Top-10% 1/3/7/0 matches independent ranking** (`JAK1/TYK2`)
- **best descriptor is max of four unflipped summary_min values** (`EGFR/HER2`) independent=clogp
- **best descriptor is max of four unflipped summary_min values** (`AChE/BChE`) independent=tpsa
- **best descriptor is max of four unflipped summary_min values** (`PIK3CA/mTOR`) independent=heavy
- **physicochemical AUROCs are not max(AUC, 1-AUC); PIK3CA TPSA/cLogP can sit below 0.5**
- **GroupKFold splits keep each Bemis-Murcko scaffold in train or test, never both**
- **16 directional incremental ΔAUROC recomputed; max |Δ|=0.023** independent=0.023400000000000087
- **Figure 3B max |Δ| matches frozen incremental tables** independent=0.023400000000000087
- **primary ECFP4+docking logistic has no StandardScaler; fold-internal scaling is the Table S5 sensitivity only**
- **LogisticRegression is fit on training indices only; AUROC uses out-of-fold scores**
- **main-panel matched-minus-mismatched Δsummary_min** (`EGFR/HER2`) independent=0.1697
- **five-pair matched-minus-mismatched Δsummary_min** (`JAK1/JAK2`) independent=-0.01904296875
- **five-pair matched-minus-mismatched Δsummary_min** (`JAK1/TYK2`) independent=-0.06451612903225806
- **main-panel matched-minus-mismatched Δsummary_min** (`PIK3CA/mTOR`) independent=0.0902
- **main-panel matched-minus-mismatched Δsummary_min** (`AChE/BChE`) independent=0.1614
- **five-pair matched-minus-mismatched Δsummary_min** (`F2/F10`) independent=-0.030745967741935443
- **five-pair matched-minus-mismatched Δsummary_min** (`PPARG/PPARA`) independent=0.029564642137096753
- **five-pair matched-minus-mismatched Δsummary_min** (`PPARA/PPARD`) independent=0.01171875
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`AChE/BChE`) independent=-0.02499999999999991
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`PIK3CA/mTOR`) independent=-0.022499999999999964
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`JAK1/JAK2`) independent=0.008333333333333304
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`JAK1/TYK2`) independent=0.024999999999999967
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`F2/F10`) independent=-0.07894736842105265
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`PPARG/PPARA`) independent=0.0060526315789474205
- **holdout matched-minus-mismatched Δsummary_min matches locked CSV** (`PPARA/PPARD`) independent=0.15000000000000002
- **EGFR/HER2 has no unused-pool holdout, matching Figure 4 dagger** (`EGFR/HER2`)
- **seven holdout Δsummary_min values span manuscript -0.079 to +0.150** independent=-0.0789 to 0.1500
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`JAK1/JAK2`)
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`JAK1/TYK2`)
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`PIK3CA/mTOR`)
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`AChE/BChE`)
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`F2/F10`)
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`PPARG/PPARA`)
- **within-pair unused-pool holdout disjoint from main ChEMBL IDs** (`PPARA/PPARD`)
- **PM110 is a superset of main PM48 (expanded panel, not a disjoint holdout)** independent=PM48=48 PM110=115 extras=67
- **PIK3CA/mTOR holdout disjoint from PM48 and PM110**
- **EGFR main analysis set is the 110-ligand panel containing the original 40**
- **BindingDB ligands marked independent should not also be in scored-panel structures** independent=0
- **Figure 2A EGFR Δ 0.378 is in plotted_values and equal-score CSV** independent=0.378
- **Figure 2A JAK1/TYK2 Δ 0.444 is in plotted_values and equal-score CSV** independent=0.444

## WARNING

- **ligand-level source document IDs are not deposited for this pair**
  - pair: `JAK1/JAK2`
  - file: `document_blocked_ligand_groups_v1.csv covers only EGFR/AChE/PIK3CA`
  - independent: `110`
- **ligand-level source document IDs are not deposited for this pair**
  - pair: `JAK1/TYK2`
  - file: `document_blocked_ligand_groups_v1.csv covers only EGFR/AChE/PIK3CA`
  - independent: `109`
- **ligand-level source document IDs are not deposited for this pair**
  - pair: `F2/F10`
  - file: `document_blocked_ligand_groups_v1.csv covers only EGFR/AChE/PIK3CA`
  - independent: `107`
- **acyclic ligands have an empty Bemis-Murcko core; production GroupKFold uses singleton NONE:ligand groups**
  - pair: `PPARG/PPARA`
  - file: `analyze_five_pair_stack_v1.py:275-277`
  - independent: `3`
  - detail: `['PGPA_029', 'PGPA_097', 'PGPA_099']`
- **ligand-level source document IDs are not deposited for this pair**
  - pair: `PPARG/PPARA`
  - file: `document_blocked_ligand_groups_v1.csv covers only EGFR/AChE/PIK3CA`
  - independent: `109`
- **acyclic ligands have an empty Bemis-Murcko core; production GroupKFold uses singleton NONE:ligand groups**
  - pair: `PPARA/PPARD`
  - file: `analyze_five_pair_stack_v1.py:275-277`
  - independent: `1`
  - detail: `['PAPD_097']`
- **ligand-level source document IDs are not deposited for this pair**
  - pair: `PPARA/PPARD`
  - file: `document_blocked_ligand_groups_v1.csv covers only EGFR/AChE/PIK3CA`
  - independent: `110`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `JAK1/JAK2`
  - independent: `8`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `JAK1/TYK2`
  - independent: `8`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `PIK3CA/mTOR`
  - independent: `2`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `AChE/BChE`
  - independent: `9`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `F2/F10`
  - independent: `12`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `PPARG/PPARA`
  - independent: `15`
- **holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)**
  - pair: `PPARA/PPARD`
  - independent: `17`
- **holdout IDs appear in a different pair's main panel (per-pair leftover by design)**
  - independent: `4`
  - detail: `JAK1/TYK2→JAK1/JAK2:3; JAK1/JAK2→JAK1/TYK2:3; PPARG/PPARA→PPARA/PPARD:3; PPARA/PPARD→PPARG/PPARA:4`

## FAIL

- None. Independent ligand-level recomputation matched the manuscript, figures, and frozen CSVs on all blocking checks.

## Notes on items that are not FAILs

- Strict 6.5/5.5 is the **sampling** rule for six pairs. On those panels it coincides with θ=6.0 because gray-zone ligands were never sampled. Primary AUROCs were recomputed from pChEMBL with θ=6.0.
- PM110 includes all PM48 IDs; that is an expanded protocol-sensitivity panel, not an unused-pool holdout.
- Cross-pair JAK/PPAR ChEMBL sharing is expected: leftover holdouts exclude only that pair's main panel.
- Unused-pool holdouts are not scaffold-split versus the same pair's main panel. Shared scaffolds are expected leftover chemistry, not GroupKFold leakage.
- Ligand-level document IDs are deposited only for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR.
- `formulation_conventional_vs_directional_v1.csv` still contains empty summary_min CI cells and a note about separate-arm CIs. Table 2 uses `primary_directional_intervals_review_v1.csv`, whose intervals were recovered here by min-inside-replicate bootstrap.

