# Dual-Compatibility Scoring with Dual-VSDS Benchmark

## Motivation

Dual-target small-molecule design is usually assessed by docking the same ligand independently to targets A and B and fusing raw scores by mean, min, or rank average. Because docking scores are incommensurable across pockets and potency is often asymmetric, this practice systematically promotes A-only or B-only molecules into the dual-hit list. Existing dual-target docking papers fuse ranks for case studies but do not treat cross-target calibration, dual-vs-single hard negatives, fused/linked architecture conditioning, and leakage-controlled paired evaluation as a first-class machine-learning problem with an open benchmark. Generative polypharmacology papers design molecules but do not solve the dual-compatibility ranking decision.

## Method

1. Build Dual-VSDS: paired activity tables for multiple target pairs with labels dual / A-only / B-only / inactive (never treating untested as inactive), TrueNegative and RandomDecoy protocols, scaffold and leave-target-pair splits, plus curated design_type (fused/linked/merged) from literature.
2. Dock each ligand independently to A and B with a fixed engine (GNINA or Vina), optional RTMScore rescoring, and PoseBusters validity gating; keep top-K poses per target.
3. Calibrate each target’s scores to \(\hat p_t(\mathrm{active})\) on that target’s labeled actives/inactives (isotonic/Platt or a small calibrator).
4. Fuse with a shortfall-sensitive dual-compatibility score, e.g. softmin of \((\hat p_A-\theta_A,\hat p_B-\theta_B)\), optionally conditioned on design_type experts.
5. Train/rank with an explicit dual-vs-single objective so dual positives outrank A-only and B-only hard negatives; evaluate with dual-vs-single pairwise accuracy, PR-AUC, EF under both decoy regimes, and calibration ECE.
6. Hold out private NLRP3/JNK1 cell-active fused/linked molecules as a time-split external ranking test (cell as L2/L3 holdout, not binding gold standard); report PK failures as exposure decoupling, not as a PK predictor claim.
7. Release data schemas, docking YAML, fusion code, and the Dual-VSDS splits as the primary community resource.
