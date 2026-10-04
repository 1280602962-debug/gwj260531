# STAGE2_EVIDENCE_NOTE

This note describes the Stage 2 tables written on 2026-10-04. It does not add a scientific estimate, a success threshold, or a revision of the paper's main conclusion. Figures under `results/jcim_stage2_figures/` plot the written PRIMARY and Stage 2 tables in the display order EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, PPARA/PPARD. CSV row order remains alphabetical. G, J, and K were not entered.

## Directions that move both ways

PRIMARY M0, on the shared-dual universe, is not uniformly above 0.5 in both directions. EGFR/HER2 dual-versus-B_only and F2/F10 dual-versus-A_only sit below 0.5; other arms sit above. The paired `summary_min` deltas for M1−M0, M1b−M1, M2−M0, and M3−M0 also change sign across pairs. PPARG/PPARA M2−M0 is a large negative point. These are the frozen paired contrasts, not a new ranking of methods.

Module B compares dual-versus-neither with dual-versus-single. Eleven of sixteen arms have a positive point (neither higher than single) and five have a negative point. Fourteen of the sixteen paired intervals contain zero. The two that do not are EGFR/HER2 D_vs_B (+0.500, interval +0.302 to +0.682) and JAK1/TYK2 D_vs_B (+0.384, interval +0.195 to +0.569). PIK3CA/mTOR keeps neither n=4 in both arms.

Module C is the wrong-pocket contrast on the same directional M0 population, with members dropped only when the opposite-pocket score is not finite. AChE D_vs_A has no such extra drop (25 to 25). Eleven points are positive, four are negative, and JAK1/JAK2 D_vs_A is exactly zero. Fourteen intervals contain zero. AChE/BChE D_vs_A (+0.166, +0.055 to +0.283) and PPARA/PPARD D_vs_A (+0.235, +0.069 to +0.404) have intervals above zero. **错口袋对照不能单独证明结合机制。** A positive pocket delta is a score contrast under the frozen pocket assignment. It does not identify a binding pose, a residue contact, or a causal mechanism.

Module D reports fixed-baseline point estimates on the chemistry OOF population (928 records). The ECFP4+M0 minus ECFP4 increment is negative on eleven arms, positive on four, and zero on PIK3CA/mTOR D_vs_B. The increments are small relative to the gap between ECFP4 and the physchem model on most arms. Supporting M0 is the same layer-3 population, not the full labeled PRIMARY population, so it is not subtracted from the primary M0. No chemistry interval was added. The increment is a point estimate on one frozen split, not a demonstrated stable gain.

## Sensitivity, not a second primary result

Module E keeps five seeds (17, 29, 42, 71, 101). All 24 pair-level summaries had five finite points, so the linear IQR was written; none fell back to NA. Seed 42 matches formal M0. AB_046 seed 17 on 4BDS remains TIMEOUT and is listed, not filled in. The figure shows the stored seed points and annotates the stored IQR. That span is not a confidence interval. **结果变化可能同时来自分数和有效成员变化。** This note does not assign each pair's seed-to-seed movement to only one of those sources.

Module F replaces only the specified side on the common complete population. The unreplaced AUROC delta is 0.0 on all eight rows, with a zero interval. That is a program-consistency check of the copy-and-replace construction. The affected-side deltas move in both directions. Seven of eight affected intervals contain zero. PIK3CA/mTOR with ALT_4JT5 on side B is the exception (+0.123, interval +0.004 to +0.266). Unreplaced exact zero is not evidence that the alternate receptor is equivalent.

## Joint list

Module I is a constructed challenge-set description on dual + A_only + B_only. Top-k uses k=ceil(0.10×n), average midrank, worst rank, and the two frozen tie-breaks. M0 Top-k sizes are k=10, 10, 10, 5, 8, 10, 10, 10 in display order. Class composition inside those lists varies by pair. The four concordance states all occur: M1 is 2/2/2/2 across improve, worsen, discordant, and no change; M1b is 3/3/1/1; M2 is 4/3/1/0; M3 is 3/2/2/1. M1b is a supporting sum of PRIMARY points on a membership-matched population. It has no new interval. The list is not a natural-library hit rate, not an independent validation, and not a causal effect.

## Limits that stay attached to every number

Activity labels come from the frozen source tables. Endpoint definitions and label aggregation differ across the eight pairs, and this round does not re-derive them. AUROC ties use 0.5. Intervals that contain zero, wide intervals, and negative points remain visible. No pair is marked successful or unsuccessful by a new cutoff. Source, endpoint, and aggregation limits are not removed by a high ECFP4 point or by a pocket interval that excludes zero.
