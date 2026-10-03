# FORBIDDEN_POSTHOC_ACTIONS

These actions remain forbidden after this freeze, including during later compute phases unless a new freeze explicitly replaces a listed item.

1. Do not add a scoring function to chase a higher AUROC.
2. Do not add a receptor to chase a better result.
3. Do not select the best seed.
4. Do not drop a failing pair (including EGFR/HER2, AChE/BChE, PPARG/PPARA).
5. Do not change primary θ = 6.0.
6. Do not change the frozen Top 10% rule.
7. Do not add Top 20% or other cutoffs and then keep the prettier one.
8. Do not change the minimax / midrank / rank-sum / ID tie formula after seeing ranking results.
9. Do not change the fingerprint after seeing chemistry results.
10. Do not search chemistry-model hyperparameters.
11. Do not drop holdout ligands for high similarity after seeing overlap.
12. Do not rebalance or resample holdout using main-set results.
13. Do not impute missing scores (0, 99, median, worst, best, or manual fill).
14. Do not re-query ChEMBL or re-adjudicate primary labels.
15. Do not invent a consensus score after seeing results.
16. Do not call untested new-library candidates dual inhibitors.
17. Do not treat MD as activity validation.
18. Do not call the joint-ranking challenge a real virtual screen of a natural library.
19. Do not call same-set ranking results an independent validation of AUROC.
20. Do not describe TYPE B new rules as historical preregistration.
21. Do not promote M1, M1b, M2, or M3 to a new primary method after seeing results.
22. Do not average five seed scores into a new official score.
23. Do not replace a primary receptor with an alternative receptor.
24. Do not use `max_saved_RTMScore` as formal M2.
25. Do not compute `AUROC(min(score_A, score_B))`.
26. Do not use independent CI overlap in place of paired delta bootstrap.
27. Do not restore `analysis_freeze.yaml` as method authority.
28. Do not use `/tmp/pr39_fiveseed` as an execution or input root.
29. Do not add EF1%, EF5%, BEDROC, MCC, F1, accuracy, or PR-AUC unless a later freeze states a new question.
30. Do not add a new compute module because “the analysis would be more complete”.
