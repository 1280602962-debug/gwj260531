# HOLDOUT_WRONG_POCKET_POTENCY_VERDICT_V1

Zero new docking. Matching copied from T0 (`|Δp|≤0.5` potency; `|Δheavy|≤2` size).
Pocket-matched: D/A uses vina_B; D/B uses vina_A. Wrong-pocket: the reverse.

## Does unused-pool sampling shift potency vs the main panel?

| pair | class | feature | holdout mean (n) | main mean (n) | Δ |
|------|-------|---------|-----------------:|--------------:|--:|
| AChE/BChE | dual | pA | 7.554 (20) | 7.713 (27) | -0.159 |
| AChE/BChE | dual | pB | 7.164 (20) | 7.622 (27) | -0.458 |
| AChE/BChE | dual | heavy | 35.1 (20) | 33.444 (27) | 1.656 |
| AChE/BChE | A_only | pA | 7.555 (20) | 7.31 (25) | 0.245 |
| AChE/BChE | A_only | pB | 4.989 (20) | 4.873 (25) | 0.116 |
| AChE/BChE | A_only | heavy | 33.95 (20) | 31.12 (25) | 2.83 |
| AChE/BChE | B_only | pA | 5.045 (20) | 5.048 (28) | -0.003 |
| AChE/BChE | B_only | pB | 7.223 (20) | 7.376 (28) | -0.153 |
| AChE/BChE | B_only | heavy | 29.45 (20) | 31.964 (28) | -2.514 |
| PIK3CA/PIK3CB | dual | pA | 8.35 (20) | 7.745 (28) | 0.605 |
| PIK3CA/PIK3CB | dual | pB | 7.544 (20) | 7.745 (28) | -0.201 |
| PIK3CA/PIK3CB | dual | heavy | 34.5 (20) | 33.464 (28) | 1.036 |
| PIK3CA/PIK3CB | A_only | pA | 7.276 (19) | 7.093 (27) | 0.183 |
| PIK3CA/PIK3CB | A_only | pB | 5.234 (19) | 5.267 (27) | -0.033 |
| PIK3CA/PIK3CB | A_only | heavy | 31.579 (19) | 30.667 (27) | 0.912 |
| PIK3CA/PIK3CB | B_only | pA | 5.131 (20) | 5.158 (28) | -0.027 |
| PIK3CA/PIK3CB | B_only | pB | 7.204 (20) | 7.349 (28) | -0.145 |
| PIK3CA/PIK3CB | B_only | heavy | 28.3 (20) | 29.893 (28) | -1.593 |
| PIK3CA/mTOR | dual | pA | 7.567 (20) | 8.639 (18) | -1.072 |
| PIK3CA/mTOR | dual | pB | 8.061 (20) | 8.399 (18) | -0.338 |
| PIK3CA/mTOR | dual | heavy | 33.5 (20) | 31.667 (18) | 1.833 |
| PIK3CA/mTOR | A_only | pA | 7.472 (20) | 8.731 (14) | -1.259 |
| PIK3CA/mTOR | A_only | pB | 5.129 (20) | 5.424 (14) | -0.295 |
| PIK3CA/mTOR | A_only | heavy | 32.25 (20) | 32.214 (14) | 0.036 |
| PIK3CA/mTOR | B_only | pA | 5.202 (20) | 5.629 (12) | -0.427 |
| PIK3CA/mTOR | B_only | pB | 7.286 (20) | 9.046 (12) | -1.76 |
| PIK3CA/mTOR | B_only | heavy | 30.95 (20) | 32.5 (12) | -1.55 |

## Matched-subset directional AUROC (holdout)

| pair | family | aggregation | D/A | D/B | summary_min | n_min | underpowered |
|------|--------|-------------|----:|----:|------------:|------:|:------------:|
| AChE/BChE | potency_matched | pocket_matched | 0.6288 | 0.5926 | 0.5926 | 18 | N |
| AChE/BChE | potency_matched | wrong_pocket | 0.6454 | 0.642 | 0.642 | 18 | N |
| AChE/BChE | size_matched | pocket_matched | 0.6923 | 0.4074 | 0.4074 | 9 | N |
| AChE/BChE | size_matched | wrong_pocket | 0.7041 | 0.4321 | 0.4321 | 9 | N |
| AChE/BChE | unmatched | pocket_matched | 0.635 | 0.6175 | 0.6175 | 20 | N |
| AChE/BChE | unmatched | wrong_pocket | 0.6425 | 0.6525 | 0.6425 | 20 | N |
| PIK3CA/PIK3CB | potency_matched | pocket_matched | 0.9008 | 0.3633 | 0.3633 | 11 | N |
| PIK3CA/PIK3CB | potency_matched | wrong_pocket | 0.7934 | 0.5625 | 0.5625 | 11 | N |
| PIK3CA/PIK3CB | size_matched | pocket_matched | 0.6786 | 0.3018 | 0.3018 | 13 | N |
| PIK3CA/PIK3CB | size_matched | wrong_pocket | 0.5969 | 0.426 | 0.426 | 13 | N |
| PIK3CA/PIK3CB | unmatched | pocket_matched | 0.7658 | 0.425 | 0.425 | 19 | N |
| PIK3CA/PIK3CB | unmatched | wrong_pocket | 0.6395 | 0.52 | 0.52 | 19 | N |
| PIK3CA/mTOR | potency_matched | pocket_matched | 0.832 | 0.7153 | 0.7153 | 12 | N |
| PIK3CA/mTOR | potency_matched | wrong_pocket | 0.7344 | 0.7361 | 0.7344 | 12 | N |
| PIK3CA/mTOR | size_matched | pocket_matched | 0.8756 | 0.7153 | 0.7153 | 12 | N |
| PIK3CA/mTOR | size_matched | wrong_pocket | 0.8178 | 0.8542 | 0.8178 | 12 | N |
| PIK3CA/mTOR | unmatched | pocket_matched | 0.86 | 0.765 | 0.765 | 20 | N |
| PIK3CA/mTOR | unmatched | wrong_pocket | 0.7875 | 0.8575 | 0.7875 | 20 | N |

## Does matching flip wrong-pocket ≥ matched?

- **AChE/BChE / unmatched**: pocket-matched 0.6175 vs wrong-pocket 0.6425 (gap matched−wrong = -0.025); wrong ≥ matched: **yes**.
- **AChE/BChE / potency_matched**: pocket-matched 0.5926 vs wrong-pocket 0.642 (gap matched−wrong = -0.049); wrong ≥ matched: **yes**.
- **AChE/BChE / size_matched**: pocket-matched 0.4074 vs wrong-pocket 0.4321 (gap matched−wrong = -0.025); wrong ≥ matched: **yes**.
- **PIK3CA/PIK3CB / unmatched**: pocket-matched 0.425 vs wrong-pocket 0.52 (gap matched−wrong = -0.095); wrong ≥ matched: **yes**.
- **PIK3CA/PIK3CB / potency_matched**: pocket-matched 0.3633 vs wrong-pocket 0.5625 (gap matched−wrong = -0.199); wrong ≥ matched: **yes**.
- **PIK3CA/PIK3CB / size_matched**: pocket-matched 0.3018 vs wrong-pocket 0.426 (gap matched−wrong = -0.124); wrong ≥ matched: **yes**.
- **PIK3CA/mTOR / unmatched**: pocket-matched 0.765 vs wrong-pocket 0.7875 (gap matched−wrong = -0.022); wrong ≥ matched: **yes**.
- **PIK3CA/mTOR / potency_matched**: pocket-matched 0.7153 vs wrong-pocket 0.7344 (gap matched−wrong = -0.019); wrong ≥ matched: **yes**.
- **PIK3CA/mTOR / size_matched**: pocket-matched 0.7153 vs wrong-pocket 0.8178 (gap matched−wrong = -0.102); wrong ≥ matched: **yes**.

### One-line verdict

**Holdout wrong-pocket ≥ pocket-matched survives potency matching on all three pairs.** Unused-pool sampling is not a sufficient explanation of the paradox. This is a diagnostic, not a protocol change: primary holdout numbers stay unmatched.

## What this is not

- Not a new docking run and not a replacement of Table S8 / HOLDOUT_VERDICT.
- Matched subsets are small (n≈20 before matching); n_min<8 cells are underpowered.
- Does not claim the main-panel pocket-matched > wrong-pocket gap is explained.

```bash
python3 Dual_Target_Docking/data/jcim_holdout_v0/scripts/wrong_pocket_potency_match_v1.py
```

