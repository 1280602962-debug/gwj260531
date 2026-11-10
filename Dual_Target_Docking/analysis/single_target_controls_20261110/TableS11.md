# Table S11 Exploratory single-target controls using frozen Vina scores

## A. Single-target-active versus neither

| Target pair | Score target | Positive class | n positive | n negative | AUROC (95% CI) |
| --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | EGFR | A-only | 37 | 12 | 0.5766 [0.3761, 0.7680] |
| EGFR/HER2 | HER2 | B-only | 32 | 12 | 0.8672 [0.7109, 0.9818] |
| JAK1/JAK2 | JAK1 | A-only | 32 | 14 | 0.6049 [0.4241, 0.7790] |
| JAK1/JAK2 | JAK2 | B-only | 32 | 14 | 0.5446 [0.3504, 0.7388] |
| JAK1/TYK2 | JAK1 | A-only | 32 | 14 | 0.7455 [0.5558, 0.9085] |
| JAK1/TYK2 | TYK2 | B-only | 32 | 14 | 0.7522 [0.5804, 0.8996] |
| PIK3CA/mTOR | PIK3CA | A-only | 32 | 14 | 0.5357 [0.3549, 0.7121] |
| PIK3CA/mTOR | mTOR | B-only | 32 | 14 | 0.6183 [0.4576, 0.7768] |
| AChE/BChE | AChE | A-only | 25 | 15 | 0.6213 [0.4427, 0.7867] |
| AChE/BChE | BChE | B-only | 28 | 15 | 0.6071 [0.4310, 0.7786] |
| F2/F10 | F2 | A-only | 32 | 14 | 0.6339 [0.4621, 0.7879] |
| F2/F10 | F10 | B-only | 32 | 14 | 0.6406 [0.4353, 0.8304] |
| PPARG/PPARA | PPARG | A-only | 32 | 14 | 0.6685 [0.4631, 0.8504] |
| PPARG/PPARA | PPARA | B-only | 32 | 14 | 0.6384 [0.4330, 0.8281] |
| PPARA/PPARD | PPARA | A-only | 32 | 14 | 0.5223 [0.3303, 0.7143] |
| PPARA/PPARD | PPARD | B-only | 32 | 14 | 0.7009 [0.5290, 0.8527] |

## B. Single-target-active versus opposite single-target-active

| Target pair | Score target | Positive class | n positive | n negative | AUROC (95% CI) |
| --- | --- | --- | --- | --- | --- |
| EGFR/HER2 | EGFR | A-only | 37 | 32 | 0.1807 [0.0845, 0.2914] |
| EGFR/HER2 | HER2 | B-only | 32 | 37 | 0.7939 [0.6748, 0.8970] |
| JAK1/JAK2 | JAK1 | A-only | 32 | 32 | 0.6250 [0.4844, 0.7588] |
| JAK1/JAK2 | JAK2 | B-only | 32 | 32 | 0.3755 [0.2417, 0.5186] |
| JAK1/TYK2 | JAK1 | A-only | 32 | 32 | 0.3125 [0.1865, 0.4512] |
| JAK1/TYK2 | TYK2 | B-only | 32 | 32 | 0.6777 [0.5410, 0.8057] |
| PIK3CA/mTOR | PIK3CA | A-only | 32 | 32 | 0.4268 [0.2900, 0.5713] |
| PIK3CA/mTOR | mTOR | B-only | 32 | 32 | 0.6348 [0.4941, 0.7686] |
| AChE/BChE | AChE | A-only | 25 | 28 | 0.7000 [0.5500, 0.8286] |
| AChE/BChE | BChE | B-only | 28 | 25 | 0.5293 [0.3629, 0.6900] |
| F2/F10 | F2 | A-only | 32 | 32 | 0.4619 [0.3213, 0.6084] |
| F2/F10 | F10 | B-only | 32 | 32 | 0.5039 [0.3584, 0.6445] |
| PPARG/PPARA | PPARG | A-only | 32 | 32 | 0.5044 [0.3584, 0.6475] |
| PPARG/PPARA | PPARA | B-only | 32 | 32 | 0.5083 [0.3623, 0.6523] |
| PPARA/PPARD | PPARA | A-only | 32 | 32 | 0.5542 [0.4102, 0.6958] |
| PPARA/PPARD | PPARD | B-only | 32 | 32 | 0.7144 [0.5801, 0.8345] |

S1: negative Vina docking energy, with higher values favorable (search seed 42). The positive class is the single-target class active at the scored target. Neither ligands are negative in panel A; the opposite single-target class is negative in panel B. Dual ligands are excluded. Counts include eligible observations with a valid score at that target. CI: 95% class-stratified percentile bootstrap, 10,000 resamples, seed 271828; linear 2.5th and 97.5th percentiles. Positive resamples are shared between the two comparisons. Score ties count as 0.5. Intervals are exploratory and unadjusted for multiplicity.
