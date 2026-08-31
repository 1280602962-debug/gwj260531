# Four-pair multi-seed Vina sensitivity

Source scores: `data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv`
Frozen seeds: 20260727 (primary, reused) + 20260811–20260814.
Protocol otherwise identical to production (receptors, boxes, exhaustiveness, modes, energy_range).

## Per-seed metrics

| pair | seed | dual_vs_A | dual_vs_B | summary_min | dual_vs_neither | gap |
|---|---:|---:|---:|---:|---:|---:|
| AChE/BChE | 20260727 | 0.6504 | 0.6058 | 0.6058 | 0.6494 | 0.0436 |
| AChE/BChE | 20260811 | 0.5988 | 0.6396 | 0.5988 | 0.659 | 0.0602 |
| AChE/BChE | 20260812 | 0.5645 | 0.6104 | 0.5645 | 0.6373 | 0.0728 |
| AChE/BChE | 20260813 | 0.6017 | 0.6104 | 0.6017 | 0.6387 | 0.0369 |
| AChE/BChE | 20260814 | 0.5527 | 0.6119 | 0.5527 | 0.6507 | 0.0979 |
| EGFR/HER2 | 20260727 | 0.6664 | 0.4297 | 0.4297 | 0.7641 | 0.3344 |
| EGFR/HER2 | 20260811 | 0.6598 | 0.3728 | 0.3728 | 0.7798 | 0.407 |
| EGFR/HER2 | 20260812 | 0.6391 | 0.3214 | 0.3214 | 0.7634 | 0.442 |
| EGFR/HER2 | 20260813 | 0.6532 | 0.3672 | 0.3672 | 0.7693 | 0.4022 |
| EGFR/HER2 | 20260814 | 0.6692 | 0.3956 | 0.3956 | 0.7768 | 0.3811 |
| PIK3CA/PIK3CB | 20260727 | 0.6905 | 0.5 | 0.5 | 0.5798 | 0.0798 |
| PIK3CA/PIK3CB | 20260811 | 0.6587 | 0.5019 | 0.5019 | 0.5926 | 0.0907 |
| PIK3CA/PIK3CB | 20260812 | 0.6634 | 0.4783 | 0.4783 | 0.5977 | 0.1193 |
| PIK3CA/PIK3CB | 20260813 | 0.6541 | 0.4681 | 0.4681 | 0.5915 | 0.1234 |
| PIK3CA/PIK3CB | 20260814 | 0.6772 | 0.4707 | 0.4707 | 0.5993 | 0.1287 |
| PIK3CA/mTOR | 20260727 | 0.7143 | 0.6921 | 0.6921 | 0.5278 | -0.1644 |
| PIK3CA/mTOR | 20260811 | 0.746 | 0.713 | 0.713 | 0.5486 | -0.1644 |
| PIK3CA/mTOR | 20260812 | 0.7262 | 0.7778 | 0.7262 | 0.5486 | -0.1776 |
| PIK3CA/mTOR | 20260813 | 0.7381 | 0.6759 | 0.6759 | 0.5347 | -0.1412 |
| PIK3CA/mTOR | 20260814 | 0.744 | 0.7037 | 0.7037 | 0.5347 | -0.169 |

## Consistency vs primary seed

- **AChE/BChE**: summary_min median 0.5988 (range 0.0531); gap-sign match 5/5; neither>summary_min order match 5/5.
- **EGFR/HER2**: summary_min median 0.3728 (range 0.1083); gap-sign match 5/5; neither>summary_min order match 5/5.
- **PIK3CA/PIK3CB**: summary_min median 0.4783 (range 0.0338); gap-sign match 5/5; neither>summary_min order match 5/5.
- **PIK3CA/mTOR**: summary_min median 0.7037 (range 0.0503); gap-sign match 5/5; neither>summary_min order match 5/5.

## Claim ceiling

- Allowed: report median/IQR/range across frozen seeds; state whether the primary qualitative pattern held.
- Forbidden: picking a favorable seed; replacing primary Table 2 with a multi-seed mean; claiming seed robustness beyond these four pairs.

