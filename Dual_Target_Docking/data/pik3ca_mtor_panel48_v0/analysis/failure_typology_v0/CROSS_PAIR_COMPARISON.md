# Cross-pair comparison — EGFR/HER2 vs PIK3CA/mTOR

| Dimension | EGFR/HER2 (panel40) | PIK3CA/mTOR (panel48) |
|-----------|---------------------|------------------------|
| Structures | 3POZ / 3RCD (TAK-285) | 4L23 / 4JT6 (PI-103/X6K) |
| Cognate exhaustiveness | E=8 sufficient | **E=16 required** |
| vina_mean AUROC (Dual vs rest/hardneg) | ~0.55 | ~0.63 |
| Best arm | rtm_min_z ~0.71 | rtm_min_z ~0.69 |
| Top10 hardneg trend after RTM | ↓ (6→3) | **mixed: B↓ A↑** |
| Rescued hardneg | EH40_18-like | WYE-132 (PM48_34) |
| Stubborn hardneg | EH40_23 | PM48_26/20/21 |
| Injured true duals | not central | Torin1, Omipalisib |
| Clash gate on stubborn | fail (0 clashes) | fail (0 clashes) |
| Simple shortfall on stubborn | limited | **pre-test negative** (PM48_26 stays #1) |
| C4 claim | partial success narrative | **not closed** — ruler needs diagnostics |

## Shared protocol lessons
1. Vina Top-1 mean is insufficient on both pairs.
2. Best-of-9 + RTM changes ranks materially (RTM-best ≠ mode1 often).
3. Geometric clash gating does not solve chemotype-homolog false duals.
4. Failure typology transfers (T1/T2); new T5 appears on PIK3CA/mTOR.

## Honest one-liner for paper
Rescoring and weak-end fusion are **necessary but not sufficient**; the transferable product is a **reproducible dual-target decision protocol with typed failure modes**, not a universal AUROC win.
