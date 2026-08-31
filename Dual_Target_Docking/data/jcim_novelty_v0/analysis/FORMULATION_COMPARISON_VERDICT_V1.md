# Formulation / chemotype / incremental verdict (frozen Vina, no new docking)

Numbers from `data/jcim_novelty_v0/tables/` after `scripts/benchmark_formulation_v1.py` (bootstrap B = 2000, seed 20260729; contrast CIs use a stable md5 offset). Table 2 directional arms are reproduced to 0.0001.

## 1. Conventional vs DualFourClass

Do **not** claim a four-pair overestimation law. EGFR/HER2 is the only pair where Dual-vs-neither looks like a successful dual-target docking result while directional `summary_min` does not. Dual-versus-neither is a **nonselectivity-controlled comparator**, not “the conventional dual-target benchmark.” The contrast is descriptive (different negatives), not a paired Δ.

| Pair | summary_min | Dual vs neither mean [95% CI] | n_neither | Dual vs all non-duals | Single-target pocket A / B |
|---|---:|---:|---:|---:|---|
| EGFR/HER2 | 0.4297 | 0.7560 [0.5625, 0.9197] | 12 | 0.5514 [0.4429, 0.6664] | 0.4254 / 0.7593 |
| AChE/BChE | 0.6058 | 0.6494 [0.4864, 0.8049] | 15 | 0.5792 | 0.6395 / 0.6368 |
| PIK3CA/PIK3CB | 0.5000 | 0.5592 [0.3638, 0.7455] | 16 | 0.5558 | 0.4804 / 0.6551 |
| PIK3CA/mTOR | 0.6921 | 0.5139 [0.2222, 0.8056] | **4 underpowered** | 0.6741 [0.5130, 0.8222] | 0.5244 / 0.6519 |

The Dual-vs-neither and directional CIs use different negatives; they are not a paired Δ. EGFR intervals barely overlap.

EGFR single-target analogue is itself asymmetric: pocket B (HER2) inhibitor vs noninhibitor AUROC 0.759, pocket A (EGFR) 0.425. That is compatible with Table 2’s weak Dual vs B-only arm (pocket A).

## 2. Chemotype-matched hard negatives

Nearest-neighbor ECFP4 Tanimoto ≥ 0.7 matching was empty. Cuts 0.3 / 0.4 / 0.5 were used.

Strongest drop: PIK3CA/PIK3CB Dual vs A-only 0.6905 → T ≥ 0.3 n_neg = 11 AUROC **0.5032**; distant T < 0.3 n_neg = 16 AUROC **0.8192**.

High-T cells with n_neg ≤ 3 are listed in the CSV and are not interpreted.

## 3. Incremental information

`GroupKFold` logistic CV AUROC. Docking feature = the pocket-matched score for that arm. Rank docking AUROC is reported alongside and is **not** the logistic number.

| Pair | Arm | ECFP4 | docking (logistic) | ECFP4+docking | rank docking |
|---|---|---:|---:|---:|---:|
| EGFR/HER2 | D vs A | 0.7453 | 0.5912 | 0.7509 | 0.6664 |
| EGFR/HER2 | D vs B | 0.8895 | 0.4252 | 0.8873 | 0.4297 |
| AChE/BChE | D vs A | 0.8948 | 0.6119 | 0.8933 | 0.6504 |
| AChE/BChE | D vs B | 0.8214 | 0.5291 | 0.8082 | 0.6058 |
| PIK3CA/PIK3CB | D vs A | 0.7817 | 0.6534 | 0.7857 | 0.6905 |
| PIK3CA/PIK3CB | D vs B | 0.7691 | 0.2449 | 0.7717 | 0.5000 |
| PIK3CA/mTOR | D vs A | 0.7619 | 0.6270 | 0.7421 | 0.7143 |
| PIK3CA/mTOR | D vs B | 0.8889 | 0.6343 | 0.8981 | 0.6921 |

Docking incremental |Δ| vs ECFP4 is at most 0.0198 and is negative on several arms. PM n_scaffolds ≈ n, so CV is nearly leave-one-scaffold.

## 4. Mixed-library enrichment

EGFR/HER2 `vina_mean` Top-10: 1 dual, 9 selectives, EF = 0.393, hard-neg fraction 0.90. EF5 = 0.655 (worse than random). This is the screening-facing form of the EGFR formulation gap.
