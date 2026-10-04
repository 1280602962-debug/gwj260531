# JOINT_RANKING_FREEZE

Status: TYPE B rules specified before any joint-ranking result is computed.
Inherited TYPE A pieces are marked explicitly.

## Name

Formal name: 双靶联合排序挑战 (dual-target joint-ranking challenge).

Forbidden names for this analysis:

- 真实虚拟筛选模拟
- natural-library virtual screening
- external validation of AUROC
- independent replication

## Scientific question

In an experimentally labeled, selectivity-enriched challenge set, when a candidate must rank well on both targets, can docking/scoring methods place dual-active ligands toward the top and reduce A-only/B-only entry into the top list?

This translates directional AUROC into candidate-list language. It does not predict hit rates in an unlabeled commercial library.

## Population

Main pool (TYPE B):

- `dual` + `A_only` + `B_only`
- `activity_eligible=1`
- `neither` is excluded from the main ranking challenge

This is a **selectivity-enriched challenge population**, not a natural chemical-library class mix.

`neither` remains available to Module B and is not used here.

## Completeness

Method-specific descriptive ranking (not for causal method comparison):

- include a ligand only if that method has finite `score_A` and finite `score_B`.

Pairwise common population for method comparison (TYPE A inheritance of `pairwise_common_population: required`):

- M0 vs M1, M0 vs M1b, M0 vs M2, M0 vs M3 each use the ligands for which **both** methods have finite score_A and finite score_B.
- Do not compare Top-10 composition across different ligand sets.

All-five-method intersection is `secondary_only` (TYPE A). It may be shown as a supporting side-by-side table. It must not replace the four pairwise-vs-M0 comparisons.

## Single-target ranks (TYPE B)

All formal scores are higher-better.

For a fixed population and method:

```
rank_A = rankdata(-score_A, method="average")
rank_B = rankdata(-score_B, method="average")
```

`rank 1` is best. Exact score ties receive the same midrank.

`global_ligand_entity_id` must not change the scientific single-target rank.

Scope note: `SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml` already has TYPE A `tiebreak: [score_desc, global_ligand_entity_id_lex]`. That rule applies only if a later analysis ranks a **single** score list into Top-K by raw score. It does **not** apply to the scientific ranks inside this joint-ranking challenge.

## Joint rank (TYPE B)

```
worst_rank = max(rank_A, rank_B)
```

Smaller `worst_rank` is better. A ligand rises in the joint list only if both pockets rank it well.

## Membership tie-breaks for the exact Top-k set (TYPE B)

If two ligands have the same `worst_rank`, prefer the smaller `rank_A + rank_B`.

If still tied, prefer `global_ligand_entity_id` lexical ascending.

The ID rule only decides which ligands occupy the fixed k slots. It is not a scientific rank.

## Top-K (TYPE A)

From `SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml`:

- `top_fraction = 0.10`
- `k = ceil(0.10 * n)`

Only Top 10% is used. Forbidden: Top 1%, Top 5%, Top 20%, Top 50, or trying several cutoffs.

`n` is the size of the population actually ranked in that table (method-specific complete, pairwise common, or secondary intersection).

## Required outputs

For each pair × method/population:

- `n_total`
- `n_dual_total`
- `n_A_only_total`
- `n_B_only_total`
- `k_top10`
- `top10_dual`
- `top10_A_only`
- `top10_B_only`
- `top10_dual_fraction`
- `top10_single_target_fraction = (top10_A_only + top10_B_only) / k_top10`
- `dual_recall_at_top10 = top10_dual / n_dual_total` (if `n_dual_total=0`, write NA)

## Random expectation (TYPE B)

Closed-form references only:

- `expected_dual_count = k * N_dual / N_total`
- `expected_dual_fraction = N_dual / N_total`

Forbidden: hypergeometric p-values, permutation tests, 10000 random shuffles, multiple-testing correction.

## Link to directional benchmark

Do not treat 8×5=40 points as independent experiments.

Within each pair, for M0→M1, M0→M1b, M0→M2, M0→M3, compare:

- `delta_summary_min = summary_min(new) − summary_min(M0)` using the already frozen pairwise common directional rule;
- `delta_top10_single_target_fraction = top10_single_target_fraction(new) − top10_single_target_fraction(M0)` on the matching pairwise common ranking population.

Expected descriptive direction: positive `delta_summary_min` with negative `delta_top10_single_target_fraction`.

Concordance coding (TYPE B, frozen before ranking results):

- `CONCORDANT` if `delta_summary_min > 0` and `delta_top10_single_target_fraction < 0`;
- `NO_CHANGE` if either delta is exactly 0 (zero is not discordant);
- `DISCORDANT` otherwise.

M1b−M0 supporting rule: first verify that the four-sided common-complete membership of M1b−M0 equals that of M1−M0 for every pair. Only after that identity check may a supporting point estimate be added. Do not add a new bootstrap CI for this supporting comparison. Do not write it into `PRIMARY_DIRECTIONAL_METRICS.csv` or `PRIMARY_METHOD_DELTA_METRICS.csv`.

Allowed later outputs: paired/linked plot and concordant / discordant / no-change counts.

Forbidden: Spearman or Pearson correlation; calling this an independent validation of AUROC.

## Interpretation boundary

Top-10 fractions describe single-target occupancy in this constructed selectivity-enriched set only. They must not be written as the false-positive rate of a real screening library.
