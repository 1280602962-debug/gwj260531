"""Post hoc single-target discrimination check (M0 / Vina).

Computed after the V4.2 result freeze; not part of the pre-specified analysis plan.
Reads the frozen primary population table and asks, for each of the 16 target arms,
whether the target's own Vina score ranks its measured actives above measured
inactives drawn from three different negative classes:

  pooled          positives = dual + own-target-only, negatives = partner-only + neither
  vs_neither      positives = own-target-only,        negatives = neither
  vs_other_only   positives = own-target-only,        negatives = partner-only

The pooled contrast shares ligand pairs with the primary directional contrast
(dual vs partner-only), so its agreement with the directional result is partly
arithmetic. The other two contrasts share no positives with it.

Bootstrap settings match the primary analysis: 10,000 stratified resamples, seed 271828.
"""

import numpy as np
import pandas as pd

# Frozen population table, carried on branch cursor/v4-2-final-audited-results-20260930.
POPULATION = "results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
N_BOOT = 10000
SEED = 271828


def auroc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    ranks = pd.Series(np.concatenate([pos, neg])).rank().values
    return (ranks[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def boot_ci(pos, neg, rng):
    reps = [auroc(rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)) for _ in range(N_BOOT)]
    return np.percentile(reps, [2.5, 97.5])


def main():
    pop = pd.read_csv(POPULATION)
    pop = pop[(pop.activity_eligible == 1) & pop["class"].notna()]
    rng = np.random.default_rng(SEED)

    rows = []
    for pair, group in pop.groupby("pair"):
        target_a, target_b = pair.split("/")
        arms = [
            (target_a, "M0_score_A", "M0_A_valid", "A_only", "B_only"),
            (target_b, "M0_score_B", "M0_B_valid", "B_only", "A_only"),
        ]
        for target, score_col, valid_col, own_only, partner_only in arms:
            scored = group[(group[valid_col] == 1) & group[score_col].notna()]
            by_class = {k: v[score_col].values for k, v in scored.groupby("class")}
            contrasts = {
                "pooled": (
                    np.concatenate([by_class.get("dual", []), by_class.get(own_only, [])]),
                    np.concatenate([by_class.get(partner_only, []), by_class.get("neither", [])]),
                ),
                "vs_neither": (by_class.get(own_only, []), by_class.get("neither", [])),
                "vs_other_only": (by_class.get(own_only, []), by_class.get(partner_only, [])),
            }
            for name, (pos, neg) in contrasts.items():
                lo, hi = boot_ci(pos, neg, rng)
                rows.append(
                    dict(
                        pair=pair,
                        target=target,
                        contrast=name,
                        n_pos=len(pos),
                        n_neg=len(neg),
                        auroc=round(auroc(pos, neg), 4),
                        ci_lo=round(lo, 4),
                        ci_hi=round(hi, 4),
                    )
                )

    result = pd.DataFrame(rows)
    result.to_csv("results/post_hoc_single_target/POST_HOC_SINGLE_TARGET_AUROC.csv", index=False)
    for name, block in result.groupby("contrast"):
        print(f"=== {name} ===")
        print(block.to_string(index=False))
        print(
            f"point>0.5: {(block.auroc > 0.5).sum()}/16   "
            f"CI entirely >0.5: {(block.ci_lo > 0.5).sum()}   "
            f"CI entirely <0.5: {(block.ci_hi < 0.5).sum()}\n"
        )


if __name__ == "__main__":
    main()
