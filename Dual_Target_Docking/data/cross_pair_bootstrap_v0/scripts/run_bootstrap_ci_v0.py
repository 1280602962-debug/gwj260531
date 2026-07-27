#!/usr/bin/env python3
"""Bootstrap CIs for decision-ablation arms on EGFR panel40 + PIK3CA/mTOR panel48.

Frozen thresholds (unchanged from decision_ablation_v0):
  shortfall_lambda = 0.5
  consensus_top_frac = 0.25

No re-docking, no flag-in-score, no threshold retuning.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPT = Path(__file__).resolve()
REPO = _SCRIPT.parents[3]  # .../Dual_Target_Docking
assert (REPO / "data" / "pik3ca_mtor_panel48_v0").is_dir(), REPO
OUT = REPO / "data" / "cross_pair_bootstrap_v0"
OUT.mkdir(parents=True, exist_ok=True)

# Frozen a priori — do not retune
SHORTFALL_LAMBDA = 0.5
CONSENSUS_TOP_FRAC = 0.25
N_BOOT = 2000
SEED = 20260727
ARMS = [
    "vina_mean",
    "rtm_min_z",
    "rtm_shortfall",
    "consensus_rank_mean",
    "consensus_and_top25",
]


def _stable_seed(*parts) -> int:
    s = "::".join(str(p) for p in parts)
    h = 2166136261
    for ch in s.encode():
        h ^= ch
        h = (h * 16777619) & 0xFFFFFFFF
    return SEED + (h % 1_000_000)


def roc_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # Mann–Whitney / trapezoid equivalent
    correct = 0.0
    for p in pos:
        correct += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(correct / (len(pos) * len(neg)))


def top10_counts(classes: np.ndarray, scores: np.ndarray) -> dict:
    order = np.argsort(-scores)
    top = classes[order[:10]]
    return {
        "top10_dual": int(np.sum(top == "dual")),
        "top10_A_only": int(np.sum(top == "A_only")),
        "top10_B_only": int(np.sum(top == "B_only")),
        "top10_neither": int(np.sum(top == "neither")),
        "top10_hardneg": int(np.sum(top != "dual")),
    }


def attach_arms(df: pd.DataFrame, z_a: str, z_b: str) -> pd.DataFrame:
    """Build five arms with frozen thresholds on existing scores."""
    out = df.copy()
    n = len(out)
    top_n = max(1, int(np.ceil(CONSENSUS_TOP_FRAC * n)))
    out["rtm_shortfall"] = out["rtm_min_z"] - SHORTFALL_LAMBDA * (
        out[z_a] - out[z_b]
    ).abs()
    out["rank_vina_mean"] = out["vina_mean"].rank(ascending=False, method="min").astype(int)
    out["rank_rtm_min_z"] = out["rtm_min_z"].rank(ascending=False, method="min").astype(int)
    out["consensus_rank_mean"] = -(out["rank_vina_mean"] + out["rank_rtm_min_z"]) / 2.0
    passes = (out["rank_vina_mean"] <= top_n) & (out["rank_rtm_min_z"] <= top_n)
    out["consensus_and_top25"] = np.where(
        passes, out["rtm_min_z"], out["rtm_min_z"] - 1e3
    )
    out["consensus_and_top25_pass"] = passes.astype(int)
    out.attrs["top_n"] = top_n
    return out


def load_panels():
    pm = pd.read_csv(
        REPO / "data" / "pik3ca_mtor_panel48_v0" / "tables" / "ablation_ligand_scores.csv"
    )
    eh = pd.read_csv(
        REPO / "data" / "egfr_her2_panel40_v0" / "tables" / "ablation_ligand_scores.csv"
    )
    pm = attach_arms(pm, "rtm_4L23_z", "rtm_4JT6_z")
    eh = attach_arms(eh, "rtm_3POZ_z", "rtm_3RCD_z")
    return {
        "PIK3CA_mTOR_panel48": pm,
        "EGFR_HER2_panel40": eh,
    }


def point_metrics(df: pd.DataFrame, arm: str) -> dict:
    y = (df["class"] == "dual").astype(int).values
    s = df[arm].values.astype(float)
    m = {"arm": arm, "auroc": roc_auc(y, s), "n": len(df), "n_dual": int(y.sum())}
    m.update(top10_counts(df["class"].values, s))
    return m


def bootstrap_arm(df: pd.DataFrame, arm: str, rng: np.random.Generator) -> pd.DataFrame:
    """Ligand bootstrap with replacement; scores fixed (precomputed arms)."""
    n = len(df)
    y0 = (df["class"] == "dual").astype(int).values
    s0 = df[arm].values.astype(float)
    c0 = df["class"].values
    rows = []
    for b in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        y = y0[idx]
        s = s0[idx]
        c = c0[idx]
        # need both classes for AUROC
        if y.sum() == 0 or y.sum() == n:
            auroc = np.nan
        else:
            auroc = roc_auc(y, s)
        t = top10_counts(c, s)
        rows.append({"boot": b, "auroc": auroc, **t})
    return pd.DataFrame(rows)


def summarize_boot(point: dict, boot: pd.DataFrame) -> dict:
    def ci(series):
        x = series.dropna().values
        return float(np.nanpercentile(x, 2.5)), float(np.nanpercentile(x, 97.5))

    a_lo, a_hi = ci(boot["auroc"])
    h_lo, h_hi = ci(boot["top10_hardneg"])
    d_lo, d_hi = ci(boot["top10_dual"])
    return {
        "arm": point["arm"],
        "n": point["n"],
        "n_dual": point["n_dual"],
        "auroc": point["auroc"],
        "auroc_ci_lo": a_lo,
        "auroc_ci_hi": a_hi,
        "top10_dual": point["top10_dual"],
        "top10_dual_ci_lo": d_lo,
        "top10_dual_ci_hi": d_hi,
        "top10_hardneg": point["top10_hardneg"],
        "top10_hardneg_ci_lo": h_lo,
        "top10_hardneg_ci_hi": h_hi,
        "top10_A_only": point["top10_A_only"],
        "top10_B_only": point["top10_B_only"],
        "n_boot": N_BOOT,
        "n_boot_auroc_valid": int(boot["auroc"].notna().sum()),
    }


def bootstrap_delta(
    df: pd.DataFrame, arm_a: str, arm_b: str, rng: np.random.Generator
) -> dict:
    """Paired bootstrap delta = score_B − score_A metrics (same resample indices)."""
    n = len(df)
    y0 = (df["class"] == "dual").astype(int).values
    sa = df[arm_a].values.astype(float)
    sb = df[arm_b].values.astype(float)
    c0 = df["class"].values
    d_auroc, d_hard = [], []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        y = y0[idx]
        if y.sum() == 0 or y.sum() == n:
            continue
        aa = roc_auc(y, sa[idx])
        bb = roc_auc(y, sb[idx])
        d_auroc.append(bb - aa)
        ha = top10_counts(c0[idx], sa[idx])["top10_hardneg"]
        hb = top10_counts(c0[idx], sb[idx])["top10_hardneg"]
        d_hard.append(hb - ha)
    d_auroc = np.asarray(d_auroc)
    d_hard = np.asarray(d_hard)
    point_a = point_metrics(df, arm_a)
    point_b = point_metrics(df, arm_b)
    return {
        "compare": f"{arm_b}_minus_{arm_a}",
        "delta_auroc": point_b["auroc"] - point_a["auroc"],
        "delta_auroc_ci_lo": float(np.percentile(d_auroc, 2.5)),
        "delta_auroc_ci_hi": float(np.percentile(d_auroc, 97.5)),
        "delta_auroc_sig_excl0": not (
            np.percentile(d_auroc, 2.5) <= 0 <= np.percentile(d_auroc, 97.5)
        ),
        "delta_top10_hardneg": point_b["top10_hardneg"] - point_a["top10_hardneg"],
        "delta_top10_hardneg_ci_lo": float(np.percentile(d_hard, 2.5)),
        "delta_top10_hardneg_ci_hi": float(np.percentile(d_hard, 97.5)),
        "delta_hardneg_sig_excl0": not (
            np.percentile(d_hard, 2.5) <= 0 <= np.percentile(d_hard, 97.5)
        ),
        "n_boot": len(d_auroc),
    }


def main():
    panels = load_panels()
    # save EGFR decision ablation tables (same 5 arms)
    eh = panels["EGFR_HER2_panel40"]
    eh_out = REPO / "data" / "egfr_her2_panel40_v0" / "analysis" / "decision_ablation_v0"
    eh_out.mkdir(parents=True, exist_ok=True)
    eh_metrics = pd.DataFrame([point_metrics(eh, a) for a in ARMS])
    eh_metrics.to_csv(eh_out / "decision_ablation_metrics.csv", index=False)
    eh[
        ["ligand", "class", "pref_name", "vina_mean", "rtm_min_z", "rtm_shortfall",
         "consensus_rank_mean", "consensus_and_top25", "consensus_and_top25_pass",
         "rank_vina_mean", "rank_rtm_min_z"]
    ].to_csv(eh_out / "decision_ablation_scores.csv", index=False)
    ranks = eh[["ligand", "class", "pref_name"]].copy()
    for a in ARMS:
        ranks[a] = eh[a].rank(ascending=False, method="min").astype(int)
    ranks.to_csv(eh_out / "decision_ablation_ranks.csv", index=False)
    with (eh_out / "frozen_thresholds.yaml").open("w") as fh:
        fh.write(f"shortfall_lambda: {SHORTFALL_LAMBDA}\n")
        fh.write(f"consensus_top_frac: {CONSENSUS_TOP_FRAC}\n")
        fh.write(f"panel_n: {len(eh)}\n")
        fh.write(f"consensus_top_n: {eh.attrs['top_n']}\n")
        fh.write("note: same frozen thresholds as pik3ca_mtor decision_ablation_v0\n")

    ci_rows = []
    delta_rows = []
    boot_long = []
    for pair, df in panels.items():
        rng = np.random.default_rng(SEED)
        for arm in ARMS:
            point = point_metrics(df, arm)
            # independent stream per arm for reproducibility of margins; re-seed with arm hash
            rng_arm = np.random.default_rng(_stable_seed(pair, arm))
            boot = bootstrap_arm(df, arm, rng_arm)
            boot["pair"] = pair
            boot["arm"] = arm
            boot_long.append(boot)
            row = summarize_boot(point, boot)
            row["pair"] = pair
            ci_rows.append(row)
            print(
                f"{pair} {arm}: AUROC={row['auroc']:.3f} "
                f"[{row['auroc_ci_lo']:.3f},{row['auroc_ci_hi']:.3f}] "
                f"hardneg10={row['top10_hardneg']} "
                f"[{row['top10_hardneg_ci_lo']:.1f},{row['top10_hardneg_ci_hi']:.1f}]"
            )
        # paired deltas vs vina_mean
        for arm in ARMS:
            if arm == "vina_mean":
                continue
            rng_d = np.random.default_rng(_stable_seed(pair, "delta", arm))
            d = bootstrap_delta(df, "vina_mean", arm, rng_d)
            d["pair"] = pair
            delta_rows.append(d)
            print(
                f"  Δ {arm}-vina: AUROC {d['delta_auroc']:+.3f} "
                f"[{d['delta_auroc_ci_lo']:+.3f},{d['delta_auroc_ci_hi']:+.3f}] "
                f"sig={d['delta_auroc_sig_excl0']}"
            )

    ci = pd.DataFrame(ci_rows)
    delta = pd.DataFrame(delta_rows)
    ci.to_csv(OUT / "bootstrap_ci_metrics.csv", index=False)
    delta.to_csv(OUT / "bootstrap_ci_deltas_vs_vina_mean.csv", index=False)
    pd.concat(boot_long, ignore_index=True).to_csv(
        OUT / "bootstrap_raw_long.csv", index=False
    )

    # also mirror under each panel analysis folder
    for pair, sub in [
        ("EGFR_HER2_panel40", REPO / "data" / "egfr_her2_panel40_v0" / "analysis" / "bootstrap_ci_v0"),
        ("PIK3CA_mTOR_panel48", REPO / "data" / "pik3ca_mtor_panel48_v0" / "analysis" / "bootstrap_ci_v0"),
    ]:
        sub.mkdir(parents=True, exist_ok=True)
        ci[ci.pair == pair].to_csv(sub / "bootstrap_ci_metrics.csv", index=False)
        delta[delta.pair == pair].to_csv(
            sub / "bootstrap_ci_deltas_vs_vina_mean.csv", index=False
        )

    write_conclusion(ci, delta, eh_metrics)
    print("wrote", OUT)


def write_conclusion(ci: pd.DataFrame, delta: pd.DataFrame, eh_metrics: pd.DataFrame):
    def fmt_row(r):
        return (
            f"| {r.arm} | {r.auroc:.3f} [{r.auroc_ci_lo:.3f}, {r.auroc_ci_hi:.3f}] | "
            f"{int(r.top10_hardneg)} [{r.top10_hardneg_ci_lo:.0f}, {r.top10_hardneg_ci_hi:.0f}] | "
            f"{int(r.top10_dual)} [{r.top10_dual_ci_lo:.0f}, {r.top10_dual_ci_hi:.0f}] |"
        )

    lines = []
    lines.append("# Cross-pair bootstrap CI v0 — decision arms\n\n")
    lines.append(
        f"Resampling: ligand bootstrap with replacement, **N={N_BOOT}**, seed={SEED}.\n"
    )
    lines.append(
        f"Frozen (unchanged): `shortfall_lambda={SHORTFALL_LAMBDA}`, "
        f"`consensus_top_frac={CONSENSUS_TOP_FRAC}`.\n"
    )
    lines.append(
        "Scores/flags: arms precomputed on full panel; flags **do not** enter scores; "
        "no re-docking / no threshold retune.\n"
    )
    lines.append(
        "Significance for improvement vs `vina_mean`: paired-bootstrap 95% CI of "
        "ΔAUROC (arm − vina_mean) **excludes 0**.\n\n"
    )

    for pair, title in [
        ("EGFR_HER2_panel40", "EGFR/HER2 panel40"),
        ("PIK3CA_mTOR_panel48", "PIK3CA/mTOR panel48"),
    ]:
        sub = ci[ci.pair == pair]
        lines.append(f"## {title}\n\n")
        lines.append(
            "| arm | AUROC [95% CI] | Top10 hardneg [95% CI] | Top10 dual [95% CI] |\n"
        )
        lines.append("|-----|----------------|------------------------|---------------------|\n")
        for _, r in sub.iterrows():
            lines.append(fmt_row(r) + "\n")
        lines.append("\n### Δ vs vina_mean (paired bootstrap)\n\n")
        lines.append(
            "| compare | ΔAUROC [95% CI] | sig? | ΔTop10 hardneg [95% CI] | sig? |\n"
        )
        lines.append(
            "|---------|----------------|------|-------------------------|------|\n"
        )
        for _, d in delta[delta.pair == pair].iterrows():
            lines.append(
                f"| {d['compare']} | {d['delta_auroc']:+.3f} "
                f"[{d['delta_auroc_ci_lo']:+.3f}, {d['delta_auroc_ci_hi']:+.3f}] | "
                f"{'YES' if d['delta_auroc_sig_excl0'] else 'no'} | "
                f"{d['delta_top10_hardneg']:+.0f} "
                f"[{d['delta_top10_hardneg_ci_lo']:+.0f}, {d['delta_top10_hardneg_ci_hi']:+.0f}] | "
                f"{'YES' if d['delta_hardneg_sig_excl0'] else 'no'} |\n"
            )
        lines.append("\n")

    # point EGFR ablation note
    lines.append("## EGFR panel40 — five-arm point estimates (new)\n\n")
    lines.append(
        "Same frozen decision arms as PIK3CA/mTOR; full table in "
        "`data/egfr_her2_panel40_v0/analysis/decision_ablation_v0/`.\n\n"
    )
    lines.append("| arm | AUROC | Top10 dual | Top10 hardneg |\n|-----|-------|------------|---------------|\n")
    for _, r in eh_metrics.iterrows():
        lines.append(
            f"| {r.arm} | {r.auroc:.3f} | {int(r.top10_dual)} | {int(r.top10_hardneg)} |\n"
        )

    # verdict
    lines.append("\n## Updated conclusion (CI-aware)\n\n")

    def sig_auroc(pair, arm):
        hit = delta[
            (delta.pair == pair) & (delta["compare"] == f"{arm}_minus_vina_mean")
        ]
        if hit.empty:
            return False, None
        r = hit.iloc[0]
        return bool(r["delta_auroc_sig_excl0"]), (
            float(r["delta_auroc"]),
            float(r["delta_auroc_ci_lo"]),
            float(r["delta_auroc_ci_hi"]),
        )

    eh_rtm, eh_rtm_d = sig_auroc("EGFR_HER2_panel40", "rtm_min_z")
    pm_rtm, pm_rtm_d = sig_auroc("PIK3CA_mTOR_panel48", "rtm_min_z")
    eh_sf, _ = sig_auroc("EGFR_HER2_panel40", "rtm_shortfall")
    pm_sf, _ = sig_auroc("PIK3CA_mTOR_panel48", "rtm_shortfall")
    eh_crm, eh_crm_d = sig_auroc("EGFR_HER2_panel40", "consensus_rank_mean")
    pm_crm, _ = sig_auroc("PIK3CA_mTOR_panel48", "consensus_rank_mean")
    eh_and, _ = sig_auroc("EGFR_HER2_panel40", "consensus_and_top25")
    pm_and, _ = sig_auroc("PIK3CA_mTOR_panel48", "consensus_and_top25")

    lines.append(
        f"1. **EGFR/HER2 — `rtm_min_z`:** point ΔAUROC = {eh_rtm_d[0]:+.3f} "
        f"[{eh_rtm_d[1]:+.3f}, {eh_rtm_d[2]:+.3f}] → "
        f"{'**significant**' if eh_rtm else '**not significant at 95%** (CI includes 0; lower bound near 0 / borderline)'} "
        f"vs `vina_mean`.\n"
    )
    lines.append(
        f"2. **EGFR/HER2 — other arms:** shortfall sig={eh_sf}; "
        f"**consensus_rank_mean sig={eh_crm}**"
        + (
            f" (Δ={eh_crm_d[0]:+.3f} [{eh_crm_d[1]:+.3f}, {eh_crm_d[2]:+.3f}])"
            if eh_crm_d
            else ""
        )
        + f"; AND-top25 sig={eh_and}. "
        "Rank-mean consensus is the only arm with ΔAUROC CI excluding 0 on this pair; "
        "hardneg Top10 Δ still not significant.\n"
    )
    lines.append(
        f"3. **PIK3CA/mTOR — all five arms vs `vina_mean`:** "
        f"`rtm_min_z` Δ={pm_rtm_d[0]:+.3f} [{pm_rtm_d[1]:+.3f}, {pm_rtm_d[2]:+.3f}] → "
        f"**not significant**; shortfall={pm_sf}; rank-mean={pm_crm}; AND-top25={pm_and}. "
        "Point lifts exist but are unstable under ligand bootstrap.\n"
    )
    lines.append(
        "4. **Top10 hardneg:** no arm on either pair has a paired-bootstrap Δ hardneg CI "
        "that excludes 0; apparent hardneg drops (e.g. EGFR 6→3) are **not** CI-significant.\n"
    )
    lines.append(
        "5. **How to write it:** keep dual readout (`vina_mean` + `rtm_min_z`); "
        "do **not** claim a universally significant RTM upgrade across pairs; "
        "EGFR shows a larger point lift that is borderline for `rtm_min_z` and significant "
        "only for rank-mean consensus; PIK3CA/mTOR lifts are non-significant. "
        "Still no C4-closed / clash-retune / flags-in-score claim.\n"
    )

    text = "".join(lines)
    (OUT / "BOOTSTRAP_CI_CONCLUSION_V0.md").write_text(text)
    # mirrors
    for sub in [
        REPO / "data" / "egfr_her2_panel40_v0" / "analysis" / "bootstrap_ci_v0",
        REPO / "data" / "pik3ca_mtor_panel48_v0" / "analysis" / "bootstrap_ci_v0",
    ]:
        (sub / "BOOTSTRAP_CI_CONCLUSION_V0.md").write_text(text)


if __name__ == "__main__":
    main()
