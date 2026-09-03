#!/usr/bin/env python3
"""Multi-seed Vina sensitivity with Table 3 Dual-versus-neither estimand.

Does not re-dock. Reads frozen `multiseed_scores_long_v1.csv`.

Primary Dual-versus-neither (matches Table 3):
    vina_mean = (S_A + S_B) / 2 per ligand, then one AUROC.

Legacy sensitivity (v1 column; do not cite as Table 3 / Table S54 primary):
    mean_marginal_pocket_auroc_D_vs_neither = mean(AUC_A, AUC_B).

Directional D vs A, D vs B, and summary_min are unchanged from v1
(pocket-matched). Formulation gap uses the vina_mean estimand.

Writes v2 tables only; v1 files remain as a dated wrong-estimand snapshot.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_multiseed_v0"
TAB = OUT / "tables"
AN = OUT / "analysis"
AN.mkdir(parents=True, exist_ok=True)

SCORE_PATH = TAB / "multiseed_scores_long_v1.csv"
FORMULATION_PATH = (
    ROOT / "data" / "jcim_novelty_v0" / "tables" / "formulation_conventional_vs_directional_v1.csv"
)
PRIMARY_SEED = 20260727
# Table 3 lock: AUC(vina_mean) Dual vs neither, 4 d.p. (EGFR stored as 0.756).
TABLE3_VINA_MEAN = {
    "EGFR/HER2": 0.7560,
    "AChE/BChE": 0.6494,
    "PIK3CA/PIK3CB": 0.5592,
    "PIK3CA/mTOR": 0.5139,
}
PAIR_ORDER = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]


def auroc(pos, neg) -> float:
    if not pos or not neg:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def r4(x: float) -> float:
    if x != x:  # NaN
        return x
    return float(f"{x:.4f}")


def load_scores():
    if not SCORE_PATH.exists():
        raise SystemExit(f"missing {SCORE_PATH}")
    rows = list(csv.DictReader(SCORE_PATH.open(encoding="utf-8")))
    return rows


def main() -> int:
    rows = load_scores()
    by = defaultdict(dict)
    for r in rows:
        if r["status"] not in {"ok", "exists", "primary_reused"}:
            continue
        if r["vina_mode1"] in ("", None):
            continue
        key = (r["pair"], int(r["seed"]), r["ligand"])
        by[key].setdefault("class", r["class"])
        by[key][r["pocket"]] = -float(r["vina_mode1"])  # higher = better

    grouped = defaultdict(list)
    for (pair, seed, lig), d in by.items():
        if "A" in d and "B" in d:
            grouped[(pair, seed)].append({"ligand": lig, **d})

    metrics = []
    for (pair, seed), recs in sorted(grouped.items()):
        dual = [r for r in recs if r["class"] == "dual"]
        aonly = [r for r in recs if r["class"] in ("A_only", "A-only")]
        bonly = [r for r in recs if r["class"] in ("B_only", "B-only")]
        neither = [r for r in recs if r["class"] == "neither"]
        d_vs_a = auroc([r["B"] for r in dual], [r["B"] for r in aonly])
        d_vs_b = auroc([r["A"] for r in dual], [r["A"] for r in bonly])
        d_vs_n_a = auroc([r["A"] for r in dual], [r["A"] for r in neither])
        d_vs_n_b = auroc([r["B"] for r in dual], [r["B"] for r in neither])
        mean_marginal = float(np.nanmean([d_vs_n_a, d_vs_n_b]))
        dual_mean = [(r["A"] + r["B"]) / 2.0 for r in dual]
        neither_mean = [(r["A"] + r["B"]) / 2.0 for r in neither]
        d_vs_n_vina_mean = auroc(dual_mean, neither_mean)
        smin = float(np.nanmin([d_vs_a, d_vs_b]))
        gap = d_vs_n_vina_mean - smin
        metrics.append(
            {
                "pair": pair,
                "seed": seed,
                "n_complete": len(recs),
                "n_dual": len(dual),
                "n_A_only": len(aonly),
                "n_B_only": len(bonly),
                "n_neither": len(neither),
                "auroc_dual_vs_A_only": r4(d_vs_a),
                "auroc_dual_vs_B_only": r4(d_vs_b),
                "summary_min": r4(smin),
                "auroc_dual_vs_neither_vina_mean": r4(d_vs_n_vina_mean),
                "mean_marginal_pocket_auroc_D_vs_neither": r4(mean_marginal),
                "formulation_gap_neither_minus_summary_min": r4(gap),
                "_raw_vina_mean": d_vs_n_vina_mean,
                "_raw_smin": smin,
                "_raw_gap": gap,
                "_raw_marginal": mean_marginal,
            }
        )

    # Gate: primary seed must recover Table 3 Dual-versus-neither.
    mismatches = []
    for pair, expected in TABLE3_VINA_MEAN.items():
        rec = next(
            (m for m in metrics if m["pair"] == pair and int(m["seed"]) == PRIMARY_SEED),
            None,
        )
        if rec is None:
            mismatches.append((pair, "missing primary seed", expected, None))
            continue
        got = rec["_raw_vina_mean"]
        if abs(got - expected) > 5e-4:
            mismatches.append((pair, "vina_mean mismatch", expected, got))
    if mismatches:
        print("FATAL: primary seed did not recover Table 3 AUC(vina_mean).", file=sys.stderr)
        for row in mismatches:
            print(f"  {row}", file=sys.stderr)
        print("Stop. Do not write manuscript numbers from this run.", file=sys.stderr)
        return 1

    out_fields = [
        "pair",
        "seed",
        "n_complete",
        "n_dual",
        "n_A_only",
        "n_B_only",
        "n_neither",
        "auroc_dual_vs_A_only",
        "auroc_dual_vs_B_only",
        "summary_min",
        "auroc_dual_vs_neither_vina_mean",
        "mean_marginal_pocket_auroc_D_vs_neither",
        "formulation_gap_neither_minus_summary_min",
    ]
    out_csv = TAB / "multiseed_auroc_by_seed_v2.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, extrasaction="ignore")
        w.writeheader()
        for m in metrics:
            w.writerow({k: m[k] for k in out_fields})

    agg_keys = [
        "auroc_dual_vs_A_only",
        "auroc_dual_vs_B_only",
        "summary_min",
        "auroc_dual_vs_neither_vina_mean",
        "mean_marginal_pocket_auroc_D_vs_neither",
        "formulation_gap_neither_minus_summary_min",
    ]
    agg = []
    by_pair = defaultdict(list)
    for m in metrics:
        by_pair[m["pair"]].append(m)
    for pair in PAIR_ORDER:
        ms = by_pair.get(pair, [])
        for key in agg_keys:
            vals = [float(m[key]) for m in ms]
            primary = next((m for m in ms if int(m["seed"]) == PRIMARY_SEED), None)
            agg.append(
                {
                    "pair": pair,
                    "metric": key,
                    "n_seeds": len(vals),
                    "primary_seed_value": primary[key] if primary else "",
                    "median": r4(float(np.median(vals))),
                    "min": r4(float(np.min(vals))),
                    "max": r4(float(np.max(vals))),
                    "iqr": r4(float(np.percentile(vals, 75) - np.percentile(vals, 25))),
                    "range": r4(float(np.max(vals) - np.min(vals))),
                }
            )

    consistency = []
    for pair in PAIR_ORDER:
        ms = by_pair.get(pair, [])
        primary = next((m for m in ms if int(m["seed"]) == PRIMARY_SEED), None)
        if primary is None:
            continue
        prim_gap = float(primary["formulation_gap_neither_minus_summary_min"])
        prim_smin = float(primary["summary_min"])
        same_gap_sign = sum(
            1
            for m in ms
            if (
                np.sign(float(m["formulation_gap_neither_minus_summary_min"])) == np.sign(prim_gap)
                or (prim_gap == 0 and float(m["formulation_gap_neither_minus_summary_min"]) == 0)
            )
        )
        n_positive_gap = sum(
            1 for m in ms if float(m["formulation_gap_neither_minus_summary_min"]) > 0
        )
        same_form_order = sum(
            1
            for m in ms
            if (
                float(m["auroc_dual_vs_neither_vina_mean"]) > float(m["summary_min"])
            )
            == (
                float(primary["auroc_dual_vs_neither_vina_mean"])
                > float(primary["summary_min"])
            )
        )
        consistency.append(
            {
                "pair": pair,
                "n_seeds": len(ms),
                "primary_summary_min": prim_smin,
                "primary_auroc_dual_vs_neither_vina_mean": primary[
                    "auroc_dual_vs_neither_vina_mean"
                ],
                "primary_formulation_gap": prim_gap,
                "n_seeds_same_gap_sign": same_gap_sign,
                "n_seeds_positive_gap": n_positive_gap,
                "n_seeds_same_neither_gt_summary_min_order": same_form_order,
                "summary_min_median": r4(
                    float(np.median([float(m["summary_min"]) for m in ms]))
                ),
                "summary_min_min": r4(
                    float(np.min([float(m["summary_min"]) for m in ms]))
                ),
                "summary_min_max": r4(
                    float(np.max([float(m["summary_min"]) for m in ms]))
                ),
                "summary_min_range": r4(
                    float(
                        np.max([float(m["summary_min"]) for m in ms])
                        - np.min([float(m["summary_min"]) for m in ms])
                    )
                ),
            }
        )

    with (TAB / "multiseed_auroc_aggregate_v2.csv").open(
        "w", encoding="utf-8", newline=""
    ) as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0].keys()) if agg else ["pair"])
        w.writeheader()
        w.writerows(agg)
    with (TAB / "multiseed_consistency_v2.csv").open(
        "w", encoding="utf-8", newline=""
    ) as fh:
        w = csv.DictWriter(
            fh, fieldnames=list(consistency[0].keys()) if consistency else ["pair"]
        )
        w.writeheader()
        w.writerows(consistency)

    lines = [
        "# Four-pair multi-seed Vina sensitivity (v2; Table 3 estimand)",
        "",
        f"Source scores: `{SCORE_PATH.relative_to(ROOT)}`",
        "Frozen seeds: 20260727 (primary, reused) + 20260811–20260814.",
        "Protocol otherwise identical to production (receptors, boxes, exhaustiveness, modes, energy_range).",
        "",
        "## Estimands",
        "",
        "- Directional Dual vs A-only / Dual vs B-only / `summary_min`: pocket-matched; unchanged from v1.",
        "- **Primary Dual-versus-neither:** per-ligand `vina_mean = (S_A+S_B)/2`, then one AUROC.",
        "  This is the Table 3 estimand. Primary seed 20260727 recovered",
        "  0.7560 / 0.6494 / 0.5592 / 0.5139.",
        "- **Sensitivity only:** `mean_marginal_pocket_auroc_D_vs_neither = mean(AUC_A, AUC_B)`",
        "  (v1 column; do not cite as Table 3 or as the Table S54 Dual-versus-neither value).",
        "- Formulation gap = Dual-versus-neither (`vina_mean`) − `summary_min`.",
        "",
        "## Per-seed metrics (primary Dual-versus-neither = vina_mean)",
        "",
        "| pair | seed | dual_vs_A | dual_vs_B | summary_min | dual_vs_neither_vina_mean | gap | mean_marginal (legacy) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for m in metrics:
        lines.append(
            f"| {m['pair']} | {m['seed']} | {m['auroc_dual_vs_A_only']} | {m['auroc_dual_vs_B_only']} | "
            f"{m['summary_min']} | {m['auroc_dual_vs_neither_vina_mean']} | "
            f"{m['formulation_gap_neither_minus_summary_min']} | "
            f"{m['mean_marginal_pocket_auroc_D_vs_neither']} |"
        )
    lines += ["", "## Consistency vs primary seed (vina_mean gap)", ""]
    for c in consistency:
        lines.append(
            f"- **{c['pair']}**: summary_min median {c['summary_min_median']} "
            f"(range {c['summary_min_min']}–{c['summary_min_max']}); "
            f"gap-sign match {c['n_seeds_same_gap_sign']}/{c['n_seeds']}; "
            f"positive gap {c['n_seeds_positive_gap']}/{c['n_seeds']}; "
            f"neither>summary_min order match {c['n_seeds_same_neither_gt_summary_min_order']}/{c['n_seeds']}."
        )
    lines += [
        "",
        "## Claim ceiling",
        "",
        "- Allowed: report median/IQR/range across frozen seeds; state whether the primary qualitative pattern held.",
        "- Forbidden: picking a favorable seed; replacing primary Table 2 with a multi-seed mean; claiming seed robustness beyond these four pairs; citing v1 Dual-versus-neither as Table 3.",
        "",
        "v1 tables remain as a dated wrong-estimand snapshot and must not be copied into the article.",
        "",
    ]
    (AN / "MULTISEED_VINA_VERDICT_V2.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out_csv}")
    print(f"wrote {TAB / 'multiseed_auroc_aggregate_v2.csv'}")
    print(f"wrote {TAB / 'multiseed_consistency_v2.csv'}")
    print(f"wrote {AN / 'MULTISEED_VINA_VERDICT_V2.md'}")
    print("primary-seed Table 3 recovery: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
