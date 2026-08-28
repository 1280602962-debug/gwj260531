#!/usr/bin/env python3
"""Summarize four-pair multi-seed Vina sensitivity (after docking completes)."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_multiseed_v0"
TAB = OUT / "tables"
AN = OUT / "analysis"
AN.mkdir(parents=True, exist_ok=True)

SCORE_PATHS = [
    TAB / "multiseed_scores_long_v1.csv",
    TAB / "multiseed_scores_long_partial_v1.csv",
]


def auroc(pos, neg) -> float:
    if not pos or not neg:
        return float("nan")
    wins = sum(p > n for p in pos for n in neg)
    ties = sum(p == n for p in pos for n in neg)
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def load_scores():
    path = next((p for p in SCORE_PATHS if p.exists()), None)
    if path is None:
        raise SystemExit("no multiseed score table yet")
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    return path, rows


def main():
    path, rows = load_scores()
    # pair,seed,ligand -> {A,B,class,status}
    by = defaultdict(dict)
    for r in rows:
        if r["status"] not in {"ok", "exists", "primary_reused"}:
            continue
        if r["vina_mode1"] in ("", None):
            continue
        key = (r["pair"], int(r["seed"]), r["ligand"])
        by[key].setdefault("class", r["class"])
        by[key][r["pocket"]] = -float(r["vina_mode1"])  # higher = better

    metrics = []
    for (pair, seed), group in defaultdict(list).items() if False else []:
        pass

    grouped = defaultdict(list)
    for (pair, seed, lig), d in by.items():
        if "A" in d and "B" in d:
            grouped[(pair, seed)].append({"ligand": lig, **d})

    for (pair, seed), recs in sorted(grouped.items()):
        dual = [r for r in recs if r["class"] == "dual"]
        aonly = [r for r in recs if r["class"] in ("A_only", "A-only")]
        bonly = [r for r in recs if r["class"] in ("B_only", "B-only")]
        neither = [r for r in recs if r["class"] == "neither"]
        # pocket-matched: D vs A uses B score; D vs B uses A score
        d_vs_a = auroc([r["B"] for r in dual], [r["B"] for r in aonly])
        d_vs_b = auroc([r["A"] for r in dual], [r["A"] for r in bonly])
        d_vs_n_a = auroc([r["A"] for r in dual], [r["A"] for r in neither])
        d_vs_n_b = auroc([r["B"] for r in dual], [r["B"] for r in neither])
        d_vs_n = float(np.nanmean([d_vs_n_a, d_vs_n_b]))
        smin = float(np.nanmin([d_vs_a, d_vs_b]))
        gap = d_vs_n - smin
        metrics.append(
            {
                "pair": pair,
                "seed": seed,
                "n_complete": len(recs),
                "n_dual": len(dual),
                "n_A_only": len(aonly),
                "n_B_only": len(bonly),
                "n_neither": len(neither),
                "auroc_dual_vs_A_only": round(d_vs_a, 4),
                "auroc_dual_vs_B_only": round(d_vs_b, 4),
                "summary_min": round(smin, 4),
                "auroc_dual_vs_neither_mean": round(d_vs_n, 4),
                "formulation_gap_neither_minus_summary_min": round(gap, 4),
            }
        )

    out_csv = TAB / "multiseed_auroc_by_seed_v1.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        fields = list(metrics[0].keys()) if metrics else [
            "pair",
            "seed",
            "n_complete",
            "summary_min",
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(metrics)

    # Aggregate across seeds
    agg = []
    by_pair = defaultdict(list)
    for m in metrics:
        by_pair[m["pair"]].append(m)
    for pair, ms in sorted(by_pair.items()):
        for key in [
            "auroc_dual_vs_A_only",
            "auroc_dual_vs_B_only",
            "summary_min",
            "auroc_dual_vs_neither_mean",
            "formulation_gap_neither_minus_summary_min",
        ]:
            vals = [float(m[key]) for m in ms]
            primary = next((m for m in ms if int(m["seed"]) == 20260727), None)
            # directional consistency vs primary for summary_min / gap sign
            agg.append(
                {
                    "pair": pair,
                    "metric": key,
                    "n_seeds": len(vals),
                    "primary_seed_value": primary[key] if primary else "",
                    "median": round(float(np.median(vals)), 4),
                    "min": round(float(np.min(vals)), 4),
                    "max": round(float(np.max(vals)), 4),
                    "iqr": round(float(np.percentile(vals, 75) - np.percentile(vals, 25)), 4),
                    "range": round(float(np.max(vals) - np.min(vals)), 4),
                }
            )

    # Consistency: how many seeds preserve EGFR formulation gap sign and weak-arm pattern
    consistency = []
    for pair, ms in sorted(by_pair.items()):
        primary = next((m for m in ms if int(m["seed"]) == 20260727), None)
        if primary is None:
            continue
        prim_gap = float(primary["formulation_gap_neither_minus_summary_min"])
        prim_smin = float(primary["summary_min"])
        same_gap_sign = sum(
            1
            for m in ms
            if np.sign(float(m["formulation_gap_neither_minus_summary_min"])) == np.sign(prim_gap)
            or (prim_gap == 0 and float(m["formulation_gap_neither_minus_summary_min"]) == 0)
        )
        # "weak arm" for EGFR: summary_min < dual-vs-neither
        same_form_order = sum(
            1
            for m in ms
            if (float(m["auroc_dual_vs_neither_mean"]) > float(m["summary_min"]))
            == (float(primary["auroc_dual_vs_neither_mean"]) > float(primary["summary_min"]))
        )
        consistency.append(
            {
                "pair": pair,
                "n_seeds": len(ms),
                "primary_summary_min": prim_smin,
                "primary_formulation_gap": prim_gap,
                "n_seeds_same_gap_sign": same_gap_sign,
                "n_seeds_same_neither_gt_summary_min_order": same_form_order,
                "summary_min_median": round(float(np.median([float(m["summary_min"]) for m in ms])), 4),
                "summary_min_range": round(
                    float(np.max([float(m["summary_min"]) for m in ms]) - np.min([float(m["summary_min"]) for m in ms])),
                    4,
                ),
            }
        )

    with (TAB / "multiseed_auroc_aggregate_v1.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0].keys()) if agg else ["pair"])
        w.writeheader()
        w.writerows(agg)
    with (TAB / "multiseed_consistency_v1.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(consistency[0].keys()) if consistency else ["pair"])
        w.writeheader()
        w.writerows(consistency)

    lines = [
        "# Four-pair multi-seed Vina sensitivity",
        "",
        f"Source scores: `{path.relative_to(ROOT)}`",
        "Frozen seeds: 20260727 (primary, reused) + 20260811–20260814.",
        "Protocol otherwise identical to production (receptors, boxes, exhaustiveness, modes, energy_range).",
        "",
        "## Per-seed metrics",
        "",
        "| pair | seed | dual_vs_A | dual_vs_B | summary_min | dual_vs_neither | gap |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for m in metrics:
        lines.append(
            f"| {m['pair']} | {m['seed']} | {m['auroc_dual_vs_A_only']} | {m['auroc_dual_vs_B_only']} | "
            f"{m['summary_min']} | {m['auroc_dual_vs_neither_mean']} | {m['formulation_gap_neither_minus_summary_min']} |"
        )
    lines += ["", "## Consistency vs primary seed", ""]
    for c in consistency:
        lines.append(
            f"- **{c['pair']}**: summary_min median {c['summary_min_median']} "
            f"(range {c['summary_min_range']}); "
            f"gap-sign match {c['n_seeds_same_gap_sign']}/{c['n_seeds']}; "
            f"neither>summary_min order match {c['n_seeds_same_neither_gt_summary_min_order']}/{c['n_seeds']}."
        )
    lines += [
        "",
        "## Claim ceiling",
        "",
        "- Allowed: report median/IQR/range across frozen seeds; state whether the primary qualitative pattern held.",
        "- Forbidden: picking a favorable seed; replacing primary Table 2 with a multi-seed mean; claiming seed robustness beyond these four pairs.",
        "",
    ]
    (AN / "MULTISEED_VINA_VERDICT_V1.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out_csv}")
    print(f"wrote {TAB / 'multiseed_auroc_aggregate_v1.csv'}")
    print(f"wrote {AN / 'MULTISEED_VINA_VERDICT_V1.md'}")


if __name__ == "__main__":
    main()
