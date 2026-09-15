#!/usr/bin/env python3
"""Unified eight-pair top-10% and AND-filter operating points from frozen scores.

Reads review_scored_membership_v1.csv. Ranking uses k = ceil(0.10 n) on the
full four-state panel so pairs of different size are compared at the same
screening fraction. AND-filter rules still match
operating_point_examples_review_v1.csv (EGFR/HER2 and JAK1/TYK2).
Does not redock.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data/jcim_novelty_v0/tables"
MEMBERSHIP = TAB / "review_scored_membership_v1.csv"
EXAMPLES = TAB / "operating_point_examples_review_v1.csv"
OUT = TAB / "eight_pair_ranking_operating_point_v1.csv"
ORDER = ("EGFR/HER2", "JAK1/JAK2", "JAK1/TYK2", "PIK3CA/mTOR",
         "AChE/BChE", "F2/F10", "PPARG/PPARA", "PPARA/PPARD")
CLASSES = ("dual", "A_only", "B_only")
TOP_FRACTION = 0.10
FIELDS = (
    "pair", "n_ranked", "n_dual", "n_A_only", "n_B_only", "n_neither",
    "top_fraction", "top_k", "top_k_fraction", "top_dual", "top_A_only",
    "top_B_only", "top_neither", "top_dual_fraction", "top_ligand_ids",
    "threshold", "n_filter_input", "retained_dual", "retained_A_only",
    "retained_B_only", "dual_recall", "dual_precision", "filter_rule",
    "tie_rule", "ranking_rule",
)


def rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def top_k(n: int) -> int:
    return max(1, math.ceil(TOP_FRACTION * n))


def operating_point(recs: list[dict]) -> dict:
    n = len(recs)
    k = top_k(n)
    counts = Counter(r["cls"] for r in recs)
    ranked = sorted(recs, key=lambda r: (-r["vina_mean"], r["ligand"]))
    top = ranked[:k]
    ct = Counter(r["cls"] for r in top)
    duals = [r for r in recs if r["cls"] == "dual"]
    threshold = float(np.median([r["vina_worst"] for r in duals]))
    use = [r for r in recs if r["cls"] in CLASSES]
    retained = [r for r in use if r["vina_worst"] >= threshold]
    cr = Counter(r["cls"] for r in retained)
    n_dual = sum(r["cls"] == "dual" for r in use)
    return dict(
        n_ranked=n,
        n_dual=counts["dual"],
        n_A_only=counts["A_only"],
        n_B_only=counts["B_only"],
        n_neither=counts["neither"],
        top_fraction=TOP_FRACTION,
        top_k=k,
        top_k_fraction=k / n,
        top_dual=ct["dual"],
        top_A_only=ct["A_only"],
        top_B_only=ct["B_only"],
        top_neither=ct["neither"],
        top_dual_fraction=ct["dual"] / k,
        top_ligand_ids=";".join(r["ligand"] for r in top),
        threshold=threshold,
        n_filter_input=len(use),
        retained_dual=cr["dual"],
        retained_A_only=cr["A_only"],
        retained_B_only=cr["B_only"],
        dual_recall=cr["dual"] / n_dual if n_dual else float("nan"),
        dual_precision=cr["dual"] / len(retained) if retained else float("nan"),
        filter_rule="score_worst >= median_dual_score_worst; neither excluded",
        tie_rule="descending score_mean then ascending ligand ID",
        ranking_rule="top 10% on full four-state panel; k=ceil(0.10 n)",
    )


def build() -> list[dict]:
    by: dict[str, list[dict]] = defaultdict(list)
    for raw in rows(MEMBERSHIP):
        sa, sb = float(raw["score_A"]), float(raw["score_B"])
        by[raw["pair"]].append({
            "pair": raw["pair"],
            "ligand": raw["ligand"],
            "cls": raw["cls"],
            "vina_mean": (sa + sb) / 2.0,
            "vina_worst": min(sa, sb),
        })
    assert set(by) == set(ORDER), set(by)
    out = []
    for pair in ORDER:
        recs = by[pair]
        assert len({r["ligand"] for r in recs}) == len(recs), pair
        row = {"pair": pair, **operating_point(recs)}
        assert row["top_k"] == top_k(row["n_ranked"]), pair
        assert row["top_dual"] + row["top_A_only"] + row["top_B_only"] + row["top_neither"] == row["top_k"], pair
        out.append(row)
    examples = {r["pair"]: r for r in rows(EXAMPLES)}
    for pair, example in examples.items():
        got = next(r for r in out if r["pair"] == pair)
        for key in ("n_ranked", "n_filter_input", "retained_dual", "retained_A_only", "retained_B_only"):
            assert int(example[key]) == int(got[key]), (pair, key, example[key], got[key])
        assert abs(float(example["threshold"]) - float(got["threshold"])) < 1e-9, pair
        assert abs(float(example["dual_precision"]) - float(got["dual_precision"])) < 1e-12, pair
    return out


def csv_text(data: list[dict]) -> str:
    handle = __import__("io").StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(data)
    return handle.getvalue()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = csv_text(build())
    if args.check:
        assert OUT.read_text(encoding="utf-8") == text, f"replay differs: {OUT}"
        print("checked", OUT.name)
    else:
        OUT.write_text(text, encoding="utf-8", newline="")
        print("wrote", OUT.name)
    print("PASS: eight-pair top-10% / AND-filter operating points")


if __name__ == "__main__":
    main()
