#!/usr/bin/env python3
"""Audit measurement-frequency imbalance in the frozen scored panels.

This script uses the deposited API-refetch table. It does not reconstruct
assay confidence, species, target form, or source-document provenance and must
not be described as a high-confidence label rebuild.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "data" / "jcim_novelty_v0" / "tables" / "assay_max_vs_median_ligand_v1.csv"
OUT_SUMMARY = ROOT / "data" / "jcim_novelty_v0" / "tables" / "measurement_frequency_by_class_v1.csv"
OUT_CORR = ROOT / "data" / "jcim_novelty_v0" / "tables" / "measurement_frequency_max_median_v1.csv"


def read_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def num(row: dict, key: str) -> float:
    return float(row[key])


def summarize(values: list[float]) -> tuple[float, float, float, float]:
    arr = np.asarray(values, dtype=float)
    return tuple(float(x) for x in (np.median(arr), np.percentile(arr, 25), np.percentile(arr, 75), np.max(arr)))


def main() -> None:
    rows = read_rows(SOURCE)
    pairs = list(dict.fromkeys(r["pair"] for r in rows))
    classes = ("dual", "A_only", "B_only", "neither")
    summary = []
    correlations = []

    for pair in pairs:
        pair_rows = [r for r in rows if r["pair"] == pair]
        for cls in classes:
            kept = [r for r in pair_rows if r["frozen_class"] == cls]
            if not kept:
                continue
            n_a = [num(r, "n_activity_A") for r in kept]
            n_b = [num(r, "n_activity_B") for r in kept]
            total = [a + b for a, b in zip(n_a, n_b)]
            med, q1, q3, maximum = summarize(total)
            summary.append(
                {
                    "pair": pair,
                    "class": cls,
                    "n_ligands": len(kept),
                    "median_activity_records_total": round(med, 3),
                    "q1_activity_records_total": round(q1, 3),
                    "q3_activity_records_total": round(q3, 3),
                    "max_activity_records_total": int(maximum),
                    "fraction_repeated_pchembl_A": round(sum(num(r, "n_pchembl_A") > 1 for r in kept) / len(kept), 4),
                    "fraction_repeated_pchembl_B": round(sum(num(r, "n_pchembl_B") > 1 for r in kept) / len(kept), 4),
                    "note": "API activity-record counts; not assay-confidence or source-document audit",
                }
            )

        for end in ("A", "B"):
            counts = np.asarray([num(r, f"n_pchembl_{end}") for r in pair_rows], dtype=float)
            gaps = np.asarray([num(r, f"api_max_{end}") - num(r, f"api_median_{end}") for r in pair_rows], dtype=float)
            rho, p = spearmanr(np.log1p(counts), gaps)
            correlations.append(
                {
                    "pair": pair,
                    "end": end,
                    "n_ligands": len(pair_rows),
                    "n_repeated": int((counts > 1).sum()),
                    "median_max_minus_median": round(float(np.median(gaps)), 4),
                    "q95_max_minus_median": round(float(np.percentile(gaps, 95)), 4),
                    "spearman_log_count_vs_gap": round(float(rho), 4),
                    "spearman_p_unadjusted": round(float(p), 6),
                    "note": "diagnostic association; includes singletons with zero aggregation gap; not causal",
                }
            )

    write_rows(OUT_SUMMARY, summary)
    write_rows(OUT_CORR, correlations)
    print(f"wrote {OUT_SUMMARY} ({len(summary)} rows)")
    print(f"wrote {OUT_CORR} ({len(correlations)} rows)")


if __name__ == "__main__":
    main()
