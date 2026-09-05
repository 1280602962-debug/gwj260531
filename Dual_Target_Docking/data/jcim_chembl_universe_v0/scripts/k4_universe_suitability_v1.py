#!/usr/bin/env python3
"""Rank frozen K=4 pairs against the ChEMBL 37 universe census tables.

Does not query SQLite. Does not dock. Does not change Table 2 / K = 4.
Reads tables written by chembl_exhaustive_pair_census_v1.py.
"""
from __future__ import annotations

import csv
from pathlib import Path

TABLES = Path(__file__).resolve().parents[1] / "tables"

K4 = [
    {
        "pair": "EGFR/HER2",
        "gene_A": "EGFR",
        "gene_B": "ERBB2",
        "in_j0": 1,
        "in_k4": 1,
        "k4_role": "case_supply_limited",
        "paper_role": "pose-gold kinase case; not a thick-supply seat",
    },
    {
        "pair": "AChE/BChE",
        "gene_A": "ACHE",
        "gene_B": "BCHE",
        "in_j0": 1,
        "in_k4": 1,
        "k4_role": "development_homolog",
        "paper_role": "literature cholinesterase dual; homologous hydrolase",
    },
    {
        "pair": "PIK3CA/PIK3CB",
        "gene_A": "PIK3CA",
        "gene_B": "PIK3CB",
        "in_j0": 1,
        "in_k4": 1,
        "k4_role": "isoform_control",
        "paper_role": "too_close_for_primary; isoform control, not a dual optimum",
    },
    {
        "pair": "PIK3CA/mTOR",
        "gene_A": "PIK3CA",
        "gene_B": "MTOR",
        "in_j0": 1,
        "in_k4": 1,
        "k4_role": "development_cross_class",
        "paper_role": "literature PI3K/mTOR dual; pose-gold PI-103",
    },
]

COMPARATORS = [
    {
        "pair": "CREBBP/BRD4",
        "gene_A": "CREBBP",
        "gene_B": "BRD4",
        "in_j0": 0,
        "in_k4": 0,
        "k4_role": "not_in_candidate_list",
        "paper_role": "thickest conventional cross-class pair in this dump; never listed in J0",
    },
    {
        "pair": "HDAC1/HDAC6",
        "gene_A": "HDAC1",
        "gene_B": "HDAC6",
        "in_j0": 1,
        "in_k4": 0,
        "k4_role": "j0_thick_metal_excluded",
        "paper_role": "J0 thick pair correctly dropped as Zn-metal + isozyme",
    },
    {
        "pair": "PIK3CG/PIK3CB",
        "gene_A": "PIK3CG",
        "gene_B": "PIK3CB",
        "in_j0": 0,
        "in_k4": 0,
        "k4_role": "thicker_isoform_not_in_j0",
        "paper_role": "thicker PI3K isoform pair than PIK3CA/PIK3CB; not a J0 candidate",
    },
    {
        "pair": "CNR1/CNR2",
        "gene_A": "CNR1",
        "gene_B": "CNR2",
        "in_j0": 0,
        "in_k4": 0,
        "k4_role": "thicker_homolog_not_in_j0",
        "paper_role": "example: several GPCR/transporter homologs beat AChE/BChE on min HN",
    },
]


def load(name: str) -> list[dict]:
    with (TABLES / name).open(newline="") as fh:
        return list(csv.DictReader(fh))


def pair_key(gene_a: str, gene_b: str) -> frozenset[str]:
    return frozenset([gene_a, gene_b])


def find_rank(rows: list[dict], gene_a: str, gene_b: str) -> tuple[int | None, dict | None]:
    want = pair_key(gene_a, gene_b)
    for i, rec in enumerate(rows, 1):
        if pair_key(rec["gene_A"], rec["gene_B"]) == want:
            return i, rec
    return None, None


def verdict_for(spec: dict, thick_rank: int | None, rec: dict | None) -> str:
    if spec["in_k4"] and spec["k4_role"] == "case_supply_limited":
        return (
            "NOT a ChEMBL-wide thick optimum. Directional at θ=6.0 but min strict "
            "hard-neg = 7 (thick gate is 50). Kept as a pose-gold / supply-limited case."
        )
    if spec["in_k4"] and spec["k4_role"] == "isoform_control":
        return (
            "Thick, but not a dual-target optimum: J0 labelled too_close_for_primary. "
            "A thicker isoform pair (PIK3CG/PIK3CB) exists outside J0."
        )
    if spec["in_k4"] and spec["k4_role"] == "development_homolog":
        return (
            "Thick homolog dual with literature hybrids. Valid under J0+scientific "
            "filters; not rank-1 among homologs (CNR, SLC6, adenosine, opioid are thicker)."
        )
    if spec["in_k4"] and spec["k4_role"] == "development_cross_class":
        return (
            "Among the scarce conventional cross-class thick pairs, and the only such "
            "pair that was both on the J0 list and pose-gold. Not unique in the dump: "
            "CREBBP/BRD4 is thicker and was never listed."
        )
    if spec["k4_role"] == "not_in_candidate_list":
        return (
            "Candidate-list incompleteness, not a Table 2 bug. Do not dock in this paper."
        )
    if spec["k4_role"] == "j0_thick_metal_excluded":
        return "Correctly excluded from K=4 (metal). Not a missed dual seat."
    if rec is None:
        return "Not recovered in this table."
    return "Comparator only; not selected and not a docking expansion."


def main() -> int:
    thick = load("universe_pairs_strict_thick_annotated_v1.csv")
    directional = load("universe_pairs_directional_n10_all.csv")
    n_thick = len(thick)
    n_dir = len(directional)
    n_cross = sum(1 for r in thick if r["supply_bucket"] == "cross_class")
    rows_out = []
    for spec in K4 + COMPARATORS:
        t_rank, t_rec = find_rank(thick, spec["gene_A"], spec["gene_B"])
        d_rank, d_rec = find_rank(directional, spec["gene_A"], spec["gene_B"])
        rec = t_rec or d_rec
        min_hn = int(rec["min_strict_hardneg"]) if rec else ""
        n_both = int(rec["n_both_measured"]) if rec else ""
        bucket = (t_rec or {}).get("supply_bucket", "")
        same_class = (rec or {}).get("same_class", "")
        metal = (rec or {}).get("metal_either", "")
        directional_n10 = (rec or {}).get("directional_n10", "")
        thick_flag = int(t_rec is not None)
        rows_out.append(
            {
                "pair": spec["pair"],
                "gene_A": spec["gene_A"],
                "gene_B": spec["gene_B"],
                "in_j0_candidate_list": spec["in_j0"],
                "in_k4": spec["in_k4"],
                "k4_role": spec["k4_role"],
                "paper_role": spec["paper_role"],
                "universe_thick": thick_flag,
                "rank_among_thick": t_rank if t_rank is not None else "",
                "n_thick": n_thick,
                "rank_among_directional": d_rank if d_rank is not None else "",
                "n_directional": n_dir,
                "n_cross_class_thick": n_cross,
                "min_strict_hardneg": min_hn,
                "n_both_measured": n_both,
                "supply_bucket": bucket,
                "same_class": same_class,
                "metal_either": metal,
                "directional_n10": directional_n10,
                "is_chembl_wide_top4": 0,
                "suitability_verdict": verdict_for(spec, t_rank, rec),
            }
        )
    fields = list(rows_out[0].keys())
    out_path = TABLES / "k4_vs_universe_suitability_v1.csv"
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows_out)
    print(f"wrote {out_path} rows={len(rows_out)}")
    for rec in rows_out:
        print(
            f"  {rec['pair']:18s} k4={rec['in_k4']} thick={rec['universe_thick']} "
            f"rank_thick={rec['rank_among_thick'] or '-':>3} "
            f"minHN={rec['min_strict_hardneg']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
