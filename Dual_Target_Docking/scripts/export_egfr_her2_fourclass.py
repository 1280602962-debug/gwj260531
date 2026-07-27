#!/usr/bin/env python3
"""Export EGFR/HER2 four-class molecule list from cached ChEMBL pChEMBL maps.

Definition (operational "true dual"):
  both targets measured on the same molecule_chembl_id, and
  max pChEMBL(EGFR) >= theta AND max pChEMBL(HER2) >= theta.

Default theta = 6.0 (≈ 1 µM). Untested ≠ inactive: A_only/B_only require
the weak end to be measured below theta.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "public_pair_selection"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", type=float, default=6.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=CACHE / "egfr_her2_fourclass_chembl_ids.csv",
    )
    args = ap.parse_args()

    egfr = {k: float(v) for k, v in json.load(open(CACHE / "mols_EGFR.json")).items()}
    her2 = {k: float(v) for k, v in json.load(open(CACHE / "mols_HER2.json")).items()}
    both = set(egfr) & set(her2)

    rows = []
    for mid in both:
        pa, pb = egfr[mid], her2[mid]
        if pa >= args.theta and pb >= args.theta:
            cls = "dual"
        elif pa >= args.theta and pb < args.theta:
            cls = "A_only"
        elif pb >= args.theta and pa < args.theta:
            cls = "B_only"
        else:
            cls = "neither"
        rows.append((mid, cls, pa, pb, min(pa, pb)))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "molecule_chembl_id",
                "class",
                "pchembl_EGFR",
                "pchembl_HER2",
                "min_pchembl",
                "theta",
            ]
        )
        for mid, cls, pa, pb, mn in sorted(rows, key=lambda x: (x[1], -x[4])):
            w.writerow([mid, cls, f"{pa:.2f}", f"{pb:.2f}", f"{mn:.2f}", args.theta])

    from collections import Counter

    c = Counter(r[1] for r in rows)
    print(f"wrote {args.out}  paired={len(rows)}  {dict(c)}")
    print(
        "Next: fetch SMILES via ChEMBL molecule API / UniChem for docking panel."
    )


if __name__ == "__main__":
    main()
