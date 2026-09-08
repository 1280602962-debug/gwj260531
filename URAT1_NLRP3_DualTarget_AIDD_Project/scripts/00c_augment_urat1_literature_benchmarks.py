#!/usr/bin/env python3
"""
Augment URAT1 training set with literature benchmark activities (Burns 2016 / Nakamura 2019).

Purpose: recover scaffold-novel drugs (lesinurad, dotinurad) dropped by ChEMBL conflict curation.
This is an SI ablation — not the primary URAT1 evidence path (docking remains primary).

Output:
  data/processed/urat1_curated_literature_augmented.csv
  data/processed/urat1_augmentation_manifest.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from utils_ml import canonicalize, murcko_scaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks_summary.csv"
CURATED = PROJECT_ROOT / "data" / "processed" / "urat1_curated.csv"

# Burns 2016 HEK-URAT1 primary references used in project benchmarks
LITERATURE_ROWS = [
    {
        "compound_name": "lesinurad",
        "canonical_smiles": "O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12",
        "pActivity": 8.45,
        "ref_pmid": "27716403",
        "source": "literature_burns_2016",
    },
    {
        "compound_name": "dotinurad",
        "canonical_smiles": "O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21",
        "pActivity": 10.43,
        "ref_pmid": "31371478",
        "source": "literature_nakamura_2019",
    },
    {
        "compound_name": "benzbromarone",
        "canonical_smiles": "CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1",
        "pActivity": 9.54,
        "ref_pmid": "27716403",
        "source": "literature_burns_2016",
    },
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=CURATED)
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "urat1_curated_literature_augmented.csv",
    )
    args = parser.parse_args()

    base = pd.read_csv(args.input)
    existing = set(base["canonical_smiles"])

    add_rows = []
    for r in LITERATURE_ROWS:
        smi = canonicalize(r["canonical_smiles"])
        if smi in existing:
            continue
        row = {c: None for c in base.columns}
        row.update(
            {
                "canonical_smiles": smi,
                "scaffold": murcko_scaffold(smi),
                "pActivity": r["pActivity"],
                "molecule_chembl_id": None,
                "molecule_name": r["compound_name"],
                "n_measurements": 1,
            }
        )
        add_rows.append(row)

    out = pd.concat([base, pd.DataFrame(add_rows)], ignore_index=True)
    out = out.drop_duplicates("canonical_smiles", keep="first")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    manifest = {
        "base_n": int(len(base)),
        "added_n": int(len(add_rows)),
        "output_n": int(len(out)),
        "added_compounds": [r["compound_name"] for r in LITERATURE_ROWS],
        "note": "Retrain with --data-dir or copy to urat1_curated.csv only for ablation; docking remains primary.",
    }
    manifest_path = args.output.parent / "urat1_augmentation_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
