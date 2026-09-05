#!/usr/bin/env python3
"""
Build mol_XXXXX -> unique_docking_pool row map from LigPrep SDF properties.
Does not rewrite geometries. Source index is 1-based into unique_docking_pool.csv.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdf", required=True, help="ligands_all.sdf from structconvert")
    ap.add_argument("--pool", required=True, help="unique_docking_pool.csv")
    ap.add_argument("--out", required=True, help="output CSV")
    args = ap.parse_args()

    pool = list(csv.DictReader(open(args.pool, newline="")))
    rows = []
    suppl = Chem.SDMolSupplier(args.sdf, removeHs=False, sanitize=False)
    for i, mol in enumerate(suppl):
        mol_id = f"mol_{i:05d}"
        if mol is None:
            rows.append(
                {
                    "idx": i,
                    "mol_id": mol_id,
                    "source_file_index": "",
                    "canonical_smiles": "",
                    "role": "",
                    "in_true": "",
                    "in_random": "",
                    "status": "sdf_null",
                }
            )
            continue
        src = ""
        for key in ("i_m_Source_File_Index", "i_m_source_file_index"):
            if mol.HasProp(key):
                src = mol.GetProp(key).strip()
                break
        smiles = ""
        role = in_true = in_random = ""
        status = "ok"
        if src.isdigit():
            j = int(src) - 1  # 1-based
            if 0 <= j < len(pool):
                smiles = pool[j]["canonical_smiles"]
                role = pool[j].get("role", "")
                in_true = pool[j].get("in_true", "")
                in_random = pool[j].get("in_random", "")
            else:
                status = "src_oob"
        else:
            status = "no_src"
            try:
                Chem.SanitizeMol(mol)
                smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
            except Exception:
                pass
        rows.append(
            {
                "idx": i,
                "mol_id": mol_id,
                "source_file_index": src,
                "canonical_smiles": smiles,
                "role": role,
                "in_true": in_true,
                "in_random": in_random,
                "title": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
                "status": status,
            }
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "idx",
        "mol_id",
        "source_file_index",
        "canonical_smiles",
        "role",
        "in_true",
        "in_random",
        "title",
        "status",
    ]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"wrote {out} n={len(rows)} ok={ok}")


if __name__ == "__main__":
    main()
