#!/usr/bin/env python3
"""Prepare ligand PDBQT files from SMILES (RDKit + Meeko) for Vina docking."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from utils_ml import canonicalize

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def smiles_to_pdbqt(smiles: str, mol_id: str) -> tuple[str | None, str | None]:
    smi = canonicalize(smiles)
    if not smi:
        return None, "invalid_smiles"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, "rdkit_parse_fail"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None, "embed_fail"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return None, "meeko_prep_fail"
    pdbqt = PDBQTWriterLegacy.write_string(setups[0])
    if isinstance(pdbqt, tuple):
        pdbqt = pdbqt[0]
    return pdbqt, None


def prepare_pool(
    input_csv: Path,
    output_dir: Path,
    smiles_col: str = "canonical_smiles",
    id_col: str = "repurposing_id",
) -> dict:
    df = pd.read_csv(input_csv, low_memory=False)
    if smiles_col not in df.columns:
        raise ValueError(f"{input_csv} missing {smiles_col}")
    if id_col not in df.columns:
        df[id_col] = [f"LIG_{i:05d}" for i in range(len(df))]

    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, row in df.iterrows():
        rid = str(row[id_col])
        smi = row[smiles_col]
        pdbqt, err = smiles_to_pdbqt(str(smi), rid)
        out_path = output_dir / f"{rid}.pdbqt"
        status = "prepared"
        if pdbqt is None:
            status = err or "fail"
            rows.append({"repurposing_id": rid, "canonical_smiles": smi, "status": status, "pdbqt": None})
            continue
        out_path.write_text(pdbqt)
        rows.append(
            {
                "repurposing_id": rid,
                "canonical_smiles": smi,
                "status": status,
                "pdbqt": str(out_path),
            }
        )

    manifest = pd.DataFrame(rows)
    manifest_path = output_dir / "ligand_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    summary = {
        "input": str(input_csv),
        "output_dir": str(output_dir),
        "n_total": int(len(manifest)),
        "n_prepared": int((manifest["status"] == "prepared").sum()),
        "n_failed": int((manifest["status"] != "prepared").sum()),
        "manifest": str(manifest_path),
    }
    (output_dir / "prepare_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare ligand PDBQT library")
    parser.add_argument("--input", type=Path, required=True, help="CSV with canonical_smiles")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--smiles-col", type=str, default="canonical_smiles")
    parser.add_argument("--id-col", type=str, default="repurposing_id")
    args = parser.parse_args()
    summary = prepare_pool(args.input, args.output_dir, args.smiles_col, args.id_col)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
