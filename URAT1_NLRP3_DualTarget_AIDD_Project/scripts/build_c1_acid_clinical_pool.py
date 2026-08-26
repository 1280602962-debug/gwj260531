#!/usr/bin/env python3
"""Build C1 Acid-track clinical acid pool from qN>=0.5 docking pool.

Gates (campaign_c1.yaml / Amendment A1):
  - carboxylate or tetrazole (acid equivalent)
  - exclude Vecabrutinib and other non-acids from URAT1 claim
  - keep qN>=0.5 pool membership (already in docking_pool_p05)
  - optional chemistry flags recorded, not used to delete yet
No percentile ranking. Output under data/campaigns/c1/07_clinical_dock/acid_pool/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors

PROJECT_ROOT = Path(__file__).resolve().parents[1]

COOH = Chem.MolFromSmarts("[CX3](=O)[OH,O-]")
TETRAZOLE = Chem.MolFromSmarts("c1nnn[nH]1")
ACYL_SULFONAMIDE = Chem.MolFromSmarts("[CX3](=O)N[SX4](=O)(=O)")


def acid_flags(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {"ok": False, "reason": "parse_fail"}
    has_cooh = bool(COOH and m.HasSubstructMatch(COOH))
    has_tet = bool(TETRAZOLE and m.HasSubstructMatch(TETRAZOLE))
    has_asu = bool(ACYL_SULFONAMIDE and m.HasSubstructMatch(ACYL_SULFONAMIDE))
    if not (has_cooh or has_tet or has_asu):
        return {"ok": False, "reason": "no_acid_equivalent", "mol": m}
    mw = Descriptors.MolWt(m)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)
    logp = Crippen.MolLogP(m)
    rot = Lipinski.NumRotatableBonds(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    veber = tpsa <= 140 and rot <= 10
    ro5_hb = hbd <= 5 and hba <= 10
    mw_ok = 200 <= mw <= 550
    return {
        "ok": True,
        "has_carboxylate": has_cooh,
        "has_tetrazole": has_tet,
        "has_acyl_sulfonamide": has_asu,
        "mw": mw,
        "hbd": hbd,
        "hba": hba,
        "logp": logp,
        "rotbonds": rot,
        "tpsa": tpsa,
        "pass_veber": veber,
        "pass_ro5_hb": ro5_hb,
        "pass_mw_200_550": mw_ok,
        "pass_chemistry_soft": veber and ro5_hb and mw_ok,
        "mol": m,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT / "data/repurposing/screening/docking_pool_p05.csv",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_pool",
    )
    ap.add_argument("--exclude-name-substr", nargs="*", default=["VECABRUTINIB"])
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    rows = []
    for _, r in df.iterrows():
        name = str(r.get("name", ""))
        if any(x.upper() in name.upper() for x in args.exclude_name_substr):
            continue
        flags = acid_flags(str(r["canonical_smiles"]))
        if not flags.get("ok"):
            continue
        flags.pop("mol", None)
        rows.append(
            {
                "repurposing_id": r["repurposing_id"],
                "chembl_id": r.get("chembl_id"),
                "name": name,
                "canonical_smiles": r["canonical_smiles"],
                "p_active_nlrp3": r.get("p_active_nlrp3"),
                "mw_pool": r.get("mw"),
                **flags,
            }
        )

    out = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "acid_clinical_pool.csv"
    out.to_csv(csv_path, index=False)

    # soft-chemistry subset for primary Acid shortlist docking first
    soft = out[out["pass_chemistry_soft"]].copy()
    soft_path = args.output_dir / "acid_clinical_pool_chemistry_pass.csv"
    soft.to_csv(soft_path, index=False)

    summary = {
        "input": str(args.input),
        "n_input": int(len(df)),
        "n_acid_equivalent": int(len(out)),
        "n_chemistry_soft_pass": int(len(soft)),
        "n_carboxylate": int(out["has_carboxylate"].sum()),
        "n_tetrazole": int(out["has_tetrazole"].sum()),
        "n_acyl_sulfonamide": int(out["has_acyl_sulfonamide"].sum()),
        "exclude_name_substr": args.exclude_name_substr,
        "acid_pool_csv": str(csv_path),
        "chemistry_pass_csv": str(soft_path),
        "note": "Acid track: no docking percentile ranking. Geometry+chemistry gates only.",
    }
    (args.output_dir / "acid_pool_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
