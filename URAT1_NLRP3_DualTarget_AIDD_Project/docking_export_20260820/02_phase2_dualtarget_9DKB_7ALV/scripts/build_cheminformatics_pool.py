#!/usr/bin/env python3
"""Build Module B/C cheminformatics tables for the dual-docked P0.5 pool.

Missing from the server upload pack; regenerated from pareto_merged_scores.csv.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _alert_flags(mol: Chem.Mol) -> dict:
    out = {"pains_any": False, "brenk": False, "nih": False}
    for key, enum in (
        ("pains_any", FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS),
        ("brenk", FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK),
        ("nih", FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH),
    ):
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(enum)
        cat = FilterCatalog.FilterCatalog(params)
        out[key] = bool(cat.HasMatch(mol))
    # Aggregation heuristic: large aromatic + high logP
    clogp = Crippen.MolLogP(mol)
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    out["aggregation_risk_heuristic"] = bool(arom >= 4 and clogp >= 5.0)
    return out


def _row_props(mol: Chem.Mol) -> dict:
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    heavy = Lipinski.HeavyAtomCount(mol)
    qed = float(QED.qed(mol))
    lipinski = (mw <= 500) and (clogp <= 5) and (hbd <= 5) and (hba <= 10)
    veber = (rotb <= 10) and (tpsa <= 140)
    # Ghose: 160<=MW<=480, -0.4<=logP<=5.6, 20<=atoms<=70, 40<=MR<=130
    mr = Crippen.MolMR(mol)
    ghose = (160 <= mw <= 480) and (-0.4 <= clogp <= 5.6) and (20 <= heavy <= 70) and (40 <= mr <= 130)
    oral_absorption_ok = veber and (mw <= 550) and (clogp <= 5.5)
    try:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        scaffold = ""
    return {
        "qed": qed,
        "mw": mw,
        "clogp": clogp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "heavy_atoms": heavy,
        "lipinski_pass": bool(lipinski),
        "veber_pass": bool(veber),
        "ghose_pass": bool(ghose),
        "oral_absorption_ok": bool(oral_absorption_ok),
        "scaffold_rdkit": scaffold,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pool",
        type=Path,
        default=PROJECT_ROOT / "data/repurposing/pareto/pareto_merged_scores.csv",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results/cheminformatics",
    )
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pool)
    filter_rows = []
    admet_rows = []
    for _, r in df.iterrows():
        name = r.get("name")
        smi = r.get("canonical_smiles")
        mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
        if mol is None:
            filter_rows.append(
                {
                    "name": name,
                    "pains_any": True,
                    "brenk": True,
                    "nih": True,
                    "aggregation_risk_heuristic": True,
                }
            )
            admet_rows.append(
                {
                    "name": name,
                    "qed": None,
                    "mw": r.get("mw"),
                    "clogp": None,
                    "tpsa": None,
                    "hbd": None,
                    "hba": None,
                    "heavy_atoms": None,
                    "lipinski_pass": False,
                    "veber_pass": False,
                    "ghose_pass": False,
                    "oral_absorption_ok": False,
                }
            )
            continue
        alerts = _alert_flags(mol)
        props = _row_props(mol)
        filter_rows.append({"name": name, **alerts})
        admet_rows.append(
            {
                "name": name,
                "qed": props["qed"],
                "mw": props["mw"],
                "clogp": props["clogp"],
                "tpsa": props["tpsa"],
                "hbd": props["hbd"],
                "hba": props["hba"],
                "heavy_atoms": props["heavy_atoms"],
                "lipinski_pass": props["lipinski_pass"],
                "veber_pass": props["veber_pass"],
                "ghose_pass": props["ghose_pass"],
                "oral_absorption_ok": props["oral_absorption_ok"],
            }
        )

    f = pd.DataFrame(filter_rows)
    a = pd.DataFrame(admet_rows)
    f.to_csv(args.output_dir / "filters_pool.csv", index=False)
    a.to_csv(args.output_dir / "admet_pool.csv", index=False)
    print(f"wrote {args.output_dir/'filters_pool.csv'} n={len(f)}")
    print(f"wrote {args.output_dir/'admet_pool.csv'} n={len(a)}")
    print(
        "alert rates:",
        f"PAINS={f.pains_any.mean():.3f}",
        f"Brenk={f.brenk.mean():.3f}",
        f"Lipinski={a.lipinski_pass.mean():.3f}",
        f"Veber={a.veber_pass.mean():.3f}",
    )


if __name__ == "__main__":
    main()
