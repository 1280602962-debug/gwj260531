#!/usr/bin/env python3
"""
Non-docking computational module C — physicochemical / drug-likeness profiling.

Computes RDKit-derived physicochemical descriptors and applies the standard
drug-likeness rule sets (Lipinski Ro5, Veber, Egan, Ghose) plus QED and a small
set of ADMET-relevant surrogate flags (Lipinski-style oral-absorption cues).
This is a pure downstream annotation of the EXISTING docking/Pareto data; it does
not re-dock or alter any docking / Pareto / ML score.

RDKit descriptors are a transparent, license-free surrogate. For submission,
these can be cross-checked against SwissADME / ADMETlab; the rule outcomes here
are deterministic and reproducible.

Inputs (read-only):
  data/repurposing/pareto/pareto_shortlist.csv
  data/repurposing/pareto/pareto_merged_scores.csv

Outputs:
  results/cheminformatics/admet_shortlist.csv
  results/cheminformatics/admet_pool.csv
  results/cheminformatics/admet_summary.json

Usage:
  python3 scripts/10_admet_druglikeness.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARETO_DIR = PROJECT_ROOT / "data" / "repurposing" / "pareto"
OUT_DIR = PROJECT_ROOT / "results" / "cheminformatics"

SMILES_COL = "canonical_smiles"
NAME_COL = "name"

KEEP_META = (
    "repurposing_id", "chembl_id", NAME_COL, SMILES_COL,
    "max_phase", "s_u_percentile", "s_n_percentile", "pareto_front",
)


def descriptors(mol) -> dict:
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    molmr = Crippen.MolMR(mol)
    n_atoms = mol.GetNumHeavyAtoms()
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    qed = Descriptors.qed(mol)

    # Rule sets
    ro5_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    lipinski_pass = ro5_violations <= 1
    veber_pass = rotb <= 10 and tpsa <= 140
    egan_pass = tpsa <= 131.6 and logp <= 5.88
    ghose_pass = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (20 <= n_atoms <= 70) and (40 <= molmr <= 130)
    # Oral-absorption surrogate (TPSA/rotB gate commonly used with Veber)
    oral_absorption_ok = tpsa <= 140 and rotb <= 10

    return {
        "mw": round(mw, 2),
        "clogp": round(logp, 2),
        "hbd": int(hbd),
        "hba": int(hba),
        "tpsa": round(tpsa, 2),
        "rot_bonds": int(rotb),
        "aromatic_rings": int(arom),
        "fsp3": round(fsp3, 3),
        "heavy_atoms": int(n_atoms),
        "molar_refractivity": round(molmr, 2),
        "qed": round(qed, 3),
        "ro5_violations": int(ro5_violations),
        "lipinski_pass": bool(lipinski_pass),
        "veber_pass": bool(veber_pass),
        "egan_pass": bool(egan_pass),
        "ghose_pass": bool(ghose_pass),
        "oral_absorption_ok": bool(oral_absorption_ok),
    }


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    from rdkit import Chem

    rows = []
    for _, r in df.iterrows():
        rec = {c: r.get(c) for c in df.columns if c in KEEP_META}
        smi = r.get(SMILES_COL)
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None:
            rec["parse_ok"] = False
            rows.append(rec)
            continue
        rec["parse_ok"] = True
        rec.update(descriptors(mol))
        rows.append(rec)
    return pd.DataFrame(rows)


def summarize(name: str, df: pd.DataFrame) -> dict:
    parsed = df[df["parse_ok"] == True]  # noqa: E712
    if not len(parsed):
        return {"set": name, "n": int(len(df)), "n_parsed": 0}
    return {
        "set": name,
        "n": int(len(df)),
        "n_parsed": int(len(parsed)),
        "median_qed": round(float(parsed["qed"].median()), 3),
        "n_lipinski_pass": int(parsed["lipinski_pass"].sum()),
        "n_veber_pass": int(parsed["veber_pass"].sum()),
        "n_egan_pass": int(parsed["egan_pass"].sum()),
        "n_ghose_pass": int(parsed["ghose_pass"].sum()),
        "n_oral_absorption_ok": int(parsed["oral_absorption_ok"].sum()),
        "n_all_rules_pass": int(
            (parsed["lipinski_pass"] & parsed["veber_pass"] & parsed["egan_pass"]).sum()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="ADMET / drug-likeness profiling (non-docking module C)")
    parser.add_argument("--shortlist", type=Path, default=PARETO_DIR / "pareto_shortlist.csv")
    parser.add_argument("--pool", type=Path, default=PARETO_DIR / "pareto_merged_scores.csv")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []

    shortlist = pd.read_csv(args.shortlist)
    ann_short = annotate(shortlist)
    ann_short.to_csv(args.output_dir / "admet_shortlist.csv", index=False)
    summaries.append(summarize("pareto_shortlist", ann_short))
    print(f"Shortlist profiled: {len(ann_short)} rows -> admet_shortlist.csv")

    pool = pd.read_csv(args.pool)
    ann_pool = annotate(pool)
    ann_pool.to_csv(args.output_dir / "admet_pool.csv", index=False)
    summaries.append(summarize("dual_docked_pool", ann_pool))
    print(f"Pool profiled: {len(ann_pool)} rows -> admet_pool.csv")

    with open(args.output_dir / "admet_summary.json", "w") as f:
        json.dump({"module": "C_admet_druglikeness", "sets": summaries}, f, indent=2)

    print("\n=== Drug-likeness summary ===")
    for s in summaries:
        if s.get("n_parsed"):
            print(
                f"  {s['set']}: n={s['n']} medQED={s['median_qed']} "
                f"Lipinski={s['n_lipinski_pass']} Veber={s['n_veber_pass']} "
                f"Egan={s['n_egan_pass']} all3={s['n_all_rules_pass']}"
            )

    print("\nShortlist detail:")
    cols = [NAME_COL, "max_phase", "mw", "clogp", "tpsa", "hbd", "hba",
            "qed", "lipinski_pass", "veber_pass", "oral_absorption_ok"]
    print(ann_short[[c for c in cols if c in ann_short.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
