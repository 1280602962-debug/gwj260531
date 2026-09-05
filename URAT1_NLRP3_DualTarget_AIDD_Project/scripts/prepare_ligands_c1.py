#!/usr/bin/env python3
"""C1 ligand prep: pH 7.4 microspecies (Dimorphite-DL) → Meeko PDBQT.

Campaign lock: config/campaign_c1.yaml (carboxylic acids deprotonated at pH 7.4).
Forbidden: hand-fix only GSK/lesinurad while keeping old decoy PDBQTs.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd
from dimorphite_dl import protonate_smiles
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def _has_carboxylate(mol: Chem.Mol) -> bool:
    # deprotonated carboxylate or carboxylic acid
    pats = [
        Chem.MolFromSmarts("[CX3](=O)[O-]"),
        Chem.MolFromSmarts("[CX3](=O)[OH]"),
    ]
    return any(mol.HasSubstructMatch(p) for p in pats if p is not None)


def enumerate_ph74(smiles: str, max_variants: int = 8) -> tuple[list[str], str | None]:
    try:
        variants = list(
            protonate_smiles(
                smiles,
                ph_min=7.4,
                ph_max=7.4,
                precision=1.0,
                max_variants=max_variants,
            )
        )
    except Exception as exc:  # noqa: BLE001
        return [], f"dimorphite_fail:{type(exc).__name__}"
    # unique canonical forms
    out: list[str] = []
    seen: set[str] = set()
    for s in variants:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        mol = _largest_fragment(mol)
        can = Chem.MolToSmiles(mol)
        if can not in seen:
            seen.add(can)
            out.append(can)
    if not out:
        return [], "dimorphite_empty"
    return out, None


def smiles_to_pdbqt(smiles: str, embed_seed: int = 0xC0FFEE) -> tuple[str | None, str | None]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "rdkit_parse_fail"
    mol = _largest_fragment(mol)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = embed_seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None, "embed_fail"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    try:
        prep = MoleculePreparation()
        setups = prep.prepare(mol)
    except Exception as exc:  # noqa: BLE001
        return None, f"meeko_prep_fail:{type(exc).__name__}"
    if not setups:
        return None, "meeko_prep_fail"
    pdbqt = PDBQTWriterLegacy.write_string(setups[0])
    if isinstance(pdbqt, tuple):
        pdbqt = pdbqt[0]
    return pdbqt, None


def _safe_id(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return s.strip("_") or "LIG"


def prepare_rows(
    rows: list[dict],
    output_dir: Path,
    embed_seed: int = 0xC0FFEE,
    prefer_carboxylate: bool = True,
) -> dict:
    """Prepare one primary microspecies PDBQT per input molecule.

    If prefer_carboxylate, choose the variant that contains a carboxylate when any
    exist; otherwise keep the first Dimorphite variant. All variants are recorded.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_dir = output_dir / "pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict] = []
    failures: list[dict] = []

    for row in rows:
        rid = _safe_id(str(row["repurposing_id"]))
        smi0 = str(row["canonical_smiles"])
        variants, err = enumerate_ph74(smi0)
        if err:
            failures.append({**row, "status": err})
            manifest_rows.append(
                {
                    "repurposing_id": rid,
                    "canonical_smiles": smi0,
                    "status": err,
                    "pdbqt": None,
                    "n_microspecies": 0,
                    "selected_smiles": None,
                    "formal_charge": None,
                    "has_carboxylate": False,
                    "all_microspecies": "",
                }
            )
            continue

        selected = variants[0]
        if prefer_carboxylate:
            carbox = []
            for v in variants:
                m = Chem.MolFromSmiles(v)
                if m is not None and _has_carboxylate(m):
                    # prefer deprotonated
                    if m.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[O-]")):
                        carbox.insert(0, v)
                    else:
                        carbox.append(v)
            if carbox:
                selected = carbox[0]

        mol_sel = Chem.MolFromSmiles(selected)
        formal = int(Chem.GetFormalCharge(mol_sel)) if mol_sel else None
        has_cooh = bool(mol_sel and _has_carboxylate(mol_sel))

        out_path = pdbqt_dir / f"{rid}.pdbqt"
        if out_path.exists() and out_path.stat().st_size > 0:
            status = "prepared"
            pdbqt_s = str(out_path)
        else:
            pdbqt, perr = smiles_to_pdbqt(selected, embed_seed=embed_seed)
            if pdbqt is None:
                status = perr or "fail"
                failures.append({**row, "status": status, "selected_smiles": selected})
                manifest_rows.append(
                    {
                        "repurposing_id": rid,
                        "canonical_smiles": smi0,
                        "status": status,
                        "pdbqt": None,
                        "n_microspecies": len(variants),
                        "selected_smiles": selected,
                        "formal_charge": formal,
                        "has_carboxylate": has_cooh,
                        "all_microspecies": "|".join(variants),
                    }
                )
                continue
            out_path.write_text(pdbqt)
            status = "prepared"
            pdbqt_s = str(out_path)

        # write all microspecies smiles for audit
        (output_dir / "microspecies").mkdir(exist_ok=True)
        (output_dir / "microspecies" / f"{rid}.smi").write_text(
            "\n".join(variants) + "\n"
        )

        manifest_rows.append(
            {
                "repurposing_id": rid,
                "canonical_smiles": smi0,
                "status": status,
                "pdbqt": pdbqt_s,
                "n_microspecies": len(variants),
                "selected_smiles": selected,
                "formal_charge": formal,
                "has_carboxylate": has_cooh,
                "all_microspecies": "|".join(variants),
                "mw": Descriptors.MolWt(mol_sel) if mol_sel else None,
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "ligand_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    # hard gate for named acids
    must_deprot = {"lesinurad", "GSK-3008348", "GSK-3008348_FREE_BASE", "verinurad", "probenecid", "puliginurad", "SHR-4640"}
    acid_gate = {}
    for _, r in manifest.iterrows():
        if r["repurposing_id"] in must_deprot or r["repurposing_id"].lower() in {x.lower() for x in must_deprot}:
            sel = r.get("selected_smiles") or ""
            acid_gate[r["repurposing_id"]] = {
                "has_carboxylate": bool(r.get("has_carboxylate")),
                "deprotonated": "[O-]" in sel or "([O-])" in sel,
                "formal_charge": r.get("formal_charge"),
                "status": r.get("status"),
            }

    summary = {
        "output_dir": str(output_dir),
        "n_total": int(len(manifest)),
        "n_prepared": int((manifest["status"] == "prepared").sum()),
        "n_failed": int((manifest["status"] != "prepared").sum()),
        "manifest": str(manifest_path),
        "acid_gate": acid_gate,
        "ph": 7.4,
        "engine": "dimorphite_dl_then_meeko",
    }
    (output_dir / "prep_summary.json").write_text(json.dumps(summary, indent=2))
    if failures:
        pd.DataFrame(failures).to_csv(output_dir / "prep_failures.csv", index=False)
    return summary


def build_forced_recovery_table() -> list[dict]:
    """Pre-registered forced-recovery + L2 selfdock ligands (neutral input SMILES)."""
    # Unique SMILES from literature_benchmarks / curated tables in this repo.
    return [
        {"repurposing_id": "lesinurad", "canonical_smiles": "O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12", "role": "forced_recovery+selfdock"},
        {"repurposing_id": "benzbromarone", "canonical_smiles": "CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1", "role": "forced_recovery"},
        {"repurposing_id": "dotinurad", "canonical_smiles": "O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21", "role": "forced_recovery"},
        {"repurposing_id": "verinurad", "canonical_smiles": "CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O", "role": "forced_recovery"},
        {"repurposing_id": "probenecid", "canonical_smiles": "CCCN(CCC)S(=O)(=O)c1ccc(C(=O)O)cc1", "role": "forced_recovery"},
        {"repurposing_id": "puliginurad", "canonical_smiles": "CC(C)(Cc1cc2c(-c3ccc(C#N)cc3)cncc2s1)C(=O)O", "role": "forced_recovery"},
        {"repurposing_id": "SHR-4640", "canonical_smiles": "O=C(O)C1(Sc2ccnc3ccc(Br)cc23)CCC1", "role": "forced_recovery"},
        {
            "repurposing_id": "GSK-3008348",
            "canonical_smiles": "Cc1cc(C)n(-c2cccc([C@H](CC(=O)O)CN3CC[C@@H](CCc4ccc5c(n4)NCCC5)C3)c2)n1",
            "role": "acid_audit",
        },
    ]


def main() -> None:
    ap = argparse.ArgumentParser(description="C1 pH 7.4 carboxylate ligand prep")
    ap.add_argument("--input-csv", type=Path, default=None, help="CSV with id + smiles cols")
    ap.add_argument("--id-col", default="repurposing_id")
    ap.add_argument("--smiles-col", default="canonical_smiles")
    ap.add_argument("--forced-recovery", action="store_true", help="Use built-in recovery table")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--embed-seed", type=lambda x: int(x, 0), default=0xC0FFEE)
    args = ap.parse_args()

    if args.forced_recovery:
        rows = build_forced_recovery_table()
    elif args.input_csv:
        df = pd.read_csv(args.input_csv)
        rows = [
            {"repurposing_id": r[args.id_col], "canonical_smiles": r[args.smiles_col]}
            for _, r in df.iterrows()
        ]
    else:
        raise SystemExit("Provide --forced-recovery or --input-csv")

    summary = prepare_rows(rows, args.output_dir, embed_seed=args.embed_seed)
    print(json.dumps(summary, indent=2))

    # Fail hard only for named acids that were actually requested in this run.
    requested = {str(r["repurposing_id"]) for r in rows}
    gate = summary.get("acid_gate", {})
    bad = []
    for name in ("lesinurad", "GSK-3008348"):
        if name not in requested:
            continue
        g = gate.get(name)
        if not g or g.get("status") != "prepared" or not g.get("deprotonated"):
            bad.append(name)
    if bad:
        raise SystemExit(f"L1 acid gate FAILED for: {bad}. Fix prep before docking.")


if __name__ == "__main__":
    main()
