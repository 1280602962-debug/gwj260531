#!/usr/bin/env python3
"""Extract NLRP3 (7ALV) key-residue map for Acid-track interaction fingerprints.

UniProt / Dekker 2021 key residues:
  Ala227, Ala228, Arg351, Met408, Tyr443, Phe575, Arg578
Mapped onto prepared receptor atoms via Cα matching to 7ALV.pdb.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]

NLRP3_KEY = {
    "ALA227": "Walker A backbone (sulfonylurea)",
    "ALA228": "Walker A polar contact (urea)",
    "ARG351": "NBD clamp on sulfonyl",
    "MET408": "HD1 hydrophobic wing",
    "TYR443": "WHD hydrophobic wing",
    "PHE575": "HD2 hydrophobic",
    "ARG578": "HD2 clamp (critical H-bond)",
}


def load_pdb_atoms(path: Path, hetatm: bool = True) -> list[dict]:
    atoms = []
    with path.open() as f:
        for line in f:
            if line.startswith("ATOM") or (hetatm and line.startswith("HETATM")):
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21]
                try:
                    resi = int(line[22:26])
                except ValueError:
                    continue
                xyz = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
                elem = (line[76:78].strip() if len(line) >= 78 else name[0]).upper()
                if not elem:
                    elem = name[0]
                atoms.append(
                    {
                        "name": name,
                        "resn": resn,
                        "chain": chain,
                        "resi": resi,
                        "xyz": xyz.tolist(),
                        "elem": elem,
                        "is_h": elem == "H" or name.startswith("H"),
                    }
                )
    return atoms


def load_pdbqt_atoms(path: Path) -> list[dict]:
    """Approximate ATOM records from AutoDock PDBQT (same columns as PDB for coords)."""
    atoms = []
    with path.open() as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21] if len(line) > 21 else "A"
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            xyz = [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]
            # PDBQT atom type at end; element from name
            elem = name[0].upper()
            if name[:2].upper() in {"CL", "BR", "FE", "ZN", "MG"}:
                elem = name[:2].upper().title()
            atoms.append(
                {
                    "name": name,
                    "resn": resn,
                    "chain": chain,
                    "resi": resi,
                    "xyz": xyz,
                    "elem": elem,
                    "is_h": elem == "H" or name.startswith("H") or name.upper().startswith("HD"),
                }
            )
    return atoms


def map_key_residues(src_atoms: list[dict], prep_atoms: list[dict]) -> dict:
    want = {227, 228, 351, 408, 443, 575, 578}
    expected = {
        227: "ALA",
        228: "ALA",
        351: "ARG",
        408: "MET",
        443: "TYR",
        575: "PHE",
        578: "ARG",
    }
    src_ca = {
        a["resi"]: a
        for a in src_atoms
        if a["name"] == "CA" and a["chain"] == "A" and a["resi"] in want
    }
    prep_ca = [a for a in prep_atoms if a["name"] == "CA"]
    mapped = {}
    for uni, src in src_ca.items():
        if src["resn"] != expected[uni]:
            continue
        best, bd = None, 1e9
        sx = np.array(src["xyz"])
        for p in prep_ca:
            d = float(np.linalg.norm(np.array(p["xyz"]) - sx))
            if d < bd:
                bd, best = d, p
        if best is None or bd >= 1.5:
            continue
        # collect all prep atoms of that residue
        res_atoms = [
            a
            for a in prep_atoms
            if a["chain"] == best["chain"] and a["resi"] == best["resi"] and not a["is_h"]
        ]
        label = f"{expected[uni]}{uni}"
        mapped[label] = {
            "uniprot": uni,
            "prep_chain": best["chain"],
            "prep_resi": best["resi"],
            "prep_resn": best["resn"],
            "ca_match_d_A": bd,
            "role": NLRP3_KEY[label],
            "heavy_atoms": [{"name": a["name"], "elem": a["elem"], "xyz": a["xyz"]} for a in res_atoms],
        }
    return mapped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--src-pdb",
        type=Path,
        default=PROJECT_ROOT / "data/structures/pdb/7ALV.pdb",
    )
    ap.add_argument(
        "--prep-pdbqt",
        type=Path,
        default=PROJECT_ROOT / "data/structures/prepared/7ALV_receptor.pdbqt",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT
        / "data/campaigns/c1/01_ligand_prep/selfdock_refs/nlrp3_key_residues.json",
    )
    args = ap.parse_args()
    src = load_pdb_atoms(args.src_pdb)
    prep = load_pdbqt_atoms(args.prep_pdbqt)
    mapped = map_key_residues(src, prep)
    out = {
        "pdb": "7ALV",
        "ligand_ref": "NP3-146/RM5",
        "n_key_expected": 7,
        "n_key_mapped": len(mapped),
        "residues": mapped,
        "note": "Contacts use prepared-receptor atom coords; UniProt labels via Cα match to 7ALV.pdb",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"n_mapped": len(mapped), "labels": list(mapped.keys())}, indent=2))
    if len(mapped) < 7:
        raise SystemExit(f"Only mapped {len(mapped)}/7 key residues")


if __name__ == "__main__":
    main()
