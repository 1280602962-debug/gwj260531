#!/usr/bin/env python3
"""Extract URAT1 (9DKB) literature key-residue map for C5 W2 IFP gate.

Literature labels (campaign_c5.yaml):
  S35, M214, F241, F360, F364, F365, D389, K393, Q437, F449, R477, Q473

Mapped onto prepared 9DKB PDBQT via Cα match to 9DKB.cif (same protocol as
extract_nlrp3_key_residues.py). Prepared numbering is typically lit-1
(e.g. Arg477 → ARG A 476).

Q437 is LEU in 9DKB chain A — recorded as unmappable, not invented.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import gemmi
import numpy as np

from extract_nlrp3_key_residues import load_pdbqt_atoms

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# literature / UniProt-style labels used in campaign docs
URAT1_KEY = {
    35: ("SER", "S35", "access path / Tan 2017"),
    214: ("MET", "M214", "pocket wall"),
    241: ("PHE", "F241", "Phe cage (Dai & Lee 2024)"),
    360: ("PHE", "F360", "Phe cage"),
    364: ("PHE", "F364", "Phe cage"),
    365: ("PHE", "F365", "Phe cage"),
    389: ("ASP", "D389", "Guo/Chen 2025 mutant"),
    393: ("LYS", "K393", "pocket electrostatics"),
    437: ("GLN", "Q437", "literature label — may not match 9DKB sequence"),
    449: ("PHE", "F449", "Phe cage"),
    473: ("GLN", "Q473", "near Arg477"),
    477: ("ARG", "R477", "acid anchor (Tan / Guo)"),
}


def cif_chain_a_residues(path: Path) -> dict[int, gemmi.Residue]:
    st = gemmi.read_structure(str(path))
    out: dict[int, gemmi.Residue] = {}
    for res in st[0]["A"]:
        try:
            resi = int(res.seqid.num)
        except Exception:
            continue
        out[resi] = res
    return out


def ca_xyz(res: gemmi.Residue) -> np.ndarray | None:
    for atom in res:
        if atom.name.strip() == "CA":
            return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
    return None


def map_key_residues(src_cif: Path, prep_pdbqt: Path) -> tuple[dict, dict]:
    src = cif_chain_a_residues(src_cif)
    prep = load_pdbqt_atoms(prep_pdbqt)
    prep_ca = [a for a in prep if a["name"] == "CA"]
    mapped: dict = {}
    unmatched: dict = {}
    for uni, (exp_resn, label, role) in URAT1_KEY.items():
        if uni not in src:
            unmatched[label] = {"reason": "missing_in_cif", "lit_resi": uni}
            continue
        res = src[uni]
        got = res.name.strip()
        if got != exp_resn:
            unmatched[label] = {
                "reason": "sequence_mismatch",
                "lit_resi": uni,
                "expected_resn": exp_resn,
                "cif_resn": got,
            }
            continue
        xyz = ca_xyz(res)
        if xyz is None:
            unmatched[label] = {"reason": "no_CA", "lit_resi": uni}
            continue
        best, bd = None, 1e9
        for a in prep_ca:
            d = float(np.linalg.norm(np.array(a["xyz"]) - xyz))
            if d < bd:
                bd, best = d, a
        if best is None or bd >= 1.5:
            unmatched[label] = {"reason": "ca_match_fail", "lit_resi": uni, "best_d": bd}
            continue
        res_atoms = [
            a
            for a in prep
            if a["chain"] == best["chain"] and a["resi"] == best["resi"] and not a["is_h"]
        ]
        # Key = RESN+lit (e.g. ARG477) so c1_nlrp3_pose_metrics IFP polarity
        # startswith("ARG"/"PHE"/...) still works; keep short label in metadata.
        key = f"{exp_resn}{uni}"
        mapped[key] = {
            "short_label": label,
            "uniprot_or_lit": uni,
            "prep_chain": best["chain"],
            "prep_resi": best["resi"],
            "prep_resn": best["resn"],
            "ca_match_d_A": bd,
            "role": role,
            "heavy_atoms": [
                {"name": a["name"], "elem": a["elem"], "xyz": a["xyz"]} for a in res_atoms
            ],
        }
    return mapped, unmatched


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--src-cif",
        type=Path,
        default=PROJECT_ROOT / "data/structures/pdb/9DKB.cif",
    )
    ap.add_argument(
        "--prep-pdbqt",
        type=Path,
        default=PROJECT_ROOT / "data/structures/prepared/9DKB_receptor.pdbqt",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT
        / "data/campaigns/c5/02_urat1_ifp/urat1_key_residues.json",
    )
    args = ap.parse_args()
    mapped, unmatched = map_key_residues(args.src_cif, args.prep_pdbqt)
    out = {
        "pdb": "9DKB",
        "n_key_expected": len(URAT1_KEY),
        "n_key_mapped": len(mapped),
        "residues": mapped,
        "unmatched": unmatched,
        "note": (
            "Contacts use prepared-receptor atom coords; literature labels via "
            "Cα match to 9DKB.cif. Q437 is LEU in this structure and is left unmatched."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"n_mapped": len(mapped), "labels": list(mapped), "unmatched": unmatched}, indent=2))
    if len(mapped) < 10:
        raise SystemExit(f"Only mapped {len(mapped)} key residues (need >=10)")


if __name__ == "__main__":
    main()
