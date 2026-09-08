#!/usr/bin/env python3
"""Extract crystal reference ligands and Arg477 guanidinium coords for C1 L2."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def extract_hetatm_pdb(pdb_path: Path, resn: str, out_pdb: Path) -> int:
    lines = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("HETATM", "ATOM")) and len(line) >= 20:
            if line[17:20].strip() == resn:
                lines.append(line)
        elif line.startswith("TER") and lines:
            break
    if not lines:
        raise SystemExit(f"No residue {resn} in {pdb_path}")
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(lines) + "\nEND\n")
    return len(lines)


def pdb_het_to_sdf(pdb_path: Path, sdf_path: Path, template_smiles: str | None = None) -> None:
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=True, sanitize=False)
    if mol is None:
        raise SystemExit(f"RDKit failed to read {pdb_path}")
    if template_smiles:
        tmpl = Chem.MolFromSmiles(template_smiles)
        if tmpl is None:
            raise SystemExit(f"bad template SMILES: {template_smiles}")
        mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    w = Chem.SDWriter(str(sdf_path))
    w.write(mol)
    w.close()


def arg477_from_cif(cif_path: Path) -> dict:
    """Return guanidinium N coords for ARG auth_seq_id 477 (paper numbering)."""
    lines = cif_path.read_text().splitlines()
    start = None
    cols: list[str] = []
    for i, line in enumerate(lines):
        if line.startswith("loop_") and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            j = i + 1
            while j < len(lines) and lines[j].startswith("_atom_site."):
                cols.append(lines[j].strip())
                j += 1
            start = j
            break
    if start is None:
        raise SystemExit("atom_site loop not found")
    idx = {c: i for i, c in enumerate(cols)}
    atoms = {}
    for line in lines[start:]:
        if line.startswith(("#", "_", "loop_")) or not line.strip():
            break
        parts = line.split()
        if len(parts) < len(cols):
            continue
        if parts[idx["_atom_site.label_comp_id"]] != "ARG":
            continue
        auth = parts[idx["_atom_site.auth_seq_id"]]
        if auth != "477":
            continue
        atom = parts[idx["_atom_site.label_atom_id"]]
        if atom in {"NE", "CZ", "NH1", "NH2"}:
            atoms[atom] = [
                float(parts[idx["_atom_site.Cartn_x"]]),
                float(parts[idx["_atom_site.Cartn_y"]]),
                float(parts[idx["_atom_site.Cartn_z"]]),
            ]
    if "NH1" not in atoms or "NH2" not in atoms:
        raise SystemExit(f"Incomplete Arg477 sidechain in CIF: {atoms.keys()}")
    return {
        "residue_label": "Arg477",
        "source": str(cif_path),
        "note": "Prepared 9DKB PDBQT renumbers this residue as ARG A 476; same guanidinium.",
        "atoms": atoms,
    }


def copy_lesinurad_crystal_ref(out_dir: Path) -> Path:
    """Prefer existing redock crystal SDF if present."""
    candidates = [
        PROJECT_ROOT / "data/redock_smoke/lesinurad_9DKB/prep/lesinurad_crystal_ref_heavy.sdf",
        PROJECT_ROOT / "data/redock_smoke/lesinurad_9DKB/prep/lesinurad_crystal_ref.sdf",
        PROJECT_ROOT / "data/structures/pdb/lesinurad_crystal_ref_heavy.sdf",
    ]
    for c in candidates:
        if c.exists():
            dest = out_dir / "lesinurad_crystal_ref.sdf"
            dest.write_bytes(c.read_bytes())
            return dest
    raise SystemExit("No lesinurad crystal reference SDF found")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs",
    )
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    # NP3-146 / RM5 from 7ALV
    pdb7 = PROJECT_ROOT / "data/structures/pdb/7ALV.pdb"
    rm5_pdb = out / "RM5_crystal.pdb"
    n = extract_hetatm_pdb(pdb7, "RM5", rm5_pdb)
    rm5_sdf = out / "NP3-146_RM5_crystal_ref.sdf"
    # Dekker 2021 / PDB RM5 (NP3-146)
    rm5_smi = "CC(C)c1cc(Cl)cc(C(C)C)c1NC(=O)NS(=O)(=O)c1cc(C(C)(C)O)co1"
    pdb_het_to_sdf(rm5_pdb, rm5_sdf, template_smiles=rm5_smi)

    les_sdf = copy_lesinurad_crystal_ref(out)
    arg = arg477_from_cif(PROJECT_ROOT / "data/structures/pdb/9DKB.cif")
    (out / "arg477_coords.json").write_text(json.dumps(arg, indent=2))

    meta = {
        "rm5_atoms": n,
        "rm5_sdf": str(rm5_sdf),
        "lesinurad_sdf": str(les_sdf),
        "arg477": arg,
    }
    (out / "refs_summary.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
