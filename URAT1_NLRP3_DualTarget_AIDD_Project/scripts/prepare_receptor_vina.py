#!/usr/bin/env python3
"""Prepare protein receptor PDBQT for AutoDock Vina from PDB/mmCIF."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import gemmi
import yaml
from openbabel import openbabel as ob

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "docking_open_source.yaml"


def load_target_config(config_path: Path, target_key: str) -> dict:
    cfg = yaml.safe_load(config_path.read_text())
    if target_key not in cfg.get("targets", {}):
        raise KeyError(f"Unknown target {target_key}; keys={list(cfg.get('targets', {}))}")
    return cfg["targets"][target_key], cfg.get("receptor_prep", {})


def read_structure(path: Path) -> gemmi.Structure:
    path = Path(path)
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return gemmi.read_structure(str(path))
    return gemmi.read_pdb(str(path))


def structure_to_protein_pdb(
    st: gemmi.Structure,
    chain_id: str = "A",
    remove_waters: bool = True,
) -> str:
    """Return PDB text for protein only (no ligand/ions/waters)."""
    water = {"HOH", "WAT", "DOD"}
    lines: list[str] = []
    atom_serial = 1
    for model in st:
        for chain in model:
            if chain_id and chain.name != chain_id:
                continue
            for residue in chain:
                resn = residue.name.strip()
                if remove_waters and resn in water:
                    continue
                if residue.het_flag == "H":
                    continue
                for atom in residue:
                    lines.append(
                        f"ATOM  {atom_serial:5d} {atom.name:^4s}{resn:>3s} "
                        f"{chain.name}{residue.seqid.num:4d}    "
                        f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                        f"  1.00  0.00           {atom.element.name:>2s}\n"
                    )
                    atom_serial += 1
        break
    lines.append("END\n")
    if atom_serial == 1:
        raise ValueError(f"No protein atoms found for chain {chain_id}")
    return "".join(lines)


def pdb_to_pdbqt(pdb_text: str, ph: float = 7.4) -> str:
    obconv = ob.OBConversion()
    obconv.SetInAndOutFormats("pdb", "pdbqt")
    # Rigid receptor: no ROOT/BRANCH torsion tree (required by Vina)
    obconv.AddOption("r", ob.OBConversion.OUTOPTIONS)
    obconv.AddOption("p", ob.OBConversion.OUTOPTIONS, str(ph))
    obconv.AddOption("h", ob.OBConversion.OUTOPTIONS)  # preserve hydrogens added below
    mol = ob.OBMol()
    if not obconv.ReadString(mol, pdb_text):
        raise RuntimeError("Open Babel failed to read protein PDB")
    mol.AddHydrogens(False, True, ph)
    out = obconv.WriteString(mol)
    if not out:
        raise RuntimeError("Open Babel failed to write receptor PDBQT")
    # Safety: strip any ligand-style tags if present
    lines = [
        ln
        for ln in out.splitlines()
        if ln.startswith(("ATOM", "HETATM", "REMARK")) or not ln.strip()
    ]
    # Vina rigid receptor must not contain ROOT/BRANCH/TORSDOF
    return "\n".join(lines) + "\n"


def prepare_receptor(
    structure_file: Path,
    output_pdbqt: Path,
    chain_id: str = "A",
    remove_waters: bool = True,
    ph: float = 7.4,
) -> dict:
    st = read_structure(structure_file)
    pdb_text = structure_to_protein_pdb(st, chain_id=chain_id, remove_waters=remove_waters)
    pdbqt = pdb_to_pdbqt(pdb_text, ph=ph)
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    output_pdbqt.write_text(pdbqt)
    return {
        "structure_file": str(structure_file),
        "output": str(output_pdbqt),
        "chain_id": chain_id,
        "n_atom_lines": pdbqt.count("ATOM") + pdbqt.count("HETATM"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare receptor PDBQT for Vina")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--target", type=str, required=True, help="e.g. urat1_9dkb or nlrp3_7alv")
    parser.add_argument("--structure", type=Path, default=None, help="Override structure file")
    parser.add_argument("--output", type=Path, default=None, help="Override output PDBQT")
    args = parser.parse_args()

    target, prep = load_target_config(args.config, args.target)
    structure = args.structure or (PROJECT_ROOT / target["structure_file"])
    output = args.output or (PROJECT_ROOT / target["prepared_receptor"])
    summary = prepare_receptor(
        structure,
        output,
        chain_id=prep.get("chain", "A"),
        remove_waters=prep.get("remove_waters", True),
        ph=float(prep.get("protonate_ph", 7.4)),
    )
    summary["pdb_id"] = target["pdb_id"]
    summary["center"] = target["center"]
    summary["size"] = target["size"]
    qc = output.with_suffix(".qc.json")
    qc.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
