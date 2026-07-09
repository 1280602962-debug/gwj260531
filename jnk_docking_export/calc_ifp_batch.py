#!/usr/bin/env python3
"""
Batch SIFt (Structural Interaction Fingerprint) for prepped complexes.

Requires Schrödinger Python:
  %SCHRODINGER%\\run.exe python3 calc_ifp_batch.py --config jobs_step1.json

Outputs (default ifp_results/):
  - ifp_all.csv
  - ifp_interactions.csv
  - ifp_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.structutils.interactionfp import StructuralInteractionFingerprintGenerator

KINASE_MAP = {
    "3ELJ": "JNK1",
    "4L7F": "JNK1",
    "3E7O": "JNK2",
    "3TTI": "JNK3",
    "4WHZ": "JNK3",
}


def sanitize(name: str) -> str:
    return re.sub(r"\s+", "_", name.strip())


def ligand_matches(name: str, include: list[str], exclude: list[str]) -> bool:
    upper = name.upper()
    for ex in exclude:
        if ex.upper().replace("-", "_") in upper.replace("-", "_"):
            return False
    if not include:
        return True
    norm = upper.replace("-", "_").replace(" ", "_")
    for inc in include:
        token = inc.upper().replace("-", "_").replace(" ", "_")
        if token in norm or norm in token:
            return True
    return False


def split_receptor_ligand(st: structure.Structure) -> tuple[structure.Structure, structure.Structure]:
    ligand_atoms = analyze.evaluate_asl(st, "ligand")
    if not ligand_atoms:
        raise ValueError("ASL 'ligand' matched no atoms; check complex structure")
    protein_atoms = analyze.evaluate_asl(st, "protein")
    if not protein_atoms:
        protein_atoms = [a.index for a in st.atom if a.index not in ligand_atoms]

    ligand_st = st.extract(ligand_atoms, copy_props=True)
    receptor_st = st.extract(protein_atoms, copy_props=True)
    return receptor_st, ligand_st


def parse_name(stem: str) -> tuple[str, str]:
    base = stem.replace("_prepped", "")
    parts = base.split("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return base, base


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch SIFt IFP on prepped complexes.")
    parser.add_argument("--config", default="jobs_step1.json")
    parser.add_argument("--in-dir", default="complexes_prepped")
    parser.add_argument("--out-dir", default="ifp_results")
    args = parser.parse_args()

    cfg = load_config(Path(args.config).resolve())
    options = cfg.get("options", {})
    ifp_cfg = cfg.get("ifp", {})
    include = options.get("ligand_include", [])
    exclude = options.get("ligand_exclude", ["JNK-IN-8"])

    root = Path(cfg.get("root", ".")).resolve()
    in_root = Path(args.in_dir)
    if not in_root.is_absolute():
        in_root = root / in_root
    out_root = Path(args.out_dir)
    if not out_root.is_absolute():
        out_root = root / out_root
    out_root.mkdir(parents=True, exist_ok=True)

    files = sorted(in_root.glob("**/*_prepped.maegz"))
    if not files:
        files = sorted(in_root.glob("**/*.maegz"))
    if not files:
        print(f"ERROR: no complexes under {in_root}", file=sys.stderr)
        return 1

    gen = StructuralInteractionFingerprintGenerator()
    summary_rows: list[dict] = []
    errors: list[str] = []
    count = 0

    print(f"Input : {in_root}")
    print(f"Output: {out_root}")
    print()

    for path in files:
        stem = path.stem.replace("_prepped", "")
        pdb_id, lig_name = parse_name(stem)
        if not ligand_matches(lig_name, include, exclude):
            print(f"SKIP {path.name}")
            continue

        kinase = KINASE_MAP.get(pdb_id.upper(), "")
        fp_id = f"{pdb_id}_{sanitize(lig_name)}"
        print(f"IFP {path.name} ({kinase})")

        try:
            st = structure.StructureReader.read(str(path))
            receptor_st, ligand_st = split_receptor_ligand(st)
            gen.setReceptorStructure(receptor_st)
            gen.generateFingerprint(
                ligand_st,
                fp_id,
                receptor_region=None,
                ligand_title=lig_name,
                nonpolar_hydrogens=bool(ifp_cfg.get("nonpolar_hydrogens", False)),
                receptor_st=receptor_st,
            )
            count += 1
            summary_rows.append(
                {
                    "pdb_id": pdb_id,
                    "kinase": kinase,
                    "ligand": lig_name,
                    "ifp_id": fp_id,
                    "source_file": str(path.relative_to(in_root)),
                    "status": "ok",
                }
            )
        except Exception as exc:  # noqa: BLE001
            msg = f"{path.name}: {exc}"
            print(f"  ERROR: {exc}")
            errors.append(msg)
            summary_rows.append(
                {
                    "pdb_id": pdb_id,
                    "kinase": kinase,
                    "ligand": lig_name,
                    "ifp_id": fp_id,
                    "source_file": str(path.relative_to(in_root)),
                    "status": f"error: {exc}",
                }
            )

    if count == 0:
        print("ERROR: no IFPs generated", file=sys.stderr)
        return 1

    ifp_all = out_root / "ifp_all.csv"
    ifp_int = out_root / "ifp_interactions.csv"
    ifp_sum = out_root / "ifp_summary.tsv"

    interacting_only = bool(ifp_cfg.get("interacting_only", True))
    gen.writeCSVFile(str(ifp_all), all_props=True, interacting_only=interacting_only)
    gen.writeInteractionsFile(str(ifp_int), all_props=True)

    with ifp_sum.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "kinase", "ligand", "ifp_id", "source_file", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print()
    print(f"Generated IFP for {count} complexes")
    print(f"  {ifp_all}")
    print(f"  {ifp_int}")
    print(f"  {ifp_sum}")
    if errors:
        err_log = out_root / "ifp_errors.log"
        err_log.write_text("\n".join(errors) + "\n", encoding="utf-8")
        print(f"  {err_log} ({len(errors)} errors)")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
