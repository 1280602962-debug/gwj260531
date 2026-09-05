#!/usr/bin/env python3
"""
Export MD-ready receptor + ligand files for the molecules selected by
scripts/select_md_candidates.py. No MD simulation is run by this script —
it only prepares input files (protein PDB, ligand SDF/PDB/SMILES) so MD can
be executed elsewhere (cloud/HPC/another workstation).

For each selected (compound, target) pair with a successful docked pose:
  - receptor.pdb : protein-only, water/hetatm-removed, protonated at the same
                    pH used for docking (reuses scripts/prepare_receptor_vina.py
                    logic so the receptor matches the frame the pose was docked into).
  - ligand.sdf   : the docked pose exactly as produced by gnina (P2 protocol).
  - ligand.pdb   : the same pose converted to PDB for quick visualization/complex assembly.
  - ligand.smi   : the ground-truth canonical SMILES (for force-field parameterization).
  - README.md    : provenance (PDB ID, docking engine/score, protonation, next steps).

Example:
  python3 scripts/export_md_ready_candidates.py \\
    --selection data/md_candidates/md_candidate_selection.csv \\
    --output-dir data/md_candidates
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from prepare_receptor_vina import load_target_config, read_structure, structure_to_protein_pdb  # noqa: E402

DEFAULT_SELECTION = PROJECT_ROOT / "data" / "md_candidates" / "md_candidate_selection.csv"
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "docking_production_p2.yaml"
DEFAULT_9DKB_DOCK_DIR = PROJECT_ROOT / "results" / "repurposing" / "docking_p2" / "9dkb"
DEFAULT_7ALV_DOCK_DIR = PROJECT_ROOT / "results" / "repurposing" / "docking_p2" / "7alv"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "md_candidates"

TARGET_KEY_BY_PDB = {"9DKB": "urat1_9dkb", "7ALV": "nlrp3_7alv"}


def protein_pdb_with_h(structure_file: Path, chain_id: str, remove_waters: bool, ph: float) -> str:
    """Protein-only PDB (no waters/hetatm), hydrogens added at the docking pH.

    Mirrors scripts/prepare_receptor_vina.py's preparation but writes plain
    PDB instead of PDBQT, so it matches the exact receptor frame used for
    docking while being directly usable by external MD prep tools
    (pdb2gmx/tleap etc. will still assign their own force-field-consistent
    protonation/topology; this file establishes the coordinate frame and a
    reasonable starting protonation state).
    """
    from openbabel import openbabel as ob

    st = read_structure(structure_file)
    pdb_text = structure_to_protein_pdb(st, chain_id=chain_id, remove_waters=remove_waters)

    obconv = ob.OBConversion()
    obconv.SetInAndOutFormats("pdb", "pdb")
    mol = ob.OBMol()
    if not obconv.ReadString(mol, pdb_text):
        raise RuntimeError(f"Open Babel failed to read protein PDB from {structure_file}")
    mol.AddHydrogens(False, True, ph)
    out = obconv.WriteString(mol)
    if not out:
        raise RuntimeError("Open Babel failed to write protonated receptor PDB")
    return out


def get_or_build_receptor(pdb_id: str, cfg: dict, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{pdb_id}_receptor.pdb"
    if out_path.exists():
        return out_path
    target_key = TARGET_KEY_BY_PDB[pdb_id]
    target, prep = load_target_config(cfg, target_key)
    structure_file = PROJECT_ROOT / target["structure_file"]
    text = protein_pdb_with_h(
        structure_file,
        chain_id=prep.get("chain", "A"),
        remove_waters=prep.get("remove_waters", True),
        ph=float(prep.get("protonate_ph", 7.4)),
    )
    out_path.write_text(text)
    return out_path


def export_ligand(pose_sdf: Path, canonical_smiles: str, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    status = {"sdf": False, "pdb": False, "smi": False, "warning": None}

    sdf_out = out_dir / "ligand.sdf"
    sdf_out.write_text(pose_sdf.read_text())
    status["sdf"] = True

    smi_out = out_dir / "ligand.smi"
    smi_out.write_text(f"{canonical_smiles}\n")
    status["smi"] = True

    try:
        suppl = Chem.SDMolSupplier(str(pose_sdf), removeHs=False, sanitize=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            status["warning"] = "Could not parse pose SDF with RDKit; ligand.pdb not written."
        else:
            try:
                Chem.SanitizeMol(mol, catchErrors=True)
            except Exception:
                pass
            Chem.MolToPDBFile(mol, str(out_dir / "ligand.pdb"))
            status["pdb"] = True
    except Exception as exc:  # pragma: no cover - defensive
        status["warning"] = f"RDKit ligand.pdb export failed: {exc}"

    return status


def write_readme(
    out_dir: Path,
    *,
    compound_name: str,
    chembl_id: str,
    repurposing_id: str,
    pdb_id: str,
    md_category: str,
    dock_score: float | None,
    docking_engine: str,
    receptor_ph: float,
    ligand_status: dict,
) -> None:
    lines = [
        f"# MD-ready bundle: {compound_name} @ {pdb_id}",
        "",
        f"- repurposing_id: {repurposing_id}",
        f"- chembl_id: {chembl_id}",
        f"- category: {md_category}",
        f"- target PDB: {pdb_id}",
        f"- docking engine / protocol: {docking_engine}",
        f"- docking score (dock_score, lower=better convention): {dock_score}",
        f"- receptor protonation pH: {receptor_ph} (Open Babel AddHydrogens, same as used for docking)",
        "",
        "## Files",
        "- `receptor.pdb` (shared per PDB ID, see `../_receptors/`) — protein only, waters/hetatm removed",
        "- `ligand.sdf` — docked pose exactly as produced by gnina (P2)",
        f"- `ligand.pdb` — {'available' if ligand_status['pdb'] else 'NOT generated (see warning below)'}",
        "- `ligand.smi` — ground-truth canonical SMILES",
        "",
        "## Status",
    ]
    if ligand_status.get("warning"):
        lines.append(f"- WARNING: {ligand_status['warning']}")
    else:
        lines.append("- Ligand files exported without errors.")
    lines += [
        "",
        "## Suggested next steps for MD (to run on a machine with sufficient compute)",
        "1. Re-check protonation states for both receptor and ligand at physiological pH "
        "(this bundle's protonation is a docking-time approximation, not MD-validated).",
        "2. Parameterize the ligand (e.g. GAFF2/OpenFF via antechamber, or CGenFF) from `ligand.sdf`.",
        "3. Build the complex (align `ligand.pdb` pose into `receptor.pdb`), solvate, add ions, "
        "minimize, equilibrate, then run production MD (target 50-100 ns).",
        "4. Report RMSD/RMSF, key residue contacts, and (optionally) MM-GBSA as RELATIVE "
        "comparisons only — do not report as absolute binding free energies.",
    ]
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export MD-ready receptor/ligand files (no MD is run)")
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--dock-dir-9dkb", type=Path, default=DEFAULT_9DKB_DOCK_DIR)
    parser.add_argument("--dock-dir-7alv", type=Path, default=DEFAULT_7ALV_DOCK_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if not args.selection.exists():
        raise FileNotFoundError(
            f"{args.selection} not found. Run scripts/select_md_candidates.py first."
        )

    sel = pd.read_csv(args.selection)
    receptor_cache = args.output_dir / "_receptors"
    manifest_rows = []

    dock_dirs = {"9DKB": args.dock_dir_9dkb, "7ALV": args.dock_dir_7alv}
    docking_engine = "gnina_cpu (Pi*=P2, CNNaffinity)"

    for _, row in sel.iterrows():
        rid = row["repurposing_id"]
        for pdb_id, has_pose_col in (("9DKB", "has_9dkb_pose"), ("7ALV", "has_7alv_pose")):
            if not bool(row.get(has_pose_col, False)):
                continue

            pose_sdf = dock_dirs[pdb_id] / "poses" / f"{rid}_out.sdf"
            if not pose_sdf.exists():
                print(f"WARNING: pose file missing for {rid} @ {pdb_id}: {pose_sdf} — skipping.")
                continue

            receptor_ph = 7.4
            try:
                target_key = TARGET_KEY_BY_PDB[pdb_id]
                _, prep = load_target_config(args.config, target_key)
                receptor_ph = float(prep.get("protonate_ph", 7.4))
                receptor_path = get_or_build_receptor(pdb_id, args.config, receptor_cache)
                receptor_ok = True
            except Exception as exc:
                print(f"WARNING: receptor prep failed for {pdb_id}: {exc}")
                receptor_path = None
                receptor_ok = False

            folder_name = f"{pdb_id}_{rid}"
            out_dir = args.output_dir / folder_name
            ligand_status = export_ligand(pose_sdf, row["canonical_smiles"], out_dir)

            if receptor_ok:
                # Symlink-by-copy so each bundle is self-contained (avoids relative-path fragility).
                (out_dir / "receptor.pdb").write_text(receptor_path.read_text())

            # Try to recover the docking score from the docking CSV for provenance.
            dock_csv = dock_dirs[pdb_id] / f"docking_{pdb_id.lower()}_gnina.csv"
            dock_score = None
            if dock_csv.exists():
                try:
                    dscores = pd.read_csv(dock_csv, low_memory=False)
                    match = dscores[dscores["repurposing_id"] == rid]
                    if not match.empty and "dock_score" in match.columns:
                        dock_score = float(match.iloc[0]["dock_score"])
                except Exception:
                    pass

            write_readme(
                out_dir,
                compound_name=row["name"],
                chembl_id=row["chembl_id"],
                repurposing_id=rid,
                pdb_id=pdb_id,
                md_category=row["md_category"],
                dock_score=dock_score,
                docking_engine=docking_engine,
                receptor_ph=receptor_ph,
                ligand_status=ligand_status,
            )

            manifest_rows.append(
                {
                    "folder": folder_name,
                    "repurposing_id": rid,
                    "chembl_id": row["chembl_id"],
                    "name": row["name"],
                    "md_category": row["md_category"],
                    "pdb_id": pdb_id,
                    "dock_score": dock_score,
                    "receptor_ready": receptor_ok,
                    "ligand_sdf_ready": ligand_status["sdf"],
                    "ligand_pdb_ready": ligand_status["pdb"],
                    "warning": ligand_status.get("warning"),
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(manifest_rows)
    manifest_path = args.output_dir / "md_ready_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    print(f"\nExported {len(manifest_rows)} MD-ready (compound, target) bundles under {args.output_dir}/")
    print(f"Manifest: {manifest_path}")
    if len(manifest_rows):
        n_warn = manifest["warning"].notna().sum()
        if n_warn:
            print(f"WARNING: {n_warn} bundle(s) have warnings — check md_ready_manifest.csv before handoff.")


if __name__ == "__main__":
    main()
