#!/usr/bin/env python3
"""P0 cognate-redocking QC for alternative PIK3CA/mTOR crystal structures.

This script deliberately performs a structure-level gate only.  It does not
redock the frozen PM panel.  A receptor PDBQT and its cognate-derived box are
promoted into ``receptors/`` only after best-of-nine heavy-atom RMSD < 2 Å.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_structure_robust_v0"
PYTHON = Path("/home/gwj/miniconda3/bin/python3")
VINA = Path("/home/gwj/miniconda3/bin/vina")
MEEKO_REC = Path("/home/gwj/miniconda3/bin/mk_prepare_receptor.py")
MEEKO_LIG = Path("/home/gwj/miniconda3/bin/mk_prepare_ligand.py")
OBABEL = Path("/home/gwj/miniconda3/envs/cadd_tools/bin/obabel")

SEED = 20260727
EXHAUSTIVENESS = 16
N_POSES = 9
PAD = 5.0
MIN_EDGE = 20.0

# Entries were selected in docs/JCIM_SUPPLEMENTARY_EXPERIMENTS_PLAN_V2.md §3.
SPECS = (
    {"pdb_id": "4JPS", "target": "PIK3CA"},
    {"pdb_id": "5DXT", "target": "PIK3CA"},
    {"pdb_id": "4JSX", "target": "mTOR"},
)

# Solvents, elemental ions, common precipitants, and crystallographic buffers.
# Any remaining HET residue is recorded in the report before the largest organic
# component is selected as the cognate-ligand candidate.
EXCLUDED_HET = {
    "HOH", "DOD", "WAT", "H2O", "CL", "BR", "IOD", "F", "NA", "K", "CA",
    "MG", "MN", "ZN", "CO", "CU", "CD", "NI", "FE", "HG", "SO4", "PO4",
    "GOL", "EDO", "PEG", "PG4", "PGE", "MPD", "DMS", "ACT", "ACE", "FMT",
    "MES", "TRS", "BME", "CIT", "TLA", "NAG", "MAN", "BGC",
}


def fetch(url: str, path: Path) -> None:
    """Download once, with an explicit User-Agent accepted by RCSB."""
    if path.exists() and path.stat().st_size:
        return
    request = urllib.request.Request(url, headers={"User-Agent": "jcim-structure-qc/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        path.write_bytes(response.read())


def api_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "jcim-structure-qc/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def field(line: str, start: int, end: int) -> str:
    return line[start:end].strip()


def atom_element(line: str) -> str:
    element = field(line, 76, 78)
    if element:
        return element.upper()
    name = field(line, 12, 16).upper()
    return name.lstrip("0123456789")[:1]


def hetero_residues(pdb: Path) -> dict[tuple[str, str, str, str], list[str]]:
    residues: dict[tuple[str, str, str, str], list[str]] = {}
    for line in pdb.read_text(errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        key = (field(line, 17, 20), field(line, 21, 22), field(line, 22, 26), field(line, 26, 27))
        residues.setdefault(key, []).append(line)
    return residues


def choose_ligand(pdb: Path) -> tuple[tuple[str, str, str, str], list[str], list[dict]]:
    """Choose the largest non-excluded, organic HET residue and retain an audit."""
    candidates = []
    for key, lines in hetero_residues(pdb).items():
        resname, chain, resid, icode = key
        elements = [atom_element(line) for line in lines]
        heavy = sum(element != "H" for element in elements)
        organic = any(element == "C" for element in elements)
        candidates.append(
            {
                "resname": resname,
                "chain": chain or "_",
                "resid": resid,
                "icode": icode or "_",
                "heavy_atoms": heavy,
                "organic": organic,
                "excluded": resname in EXCLUDED_HET,
                "_lines": lines,
            }
        )
    eligible = [
        item for item in candidates
        if not item["excluded"] and item["organic"] and 6 <= item["heavy_atoms"] <= 100
    ]
    if not eligible:
        raise RuntimeError(f"No non-solvent organic HET ligand found in {pdb.name}")
    selected = max(eligible, key=lambda item: item["heavy_atoms"])
    audit = [{k: v for k, v in item.items() if k != "_lines"} for item in candidates]
    return (
        (selected["resname"], selected["chain"], selected["resid"], selected["icode"]),
        selected["_lines"],
        audit,
    )


def write_pdb(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\nEND\n")


def write_protein(pdb: Path, out: Path) -> None:
    # Match the frozen receptor preparation style: standard polymer ATOM records
    # only, retaining TER boundaries and avoiding crystallographic additives.
    lines = [line for line in pdb.read_text(errors="replace").splitlines() if line.startswith(("ATOM  ", "TER"))]
    if not any(line.startswith("ATOM  ") for line in lines):
        raise RuntimeError(f"No ATOM records in {pdb.name}")
    write_pdb(out, lines)


def heavy_xyz(lines: list[str]) -> np.ndarray:
    xyz = []
    for line in lines:
        if atom_element(line) == "H":
            continue
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not xyz:
        raise RuntimeError("Cognate ligand has no heavy atoms")
    return np.asarray(xyz, dtype=float)


def cognate_box(xyz: np.ndarray) -> dict:
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    center = (lo + hi) / 2.0
    size = np.maximum(hi - lo + 2 * PAD, MIN_EDGE)
    return {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": round(float(size[0]), 3),
        "size_y": round(float(size[1]), 3),
        "size_z": round(float(size[2]), 3),
        "n_heavy_atoms": int(len(xyz)),
        "definition": "cognate heavy-atom AABB center; AABB extent + 5 Å each side; 20 Å minimum edge",
    }


def checked_run(command: list[str], label: str) -> subprocess.CompletedProcess:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        detail = (result.stderr or result.stdout)[-1200:]
        raise RuntimeError(f"{label} failed (rc={result.returncode}):\n{detail}")
    return result


def prepare_receptor(protein: Path, receptor: Path) -> None:
    base = receptor.with_suffix("")
    checked_run(
        [
            str(PYTHON), str(MEEKO_REC), "--read_pdb", str(protein), "-o", str(base),
            "-p", "-a", "--default_altloc", "A",
        ],
        "Meeko receptor preparation",
    )
    produced = base.with_suffix(".pdbqt")
    if not produced.exists():
        matches = list(receptor.parent.glob(base.name + "*.pdbqt"))
        if not matches:
            raise RuntimeError("Meeko receptor preparation completed without a PDBQT")
        produced = matches[0]
    if produced != receptor:
        shutil.copyfile(produced, receptor)


def prepare_ligand(ligand_pdb: Path, ligand_sdf: Path, ligand_pdbqt: Path) -> None:
    # Open Babel is only used to infer the PDB residue bond graph into SDF;
    # PDBQT preparation itself is exclusively Meeko 0.7.1.
    checked_run([str(OBABEL), "-ipdb", str(ligand_pdb), "-osdf", "-O", str(ligand_sdf), "-h"], "Open Babel PDB→SDF")
    checked_run(
        [str(PYTHON), str(MEEKO_LIG), "-i", str(ligand_sdf), "-o", str(ligand_pdbqt)],
        "Meeko ligand preparation",
    )
    if not ligand_pdbqt.exists():
        raise RuntimeError("Meeko ligand preparation completed without a PDBQT")


def vina_redock(receptor: Path, ligand: Path, box: dict, out: Path, log: Path) -> None:
    config = out.with_suffix(".vina.txt")
    config.write_text(
        "\n".join(
            [
                f"receptor = {receptor}", f"ligand = {ligand}", f"out = {out}",
                f"center_x = {box['center_x']}", f"center_y = {box['center_y']}", f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}", f"size_y = {box['size_y']}", f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUSTIVENESS}", f"num_modes = {N_POSES}",
                "energy_range = 3", "cpu = 1", f"seed = {SEED}",
            ]
        ) + "\n"
    )
    result = subprocess.run([str(VINA), "--config", str(config)], capture_output=True, text=True)
    log.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode or not out.exists():
        raise RuntimeError(f"Vina failed (rc={result.returncode}):\n{(result.stderr or result.stdout)[-1200:]}")
    if f"random seed: {SEED}" not in (result.stdout + result.stderr):
        raise RuntimeError("Vina log does not confirm the requested seed")


def meeko_smiles_and_map(pdbqt: Path) -> tuple[str, list[tuple[int, int]]]:
    smiles = None
    mapped: list[int] = []
    for line in pdbqt.read_text(errors="replace").splitlines():
        if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
            smiles = line.removeprefix("REMARK SMILES ").strip()
        elif line.startswith("REMARK SMILES IDX"):
            mapped.extend(int(value) for value in line.split()[3:])
    if not smiles or len(mapped) % 2:
        raise RuntimeError(f"Could not recover Meeko SMILES atom map from {pdbqt.name}")
    return smiles, list(zip(mapped[0::2], mapped[1::2]))


def pose_molecule(pdbqt: Path, template: Chem.Mol, mapping: list[tuple[int, int]]) -> Chem.Mol:
    atoms = {}
    for line in pdbqt.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atoms[int(line[6:11])] = (atom_element(line), (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    molecule = Chem.Mol(template)
    conf = Chem.Conformer(molecule.GetNumAtoms())
    n_heavy = 0
    for smile_index, serial in mapping:
        element, xyz = atoms.get(serial, (None, None))
        if element in {"H", "HD", "HS"}:
            continue
        if xyz is None:
            raise RuntimeError(f"Mapped atom serial {serial} absent from {pdbqt.name}")
        conf.SetAtomPosition(smile_index - 1, xyz)
        n_heavy += 1
    if n_heavy != molecule.GetNumAtoms():
        raise RuntimeError(f"Mapped {n_heavy}/{molecule.GetNumAtoms()} heavy atoms in {pdbqt.name}")
    molecule.RemoveAllConformers()
    molecule.AddConformer(conf, assignId=True)
    return molecule


def reference_molecules(sdf: Path, template: Chem.Mol) -> list[Chem.Mol]:
    reference = Chem.SDMolSupplier(str(sdf), removeHs=True, sanitize=True)[0]
    if reference is None:
        raise RuntimeError(f"RDKit could not read {sdf.name}")
    matches = reference.GetSubstructMatches(template, uniquify=False)
    if not matches:
        raise RuntimeError("Crystal SDF does not match the Meeko chemical template")
    source = reference.GetConformer()
    refs = []
    for match in matches:
        mol = Chem.Mol(template)
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for template_index, reference_index in enumerate(match):
            point = source.GetAtomPosition(reference_index)
            conformer.SetAtomPosition(template_index, point)
        mol.RemoveAllConformers()
        mol.AddConformer(conformer, assignId=True)
        refs.append(mol)
    return refs


def rmsd_no_superposition(reference: Chem.Mol, pose: Chem.Mol, automorphisms: tuple[tuple[int, ...], ...]) -> float:
    return min(
        float(rdMolAlign.CalcRMS(pose, reference, map=[list(enumerate(permutation))]))
        for permutation in automorphisms
    )


def entity_descriptions(pdb_id: str) -> list[str]:
    entry = api_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    descriptions = []
    for entity_id in entity_ids:
        entity = api_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
        description = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        descriptions.append(f"entity {entity_id}: {description}")
    return descriptions


def process(spec: dict, force: bool) -> dict:
    pdb_id, target = spec["pdb_id"], spec["target"]
    work = OUT / "cognate_qc" / pdb_id
    work.mkdir(parents=True, exist_ok=True)
    pdb = work / f"{pdb_id}.pdb"
    if force and pdb.exists():
        pdb.unlink()
    fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb", pdb)
    entities = entity_descriptions(pdb_id)

    ligand_key, ligand_lines, ligand_audit = choose_ligand(pdb)
    ligand_code, chain, resid, icode = ligand_key
    protein = work / f"{pdb_id}_protein.pdb"
    ligand_pdb = work / f"{pdb_id}_{ligand_code}_crystal.pdb"
    ligand_sdf = work / f"{pdb_id}_{ligand_code}_crystal.sdf"
    receptor = work / f"{pdb_id}_receptor.pdbqt"
    ligand_pdbqt = work / f"{pdb_id}_{ligand_code}.pdbqt"
    docked = work / f"{pdb_id}_cognate_all_modes.pdbqt"
    log = work / f"{pdb_id}_cognate_vina.log"

    write_protein(pdb, protein)
    write_pdb(ligand_pdb, ligand_lines)
    box = cognate_box(heavy_xyz(ligand_lines))
    box.update({"pdb_id": pdb_id, "target": target, "ligand_resname": ligand_code, "chain": chain or "_", "resid": resid, "icode": icode or "_"})
    (work / f"{pdb_id}_box.json").write_text(json.dumps(box, indent=2) + "\n")
    prepare_receptor(protein, receptor)
    prepare_ligand(ligand_pdb, ligand_sdf, ligand_pdbqt)
    vina_redock(receptor, ligand_pdbqt, box, docked, log)

    smiles, atom_map = meeko_smiles_and_map(ligand_pdbqt)
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise RuntimeError("RDKit could not parse Meeko's ligand SMILES")
    references = reference_molecules(ligand_sdf, template)
    automorphisms = template.GetSubstructMatches(template, uniquify=False)
    mode_rmsds = []
    text = docked.read_text(errors="replace")
    models = []
    current = []
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            current = [line]
        elif line.startswith("ENDMDL"):
            current.append(line)
            models.append(current)
            current = []
        elif current:
            current.append(line)
    if len(models) != N_POSES:
        raise RuntimeError(f"Expected {N_POSES} docking poses, found {len(models)}")
    for index, model in enumerate(models, 1):
        pose_file = work / f"mode_{index:02d}.pdbqt"
        pose_file.write_text("".join(model))
        pose = pose_molecule(pose_file, template, atom_map)
        rmsd = min(rmsd_no_superposition(reference, pose, automorphisms) for reference in references)
        mode_rmsds.append(round(rmsd, 3))
    best_mode, best = min(enumerate(mode_rmsds, 1), key=lambda value: value[1])
    passed = best < 2.0

    if passed:
        receptors = OUT / "receptors"
        shutil.copyfile(protein, receptors / f"{pdb_id}_protein.pdb")
        shutil.copyfile(receptor, receptors / f"{pdb_id}_receptor.pdbqt")
        (receptors / f"{pdb_id}_box.json").write_text(json.dumps(box, indent=2) + "\n")
    return {
        "pdb_id": pdb_id, "target": target, "polymer_entity_descriptions": entities,
        "selected_ligand": {"resname": ligand_code, "chain": chain or "_", "resid": resid, "icode": icode or "_"},
        "hetero_residue_audit": ligand_audit, "box": box,
        "vina_version": "1.2.7", "meeko_version": "0.7.1", "seed": SEED,
        "exhaustiveness": EXHAUSTIVENESS, "n_poses": N_POSES,
        "rmsd_definition": "heavy atom; Meeko SMILES-atom mapping; template automorphism minimum; no protein or ligand superposition",
        "rmsd_mode1_angstrom": mode_rmsds[0], "rmsd_best_of_9_angstrom": best,
        "best_mode": best_mode, "pass_best_of_9_lt_2A": passed,
        "status": "PASS" if passed else "FAIL_RMSD",
    }


def write_report(rows: list[dict]) -> None:
    report = OUT / "analysis" / "STRUCTURE_ROBUSTNESS_QC_V1.md"
    lines = [
        "# Structure robustness cognate QC v1",
        "",
        "P0 alternative-receptor gate for the PIK3CA/mTOR PM study. This is cognate redocking QC only; no frozen-panel redocking was started.",
        "",
        "## Frozen protocol",
        "",
        f"- AutoDock Vina 1.2.7; Meeko 0.7.1; seed `{SEED}`; exhaustiveness `{EXHAUSTIVENESS}`; `n_poses={N_POSES}`; `cpu=1`.",
        f"- Box: cognate-ligand heavy-atom AABB center, `{PAD:.0f}` Å padding on each side, `{MIN_EDGE:.0f}` Å minimum edge (the same construction used for the frozen PM boxes).",
        "- RMSD: heavy atoms only, Meeko SMILES-index mapping with all template automorphisms, no superposition. Pass gate: best-of-9 < 2.0 Å.",
        "- Receptor preparation: Meeko `mk_prepare_receptor.py --read_pdb -p -a --default_altloc A`; ligand PDBQT: Meeko. Open Babel only converts the crystal PDB residue to SDF to retain its bond graph.",
        "",
        "## Results",
        "",
        "| target | PDB | cognate ligand | RMSD mode 1 (Å) | best-of-9 RMSD (Å) | best mode | verdict |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for row in rows:
        ligand = row.get("selected_ligand", {}).get("resname", "")
        lines.append(
            f"| {row.get('target', '')} | {row.get('pdb_id', '')} | {ligand} | "
            f"{row.get('rmsd_mode1_angstrom', '')} | {row.get('rmsd_best_of_9_angstrom', '')} | "
            f"{row.get('best_mode', '')} | {row.get('status', 'ERROR')} |"
        )
    lines += ["", "## Polymer-entity identity audit", ""]
    for row in rows:
        lines.append(f"### {row.get('pdb_id')} — {row.get('target')}")
        for description in row.get("polymer_entity_descriptions", []):
            lines.append(f"- {description}")
        if row.get("status") == "ERROR":
            lines.append(f"- QC error: `{row.get('error')}`")
        else:
            selected = row["selected_ligand"]
            lines.append(
                f"- Selected cognate: `{selected['resname']}` chain `{selected['chain']}`, residue `{selected['resid']}{selected['icode']}`. "
                "Selection was the largest non-solvent organic HET residue; full candidate audit is in the JSON record."
            )
            if row["pass_best_of_9_lt_2A"]:
                lines.append(f"- PASS: promoted receptor PDBQT and box to `data/jcim_structure_robust_v0/receptors/{row['pdb_id']}_*`.")
            else:
                lines.append("- FAIL: retained QC artifacts only; this receptor is not eligible for panel redocking.")
        lines.append("")
    lines += [
        "## Interpretation",
        "",
        "A PASS establishes that the frozen docking preparation and cognate-centered box can recover the deposited ligand for that alternative crystal structure. It does not itself establish panel-level structural robustness; only PASS receptors may be used in the later, separate frozen-panel redocking step.",
        "",
    ]
    report.write_text("\n".join(lines))
    (OUT / "analysis" / "structure_robustness_qc_v1.json").write_text(json.dumps(rows, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", choices=[spec["pdb_id"] for spec in SPECS], action="append", help="Run only one candidate; repeatable.")
    parser.add_argument("--force-download", action="store_true", help="Replace cached RCSB PDB downloads.")
    args = parser.parse_args()
    selected = [spec for spec in SPECS if not args.pdb or spec["pdb_id"] in args.pdb]
    rows = []
    for spec in selected:
        print(f"QC {spec['target']} {spec['pdb_id']}", flush=True)
        try:
            row = process(spec, args.force_download)
        except Exception as exc:  # retain completed candidates and make the blocker auditable
            row = {"pdb_id": spec["pdb_id"], "target": spec["target"], "status": "ERROR", "error": str(exc)}
        rows.append(row)
        print(json.dumps(row, indent=2), flush=True)
    write_report(rows)
    return 1 if any(row["status"] == "ERROR" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
