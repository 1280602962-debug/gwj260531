#!/usr/bin/env python3
"""Phase A: maegz -> SDF/PDBQT, boxes from X6K, receptor PDBQT."""
from __future__ import annotations

import csv
import gzip
import json
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0")
MAEGZ = Path(
    "/mnt/d/CADD paper exercise/dual target docking/Maestro doc/"
    "pik3ca_mtor_panel48_v0_ligprep/pik3ca_mtor_panel48_v0_ligprep-out.maegz"
)
MK_PREP = Path("/home/gwj/miniconda3/bin/mk_prepare_receptor.py")
PAD = 5.0
MIN_EDGE = 20.0


def load_maegz_mols():
    text = gzip.open(MAEGZ, "rt", errors="replace").read().replace("\\_", "_")
    tmp = Path("/tmp/pik3ca_mtor_panel48_v0_ligprep-out.fixed.mae")
    tmp.write_text(text)
    suppl = Chem.MaeMolSupplier(str(tmp), removeHs=False)
    by_id = {}
    for mol in suppl:
        if mol is None:
            continue
        pid = mol.GetProp("s_canvas_panel_id") if mol.HasProp("s_canvas_panel_id") else None
        if not pid:
            raise RuntimeError("missing s_canvas_panel_id")
        by_id.setdefault(pid, []).append(mol)
    return by_id


def write_ligands(by_id: dict) -> list[dict]:
    preparator = MoleculePreparation()
    rows = []
    for i in range(1, 49):
        pid = f"PM48_{i:02d}"
        if pid not in by_id:
            raise RuntimeError(f"missing {pid} in maegz")
        confs = by_id[pid]
        # one conf only: take first
        mol = confs[0]
        chembl = mol.GetProp("s_canvas_molecule_chembl_id") if mol.HasProp("s_canvas_molecule_chembl_id") else ""
        pref = mol.GetProp("s_canvas_pref_name") if mol.HasProp("s_canvas_pref_name") else ""
        sdf = ROOT / "ligands_sdf" / f"{pid}.sdf"
        pdbqt = ROOT / "ligands_pdbqt" / f"{pid}.pdbqt"
        w = Chem.SDWriter(str(sdf))
        w.write(mol)
        w.close()
        setups = preparator.prepare(mol)
        if not setups:
            raise RuntimeError(f"meeko failed for {pid}")
        pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(setups[0])
        if not success:
            raise RuntimeError(f"pdbqt write failed {pid}: {error_msg}")
        pdbqt.write_text(pdbqt_string)
        note = "first_conf" if len(confs) == 1 else f"n_confs={len(confs)}; chose index 0"
        rows.append(
            {
                "panel_id": pid,
                "chembl_id": chembl,
                "pref_name": pref,
                "n_confs_in_maegz": len(confs),
                "chosen_conf_index": 0,
                "sdf_path": str(sdf),
                "pdbqt_path": str(pdbqt),
                "notes": note,
            }
        )
    with (ROOT / "tables" / "ligand_input_manifest.csv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "panel_id",
                "chembl_id",
                "pref_name",
                "n_confs_in_maegz",
                "chosen_conf_index",
                "sdf_path",
                "pdbqt_path",
                "notes",
            ],
        )
        w.writeheader()
        w.writerows(rows)
    return rows


def extract_x6k_and_box(tag: str, prepared_name: str) -> dict:
    src = ROOT / "receptors" / prepared_name
    lines = []
    coords = []
    for line in src.read_text(errors="ignore").splitlines():
        if line.startswith(("HETATM", "ATOM  ")) and line[17:20].strip() == "X6K":
            lines.append(line)
            elem = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if elem.startswith("H"):
                continue
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not coords:
        raise RuntimeError(f"no X6K heavy atoms in {src}")
    out_pdb = ROOT / "tables" / f"{tag}_cocrystal_X6K.pdb"
    out_pdb.write_text("\n".join(lines) + "\n")
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    size = [
        max(max(xs) - min(xs) + 2 * PAD, MIN_EDGE),
        max(max(ys) - min(ys) + 2 * PAD, MIN_EDGE),
        max(max(zs) - min(zs) + 2 * PAD, MIN_EDGE),
    ]
    center = [(min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0, (min(zs) + max(zs)) / 2.0]
    box = {
        "center_x": round(center[0], 3),
        "center_y": round(center[1], 3),
        "center_z": round(center[2], 3),
        "size_x": round(size[0], 3),
        "size_y": round(size[1], 3),
        "size_z": round(size[2], 3),
        "n_ligand_atoms": len(lines),
        "n_heavy_atoms": len(coords),
    }
    (ROOT / "boxes" / f"{tag}_box.json").write_text(json.dumps(box, indent=2))
    (ROOT / "boxes" / f"{tag}_box.txt").write_text(
        "\n".join(f"{k} = {v}" for k, v in box.items() if k.startswith(("center", "size"))) + "\n"
    )
    return box


def prepare_receptors():
    for tag in ["4L23", "4JT6"]:
        protein = ROOT / "receptors" / f"{tag}_protein.pdb"
        out = ROOT / "receptors" / f"{tag}_receptor.pdbqt"
        cmd = [
            str(MK_PREP),
            "--read_pdb",
            str(protein),
            "-o",
            str(out.with_suffix("")),  # meeko adds .pdbqt
            "-p",
            str(out),
        ]
        # mk_prepare_receptor API varies; try common forms
        trials = [
            [str(MK_PREP), "--read_pdb", str(protein), "-o", str(out)],
            [str(MK_PREP), "--read_pdb", str(protein), "-o", str(out.with_suffix(""))],
        ]
        ok = False
        last = None
        for c in trials:
            r = subprocess.run(c, capture_output=True, text=True)
            last = r
            if r.returncode == 0 and out.exists():
                ok = True
                break
            # sometimes writes without extension handling
            alt = out.with_suffix(".pdbqt")
            if alt.exists():
                ok = True
                break
        if not ok:
            raise RuntimeError(f"mk_prepare_receptor failed for {tag}: {last.stdout}\n{last.stderr}")
        print("receptor", tag, out.exists(), out.stat().st_size if out.exists() else 0)


def write_protocol(boxes: dict):
    yaml = f"""freeze_id: pik3ca_mtor_panel48_v0
engine: AutoDock_Vina
vina_version: "1.2.7"
seed_policy: fixed_global
seed_fixed_global: 20260727
exhaustiveness: 8
n_modes: 9
energy_range: 3
cpu_per_job: 1
box_definition: "AABB(X6K) + 5A padding; min edge 20A"
cognate_qc_gate: "both ends best_of_9 heavy-atom RMSD < 2.0 A"
targets:
  PIK3CA:
    pdb: 4L23
    prepared_protein: "D:/CADD paper exercise/dual target docking/Maestro doc/4L23_PIK3CA_prepared.pdb"
    prepared_protein_in_freeze: receptors/4L23_PIK3CA_prepared.pdb
    receptor_pdbqt: receptors/4L23_receptor.pdbqt
    ligand_cocrystal_resname: X6K
    pose_gold_panel_id: PM48_01
    box: {json.dumps(boxes['4L23'])}
  MTOR:
    pdb: 4JT6
    prepared_protein: "D:/CADD paper exercise/dual target docking/Maestro doc/4JT6_mTOR_prepared.pdb"
    prepared_protein_in_freeze: receptors/4JT6_mTOR_prepared.pdb
    receptor_pdbqt: receptors/4JT6_receptor.pdbqt
    ligand_cocrystal_resname: X6K
    pose_gold_panel_id: PM48_01
    box: {json.dumps(boxes['4JT6'])}
ligand_prep:
  tool: Schrodinger_LigPrep
  maegz: "D:/CADD paper exercise/dual target docking/Maestro doc/pik3ca_mtor_panel48_v0_ligprep/pik3ca_mtor_panel48_v0_ligprep-out.maegz"
  conversion: "RDKit MaeMolSupplier (escaped underscores fixed) -> SDF; meeko -> PDBQT; one conf per panel_id"
"""
    (ROOT / "protocol" / "protocol.yaml").write_text(yaml)
    (ROOT / "protocol" / "SEED_POLICY.md").write_text(
        """# SEED POLICY — pik3ca_mtor_panel48_v0

- seed_policy: fixed_global
- seed_fixed_global: 20260727
- exhaustiveness: 8 (from EGFR/HER2 exhaustiveness_v0_1)
- n_modes: 9
- Do not raise exhaustiveness for score cosmetics.
"""
    )


def main():
    print("Loading maegz...")
    by_id = load_maegz_mols()
    print("panel ids in maegz", len(by_id))
    rows = write_ligands(by_id)
    print("wrote ligands", len(rows))
    boxes = {
        "4L23": extract_x6k_and_box("4L23", "4L23_PIK3CA_prepared.pdb"),
        "4JT6": extract_x6k_and_box("4JT6", "4JT6_mTOR_prepared.pdb"),
    }
    (ROOT / "boxes" / "all_boxes.json").write_text(json.dumps(boxes, indent=2))
    print("boxes", boxes)
    print("Preparing receptors...")
    prepare_receptors()
    write_protocol(boxes)
    print("Phase A done")


if __name__ == "__main__":
    main()
