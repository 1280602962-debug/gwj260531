#!/usr/bin/env python3
"""Recompute C5 W1 docking-box centers from deposited coordinates.

9DKB and 7ALV keep the locked production boxes in config/docking_c1.yaml.
9DKA / 9DKC use the crystal-ligand heavy-atom centroid in the native frame.
9DK9 (apo) receives the 9DKB locked center after chain-A CA Kabsch superposition.

This script does not invent a new 9DKB box.
"""
from __future__ import annotations

import json
from pathlib import Path

import gemmi
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = PROJECT_ROOT / "data" / "structures" / "pdb"
OUT = PROJECT_ROOT / "data" / "campaigns" / "c5" / "00_verification" / "w1_box_centers.json"

LOCKED_9DKB = np.array([99.966, 102.967, 105.699], dtype=float)
LOCKED_7ALV = np.array([16.756, 35.449, 125.714], dtype=float)


def load_cif(pdb_id: str) -> gemmi.Structure:
    return gemmi.read_structure(str(PDB_DIR / f"{pdb_id}.cif"))


def ligand_heavy_atoms(st: gemmi.Structure, resname: str) -> list[dict]:
    atoms: list[dict] = []
    for chain in st[0]:
        for res in chain:
            if res.name.strip() != resname:
                continue
            for atom in res:
                if atom.is_hydrogen():
                    continue
                atoms.append(
                    {
                        "chain": chain.name,
                        "seqid": res.seqid.num,
                        "atom": atom.name.strip(),
                        "xyz": [atom.pos.x, atom.pos.y, atom.pos.z],
                    }
                )
    if not atoms:
        raise ValueError(f"no heavy atoms for residue {resname}")
    return atoms


def centroid(atoms: list[dict]) -> np.ndarray:
    return np.array([a["xyz"] for a in atoms], dtype=float).mean(axis=0)


def ca_map(st: gemmi.Structure, chain_id: str = "A") -> dict[tuple[str, int], np.ndarray]:
    out: dict[tuple[str, int], np.ndarray] = {}
    for chain in st[0]:
        if chain.name != chain_id:
            continue
        for res in chain:
            if res.het_flag == "H":
                continue
            ca = res.get_ca()
            if ca is None:
                continue
            out[(res.name.strip(), int(res.seqid.num))] = np.array(
                [ca.pos.x, ca.pos.y, ca.pos.z], dtype=float
            )
    return out


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    c_m = mobile.mean(axis=0)
    c_t = target.mean(axis=0)
    x = mobile - c_m
    y = target - c_t
    u, _, vt = np.linalg.svd(x.T @ y)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    t = c_t - r @ c_m
    fitted = (r @ mobile.T).T + t
    rmsd = float(np.sqrt(((fitted - target) ** 2).sum(axis=1).mean()))
    return r, t, rmsd


def main() -> None:
    report: dict = {
        "purpose": "W1 docking-box centers from crystal ligand centroids / CA superposition",
        "method": "gemmi parse + Kabsch CA superposition (chain A)",
        "locked_production_boxes": {
            "9DKB": {"center": LOCKED_9DKB.tolist(), "size": [22, 22, 22], "source": "config/docking_c1.yaml"},
            "7ALV": {"center": LOCKED_7ALV.tolist(), "size": [20, 20, 20], "source": "config/docking_c1.yaml"},
        },
        "receptors": {},
        "structure_files": {
            "9DK9": "data/structures/pdb/9DK9.cif",
            "9DKA": "data/structures/pdb/9DKA.cif",
            "9DKB": "data/structures/pdb/9DKB.cif",
            "9DKC": "data/structures/pdb/9DKC.cif",
            "7ALV": "data/structures/pdb/7ALV.pdb",
            "9DKC_format": "mmCIF only; legacy .pdb download is HTTP 404",
        },
    }

    for pdb_id, lig in (("9DKB", "A1AIL"), ("9DKA", "R75"), ("9DKC", "A1A45")):
        atoms = ligand_heavy_atoms(load_cif(pdb_id), lig)
        cen = centroid(atoms)
        report["receptors"][pdb_id] = {
            "ligand_ccd": lig,
            "ligand_chains": sorted({a["chain"] for a in atoms}),
            "n_heavy_atoms": len(atoms),
            "crystal_ligand_centroid": [round(float(x), 3) for x in cen],
            "size": [22, 22, 22],
            "box_rule": "crystal_ligand_heavy_atom_centroid",
        }

    comp = np.array(report["receptors"]["9DKB"]["crystal_ligand_centroid"])
    report["receptors"]["9DKB"]["vs_locked_center_abs_delta"] = [
        round(float(x), 3) for x in np.abs(comp - LOCKED_9DKB)
    ]
    report["receptors"]["9DKB"]["use"] = "locked production center from docking_c1.yaml"
    report["receptors"]["9DKB"]["note"] = (
        "Fresh heavy-atom centroid differs from the locked production box by "
        "0.076 / 0.688 / 0.229 A. Keep the locked docking_c1.yaml center for all "
        "9DKB jobs so C1 self-dock / acid-track poses stay in the same box. "
        "The 0.7 A offset is well inside the 22 A cube."
    )

    stb = load_cif("9DKB")
    tb = ca_map(stb)
    st9 = load_cif("9DK9")
    m9 = ca_map(st9)
    keys = sorted(set(tb) & set(m9))
    r, t, rmsd = kabsch(np.vstack([tb[k] for k in keys]), np.vstack([m9[k] for k in keys]))
    center_9dk9 = r @ LOCKED_9DKB + t
    report["receptors"]["9DK9"] = {
        "ligand_ccd": None,
        "state": "apo",
        "box_rule": "transfer_9DKB_locked_center_after_chainA_CA_superposition_9DKB_onto_9DK9",
        "center": [round(float(x), 3) for x in center_9dk9],
        "size": [22, 22, 22],
        "superposition_rmsd_A": round(rmsd, 3),
        "n_ca_pairs": len(keys),
        "do_not": "copy 9DKB Cartesian center without superposition (different coordinate frame)",
    }

    for pdb_id in ("9DKA", "9DKC", "9DK9"):
        mm = ca_map(load_cif(pdb_id))
        keys = sorted(set(tb) & set(mm))
        r, t, rmsd = kabsch(np.vstack([mm[k] for k in keys]), np.vstack([tb[k] for k in keys]))
        report["receptors"][pdb_id]["superpose_onto_9DKB_ca_rmsd_A"] = round(rmsd, 3)
        report["receptors"][pdb_id]["n_ca_pairs_vs_9DKB"] = len(keys)
        if pdb_id != "9DK9":
            native = np.array(report["receptors"][pdb_id]["crystal_ligand_centroid"])
            mapped = r @ native + t
            report["receptors"][pdb_id]["centroid_mapped_into_9DKB_frame"] = [
                round(float(x), 3) for x in mapped
            ]
            report["receptors"][pdb_id]["native_centroid_vs_9DKB_locked_abs_delta"] = [
                round(float(x), 3) for x in np.abs(native - LOCKED_9DKB)
            ]
            report["receptors"][pdb_id]["use_native_crystal_centroid"] = True
            report["receptors"][pdb_id]["do_not"] = (
                "do not copy the 9DKB Cartesian center onto 9DKA/9DKC"
            )

    st7 = gemmi.read_structure(str(PDB_DIR / "7ALV.pdb"))
    atoms = ligand_heavy_atoms(st7, "RM5")
    cen = centroid(atoms)
    report["7ALV_RM5"] = {
        "n_heavy": len(atoms),
        "centroid": [round(float(x), 3) for x in cen],
        "vs_locked_abs_delta": [round(float(x), 3) for x in np.abs(cen - LOCKED_7ALV)],
        "use": "locked production center from docking_c1.yaml",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
