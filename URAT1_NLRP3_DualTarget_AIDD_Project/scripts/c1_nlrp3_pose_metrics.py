#!/usr/bin/env python3
"""NLRP3 (7ALV) structural pose metrics vs NP3-146/RM5.

Adds (beyond loose COM gate):
  1. pocket volume overlap with NP3-146 heavy atoms
  2. interaction fingerprint (IFP) Jaccard vs crystal NP3-146
  3. key-residue contact recovery (Ala227/228, Arg351, Met408, Tyr443, Phe575, Arg578)
  4. helpers for multi-seed robustness (caller aggregates)

Claim language: pocket / IFP compatibility — NOT direct NLRP3 binding proof.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem

from parse_c1_sdf_readouts import _fprop, heavy_centroid, load_poses, load_ref_centroid

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFS = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
KEY_JSON = REFS / "nlrp3_key_residues.json"
REF_SDF = REFS / "NP3-146_RM5_crystal_ref.sdf"
REC_PDBQT = PROJECT_ROOT / "data/structures/prepared/7ALV_receptor.pdbqt"

CENTROID_MAX_A = 6.0
NLRP3_CNN_MIN = 0.5
OVERLAP_CUTOFF_A = 2.5
CONTACT_CUTOFF_A = 4.5
HBOND_CUTOFF_A = 3.5
CLASH_CUTOFF_A = 2.2
# Structural gate (NLRP3 upgrade 2026-08-27): stricter than loose COM+CNNscore.
# Anchored below NP3-146 self-dock (overlap≈1.0, IFP≈0.84–1.0, key 6–7/7).
MIN_OVERLAP = 0.50
MIN_IFP_JACCARD = 0.50
MIN_KEY_RECOVERY = 5 / 7  # at least 5 of 7 key residues contacted


def _mol_heavy_xyz(mol: Chem.Mol) -> tuple[np.ndarray, list[str], np.ndarray]:
    conf = mol.GetConformer()
    xyz, elems, is_h = [], [], []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        e = atom.GetSymbol()
        xyz.append([p.x, p.y, p.z])
        elems.append(e)
        is_h.append(e == "H")
    return np.array(xyz, dtype=float), elems, np.array(is_h, dtype=bool)


def _pairwise_min(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2).min(axis=1)


def pocket_overlap_frac(lig_heavy: np.ndarray, ref_heavy: np.ndarray, cutoff: float = OVERLAP_CUTOFF_A) -> float:
    if len(lig_heavy) == 0 or len(ref_heavy) == 0:
        return float("nan")
    return float((_pairwise_min(lig_heavy, ref_heavy) < cutoff).mean())


def load_key_map(path: Path = KEY_JSON) -> dict:
    return json.loads(path.read_text())["residues"]


def load_receptor_heavy(path: Path = REC_PDBQT) -> np.ndarray:
    pts = []
    with path.open() as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            if name.startswith("H") or name.upper().startswith("HD"):
                continue
            pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(pts, dtype=float)


def residue_min_dist(lig_heavy: np.ndarray, res_atoms: list[dict]) -> float:
    coords = np.array([a["xyz"] for a in res_atoms], dtype=float)
    if len(lig_heavy) == 0 or len(coords) == 0:
        return float("inf")
    return float(_pairwise_min(lig_heavy, coords).min())


def key_residue_contacts(lig_heavy: np.ndarray, key_map: dict, cutoff: float = CONTACT_CUTOFF_A) -> dict:
    hits = {}
    for label, info in key_map.items():
        d = residue_min_dist(lig_heavy, info["heavy_atoms"])
        hits[label] = {"min_dist_A": d, "contact": bool(d <= cutoff)}
    return hits


def interaction_fingerprint(
    lig_heavy: np.ndarray,
    lig_elems: list[str],
    lig_is_h: np.ndarray,
    key_map: dict,
    receptor_heavy: np.ndarray,
) -> set[str]:
    """Compact structural IFP bits over key residues + global clash/hbond proxies."""
    bits: set[str] = set()
    heavy = lig_heavy[~lig_is_h] if lig_is_h.any() else lig_heavy
    elems = [e for e, h in zip(lig_elems, lig_is_h) if not h]
    contacts = key_residue_contacts(heavy, key_map)
    for label, info in contacts.items():
        if info["contact"]:
            bits.add(f"CONTACT:{label}")
            # polarity class of residue
            if label.startswith("ARG"):
                bits.add(f"POLAR:{label}")
            elif label.startswith(("PHE", "TYR", "MET")):
                bits.add(f"HYDROPHOBIC:{label}")
            elif label.startswith("ALA"):
                bits.add(f"WALKER:{label}")
            if label.startswith("TYR") or label.startswith("PHE"):
                bits.add(f"AROMATIC:{label}")
    # H-bond proxy: N/O of ligand within 3.5 Å of key residue N/O
    lig_no = np.array([heavy[i] for i, e in enumerate(elems) if e in {"N", "O"}], dtype=float)
    for label, info in key_map.items():
        no_atoms = [a for a in info["heavy_atoms"] if a["elem"] in {"N", "O"}]
        if len(lig_no) == 0 or not no_atoms:
            continue
        coords = np.array([a["xyz"] for a in no_atoms], dtype=float)
        if float(_pairwise_min(lig_no, coords).min()) <= HBOND_CUTOFF_A:
            bits.add(f"HBOND:{label}")
    # clash bit (global)
    if len(heavy) and len(receptor_heavy):
        if float((_pairwise_min(heavy, receptor_heavy) < CLASH_CUTOFF_A).sum()) > 0:
            bits.add("CLASH:ANY")
    return bits


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return float("nan")
    inter = len(a & b)
    union = len(a | b)
    return float(inter / union) if union else float("nan")


def evaluate_nlrp3_structural(
    sdf: Path,
    ligand_id: str,
    seed: int,
    key_map: dict | None = None,
    ref_heavy: np.ndarray | None = None,
    ref_com: tuple[float, float, float] | None = None,
    receptor_heavy: np.ndarray | None = None,
    ref_ifp: set[str] | None = None,
) -> dict:
    key_map = key_map or load_key_map()
    if ref_heavy is None:
        ref_mol = load_poses(REF_SDF)[0]
        ref_xyz, _, ref_h = _mol_heavy_xyz(ref_mol)
        ref_heavy = ref_xyz[~ref_h]
    if ref_com is None:
        ref_com = load_ref_centroid(REF_SDF)
    if receptor_heavy is None:
        receptor_heavy = load_receptor_heavy()

    poses = load_poses(sdf)
    if not poses:
        return {
            "ligand_id": ligand_id,
            "target": "nlrp3_7alv",
            "seed": seed,
            "error": "no_poses",
            "keep_nlrp3_pose": False,
            "keep_nlrp3_structural": False,
        }
    i_star = max(range(len(poses)), key=lambda j: _fprop(poses[j], "CNNscore") or -1.0)
    pose = poses[i_star]
    xyz, elems, is_h = _mol_heavy_xyz(pose)
    heavy = xyz[~is_h]
    pose_com = heavy_centroid(pose)
    d_com = float(
        math.sqrt(sum((pose_com[i] - ref_com[i]) ** 2 for i in range(3)))
    )
    in_pocket = d_com <= CENTROID_MAX_A
    cnn = _fprop(pose, "CNNscore")
    cnna = _fprop(pose, "CNNaffinity")

    overlap = pocket_overlap_frac(heavy, ref_heavy)
    contacts = key_residue_contacts(heavy, key_map)
    n_key = sum(1 for v in contacts.values() if v["contact"])
    key_recovery = n_key / max(len(key_map), 1)
    ifp = interaction_fingerprint(xyz, elems, is_h, key_map, receptor_heavy)
    if ref_ifp is None:
        # build crystal IFP once
        ref_mol = load_poses(REF_SDF)[0]
        rxyz, relems, rh = _mol_heavy_xyz(ref_mol)
        ref_ifp = interaction_fingerprint(rxyz, relems, rh, key_map, receptor_heavy)
    ifp_jac = jaccard(ifp, ref_ifp)

    keep_loose = bool(in_pocket and (cnn or 0) >= NLRP3_CNN_MIN)
    keep_structural = bool(
        keep_loose
        and (overlap >= MIN_OVERLAP)
        and (ifp_jac >= MIN_IFP_JACCARD)
        and (key_recovery >= MIN_KEY_RECOVERY)
        and ("CLASH:ANY" not in ifp)
    )

    row = {
        "ligand_id": ligand_id,
        "target": "nlrp3_7alv",
        "seed": seed,
        "n_poses": len(poses),
        "selected_mode": i_star + 1,
        "CNNscore": cnn,
        "CNNaffinity": cnna,
        "centroid_to_crystal_lig_A": d_com,
        "pass_pocket_centroid": in_pocket,
        "pocket_overlap_frac": overlap,
        "n_key_contacts": n_key,
        "key_recovery_frac": key_recovery,
        "ifp_jaccard_vs_np3146": ifp_jac,
        "ifp_bits": ";".join(sorted(ifp)),
        "n_ifp_bits": len(ifp),
        "keep_nlrp3_pose": keep_loose,  # legacy A1/A2 loose gate
        "keep_nlrp3_structural": keep_structural,
        "sdf": str(sdf),
    }
    for label, info in contacts.items():
        row[f"d_{label}"] = info["min_dist_A"]
        row[f"contact_{label}"] = info["contact"]
    return row


def crystal_reference_ifp() -> tuple[set[str], np.ndarray, tuple[float, float, float]]:
    key_map = load_key_map()
    receptor_heavy = load_receptor_heavy()
    ref_mol = load_poses(REF_SDF)[0]
    rxyz, relems, rh = _mol_heavy_xyz(ref_mol)
    ref_heavy = rxyz[~rh]
    ref_ifp = interaction_fingerprint(rxyz, relems, rh, key_map, receptor_heavy)
    ref_com = load_ref_centroid(REF_SDF)
    return ref_ifp, ref_heavy, ref_com
