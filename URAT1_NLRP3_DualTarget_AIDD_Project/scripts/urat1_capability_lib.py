"""Shared inventory, RMSD, and search/selection labels for URAT1 capability docking.

D1 is the 9DKA/B/C × 3 ligands × 3 seeds matrix. 9DK9 apo is not in D1.
Selection failure is only assigned when a native-like pose exists in the nine modes.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import yaml
from rdkit import Chem

from parse_c1_sdf_readouts import (
    carboxylate_oxygens,
    load_poses,
    pose_rmsd,
)
from run_c5_w2_urat1_ifp_gate import ca_alignment, transform_xyz

PROJECT_ROOT = Path(__file__).resolve().parents[1]
NATIVE_LIKE = 2.0
SEEDS = (42, 43, 44)
D1_LIGANDS = ("benzbromarone", "lesinurad", "TD-3")
D1_TARGETS = ("urat1_9dka", "urat1_9dkb", "urat1_9dkc")
TARGET_PDB = {
    "urat1_9dka": "9DKA",
    "urat1_9dkb": "9DKB",
    "urat1_9dkc": "9DKC",
}
LIGAND_NATIVE_PDB = {
    "benzbromarone": "9DKA",
    "lesinurad": "9DKB",
    "TD-3": "9DKC",
}
LIGAND_PDBQT = {
    "lesinurad": PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/forced_recovery/pdbqt/lesinurad.pdbqt",
    "benzbromarone": PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/forced_recovery/pdbqt/benzbromarone.pdbqt",
    "TD-3": PROJECT_ROOT / "data/campaigns/c5/01_ligand_prep/w1_refs/pdbqt/TD-3.pdbqt",
    "verinurad": PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/forced_recovery/pdbqt/verinurad.pdbqt",
    "dotinurad": PROJECT_ROOT
    / "data/campaigns/c1/01_ligand_prep/forced_recovery/pdbqt/dotinurad.pdbqt",
}
CRYSTAL_REFS = {
    "lesinurad": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/lesinurad_A1AIL_crystal_ref.sdf",
    "benzbromarone": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/benzbromarone_R75_crystal_ref.sdf",
    "TD-3": PROJECT_ROOT
    / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/TD3_A1A45_crystal_ref.sdf",
}
STRUCTURE_CIF = {
    "9DKA": PROJECT_ROOT / "data/structures/pdb/9DKA.cif",
    "9DKB": PROJECT_ROOT / "data/structures/pdb/9DKB.cif",
    "9DKC": PROJECT_ROOT / "data/structures/pdb/9DKC.cif",
}
REUSE_SOURCES = {
    ("lesinurad", "urat1_9dkb", 42): PROJECT_ROOT
    / "data/campaigns/c1/02_selfdock/urat1_9dkb/seed42/lesinurad_out.sdf",
    ("lesinurad", "urat1_9dkb", 43): PROJECT_ROOT
    / "data/campaigns/c1/02_selfdock/urat1_9dkb/seed43/lesinurad_out.sdf",
    ("lesinurad", "urat1_9dkb", 44): PROJECT_ROOT
    / "data/campaigns/c1/02_selfdock/urat1_9dkb/seed44/lesinurad_out.sdf",
    ("benzbromarone", "urat1_9dka", 42): PROJECT_ROOT
    / "data/campaigns/c5/01_crossdock/9dka/seed42/benzbromarone_out.sdf",
    ("benzbromarone", "urat1_9dka", 43): PROJECT_ROOT
    / "data/campaigns/c5/01_crossdock/9dka/seed43/benzbromarone_out.sdf",
    ("benzbromarone", "urat1_9dka", 44): PROJECT_ROOT
    / "data/campaigns/c5/01_crossdock/9dka/seed44/benzbromarone_out.sdf",
    ("benzbromarone", "urat1_9dkb", 42): PROJECT_ROOT
    / "data/campaigns/c1/03_forced_recovery/urat1_9dkb/seed42/benzbromarone_out.sdf",
}
D1_OUT = PROJECT_ROOT / "data/campaigns/c5/06_capability_dock/d1"
REF_CACHE = PROJECT_ROOT / "data/campaigns/c5/06_capability_dock/refs_transformed"


def load_capability_config() -> dict:
    return yaml.safe_load((PROJECT_ROOT / "config/docking_capability.yaml").read_text())


def load_w1_engine_config() -> dict:
    return yaml.safe_load((PROJECT_ROOT / "config/docking_c5_w1.yaml").read_text())


def d1_jobs() -> list[tuple[str, str, int]]:
    return [
        (lig, tgt, seed)
        for lig in D1_LIGANDS
        for tgt in D1_TARGETS
        for seed in SEEDS
    ]


def d1_new_jobs() -> list[tuple[str, str, int]]:
    return [j for j in d1_jobs() if j not in REUSE_SOURCES]


def d1_selfdock(ligand: str, target: str) -> bool:
    return LIGAND_NATIVE_PDB[ligand] == TARGET_PDB[target]


def cell_sdf(ligand: str, target: str, seed: int) -> Path:
    pdb = TARGET_PDB[target].lower()
    return D1_OUT / pdb / f"seed{seed}" / f"{ligand}_out.sdf"


def classify_search_selection(
    top1_rmsd: float | None,
    best_of_9: float | None,
    cutoff: float = NATIVE_LIKE,
) -> str:
    if top1_rmsd is None or best_of_9 is None:
        return "unevaluable"
    if top1_rmsd <= cutoff:
        return "selection_ok"
    if best_of_9 <= cutoff:
        return "search_ok_selection_fail"
    return "search_incomplete"


def success_at_k(rmsds_by_cnn_rank: list[float], k: int, cutoff: float = NATIVE_LIKE) -> bool:
    if not rmsds_by_cnn_rank:
        return False
    return any(r <= cutoff for r in rmsds_by_cnn_rank[:k])


def cnnscore_of(mol: Chem.Mol) -> float:
    if mol.HasProp("CNNscore"):
        return float(mol.GetProp("CNNscore"))
    raise KeyError("CNNscore")


def acid_proxy_oxygens(mol: Chem.Mol) -> list[np.ndarray]:
    oxy = carboxylate_oxygens(mol)
    if oxy:
        return [np.array(p, dtype=float) for p in oxy]
    conf = mol.GetConformer()
    out: list[np.ndarray] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O" or atom.GetFormalCharge() >= 0:
            continue
        p = conf.GetAtomPosition(atom.GetIdx())
        out.append(np.array([p.x, p.y, p.z], dtype=float))
    return out


def arg_guanidinium(cif: Path, resi: int = 477, chain: str = "A") -> dict[str, list[float]]:
    import gemmi

    st = gemmi.read_structure(str(cif))
    for res in st[0][chain]:
        try:
            num = int(res.seqid.num)
        except Exception:
            continue
        if num != resi or res.name.strip() != "ARG":
            continue
        atoms: dict[str, list[float]] = {}
        for atom in res:
            name = atom.name.strip()
            if name in {"NE", "NH1", "NH2"}:
                atoms[name] = [float(atom.pos.x), float(atom.pos.y), float(atom.pos.z)]
        if {"NE", "NH1", "NH2"} <= set(atoms):
            return atoms
    raise KeyError(f"ARG {chain}{resi} not found in {cif}")


def min_acid_arg(mol: Chem.Mol, arg_atoms: dict[str, list[float]]) -> float | None:
    oxys = acid_proxy_oxygens(mol)
    if not oxys:
        return None
    nitrogens = np.array([arg_atoms[k] for k in ("NE", "NH1", "NH2")], dtype=float)
    return float(min(np.linalg.norm(o - nitrogens, axis=1).min() for o in oxys))


def _copy_mol(mol: Chem.Mol) -> Chem.Mol:
    return Chem.Mol(mol)


def apply_rt_to_mol(mol: Chem.Mol, R: np.ndarray, t: np.ndarray) -> Chem.Mol:
    out = _copy_mol(mol)
    conf = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xyz = transform_xyz(np.array([[p.x, p.y, p.z]], dtype=float), R, t)[0]
        conf.SetAtomPosition(i, [float(xyz[0]), float(xyz[1]), float(xyz[2])])
    return out


def reference_in_receptor_frame(ligand: str, target: str) -> tuple[Path, dict]:
    """Crystal reference coordinates expressed in the docking-receptor frame."""
    native_pdb = LIGAND_NATIVE_PDB[ligand]
    dock_pdb = TARGET_PDB[target]
    src = CRYSTAL_REFS[ligand]
    meta = {
        "ligand": ligand,
        "native_pdb": native_pdb,
        "dock_pdb": dock_pdb,
        "self_dock": native_pdb == dock_pdb,
        "source_ref": str(src.relative_to(PROJECT_ROOT)),
    }
    if native_pdb == dock_pdb:
        meta["transform"] = "identity"
        return src, meta
    REF_CACHE.mkdir(parents=True, exist_ok=True)
    cached = REF_CACHE / f"{ligand}_onto_{dock_pdb}.sdf"
    R, t, ca_rmsd, n_ca = ca_alignment(STRUCTURE_CIF[native_pdb], STRUCTURE_CIF[dock_pdb])
    refs = load_poses(src)
    if not refs:
        raise ValueError(f"empty crystal ref {src}")
    moved = apply_rt_to_mol(refs[0], R, t)
    writer = Chem.SDWriter(str(cached))
    writer.write(moved)
    writer.close()
    meta.update(
        {
            "transform": "chainA_CA_Kabsch_native_onto_dock_receptor",
            "ca_rmsd_A": round(float(ca_rmsd), 4),
            "n_ca": int(n_ca),
            "transformed_ref": str(cached.relative_to(PROJECT_ROOT)),
        }
    )
    (cached.with_suffix(".json")).write_text(json.dumps(meta, indent=2) + "\n")
    return cached, meta


def evaluate_sdf(ligand: str, target: str, seed: int, sdf: Path) -> dict:
    row = {
        "layer": "d1",
        "ligand": ligand,
        "target": target,
        "pdb": TARGET_PDB[target],
        "seed": seed,
        "self_dock": d1_selfdock(ligand, target),
        "sdf": str(sdf.relative_to(PROJECT_ROOT)) if sdf.exists() else str(sdf),
        "status": "ok" if sdf.exists() and sdf.stat().st_size > 0 else "missing",
        "rmsd_definition": (
            "symmetry-aware heavy-atom pose RMSD in docking-receptor frame; "
            "no ligand superposition"
        ),
        "native_like_cutoff_A": NATIVE_LIKE,
    }
    if row["status"] != "ok":
        row["capability_class"] = "unevaluable"
        return row
    poses = load_poses(sdf)
    if not poses:
        row["status"] = "no_poses"
        row["capability_class"] = "unevaluable"
        return row
    scored: list[tuple[float, int, Chem.Mol]] = []
    for i, mol in enumerate(poses):
        try:
            scored.append((cnnscore_of(mol), i, mol))
        except Exception:
            continue
    if not scored:
        row["status"] = "no_CNNscore"
        row["capability_class"] = "unevaluable"
        return row
    scored.sort(key=lambda x: x[0], reverse=True)
    ref_path, ref_meta = reference_in_receptor_frame(ligand, target)
    row["reference"] = ref_meta
    ref_mols = load_poses(ref_path)
    if not ref_mols:
        row["status"] = "ref_empty"
        row["capability_class"] = "unevaluable"
        return row
    ref = ref_mols[0]
    rmsds: list[float] = []
    for _score, _idx, mol in scored:
        try:
            rmsds.append(float(pose_rmsd(mol, ref)))
        except Exception:
            rmsds.append(float("nan"))
    finite = [(rank, r) for rank, r in enumerate(rmsds, start=1) if math.isfinite(r)]
    if not finite:
        row["status"] = "rmsd_failed"
        row["capability_class"] = "unevaluable"
        return row
    best_rank, best_rmsd = min(finite, key=lambda x: x[1])
    top1 = finite[0][1] if finite[0][0] == 1 else rmsds[0]
    row.update(
        {
            "n_poses": len(poses),
            "n_scored_poses": len(scored),
            "CNNscore_top1": scored[0][0],
            "top1_pose_rmsd_A": float(top1),
            "best_of_9_pose_rmsd_A": float(best_rmsd),
            "best_rmsd_cnnscore_rank": int(best_rank),
            "success_at_1": bool(success_at_k(rmsds, 1)),
            "success_at_3": bool(success_at_k(rmsds, 3)),
            "success_at_5": bool(success_at_k(rmsds, 5)),
            "success_at_9": bool(success_at_k(rmsds, 9)),
            "capability_class": classify_search_selection(float(top1), float(best_rmsd)),
        }
    )
    try:
        arg = arg_guanidinium(STRUCTURE_CIF[TARGET_PDB[target]])
        row["acid_arg_top1_A"] = min_acid_arg(scored[0][2], arg)
        best_mol = scored[best_rank - 1][2]
        row["acid_arg_best_rmsd_A"] = min_acid_arg(best_mol, arg)
    except Exception as exc:
        row["acid_arg_error"] = str(exc)
    row["per_mode"] = [
        {
            "cnnscore_rank": rank,
            "sdf_index": idx,
            "CNNscore": score,
            "pose_rmsd_A": rmsds[rank - 1],
        }
        for rank, (score, idx, _mol) in enumerate(scored, start=1)
    ]
    return row


def inventory_rows() -> list[dict]:
    rows = []
    for lig, tgt, seed in d1_jobs():
        reuse = (lig, tgt, seed) in REUSE_SOURCES
        rows.append(
            {
                "layer": "d1",
                "ligand": lig,
                "target": tgt,
                "pdb": TARGET_PDB[tgt],
                "seed": seed,
                "self_dock": d1_selfdock(lig, tgt),
                "action": "reuse" if reuse else "new_gnina",
                "reuse_source": (
                    str(REUSE_SOURCES[(lig, tgt, seed)].relative_to(PROJECT_ROOT))
                    if reuse
                    else ""
                ),
                "output_sdf": str(cell_sdf(lig, tgt, seed).relative_to(PROJECT_ROOT)),
            }
        )
    return rows
