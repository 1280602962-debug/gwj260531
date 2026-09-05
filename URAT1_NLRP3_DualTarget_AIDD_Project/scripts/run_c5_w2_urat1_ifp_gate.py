#!/usr/bin/env python3
"""C5 W2: URAT1 acid-anchor IFP gate (0 new docks).

Re-score Phase I 9DKB 9-mode SDFs. Thresholds are anchored on deposited crystal
ligand coordinates (R75 / A1AIL / A1A45), never on CNNscore-selected re-docks.

Protocol (isomorphic to NLRP3 structural gate):
  acid_O–Arg477_N_min + pocket overlap + key-residue IFP Jaccard + no clash

Pose selection (primary, NLRP3-symmetric): CNNscore Top-1, then structural checks.
Diagnostics: best-of-9 IFP Jaccard among modes.

Pass rule (pre-registered): OR 95% CI lower bound > 1 on 228 vs 64 once.
On fail: fall back to A1 ∩ A2; do not invent a new metric.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import fisher_exact

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from c1_acid_pose_selection import (  # noqa: E402
    ARG_THRESH_A,
    evaluate_urat1_acid_sdf,
)
from c1_nlrp3_pose_metrics import (  # noqa: E402
    CLASH_CUTOFF_A,
    _mol_heavy_xyz,
    _pairwise_min,
    interaction_fingerprint,
    jaccard,
    key_residue_contacts,
    load_receptor_heavy,
    pocket_overlap_frac,
)
from parse_c1_sdf_readouts import (  # noqa: E402
    _fprop,
    carboxylate_oxygens,
    heavy_centroid,
    load_poses,
    min_acid_arg_dist,
)
from run_acid_gate_benchmark import bootstrap_or_ci, has_acid, metrics_from_pass  # noqa: E402

OUT = PROJECT_ROOT / "data/campaigns/c5/02_urat1_ifp"
KEY_JSON = OUT / "urat1_key_residues.json"
POSE_ROOT = (
    PROJECT_ROOT
    / "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/poses/gnina_sdf"
)
MAP_CSV = (
    PROJECT_ROOT
    / "docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/meta/mol_index_map.csv"
)
ARG_JSON = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs/arg477_coords.json"
REC_PDBQT = PROJECT_ROOT / "data/structures/prepared/9DKB_receptor.pdbqt"

# Crystal refs (native frames). A1AIL already in 9DKB; R75/A1A45 transformed.
CRYSTAL_REFS = {
    "A1AIL": {
        "sdf": PROJECT_ROOT
        / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/lesinurad_A1AIL_crystal_ref.sdf",
        "src_struct": PROJECT_ROOT / "data/structures/pdb/9DKB.cif",
        "need_transform": False,
    },
    "R75": {
        "sdf": PROJECT_ROOT
        / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/benzbromarone_R75_crystal_ref.sdf",
        "src_struct": PROJECT_ROOT / "data/structures/pdb/9DKA.cif",
        "need_transform": True,
    },
    "A1A45": {
        "sdf": PROJECT_ROOT
        / "data/campaigns/c5/01_ligand_prep/w1_crystal_refs/TD3_A1A45_crystal_ref.sdf",
        "src_struct": PROJECT_ROOT / "data/structures/pdb/9DKC.cif",
        "need_transform": True,
    },
}
TARGET_STRUCT = PROJECT_ROOT / "data/structures/pdb/9DKB.cif"

# Floors isomorphic to NLRP3 (below crystal self-metrics of 1.0 / full contacts).
MIN_OVERLAP = 0.50
MIN_IFP_JACCARD = 0.50
MIN_KEY_FRAC = 0.50
CLASH_FORBIDDEN = True


def chain_a_ca(path: Path) -> dict[int, np.ndarray]:
    st = gemmi.read_structure(str(path))
    out: dict[int, np.ndarray] = {}
    for res in st[0]["A"]:
        try:
            resi = int(res.seqid.num)
        except Exception:
            continue
        for atom in res:
            if atom.name.strip() == "CA":
                out[resi] = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
                break
    return out


def kabsch_R_t(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return R, t such that P @ R.T + t ≈ Q (map P-frame → Q-frame)."""
    pc, qc = P.mean(0), Q.mean(0)
    P0, Q0 = P - pc, Q - qc
    H = P0.T @ Q0
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = qc - pc @ R.T
    return R, t


def transform_xyz(xyz: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return xyz @ R.T + t


def ca_alignment(src: Path, dst: Path) -> tuple[np.ndarray, np.ndarray, float, int]:
    a, b = chain_a_ca(src), chain_a_ca(dst)
    common = sorted(set(a) & set(b))
    P = np.vstack([a[i] for i in common])
    Q = np.vstack([b[i] for i in common])
    R, t = kabsch_R_t(P, Q)
    aligned = transform_xyz(P, R, t)
    rmsd = float(np.sqrt(((aligned - Q) ** 2).sum(axis=1).mean()))
    return R, t, rmsd, len(common)


def phenolate_or_carboxylate_oxygens(mol: Chem.Mol) -> list:
    """Carboxylate O preferred; else aromatic/phenolic O- as acid proxy (benzbromarone)."""
    oxys = carboxylate_oxygens(mol)
    if oxys:
        return oxys
    conf = mol.GetConformer()
    out = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        if atom.GetFormalCharge() < 0 or atom.GetTotalNumHs() == 0 and atom.GetDegree() == 1:
            # charged O or keto/phenolate-like terminal O
            if atom.GetFormalCharge() < 0:
                p = conf.GetAtomPosition(atom.GetIdx())
                out.append(np.array([p.x, p.y, p.z], dtype=float))
    return out


def acid_arg_min(mol: Chem.Mol, arg_atoms: dict) -> float | None:
    # Prefer standard carboxylate path; fall back to phenolate O
    d = min_acid_arg_dist(mol, arg_atoms)
    if d is not None:
        return d
    oxys = phenolate_or_carboxylate_oxygens(mol)
    if not oxys:
        return None
    arg_xyz = np.array(list(arg_atoms.values()), dtype=float)
    return float(min(np.linalg.norm(o - arg_xyz, axis=1).min() for o in oxys))


def load_key_map() -> dict:
    return json.loads(KEY_JSON.read_text())["residues"]


def mol_from_sdf_transformed(sdf: Path, R: np.ndarray | None, t: np.ndarray | None) -> Chem.Mol:
    mol = load_poses(sdf)[0]
    if R is None:
        return mol
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        v = transform_xyz(np.array([[p.x, p.y, p.z]]), R, t)[0]
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(v[0]), float(v[1]), float(v[2])))
    return mol


def build_crystal_anchors(key_map: dict, receptor_heavy: np.ndarray, arg_atoms: dict) -> dict:
    refs = {}
    transforms = {}
    for name, meta in CRYSTAL_REFS.items():
        R = t = None
        rmsd = 0.0
        n_ca = 0
        if meta["need_transform"]:
            R, t, rmsd, n_ca = ca_alignment(meta["src_struct"], TARGET_STRUCT)
            transforms[name] = {"ca_rmsd_A": rmsd, "n_ca": n_ca, "src": str(meta["src_struct"])}
        mol = mol_from_sdf_transformed(meta["sdf"], R, t)
        xyz, elems, is_h = _mol_heavy_xyz(mol)
        heavy = xyz[~is_h]
        ifp = interaction_fingerprint(xyz, elems, is_h, key_map, receptor_heavy)
        contacts = key_residue_contacts(heavy, key_map)
        n_key = sum(1 for v in contacts.values() if v["contact"])
        d_arg = acid_arg_min(mol, arg_atoms)
        refs[name] = {
            "n_heavy": int(len(heavy)),
            "n_key_contacts": n_key,
            "key_frac": n_key / max(len(key_map), 1),
            "acid_arg477_min_A": d_arg,
            "ifp_bits": sorted(ifp),
            "n_ifp_bits": len(ifp),
            "com": heavy.mean(axis=0).tolist(),
            "heavy": heavy,
            "ifp": ifp,
            "transformed": bool(meta["need_transform"]),
            "ca_rmsd_A": rmsd,
        }
    # Reference IFP = union of crystal bits (covers multi-chemotype pocket modes).
    ref_ifp: set[str] = set()
    for r in refs.values():
        ref_ifp |= set(r["ifp_bits"])
    # Overlap reference heavy = concatenate of three (pocket envelope).
    ref_heavy = np.vstack([r["heavy"] for r in refs.values()])
    # Arg gate for Phase-I@9DKB: lock to pre-registered crystal-relative threshold
    # (lesinurad A1AIL min 6.70 + 1.0 = 7.7027). Do not loosen with Kabsch-transformed
    # A1A45/R75 O–Arg distances (different holo frames / phenolate).
    arg_max = float(ARG_THRESH_A)
    thresholds = {
        "min_overlap": MIN_OVERLAP,
        "min_ifp_jaccard": MIN_IFP_JACCARD,
        "min_key_frac": MIN_KEY_FRAC,
        "min_key_contacts": int(math.ceil(MIN_KEY_FRAC * len(key_map))),
        "arg_max_A": arg_max,
        "clash_cutoff_A": CLASH_CUTOFF_A,
        "anchor_method": "deposited_crystal_R75_A1AIL_A1A45_kabsch_to_9DKB",
        "crystal_key_fracs": {k: refs[k]["key_frac"] for k in refs},
        "crystal_arg_A": {k: refs[k]["acid_arg477_min_A"] for k in refs},
        "arg_anchor": "ARG_THRESH_A_7.7027_from_A1AIL_crystal_plus_1A",
        "arg_anchor_note": (
            "R75 phenolate and Kabsch-A1A45 O–Arg reported for audit only; "
            "do not inflate 9DKB arg_max"
        ),
        "note": "Thresholds fixed from protocol floors below crystal self-metrics; not tuned on 228vs64.",
    }
    # strip heavy/ifp objects for JSON
    refs_json = {
        k: {kk: vv for kk, vv in v.items() if kk not in {"heavy", "ifp"}}
        for k, v in refs.items()
    }
    return {
        "refs": refs,
        "refs_json": refs_json,
        "transforms": transforms,
        "ref_ifp": ref_ifp,
        "ref_heavy": ref_heavy,
        "thresholds": thresholds,
    }


def score_pose(
    pose: Chem.Mol,
    key_map: dict,
    receptor_heavy: np.ndarray,
    ref_heavy: np.ndarray,
    ref_ifp: set[str],
    arg_atoms: dict,
    thresholds: dict,
) -> dict:
    xyz, elems, is_h = _mol_heavy_xyz(pose)
    heavy = xyz[~is_h]
    overlap = pocket_overlap_frac(heavy, ref_heavy)
    contacts = key_residue_contacts(heavy, key_map)
    n_key = sum(1 for v in contacts.values() if v["contact"])
    key_frac = n_key / max(len(key_map), 1)
    ifp = interaction_fingerprint(xyz, elems, is_h, key_map, receptor_heavy)
    ifp_jac = jaccard(ifp, ref_ifp)
    d_arg = acid_arg_min(pose, arg_atoms)
    clash = "CLASH:ANY" in ifp
    pass_arg = d_arg is not None and d_arg <= thresholds["arg_max_A"]
    pass_overlap = overlap >= thresholds["min_overlap"]
    pass_ifp = ifp_jac >= thresholds["min_ifp_jaccard"]
    pass_key = key_frac >= thresholds["min_key_frac"]
    pass_clash = (not clash) if CLASH_FORBIDDEN else True
    keep = bool(pass_arg and pass_overlap and pass_ifp and pass_key and pass_clash)
    return {
        "acid_arg477_min_A": d_arg,
        "pocket_overlap_frac": overlap,
        "ifp_jaccard_vs_crystal_union": ifp_jac,
        "n_key_contacts": n_key,
        "key_frac": key_frac,
        "clash": clash,
        "pass_arg": pass_arg,
        "pass_overlap": pass_overlap,
        "pass_ifp": pass_ifp,
        "pass_key": pass_key,
        "pass_clash": pass_clash,
        "keep_urat1_ifp": keep,
        "CNNscore": _fprop(pose, "CNNscore"),
        "CNNaffinity": _fprop(pose, "CNNaffinity"),
        "ifp_bits": ";".join(sorted(ifp)),
    }


def evaluate_sdf(
    sdf: Path,
    ligand_id: str,
    key_map: dict,
    receptor_heavy: np.ndarray,
    ref_heavy: np.ndarray,
    ref_ifp: set[str],
    arg_atoms: dict,
    thresholds: dict,
) -> dict:
    poses = load_poses(sdf)
    row = {
        "ligand_id": ligand_id,
        "sdf": str(sdf),
        "n_poses": len(poses),
        "keep_urat1_ifp": False,
        "keep_urat1_ifp_best_of_9": False,
    }
    if not poses:
        row["error"] = "no_poses"
        return row
    i_star = max(range(len(poses)), key=lambda j: _fprop(poses[j], "CNNscore") or -1.0)
    top = score_pose(
        poses[i_star], key_map, receptor_heavy, ref_heavy, ref_ifp, arg_atoms, thresholds
    )
    row["selected_mode_cnn"] = i_star + 1
    for k, v in top.items():
        row[k] = v
    # best-of-9 by IFP Jaccard (diagnostic)
    best_j, best_i, best_sc = -1.0, None, None
    for j, pose in enumerate(poses):
        sc = score_pose(pose, key_map, receptor_heavy, ref_heavy, ref_ifp, arg_atoms, thresholds)
        jac = sc["ifp_jaccard_vs_crystal_union"]
        if jac == jac and jac > best_j:  # not nan
            best_j, best_i, best_sc = jac, j, sc
    if best_sc is not None:
        row["best_of_9_mode"] = best_i + 1
        row["best_of_9_ifp_jaccard"] = best_sc["ifp_jaccard_vs_crystal_union"]
        row["keep_urat1_ifp_best_of_9"] = bool(best_sc["keep_urat1_ifp"])
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # 1) key residues
    if not KEY_JSON.exists():
        import subprocess

        subprocess.check_call(
            [sys.executable, str(PROJECT_ROOT / "scripts/extract_urat1_key_residues.py")]
        )
    key_map = load_key_map()
    receptor_heavy = load_receptor_heavy(REC_PDBQT)
    arg_atoms = json.loads(ARG_JSON.read_text())["atoms"]

    anchors = build_crystal_anchors(key_map, receptor_heavy, arg_atoms)
    thresholds = anchors["thresholds"]
    (OUT / "w2_crystal_anchor_metrics.json").write_text(
        json.dumps(
            {
                "thresholds": thresholds,
                "crystal_refs": anchors["refs_json"],
                "transforms": anchors["transforms"],
                "ref_ifp_bits": sorted(anchors["ref_ifp"]),
                "n_key_mapped": len(key_map),
                "key_labels": list(key_map),
            },
            indent=2,
        )
        + "\n"
    )

    # 2) pool 228 vs 64
    act = pd.read_csv(PROJECT_ROOT / "data/benchmarks/urat1_true_decoy/actives.csv")
    dec = pd.read_csv(PROJECT_ROOT / "data/benchmarks/urat1_true_decoy/true_decoys.csv")
    act = act[act.canonical_smiles.map(has_acid)].copy()
    dec = dec[dec.canonical_smiles.map(has_acid)].copy()
    act["label"] = 1
    dec["label"] = 0
    pool = pd.concat([act, dec], ignore_index=True).drop_duplicates(subset=["canonical_smiles"])
    mmap = pd.read_csv(MAP_CSV)
    smi2mol = dict(zip(mmap.canonical_smiles, mmap.mol_id))
    pool["mol_id"] = pool.canonical_smiles.map(smi2mol)
    pool = pool[pool.mol_id.notna()].copy()
    if args.limit:
        pool = pool.head(args.limit)

    rows = []
    for i, r in pool.iterrows():
        sdf = POSE_ROOT / f"{r['mol_id']}_out.sdf"
        if not sdf.exists():
            continue
        ev = evaluate_sdf(
            sdf,
            str(r["mol_id"]),
            key_map,
            receptor_heavy,
            anchors["ref_heavy"],
            anchors["ref_ifp"],
            arg_atoms,
            thresholds,
        )
        ev["label"] = int(r["label"])
        ev["canonical_smiles"] = r["canonical_smiles"]
        ev["set_role"] = "active" if r["label"] == 1 else "true_decoy"
        rows.append(ev)
        if len(rows) % 50 == 0:
            print(f"... {len(rows)} molecules", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "w2_ifp_gate_per_mol.csv", index=False)

    # 3) also A1/A2 for fallback intersection
    ref_com = heavy_centroid(load_poses(CRYSTAL_REFS["A1AIL"]["sdf"])[0])
    a_rows = []
    for _, r in pool.iterrows():
        sdf = POSE_ROOT / f"{r['mol_id']}_out.sdf"
        if not sdf.exists():
            continue
        for rule in ("a1", "a2"):
            ev = evaluate_urat1_acid_sdf(sdf, ARG_JSON, ref_com, str(r["mol_id"]), 42, rule=rule)
            ev["label"] = int(r["label"])
            a_rows.append(ev)
    adf = pd.DataFrame(a_rows)
    a1 = adf[adf.pose_selection_rule == "a1"].set_index("ligand_id")
    a2 = adf[adf.pose_selection_rule == "a2"].set_index("ligand_id")
    df["keep_a1"] = df.ligand_id.map(lambda x: bool(a1.loc[x, "keep_urat1_acid"]) if x in a1.index else False)
    df["keep_a2"] = df.ligand_id.map(lambda x: bool(a2.loc[x, "keep_urat1_acid"]) if x in a2.index else False)
    df["keep_a1_and_a2"] = df["keep_a1"] & df["keep_a2"]
    df.to_csv(OUT / "w2_ifp_gate_per_mol.csv", index=False)

    summary: dict = {
        "subset": "carboxylate actives vs true decoys (Phase I 9DKB SDFs)",
        "pose_source": str(POSE_ROOT),
        "anchor": "deposited_crystal_R75_A1AIL_A1A45",
        "n_key_mapped": len(key_map),
        "thresholds": thresholds,
        "elapsed_sec": time.time() - t0,
    }
    for name, col in [
        ("ifp_cnn_top1", "keep_urat1_ifp"),
        ("ifp_best_of_9", "keep_urat1_ifp_best_of_9"),
        ("a1", "keep_a1"),
        ("a2", "keep_a2"),
        ("a1_and_a2_fallback", "keep_a1_and_a2"),
    ]:
        summary[name] = metrics_from_pass(df.label.tolist(), df[col].tolist())

    ifp = summary["ifp_cnn_top1"]
    passed = bool(ifp["or_bootstrap_ci95_lo"] > 1.0)
    summary["pass_rule"] = "or_bootstrap_ci95_lo > 1.0 on ifp_cnn_top1"
    summary["gate_pass"] = passed
    summary["on_fail_action"] = None if passed else "fallback_to_a1_and_a2_intersection"
    if not passed:
        summary["fallback_metrics"] = summary["a1_and_a2_fallback"]

    (OUT / "w2_ifp_gate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    # short markdown report
    md = [
        "# C5 W2 URAT1 IFP gate result",
        "",
        f"- Anchor: deposited R75 / A1AIL / A1A45 (Kabsch→9DKB for R75 & A1A45).",
        f"- Key residues mapped: {len(key_map)}/12 (Q437 unmatched: LEU in 9DKB).",
        f"- Primary gate (CNNscore Top-1 + IFP structural): "
        f"OR={ifp['odds_ratio']:.3f}, CI95=[{ifp['or_bootstrap_ci95_lo']:.3f}, {ifp['or_bootstrap_ci95_hi']:.3f}], "
        f"Fisher p={ifp['fisher_exact_p']:.3g}.",
        f"- **gate_pass = {passed}** (need CI lower > 1).",
    ]
    if not passed:
        fb = summary["a1_and_a2_fallback"]
        md.append(
            f"- Fallback A1∩A2: OR={fb['odds_ratio']:.3f}, "
            f"CI95=[{fb['or_bootstrap_ci95_lo']:.3f}, {fb['or_bootstrap_ci95_hi']:.3f}]."
        )
    md.append("")
    (OUT / "W2_IFP_GATE_REPORT.md").write_text("\n".join(md) + "\n")
    print(json.dumps(summary, indent=2))
    print("WROTE", OUT / "w2_ifp_gate_summary.json")


if __name__ == "__main__":
    main()
