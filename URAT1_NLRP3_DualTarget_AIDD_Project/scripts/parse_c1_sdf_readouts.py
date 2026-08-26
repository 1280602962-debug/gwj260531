#!/usr/bin/env python3
"""Parse multi-pose GNINA SDF into C1 readouts and L2 self-dock metrics.

C1_P2star = CNNaffinity of pose with max CNNscore
C1_VS     = CNN_VS of that pose (or CNNscore*CNNaffinity)
C1_P2max  = max CNNaffinity
C1_P0     = max CNNscore

Batch CSV from run_gnina_batch.py is C1_P2max only — do not use it for Rank pass.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _fprop(mol: Chem.Mol, *keys: str) -> float | None:
    for k in keys:
        if mol.HasProp(k):
            try:
                return float(mol.GetProp(k))
            except Exception:
                continue
    return None


def load_poses(sdf: Path) -> list[Chem.Mol]:
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    return [m for m in suppl if m is not None]


def _neutralize_for_match(mol: Chem.Mol) -> Chem.Mol:
    """Drop formal charges so protonation microspecies can still graph-match."""
    m = Chem.RemoveHs(Chem.Mol(mol))
    for atom in m.GetAtoms():
        atom.SetFormalCharge(0)
        atom.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    return m


def pose_rmsd(docked: Chem.Mol, ref: Chem.Mol) -> float:
    d = _neutralize_for_match(docked)
    r = _neutralize_for_match(ref)
    matches = r.GetSubstructMatches(d, uniquify=False)
    if not matches:
        matches = d.GetSubstructMatches(r, uniquify=False)
        inv = []
        for m in matches:
            mapping = [None] * len(m)
            for qi, di in enumerate(m):
                mapping[di] = qi
            inv.append(tuple(mapping))
        matches = inv
    best = float("inf")
    cd, cr = d.GetConformer(), r.GetConformer()
    n = d.GetNumAtoms()
    for match in matches:
        if any(x is None for x in match) or len(match) != n:
            continue
        s = 0.0
        for i, j in enumerate(match):
            a, b = cd.GetAtomPosition(i), cr.GetAtomPosition(j)
            s += (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
        best = min(best, (s / n) ** 0.5)
    if not math.isfinite(best):
        raise ValueError("no RMSD match")
    return float(best)


def carboxylate_oxygens(mol: Chem.Mol) -> list[list[float]]:
    """Coords of oxygens in carboxylate / COOH."""
    pats = [
        Chem.MolFromSmarts("[CX3](=O)[O-]"),
        Chem.MolFromSmarts("[CX3](=O)[OH]"),
    ]
    conf = mol.GetConformer()
    coords: list[list[float]] = []
    seen: set[int] = set()
    for pat in pats:
        if pat is None:
            continue
        for match in mol.GetSubstructMatches(pat):
            # match: C, =O, O-
            for idx in match[1:]:
                if idx in seen:
                    continue
                seen.add(idx)
                p = conf.GetAtomPosition(idx)
                coords.append([p.x, p.y, p.z])
    return coords


def min_acid_arg_dist(mol: Chem.Mol, arg_atoms: dict) -> float | None:
    oxy = carboxylate_oxygens(mol)
    if not oxy:
        return None
    nitrogens = [arg_atoms[k] for k in ("NE", "NH1", "NH2") if k in arg_atoms]
    best = float("inf")
    for o in oxy:
        for n in nitrogens:
            d = math.sqrt(sum((o[i] - n[i]) ** 2 for i in range(3)))
            best = min(best, d)
    return None if not math.isfinite(best) else float(best)


def parse_sdf_readouts(sdf: Path) -> dict:
    poses = load_poses(sdf)
    if not poses:
        return {"n_poses": 0, "error": "no_poses"}
    rows = []
    for i, m in enumerate(poses):
        cnnscore = _fprop(m, "CNNscore")
        cnnaff = _fprop(m, "CNNaffinity")
        cnn_vs = _fprop(m, "CNN_VS")
        if cnn_vs is None and cnnscore is not None and cnnaff is not None:
            cnn_vs = cnnscore * cnnaff
        rows.append(
            {
                "mode": i + 1,
                "CNNscore": cnnscore,
                "CNNaffinity": cnnaff,
                "CNN_VS": cnn_vs,
                "minimizedAffinity": _fprop(m, "minimizedAffinity"),
            }
        )
    df = pd.DataFrame(rows)
    # select by CNNscore
    valid = df.dropna(subset=["CNNscore"])
    if valid.empty:
        return {"n_poses": len(poses), "error": "no_CNNscore", "modes": rows}
    i_star = int(valid["CNNscore"].idxmax())
    i_maxaff = int(df.dropna(subset=["CNNaffinity"])["CNNaffinity"].idxmax())
    i_p0 = i_star
    return {
        "n_poses": len(poses),
        "C1_P2star": float(df.loc[i_star, "CNNaffinity"]) if pd.notna(df.loc[i_star, "CNNaffinity"]) else None,
        "C1_VS": float(df.loc[i_star, "CNN_VS"]) if pd.notna(df.loc[i_star, "CNN_VS"]) else None,
        "C1_P2max": float(df.loc[i_maxaff, "CNNaffinity"]),
        "C1_P0": float(df.loc[i_p0, "CNNscore"]),
        "selected_mode_cnnscore": int(df.loc[i_star, "mode"]),
        "selected_mode_p2max": int(df.loc[i_maxaff, "mode"]),
        "modes": rows,
        "_pose_index_star": i_star,
        "_poses": poses,
    }


def evaluate_selfdock(
    sdf: Path,
    ref_sdf: Path | None,
    arg477_json: Path | None,
    ligand_id: str,
    seed: int,
    target: str,
) -> dict:
    parsed = parse_sdf_readouts(sdf)
    out = {
        "ligand_id": ligand_id,
        "seed": seed,
        "target": target,
        "sdf": str(sdf),
        "n_poses": parsed.get("n_poses"),
        "C1_P2star": parsed.get("C1_P2star"),
        "C1_VS": parsed.get("C1_VS"),
        "C1_P2max": parsed.get("C1_P2max"),
        "C1_P0": parsed.get("C1_P0"),
        "selected_mode_cnnscore": parsed.get("selected_mode_cnnscore"),
    }
    if parsed.get("error"):
        out["error"] = parsed["error"]
        out["pass"] = False
        return out

    poses = parsed["_poses"]
    i_star = parsed["_pose_index_star"]
    star_pose = poses[i_star]

    if ref_sdf and ref_sdf.exists():
        ref = Chem.SDMolSupplier(str(ref_sdf), removeHs=False)[0]
        if ref is None:
            out["rmsd_error"] = "ref_unreadable"
        else:
            try:
                rmsd_star = pose_rmsd(star_pose, ref)
                rmsds = []
                for p in poses:
                    try:
                        rmsds.append(pose_rmsd(p, ref))
                    except Exception:
                        rmsds.append(float("nan"))
                out["rmsd_cnnscore_selected"] = rmsd_star
                out["rmsd_top1_by_order"] = rmsds[0] if rmsds else None
                out["rmsd_best_of_n"] = float(min(x for x in rmsds if math.isfinite(x))) if any(math.isfinite(x) for x in rmsds) else None
            except Exception as exc:  # noqa: BLE001
                out["rmsd_error"] = str(exc)

    if arg477_json and arg477_json.exists():
        arg = json.loads(arg477_json.read_text())
        d = min_acid_arg_dist(star_pose, arg["atoms"])
        out["acid_arg477_min_A"] = d
        # also best-of-N
        best_d = None
        for p in poses:
            dd = min_acid_arg_dist(p, arg["atoms"])
            if dd is not None:
                best_d = dd if best_d is None else min(best_d, dd)
        out["acid_arg477_best_of_n_A"] = best_d

    # Pass criteria from campaign_c1.yaml
    if target == "urat1_9dkb" and ligand_id.lower() == "lesinurad":
        rmsd_ok = out.get("rmsd_cnnscore_selected") is not None and out["rmsd_cnnscore_selected"] <= 2.0
        acid_ok = out.get("acid_arg477_min_A") is not None and out["acid_arg477_min_A"] <= 4.0
        out["pass_rmsd"] = bool(rmsd_ok)
        out["pass_acid_arg477"] = bool(acid_ok)
        out["pass"] = bool(rmsd_ok and acid_ok)
    elif target == "nlrp3_7alv" and ligand_id.upper().startswith(("NP3", "RM5")):
        rmsd_ok = out.get("rmsd_cnnscore_selected") is not None and out["rmsd_cnnscore_selected"] <= 2.0
        out["pass_rmsd"] = bool(rmsd_ok)
        out["pass"] = bool(rmsd_ok)
    else:
        out["pass"] = None  # not a gate molecule
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdf", type=Path, required=True)
    ap.add_argument("--ref-sdf", type=Path, default=None)
    ap.add_argument("--arg477-json", type=Path, default=None)
    ap.add_argument("--ligand-id", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out-json", type=Path, required=True)
    args = ap.parse_args()

    result = evaluate_selfdock(
        args.sdf, args.ref_sdf, args.arg477_json, args.ligand_id, args.seed, args.target
    )
    # drop non-serializable
    result.pop("_poses", None)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
