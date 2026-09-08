"""C1 Acid-track pose selection (A1 exploratory vs A2 geometry-first).

A2 (Amendment accepted 2026-08-27):
  9 poses -> geometry-compatible subset -> highest CNNscore within subset.
A1 (frozen exploratory):
  highest CNNscore pose -> geometry check.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Literal

from parse_c1_sdf_readouts import (
    _fprop,
    carboxylate_oxygens,
    heavy_centroid,
    load_poses,
    load_ref_centroid,
    min_acid_arg_dist,
)

ARG_THRESH_A = 7.7027
CENTROID_MAX_A = 6.0
NLRP3_CNN_MIN = 0.5

PoseRule = Literal["a1", "a2"]


def com_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return float(math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3))))


def urat1_pose_geometry(
    pose,
    arg_atoms: dict,
    ref_com: tuple[float, float, float],
    arg_thresh: float = ARG_THRESH_A,
    centroid_max: float = CENTROID_MAX_A,
) -> dict:
    d_arg = min_acid_arg_dist(pose, arg_atoms)
    has_acid = bool(carboxylate_oxygens(pose))
    d_com = com_distance(heavy_centroid(pose), ref_com)
    in_pocket = d_com <= centroid_max
    geom_pass = bool(
        has_acid
        and d_arg is not None
        and d_arg <= arg_thresh
        and in_pocket
    )
    return {
        "acid_arg477_min_A": d_arg,
        "has_carboxylate_oxy_in_pose": has_acid,
        "centroid_to_crystal_lig_A": d_com,
        "pass_pocket_centroid": in_pocket,
        "pass_arg_A1": (d_arg is not None and d_arg <= arg_thresh),
        "geometry_pass": geom_pass,
        "CNNscore": _fprop(pose, "CNNscore"),
        "CNNaffinity": _fprop(pose, "CNNaffinity"),
    }


def select_urat1_pose_index(
    poses: list,
    arg_atoms: dict,
    ref_com: tuple[float, float, float],
    rule: PoseRule = "a2",
    arg_thresh: float = ARG_THRESH_A,
    centroid_max: float = CENTROID_MAX_A,
) -> tuple[int | None, str, list[int]]:
    """Return (selected_index, selection_status, geometry_passing_indices)."""
    if not poses:
        return None, "no_poses", []
    geom_indices: list[int] = []
    for j, pose in enumerate(poses):
        g = urat1_pose_geometry(pose, arg_atoms, ref_com, arg_thresh, centroid_max)
        if g["geometry_pass"]:
            geom_indices.append(j)
    if rule == "a1":
        i_star = max(range(len(poses)), key=lambda j: _fprop(poses[j], "CNNscore") or -1.0)
        status = "a1_cnnscore_top1"
        return i_star, status, geom_indices
    if not geom_indices:
        return None, "a2_no_geometry_pass", []
    i_star = max(geom_indices, key=lambda j: _fprop(poses[j], "CNNscore") or -1.0)
    return i_star, "a2_geometry_then_cnnscore", geom_indices


def evaluate_urat1_acid_sdf(
    sdf: Path,
    arg_json: Path,
    ref_com: tuple[float, float, float],
    ligand_id: str,
    seed: int,
    rule: PoseRule = "a2",
) -> dict:
    poses = load_poses(sdf)
    arg = json.loads(arg_json.read_text())
    i_star, status, geom_indices = select_urat1_pose_index(
        poses, arg["atoms"], ref_com, rule=rule
    )
    row: dict = {
        "ligand_id": ligand_id,
        "target": "urat1_9dkb",
        "seed": seed,
        "pose_selection_rule": rule,
        "pose_selection_status": status,
        "n_poses": len(poses),
        "n_geometry_pass_modes": len(geom_indices),
        "sdf": str(sdf),
    }
    if i_star is None:
        row.update(
            {
                "selected_mode": None,
                "CNNscore": None,
                "CNNaffinity": None,
                "acid_arg477_min_A": None,
                "has_carboxylate_oxy_in_pose": False,
                "pass_arg_A1": False,
                "centroid_to_crystal_lig_A": None,
                "pass_pocket_centroid": False,
                "keep_urat1_acid": False,
                "error": status,
            }
        )
        return row
    pose = poses[i_star]
    g = urat1_pose_geometry(pose, arg["atoms"], ref_com)
    row.update(
        {
            "selected_mode": i_star + 1,
            "CNNscore": g["CNNscore"],
            "CNNaffinity": g["CNNaffinity"],
            "C1_P2star": g["CNNaffinity"],
            "C1_VS": (g["CNNscore"] or 0) * (g["CNNaffinity"] or 0),
            "acid_arg477_min_A": g["acid_arg477_min_A"],
            "has_carboxylate_oxy_in_pose": g["has_carboxylate_oxy_in_pose"],
            "pass_arg_A1": g["pass_arg_A1"],
            "centroid_to_crystal_lig_A": g["centroid_to_crystal_lig_A"],
            "pass_pocket_centroid": g["pass_pocket_centroid"],
            "keep_urat1_acid": bool(g["geometry_pass"]) if rule == "a2" else bool(
                g["pass_arg_A1"] and g["has_carboxylate_oxy_in_pose"] and g["pass_pocket_centroid"]
            ),
        }
    )
    if rule == "a1":
        row["keep_urat1_acid"] = bool(
            row["pass_arg_A1"] and row["has_carboxylate_oxy_in_pose"] and row["pass_pocket_centroid"]
        )
    return row


def evaluate_nlrp3_pose_sdf(
    sdf: Path,
    ref_com: tuple[float, float, float],
    ligand_id: str,
    seed: int,
) -> dict:
    poses = load_poses(sdf)
    if not poses:
        return {
            "ligand_id": ligand_id,
            "target": "nlrp3_7alv",
            "seed": seed,
            "error": "no_poses",
            "keep_nlrp3_pose": False,
        }
    i_star = max(range(len(poses)), key=lambda j: _fprop(poses[j], "CNNscore") or -1.0)
    pose = poses[i_star]
    d_com = com_distance(heavy_centroid(pose), ref_com)
    in_pocket = d_com <= CENTROID_MAX_A
    cnn = _fprop(pose, "CNNscore")
    row = {
        "ligand_id": ligand_id,
        "target": "nlrp3_7alv",
        "seed": seed,
        "pose_selection_rule": "a1_nlrp3_unchanged",
        "n_poses": len(poses),
        "selected_mode": i_star + 1,
        "CNNscore": cnn,
        "CNNaffinity": _fprop(pose, "CNNaffinity"),
        "centroid_to_crystal_lig_A": d_com,
        "pass_pocket_centroid": in_pocket,
        "keep_nlrp3_pose": bool(in_pocket and (cnn or 0) >= NLRP3_CNN_MIN),
        "sdf": str(sdf),
    }
    return row
