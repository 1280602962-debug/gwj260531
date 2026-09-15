#!/usr/bin/env python3
"""Phase 0 freeze + Phase 1 corrected heavy-atom boxes.

Does not overwrite official result CSVs, figures, or manuscript files.
Writes only Dual_Target_Docking/remediation_outputs/.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs"
P0 = OUT / "phase0_freeze"
P1 = OUT / "phase1_boxes"
LOCAL120 = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0"
)
EXPAND = 5.0
MIN_EDGE = 20.0

CORE_PRE_FIX = [
    # pair, metric, value, ci_lo, ci_hi, note
    ("EGFR/HER2", "AUROC_D_vs_A_pocketB", 0.6663533834586466, 0.5237584218584672, 0.793095238095238, "n_pos=28 n_neg=38"),
    ("EGFR/HER2", "AUROC_D_vs_B_pocketA", 0.4296875, 0.2817924317062248, 0.579011972770499, "historical ~0.430; n_pos=28 n_neg=32"),
    ("EGFR/HER2", "summary_min", 0.4296875, 0.2817924317062248, 0.5775046404364821, "min of two directional arms"),
    ("EGFR/HER2", "fixed_score_delta_pocketA", 0.3783482142857143, 0.20161830357142851, 0.5524553571428571, "historical Δ=0.378; D-vs-N minus D-vs-B using S_A"),
    ("EGFR/HER2", "AUROC_mean_D_vs_neither", 0.7559523809523809, 0.5625, 0.9197, "two-pocket mean; n_neither=12"),
    ("EGFR/HER2", "delta_summary_min_matched_minus_mismatched", 0.16964285714285715, 0.059145422149122794, 0.2877330714894653, "CI excludes 0"),
    ("EGFR/HER2", "Top10_D_A_B_N", "1/5/5/0", "", "", "k=11 of n=110"),
    ("EGFR/HER2", "EF_dual_10", 0.3571428571428572, "", "", "1 dual / 11 top vs 28/110"),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "UNKNOWN"


def box_from_xyz(xyz) -> dict:
    xyz = np.asarray(xyz, dtype=float)
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    span = hi - lo
    size = np.maximum(span + 2.0 * EXPAND, MIN_EDGE)
    # Symmetric extension of any short edge keeps the heavy-atom midpoint.
    center = (lo + hi) / 2.0
    return {
        "n": int(len(xyz)),
        "center_x": float(center[0]),
        "center_y": float(center[1]),
        "center_z": float(center[2]),
        "size_x": float(size[0]),
        "size_y": float(size[1]),
        "size_z": float(size[2]),
        "span_x": float(span[0]),
        "span_y": float(span[1]),
        "span_z": float(span[2]),
    }


def round_box(box: dict, nd: int = 3) -> dict:
    out = dict(box)
    for k in (
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "span_x",
        "span_y",
        "span_z",
    ):
        if k in out:
            out[k] = round(float(out[k]), nd)
    return out


def sdf_xyz(path: Path, heavy: bool):
    mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
    xyz = []
    n_h = 0
    n_heavy = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            n_h += 1
            if heavy:
                continue
        else:
            n_heavy += 1
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz.append([pos.x, pos.y, pos.z])
    return xyz, n_heavy, n_h, mol.GetNumAtoms()


def pdb_xyz(path: Path, heavy: bool):
    xyz = []
    n_h = 0
    n_heavy = 0
    for line in path.read_text().splitlines():
        if not (line.startswith("HETATM") or line.startswith("ATOM")):
            continue
        el = line[76:78].strip() if len(line) >= 78 else ""
        if not el:
            name = line[12:16].strip()
            el = "".join(c for c in name if c.isalpha())[:1]
        is_h = el.upper() == "H"
        if is_h:
            n_h += 1
            if heavy:
                continue
        else:
            n_heavy += 1
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return xyz, n_heavy, n_h, n_heavy + n_h


def load_json_box(path: Path) -> dict:
    return json.loads(path.read_text())


def freeze_phase0() -> str:
    P0.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    commit = git_commit()
    (P0 / "git_commit.txt").write_text(
        f"commit={commit}\n"
        f"generated={datetime.now(timezone.utc).isoformat()}\n"
        f"worktree={ROOT}\n",
        encoding="utf-8",
    )

    files = []
    candidates = [
        ROOT / "data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv",
        ROOT / "data/jcim_novelty_v0/tables/primary_directional_intervals_review_v1.csv",
        ROOT / "audit_outputs/audit_master_metrics.csv",
        ROOT / "audit_outputs/canonical_ligand_table.csv",
        ROOT / "audit_outputs/candidate_ranking_audit.csv",
        ROOT / "docs/MANUSCRIPT_JCIM_EN.md",
        ROOT / "docs/MANUSCRIPT_JCIM_ZH.md",
        ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md",
        ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_ZH.md",
        ROOT / "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
        ROOT / "data/egfr_her2_panel40_v0/boxes/3POZ_box.json",
        ROOT / "data/egfr_her2_panel40_v0/boxes/3RCD_box.json",
        ROOT / "data/egfr_her2_panel40_v0/boxes/3POZ_box.txt",
        ROOT / "data/egfr_her2_panel40_v0/boxes/3RCD_box.txt",
        ROOT / "data/egfr_her2_panel40_v0/tables/panel_v0_40.csv",
        ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        ROOT / "data/egfr_her2_panel120_v0/tables/scores_vina.csv",
        ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
        ROOT / "data/ache_bche_panel_v0/scripts/build_strict_panels.py",
        ROOT / "figures/jcim_article/scripts/audit_figures_pr32.py",
        ROOT / "figures/jcim_article/scripts/plot_jcim_article_figures_v3.py",
        ROOT / "figures/jcim_article/plotted_values.json",
    ]
    for folder, patterns in (
        (ROOT / "figures/jcim_article", ("Fig1*", "Fig2*", "Fig3*", "Fig4*", "Fig5*", "Fig6*", "FigS1*", "FigS2*", "FigS3*", "FigS4*")),
        (ROOT / "submission_pack/tables", ("Table*", "table*")),
        (ROOT / "data/jcim_novelty_v0/tables", ("*table*", "*Table*")),
    ):
        if folder.exists():
            for pat in patterns:
                candidates.extend(sorted(folder.glob(pat)))

    seen = set()
    rows = []
    for path in candidates:
        if not path.exists() or not path.is_file():
            continue
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        rel = str(path.relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)
        rows.append(
            {
                "path": rel,
                "present": 1,
                "nbytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "role": "pre_fix_official",
            }
        )
        files.append(rel)

    man_path = OUT / "pre_fix_manifest.csv"
    with man_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "present", "nbytes", "sha256", "role"])
        w.writeheader()
        w.writerows(rows)

    core_path = P0 / "pre_fix_core_eight_metrics.csv"
    with core_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pair", "metric", "pre_fix", "ci_lo", "ci_hi", "note", "git_commit"])
        for rec in CORE_PRE_FIX:
            w.writerow([*rec, commit])

    # post-fix placeholders until docking finishes
    post = OUT / "post_fix_manifest.csv"
    with post.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["path", "present", "nbytes", "sha256", "role", "status"]
        )
        w.writeheader()
        w.writerow(
            {
                "path": "PENDING_PHASE1_DOCKING",
                "present": 0,
                "nbytes": 0,
                "sha256": "",
                "role": "post_fix",
                "status": "not_yet_generated",
            }
        )

    for name in (
        "metric_diff_pre_vs_post.csv",
        "figure_diff_pre_vs_post.csv",
        "sample_diff_pre_vs_post.csv",
        "provenance_post_fix.csv",
    ):
        p = OUT / name
        if not p.exists():
            with p.open("w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["status", "note"])
                w.writerow(["PENDING", "filled after Phase 1/2 recomputation"])

    print(f"Phase 0 freeze: commit={commit} n_files={len(rows)} -> {man_path}")
    return commit


def compare_row(pdb: str, old: dict, new: dict, n_heavy: int, n_h: int, n_all: int) -> dict:
    def g(d, k):
        return float(d[k])

    rec = {
        "pdb": pdb,
        "old_center_x": g(old, "center_x"),
        "old_center_y": g(old, "center_y"),
        "old_center_z": g(old, "center_z"),
        "old_size_x": g(old, "size_x"),
        "old_size_y": g(old, "size_y"),
        "old_size_z": g(old, "size_z"),
        "new_center_x": round(new["center_x"], 3),
        "new_center_y": round(new["center_y"], 3),
        "new_center_z": round(new["center_z"], 3),
        "new_size_x": round(new["size_x"], 3),
        "new_size_y": round(new["size_y"], 3),
        "new_size_z": round(new["size_z"], 3),
        "abs_diff_center_x": abs(g(old, "center_x") - new["center_x"]),
        "abs_diff_center_y": abs(g(old, "center_y") - new["center_y"]),
        "abs_diff_center_z": abs(g(old, "center_z") - new["center_z"]),
        "abs_diff_size_x": abs(g(old, "size_x") - new["size_x"]),
        "abs_diff_size_y": abs(g(old, "size_y") - new["size_y"]),
        "abs_diff_size_z": abs(g(old, "size_z") - new["size_z"]),
        "heavy_atom_count": n_heavy,
        "hydrogen_count": n_h,
        "n_ligand_atoms_all": n_all,
        "old_n_ligand_atoms_field": old.get("n_ligand_atoms", old.get("n_heavy_atoms", "")),
    }
    rec["max_abs_diff_any_edge_or_center"] = max(
        rec[k]
        for k in rec
        if k.startswith("abs_diff_")
    )
    return rec


def phase1_boxes() -> dict:
    P1.mkdir(parents=True, exist_ok=True)
    specs = {
        "3POZ": (
            ROOT / "data/egfr_her2_panel40_v0/cognate_qc/3POZ_03P_crystal.sdf",
            ROOT / "data/egfr_her2_panel40_v0/boxes/3POZ_box.json",
            "sdf",
        ),
        "3RCD": (
            ROOT / "data/egfr_her2_panel40_v0/cognate_qc/3RCD_03P_crystal.sdf",
            ROOT / "data/egfr_her2_panel40_v0/boxes/3RCD_box.json",
            "sdf",
        ),
    }
    out_rows = {}
    for pdb, (src, old_path, kind) in specs.items():
        xyz_h, n_heavy, n_h, n_all = sdf_xyz(src, heavy=True)
        xyz_all, _, _, _ = sdf_xyz(src, heavy=False)
        new_h = box_from_xyz(xyz_h)
        new_all = box_from_xyz(xyz_all)
        old = load_json_box(old_path)
        rec = compare_row(pdb, old, new_h, n_heavy, n_h, n_all)
        rec["all_atom_center_x"] = round(new_all["center_x"], 3)
        rec["all_atom_center_y"] = round(new_all["center_y"], 3)
        rec["all_atom_center_z"] = round(new_all["center_z"], 3)
        rec["all_atom_size_x"] = round(new_all["size_x"], 3)
        rec["all_atom_size_y"] = round(new_all["size_y"], 3)
        rec["all_atom_size_z"] = round(new_all["size_z"], 3)
        rec["source_cognate"] = str(src.relative_to(ROOT))
        rec["rule"] = "heavy-atom AABB + 5 Å/side, min edge 20 Å, center = heavy midpoint"
        out_rows[pdb] = rec
        rb = round_box(new_h)
        payload = {
            "center_x": rb["center_x"],
            "center_y": rb["center_y"],
            "center_z": rb["center_z"],
            "size_x": rb["size_x"],
            "size_y": rb["size_y"],
            "size_z": rb["size_z"],
            "n_heavy_atoms": n_heavy,
            "n_hydrogen_atoms": n_h,
            "n_ligand_atoms": n_heavy,
            "pdb": pdb,
            "ligand": "03P",
            "construction": "cognate_heavy_atom_AABB_plus_5A_min20",
            "status": "corrected_not_yet_canonical",
        }
        (P1 / f"{pdb}_box_corrected.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        (P1 / f"{pdb}_box_corrected.txt").write_text(
            "\n".join(
                [
                    f"center_x = {payload['center_x']}",
                    f"center_y = {payload['center_y']}",
                    f"center_z = {payload['center_z']}",
                    f"size_x = {payload['size_x']}",
                    f"size_y = {payload['size_y']}",
                    f"size_z = {payload['size_z']}",
                    f"n_heavy_atoms = {n_heavy}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    def write_cmp(name: str, rec: dict):
        path = OUT / name
        fields = list(rec.keys())
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerow(rec)
        print(f"wrote {path}")

    write_cmp("egfr_3POZ_box_old_vs_corrected.csv", out_rows["3POZ"])
    write_cmp("her2_3RCD_box_old_vs_corrected.csv", out_rows["3RCD"])

    # Consistency vs AChE / PIK3CA deposited boxes
    checks = [
        (
            "4EY7",
            ROOT / "data/ache_bche_panel_v0/cognate_qc/4EY7_E20_crystal.sdf",
            ROOT / "data/ache_bche_panel_v0/boxes/4EY7_box.json",
            "sdf",
        ),
        (
            "4BDS",
            ROOT / "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_crystal.sdf",
            ROOT / "data/ache_bche_panel_v0/boxes/4BDS_box.json",
            "sdf",
        ),
        (
            "4L23",
            ROOT / "data/pik3ca_mtor_panel48_v0/tables/4L23_cocrystal_X6K.pdb",
            ROOT / "data/pik3ca_mtor_panel48_v0/boxes/4L23_box.json",
            "pdb",
        ),
        (
            "4JT6",
            ROOT / "data/pik3ca_mtor_panel48_v0/tables/4JT6_cocrystal_X6K.pdb",
            ROOT / "data/pik3ca_mtor_panel48_v0/boxes/4JT6_box.json",
            "pdb",
        ),
    ]
    cons = []
    for pdb, src, old_path, kind in checks:
        if not src.exists() or not old_path.exists():
            cons.append(
                {
                    "pdb": pdb,
                    "status": "MISSING_INPUT",
                    "source": str(src),
                    "old_box": str(old_path),
                }
            )
            continue
        if kind == "sdf":
            xyz_h, n_heavy, n_h, n_all = sdf_xyz(src, heavy=True)
        else:
            xyz_h, n_heavy, n_h, n_all = pdb_xyz(src, heavy=True)
        new_h = round_box(box_from_xyz(xyz_h))
        old = load_json_box(old_path)
        diffs = [
            abs(float(old[k]) - new_h[k])
            for k in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
        ]
        cons.append(
            {
                "pdb": pdb,
                "status": "MATCH" if max(diffs) < 0.0015 else "MISMATCH",
                "max_abs_diff": max(diffs),
                "n_heavy": n_heavy,
                "n_hydrogen": n_h,
                "old_n_heavy_field": old.get("n_heavy_atoms", old.get("n_ligand_atoms", "")),
                "recomputed_center_xyz": f"{new_h['center_x']},{new_h['center_y']},{new_h['center_z']}",
                "deposited_center_xyz": f"{old['center_x']},{old['center_y']},{old['center_z']}",
                "recomputed_size_xyz": f"{new_h['size_x']},{new_h['size_y']},{new_h['size_z']}",
                "deposited_size_xyz": f"{old['size_x']},{old['size_y']},{old['size_z']}",
                "source": str(src.relative_to(ROOT)),
            }
        )
    cons_path = P1 / "box_rule_consistency_ache_pik3ca.csv"
    fields = sorted({k for r in cons for k in r})
    with cons_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(cons)
    print(f"wrote {cons_path}")
    return out_rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    commit = freeze_phase0()
    boxes = phase1_boxes()
    print("3POZ max|Δ|", boxes["3POZ"]["max_abs_diff_any_edge_or_center"])
    print("3RCD max|Δ|", boxes["3RCD"]["max_abs_diff_any_edge_or_center"])
    print("commit", commit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
