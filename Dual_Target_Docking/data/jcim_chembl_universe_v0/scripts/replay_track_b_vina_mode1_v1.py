#!/usr/bin/env python3
"""Replay production Vina mode-1 affinities from saved pose PDBQT files.

Does not redock. Looks for mode_01.pdbqt under --pose-root (default:
local_track_b_v0/poses). Writes a per-job mapping even when the pose is
missing, so absence is explicit rather than a silent PASS.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
SCORE = LOCAL / "tables" / "scores_vina_mode1_v1.csv"
STATUS = LOCAL / "tables" / "job_status.csv"


def parse_mode1(path: Path):
    in_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            in_model = True
            continue
        if in_model and "VINA RESULT" in line:
            m = re.search(r"VINA RESULT:\s+(-?\d+\.?\d*)", line)
            return float(m.group(1)) if m else None
        if line.startswith("ENDMDL") and in_model:
            break
    return None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pose-root", type=Path, default=LOCAL / "poses")
    args = ap.parse_args()
    pose_root = args.pose_root
    scores = {(r["pair"], r["target"], r["ligand"]): r for r in csv.DictReader(SCORE.open())}
    jobs = list(csv.DictReader(STATUS.open()))
    rows = []
    n_match = n_miss = n_mis = n_skip = 0
    for job in jobs:
        key = (job["pair"], job["target"], job["ligand"])
        pose = pose_root / job["target"] / job["ligand"] / "mode_01.pdbqt"
        rec = {
            "pair": job["pair"],
            "target": job["target"],
            "ligand": job["ligand"],
            "job_status": job["status"],
            "csv_mode1_energy": scores.get(key, {}).get("mode1_energy", ""),
            "pose_path": str(pose) if pose.exists() else "",
            "pose_mode1_energy": "",
            "sha256_mode01": "",
            "replay": "",
        }
        if job["status"] != "success":
            rec["replay"] = "not_applicable_skipped"
            n_skip += 1
        elif not pose.exists():
            rec["replay"] = "pose_missing"
            n_miss += 1
        else:
            e = parse_mode1(pose)
            rec["pose_mode1_energy"] = "" if e is None else f"{e:.6g}"
            rec["sha256_mode01"] = sha256(pose)
            csv_e = scores.get(key, {}).get("mode1_energy")
            if e is None or csv_e in (None, ""):
                rec["replay"] = "parse_fail"
                n_mis += 1
            elif abs(float(e) - float(csv_e)) <= 1e-6:
                rec["replay"] = "match"
                n_match += 1
            else:
                rec["replay"] = "mismatch"
                n_mis += 1
        rows.append(rec)

    dest = LOCAL / "tables" / "scores_vina_mode1_pose_replay_v1.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"jobs={len(rows)} match={n_match} missing={n_miss} mismatch_or_parse={n_mis} skipped={n_skip}")
    print("wrote", dest)
    return 0 if n_mis == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
