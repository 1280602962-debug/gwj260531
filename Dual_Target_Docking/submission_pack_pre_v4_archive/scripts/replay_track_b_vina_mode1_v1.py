#!/usr/bin/env python3
"""Replay production Vina mode-1 affinities from saved pose PDBQT files.

Does not redock. Looks for mode_01.pdbqt under --pose-root (default:
local_track_b_v0/poses). Writes a per-job mapping even when the pose is
missing, so absence is explicit rather than a silent PASS.

Exit codes:
  0  recorded skips only, or every success job matched
  2  success-job pose missing, parse failure, score mismatch, duplicate
     keys, or score rows without a job
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
SCORE = LOCAL / "tables" / "scores_vina_mode1_v1.csv"
STATUS = LOCAL / "tables" / "job_status.csv"
DEFAULT_OUT = LOCAL / "tables" / "scores_vina_mode1_pose_replay_v1.csv"
FIELDS = [
    "pair",
    "target",
    "ligand",
    "job_status",
    "csv_mode1_energy",
    "pose_path",
    "pose_mode1_energy",
    "sha256_mode01",
    "replay",
    "pose_location",
]


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


def load_unique_csv(path: Path, key_fields: tuple[str, ...], kind: str) -> dict[tuple, dict]:
    seen: dict[tuple, dict] = {}
    dups = []
    for r in csv.DictReader(path.open(encoding="utf-8", newline="")):
        key = tuple(r[k] for k in key_fields)
        if key in seen:
            dups.append(key)
            continue
        seen[key] = r
    if dups:
        raise SystemExit(f"duplicate {kind} keys: {dups[:8]!r} (n={len(dups)})")
    return seen


def rel_pose_path(pose: Path, pose_root: Path) -> str:
    try:
        return str(Path("poses") / pose.relative_to(pose_root))
    except ValueError:
        return str(Path("poses") / pose.name)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pose-root", type=Path, default=LOCAL / "poses")
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="CSV destination. Default is the committed replay table; "
        "verification runs must pass a temporary path.",
    )
    ap.add_argument(
        "--allow-overwrite-committed",
        action="store_true",
        help="Permit writing the committed replay evidence table.",
    )
    args = ap.parse_args(argv)
    dest = args.output or DEFAULT_OUT
    if dest.resolve() == DEFAULT_OUT.resolve() and not args.allow_overwrite_committed:
        print(
            "refusing to overwrite committed replay table; "
            "pass --output PATH or --allow-overwrite-committed",
            file=sys.stderr,
        )
        return 2

    scores = load_unique_csv(SCORE, ("pair", "target", "ligand"), "score")
    jobs = load_unique_csv(STATUS, ("pair", "target", "ligand"), "job")
    orphan_scores = sorted(set(scores) - set(jobs))
    if orphan_scores:
        print(f"score rows without a job: {orphan_scores[:8]!r} n={len(orphan_scores)}", file=sys.stderr)
        return 2

    pose_root = args.pose_root
    rows = []
    n_match = n_miss = n_mis = n_skip = 0
    hard_fail = bool(orphan_scores)
    for key, job in jobs.items():
        pose = pose_root / job["target"] / job["ligand"] / "mode_01.pdbqt"
        rec = {
            "pair": job["pair"],
            "target": job["target"],
            "ligand": job["ligand"],
            "job_status": job["status"],
            "csv_mode1_energy": scores.get(key, {}).get("mode1_energy", ""),
            "pose_path": rel_pose_path(pose, pose_root) if pose.exists() else "",
            "pose_mode1_energy": "",
            "sha256_mode01": "",
            "replay": "",
            "pose_location": "",
        }
        if job["status"] != "success":
            rec["replay"] = "not_applicable_skipped"
            rec["pose_location"] = "not_applicable_skipped"
            n_skip += 1
        elif not pose.exists():
            rec["replay"] = "pose_missing"
            rec["pose_location"] = "missing"
            n_miss += 1
            hard_fail = True
        else:
            e = parse_mode1(pose)
            rec["pose_mode1_energy"] = "" if e is None else f"{e:.6g}"
            rec["sha256_mode01"] = sha256(pose)
            rec["pose_location"] = "local_disk_gitignored"
            csv_e = scores.get(key, {}).get("mode1_energy")
            if e is None or csv_e in (None, ""):
                rec["replay"] = "parse_fail"
                n_mis += 1
                hard_fail = True
            elif abs(float(e) - float(csv_e)) <= 1e-6:
                rec["replay"] = "match"
                n_match += 1
            else:
                rec["replay"] = "mismatch"
                n_mis += 1
                hard_fail = True
        rows.append(rec)

    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(
        f"jobs={len(rows)} match={n_match} missing={n_miss} "
        f"mismatch_or_parse={n_mis} skipped={n_skip}"
    )
    print("wrote", dest)
    return 2 if hard_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
