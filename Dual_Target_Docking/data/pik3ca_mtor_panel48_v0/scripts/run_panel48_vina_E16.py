#!/usr/bin/env python3
"""Full panel 48×2 @ E=16 seed=20260727. Only after COGNATE_QC_VERDICT_E16 Go."""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0")
VINA = Path("/home/gwj/miniconda3/bin/vina")
SEED = 20260727
EXHAUST = 16
N_MODES = 9
ENERGY_RANGE = 3
TARGETS = ["4L23", "4JT6"]


def assert_go() -> None:
    verdict = (ROOT / "analysis" / "cognate_redock_v0" / "COGNATE_QC_VERDICT_E16.md").read_text()
    if "Verdict: Go" not in verdict:
        raise SystemExit("Gate blocked: COGNATE_QC_VERDICT_E16.md is not Go")
    proto = (ROOT / "protocol" / "protocol.yaml").read_text()
    if "exhaustiveness: 16" not in proto and "exhaustiveness_v0_1: 16" not in proto:
        raise SystemExit("Gate blocked: protocol.yaml not frozen at E=16")
    n_lig = len(list((ROOT / "ligands_pdbqt").glob("PM48_*.pdbqt")))
    if n_lig != 48:
        raise SystemExit(f"Gate blocked: expected 48 ligands, found {n_lig}")


def split_models(all_pose: Path, pose_dir: Path) -> int:
    text = all_pose.read_text(errors="ignore")
    blocks, cur, inn = [], [], False
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            if cur:
                blocks.append(cur)
            cur = [line]
            inn = True
        elif line.startswith("ENDMDL"):
            cur.append(line)
            blocks.append(cur)
            cur, inn = [], False
        elif inn:
            cur.append(line)
    if cur:
        blocks.append(cur)
    if not blocks and text.strip():
        blocks = [[text]]
    pose_dir.mkdir(parents=True, exist_ok=True)
    for i, b in enumerate(blocks, 1):
        (pose_dir / f"mode_{i:02d}.pdbqt").write_text("".join(b))
    return len(blocks)


def parse_mode1_score(path: Path) -> float | None:
    if not path.exists():
        return None
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def run_one(target: str, lig: str) -> dict:
    boxes = json.loads((ROOT / "boxes" / "all_boxes.json").read_text())
    box = boxes[target]
    pose_dir = ROOT / "poses" / target / lig
    pose_dir.mkdir(parents=True, exist_ok=True)
    all_pose = pose_dir / f"{lig}_all_modes.pdbqt"
    log = ROOT / "logs" / "vina" / f"{target}_{lig}.log"
    conf = ROOT / "logs" / "vina_confs" / f"{target}_{lig}.txt"
    conf.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    conf.write_text(
        "\n".join(
            [
                f"receptor = {ROOT / 'receptors' / f'{target}_receptor.pdbqt'}",
                f"ligand = {ROOT / 'ligands_pdbqt' / f'{lig}.pdbqt'}",
                f"out = {all_pose}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                f"energy_range = {ENERGY_RANGE}",
                "cpu = 1",
                f"seed = {SEED}",
            ]
        )
        + "\n"
    )
    with log.open("w") as fh:
        proc = subprocess.run([str(VINA), "--config", str(conf)], stdout=fh, stderr=subprocess.STDOUT)
    status = "success"
    reason = ""
    n_modes = 0
    if proc.returncode != 0 or not all_pose.exists():
        status = "fail"
        reason = f"vina_rc={proc.returncode}"
    else:
        n_modes = split_models(all_pose, pose_dir)
        log_txt = log.read_text(errors="ignore")
        m = re.search(r"random seed:\s*(-?\d+)", log_txt)
        if not m or int(m.group(1)) != SEED:
            status = "fail"
            reason = f"seed_mismatch:{m.group(1) if m else None}"
        elif n_modes != 9:
            status = "fail"
            reason = f"n_modes={n_modes}"
    return {
        "target": target,
        "ligand_id": lig,
        "seed": SEED,
        "exhaustiveness": EXHAUST,
        "status": status,
        "reason": reason,
        "n_modes_written": n_modes,
        "log_path": str(log),
        "pose_dir": str(pose_dir),
        "vina_mode1": parse_mode1_score(pose_dir / "mode_01.pdbqt"),
    }


def write_score_tables(job_rows: list[dict]) -> None:
    panel = {r["panel_id"]: r for r in csv.DictReader((ROOT / "tables" / "panel_v0_48.csv").open())}
    long_rows = []
    for target in TARGETS:
        for i in range(1, 49):
            lig = f"PM48_{i:02d}"
            for mode in range(1, 10):
                pose = ROOT / "poses" / target / lig / f"mode_{mode:02d}.pdbqt"
                long_rows.append(
                    {
                        "panel_id": lig,
                        "target": target,
                        "mode": mode,
                        "vina_score": parse_mode1_score(pose),
                        "seed": SEED,
                        "exhaustiveness": EXHAUST,
                    }
                )
    with (ROOT / "tables" / "scores_vina_long.csv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["panel_id", "target", "mode", "vina_score", "seed", "exhaustiveness"]
        )
        w.writeheader()
        w.writerows(long_rows)

    summary = []
    for i in range(1, 49):
        lig = f"PM48_{i:02d}"
        s23 = parse_mode1_score(ROOT / "poses" / "4L23" / lig / "mode_01.pdbqt")
        s46 = parse_mode1_score(ROOT / "poses" / "4JT6" / lig / "mode_01.pdbqt")
        if s23 is None or s46 is None:
            continue
        mean = (s23 + s46) / 2.0
        vmin = max(s23, s46)  # less negative = worse for vina
        # keep native: vina_min as the worse (higher) affinity
        summary.append(
            {
                "panel_id": lig,
                "class": panel.get(lig, {}).get("class", ""),
                "vina_4L23_mode1": s23,
                "vina_4JT6_mode1": s46,
                "vina_mean": mean,
                "vina_min": vmin,
                "vina_delta": s23 - s46,
                "best_mode_4L23": 1,
                "best_mode_4JT6": 1,
            }
        )
    with (ROOT / "tables" / "scores_vina.csv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "panel_id",
                "class",
                "vina_4L23_mode1",
                "vina_4JT6_mode1",
                "vina_mean",
                "vina_min",
                "vina_delta",
                "best_mode_4L23",
                "best_mode_4JT6",
            ],
        )
        w.writeheader()
        w.writerows(summary)


def write_manifest() -> None:
    path = ROOT / "MANIFEST.md"
    text = f"""# MANIFEST — pik3ca_mtor_panel48_v0

- exhaustiveness: **16** (`exhaustiveness_v0_1`)
- seed_fixed_global: 20260727
- cognate QC E=8: No-Go (archived) — `analysis/cognate_redock_v0/COGNATE_QC_VERDICT.md`
- cognate QC E=16: **Go** — `analysis/cognate_redock_v0/COGNATE_QC_VERDICT_E16.md`
- E=8 cognate poses archived at: `poses/cognate_E8_archive/`
- cognate E16 poses: `poses/cognate_E16/`
- full panel poses: `poses/4L23|4JT6/<PM48_XX>/`
- note: 4JT6 mode1 may fail redock RMSD while best_of_9 passes → keep 9 modes; plan RTM best-of-9
"""
    path.write_text(text)


def main() -> int:
    assert_go()
    ligands = [f"PM48_{i:02d}" for i in range(1, 49)]
    jobs = [(t, lig) for t in TARGETS for lig in ligands]
    workers = min(8, max(1, os.cpu_count() or 1))
    print(f"Launching {len(jobs)} jobs with {workers} workers", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, t, lig): (t, lig) for t, lig in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            print(
                f"[{i}/{len(jobs)}] {r['target']} {r['ligand_id']} -> {r['status']} {r['reason']}",
                flush=True,
            )
    results.sort(key=lambda r: (r["target"], r["ligand_id"]))
    with (ROOT / "tables" / "job_status.csv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "target",
                "ligand_id",
                "seed",
                "exhaustiveness",
                "status",
                "reason",
                "n_modes_written",
                "log_path",
                "pose_dir",
                "vina_mode1",
            ],
        )
        w.writeheader()
        w.writerows(results)
    n_ok = sum(1 for r in results if r["status"] == "success")
    write_score_tables(results)
    write_manifest()
    print(json.dumps({"n_jobs": len(results), "n_success": n_ok, "n_fail": len(results) - n_ok}), flush=True)
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
