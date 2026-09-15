#!/usr/bin/env python3
"""Independent GNINA docking for EGFR/HER2 using canonical heavy-atom boxes.

Does not reuse old-box poses. Writes to remediation_outputs/phase_gnina_independent/.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
OUT = ROOT / "remediation_outputs/phase_gnina_independent"
LIGDIR = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0/ligands_pdbqt"
)
RECDIR = ROOT / "data/egfr_her2_panel40_v0/receptors"
BOXDIR = ROOT / "data/egfr_her2_panel40_v0/boxes"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
GNINA = Path("/mnt/d/CADD paper exercise/gnina/bin/gnina")
GNINA_LIB = Path("/mnt/d/CADD paper exercise/gnina/conda_env/lib")
SEED = 20260727
WORKERS = int(os.environ.get("GNINA_WORKERS", "2"))


def parse_aff(path: Path) -> float | None:
    in_model = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("MODEL 1"):
            in_model = True
        elif in_model and line.startswith("REMARK minimizedAffinity"):
            return float(line.split()[2])
        elif in_model and line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def dock(target: str, lig: str, box: dict) -> dict:
    t0 = time.time()
    rec = RECDIR / f"{target}_protein.pdb"
    if not rec.exists():
        rec = RECDIR / f"{target}_receptor.pdbqt"
    ligand = LIGDIR / f"{lig}.pdbqt"
    out_dir = OUT / "out"
    log_dir = OUT / "logs"
    pose_dir = OUT / "poses" / target / lig
    for d in (out_dir, log_dir, pose_dir):
        d.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{target}_{lig}_gnina_out.pdbqt"
    log = log_dir / f"{target}_{lig}.log"
    mode1 = pose_dir / "mode_01.pdbqt"
    if mode1.exists() or (out.exists() and parse_aff(out) is not None):
        aff = parse_aff(mode1 if mode1.exists() else out)
        return {"target": target, "ligand": lig, "status": "cached", "gnina_mode1": aff, "seconds": 0}
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{env.get('LD_LIBRARY_PATH', '')}"
    cmd = [
        str(GNINA),
        "--no_gpu",
        "-r", str(rec),
        "-l", str(ligand),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", "8",
        "--num_modes", "9",
        "--seed", str(SEED),
        "--cpu", "1",
        "-o", str(out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    aff = parse_aff(out) if out.exists() else None
    if aff is not None and out.exists():
        # split MODEL 1 only for cache marker
        text = out.read_text(errors="replace").splitlines()
        cur, n = [], 0
        for line in text:
            if line.startswith("MODEL"):
                cur = [line]
            elif line.startswith("ENDMDL"):
                cur.append(line)
                n += 1
                (pose_dir / f"mode_{n:02d}.pdbqt").write_text("\n".join(cur) + "\n")
                cur = []
            elif cur:
                cur.append(line)
    return {
        "target": target,
        "ligand": lig,
        "status": "ok" if aff is not None else f"FAIL rc={proc.returncode}",
        "gnina_mode1": aff if aff is not None else "",
        "seconds": round(time.time() - t0, 1),
    }


def main() -> int:
    if not GNINA.exists():
        raise SystemExit(f"missing {GNINA}")
    OUT.mkdir(parents=True, exist_ok=True)
    boxes = {pdb: json.loads((BOXDIR / f"{pdb}_box.json").read_text()) for pdb in ("3POZ", "3RCD")}
    ligs = [r["panel_id"] for r in csv.DictReader(PANEL.open())]
    jobs = [(pdb, lig) for lig in ligs for pdb in ("3POZ", "3RCD")]
    print(f"gnina independent jobs={len(jobs)} workers={WORKERS}", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(dock, pdb, lig, boxes[pdb]) for pdb, lig in jobs]
        done = 0
        for fut in as_completed(futs):
            rec = fut.result()
            results.append(rec)
            done += 1
            print(f"[{done}/{len(jobs)}] {rec['target']} {rec['ligand']} {rec['status']} {rec.get('gnina_mode1','')} {rec.get('seconds')}s", flush=True)
            if done % 10 == 0:
                dest = OUT / "gnina_dock_scores_EGFR_HER2.csv"
                with dest.open("w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
                    w.writeheader()
                    w.writerows(results)
    dest = OUT / "gnina_dock_scores_EGFR_HER2.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(sorted(results, key=lambda r: (r["ligand"], r["target"])))
    n_ok = sum(1 for r in results if r["status"] in {"ok", "cached"})
    print(f"done {n_ok}/{len(results)} -> {dest}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
