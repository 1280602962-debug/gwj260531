#!/usr/bin/env python3
"""Independent GNINA search on the uniform EGFR/HER2 RDKit/Meeko ligands.

Same frozen receptors and corrected boxes as the uniform Vina rebuild.
Does not redock PIK3CA/mTOR or JAK1/TYK2.
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/egfr_her2_uniform_rdkit_v1"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
BOX_DIR = ROOT / "data/egfr_her2_panel120_v0/boxes"
REC_DIR = OUT / "receptors"
LIG_DIR = OUT / "ligands_pdbqt"
WORK = OUT / "gnina_work"
GNINA = Path("/mnt/d/CADD paper exercise/gnina/bin/gnina")
GNINA_LIB = Path("/mnt/d/CADD paper exercise/gnina/conda_env/lib")
SEED = 20260727
TARGETS = ("3POZ", "3RCD")
TIMEOUT_S = int(os.environ.get("GNINA_TIMEOUT_S", "600"))
AFF_RE = re.compile(r"REMARK minimizedAffinity\s+([-\d.]+)")
VINA_RE = re.compile(r"REMARK VINA RESULT:\s+([-\d.]+)")


def load_box(pdb: str) -> dict:
    return json.loads((BOX_DIR / f"{pdb}_box_corrected.json").read_text())


def parse_affinity(pdbqt: Path) -> float | None:
    if not pdbqt.is_file():
        return None
    text = pdbqt.read_text(errors="replace")
    match = AFF_RE.search(text) or VINA_RE.search(text)
    return float(match.group(1)) if match else None


def dock_one(job: dict) -> dict:
    t0 = time.time()
    out_dir = WORK / job["pdb"] / job["ligand"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = out_dir / "gnina_out.pdbqt"
    log_path = out_dir / "gnina.log"
    rec = {
        "target": job["pdb"],
        "ligand": job["ligand"],
        "status": "fail",
        "gnina_mode1": "",
        "seconds": "",
        "reason": "",
    }
    lig = LIG_DIR / f"{job['ligand']}.pdbqt"
    recp = REC_DIR / f"{job['pdb']}_receptor.pdbqt"
    if not lig.is_file():
        rec["reason"] = "missing_ligand_pdbqt"
        rec["seconds"] = round(time.time() - t0, 1)
        return rec
    if out_pdbqt.is_file() and out_pdbqt.stat().st_size > 0:
        aff = parse_affinity(out_pdbqt)
        if aff is not None:
            rec["gnina_mode1"] = aff
            rec["status"] = "ok"
            rec["seconds"] = 0.0
            rec["reason"] = "cached"
            return rec
    box = load_box(job["pdb"])
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{env.get('LD_LIBRARY_PATH', '')}"
    cmd = [
        str(GNINA),
        "--no_gpu",
        "-r",
        str(recp),
        "-l",
        str(lig),
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--exhaustiveness",
        "8",
        "--num_modes",
        "9",
        "--seed",
        str(SEED),
        "--cpu",
        "1",
        "-o",
        str(out_pdbqt),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        rec["seconds"] = round(time.time() - t0, 1)
        rec["status"] = "timeout_skipped"
        rec["reason"] = f"timeout_{TIMEOUT_S}s"
        log_path.write_text(f"TIMEOUT {TIMEOUT_S}s\n")
        if out_pdbqt.exists() and out_pdbqt.stat().st_size == 0:
            out_pdbqt.unlink()
        return rec
    log_path.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    aff = parse_affinity(out_pdbqt)
    rec["seconds"] = round(time.time() - t0, 1)
    if proc.returncode != 0 or aff is None:
        rec["reason"] = (proc.stderr or proc.stdout or "gnina_fail")[-300:]
        return rec
    rec["gnina_mode1"] = aff
    rec["status"] = "ok"
    return rec


def main() -> int:
    workers = int(os.environ.get("GNINA_WORKERS", "4"))
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8-sig", newline="")))
    jobs = [{"ligand": r["panel_id"], "pdb": pdb} for r in panel for pdb in TARGETS]
    print(f"gnina jobs={len(jobs)} workers={workers} timeout_s={TIMEOUT_S}", flush=True)
    rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(dock_one, job) for job in jobs]
        for fut in as_completed(futs):
            rec = fut.result()
            rows.append(rec)
            done += 1
            if done % 10 == 0 or rec["status"] != "ok":
                n_ok = sum(1 for r in rows if r["status"] == "ok")
                n_skip = sum(1 for r in rows if r["status"] == "timeout_skipped")
                print(
                    f"[{done}/{len(jobs)}] {rec['status']} {rec['target']} {rec['ligand']} "
                    f"aff={rec['gnina_mode1']} ok={n_ok} timeout_skipped={n_skip}",
                    flush=True,
                )
    rows.sort(key=lambda r: (r["ligand"], r["target"]))
    dest = OUT / "tables"
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / "gnina_dock_scores_EGFR_HER2.csv"
    fields = ["target", "ligand", "status", "gnina_mode1", "seconds", "reason"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_skip = sum(1 for r in rows if r["status"] == "timeout_skipped")
    n_fail = len(rows) - n_ok - n_skip
    print(f"wrote {path} ok={n_ok} timeout_skipped={n_skip} fail={n_fail}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
