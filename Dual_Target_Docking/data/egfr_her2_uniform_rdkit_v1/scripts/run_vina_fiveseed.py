#!/usr/bin/env python3
"""EGFR/HER2 five-seed Vina on uniform RDKit/Meeko ligands and frozen receptors.

110 ligands × 2 receptors × 5 seeds. Failures are recorded; seed/parameters are
not swapped.
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
WORK = OUT / "work" / "vina"
SEEDS = (20260727, 20260811, 20260812, 20260813, 20260814)
TARGETS = ("3POZ", "3RCD")
TIMEOUT_S = int(os.environ.get("VINA_TIMEOUT_S", "600"))
AFF_RE = re.compile(r"REMARK VINA RESULT:\s+([-\d.]+)")


def load_box(pdb: str) -> dict:
    box = json.loads((BOX_DIR / f"{pdb}_box_corrected.json").read_text())
    return box


def parse_affinity(pdbqt: Path, log_text: str) -> float | None:
    if pdbqt.is_file():
        text = pdbqt.read_text(errors="replace")
        match = AFF_RE.search(text)
        if match:
            return float(match.group(1))
        for line in text.splitlines():
            if line.startswith("REMARK minimizedAffinity"):
                return float(line.split()[2])
    for line in log_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("1 ") or re.match(r"^1\s+-", stripped):
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    pass
    return None


def dock_one(job: dict) -> dict:
    t0 = time.time()
    out_dir = WORK / str(job["seed"]) / job["pdb"] / job["ligand"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = out_dir / "out.pdbqt"
    log_path = out_dir / "vina.log"
    rec = dict(job)
    rec["vina_mode1"] = ""
    rec["status"] = "fail"
    rec["reason"] = ""
    rec["seconds"] = ""
    lig = LIG_DIR / f"{job['ligand']}.pdbqt"
    recp = REC_DIR / f"{job['pdb']}_receptor.pdbqt"
    if not lig.is_file():
        rec["reason"] = "missing_ligand_pdbqt"
        rec["seconds"] = round(time.time() - t0, 1)
        return rec
    if out_pdbqt.is_file() and out_pdbqt.stat().st_size > 0:
        aff = parse_affinity(out_pdbqt, log_path.read_text(errors="replace") if log_path.is_file() else "")
        if aff is not None:
            rec["vina_mode1"] = aff
            rec["status"] = "ok"
            rec["seconds"] = 0.0
            rec["reason"] = "cached"
            return rec
    box = load_box(job["pdb"])
    cmd = [
        "vina",
        "--receptor",
        str(recp),
        "--ligand",
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
        "--energy_range",
        "3",
        "--cpu",
        "1",
        "--seed",
        str(job["seed"]),
        "--out",
        str(out_pdbqt),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired as exc:
        rec["seconds"] = round(time.time() - t0, 1)
        rec["status"] = "timeout_skipped"
        rec["reason"] = f"timeout_{TIMEOUT_S}s"
        out_txt = ""
        if exc.stdout:
            out_txt += exc.stdout if isinstance(exc.stdout, str) else exc.stdout.decode("utf-8", "replace")
        if exc.stderr:
            out_txt += "\n" + (exc.stderr if isinstance(exc.stderr, str) else exc.stderr.decode("utf-8", "replace"))
        log_path.write_text(out_txt + f"\nTIMEOUT {TIMEOUT_S}s\n")
        if out_pdbqt.exists() and out_pdbqt.stat().st_size == 0:
            out_pdbqt.unlink()
        return rec
    log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    log_path.write_text(log_text)
    aff = parse_affinity(out_pdbqt, log_text)
    rec["seconds"] = round(time.time() - t0, 1)
    if proc.returncode != 0 or aff is None:
        rec["reason"] = (proc.stderr or proc.stdout or "vina_fail")[-300:]
        rec["status"] = "fail"
        return rec
    rec["vina_mode1"] = aff
    rec["status"] = "ok"
    return rec


def main() -> int:
    workers = int(os.environ.get("VINA_WORKERS", "7"))
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8-sig", newline="")))
    jobs = []
    for row in panel:
        for seed in SEEDS:
            for pdb in TARGETS:
                jobs.append(
                    {
                        "pair": "EGFR/HER2",
                        "seed": seed,
                        "ligand": row["panel_id"],
                        "class": row["class"],
                        "pdb": pdb,
                    }
                )
    print(f"vina jobs={len(jobs)} workers={workers} timeout_s={TIMEOUT_S}", flush=True)
    rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(dock_one, job) for job in jobs]
        for fut in as_completed(futs):
            rec = fut.result()
            rows.append(rec)
            done += 1
            if done % 20 == 0 or rec["status"] != "ok":
                n_ok = sum(1 for r in rows if r["status"] == "ok")
                n_skip = sum(1 for r in rows if r["status"] == "timeout_skipped")
                print(
                    f"[{done}/{len(jobs)}] {rec['status']} seed={rec['seed']} "
                    f"{rec['pdb']} {rec['ligand']} aff={rec['vina_mode1']} "
                    f"ok={n_ok} timeout_skipped={n_skip}",
                    flush=True,
                )
    rows.sort(key=lambda r: (int(r["seed"]), r["pdb"], r["ligand"]))
    dest = OUT / "tables"
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / "scores_vina_mode1_fiveseed.csv"
    fields = ["pair", "seed", "ligand", "class", "pdb", "status", "vina_mode1", "seconds", "reason"]
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
