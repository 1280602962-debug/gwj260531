#!/usr/bin/env python3
"""Seed42-only production docking on the eight alternative receptors.

Reuses frozen V4.2 ligand pdbqt. Does not redock unreplaced primaries.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.alt_receptor_config import (  # noqa: E402
    ALTS,
    FORMAL_PAIRS,
    JOB_TIMEOUT_S,
    PRODUCTION_SEED,
    RUN,
    VINA,
    WORKERS,
)
from rerun_v4_2.phase2_prepare_receptors import assert_v42  # noqa: E402
from rerun_v4_2.phase7_fiveseed_production import parse_vina_out  # noqa: E402

REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
LIG = RUN / "02_ligands" / "pdbqt"
QA = RUN / "13_qa"
PROTO = RUN / "00_protocol"
MAP = PROTO / "pair_ligand_mapping.csv"
JOBS = RUN / "08_alt_vina_seed42" / "jobs"
LEDGER = QA / "alt_production_seed42_ledger.jsonl"
MASTER = QA / "alt_production_seed42_master.csv"
STATUS = QA / "alt_production_status.json"
DONE = {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"}


def pair_safe(pair: str) -> str:
    return pair.replace("/", "_")


def job_id(pair: str, pdb: str, eid: str) -> str:
    return f"{pair_safe(pair)}__{pdb}__{eid}__seed{PRODUCTION_SEED}"


def build_jobs() -> list[dict]:
    rows = [
        r for r in csv.DictReader(MAP.open(encoding="utf-8-sig", newline=""))
        if r["pair"] in FORMAL_PAIRS
        and r.get("parent_independent") == "1"
        and r.get("prepare_3d") == "1"
    ]
    jobs = []
    for r in rows:
        pair = r["pair"]
        eid = r["global_ligand_entity_id"]
        for pdb, cfg in ALTS.items():
            if cfg["pair"] != pair:
                continue
            jobs.append({
                "job_id": job_id(pair, pdb, eid),
                "pair": pair,
                "pdb_id": pdb,
                "alt_id": cfg["id"],
                "side": cfg["side"],
                "design": cfg["design"],
                "global_ligand_entity_id": eid,
                "canonical_ligand_id": r["canonical_ligand_id"],
                "seed": PRODUCTION_SEED,
                "primary_score": 1,
            })
    return jobs


def load_done() -> dict[str, dict]:
    done: dict[str, dict] = {}
    if LEDGER.is_file():
        for line in LEDGER.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec.get("status") in DONE:
                done[rec["job_id"]] = rec
    for st in JOBS.glob("*/status.json"):
        rec = json.loads(st.read_text(encoding="utf-8"))
        if rec.get("status") in DONE:
            done[rec["job_id"]] = rec
    return done


def append_ledger(rec: dict) -> None:
    assert_v42(LEDGER)
    with LEDGER.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(rec, sort_keys=True) + "\n")


def write_progress(n_total: int, done: dict[str, dict], extra: dict | None = None) -> None:
    counts = {"SUCCESS": 0, "TIMEOUT": 0, "PERMANENT_FAIL": 0}
    for rec in done.values():
        counts[rec["status"]] = counts.get(rec["status"], 0) + 1
    payload = {
        "jobs_expected": n_total,
        "jobs_accounted": sum(counts.values()),
        **counts,
        "remaining": n_total - sum(counts.values()),
        "AUROC_computed": "NO",
        "unreplaced_primaries_redocked": "NO",
        "JOB_TIMEOUT_s": JOB_TIMEOUT_S,
        "max_concurrent": WORKERS,
        "updated_unix": time.time(),
    }
    if extra:
        payload.update(extra)
    STATUS.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_one(job: dict) -> dict:
    pdb = job["pdb_id"]
    eid = job["global_ligand_entity_id"]
    jobdir = JOBS / job["job_id"]
    jobdir.mkdir(parents=True, exist_ok=True)
    assert_v42(jobdir)
    box = json.loads((BOX / f"{pdb}_box.json").read_text())
    rec_p = REC / f"{pdb}_receptor.pdbqt"
    lig_p = LIG / f"{eid}.pdbqt"
    out = jobdir / "out.pdbqt"
    rec = {
        **job,
        "scoring": "vina", "exhaustiveness": 16, "num_modes": 9, "energy_range": 6, "cpu": 1,
        "timeout_s": JOB_TIMEOUT_S, "timeout_frozen": True, "n_retries": 0,
        "status": "", "runtime_s": "", "vina_rc": "", "mode1_affinity": "",
        "n_modes_returned": "", "reason": "", "rescue_params": 0,
    }
    if not lig_p.is_file():
        rec.update(status="PERMANENT_FAIL", reason="missing_ligand_pdbqt")
        (jobdir / "status.json").write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        return rec
    cmd = [
        str(VINA), "--receptor", str(rec_p), "--ligand", str(lig_p),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--scoring", "vina", "--exhaustiveness", "16", "--num_modes", "9",
        "--energy_range", "6", "--cpu", "1", "--seed", str(PRODUCTION_SEED),
        "--out", str(out),
    ]

    def once() -> tuple[int | str, float]:
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=JOB_TIMEOUT_S)
            dt = time.time() - t0
            (jobdir / "vina.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
            (jobdir / "vina.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
            return proc.returncode, dt
        except subprocess.TimeoutExpired:
            (jobdir / "vina.stderr.txt").write_text("JOB_TIMEOUT\n", encoding="utf-8")
            return "TIMEOUT", time.time() - t0

    rc, dt = once()
    if rc != 0:
        rec["n_retries"] = 1
        rc, dt = once()
    rec["runtime_s"] = round(dt, 3)
    rec["vina_rc"] = rc
    n_modes, mode1 = parse_vina_out(out)
    rec["n_modes_returned"] = n_modes
    rec["mode1_affinity"] = mode1
    if rc == "TIMEOUT":
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
    elif rc != 0 or not out.is_file() or out.stat().st_size == 0:
        rec.update(status="PERMANENT_FAIL", reason=f"vina_rc={rc}")
    else:
        rec.update(status="SUCCESS")
    (jobdir / "status.json").write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    return rec


def rewrite_master(done: dict[str, dict]) -> None:
    rows = [done[k] for k in sorted(done)]
    if not rows:
        return
    fields = [
        "job_id", "pair", "pdb_id", "alt_id", "side", "design",
        "global_ligand_entity_id", "canonical_ligand_id", "seed", "primary_score",
        "status", "runtime_s", "vina_rc", "n_retries", "mode1_affinity",
        "n_modes_returned", "reason", "timeout_s",
    ]
    with MASTER.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    JOBS.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    assert_v42(JOBS)
    jobs = build_jobs()
    done = load_done()
    pending = [j for j in jobs if j["job_id"] not in done]
    write_progress(len(jobs), done, {"pending_at_start": len(pending)})
    print(f"alt_prod jobs={len(jobs)} done={len(done)} pending={len(pending)}", flush=True)
    if not pending:
        rewrite_master(done)
        write_progress(len(jobs), done, {"PRODUCTION_COMPLETE": True})
        print("ALT_PRODUCTION_ALREADY_COMPLETE")
        return 0
    finished = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(run_one, job): job for job in pending}
        for fut in as_completed(futs):
            rec = fut.result()
            done[rec["job_id"]] = rec
            append_ledger(rec)
            finished += 1
            if finished % 10 == 0 or rec["status"] != "SUCCESS":
                write_progress(len(jobs), done)
                if finished % 25 == 0:
                    rewrite_master(done)
            print(
                f"{len(done)}/{len(jobs)} {rec['status']} {rec['pair']} {rec['pdb_id']} "
                f"{rec['global_ligand_entity_id']} {rec['runtime_s']}s",
                flush=True,
            )
    rewrite_master(done)
    write_progress(len(jobs), done, {"PRODUCTION_COMPLETE": True})
    print("ALT_PRODUCTION_DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
