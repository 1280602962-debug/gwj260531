#!/usr/bin/env python3
"""V4.2 Phase 7: five-seed production Vina. No AUROC / summary_min / Top-K / LEGACY."""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.uniform_protocol_config import RECEPTORS  # noqa: E402

RUN_ID = "UNIFORM_RERUN_V4_2_20260921"
RUN = ROOT / "reruns" / RUN_ID
REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
LIG = RUN / "02_ligands" / "pdbqt"
QA = RUN / "13_qa"
JOBS = RUN / "06_vina_fiveseed" / "jobs"
PROTO = RUN / "00_protocol"
MAP = PROTO / "pair_ligand_mapping.csv"
VINA = Path("/home/gwj/miniconda3/bin/vina")
SEEDS = (17, 29, 42, 71, 101)
PRIMARY_SEED = 42
WORKERS = 6
LEDGER = QA / "phase7_fiveseed_ledger.jsonl"
MASTER = QA / "phase7_fiveseed_master.csv"
STATUS = QA / "phase7_status.json"
PROGRESS = QA / "phase7_progress.json"
DONE = {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"}


def assert_v42(path: Path) -> None:
    if RUN_ID not in str(path.resolve()):
        raise SystemExit(f"refusing write outside V4.2: {path}")


def load_timeout() -> int:
    freeze = json.loads((QA / "JOB_TIMEOUT_FREEZE.json").read_text())
    if not freeze.get("timeout_frozen"):
        raise SystemExit("JOB_TIMEOUT not frozen")
    return int(freeze["JOB_TIMEOUT_s"])


def pair_safe(pair: str) -> str:
    return pair.replace("/", "_")


def job_id(pair: str, pdb: str, eid: str, seed: int) -> str:
    return f"{pair_safe(pair)}__{pdb}__{eid}__seed{seed}"


def parse_vina_out(path: Path) -> tuple[int, str]:
    n = 0
    mode1 = ""
    if not path.is_file():
        return 0, ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "VINA RESULT" in line:
            n += 1
            if n == 1:
                parts = line.split()
                for i, tok in enumerate(parts):
                    if tok == "RESULT:" and i + 1 < len(parts):
                        mode1 = parts[i + 1]
                        break
    return n, mode1


def build_jobs() -> list[dict]:
    rows = [r for r in csv.DictReader(MAP.open(encoding="utf-8-sig", newline=""))
            if r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"]
    jobs = []
    for r in rows:
        pair = r["pair"]
        if pair not in RECEPTORS:
            continue
        eid = r["global_ligand_entity_id"]
        for pdb in RECEPTORS[pair]:
            for seed in SEEDS:
                jobs.append({
                    "job_id": job_id(pair, pdb, eid, seed),
                    "pair": pair,
                    "pdb_id": pdb,
                    "global_ligand_entity_id": eid,
                    "canonical_ligand_id": r["canonical_ligand_id"],
                    "seed": seed,
                    "primary_score": int(seed == PRIMARY_SEED),
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
        "RUN_ID": RUN_ID,
        "PHASE": 7,
        "jobs_expected": n_total,
        "jobs_accounted": sum(counts.values()),
        **counts,
        "remaining": n_total - sum(counts.values()),
        "AUROC_computed": "NO",
        "summary_min_computed": "NO",
        "topk_computed": "NO",
        "LEGACY_compared": "NO",
        "seeds": list(SEEDS),
        "primary": "seed42_mode1",
        "max_concurrent": WORKERS,
        "updated_unix": time.time(),
    }
    if extra:
        payload.update(extra)
    STATUS.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    PROGRESS.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_one(job: dict, timeout_s: int) -> dict:
    pdb = job["pdb_id"]
    eid = job["global_ligand_entity_id"]
    seed = job["seed"]
    jobdir = JOBS / job["job_id"]
    jobdir.mkdir(parents=True, exist_ok=True)
    assert_v42(jobdir)
    box = json.loads((BOX / f"{pdb}_box.json").read_text())
    rec_p = REC / f"{pdb}_receptor.pdbqt"
    lig_p = LIG / f"{eid}.pdbqt"
    out = jobdir / "out.pdbqt"
    rec = {
        **job,
        "scoring": "vina",
        "exhaustiveness": 16,
        "num_modes": 9,
        "energy_range": 6,
        "cpu": 1,
        "timeout_s": timeout_s,
        "timeout_frozen": True,
        "n_retries": 0,
        "status": "",
        "runtime_s": "",
        "vina_rc": "",
        "mode1_affinity": "",
        "n_modes_returned": "",
        "reason": "",
        "rescue_params": 0,
        "primary_seed": PRIMARY_SEED,
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
        "--energy_range", "6", "--cpu", "1", "--seed", str(seed),
        "--out", str(out),
    ]

    def once() -> tuple[int | str, float]:
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
            dt = time.time() - t0
            (jobdir / "vina.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
            (jobdir / "vina.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
            (jobdir / "vina.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
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
        "job_id", "pair", "pdb_id", "global_ligand_entity_id", "canonical_ligand_id",
        "seed", "primary_score", "status", "runtime_s", "vina_rc", "n_retries",
        "mode1_affinity", "n_modes_returned", "reason", "timeout_s",
    ]
    with MASTER.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    JOBS.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    assert_v42(JOBS)
    timeout_s = load_timeout()
    jobs = build_jobs()
    done = load_done()
    pending = [j for j in jobs if j["job_id"] not in done]
    write_progress(len(jobs), done, {"JOB_TIMEOUT_s": timeout_s, "pending_at_start": len(pending)})
    print(f"phase7 jobs={len(jobs)} done={len(done)} pending={len(pending)} timeout_s={timeout_s}", flush=True)
    if not pending:
        rewrite_master(done)
        write_progress(len(jobs), done, {"PHASE_7_COMPLETE": True})
        print("PHASE7_ALREADY_COMPLETE")
        return 0
    finished = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(run_one, job, timeout_s): job for job in pending}
        for fut in as_completed(futs):
            rec = fut.result()
            done[rec["job_id"]] = rec
            append_ledger(rec)
            finished += 1
            if finished % 10 == 0 or rec["status"] != "SUCCESS":
                write_progress(len(jobs), done, {"JOB_TIMEOUT_s": timeout_s})
                if finished % 50 == 0:
                    rewrite_master(done)
            print(
                f"{len(done)}/{len(jobs)} {rec['status']} {rec['pair']} {rec['pdb_id']} "
                f"{rec['global_ligand_entity_id']} seed{rec['seed']} {rec['runtime_s']}s",
                flush=True,
            )
    rewrite_master(done)
    write_progress(len(jobs), done, {"JOB_TIMEOUT_s": timeout_s, "PHASE_7_COMPLETE": True})
    print("PHASE7_DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
