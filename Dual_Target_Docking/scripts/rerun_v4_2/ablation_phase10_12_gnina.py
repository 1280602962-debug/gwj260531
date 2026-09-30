#!/usr/bin/env python3
"""Phase 10-12: GNINA CPU calibration, cognate QC, 1592 seed42 independent jobs.

Does not compute AUROC. Starts only when invoked (watcher waits for alt production).
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    BOX,
    GIND,
    GNINA_BIN,
    GNINA_LIB,
    LIG_GNINA_PDBQT,
    LIG_PDBQT,
    LIG_SDF,
    QA,
    REC,
    UNIQUE_14,
)
try:
    from rerun_v4_2.phase5_cognate_redock import calc_pose_rmsds, crystal_mol  # noqa: E402
except ImportError:
    calc_pose_rmsds = None
    crystal_mol = None

MANIFEST = GIND / "gnina_seed42_job_manifest.csv"
CALIB = GIND / "GNINA_RUNTIME_CALIBRATION.csv"
QC = GIND / "GNINA_COGNATE_REDOCKING_QC.csv"
LEDGER = GIND / "GNINA_SEED42_PRODUCTION_LEDGER.csv"
MASTER = GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv"
JOBS = GIND / "jobs"
RUNTIME_FREEZE = QA / "GNINA_RUNTIME_FREEZE.json"
COGNATE_LIG = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/04_cognate_redocking/ligands"

import re
CNN_RE = {
    "CNNscore": re.compile(r"CNNscore[:\s]+([0-9eE.+-]+)"),
    "CNNaffinity": re.compile(r"CNNaffinity[:\s]+([0-9eE.+-]+)"),
    "CNN_VS": re.compile(r"CNN_?VS[:\s]+([0-9eE.+-]+)"),
    "affinity": re.compile(r"(?:Affinity|VINA RESULT:|minimizedAffinity)[:\s]+([0-9eE.+-]+)"),
}


def env() -> dict:
    e = os.environ.copy()
    e["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{e.get('LD_LIBRARY_PATH', '')}"
    return e


def parse_scores(text: str) -> dict:
    out = {"CNNscore": "", "CNNaffinity": "", "CNN_VS": "", "empirical_affinity": ""}
    for k, rx in CNN_RE.items():
        m = rx.search(text)
        if not m:
            continue
        if k == "affinity":
            out["empirical_affinity"] = m.group(1)
        else:
            out[k if k != "CNN_VS" else "CNN_VS"] = m.group(1)
    return out


def ligand_complexity(eid: str) -> tuple[int, int]:
    sdf = LIG_SDF / f"{eid}.sdf"
    mol = Chem.MolFromMolFile(str(sdf), sanitize=True) if sdf.is_file() else None
    if mol is None:
        return (9999, 9999)
    heavy = mol.GetNumHeavyAtoms()
    rot = int(rdMolDescriptors.CalcNumRotatableBonds(mol))
    return heavy, rot


def pick_calib_ligands(rows: list[dict]) -> list[tuple[str, str, str]]:
    by_rec = {}
    for r in rows:
        by_rec.setdefault(r["pdb_id"], set()).add(r["global_ligand_entity_id"])
    picked = []
    for pdb in UNIQUE_14:
        eids = sorted(by_rec.get(pdb, []))
        ranked = sorted(eids, key=lambda e: (*ligand_complexity(e), e))
        if not ranked:
            continue
        lo, hi = ranked[0], ranked[-1]
        picked.append((pdb, lo, "low"))
        if hi != lo:
            picked.append((pdb, hi, "high"))
        else:
            picked.append((pdb, lo, "high_same"))
    return picked


# Meeko 0.7.1 extra types rejected by GNINA 1.3.2. Map/drop only; coords of real atoms unchanged.
_GNINA_TYPE_MAP = {"CG0": "C"}
_GNINA_DROP_TYPES = {"G0"}


def write_gnina_ad4_ligand(src: Path, dest: Path) -> dict:
    n_map = n_drop = 0
    out_lines = []
    for ln in src.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            typ = ln[77:].strip() if len(ln) >= 78 else ln.split()[-1]
            if typ in _GNINA_DROP_TYPES:
                n_drop += 1
                continue
            if typ in _GNINA_TYPE_MAP:
                n_map += 1
                mapped = _GNINA_TYPE_MAP[typ]
                if len(ln) >= 78:
                    ln = ln[:77] + f"{mapped:<3}"
                else:
                    parts = ln.rsplit(None, 1)
                    ln = parts[0] + f" {mapped}"
        out_lines.append(ln)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return {"src": src.name, "n_type_mapped": n_map, "n_dummy_dropped": n_drop}


def materialize_gnina_ligands() -> int:
    LIG_GNINA_PDBQT.mkdir(parents=True, exist_ok=True)
    n = 0
    changed = []
    for src in sorted(LIG_PDBQT.glob("*.pdbqt")):
        dest = LIG_GNINA_PDBQT / src.name
        rec = write_gnina_ad4_ligand(src, dest)
        n += 1
        if rec["n_type_mapped"] or rec["n_dummy_dropped"]:
            changed.append(rec)
    note = {
        "role": "GNINA_reader_format_conversion_only",
        "frozen_ligands_untouched": str(LIG_PDBQT),
        "n_copied": n,
        "n_with_cg0_g0": len(changed),
        "map": _GNINA_TYPE_MAP,
        "drop_dummy": sorted(_GNINA_DROP_TYPES),
        "coordinates_of_real_atoms": "unchanged",
        "chemistry_reprepared": False,
    }
    (GIND / "GNINA_LIGAND_AD4_TYPE_NOTE.json").write_text(json.dumps(note, indent=2) + "\n")
    print("gnina_ligand_compat", n, "changed", len(changed), flush=True)
    return n


def gnina_cmd(receptor: Path, ligand: Path, box: dict, out: Path, seed: int = 42) -> list[str]:
    return [
        str(GNINA_BIN), "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--no_gpu", "--scoring", "vina", "--cnn_scoring", "rescore",
        "--pose_sort_order", "CNNscore", "--seed", str(seed),
        "--exhaustiveness", "16", "--num_modes", "9", "--cpu", "1",
        "--out", str(out),
    ]


def run_gnina(cmd: list[str], log: Path, timeout_s: int) -> tuple[str | int, float, dict]:
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env(), timeout=timeout_s)
        dt = time.time() - t0
        log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
        scores = parse_scores((proc.stdout or "") + "\n" + (proc.stderr or ""))
        if out := Path(cmd[cmd.index("--out") + 1]):
            if out.is_file():
                extra = parse_scores(out.read_text(errors="replace"))
                for k, v in extra.items():
                    if not scores.get(k) and v:
                        scores[k] = v
        return proc.returncode, dt, scores
    except subprocess.TimeoutExpired:
        log.write_text("JOB_TIMEOUT\n", encoding="utf-8")
        return "TIMEOUT", time.time() - t0, {}


def phase10(rows: list[dict]) -> dict:
    picked = pick_calib_ligands(rows)
    recs = []
    # conservative timeout for calibration itself: 30 min
    for pdb, eid, kind in picked:
        box = json.loads((BOX / f"{pdb}_box.json").read_text())
        jobdir = GIND / "calibration" / f"{pdb}__{eid}"
        jobdir.mkdir(parents=True, exist_ok=True)
        out = jobdir / "out.pdbqt"
        log = jobdir / "gnina.log"
        cmd = gnina_cmd(REC / f"{pdb}_receptor.pdbqt", LIG_GNINA_PDBQT / f"{eid}.pdbqt", box, out)
        rc, dt, sc = run_gnina(cmd, log, 1800)
        if rc != 0 and rc != "TIMEOUT":
            rc, dt, sc = run_gnina(cmd, log, 1800)
        recs.append({
            "pdb_id": pdb, "global_ligand_entity_id": eid, "complexity": kind,
            "heavy_atoms": ligand_complexity(eid)[0], "rotatable_bonds": ligand_complexity(eid)[1],
            "runtime_s": round(dt, 3), "status": "TIMEOUT" if rc == "TIMEOUT" else ("SUCCESS" if rc == 0 else "FAIL"),
            "CNNscore": sc.get("CNNscore", ""),
        })
        print("calib", pdb, eid, recs[-1]["status"], recs[-1]["runtime_s"], flush=True)
    with CALIB.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(recs[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(recs)
    ok_rt = [r["runtime_s"] for r in recs if r["status"] == "SUCCESS"]
    if not ok_rt:
        raise SystemExit("calibration produced no SUCCESS runtimes")
    ok_rt.sort()
    p95 = ok_rt[max(0, math.ceil(0.95 * len(ok_rt)) - 1)]
    timeout = int(math.ceil((3 * p95) / 300.0) * 300)
    # concurrency: 8 logical CPUs, reserve 2, each job cpu=1; freeze 2 for CPU CNN memory
    conc = 2
    freeze = {
        "n_calib": len(recs),
        "n_success": len(ok_rt),
        "P95_s": p95,
        "GNINA_JOB_TIMEOUT_s": timeout,
        "GNINA_GLOBAL_CONCURRENCY": conc,
        "formula": "ceil_to_5min(3*P95)",
        "AUROC_used": False,
        "score_used": False,
    }
    RUNTIME_FREEZE.write_text(json.dumps(freeze, indent=2) + "\n")
    print(json.dumps(freeze, indent=2), flush=True)
    return freeze


def phase11() -> None:
    from rerun_v4_2.ablation_config import RUN
    ccd = {}
    # use already saved cognate crystals + primary cognate pdbqt
    rows = []
    for pdb in UNIQUE_14:
        box = json.loads((BOX / f"{pdb}_box.json").read_text())
        lig = COGNATE_LIG / f"{pdb}_cognate.pdbqt"
        jobdir = GIND / "cognate_qc" / pdb
        jobdir.mkdir(parents=True, exist_ok=True)
        out = jobdir / "out.pdbqt"
        cmd = gnina_cmd(REC / f"{pdb}_receptor.pdbqt", lig, box, out)
        rc, dt, sc = run_gnina(cmd, jobdir / "gnina.log", 1800)
        if rc != 0 and rc != "TIMEOUT":
            rc, dt, sc = run_gnina(cmd, jobdir / "gnina.log", 1800)
        rec = {
            "pdb_id": pdb, "seed": 42, "status": "TIMEOUT" if rc == "TIMEOUT" else ("SUCCESS" if rc == 0 and out.is_file() else "FAIL"),
            "runtime_s": round(dt, 3), "CNNscore": sc.get("CNNscore", ""),
            "top1_rmsd_A": "", "min_saved_pose_rmsd_A": "", "n_modes": "",
            "receptor_modified_from_rmsd": 0, "box_modified_from_rmsd": 0,
        }
        # RMSD if we can build crystal mol from cognate json + any smiles in crystal prep sdf
        sdf = RUN / "04_cognate_redocking" / "ligands" / f"{pdb}_cognate.sdf"
        if rec["status"] == "SUCCESS" and sdf.is_file() and crystal_mol and calc_pose_rmsds:
            suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
            mol0 = next((m for m in suppl if m is not None), None)
            if mol0 is not None:
                smi = Chem.MolToSmiles(Chem.RemoveHs(mol0))
                ref, why = crystal_mol(pdb, smi)
                if ref is not None:
                    rms = calc_pose_rmsds(out, ref)
                    rec["top1_rmsd_A"] = rms.get("top1_rmsd_A", "")
                    rec["min_saved_pose_rmsd_A"] = rms.get("min_saved_pose_rmsd_A", "")
                    rec["n_modes"] = rms.get("n_modes_returned", "")
                    rec["crystal_ref"] = why
        rows.append(rec)
        print("qc", pdb, rec["status"], rec.get("top1_rmsd_A"), flush=True)
    with QC.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0].keys()), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"] == "SUCCESS")
    if n_ok < 14:
        raise SystemExit(f"cognate QC SUCCESS {n_ok}/14")


def phase12(rows: list[dict], freeze: dict) -> None:
    timeout = int(freeze["GNINA_JOB_TIMEOUT_s"])
    workers = int(freeze["GNINA_GLOBAL_CONCURRENCY"])
    JOBS.mkdir(parents=True, exist_ok=True)
    pending = []
    done = {}
    for r in rows:
        jid = f"{r['pdb_id']}__{r['global_ligand_entity_id']}__seed42"
        st = JOBS / jid / "status.json"
        if st.is_file():
            prev = json.loads(st.read_text())
            if prev.get("status") in {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"}:
                done[jid] = prev
                continue
        pending.append({**r, "job_id": jid})
    print(f"gnina_prod total={len(rows)} cached={len(done)} pending={len(pending)} workers={workers} timeout={timeout}", flush=True)

    def work(job: dict) -> dict:
        jobdir = JOBS / job["job_id"]
        jobdir.mkdir(parents=True, exist_ok=True)
        box = json.loads((BOX / f"{job['pdb_id']}_box.json").read_text())
        lig = LIG_GNINA_PDBQT / f"{job['global_ligand_entity_id']}.pdbqt"
        out = jobdir / "out.pdbqt"
        rec = {
            "job_id": job["job_id"], "pdb_id": job["pdb_id"],
            "global_ligand_entity_id": job["global_ligand_entity_id"],
            "pairs": job["pairs"], "seed": 42, "status": "",
            "CNNscore": "", "CNNaffinity": "", "CNN_VS": "",
            "empirical_affinity": "", "M3_empirical_higher_better": "",
            "runtime_s": "", "n_retries": 0, "reason": "",
        }
        if not lig.is_file():
            rec.update(status="PERMANENT_FAIL", reason="missing_ligand")
            (jobdir / "status.json").write_text(json.dumps(rec, indent=2) + "\n")
            return rec
        cmd = gnina_cmd(REC / f"{job['pdb_id']}_receptor.pdbqt", lig, box, out)
        rc, dt, sc = run_gnina(cmd, jobdir / "gnina.log", timeout)
        if rc != 0:
            rec["n_retries"] = 1
            rc, dt, sc = run_gnina(cmd, jobdir / "gnina.log", timeout)
        rec["runtime_s"] = round(dt, 3)
        rec.update(sc)
        if sc.get("empirical_affinity"):
            rec["M3_empirical_higher_better"] = f"{-float(sc['empirical_affinity']):.6f}"
        if rc == "TIMEOUT":
            rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
        elif rc != 0 or not out.is_file() or not sc.get("CNNscore"):
            rec.update(status="PERMANENT_FAIL", reason=f"rc={rc}")
        else:
            rec.update(status="SUCCESS")
        (jobdir / "status.json").write_text(json.dumps(rec, indent=2) + "\n")
        return rec

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work, j): j for j in pending}
        for fut in as_completed(futs):
            rec = fut.result()
            done[rec["job_id"]] = rec
            if len(done) % 10 == 0 or rec["status"] != "SUCCESS":
                print(f"{len(done)}/{len(rows)} {rec['status']} {rec['pdb_id']} {rec['global_ligand_entity_id']} {rec['runtime_s']}", flush=True)

    rows_out = [done[k] for k in sorted(done)]
    fields = [
        "job_id", "pdb_id", "global_ligand_entity_id", "pairs", "seed", "status",
        "CNNscore", "CNNaffinity", "CNN_VS", "empirical_affinity",
        "M3_empirical_higher_better", "runtime_s", "n_retries", "reason",
    ]
    with LEDGER.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows_out)
    with MASTER.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows_out)
    counts = {}
    for r in rows_out:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    (QA / "GNINA_INDEPENDENT_COMPLETENESS_AUDIT.md").write_text(
        "# GNINA INDEPENDENT COMPLETENESS AUDIT\n\n"
        + "\n".join(f"- {k}: **{v}**" for k, v in {"expected": 1592, "accounted": len(rows_out), **counts, "AUROC": "NOT_COMPUTED"}.items())
        + "\n",
        encoding="utf-8",
    )
    print("PHASE12_DONE", counts, flush=True)


def main() -> int:
    materialize_gnina_ligands()
    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig", newline="")))
    if RUNTIME_FREEZE.is_file():
        freeze = json.loads(RUNTIME_FREEZE.read_text())
        print("reuse_runtime_freeze", freeze, flush=True)
    else:
        freeze = phase10(rows)
    if not (QC.is_file() and sum(1 for r in csv.DictReader(QC.open()) if r.get("status") == "SUCCESS") == 14):
        phase11()
    phase12(rows, freeze)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
