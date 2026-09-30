#!/usr/bin/env python3
"""V4.2 Phase 6: runtime calibration only. No AUROC. Then freeze JOB_TIMEOUT."""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Lipinski

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.uniform_protocol_config import RECEPTORS, UNIQUE_14  # noqa: E402

RUN_ID = "UNIFORM_RERUN_V4_2_20260921"
RUN = ROOT / "reruns" / RUN_ID
REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
LIG = RUN / "02_ligands" / "pdbqt"
QA = RUN / "13_qa"
CAL = RUN / "05_vina" / "calibration"
PROTO = RUN / "00_protocol"
MAP = PROTO / "pair_ligand_mapping.csv"
VINA = Path("/home/gwj/miniconda3/bin/vina")
CAL_SEED = 42  # frozen primary seed; runtime only
CAL_TIMEOUT_S = 1800  # provisional until P95 freeze
WORKERS = 6
PDB_TO_PAIRS = defaultdict(list)
for pair, pdbs in RECEPTORS.items():
    for pdb in pdbs:
        PDB_TO_PAIRS[pdb].append(pair)


def assert_v42(path: Path) -> None:
    if RUN_ID not in str(path.resolve()):
        raise SystemExit(f"refusing write outside V4.2: {path}")


def complexity(smiles: str) -> tuple[int, int, int]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (-1, -1, -1)
    n_heavy = int(mol.GetNumHeavyAtoms())
    n_rot = int(Lipinski.NumRotatableBonds(mol))
    return (n_heavy + n_rot, n_heavy, n_rot)


def select_ligands() -> list[dict]:
    rows = [r for r in csv.DictReader(MAP.open(encoding="utf-8-sig", newline=""))
            if r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"]
    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        smi = r.get("isomeric_smiles") or r.get("canonical_smiles")
        score, n_heavy, n_rot = complexity(smi)
        if score < 0:
            continue
        rec = dict(r)
        rec["complexity_score"] = score
        rec["n_heavy"] = n_heavy
        rec["n_rotatable"] = n_rot
        rec["smiles_used"] = smi
        by_pair[r["pair"]].append(rec)
    chosen = []
    for pdb in UNIQUE_14:
        pool = []
        seen = set()
        for pair in PDB_TO_PAIRS[pdb]:
            for r in by_pair[pair]:
                eid = r["global_ligand_entity_id"]
                if eid in seen:
                    continue
                seen.add(eid)
                pool.append(r)
        pool.sort(key=lambda r: (r["complexity_score"], r["n_heavy"], r["n_rotatable"], r["global_ligand_entity_id"]))
        if len(pool) < 2:
            raise SystemExit(f"{pdb}: fewer than 2 dockable ligands")
        mid = pool[len(pool) // 2]
        hi = pool[min(len(pool) - 1, max(len(pool) // 2 + 1, int(round(0.90 * (len(pool) - 1)))))]
        if hi["global_ligand_entity_id"] == mid["global_ligand_entity_id"]:
            hi = pool[-1]
        for role, lig in (("medium", mid), ("higher", hi)):
            chosen.append({
                "pdb_id": pdb,
                "complexity_role": role,
                "pair": lig["pair"],
                "global_ligand_entity_id": lig["global_ligand_entity_id"],
                "canonical_ligand_id": lig["canonical_ligand_id"],
                "n_heavy": lig["n_heavy"],
                "n_rotatable": lig["n_rotatable"],
                "complexity_score": lig["complexity_score"],
                "n_pool": len(pool),
                "selection": "median_score" if role == "medium" else "p90_or_max_distinct",
            })
    return chosen


def run_job(sel: dict) -> dict:
    pdb = sel["pdb_id"]
    eid = sel["global_ligand_entity_id"]
    box = json.loads((BOX / f"{pdb}_box.json").read_text())
    rec_p = REC / f"{pdb}_receptor.pdbqt"
    lig_p = LIG / f"{eid}.pdbqt"
    jobdir = CAL / f"{pdb}_{sel['complexity_role']}_{eid}"
    jobdir.mkdir(parents=True, exist_ok=True)
    assert_v42(jobdir)
    out = jobdir / "out.pdbqt"
    rec = {
        **sel,
        "seed": CAL_SEED,
        "scoring": "vina",
        "exhaustiveness": 16,
        "num_modes": 9,
        "energy_range": 6,
        "cpu": 1,
        "status": "",
        "runtime_s": "",
        "vina_rc": "",
        "n_retries": 0,
        "reason": "",
        "auroc_computed": 0,
    }
    if not lig_p.is_file():
        rec.update(status="LIGAND_MISSING", reason=str(lig_p))
        return rec
    cmd = [
        str(VINA), "--receptor", str(rec_p), "--ligand", str(lig_p),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--scoring", "vina", "--exhaustiveness", "16", "--num_modes", "9",
        "--energy_range", "6", "--cpu", "1", "--seed", str(CAL_SEED),
        "--out", str(out),
    ]

    def once():
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=CAL_TIMEOUT_S)
            dt = time.time() - t0
            (jobdir / "vina.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
            (jobdir / "vina.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
            return proc.returncode, dt, ""
        except subprocess.TimeoutExpired:
            return "TIMEOUT", time.time() - t0, "JOB_TIMEOUT"

    rc, dt, err = once()
    if rc != 0:
        rec["n_retries"] = 1
        rc, dt, err = once()
    rec["runtime_s"] = round(dt, 3)
    rec["vina_rc"] = rc
    if rc == "TIMEOUT":
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
    elif rc != 0 or not out.is_file() or out.stat().st_size == 0:
        rec.update(status="FAIL", reason=err or f"vina_rc={rc}")
    else:
        rec.update(status="OK")
    (jobdir / "status.json").write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    return rec


def percentile(values: list[float], p: float) -> float:
    xs = sorted(values)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p / 100.0
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return xs[int(k)]
    return xs[f] * (c - k) + xs[c] * (k - f)


def freeze_timeout(ok_runtimes: list[float]) -> dict:
    p95 = percentile(ok_runtimes, 95)
    raw = 3.0 * p95
    timeout_s = int(round(min(30 * 60, max(15 * 60, raw))))
    rec = {
        "n_ok_runtimes": len(ok_runtimes),
        "p50_s": round(percentile(ok_runtimes, 50), 3),
        "p95_s": round(p95, 3),
        "three_p95_s": round(raw, 3),
        "JOB_TIMEOUT_s": timeout_s,
        "JOB_TIMEOUT_min": round(timeout_s / 60.0, 2),
        "formula": "clamp(3*P95, 15min, 30min)",
        "max_concurrent": WORKERS,
        "concurrency_reduced_for_memory": False,
        "timeout_frozen": True,
        "auroc_used": False,
    }
    rt = PROTO / "runtime_config.yaml"
    text = rt.read_text(encoding="utf-8")
    text = text.replace("timeout_frozen: false", "timeout_frozen: true")
    if "JOB_TIMEOUT_s:" not in text:
        text += (
            f"\nJOB_TIMEOUT_s: {timeout_s}\n"
            f"JOB_TIMEOUT_min: {rec['JOB_TIMEOUT_min']}\n"
            f"calibration_p95_s: {rec['p95_s']}\n"
            f"max_concurrent_jobs: {WORKERS}\n"
        )
    else:
        # rewrite later via append-only freeze file
        pass
    assert_v42(rt)
    rt.write_text(text, encoding="utf-8")
    af = PROTO / "analysis_freeze.yaml"
    aft = af.read_text(encoding="utf-8")
    aft = aft.replace("  timeout_calibration_required: true", "  timeout_calibration_required: false")
    aft = aft.replace("  timeout_frozen: false", "  timeout_frozen: true")
    if "JOB_TIMEOUT_s:" not in aft:
        aft = aft.replace(
            "  suggested_max_concurrent_jobs: 6\n",
            f"  suggested_max_concurrent_jobs: 6\n  JOB_TIMEOUT_s: {timeout_s}\n  calibration_p95_s: {rec['p95_s']}\n",
        )
    assert_v42(af)
    af.write_text(aft, encoding="utf-8")
    (QA / "JOB_TIMEOUT_FREEZE.json").write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    return rec


def main() -> int:
    CAL.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    chosen = select_ligands()
    sel_path = QA / "phase6_calibration_ligand_selection.csv"
    with sel_path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(chosen[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(chosen)
    print(f"selected {len(chosen)} calibration ligands", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(run_job, sel) for sel in chosen]
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            rows.append(rec)
            print(i, rec["pdb_id"], rec["complexity_role"], rec["status"], rec["runtime_s"], flush=True)
    fields = sorted({k for r in rows for k in r})
    with (QA / "phase6_calibration_runtimes.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["pdb_id"], r["complexity_role"])))
    ok = [float(r["runtime_s"]) for r in rows if r["status"] == "OK" and r["runtime_s"] != ""]
    n_ok = len(ok)
    n_fail = sum(1 for r in rows if r["status"] not in {"OK"})
    systematic = n_ok < 24 or n_fail >= 5
    freeze = freeze_timeout(ok) if n_ok else {}
    answers = {
        "RUN_ID": RUN_ID,
        "PHASE": 6,
        "jobs": len(rows),
        "OK": n_ok,
        "FAIL_OR_TIMEOUT": n_fail,
        "systematic_error": systematic,
        "AUROC_computed": "NO",
        "receptor_box_ptr_his_exhaustiveness_changed": "NO",
        "JOB_TIMEOUT_FREEZE": freeze,
        "PHASE_7_GO": "YES" if not systematic and n_ok >= 24 else "NO",
        "STOP_IF_SYSTEMATIC": systematic,
    }
    (QA / "phase6_status.json").write_text(json.dumps(answers, indent=2) + "\n", encoding="utf-8")
    md = ["# PHASE 6 RUNTIME CALIBRATION", "", f"RUN_ID = {RUN_ID}", ""]
    for k, v in answers.items():
        md.append(f"- {k}: **{v}**")
    (QA / "PHASE6_RUNTIME_CALIBRATION.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(answers, indent=2))
    print("PHASE6_DONE")
    return 1 if systematic else 0


if __name__ == "__main__":
    raise SystemExit(main())
