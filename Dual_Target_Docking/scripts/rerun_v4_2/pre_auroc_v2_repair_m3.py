#!/usr/bin/env python3
"""M3 topology repair for pre-frozen FAIL set only. No score cherry-pick. No AUROC."""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    BOX,
    GIND,
    GNINA_BIN,
    GNINA_LIB,
    LIG_SDF,
    QA,
    REC,
    SEC,
)
from rerun_v4_2.pre_auroc_v2_lib import (  # noqa: E402
    finite_float,
    load_frozen_sdf,
    parse_model1_remarks,
    parse_pdbqt_atoms,
    sdf_heavy,
    special_types,
)
from rerun_v4_2.pre_auroc_v2_rebuild_m1 import meeko_prepare, validate_meeko_output  # noqa: E402

MEEKO = Path("/home/gwj/miniconda3/bin/mk_prepare_ligand.py")
TOPO = QA / "M3_MACROCYCLE_TOPOLOGY_DECISION.csv"
AUDIT = QA / "GNINA_LIGAND_REPRESENTATION_AUDIT.csv"
REPARSED = GIND / "GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv"
OLD_MASTER = GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv"
LIG_OUT = GIND / "ligands_gnina_rigid_macrocycles_v2"
JOB_OUT = GIND / "jobs_topology_repair_v2"
VERIFIED = GIND / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv"
TIMEOUT_S = 2400
WORKERS = 2
MANIFEST = GIND / "gnina_seed42_job_manifest.csv"


def env() -> dict:
    e = os.environ.copy()
    e["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{e.get('LD_LIBRARY_PATH', '')}"
    return e


def prepare_ligand(eid: str) -> dict:
    dest = LIG_OUT / f"{eid}.pdbqt"
    sdf = load_frozen_sdf(LIG_SDF / f"{eid}.sdf")
    rec = {"global_ligand_entity_id": eid, "status": "", "reason": "", "n_heavy": "", "special": ""}
    if sdf is None:
        rec.update(status="FAIL", reason="frozen_sdf_unreadable")
        return rec
    heavy = sdf_heavy(sdf)
    # Frozen SDF already has the official 3D + protonation. Materialize Hs if needed only.
    mol = sdf
    if any(a.GetAtomicNum() == 1 for a in sdf.GetAtoms()):
        mol = sdf
    else:
        mol = Chem.AddHs(sdf, addCoords=True)
    st, log = meeko_prepare(mol if mol.GetNumAtoms() == heavy.GetNumAtoms() else mol, dest)
    # meeko_prepare AddHs again internally; pass heavy-with-coords if SDF already explicit
    # Re-run on frozen SDF file directly to keep official starting coords of ALL atoms including H.
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(MEEKO), "-i", str(LIG_SDF / f"{eid}.sdf"), "-o", str(dest), "--rigid_macrocycles", "--add_index_map"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log2 = (proc.stdout or "") + "\n" + (proc.stderr or "")
    (LIG_OUT / f"{eid}.meeko.log").write_text(log2, encoding="utf-8")
    if proc.returncode != 0 or not dest.is_file():
        rec.update(status="FAIL", reason="meeko_fail")
        return rec
    atoms = parse_pdbqt_atoms(dest.read_text(errors="replace"))
    spec = special_types(atoms)
    rec["special"] = ";".join(spec)
    rec["n_heavy"] = str(sum(1 for a in atoms if a["type"] not in {"H", "HD", "HS"} and not a["is_glue"]))
    if spec:
        rec.update(status="FAIL", reason="special_type_remain")
        return rec
    if int(rec["n_heavy"]) != heavy.GetNumAtoms():
        rec.update(status="FAIL", reason=f"heavy_count {rec['n_heavy']}!={heavy.GetNumAtoms()}")
        return rec
    rec["status"] = "OK"
    return rec


def gnina_cmd(pdb: str, eid: str, out: Path) -> list[str]:
    box = json.loads((BOX / f"{pdb}_box.json").read_text())
    return [
        str(GNINA_BIN), "--receptor", str(REC / f"{pdb}_receptor.pdbqt"),
        "--ligand", str(LIG_OUT / f"{eid}.pdbqt"),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--no_gpu", "--scoring", "vina", "--cnn_scoring", "rescore",
        "--pose_sort_order", "CNNscore", "--seed", "42",
        "--exhaustiveness", "16", "--num_modes", "9", "--cpu", "1",
        "--out", str(out),
    ]


def run_job(pdb: str, eid: str, pairs: str) -> dict:
    job_id = f"{pdb}__{eid}__seed42"
    jobdir = JOB_OUT / job_id
    jobdir.mkdir(parents=True, exist_ok=True)
    st_p = jobdir / "status.json"
    if st_p.is_file():
        prev = json.loads(st_p.read_text())
        if prev.get("status") in {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"}:
            return prev
    out = jobdir / "out.pdbqt"
    log = jobdir / "gnina.log"
    cmd = gnina_cmd(pdb, eid, out)
    rec = {
        "job_id": job_id, "pdb_id": pdb, "global_ligand_entity_id": eid,
        "pairs": pairs, "seed": 42, "status": "", "CNNscore": "",
        "CNNaffinity": "", "CNN_VS": "", "empirical_affinity": "",
        "M3_empirical_higher_better": "", "runtime_s": "", "n_retries": 0,
        "reason": "", "source": "topology_repaired_rigid_macrocycles",
    }

    def once() -> tuple[str | int, float]:
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, env=env(), timeout=TIMEOUT_S)
            dt = time.time() - t0
            log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
            return proc.returncode, dt
        except subprocess.TimeoutExpired:
            log.write_text("JOB_TIMEOUT\n", encoding="utf-8")
            return "TIMEOUT", time.time() - t0

    rc, dt = once()
    if rc != 0 and rc != "TIMEOUT":
        rec["n_retries"] = 1
        rc, dt = once()
    rec["runtime_s"] = round(dt, 3)
    if rc == "TIMEOUT":
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
    elif rc != 0 or not out.is_file():
        rec.update(status="PERMANENT_FAIL", reason=f"rc={rc}")
    else:
        sc = parse_model1_remarks(out.read_text(errors="replace"))
        rec["CNNscore"] = sc["CNNscore"]
        rec["CNNaffinity"] = sc["CNNaffinity"]
        rec["CNN_VS"] = sc["CNN_VS"]
        rec["empirical_affinity"] = sc["minimizedAffinity"]
        if finite_float(sc["minimizedAffinity"]):
            rec["M3_empirical_higher_better"] = f"{-float(sc['minimizedAffinity']):.6f}"
        rec["status"] = "SUCCESS" if finite_float(sc["CNNscore"]) else "PERMANENT_FAIL"
        rec["reason"] = "" if rec["status"] == "SUCCESS" else "no_finite_model1_CNNscore"
    st_p.write_text(json.dumps(rec, indent=2) + "\n")
    return rec


def main() -> int:
    LIG_OUT.mkdir(parents=True, exist_ok=True)
    JOB_OUT.mkdir(parents=True, exist_ok=True)
    topo = {r["global_ligand_entity_id"]: r for r in csv.DictReader(TOPO.open(encoding="utf-8-sig"))}
    fail_eids = sorted(eid for eid, r in topo.items() if r["decision"] == "M3_INPUT_TOPOLOGY_FAIL")
    audit = {r["global_ligand_entity_id"]: r for r in csv.DictReader(AUDIT.open(encoding="utf-8-sig"))}
    prep_rows = []
    for eid in fail_eids:
        rec = prepare_ligand(eid)
        prep_rows.append(rec)
        print("prep", eid, rec["status"], rec["reason"], rec.get("n_heavy"), flush=True)
        if rec["status"] != "OK":
            print("STOP ligand prep fail", eid, flush=True)
            return 2
    with (QA / "M3_TOPOLOGY_REPAIR_LIGAND_PREP.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(prep_rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(prep_rows)

    man = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in csv.DictReader(MANIFEST.open(encoding="utf-8-sig"))}
    jobs = []
    for eid in fail_eids:
        for jid in audit[eid]["affected_physical_jobs_M3"].split(";"):
            if not jid:
                continue
            pdb, rest = jid.split("__", 1)
            assert rest == eid
            jobs.append((pdb, eid, man.get((pdb, eid), {}).get("pairs", "")))
    print(f"m3_repair jobs={len(jobs)} workers={WORKERS}", flush=True)
    done = {}
    pending = []
    for pdb, eid, pairs in jobs:
        st = JOB_OUT / f"{pdb}__{eid}__seed42" / "status.json"
        if st.is_file():
            prev = json.loads(st.read_text())
            if prev.get("status") in {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"}:
                done[(pdb, eid)] = prev
                continue
        pending.append((pdb, eid, pairs))
    if pending:
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futs = {pool.submit(run_job, pdb, eid, pairs): (pdb, eid) for pdb, eid, pairs in pending}
            for i, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                done[(rec["pdb_id"], rec["global_ligand_entity_id"])] = rec
                print(f"{len(done)}/{len(jobs)} {rec['status']} {rec['job_id']}", flush=True)

    repaired = {f"{pdb}__{eid}__seed42": rec for (pdb, eid), rec in done.items()}
    old = list(csv.DictReader(REPARSED.open(encoding="utf-8-sig")))
    out_rows = []
    for r in old:
        eid = r["global_ligand_entity_id"]
        if eid in fail_eids:
            rec = repaired.get(r["job_id"])
            if rec is None:
                row = dict(r)
                row["status"] = "INVALID_INPUT_REPRESENTATION"
                row["score_source"] = "missing_repair"
                row["official_m3"] = "NO"
                out_rows.append(row)
                continue
            row = {
                "job_id": rec["job_id"], "pdb_id": rec["pdb_id"],
                "global_ligand_entity_id": rec["global_ligand_entity_id"],
                "pairs": rec.get("pairs") or r.get("pairs", ""),
                "seed": 42, "status": rec["status"],
                "CNNscore": rec.get("CNNscore", ""),
                "CNNaffinity": rec.get("CNNaffinity", ""),
                "CNN_VS": rec.get("CNN_VS", ""),
                "empirical_affinity": rec.get("empirical_affinity", ""),
                "M3_empirical_higher_better": rec.get("M3_empirical_higher_better", ""),
                "runtime_s": rec.get("runtime_s", ""),
                "n_retries": rec.get("n_retries", ""),
                "reason": rec.get("reason", ""),
                "score_source": "topology_repaired_rigid_macrocycles",
                "official_m3": "YES" if rec.get("status") == "SUCCESS" else "NO",
                "original_invalid_input_retained_for_provenance": "YES",
            }
            out_rows.append(row)
        else:
            row = dict(r)
            row["score_source"] = "original_reparsed_model1"
            row["official_m3"] = "YES" if r.get("status") == "SUCCESS" else "NO"
            if r.get("status") == "TIMEOUT":
                row["score_source"] = "original_timeout_kept_no_rescue"
            out_rows.append(row)
    fields: list[str] = []
    seen: set[str] = set()
    for row in out_rows:
        for k in row:
            if k not in seen:
                seen.add(k)
                fields.append(k)
    with VERIFIED.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(json.dumps({
        "n_verified": len(out_rows),
        "n_repaired_jobs": len(repaired),
        "repaired_status": {s: sum(1 for r in repaired.values() if r.get("status")==s) for s in {x.get("status") for x in repaired.values()}},
        "AUROC": "NOT_COMPUTED",
        "old_master_preserved": str(OLD_MASTER),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
