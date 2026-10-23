#!/usr/bin/env python3
"""Q110 V4.2 production.

Writes only under this directory. Does not replace members or edit formal assets.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
RUN = PROJECT / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
MEMBERSHIP_SRC = Path(
    "/mnt/d/CADD paper exercise/dual target docking/repo_jcim_holdout_audit/"
    "Dual_Target_Docking/extension_jcim_validation/panel_construction_audit_20261022/"
    "PIK3CA_MTOR_Q110_MEMBERSHIP.csv"
)
MEMBERSHIP_SHA = "1e85e8c25af4dbe14e756f7042993ae86ac7f89ae06f810379372d739246b550"
PAIR = "PIK3CA/mTOR"
RECEPTORS = (("A", "4L23"), ("B", "4JT6"))
SEEDS = (17, 29, 42, 71, 101)
EMBED_SEED = 20260727
VINA_TIMEOUT = 1586
GNINA_TIMEOUT = 2400
GNINA_SCORE_TIMEOUT = 300
RTM_TIMEOUT = 180
VINA_WORKERS = 6
GNINA_WORKERS = 2
VINA = Path("/home/gwj/miniconda3/bin/vina")
MEEKO = Path("/home/gwj/miniconda3/bin/mk_prepare_ligand.py")
GNINA_BIN = Path("/mnt/d/CADD paper exercise/gnina/bin/gnina")
GNINA_LIB = Path("/mnt/d/CADD paper exercise/gnina/conda_env/lib")
RTM_PY = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_SCRIPT = RTM_ROOT / "example" / "rtmscore.py"
RTM_CKPT = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"

LIG_SDF = HERE / "ligands" / "sdf"
LIG_PDBQT = HERE / "ligands" / "pdbqt"
LIG_GNINA = HERE / "ligands" / "pdbqt_gnina_ad4"
VINA_ROOT = HERE / "vina"
GNINA_RESCORE = HERE / "gnina_rescore"
RTM_ROOT_OUT = HERE / "rtmscore"
GNINA_DOCK = HERE / "gnina_independent"


def assert_inside(path: Path) -> None:
    path.resolve().relative_to(HERE.resolve())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    assert_inside(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fields:
        path.write_text("", encoding="utf-8")
        return
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_members() -> list[dict]:
    dest = HERE / "PIK3CA_MTOR_Q110_MEMBERSHIP.csv"
    if not dest.is_file():
        shutil.copy2(MEMBERSHIP_SRC, dest)
    digest = sha256(dest)
    if digest != MEMBERSHIP_SHA or sha256(MEMBERSHIP_SRC) != MEMBERSHIP_SHA:
        raise SystemExit(f"membership hash changed: {digest}")
    rows = list(csv.DictReader(dest.open(encoding="utf-8")))
    if len(rows) != 110:
        raise SystemExit(f"expected 110 members, found {len(rows)}")
    return rows


def assign_folds(rows: list[dict]) -> None:
    scaffolds = sorted({r["murcko_scaffold"] for r in rows})
    fold_of = {scaffold: index % 5 for index, scaffold in enumerate(scaffolds)}
    out = []
    for row in rows:
        fold = fold_of[row["murcko_scaffold"]]
        for arm in ("D_vs_A", "D_vs_B"):
            out.append({
                "pair": PAIR,
                "arm": arm,
                "ligand_id": row["set_id"],
                "scaffold": row["murcko_scaffold"],
                "fold_id": fold,
                "class": row["class"],
                "n_splits": 5,
            })
    write_csv(HERE / "Q110_FOLD_ASSIGNMENTS.csv", out)
    (HERE / "Q110_FOLD_ASSIGNMENT_RULE.json").write_text(json.dumps({
        "rule": "sorted unique frozen murcko_scaffold, round-robin folds 0-4",
        "same_scaffold_same_fold": True,
        "written_before_reading_q110_scores": True,
        "n_scaffolds": len(scaffolds),
        "n_molecules": len(rows),
    }, indent=2) + "\n", encoding="utf-8")


def prepare_one(row: dict) -> dict:
    eid = row["parent_inchikey"]
    smiles = row["canonical_smiles"]
    rec = {
        "set_id": row["set_id"],
        "global_ligand_entity_id": eid,
        "membership": row["membership"],
        "status": "",
        "reason": "",
        "rdkit": Chem.rdBase.rdkitVersion,
        "embed_seed": EMBED_SEED,
        "mmff_max_its": 200,
        "meeko": "0.7.1",
    }
    sdf = LIG_SDF / f"{eid}.sdf"
    pdbqt = LIG_PDBQT / f"{eid}.pdbqt"
    assert_inside(sdf)
    if row["in_pm48"] == "1":
        src_sdf = RUN / "02_ligands" / "sdf" / f"{eid}.sdf"
        src_pdbqt = RUN / "02_ligands" / "pdbqt" / f"{eid}.pdbqt"
        if not src_sdf.is_file() or not src_pdbqt.is_file():
            rec.update(status="REUSE_MISSING", reason="formal ligand asset missing")
            return rec
        LIG_SDF.mkdir(parents=True, exist_ok=True)
        LIG_PDBQT.mkdir(parents=True, exist_ok=True)
        if not sdf.is_file():
            shutil.copy2(src_sdf, sdf)
        if not pdbqt.is_file():
            shutil.copy2(src_pdbqt, pdbqt)
        rec.update(status="REUSED_FORMAL", reason="copied_not_regenerated")
        return rec
    if sdf.is_file() and pdbqt.is_file() and pdbqt.stat().st_size > 0:
        rec.update(status="SUCCESS", reason="already_prepared")
        return rec
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        rec.update(status="PERMANENT_FAIL", reason="bad_smiles")
        return rec
    got = Chem.MolToInchiKey(ref)
    if got != eid:
        rec.update(status="PERMANENT_FAIL", reason=f"inchikey_mismatch:{got}")
        return rec
    mol = Chem.AddHs(Chem.Mol(ref))
    params = AllChem.ETKDGv3()
    params.randomSeed = EMBED_SEED
    params.numThreads = 1
    if AllChem.EmbedMolecule(mol, params) < 0:
        rec.update(status="PERMANENT_FAIL", reason="LIGAND_GEOMETRY_FAILED")
        return rec
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        rec.update(status="PERMANENT_FAIL", reason="LIGAND_GEOMETRY_FAILED_no_mmff")
        return rec
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    ff.Minimize(maxIts=200)
    heavy = Chem.RemoveHs(mol)
    if Chem.MolToSmiles(heavy, isomericSmiles=False) != Chem.MolToSmiles(ref, isomericSmiles=False):
        rec.update(status="PERMANENT_FAIL", reason="graph_changed")
        return rec
    Chem.AssignStereochemistryFrom3D(heavy, replaceExistingTags=False)
    specified = {idx: tag for idx, tag in Chem.FindMolChiralCenters(ref, includeUnassigned=True) if tag != "?"}
    produced = dict(Chem.FindMolChiralCenters(heavy, includeUnassigned=True))
    if any(produced.get(idx) != tag for idx, tag in specified.items()):
        rec.update(status="PERMANENT_FAIL", reason="specified_stereo_not_preserved")
        return rec
    LIG_SDF.mkdir(parents=True, exist_ok=True)
    LIG_PDBQT.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", eid)
    writer.write(mol)
    writer.close()
    proc = subprocess.run([str(MEEKO), "-i", str(sdf), "-o", str(pdbqt)], capture_output=True, text=True, timeout=180)
    if proc.returncode != 0 or not pdbqt.is_file() or pdbqt.stat().st_size == 0:
        rec.update(status="PERMANENT_FAIL", reason=(proc.stderr or proc.stdout or "meeko")[-240:])
        return rec
    rec.update(status="SUCCESS", reason="etkdgv3_mmff200_meeko071")
    return rec


def prepare_all(rows: list[dict]) -> list[dict]:
    new = [r for r in rows if r["in_pm48"] != "1"]
    if len(new) != 101:
        raise SystemExit(f"expected 101 new ligands, found {len(new)}")
    done = []
    with ThreadPoolExecutor(max_workers=VINA_WORKERS) as pool:
        futures = [pool.submit(prepare_one, row) for row in rows]
        for fut in as_completed(futures):
            rec = fut.result()
            done.append(rec)
            print("prep", rec["set_id"], rec["status"], flush=True)
    done.sort(key=lambda r: r["set_id"])
    write_csv(HERE / "ligand_prep_master.csv", done)
    return done


def vina_env_box(pdb: str) -> dict:
    return json.loads((RUN / "03_boxes" / f"{pdb}_box.json").read_text(encoding="utf-8"))


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


def run_vina_job(row: dict, side: str, pdb: str, seed: int) -> dict:
    eid = row["parent_inchikey"]
    job_id = f"PIK3CA_mTOR__{pdb}__{eid}__seed{seed}"
    jobdir = VINA_ROOT / f"seed{seed}" / "jobs" / job_id
    assert_inside(jobdir)
    jobdir.mkdir(parents=True, exist_ok=True)
    status_p = jobdir / "status.json"
    if status_p.is_file():
        prev = json.loads(status_p.read_text(encoding="utf-8"))
        if prev.get("status") in {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"}:
            return prev
    lig = LIG_PDBQT / f"{eid}.pdbqt"
    out = jobdir / "out.pdbqt"
    box = vina_env_box(pdb)
    rec = {
        "job_id": job_id,
        "pair": PAIR,
        "pdb_id": pdb,
        "target_side": side,
        "global_ligand_entity_id": eid,
        "canonical_ligand_id": row["set_id"],
        "seed": seed,
        "scoring": "vina",
        "exhaustiveness": 16,
        "num_modes": 9,
        "energy_range": 6,
        "cpu": 1,
        "timeout_s": VINA_TIMEOUT,
        "n_retries": 0,
        "status": "",
        "runtime_s": "",
        "vina_rc": "",
        "mode1_affinity": "",
        "n_modes_returned": "",
        "reason": "",
        "source": "q110_new",
    }
    if not lig.is_file():
        rec.update(status="PERMANENT_FAIL", reason="missing_ligand_pdbqt")
        status_p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        return rec
    cmd = [
        str(VINA), "--receptor", str(RUN / "01_receptors" / f"{pdb}_receptor.pdbqt"),
        "--ligand", str(lig),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]), "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--scoring", "vina", "--exhaustiveness", "16", "--num_modes", "9", "--energy_range", "6",
        "--cpu", "1", "--seed", str(seed), "--out", str(out),
    ]

    def once() -> tuple[int | str, float]:
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=VINA_TIMEOUT)
            (jobdir / "vina.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
            return proc.returncode, time.time() - t0
        except subprocess.TimeoutExpired:
            (jobdir / "vina.log").write_text("JOB_TIMEOUT\n", encoding="utf-8")
            return "TIMEOUT", time.time() - t0

    rc, dt = once()
    if rc != 0:
        rec["n_retries"] = 1
        rc, dt = once()
    n_modes, mode1 = parse_vina_out(out)
    rec.update(runtime_s=round(dt, 3), vina_rc=rc, n_modes_returned=n_modes, mode1_affinity=mode1)
    if rc == "TIMEOUT":
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
    elif rc != 0 or not out.is_file() or out.stat().st_size == 0 or not mode1:
        rec.update(status="PERMANENT_FAIL", reason=f"vina_rc={rc}")
    else:
        rec.update(status="SUCCESS", reason="")
    status_p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    print("vina", job_id, rec["status"], rec["mode1_affinity"], flush=True)
    return rec


def reuse_formal_vina(rows: list[dict]) -> list[dict]:
    out = []
    for row in rows:
        if row["in_pm48"] != "1":
            continue
        eid = row["parent_inchikey"]
        for side, pdb in RECEPTORS:
            for seed in SEEDS:
                job_id = f"PIK3CA_mTOR__{pdb}__{eid}__seed{seed}"
                formal = RUN / "06_vina_fiveseed" / "jobs" / job_id / "out.pdbqt"
                n_modes, mode1 = parse_vina_out(formal)
                status = "SUCCESS" if mode1 else "REUSE_MISSING"
                out.append({
                    "job_id": job_id,
                    "pair": PAIR,
                    "pdb_id": pdb,
                    "target_side": side,
                    "global_ligand_entity_id": eid,
                    "canonical_ligand_id": row["set_id"],
                    "seed": seed,
                    "status": status,
                    "mode1_affinity": mode1,
                    "n_modes_returned": n_modes,
                    "reason": "reused_formal_out_pdbqt" if mode1 else "formal_pose_missing",
                    "source": "formal_v42",
                    "exhaustiveness": 16,
                    "timeout_s": VINA_TIMEOUT,
                    "n_retries": "",
                    "runtime_s": "",
                    "vina_rc": 0 if mode1 else "",
                })
    return out


def run_vina(rows: list[dict], seeds: tuple[int, ...]) -> list[dict]:
    new = [r for r in rows if r["in_pm48"] != "1"]
    jobs = [(row, side, pdb, seed) for row in new for side, pdb in RECEPTORS for seed in seeds]
    done = []
    with ThreadPoolExecutor(max_workers=VINA_WORKERS) as pool:
        futures = [pool.submit(run_vina_job, *job) for job in jobs]
        for fut in as_completed(futures):
            done.append(fut.result())
    return done


def collect_vina(rows: list[dict]) -> list[dict]:
    records = reuse_formal_vina(rows)
    for seed in SEEDS:
        root = VINA_ROOT / f"seed{seed}" / "jobs"
        if not root.is_dir():
            continue
        for status in root.glob("*/status.json"):
            records.append(json.loads(status.read_text(encoding="utf-8")))
    records.sort(key=lambda r: (int(r["seed"]), r["canonical_ligand_id"], r["pdb_id"]))
    write_csv(HERE / "q110_vina_master.csv", records)
    return records


def accept_seed42(records: list[dict]) -> dict:
    new = [r for r in records if r.get("source") == "q110_new" and int(r["seed"]) == 42]
    counts = {}
    for row in new:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    bad_params = [
        r["job_id"] for r in new
        if int(r.get("exhaustiveness", 0)) != 16 or str(r.get("seed")) != "42"
    ]
    report = {
        "expected_new_jobs": 202,
        "observed_new_jobs": len(new),
        "counts": counts,
        "bad_parameter_jobs": bad_params,
        "accepted": len(new) == 202 and not bad_params,
    }
    (HERE / "VINA_SEED42_ACCEPTANCE.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("seed42 acceptance", json.dumps(report), flush=True)
    return report


def gnina_env() -> dict:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{env.get('LD_LIBRARY_PATH', '')}"
    return env


def write_gnina_ligand(src: Path, dest: Path) -> None:
    assert_inside(dest)
    out = []
    for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            typ = line[77:].strip() if len(line) >= 78 else line.split()[-1]
            if typ == "G0":
                continue
            if typ == "CG0" and len(line) >= 78:
                line = line[:77] + f"{'C':<3}"
        out.append(line)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(out) + "\n", encoding="utf-8")


def materialize_gnina_ligands(rows: list[dict]) -> None:
    for row in rows:
        if row["in_pm48"] == "1":
            continue
        src = LIG_PDBQT / f"{row['parent_inchikey']}.pdbqt"
        if src.is_file():
            write_gnina_ligand(src, LIG_GNINA / src.name)


def main() -> int:
    rows = load_members()
    assign_folds(rows)
    (HERE / "Q110_STATUS.json").write_text(json.dumps({
        "state": "proposed_not_current",
        "membership_sha256": MEMBERSHIP_SHA,
        "rdkit": Chem.rdBase.rdkitVersion,
        "meeko": "0.7.1",
        "vina": "1.2.7",
        "note": "RDKit in this environment is the installed build. Embed parameters remain ETKDGv3 seed 20260727 and MMFF max 200.",
    }, indent=2) + "\n", encoding="utf-8")
    prepare_all(rows)
    run_vina(rows, (42,))
    records = collect_vina(rows)
    report = accept_seed42(records)
    if not report["accepted"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
