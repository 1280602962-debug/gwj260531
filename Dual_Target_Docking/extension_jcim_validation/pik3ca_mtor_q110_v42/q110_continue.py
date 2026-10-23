#!/usr/bin/env python3
"""Stages after an accepted Q110 Vina seed42 run. Writes only in this directory."""
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
import q110_run as run  # noqa: E402

RUN = run.RUN
PAIR = run.PAIR
SEEDS = run.SEEDS
CNN_RE = {
    "CNNscore": re.compile(r"CNNscore[:\s]+([0-9eE.+-]+)"),
    "CNNaffinity": re.compile(r"CNNaffinity[:\s]+([0-9eE.+-]+)"),
    "CNN_VS": re.compile(r"CNN_?VS[:\s]+([0-9eE.+-]+)"),
    "affinity": re.compile(r"(?:Affinity|VINA RESULT:|minimizedAffinity)[:\s]+([0-9eE.+-]+)"),
}


def mode1_cnnscore(text: str) -> str:
    """Mode-1 CNNscore from a pose remark or the GNINA result table, not the command line."""
    for line in text.splitlines():
        if line.startswith("REMARK CNNscore "):
            token = line.split()[-1]
            try:
                if math.isfinite(float(token)):
                    return token
            except ValueError:
                continue
            break
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0] == "1":
            try:
                value = float(parts[3])
            except ValueError:
                continue
            if math.isfinite(value) and 0.0 <= value <= 1.0:
                return parts[3]
    return ""


def parse_scores(text: str) -> dict:
    out = {"CNNscore": "", "CNNaffinity": "", "CNN_VS": "", "empirical_affinity": ""}
    for key, rx in CNN_RE.items():
        found = rx.findall(text)
        if not found:
            continue
        if key == "affinity":
            out["empirical_affinity"] = found[-1]
        elif key == "CNN_VS":
            out["CNN_VS"] = found[0]
        else:
            out[key] = found[0]
    return out


def split_models(text: str) -> list[str]:
    chunks, cur = [], []
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL") and cur:
            chunks.append("".join(cur))
            cur = [line]
        else:
            cur.append(line)
    if cur:
        chunks.append("".join(cur))
    cleaned = []
    for chunk in chunks:
        body = "".join(line for line in chunk.splitlines(keepends=True) if not line.startswith(("MODEL", "ENDMDL")))
        if "ATOM" in body or "HETATM" in body:
            cleaned.append(body)
    return cleaned


def gnina_rescore_one(row: dict, side: str, pdb: str) -> dict:
    eid = row["parent_inchikey"]
    jobdir = run.GNINA_RESCORE / f"{pdb}__{eid}"
    run.assert_inside(jobdir)
    status_p = jobdir / "status.json"
    if status_p.is_file():
        prev = json.loads(status_p.read_text(encoding="utf-8"))
        if prev.get("job_status") in {"SUCCESS", "MISSING_NO_VINA_POSE", "TECHNICAL_FAIL"}:
            return prev
    vina_status = run.VINA_ROOT / "seed42" / "jobs" / f"PIK3CA_mTOR__{pdb}__{eid}__seed42" / "status.json"
    pose = vina_status.parent / "out.pdbqt"
    rec = {
        "pdb_id": pdb, "target_side": side, "global_ligand_entity_id": eid,
        "canonical_ligand_id": row["set_id"], "job_status": "", "n_poses": 0,
        "M1_CNNscore": "", "M1b_CNNscore": "", "reason": "",
    }
    if not vina_status.is_file() or json.loads(vina_status.read_text())["status"] != "SUCCESS" or not pose.is_file():
        rec.update(job_status="MISSING_NO_VINA_POSE", reason="no_seed42_pose")
        jobdir.mkdir(parents=True, exist_ok=True)
        status_p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        return rec
    poses = split_models(pose.read_text(errors="replace"))
    rec["n_poses"] = len(poses)
    jobdir.mkdir(parents=True, exist_ok=True)
    receptor = RUN / "01_receptors" / f"{pdb}_receptor.pdbqt"
    scored = []
    for i, block in enumerate(poses, 1):
        pose_p = jobdir / f"vina_mode{i}.pdbqt"
        pose_p.write_text(block, encoding="utf-8")
        cmd = [str(run.GNINA_BIN), "--receptor", str(receptor), "--ligand", str(pose_p),
               "--score_only", "--no_gpu", "--scoring", "vina", "--cnn_scoring", "rescore"]
        proc = subprocess.run(cmd, capture_output=True, text=True, env=run.gnina_env(), timeout=run.GNINA_SCORE_TIMEOUT)
        (jobdir / f"vina_mode{i}.gnina.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
        sc = parse_scores((proc.stdout or "") + "\n" + (proc.stderr or ""))
        sc["rc"] = proc.returncode
        scored.append(sc)
    if any(not s.get("CNNscore") for s in scored):
        rec.update(job_status="TECHNICAL_FAIL", reason="parse_or_rc")
    else:
        vals = [float(s["CNNscore"]) for s in scored]
        rec.update(job_status="SUCCESS", M1_CNNscore=scored[0]["CNNscore"], M1b_CNNscore=f"{max(vals):.6f}")
    status_p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    print("gnina_rescore", pdb, row["set_id"], rec["job_status"], flush=True)
    return rec


def copy_pockets() -> None:
    dest = run.RTM_ROOT_OUT / "pockets"
    dest.mkdir(parents=True, exist_ok=True)
    src = RUN / "09_secondary_scoring" / "rtmscore" / "pockets"
    for pdb in ("4L23", "4JT6"):
        target = dest / f"{pdb}_rtmscore_pocket10.pdb"
        run.assert_inside(target)
        if not target.is_file():
            shutil.copy2(src / target.name, target)


def map_pose(eid: str, pose_text: str):
    from rdkit import Chem
    from rdkit.Chem import rdmolops
    from rerun_v4_2.ablation_phase578_rtmscore import parse_smiles_idx, pdbqt_coords_by_serial
    lig_p = run.LIG_PDBQT / f"{eid}.pdbqt"
    sdf_p = run.LIG_SDF / f"{eid}.sdf"
    smiles, idx = parse_smiles_idx(lig_p.read_text(errors="replace"))
    if not smiles or not idx:
        return None, "ATOM_MAPPING_FAILURE_no_smiles_idx"
    sdf = Chem.MolFromMolFile(str(sdf_p), removeHs=False, sanitize=True)
    if sdf is None:
        return None, "ATOM_MAPPING_FAILURE_sdf_parse"
    heavy = Chem.Mol(rdmolops.RemoveHs(sdf))
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None or heavy.GetNumAtoms() != tmpl.GetNumAtoms():
        return None, "ATOM_MAPPING_FAILURE_atom_count"
    order_ok = all(heavy.GetAtomWithIdx(i).GetAtomicNum() == tmpl.GetAtomWithIdx(i).GetAtomicNum() for i in range(tmpl.GetNumAtoms()))
    match = tuple(range(tmpl.GetNumAtoms())) if order_ok else None
    if match is None:
        matches = heavy.GetSubstructMatches(tmpl, uniquify=True)
        if len(matches) != 1:
            return None, "ATOM_MAPPING_FAILURE_ambiguous_isomorphism"
        match = matches[0]
    coords = pdbqt_coords_by_serial(pose_text)
    conf = Chem.Conformer(heavy.GetNumAtoms())
    for smi_idx_1, serial in idx.items():
        if serial not in coords:
            return None, "ATOM_MAPPING_FAILURE_missing_vina_serial"
        x, y, z = coords[serial]
        conf.SetAtomPosition(match[smi_idx_1 - 1], (x, y, z))
    if len(idx) != tmpl.GetNumAtoms():
        return None, "ATOM_MAPPING_FAILURE_incomplete_idx"
    heavy.RemoveAllConformers()
    heavy.AddConformer(conf, assignId=True)
    return heavy, "ok"


def rtm_one(row: dict, side: str, pdb: str) -> dict:
    from rdkit import Chem
    eid = row["parent_inchikey"]
    work = run.RTM_ROOT_OUT / "jobs"
    prefix = work / f"{pdb}__{eid}__mode1"
    run.assert_inside(prefix)
    cached = prefix.with_suffix(".csv")
    rec = {"pdb_id": pdb, "target_side": side, "global_ligand_entity_id": eid,
           "canonical_ligand_id": row["set_id"], "job_status": "", "M2_RTMScore": "", "reason": ""}
    if cached.is_file():
        rows = list(csv.DictReader(cached.open()))
        if rows and rows[0].get("score"):
            rec.update(job_status="SUCCESS", M2_RTMScore=rows[0]["score"])
            return rec
    pose = run.GNINA_RESCORE / f"{pdb}__{eid}" / "vina_mode1.pdbqt"
    if not pose.is_file():
        rec.update(job_status="MISSING_NO_VINA_POSE", reason="no_mode1")
        return rec
    mol, why = map_pose(eid, pose.read_text(errors="replace"))
    if mol is None:
        rec.update(job_status="ATOM_MAPPING_FAILURE", reason=why)
        return rec
    sdf_out = run.RTM_ROOT_OUT / "mapped" / f"{eid}__{pdb}__mode1.sdf"
    run.assert_inside(sdf_out)
    sdf_out.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(sdf_out))
    writer.write(mol)
    writer.close()
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{run.RTM_ROOT}:{env.get('PYTHONPATH', '')}"
    cmd = [str(run.RTM_PY), str(run.RTM_SCRIPT), "-p", str(run.RTM_ROOT_OUT / "pockets" / f"{pdb}_rtmscore_pocket10.pdb"),
           "-l", str(sdf_out), "-m", str(run.RTM_CKPT), "-o", str(prefix)]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(run.RTM_ROOT / "example"), timeout=run.RTM_TIMEOUT)
    (work / f"{pdb}__{eid}__mode1.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
    if proc.returncode == 0 and cached.is_file():
        score = list(csv.DictReader(cached.open()))[0].get("score", "")
        rec.update(job_status="SUCCESS" if score else "TECHNICAL_FAIL", M2_RTMScore=score, reason="" if score else "empty_score")
    else:
        rec.update(job_status="TECHNICAL_FAIL", reason=f"rtm_rc={proc.returncode}")
    print("rtm", pdb, row["set_id"], rec["job_status"], flush=True)
    return rec


def gnina_dock_one(row: dict, side: str, pdb: str) -> dict:
    eid = row["parent_inchikey"]
    jobdir = run.GNINA_DOCK / "jobs" / f"{pdb}__{eid}__seed42"
    run.assert_inside(jobdir)
    status_p = jobdir / "status.json"
    if status_p.is_file():
        prev = json.loads(status_p.read_text(encoding="utf-8"))
        if prev.get("status") in {"SUCCESS", "TIMEOUT", "FAIL"}:
            return prev
    lig = run.LIG_GNINA / f"{eid}.pdbqt"
    out = jobdir / "out.pdbqt"
    box = run.vina_env_box(pdb)
    jobdir.mkdir(parents=True, exist_ok=True)
    rec = {"pdb_id": pdb, "target_side": side, "global_ligand_entity_id": eid,
           "canonical_ligand_id": row["set_id"], "seed": 42, "status": "", "CNNscore": "", "reason": ""}
    if not lig.is_file():
        rec.update(status="FAIL", reason="missing_gnina_ligand")
        status_p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        return rec
    cmd = [str(run.GNINA_BIN), "--receptor", str(RUN / "01_receptors" / f"{pdb}_receptor.pdbqt"), "--ligand", str(lig),
           "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]), "--center_z", str(box["center_z"]),
           "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
           "--no_gpu", "--scoring", "vina", "--cnn_scoring", "rescore", "--pose_sort_order", "CNNscore",
           "--seed", "42", "--exhaustiveness", "16", "--num_modes", "9", "--cpu", "1", "--out", str(out)]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=run.gnina_env(), timeout=run.GNINA_TIMEOUT)
        rc, timed = proc.returncode, False
        (jobdir / "gnina.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
        text = (proc.stdout or "") + "\n" + (proc.stderr or "")
        if out.is_file():
            text += "\n" + out.read_text(errors="replace")
        scores = parse_scores(text)
        cnn = mode1_cnnscore(text)
    except subprocess.TimeoutExpired:
        rc, timed, scores, cnn = "TIMEOUT", True, {}, ""
        (jobdir / "gnina.log").write_text("JOB_TIMEOUT\n", encoding="utf-8")
    rec["runtime_s"] = round(time.time() - t0, 3)
    if timed:
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
    elif rc == 0 and cnn:
        rec.update(status="SUCCESS", CNNscore=cnn)
    else:
        rec.update(status="FAIL", reason=f"gnina_rc={rc}")
    status_p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    print("gnina_dock", pdb, row["set_id"], rec["status"], flush=True)
    return rec


def formal_lookup() -> dict:
    base = RUN
    gnina = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in csv.DictReader((base / "09_secondary_scoring" / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv").open())}
    rtm = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in csv.DictReader((base / "09_secondary_scoring" / "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv").open())}
    dock = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in csv.DictReader((base / "10_gnina_independent_docking" / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv").open())}
    return {"gnina": gnina, "rtm": rtm, "dock": dock}


def load_json_dir(root: Path) -> dict:
    out = {}
    if not root.is_dir():
        return out
    for path in root.glob("*/status.json"):
        rec = json.loads(path.read_text(encoding="utf-8"))
        out[(rec.get("pdb_id"), rec.get("global_ligand_entity_id"))] = rec
    for path in root.glob("*/*/status.json"):
        rec = json.loads(path.read_text(encoding="utf-8"))
        out[(rec.get("pdb_id"), rec.get("global_ligand_entity_id"))] = rec
    return out


def build_population(rows: list[dict]) -> list[dict]:
    vina = {}
    for seed_dir in (run.VINA_ROOT / "seed42" / "jobs").glob("*/status.json"):
        rec = json.loads(seed_dir.read_text(encoding="utf-8"))
        vina[(rec["pdb_id"], rec["global_ligand_entity_id"])] = rec
    reused = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in run.reuse_formal_vina(rows) if int(r["seed"]) == 42}
    gnina_new = load_json_dir(run.GNINA_RESCORE)
    rtm_rows = {}
    for path in (run.RTM_ROOT_OUT / "jobs").glob("*__mode1.csv"):
        pdb, eid, _ = path.name.split("__")
        score = list(csv.DictReader(path.open()))[0].get("score", "")
        rtm_rows[(pdb, eid)] = score
    dock_new = load_json_dir(run.GNINA_DOCK / "jobs")
    formal = formal_lookup()
    population = []
    for row in rows:
        eid = row["parent_inchikey"]
        item = {
            "pair": PAIR, "global_ligand_entity_id": eid, "canonical_ligand_id": row["set_id"],
            "panel_id": row["set_id"], "class": row["class"], "activity_eligible": "1",
            "analysis_set": "q110_primary", "parent_independent": "1",
            "pdb_A": "4L23", "pdb_B": "4JT6", "level": "pair_observation",
        }
        for side, pdb in run.RECEPTORS:
            if row["in_pm48"] == "1":
                v = reused.get((pdb, eid), {})
                g = formal["gnina"].get((pdb, eid), {})
                t = formal["rtm"].get((pdb, eid), {})
                d = formal["dock"].get((pdb, eid), {})
                m0_ok = v.get("status") == "SUCCESS" and v.get("mode1_affinity") not in ("", None)
                m1_ok = g.get("job_status") == "SUCCESS" and g.get("M1_CNNscore") not in ("", None)
                m2_ok = t.get("job_status") == "SUCCESS" and t.get("M2_RTMScore") not in ("", None)
                m3_ok = d.get("status") == "SUCCESS" and d.get("CNNscore") not in ("", None)
                m0 = str(-float(v["mode1_affinity"])) if m0_ok else ""
                m1 = g.get("M1_CNNscore", "") if m1_ok else ""
                m1b = g.get("M1b_CNNscore", "") if m1_ok else ""
                m2 = t.get("M2_RTMScore", "") if m2_ok else ""
                m3 = d.get("CNNscore", "") if m3_ok else ""
            else:
                v = vina.get((pdb, eid), {})
                g = gnina_new.get((pdb, eid), {})
                m0_ok = v.get("status") == "SUCCESS" and v.get("mode1_affinity") not in ("", None)
                m1_ok = g.get("job_status") == "SUCCESS" and g.get("M1_CNNscore") not in ("", None)
                m2 = rtm_rows.get((pdb, eid), "")
                m2_ok = m2 not in ("", None)
                d = dock_new.get((pdb, eid), {})
                m3_ok = d.get("status") == "SUCCESS" and d.get("CNNscore") not in ("", None)
                m0 = str(-float(v["mode1_affinity"])) if m0_ok else ""
                m1 = g.get("M1_CNNscore", "") if m1_ok else ""
                m1b = g.get("M1b_CNNscore", "") if m1_ok else ""
                m3 = d.get("CNNscore", "") if m3_ok else ""
            item[f"M0_{side}_valid"] = "1" if m0_ok else "0"
            item[f"M0_score_{side}"] = m0
            item[f"M0_{side}_missing"] = "" if m0_ok else "missing"
            item[f"M1_{side}_valid"] = "1" if m1_ok else "0"
            item[f"M1_score_{side}"] = m1
            item[f"M1_{side}_missing"] = "" if m1_ok else "missing"
            item[f"M1b_{side}_valid"] = "1" if m1_ok and m1b not in ("", None) else "0"
            item[f"M1b_score_{side}"] = m1b if m1_ok else ""
            item[f"M1b_{side}_missing"] = "" if item[f"M1b_{side}_valid"] == "1" else "missing"
            item[f"M2_{side}_valid"] = "1" if m2_ok else "0"
            item[f"M2_score_{side}"] = m2 if m2_ok else ""
            item[f"M2_{side}_missing"] = "" if m2_ok else "missing"
            item[f"M3_{side}_valid"] = "1" if m3_ok else "0"
            item[f"M3_score_{side}"] = m3 if m3_ok else ""
            item[f"M3_{side}_missing"] = "" if m3_ok else "missing"
        population.append(item)
    run.write_csv(HERE / "Q110_POPULATION.csv", population)
    return population


def run_stats(rows: list[dict], population: list[dict]) -> None:
    from q110_stats import primary_rows
    from jcim_stage2_compute import module_b, module_c, module_d, module_d1, module_e, module_h, module_i
    metrics, deltas = primary_rows(population)
    run.write_csv(HERE / "Q110_PRIMARY_METRICS.csv", metrics)
    run.write_csv(HERE / "Q110_PRIMARY_DELTAS.csv", deltas)
    vina_rows = list(csv.DictReader((HERE / "q110_vina_master.csv").open()))
    tables = {
        "population": population,
        "fiveseed": vina_rows,
        "m0_master": [
            {"pair": PAIR, "pdb_id": "4L23", "target_side": "A"},
            {"pair": PAIR, "pdb_id": "4JT6", "target_side": "B"},
        ],
        "folds": list(csv.DictReader((HERE / "Q110_FOLD_ASSIGNMENTS.csv").open())),
        "mapping": [{"pair": PAIR, "canonical_ligand_id": r["set_id"], "canonical_smiles": r["canonical_smiles"], "parent_independent": "1"} for r in rows],
        "alias": [],
    }
    out = HERE / "q110_stage2"
    out.mkdir(exist_ok=True)
    run.write_csv(out / "module_B.csv", module_b(population))
    run.write_csv(out / "module_C.csv", module_c(population))
    units, oof = module_d(tables)
    run.write_csv(out / "module_D_units.csv", units)
    run.write_csv(out / "module_D_oof.csv", oof)
    run.write_csv(out / "module_D1.csv", module_d1(tables))
    e_seed, e_sum = module_e(tables)
    run.write_csv(out / "module_E_seed.csv", e_seed)
    run.write_csv(out / "module_E_summary.csv", e_sum)
    run.write_csv(out / "module_H.csv", module_h(population))
    desc, ligs, cmps = module_i(population, deltas)
    run.write_csv(out / "module_I_desc.csv", desc)
    run.write_csv(out / "module_I_lig.csv", ligs)
    run.write_csv(out / "module_I_cmp.csv", cmps)
    strict_path = Path("/mnt/d/CADD paper exercise/dual target docking/repo_jcim_holdout_audit/Dual_Target_Docking/extension_jcim_validation/panel_construction_audit_20261022/PIK3CA_MTOR_Q110_STRICT_STATUS.csv")
    strict = {r["set_id"]: r for r in csv.DictReader(strict_path.open())}
    subset = []
    for row in population:
        src = strict[row["canonical_ligand_id"]]
        if src["strict_same_class"] == "1":
            copied = dict(row)
            copied["class"] = src["strict_class"]
            subset.append(copied)
    s_metrics, _ = primary_rows(subset)
    run.write_csv(HERE / "Q110_STRICT_SUBSET_METRICS.csv", s_metrics)
    (HERE / "Q110_STRICT_SUBSET_NOTE.json").write_text(json.dumps({
        "role": "supporting_information",
        "replaces_primary": False,
        "n": len(subset),
        "classes": {k: sum(1 for r in subset if r["class"] == k) for k in ("dual", "A_only", "B_only", "neither")},
    }, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    acceptance = json.loads((HERE / "VINA_SEED42_ACCEPTANCE.json").read_text(encoding="utf-8"))
    if not acceptance.get("accepted"):
        raise SystemExit("seed42 acceptance is not passed")
    rows = run.load_members()
    new = [r for r in rows if r["in_pm48"] != "1"]
    print("other vina seeds", flush=True)
    run.run_vina(rows, (17, 29, 71, 101))
    run.collect_vina(rows)
    print("gnina rescore", flush=True)
    with ThreadPoolExecutor(max_workers=run.GNINA_WORKERS) as pool:
        futures = [pool.submit(gnina_rescore_one, row, side, pdb) for row in new for side, pdb in run.RECEPTORS]
        for fut in as_completed(futures):
            fut.result()
    copy_pockets()
    print("rtmscore mode1", flush=True)
    for row in new:
        for side, pdb in run.RECEPTORS:
            rtm_one(row, side, pdb)
    run.materialize_gnina_ligands(rows)
    print("gnina independent", flush=True)
    with ThreadPoolExecutor(max_workers=run.GNINA_WORKERS) as pool:
        futures = [pool.submit(gnina_dock_one, row, side, pdb) for row in new for side, pdb in run.RECEPTORS]
        for fut in as_completed(futures):
            fut.result()
    population = build_population(rows)
    run_stats(rows, population)
    (HERE / "Q110_STATUS.json").write_text(json.dumps({
        "state": "accepted_current_primary",
        "membership_sha256": run.MEMBERSHIP_SHA,
        "scope": "PIK3CA/mTOR primary population replaced by frozen Q110",
        "other_seven_pairs": "unchanged formal V4.2",
        "alt_receptor": "not_authorized_this_round",
        "baseline_sha": "7777e42e8f7f204d4ba275c7f1cbefc5f8f6bd5b",
    }, indent=2) + "\n", encoding="utf-8")
    print("Q110_CONTINUE_DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
