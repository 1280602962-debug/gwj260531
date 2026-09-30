#!/usr/bin/env python3
"""V4.2 Phase 5: 70 cognate redocking jobs. STOP. Do not start Phase 6/7."""
from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, rdmolops
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.uniform_protocol_config import BOX_DEFINING_LIGAND, UNIQUE_14  # noqa: E402
from rerun_v4_1.phase3_prepare_ligands import graph_sig, resolve_mk  # noqa: E402

RUN_ID = "UNIFORM_RERUN_V4_2_20260921"
RUN = ROOT / "reruns" / RUN_ID
REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
QA = RUN / "13_qa"
DOCK = RUN / "04_cognate_redocking"
LIG = DOCK / "ligands"
JOBS = DOCK / "jobs"
SEEDS = (17, 29, 42, 71, 101)
VINA = Path("/home/gwj/miniconda3/bin/vina")
TIMEOUT_S = 1800  # provisional clamp upper bound; timeout_frozen=false
WORKERS = 6
NCONFS = 20
EMBED_SEED = 42
COGNATE_ID = ROOT / "qa/final_input_audit/cognate_identity_master.csv"


def assert_v42(path: Path) -> None:
    if RUN_ID not in str(path.resolve()):
        raise SystemExit(f"refusing write outside V4.2: {path}")


def load_cognate_table() -> dict[str, dict]:
    rows = list(csv.DictReader(COGNATE_ID.open(encoding="utf-8-sig", newline="")))
    return {r["pdb_id"]: r for r in rows}


def parse_pdb_atoms(text: str) -> list[dict]:
    atoms = []
    for ln in text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        atoms.append({
            "name": ln[12:16].strip(),
            "comp": ln[17:20].strip(),
            "ch": ln[21].strip(),
            "seq": int(ln[22:26]),
            "x": float(ln[30:38]),
            "y": float(ln[38:46]),
            "z": float(ln[46:54]),
        })
    return atoms


def crystal_mol(pdb: str, smiles: str):
    data = json.loads((REC / f"{pdb}_cognate_crystal.json").read_text())
    atoms = [a for a in data["atoms"] if a.get("elem") != "H"]
    block = ["COMPND    %s\n" % pdb]
    for i, a in enumerate(atoms, 1):
        name = a["name"]
        name4 = f"{name:>4s}" if len(name) < 4 else name[:4]
        block.append(
            f"HETATM{i:5d} {name4} {a['comp']:>3s} {a['auth_asym']}{int(a['auth_seq']):4d}    "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  1.00          {a['elem']:>2s}\n"
        )
    block.append("END\n")
    pdbmol = Chem.MolFromPDBBlock("".join(block), sanitize=False, removeHs=True)
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None or pdbmol is None:
        return None, "crystal_or_smiles_parse_fail"
    tmpl = Chem.RemoveHs(tmpl)
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(tmpl, pdbmol)
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return None, f"assign_bond_orders:{exc}"
    return Chem.RemoveHs(mol), "ok"


def prep_cognate(pdb: str, smiles: str, charge: int, mk: str) -> dict:
    rec = {
        "pdb_id": pdb,
        "ccd": BOX_DEFINING_LIGAND[pdb]["ccd"],
        "smiles": smiles,
        "formal_charge": charge,
        "embed_status": "",
        "ff": "",
        "n_confs": 0,
        "identity_status": "",
        "pdbqt_status": "",
        "reason": "",
    }
    sdf = LIG / f"{pdb}_cognate.sdf"
    pdbqt = LIG / f"{pdb}_cognate.pdbqt"
    assert_v42(sdf)
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        rec.update(embed_status="fail", pdbqt_status="fail", identity_status="FAIL", reason="bad_smiles")
        return rec
    if int(rdmolops.GetFormalCharge(ref)) != int(charge):
        rec.update(identity_status="FAIL", pdbqt_status="fail", reason="formal_charge_mismatch")
        return rec
    mol = Chem.AddHs(Chem.Mol(ref))
    params = AllChem.ETKDGv3()
    params.randomSeed = EMBED_SEED
    params.numThreads = 1
    params.pruneRmsThresh = -1.0
    AllChem.EmbedMultipleConfs(mol, numConfs=NCONFS, params=params)
    n = mol.GetNumConformers()
    if n < 1:
        rec.update(embed_status="fail", pdbqt_status="fail", reason="LIGAND_GEOMETRY_FAILED")
        return rec
    rec["n_confs"] = int(n)
    rec["embed_status"] = "etkdgv3_ok"
    energies = []
    used_uff = False
    for cid in range(n):
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            if props is None:
                raise ValueError("no mmff")
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            ff.Minimize()
            energies.append((ff.CalcEnergy(), cid, "MMFF94s"))
        except Exception:
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                ff.Minimize()
                energies.append((ff.CalcEnergy(), cid, "UFF"))
                used_uff = True
            except Exception:
                continue
    if not energies:
        rec.update(ff="fail", pdbqt_status="fail", reason="LIGAND_GEOMETRY_FAILED")
        return rec
    energies.sort()
    _, best_cid, ffname = energies[0]
    rec["ff"] = ffname
    rec["uff_fallback"] = int(used_uff or ffname == "UFF")
    best = Chem.Mol(mol)
    best.RemoveAllConformers()
    best.AddConformer(mol.GetConformer(best_cid), assignId=True)
    if graph_sig(Chem.RemoveHs(best)) != graph_sig(ref):
        rec.update(identity_status="FAIL", pdbqt_status="fail", reason="graph_changed")
        return rec
    rec["identity_status"] = "PASS"
    LIG.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(sdf))
    best.SetProp("_Name", f"{pdb}_{BOX_DEFINING_LIGAND[pdb]['ccd']}")
    w.write(best)
    w.close()
    proc = subprocess.run([mk, "-i", str(sdf), "-o", str(pdbqt)], capture_output=True, text=True, timeout=180)
    if proc.returncode == 0 and pdbqt.is_file() and pdbqt.stat().st_size > 0:
        rec.update(meeko_status="ok", pdbqt_status="ok", reason="ok")
    else:
        rec.update(meeko_status="fail", pdbqt_status="fail", reason=(proc.stderr or proc.stdout or "meeko")[-200:])
    return rec


def calc_pose_rmsds(out_pdbqt: Path, ref_mol) -> dict:
    if ref_mol is None:
        return {"mapping_status": "NO_CRYSTAL_MOL", "n_modes_returned": 0}
    try:
        pq = PDBQTMolecule.from_file(str(out_pdbqt))
        mols = RDKitMolCreate.from_pdbqt_mol(pq)
    except Exception as exc:
        return {"mapping_status": f"pdbqt_parse:{exc}", "n_modes_returned": 0}
    if not mols or mols[0] is None:
        return {"mapping_status": "meeko_rdkit_none", "n_modes_returned": 0}
    docked = Chem.RemoveHs(mols[0])
    ref = Chem.RemoveHs(ref_mol)
    if docked.GetNumHeavyAtoms() != ref.GetNumHeavyAtoms():
        return {
            "mapping_status": f"heavy_count_mismatch_{docked.GetNumHeavyAtoms()}_{ref.GetNumHeavyAtoms()}",
            "n_modes_returned": docked.GetNumConformers(),
        }
    rmsds = []
    for cid in range(docked.GetNumConformers()):
        try:
            rms = float(rdMolAlign.CalcRMS(docked, ref, prbId=cid, refId=0, maxMatches=10000))
            rmsds.append(round(rms, 3))
        except Exception as exc:
            return {"mapping_status": f"calcrrms:{exc}", "n_modes_returned": docked.GetNumConformers()}
    energies = list(pq._pose_data.get("free_energies") or [])
    return {
        "mapping_status": "CalcRMS_symmetry_aware_in_place",
        "n_modes_returned": len(rmsds),
        "mode1_affinity": energies[0] if energies else "",
        "top1_rmsd_A": rmsds[0] if rmsds else "",
        "min_saved_pose_rmsd_A": min(rmsds) if rmsds else "",
        "min_rmsd_mode": (1 + min(range(len(rmsds)), key=lambda i: rmsds[i])) if rmsds else "",
        "all_rmsds": rmsds,
        "all_affinities": energies,
    }


def run_vina(pdb: str, seed: int) -> dict:
    box = json.loads((BOX / f"{pdb}_box.json").read_text())
    rec_pdbqt = REC / f"{pdb}_receptor.pdbqt"
    lig_pdbqt = LIG / f"{pdb}_cognate.pdbqt"
    jobdir = JOBS / f"{pdb}_seed{seed}"
    jobdir.mkdir(parents=True, exist_ok=True)
    assert_v42(jobdir)
    out = jobdir / "out.pdbqt"
    log = jobdir / "vina.log"
    cmd = [
        str(VINA),
        "--receptor", str(rec_pdbqt),
        "--ligand", str(lig_pdbqt),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--scoring", "vina",
        "--exhaustiveness", "16",
        "--num_modes", "9",
        "--energy_range", "6",
        "--cpu", "1",
        "--seed", str(seed),
        "--out", str(out),
    ]
    rec = {
        "pdb_id": pdb,
        "seed": seed,
        "ccd": BOX_DEFINING_LIGAND[pdb]["ccd"],
        "box_instance": box.get("reference_instance"),
        "scoring": "vina",
        "exhaustiveness": 16,
        "num_modes": 9,
        "energy_range": 6,
        "cpu": 1,
        "timeout_s": TIMEOUT_S,
        "timeout_frozen": False,
        "n_retries": 0,
        "status": "",
        "runtime_s": "",
        "vina_rc": "",
        "mode1_affinity": "",
        "n_modes_returned": "",
        "top1_rmsd_A": "",
        "min_saved_pose_rmsd_A": "",
        "min_rmsd_mode": "",
        "mapping_status": "",
        "reason": "",
        "receptor_modified_from_rmsd": 0,
        "box_modified_from_rmsd": 0,
        "his_modified_from_rmsd": 0,
        "ptr_modified_from_rmsd": 0,
        "exhaustiveness_modified_from_rmsd": 0,
    }
    if not lig_pdbqt.is_file():
        rec.update(status="LIGAND_PREP_FAIL", reason="missing_cognate_pdbqt")
        return rec

    def once() -> tuple[int | str, float, str]:
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
            runtime = time.time() - t0
            (jobdir / "vina.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
            (jobdir / "vina.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
            log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
            return proc.returncode, runtime, ""
        except subprocess.TimeoutExpired:
            runtime = time.time() - t0
            (jobdir / "vina.stderr.txt").write_text("TIMEOUT\n", encoding="utf-8")
            return "TIMEOUT", runtime, "JOB_TIMEOUT"

    rc, runtime, err = once()
    if rc != 0:
        rec["n_retries"] = 1
        rc, runtime, err = once()
    rec["runtime_s"] = round(runtime, 3)
    rec["vina_rc"] = rc
    if rc == "TIMEOUT":
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
        return rec
    if rc != 0 or not out.is_file() or out.stat().st_size == 0:
        rec.update(status="FAIL", reason=err or f"vina_rc={rc}")
        return rec
    rec["status"] = "OK"
    return rec


def write_master(rows: list[dict]) -> None:
    path = QA / "redocking_master.csv"
    assert_v42(path)
    fields = [
        "pdb_id", "seed", "ccd", "box_instance", "status", "runtime_s", "vina_rc",
        "n_retries", "mode1_affinity", "n_modes_returned", "top1_rmsd_A",
        "min_saved_pose_rmsd_A", "min_rmsd_mode", "mapping_status", "reason",
        "scoring", "exhaustiveness", "num_modes", "energy_range", "cpu",
        "timeout_s", "timeout_frozen",
        "receptor_modified_from_rmsd", "box_modified_from_rmsd",
        "his_modified_from_rmsd", "ptr_modified_from_rmsd",
        "exhaustiveness_modified_from_rmsd",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["pdb_id"], int(r["seed"]))))


def write_audit(rows: list[dict], lig_rows: list[dict], gate: dict) -> dict:
    n = len(rows)
    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_to = sum(1 for r in rows if r["status"] == "TIMEOUT")
    n_fail = sum(1 for r in rows if r["status"] == "FAIL")
    n_hold = sum(1 for r in rows if r["status"] == "HOLD_REVIEW_REQUIRED")
    n_ligfail = sum(1 for r in rows if r["status"] == "LIGAND_PREP_FAIL")
    phase6 = (RUN / "05_vina").exists() and any((RUN / "05_vina").rglob("*"))
    phase7 = (RUN / "06_vina_fiveseed").exists() and any((RUN / "06_vina_fiveseed").rglob("*"))
    answers = {
        "RUN_ID": RUN_ID,
        "PHASE": 5,
        "jobs_accounted": n,
        "jobs_expected": 70,
        "OK": n_ok,
        "TIMEOUT": n_to,
        "FAIL": n_fail,
        "HOLD_REVIEW_REQUIRED": n_hold,
        "LIGAND_PREP_FAIL": n_ligfail,
        "cognate_ligands_ok": sum(1 for r in lig_rows if r.get("pdbqt_status") == "ok"),
        "any_hold_review_required": bool(gate.get("any_hold_review_required")),
        "REDOCKING_TECHNICAL_INTEGRITY": "YES" if n == 70 else "NO",
        "parameters_changed_from_RMSD": "NO",
        "PHASE_6_STARTED": "YES" if phase6 else "NO",
        "PHASE_7_STARTED": "YES" if phase7 else "NO",
        "Production_Vina": "NOT STARTED",
        "STOP_AFTER_PHASE_5": "YES",
    }
    md = ["# REDOCKING TECHNICAL INTEGRITY AUDIT", "", f"RUN_ID = {RUN_ID}", ""]
    for k, v in answers.items():
        md.append(f"- {k}: **{v}**")
    md.append("")
    md.append("Vina: scoring=vina exhaustiveness=16 num_modes=9 energy_range=6 cpu=1 seeds=17/29/42/71/101.")
    md.append("RMSD: RDKit CalcRMS, symmetry-aware, no superposition.")
    md.append("Timeout 1800 s is provisional (timeout_frozen=false). One same-parameter technical retry.")
    md.append("No receptor/box/His/PTR/exhaustiveness change from RMSD. Phase 6/7 not started.")
    md.append("")
    md.append("## Per receptor")
    by = {}
    for r in rows:
        by.setdefault(r["pdb_id"], []).append(r)
    for pdb in UNIQUE_14:
        recs = by.get(pdb, [])
        ok = sum(1 for r in recs if r["status"] == "OK")
        rms = [r.get("top1_rmsd_A") for r in recs if r.get("top1_rmsd_A") != ""]
        md.append(f"- {pdb}: {ok}/5 OK; mode1 RMSDs={rms}")
    (QA / "REDOCKING_TECHNICAL_INTEGRITY_AUDIT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (QA / "phase5_status.json").write_text(json.dumps(answers, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(answers, indent=2))
    return answers


def main() -> int:
    if not VINA.is_file():
        raise SystemExit("vina 1.2.7 not found")
    ver = subprocess.run([str(VINA), "--version"], capture_output=True, text=True).stdout
    if "1.2.7" not in ver:
        raise SystemExit(f"unexpected vina: {ver}")
    DOCK.mkdir(parents=True, exist_ok=True)
    LIG.mkdir(parents=True, exist_ok=True)
    JOBS.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    gate = json.loads((QA / "PREPARATION_EXCEPTION_COGNATE_DISTANCE.json").read_text())
    hold_pdbs = {e["pdb_id"] for e in gate["exceptions"] if e.get("gate") == "HOLD_REVIEW_REQUIRED"}

    table = load_cognate_table()
    mk = resolve_mk()
    lig_rows = []
    refs = {}
    for pdb in UNIQUE_14:
        row = table[pdb]
        smiles = row["smiles"]
        charge = int(float(row["formal_charge"] or 0))
        lig_rows.append(prep_cognate(pdb, smiles, charge, mk))
        ref, why = crystal_mol(pdb, smiles)
        refs[pdb] = (ref, why)
        print("ligand", pdb, lig_rows[-1]["pdbqt_status"], "crystal_ref", why, flush=True)
    with (QA / "cognate_ligand_prep_status.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=sorted({k for r in lig_rows for k in r}), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(lig_rows)

    jobs = [(pdb, seed) for pdb in UNIQUE_14 for seed in SEEDS]
    rows: list[dict] = []

    def work(item: tuple[str, int]) -> dict:
        pdb, seed = item
        if pdb in hold_pdbs:
            return {
                "pdb_id": pdb, "seed": seed, "ccd": BOX_DEFINING_LIGAND[pdb]["ccd"],
                "status": "HOLD_REVIEW_REQUIRED", "reason": "incomplete_residue_within_5A_of_cognate",
                "runtime_s": "", "vina_rc": "", "n_retries": 0, "mode1_affinity": "",
                "n_modes_returned": "", "top1_rmsd_A": "", "min_saved_pose_rmsd_A": "",
                "min_rmsd_mode": "", "mapping_status": "", "box_instance": "",
                "scoring": "vina", "exhaustiveness": 16, "num_modes": 9, "energy_range": 6,
                "cpu": 1, "timeout_s": TIMEOUT_S, "timeout_frozen": False,
                "receptor_modified_from_rmsd": 0, "box_modified_from_rmsd": 0,
                "his_modified_from_rmsd": 0, "ptr_modified_from_rmsd": 0,
                "exhaustiveness_modified_from_rmsd": 0,
            }
        rec = run_vina(pdb, seed)
        if rec["status"] == "OK":
            rms = calc_pose_rmsds(JOBS / f"{pdb}_seed{seed}" / "out.pdbqt", refs[pdb][0])
            rec.update({k: rms.get(k, rec.get(k, "")) for k in (
                "mapping_status", "n_modes_returned", "mode1_affinity",
                "top1_rmsd_A", "min_saved_pose_rmsd_A", "min_rmsd_mode",
            )})
            (JOBS / f"{pdb}_seed{seed}" / "rmsd.json").write_text(json.dumps(rms, indent=2) + "\n", encoding="utf-8")
        rec["crystal_ref_status"] = refs[pdb][1]
        (JOBS / f"{pdb}_seed{seed}" / "status.json").write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        return rec

    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(work, job): job for job in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            rows.append(rec)
            write_master(rows)
            print(f"{i}/70", rec["pdb_id"], rec["seed"], rec["status"], rec.get("top1_rmsd_A"), rec.get("runtime_s"), flush=True)

    write_master(rows)
    write_audit(rows, lig_rows, gate)
    # refuse Phase 6/7
    for banned in (RUN / "05_vina", RUN / "06_vina_fiveseed", RUN / "07_gnina"):
        if banned.exists() and any(banned.rglob("*")):
            print("WARNING: unexpected Phase 6/7 outputs at", banned)
    print("PHASE5_DONE PHASE6_NOT_STARTED PHASE7_NOT_STARTED")
    return 0 if len(rows) == 70 else 1


if __name__ == "__main__":
    raise SystemExit(main())
