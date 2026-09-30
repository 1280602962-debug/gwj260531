#!/usr/bin/env python3
"""Eight-alternative five-seed cognate redocking QC. RMSD is diagnostic only."""
from __future__ import annotations

import csv
import json
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
from rerun_v4_1.phase3_prepare_ligands import graph_sig, resolve_mk  # noqa: E402
from rerun_v4_2.alt_receptor_config import (  # noqa: E402
    ALTS,
    JOB_TIMEOUT_S,
    RUN,
    SEEDS_REDOCK,
    VINA,
    WORKERS,
)
from rerun_v4_2.phase2_prepare_receptors import assert_v42  # noqa: E402
from rerun_v4_2.phase5_cognate_redock import calc_pose_rmsds  # noqa: E402

REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
QA = RUN / "13_qa"
DOCK = RUN / "07_alt_cognate_redocking"
LIG = DOCK / "ligands"
JOBS = DOCK / "jobs"
CCD = QA / "alt_cognate_ccd.json"
NCONFS = 20
EMBED_SEED = 42


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
    cfg = ALTS[pdb]
    rec = {
        "pdb_id": pdb, "ccd": cfg["box"]["ccd"], "smiles": smiles,
        "formal_charge": charge, "embed_status": "", "ff": "", "n_confs": 0,
        "identity_status": "", "pdbqt_status": "", "reason": "",
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
    best.SetProp("_Name", f"{pdb}_{cfg['box']['ccd']}")
    w.write(best)
    w.close()
    proc = subprocess.run([mk, "-i", str(sdf), "-o", str(pdbqt)], capture_output=True, text=True, timeout=180)
    if proc.returncode == 0 and pdbqt.is_file() and pdbqt.stat().st_size > 0:
        rec.update(meeko_status="ok", pdbqt_status="ok", reason="ok")
    else:
        rec.update(meeko_status="fail", pdbqt_status="fail", reason=(proc.stderr or proc.stdout or "meeko")[-200:])
    return rec


def run_vina(pdb: str, seed: int) -> dict:
    box = json.loads((BOX / f"{pdb}_box.json").read_text())
    rec_pdbqt = REC / f"{pdb}_receptor.pdbqt"
    lig_pdbqt = LIG / f"{pdb}_cognate.pdbqt"
    jobdir = JOBS / f"{pdb}_seed{seed}"
    jobdir.mkdir(parents=True, exist_ok=True)
    assert_v42(jobdir)
    out = jobdir / "out.pdbqt"
    cmd = [
        str(VINA), "--receptor", str(rec_pdbqt), "--ligand", str(lig_pdbqt),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--scoring", "vina", "--exhaustiveness", "16", "--num_modes", "9",
        "--energy_range", "6", "--cpu", "1", "--seed", str(seed), "--out", str(out),
    ]
    rec = {
        "pdb_id": pdb, "seed": seed, "ccd": ALTS[pdb]["box"]["ccd"],
        "box_instance": box.get("reference_instance"), "scoring": "vina",
        "exhaustiveness": 16, "num_modes": 9, "energy_range": 6, "cpu": 1,
        "timeout_s": JOB_TIMEOUT_S, "timeout_frozen": True, "n_retries": 0,
        "status": "", "runtime_s": "", "vina_rc": "", "mode1_affinity": "",
        "n_modes_returned": "", "top1_rmsd_A": "", "min_saved_pose_rmsd_A": "",
        "min_rmsd_mode": "", "mapping_status": "", "reason": "",
        "receptor_modified_from_rmsd": 0, "box_modified_from_rmsd": 0,
        "his_modified_from_rmsd": 0, "ptr_modified_from_rmsd": 0,
        "exhaustiveness_modified_from_rmsd": 0,
    }

    def once() -> tuple[int | str, float]:
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=JOB_TIMEOUT_S)
            (jobdir / "vina.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
            (jobdir / "vina.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
            (jobdir / "vina.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
            return proc.returncode, time.time() - t0
        except subprocess.TimeoutExpired:
            (jobdir / "vina.stderr.txt").write_text("JOB_TIMEOUT\n", encoding="utf-8")
            return "TIMEOUT", time.time() - t0

    rc, runtime = once()
    if rc != 0:
        rec["n_retries"] = 1
        rc, runtime = once()
    rec["runtime_s"] = round(runtime, 3)
    rec["vina_rc"] = rc
    if rc == "TIMEOUT":
        rec.update(status="TIMEOUT", reason="JOB_TIMEOUT")
        return rec
    if rc != 0 or not out.is_file() or out.stat().st_size == 0:
        rec.update(status="FAIL", reason=f"vina_rc={rc}")
        return rec
    rec["status"] = "OK"
    return rec


def main() -> int:
    if not VINA.is_file():
        raise SystemExit("vina not found")
    ver = subprocess.run([str(VINA), "--version"], capture_output=True, text=True).stdout
    if "1.2.7" not in ver:
        raise SystemExit(f"unexpected vina: {ver}")
    ccd = json.loads(CCD.read_text())
    DOCK.mkdir(parents=True, exist_ok=True)
    LIG.mkdir(parents=True, exist_ok=True)
    JOBS.mkdir(parents=True, exist_ok=True)
    mk = resolve_mk()
    lig_rows = []
    refs = {}
    for pdb in ALTS:
        info = ccd[ALTS[pdb]["box"]["ccd"]]
        smiles = info["smiles"]
        charge = int(info["charge"] or 0)
        lig_rows.append(prep_cognate(pdb, smiles, charge, mk))
        ref, why = crystal_mol(pdb, smiles)
        refs[pdb] = (ref, why)
        print("ligand", pdb, lig_rows[-1]["pdbqt_status"], "crystal_ref", why, flush=True)
    with (QA / "alt_cognate_ligand_prep_status.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=sorted({k for r in lig_rows for k in r}), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(lig_rows)

    jobs = [(pdb, seed) for pdb in ALTS for seed in SEEDS_REDOCK]
    rows = []

    def work(item):
        pdb, seed = item
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
            print(f"{i}/40", rec["pdb_id"], rec["seed"], rec["status"], rec.get("top1_rmsd_A"), rec.get("runtime_s"), flush=True)

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
    with (QA / "alt_redocking_master.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["pdb_id"], int(r["seed"]))))

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_to = sum(1 for r in rows if r["status"] == "TIMEOUT")
    n_fail = sum(1 for r in rows if r["status"] in {"FAIL", "LIGAND_PREP_FAIL"})
    gate = {
        "jobs_expected": 40,
        "jobs_accounted": len(rows),
        "OK": n_ok,
        "TIMEOUT": n_to,
        "FAIL": n_fail,
        "REDOCKING_TECHNICAL_INTEGRITY": "YES" if len(rows) == 40 and n_fail == 0 else "NO",
        "parameters_changed_from_RMSD": "NO",
        "PRODUCTION_STARTED": "NO",
    }
    (QA / "alt_redocking_technical_integrity.json").write_text(json.dumps(gate, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(gate, indent=2))
    return 0 if gate["REDOCKING_TECHNICAL_INTEGRITY"] == "YES" else 1


if __name__ == "__main__":
    raise SystemExit(main())
