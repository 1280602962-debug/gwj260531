#!/usr/bin/env python3
"""Phase 5-8: RTM pockets, frozen-SDF+Vina coordinate mapping, rescore, audit. No AUROC."""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdmolops

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    LIG_PDBQT,
    LIG_SDF,
    POCKET_CUTOFF,
    QA,
    REC,
    RTM_CKPT,
    RTM_PY,
    RTM_ROOT,
    RTM_SCRIPT,
    SEC,
    UNIQUE_14,
)

POCKET_DIR = SEC / "rtmscore" / "pockets"
MAPPED = SEC / "rtmscore" / "mapped_poses"
WORK = SEC / "rtmscore" / "jobs"
MANIFEST = SEC / "rescoring_job_manifest.csv"


def parse_pdb_atoms(text: str) -> list[dict]:
    atoms = []
    for ln in text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        atoms.append({
            "ln": ln,
            "name": ln[12:16].strip(),
            "comp": ln[17:20].strip(),
            "ch": ln[21],
            "seq": ln[22:26],
            "x": float(ln[30:38]),
            "y": float(ln[38:46]),
            "z": float(ln[46:54]),
            "serial": int(ln[6:11]),
        })
    return atoms


def write_pocket(pdb: str) -> dict:
    prot = parse_pdb_atoms((REC / f"{pdb}_protein_prepared.pdb").read_text(errors="replace"))
    cog = json.loads((REC / f"{pdb}_cognate_crystal.json").read_text())
    cref = [(a["x"], a["y"], a["z"]) for a in cog["atoms"] if a.get("elem") != "H"]
    keep_res = set()
    for a in prot:
        for x, y, z in cref:
            if math.dist((a["x"], a["y"], a["z"]), (x, y, z)) <= POCKET_CUTOFF:
                keep_res.add((a["ch"], a["seq"], a["comp"]))
                break
    kept = [a for a in prot if (a["ch"], a["seq"], a["comp"]) in keep_res]
    dest = POCKET_DIR / f"{pdb}_rtmscore_pocket10.pdb"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("".join(a["ln"] + ("\n" if not a["ln"].endswith("\n") else "") for a in kept) + "END\n", encoding="utf-8")
    return {"pdb_id": pdb, "n_pocket_atoms": len(kept), "n_residues": len(keep_res), "path": str(dest)}


def parse_smiles_idx(pdbqt_text: str) -> tuple[str, dict[int, int]]:
    smiles = ""
    pairs: list[int] = []
    for ln in pdbqt_text.splitlines():
        if ln.startswith("REMARK SMILES IDX"):
            pairs.extend(int(x) for x in ln.split()[3:])
        elif ln.startswith("REMARK SMILES") and "IDX" not in ln:
            smiles = ln[len("REMARK SMILES "):].strip()
    if len(pairs) % 2:
        return smiles, {}
    mapping = {}
    for i in range(0, len(pairs), 2):
        mapping[pairs[i]] = pairs[i + 1]  # smiles_idx (1-based) -> pdbqt serial
    return smiles, mapping


def pdbqt_coords_by_serial(text: str) -> dict[int, tuple[float, float, float]]:
    out = {}
    for ln in text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        serial = int(ln[6:11])
        out[serial] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return out


def map_vina_to_frozen_sdf(eid: str, vina_pose_pdbqt: str) -> tuple[Chem.Mol | None, str]:
    lig_p = LIG_PDBQT / f"{eid}.pdbqt"
    sdf_p = LIG_SDF / f"{eid}.sdf"
    if not lig_p.is_file() or not sdf_p.is_file():
        return None, "missing_frozen_ligand_asset"
    smiles, idx = parse_smiles_idx(lig_p.read_text(errors="replace"))
    if not smiles or not idx:
        return None, "ATOM_MAPPING_FAILURE_no_smiles_idx"
    sdf = Chem.MolFromMolFile(str(sdf_p), removeHs=False, sanitize=True)
    if sdf is None:
        return None, "ATOM_MAPPING_FAILURE_sdf_parse"
    heavy = Chem.Mol(rdmolops.RemoveHs(sdf))
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        return None, "ATOM_MAPPING_FAILURE_remark_smiles"
    if heavy.GetNumAtoms() != tmpl.GetNumAtoms():
        return None, "ATOM_MAPPING_FAILURE_atom_count"
    order_ok = all(
        heavy.GetAtomWithIdx(i).GetAtomicNum() == tmpl.GetAtomWithIdx(i).GetAtomicNum()
        for i in range(tmpl.GetNumAtoms())
    )
    if order_ok:
        match = tuple(range(tmpl.GetNumAtoms()))
    else:
        matches = heavy.GetSubstructMatches(tmpl, uniquify=True)
        if len(matches) != 1:
            return None, "ATOM_MAPPING_FAILURE_ambiguous_isomorphism"
        match = matches[0]
    coords = pdbqt_coords_by_serial(vina_pose_pdbqt)
    conf = Chem.Conformer(heavy.GetNumAtoms())
    for smi_idx_1, pdbqt_serial in idx.items():
        if smi_idx_1 < 1 or smi_idx_1 > tmpl.GetNumAtoms():
            return None, "ATOM_MAPPING_FAILURE_idx_range"
        if pdbqt_serial not in coords:
            return None, "ATOM_MAPPING_FAILURE_missing_vina_serial"
        sdf_idx = match[smi_idx_1 - 1]
        x, y, z = coords[pdbqt_serial]
        conf.SetAtomPosition(sdf_idx, (x, y, z))
        if heavy.GetAtomWithIdx(sdf_idx).GetAtomicNum() != tmpl.GetAtomWithIdx(smi_idx_1 - 1).GetAtomicNum():
            return None, "ATOM_MAPPING_FAILURE_element_mismatch"
    if len(idx) != tmpl.GetNumAtoms():
        return None, "ATOM_MAPPING_FAILURE_incomplete_idx"
    heavy.RemoveAllConformers()
    heavy.AddConformer(conf, assignId=True)
    return heavy, "ok"


def read_rtm_csv(out_prefix: Path) -> tuple[str, str]:
    csv_p = Path(str(out_prefix) + ".csv")
    if not csv_p.is_file():
        return "", "rtm_missing_csv"
    rows = list(csv.DictReader(csv_p.open(encoding="utf-8-sig")))
    if not rows:
        return "", "rtm_empty"
    sc = rows[0].get("score", "")
    if sc in ("", None):
        return "", "rtm_empty"
    return sc, "ok"


def run_rtm(pocket: Path, lig_sdf: Path, out_prefix: Path) -> tuple[str, str]:
    cached = read_rtm_csv(out_prefix)
    if cached[1] == "ok":
        return cached
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{RTM_ROOT}:{env.get('PYTHONPATH', '')}"
    cmd = [str(RTM_PY), str(RTM_SCRIPT), "-p", str(pocket), "-l", str(lig_sdf), "-m", str(RTM_CKPT), "-o", str(out_prefix)]
    log_p = out_prefix.parent / (out_prefix.name + ".log")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(RTM_ROOT / "example"), timeout=180)
    except subprocess.TimeoutExpired as exc:
        log_p.write_text((exc.stdout or "") + "\n" + (exc.stderr or "") + "\nJOB_TIMEOUT 180s\n", encoding="utf-8")
        return "", "rtm_timeout"
    log_p.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
    if proc.returncode != 0:
        return "", f"rtm_rc={proc.returncode}"
    return read_rtm_csv(out_prefix)


def main() -> int:
    POCKET_DIR.mkdir(parents=True, exist_ok=True)
    MAPPED.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    pockets = []
    for pdb in UNIQUE_14:
        pockets.append(write_pocket(pdb))
        print("pocket", pdb, pockets[-1]["n_residues"], flush=True)
    with (SEC / "rtmscore" / "pocket_manifest.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(pockets[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(pockets)

    jobs = list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig", newline="")))
    long_rows = []
    masters = []
    gnina_work = SEC / "gnina_rescore"
    for n, row in enumerate(jobs, 1):
        pdb, eid = row["pdb_id"], row["global_ligand_entity_id"]
        rec = {
            "pdb_id": pdb, "global_ligand_entity_id": eid, "pairs": row["pairs"],
            "job_status": "", "M2_RTMScore": "", "max_saved_RTMScore": "",
            "n_poses": 0, "reason": "", "n_retries": 0,
        }
        if row.get("rescoring_eligible") != "1":
            rec.update(job_status="MISSING_NO_VINA_POSE", reason="no_vina_pose")
            masters.append(rec)
            continue
        pose_dir = gnina_work / f"{pdb}__{eid}"
        mode_files = sorted(pose_dir.glob("vina_mode*.pdbqt"), key=lambda p: int(p.stem.replace("vina_mode", "")))
        if not mode_files:
            rec.update(job_status="TECHNICAL_FAIL", reason="missing_split_vina_poses")
            masters.append(rec)
            continue
        scores = []
        mapped_ok = True
        for mf in mode_files:
            mode = int(mf.stem.replace("vina_mode", ""))
            mol, why = map_vina_to_frozen_sdf(eid, mf.read_text(errors="replace"))
            if mol is None:
                rec.update(job_status="ATOM_MAPPING_FAILURE", reason=why)
                mapped_ok = False
                break
            sdf_out = MAPPED / f"{eid}__{pdb}__vina_mode{mode}_for_rtmscore.sdf"
            w = Chem.SDWriter(str(sdf_out))
            mol.SetProp("_Name", f"{eid}_{pdb}_m{mode}")
            w.write(mol)
            w.close()
            prefix = WORK / f"{pdb}__{eid}__mode{mode}"
            sc, st = run_rtm(POCKET_DIR / f"{pdb}_rtmscore_pocket10.pdb", sdf_out, prefix)
            if st != "ok":
                sc, st = run_rtm(POCKET_DIR / f"{pdb}_rtmscore_pocket10.pdb", sdf_out, prefix)
                rec["n_retries"] = 1
            long_rows.append({
                "pdb_id": pdb, "global_ligand_entity_id": eid, "pairs": row["pairs"],
                "vina_seed": 42, "vina_mode": mode, "is_vina_mode1": int(mode == 1),
                "RTMScore": sc, "status": st,
            })
            if st != "ok":
                rec.update(job_status="TECHNICAL_FAIL", reason=st)
                mapped_ok = False
                break
            scores.append((mode, float(sc)))
        if mapped_ok and scores:
            s1 = next(s for m, s in scores if m == 1)
            rec.update(
                job_status="SUCCESS", n_poses=len(scores),
                M2_RTMScore=f"{s1:.6f}",
                max_saved_RTMScore=f"{max(s for _, s in scores):.6f}",
            )
        masters.append(rec)
        if n % 25 == 0 or rec["job_status"] != "SUCCESS":
            print(f"rtm {n}/{len(jobs)} last={rec['job_status']} {pdb} {eid} {rec.get('reason','')}", flush=True)

    with (SEC / "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(masters[0].keys()), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(masters)
    if long_rows:
        with (SEC / "RTMSCORE_VINA_POSE_RESCORE_POSE_LONG.csv").open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(long_rows[0].keys()), extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows(long_rows)
    counts = {}
    for r in masters:
        counts[r["job_status"]] = counts.get(r["job_status"], 0) + 1
    audit = {
        "expected_keys": 1592,
        "accounted": len(masters),
        **counts,
        "AUROC": "NOT_COMPUTED",
    }
    (QA / "RESCORING_COMPLETENESS_AUDIT.md").write_text(
        "# RESCORING COMPLETENESS AUDIT\n\n" + "\n".join(f"- {k}: **{v}**" for k, v in audit.items()) + "\n",
        encoding="utf-8",
    )
    (QA / "rescoring_completeness_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))
    return 0 if len(masters) == 1592 else 1


if __name__ == "__main__":
    raise SystemExit(main())
