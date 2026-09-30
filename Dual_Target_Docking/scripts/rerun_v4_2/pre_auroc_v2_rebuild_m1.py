#!/usr/bin/env python3
"""Rebuild unified GNINA-compatible Vina-pose representations and score_only M1/M1b.
Does not overwrite gnina_rescore/ or original masters. No AUROC.
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdmolops
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    GNINA_BIN,
    GNINA_LIB,
    LIG_PDBQT,
    LIG_SDF,
    QA,
    REC,
    SEC,
)
from rerun_v4_2.pre_auroc_v2_lib import (  # noqa: E402
    finite_float,
    formal_charge,
    is_glue_type,
    is_real_macrocycle_type,
    load_frozen_sdf,
    map_serial_to_sdf_heavy,
    parse_pdbqt_atoms,
    parse_smiles_idx,
    ring_count,
    sdf_heavy,
    special_types,
    split_models,
)

MEEKO = Path("/home/gwj/miniconda3/bin/mk_prepare_ligand.py")
MANIFEST = SEC / "rescoring_job_manifest.csv"
WORK = SEC / "gnina_rescore_representation_v2"
LONG_OUT = SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv"
MASTER_OUT = SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv"
RMS_MAX = 0.001
WORKERS = 8


def pdbqt_coords_by_serial(text: str) -> dict[int, tuple[float, float, float]]:
    out = {}
    for ln in text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        out[int(ln[6:11])] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return out


def chiral_tuple(mol: Chem.Mol) -> tuple:
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    cents = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    return tuple(sorted(cents))


def ring_bond_set(mol: Chem.Mol) -> set[tuple[int, int]]:
    bonds = set()
    ri = mol.GetRingInfo()
    for ring in ri.BondRings():
        for bidx in ring:
            b = mol.GetBondWithIdx(bidx)
            a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            bonds.add((min(a, c), max(a, c)))
    return bonds


def graph_signature(mol: Chem.Mol) -> dict:
    heavy = sdf_heavy(mol) if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) else mol
    elems = tuple(heavy.GetAtomWithIdx(i).GetAtomicNum() for i in range(heavy.GetNumAtoms()))
    charges = tuple(heavy.GetAtomWithIdx(i).GetFormalCharge() for i in range(heavy.GetNumAtoms()))
    bonds = tuple(sorted(
        (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
         max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
         int(b.GetBondTypeAsDouble() * 10),
         int(b.GetIsAromatic()))
        for b in heavy.GetBonds()
    ))
    return {
        "n_heavy": heavy.GetNumAtoms(),
        "formal_charge": formal_charge(heavy),
        "elems": elems,
        "charges": charges,
        "bonds": bonds,
        "rings": tuple(sorted(ring_bond_set(heavy))),
        "stereo": chiral_tuple(heavy),
        "n_rings": ring_count(heavy),
    }


def build_mapped_heavy(eid: str, vina_pose_text: str) -> tuple[Chem.Mol | None, str, dict]:
    sdf_p = LIG_SDF / f"{eid}.sdf"
    lig_p = LIG_PDBQT / f"{eid}.pdbqt"
    meta = {}
    if not sdf_p.is_file() or not lig_p.is_file():
        return None, "ATOM_MAPPING_FAILURE_missing_frozen_asset", meta
    sdf = load_frozen_sdf(sdf_p)
    if sdf is None:
        return None, "ATOM_MAPPING_FAILURE_sdf_parse", meta
    frozen_text = lig_p.read_text(errors="replace")
    serial_to_sdf, why = map_serial_to_sdf_heavy(frozen_text, sdf)
    if not serial_to_sdf:
        return None, why, meta
    sdf_to_serial: dict[int, int] = {}
    for serial, sdf_idx in serial_to_sdf.items():
        if sdf_idx in sdf_to_serial:
            return None, "ATOM_MAPPING_FAILURE_nonunique_sdf_index", meta
        sdf_to_serial[sdf_idx] = serial
    heavy = sdf_heavy(sdf)
    if len(sdf_to_serial) != heavy.GetNumAtoms():
        return None, "ATOM_MAPPING_FAILURE_incomplete_idx", meta
    frozen_atoms = parse_pdbqt_atoms(frozen_text)
    glue = {a["serial"] for a in frozen_atoms if a["is_glue"]}
    cgn = {a["serial"] for a in frozen_atoms if a["is_real_macro"]}
    used = set(sdf_to_serial.values())
    if used & glue:
        return None, "ATOM_MAPPING_FAILURE_used_glue", meta
    if not cgn.issubset(used):
        return None, "ATOM_MAPPING_FAILURE_dropped_CGn", meta
    coords = pdbqt_coords_by_serial(vina_pose_text)
    conf = Chem.Conformer(heavy.GetNumAtoms())
    for sdf_idx, serial in sdf_to_serial.items():
        if serial not in coords:
            return None, "ATOM_MAPPING_FAILURE_missing_vina_serial", meta
        x, y, z = coords[serial]
        conf.SetAtomPosition(sdf_idx, Point3D(float(x), float(y), float(z)))
        if heavy.GetAtomWithIdx(sdf_idx).GetAtomicNum() == 1:
            return None, "ATOM_MAPPING_FAILURE_hydrogen_in_heavy_map", meta
    before = graph_signature(sdf)
    heavy.RemoveAllConformers()
    heavy.AddConformer(conf, assignId=True)
    after = graph_signature(heavy)
    if before["n_heavy"] != after["n_heavy"] or before["bonds"] != after["bonds"] or before["charges"] != after["charges"]:
        return None, "ATOM_MAPPING_FAILURE_graph_changed", meta
    if before["rings"] != after["rings"]:
        return None, "ATOM_MAPPING_FAILURE_ring_changed", meta
    if before["formal_charge"] != after["formal_charge"]:
        return None, "ATOM_MAPPING_FAILURE_charge_changed", meta
    if before["stereo"] != after["stereo"]:
        return None, "ATOM_MAPPING_FAILURE_stereo_changed", meta
    meta = {
        "n_heavy": after["n_heavy"],
        "n_cgn_mapped": len(cgn),
        "n_glue_excluded": len(glue),
        "graph_ok": True,
        "vina_serials": sdf_to_serial,
    }
    return heavy, "ok", meta


def parse_index_map(pdbqt_text: str) -> dict[int, int]:
    """input_atom_index (1-based if Meeko) -> pdbqt serial from REMARK INDEX MAP if present."""
    pairs: list[int] = []
    for ln in pdbqt_text.splitlines():
        if "INDEX MAP" in ln and ln.startswith("REMARK"):
            pairs.extend(int(x) for x in ln.split() if x.lstrip("-").isdigit())
    if len(pairs) >= 2 and len(pairs) % 2 == 0:
        return {pairs[i]: pairs[i + 1] for i in range(0, len(pairs), 2)}
    return {}


def meeko_prepare(sdf_mol: Chem.Mol, dest_pdbqt: Path) -> tuple[str, str]:
    dest_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    # Materialize frozen-graph implicit Hs only. Do not ETKDG/MMFF or change protonation.
    mol_h = Chem.AddHs(sdf_mol, addCoords=True)
    with tempfile.TemporaryDirectory(prefix="m1rep_") as td:
        sdf_p = Path(td) / "pose.sdf"
        w = Chem.SDWriter(str(sdf_p))
        w.write(mol_h)
        w.close()
        tmp_out = Path(td) / "pose.pdbqt"
        cmd = [
            str(MEEKO), "-i", str(sdf_p), "-o", str(tmp_out),
            "--rigid_macrocycles", "--add_index_map",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        log = (proc.stdout or "") + "\n" + (proc.stderr or "")
        if proc.returncode != 0 or not tmp_out.is_file():
            return "MEEKO_FAIL", log[-2000:]
        dest_pdbqt.write_text(tmp_out.read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
        return "ok", log[-500:]


def validate_meeko_output(eid: str, mapped: Chem.Mol, pdbqt_p: Path) -> tuple[str, dict]:
    text = pdbqt_p.read_text(errors="replace")
    atoms = parse_pdbqt_atoms(text)
    spec = special_types(atoms)
    glue = [a for a in atoms if a["is_glue"] or is_glue_type(a["type"])]
    macro = [a for a in atoms if a["is_real_macro"] or is_real_macrocycle_type(a["type"])]
    if spec or glue or macro:
        return "ATOM_MAPPING_FAILURE_special_type_in_output", {"special": spec}
    smiles, idx = parse_smiles_idx(text)
    heavy = sdf_heavy(mapped)
    if not idx or len(idx) != heavy.GetNumAtoms():
        return "ATOM_MAPPING_FAILURE_meeko_idx_count", {"n_idx": len(idx), "n_heavy": heavy.GetNumAtoms()}
    coords = pdbqt_coords_by_serial(text)
    # SMILES IDX: smiles_idx -> serial. If sequential match, smiles atom i maps to heavy i.
    tmpl = Chem.MolFromSmiles(smiles) if smiles else None
    if tmpl is None or tmpl.GetNumAtoms() != heavy.GetNumAtoms():
        return "ATOM_MAPPING_FAILURE_meeko_smiles", {}
    order_ok = all(
        heavy.GetAtomWithIdx(i).GetAtomicNum() == tmpl.GetAtomWithIdx(i).GetAtomicNum()
        for i in range(heavy.GetNumAtoms())
    )
    if order_ok:
        match = tuple(range(heavy.GetNumAtoms()))
    else:
        matches = heavy.GetSubstructMatches(tmpl, uniquify=True)
        if len(matches) != 1:
            return "ATOM_MAPPING_FAILURE_meeko_isomorphism", {}
        match = matches[0]
    d2 = []
    conf = mapped.GetConformer()
    for smi_1, serial in idx.items():
        sdf_idx = match[smi_1 - 1]
        if serial not in coords:
            return "ATOM_MAPPING_FAILURE_meeko_missing_serial", {}
        x, y, z = coords[serial]
        p = conf.GetAtomPosition(sdf_idx)
        d2.append((x - p.x) ** 2 + (y - p.y) ** 2 + (z - p.z) ** 2)
        if heavy.GetAtomWithIdx(sdf_idx).GetAtomicNum() != tmpl.GetAtomWithIdx(smi_1 - 1).GetAtomicNum():
            return "ATOM_MAPPING_FAILURE_meeko_element", {}
    rms = math.sqrt(sum(d2) / len(d2)) if d2 else 999.0
    if rms > RMS_MAX:
        return "ATOM_MAPPING_FAILURE_coord_rms", {"rms": rms}
    return "ok", {"rms": rms, "n_heavy": heavy.GetNumAtoms(), "special": []}


def gnina_env() -> dict:
    e = os.environ.copy()
    e["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{e.get('LD_LIBRARY_PATH', '')}"
    return e


def score_one(receptor: Path, pose: Path, log: Path) -> dict:
    cmd = [
        str(GNINA_BIN), "--receptor", str(receptor), "--ligand", str(pose),
        "--score_only", "--no_gpu", "--scoring", "vina", "--cnn_scoring", "rescore",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=gnina_env(), timeout=300)
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    log.write_text(text, encoding="utf-8")
    out = {"rc": proc.returncode, "CNNscore": "", "CNNaffinity": "", "CNN_VS": "", "empirical_affinity": ""}
    for key, pat in (
        ("CNNscore", r"CNNscore:\s+([0-9eE.+-]+)"),
        ("CNNaffinity", r"CNNaffinity:\s+([0-9eE.+-]+)"),
        ("CNN_VS", r"CNN_?VS:\s+([0-9eE.+-]+)"),
        ("empirical_affinity", r"(?:Affinity:|REMARK minimizedAffinity)\s+([0-9eE.+-]+)"),
    ):
        import re
        m = re.search(pat, text)
        if m and m.group(1) != "--":
            out[key] = m.group(1)
    return out


def process_job(row: dict) -> dict:
    pdb, eid = row["pdb_id"], row["global_ligand_entity_id"]
    jobdir = WORK / f"{pdb}__{eid}"
    jobdir.mkdir(parents=True, exist_ok=True)
    status_p = jobdir / "status.json"
    if status_p.is_file():
        prev = json.loads(status_p.read_text())
        if prev.get("job_status") in {"SUCCESS", "MISSING_NO_VINA_POSE", "ATOM_MAPPING_FAILURE"}:
            return prev
    rec = {
        "pdb_id": pdb,
        "global_ligand_entity_id": eid,
        "pairs": row.get("pairs", ""),
        "vina_status": row.get("vina_status", ""),
        "job_status": "",
        "n_poses": 0,
        "M1_CNNscore": "",
        "M1b_CNNscore": "",
        "CNNaffinity_mode1": "",
        "CNN_VS_mode1": "",
        "empirical_affinity_mode1": "",
        "reason": "",
        "n_retries": 0,
        "coord_rms_max": "",
    }
    if row.get("rescoring_eligible") != "1" or not row.get("pose_path"):
        rec.update(job_status="MISSING_NO_VINA_POSE", reason="no_vina_pose")
        status_p.write_text(json.dumps(rec, indent=2) + "\n")
        rec["_long"] = []
        return rec
    poses = split_models(Path(row["pose_path"]).read_text(errors="replace"))
    rec["n_poses"] = len(poses)
    receptor = REC / f"{pdb}_receptor.pdbqt"
    long_rows = []
    rms_max = 0.0
    scored = []
    for i, block in enumerate(poses, 1):
        mapped, why, meta = build_mapped_heavy(eid, block)
        if mapped is None:
            rec.update(job_status="ATOM_MAPPING_FAILURE", reason=why)
            status_p.write_text(json.dumps({k: v for k, v in rec.items() if k != "_long"}, indent=2) + "\n")
            rec["_long"] = long_rows
            return rec
        sdf_p = jobdir / f"vina_mode{i}_mapped.sdf"
        w = Chem.SDWriter(str(sdf_p))
        w.write(mapped)
        w.close()
        pdbqt_p = jobdir / f"vina_mode{i}_gnina.pdbqt"
        mk_st, mk_log = meeko_prepare(mapped, pdbqt_p)
        (jobdir / f"vina_mode{i}_meeko.log").write_text(mk_log, encoding="utf-8")
        if mk_st != "ok":
            rec.update(job_status="ATOM_MAPPING_FAILURE", reason=mk_st)
            status_p.write_text(json.dumps({k: v for k, v in rec.items() if k != "_long"}, indent=2) + "\n")
            rec["_long"] = long_rows
            return rec
        val, vmeta = validate_meeko_output(eid, mapped, pdbqt_p)
        if val != "ok":
            rec.update(job_status="ATOM_MAPPING_FAILURE", reason=val)
            rec["coord_rms_max"] = str(vmeta.get("rms", ""))
            status_p.write_text(json.dumps({k: v for k, v in rec.items() if k != "_long"}, indent=2) + "\n")
            rec["_long"] = long_rows
            return rec
        rms_max = max(rms_max, float(vmeta["rms"]))
        log_p = jobdir / f"vina_mode{i}.gnina.log"
        sc = score_one(receptor, pdbqt_p, log_p)
        if sc["rc"] != 0 or not finite_float(sc.get("CNNscore", "")):
            rec["n_retries"] = 1
            sc = score_one(receptor, pdbqt_p, log_p)
        scored.append({"mode": i, **sc})
        long_rows.append({
            "pdb_id": pdb,
            "global_ligand_entity_id": eid,
            "pairs": row.get("pairs", ""),
            "vina_seed": 42,
            "vina_mode": i,
            "is_vina_mode1": int(i == 1),
            "CNNscore": sc.get("CNNscore", ""),
            "CNNaffinity": sc.get("CNNaffinity", ""),
            "CNN_VS": sc.get("CNN_VS", ""),
            "empirical_affinity": sc.get("empirical_affinity", ""),
            "gnina_rc": sc.get("rc", ""),
            "representation": "frozen_sdf_graph_plus_vina_coords_meeko_rigid_macrocycles",
        })
    rec["coord_rms_max"] = f"{rms_max:.8f}"
    ok = scored and all(finite_float(s.get("CNNscore", "")) for s in scored)
    if not ok:
        rec.update(job_status="GNINA_TECHNICAL_FAIL", reason="score_only_failed_after_retry")
    else:
        cnn = [float(s["CNNscore"]) for s in scored]
        rec.update(
            job_status="SUCCESS",
            M1_CNNscore=scored[0]["CNNscore"],
            M1b_CNNscore=f"{max(cnn):.6f}",
            CNNaffinity_mode1=scored[0].get("CNNaffinity", ""),
            CNN_VS_mode1=scored[0].get("CNN_VS", ""),
            empirical_affinity_mode1=scored[0].get("empirical_affinity", ""),
        )
        if max(cnn) + 1e-12 < cnn[0]:
            rec.update(job_status="GNINA_TECHNICAL_FAIL", reason="M1b_lt_M1_STOP")
    slim = {k: v for k, v in rec.items() if k != "_long"}
    status_p.write_text(json.dumps(slim, indent=2) + "\n")
    with (jobdir / "long.csv").open("w", encoding="utf-8", newline="") as h:
        if long_rows:
            w = csv.DictWriter(h, fieldnames=list(long_rows[0].keys()), lineterminator="\n")
            w.writeheader()
            w.writerows(long_rows)
    rec["_long"] = long_rows
    return rec


def write_tables(masters: dict, long_all: list[dict]) -> None:
    fields = [
        "pdb_id", "global_ligand_entity_id", "pairs", "vina_status", "job_status",
        "n_poses", "M1_CNNscore", "M1b_CNNscore", "CNNaffinity_mode1",
        "CNN_VS_mode1", "empirical_affinity_mode1", "reason", "n_retries", "coord_rms_max",
    ]
    rows = [masters[k] for k in sorted(masters)]
    with MASTER_OUT.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    if long_all:
        with LONG_OUT.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(long_all[0].keys()), extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows(long_all)


def main() -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig", newline="")))
    pending, done = [], []
    for r in rows:
        st = WORK / f"{r['pdb_id']}__{r['global_ligand_entity_id']}" / "status.json"
        if st.is_file():
            prev = json.loads(st.read_text())
            if prev.get("job_status") in {"SUCCESS", "MISSING_NO_VINA_POSE", "ATOM_MAPPING_FAILURE"}:
                done.append(prev)
                continue
        pending.append(r)
    print(f"m1_rebuild total={len(rows)} cached={len(done)} pending={len(pending)} workers={WORKERS}", flush=True)
    long_all = []
    masters = {(d["pdb_id"], d["global_ligand_entity_id"]): d for d in done if "pdb_id" in d}
    for d in done:
        lp = WORK / f"{d['pdb_id']}__{d['global_ligand_entity_id']}" / "long.csv"
        if lp.is_file():
            long_all.extend(csv.DictReader(lp.open(encoding="utf-8-sig")))
    if pending:
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futs = {pool.submit(process_job, r): r for r in pending}
            for i, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                long_rows = rec.pop("_long", [])
                long_all.extend(long_rows)
                masters[(rec["pdb_id"], rec["global_ligand_entity_id"])] = rec
                if rec.get("job_status") == "GNINA_TECHNICAL_FAIL" and rec.get("reason") == "M1b_lt_M1_STOP":
                    write_tables(masters, long_all)
                    print("STOP M1b < M1", rec["pdb_id"], rec["global_ligand_entity_id"], flush=True)
                    return 3
                if i % 5 == 0 or rec.get("job_status") != "SUCCESS":
                    print(f"{len(masters)}/{len(rows)} {rec.get('job_status')} {rec.get('pdb_id')} {rec.get('global_ligand_entity_id')} {rec.get('reason','')}", flush=True)
    write_tables(masters, long_all)
    counts = {}
    for r in masters.values():
        counts[r.get("job_status", "?")] = counts.get(r.get("job_status", "?"), 0) + 1
    print(json.dumps({"n": len(masters), "long": len(long_all), **counts, "AUROC": "NOT_COMPUTED"}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
