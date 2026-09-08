#!/usr/bin/env python3
"""Optional D2 native self-docks: dotinurad@9B1G, verinurad@9IRY, lesinurad@9B1H.

These numbers are not pooled with the D1 9DK matrix.
TD-3@9DKC is the D1 diagonal and is not rerun here.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from parse_c1_sdf_readouts import load_poses, pose_rmsd  # noqa: E402
from prepare_receptor_vina import prepare_receptor  # noqa: E402
from run_urat1_capability_dock import find_gnina, gnina_version, run_gnina  # noqa: E402
from urat1_capability_lib import (  # noqa: E402
    LIGAND_PDBQT,
    NATIVE_LIKE,
    SEEDS,
    classify_search_selection,
    cnnscore_of,
    success_at_k,
)

SKIP_RES = {
    "HOH",
    "WAT",
    "DOD",
    "NAG",
    "BMA",
    "MAN",
    "GLC",
    "BGC",
    "SO4",
    "PO4",
    "GOL",
    "EDO",
    "PEG",
    "NA",
    "CL",
    "ZN",
    "CA",
    "MG",
    "K",
}
D2_JOBS = (
    {"ligand": "dotinurad", "pdb": "9B1G"},
    {"ligand": "verinurad", "pdb": "9IRY"},
    {"ligand": "lesinurad", "pdb": "9B1H"},
)
D2_OUT = PROJECT_ROOT / "data/campaigns/c5/06_capability_dock/d2"
PDB_DIR = PROJECT_ROOT / "data/structures/pdb"


def fetch_cif(pdb_id: str) -> Path:
    dest = PDB_DIR / f"{pdb_id}.cif"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"FETCH {url}", flush=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def largest_organic_het(st: gemmi.Structure) -> tuple[str, str, int]:
    best = None
    for chain in st[0]:
        for res in chain:
            name = res.name.strip()
            if name in SKIP_RES:
                continue
            if res.het_flag != "H":
                continue
            n = sum(1 for atom in res if not atom.is_hydrogen())
            if best is None or n > best[0]:
                best = (n, chain.name, name, int(res.seqid.num))
    if best is None:
        raise ValueError("no organic HETATM ligand found")
    return best[1], best[2], best[3]


def ligand_heavy_xyz(st: gemmi.Structure, chain: str, resn: str, seqid: int) -> np.ndarray:
    xyz = []
    for res in st[0][chain]:
        if res.name.strip() != resn or int(res.seqid.num) != seqid:
            continue
        for atom in res:
            if atom.is_hydrogen():
                continue
            xyz.append([atom.pos.x, atom.pos.y, atom.pos.z])
    if not xyz:
        raise ValueError(f"no heavy atoms for {resn} {chain}{seqid}")
    return np.array(xyz, dtype=float)


def write_ligand_sdf(st: gemmi.Structure, chain: str, resn: str, seqid: int, dest: Path) -> Path:
    model = st[0]
    pdb_lines = []
    serial = 1
    for res in model[chain]:
        if res.name.strip() != resn or int(res.seqid.num) != seqid:
            continue
        for atom in res:
            if atom.is_hydrogen():
                continue
            pdb_lines.append(
                f"HETATM{serial:5d} {atom.name:^4s}{resn:>3s} "
                f"{chain}{seqid:4d}    "
                f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                f"  1.00  0.00           {atom.element.name:>2s}\n"
            )
            serial += 1
    pdb_lines.append("END\n")
    dest.parent.mkdir(parents=True, exist_ok=True)
    pdb_path = dest.with_suffix(".pdb")
    pdb_path.write_text("".join(pdb_lines))
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError(f"RDKit failed to read {pdb_path}")
    mol.UpdatePropertyCache(strict=False)
    writer = Chem.SDWriter(str(dest))
    writer.write(mol)
    writer.close()
    return dest


def polymer_chain(st: gemmi.Structure) -> str:
    counts: dict[str, int] = {}
    for chain in st[0]:
        n = sum(1 for res in chain if res.het_flag != "H")
        if n:
            counts[chain.name] = n
    if not counts:
        return "A"
    return max(counts, key=counts.get)


def evaluate_selfdock(ligand: str, pdb_id: str, seed: int, sdf: Path, ref: Path) -> dict:
    row = {
        "layer": "d2",
        "ligand": ligand,
        "pdb": pdb_id,
        "seed": seed,
        "self_dock": True,
        "sdf": str(sdf.relative_to(PROJECT_ROOT)) if sdf.exists() else str(sdf),
        "status": "ok" if sdf.exists() and sdf.stat().st_size > 0 else "missing",
        "native_like_cutoff_A": NATIVE_LIKE,
    }
    if row["status"] != "ok":
        row["capability_class"] = "unevaluable"
        return row
    poses = load_poses(sdf)
    refs = load_poses(ref)
    if not poses or not refs:
        row["status"] = "no_poses"
        row["capability_class"] = "unevaluable"
        return row
    scored = []
    for i, mol in enumerate(poses):
        try:
            scored.append((cnnscore_of(mol), i, mol))
        except Exception:
            continue
    scored.sort(key=lambda x: x[0], reverse=True)
    rmsds = []
    for _score, _idx, mol in scored:
        try:
            rmsds.append(float(pose_rmsd(mol, refs[0])))
        except Exception as exc:
            rmsds.append(float("nan"))
            row.setdefault("rmsd_errors", []).append(str(exc))
    finite = [(rank, r) for rank, r in enumerate(rmsds, start=1) if r == r]
    if not finite:
        row["capability_class"] = "unevaluable"
        return row
    best_rank, best = min(finite, key=lambda x: x[1])
    top1 = rmsds[0]
    row.update(
        {
            "n_poses": len(poses),
            "CNNscore_top1": scored[0][0],
            "top1_pose_rmsd_A": float(top1),
            "best_of_9_pose_rmsd_A": float(best),
            "best_rmsd_cnnscore_rank": int(best_rank),
            "success_at_1": bool(success_at_k(rmsds, 1)),
            "success_at_3": bool(success_at_k(rmsds, 3)),
            "success_at_5": bool(success_at_k(rmsds, 5)),
            "success_at_9": bool(success_at_k(rmsds, 9)),
            "capability_class": classify_search_selection(float(top1), float(best)),
        }
    )
    return row


def prepare_d2(fetch: bool) -> list[dict]:
    rows = []
    for job in D2_JOBS:
        pdb_id = job["pdb"]
        cif = PDB_DIR / f"{pdb_id}.cif"
        if fetch or not cif.exists():
            cif = fetch_cif(pdb_id)
        if not cif.exists():
            raise SystemExit(f"missing {cif}; rerun with --fetch")
        st = gemmi.read_structure(str(cif))
        chain, resn, seqid = largest_organic_het(st)
        xyz = ligand_heavy_xyz(st, chain, resn, seqid)
        center = xyz.mean(axis=0).tolist()
        out_dir = D2_OUT / pdb_id.lower()
        ref = write_ligand_sdf(st, chain, resn, seqid, out_dir / f"{job['ligand']}_crystal_ref.sdf")
        receptor = PROJECT_ROOT / "data/structures/prepared" / f"{pdb_id}_receptor.pdbqt"
        prep_chain = polymer_chain(st)
        prepare_receptor(cif, receptor, chain_id=prep_chain, remove_waters=True, ph=7.4)
        meta = {
            "ligand": job["ligand"],
            "pdb": pdb_id,
            "cif": str(cif.relative_to(PROJECT_ROOT)),
            "ligand_chain": chain,
            "ligand_ccd": resn,
            "ligand_seqid": seqid,
            "n_heavy": int(len(xyz)),
            "center": [round(float(x), 3) for x in center],
            "size": [22, 22, 22],
            "receptor": str(receptor.relative_to(PROJECT_ROOT)),
            "prep_chain": prep_chain,
            "crystal_ref": str(ref.relative_to(PROJECT_ROOT)),
            "pdbqt": str(LIGAND_PDBQT[job["ligand"]].relative_to(PROJECT_ROOT)),
        }
        (out_dir / "prep.json").write_text(json.dumps(meta, indent=2) + "\n")
        rows.append(meta)
        print(json.dumps(meta), flush=True)
    (D2_OUT / "d2_prep.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--fetch", action="store_true", help="Download missing CIFs from RCSB")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--score-only", action="store_true")
    ap.add_argument("--gnina", default=None)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--no-gpu", action="store_true")
    args = ap.parse_args()

    D2_OUT.mkdir(parents=True, exist_ok=True)
    if args.prepare_only or args.dry_run:
        if args.dry_run and not all((PDB_DIR / f"{j['pdb']}.cif").exists() for j in D2_JOBS):
            print("D2 dry-run without local CIFs: 9 new self-dock jobs (9B1G, 9IRY, 9B1H) × seeds 42/43/44")
            return
        prepare_d2(fetch=args.fetch)
        if args.prepare_only or args.dry_run:
            return

    prep_path = D2_OUT / "d2_prep.json"
    if not prep_path.exists():
        prepare_d2(fetch=args.fetch)
    preps = {row["ligand"]: row for row in json.loads(prep_path.read_text())}

    jobs = [(job["ligand"], seed) for job in D2_JOBS for seed in SEEDS]
    print(f"D2 new GNINA jobs: {len(jobs)}", flush=True)
    if args.score_only:
        binary = None
    else:
        binary = find_gnina(args.gnina)
        (D2_OUT / "run_meta.json").write_text(
            json.dumps(
                {
                    "layer": "d2",
                    "gnina": str(binary),
                    "gnina_version": gnina_version(binary),
                    "no_gpu": bool(args.no_gpu),
                    "note": "native self-dock only; do not merge with D1 RMSD table",
                },
                indent=2,
            )
            + "\n"
        )

    rows = []
    for ligand, seed in jobs:
        meta = preps[ligand]
        sdf = D2_OUT / meta["pdb"].lower() / f"seed{seed}" / f"{ligand}_out.sdf"
        if not args.score_only:
            run_gnina(
                binary,
                PROJECT_ROOT / meta["receptor"],
                LIGAND_PDBQT[ligand],
                meta["center"],
                meta["size"],
                sdf,
                seed,
                args.cpu,
                args.timeout,
                args.no_gpu,
            )
        row = evaluate_selfdock(ligand, meta["pdb"], seed, sdf, PROJECT_ROOT / meta["crystal_ref"])
        rows.append(row)
        print(
            f"EVAL {ligand}@{meta['pdb']} seed{seed}: "
            f"top1={row.get('top1_pose_rmsd_A')} best9={row.get('best_of_9_pose_rmsd_A')} "
            f"class={row.get('capability_class')}",
            flush=True,
        )
    pd.DataFrame(rows).to_csv(D2_OUT / "d2_capability_summary.csv", index=False)
    (D2_OUT / "d2_capability_summary.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
