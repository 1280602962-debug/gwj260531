#!/usr/bin/env python3
"""
Split LigPrep multi-mol SDF into per-molecule SDF + PDBQT (Meeko).
Keeps LigPrep 3D coordinates (no re-embedding).
"""
from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def _safe_props(mol: Chem.Mol) -> dict:
    out = {}
    for k in mol.GetPropNames():
        try:
            out[k] = mol.GetProp(k)
        except Exception:
            pass
    if mol.HasProp("_Name"):
        out["_Name"] = mol.GetProp("_Name")
    return out


def process_one(args):
    idx, sdf_block, out_sdf, out_pdbqt, write_pdbqt = args
    mol_id = f"mol_{idx:05d}"
    try:
        mol = Chem.MolFromMolBlock(sdf_block, removeHs=False, sanitize=False)
        if mol is None:
            return idx, mol_id, "", "parse_fail", ""
        # light sanitize for meeko; keep coords
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            # still try write SDF unsanitized
            pass

        props = _safe_props(mol)
        src_idx = props.get("i_m_Source_File_Index") or props.get("i_m_source_file_index") or ""
        name = props.get("s_m_entry_name") or props.get("_Name") or mol_id
        smiles = ""
        try:
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        except Exception:
            try:
                smiles = Chem.MolToSmiles(mol, canonical=False)
            except Exception:
                smiles = ""

        Path(out_sdf).parent.mkdir(parents=True, exist_ok=True)
        w = Chem.SDWriter(out_sdf)
        mol.SetProp("_Name", mol_id)
        for k, v in props.items():
            if k != "_Name":
                try:
                    mol.SetProp(k, str(v))
                except Exception:
                    pass
        w.write(mol)
        w.close()

        pdbqt_status = "skip"
        if write_pdbqt:
            try:
                from meeko import MoleculePreparation, PDBQTWriterLegacy

                preparator = MoleculePreparation()
                setups = preparator.prepare(mol)
                setup = setups[0] if isinstance(setups, (list, tuple)) else setups
                pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(setup)
                if not is_ok:
                    return idx, mol_id, smiles, f"pdbqt_fail:{err}", src_idx
                Path(out_pdbqt).write_text(pdbqt_string)
                pdbqt_status = "ok"
            except Exception as e:
                return idx, mol_id, smiles, f"pdbqt_exc:{e}", src_idx

        return idx, mol_id, smiles, f"ok:{pdbqt_status}", src_idx
    except Exception as e:
        return idx, mol_id, "", f"exc:{e}", ""


def iter_sdf_blocks(path: Path):
    """Yield (0-based index, molblock string) without full RDKit sanitize."""
    buf = []
    idx = 0
    with open(path, "r", errors="ignore") as f:
        for line in f:
            buf.append(line)
            if line.strip() == "$$$$":
                yield idx, "".join(buf)
                idx += 1
                buf = []
    if buf and any(x.strip() for x in buf):
        yield idx, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--nproc", type=int, default=max(1, mp.cpu_count() - 1))
    ap.add_argument("--write-pdbqt", action="store_true")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=-1, help="exclusive; -1=all")
    args = ap.parse_args()

    out = Path(args.outdir)
    sdf_dir = out / "ligands_sdf"
    pdbqt_dir = out / "ligands_pdbqt"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    if args.write_pdbqt:
        pdbqt_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for i, block in iter_sdf_blocks(Path(args.sdf)):
        if i < args.start:
            continue
        if args.end >= 0 and i >= args.end:
            break
        out_sdf = str(sdf_dir / f"mol_{i:05d}.sdf")
        out_pdbqt = str(pdbqt_dir / f"mol_{i:05d}.pdbqt")
        if Path(out_sdf).exists() and (not args.write_pdbqt or Path(out_pdbqt).exists()):
            continue
        jobs.append((i, block, out_sdf, out_pdbqt, args.write_pdbqt))

    print(f"jobs_to_run={len(jobs)} nproc={args.nproc}", flush=True)
    manifest = out / "ligand_manifest.csv"
    write_header = not manifest.exists()
    with open(manifest, "a", newline="") as mf:
        w = csv.DictWriter(
            mf,
            fieldnames=["idx", "mol_id", "smiles", "status", "source_file_index"],
        )
        if write_header:
            w.writeheader()
        if not jobs:
            print("nothing to do")
            return
        with mp.Pool(args.nproc) as pool:
            for idx, mol_id, smiles, status, src_idx in pool.imap_unordered(
                process_one, jobs, chunksize=8
            ):
                w.writerow(
                    {
                        "idx": idx,
                        "mol_id": mol_id,
                        "smiles": smiles,
                        "status": status,
                        "source_file_index": src_idx,
                    }
                )
                if idx % 200 == 0:
                    mf.flush()
                    print(f"done idx={idx} status={status}", flush=True)
    print("finished", flush=True)


if __name__ == "__main__":
    main()
