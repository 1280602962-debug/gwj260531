#!/usr/bin/env python3
"""C2 open-source fallback: AutoDock Vina multi-seed redock for 690/2157.

Uses cognate-ligand centroid boxes from RCSB PDBs.
Does NOT replace archived Glide ranks; reports pose consensus only.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
except ImportError as e:
    raise SystemExit(f"meeko required: {e}")

ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "data" / "structures" / "pdb"
OUT = ROOT / "results" / "pose_consensus"
OUT.mkdir(parents=True, exist_ok=True)

VINA = Path("/tmp/vina")
if not VINA.exists():
    VINA = Path("vina")

LIGANDS = {
    "690": "Cc1ccc(NC(=O)Nc2ncc3c(n2)CC(C)(C)CC3=O)cc1",
    "2231": "COc1nc(NCc2ccccc2CN2CCCC2=O)ncc1F",
}

# Primary receptors used in MD shortlist
RECEPTORS = {
    "JNK1": {"pdb": "3ELJ", "chain": "A"},
    "JNK2": {"pdb": "3E7O", "chain": "A"},
    "JNK3": {"pdb": "3TTI", "chain": "A"},
}

SEEDS = [1, 2, 3]
BOX_PAD = 8.0  # Å half-size around cognate centroid -> size 2*pad if small ligand extent


def load_pdb_lines(pdb_id: str):
    path = PDB_DIR / f"{pdb_id}.pdb"
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def cognate_coords(lines, chain: str):
    """Heavy-atom coords of first non-water HETATM ligand on chain."""
    skip = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "DMSO", "CL", "NA", "MG", "ZN", "CA"}
    by_res = {}
    for ln in lines:
        if not ln.startswith("HETATM"):
            continue
        resn = ln[17:20].strip()
        ch = ln[21]
        if ch != chain or resn in skip:
            continue
        # skip covalently modified AA residues sometimes labeled HETATM
        if resn in {"MSE", "CSO", "PTR", "TPO", "SEP"}:
            continue
        resi = ln[22:26].strip()
        key = (resn, resi)
        elem = ln[76:78].strip() or ln[12:16].strip()[0]
        if elem == "H":
            continue
        x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        by_res.setdefault(key, []).append((x, y, z))
    if not by_res:
        raise RuntimeError("No cognate ligand found")
    # pick largest hetero residue
    key = max(by_res, key=lambda k: len(by_res[k]))
    coords = np.array(by_res[key], dtype=float)
    return key, coords


def prepare_receptor_meeko(pdb_id: str, chain: str, out_pdbqt: Path, box_center, box_size):
    """Write chain-only PDB and prepare receptor PDBQT via meeko CLI."""
    lines = load_pdb_lines(pdb_id)
    prot = [ln for ln in lines if ln.startswith("ATOM") and ln[21] == chain]
    prot.append("END")
    tmp_pdb = out_pdbqt.with_suffix(".chain.pdb")
    tmp_pdb.write_text("\n".join(prot) + "\n", encoding="utf-8")

    mk = Path.home() / ".local/bin/mk_prepare_receptor.py"
    if not mk.exists():
        raise RuntimeError("mk_prepare_receptor.py not found; install meeko")
    basename = out_pdbqt.with_suffix("")
    cmd = [
        "python3",
        str(mk),
        "--read_pdb",
        str(tmp_pdb),
        "-o",
        str(basename),
        "-p",
        "--default_altloc",
        "A",
        "-a",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # meeko writes basename.pdbqt
    produced = Path(str(basename) + ".pdbqt")
    if not produced.exists():
        raise RuntimeError(f"receptor prep failed: {proc.stderr or proc.stdout}")
    if produced != out_pdbqt:
        out_pdbqt.write_bytes(produced.read_bytes())


def prepare_ligand_pdbqt(smiles: str, out_pdbqt: Path, seed: int):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()
    if isinstance(pdbqt_string, tuple):
        pdbqt_string, is_ok, error_msg = pdbqt_string
        if not is_ok:
            raise RuntimeError(error_msg or "ligand pdbqt failed")
    if not pdbqt_string:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setups[0])
        if not is_ok:
            raise RuntimeError(error_msg)
    out_pdbqt.write_text(pdbqt_string, encoding="utf-8")
    return mol


def rmsd_ignore_h(mol_ref, mol_prb):
    """Best-effort heavy-atom RMSD after RDKit alignment of docked conformer graphs.
    For Vina poses we rebuild from pdbqt via meeko if possible; else return nan.
    """
    try:
        ref = Chem.RemoveHs(mol_ref)
        prb = Chem.RemoveHs(mol_prb)
        if ref.GetNumAtoms() != prb.GetNumAtoms():
            return float("nan")
        # map by atomic numbers order if isomorphic
        AllChem.AlignMol(prb, ref)
        conf1 = ref.GetConformer()
        conf2 = prb.GetConformer()
        d = []
        for i in range(ref.GetNumAtoms()):
            p1 = conf1.GetAtomPosition(i)
            p2 = conf2.GetAtomPosition(i)
            d.append((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
        return float(np.sqrt(np.mean(d)))
    except Exception:
        return float("nan")


def parse_vina_pdbqt_poses(pdbqt_path: Path):
    """Return list of (score, pdbqt_block)."""
    text = pdbqt_path.read_text(encoding="utf-8", errors="replace")
    blocks = []
    cur = []
    score = None
    for ln in text.splitlines():
        if ln.startswith("MODEL"):
            cur = [ln]
            score = None
        elif ln.startswith("REMARK VINA RESULT:"):
            score = float(ln.split()[3])
            cur.append(ln)
        elif ln.startswith("ENDMDL"):
            cur.append(ln)
            blocks.append((score, "\n".join(cur) + "\n"))
            cur = []
        else:
            if cur:
                cur.append(ln)
    if not blocks and "REMARK VINA RESULT" in text:
        # single pose
        for ln in text.splitlines():
            if ln.startswith("REMARK VINA RESULT:"):
                score = float(ln.split()[3])
        blocks = [(score, text)]
    return blocks


def pdbqt_pose_to_rdkit(block: str):
    try:
        pdbqt_mol = PDBQTMolecule(block, is_dlg=False, skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        return mols[0] if mols else None
    except Exception:
        return None


def run_vina(receptor, ligand, out_pdbqt, center, size, seed, exhaustiveness=16):
    cmd = [
        str(VINA),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--out",
        str(out_pdbqt),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        "5",
        "--seed",
        str(seed),
        "--cpu",
        "2",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout)
    return proc.stdout


def main():
    if not Path(VINA).exists() and subprocess.run(["which", "vina"], capture_output=True).returncode != 0:
        raise SystemExit("AutoDock Vina binary not found")

    work = OUT / "vina_work"
    work.mkdir(exist_ok=True)

    # prepare receptors + boxes
    boxes = {}
    for iso, meta in RECEPTORS.items():
        lines = load_pdb_lines(meta["pdb"])
        key, coords = cognate_coords(lines, meta["chain"])
        center = coords.mean(axis=0)
        extent = coords.max(axis=0) - coords.min(axis=0)
        size = np.maximum(extent + 2 * BOX_PAD, 18.0)
        rec_pdbqt = work / f"{meta['pdb']}_receptor.pdbqt"
        prepare_receptor_meeko(meta["pdb"], meta["chain"], rec_pdbqt, center, size)
        boxes[iso] = {
            "pdb": meta["pdb"],
            "cognate": f"{key[0]}:{key[1]}",
            "center": center.tolist(),
            "size": size.tolist(),
            "receptor_pdbqt": str(rec_pdbqt),
        }

    score_rows = []
    pose_mols = {}  # (cid, iso, seed) -> rdkit mol top1

    for cid, smi in LIGANDS.items():
        for iso, box in boxes.items():
            for seed in SEEDS:
                lig_q = work / f"{cid}_{iso}_seed{seed}_lig.pdbqt"
                out_q = work / f"{cid}_{iso}_seed{seed}_out.pdbqt"
                prepare_ligand_pdbqt(smi, lig_q, seed=seed)
                log = run_vina(
                    box["receptor_pdbqt"],
                    lig_q,
                    out_q,
                    box["center"],
                    box["size"],
                    seed=seed,
                )
                poses = parse_vina_pdbqt_poses(out_q)
                if not poses:
                    score_rows.append(
                        {
                            "compound_id": cid,
                            "isoform": iso,
                            "pdb": box["pdb"],
                            "seed": seed,
                            "vina_score": np.nan,
                            "status": "no_pose",
                        }
                    )
                    continue
                best_score, best_block = poses[0]
                mol = pdbqt_pose_to_rdkit(best_block)
                pose_mols[(cid, iso, seed)] = mol
                score_rows.append(
                    {
                        "compound_id": cid,
                        "isoform": iso,
                        "pdb": box["pdb"],
                        "seed": seed,
                        "vina_score": best_score,
                        "status": "ok",
                        "n_modes": len(poses),
                    }
                )

    scores = pd.DataFrame(score_rows)
    scores.to_csv(OUT / "c2_vina_scores_by_seed.csv", index=False)

    # pairwise RMSD among seeds for same compound/isoform
    rmsd_rows = []
    for cid in LIGANDS:
        for iso in RECEPTORS:
            seeds_present = [s for s in SEEDS if pose_mols.get((cid, iso, s)) is not None]
            for i, s1 in enumerate(seeds_present):
                for s2 in seeds_present[i + 1 :]:
                    r = rmsd_ignore_h(pose_mols[(cid, iso, s1)], pose_mols[(cid, iso, s2)])
                    rmsd_rows.append(
                        {
                            "compound_id": cid,
                            "isoform": iso,
                            "seed_a": s1,
                            "seed_b": s2,
                            "heavy_atom_rmsd_A": r,
                            "consensus_le_2A": (r <= 2.0) if r == r else False,
                        }
                    )
    rmsd_df = pd.DataFrame(rmsd_rows)
    rmsd_df.to_csv(OUT / "c2_pose_rmsd_matrix.csv", index=False)

    # summary
    summary = []
    for cid in LIGANDS:
        for iso in RECEPTORS:
            sub_s = scores[(scores.compound_id == cid) & (scores.isoform == iso)]
            sub_r = rmsd_df[(rmsd_df.compound_id == cid) & (rmsd_df.isoform == iso)]
            frac = float(sub_r["consensus_le_2A"].mean()) if len(sub_r) else float("nan")
            summary.append(
                {
                    "compound_id": cid,
                    "isoform": iso,
                    "mean_vina_score": float(sub_s["vina_score"].mean()),
                    "std_vina_score": float(sub_s["vina_score"].std(ddof=0)) if len(sub_s) else np.nan,
                    "pairwise_rmsd_mean": float(sub_r["heavy_atom_rmsd_A"].mean()) if len(sub_r) else np.nan,
                    "fraction_pairs_rmsd_le_2A": frac,
                    "pose_consensus_pass": bool(frac >= 0.66) if frac == frac else False,
                }
            )
    summ = pd.DataFrame(summary)
    summ.to_csv(OUT / "c2_pose_consensus_summary.csv", index=False)

    meta = {
        "engine": "AutoDock Vina 1.2.5",
        "role": "open-source C2 fallback; does not replace Glide shortlist ranks",
        "seeds": SEEDS,
        "exhaustiveness": 16,
        "boxes": boxes,
        "note": "Receptor prep is chain-ATOM pdbqt without full Protonation Wizard; interpret qualitatively.",
    }
    (OUT / "c2_vina_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    md = [
        "# C2 Pose Consensus — AutoDock Vina Multi-Seed (Open Fallback)",
        "",
        "**Scope:** geometry consensus across seeds, not Glide score validation.",
        "",
        "## Box definitions (cognate ligand)",
        "",
        "```json",
        json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "receptor_pdbqt"} for k, v in boxes.items()}, indent=2),
        "```",
        "",
        "## Consensus summary",
        "",
        summ.to_markdown(index=False),
        "",
        "## Pairwise seed RMSD",
        "",
        rmsd_df.to_markdown(index=False) if len(rmsd_df) else "_no RMSD (pose parse failed)_",
        "",
        "## Scores",
        "",
        scores.to_markdown(index=False),
        "",
    ]
    (OUT / "C2_POSE_CONSENSUS_REPORT.md").write_text("\n".join(md), encoding="utf-8")
    print(summ.to_string(index=False))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
