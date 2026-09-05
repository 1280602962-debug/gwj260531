#!/usr/bin/env python3
"""C1 L2b: Arg477 / crystal-pose constrained lesinurad self-dock.

Triggered because free docking (L2) failed CNNscore-selected RMSD<=2 & Arg<=4.
Per LOCAL_C1_CANDIDATE_CAMPAIGN.md: do constrained self-dock before any L3;
if this path is used for candidates, article language = crystal acid-pose matching,
NOT activity retrieval.

Protocol (pre-registered here):
  A) Build deprotonated lesinurad on crystal heavy-atom coords → Meeko PDBQT
  B) gnina --local_only minimize (seed 42/43/44)
  C) Optional tight-box search (8 Å cube on crystal COM) as sensitivity
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from parse_c1_sdf_readouts import evaluate_selfdock  # noqa: E402


def crystal_carboxylate_pdbqt(crystal_sdf: Path, out_pdbqt: Path) -> str:
    ref = Chem.SDMolSupplier(str(crystal_sdf), removeHs=False)[0]
    if ref is None:
        raise SystemExit(f"cannot read {crystal_sdf}")
    # deprotonated template
    smi = "O=C([O-])CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12"
    tmpl = Chem.MolFromSmiles(smi)
    ref_h = Chem.RemoveHs(ref)
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, ref_h)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    pdbqt = PDBQTWriterLegacy.write_string(setups[0])
    if isinstance(pdbqt, tuple):
        pdbqt = pdbqt[0]
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.write_text(pdbqt)
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def run_gnina_local(gnina: Path, receptor: Path, ligand: Path, center, size, out_sdf: Path, seed: int, local_only: bool, cpu: int = 4) -> None:
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(gnina), "-r", str(receptor), "-l", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", "32", "--num_modes", "9", "--cpu", str(cpu),
        "--cnn_scoring", "rescore", "--seed", str(seed),
        "-o", str(out_sdf), "--log", str(log), "--no_gpu",
    ]
    if local_only:
        cmd.append("--local_only")
    print("RUN", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    (out_sdf.parent / (out_sdf.stem + "_stdout.txt")).write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0 and not out_sdf.exists():
        raise RuntimeError(proc.stderr[:500])


def main() -> None:
    refs = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
    out_root = PROJECT_ROOT / "data/campaigns/c1/02_selfdock/constrained_arg477"
    out_root.mkdir(parents=True, exist_ok=True)
    crystal_sdf = refs / "lesinurad_crystal_ref.sdf"
    lig_pdbqt = out_root / "lesinurad_crystal_carboxylate.pdbqt"
    smi = crystal_carboxylate_pdbqt(crystal_sdf, lig_pdbqt)
    (out_root / "ligand_smiles.txt").write_text(smi + "\n")

    # crystal COM for tight box
    mol = Chem.SDMolSupplier(str(crystal_sdf), removeHs=True)[0]
    conf = mol.GetConformer()
    import numpy as np
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    com = xyz.mean(axis=0).tolist()
    (out_root / "crystal_com.json").write_text(json.dumps({"com": com, "tight_box": 8.0}, indent=2))

    gnina = PROJECT_ROOT / "tools" / "gnina"
    receptor = PROJECT_ROOT / "data/structures/prepared/9DKB_receptor.pdbqt"
    # production box still used for local_only (ignored for search extent mostly)
    prod_center = [99.966, 102.967, 105.699]
    prod_size = [22, 22, 22]
    tight_size = [8, 8, 8]

    results = []
    for seed in [42, 43, 44]:
        # A: local_only from crystal carboxylate
        sdf_a = out_root / f"local_only/seed{seed}/lesinurad_out.sdf"
        if not (sdf_a.exists() and sdf_a.stat().st_size > 0):
            run_gnina_local(gnina, receptor, lig_pdbqt, prod_center, prod_size, sdf_a, seed, local_only=True)
        m = evaluate_selfdock(sdf_a, crystal_sdf, refs / "arg477_coords.json", "lesinurad", seed, "urat1_9dkb")
        m.pop("_poses", None)
        m["protocol"] = "local_only_crystal_carboxylate"
        results.append(m)
        print(json.dumps(m, indent=2), flush=True)

        # B: tight-box free search (sensitivity)
        sdf_b = out_root / f"tight_box8/seed{seed}/lesinurad_out.sdf"
        if not (sdf_b.exists() and sdf_b.stat().st_size > 0):
            run_gnina_local(gnina, receptor, lig_pdbqt, com, tight_size, sdf_b, seed, local_only=False)
        m2 = evaluate_selfdock(sdf_b, crystal_sdf, refs / "arg477_coords.json", "lesinurad", seed, "urat1_9dkb")
        m2.pop("_poses", None)
        m2["protocol"] = "tight_box_8A_crystal_com"
        results.append(m2)
        print(json.dumps(m2, indent=2), flush=True)

    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(out_root / "l2b_constrained_metrics.csv", index=False)
    summary = {
        "protocols": ["local_only_crystal_carboxylate", "tight_box_8A_crystal_com"],
        "local_only_pass": {
            r["seed"]: r.get("pass") for r in results if r["protocol"] == "local_only_crystal_carboxylate"
        },
        "tight_box_pass": {
            r["seed"]: r.get("pass") for r in results if r["protocol"] == "tight_box_8A_crystal_com"
        },
        "claim_language": "crystal acid-pose matching / constrained geometry — NOT docking rank activity retrieval",
        "allow_L3_rank_track": False,
    }
    # L3 Rank still requires free-dock CNNscore gate; constrained success only enables Acid-track geometry language
    (out_root / "l2b_gate_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
