#!/usr/bin/env python3
"""T1.3 B2: Single-target enrichment sanity for 4L23 (PIK3CA) and 4JT6 (mTOR).

Actives pChEMBL>=6.5; property-matched decoys pChEMBL<=5.5.
Dock E=16, report AUROC / EF1% / EF5%.
"""
from __future__ import annotations

import csv
import json
import random
import re
import subprocess
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "data/jcim_strengthen_t0t1_v0"
WORK = Path("/mnt/d/CADD paper exercise/dual target docking/results/jcim_strengthen_enrichment_v0")
PM48_ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0")

TARGETS = {
    "4L23": {"mols": "mols_PIK3CA.json", "chembl_target": "PIK3CA"},
    "4JT6": {"mols": "mols_MTOR.json", "chembl_target": "MTOR"},
}
VINA = "/home/gwj/miniconda3/bin/vina"
MEKO_PY = "/home/gwj/miniconda3/bin/python"
SEED = 20260727
EXHAUST = 16
N_MODES = 9
N_ACTIVE = 80
N_DECOY = 300
MAX_WORKERS = 6
CACHE = OUT / "tables" / "enrichment_smiles_cache.json"


def fetch_smiles(cid: str, cache: dict) -> str | None:
    if cid in cache:
        return cache[cid]
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
    req = urllib.request.Request(url, headers={"User-Agent": "dualfourclass-jcim/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.load(resp)
        smi = (data.get("molecule_structures") or {}).get("canonical_smiles")
        cache[cid] = smi
        time.sleep(0.03)
        return smi
    except Exception:
        cache[cid] = None
        return None


def phys(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "heavy": mol.GetNumHeavyAtoms(),
    }


def match_decoys(actives, pool, n_decoy, rng):
    """Greedy property match each active to a decoy (MW±50, logP±1, TPSA±20)."""
    chosen = []
    used = set()
    for a in actives:
        cands = []
        for d in pool:
            if d["chembl_id"] in used:
                continue
            if abs(d["mw"] - a["mw"]) > 50 or abs(d["logp"] - a["logp"]) > 1.5 or abs(d["tpsa"] - a["tpsa"]) > 25:
                continue
            dist = abs(d["mw"] - a["mw"]) + 2 * abs(d["logp"] - a["logp"]) + 0.1 * abs(d["tpsa"] - a["tpsa"])
            cands.append((dist, d))
        if cands:
            cands.sort(key=lambda x: x[0])
            d = cands[0][1]
            used.add(d["chembl_id"])
            chosen.append(d)
        if len(chosen) >= n_decoy:
            break
    if len(chosen) < n_decoy:
        rest = [d for d in pool if d["chembl_id"] not in used]
        rng.shuffle(rest)
        chosen.extend(rest[: n_decoy - len(chosen)])
    return chosen[:n_decoy]


def prep_pdbqt(mol_id, smi, root: Path) -> Path | None:
    from rdkit.Chem import AllChem

    pdbqt = root / "ligands_pdbqt" / f"{mol_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=SEED)
    sdf = root / "ligands_sdf" / f"{mol_id}.sdf"
    sdf.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(sdf))
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [MEKO_PY, "/home/gwj/miniconda3/bin/mk_prepare_ligand.py", "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True, text=True,
    )
    return pdbqt if proc.returncode == 0 and pdbqt.exists() else None


def parse_best(log_text):
    affs = []
    for line in log_text.splitlines():
        m = re.search(r"^\s*1\s+(-?\d+\.\d+)", line)
        if m:
            affs.append(float(m.group(1)))
    return min(affs) if affs else None


def dock_one(target, mol_id, pdbqt, root, box, rec):
    log = root / "logs" / target / f"{mol_id}.log"
    out = root / "logs" / target / f"{mol_id}_out.pdbqt"
    log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        VINA, "--receptor", str(rec), "--ligand", str(pdbqt),
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]), "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(EXHAUST), "--num_modes", str(N_MODES), "--seed", str(SEED), "--cpu", "1", "--out", str(out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.write_text(proc.stdout + proc.stderr)
    best = parse_best(proc.stdout + proc.stderr)
    return mol_id, best


def auroc(scores, labels):
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return float("nan")
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def enrichment_factor(scores, labels, pct):
    order = np.argsort(-np.array(scores, float))
    labs = np.array(labels)[order]
    top = max(1, int(len(labs) * pct / 100))
    return float(labs[:top].sum() / (sum(labels) * pct / 100)) if sum(labels) else float("nan")


def main():
    rng = random.Random(SEED)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    WORK.mkdir(parents=True, exist_ok=True)
    # symlink receptors/boxes from PM48
    for sub in ("receptors", "boxes"):
        dst = WORK / sub
        if not dst.exists():
            dst.symlink_to(PM48_ROOT / sub)

    results = []
    for target, cfg in TARGETS.items():
        mols_path = REPO / f"data/public_pair_selection/{cfg['mols']}"
        mols = json.loads(mols_path.read_text())
        act_ids = [k for k, v in mols.items() if v >= 6.5]
        dec_ids = [k for k, v in mols.items() if v <= 5.5]
        rng.shuffle(act_ids)
        rng.shuffle(dec_ids)
        act_ids = act_ids[:N_ACTIVE]
        print(f"{target}: fetching smiles for {len(act_ids)} actives...", flush=True)

        actives = []
        for cid in act_ids:
            smi = fetch_smiles(cid, cache)
            if not smi:
                continue
            p = phys(smi)
            if p:
                actives.append({"chembl_id": cid, "smiles": smi, "label": 1, **p})
        print(f"{target}: {len(actives)} actives with SMILES", flush=True)

        decoy_pool = []
        for cid in dec_ids[:2000]:
            smi = fetch_smiles(cid, cache)
            if not smi:
                continue
            p = phys(smi)
            if p and 150 < p["mw"] < 800:
                decoy_pool.append({"chembl_id": cid, "smiles": smi, "label": 0, **p})
        decoys = match_decoys(actives, decoy_pool, N_DECOY, rng)
        panel = actives + decoys
        print(f"{target}: panel n={len(panel)} act={len(actives)} dec={len(decoys)}", flush=True)

        CACHE.parent.mkdir(parents=True, exist_ok=True)
        CACHE.write_text(json.dumps(cache, indent=2))

        box = json.loads((WORK / "boxes" / f"{target}_box.json").read_text())
        rec = WORK / "receptors" / f"{target}_receptor.pdbqt"

        jobs = []
        for m in panel:
            mid = f"{m['chembl_id']}"
            pdbqt = prep_pdbqt(mid, m["smiles"], WORK)
            if pdbqt:
                jobs.append((mid, pdbqt, m["label"]))

        scores_map = {}
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(dock_one, target, mid, pq, WORK, box, rec): mid for mid, pq, _ in jobs}
            for i, fut in enumerate(as_completed(futs), 1):
                mid, best = fut.result()
                if best is not None:
                    scores_map[mid] = -best  # higher=better
                if i % 20 == 0:
                    print(f"  [{target}] docked {i}/{len(futs)}", flush=True)

        scored = [(scores_map[j[0]], j[2]) for j in jobs if j[0] in scores_map]
        if len(scored) < 20:
            continue
        sc, lab = zip(*scored)
        auc = auroc(sc, lab)
        ef1 = enrichment_factor(sc, lab, 1)
        ef5 = enrichment_factor(sc, lab, 5)
        results.append({
            "receptor": target,
            "target_name": cfg["chembl_target"],
            "n_active": sum(lab),
            "n_decoy": len(lab) - sum(lab),
            "n_docked": len(scored),
            "auroc": round(auc, 4),
            "EF_1pct": round(ef1, 4),
            "EF_5pct": round(ef5, 4),
            "exhaustiveness": EXHAUST,
            "seed": SEED,
        })
        print(f"{target} AUROC={auc:.3f} EF1%={ef1:.2f} EF5%={ef5:.2f}", flush=True)

    write_path = OUT / "tables" / "single_target_enrichment_v1.csv"
    with write_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print("Wrote", write_path)


if __name__ == "__main__":
    main()
