#!/usr/bin/env python3
"""Compute redock metrics for vina/gnina pose ensembles."""
from __future__ import annotations
import argparse, json, re, csv
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def pose_rmsd(docked: Chem.Mol, ref: Chem.Mol) -> float:
    """Heavy-atom pose RMSD in place (no Kabsch), best symmetry match."""
    d = Chem.RemoveHs(Chem.Mol(docked))
    r = Chem.RemoveHs(Chem.Mol(ref))
    if d.GetNumAtoms() != r.GetNumAtoms():
        raise ValueError(f"atom count mismatch {d.GetNumAtoms()} vs {r.GetNumAtoms()}")
    matches = r.GetSubstructMatches(d, uniquify=False)
    if not matches:
        matches = d.GetSubstructMatches(r, uniquify=False)
        # match[i]=index in d of query atom i in r -> invert
        inv = []
        for m in matches:
            mapping = [None]*len(m)
            for qi, di in enumerate(m):
                mapping[di] = qi
            inv.append(tuple(mapping))
        matches = inv
        # now match[i] should be index in r for atom i in d
        # after invert: mapping[di]=qi means d atom di <-> r atom qi, so match_for_d[di]=qi
    best = float("inf")
    cd, cr = d.GetConformer(), r.GetConformer()
    n = d.GetNumAtoms()
    for match in matches:
        if any(x is None for x in match) or len(match) != n:
            continue
        s = 0.0
        for i, j in enumerate(match):
            a, b = cd.GetAtomPosition(i), cr.GetAtomPosition(j)
            s += (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2
        best = min(best, (s/n)**0.5)
    if not np.isfinite(best):
        raise ValueError("no isomorphism match for RMSD")
    return float(best)

def load_sdf_poses(path: Path):
    suppl = Chem.SDMolSupplier(str(path), removeHs=False)
    mols = []
    for m in suppl:
        if m is not None:
            mols.append(m)
    return mols

def pdbqt_to_mols(path: Path):
    """Split multi-model PDBQT and convert via meeko or RDKit PDB block (lossy). Prefer SDF."""
    text = path.read_text()
    models = re.split(r"MODEL\s+\d+", text)
    # first chunk may be without MODEL header
    chunks = []
    if "MODEL" in text:
        parts = re.split(r"\nMODEL\s+\d+\n", text)
        for p in parts:
            if "ATOM" in p or "HETATM" in p:
                chunks.append(p)
    else:
        chunks = [text]
    # Use openbabel via temporary files is easier - caller should provide SDF
    raise NotImplementedError("provide SDF poses")

def parse_vina_scores_from_pdbqt(path: Path):
    scores = []
    for line in path.read_text().splitlines():
        if "REMARK VINA RESULT:" in line:
            # REMARK VINA RESULT:    -8.123    0.000    0.000
            parts = line.split()
            scores.append(float(parts[3]))
    return scores

def gnina_props(m):
    # RDKit stores SDF props
    def getf(*keys, default=None):
        for k in keys:
            if m.HasProp(k):
                try:
                    return float(m.GetProp(k))
                except Exception:
                    pass
        return default
    return {
        "affinity": getf("minimizedAffinity", "CNNaffinity", "affinity"),  # gnina uses minimizedAffinity
        "CNNscore": getf("CNNscore"),
        "CNNaffinity": getf("CNNaffinity"),
        "vina_affinity": getf("minimizedAffinity"),
    }

def evaluate(poses, scores, higher_better, ref):
    # scores aligned with poses; rank by score
    order = sorted(range(len(poses)), key=lambda i: scores[i], reverse=higher_better)
    rmsds = [pose_rmsd(poses[i], ref) for i in range(len(poses))]
    top1 = order[0]
    return {
        "top1_score": scores[top1],
        "top1_rmsd_A": rmsds[top1],
        "best_ensemble_rmsd_A": min(rmsds),
        "sampling_ok_le2": min(rmsds) <= 2.0,
        "scoring_ok_top1_le2": rmsds[top1] <= 2.0,
        "all_rmsds": rmsds,
        "ranked_indices": order,
        "ranked_rmsds": [rmsds[i] for i in order],
        "ranked_scores": [scores[i] for i in order],
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--poses-sdf", required=True)
    ap.add_argument("--score-field", required=True, help="vina|CNNaffinity|CNNscore|affinity|rtmscore|propname")
    ap.add_argument("--higher-better", action="store_true")
    ap.add_argument("--scores-json", help="optional external scores list JSON")
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()
    ref = Chem.SDMolSupplier(args.ref, removeHs=False)[0]
    poses = load_sdf_poses(Path(args.poses_sdf))
    if args.scores_json:
        scores = json.loads(Path(args.scores_json).read_text())
    else:
        scores = []
        for m in poses:
            if args.score_field == "vina":
                # minimizedAffinity or SD prop from conversion
                for k in ("minimizedAffinity", "affinity", "VINA_SCORE"):
                    if m.HasProp(k):
                        scores.append(float(m.GetProp(k))); break
                else:
                    raise SystemExit("no vina score prop on pose")
            else:
                if not m.HasProp(args.score_field):
                    raise SystemExit(f"missing prop {args.score_field}; have {list(m.GetPropNames())}")
                scores.append(float(m.GetProp(args.score_field)))
    higher = args.higher_better
    if args.score_field in ("CNNaffinity", "CNNscore", "rtmscore"):
        higher = True
    if args.score_field in ("vina", "affinity", "minimizedAffinity"):
        higher = False
    res = evaluate(poses, scores, higher, ref)
    res["n_poses"] = len(poses)
    res["score_field"] = args.score_field
    res["higher_better"] = higher
    Path(args.out_json).write_text(json.dumps(res, indent=2))
    print(json.dumps({k: res[k] for k in ("top1_score","top1_rmsd_A","best_ensemble_rmsd_A","sampling_ok_le2","scoring_ok_top1_le2")}, indent=2))

if __name__ == "__main__":
    main()
