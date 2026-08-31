#!/usr/bin/env python3
"""Pose-level diagnostics for PIK3CA/mTOR failure typology v0.

Reads RTM-best export pack + panel/score tables; writes interaction,
chemotype, and asymmetry CSVs under analysis/failure_typology_v0/tables/.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS

# scripts/ → failure_typology_v0/ → analysis/ → pik3ca_mtor_panel48_v0/
ROOT = Path(__file__).resolve().parents[3]
EXPORT = ROOT / "analysis" / "rtm_best_pose_export_v1"
OUT = ROOT / "analysis" / "failure_typology_v0"
TABLES = OUT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

HINGE = {
    "4L23": {"resi": "851", "name": "VAL851"},
    "4JT6": {"resi": "2240", "name": "VAL2240"},
}
RECEPTORS = {
    "4L23": EXPORT / "receptors" / "4L23_PIK3CA_prepared.pdb",
    "4JT6": EXPORT / "receptors" / "4JT6_mTOR_prepared.pdb",
}
CRYSTAL = {
    "4L23": ROOT / "tables" / "4L23_cocrystal_X6K.pdb",
    "4JT6": ROOT / "tables" / "4JT6_cocrystal_X6K.pdb",
}
CASE_IDS = [
    "PM48_01",
    "PM48_02",
    "PM48_10",
    "PM48_20",
    "PM48_21",
    "PM48_26",
    "PM48_34",
]


def parse_pdb(path: Path, atom_only: bool = False, hetatm_only: bool = False):
    rows = []
    for line in path.read_text().splitlines():
        if atom_only and not line.startswith("ATOM"):
            continue
        if hetatm_only and not line.startswith("HETATM"):
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        resi = line[22:26].strip()
        elem = (line[76:78].strip() or name[:1]).upper()
        if elem.startswith("H"):
            continue
        xyz = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
        rows.append({"name": name, "resi": resi, "elem": elem, "xyz": xyz})
    return rows


def lig_from_sdf(path: Path):
    mol = Chem.MolFromMolFile(str(path), removeHs=True)
    conf = mol.GetConformer()
    xyz, elems = [], []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        xyz.append([p.x, p.y, p.z])
        elems.append(atom.GetSymbol())
    return np.asarray(xyz), elems, mol


def hinge_bb(receptor: Path, resi: str):
    return [
        r
        for r in parse_pdb(receptor)
        if r["resi"] == resi and r["name"] in ("N", "O", "C", "CA")
    ]


def hinge_hbond(lig_xyz, elems, hinge_atoms, cutoff=3.5):
    lig_no = np.array([x for x, e in zip(lig_xyz, elems) if e in ("N", "O")])
    bb_no = np.array([a["xyz"] for a in hinge_atoms if a["name"] in ("N", "O")])
    if len(lig_no) == 0 or len(bb_no) == 0:
        return False, float("nan")
    d = np.linalg.norm(lig_no[:, None, :] - bb_no[None, :, :], axis=2)
    mind = float(d.min())
    return mind < cutoff, mind


def clash_count(lig_xyz, prot_xyz, cutoff=2.2):
    n = 0
    for i in range(0, len(lig_xyz), 25):
        L = lig_xyz[i : i + 25]
        d = np.linalg.norm(L[:, None, :] - prot_xyz[None, :, :], axis=2)
        n += int((d < cutoff).sum())
    return n


def occupancy(lig_xyz, crystal_xyz, cutoff=2.0):
    d = np.linalg.norm(lig_xyz[:, None, :] - crystal_xyz[None, :, :], axis=2)
    frac = float((d.min(axis=1) < cutoff).mean())
    cent = float(np.linalg.norm(lig_xyz.mean(0) - crystal_xyz.mean(0)))
    return frac, cent


def mcs_rmsd(ref_mol, pose_mol):
    try:
        return float(AllChem.GetBestRMS(Chem.Mol(ref_mol), Chem.Mol(pose_mol)))
    except Exception:
        pass
    mcs = rdFMCS.FindMCS(
        [ref_mol, pose_mol],
        timeout=30,
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
    )
    if mcs.numAtoms < 5:
        return float("nan")
    patt = Chem.MolFromSmarts(mcs.smartsString)
    rm = ref_mol.GetSubstructMatch(patt)
    pm = pose_mol.GetSubstructMatch(patt)
    if not rm or not pm:
        return float("nan")
    cr, cp = ref_mol.GetConformer(), pose_mol.GetConformer()
    sq = 0.0
    for i, j in zip(rm, pm):
        a, b = cr.GetAtomPosition(i), cp.GetAtomPosition(j)
        sq += (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
    return math.sqrt(sq / len(rm))


def main() -> None:
    inv = list(csv.DictReader(open(EXPORT / "pose_inventory.csv")))
    panel = {
        r["panel_id"]: r
        for r in csv.DictReader(open(ROOT / "tables" / "panel_v0_48.csv"))
    }
    scores = {
        r["ligand"]: r
        for r in csv.DictReader(
            open(ROOT / "tables" / "ablation_ligand_scores.csv")
        )
    }
    ranks = {
        r["ligand"]: r
        for r in csv.DictReader(open(ROOT / "tables" / "ablation_ranks.csv"))
    }

    prot = {
        t: np.array([a["xyz"] for a in parse_pdb(p, atom_only=True)])
        for t, p in RECEPTORS.items()
    }
    crystal_xyz = {
        t: np.array([a["xyz"] for a in parse_pdb(p, hetatm_only=True)])
        for t, p in CRYSTAL.items()
    }
    crystal_mol = {}
    for t, p in CRYSTAL.items():
        m = Chem.MolFromPDBFile(str(p), removeHs=True, sanitize=False)
        if m is not None:
            try:
                Chem.SanitizeMol(m)
            except Exception:
                pass
        crystal_mol[t] = m
    hinge_atoms = {
        t: hinge_bb(RECEPTORS[t], HINGE[t]["resi"]) for t in HINGE
    }

    interaction = []
    for row in inv:
        lig, tgt = row["ligand"], row["target"]
        xyz, elems, mol = lig_from_sdf(EXPORT / row["sdf"])
        hb, hb_d = hinge_hbond(xyz, elems, hinge_atoms[tgt])
        frac, cent = occupancy(xyz, crystal_xyz[tgt])
        rms = (
            mcs_rmsd(crystal_mol[tgt], mol)
            if crystal_mol[tgt] is not None
            else float("nan")
        )
        interaction.append(
            {
                "ligand": lig,
                "role": row["role"],
                "target": tgt,
                "rtm_best_mode": row["rtm_best_mode"],
                "vina_score_mode": row["vina_score"],
                "hinge_residue": HINGE[tgt]["name"],
                "hinge_hbond": "yes" if hb else "no",
                "hinge_min_NO_to_bb_A": round(hb_d, 3) if hb_d == hb_d else "",
                "n_clash_lt_2.2A": clash_count(xyz, prot[tgt]),
                "crystal_occupancy_frac_lt2A": round(frac, 3),
                "centroid_dist_to_X6K_A": round(cent, 3),
                "mcs_rmsd_vs_X6K_A": round(rms, 3) if rms == rms else "",
                "n_heavy": len(xyz),
                "rtmscore": scores[lig][f"rtm_{tgt}"],
                "rtm_z": scores[lig][f"rtm_{tgt}_z"],
            }
        )
    with open(TABLES / "interaction_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(interaction[0].keys()))
        w.writeheader()
        w.writerows(interaction)

    fps, mols = {}, {}
    for lid in CASE_IDS:
        m = Chem.MolFromSmiles(panel[lid]["smiles"])
        mols[lid] = m
        fps[lid] = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)

    sim_rows = []
    for i, a in enumerate(CASE_IDS):
        for b in CASE_IDS[i + 1 :]:
            mcs = rdFMCS.FindMCS(
                [mols[a], mols[b]], timeout=20, completeRingsOnly=True
            )
            sim_rows.append(
                {
                    "ligand_a": a,
                    "name_a": panel[a]["pref_name"] or a,
                    "ligand_b": b,
                    "name_b": panel[b]["pref_name"] or b,
                    "tanimoto_morgan2": round(
                        DataStructs.TanimotoSimilarity(fps[a], fps[b]), 3
                    ),
                    "mcs_atoms": mcs.numAtoms,
                }
            )
    with open(TABLES / "chem_similarity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sim_rows[0].keys()))
        w.writeheader()
        w.writerows(sim_rows)

    pi_rows = []
    for a in CASE_IDS:
        pi_rows.append(
            {
                "ligand": a,
                "pref_name": panel[a]["pref_name"],
                "class": panel[a]["class"],
                "tanimoto_to_PI103": round(
                    DataStructs.TanimotoSimilarity(fps[a], fps["PM48_01"]), 3
                ),
                "scaffold": panel[a]["murcko_scaffold"][:100],
            }
        )
    with open(TABLES / "chem_similarity_to_PI103.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pi_rows[0].keys()))
        w.writeheader()
        w.writerows(pi_rows)

    smarts = {
        "morpholine": "C1COCCN1",
        "amino_triazine_loose": "c1ncnc(N)n1",
    }
    flags = []
    for lid in CASE_IDS:
        m = Chem.MolFromSmiles(panel[lid]["smiles"])
        rec = {
            "ligand": lid,
            "pref_name": panel[lid]["pref_name"],
            "class": panel[lid]["class"],
        }
        for name, sma in smarts.items():
            q = Chem.MolFromSmarts(sma)
            rec[name] = bool(m.GetSubstructMatch(q)) if q else False
        flags.append(rec)
    with open(TABLES / "scaffold_flags.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(flags[0].keys()))
        w.writeheader()
        w.writerows(flags)

    asym = []
    for lid in CASE_IDS:
        s, r, p = scores[lid], ranks[lid], panel[lid]
        rtm_a, rtm_b = float(s["rtm_4L23"]), float(s["rtm_4JT6"])
        asym.append(
            {
                "ligand": lid,
                "pref_name": p["pref_name"],
                "class": p["class"],
                "pchembl_PIK3CA": p["pchembl_PIK3CA"],
                "pchembl_MTOR": p["pchembl_MTOR"],
                "vina_mean_rank": r["vina_mean"],
                "rtm_min_z_rank": r["rtm_min_z"],
                "rtm_4L23": round(rtm_a, 2),
                "rtm_4JT6": round(rtm_b, 2),
                "rtm_shortfall_abs": round(abs(rtm_a - rtm_b), 2),
                "weak_end": "4L23" if rtm_a < rtm_b else "4JT6",
                "rtm_min_z": round(float(s["rtm_min_z"]), 3),
                "clash_fail": s["clash_fail"],
                "role": next(x["role"] for x in inv if x["ligand"] == lid),
            }
        )
    with open(TABLES / "score_asymmetry.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asym[0].keys()))
        w.writeheader()
        w.writerows(asym)

    # shortfall pre-test summary
    rows = []
    for lig, s in scores.items():
        za, zb = float(s["rtm_4L23_z"]), float(s["rtm_4JT6_z"])
        for lam in (0.0, 0.25, 0.5, 1.0):
            rows.append(
                {
                    "ligand": lig,
                    "class": s["class"],
                    "lambda": lam,
                    "score": min(za, zb) - lam * abs(za - zb),
                }
            )
    with open(TABLES / "shortfall_pretest_scores.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = []
    for lam in (0.0, 0.25, 0.5, 1.0):
        sub = [r for r in rows if r["lambda"] == lam]
        sub = sorted(sub, key=lambda r: -r["score"])
        top = sub[:10]
        duals = [r["score"] for r in sub if r["class"] == "dual"]
        rest = [r["score"] for r in sub if r["class"] != "dual"]
        correct = 0.0
        for d in duals:
            correct += sum(1 for x in rest if d > x) + 0.5 * sum(
                1 for x in rest if d == x
            )
        auroc = correct / (len(duals) * len(rest))
        summary.append(
            {
                "lambda": lam,
                "auroc_dual_vs_rest": round(auroc, 4),
                "top10_dual": sum(1 for r in top if r["class"] == "dual"),
                "top10_A_only": sum(1 for r in top if r["class"] == "A_only"),
                "top10_B_only": sum(1 for r in top if r["class"] == "B_only"),
                "top1": top[0]["ligand"],
            }
        )
    with open(TABLES / "shortfall_pretest_metrics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    print("Wrote tables to", TABLES)
    for s in summary:
        print(s)


if __name__ == "__main__":
    main()
