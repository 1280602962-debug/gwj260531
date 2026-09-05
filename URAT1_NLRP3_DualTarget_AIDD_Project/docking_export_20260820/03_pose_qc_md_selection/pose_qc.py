#!/usr/bin/env python3
"""Pose reasonableness QC for 7 dual-target nominees vs crystal references."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path("/home/hww/gwj/NLRP3_URAT1")
P2 = ROOT / "server_p2_gpu_upload"
OUT = ROOT / "md_dual_target" / "pose_qc"
OUT.mkdir(parents=True, exist_ok=True)

MOLS = [
    ("VECABRUTINIB", "REP_06421"),
    ("ZELENIRSTAT", "REP_05621"),
    ("DEUCRICTIBANT", "REP_06841"),
    ("PRALICIGUAT", "REP_06295"),
    ("GSK-3008348", "REP_07907"),
    ("MLN-0415", "REP_08167"),
    ("BI 653048", "REP_07792"),
]
CONTROLS = [
    ("LESINURAD", "REP_00207"),
    ("VERINURAD", "REP_05846"),
]

# Literature / UniProt numbers
URAT1_KEY = {
    "SER35": "H-bond (lesinurad triazole / verinurad pyridine)",
    "MET214": "hydrophobic (lesinurad)",
    "PHE241": "pi-stack naphthalene cage",
    "PHE360": "T-stack cage",
    "PHE364": "hydrophobic cage",
    "PHE365": "critical cage (F365Y kills potency)",
    "LYS393": "basic / verinurad H-bond",
    "PHE449": "hydrophobic cage",
    "ARG477": "gating / carboxylate electrostatics",
}
NLRP3_KEY = {
    "ALA227": "Walker A backbone (sulfonylurea)",
    "ALA228": "Walker A polar contact (urea)",
    "ARG351": "NBD clamp on sulfonyl",
    "MET408": "HD1 hydrophobic wing",
    "TYR443": "WHD hydrophobic wing",
    "PHE575": "HD2 hydrophobic",
    "ARG578": "HD2 clamp (critical H-bond)",
}


def load_pdb_atoms(path: Path, hetatm=True):
    atoms = []
    with path.open() as f:
        for line in f:
            if line.startswith("ATOM") or (hetatm and line.startswith("HETATM")):
                name = line[12:16].strip()
                resn = line[17:21].strip()
                chain = line[21]
                resi = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                elem = (line[76:78].strip() if len(line) >= 78 else name[0]).upper()
                if not elem:
                    elem = name[0]
                atoms.append(
                    {
                        "name": name,
                        "resn": resn,
                        "chain": chain,
                        "resi": resi,
                        "xyz": xyz,
                        "elem": elem,
                        "is_h": elem == "H" or name.startswith("H"),
                    }
                )
    return atoms


def sdf_heavy_xyz(path: Path):
    suppl = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise RuntimeError(f"unreadable sdf: {path}")
    conf = mol.GetConformer()
    xyz, elems, is_h = [], [], []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        e = atom.GetSymbol()
        xyz.append([p.x, p.y, p.z])
        elems.append(e)
        is_h.append(e == "H")
    return mol, np.array(xyz), elems, np.array(is_h)


def pdb_heavy_xyz(path: Path):
    atoms = [a for a in load_pdb_atoms(path, hetatm=True) if not a["is_h"]]
    return np.array([a["xyz"] for a in atoms]), [a["elem"] for a in atoms]


def com(xyz, mask=None):
    pts = xyz if mask is None else xyz[mask]
    return pts.mean(axis=0)


def pairwise_min(a, b):
    # a: (n,3) b: (m,3) -> min dist per a
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return d.min(axis=1), d.min(axis=0), d


def hbond_count(lig_xyz, lig_elem, rec_atoms, cutoff=3.5):
    donors_acc = {"N", "O"}
    lig_idx = [i for i, e in enumerate(lig_elem) if e in donors_acc]
    rec = [a for a in rec_atoms if a["elem"] in donors_acc and not a["is_h"]]
    n = 0
    pairs = []
    for i in lig_idx:
        for a in rec:
            d = float(np.linalg.norm(lig_xyz[i] - a["xyz"]))
            if d <= cutoff:
                n += 1
                pairs.append((lig_elem[i], f"{a['resn']}{a['resi']}", a["name"], d))
    return n, pairs


def clash_count(lig_xyz, lig_h, rec_atoms, cutoff=2.2):
    rec = np.array([a["xyz"] for a in rec_atoms if not a["is_h"]])
    lig = lig_xyz[~lig_h]
    if len(lig) == 0 or len(rec) == 0:
        return 0
    dmin, _, _ = pairwise_min(lig, rec)
    return int((dmin < cutoff).sum())


def contacts_by_res(lig_xyz, lig_h, rec_atoms, cutoff=4.5):
    lig = lig_xyz[~lig_h]
    out = {}
    for a in rec_atoms:
        if a["is_h"]:
            continue
        key = (a["chain"], a["resi"], a["resn"])
        d = float(np.min(np.linalg.norm(lig - a["xyz"], axis=1)))
        if key not in out or d < out[key]:
            out[key] = d
    return {k: v for k, v in out.items() if v <= cutoff}


def map_urat1_uniprot(prep_atoms):
    """prep numbering is UniProt-1 (ARG476 == ARG477)."""
    want = {35, 214, 241, 360, 364, 365, 393, 449, 477}
    mapped = {}
    for a in prep_atoms:
        if a["name"] != "CA":
            continue
        uni = a["resi"] + 1
        label = f"{a['resn']}{uni}"
        # expected names
        expected = {
            35: "SER",
            214: "MET",
            241: "PHE",
            360: "PHE",
            364: "PHE",
            365: "PHE",
            393: "LYS",
            449: "PHE",
            477: "ARG",
        }
        if uni in want and a["resn"] == expected[uni]:
            mapped[label] = a
    return mapped


def map_nlrp3_uniprot(prep_atoms, src_atoms):
    src_ca = {
        a["resi"]: a
        for a in src_atoms
        if a["name"] == "CA" and a["chain"] == "A"
    }
    prep_ca = [a for a in prep_atoms if a["name"] == "CA"]
    mapped = {}
    for uni, src in src_ca.items():
        best, bd = None, 1e9
        for p in prep_ca:
            d = float(np.linalg.norm(p["xyz"] - src["xyz"]))
            if d < bd:
                bd, best = d, p
        if best is not None and bd < 1.0:
            mapped[f"{src['resn']}{uni}"] = {**best, "uniprot": uni, "match_d": bd}
    return mapped


def acid_oxygens(mol, xyz):
    """Carboxylate / acid oxygens if present."""
    outs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        # C(=O)O pattern: oxygen attached to carbon that has another O
        neigh = [n for n in atom.GetNeighbors()]
        if any(n.GetSymbol() == "C" for n in neigh):
            c = next(n for n in neigh if n.GetSymbol() == "C")
            o_on_c = [n for n in c.GetNeighbors() if n.GetSymbol() == "O"]
            if len(o_on_c) >= 2:
                outs.append(xyz[atom.GetIdx()])
    return np.array(outs) if outs else np.zeros((0, 3))


def sulfo_atoms(mol, xyz):
    outs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "S":
            o_n = [n for n in atom.GetNeighbors() if n.GetSymbol() in ("O", "N")]
            if len(o_n) >= 2:
                outs.append(xyz[atom.GetIdx()])
                for n in o_n:
                    outs.append(xyz[n.GetIdx()])
    return np.array(outs) if outs else np.zeros((0, 3))


def score_pose(row):
    """Higher is more MD-ready."""
    s = 0.0
    s += max(0.0, 8.0 - row["com_to_ref"]) * 1.2
    s += min(row["n_key_contacts"], 7) * 1.5
    s += min(row["n_hbonds"], 4) * 0.8
    s -= row["n_clashes"] * 2.5
    s += 2.0 if row["in_pocket"] else -4.0
    if row["target"] == "URAT1":
        if np.isfinite(row.get("acid_arg477", np.nan)):
            s += max(0.0, 8.0 - row["acid_arg477"])
        s += min(row.get("n_phe_cage", 0), 5) * 0.8
    else:
        if np.isfinite(row.get("d_arg578", np.nan)):
            s += max(0.0, 6.0 - row["d_arg578"]) * 1.2
        if np.isfinite(row.get("d_ala228", np.nan)):
            s += max(0.0, 6.0 - row["d_ala228"])
        if np.isfinite(row.get("d_arg351", np.nan)):
            s += max(0.0, 6.0 - row["d_arg351"])
    return s


def analyze_target(target, rec_path, ref_xyz, key_map, key_labels, pose_dir, extra_fn=None):
    rec = load_pdb_atoms(rec_path)
    rec_heavy = np.array([a["xyz"] for a in rec if not a["is_h"]])
    rows = []
    items = MOLS + (CONTROLS if target == "URAT1" else [])
    for name, rid in items:
        sdf = pose_dir / f"{rid}_out.sdf"
        mol, xyz, elems, is_h = sdf_heavy_xyz(sdf)
        lig = xyz[~is_h]
        lig_elem = [e for e, h in zip(elems, is_h) if not h]
        c = com(lig)
        dcom = float(np.linalg.norm(c - com(ref_xyz)))
        # overlap with reference: fraction of lig atoms within 2.5 A of any ref atom
        dmin_ref, _, _ = pairwise_min(lig, ref_xyz)
        frac_overlap = float((dmin_ref < 2.5).mean())
        nclash = clash_count(xyz, is_h, rec)
        nhb, hb_pairs = hbond_count(xyz[~is_h], lig_elem, rec)
        contacts = contacts_by_res(xyz, is_h, rec, 4.5)
        # key contacts using mapped CA + any sidechain atom of that residue
        key_hits = {}
        for label, rec_ca in key_map.items():
            if label not in key_labels:
                continue
            # all atoms of that mapped residue
            atoms = [
                a
                for a in rec
                if a["chain"] == rec_ca["chain"]
                and a["resi"] == rec_ca["resi"]
                and not a["is_h"]
            ]
            if not atoms:
                continue
            rax = np.array([a["xyz"] for a in atoms])
            d = float(np.min(np.linalg.norm(lig[:, None, :] - rax[None, :, :], axis=2)))
            key_hits[label] = d
        n_key = sum(1 for d in key_hits.values() if d <= 4.5)
        in_pocket = dcom <= 6.0 or n_key >= 3
        row = {
            "name": name,
            "rep_id": rid,
            "target": target,
            "n_heavy": int(len(lig)),
            "com_to_ref": dcom,
            "frac_ref_overlap_2p5": frac_overlap,
            "n_clashes": nclash,
            "n_hbonds": nhb,
            "n_key_contacts": n_key,
            "in_pocket": bool(in_pocket),
            "hbonds": "; ".join(f"{a}-{b}.{c}:{d:.2f}" for a, b, c, d in hb_pairs[:8]),
        }
        for lab in key_labels:
            row[f"d_{lab.lower()}"] = key_hits.get(lab, np.nan)
        if extra_fn:
            row.update(extra_fn(mol, xyz, elems, is_h, key_map, rec))
        row["pose_score"] = 0.0
        rows.append(row)
        # save stripped pose pdb for MD
        wdir = OUT / "poses" / target / rid
        wdir.mkdir(parents=True, exist_ok=True)
        Chem.MolToPDBFile(mol, str(wdir / "ligand.pdb"))
    df = pd.DataFrame(rows)
    df["pose_score"] = df.apply(score_pose, axis=1)
    return df


def urat1_extra(mol, xyz, elems, is_h, key_map, rec):
    out = {"acid_arg477": np.nan, "n_phe_cage": 0}
    acids = acid_oxygens(mol, xyz)
    if "ARG477" in key_map:
        arg_atoms = [
            a
            for a in rec
            if a["chain"] == key_map["ARG477"]["chain"]
            and a["resi"] == key_map["ARG477"]["resi"]
            and a["name"] in ("CZ", "NH1", "NH2", "NE")
        ]
        if len(acids) and arg_atoms:
            rax = np.array([a["xyz"] for a in arg_atoms])
            out["acid_arg477"] = float(np.min(np.linalg.norm(acids[:, None, :] - rax[None, :, :], axis=2)))
    cage = []
    lig = xyz[~is_h]
    for lab in ("PHE241", "PHE360", "PHE364", "PHE365", "PHE449"):
        if lab not in key_map:
            continue
        atoms = [
            a
            for a in rec
            if a["chain"] == key_map[lab]["chain"]
            and a["resi"] == key_map[lab]["resi"]
            and not a["is_h"]
        ]
        if not atoms:
            continue
        rax = np.array([a["xyz"] for a in atoms])
        d = float(np.min(np.linalg.norm(lig[:, None, :] - rax[None, :, :], axis=2)))
        if d <= 4.5:
            cage.append(lab)
    out["n_phe_cage"] = len(cage)
    out["phe_cage"] = ",".join(cage)
    return out


def nlrp3_extra(mol, xyz, elems, is_h, key_map, rec):
    out = {}
    sulfo = sulfo_atoms(mol, xyz)
    lig = xyz[~is_h]
    for lab, key in (("ARG578", "d_polar_arg578"), ("ALA228", "d_polar_ala228"), ("ARG351", "d_polar_arg351")):
        if lab not in key_map:
            out[key] = np.nan
            continue
        atoms = [
            a
            for a in rec
            if a["chain"] == key_map[lab]["chain"]
            and a["resi"] == key_map[lab]["resi"]
            and a["elem"] in ("N", "O")
            and not a["is_h"]
        ]
        if not atoms:
            out[key] = np.nan
            continue
        rax = np.array([a["xyz"] for a in atoms])
        src = sulfo if len(sulfo) else lig
        out[key] = float(np.min(np.linalg.norm(src[:, None, :] - rax[None, :, :], axis=2)))
    return out


def main():
    rec_u = P2 / "data/structures/prepared/9DKB_receptor.pdbqt"
    rec_n = P2 / "data/structures/prepared/7ALV_receptor.pdbqt"
    src_n = Path(
        "/home/hww/gwj/NLRP3_PLK1/PLK1_NLRP3_opensource_docking_server_pack/inputs/7ALV_receptor_src.pdb"
    )
    les_sdf = Path("/home/hww/gwj/NLRP3_URAT1/server_dock_maestro_prep/inputs/lesinurad_crystal_ref.sdf")
    rm5_pdb = Path(
        "/home/hww/gwj/NLRP3_PLK1/PLK1_NLRP3_opensource_docking_server_pack/inputs/7ALV_crystal_RM5.pdb"
    )

    u_atoms = load_pdb_atoms(rec_u)
    n_atoms = load_pdb_atoms(rec_n)
    src_atoms = load_pdb_atoms(src_n)
    u_map = map_urat1_uniprot(u_atoms)
    n_map = map_nlrp3_uniprot(n_atoms, src_atoms)
    n_map = {k: v for k, v in n_map.items() if k in NLRP3_KEY}

    print("URAT1 mapped", sorted(u_map))
    print("NLRP3 mapped", {k: (v["chain"], v["resi"], v["resn"], round(v["match_d"], 3)) for k, v in n_map.items()})

    _, les_xyz, _, les_h = sdf_heavy_xyz(les_sdf)
    les_xyz = les_xyz[~les_h]
    rm5_xyz, _ = pdb_heavy_xyz(rm5_pdb)

    df_u = analyze_target(
        "URAT1",
        rec_u,
        les_xyz,
        u_map,
        list(URAT1_KEY),
        P2 / "results/repurposing/docking_p2/9dkb/poses",
        urat1_extra,
    )
    df_n = analyze_target(
        "NLRP3",
        rec_n,
        rm5_xyz,
        n_map,
        list(NLRP3_KEY),
        P2 / "results/repurposing/docking_p2/7alv/poses",
        nlrp3_extra,
    )
    df = pd.concat([df_u, df_n], ignore_index=True)
    df.to_csv(OUT / "pose_qc_table.csv", index=False)

    # dual summary
    dual = []
    for name, rid in MOLS:
        u = df_u[df_u.rep_id == rid].iloc[0]
        n = df_n[df_n.rep_id == rid].iloc[0]
        dual.append(
            {
                "name": name,
                "rep_id": rid,
                "urat1_com_to_lesinurad": u.com_to_ref,
                "urat1_key_contacts": u.n_key_contacts,
                "urat1_clashes": u.n_clashes,
                "urat1_hbonds": u.n_hbonds,
                "urat1_acid_arg477": u.acid_arg477,
                "urat1_phe_cage": u.n_phe_cage,
                "urat1_in_pocket": u.in_pocket,
                "urat1_pose_score": u.pose_score,
                "nlrp3_com_to_rm5": n.com_to_ref,
                "nlrp3_key_contacts": n.n_key_contacts,
                "nlrp3_clashes": n.n_clashes,
                "nlrp3_hbonds": n.n_hbonds,
                "nlrp3_d_arg578": n.d_arg578,
                "nlrp3_d_ala228": n.d_ala228,
                "nlrp3_d_arg351": n.d_arg351,
                "nlrp3_in_pocket": n.in_pocket,
                "nlrp3_pose_score": n.pose_score,
                "dual_pose_score": 0.5 * (u.pose_score + n.pose_score),
                "both_in_pocket": bool(u.in_pocket and n.in_pocket),
            }
        )
    ddf = pd.DataFrame(dual).sort_values("dual_pose_score", ascending=False)
    ddf.to_csv(OUT / "pose_qc_dual.csv", index=False)
    print("\n==== URAT1 ====")
    cols = [
        "name",
        "com_to_ref",
        "n_key_contacts",
        "n_clashes",
        "n_hbonds",
        "acid_arg477",
        "n_phe_cage",
        "in_pocket",
        "pose_score",
    ]
    print(df_u[cols].round(2).to_string(index=False))
    print("\n==== NLRP3 ====")
    cols = [
        "name",
        "com_to_ref",
        "n_key_contacts",
        "n_clashes",
        "n_hbonds",
        "d_arg578",
        "d_ala228",
        "d_arg351",
        "in_pocket",
        "pose_score",
    ]
    print(df_n[cols].round(2).to_string(index=False))
    print("\n==== DUAL ====")
    print(ddf.round(2).to_string(index=False))

    summary = {
        "urat1_key_map": {k: f"{v['resn']}{v['resi']}" for k, v in u_map.items()},
        "nlrp3_key_map": {
            k: f"{v['chain']}:{v['resn']}{v['resi']}" for k, v in n_map.items()
        },
        "note": "9DKB prep numbering = UniProt-1; NLRP3 mapped via CA to 7ALV src.",
    }
    (OUT / "residue_map.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
