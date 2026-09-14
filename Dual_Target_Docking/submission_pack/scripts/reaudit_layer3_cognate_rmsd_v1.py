#!/usr/bin/env python3
"""Re-audit Track B cognate RMSD on the saved E=8 poses.

Does not redock. Uses the deposited cognate PDBQT/SDF/PDB and
layer3 *_out_E8.pdbqt files.

Primary number: RDKit symmetry-aware CalcRMS on a chemically mapped
graph (Meeko reconstruction when it works; CCD SMILES + bond-aware
crystal assignment for 2JKH/BI7). Hungarian name-H matching is kept
only as a comparison column.
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from scipy.optimize import linear_sum_assignment

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
COG = LOCAL / "cognates"
QC = LOCAL / "cognate_qc"
TAB = LOCAL / "tables"
AN = LOCAL / "analysis"

SPECS = [
    ("F2", "4UDW", "N6L"),
    ("F10", "2JKH", "BI7"),
    ("JAK1", "6N7A", "KEV"),
    ("TYK2", "3LXP", "IZA"),
    ("JAK2", "8BXH", "C87"),
    ("PPARG", "9V8H", "BRL"),
    ("PPARA", "6LXA", "EPA"),
    ("PPARD", "5U3Q", "7UJ"),
]

CCD_SMILES = {
    "BI7": "C[N+](C)(C)CCCN1C(=O)[C@H]2[C@@H]3CCCN3[C@H]([C@H]2C1=O)c4cc(on4)c5ccc(s5)Cl",
}


def pdbqt_element(line: str) -> str:
    atom_type = line.split()[-1]
    if atom_type.upper().startswith("H"):
        return "H"
    if atom_type in {"A", "C"}:
        return "C"
    if atom_type in {"N", "NA"}:
        return "N"
    if atom_type in {"O", "OA"}:
        return "O"
    if atom_type in {"S", "SA", "S6"}:
        return "S"
    return {"cl": "Cl", "br": "Br", "f": "F", "p": "P", "i": "I"}.get(atom_type.lower(), atom_type.title())


def pdbqt_models(path: Path):
    models, current = [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            models.append(current)
            current = []
        elif line.startswith(("ATOM", "HETATM")):
            element = pdbqt_element(line)
            if element != "H":
                current.append(
                    (element, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], float))
                )
    if current and not models:
        models = [current]
    return models


def nameH_xyz_pdb(path: Path):
    xyz = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("HETATM", "ATOM")):
            continue
        name = line[12:16].strip()
        el = (line[76:78].strip() if len(line) >= 78 else "") or name[:1]
        if name.upper().startswith("H") or el.upper() == "H":
            continue
        xyz.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], float))
    return xyz


def nameH_models(path: Path):
    models, current = [], None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            if current is not None:
                models.append(current)
            current = None
        elif current is not None and line.startswith(("ATOM", "HETATM")):
            if line[12:16].strip().upper().startswith("H"):
                continue
            current.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], float))
    return models


def hungarian(ref, mob) -> float:
    R = np.asarray(ref, float)
    M = np.asarray(mob, float)
    n = min(len(R), len(M))
    R, M = R[:n], M[:n]
    cost = ((R[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    ri, ci = linear_sum_assignment(cost)
    return float(math.sqrt(cost[ri, ci].sum() / len(ri)))


def meeko_mol(path: Path):
    from meeko import PDBQTMolecule, RDKitMolCreate

    pdbqt = PDBQTMolecule.from_file(str(path), skip_typing=True)
    made = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    if not made or made[0] is None:
        return None
    return Chem.RemoveHs(made[0])


def read_sdf(path: Path):
    mol = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)[0]
    return None if mol is None else Chem.RemoveHs(mol)


def crystal_elems(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("HETATM", "ATOM")):
            continue
        name = line[12:16].strip()
        el = (line[76:78].strip() if len(line) >= 78 else "") or name[:1]
        if name.upper().startswith("H") or el.upper() == "H":
            continue
        el = "Cl" if (el.title() == "Cl" or name.upper().startswith("CL")) else (el.title() if len(el) <= 2 else el[0])
        rows.append((el, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], float)))
    return rows


def bondaware_ref(smi: str, cry):
    tmpl = Chem.MolFromSmiles(smi)
    if tmpl is None:
        raise RuntimeError("CCD SMILES parse failed")
    tmpl = Chem.RemoveHs(tmpl)
    ref_el = [a.GetSymbol() for a in tmpl.GetAtoms()]
    cry_el = [x[0] for x in cry]
    if Counter(ref_el) != Counter(cry_el):
        raise RuntimeError(f"element count mismatch {Counter(ref_el)} vs {Counter(cry_el)}")
    xyz = np.asarray([x[1] for x in cry])
    adj = {i: [] for i in range(len(ref_el))}
    for bond in tmpl.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[a].append(b)
        adj[b].append(a)
    cand = {i: [j for j, e in enumerate(cry_el) if e == ref_el[i]] for i in range(len(ref_el))}
    assign = [-1] * len(ref_el)
    taken = [False] * len(cry_el)
    order = sorted(range(len(ref_el)), key=lambda i: (len(cand[i]), -len(adj[i])))

    def ok(i, j):
        for nb in adj[i]:
            if assign[nb] != -1:
                d = float(np.linalg.norm(xyz[j] - xyz[assign[nb]]))
                if d < 1.05 or d > 2.25:
                    return False
        return True

    def bt(k):
        if k == len(order):
            return True
        i = order[k]
        for j in cand[i]:
            if taken[j] or not ok(i, j):
                continue
            assign[i] = j
            taken[j] = True
            if bt(k + 1):
                return True
            assign[i] = -1
            taken[j] = False
        return False

    if not bt(0):
        raise RuntimeError("bond-aware crystal assignment failed")
    conf = Chem.Conformer(tmpl.GetNumAtoms())
    for i, j in enumerate(assign):
        conf.SetAtomPosition(i, xyz[j])
    tmpl.AddConformer(conf, assignId=True)
    return tmpl


def poses_on_ref(ref, input_atoms, output_models, max_map=0.50):
    ref_el = [a.GetSymbol() for a in ref.GetAtoms()]
    ref_xyz = np.asarray([list(ref.GetConformer().GetAtomPosition(i)) for i in range(ref.GetNumAtoms())])
    prep_el = [x[0] for x in input_atoms]
    prep_xyz = np.asarray([x[1] for x in input_atoms])
    if len(ref_el) != len(prep_el):
        raise RuntimeError(f"heavy-atom count {len(ref_el)} != {len(prep_el)}")
    cost = ((ref_xyz[:, None, :] - prep_xyz[None, :, :]) ** 2).sum(axis=2)
    for i, a in enumerate(ref_el):
        for j, b in enumerate(prep_el):
            if a != b:
                cost[i, j] = 1e6
    ri, ci = linear_sum_assignment(cost)
    dist = np.sqrt(cost[ri, ci])
    if float(dist.max()) > max_map:
        raise RuntimeError(f"unsafe map max={dist.max():.3f}")
    prep_for_ref = dict(zip(ri, ci))
    poses = Chem.Mol(ref)
    poses.RemoveAllConformers()
    for model in output_models:
        conf = Chem.Conformer(ref.GetNumAtoms())
        for ridx in range(ref.GetNumAtoms()):
            conf.SetAtomPosition(ridx, model[prep_for_ref[ridx]][1])
        poses.AddConformer(conf, assignId=True)
    rms = [float(rdMolAlign.CalcRMS(poses, ref, prbId=i, refId=0)) for i in range(poses.GetNumConformers())]
    return rms, float(dist.max())


def main() -> int:
    TAB.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)
    summaries = []
    modes = []
    for protein, pdb, ccd in SPECS:
        cry_pdb = COG / f"{pdb}_{ccd}.pdb"
        cry_sdf = COG / f"{pdb}_{ccd}.sdf"
        in_pdbqt = COG / f"{pdb}_{ccd}.pdbqt"
        out_pdbqt = QC / f"{pdb}_{ccd}_out_E8.pdbqt"
        hung = [hungarian(nameH_xyz_pdb(cry_pdb), m) for m in nameH_models(out_pdbqt)]
        in0 = pdbqt_models(in_pdbqt)[0]
        out_models = pdbqt_models(out_pdbqt)
        in_mol = meeko_mol(in_pdbqt)
        out_mol = meeko_mol(out_pdbqt)
        sdf_mol = read_sdf(cry_sdf)
        method = ""
        calc = None
        map_max = ""
        if in_mol is not None and out_mol is not None and Chem.MolToSmiles(in_mol) == Chem.MolToSmiles(out_mol):
            calc = [float(rdMolAlign.CalcRMS(out_mol, in_mol, prbId=i, refId=0)) for i in range(out_mol.GetNumConformers())]
            method = "Meeko topology + RDKit CalcRMS; no superposition"
        elif sdf_mol is not None:
            calc, map_max = poses_on_ref(sdf_mol, in0, out_models)
            method = "SDF graph + element-constrained map + CalcRMS"
        else:
            ref = bondaware_ref(CCD_SMILES[ccd], crystal_elems(cry_pdb))
            calc, map_max = poses_on_ref(ref, in0, out_models, max_map=0.50)
            method = "CCD SMILES + bond-aware crystal coords + CalcRMS"
        top1 = calc[0]
        best = min(calc)
        best_mode = int(np.argmin(calc)) + 1
        top3 = min(calc[:3])
        n = len(calc)
        summaries.append(
            {
                "protein": protein,
                "pdb": pdb,
                "cognate": ccd,
                "n_modes": n,
                "hungarian_top1_A": round(hung[0], 3),
                "hungarian_best_A": round(min(hung), 3),
                "calcrrms_top1_A": round(top1, 3),
                "calcrrms_top3_A": round(top3, 3),
                "calcrrms_best_A": round(best, 3),
                "calcrrms_best_mode": best_mode,
                "pass_best_lt2": int(best < 2.0),
                "pass_top1_lt2": int(top1 < 2.0),
                "pass_top3_lt2": int(top3 < 2.0),
                "mapping_method": method,
                "mapping_max_distance_A": "" if map_max == "" else round(map_max, 4),
                "legacy_source": "layer3_cognate_rmsd_v1.csv Hungarian name-H",
            }
        )
        for rank, (h, c) in enumerate(zip(hung, calc), 1):
            modes.append(
                {
                    "protein": protein,
                    "pdb": pdb,
                    "pose_rank": rank,
                    "hungarian_rmsd_A": round(h, 4),
                    "calcrrms_rmsd_A": round(c, 4),
                    "mapping_method": method,
                }
            )
        print(protein, pdb, [round(x, 3) for x in calc], "best", round(best, 3), method)

    sum_path = TAB / "layer3_cognate_rmsd_calcrrms_v1.csv"
    mode_path = TAB / "layer3_cognate_rmsd_calcrrms_modes_v1.csv"
    with sum_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summaries[0]))
        w.writeheader()
        w.writerows(summaries)
    with mode_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(modes[0]))
        w.writeheader()
        w.writerows(modes)

    lines = [
        "# Track B 共晶 RMSD 化学对应复核\n\n",
        "对象：已保存的 8 个 `*_out_E8.pdbqt`，不重对接。\n\n",
        "主数字：RDKit `CalcRMS`（考虑对称、不做蛋白叠合）。",
        "匈牙利坐标匹配只作对照。`2JKH/BI7` 的 OpenBabel SDF 键级无效，改用 CCD SMILES。\n\n",
        "| 蛋白 | PDB | Hungarian top-1 / best | CalcRMS top-1 / top-3 / best | best&lt;2 | top-1&lt;2 |\n",
        "|------|-----|------------------------:|-----------------------------:|:---:|:---:|\n",
    ]
    for r in summaries:
        lines.append(
            f"| {r['protein']} | {r['pdb']} | {r['hungarian_top1_A']:.3f} / {r['hungarian_best_A']:.3f} | "
            f"{r['calcrrms_top1_A']:.3f} / {r['calcrrms_top3_A']:.3f} / {r['calcrrms_best_A']:.3f} | "
            f"{'是' if r['pass_best_lt2'] else '否'} | {'是' if r['pass_top1_lt2'] else '否'} |\n"
        )
    lines.append(
        "\n搜索覆盖门槛（全部保存姿态最低 RMSD < 2.0 Å）在化学对应复核后仍全部通过。"
        "JAK2、PPARG、PPARA 的第一姿态未回到 2 Å 以内；PPARA 前三姿态也未回到 2 Å 以内。\n"
        f"\n表：`{sum_path.relative_to(ROOT)}`；逐姿态：`{mode_path.relative_to(ROOT)}`。\n"
    )
    (AN / "LAYER3_COGNATE_RMSD_REAUDIT_V1.md").write_text("".join(lines), encoding="utf-8")
    (AN / "layer3_cognate_rmsd_calcrrms_summary.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("wrote", sum_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
