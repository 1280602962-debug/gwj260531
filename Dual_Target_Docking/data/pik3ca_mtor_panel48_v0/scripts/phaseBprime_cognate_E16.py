#!/usr/bin/env python3
"""Phase B′: cognate E=16 for PM48_01 (4 jobs) + RMSD + verdict."""
from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0")
VINA = Path("/home/gwj/miniconda3/bin/vina")
LIG = "PM48_01"
EXHAUST = 16
N_MODES = 9
ENERGY_RANGE = 3
JOBS = [
    ("B16_A", "4L23", 20260727),
    ("B16_B", "4JT6", 20260727),
    ("B16_C", "4JT6", 7),
    ("B16_D", "4JT6", 42),
]
TPL = Chem.MolFromSmiles("Oc1cccc(-c2nc(N3CCOCC3)c3oc4ncccc4c3n2)c1")
SYMS = TPL.GetSubstructMatches(TPL, uniquify=False)


def split_models(all_pose: Path, pose_dir: Path) -> int:
    text = all_pose.read_text(errors="ignore")
    blocks, cur, inn = [], [], False
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            if cur:
                blocks.append(cur)
            cur = [line]
            inn = True
        elif line.startswith("ENDMDL"):
            cur.append(line)
            blocks.append(cur)
            cur, inn = [], False
        elif inn:
            cur.append(line)
    if cur:
        blocks.append(cur)
    if not blocks and text.strip():
        blocks = [[text]]
    pose_dir.mkdir(parents=True, exist_ok=True)
    for i, b in enumerate(blocks, 1):
        (pose_dir / f"mode_{i:02d}.pdbqt").write_text("".join(b))
    return len(blocks)


def run_job(job_id: str, target: str, seed: int) -> Path:
    boxes = json.loads((ROOT / "boxes" / "all_boxes.json").read_text())
    box = boxes[target]
    pose_dir = ROOT / "poses" / "cognate_E16" / target / f"{LIG}_seed{seed}"
    pose_dir.mkdir(parents=True, exist_ok=True)
    all_pose = pose_dir / f"{LIG}_all_modes.pdbqt"
    conf_dir = ROOT / "logs" / "cognate_E16" / "confs"
    conf_dir.mkdir(parents=True, exist_ok=True)
    log_path = ROOT / "logs" / "cognate_E16" / f"{job_id}.log"
    conf = conf_dir / f"{job_id}.txt"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {ROOT / 'receptors' / f'{target}_receptor.pdbqt'}",
                f"ligand = {ROOT / 'ligands_pdbqt' / f'{LIG}.pdbqt'}",
                f"out = {all_pose}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                f"energy_range = {ENERGY_RANGE}",
                "cpu = 1",
                f"seed = {seed}",
            ]
        )
        + "\n"
    )
    with log_path.open("w") as fh:
        proc = subprocess.run([str(VINA), "--config", str(conf)], stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not all_pose.exists():
        raise RuntimeError(f"{job_id} failed rc={proc.returncode}\n{log_path.read_text()[-800:]}")
    n = split_models(all_pose, pose_dir)
    log_txt = log_path.read_text(errors="ignore")
    m = re.search(r"random seed:\s*(-?\d+)", log_txt)
    if not m or int(m.group(1)) != seed:
        raise RuntimeError(f"{job_id} seed mismatch: expected {seed}, got {m.group(1) if m else None}")
    if n != 9:
        raise RuntimeError(f"{job_id} expected 9 modes, got {n}")
    print(f"{job_id} OK target={target} seed={seed} modes={n}", flush=True)
    return pose_dir


def parse_smiles_idx(path: Path):
    nums = []
    for line in path.read_text().splitlines():
        if line.startswith("REMARK SMILES IDX"):
            nums.extend(int(x) for x in line.split()[3:])
    return list(zip(nums[0::2], nums[1::2]))


def parse_coords(path: Path):
    atoms = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms[int(line[6:11])] = (
            line.split()[-1],
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        )
    return atoms


def pose_mol(path: Path) -> Chem.Mol:
    pairs = parse_smiles_idx(path)
    atoms = parse_coords(path)
    m = Chem.Mol(TPL)
    conf = Chem.Conformer(m.GetNumAtoms())
    n = 0
    for s1, p1 in pairs:
        e, xyz = atoms[p1]
        if e in ("H", "HD", "HS"):
            continue
        conf.SetAtomPosition(s1 - 1, xyz.tolist())
        n += 1
    if n != TPL.GetNumAtoms():
        raise RuntimeError(f"mapped {n} heavy atoms in {path}")
    m.RemoveAllConformers()
    m.AddConformer(conf, assignId=True)
    return m


def crystal_mols(pdb: Path):
    lines = [
        l
        for l in pdb.read_text().splitlines()
        if l.startswith(("HETATM", "ATOM  ")) and l[17:20].strip() == "X6K"
    ]
    mol = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", removeHs=True, sanitize=False)
    mol = AllChem.AssignBondOrdersFromTemplate(TPL, mol)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    out = []
    src = mol.GetConformer()
    for match in mol.GetSubstructMatches(TPL, uniquify=False):
        m = Chem.Mol(TPL)
        conf = Chem.Conformer(m.GetNumAtoms())
        for i, j in enumerate(match):
            p = src.GetAtomPosition(j)
            conf.SetAtomPosition(i, [p.x, p.y, p.z])
        m.RemoveAllConformers()
        m.AddConformer(conf, assignId=True)
        out.append(m)
    return out


def rmsd_no_align(ref: Chem.Mol, pose: Chem.Mol) -> float:
    best = None
    for perm in SYMS:
        amap = [(i, perm[i]) for i in range(len(perm))]
        r = float(rdMolAlign.CalcRMS(pose, ref, map=[amap]))
        if best is None or r < best:
            best = r
    return float(best)


def compute_rmsd_table(pose_dirs: dict) -> list[dict]:
    rows = []
    refs = {
        "4L23": crystal_mols(ROOT / "tables" / "4L23_cocrystal_X6K.pdb"),
        "4JT6": crystal_mols(ROOT / "tables" / "4JT6_cocrystal_X6K.pdb"),
    }
    for (job_id, target, seed), pose_dir in pose_dirs.items():
        vals = []
        for mode in range(1, 10):
            pose = pose_mol(pose_dir / f"mode_{mode:02d}.pdbqt")
            best = min(rmsd_no_align(ref, pose) for ref in refs[target])
            vals.append((mode, best))
        mode1 = next(r for m, r in vals if m == 1)
        bm, br = min(vals, key=lambda x: x[1])
        rows.append(
            {
                "target": target,
                "seed": seed,
                "exhaustiveness": EXHAUST,
                "rmsd_mode1": round(mode1, 3),
                "rmsd_best_of_9": round(br, 3),
                "best_of_9_mode": bm,
                "pass_mode1_lt2": mode1 < 2.0,
                "pass_best_of_9_lt2": br < 2.0,
                "job_id": job_id,
            }
        )
        print(rows[-1], flush=True)
    out = ROOT / "analysis" / "cognate_redock_v0" / "tables" / "pm48_01_rmsd_E16.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "target",
                "seed",
                "exhaustiveness",
                "rmsd_mode1",
                "rmsd_best_of_9",
                "best_of_9_mode",
                "pass_mode1_lt2",
                "pass_best_of_9_lt2",
                "job_id",
            ],
        )
        w.writeheader()
        w.writerows(rows)
    return rows


def write_verdict(rows: list[dict]) -> str:
    r_l23 = next(r for r in rows if r["target"] == "4L23" and r["seed"] == 20260727)
    r_jt6_main = next(r for r in rows if r["target"] == "4JT6" and r["seed"] == 20260727)
    noise = [r for r in rows if r["target"] == "4JT6"]
    n_pass = sum(1 for r in noise if r["pass_best_of_9_lt2"])
    go = (
        r_l23["pass_best_of_9_lt2"]
        and r_jt6_main["pass_best_of_9_lt2"]
        and n_pass >= 2
    )
    verdict = "Go" if go else "No-Go"
    notes = []
    if go and (not r_l23["pass_mode1_lt2"] or not r_jt6_main["pass_mode1_lt2"]):
        notes.append(
            "mode1 失败但 best_of_9 成功：全面板必须输出 9 modes，并计划 RTM best-of-9。"
        )
    if go and n_pass == 2:
        notes.append("噪声层恰好 2/3 通过（含主 seed）。")
    if not go:
        notes.append("主 seed 或噪声条件未满足；禁止开全面板。")

    lines = [
        "# COGNATE QC VERDICT E16 — PM48_01 (PI-103 / X6K)",
        "",
        f"**Verdict: {verdict}**",
        "",
        "## Go 条件核对",
        f"1. 4L23 @ 20260727 best9={r_l23['rmsd_best_of_9']} (<2? {r_l23['pass_best_of_9_lt2']})",
        f"2. 4JT6 @ 20260727 best9={r_jt6_main['rmsd_best_of_9']} (<2? {r_jt6_main['pass_best_of_9_lt2']})",
        f"3. 4JT6 噪声 seeds{{20260727,7,42}} pass count = {n_pass}/3（需≥2）",
        "",
        "## 全表",
        "",
        "| target | seed | rmsd_mode1 | rmsd_best_of_9 | best_mode | mode1<2 | best9<2 |",
        "|--------|------|------------|----------------|-----------|---------|---------|",
    ]
    for r in rows:
        lines.append(
            f"| {r['target']} | {r['seed']} | {r['rmsd_mode1']} | {r['rmsd_best_of_9']} | {r['best_of_9_mode']} | {r['pass_mode1_lt2']} | {r['pass_best_of_9_lt2']} |"
        )
    lines += ["", "## 备注", ""]
    lines.extend(f"- {n}" for n in notes or ["无额外备注。"])
    lines += [
        "",
        "## 参数",
        "- exhaustiveness: 16",
        "- n_modes: 9",
        "- RMSD: heavy-atom, meeko SMILES IDX map, template automorphism min CalcRMS, no superposition",
        "- poses: `poses/cognate_E16/`（未覆盖 E=8）",
        "- table: `analysis/cognate_redock_v0/tables/pm48_01_rmsd_E16.csv`",
        "",
    ]
    if go:
        lines.append("**允许启动全面板 48×2 @ E=16 @ seed=20260727。**")
    else:
        lines.append("**禁止启动全面板。**")
    (ROOT / "analysis" / "cognate_redock_v0" / "COGNATE_QC_VERDICT_E16.md").write_text(
        "\n".join(lines) + "\n"
    )
    return verdict


def main():
    (ROOT / "poses" / "cognate_E16").mkdir(parents=True, exist_ok=True)
    (ROOT / "logs" / "cognate_E16").mkdir(parents=True, exist_ok=True)
    pose_dirs = {}
    for job_id, target, seed in JOBS:
        pose_dirs[(job_id, target, seed)] = run_job(job_id, target, seed)
    rows = compute_rmsd_table(pose_dirs)
    verdict = write_verdict(rows)
    print("VERDICT", verdict, flush=True)


if __name__ == "__main__":
    main()
