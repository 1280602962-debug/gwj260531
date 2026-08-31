#!/usr/bin/env python3
"""Phase B: PM48_01 cognate redock on 4L23 and 4JT6 + RMSD QC."""
from __future__ import annotations

import csv
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0")
VINA = Path("/home/gwj/miniconda3/bin/vina")
SEED = 20260727
EXHAUST = 8
N_MODES = 9
ENERGY_RANGE = 3
TARGETS = ["4L23", "4JT6"]
LIG = "PM48_01"


def split_models(all_pose: Path, pose_dir: Path) -> int:
    text = all_pose.read_text(errors="ignore")
    blocks = []
    current = []
    in_model = False
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            if current:
                blocks.append(current)
            current = [line]
            in_model = True
        elif line.startswith("ENDMDL"):
            current.append(line)
            blocks.append(current)
            current = []
            in_model = False
        elif in_model:
            current.append(line)
    if current:
        blocks.append(current)
    if not blocks and text.strip():
        blocks = [[text]]
    pose_dir.mkdir(parents=True, exist_ok=True)
    for i, block in enumerate(blocks, 1):
        (pose_dir / f"mode_{i:02d}.pdbqt").write_text("".join(block))
    return len(blocks)


def run_vina(target: str) -> None:
    boxes = json.loads((ROOT / "boxes" / "all_boxes.json").read_text())
    box = boxes[target]
    pose_dir = ROOT / "poses" / target / LIG
    pose_dir.mkdir(parents=True, exist_ok=True)
    all_pose = pose_dir / f"{LIG}_all_modes.pdbqt"
    log = ROOT / "logs" / "vina" / f"{target}_{LIG}.log"
    conf = ROOT / "logs" / "vina_confs" / f"{target}_{LIG}.txt"
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
                f"seed = {SEED}",
                "cpu = 1",
            ]
        )
        + "\n"
    )
    cmd = [
        str(VINA),
        "--receptor",
        str(ROOT / "receptors" / f"{target}_receptor.pdbqt"),
        "--ligand",
        str(ROOT / "ligands_pdbqt" / f"{LIG}.pdbqt"),
        "--out",
        str(all_pose),
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--exhaustiveness",
        str(EXHAUST),
        "--num_modes",
        str(N_MODES),
        "--energy_range",
        str(ENERGY_RANGE),
        "--seed",
        str(SEED),
        "--cpu",
        "1",
    ]
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not all_pose.exists():
        raise RuntimeError(f"vina failed {target}: rc={proc.returncode}")
    n = split_models(all_pose, pose_dir)
    log_txt = log.read_text(errors="ignore")
    if "20260727" not in log_txt:
        raise RuntimeError(f"seed 20260727 not found in {log}")
    print(f"{target}: ok modes={n} seed_ok")


def crystal_mol(pdb_path: Path, tmpl: Chem.Mol) -> Chem.Mol:
    lines = [
        line
        for line in pdb_path.read_text(errors="ignore").splitlines()
        if line.startswith(("HETATM", "ATOM  ")) and line[17:20].strip() == "X6K"
    ]
    mol = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", removeHs=True, sanitize=False)
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol)


def pose_mol(pdbqt_path: Path, tmpl: Chem.Mol) -> Chem.Mol:
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tf:
        out = tf.name
    subprocess.check_call(
        ["bash", "-lc", f"source /home/gwj/miniconda3/etc/profile.d/conda.sh && conda activate cadd_tools && obabel -ipdbqt '{pdbqt_path}' -osdf -O '{out}' >/dev/null 2>&1"]
    )
    mol = Chem.SDMolSupplier(out, removeHs=False, sanitize=False)[0]
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, Chem.RemoveHs(mol))
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol)


def symmetry_corrected_rmsd(ref: Chem.Mol, pose: Chem.Mol, tmpl: Chem.Mol) -> float:
    matches_ref = ref.GetSubstructMatches(tmpl, uniquify=False)
    matches_pose = pose.GetSubstructMatches(tmpl, uniquify=False)
    if not matches_ref or not matches_pose:
        raise RuntimeError("no template substructure match")
    best = None
    for a in matches_ref:
        for b in matches_pose:
            rms = float(rdMolAlign.CalcRMS(pose, ref, map=[list(zip(b, a))]))
            if best is None or rms < best:
                best = rms
    return float(best)


def compute_rmsd_and_verdict() -> None:
    tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(ROOT / "ligands_sdf" / f"{LIG}.sdf"), removeHs=False)[0])
    rows = []
    for target in TARGETS:
        ref = crystal_mol(ROOT / "tables" / f"{target}_cocrystal_X6K.pdb", tmpl)
        rmsds = []
        for mode in range(1, 10):
            pose = pose_mol(ROOT / "poses" / target / LIG / f"mode_{mode:02d}.pdbqt", tmpl)
            rms = symmetry_corrected_rmsd(ref, pose, tmpl)
            rmsds.append((mode, rms))
        mode1 = next(r for m, r in rmsds if m == 1)
        best_mode, best_rmsd = min(rmsds, key=lambda x: x[1])
        rows.append(
            {
                "target": target,
                "seed": SEED,
                "exhaustiveness": EXHAUST,
                "rmsd_mode1": round(mode1, 3),
                "rmsd_best_of_9": round(best_rmsd, 3),
                "best_of_9_mode": best_mode,
                "pass_mode1_lt2": mode1 < 2.0,
                "pass_best_of_9_lt2": best_rmsd < 2.0,
                "rmsd_rtm_best": "",
                "rtm_best_mode": "",
            }
        )
        print(target, "mode1", mode1, "best", best_rmsd, "mode", best_mode)

    out_csv = ROOT / "analysis" / "cognate_redock_v0" / "tables" / "pm48_01_rmsd.csv"
    with out_csv.open("w", newline="") as fh:
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
                "rmsd_rtm_best",
                "rtm_best_mode",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    both_best = all(r["pass_best_of_9_lt2"] for r in rows)
    both_mode1 = all(r["pass_mode1_lt2"] for r in rows)
    if both_best:
        decision = "Go"
        detail = "两端 best_of_9 < 2.0 Å；采样 QC 通过，可进入 Phase C 全面板。"
        if not both_mode1:
            detail += " 注意：至少一端 mode1 失败但 best_of_9 成功（与 EGFR/TAK-285 在 3POZ 上类似）；全面板可开，需在报告写明并计划 RTM best-of-9。"
    else:
        failed = [r["target"] for r in rows if not r["pass_best_of_9_lt2"]]
        # check if any near-miss
        near = [r for r in rows if (not r["pass_best_of_9_lt2"]) and r["rmsd_best_of_9"] < 3.0]
        if near:
            decision = "No-Go（暂缓）"
            detail = f"端 {failed} best_of_9 ≥ 2.0 但接近失败；先查盒子/质子化/配体映射；可仅对 PM48_01 试 E=16 诊断，不自动全面板升 E。"
        else:
            decision = "No-Go"
            detail = f"端 {failed} 采样失败（mode1 与 best_of_9 均未过或远超阈值）；停全面板，修蛋白/盒子/配体后重做 Phase B。"

    md = ROOT / "analysis" / "cognate_redock_v0" / "COGNATE_QC_VERDICT.md"
    lines = [
        "# COGNATE QC VERDICT — PM48_01 (PI-103 / X6K)",
        "",
        f"**结论：{decision}**",
        "",
        detail,
        "",
        "## 参数",
        f"- seed: `{SEED}`",
        f"- exhaustiveness: `{EXHAUST}`",
        f"- n_modes: `{N_MODES}`",
        f"- RMSD: heavy-atom, template-constrained min CalcRMS vs `tables/<target>_cocrystal_X6K.pdb`",
        "",
        "## 结果",
        "",
        "| target | rmsd_mode1 | rmsd_best_of_9 | best_mode | mode1<2 | best9<2 |",
        "|--------|------------|----------------|-----------|---------|---------|",
    ]
    for r in rows:
        lines.append(
            f"| {r['target']} | {r['rmsd_mode1']} | {r['rmsd_best_of_9']} | {r['best_of_9_mode']} | {r['pass_mode1_lt2']} | {r['pass_best_of_9_lt2']} |"
        )
    lines += ["", f"表文件：`analysis/cognate_redock_v0/tables/pm48_01_rmsd.csv`", ""]
    md.write_text("\n".join(lines))
    print("VERDICT", decision)
    (ROOT / "analysis" / "cognate_redock_v0" / "tables" / "rmsd_definition.md").write_text(
        """# RMSD definition — cognate redock v0\n\n- Reference: `tables/4L23_cocrystal_X6K.pdb`, `tables/4JT6_cocrystal_X6K.pdb`\n- Atoms: heavy atoms only\n- Symmetry: template-constrained min CalcRMS using PM48_01 SDF as chemistry template (AssignBondOrdersFromTemplate)\n- No protein superposition; docking frame coordinates\n- Threshold: best_of_9 < 2.0 Å on BOTH ends for Go\n"""
    )


def main():
    for t in TARGETS:
        print("Docking", t, LIG)
        run_vina(t)
    print("Computing RMSD...")
    compute_rmsd_and_verdict()


if __name__ == "__main__":
    main()
