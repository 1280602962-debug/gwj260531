#!/usr/bin/env python3
"""Split corrected-box Vina 9-mode outputs and run RTMScore best-of-9.

Does not reuse old-box poses. Writes under remediation_outputs/phase_rtm/.
"""
from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
VINA_OUT = ROOT / "remediation_outputs/phase1_vina/out"
POSE = ROOT / "remediation_outputs/phase_rtm/poses"
LOGS = ROOT / "remediation_outputs/phase_rtm/logs/rtmscore"
PANEL = ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"
LIG_SDF = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0/ligands_sdf"
)
RECDIR = ROOT / "data/egfr_her2_panel40_v0/receptors"
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
TARGETS = ["3POZ", "3RCD"]


def split_modes(out_pdbqt: Path, dest_dir: Path) -> int:
    text = out_pdbqt.read_text(errors="ignore")
    blocks = re.split(r"(?=^MODEL\s+\d+)", text, flags=re.M)
    dest_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for block in blocks:
        m = re.match(r"MODEL\s+(\d+)", block.strip())
        if not m:
            continue
        mode = int(m.group(1))
        path = dest_dir / f"mode_{mode:02d}.pdbqt"
        body = block if block.rstrip().endswith("ENDMDL") else block.rstrip() + "\nENDMDL\n"
        path.write_text(body)
        n += 1
    return n


def pdbqt_xyz(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return xyz


def smiles_idx_pairs(path: Path):
    nums = []
    for line in path.read_text().splitlines():
        if line.startswith("REMARK SMILES IDX"):
            nums.extend(int(x) for x in line.split()[3:])
    if not nums:
        return None
    return list(zip(nums[0::2], nums[1::2]))


def split_all(panel: pd.DataFrame) -> None:
    n_ok = 0
    for lig in panel["panel_id"]:
        for t in TARGETS:
            dest = POSE / t / lig
            if (dest / "mode_01.pdbqt").exists():
                n_ok += 1
                continue
            src = VINA_OUT / f"{t}_{lig}_out.pdbqt"
            if not src.exists():
                print("MISSING", src, flush=True)
                continue
            n = split_modes(src, dest)
            n_ok += 1
            if n == 0:
                print("empty split", src, flush=True)
    print(f"split ligand-target pairs with mode_01: {n_ok}", flush=True)


def write_sdfs(panel: pd.DataFrame) -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    for target in TARGETS:
        out_sdf = LOGS / f"{target}_poses.sdf"
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        skipped = 0
        for lig in panel["panel_id"]:
            sdf_lig = LIG_SDF / f"{lig}.sdf"
            pose_dir = POSE / target / lig
            modes = sorted(pose_dir.glob("mode_*.pdbqt"))
            if not modes or not sdf_lig.exists():
                skipped += 1
                continue
            tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(sdf_lig), removeHs=False)[0])
            if tmpl is None:
                skipped += 1
                continue
            for mp in modes:
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                xyz = pdbqt_xyz(mp)
                pairs = smiles_idx_pairs(mp)
                if not pairs or len(pairs) != tmpl.GetNumAtoms():
                    skipped += 1
                    continue
                mol = Chem.Mol(tmpl)
                conf = Chem.Conformer(mol.GetNumAtoms())
                for s_idx, p_idx in pairs:
                    x, y, z = xyz[p_idx - 1]
                    conf.SetAtomPosition(s_idx - 1, (x, y, z))
                mol.RemoveAllConformers()
                mol.AddConformer(conf, assignId=True)
                mol.SetProp("_Name", f"{lig}_mode{mode}")
                w.write(mol)
                n += 1
        w.close()
        print(f"wrote {out_sdf} n={n} skipped={skipped}", flush=True)


def run_rtm() -> None:
    for target in TARGETS:
        sdf = LOGS / f"{target}_poses.sdf"
        pocket = RECDIR / f"{target}_pocket_10.0.pdb"
        out_prefix = LOGS / f"{target}_rtmscore"
        log = LOGS / f"{target}_rtmscore.log"
        print("RTM", target, flush=True)
        with log.open("w") as fh:
            proc = subprocess.run(
                [
                    str(RTM_PYTHON),
                    str(RTM_PY),
                    "-p",
                    str(pocket),
                    "-l",
                    str(sdf),
                    "-m",
                    str(MODEL),
                    "-o",
                    str(out_prefix),
                ],
                cwd=str(RTM_ROOT / "example"),
                stdout=fh,
                stderr=subprocess.STDOUT,
            )
        if proc.returncode != 0:
            raise RuntimeError(f"RTM failed {target}; see {log}")
        csv_path = Path(f"{out_prefix}.csv")
        if not csv_path.exists():
            alt = RTM_ROOT / "example" / f"{out_prefix.name}.csv"
            if alt.exists():
                alt.rename(csv_path)
        print("OK", csv_path, flush=True)


def best_of_9() -> Path:
    rows = []
    for target in TARGETS:
        path = LOGS / f"{target}_rtmscore.csv"
        df = pd.read_csv(path)
        sc_col = "rtmscore" if "rtmscore" in df.columns else df.columns[1]
        name_col = df.columns[0]
        for _, r in df.iterrows():
            name = str(r[name_col])
            lig, _, mode = name.rpartition("_mode")
            rows.append(
                {
                    "ligand": lig,
                    "target": target,
                    "vina_mode": int(mode) if str(mode).isdigit() else mode,
                    "rtmscore": float(r[sc_col]),
                }
            )
    long = pd.DataFrame(rows)
    dest_long = ROOT / "remediation_outputs/phase_rtm/rtmscore_long.csv"
    dest_long.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(dest_long, index=False)
    best = long.sort_values(["ligand", "target", "rtmscore"], ascending=[True, True, False]).drop_duplicates(
        ["ligand", "target"]
    )
    a = best[best.target == "3POZ"][["ligand", "rtmscore", "vina_mode"]].rename(
        columns={"rtmscore": "rtm_3POZ", "vina_mode": "rtm_mode_3POZ"}
    )
    b = best[best.target == "3RCD"][["ligand", "rtmscore", "vina_mode"]].rename(
        columns={"rtmscore": "rtm_3RCD", "vina_mode": "rtm_mode_3RCD"}
    )
    wide = a.merge(b, on="ligand", how="outer")
    dest = ROOT / "remediation_outputs/phase_rtm/rtmscore_best.csv"
    wide.to_csv(dest, index=False)
    print("wrote", dest, "n=", len(wide), flush=True)
    return dest


def main() -> int:
    panel = pd.read_csv(PANEL)
    split_all(panel)
    write_sdfs(panel)
    run_rtm()
    best_of_9()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
