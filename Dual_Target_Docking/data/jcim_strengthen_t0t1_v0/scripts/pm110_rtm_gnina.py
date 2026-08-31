#!/usr/bin/env python3
"""PM110: split poses → RTM (best-of-9) → GNINA mode_01 → ablation tables."""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel110_rdkit_v0")
R48 = Path("/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0")
REPO = Path("/home/gwj/repos/gwj260531/Dual_Target_Docking")
REPO_OUT = REPO / "data/pik3ca_mtor_panel110_rdkit_v0"
PANEL = REPO_OUT / "tables/panel_v0_110.csv"

RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
GNINA_SCRIPT = REPO / "data/jcim_bench_v0/scripts/gnina_rescore_panel.py"
TARGETS = ["4L23", "4JT6"]


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
        path.write_text(block if block.rstrip().endswith("ENDMDL") else block.rstrip() + "\nENDMDL\n")
        n += 1
    return n


def ensure_poses(panel: pd.DataFrame):
    for _, r in panel.iterrows():
        lig = r["panel_id"]
        for t in TARGETS:
            dest = ROOT / "poses" / t / lig
            if (dest / "mode_01.pdbqt").exists():
                continue
            if lig.startswith("PM48_"):
                src = R48 / "poses" / t / lig
                if src.exists():
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    if dest.exists() or dest.is_symlink():
                        if dest.is_symlink() or dest.is_file():
                            dest.unlink()
                        else:
                            shutil.rmtree(dest)
                    dest.symlink_to(src)
                    continue
            out = ROOT / "logs" / "vina" / f"{t}_{lig}_out.pdbqt"
            if not out.exists():
                print("MISSING out", out, flush=True)
                continue
            n = split_modes(out, dest)
            print(f"split {t}/{lig} modes={n}", flush=True)


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


def vina_from_pdbqt(path: Path):
    for line in path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def write_sdfs(panel: pd.DataFrame):
    logs = ROOT / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in TARGETS:
        out_sdf = logs / f"{target}_poses.sdf"
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        skipped = 0
        for lig in panel["panel_id"]:
            sdf_path = ROOT / "ligands_sdf" / f"{lig}.sdf"
            if not sdf_path.exists() and lig.startswith("PM48_"):
                sdf_path = R48 / "ligands_sdf" / f"{lig}.sdf"
            if not sdf_path.exists():
                skipped += 1
                continue
            tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0])
            pose_dir = ROOT / "poses" / target / lig
            for mp in sorted(pose_dir.glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                xyz = pdbqt_xyz(mp)
                pairs = smiles_idx_pairs(mp)
                if not pairs or len(pairs) != tmpl.GetNumAtoms():
                    print(f"SKIP pairs {mp}", flush=True)
                    skipped += 1
                    continue
                mol = Chem.Mol(tmpl)
                conf = Chem.Conformer(mol.GetNumAtoms())
                for s_idx, p_idx in pairs:
                    conf.SetAtomPosition(s_idx - 1, xyz[p_idx - 1])
                mol.RemoveAllConformers()
                mol.AddConformer(conf, assignId=True)
                mol.SetProp("_Name", f"{lig}_mode{mode}")
                w.write(mol)
                n += 1
        w.close()
        print(f"wrote {out_sdf} n={n} skipped={skipped}", flush=True)


def run_rtm():
    logs = ROOT / "logs" / "rtmscore"
    for target in TARGETS:
        sdf = logs / f"{target}_poses.sdf"
        pocket = ROOT / "receptors" / f"{target}_pocket_10.0.pdb"
        out_prefix = logs / f"{target}_rtmscore"
        log = logs / f"{target}_rtmscore.log"
        print("RTM", target, flush=True)
        with log.open("w") as fh:
            proc = subprocess.run(
                [
                    str(RTM_PYTHON), str(RTM_PY),
                    "-p", str(pocket), "-l", str(sdf),
                    "-m", str(MODEL), "-o", str(out_prefix),
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


def run_gnina():
    cmd = [
        sys.executable, str(GNINA_SCRIPT),
        "--root", str(ROOT),
        "--targets", *TARGETS,
        "--receptor-map",
        f"4L23={ROOT}/receptors/4L23_protein.pdb",
        f"4JT6={ROOT}/receptors/4JT6_protein.pdb",
        "--workers", "4",
        "--modes", "mode_01",
        "--timeout", "300",
    ]
    print("GNINA", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError("GNINA failed")


def parse_id(s):
    m = re.search(r"(PM(?:48|110)_\d+)_mode(\d+)", str(s))
    if not m:
        raise ValueError(s)
    return m.group(1), int(m.group(2))


def build_tables(panel: pd.DataFrame):
    ids = panel["panel_id"].tolist()
    vina_rows = []
    for lig in ids:
        for t in TARGETS:
            for mp in sorted((ROOT / "poses" / t / lig).glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                vina_rows.append({
                    "ligand": lig, "target": t, "vina_mode": mode,
                    "vina_score": vina_from_pdbqt(mp),
                })
    vina_long = pd.DataFrame(vina_rows)
    vina_long.to_csv(ROOT / "tables" / "scores_vina_long.csv", index=False)
    m1 = vina_long.loc[vina_long["vina_mode"] == 1].pivot(
        index="ligand", columns="target", values="vina_score"
    )

    rtm_rows = []
    for t in TARGETS:
        d = pd.read_csv(ROOT / "logs" / "rtmscore" / f"{t}_rtmscore.csv")
        id_col = "id" if "id" in d.columns else d.columns[0]
        sc_col = "score" if "score" in d.columns else d.columns[1]
        for _, r in d.iterrows():
            lig, mode = parse_id(r[id_col])
            rtm_rows.append({
                "ligand": lig, "target": t, "vina_mode": mode,
                "rtmscore": float(r[sc_col]),
            })
    rtm_long = pd.DataFrame(rtm_rows)
    rtm_long.to_csv(ROOT / "tables" / "scores_rtm_all_poses.csv", index=False)
    best = (
        rtm_long.sort_values(["ligand", "target", "rtmscore"], ascending=[True, True, False])
        .groupby(["ligand", "target"], as_index=False).first()
    )
    a = best[best.target == "4L23"][["ligand", "rtmscore"]].rename(columns={"rtmscore": "rtm_4L23"})
    b = best[best.target == "4JT6"][["ligand", "rtmscore"]].rename(columns={"rtmscore": "rtm_4JT6"})

    df = panel.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "molecule_chembl_id", "smiles", "pchembl_PIK3CA", "pchembl_MTOR"]
    ]
    df = df.merge(m1.reset_index(), on="ligand").merge(a, on="ligand").merge(b, on="ligand")
    df = df.rename(columns={"4L23": "4L23_affinity", "4JT6": "4JT6_affinity"})

    gnina_best = ROOT / "tables" / "scores_gnina_best.csv"
    if gnina_best.exists():
        g = pd.read_csv(gnina_best)
        df = df.merge(g, on="ligand", how="left")

    out = ROOT / "tables" / "ablation_ligand_scores.csv"
    df.to_csv(out, index=False)
    REPO_OUT.mkdir(parents=True, exist_ok=True)
    for name in (
        "ablation_ligand_scores.csv",
        "scores_vina_long.csv",
        "scores_rtm_all_poses.csv",
        "scores_gnina_best.csv",
        "scores_gnina_long.csv",
        "scores_vina_best.csv",
    ):
        src = ROOT / "tables" / name
        if src.exists():
            shutil.copy2(src, REPO_OUT / "tables" / name)
    print("wrote", out, "n=", len(df), flush=True)


def main():
    panel = pd.read_csv(PANEL)
    print("panel", len(panel), flush=True)
    ensure_poses(panel)
    write_sdfs(panel)
    run_rtm()
    run_gnina()
    build_tables(panel)
    return 0


if __name__ == "__main__":
    sys.exit(main())
