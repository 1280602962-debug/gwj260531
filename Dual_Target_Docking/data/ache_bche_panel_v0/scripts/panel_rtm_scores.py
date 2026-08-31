#!/usr/bin/env python3
"""Panel RTM best-of-K + Vina tables (PM48-style SMILES-IDX SDF rebuild)."""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")


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


def auroc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = 0.0
    for p in pos:
        wins += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(wins / (len(pos) * len(neg)))


def write_sdfs(root: Path, targets: list[str], panel: pd.DataFrame):
    logs = root / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in targets:
        out_sdf = logs / f"{target}_poses.sdf"
        if out_sdf.exists() and out_sdf.stat().st_size > 1000:
            print("reuse", out_sdf, flush=True)
            continue
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        for lig in panel["panel_id"]:
            sdf_lig = root / "ligands_sdf" / f"{lig}.sdf"
            pose_dir = root / "poses" / target / lig
            modes = sorted(pose_dir.glob("mode_*.pdbqt"))
            if not modes or not sdf_lig.exists():
                continue
            tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(sdf_lig), removeHs=False)[0])
            if tmpl is None:
                continue
            for mp in modes:
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                xyz = pdbqt_xyz(mp)
                pairs = smiles_idx_pairs(mp)
                if not pairs or len(pairs) != tmpl.GetNumAtoms():
                    print(f"skip pairs {mp}", flush=True)
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
        print("wrote", out_sdf, n, flush=True)


def ensure_pocket(root: Path, target: str, protein: Path) -> Path:
    pocket = root / "receptors" / f"{target}_pocket_10.0.pdb"
    if pocket.exists() and pocket.stat().st_size > 100:
        return pocket
    # RTM can generate pocket; for speed use protein as -p with -gen_pocket
    return protein


def run_rtm(root: Path, target: str, protein: Path, reflig: Path | None = None):
    logs = root / "logs" / "rtmscore"
    sdf = logs / f"{target}_poses.sdf"
    out_prefix = logs / f"{target}_rtmscore"
    csv_path = Path(f"{out_prefix}.csv")
    if csv_path.exists() and csv_path.stat().st_size > 100:
        print("reuse", csv_path, flush=True)
        return csv_path
    pocket = root / "receptors" / f"{target}_pocket_10.0.pdb"
    log = logs / f"{target}_rtmscore.log"
    print("RTM", target, flush=True)
    if pocket.exists() and pocket.stat().st_size > 100:
        cmd = [
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
        ]
    else:
        if reflig is None or not Path(reflig).exists():
            raise RuntimeError(f"need pocket or reflig for {target}")
        cmd = [
            str(RTM_PYTHON),
            str(RTM_PY),
            "-p",
            str(protein),
            "-l",
            str(sdf),
            "-m",
            str(MODEL),
            "-o",
            str(out_prefix),
            "-gen_pocket",
            "-c",
            "10.0",
            "-rl",
            str(reflig),
        ]
    with log.open("w") as fh:
        proc = subprocess.run(
            cmd, cwd=str(RTM_ROOT / "example"), stdout=fh, stderr=subprocess.STDOUT
        )
    if proc.returncode != 0:
        raise RuntimeError(f"RTM failed {target}; see {log}")
    # keep generated pocket if produced next to protein
    gen = protein.with_name(protein.stem + "_pocket_10.0.pdb")
    if gen.exists() and not pocket.exists():
        pocket.write_bytes(gen.read_bytes())
    if not csv_path.exists():
        alt = RTM_ROOT / "example" / f"{out_prefix.name}.csv"
        if alt.exists():
            alt.rename(csv_path)
    print("OK", csv_path, flush=True)
    return csv_path


def parse_id(s):
    m = re.search(r"([A-Z]+_\d+)_mode(\d+)", str(s))
    if not m:
        raise ValueError(s)
    return m.group(1), int(m.group(2))


def build_tables(root: Path, targets: list[str], panel: pd.DataFrame, repo: Path | None):
    ids = panel["panel_id"].tolist()
    a, b = targets
    vina_rows = []
    for lig in ids:
        for t in targets:
            for mp in sorted((root / "poses" / t / lig).glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                vina_rows.append(
                    {
                        "ligand": lig,
                        "target": t,
                        "vina_mode": mode,
                        "vina_score": vina_from_pdbqt(mp),
                    }
                )
    vina_long = pd.DataFrame(vina_rows)
    vina_long.to_csv(root / "tables" / "scores_vina_long.csv", index=False)
    m1 = (
        vina_long.loc[vina_long["vina_mode"] == 1]
        .pivot(index="ligand", columns="target", values="vina_score")
        if len(vina_long)
        else pd.DataFrame()
    )

    rtm_rows = []
    for t in targets:
        d = pd.read_csv(root / "logs" / "rtmscore" / f"{t}_rtmscore.csv")
        id_col = "id" if "id" in d.columns else d.columns[0]
        sc_col = "score" if "score" in d.columns else d.columns[1]
        for _, r in d.iterrows():
            lig, mode = parse_id(r[id_col])
            rtm_rows.append(
                {"ligand": lig, "target": t, "vina_mode": mode, "rtmscore": float(r[sc_col])}
            )
    rtm_long = pd.DataFrame(rtm_rows)
    rtm_long.to_csv(root / "tables" / "scores_rtm_all_poses.csv", index=False)
    best = (
        rtm_long.sort_values(["ligand", "target", "rtmscore"], ascending=[True, True, False])
        .groupby(["ligand", "target"], as_index=False)
        .first()
    )
    ra = best[best.target == a][["ligand", "rtmscore"]].rename(columns={"rtmscore": f"rtm_{a}"})
    rb = best[best.target == b][["ligand", "rtmscore"]].rename(columns={"rtmscore": f"rtm_{b}"})

    df = panel.rename(columns={"panel_id": "ligand"}).copy()
    if len(m1):
        df = df.merge(m1.reset_index(), on="ligand", how="left")
        df = df.rename(columns={a: f"vina_{a}", b: f"vina_{b}"})
    df = df.merge(ra, on="ligand", how="left").merge(rb, on="ligand", how="left")
    df[f"vina_{a}_hb"] = -df[f"vina_{a}"].astype(float)
    df[f"vina_{b}_hb"] = -df[f"vina_{b}"].astype(float)
    df["vina_mean"] = (df[f"vina_{a}_hb"] + df[f"vina_{b}_hb"]) / 2
    df["vina_min"] = df[[f"vina_{a}_hb", f"vina_{b}_hb"]].min(axis=1)
    df["rtm_mean"] = (df[f"rtm_{a}"] + df[f"rtm_{b}"]) / 2
    df["rtm_min"] = df[[f"rtm_{a}", f"rtm_{b}"]].min(axis=1)
    for col, zcol in [(f"rtm_{a}", f"rtm_{a}_z"), (f"rtm_{b}", f"rtm_{b}_z")]:
        mu, sd = df[col].mean(), df[col].std(ddof=0)
        df[zcol] = (df[col] - mu) / (sd if sd and sd > 0 else 1.0)
    df["rtm_min_z"] = df[[f"rtm_{a}_z", f"rtm_{b}_z"]].min(axis=1)
    df["prep"] = "rdkit_meeko"
    out = root / "tables" / "ablation_ligand_scores.csv"
    df.to_csv(out, index=False)
    print("wrote", out, flush=True)

    rows = []
    for arm, col in (("vina_mean", "vina_mean"), ("rtm_min_z", "rtm_min_z")):
        d = df.loc[df["class"] == "dual", col].astype(float)
        ao = df.loc[df["class"] == "A_only", col].astype(float)
        bo = df.loc[df["class"] == "B_only", col].astype(float)
        rows.append(
            {
                "arm": arm,
                "auroc_D_vs_A": auroc(d, ao),
                "auroc_D_vs_B": auroc(d, bo),
                "summary_min": float(np.nanmin([auroc(d, ao), auroc(d, bo)])),
                "n_scored": int(df[col].notna().sum()),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(root / "tables" / "directional_summary.csv", index=False)
    print(summary.to_string(index=False), flush=True)

    if repo:
        (repo / "tables").mkdir(parents=True, exist_ok=True)
        for name in (
            "ablation_ligand_scores.csv",
            "directional_summary.csv",
            "scores_vina_long.csv",
            "scores_rtm_all_poses.csv",
            "job_status.csv",
            "panel_v0_strict_with_smiles.csv",
        ):
            src = root / "tables" / name
            if src.exists():
                (repo / "tables" / name).write_bytes(src.read_bytes())
        (repo / "analysis").mkdir(exist_ok=True)
        (repo / "analysis" / "DIRECTIONAL.md").write_text(
            "# Directional AUROC\n\n```\n" + summary.to_string(index=False) + "\n```\n"
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--targets", nargs=2, required=True)
    ap.add_argument(
        "--proteins",
        nargs=2,
        required=True,
        help="protein PDB paths matching --targets order",
    )
    ap.add_argument(
        "--refligs",
        nargs=2,
        default=None,
        help="cognate SDF refs for pocket gen when no precomputed pocket",
    )
    ap.add_argument("--repo", default="")
    args = ap.parse_args()
    root = Path(args.root)
    panel = pd.read_csv(args.panel)
    write_sdfs(root, args.targets, panel)
    refligs = [Path(x) for x in args.refligs] if args.refligs else [None, None]
    for t, prot, rl in zip(args.targets, args.proteins, refligs):
        run_rtm(root, t, Path(prot), rl)
    build_tables(root, args.targets, panel, Path(args.repo) if args.repo else None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
