#!/usr/bin/env python3
"""PM48 RDKit postprocess: RTM + LigPrep delta + directional tables."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_rdkit_v0"
)
OLD = Path(
    "/home/gwj/repos/gwj260531/Dual_Target_Docking/data/pik3ca_mtor_panel48_v0/tables"
)
REPO = Path(
    "/home/gwj/repos/gwj260531/Dual_Target_Docking/data/pik3ca_mtor_panel48_rdkit_v0"
)
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
TARGETS = ["4L23", "4JT6"]
SEED = 20260727


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


def directional(df, arm):
    d = df.loc[df["class"] == "dual", arm].astype(float)
    a = df.loc[df["class"] == "A_only", arm].astype(float)
    b = df.loc[df["class"] == "B_only", arm].astype(float)
    return {
        "auroc_D_vs_A": auroc(d, a),
        "auroc_D_vs_B": auroc(d, b),
        "summary_min": float(np.nanmin([auroc(d, a), auroc(d, b)])),
    }


def write_sdfs():
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_48.csv")
    logs = ROOT / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in TARGETS:
        out_sdf = logs / f"{target}_poses.sdf"
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        for lig in panel["panel_id"]:
            tmpl = Chem.RemoveHs(
                Chem.SDMolSupplier(
                    str(ROOT / "ligands_sdf" / f"{lig}.sdf"), removeHs=False
                )[0]
            )
            for mp in sorted((ROOT / "poses" / target / lig).glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                xyz = pdbqt_xyz(mp)
                pairs = smiles_idx_pairs(mp)
                if not pairs or len(pairs) != tmpl.GetNumAtoms():
                    raise RuntimeError(f"{mp} pairs mismatch")
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
        print("wrote", out_sdf, n)


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
        print("OK", csv_path)


def parse_id(s):
    m = re.search(r"(PM48_\d+)_mode(\d+)", str(s))
    return m.group(1), int(m.group(2))


def build_and_compare():
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_48.csv")
    ids = panel["panel_id"].tolist()
    vina_rows = []
    for lig in ids:
        for t in TARGETS:
            for mp in sorted((ROOT / "poses" / t / lig).glob("mode_*.pdbqt")):
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
            rtm_rows.append(
                {
                    "ligand": lig,
                    "target": t,
                    "vina_mode": mode,
                    "rtmscore": float(r[sc_col]),
                }
            )
    rtm_long = pd.DataFrame(rtm_rows)
    rtm_long.to_csv(ROOT / "tables" / "scores_rtm_all_poses.csv", index=False)
    best = (
        rtm_long.sort_values(
            ["ligand", "target", "rtmscore"], ascending=[True, True, False]
        )
        .groupby(["ligand", "target"], as_index=False)
        .first()
    )
    a = best[best.target == "4L23"][["ligand", "rtmscore"]].rename(
        columns={"rtmscore": "rtm_4L23"}
    )
    b = best[best.target == "4JT6"][["ligand", "rtmscore"]].rename(
        columns={"rtmscore": "rtm_4JT6"}
    )
    df = panel.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "pref_name", "molecule_chembl_id"]
    ]
    df = df.merge(m1.reset_index(), on="ligand").merge(a, on="ligand").merge(b, on="ligand")
    df = df.rename(columns={"4L23": "4L23_affinity", "4JT6": "4JT6_affinity"})
    df["vina_4L23_hb"] = -df["4L23_affinity"]
    df["vina_4JT6_hb"] = -df["4JT6_affinity"]
    df["vina_mean"] = (df["vina_4L23_hb"] + df["vina_4JT6_hb"]) / 2
    df["vina_min"] = df[["vina_4L23_hb", "vina_4JT6_hb"]].min(axis=1)
    df["rtm_mean"] = (df["rtm_4L23"] + df["rtm_4JT6"]) / 2
    df["rtm_min"] = df[["rtm_4L23", "rtm_4JT6"]].min(axis=1)
    for col, zcol in [("rtm_4L23", "rtm_4L23_z"), ("rtm_4JT6", "rtm_4JT6_z")]:
        mu, sd = df[col].mean(), df[col].std(ddof=0)
        df[zcol] = (df[col] - mu) / (sd if sd > 0 else 1.0)
    df["rtm_min_z"] = df[["rtm_4L23_z", "rtm_4JT6_z"]].min(axis=1)
    df["prep"] = "rdkit_meeko"
    df.to_csv(ROOT / "tables" / "ablation_ligand_scores.csv", index=False)

    old = pd.read_csv(OLD / "ablation_ligand_scores.csv").rename(
        columns={
            "vina_mean": "vina_mean_ligprep",
            "rtm_min_z": "rtm_min_z_ligprep",
            "rtm_mean": "rtm_mean_ligprep",
        }
    )
    m = df.merge(
        old[["ligand", "vina_mean_ligprep", "rtm_min_z_ligprep", "rtm_mean_ligprep"]],
        on="ligand",
        how="left",
    )
    m["delta_vina_mean"] = m["vina_mean"] - m["vina_mean_ligprep"]
    m["delta_rtm_min_z"] = m["rtm_min_z"] - m["rtm_min_z_ligprep"]
    m.to_csv(ROOT / "tables" / "prep_delta_vs_ligprep.csv", index=False)

    rows = []
    for label, vcol, rcol in (
        ("ligprep_old", "vina_mean_ligprep", "rtm_min_z_ligprep"),
        ("rdkit_new", "vina_mean", "rtm_min_z"),
    ):
        tmp = m.copy()
        for arm, col in (("vina_mean", vcol), ("rtm_min_z", rcol)):
            dmet = directional(tmp.assign(**{arm: tmp[col]}), arm)
            rows.append({"prep": label, "arm": arm, **dmet})
    dir_cmp = pd.DataFrame(rows)
    dir_cmp.to_csv(ROOT / "tables" / "directional_by_prep.csv", index=False)

    REPO.mkdir(parents=True, exist_ok=True)
    (REPO / "tables").mkdir(exist_ok=True)
    (REPO / "analysis").mkdir(exist_ok=True)
    for name in (
        "ablation_ligand_scores.csv",
        "prep_delta_vs_ligprep.csv",
        "directional_by_prep.csv",
        "job_status.csv",
        "scores_vina_long.csv",
        "scores_rtm_all_poses.csv",
        "panel_v0_48.csv",
    ):
        src = ROOT / "tables" / name
        if src.exists():
            (REPO / "tables" / name).write_bytes(src.read_bytes())

    md = [
        "# PM48 RDKit reprep — prep sensitivity\n\n",
        f"- Protocol: 4L23/4JT6, E=16, seed={SEED}, n_modes=9, RTM best-of-9\n",
        "- Prep: RDKit ETKDG + meeko\n\n",
        "## Directional by prep\n\n```\n",
        dir_cmp.to_string(index=False) + "\n```\n\n",
        f"- |Δvina_mean| median={m['delta_vina_mean'].abs().median():.3f}\n",
        f"- |Δrtm_min_z| median={m['delta_rtm_min_z'].abs().median():.3f}\n",
    ]
    text = "".join(md)
    (ROOT / "analysis" / "PREP_DELTA.md").write_text(text)
    (REPO / "analysis" / "PREP_DELTA.md").write_text(text)
    print(dir_cmp.to_string(index=False))


def main():
    missing = []
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_48.csv")
    for lig in panel["panel_id"]:
        for t in TARGETS:
            if not (ROOT / "poses" / t / lig / "mode_01.pdbqt").exists():
                missing.append((t, lig))
    if missing:
        raise SystemExit(f"missing poses n={len(missing)} e.g. {missing[:5]}")
    write_sdfs()
    run_rtm()
    build_and_compare()


if __name__ == "__main__":
    main()
