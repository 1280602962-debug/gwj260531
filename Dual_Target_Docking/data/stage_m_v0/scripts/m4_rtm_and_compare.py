#!/usr/bin/env python3
"""M4-min postprocess: RTM on re-prepped EH40 + LigPrep vs RDKit delta table."""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(
    "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_reprep_rdkit_v0"
)
OLD = Path(
    "/home/gwj/repos/gwj260531/Dual_Target_Docking/data/egfr_her2_panel120_v0/tables"
)
PANEL = OLD / "panel_v0_120.csv"
STAGE_M = Path(
    "/home/gwj/repos/gwj260531/Dual_Target_Docking/data/stage_m_v0"
)
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
TARGETS = ["3POZ", "3RCD"]
SEED = 20260727


def pdbqt_xyz(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
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
        "auroc_pooled": auroc(d, list(a) + list(b)),
        "summary_min": float(
            np.nanmin([auroc(d, a), auroc(d, b)])
        ),
    }


def write_pose_sdfs():
    panel = pd.read_csv(PANEL)
    ids = panel.loc[panel.from_panel40 == "yes", "panel_id"].tolist()
    logs = ROOT / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in TARGETS:
        out_sdf = logs / f"{target}_poses.sdf"
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        for lig in ids:
            tmpl = Chem.RemoveHs(
                Chem.SDMolSupplier(str(ROOT / "ligands_sdf" / f"{lig}.sdf"), removeHs=False)[0]
            )
            for mode_path in sorted((ROOT / "poses" / target / lig).glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mode_path.name).group(1))
                xyz = pdbqt_xyz(mode_path)
                pairs = smiles_idx_pairs(mode_path)
                if not pairs or len(pairs) != tmpl.GetNumAtoms():
                    raise RuntimeError(
                        f"{mode_path}: pairs={0 if not pairs else len(pairs)} tmpl={tmpl.GetNumAtoms()}"
                    )
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


def parse_rtm_id(s: str):
    m = re.search(r"(EH40_\d+)_mode(\d+)", str(s))
    if not m:
        raise ValueError(s)
    return m.group(1), int(m.group(2))


def build_scores():
    panel = pd.read_csv(PANEL)
    old40 = panel[panel.from_panel40 == "yes"].copy()
    ids = old40["panel_id"].tolist()

    # vina from new poses
    vina_rows = []
    for lig in ids:
        for target in TARGETS:
            for mp in sorted((ROOT / "poses" / target / lig).glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                vina_rows.append(
                    {
                        "ligand": lig,
                        "target": target,
                        "vina_mode": mode,
                        "vina_score": vina_from_pdbqt(mp),
                    }
                )
    vina_long = pd.DataFrame(vina_rows)
    vina_long.to_csv(ROOT / "tables" / "scores_vina_long.csv", index=False)
    m1 = vina_long.loc[vina_long["vina_mode"] == 1].pivot(
        index="ligand", columns="target", values="vina_score"
    )
    # rtm
    rtm_rows = []
    for target in TARGETS:
        d = pd.read_csv(ROOT / "logs" / "rtmscore" / f"{target}_rtmscore.csv")
        id_col = "id" if "id" in d.columns else d.columns[0]
        sc_col = "score" if "score" in d.columns else (
            "rtmscore" if "rtmscore" in d.columns else d.columns[1]
        )
        for _, r in d.iterrows():
            lig, mode = parse_rtm_id(r[id_col])
            rtm_rows.append(
                {
                    "ligand": lig,
                    "target": target,
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
    a = best[best.target == "3POZ"][["ligand", "rtmscore", "vina_mode"]].rename(
        columns={"rtmscore": "rtm_3POZ", "vina_mode": "rtm_mode_3POZ"}
    )
    b = best[best.target == "3RCD"][["ligand", "rtmscore", "vina_mode"]].rename(
        columns={"rtmscore": "rtm_3RCD", "vina_mode": "rtm_mode_3RCD"}
    )
    df = old40.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "pref_name", "molecule_chembl_id"]
    ]
    df = df.merge(m1.reset_index(), on="ligand", how="left")
    df = df.rename(columns={"3POZ": "3POZ_affinity", "3RCD": "3RCD_affinity"})
    df = df.merge(a, on="ligand").merge(b, on="ligand")
    df["vina_3POZ_hb"] = -df["3POZ_affinity"]
    df["vina_3RCD_hb"] = -df["3RCD_affinity"]
    df["vina_mean"] = (df["vina_3POZ_hb"] + df["vina_3RCD_hb"]) / 2
    df["vina_min"] = df[["vina_3POZ_hb", "vina_3RCD_hb"]].min(axis=1)
    df["rtm_mean"] = (df["rtm_3POZ"] + df["rtm_3RCD"]) / 2
    df["rtm_min"] = df[["rtm_3POZ", "rtm_3RCD"]].min(axis=1)
    for col, zcol in [("rtm_3POZ", "rtm_3POZ_z"), ("rtm_3RCD", "rtm_3RCD_z")]:
        mu, sd = df[col].mean(), df[col].std(ddof=0)
        df[zcol] = (df[col] - mu) / (sd if sd > 0 else 1.0)
    df["rtm_min_z"] = df[["rtm_3POZ_z", "rtm_3RCD_z"]].min(axis=1)
    df["prep"] = "rdkit_meeko"
    df.to_csv(ROOT / "tables" / "ablation_ligand_scores.csv", index=False)
    return df


def compare_to_ligprep(new_df: pd.DataFrame):
    old = pd.read_csv(OLD / "ablation_ligand_scores.csv")
    old = old[old.from_panel40.astype(str) == "yes"].copy()
    old = old.rename(
        columns={
            "vina_mean": "vina_mean_ligprep",
            "rtm_min_z": "rtm_min_z_ligprep",
            "rtm_mean": "rtm_mean_ligprep",
            "vina_min": "vina_min_ligprep",
        }
    )
    m = new_df.merge(
        old[
            [
                "ligand",
                "class",
                "vina_mean_ligprep",
                "vina_min_ligprep",
                "rtm_mean_ligprep",
                "rtm_min_z_ligprep",
            ]
        ],
        on=["ligand", "class"],
        how="left",
    )
    m["delta_vina_mean"] = m["vina_mean"] - m["vina_mean_ligprep"]
    m["delta_rtm_min_z"] = m["rtm_min_z"] - m["rtm_min_z_ligprep"]
    m["abs_delta_vina_mean"] = m["delta_vina_mean"].abs()
    m["abs_delta_rtm_min_z"] = m["delta_rtm_min_z"].abs()

    rows = []
    for prep_label, frame, vcol, rcol in (
        ("ligprep_old", m, "vina_mean_ligprep", "rtm_min_z_ligprep"),
        ("rdkit_new", m, "vina_mean", "rtm_min_z"),
    ):
        tmp = m.copy()
        tmp["_v"] = tmp[vcol]
        tmp["_r"] = tmp[rcol]
        for arm, col in (("vina_mean", "_v"), ("rtm_min_z", "_r")):
            dmet = directional(tmp.assign(**{arm: tmp[col]}), arm)
            rows.append({"prep": prep_label, "arm": arm, **dmet})

    dir_cmp = pd.DataFrame(rows)
    for c in ("auroc_D_vs_A", "auroc_D_vs_B", "auroc_pooled", "summary_min"):
        dir_cmp[c] = dir_cmp[c].round(4)

    out_delta = STAGE_M / "tables" / "m4_old40_prep_delta.csv"
    out_delta.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(out_delta, index=False)
    dir_cmp.to_csv(STAGE_M / "tables" / "m4_directional_by_prep.csv", index=False)
    m.to_csv(ROOT / "tables" / "m4_old40_prep_delta.csv", index=False)
    dir_cmp.to_csv(ROOT / "tables" / "m4_directional_by_prep.csv", index=False)

    # conclusion: does new70-style RTM collapse persist on unified old40?
    # Compare RDKit old40 rtm vs vina directional summaries
    rd = dir_cmp[dir_cmp.prep == "rdkit_new"]
    lg = dir_cmp[dir_cmp.prep == "ligprep_old"]
    rtm_rd = rd[rd.arm == "rtm_min_z"].iloc[0]
    vina_rd = rd[rd.arm == "vina_mean"].iloc[0]
    rtm_lg = lg[lg.arm == "rtm_min_z"].iloc[0]
    vina_lg = lg[lg.arm == "vina_mean"].iloc[0]

    # "RTM worse than vina on summary_min" under each prep
    rtm_worse_rd = rtm_rd.summary_min < vina_rd.summary_min - 0.02
    rtm_worse_lg = rtm_lg.summary_min < vina_lg.summary_min - 0.02
    # also check if LigPrep advantage of RTM disappears after RDKit
    ligprep_rtm_advantage = rtm_lg.summary_min > vina_lg.summary_min + 0.02
    advantage_gone = ligprep_rtm_advantage and (
        rtm_rd.summary_min <= vina_rd.summary_min + 0.02
    )

    md = []
    md.append("# M4 — Unified prep (M4-min: EH40 RDKit re-dock)\n\n")
    md.append("## STATUS: **Go** (M4-min completed)\n\n")
    md.append(
        "- Protocol: 3POZ/3RCD, E=8, seed=20260727, n_modes=9, RTM best-of-9\n"
        "- Prep: **RDKit ETKDG + meeko** (same as panel120 new70)\n"
        "- Scope: EH40 only (from_panel40=yes); M4-full not run\n\n"
    )
    md.append("## Directional AUROC by prep (same 40 ligands)\n\n")
    try:
        md.append(dir_cmp.to_markdown(index=False) + "\n\n")
    except ImportError:
        md.append("```\n" + dir_cmp.to_string(index=False) + "\n```\n\n")
    md.append("## Per-ligand score deltas (RDKit − LigPrep)\n\n")
    md.append(
        f"- |Δvina_mean| median={m.abs_delta_vina_mean.median():.3f}; "
        f"p90={m.abs_delta_vina_mean.quantile(0.9):.3f}\n"
        f"- |Δrtm_min_z| median={m.abs_delta_rtm_min_z.median():.3f}; "
        f"p90={m.abs_delta_rtm_min_z.quantile(0.9):.3f}\n\n"
    )
    md.append("## Interpretation\n\n")
    if advantage_gone:
        md.append(
            "- LigPrep-era RTM > vina advantage **disappears** under unified RDKit → "
            "old40/new70 RTM split is largely a **prep confound**, not a confirmed method claim.\n"
        )
    elif rtm_worse_rd and not rtm_worse_lg:
        md.append(
            "- Under RDKit, RTM is worse than vina on summary_min while LigPrep showed RTM≥vina → prep-sensitive.\n"
        )
    else:
        md.append(
            f"- RDKit: vina summary_min={vina_rd.summary_min:.3f}, "
            f"rtm_min_z={rtm_rd.summary_min:.3f}; "
            f"LigPrep: vina={vina_lg.summary_min:.3f}, rtm={rtm_lg.summary_min:.3f}.\n"
            f"- rtm_worse_rd={rtm_worse_rd}, rtm_worse_lg={rtm_worse_lg}, "
            f"advantage_gone={advantage_gone}.\n"
        )
    md.append(
        "\nDo **not** treat mixed-prep panel120 RTM conclusions as method-validated.\n"
    )
    (STAGE_M / "analysis" / "M4_UNIFIED_PREP.md").write_text("".join(md))
    (ROOT / "analysis" / "M4_UNIFIED_PREP.md").write_text("".join(md))
    print(dir_cmp.to_string(index=False))
    print("advantage_gone", advantage_gone, "wrote", STAGE_M / "analysis" / "M4_UNIFIED_PREP.md")
    return dir_cmp


def main():
    missing = []
    panel = pd.read_csv(PANEL)
    for lig in panel.loc[panel.from_panel40 == "yes", "panel_id"]:
        for t in TARGETS:
            if not (ROOT / "poses" / t / lig / "mode_01.pdbqt").exists():
                missing.append((t, lig))
    if missing:
        raise SystemExit(f"missing poses n={len(missing)} e.g. {missing[:5]}")
    write_pose_sdfs()
    run_rtm()
    df = build_scores()
    compare_to_ligprep(df)


if __name__ == "__main__":
    main()
