#!/usr/bin/env python3
"""Stage-1 postprocess: RTM (new ligands) + merge EH40 scores + bootstrap S1 verdict."""
from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0")
SRC40 = Path("/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0")
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
TARGETS = ["3POZ", "3RCD"]
N_BOOT = 2000
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
        # fallback: sequential heavy atoms (RDKit meeko usually writes IDX)
        return None
    return list(zip(nums[0::2], nums[1::2]))


def vina_from_pdbqt(path: Path):
    for line in path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def write_pose_sdfs_for_new():
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_120.csv")
    new_ids = panel.loc[panel.from_panel40 == "no", "panel_id"].tolist()
    logs = ROOT / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in TARGETS:
        out_sdf = logs / f"{target}_new_poses.sdf"
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        for lig in new_ids:
            tmpl_path = ROOT / "ligands_sdf" / f"{lig}.sdf"
            tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(tmpl_path), removeHs=False)[0])
            for mode_path in sorted((ROOT / "poses" / target / lig).glob("mode_*.pdbqt")):
                m = re.search(r"mode_(\d+)", mode_path.name)
                mode = int(m.group(1))
                xyz = pdbqt_xyz(mode_path)
                pairs = smiles_idx_pairs(mode_path)
                mol = Chem.Mol(tmpl)
                conf = Chem.Conformer(mol.GetNumAtoms())
                if not pairs or len(pairs) != mol.GetNumAtoms():
                    raise RuntimeError(
                        f"{mode_path}: pairs={0 if not pairs else len(pairs)} "
                        f"tmpl={mol.GetNumAtoms()}"
                    )
                for s_idx, p_idx in pairs:
                    conf.SetAtomPosition(s_idx - 1, xyz[p_idx - 1])
                mol.RemoveAllConformers()
                mol.AddConformer(conf, assignId=True)
                mol.SetProp("_Name", f"{lig}_mode{mode}")
                w.write(mol)
                n += 1
        w.close()
        print("wrote", out_sdf, "n=", n)


def run_rtm_new():
    logs = ROOT / "logs" / "rtmscore"
    for target in TARGETS:
        sdf = logs / f"{target}_new_poses.sdf"
        pocket = ROOT / "receptors" / f"{target}_pocket_10.0.pdb"
        out_prefix = logs / f"{target}_rtmscore_new"
        log = logs / f"{target}_rtmscore_new.log"
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
    m = re.search(r"(EH(?:40|120)_\d+)_mode(\d+)", str(s))
    if not m:
        raise ValueError(s)
    return m.group(1), int(m.group(2))


def collect_vina_long():
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_120.csv")
    rows = []
    for _, r in panel.iterrows():
        lig = r["panel_id"]
        for target in TARGETS:
            for mode_path in sorted((ROOT / "poses" / target / lig).glob("mode_*.pdbqt")):
                mode = int(re.search(r"mode_(\d+)", mode_path.name).group(1))
                rows.append(
                    {
                        "ligand": lig,
                        "target": target,
                        "vina_mode": mode,
                        "vina_score": vina_from_pdbqt(mode_path),
                    }
                )
    vdf = pd.DataFrame(rows)
    vdf.to_csv(ROOT / "tables" / "scores_vina_long.csv", index=False)
    # mode1 summary — use .loc; avoid column name "mode" (pandas DataFrame.mode)
    m1 = vdf.loc[vdf["vina_mode"] == 1].pivot(
        index="ligand", columns="target", values="vina_score"
    )
    out = panel[["panel_id", "class", "pref_name", "molecule_chembl_id"]].rename(
        columns={"panel_id": "ligand"}
    )
    out = out.merge(m1.reset_index(), on="ligand", how="left")
    out = out.rename(columns={"3POZ": "3POZ_affinity", "3RCD": "3RCD_affinity"})
    out["vina_mean_raw"] = (out["3POZ_affinity"] + out["3RCD_affinity"]) / 2
    out.to_csv(ROOT / "tables" / "scores_vina.csv", index=False)
    return out


def collect_rtm_all():
    old = pd.read_csv(SRC40 / "tables" / "scores_rtm_all_poses.csv")
    # panel40 format: id,rtmscore,target,ligand,vina_mode,vina_score
    old = old[["ligand", "target", "vina_mode", "rtmscore"]].copy()
    old = old[old.ligand.str.startswith("EH40_")]

    new_rows = []
    for target in TARGETS:
        d = pd.read_csv(ROOT / "logs" / "rtmscore" / f"{target}_rtmscore_new.csv")
        id_col = "id" if "id" in d.columns else d.columns[0]
        sc_col = "score" if "score" in d.columns else (
            "rtmscore" if "rtmscore" in d.columns else d.columns[1]
        )
        for _, r in d.iterrows():
            lig, mode = parse_rtm_id(r[id_col])
            new_rows.append(
                {
                    "ligand": lig,
                    "target": target,
                    "vina_mode": mode,
                    "rtmscore": float(r[sc_col]),
                }
            )
    new = pd.DataFrame(new_rows)
    all_rtm = pd.concat([old, new], ignore_index=True)
    all_rtm.to_csv(ROOT / "tables" / "scores_rtm_all_poses.csv", index=False)
    return all_rtm


def build_ablation(vina_sum: pd.DataFrame, rtm_long: pd.DataFrame):
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_120.csv")
    best = (
        rtm_long.sort_values(
            ["ligand", "target", "rtmscore"], ascending=[True, True, False]
        )
        .groupby(["ligand", "target"], as_index=False)
        .first()
        .rename(columns={"vina_mode": "best_rtm_mode", "rtmscore": "rtmscore"})
    )
    a = best[best.target == "3POZ"][["ligand", "rtmscore", "best_rtm_mode"]].rename(
        columns={"rtmscore": "rtm_3POZ", "best_rtm_mode": "rtm_mode_3POZ"}
    )
    b = best[best.target == "3RCD"][["ligand", "rtmscore", "best_rtm_mode"]].rename(
        columns={"rtmscore": "rtm_3RCD", "best_rtm_mode": "rtm_mode_3RCD"}
    )
    df = panel.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "pref_name", "molecule_chembl_id", "from_panel40"]
    ]
    df = df.merge(
        vina_sum[["ligand", "3POZ_affinity", "3RCD_affinity"]], on="ligand", how="left"
    )
    df = df.merge(a, on="ligand", how="left").merge(b, on="ligand", how="left")
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
    df.to_csv(ROOT / "tables" / "ablation_ligand_scores.csv", index=False)
    return df


def roc_auc(y, s):
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    c = 0.0
    for p in pos:
        c += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(c / (len(pos) * len(neg)))


def bootstrap_delta(df, n_boot=N_BOOT, seed=SEED):
    y0 = (df["class"] == "dual").astype(int).values
    sa = df["vina_mean"].values.astype(float)
    sb = df["rtm_min_z"].values.astype(float)
    n = len(df)
    rng = np.random.default_rng(seed)
    d = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        y = y0[idx]
        if y.sum() == 0 or y.sum() == n:
            continue
        d.append(roc_auc(y, sb[idx]) - roc_auc(y, sa[idx]))
    d = np.asarray(d)
    point = roc_auc(y0, sb) - roc_auc(y0, sa)
    lo, hi = np.percentile(d, 2.5), np.percentile(d, 97.5)
    return {
        "n": n,
        "n_dual": int(y0.sum()),
        "auroc_vina_mean": roc_auc(y0, sa),
        "auroc_rtm_min_z": roc_auc(y0, sb),
        "delta_auroc": point,
        "delta_auroc_ci_lo": float(lo),
        "delta_auroc_ci_hi": float(hi),
        "significant_excl0": not (lo <= 0 <= hi),
        "n_boot": len(d),
    }


def top10_stats(df, arm):
    order = np.argsort(-df[arm].values)
    top = df.iloc[order[:10]]
    return {
        "arm": arm,
        "top10_dual": int((top["class"] == "dual").sum()),
        "top10_A_only": int((top["class"] == "A_only").sum()),
        "top10_B_only": int((top["class"] == "B_only").sum()),
        "top10_hardneg": int((top["class"] != "dual").sum()),
    }


def write_verdict(df, boot):
    out = ROOT / "analysis" / "STAGE1_VERDICT.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    go = bool(boot["significant_excl0"])
    t_vina = top10_stats(df, "vina_mean")
    t_rtm = top10_stats(df, "rtm_min_z")
    lines = []
    lines.append("# STAGE1 VERDICT — EGFR/HER2 panel120 expand (S1 gate)\n\n")
    lines.append(f"- N = **{boot['n']}** (dual={boot['n_dual']})\n")
    lines.append(f"- Protocol: 3POZ/3RCD, E=8, seed={SEED}, n_modes=9, RTM best-of-9\n")
    lines.append(
        f"- AUROC vina_mean = **{boot['auroc_vina_mean']:.3f}**; "
        f"rtm_min_z = **{boot['auroc_rtm_min_z']:.3f}**\n"
    )
    lines.append(
        f"- ΔAUROC(rtm_min_z − vina_mean) = **{boot['delta_auroc']:+.3f}** "
        f"95% CI [{boot['delta_auroc_ci_lo']:+.3f}, {boot['delta_auroc_ci_hi']:+.3f}] "
        f"(B={boot['n_boot']})\n"
    )
    lines.append(
        f"- Top10 hardneg: vina={t_vina['top10_hardneg']} "
        f"(A={t_vina['top10_A_only']}, B={t_vina['top10_B_only']}); "
        f"rtm={t_rtm['top10_hardneg']} "
        f"(A={t_rtm['top10_A_only']}, B={t_rtm['top10_B_only']})\n\n"
    )
    if go:
        lines.append("## Verdict: **Go (S1)**\n\n")
        lines.append("ΔAUROC CI excludes 0. Proceed to Stage-2 multi-pair planning.\n")
    else:
        lines.append("## Verdict: **No-Go (S1)**\n\n")
        lines.append(
            "ΔAUROC CI includes 0 at N≈110. Do **not** stack complex methods yet. "
            "Recommended downgrade: keep DualFourClass-Bench + failure typology as the "
            "publishable claim; optional further EGFR expansion only if effect estimate remains "
            "large; otherwise accept diagnosis-paper route (JCIM gap scenario 甲).\n"
        )
    lines.append(
        "\n### Notes\n"
        "- Panel40 poses reused; new ligands RDKit+meeko (documented).\n"
        "- Architecture not used as selection filter; flags not gated into score.\n"
        "- No clash/shortfall retune for this gate.\n"
    )
    out.write_text("".join(lines))
    pd.DataFrame([boot]).to_csv(ROOT / "tables" / "stage1_bootstrap_delta.csv", index=False)
    pd.DataFrame([t_vina, t_rtm]).to_csv(ROOT / "tables" / "stage1_top10.csv", index=False)
    print("".join(lines))


def main():
    # ensure all new poses present
    panel = pd.read_csv(ROOT / "tables" / "panel_v0_120.csv")
    new_ids = panel.loc[panel.from_panel40 == "no", "panel_id"].tolist()
    missing = []
    for lig in new_ids:
        for t in TARGETS:
            if not (ROOT / "poses" / t / lig / "mode_01.pdbqt").exists():
                missing.append((t, lig))
    if missing:
        raise SystemExit(f"missing poses: {missing[:10]} ... n={len(missing)}")

    write_pose_sdfs_for_new()
    run_rtm_new()
    vina_sum = collect_vina_long()
    # inspect old rtm format
    old = pd.read_csv(SRC40 / "tables" / "scores_rtm_all_poses.csv")
    print("old rtm cols", list(old.columns)[:10], "n", len(old))
    rtm_long = collect_rtm_all()
    df = build_ablation(vina_sum, rtm_long)
    print("ablation n", len(df), "null rtm", df["rtm_min_z"].isna().sum())
    boot = bootstrap_delta(df)
    write_verdict(df, boot)


if __name__ == "__main__":
    main()
