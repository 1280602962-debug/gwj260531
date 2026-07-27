#!/usr/bin/env python3
"""Rebuild ablation tables from freeze pack (egfr_her2_panel40_v0).

Run from anywhere:
  python results/egfr_her2_panel40_v0/scripts/build_ablation_table.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

OUT = Path(__file__).resolve().parents[1]
CLASH_CUTOFF = 2.2
CLASH_FAIL_N = 3
ARMS = ["vina_mean", "vina_min", "rtm_mean", "rtm_min", "rtm_min_z", "gated_rtm_min"]
KEY_LIGS = ["EH40_01", "EH40_02", "EH40_05", "EH40_18", "EH40_23", "EH40_28", "EH40_33"]


def roc_auc_score(y_true, y_score):
    try:
        from sklearn.metrics import roc_auc_score as sk_auc
        return float(sk_auc(y_true, y_score))
    except Exception:
        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)
        pos = y_score[y_true == 1]
        neg = y_score[y_true == 0]
        if len(pos) == 0 or len(neg) == 0:
            return float("nan")
        correct = 0.0
        for p in pos:
            correct += np.sum(p > neg) + 0.5 * np.sum(p == neg)
        return correct / (len(pos) * len(neg))


def parse_pdbqt_heavy(path: Path):
    coords = []
    if not path.exists():
        return np.zeros((0, 3))
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line.split()[-1] in ("H", "HD", "HS"):
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords) if coords else np.zeros((0, 3))


def parse_protein_heavy(pdb_path: Path):
    coords = []
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        elem = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        if elem.startswith("H"):
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords)


def clash_count(lig_xyz, prot_xyz, cutoff=CLASH_CUTOFF):
    if len(lig_xyz) == 0 or len(prot_xyz) == 0:
        return 0
    n = 0
    for i in range(0, len(lig_xyz), 20):
        L = lig_xyz[i : i + 20]
        d = np.linalg.norm(L[:, None, :] - prot_xyz[None, :, :], axis=2)
        n += int((d < cutoff).sum())
    return n


def main():
    panel = pd.read_csv(OUT / "tables" / "panel_v0_40.csv")
    vina = pd.read_csv(OUT / "tables" / "scores_vina.csv")
    rtm_poz = pd.read_csv(OUT / "tables" / "scores_rtm_3POZ.csv")
    rtm_rcd = pd.read_csv(OUT / "tables" / "scores_rtm_3RCD.csv")

    df = panel.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "pref_name", "molecule_chembl_id"]
    ].copy()
    df = df.merge(
        vina[["ligand", "3POZ_affinity", "3RCD_affinity", "mean_affinity"]],
        on="ligand",
        how="left",
    )
    df = df.merge(
        rtm_poz[["ligand", "rtmscore", "best_rtm_mode"]].rename(
            columns={"rtmscore": "rtm_3POZ", "best_rtm_mode": "rtm_mode_3POZ"}
        ),
        on="ligand",
        how="left",
    )
    df = df.merge(
        rtm_rcd[["ligand", "rtmscore", "best_rtm_mode"]].rename(
            columns={"rtmscore": "rtm_3RCD", "best_rtm_mode": "rtm_mode_3RCD"}
        ),
        on="ligand",
        how="left",
    )

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

    prot = {
        "3POZ": parse_protein_heavy(OUT / "receptors" / "3POZ_protein.pdb"),
        "3RCD": parse_protein_heavy(OUT / "receptors" / "3RCD_protein.pdb"),
    }
    rows = []
    for _, r in df.iterrows():
        lig = r["ligand"]
        rec = {"ligand": lig}
        for tag, mcol in [("3POZ", "rtm_mode_3POZ"), ("3RCD", "rtm_mode_3RCD")]:
            mode = int(r[mcol])
            pose = OUT / "poses" / tag / lig / f"mode_{mode:02d}.pdbqt"
            n = clash_count(parse_pdbqt_heavy(pose), prot[tag])
            rec[f"clash_{tag}_rtm_best"] = n
            rec[f"pb_like_fail_{tag}"] = int(n >= CLASH_FAIL_N)
        rows.append(rec)
    df = df.merge(pd.DataFrame(rows), on="ligand", how="left")
    p5_poz = float(np.percentile(df["rtm_3POZ"], 5))
    p5_rcd = float(np.percentile(df["rtm_3RCD"], 5))
    df["rtm_3POZ_gated"] = np.where(df["pb_like_fail_3POZ"] == 1, p5_poz, df["rtm_3POZ"])
    df["rtm_3RCD_gated"] = np.where(df["pb_like_fail_3RCD"] == 1, p5_rcd, df["rtm_3RCD"])
    df["gated_rtm_min"] = df[["rtm_3POZ_gated", "rtm_3RCD_gated"]].min(axis=1)

    eval_mask = df["class"].isin(["dual", "A_only", "B_only"])
    eval_df = df[eval_mask]
    y = (eval_df["class"] == "dual").astype(int).values

    metrics, ranks = [], []
    for arm in ARMS:
        auroc = roc_auc_score(y, eval_df[arm].values)
        order = df.sort_values(arm, ascending=False).reset_index(drop=True)
        order["rank"] = np.arange(1, len(order) + 1)
        top10 = order.head(10)
        n_a = int((top10["class"] == "A_only").sum())
        n_b = int((top10["class"] == "B_only").sum())
        metrics.append(
            {
                "arm": arm,
                "AUROC_dual_vs_AorB": round(float(auroc), 4),
                "top10_dual_count": int((top10["class"] == "dual").sum()),
                "top10_A_only_count": n_a,
                "top10_B_only_count": n_b,
                "top10_hardneg_count": n_a + n_b,
                "top10_hardneg_fraction": round((n_a + n_b) / 10.0, 3),
                "top10_neither_count": int((top10["class"] == "neither").sum()),
                "n_eval_AUROC": int(eval_mask.sum()),
                "n_dual_eval": int((eval_df["class"] == "dual").sum()),
                "n_hardneg_eval": int(eval_df["class"].isin(["A_only", "B_only"]).sum()),
            }
        )
        for lig in KEY_LIGS:
            rr = order[order["ligand"] == lig].iloc[0]
            ranks.append(
                {
                    "arm": arm,
                    "ligand": lig,
                    "rank": int(rr["rank"]),
                    "score": float(rr[arm]),
                    "class": rr["class"],
                    "pref_name": rr.get("pref_name", ""),
                }
            )

    pd.DataFrame(metrics).to_csv(OUT / "tables" / "ablation_metrics.csv", index=False)
    ranks_df = pd.DataFrame(ranks)
    ranks_df.to_csv(OUT / "tables" / "ablation_ranks.csv", index=False)
    ranks_df.pivot(index="ligand", columns="arm", values="rank").to_csv(
        OUT / "tables" / "ablation_ranks_wide.csv"
    )
    df.to_csv(OUT / "tables" / "ablation_ligand_scores.csv", index=False)
    print("Wrote ablation tables under", OUT / "tables")


if __name__ == "__main__":
    main()
