#!/usr/bin/env python3
"""P1: chemotype warning flags + frozen-threshold decision ablation.

Does NOT modify docking scores for ranking via flags.
Does NOT retune clash. Uses existing ablation_ligand_scores.csv only.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS

REPO = Path(__file__).resolve().parents[5]  # Dual_Target_Docking/
PM = REPO / "data" / "pik3ca_mtor_panel48_v0"
EH = REPO / "data" / "egfr_her2_panel40_v0"
OUT_PM = PM / "analysis" / "decision_ablation_v0"
OUT_FLAGS_PM = PM / "tables"
OUT_FLAGS_EH = EH / "tables"
OUT_PM.mkdir(parents=True, exist_ok=True)

# ---- FROZEN before inspecting case-level wins/losses ----
FROZEN = {
    "shortfall_lambda": 0.5,  # from typology pretest grid; not retuned here
    "consensus_top_frac": 0.25,  # AND gate: both arms in top quartile
    "mcs_atoms_high": 12,  # vs cognate heavy scaffold (~half of PI-103/TAK-285 core)
    "tanimoto_cognate_warn": 0.20,
}
KEY_LIGS_PM = ["PM48_01", "PM48_02", "PM48_10", "PM48_20", "PM48_21", "PM48_26", "PM48_34"]
KEY_LIGS_EH = ["EH40_01", "EH40_02", "EH40_18", "EH40_23"]

SMARTS = {
    "morpholine": "C1COCCN1",
    "amino_triazine_like": "c1ncnc(N)n1",  # loose; matches typology scaffold_flags
    "anilinoquinazoline": "c1ccc(Nc2ncnc3ccccc23)cc1",
    "warhead_acrylamide": "[C;H2,H1]=[C;H1]C(=O)N",
    "warhead_alkynamide": "C#CC(=O)N",
}


def roc_auc(y_true, y_score):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score, dtype=float)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    correct = 0.0
    for p in pos:
        correct += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return correct / (len(pos) * len(neg))


def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def mcs_atoms(mol_a, mol_b) -> int:
    res = rdFMCS.FindMCS(
        [mol_a, mol_b],
        timeout=8,
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
    )
    if res.canceled or not res.smartsString:
        return 0
    m = Chem.MolFromSmarts(res.smartsString)
    return 0 if m is None else m.GetNumAtoms()


def has_smarts(mol, smarts: str) -> bool:
    q = Chem.MolFromSmarts(smarts)
    return bool(q and mol.HasSubstructMatch(q))


def build_warning_flags_panel48():
    panel = pd.read_csv(PM / "tables" / "panel_v0_48.csv")
    cognate = Chem.MolFromSmiles(panel.loc[panel.panel_id == "PM48_01", "smiles"].iloc[0])
    cognate_fp = fp(cognate)
    rows = []
    for _, r in panel.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            raise RuntimeError(r["panel_id"])
        tan = DataStructs.TanimotoSimilarity(fp(mol), cognate_fp)
        mcs = mcs_atoms(mol, cognate)
        morph = has_smarts(mol, SMARTS["morpholine"])
        triaz = has_smarts(mol, SMARTS["amino_triazine_like"])
        high_mcs = mcs >= FROZEN["mcs_atoms_high"]
        high_tan = tan >= FROZEN["tanimoto_cognate_warn"]
        morph_atp = bool(morph and (high_mcs or high_tan))
        flags = []
        if triaz:
            flags.append("amino_triazine_like")
        if morph:
            flags.append("morpholine")
        if morph_atp:
            flags.append("morpholine_ATP_like")
        if high_mcs:
            flags.append("high_MCS_to_cognate")
        if high_tan and r["panel_id"] != "PM48_01":
            flags.append("tanimoto_to_cognate_ge_0.20")
        # typology anchors (diagnostic labels, still not score mods)
        if r["panel_id"] in ("PM48_26", "PM48_20", "PM48_21"):
            flags.append("T2_anchor")
        if r["panel_id"] in ("PM48_10", "PM48_02"):
            flags.append("T5_anchor")
        if r["panel_id"] == "PM48_34":
            flags.append("T1_anchor")
            flags.append("vina_modes_lt_9")
        rows.append(
            {
                "panel_id": r["panel_id"],
                "pair": "PIK3CA_mTOR",
                "class": r["class"],
                "pref_name": r.get("pref_name", ""),
                "amino_triazine_like": int(triaz),
                "morpholine": int(morph),
                "morpholine_ATP_like": int(morph_atp),
                "high_MCS_to_cognate": int(high_mcs),
                "mcs_atoms_to_cognate": mcs,
                "tanimoto_to_cognate": round(float(tan), 4),
                "anilinoquinazoline": 0,
                "warhead_covalent": 0,
                "warning_flags": "|".join(flags) if flags else "",
                "flags_enter_score": 0,
                "note": "diagnostic only; do not gate score in v0",
            }
        )
    df = pd.DataFrame(rows)
    path = OUT_FLAGS_PM / "warning_flags.csv"
    df.to_csv(path, index=False)
    print("wrote", path, "flagged", int((df.warning_flags != "").sum()))
    return df


def build_warning_flags_panel40():
    panel = pd.read_csv(EH / "tables" / "panel_v0_40.csv")
    cognate = Chem.MolFromSmiles(panel.loc[panel.panel_id == "EH40_01", "smiles"].iloc[0])
    cognate_fp = fp(cognate)
    rows = []
    for _, r in panel.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            raise RuntimeError(r["panel_id"])
        tan = DataStructs.TanimotoSimilarity(fp(mol), cognate_fp)
        mcs = mcs_atoms(mol, cognate)
        anil = has_smarts(mol, SMARTS["anilinoquinazoline"])
        war = has_smarts(mol, SMARTS["warhead_acrylamide"]) or has_smarts(
            mol, SMARTS["warhead_alkynamide"]
        )
        morph = has_smarts(mol, SMARTS["morpholine"])
        triaz = has_smarts(mol, SMARTS["amino_triazine_like"])
        high_mcs = mcs >= FROZEN["mcs_atoms_high"]
        flags = []
        if anil:
            flags.append("anilinoquinazoline")
        if war:
            flags.append("warhead_covalent")
        if triaz:
            flags.append("amino_triazine_like")
        if morph:
            flags.append("morpholine")
        if high_mcs:
            flags.append("high_MCS_to_cognate")
        if r["panel_id"] == "EH40_23":
            flags.append("T2_anchor")
        if r["panel_id"] == "EH40_18":
            flags.append("T1_anchor")
        rows.append(
            {
                "panel_id": r["panel_id"],
                "pair": "EGFR_HER2",
                "class": r["class"],
                "pref_name": r.get("pref_name", ""),
                "amino_triazine_like": int(triaz),
                "morpholine": int(morph),
                "morpholine_ATP_like": 0,
                "high_MCS_to_cognate": int(high_mcs),
                "mcs_atoms_to_cognate": mcs,
                "tanimoto_to_cognate": round(float(tan), 4),
                "anilinoquinazoline": int(anil),
                "warhead_covalent": int(war),
                "warning_flags": "|".join(flags) if flags else "",
                "flags_enter_score": 0,
                "note": "diagnostic only; do not gate score in v0",
            }
        )
    df = pd.DataFrame(rows)
    path = OUT_FLAGS_EH / "warning_flags.csv"
    df.to_csv(path, index=False)
    print("wrote", path, "flagged", int((df.warning_flags != "").sum()))
    return df


def metrics_block(df, score_col, y):
    score = df[score_col].values.astype(float)
    order = np.argsort(-score)
    top10 = df.iloc[order[:10]]
    return {
        "arm": score_col,
        "auroc_dual_vs_rest": roc_auc(y, score),
        "top10_dual": int((top10["class"] == "dual").sum()),
        "top10_A_only": int((top10["class"] == "A_only").sum()),
        "top10_B_only": int((top10["class"] == "B_only").sum()),
        "top10_neither": int((top10["class"] == "neither").sum()),
        "top10_hardneg": int((top10["class"] != "dual").sum()),
        "top1": df.iloc[order[0]]["ligand"],
    }


def df_to_md(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def decision_ablation_panel48():
    """Frozen arms A–D on existing scores only."""
    df = pd.read_csv(PM / "tables" / "ablation_ligand_scores.csv")
    n = len(df)
    top_n = max(1, int(np.ceil(FROZEN["consensus_top_frac"] * n)))
    lam = FROZEN["shortfall_lambda"]

    # A, B already present
    # C: shortfall
    df["rtm_shortfall"] = df["rtm_min_z"] - lam * (
        df["rtm_4L23_z"] - df["rtm_4JT6_z"]
    ).abs()

    # ranks (1=best) for consensus
    df["rank_vina_mean"] = df["vina_mean"].rank(ascending=False, method="min").astype(int)
    df["rank_rtm_min_z"] = df["rtm_min_z"].rank(ascending=False, method="min").astype(int)

    # D1 primary consensus: mean of ranks → higher-better score
    df["consensus_rank_mean"] = -(df["rank_vina_mean"] + df["rank_rtm_min_z"]) / 2.0

    # D2 AND top-quartile gate (frozen frac); among passers use rtm_min_z, else deep penalty
    passes = (df["rank_vina_mean"] <= top_n) & (df["rank_rtm_min_z"] <= top_n)
    df["consensus_and_top25"] = np.where(passes, df["rtm_min_z"], df["rtm_min_z"] - 1e3)
    df["consensus_and_top25_pass"] = passes.astype(int)

    y = (df["class"] == "dual").astype(int).values
    arms = [
        "vina_mean",
        "rtm_min_z",
        "rtm_shortfall",
        "consensus_rank_mean",
        "consensus_and_top25",
    ]
    metrics = [metrics_block(df, a, y) for a in arms]
    met = pd.DataFrame(metrics)
    met["frozen_shortfall_lambda"] = lam
    met["frozen_consensus_top_frac"] = FROZEN["consensus_top_frac"]
    met["frozen_top_n"] = top_n
    met.to_csv(OUT_PM / "decision_ablation_metrics.csv", index=False)

    # ranks table
    rank_cols = {a: df[a].rank(ascending=False, method="min").astype(int) for a in arms}
    ranks = df[["ligand", "class", "pref_name"]].copy()
    for a, s in rank_cols.items():
        ranks[a] = s
    ranks.to_csv(OUT_PM / "decision_ablation_ranks.csv", index=False)

    scores_out = df[
        [
            "ligand",
            "class",
            "pref_name",
            "vina_mean",
            "rtm_min_z",
            "rtm_shortfall",
            "consensus_rank_mean",
            "consensus_and_top25",
            "consensus_and_top25_pass",
            "rank_vina_mean",
            "rank_rtm_min_z",
        ]
    ].copy()
    scores_out.to_csv(OUT_PM / "decision_ablation_scores.csv", index=False)

    key = ranks[ranks.ligand.isin(KEY_LIGS_PM)].copy()
    key.to_csv(OUT_PM / "decision_ablation_key_ranks.csv", index=False)

    # success vs rtm_min_z baseline
    base = met[met.arm == "rtm_min_z"].iloc[0]
    lines = []
    lines.append("# Decision ablation v0 — frozen thresholds\n")
    lines.append("## Frozen a priori (before case inspection)\n")
    lines.append(f"- `shortfall_lambda` = **{lam}**\n")
    lines.append(
        f"- `consensus_top_frac` = **{FROZEN['consensus_top_frac']}** → top_n = **{top_n}** / {n}\n"
    )
    lines.append(
        "- Flags/chemotype warnings **do not** enter gated scores (`flags_enter_score=0`).\n"
    )
    lines.append("- No clash retune; no new docking.\n")
    lines.append("\n## Arms\n")
    lines.append("| arm | definition |\n|-----|------------|\n")
    lines.append("| A vina_mean | (−aff_A − aff_B)/2 |\n")
    lines.append("| B rtm_min_z | min of per-target RTM z |\n")
    lines.append("| C rtm_shortfall | rtm_min_z − λ\\|zA−zB\\| |\n")
    lines.append("| D1 consensus_rank_mean | −mean(rank_vina, rank_rtm) |\n")
    lines.append(
        "| D2 consensus_and_top25 | pass if both ranks ≤ top_n; else deep penalty |\n"
    )
    lines.append("\n## Metrics (Dual vs rest)\n\n")
    lines.append(df_to_md(met))
    lines.append("\n\n## Key ligand ranks (1=best)\n\n")
    lines.append(df_to_md(key))
    lines.append("\n\n## Verdict vs success criteria\n")
    lines.append(
        "Success = hardneg Top10 ↓ vs B **and** Torin1/Omipalisib (PM48_10/02) not clearly worse than under B.\n\n"
    )

    def hardneg(mrow):
        return int(mrow["top10_A_only"]) + int(mrow["top10_B_only"]) + int(mrow["top10_neither"])

    base_h = hardneg(base)
    t10 = int(key.loc[key.ligand == "PM48_10", "rtm_min_z"].iloc[0])
    t02 = int(key.loc[key.ligand == "PM48_02", "rtm_min_z"].iloc[0])

    verdict_bits = []
    for arm in ["rtm_shortfall", "consensus_rank_mean", "consensus_and_top25"]:
        row = met[met.arm == arm].iloc[0]
        h = hardneg(row)
        r10 = int(key.loc[key.ligand == "PM48_10", arm].iloc[0])
        r02 = int(key.loc[key.ligand == "PM48_02", arm].iloc[0])
        hardneg_improved = h < base_h
        # "not clearly worse": rank number not increased by >5 places (lower rank# = better)
        injured_ok = (r10 <= t10 + 5) and (r02 <= t02 + 5)
        ok = hardneg_improved and injured_ok
        verdict_bits.append((arm, h, base_h, r10, t10, r02, t02, ok))
        lines.append(
            f"- **{arm}**: hardneg Top10 {h} vs baseline {base_h}; "
            f"PM48_10 rank {r10} (B={t10}); PM48_02 rank {r02} (B={t02}) → "
            f"{'PASS' if ok else 'FAIL'}\n"
        )

    any_pass = any(v[-1] for v in verdict_bits)
    lines.append("\n### Bottom line\n\n")
    if any_pass:
        lines.append(
            "At least one consensus/shortfall arm meets the joint criterion under frozen thresholds.\n"
        )
    else:
        lines.append(
            "**无法同时满足**：在冻结阈值下，shortfall / consensus 未能在降低硬负 Top10 的同时"
            "保护 Torin1/Omipalisib；或硬负 Top10 未下降。"
            "主文应并列报告 `vina_mean` 与 `rtm_min_z`，并用化学型警告层标注 T2，"
            "而不是宣称决策规则已闭环。\n"
        )
    lines.append(
        "\nDo **not** claim C4 extrapolation success. Do **not** retune clash to drop PM48_26.\n"
    )

    (OUT_PM / "DECISION_ABLATION_V0.md").write_text("".join(lines))
    # also dump frozen params
    with (OUT_PM / "frozen_thresholds.yaml").open("w") as fh:
        for k, v in FROZEN.items():
            fh.write(f"{k}: {v}\n")
        fh.write(f"panel_n: {n}\n")
        fh.write(f"consensus_top_n: {top_n}\n")
    print(met.to_string(index=False))
    print("key ranks:\n", key.to_string(index=False))
    print("wrote", OUT_PM)


def main():
    build_warning_flags_panel48()
    build_warning_flags_panel40()
    decision_ablation_panel48()


if __name__ == "__main__":
    main()
