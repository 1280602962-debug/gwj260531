#!/usr/bin/env python3
"""Shared helpers for Stage M measurement audit (M1–M3, M5)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")

STAGE_M = Path(__file__).resolve().parents[1]
DUAL_ROOT = Path(__file__).resolve().parents[3]
TABLES = STAGE_M / "tables"
ANALYSIS = STAGE_M / "analysis"
SEED = 20260728

PAIRS = {
    "EGFR_HER2": {
        "panel": DUAL_ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        "scores": DUAL_ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        "id_col": "panel_id",
        "ends": ("pchembl_EGFR", "pchembl_HER2"),
        "score_arms": ["vina_mean", "vina_min", "rtm_mean", "rtm_min", "rtm_min_z"],
        "subset_col": "from_panel40",
    },
    "PIK3CA_mTOR": {
        "panel": DUAL_ROOT / "data/pik3ca_mtor_panel48_v0/tables/panel_v0_48.csv",
        "scores": DUAL_ROOT / "data/pik3ca_mtor_panel48_v0/tables/ablation_ligand_scores.csv",
        "id_col": "panel_id",
        "ends": ("pchembl_PIK3CA", "pchembl_MTOR"),
        "score_arms": ["vina_mean", "vina_min", "rtm_mean", "rtm_min", "rtm_min_z"],
        "subset_col": None,
    },
}

BASELINE_ARMS = ["heavy_atoms", "MW", "cLogP", "TPSA", "morgan_dual_medsim"]
DOCK_ARMS = ["vina_mean", "vina_min", "rtm_mean", "rtm_min", "rtm_min_z"]


def auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = 0.0
    for p in pos:
        wins += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(wins / (len(pos) * len(neg)))


def parse_float(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "na"}:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def descriptors_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MW": float(Descriptors.MolWt(mol)),
        "cLogP": float(Descriptors.MolLogP(mol)),
        "heavy_atoms": float(mol.GetNumHeavyAtoms()),
        "TPSA": float(Descriptors.TPSA(mol)),
        "_mol": mol,
    }


def morgan_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def tanimoto(fp_a, fp_b) -> float:
    return float(AllChem.DataStructs.TanimotoSimilarity(fp_a, fp_b))


def load_merged(pair: str) -> pd.DataFrame:
    cfg = PAIRS[pair]
    panel = pd.read_csv(cfg["panel"])
    scores = pd.read_csv(cfg["scores"])
    panel = panel.rename(columns={cfg["id_col"]: "ligand"})
    df = scores.merge(panel, on="ligand", how="left", suffixes=("", "_panel"))
    if "class" not in df.columns and "class_panel" in df.columns:
        df["class"] = df["class_panel"]
    # resolve class column if duplicated
    if "class_x" in df.columns:
        df["class"] = df["class_x"].fillna(df.get("class_y"))

    rows = []
    dual_fps = {}
    for _, r in df.iterrows():
        smiles = r.get("smiles")
        desc = descriptors_from_smiles(smiles) if pd.notna(smiles) else None
        if desc is None:
            continue
        item = r.to_dict()
        item.update({k: v for k, v in desc.items() if k != "_mol"})
        item["_mol"] = desc["_mol"]
        item["pA"] = parse_float(r.get(cfg["ends"][0]))
        item["pB"] = parse_float(r.get(cfg["ends"][1]))
        rows.append(item)
        if item.get("class") == "dual":
            dual_fps[item["ligand"]] = morgan_fp(desc["_mol"])

    # leave-one-out median Tanimoto to other duals (or all duals if singleton)
    for item in rows:
        fp = morgan_fp(item["_mol"])
        others = [dual_fps[k] for k in dual_fps if k != item["ligand"]]
        if not others:
            others = list(dual_fps.values())
        if others:
            sims = [tanimoto(fp, o) for o in others]
            item["morgan_dual_medsim"] = float(np.median(sims))
        else:
            item["morgan_dual_medsim"] = float("nan")
        del item["_mol"]

    out = pd.DataFrame(rows)
    out["pair"] = pair
    return out


def directional_metrics(df: pd.DataFrame, arm: str) -> dict:
    d = df.loc[df["class"] == "dual", arm].astype(float).tolist()
    a = df.loc[df["class"] == "A_only", arm].astype(float).tolist()
    b = df.loc[df["class"] == "B_only", arm].astype(float).tolist()
    # top10 over dual+A+B (exclude neither from ranking pool for hardneg report)
    pool = df[df["class"].isin(["dual", "A_only", "B_only"])].copy()
    pool = pool.sort_values(arm, ascending=False)
    top = pool.head(10)
    return {
        "arm": arm,
        "n_dual": len(d),
        "n_A_only": len(a),
        "n_B_only": len(b),
        "auroc_D_vs_A": auroc(d, a),
        "auroc_D_vs_B": auroc(d, b),
        "auroc_pooled": auroc(d, a + b),
        "top10_A_only": int((top["class"] == "A_only").sum()),
        "top10_B_only": int((top["class"] == "B_only").sum()),
        "top10_dual": int((top["class"] == "dual").sum()),
        "summary_min": float(np.nanmin([auroc(d, a), auroc(d, b)])),
        "summary_mean": float(np.nanmean([auroc(d, a), auroc(d, b)])),
    }


def assign_fourclass(pA, pB, cutoff: float, measured_only=True):
    """四类规则：未测 ≠ 阴。两端都测过才分 dual/A_only/B_only/neither；否则 None。"""
    if pA is None or pB is None:
        return None  # incomplete → exclude from four-class (not negative)
    a_pos = pA >= cutoff
    b_pos = pB >= cutoff
    if a_pos and b_pos:
        return "dual"
    if a_pos and not b_pos:
        return "A_only"
    if b_pos and not a_pos:
        return "B_only"
    return "neither"


def assign_margin_label(pA, pB):
    """Strict margin labels; gray if both measured but not strict."""
    if pA is None or pB is None:
        return "incomplete"
    if pA >= 6.5 and pB >= 6.5:
        return "dual_strict"
    if pA >= 6.5 and pB <= 5.5:
        return "A_only_strict"
    if pB >= 6.5 and pA <= 5.5:
        return "B_only_strict"
    if pA <= 5.5 and pB <= 5.5:
        return "neither_strict"
    return "gray"


def map_strict_to_class(label: str):
    return {
        "dual_strict": "dual",
        "A_only_strict": "A_only",
        "B_only_strict": "B_only",
        "neither_strict": "neither",
    }.get(label)
