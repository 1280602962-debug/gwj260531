#!/usr/bin/env python3
"""
Join protocol scores onto True/Random benchmarks; compute EF@1%, EF@5%, AUC;
select Π* by pre-specified rule.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.metrics import roc_auc_score


PROTOCOLS = [
    ("P1", "P1_vina_affinity", False),
    ("P2", "P2_CNNaffinity", True),
    ("P3", "P3_gnina_affinity", False),
    ("P0", "P0_CNNscore", True),
    ("P4", "P4_RTMScore", True),
    ("P5", "P5_RTMScore", True),
]

DRUGS = {
    "lesinurad": "O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12",
    "benzbromarone": "CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1",
    "verinurad": "CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O",
    "dotinurad": "O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21",
}


def canon(smi: str) -> str | None:
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m, canonical=True) if m else None


def enrichment_factor(y_true: np.ndarray, scores: np.ndarray, frac: float, higher_better: bool) -> float:
    """EF at top frac of library."""
    n = len(y_true)
    n_act = int(y_true.sum())
    if n_act == 0 or n == 0:
        return float("nan")
    k = max(1, int(round(n * frac)))
    order = np.argsort(-scores if higher_better else scores)
    top = y_true[order[:k]]
    hit_rate = top.mean()
    base = n_act / n
    return float(hit_rate / base) if base > 0 else float("nan")


def auc_safe(y_true: np.ndarray, scores: np.ndarray, higher_better: bool) -> float:
    if len(np.unique(y_true)) < 2:
        return float("nan")
    s = scores if higher_better else -scores
    try:
        return float(roc_auc_score(y_true, s))
    except Exception:
        return float("nan")


def eval_bench(bench: pd.DataFrame, scores: pd.DataFrame, score_col: str, higher_better: bool) -> dict:
    df = bench.merge(scores[["canonical_smiles", score_col]], on="canonical_smiles", how="left")
    df = df.dropna(subset=[score_col, "label"])
    y = df["label"].astype(int).to_numpy()
    s = df[score_col].astype(float).to_numpy()
    return {
        "n": int(len(df)),
        "n_act": int(y.sum()),
        "EF1": enrichment_factor(y, s, 0.01, higher_better),
        "EF5": enrichment_factor(y, s, 0.05, higher_better),
        "AUC": auc_safe(y, s, higher_better),
    }


def drug_percentiles(scores: pd.DataFrame, score_col: str, higher_better: bool) -> dict:
    """Percentile rank in full scored pool (0-100, higher = better rank)."""
    d = scores.dropna(subset=[score_col]).copy()
    if d.empty:
        return {k: None for k in DRUGS}
    vals = d[score_col].astype(float).to_numpy()
    out = {}
    smi_map = {canon(v): k for k, v in DRUGS.items()}
    for smi, row in d.set_index("canonical_smiles")[score_col].items():
        c = canon(str(smi)) or str(smi)
        if c not in smi_map:
            # try exact
            if str(smi) not in {canon(x) or x for x in DRUGS.values()}:
                continue
        name = smi_map.get(c)
        if not name:
            for dn, ds in DRUGS.items():
                if c == (canon(ds) or ds):
                    name = dn
                    break
        if not name:
            continue
        v = float(row)
        if higher_better:
            pct = float((vals <= v).mean() * 100.0)
        else:
            pct = float((vals >= v).mean() * 100.0)
        out[name] = pct
    for k in DRUGS:
        out.setdefault(k, None)
    return out


def select_pi(metrics: pd.DataFrame, random_veto: bool = True) -> dict:
    """
    TrueDecoy EF1 -> EF5 -> AUC; RandomDecoy veto if EF1 << cohort median;
    tie-break: mean four-drug percentile.
    """
    m = metrics.copy()
    # Random veto: EF1 more than 0.5 below median of protocols with finite EF1
    if random_veto and "random_EF1" in m.columns:
        med = m["random_EF1"].median(skipna=True)
        m["random_veto"] = m["random_EF1"] < (med - 0.5)
    else:
        m["random_veto"] = False

    cand = m[~m["random_veto"]].copy()
    if cand.empty:
        cand = m.copy()

    cand = cand.sort_values(
        by=["true_EF1", "true_EF5", "true_AUC", "drug_pct_mean"],
        ascending=[False, False, False, False],
        na_position="last",
    )
    best = cand.iloc[0]
    return {
        "selected_protocol": best["protocol"],
        "score_col": best["score_col"],
        "true_EF1": float(best["true_EF1"]) if pd.notna(best["true_EF1"]) else None,
        "true_EF5": float(best["true_EF5"]) if pd.notna(best["true_EF5"]) else None,
        "true_AUC": float(best["true_AUC"]) if pd.notna(best["true_AUC"]) else None,
        "random_EF1": float(best["random_EF1"]) if pd.notna(best["random_EF1"]) else None,
        "drug_pct_mean": float(best["drug_pct_mean"]) if pd.notna(best["drug_pct_mean"]) else None,
        "random_veto_applied": bool(best["random_veto"]),
        "rule": "True EF@1% > EF@5% > AUC; Random EF@1% veto vs median-0.5; tie -> mean four-drug percentile",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--true-bench", required=True)
    ap.add_argument("--random-bench", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--no-random-veto", action="store_true")
    args = ap.parse_args()

    work = Path(args.work)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    score_path = work / "scores" / "mol_protocol_scores.csv"
    if not score_path.exists():
        raise SystemExit(
            f"Missing {score_path}. Run: python collect_dock_scores.py --work ... --pool ... "
            "or bash 05_join_and_select_pi.sh"
        )

    scores = pd.read_csv(score_path)
    # canonicalize join keys
    scores["canonical_smiles"] = scores["canonical_smiles"].map(lambda s: canon(str(s)) or s)
    true_b = pd.read_csv(args.true_bench)
    rand_b = pd.read_csv(args.random_bench)
    true_b["canonical_smiles"] = true_b["canonical_smiles"].map(lambda s: canon(str(s)) or s)
    rand_b["canonical_smiles"] = rand_b["canonical_smiles"].map(lambda s: canon(str(s)) or s)

    rows = []
    for pid, col, higher in PROTOCOLS:
        if col not in scores.columns:
            continue
        t = eval_bench(true_b, scores, col, higher)
        r = eval_bench(rand_b, scores, col, higher)
        dp = drug_percentiles(scores, col, higher)
        dp_vals = [v for v in dp.values() if v is not None]
        rows.append(
            {
                "protocol": pid,
                "score_col": col,
                "higher_better": higher,
                "true_n": t["n"],
                "true_EF1": t["EF1"],
                "true_EF5": t["EF5"],
                "true_AUC": t["AUC"],
                "random_n": r["n"],
                "random_EF1": r["EF1"],
                "random_EF5": r["EF5"],
                "random_AUC": r["AUC"],
                "drug_pct_lesinurad": dp.get("lesinurad"),
                "drug_pct_benzbromarone": dp.get("benzbromarone"),
                "drug_pct_verinurad": dp.get("verinurad"),
                "drug_pct_dotinurad": dp.get("dotinurad"),
                "drug_pct_mean": float(np.mean(dp_vals)) if dp_vals else np.nan,
            }
        )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(outdir / "protocol_metrics.csv", index=False)

    selection = select_pi(metrics, random_veto=not args.no_random_veto)
    (outdir / "selected_pi.json").write_text(json.dumps(selection, indent=2))
    print(metrics.to_string(index=False))
    print("Selected:", json.dumps(selection, indent=2))


if __name__ == "__main__":
    main()
