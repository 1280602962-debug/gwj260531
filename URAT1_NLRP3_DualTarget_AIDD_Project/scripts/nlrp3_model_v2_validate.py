#!/usr/bin/env python3
"""NLRP3 classifier v2: nested scaffold CV with deployment-consistent evaluation.

Rebuilds the NLRP3 functional-activity classifier used to shrink the clinical
library, fixing two defects of the historical pipeline (scripts/02_train_asymmetric_models.py):

  1. ``fit_nlrp3_final`` fitted the base XGBoost on every record and then fitted
     the isotonic calibrator on a subset of those same records (in-sample
     calibration); the StandardScaler was also fitted on all records.
  2. Scaffold CV scored each test record with its *true* assay ChEMBL ID, while
     the clinical library is scored as max over five fixed assay contexts, so the
     reported AUROC did not evaluate the deployed prediction rule.

Here every outer-test molecule is scored exactly the way an unseen clinical
compound would be, calibration is fitted on inner scaffolds the base model never
saw, and a structure-only model is compared against the assay-conditioned model
on identical outer folds.

Nothing in this script touches the frozen production artefacts; all output goes
to results/nlrp3_model_v2/.
"""
from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    brier_score_loss,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit

from nlrp3_v2_lib import (
    N_DEPLOY_ASSAYS,
    AssayConditionedModel,
    StructureOnlyModel,
    apply_calibrator,
    fit_calibrators,
    xgb_classifier,
)
from utils_ml import enrichment_factor, scaffold_cv_indices

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
OUT_DIR = PROJECT_ROOT / "results" / "nlrp3_model_v2"

CAL_FRACTION = 0.25
PRESPECIFIED_THRESHOLDS = (0.3, 0.4, 0.5, 0.6, 0.7)


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------


def load_records() -> pd.DataFrame:
    df = pd.read_csv(PROCESSED / "nlrp3_records.csv", low_memory=False)
    df["canonical_smiles"] = df["canonical_smiles"].astype(str)
    df["Assay ChEMBL ID"] = df["Assay ChEMBL ID"].astype(str)
    return df


def molecule_table(df: pd.DataFrame) -> pd.DataFrame:
    """Molecule-level labels. Historical definition kept: active = max over records."""
    mol = (
        df.groupby("canonical_smiles", as_index=False)
        .agg(
            active=("active", "max"),
            scaffold=("scaffold", "first"),
            n_records=("active", "size"),
            max_pactivity=("pActivity", "max"),
            min_pactivity=("pActivity", "min"),
        )
        .reset_index(drop=True)
    )
    return mol


def audit_data(df: pd.DataFrame, mol: pd.DataFrame) -> dict:
    per_assay = []
    for aid, g in df.groupby("Assay ChEMBL ID"):
        per_assay.append(
            {
                "assay_chembl_id": aid,
                "n_records": int(len(g)),
                "n_compounds": int(g["canonical_smiles"].nunique()),
                "n_active": int(g["active"].sum()),
                "n_inactive": int((1 - g["active"]).sum()),
                "active_rate": float(g["active"].mean()),
                "cell_type": str(g["assay_cell_type"].iloc[0]),
            }
        )
    per_assay.sort(key=lambda r: -r["n_records"])

    multi = df.groupby("canonical_smiles")["active"].agg(["min", "max", "size"])
    discordant = multi[(multi["size"] > 1) & (multi["min"] != multi["max"])]

    scaff_sizes = mol.groupby("scaffold").size()
    return {
        "n_records": int(len(df)),
        "n_unique_compounds": int(df["canonical_smiles"].nunique()),
        "n_unique_assays": int(df["Assay ChEMBL ID"].nunique()),
        "n_unique_scaffolds": int(mol["scaffold"].nunique()),
        "endpoint": "pChEMBL Value (IL-1beta release assays, Assay Type B, relation '=')",
        "label_rule": "record active = pActivity >= 6.0; molecule active = max over records (historical rule, unchanged)",
        "record_active_rate": float(df["active"].mean()),
        "molecule_active_rate": float(mol["active"].mean()),
        "n_molecules_active": int(mol["active"].sum()),
        "n_molecules_inactive": int((1 - mol["active"]).sum()),
        "records_per_molecule": {
            "mean": float(df.groupby("canonical_smiles").size().mean()),
            "max": int(df.groupby("canonical_smiles").size().max()),
            "n_molecules_with_multiple_records": int((multi["size"] > 1).sum()),
        },
        "n_molecules_with_discordant_records": int(len(discordant)),
        "discordant_molecules": [
            {
                "canonical_smiles": smi,
                "n_records": int(r["size"]),
                "labels": df.loc[df["canonical_smiles"] == smi, "active"].tolist(),
                "pactivity": df.loc[df["canonical_smiles"] == smi, "pActivity"].round(2).tolist(),
            }
            for smi, r in discordant.iterrows()
        ],
        "scaffold_sizes": {
            "n_singleton_scaffolds": int((scaff_sizes == 1).sum()),
            "largest_scaffold_n_molecules": int(scaff_sizes.max()),
            "scaffold_definition": "Bemis-Murcko (RDKit MurckoScaffoldSmiles, includeChirality=False)",
        },
        "per_assay": per_assay,
        "standardisation": "RDKit MolFromSmiles -> SanitizeMol -> MolToSmiles; duplicates collapsed per canonical SMILES",
        "features": "Morgan radius-2 2048 bits + 12 RDKit descriptors (+ top-25 assay one-hot for the assay-conditioned model)",
    }


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------


def specificity_score(y: np.ndarray, pred: np.ndarray) -> float:
    neg = y == 0
    if neg.sum() == 0:
        return float("nan")
    return float((pred[neg] == 0).mean())


def classification_report(y: np.ndarray, p: np.ndarray, threshold: float = 0.5) -> dict:
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    pred = (p >= threshold).astype(int)
    two_class = len(np.unique(y)) == 2
    return {
        "n": int(len(y)),
        "n_active": int(y.sum()),
        "active_rate": float(y.mean()),
        "auroc": float(roc_auc_score(y, p)) if two_class else float("nan"),
        "auprc": float(average_precision_score(y, p)) if two_class else float("nan"),
        "accuracy": float(accuracy_score(y, pred)),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "recall": float(recall_score(y, pred, zero_division=0)),
        "specificity": specificity_score(y, pred),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "ef_10pct": float(enrichment_factor(y.astype(float), p, fraction=0.1)),
        "brier": float(brier_score_loss(y, np.clip(p, 0, 1))),
        "threshold": float(threshold),
    }


def reliability_curve(y: np.ndarray, p: np.ndarray, n_bins: int = 10) -> list[dict]:
    edges = np.linspace(0, 1, n_bins + 1)
    out = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        m = (p >= lo) & (p < hi) if i < n_bins - 1 else (p >= lo) & (p <= hi)
        if m.sum() == 0:
            out.append({"bin_lo": float(lo), "bin_hi": float(hi), "n": 0,
                        "mean_pred": None, "observed_rate": None})
        else:
            out.append(
                {
                    "bin_lo": float(lo),
                    "bin_hi": float(hi),
                    "n": int(m.sum()),
                    "mean_pred": float(p[m].mean()),
                    "observed_rate": float(np.asarray(y)[m].mean()),
                }
            )
    return out


def paired_bootstrap_delta_auroc(y, p_a, p_b, n_boot=2000, seed=0) -> dict:
    """Paired bootstrap over molecules for AUROC(B) - AUROC(A)."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    deltas = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        deltas.append(roc_auc_score(y[idx], p_b[idx]) - roc_auc_score(y[idx], p_a[idx]))
    deltas = np.asarray(deltas)
    observed = roc_auc_score(y, p_b) - roc_auc_score(y, p_a)
    return {
        "delta_auroc_observed": float(observed),
        "ci95_low": float(np.percentile(deltas, 2.5)),
        "ci95_high": float(np.percentile(deltas, 97.5)),
        "p_two_sided": float(2 * min((deltas <= 0).mean(), (deltas >= 0).mean())),
        "n_bootstrap": int(len(deltas)),
    }


# --------------------------------------------------------------------------
# nested CV
# --------------------------------------------------------------------------


def inner_split(mol_train: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    gss = GroupShuffleSplit(n_splits=1, test_size=CAL_FRACTION, random_state=seed)
    fit_rel, cal_rel = next(
        gss.split(np.arange(len(mol_train)), groups=mol_train["scaffold"].values)
    )
    return mol_train.iloc[fit_rel].copy(), mol_train.iloc[cal_rel].copy()


def run_nested_cv(df: pd.DataFrame, mol: pd.DataFrame, n_splits: int, seed: int,
                  inner_repeats: int) -> dict:
    smiles_all = mol["canonical_smiles"].tolist()
    folds = scaffold_cv_indices(smiles_all, n_splits=n_splits)

    prediction_rows = []
    split_rows = []
    leakage_checks = []
    inner_sensitivity = []

    for fold_i, (tr_idx, te_idx) in enumerate(folds):
        mol_tr = mol.iloc[tr_idx].copy()
        mol_te = mol.iloc[te_idx].copy()
        te_smiles = mol_te["canonical_smiles"].tolist()
        y_te = mol_te["active"].values.astype(int)
        rec_te = df[df["canonical_smiles"].isin(set(te_smiles))].copy()

        for s in mol_tr["canonical_smiles"]:
            split_rows.append({"fold": fold_i, "canonical_smiles": s, "outer_role": "train"})
        for s in te_smiles:
            split_rows.append({"fold": fold_i, "canonical_smiles": s, "outer_role": "test"})

        for rep in range(inner_repeats):
            inner_seed = seed + 100 * rep + fold_i
            mol_fit, mol_cal = inner_split(mol_tr, inner_seed)
            fit_smiles = set(mol_fit["canonical_smiles"])
            cal_smiles = set(mol_cal["canonical_smiles"])
            rec_fit = df[df["canonical_smiles"].isin(fit_smiles)].copy()
            rec_cal = df[df["canonical_smiles"].isin(cal_smiles)].copy()

            leakage_checks.append(
                {
                    "fold": fold_i,
                    "inner_repeat": rep,
                    "n_fit_molecules": int(len(mol_fit)),
                    "n_cal_molecules": int(len(mol_cal)),
                    "n_test_molecules": int(len(mol_te)),
                    "n_cal_active": int(mol_cal["active"].sum()),
                    "n_cal_inactive": int((1 - mol_cal["active"]).sum()),
                    "smiles_overlap_fit_cal": len(fit_smiles & cal_smiles),
                    "smiles_overlap_fit_test": len(fit_smiles & set(te_smiles)),
                    "smiles_overlap_cal_test": len(cal_smiles & set(te_smiles)),
                    "scaffold_overlap_fit_cal": len(
                        set(mol_fit["scaffold"]) & set(mol_cal["scaffold"])
                    ),
                    "scaffold_overlap_train_test": len(
                        set(mol_tr["scaffold"]) & set(mol_te["scaffold"])
                    ),
                    "calibration_usable": bool(
                        len(mol_cal) >= 30 and mol_cal["active"].nunique() == 2
                    ),
                }
            )

            model_a = StructureOnlyModel().fit(mol_fit, seed=seed + fold_i)
            model_b = AssayConditionedModel().fit(rec_fit, seed=seed + fold_i)

            variants = {
                "A_structure_only": {
                    "cal_raw": model_a.predict_molecules(mol_cal["canonical_smiles"].tolist()),
                    "cal_y": mol_cal["active"].values.astype(int),
                    "test_raw": model_a.predict_molecules(te_smiles),
                },
                "B_assay_deployment": {
                    "cal_raw": model_b.predict_molecules(mol_cal["canonical_smiles"].tolist()),
                    "cal_y": mol_cal["active"].values.astype(int),
                    "test_raw": model_b.predict_molecules(te_smiles),
                },
                "B_assay_legacy_true_assay": {
                    "cal_raw": model_b.predict_molecules_legacy(
                        rec_cal, mol_cal["canonical_smiles"].tolist()
                    ),
                    "cal_y": mol_cal["active"].values.astype(int),
                    "test_raw": model_b.predict_molecules_legacy(rec_te, te_smiles),
                },
            }

            for vname, v in variants.items():
                cals = fit_calibrators(v["cal_raw"], v["cal_y"])
                for cname, cal in cals.items():
                    p_test = apply_calibrator(cal, v["test_raw"])
                    rep_metrics = classification_report(y_te, p_test)
                    inner_sensitivity.append(
                        {"fold": fold_i, "inner_repeat": rep, "variant": vname,
                         "calibration": cname, **rep_metrics}
                    )
                    if rep == 0:
                        for j, smi in enumerate(te_smiles):
                            prediction_rows.append(
                                {
                                    "fold": fold_i,
                                    "canonical_smiles": smi,
                                    "scaffold": mol_te["scaffold"].values[j],
                                    "y_true": int(y_te[j]),
                                    "variant": vname,
                                    "calibration": cname,
                                    "p_raw": float(v["test_raw"][j]),
                                    "p_final": float(p_test[j]),
                                }
                            )

    return {
        "predictions": pd.DataFrame(prediction_rows),
        "splits": pd.DataFrame(split_rows),
        "leakage_checks": leakage_checks,
        "inner_sensitivity": pd.DataFrame(inner_sensitivity),
        "n_folds": n_splits,
    }


def summarise(pred: pd.DataFrame) -> dict:
    out = {}
    for (variant, cal), g in pred.groupby(["variant", "calibration"]):
        key = f"{variant}|{cal}"
        folds = []
        for f, gf in g.groupby("fold"):
            m = classification_report(gf["y_true"].values, gf["p_final"].values)
            m["fold"] = int(f)
            folds.append(m)
        pooled = classification_report(g["y_true"].values, g["p_final"].values)
        aurocs = [f["auroc"] for f in folds]
        out[key] = {
            "variant": variant,
            "calibration": cal,
            "pooled_outer_test": pooled,
            "fold_metrics": folds,
            "mean_fold_auroc": float(np.mean(aurocs)),
            "sd_fold_auroc": float(np.std(aurocs, ddof=1)),
            "min_fold_auroc": float(np.min(aurocs)),
            "max_fold_auroc": float(np.max(aurocs)),
            "threshold_sweep": {
                str(t): classification_report(g["y_true"].values, g["p_final"].values, t)
                for t in PRESPECIFIED_THRESHOLDS
            },
            "reliability_curve": reliability_curve(g["y_true"].values, g["p_final"].values),
            "prob_distribution": {
                "min": float(g["p_final"].min()),
                "p05": float(g["p_final"].quantile(0.05)),
                "median": float(g["p_final"].median()),
                "p95": float(g["p_final"].quantile(0.95)),
                "max": float(g["p_final"].max()),
                "frac_below_0_02": float((g["p_final"] < 0.02).mean()),
                "frac_above_0_98": float((g["p_final"] > 0.98).mean()),
                "n_distinct_values": int(g["p_final"].round(6).nunique()),
            },
        }
    return out


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------


def make_figures(pred: pd.DataFrame, summary: dict, out_dir: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = []
    variants = ["A_structure_only", "B_assay_deployment", "B_assay_legacy_true_assay"]
    cals = ["none", "platt", "isotonic"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True, sharey=True)
    for ax, variant in zip(axes, variants):
        ax.plot([0, 1], [0, 1], "k--", lw=1, label="perfect")
        for cal in cals:
            key = f"{variant}|{cal}"
            if key not in summary:
                continue
            pts = [(b["mean_pred"], b["observed_rate"]) for b in summary[key]["reliability_curve"]
                   if b["n"] >= 5]
            if pts:
                ax.plot(*zip(*pts), "o-", label=f"{cal} (Brier={summary[key]['pooled_outer_test']['brier']:.3f})")
        ax.set_title(variant, fontsize=10)
        ax.set_xlabel("mean predicted probability")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("observed active fraction")
    fig.suptitle("Pooled outer-test reliability (nested scaffold CV)")
    fig.tight_layout()
    p = out_dir / "calibration_curves.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    files.append(p.name)

    fig, axes = plt.subplots(len(cals), len(variants), figsize=(15, 9), sharex=True)
    for i, cal in enumerate(cals):
        for j, variant in enumerate(variants):
            ax = axes[i, j]
            g = pred[(pred["variant"] == variant) & (pred["calibration"] == cal)]
            if len(g):
                ax.hist(g.loc[g["y_true"] == 0, "p_final"], bins=20, range=(0, 1),
                        alpha=0.6, label="inactive")
                ax.hist(g.loc[g["y_true"] == 1, "p_final"], bins=20, range=(0, 1),
                        alpha=0.6, label="active")
                ax.axvline(0.5, color="k", ls="--", lw=1)
            if i == 0:
                ax.set_title(variant, fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{cal}\ncount", fontsize=9)
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Outer-test predicted probability distributions")
    fig.tight_layout()
    p = out_dir / "probability_distributions.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    files.append(p.name)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    width = 0.25
    xs = np.arange(5)
    for k, variant in enumerate(variants):
        key = f"{variant}|none"
        if key not in summary:
            continue
        vals = [f["auroc"] for f in sorted(summary[key]["fold_metrics"], key=lambda d: d["fold"])]
        ax.bar(xs + (k - 1) * width, vals, width, label=variant)
    ax.axhline(0.8934, color="r", ls="--", lw=1, label="historical mean AUROC 0.8934")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"fold {i}" for i in xs])
    ax.set_ylabel("outer-test AUROC")
    ax.set_ylim(0.5, 1.0)
    ax.legend(fontsize=7)
    ax.set_title("Per-fold outer-test AUROC (uncalibrated scores)")
    fig.tight_layout()
    p = out_dir / "fold_auroc_comparison.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    files.append(p.name)
    return files


# --------------------------------------------------------------------------
# final models
# --------------------------------------------------------------------------


def build_final_models(df: pd.DataFrame, mol: pd.DataFrame, seed: int) -> dict:
    mol_fit, mol_cal = inner_split(mol, seed)
    rec_fit = df[df["canonical_smiles"].isin(set(mol_fit["canonical_smiles"]))].copy()

    model_a = StructureOnlyModel().fit(mol_fit, seed=seed)
    model_b = AssayConditionedModel().fit(rec_fit, seed=seed)

    y_cal = mol_cal["active"].values.astype(int)
    cal_a = fit_calibrators(model_a.predict_molecules(mol_cal["canonical_smiles"].tolist()), y_cal)
    cal_b = fit_calibrators(model_b.predict_molecules(mol_cal["canonical_smiles"].tolist()), y_cal)

    bundles = {
        "structure_only": {
            "model": model_a,
            "calibrators": cal_a,
            "prediction_rule": "single molecule-level prediction from Morgan2048+12 RDKit descriptors",
        },
        "assay_conditioned": {
            "model": model_b,
            "calibrators": cal_b,
            "prediction_rule": f"max over the {N_DEPLOY_ASSAYS} most frequent training assay contexts",
        },
    }
    meta = {
        "n_fit_molecules": int(len(mol_fit)),
        "n_calibration_molecules": int(len(mol_cal)),
        "n_calibration_active": int(y_cal.sum()),
        "n_fit_records_assay_model": int(len(rec_fit)),
        "scaffold_overlap_fit_cal": len(set(mol_fit["scaffold"]) & set(mol_cal["scaffold"])),
        "smiles_overlap_fit_cal": len(
            set(mol_fit["canonical_smiles"]) & set(mol_cal["canonical_smiles"])
        ),
        "deploy_assay_ids": model_b.deploy_assay_ids,
        "top_assay_ids": model_b.top_assay_ids,
    }
    membership = pd.concat(
        [
            mol_fit[["canonical_smiles", "scaffold", "active"]].assign(final_role="base_fit"),
            mol_cal[["canonical_smiles", "scaffold", "active"]].assign(final_role="calibration"),
        ]
    )
    return bundles, meta, membership


# --------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUT_DIR)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--inner-repeats", type=int, default=3)
    ap.add_argument("--n-bootstrap", type=int, default=2000)
    args = ap.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    df = load_records()
    mol = molecule_table(df)
    audit = audit_data(df, mol)
    with open(out / "data_label_audit.json", "w") as f:
        json.dump(audit, f, indent=2)
    print(f"records={audit['n_records']} molecules={audit['n_unique_compounds']} "
          f"assays={audit['n_unique_assays']} scaffolds={audit['n_unique_scaffolds']} "
          f"molecule active rate={audit['molecule_active_rate']:.3f}")

    print("Running nested scaffold CV ...")
    cv = run_nested_cv(df, mol, args.n_splits, args.seed, args.inner_repeats)
    pred = cv["predictions"]
    pred.to_csv(out / "outer_test_predictions.csv", index=False)
    cv["splits"].to_csv(out / "fold_assignments.csv", index=False)
    cv["inner_sensitivity"].to_csv(out / "inner_repeat_sensitivity.csv", index=False)

    summary = summarise(pred)

    comparisons = {}
    for cal in ["none", "platt", "isotonic"]:
        a = pred[(pred["variant"] == "A_structure_only") & (pred["calibration"] == cal)]
        b = pred[(pred["variant"] == "B_assay_deployment") & (pred["calibration"] == cal)]
        l = pred[(pred["variant"] == "B_assay_legacy_true_assay") & (pred["calibration"] == cal)]
        if not len(a) or not len(b):
            continue
        a = a.set_index("canonical_smiles").loc[a["canonical_smiles"]]
        key = a.index
        b = b.set_index("canonical_smiles").reindex(key)
        l = l.set_index("canonical_smiles").reindex(key) if len(l) else None
        comparisons[f"B_minus_A|{cal}"] = paired_bootstrap_delta_auroc(
            a["y_true"].values, a["p_final"].values, b["p_final"].values, args.n_bootstrap
        )
        if l is not None:
            comparisons[f"legacy_minus_deployment|{cal}"] = paired_bootstrap_delta_auroc(
                a["y_true"].values, b["p_final"].values, l["p_final"].values, args.n_bootstrap
            )

    inner = cv["inner_sensitivity"]
    inner_var = (
        inner.groupby(["variant", "calibration"])[["auroc", "brier", "f1"]]
        .agg(["mean", "std"])
        .round(4)
    )
    inner_var.columns = ["_".join(c) for c in inner_var.columns]
    inner_var.reset_index().to_csv(out / "inner_repeat_summary.csv", index=False)

    figures = make_figures(pred, summary, out)

    print("Fitting final models ...")
    bundles, final_meta, final_membership = build_final_models(df, mol, args.seed)
    for name, bundle in bundles.items():
        joblib.dump(bundle, out / f"nlrp3_v2_{name}.joblib")
    final_membership.to_csv(out / "final_model_split_membership.csv", index=False)

    historical = {
        "source": "data/models/training_report.json (nlrp3) and results/model_validation/ml_rigor_summary.json",
        "auroc_mean_of_folds": 0.8934017784193573,
        "auprc_mean_of_folds": 0.9143469018214805,
        "ef_10pct": 1.5706359287488578,
        "brier_score": 0.1284,
        "comparability_note": (
            "The historical AUROC/AUPRC are means over folds of a CV in which each test record carried its "
            "true assay ChEMBL ID, so they are comparable only with the B_assay_legacy_true_assay variant. "
            "The historical Brier score was computed from the production model whose isotonic calibrator was "
            "fitted on records the base model had already been trained on, so it has no honest counterpart here."
        ),
    }

    report = {
        "generated_by": "scripts/nlrp3_model_v2_validate.py",
        "design": {
            "outer": f"{args.n_splits}-fold GroupKFold on Bemis-Murcko scaffold, molecule level",
            "inner": f"GroupShuffleSplit by scaffold, {CAL_FRACTION:.0%} of outer-train molecules reserved for calibration",
            "inner_repeats": args.inner_repeats,
            "label": audit["label_rule"],
            "variants": {
                "A_structure_only": "Morgan2048 + 12 RDKit descriptors, one row per molecule",
                "B_assay_deployment": "assay-conditioned model scored as max over top-5 training assay contexts (matches clinical-library screening)",
                "B_assay_legacy_true_assay": "assay-conditioned model scored with the test molecule's true assay IDs (historical evaluation, not achievable at deployment)",
            },
            "calibration_methods": ["none", "platt", "isotonic"],
            "prespecified_thresholds": list(PRESPECIFIED_THRESHOLDS),
        },
        "data_audit": audit,
        "leakage_checks": cv["leakage_checks"],
        "results": summary,
        "model_comparisons_paired_bootstrap": comparisons,
        "final_models": final_meta,
        "historical_reference": historical,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "xgboost": xgb.__version__,
            "xgb_params": xgb_classifier(args.seed).get_params(),
        },
        "figures": figures,
        "seed": args.seed,
    }
    with open(out / "validation_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    rows = []
    for key, r in summary.items():
        p = r["pooled_outer_test"]
        rows.append(
            {
                "variant": r["variant"],
                "calibration": r["calibration"],
                **{k: round(p[k], 4) for k in
                   ["auroc", "auprc", "accuracy", "precision", "recall", "specificity",
                    "f1", "ef_10pct", "brier", "active_rate"]},
                "mean_fold_auroc": round(r["mean_fold_auroc"], 4),
                "sd_fold_auroc": round(r["sd_fold_auroc"], 4),
                "min_fold_auroc": round(r["min_fold_auroc"], 4),
            }
        )
    tbl = pd.DataFrame(rows).sort_values(["variant", "calibration"])
    tbl.to_csv(out / "performance_summary.csv", index=False)

    recommended = "A_structure_only|none"
    model_card = {
        "model_id": "nlrp3_v2_structure_only",
        "artifact": "results/nlrp3_model_v2/nlrp3_v2_structure_only.joblib",
        "task": "binary classification of NLRP3-related functional activity (IL-1beta release assays)",
        "label": audit["label_rule"],
        "features": "Morgan radius-2 2048 bits + 12 RDKit descriptors, StandardScaler fitted on base-fit molecules only",
        "algorithm": "XGBoost classifier",
        "hyperparameters": xgb_classifier(args.seed).get_params(),
        "training_data": "data/processed/nlrp3_records.csv collapsed to molecules; "
                         f"{final_meta['n_fit_molecules']} molecules for the base model, "
                         f"{final_meta['n_calibration_molecules']} scaffold-disjoint molecules reserved for calibration",
        "split_membership": "results/nlrp3_model_v2/final_model_split_membership.csv",
        "prediction_rule": bundles["structure_only"]["prediction_rule"],
        "recommended_calibration": "none",
        "recommended_threshold": 0.5,
        "held_out_performance": summary[recommended]["pooled_outer_test"],
        "available_calibrators": sorted(bundles["structure_only"]["calibrators"]),
        "software": report["software"],
        "seed": args.seed,
        "alternative_model": {
            "model_id": "nlrp3_v2_assay_conditioned",
            "artifact": "results/nlrp3_model_v2/nlrp3_v2_assay_conditioned.joblib",
            "prediction_rule": bundles["assay_conditioned"]["prediction_rule"],
            "held_out_performance": summary["B_assay_deployment|none"]["pooled_outer_test"],
        },
    }
    with open(out / "model_card.json", "w") as f:
        json.dump(model_card, f, indent=2, default=str)

    print(tbl.to_string(index=False))
    print(f"\nWrote results to {out}")


if __name__ == "__main__":
    main()
