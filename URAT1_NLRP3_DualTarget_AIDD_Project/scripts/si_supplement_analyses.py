#!/usr/bin/env python3
"""SI analyses that do not re-dock the clinical pool or re-lock Π*.

1. NLRP3 assay-context top-1 / top-3 / top-5 shrink-set overlap (ML only).
2. Protocol-table EF 95% CIs from published top-slice hit counts
   (not ranking bootstrap; per-molecule P0–P5 scores are not archived).
3. Optional MCC950@7ALV analog dock if gnina is on PATH / tools/gnina.

Production docking pool remains data/repurposing/screening/docking_pool_p05.csv.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import binomtest, hypergeom

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from prepare_ligands_vina import smiles_to_pdbqt  # noqa: E402
from utils_ml import assay_one_hot_matrix, canonicalize, featurize_smiles  # noqa: E402

OUT_DIR = PROJECT_ROOT / "data" / "si"


def _ensure_gnina_runtime() -> None:
    """Prepend CUDA 12 redistributable lib dirs so the gnina binary can start on CPU."""
    extra = []
    nvidia = Path.home() / ".local" / "lib" / "python3.12" / "site-packages" / "nvidia"
    if nvidia.exists():
        extra.extend(sorted({str(p.parent) for p in nvidia.glob("*/lib/lib*.so*")}))
    nvtx = Path.home() / ".local" / "lib" / "nvtx"
    if nvtx.exists():
        extra.append(str(nvtx))
    if not extra:
        return
    current = os.environ.get("LD_LIBRARY_PATH", "")
    os.environ["LD_LIBRARY_PATH"] = ":".join(extra + ([current] if current else []))


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    return float(len(a & b) / len(union)) if union else float("nan")



def assay_shrink_overlap(model_path: Path) -> dict:
    scores = pd.read_csv(
        PROJECT_ROOT / "data" / "repurposing" / "screening" / "nlrp3_ml_scores_clinical_all.csv",
        low_memory=False,
    )
    prod_pool = pd.read_csv(
        PROJECT_ROOT / "data" / "repurposing" / "screening" / "docking_pool_p05.csv",
        low_memory=False,
    )
    prod_set = set(prod_pool["canonical_smiles"].dropna().astype(str))
    bundle = joblib.load(model_path)
    assay_ids = list(bundle["top_assay_ids"][:5])
    smiles = scores["canonical_smiles"].astype(str).tolist()
    x_mol = featurize_smiles(smiles)
    per_assay = {}
    for aid in assay_ids:
        assay_col = pd.Series([aid] * len(smiles))
        x_assay = assay_one_hot_matrix(assay_col, bundle["top_assay_ids"])
        x_s = bundle["scaler"].transform(np.hstack([x_mol, x_assay]))
        raw = bundle["model"].predict_proba(x_s)[:, 1]
        if bundle.get("calibrator") is not None:
            raw = bundle["calibrator"].predict(raw)
        per_assay[aid] = np.asarray(raw, dtype=float)

    rows = []
    sets: dict[int, set[str]] = {}
    score_frame = scores[["canonical_smiles", "name", "chembl_id", "p_active_nlrp3"]].copy()
    for n in (1, 3, 5):
        stacked = np.vstack([per_assay[aid] for aid in assay_ids[:n]])
        q = stacked.max(axis=0)
        col = f"q_n_top{n}"
        score_frame[col] = q
        mask = q >= 0.5
        selected = set(score_frame.loc[mask, "canonical_smiles"].astype(str))
        sets[n] = selected
        rows.append(
            {
                "n_ensemble_assays": n,
                "assay_ids": ";".join(assay_ids[:n]),
                "n_q_ge_0.5": int(mask.sum()),
                "n_overlap_production_1588": int(len(selected & prod_set)),
                "n_only_this_set": int(len(selected - prod_set)),
                "n_only_production_1588": int(len(prod_set - selected)),
                "jaccard_vs_production_1588": round(_jaccard(selected, prod_set), 4),
                "frac_production_recovered": round(len(selected & prod_set) / max(len(prod_set), 1), 4),
            }
        )

    pairwise = [
        {
            "pair": "top1_vs_top3",
            "jaccard": round(_jaccard(sets[1], sets[3]), 4),
            "n_intersection": int(len(sets[1] & sets[3])),
            "n_union": int(len(sets[1] | sets[3])),
        },
        {
            "pair": "top1_vs_top5",
            "jaccard": round(_jaccard(sets[1], sets[5]), 4),
            "n_intersection": int(len(sets[1] & sets[5])),
            "n_union": int(len(sets[1] | sets[5])),
        },
        {
            "pair": "top3_vs_top5",
            "jaccard": round(_jaccard(sets[3], sets[5]), 4),
            "n_intersection": int(len(sets[3] & sets[5])),
            "n_union": int(len(sets[3] | sets[5])),
        },
    ]

    abs_diff = np.abs(score_frame["q_n_top5"].to_numpy() - scores["p_active_nlrp3"].to_numpy())
    summary = {
        "n_scored": int(len(scores)),
        "production_n": int(len(prod_set)),
        "threshold": 0.5,
        "ensemble_assays_top5": assay_ids,
        "max_abs_q_top5_vs_archived_p_active": float(np.nanmax(abs_diff)),
        "mean_abs_q_top5_vs_archived_p_active": float(np.nanmean(abs_diff)),
        "n_q_top5_disagrees_with_archived_ge_0.5": int(
            ((score_frame["q_n_top5"] >= 0.5) != (scores["p_active_nlrp3"] >= 0.5)).sum()
        ),
        "sets": rows,
        "pairwise_rescore": pairwise,
        "decision": (
            "Keep production docking_pool_p05.csv (n=1588) unchanged. "
            "Assay-count variants are SI only; do not re-dock."
        ),
    }

    out = OUT_DIR / "assay_shrink_overlap"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "assay_top_n_vs_production.csv", index=False)
    pd.DataFrame(pairwise).to_csv(out / "assay_top_n_pairwise.csv", index=False)
    # Compact membership table for SI (not full 8319 score dump)
    membership = pd.DataFrame(
        {
            "canonical_smiles": score_frame["canonical_smiles"],
            "name": score_frame["name"],
            "chembl_id": score_frame["chembl_id"],
            "in_production_1588": score_frame["canonical_smiles"].astype(str).isin(prod_set),
            "q_ge_0.5_top1": score_frame["q_n_top1"] >= 0.5,
            "q_ge_0.5_top3": score_frame["q_n_top3"] >= 0.5,
            "q_ge_0.5_top5": score_frame["q_n_top5"] >= 0.5,
        }
    )
    membership.to_csv(out / "assay_membership_flags.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary



def _ef_ci(k: int, n: int, n_actives: int, n_total: int) -> dict:
    prev = n_actives / n_total
    point = (k / n) / prev if n and prev else float("nan")
    ci = binomtest(k, n).proportion_ci(confidence_level=0.95)
    hyper_p = float(hypergeom.sf(k - 1, n_total, n_actives, n)) if k > 0 else 1.0
    return {
        "hits": k,
        "n_top": n,
        "n_actives": n_actives,
        "n_total": n_total,
        "prevalence": round(prev, 6),
        "ef": round(point, 4),
        "ef_ci95_low": round(ci.low / prev, 4),
        "ef_ci95_high": round(ci.high / prev, 4),
        "hitrate_ci95_low": round(ci.low, 4),
        "hitrate_ci95_high": round(ci.high, 4),
        "hypergeom_p_right": round(hyper_p, 6),
        "interval_method": "Clopper-Pearson 95% CI on published k/n, divided by prevalence; not ranking bootstrap",
    }


def protocol_enrichment_ci() -> dict:
    archived_sum = OUT_DIR / "protocol_enrichment_ci" / "summary.json"
    if archived_sum.exists() and "ranking_files_present" in archived_sum.read_text():
        summary = json.loads(archived_sum.read_text())
        if summary.get("ranking_files_present"):
            return {"summary": summary, "n_rows": 12, "source": "archived bootstrap"}
    n_actives = 469
    n_true = 469 + 4690
    n_rand = 469 + 4690
    # Published hits@1% from docs/PROTOCOL_SELECTION_RESULT.md; EF@5% k reconstructed.
    table = [
        {
            "protocol": "P5",
            "readout": "RTMScore (gnina pose)",
            "true_hits_1pct": (13, 51),
            "true_ef_1pct_published": 2.80,
            "true_ef_5pct_published": 2.44,
            "true_auc_published": 0.590,
            "random_ef_1pct_published": 0.00,
            "random_ef_5pct_published": 1.20,
            "random_auc_published": 0.553,
            "random_hits_1pct": (0, 51),
        },
        {
            "protocol": "P2",
            "readout": "gnina CNNaffinity",
            "true_hits_1pct": (12, 52),
            "true_ef_1pct_published": 2.54,
            "true_ef_5pct_published": 1.88,
            "true_auc_published": 0.580,
            "random_ef_1pct_published": 0.21,
            "random_ef_5pct_published": 0.77,
            "random_auc_published": 0.540,
            "random_hits_1pct": (1, 52),
        },
        {
            "protocol": "P0",
            "readout": "gnina CNNscore (negative control)",
            "true_hits_1pct": (9, 52),
            "true_ef_1pct_published": 1.90,
            "true_ef_5pct_published": 2.39,
            "true_auc_published": 0.647,
            "random_ef_1pct_published": 1.90,
            "random_ef_5pct_published": 2.09,
            "random_auc_published": 0.631,
            "random_hits_1pct": (9, 52),
        },
        {
            "protocol": "P4",
            "readout": "RTMScore (Vina pose)",
            "true_hits_1pct": (3, 50),
            "true_ef_1pct_published": 0.65,
            "true_ef_5pct_published": 1.03,
            "true_auc_published": 0.625,
            "random_ef_1pct_published": 0.00,
            "random_ef_5pct_published": 0.21,
            "random_auc_published": 0.562,
            "random_hits_1pct": (0, 50),
        },
        {
            "protocol": "P1",
            "readout": "Vina affinity",
            "true_hits_1pct": (2, 51),
            "true_ef_1pct_published": 0.43,
            "true_ef_5pct_published": 0.51,
            "true_auc_published": 0.531,
            "random_ef_1pct_published": 0.85,
            "random_ef_5pct_published": 0.81,
            "random_auc_published": 0.630,
            "random_hits_1pct": (4, 51),
        },
        {
            "protocol": "P3",
            "readout": "gnina minimizedAffinity",
            "true_hits_1pct": (2, 52),
            "true_ef_1pct_published": 0.42,
            "true_ef_5pct_published": 0.51,
            "true_auc_published": 0.503,
            "random_ef_1pct_published": 0.63,
            "random_ef_5pct_published": 0.64,
            "random_auc_published": 0.564,
            "random_hits_1pct": (3, 52),
        },
    ]

    rows = []
    for rec in table:
        k_t, n_t = rec["true_hits_1pct"]
        k_r, n_r = rec["random_hits_1pct"]
        true_ci = _ef_ci(k_t, n_t, n_actives, n_true)
        rand_ci = _ef_ci(k_r, n_r, n_actives, n_rand)
        n5 = max(1, math.floor(0.05 * n_true))
        k5_true = int(round(rec["true_ef_5pct_published"] * n5 * (n_actives / n_true)))
        k5_rand = int(round(rec["random_ef_5pct_published"] * n5 * (n_actives / n_rand)))
        true_ci5 = _ef_ci(k5_true, n5, n_actives, n_true)
        rand_ci5 = _ef_ci(k5_rand, n5, n_actives, n_rand)
        rows.append(
            {
                "protocol": rec["protocol"],
                "readout": rec["readout"],
                "true_hits_at_1pct": f"{k_t}/{n_t}",
                "true_EF1pct_published": rec["true_ef_1pct_published"],
                "true_EF1pct_from_counts": true_ci["ef"],
                "true_EF1pct_ci95": f"{true_ci['ef_ci95_low']:.2f}–{true_ci['ef_ci95_high']:.2f}",
                "true_EF1pct_ci95_low": true_ci["ef_ci95_low"],
                "true_EF1pct_ci95_high": true_ci["ef_ci95_high"],
                "true_EF5pct_published": rec["true_ef_5pct_published"],
                "true_EF5pct_k_reconstructed": k5_true,
                "true_EF5pct_n": n5,
                "true_EF5pct_ci95": f"{true_ci5['ef_ci95_low']:.2f}–{true_ci5['ef_ci95_high']:.2f}",
                "true_AUC_published": rec["true_auc_published"],
                "true_AUC_ci95": "unavailable",
                "random_hits_at_1pct": f"{k_r}/{n_r}",
                "random_EF1pct_published": rec["random_ef_1pct_published"],
                "random_EF1pct_from_counts": rand_ci["ef"],
                "random_EF1pct_ci95": f"{rand_ci['ef_ci95_low']:.2f}–{rand_ci['ef_ci95_high']:.2f}",
                "random_EF1pct_ci95_low": rand_ci["ef_ci95_low"],
                "random_EF1pct_ci95_high": rand_ci["ef_ci95_high"],
                "random_EF5pct_published": rec["random_ef_5pct_published"],
                "random_AUC_published": rec["random_auc_published"],
                "random_AUC_ci95": "unavailable",
                "true_hypergeom_p_1pct": true_ci["hypergeom_p_right"],
            }
        )

    summary = {
        "n_actives": n_actives,
        "n_true_benchmark": n_true,
        "n_random_benchmark": n_rand,
        "ranking_files_present": False,
        "auc_ci": "unavailable — per-molecule P0–P5 scores were not archived",
        "ef_ci_method": (
            "Clopper–Pearson 95% CI on the published hits@1% counts, converted to EF "
            "by dividing by prevalence 469/5159. This is not a bootstrap of docking ranks. "
            "EF@5% intervals use k reconstructed as round(EF_published × n_5% × prevalence) "
            "and are labeled as such. AUC confidence intervals cannot be recovered from counts."
        ),
        "random_hits_1pct_source": (
            "P2/P5/P4 taken from published EF@1% (0.21 → 1/52; 0 → 0/n). "
            "P0 uses the published 1.90 matching True (9/52). "
            "P1/P3 k reconstructed from published Random EF@1%."
        ),
    }
    out = OUT_DIR / "protocol_enrichment_ci"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "protocol_ef_ci.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return {"summary": summary, "n_rows": len(rows)}


def try_mcc950_dock() -> dict:
    _ensure_gnina_runtime()
    gnina = PROJECT_ROOT / "tools" / "gnina"
    which = None
    from shutil import which as _which

    if gnina.exists() and gnina.stat().st_mode & 0o111:
        which = str(gnina)
    elif _which("gnina"):
        which = _which("gnina")

    mcc_smi = "CC(C)(O)c1coc(S(=O)(=O)NC(=O)Nc2c3c(cc4c2CCC4)CCC3)c1"
    rec = PROJECT_ROOT / "data" / "structures" / "prepared" / "7ALV_receptor.pdbqt"
    out = OUT_DIR / "mcc950_7alv"
    out.mkdir(parents=True, exist_ok=True)
    note = {
        "ligand": "MCC950",
        "chembl_id": "CHEMBL230208",
        "smiles": mcc_smi,
        "receptor": "7ALV",
        "co_crystal_ligand": "NP3-146 / RM5 (not MCC950)",
        "interpretation": (
            "Pharmacological analog/control dock, not a true self-dock. "
            "RMSD versus RM5 is not a self-dock pass/fail. Does not re-lock Π*."
        ),
        "in_production_1588": False,
        "in_nlrp3_train_chembl230208": False,
        "gnina_binary": which,
        "status": "blocked_no_gnina" if which is None else "ready",
    }
    nl = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "nlrp3_records.csv", low_memory=False)
    note["in_nlrp3_train_chembl230208"] = bool((nl["Molecule ChEMBL ID"].astype(str) == "CHEMBL230208").any())
    can = canonicalize(mcc_smi)
    if can and "canonical_smiles" in nl.columns:
        note["in_nlrp3_train_canonical_smiles"] = bool((nl["canonical_smiles"].astype(str) == can).any())

    if which is None or not rec.exists():
        note["status"] = "blocked_no_gnina" if which is None else "blocked_no_receptor"
        (out / "status.json").write_text(json.dumps(note, indent=2))
        return note

    from run_gnina_batch import run_batch
    import yaml

    lig_dir = out / "ligand"
    lig_dir.mkdir(exist_ok=True)
    pdbqt, err = smiles_to_pdbqt(mcc_smi, "MCC950")
    if pdbqt is None:
        note["status"] = f"ligand_prep_fail:{err}"
        (out / "status.json").write_text(json.dumps(note, indent=2))
        return note
    lig_path = lig_dir / "MCC950.pdbqt"
    lig_path.write_text(pdbqt)
    pd.DataFrame(
        [
            {
                "repurposing_id": "MCC950",
                "canonical_smiles": can or mcc_smi,
                "status": "prepared",
                "pdbqt": str(lig_path),
            }
        ]
    ).to_csv(lig_dir / "ligand_manifest.csv", index=False)

    cfg = yaml.safe_load((PROJECT_ROOT / "config" / "docking_production_p2.yaml").read_text())
    gnina_cfg = dict(cfg["gnina"])
    gnina_cfg["binary"] = which
    target = cfg["targets"]["nlrp3_7alv"]
    dock_dir = out / "dock"
    run_batch(
        lig_dir / "ligand_manifest.csv",
        rec,
        target["center"],
        target["size"],
        dock_dir,
        "7ALV",
        gnina_cfg,
        jobs=1,
        limit=1,
    )
    csv_path = dock_dir / "docking_7alv_gnina.csv"
    note["status"] = "docked" if csv_path.exists() else "dock_failed"
    note["output_csv"] = str(csv_path) if csv_path.exists() else None
    if csv_path.exists():
        note["dock_row"] = pd.read_csv(csv_path).to_dict(orient="records")[0]
    (out / "status.json").write_text(json.dumps(note, indent=2, default=str))
    return note


def main() -> None:
    parser = argparse.ArgumentParser(description="SI analyses without re-locking Π* or replacing 1588")
    parser.add_argument("--model", type=Path, default=PROJECT_ROOT / "data" / "models" / "nlrp3_model.joblib")
    parser.add_argument("--skip-dock", action="store_true")
    parser.add_argument("--only-mcc950", action="store_true", help="Skip ML/CI tables; only try MCC950@7ALV")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.only_mcc950:
        report = {"mcc950_7alv": try_mcc950_dock()}
        (OUT_DIR / "mcc950_7alv" / "status.json").write_text(json.dumps(report["mcc950_7alv"], indent=2, default=str))
        print(json.dumps(report, indent=2, default=str))
        return

    report = {
        "assay_shrink_overlap": assay_shrink_overlap(args.model),
        "protocol_enrichment_ci": protocol_enrichment_ci(),
    }
    if not args.skip_dock:
        report["mcc950_7alv"] = try_mcc950_dock()
    (OUT_DIR / "si_supplement_summary.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps({k: (v if k != "assay_shrink_overlap" else {kk: v[kk] for kk in v if kk != "sets"}) for k, v in report.items()}, indent=2, default=str))


if __name__ == "__main__":
    main()
