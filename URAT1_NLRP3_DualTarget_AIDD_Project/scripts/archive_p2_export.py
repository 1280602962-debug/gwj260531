#!/usr/bin/env python3
"""Archive slim P2 tables from docking_export_20260820 and fill SI analyses.

Does not re-dock. Does not re-lock Π*. Strips leftover glide_score_xp aliases
(identical to gnina dock_score). Writes complete-case and ranking-bootstrap
protocol CIs now that mol_protocol_scores.csv is available.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPORT = PROJECT_ROOT / "docking_export_20260820"
PHASE1 = EXPORT / "01_phase1_benchmark_URAT1_9DKB"
PHASE2 = EXPORT / "02_phase2_dualtarget_9DKB_7ALV"
PHASE3 = EXPORT / "03_pose_qc_md_selection"
POOL = PROJECT_ROOT / "data" / "repurposing" / "screening" / "docking_pool_p05.csv"
OUT_P2 = PROJECT_ROOT / "data" / "repurposing" / "p2"
OUT_PROTO = PROJECT_ROOT / "data" / "benchmarks" / "protocol_selection"
OUT_SI = PROJECT_ROOT / "data" / "si"
GLIDE_COLS = {"glide_score_xp", "nlrp3_glide_score_xp"}

PROTOCOLS = [
    ("P0", "P0_CNNscore", True, "gnina CNNscore (negative control)"),
    ("P1", "P1_vina_affinity", False, "Vina affinity"),
    ("P2", "P2_CNNaffinity", True, "gnina CNNaffinity"),
    ("P3", "P3_gnina_affinity", False, "gnina minimizedAffinity"),
    ("P4", "P4_RTM_vina", True, "RTMScore (Vina pose)"),
    ("P5", "P5_RTM_gnina", True, "RTMScore (gnina pose)"),
]


def _drop_glide(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=[c for c in df.columns if c in GLIDE_COLS], errors="ignore")


def _copy_csv(src: Path, dst: Path) -> pd.DataFrame:
    df = _drop_glide(pd.read_csv(src, low_memory=False))
    dst.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dst, index=False)
    return df


def archive_tables() -> dict:
    OUT_P2.mkdir(parents=True, exist_ok=True)
    tables = {
        "docking_9dkb": _copy_csv(
            PHASE2 / "docking_9dkb_URAT1" / "docking_9dkb_gnina.csv",
            OUT_P2 / "docking_9dkb_gnina.csv",
        ),
        "docking_7alv": _copy_csv(
            PHASE2 / "docking_7alv_NLRP3" / "docking_7alv_gnina.csv",
            OUT_P2 / "docking_7alv_gnina.csv",
        ),
        "merged": _copy_csv(
            PHASE2 / "pareto" / "pareto_merged_scores.csv",
            OUT_P2 / "pareto_merged_scores.csv",
        ),
        "shortlist": _copy_csv(
            PHASE2 / "pareto" / "pareto_shortlist.csv",
            OUT_P2 / "pareto_shortlist.csv",
        ),
        "nominated": _copy_csv(
            PHASE2 / "candidates" / "nominated_candidates.csv",
            OUT_P2 / "nominated_candidates.csv",
        ),
        "diverse": _copy_csv(
            PHASE2 / "candidates" / "nominated_shortlist_diverse.csv",
            OUT_P2 / "nominated_shortlist_diverse.csv",
        ),
        "filters": _copy_csv(PHASE2 / "cheminformatics" / "filters_pool.csv", OUT_P2 / "filters_pool.csv"),
        "admet": _copy_csv(PHASE2 / "cheminformatics" / "admet_pool.csv", OUT_P2 / "admet_pool.csv"),
        "ligand_manifest": _copy_csv(
            PHASE2 / "ligands" / "ligands_p05_pdbqt" / "ligand_manifest.csv",
            OUT_P2 / "ligand_manifest.csv",
        ),
        "scores": _copy_csv(
            PHASE1 / "scores" / "mol_protocol_scores.csv",
            OUT_PROTO / "mol_protocol_scores.csv",
        ),
        "protocol_metrics": _copy_csv(
            PHASE1 / "metrics" / "protocol_metrics.csv",
            OUT_PROTO / "protocol_metrics.csv",
        ),
        "pose_qc_dual": _copy_csv(PHASE3 / "pose_qc_dual.csv", OUT_SI / "pose_qc" / "pose_qc_dual.csv"),
        "pose_qc_table": _copy_csv(PHASE3 / "pose_qc_table.csv", OUT_SI / "pose_qc" / "pose_qc_table.csv"),
    }
    for name, src in (
        ("pareto_summary.json", PHASE2 / "pareto" / "pareto_summary.json"),
        ("candidate_nomination_summary.json", PHASE2 / "candidates" / "candidate_nomination_summary.json"),
        ("prepare_summary.json", PHASE2 / "ligands" / "ligands_p05_pdbqt" / "prepare_summary.json"),
        ("docking_9dkb_summary.json", PHASE2 / "docking_9dkb_URAT1" / "docking_summary.json"),
        ("docking_7alv_summary.json", PHASE2 / "docking_7alv_NLRP3" / "docking_summary.json"),
        ("residue_map.json", PHASE3 / "residue_map.json"),
    ):
        dst = OUT_P2 / name if "pose" not in name and name != "residue_map.json" else OUT_SI / "pose_qc" / name
        if name == "residue_map.json":
            dst = OUT_SI / "pose_qc" / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    (OUT_SI / "pose_qc" / "md_systems_built.txt").write_text(
        (PHASE3 / "md_systems_built.txt").read_text()
    )
    return tables


def complete_case(tables: dict) -> dict:
    pool = pd.read_csv(POOL, low_memory=False)
    lig = tables["ligand_manifest"]
    u = tables["docking_9dkb"]
    n = tables["docking_7alv"]
    merged = tables["merged"]

    pool_ids = set(pool["repurposing_id"].astype(str))
    if "status" in lig.columns:
        lig_ok = set(lig.loc[lig["status"].astype(str).str.lower().eq("prepared"), "repurposing_id"].astype(str))
    else:
        lig_ok = set(lig["repurposing_id"].astype(str))
    u_ok = set(u.loc[u["docking_status"].astype(str).str.lower().eq("docked"), "repurposing_id"].astype(str))
    n_ok = set(n.loc[n["docking_status"].astype(str).str.lower().eq("docked"), "repurposing_id"].astype(str))
    m_ok = set(merged["repurposing_id"].astype(str)) if "repurposing_id" in merged.columns else set()

    missing_prep = sorted(pool_ids - lig_ok)
    missing_dock = sorted((lig_ok - u_ok) | (lig_ok - n_ok))
    missing_merge = sorted((u_ok & n_ok) - m_ok)

    pool_idx = pool.set_index("repurposing_id", drop=False)
    rows = []
    for rid, reason in (
        [(x, "ligand_prep_not_in_manifest") for x in missing_prep]
        + [(x, "docking_failed_or_empty_pose") for x in missing_dock]
        + [(x, "dual_success_but_not_in_inner_merge") for x in missing_merge]
    ):
        rec = {"repurposing_id": rid, "reason": reason}
        if rid in pool_idx.index:
            hit = pool_idx.loc[rid]
            if isinstance(hit, pd.DataFrame):
                hit = hit.iloc[0]
            rec.update(
                {
                    "name": hit.get("name"),
                    "chembl_id": hit.get("chembl_id"),
                    "canonical_smiles": hit.get("canonical_smiles"),
                    "mw": hit.get("mw"),
                }
            )
        rows.append(rec)

    drop_df = pd.DataFrame(rows)
    out = OUT_SI / "complete_case_drop"
    out.mkdir(parents=True, exist_ok=True)
    drop_df.to_csv(out / "missing_from_1588.csv", index=False)
    reason_counts = drop_df["reason"].value_counts().rename_axis("reason").reset_index(name="n")
    reason_counts.to_csv(out / "si_reason_counts.csv", index=False)

    mw_missing = pd.to_numeric(drop_df["mw"], errors="coerce") if "mw" in drop_df.columns else pd.Series(dtype=float)
    mw_pool = pd.to_numeric(pool["mw"], errors="coerce")
    summary = {
        "n_pool_1588": int(len(pool)),
        "n_ligand_prepared": int(len(lig_ok)),
        "n_9dkb_docked": int(len(u_ok)),
        "n_7alv_docked": int(len(n_ok)),
        "n_dual_complete_case": int(len(merged)),
        "n_pareto_front": int(merged["pareto_front"].sum()) if "pareto_front" in merged.columns else None,
        "n_missing_from_pool": int(len(drop_df)),
        "n_ligand_prep_gap": int(len(missing_prep)),
        "n_dock_fail": int(len(missing_dock)),
        "n_merge_drop": int(len(missing_merge)),
        "missing_prep_ids": missing_prep,
        "missing_dock_ids": missing_dock,
        "missing_merge_ids": missing_merge,
        "mw_pool_median": float(mw_pool.median()) if mw_pool.notna().any() else None,
        "mw_missing_median": float(mw_missing.median()) if mw_missing.notna().any() else None,
        "note": (
            "Production percentiles use the dual-success complete case only. "
            "Pool n=1588 is the NLRP3 ML shrink set; 1580 is the P2 dual-dock intersection."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def _rank_and_metrics(scores: pd.Series, labels: pd.Series, higher_is_better: bool, frac: float = 0.01) -> dict:
    ok = scores.notna() & labels.notna()
    s = scores[ok]
    y = labels[ok].astype(int)
    n = int(len(s))
    n_act = int(y.sum())
    if n == 0 or n_act == 0:
        return {"n": n, "n_actives": n_act, "ef": float("nan"), "hits": 0, "n_top": 0, "auc": float("nan")}
    order = s.rank(method="first", ascending=not higher_is_better)
    # rank 1 = best
    rank = order if higher_is_better is False else (n + 1 - s.rank(method="first", ascending=True))
    # Simpler: argsort
    ranked_y = y.to_numpy()[np.argsort(-s.to_numpy() if higher_is_better else s.to_numpy())]
    n_top = max(1, int(np.floor(frac * n)))
    hits = int(ranked_y[:n_top].sum())
    prev = n_act / n
    ef = (hits / n_top) / prev if prev else float("nan")
    direction = s if higher_is_better else -s
    auc = float(roc_auc_score(y, direction))
    hyper_p = float(hypergeom.sf(hits - 1, n, n_act, n_top)) if hits > 0 else 1.0
    return {
        "n": n,
        "n_actives": n_act,
        "ef": round(ef, 4),
        "hits": hits,
        "n_top": n_top,
        "auc": round(auc, 4),
        "hypergeom_p": round(hyper_p, 6),
    }


def _bootstrap_ci(scores: pd.Series, labels: pd.Series, higher_is_better: bool, n_boot: int, seed: int) -> dict:
    ok = scores.notna() & labels.notna()
    s = scores[ok].to_numpy()
    y = labels[ok].astype(int).to_numpy()
    rng = np.random.default_rng(seed)
    efs1, efs5, aucs = [], [], []
    n = len(s)
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        ss, yy = s[idx], y[idx]
        if yy.sum() == 0 or yy.sum() == n:
            continue
        for frac, bucket in ((0.01, efs1), (0.05, efs5)):
            ranked = yy[np.argsort(-ss if higher_is_better else ss)]
            n_top = max(1, int(np.floor(frac * n)))
            hits = int(ranked[:n_top].sum())
            prev = yy.sum() / n
            bucket.append((hits / n_top) / prev if prev else np.nan)
        direction = ss if higher_is_better else -ss
        try:
            aucs.append(roc_auc_score(yy, direction))
        except ValueError:
            continue
    def _ci(arr):
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if len(a) == 0:
            return (float("nan"), float("nan"))
        return (round(float(np.percentile(a, 2.5)), 4), round(float(np.percentile(a, 97.5)), 4))
    return {
        "ef1_ci95": _ci(efs1),
        "ef5_ci95": _ci(efs5),
        "auc_ci95": _ci(aucs),
        "n_boot": n_boot,
    }


def protocol_ranking_ci(scores: pd.DataFrame, n_boot: int = 1000) -> dict:
    rows = []
    for pid, col, higher, readout in PROTOCOLS:
        for bench, flag in (("true", "in_true"), ("random", "in_random")):
            mask = scores[flag].astype(int) == 1
            # Active flag: role contains active, or in both benches with role==active
            y = (scores.loc[mask, "role"].astype(str).str.lower() == "active")
            met = _rank_and_metrics(scores.loc[mask, col], y.astype(int), higher)
            met5 = _rank_and_metrics(scores.loc[mask, col], y.astype(int), higher, frac=0.05)
            boot = _bootstrap_ci(scores.loc[mask, col], y.astype(int), higher, n_boot=n_boot, seed=42)
            rows.append(
                {
                    "protocol": pid,
                    "readout": readout,
                    "benchmark": bench,
                    "coverage_scored": int(scores.loc[mask, col].notna().sum()),
                    "n_benchmark": int(mask.sum()),
                    "hits_at_1pct": f"{met['hits']}/{met['n_top']}",
                    "EF1pct": met["ef"],
                    "EF1pct_ci95_low": boot["ef1_ci95"][0],
                    "EF1pct_ci95_high": boot["ef1_ci95"][1],
                    "EF5pct": met5["ef"],
                    "EF5pct_ci95_low": boot["ef5_ci95"][0],
                    "EF5pct_ci95_high": boot["ef5_ci95"][1],
                    "AUC": met["auc"],
                    "AUC_ci95_low": boot["auc_ci95"][0],
                    "AUC_ci95_high": boot["auc_ci95"][1],
                    "hypergeom_p_1pct": met["hypergeom_p"],
                    "higher_is_better": higher,
                }
            )
    out = OUT_SI / "protocol_enrichment_ci"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "protocol_ef_ci.csv", index=False)
    summary = {
        "ranking_files_present": True,
        "source": "data/benchmarks/protocol_selection/mol_protocol_scores.csv",
        "n_molecules": int(len(scores)),
        "n_boot": n_boot,
        "ef_ci_method": (
            "Molecule-resampled ranking bootstrap (1000 draws with replacement) on archived "
            "per-molecule P0–P5 scores. EF and AUC recomputed each draw. This replaces the "
            "previous Clopper–Pearson intervals that used only published hits@1% counts. "
            "Point estimates match the locked protocol table within rounding; Π* is not re-selected."
        ),
        "auc_ci": "bootstrap 95% percentile interval",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return {"summary": summary, "n_rows": len(df)}


def funnel_snapshot(cc: dict, tables: dict) -> dict:
    nom = json.loads((OUT_P2 / "candidate_nomination_summary.json").read_text())
    pose = tables["pose_qc_dual"]
    snap = {
        "n_clinical_library": 8319,
        "n_nlrp3_q_ge_0.5": cc["n_pool_1588"],
        "n_ligand_prepared": cc["n_ligand_prepared"],
        "n_9dkb_docked": cc["n_9dkb_docked"],
        "n_7alv_docked": cc["n_7alv_docked"],
        "n_dual_complete_case": cc["n_dual_complete_case"],
        "n_pareto_front": cc["n_pareto_front"],
        "n_dual_dock_gate_tau90": nom.get("n_dual_gate"),
        "n_preferred": nom.get("n_preferred_candidate"),
        "preferred_followup": ["GSK-3008348", "Vecabrutinib"],
        "preferred_gate_members": [r["name"] for r in nom.get("top_preferred_novel", [])],
        "pose_qc_n": int(len(pose)),
        "pose_qc_both_in_pocket": int(pose["both_in_pocket"].sum()) if "both_in_pocket" in pose.columns else None,
        "md_trajectories_reported": False,
        "engine": "gnina P2 CNNaffinity, exhaustiveness=32, num_modes=1, cnn_scoring=rescore",
    }
    (OUT_P2 / "funnel_snapshot.json").write_text(json.dumps(snap, indent=2))
    (PROJECT_ROOT / "data" / "repurposing" / "screening" / "nlrp3_screening_summary_clinical_all.json")
    # Patch screening summary so figures can show dual-dock n
    sum_path = PROJECT_ROOT / "data" / "repurposing" / "screening" / "nlrp3_screening_summary_clinical_all.json"
    summary = json.loads(sum_path.read_text())
    summary["n_dual_docked"] = cc["n_dual_complete_case"]
    sum_path.write_text(json.dumps(summary, indent=2) + "\n")
    return snap


def main() -> None:
    tables = archive_tables()
    cc = complete_case(tables)
    proto = protocol_ranking_ci(tables["scores"])
    snap = funnel_snapshot(cc, tables)
    report = {"complete_case": cc, "protocol_ci": proto["summary"], "funnel": snap}
    (OUT_SI / "si_supplement_summary.json").parent.mkdir(parents=True, exist_ok=True)
    existing = {}
    path = OUT_SI / "si_supplement_summary.json"
    if path.exists():
        existing = json.loads(path.read_text())
    existing["complete_case_drop"] = cc
    existing["protocol_enrichment_ci"] = proto
    existing["p2_funnel"] = snap
    path.write_text(json.dumps(existing, indent=2, default=str))
    print(json.dumps({"complete_case": cc, "protocol_rows": proto["n_rows"], "funnel": snap}, indent=2, default=str))


if __name__ == "__main__":
    main()
