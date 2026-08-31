#!/usr/bin/env python3
"""Fill remaining P0/P1 analysis tables from frozen scores (no new docking).

Writes, under data/jcim_strengthen_t0t1_v0/:
  tables/endpoint_hierarchy_v1.csv
  tables/frozen_vs_holdout_v1.csv
  tables/wrong_pocket_paired_delta_bootstrap_v1.csv
  tables/pocket_matched_vs_best_descriptor_delta_v1.csv
  tables/pose_fairness_channels_v1.csv
  tables/ligand_ml_scaffold_vs_random_v1.csv
  tables/ranking_top10_vina_mean_exploratory_v1.csv
  analysis/P0_MISSING_CONTENT_VERDICT_V1.md

Paired Δ uses the same ligand-bootstrap protocol as the rest of the paper
(B=2000, seed 20260729). Point estimates are checksummed against the frozen
pocket-matched / wrong-pocket / descriptor CSVs; the script exits nonzero if
they disagree.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]

N_BOOT = 2000
SEED = 20260729

# pair -> (score table, pocket-A vina col, pocket-B vina col)
# Affinity columns are Vina kcal/mol (more negative = better); flipped to S = -E.
SPEC = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        cls="class",
        lig="ligand",
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        cls="class",
        lig="ligand",
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        cls="class",
        lig="ligand",
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        cls="class",
        lig="ligand",
    ),
}


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def assemble(pair: str, cfg: dict) -> list[dict]:
    rows = load_csv(ROOT / cfg["scores"])
    out = []
    for r in rows:
        a, b = fnum(r.get(cfg["vina_a"])), fnum(r.get(cfg["vina_b"]))
        if a is None or b is None:
            continue
        out.append(
            {
                "pair": pair,
                "ligand": r.get(cfg["lig"]),
                "cls": r.get(cfg["cls"]),
                "vina_A": -a,
                "vina_B": -b,
            }
        )
    return out

OUT = ROOT / "data" / "jcim_strengthen_t0t1_v0"
TAB = OUT / "tables"
AN = OUT / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
AN.mkdir(parents=True, exist_ok=True)

PAIR_ORDER = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
HOLD_PAIRS = ["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
BEST_DESC = {
    "EGFR/HER2": "clogp",
    "AChE/BChE": "tpsa",
    "PIK3CA/PIK3CB": "heavy",
    "PIK3CA/mTOR": "heavy",
}
ASSEMBLED = {
    "EGFR/HER2": ROOT / "data/jcim_bench_v0/tables/assembled_EGFR_HER2.csv",
    "AChE/BChE": ROOT / "data/jcim_bench_v0/tables/assembled_AChE_BChE.csv",
    "PIK3CA/PIK3CB": ROOT / "data/jcim_bench_v0/tables/assembled_PIK3CA_PIK3CB.csv",
    "PIK3CA/mTOR": ROOT / "data/jcim_bench_v0/tables/assembled_PIK3CA_mTOR.csv",
}

# Frozen checksums (unrounded CSV values used in Table 2 / Fig 6).
EXPECTED_MATCHED = {
    "EGFR/HER2": (0.6664, 0.4297, 0.4297),
    "AChE/BChE": (0.6504, 0.6058, 0.6058),
    "PIK3CA/PIK3CB": (0.6905, 0.5, 0.5),
    "PIK3CA/mTOR": (0.7143, 0.6921, 0.6921),
}
EXPECTED_WRONG = {
    "EGFR/HER2": (0.6983, 0.26, 0.26),
    "AChE/BChE": (0.4444, 0.5582, 0.4444),
    "PIK3CA/PIK3CB": (0.6442, 0.3489, 0.3489),
    "PIK3CA/mTOR": (0.7103, 0.6019, 0.6019),
}
EXPECTED_HOLD_MATCHED = {
    "AChE/BChE": (0.635, 0.6175, 0.6175),
    "PIK3CA/PIK3CB": (0.7658, 0.425, 0.425),
    "PIK3CA/mTOR": (0.86, 0.765, 0.765),
}
EXPECTED_HOLD_WRONG = {
    "AChE/BChE": (0.6425, 0.6525, 0.6425),
    "PIK3CA/PIK3CB": (0.6395, 0.52, 0.52),
    "PIK3CA/mTOR": (0.7875, 0.8575, 0.7875),
}
EXPECTED_DESC = {
    "EGFR/HER2": 0.4821,
    "AChE/BChE": 0.7333,
    "PIK3CA/PIK3CB": 0.6217,
    "PIK3CA/mTOR": 0.463,
}


def load_csv(p: Path) -> list[dict]:
    with p.open() as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    fields = fields or list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def r4(x) -> float:
    return round(float(x), 4)


def eq4(a, b, msg: str, errors: list[str], tol: float = 5.5e-5) -> None:
    if abs(float(a) - float(b)) > tol:
        errors.append(f"{msg}: {a} != {b}")


def directional(recs: list[dict], key_da: str, key_db: str):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc([r[key_da] for r in D], [r[key_da] for r in A])
    db = auroc([r[key_db] for r in D], [r[key_db] for r in B])
    return da, db, min(da, db), len(D), len(A), len(B)


def boot_paired(
    recs: list[dict],
    key_m_da: str,
    key_m_db: str,
    key_w_da: str,
    key_w_db: str,
    n_boot: int = N_BOOT,
    seed: int = SEED,
):
    usable = [
        r
        for r in recs
        if r["cls"] in ("dual", "A_only", "B_only")
        and key_m_da in r
        and key_m_db in r
        and key_w_da in r
        and key_w_db in r
        and r[key_m_da] is not None
        and r[key_m_db] is not None
    ]
    da_m, db_m, sm_m, nD, nA, nB = directional(usable, key_m_da, key_m_db)
    da_w, db_w, sm_w, _, _, _ = directional(usable, key_w_da, key_w_db)
    # Report Δ from the same 4-decimal AUROCs used in Table 2 / Fig 6 so the
    # figure arithmetic matches the frozen CSVs (avoid round(a-b) ≠ round(a)-round(b)).
    point = r4(r4(sm_m) - r4(sm_w))
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    deltas, sm_ms, sm_ws = [], [], []
    for _ in range(n_boot):
        ii = rng.choice(idx, size=len(idx), replace=True)
        sub = [usable[int(i)] for i in ii]
        _, _, m, *_ = directional(sub, key_m_da, key_m_db)
        _, _, w, *_ = directional(sub, key_w_da, key_w_db)
        if m != m or w != w:
            continue
        sm_ms.append(m)
        sm_ws.append(w)
        deltas.append(m - w)
    if len(deltas) < n_boot // 2:
        raise RuntimeError("too few successful bootstrap replicates")
    dlo, dhi = np.percentile(deltas, [2.5, 97.5])
    mlo, mhi = np.percentile(sm_ms, [2.5, 97.5])
    wlo, whi = np.percentile(sm_ws, [2.5, 97.5])
    return {
        "n_dual": nD,
        "n_A_only": nA,
        "n_B_only": nB,
        "n_boot_ok": len(deltas),
        "matched_D_vs_A": r4(da_m),
        "matched_D_vs_B": r4(db_m),
        "matched_summary_min": r4(sm_m),
        "matched_ci_lo": r4(mlo),
        "matched_ci_hi": r4(mhi),
        "wrong_D_vs_A": r4(da_w),
        "wrong_D_vs_B": r4(db_w),
        "wrong_summary_min": r4(sm_w),
        "wrong_ci_lo": r4(wlo),
        "wrong_ci_hi": r4(whi),
        "delta_matched_minus_wrong": r4(point),
        "delta_boot_mean": r4(float(np.mean(deltas))),
        "delta_ci_lo": r4(dlo),
        "delta_ci_hi": r4(dhi),
        "ci_excludes_zero": bool(dlo > 0 or dhi < 0),
        "point_matched_gt_wrong": bool(sm_m > sm_w),
    }


def load_holdout() -> dict[str, list[dict]]:
    rows = load_csv(ROOT / "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv")
    packs: dict[str, list[dict]] = {p: [] for p in HOLD_PAIRS}
    for r in rows:
        pair = r["pair"]
        if pair not in packs:
            continue
        a, b = r.get("vina_A"), r.get("vina_B")
        if a in ("", None) or b in ("", None):
            continue
        packs[pair].append(
            {
                "ligand": r["ligand"],
                "cls": r["cls"],
                "vina_A": float(a),
                "vina_B": float(b),
                "heavy": float(r["heavy"]) if r.get("heavy") else None,
                "mw": float(r["mw"]) if r.get("mw") else None,
                "clogp": float(r["clogp"]) if r.get("clogp") else None,
                "tpsa": float(r["tpsa"]) if r.get("tpsa") else None,
            }
        )
    return packs


def attach_descriptors(packs: dict[str, list[dict]]) -> dict[str, list[dict]]:
    out = {}
    for pair, recs in packs.items():
        desc = {r["ligand"]: r for r in load_csv(ASSEMBLED[pair])}
        merged = []
        for r in recs:
            d = desc.get(r["ligand"])
            if not d:
                continue
            rec = dict(r)
            for k in ("heavy", "mw", "clogp", "tpsa"):
                rec[k] = float(d[k]) if d.get(k) not in ("", None) else None
            merged.append(rec)
        out[pair] = merged
    return out


def write_endpoint_hierarchy() -> list[dict]:
    rows = [
        {
            "role": "primary",
            "endpoint": "pocket-matched Vina summary_min at unified θ=6.0",
            "definition": "min(AUROC(D vs A_only; S_B), AUROC(D vs B_only; S_A)); S=-E_Vina",
            "panel": "K=4 main panels (PM48 for PIK3CA/mTOR)",
            "reported_in": "Table 2; Figure 3",
            "note": "Single primary. PM48 is the main PIK3CA/mTOR panel, not a secondary metric.",
        },
        {
            "role": "pre-specified secondary",
            "endpoint": "directional arms D/A and D/B",
            "definition": "same pocket-matched protocol, reported separately",
            "panel": "K=4 main panels",
            "reported_in": "Table 2; Figure 4A",
            "note": "Weak-arm inspection; not a second primary.",
        },
        {
            "role": "pre-specified secondary",
            "endpoint": "RTMScore pocket-matched summary_min (best-of-9)",
            "definition": "same pocket-matched protocol on max RTM over the 9 Vina poses",
            "panel": "K=4 main panels",
            "reported_in": "Figure 3; Table S6",
            "note": "Qualitative channel contrast, not a three-engine competition.",
        },
        {
            "role": "pre-specified secondary",
            "endpoint": "GNINA CNN pocket-matched summary_min (best-of-9)",
            "definition": "same protocol on max CNNscore over the 9 Vina poses",
            "panel": "K=4 main panels",
            "reported_in": "Figure 3; Table S15",
            "note": "Fair pose coverage vs RTM. mode-1 retained as historical contrast.",
        },
        {
            "role": "pre-specified secondary",
            "endpoint": "strongest trivial descriptor summary_min",
            "definition": "same AUROC protocol with heavy/mw/clogp/tpsa replacing S",
            "panel": "K=4 main panels",
            "reported_in": "Figure 3; Figure 4B",
            "note": "Paired Δ(Vina−descriptor) CIs in Table S19 / Figure S3C.",
        },
        {
            "role": "robustness",
            "endpoint": "unified θ grid {5.5, 6.0, 6.5, strict 6.5/5.5}",
            "definition": "relabel then recompute pocket-matched Vina summary_min",
            "panel": "K=4 main panels",
            "reported_in": "Table S4; Figure S1A",
            "note": "Not a second primary. Underpowered cells are marked, not ranked.",
        },
        {
            "role": "robustness",
            "endpoint": "PM110 expansion",
            "definition": "same protocol on the PM48 superset",
            "panel": "PIK3CA/mTOR only",
            "reported_in": "Figure S1C",
            "note": "Stability check, not independent validation.",
        },
        {
            "role": "robustness",
            "endpoint": "exhaustiveness 8 vs 16 on PM48",
            "definition": "same pocket-matched Vina on frozen E=8 scores",
            "panel": "PIK3CA/mTOR PM48",
            "reported_in": "Figure S1D",
            "note": "E=16 remains primary for this pair.",
        },
        {
            "role": "robustness",
            "endpoint": "unused-pool holdout pocket-matched Vina",
            "definition": "same protocol; dual/A_only/B_only = 20/20/20; seed 20260731",
            "panel": "AChE/BChE, PIK3CA/PIK3CB, PIK3CA/mTOR only",
            "reported_in": "Table S8 / Table S16; Figure 5A",
            "note": "EGFR/HER2 not eligible. Same ChEMBL batch, not cross-database validation.",
        },
        {
            "role": "robustness",
            "endpoint": "crystal swap on PM48 (4JPS, 5DXT, 4JSX)",
            "definition": "one pocket replaced; other pocket kept at frozen scores",
            "panel": "PIK3CA/mTOR PM48",
            "reported_in": "Table S9; Figure 5B",
            "note": "Pre-specified: true gene, non-chimera, cognate best-of-9 < 2 Å. 3T8M excluded.",
        },
        {
            "role": "robustness",
            "endpoint": "wrong-pocket control and paired Δ(matched−wrong)",
            "definition": "swap pockets then same summary_min; Δ on the same bootstrap sample",
            "panel": "K=4 main + three holdout pairs",
            "reported_in": "Figure 6; Table S17; Figure S3A–B",
            "note": "Point AUROCs already in Table S6 / S8; this table adds paired CIs.",
        },
        {
            "role": "exploratory",
            "endpoint": "ECFP4 logistic GroupKFold vs random KFold",
            "definition": "chemotype–label association, not pocket physics",
            "panel": "K=4 main panels",
            "reported_in": "Figure 7A; Table S20",
            "note": "Scaffold split is primary ML readout. Random split is leakage check only.",
        },
        {
            "role": "exploratory",
            "endpoint": "contact_count (not PLIF)",
            "definition": "ligand heavy atoms within 4.0 Å of receptor heavy atoms, mode-1 poses",
            "panel": "holdout only",
            "reported_in": "Table S11; Figure 6D",
            "note": "Scoring-free burial proxy. Does not claim a residue-level mechanism.",
        },
        {
            "role": "exploratory",
            "endpoint": "top-10 hard-negative count on vina_mean",
            "definition": "count of A_only+B_only among the 10 highest vina_mean ligands",
            "panel": "K=4 main panels",
            "reported_in": "Table S21",
            "note": "Ranking readout on pooled vina_mean, NOT the Table 2 pocket-matched metric.",
        },
        {
            "role": "not used as primary",
            "endpoint": "pooled vina_mean directional AUROC",
            "definition": "legacy forest arm; EGFR 0.2824 is not Table 2 0.4297",
            "panel": "K=4 main panels",
            "reported_in": "Table S6 (pooled_mean row)",
            "note": "Kept as aggregation contrast only.",
        },
    ]
    write_csv(TAB / "endpoint_hierarchy_v1.csv", rows)
    return rows


def write_pose_fairness() -> list[dict]:
    rows = [
        {
            "channel": "Vina 1.2.7",
            "receptor_box_ligand": "frozen (same for all channels)",
            "poses_generated": 9,
            "poses_scored": "Vina best pose (mode 1 = most negative E)",
            "per_pocket_aggregation": "S = -E_Vina",
            "used_in_Table2": "yes (primary)",
            "note": "n_modes=9, energy_range=3, seed 20260727",
        },
        {
            "channel": "RTMScore rtmscore_model1",
            "receptor_box_ligand": "frozen (same for all channels)",
            "poses_generated": 9,
            "poses_scored": "all 9 Vina poses",
            "per_pocket_aggregation": "max RTM per pocket",
            "used_in_Table2": "no (secondary)",
            "note": "Fair pose coverage by construction.",
        },
        {
            "channel": "GNINA 1.3.2 CNN mode-1",
            "receptor_box_ligand": "frozen (same for all channels)",
            "poses_generated": 9,
            "poses_scored": "Vina mode 1 only",
            "per_pocket_aggregation": "CNNscore after --minimize",
            "used_in_Table2": "no (historical)",
            "note": "Asymmetric vs RTM. Retained as Table S14/S15 contrast.",
        },
        {
            "channel": "GNINA 1.3.2 CNN best-of-9",
            "receptor_box_ligand": "frozen (same for all channels)",
            "poses_generated": 9,
            "poses_scored": "all 9 Vina poses",
            "per_pocket_aggregation": "max CNNscore per pocket",
            "used_in_Table2": "no (secondary; pose-fair vs RTM)",
            "note": "2026-08-24 rescore of frozen poses. best9−mode01 is −0.04 to +0.08, not vs Vina.",
        },
    ]
    write_csv(TAB / "pose_fairness_channels_v1.csv", rows)
    return rows


def write_frozen_vs_holdout(errors: list[str]) -> list[dict]:
    theta = {
        r["pair"]: r
        for r in load_csv(ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv")
        if r["label_rule"] == "theta_6.0"
    }
    hold = {
        r["pair"]: r
        for r in load_csv(ROOT / "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv")
        if r["variant"] == "pocket_matched_vina"
    }
    rows = []
    for pair in PAIR_ORDER:
        t = theta[pair]
        da, db, sm = float(t["auroc_D_vs_A"]), float(t["auroc_D_vs_B"]), float(t["pocket_matched_summary_min"])
        eq4(da, EXPECTED_MATCHED[pair][0], f"frozen {pair} DA", errors)
        eq4(db, EXPECTED_MATCHED[pair][1], f"frozen {pair} DB", errors)
        eq4(sm, EXPECTED_MATCHED[pair][2], f"frozen {pair} summary_min", errors)
        rows.append(
            {
                "pair": pair,
                "set": "main_panel",
                "eligible": "yes",
                "n_dual": t["n_dual"],
                "n_A_only": t["n_A_only"],
                "n_B_only": t["n_B_only"],
                "auroc_D_vs_A": t["auroc_D_vs_A"],
                "auroc_D_vs_B": t["auroc_D_vs_B"],
                "summary_min": t["pocket_matched_summary_min"],
                "summary_min_ci_lo": t["ci_lo"],
                "summary_min_ci_hi": t["ci_hi"],
                "delta_holdout_minus_main": "",
                "note": "Table 2 primary (unified θ=6.0 Vina)",
            }
        )
        if pair == "EGFR/HER2":
            rows.append(
                {
                    "pair": pair,
                    "set": "unused_pool_holdout",
                    "eligible": "no",
                    "n_dual": "",
                    "n_A_only": "",
                    "n_B_only": "",
                    "auroc_D_vs_A": "",
                    "auroc_D_vs_B": "",
                    "summary_min": "",
                    "summary_min_ci_lo": "",
                    "summary_min_ci_hi": "",
                    "delta_holdout_minus_main": "",
                    "note": "Not eligible: strict B_only=7; leftover pool cannot support 20/20/20",
                }
            )
            continue
        h = hold[pair]
        hsm = float(h["summary_min"])
        eq4(hsm, EXPECTED_HOLD_MATCHED[pair][2], f"holdout {pair} summary_min", errors)
        rows.append(
            {
                "pair": pair,
                "set": "unused_pool_holdout",
                "eligible": "yes",
                "n_dual": h["n_dual"],
                "n_A_only": h["n_A_only"],
                "n_B_only": h["n_B_only"],
                "auroc_D_vs_A": h["auroc_D_vs_A"],
                "auroc_D_vs_B": h["auroc_D_vs_B"],
                "summary_min": h["summary_min"],
                "summary_min_ci_lo": h["summary_min_ci_lo"],
                "summary_min_ci_hi": h["summary_min_ci_hi"],
                "delta_holdout_minus_main": r4(hsm - sm),
                "note": h.get("note", ""),
            }
        )
    write_csv(TAB / "frozen_vs_holdout_v1.csv", rows)
    return rows


def write_ml_merge() -> list[dict]:
    scaf = load_csv(TAB / "ligand_ml_baseline_scaffold_cv_v1.csv")
    rand = load_csv(TAB / "ligand_ml_baseline_random_cv_v1.csv")
    rand_by = {(r["pair"], r["contrast"]): r for r in rand}
    rows = []
    for s in scaf:
        r = rand_by[(s["pair"], s["contrast"])]
        d = r4(float(r["auroc_ml"]) - float(s["auroc_ml"]))
        rows.append(
            {
                "pair": s["pair"],
                "contrast": s["contrast"],
                "n": s["n"],
                "n_scaffolds": s["n_scaffolds"],
                "n_splits": s["n_splits"],
                "auroc_scaffold_GroupKFold": s["auroc_ml"],
                "auroc_random_StratifiedKFold": r["auroc_ml"],
                "delta_random_minus_scaffold": d,
                "auroc_dock_pocket_matched": s["auroc_dock_pocket_matched"],
                "delta_scaffold_ml_minus_dock": s["delta_ml_minus_dock"],
                "note": "scaffold is the primary ML readout; random is a leakage check, not a hunt for a larger gap",
            }
        )
    write_csv(TAB / "ligand_ml_scaffold_vs_random_v1.csv", rows)
    return rows


def write_top10() -> list[dict]:
    src = load_csv(ROOT / "data/jcim_bench_v0/tables/top10_hardneg_bootstrap_v1.csv")
    rows = []
    for r in src:
        if r["arm"] != "vina_mean":
            continue
        rows.append(
            {
                **r,
                "metric_family": "pooled vina_mean (NOT pocket-matched Table 2)",
                "note": (
                    f"Exploratory ranking readout on vina_mean. "
                    f"{int(r['n_top10_hardneg'])}/10 hard-negatives in the top 10."
                ),
            }
        )
    write_csv(TAB / "ranking_top10_vina_mean_exploratory_v1.csv", rows)
    return rows


def write_verdict(
    delta_rows: list[dict],
    desc_rows: list[dict],
    hold_rows: list[dict],
    ml_rows: list[dict],
) -> None:
    lines = [
        "# P0 missing-content fill (frozen scores only)",
        "",
        "No new docking. Bootstrap B=2000, seed 20260729, ligand units.",
        "Point AUROCs are checksummed against `unified_threshold_sensitivity_v2.csv` (θ=6.0),",
        "`pocket_matched_directional_v1.csv`, `holdout_pocket_matched_v1.csv`, and",
        "`forest_summary_min_ci_v1.csv`.",
        "",
        "## What was filled",
        "",
        "1. Endpoint hierarchy (primary / pre-specified secondary / robustness / exploratory).",
        "2. Formal frozen vs holdout table; EGFR/HER2 holdout = not eligible.",
        "3. Paired Δ = summary_min(matched) − summary_min(wrong) on the **same** bootstrap sample.",
        "4. Pose-fairness table (9 Vina poses; RTM and GNINA best-of-9; GNINA mode-1 historical).",
        "5. Crystal-swap criteria remain pre-specified (true gene, non-chimera, cognate best-of-9 < 2 Å; 3T8M excluded).",
        "6. Scaffold vs random ECFP4 merged; leakage is small — not a hunt for a leakier split.",
        "7. Pocket-matched (not vina_mean) paired Δ vs the strongest trivial descriptor.",
        "8. Exploratory top-10 hard-negative counts on **vina_mean** (not Table 2).",
        "",
        "## Paired Δ (matched − wrong)",
        "",
        "| set | pair | matched | wrong | Δ | 95% CI | excludes 0? |",
        "|-----|------|--------:|------:|--:|--------|-------------|",
    ]
    for r in delta_rows:
        lines.append(
            f"| {r['set']} | {r['pair']} | {r['matched_summary_min']} | {r['wrong_summary_min']} | "
            f"{r['delta_matched_minus_wrong']} | [{r['delta_ci_lo']}, {r['delta_ci_hi']}] | "
            f"{r['ci_excludes_zero']} |"
        )
    lines += [
        "",
        "Main-panel point Δ is positive on all four pairs. Holdout point Δ is negative on all three",
        "eligible pairs (wrong ≥ matched). Whether the holdout Δ CI excludes 0 is an empirical",
        "result of this bootstrap, not a claim that the paradox is explained.",
        "",
        "## Pocket-matched Δ (Vina − best descriptor)",
        "",
        "| pair | descriptor | Vina | descriptor | Δ | 95% CI | excludes 0? |",
        "|------|------------|-----:|-----------:|--:|--------|-------------|",
    ]
    for r in desc_rows:
        lines.append(
            f"| {r['pair']} | {r['best_descriptor']} | {r['vina_summary_min']} | {r['descriptor_summary_min']} | "
            f"{r['delta_vina_minus_descriptor']} | [{r['delta_ci_lo']}, {r['delta_ci_hi']}] | "
            f"{r['ci_excludes_zero']} |"
        )
    lines += [
        "",
        "This is **not** `baseline_gate_bootstrap_v1.csv` (that file uses pooled `vina_mean`;",
        "EGFR 0.2824 ≠ Table 2 0.4297).",
        "",
        "## Frozen vs holdout",
        "",
        "See `tables/frozen_vs_holdout_v1.csv`. EGFR/HER2 holdout remains not eligible.",
        "",
        "## ML leakage check",
        "",
        "Mean (random − scaffold) on the eight directional contrasts:",
        f" {r4(float(np.mean([float(r['delta_random_minus_scaffold']) for r in ml_rows])))}.",
        "Do not hunt a split that inflates the gap.",
        "",
        "## Not done (and not invented)",
        "",
        "- EGFR unused-pool holdout",
        "- BindingDB/`as_is` docking panel",
        "- 1000 independent panels",
        "- LigPrep, PLIF, B=10000 as if it changed conclusions",
        "- 4-class accuracy as primary",
        "- receptor-only scorer as hard P0",
        "",
        f"Holdout rows written: {len(hold_rows)}.",
        "",
    ]
    (AN / "P0_MISSING_CONTENT_VERDICT_V1.md").write_text("\n".join(lines))


def main() -> None:
    errors: list[str] = []
    write_endpoint_hierarchy()
    write_pose_fairness()
    hold_rows = write_frozen_vs_holdout(errors)
    ml_rows = write_ml_merge()
    write_top10()

    main_packs = attach_descriptors({pair: assemble(pair, SPEC[pair]) for pair in PAIR_ORDER})
    hold_packs = load_holdout()

    delta_rows = []
    for pair in PAIR_ORDER:
        recs = [r for r in main_packs[pair] if r["cls"] in ("dual", "A_only", "B_only")]
        stats = boot_paired(recs, "vina_B", "vina_A", "vina_A", "vina_B")
        exp_m, exp_w = EXPECTED_MATCHED[pair], EXPECTED_WRONG[pair]
        eq4(stats["matched_D_vs_A"], exp_m[0], f"main {pair} matched DA", errors)
        eq4(stats["matched_D_vs_B"], exp_m[1], f"main {pair} matched DB", errors)
        eq4(stats["matched_summary_min"], exp_m[2], f"main {pair} matched min", errors)
        eq4(stats["wrong_D_vs_A"], exp_w[0], f"main {pair} wrong DA", errors)
        eq4(stats["wrong_D_vs_B"], exp_w[1], f"main {pair} wrong DB", errors)
        eq4(stats["wrong_summary_min"], exp_w[2], f"main {pair} wrong min", errors)
        eq4(
            stats["delta_matched_minus_wrong"],
            r4(r4(exp_m[2]) - r4(exp_w[2])),
            f"main {pair} delta checksum",
            errors,
        )
        delta_rows.append({"pair": pair, "set": "main_panel", **stats})

    for pair in HOLD_PAIRS:
        recs = [r for r in hold_packs[pair] if r["cls"] in ("dual", "A_only", "B_only")]
        stats = boot_paired(recs, "vina_B", "vina_A", "vina_A", "vina_B")
        exp_m, exp_w = EXPECTED_HOLD_MATCHED[pair], EXPECTED_HOLD_WRONG[pair]
        eq4(stats["matched_D_vs_A"], exp_m[0], f"hold {pair} matched DA", errors)
        eq4(stats["matched_D_vs_B"], exp_m[1], f"hold {pair} matched DB", errors)
        eq4(stats["matched_summary_min"], exp_m[2], f"hold {pair} matched min", errors)
        eq4(stats["wrong_D_vs_A"], exp_w[0], f"hold {pair} wrong DA", errors)
        eq4(stats["wrong_D_vs_B"], exp_w[1], f"hold {pair} wrong DB", errors)
        eq4(stats["wrong_summary_min"], exp_w[2], f"hold {pair} wrong min", errors)
        eq4(
            stats["delta_matched_minus_wrong"],
            r4(r4(exp_m[2]) - r4(exp_w[2])),
            f"hold {pair} delta checksum",
            errors,
        )
        delta_rows.append({"pair": pair, "set": "unused_pool_holdout", **stats})

    delta_fields = [
        "pair",
        "set",
        "n_dual",
        "n_A_only",
        "n_B_only",
        "n_boot_ok",
        "matched_D_vs_A",
        "matched_D_vs_B",
        "matched_summary_min",
        "matched_ci_lo",
        "matched_ci_hi",
        "wrong_D_vs_A",
        "wrong_D_vs_B",
        "wrong_summary_min",
        "wrong_ci_lo",
        "wrong_ci_hi",
        "delta_matched_minus_wrong",
        "delta_boot_mean",
        "delta_ci_lo",
        "delta_ci_hi",
        "ci_excludes_zero",
        "point_matched_gt_wrong",
    ]
    write_csv(TAB / "wrong_pocket_paired_delta_bootstrap_v1.csv", delta_rows, delta_fields)

    desc_rows = []
    for pair in PAIR_ORDER:
        arm = BEST_DESC[pair]
        recs = [
            r
            for r in main_packs[pair]
            if r["cls"] in ("dual", "A_only", "B_only") and r.get(arm) is not None
        ]
        stats = boot_paired(recs, "vina_B", "vina_A", arm, arm)
        eq4(stats["matched_summary_min"], EXPECTED_MATCHED[pair][2], f"desc {pair} vina", errors)
        eq4(stats["wrong_summary_min"], EXPECTED_DESC[pair], f"desc {pair} {arm}", errors)
        eq4(
            stats["delta_matched_minus_wrong"],
            r4(r4(EXPECTED_MATCHED[pair][2]) - r4(EXPECTED_DESC[pair])),
            f"desc {pair} delta checksum",
            errors,
        )
        desc_rows.append(
            {
                "pair": pair,
                "best_descriptor": arm,
                "n_dual": stats["n_dual"],
                "n_A_only": stats["n_A_only"],
                "n_B_only": stats["n_B_only"],
                "n_boot_ok": stats["n_boot_ok"],
                "vina_D_vs_A": stats["matched_D_vs_A"],
                "vina_D_vs_B": stats["matched_D_vs_B"],
                "vina_summary_min": stats["matched_summary_min"],
                "descriptor_D_vs_A": stats["wrong_D_vs_A"],
                "descriptor_D_vs_B": stats["wrong_D_vs_B"],
                "descriptor_summary_min": stats["wrong_summary_min"],
                "delta_vina_minus_descriptor": stats["delta_matched_minus_wrong"],
                "delta_boot_mean": stats["delta_boot_mean"],
                "delta_ci_lo": stats["delta_ci_lo"],
                "delta_ci_hi": stats["delta_ci_hi"],
                "ci_excludes_zero": stats["ci_excludes_zero"],
                "note": "pocket-matched Vina vs ligand-level descriptor; NOT pooled vina_mean",
            }
        )
    write_csv(TAB / "pocket_matched_vs_best_descriptor_delta_v1.csv", desc_rows)

    write_verdict(delta_rows, desc_rows, hold_rows, ml_rows)

    if errors:
        raise SystemExit("CHECKSUM FAILED:\n" + "\n".join(errors))
    print("wrote", TAB)
    print("checksum OK")


if __name__ == "__main__":
    main()
