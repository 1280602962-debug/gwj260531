#!/usr/bin/env python3
"""JCIM eight-row submission figures from current canonical CSVs.

Official plotter is update_figures_pr32.py. Analysis-derived values are read
from results/canonical. Frozen experimental/eligibility inputs are used only
when they cannot be recomputed upward.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jcim_figure_style import (  # noqa: E402
    C,
    COMPARABLE_THETA6_PAIRS,
    DESC_LABEL,
    FS_ANNO,
    FS_AXIS,
    HOLDOUT_PAIRS,
    UNIFIED_THRESHOLD_PAIRS,
    OUT,
    PAIR_COLOR,
    PAIR_SHORT,
    PRIMARY_PAIRS,
    ROOT,
    apply_style,
    panel_label,
    save_all,
)

DATA = ROOT / "data"
CANON = ROOT / "results" / "canonical"
PROVENANCE: dict = {"source_files": {}, "plotted": {}}

S34_CONTRAST = "D_vs_B_or_neither_pocketA"
GNINA_INDEP_PAIRS = ["EGFR/HER2", "JAK1/TYK2", "PIK3CA/mTOR"]


def _read(path: Path) -> list[dict]:
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    try:
        rel = str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        rel = str(path)
    PROVENANCE["source_files"].setdefault(rel, len(rows))
    return rows


def fnum(x) -> float:
    return float(x)


def _eq(errors: list[str], a, b, tol: float, msg: str) -> None:
    if abs(float(a) - float(b)) > tol:
        errors.append(f"{msg}: plotted {a} != source {b}")


def load() -> dict:
    canon = ROOT / "results" / "canonical"
    smin = {r["pair"]: r for r in _read(canon / "primary_summary_min.csv")}
    direc = {(r["pair"], r["estimand"]): r for r in _read(canon / "primary_directional_auroc.csv")}
    two = {r["pair"]: r for r in _read(canon / "two_pocket_mean_ranking.csv")}
    counts = {r["pair"]: r for r in _read(canon / "class_counts.csv")}
    fixed = {(r["pair"], r["contrast"]): r for r in _read(canon / "fixed_score_negative_class_delta.csv")}
    ranking = {r["pair"]: r for r in _read(canon / "top10_operating_points.csv")}
    and_filter = {r["pair"]: r for r in _read(canon / "and_filter_operating_points.csv")}
    matched = {r["pair"]: r for r in _read(canon / "matched_mismatched_pocket.csv")}
    hold = {r["pair"]: r for r in _read(canon / "holdout_metrics.csv")}
    canon_ecfp = {(r["pair"], r["contrast"]): r for r in _read(canon / "ecfp4_incremental_information.csv")}
    canon_desc = {r["pair"]: r for r in _read(canon / "descriptor_baselines.csv")}
    robust = {(r["pair"], r["engine"]): r for r in _read(canon / "computational_robustness.csv")}
    rec_sub = {r["replacement"]: r for r in _read(canon / "receptor_substitution.csv")}
    labels = {(r["pair"], r["label_rule"]): r for r in _read(canon / "label_aggregation_sensitivity.csv")}
    cluster = [r for r in _read(canon / "cluster_bootstrap_sensitivity.csv")]
    rmsd = _read(canon / "cognate_rmsd.csv")
    seeds = _read(canon / "five_seed_summary_min.csv")
    protocol = {(r["panel"], r["setting"]): r for r in _read(canon / "protocol_sensitivity.csv")}
    native = _read(canon / "external_eligibility.csv")
    census = _read(ROOT / "data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv")
    master = _read(canon / "current_score_master.csv")
    tpsa = defaultdict(list)
    for r in master:
        if r.get("pair") != "AChE/BChE":
            continue
        if r.get("analysis_set") != "main" or r.get("complete_case") not in ("1", "True"):
            continue
        if str(r.get("activity_eligible", "1")) not in ("1", "True"):
            continue
        cls = r.get("primary_class_theta6")
        if cls not in ("dual", "A_only", "B_only"):
            continue
        mol = Chem.MolFromSmiles(r.get("smiles") or "")
        if mol is None:
            continue
        tpsa[cls].append(float(Descriptors.TPSA(mol)))
    return {
        "smin": smin,
        "direc": direc,
        "two": two,
        "counts": counts,
        "fixed": fixed,
        "ranking": ranking,
        "and_filter": and_filter,
        "matched": matched,
        "hold": hold,
        "canon_ecfp": canon_ecfp,
        "canon_desc": canon_desc,
        "robust": robust,
        "rec_sub": rec_sub,
        "labels": labels,
        "cluster": cluster,
        "rmsd": rmsd,
        "seeds": seeds,
        "protocol": protocol,
        "native": native,
        "census": census,
        "tpsa": tpsa,
        "jps": rec_sub["4JPS"],
        "dxt": rec_sub["5DXT"],
        "jsx": rec_sub["4JSX"],
        "equal": fixed,
        "theta6": smin,
        "theta_all": [
            {"pair": pair, "label_rule": rule, **labels[(pair, rule)]}
            for pair in PRIMARY_PAIRS
            for rule in ("theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5")
            if (pair, rule) in labels
        ],
    }


def primary_row(D: dict, pair: str) -> dict:
    s = D["smin"][pair]
    da = D["direc"][(pair, "AUROC_D_vs_A_pocketB")]
    db = D["direc"][(pair, "AUROC_D_vs_B_pocketA")]
    two = D["two"][pair]
    n = D["counts"][pair]
    return {
        "da": fnum(da["point"]),
        "db": fnum(db["point"]),
        "smin": fnum(s["summary_min"]),
        "lo": fnum(s["ci_lo"]),
        "hi": fnum(s["ci_hi"]),
        "nei": fnum(two["two_pocket_mean_D_vs_neither"]),
        "nei_lo": fnum(two["ci_lo"]),
        "nei_hi": fnum(two["ci_hi"]),
        "n_neg": int(n["n_neither"]),
    }


def s34_row(D: dict, pair: str) -> dict:
    r = D["fixed"][(pair, S34_CONTRAST)]
    return {
        "delta": fnum(r["delta_neither_minus_selective"]),
        "lo": fnum(r["delta_ci_lo"]),
        "hi": fnum(r["delta_ci_hi"]),
        "under": str(r.get("ci_excludes_zero", "0")) in {"0", "False"},
    }


def ecfp_row(D: dict, pair: str, contrast: str) -> dict:
    r = D["canon_ecfp"][(pair, contrast)]
    return {
        "vina": fnum(r["rank_auroc_docking"]),
        "ecfp": fnum(r["cv_auroc_ECFP4"]),
        "ecfp_plus": fnum(r["cv_auroc_ECFP4_docking"]),
        "ecfp_base": fnum(r["cv_auroc_ECFP4"]),
    }


def best_desc(D: dict, pair: str) -> tuple[str, float]:
    r = D["canon_desc"][pair]
    return r["best_descriptor"], fnum(r["best_descriptor_summary_min"])


def five_seed_range(D: dict, pair: str) -> dict:
    """Range of deposited five-seed summary_min values.

    Overlay the Table 2 primary seed when same_protocol_as_primary=1.
    same_membership_as_primary is recorded separately; AChE is 0.
    """
    rows = [r for r in D["seeds"] if r["pair"] == pair]
    same_protocol = True
    same_membership = True
    status = "current_or_compatible"
    vals = []
    for r in rows:
        sm = fnum(r.get("summary_min"))
        if sm is not None:
            vals.append(sm)
        if str(r.get("same_protocol_as_primary", "1")).strip() in {"0", "False", "false"}:
            same_protocol = False
        if str(r.get("same_membership_as_primary", "1")).strip() in {"0", "False", "false"}:
            same_membership = False
        if r.get("realization_status"):
            status = r["realization_status"]
    primary = None
    if same_protocol:
        if rows:
            primary = fnum(rows[0].get("primary_summary_min"))
        if primary is None:
            primary = primary_row(D, pair)["smin"]
    return {
        "primary": primary,
        "median": float(np.median(vals)) if vals else None,
        "min": float(np.min(vals)) if vals else None,
        "max": float(np.max(vals)) if vals else None,
        "n": len(vals),
        "same_protocol_as_primary": int(same_protocol),
        "same_membership_as_primary": int(same_membership),
        "realization_status": status,
    }


def wp_main(D: dict, pair: str) -> dict:
    r = D["matched"][pair]
    lo, hi = fnum(r["delta_ci_lo"]), fnum(r["delta_ci_hi"])
    return {
        "delta": fnum(r["delta"]),
        "lo": lo,
        "hi": hi,
        "excl": not (lo <= 0 <= hi),
    }


def wp_hold(D: dict, pair: str) -> dict:
    r = D["hold"][pair]
    lo, hi = fnum(r["mm_ci_lo"]), fnum(r["mm_ci_hi"])
    return {
        "delta": fnum(r["matched_minus_mismatched"]),
        "lo": lo,
        "hi": hi,
        "excl": not (lo <= 0 <= hi),
        "smin": fnum(r["summary_min"]),
        "smin_lo": fnum(r["ci_lo"]),
        "smin_hi": fnum(r["ci_hi"]),
    }


def holdout_smin(D: dict, pair: str) -> dict:
    r = D["hold"][pair]
    return {"y": fnum(r["summary_min"]), "lo": fnum(r["ci_lo"]), "hi": fnum(r["ci_hi"])}


def theta_grid_smin(D: dict, pair: str, rule: str) -> float:
    return fnum(D["labels"][(pair, rule)]["summary_min"])


def theta_grid_record(D: dict, pair: str, rule: str) -> dict:
    r = D["labels"][(pair, rule)]
    value = fnum(r["summary_min"])
    under = min(int(r["n_dual"]), int(r["n_A_only"]), int(r["n_B_only"])) < 10
    return {"value": value, "under": under}


def gnina_indep(D: dict, pair: str) -> dict:
    r = D["robust"][(pair, "gnina_dock_mode1")]
    return {
        "smin": fnum(r["summary_min"]),
        "smin_lo": fnum(r["summary_min_ci_lo"]),
        "smin_hi": fnum(r["summary_min_ci_hi"]),
        "nei": fnum(r["auroc_D_vs_neither_mean"]),
        "nei_lo": fnum(r.get("d_vs_neither_ci_lo") or 0),
        "nei_hi": fnum(r.get("d_vs_neither_ci_hi") or 0),
    }


def forest_pairs(ax, pairs, records, title, xlabel, xlim):
    y = np.arange(len(pairs))
    for i, rec in enumerate(records):
        col = C["vina"] if rec["excl"] else "#888888"
        ax.plot([rec["lo"], rec["hi"]], [i, i], color=col, lw=1.45, zorder=3)
        ax.plot(rec["y"], i, "o", color=col, markersize=5.2, zorder=4)
    ax.axvline(0 if "Δ" in xlabel or "delta" in xlabel.lower() or xlabel.startswith("Δ") else -99,
               color=C["ink"] if "Δ" in xlabel or xlabel.startswith("Δ") else "none",
               ls="--", lw=0.85, zorder=1)
    if not (xlabel.startswith("Δ") or "Δ" in xlabel):
        ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in pairs], fontsize=6.4)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=FS_AXIS, pad=3)
    ax.invert_yaxis()
    ax.set_xlim(*xlim)


def fig1_framework(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig2_formulation(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig3_chemistry(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig4_realization(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig5_mismatched(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig6_boundary(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig_s7_diagnostics(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def fig_s8_bindingdb(D: dict) -> None:
    raise SystemExit("archived unused builder; official artwork is update_figures_pr32.py")


def verify(D: dict) -> None:
    raise SystemExit("archived unused verifier; official artwork is update_figures_pr32.py")


def write_lock_and_captions() -> None:
    """Official lock and captions live in docs/; this plotter does not rewrite them."""
    return


def main() -> None:
    print("Official artwork writer is figures/jcim_article/scripts/update_figures_pr32.py")
    print("This module is a helper library (load/helpers). It does not write plotted_values.json.")
    raise SystemExit(0)


if __name__ == "__main__":
    main()
