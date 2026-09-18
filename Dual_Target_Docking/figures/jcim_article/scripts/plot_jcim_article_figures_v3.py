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
    vals = [fnum(r["summary_min"]) for r in D["seeds"] if r["pair"] == pair]
    prim = primary_row(D, pair)["smin"]
    return {
        "primary": prim,
        "median": float(np.median(vals)),
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
        "n": len(vals),
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
    fig = plt.figure(figsize=(7.0, 5.10))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.72, 1.28], hspace=0.28, wspace=0.18)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    panel_label(ax, "A", x=-0.02, y=1.04)
    ax.text(5.0, 9.55, "Four experimental states", ha="center", fontsize=FS_AXIS, fontweight="bold")
    ax.text(4.55, 8.55, "B active", ha="center", fontsize=FS_ANNO, fontweight="bold")
    ax.text(7.55, 8.55, "B low activity", ha="center", fontsize=FS_ANNO, fontweight="bold")
    ax.text(1.65, 6.80, "A active", ha="center", va="center", fontsize=FS_ANNO, fontweight="bold")
    ax.text(1.65, 3.75, "A low activity", ha="center", va="center", fontsize=FS_ANNO, fontweight="bold")
    cells = [
        (3.05, 5.45, "dual", C["dual"]),
        (6.05, 5.45, "A-only", C["a_only"]),
        (3.05, 2.40, "B-only", C["b_only"]),
        (6.05, 2.40, "neither", "#D9D9D9"),
    ]
    for x, y0, name, col in cells:
        ax.add_patch(FancyBboxPatch(
            (x, y0), 3.0, 2.35, boxstyle="round,pad=0.04,rounding_size=0.20",
            facecolor=col, alpha=0.20 if name != "neither" else 0.55,
            edgecolor=col if name != "neither" else "#A5A5A5", lw=1.0,
        ))
        ax.text(x + 1.50, y0 + 1.35, name, ha="center", va="center", fontsize=FS_AXIS, fontweight="bold")
        if name in {"A-only", "B-only"}:
            ax.text(x + 1.50, y0 + 0.62, "selective control", ha="center", va="center",
                    fontsize=6.0, color="#555555")

    ax = fig.add_subplot(gs[0, 1])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    panel_label(ax, "B", x=-0.02, y=1.04)
    ax.text(5.0, 9.55, "Pocket-matched directional tasks", ha="center", fontsize=FS_AXIS, fontweight="bold")
    boxes = [
        (0.4, 5.55, "Dual vs A-only", "scored in pocket B"),
        (0.4, 2.35, "Dual vs B-only", "scored in pocket A"),
    ]
    for x, y, head, mid in boxes:
        ax.add_patch(FancyBboxPatch(
            (x, y), 9.2, 2.55, boxstyle="round,pad=0.05,rounding_size=0.22",
            facecolor="#F7F7F7", edgecolor="#CCCCCC", lw=0.8, clip_on=False,
        ))
        ax.text(x + 0.35, y + 1.55, head, ha="left", fontsize=FS_AXIS, fontweight="bold")
        ax.text(x + 0.35, y + 0.70, mid, ha="left", va="center", fontsize=FS_ANNO, color=C["vina"])
    ax.text(5.0, 1.15, r"summary$_{\mathrm{min}}$ = min(AUROC$_{D/A}$, AUROC$_{D/B}$)",
            ha="center", fontsize=FS_ANNO)
    ax.text(5.0, 0.40, "Descriptive worst-arm summary; both arms are reported.",
            ha="center", fontsize=6.5, color="#555555")

    gs_c = gs[1, :].subgridspec(1, 2, width_ratios=[1.05, 1.20], wspace=0.22)
    ax = fig.add_subplot(gs_c[0, 0])
    panel_label(ax, "C", x=-0.10, y=1.08)
    rows = [r for r in D["j0"] if (r.get("min_strict_hardneg") or "").strip()]
    n_pairs = len(rows)
    n_thick = sum(1 for r in rows if fnum(r["min_strict_hardneg"]) >= 50)
    stages = [
        (n_pairs, "J0 scrape"),
        (n_thick, "min hard-neg ≥50"),
        (3, "after HDAC metal exclusion"),
        (4, "historically docked"),
        (3, "after PIK3CB withdrawal"),
        (8, "primary rows after census"),
    ]
    vals = [s[0] for s in stages]
    stage_labels = [
        "J0 candidates", "hard-negative gate", "after HDAC exclusion",
        "plus supply-limited EGFR/HER2", "after PIK3CB withdrawal", "plus five census pairs",
    ]
    cols = [C["other"], C["thick"], C["metal"], C["egfr"], C["vina"], C["dual"]]
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.45, 5.55)
    ax.axis("off")
    for i, (value, label, col) in enumerate(zip(vals, stage_labels, cols)):
        yy0 = 5 - i
        ax.add_patch(FancyBboxPatch(
            (0.09, yy0 - 0.31), 0.82, 0.62, boxstyle="round,pad=0.015,rounding_size=0.035",
            facecolor="white", edgecolor=col, lw=1.0,
        ))
        ax.text(0.18, yy0, str(value), ha="center", va="center", fontsize=7.0,
                fontweight="bold", color=col)
        ax.text(0.28, yy0, label, ha="left", va="center", fontsize=6.1)
        if i < len(vals) - 1:
            ax.annotate("", xy=(0.50, yy0 - 0.50), xytext=(0.50, yy0 - 0.32),
                        arrowprops={"arrowstyle": "-|>", "color": "#999999", "lw": 0.75})
    ax.set_title("J0 scrape to eight-row primary set", fontsize=FS_AXIS, pad=4)

    bx = fig.add_subplot(gs_c[0, 1])
    j0_keys = ["HDAC1/HDAC6", "PIK3CA/MTOR", "ACHE/BCHE", "PIK3CA/PIK3CB", "EGFR/HER2"]
    j0_lab = ["HDAC1/\nHDAC6", "PIK3CA/\nmTOR", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "EGFR/\nHER2"]
    j0_map = {r["pair"]: fnum(r["min_strict_hardneg"]) for r in rows}
    j0_vals = [j0_map[k] for k in j0_keys]
    j0_cols = [C["metal"], C["thick"], C["thick"], C["thick"], C["egfr"]]
    dump_vals = [fnum(D["five_xdb"][(p, "ChEMBL37_dump", "pChEMBL_STANDARD_OK")]["min_strict_hardneg"])
                 for p in COMPARABLE_THETA6_PAIRS]
    dump_lab = ["F2/\nF10", "JAK1/\nTYK2", "JAK1/\nJAK2", "PPARG/\nPPARA", "PPARA/\nPPARD"]
    labels = [s.replace("\n", "/") for s in j0_lab + dump_lab]
    vals = j0_vals + dump_vals
    yy = np.arange(len(vals))
    bx.barh(yy, vals, color=j0_cols + [C["vina"]] * 5, height=0.66, zorder=3)
    bx.axvline(50, color=C["ink"], ls="--", lw=0.8, zorder=2)
    bx.set_yticks(yy)
    bx.set_yticklabels(labels, fontsize=6.0)
    bx.invert_yaxis()
    bx.set_xlabel("Minimum strict hard-negative count")
    bx.set_xlim(0, 140)
    bx.set_title("Hard-negative supply (two audits)", fontsize=FS_AXIS, pad=4)
    bx.axhline(4.5, color="#DDDDDD", lw=0.8, zorder=1)
    bx.text(50, -0.43, "gate ≥50", ha="center", va="center", fontsize=6.0,
            color="#555555", zorder=4,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})

    fracs = [fnum(D["overlap"][p]["fraction_union_measured_both"]) for p in
             ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]]
    PROVENANCE["plotted"]["fig1C"] = {
        "n_pairs": n_pairs,
        "n_thick": n_thick,
        "stages": {lab: v for v, lab in stages},
        "j0_hardneg": {k: j0_map[k] for k in j0_keys},
        "dump_hardneg": {p: dump_vals[i] for i, p in enumerate(COMPARABLE_THETA6_PAIRS)},
        "complete_case_min": min(fracs),
        "complete_case_max": max(fracs),
    }
    fig.subplots_adjust(left=0.13, right=0.98, top=0.94, bottom=0.10)
    save_all(fig, "Fig1_four_state_and_supply")
    plt.close(fig)


def fig2_formulation(D: dict) -> None:
    fig = plt.figure(figsize=(7.0, 5.95))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.00, 1.12, 1.08], hspace=0.42)
    y = np.arange(len(PRIMARY_PAIRS))
    off = 0.12

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A", x=-0.16, y=1.04)
    da = [primary_row(D, p)["da"] for p in PRIMARY_PAIRS]
    db = [primary_row(D, p)["db"] for p in PRIMARY_PAIRS]
    for i, (a, b) in enumerate(zip(da, db)):
        ax.plot([a, b], [i, i], color="#B8B8B8", lw=0.85, zorder=2)
        ax.plot(a, i, "o", color=C["vina"], markersize=4.8, zorder=4)
        ax.plot(b, i, "s", color=C["a_only"], markersize=4.5, zorder=4)
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.25, 0.85)
    ax.set_title("Pocket-matched directional AUROC", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5,
               label="Dual vs A-only, pocket B"),
        Line2D([0], [0], marker="s", color=C["a_only"], ls="none", ms=4.7,
               label="Dual vs B-only, pocket A"),
    ], loc="lower right", fontsize=6.0, frameon=False)
    PROVENANCE["plotted"]["fig2A"] = {"DA": da, "DB": db, "pairs": list(PRIMARY_PAIRS)}

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B", x=-0.16, y=1.04)
    plotted = {}
    for i, p in enumerate(PRIMARY_PAIRS):
        r = primary_row(D, p)
        plotted[p] = r
        dir_err = np.array([[r["smin"] - r["lo"]], [r["hi"] - r["smin"]]])
        nei_err = np.array([[r["nei"] - r["nei_lo"]], [r["nei_hi"] - r["nei"]]])
        ax.errorbar(r["smin"], i + off, xerr=dir_err, fmt="o", color=C["vina"],
                    ecolor=C["vina"], capsize=1.8, elinewidth=0.9, markersize=4.8, zorder=4)
        marker = "D" if p == "PIK3CA/mTOR" else "s"
        ax.errorbar(r["nei"], i - off, xerr=nei_err, fmt=marker, color=C["desc"],
                    ecolor=C["desc"], capsize=1.8, elinewidth=0.9, markersize=4.5, zorder=4)
        if p == "PIK3CA/mTOR":
            ax.text(r["nei"] + 0.025, i - off, f"n={r['n_neg']}", va="center", fontsize=6.0,
                    color=C["desc"], fontweight="bold")
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.15, 1.02)
    ax.set_title("Directional summary and Dual vs neither", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5,
               label=r"directional summary$_{\mathrm{min}}$"),
        Line2D([0], [0], marker="s", color=C["desc"], ls="none", ms=4.8,
               label=r"Dual vs neither, Vina$_{\mathrm{mean}}$"),
        Line2D([0], [0], marker="D", color=C["desc"], ls="none", ms=4.8,
               label=f"neither n={primary_row(D, 'PIK3CA/mTOR')['n_neg']}"),
    ], loc="lower right", fontsize=5.9, frameon=False)
    PROVENANCE["plotted"]["fig2B"] = {
        p: {"directional": {"y": plotted[p]["smin"], "lo": plotted[p]["lo"], "hi": plotted[p]["hi"]},
            "neither": {"y": plotted[p]["nei"], "lo": plotted[p]["nei_lo"], "hi": plotted[p]["nei_hi"],
                        "n_neg": plotted[p]["n_neg"]}}
        for p in PRIMARY_PAIRS
    }

    ax = fig.add_subplot(gs[2, 0])
    panel_label(ax, "C", x=-0.16, y=1.04)
    recs = {}
    for i, p in enumerate(PRIMARY_PAIRS):
        r = s34_row(D, p)
        recs[p] = r
        col = C["egfr"] if p in ("EGFR/HER2", "JAK1/TYK2") else C["vina"]
        ax.plot([r["lo"], r["hi"]], [i, i], color=col, lw=1.45, zorder=3)
        m = "D" if r["under"] else "o"
        ax.plot(r["delta"], i, m, color=col, markersize=6.0 if p in ("EGFR/HER2", "JAK1/TYK2") else 5.2, zorder=4)
        if p in ("EGFR/HER2", "JAK1/TYK2"):
            ax.text(r["hi"] + 0.015, i, f"{r['delta']:.3f}", va="center", fontsize=6.0, color=col)
    ax.axvline(0, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("ΔAUROC (neither − B-only), pocket A score fixed")
    ax.set_xlim(-0.62, 0.78)
    ax.set_title("Fixed score; swap negative class", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["egfr"], ls="none", ms=6, label="EGFR/HER2 and JAK1/TYK2"),
        Line2D([0], [0], marker="D", color=C["vina"], ls="none", ms=6, label="underpowered neither"),
    ], loc="lower right", fontsize=6.0, frameon=False)
    PROVENANCE["plotted"]["fig2C"] = recs

    fig.subplots_adjust(left=0.18, right=0.97, top=0.96, bottom=0.075)
    save_all(fig, "Fig2_negative_class_formulation")
    plt.close(fig)


def fig3_chemistry(D: dict) -> None:
    fig = plt.figure(figsize=(7.0, 5.95))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 1.00], hspace=0.36, wspace=0.32)

    ax = fig.add_subplot(gs[0, :])
    panel_label(ax, "A", x=-0.08, y=1.14)
    y = np.arange(len(PRIMARY_PAIRS))
    off = {"vina_da": 0.30, "vina_db": 0.10, "ecfp_da": -0.10, "ecfp_db": -0.30}
    cols = {"vina_da": C["vina"], "vina_db": C["gnina"], "ecfp_da": C["desc"], "ecfp_db": C["a_only"]}
    marks = {"vina_da": "o", "vina_db": "s", "ecfp_da": "^", "ecfp_db": "D"}
    plotted = {k: [] for k in off}
    for i, p in enumerate(PRIMARY_PAIRS):
        a = ecfp_row(D, p, "D_vs_A")
        b = ecfp_row(D, p, "D_vs_B")
        vals = {"vina_da": a["vina"], "vina_db": b["vina"], "ecfp_da": a["ecfp"], "ecfp_db": b["ecfp"]}
        for k, v in vals.items():
            plotted[k].append(v)
            ax.plot(v, i + off[k], marks[k], color=cols[k], markersize=5.0, zorder=4)
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_ylim(7.60, -2.15)  # headroom so the legend sits above EGFR/HER2
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.20, 1.05)
    ax.set_title("Vina rank AUROC versus ECFP4 scaffold GroupKFold", fontsize=FS_AXIS, pad=8)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5.5, label="Vina D/A"),
        Line2D([0], [0], marker="s", color=C["gnina"], ls="none", ms=5.0, label="Vina D/B"),
        Line2D([0], [0], marker="^", color=C["desc"], ls="none", ms=5.5, label="ECFP4 D/A"),
        Line2D([0], [0], marker="D", color=C["a_only"], ls="none", ms=5.0, label="ECFP4 D/B"),
    ], loc="upper center", ncol=4, fontsize=6.0, frameon=True,
       columnspacing=0.9, handletextpad=0.35, borderaxespad=0.1,
       fancybox=False, edgecolor="none", facecolor="white", framealpha=0.92)
    PROVENANCE["plotted"]["fig3A"] = plotted

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B", x=-0.22, y=1.06)
    deltas, labels, arm_tags = [], [], []
    for p in PRIMARY_PAIRS:
        for contrast, tag in (("D_vs_A", "D/A"), ("D_vs_B", "D/B")):
            r = ecfp_row(D, p, contrast)
            dlt = r["ecfp_plus"] - r["ecfp_base"]
            deltas.append(dlt)
            labels.append(PAIR_SHORT[p])
            arm_tags.append(tag)
    yy = np.arange(len(PRIMARY_PAIRS))
    for i, p in enumerate(PRIMARY_PAIRS):
        d_a, d_b = deltas[2 * i], deltas[2 * i + 1]
        ax.plot([0, d_a], [i + 0.10, i + 0.10], color=C["vina"], lw=0.9, zorder=2)
        ax.plot([0, d_b], [i - 0.10, i - 0.10], color=C["a_only"], lw=0.9, zorder=2)
        ax.plot(d_a, i + 0.10, "o", color=C["vina"], markersize=4.2, zorder=4)
        ax.plot(d_b, i - 0.10, "s", color=C["a_only"], markersize=4.0, zorder=4)
    ax.axvline(0, color=C["ink"], lw=0.8, zorder=2)
    ax.set_yticks(yy)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("ΔAUROC (ECFP4+Vina − ECFP4)")
    ax.set_xlim(-0.028, 0.028)
    ax.set_title("Increment from adding Vina", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=4.5, label="D/A"),
        Line2D([0], [0], marker="s", color=C["a_only"], ls="none", ms=4.3, label="D/B"),
    ], loc="lower right", bbox_to_anchor=(0.98, 0.06), fontsize=5.9, frameon=True, fancybox=False, edgecolor="none",
       facecolor="white", framealpha=0.92, borderpad=0.35)
    PROVENANCE["plotted"]["fig3B_deltas"] = deltas
    PROVENANCE["plotted"]["fig3B_max_abs"] = float(max(abs(d) for d in deltas))

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "C", x=-0.22, y=1.06)
    rng = np.random.default_rng(20260729)
    data = [D["tpsa"]["dual"], D["tpsa"]["A_only"], D["tpsa"]["B_only"]]
    colors = [C["dual"], C["a_only"], C["b_only"]]
    ns = [len(d) for d in data]
    for i, (vals, col) in enumerate(zip(data, colors), start=1):
        arr = np.asarray(vals, float)
        jitter = rng.uniform(-0.12, 0.12, size=len(arr))
        ax.scatter(np.full(len(arr), i) + jitter, arr, s=10, color=col, alpha=0.75,
                   edgecolors=C["ink"], linewidths=0.25, zorder=3)
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        ax.plot([i - 0.22, i + 0.22], [med, med], color=C["ink"], lw=1.3, zorder=4)
        ax.plot([i, i], [q1, q3], color=C["ink"], lw=1.0, zorder=4)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([f"dual\nn={ns[0]}", f"A-only\nn={ns[1]}", f"B-only\nn={ns[2]}"], fontsize=6.2)
    ax.set_ylabel(r"TPSA ($\mathrm{\AA}^2$)")
    ax.set_title("AChE/BChE TPSA by activity class", fontsize=FS_AXIS, pad=3)
    PROVENANCE["plotted"]["fig3C"] = {
        "n": ns,
        "mean": [float(np.mean(d)) for d in data],
        "median": [float(np.median(d)) for d in data],
    }

    fig.subplots_adjust(left=0.16, right=0.98, top=0.94, bottom=0.09)
    save_all(fig, "Fig3_ligand_chemistry")
    plt.close(fig)


def fig4_realization(D: dict) -> None:
    fig = plt.figure(figsize=(7.0, 3.55))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.28, 0.94, 0.28, 1.16], wspace=0.18)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 3])]

    ax = axes[0]
    panel_label(ax, "A", x=-0.20, y=1.06)
    vina_smin, vina_nei, g_smin, g_nei = [], [], [], []
    for p in GNINA_INDEP_PAIRS:
        pr = primary_row(D, p)
        g = gnina_indep(D, p)
        vina_smin.append(pr["smin"])
        vina_nei.append(pr["nei"])
        g_smin.append(g["smin"])
        g_nei.append(g["nei"])
    yy = np.arange(len(GNINA_INDEP_PAIRS))
    for i in range(len(yy)):
        ax.plot([vina_smin[i], g_smin[i]], [i - 0.12, i - 0.12], color="#B8B8B8", lw=0.9, zorder=2)
        ax.plot([vina_nei[i], g_nei[i]], [i + 0.12, i + 0.12], color="#B8B8B8", lw=0.9, zorder=2)
    ax.scatter(vina_smin, yy - 0.12, s=27, marker="o", color=C["vina"], zorder=4, label="Vina directional")
    ax.scatter(g_smin, yy - 0.12, s=27, marker="s", color=C["vina"], zorder=4, label="GNINA directional")
    ax.scatter(vina_nei, yy + 0.12, s=27, marker="o", facecolors="white", edgecolors=C["desc"], linewidths=1.0,
               zorder=4, label="Vina Dual vs neither")
    ax.scatter(g_nei, yy + 0.12, s=27, marker="s", facecolors="white", edgecolors=C["desc"], linewidths=1.0,
               zorder=4, label="GNINA Dual vs neither")
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(yy)
    ax.set_yticklabels([PAIR_SHORT[p] for p in GNINA_INDEP_PAIRS], fontsize=6.2)
    ax.invert_yaxis()
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.12, 0.92)
    ax.set_title("Independent pose generation", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=5.7, frameon=False,
              columnspacing=0.7, handletextpad=0.35)
    PROVENANCE["plotted"]["fig4A"] = {
        "pairs": list(GNINA_INDEP_PAIRS),
        "vina_smin": vina_smin, "gnina_smin": g_smin,
        "vina_neither": vina_nei, "gnina_neither": g_nei,
    }

    ax = axes[1]
    panel_label(ax, "B", x=-0.20, y=1.06)
    crystals = ["4L23", "4JPS", "5DXT", "4JSX"]
    pm = [
        (primary_row(D, "PIK3CA/mTOR")["smin"], primary_row(D, "PIK3CA/mTOR")["lo"], primary_row(D, "PIK3CA/mTOR")["hi"]),
        (fnum(D["jps"]["summary_min"]), fnum(D["jps"]["summary_min_ci_lo"]), fnum(D["jps"]["summary_min_ci_hi"])),
        (fnum(D["dxt"]["summary_min"]), fnum(D["dxt"]["summary_min_ci_lo"]), fnum(D["dxt"]["summary_min_ci_hi"])),
        (fnum(D["jsx"]["summary_min"]), fnum(D["jsx"]["summary_min_ci_lo"]), fnum(D["jsx"]["summary_min_ci_hi"])),
    ]
    fig4b = []
    cols = [C["vina"], C["holdout"], C["a_only"], C["gnina"]]
    for i, ((y, lo, hi), col, cr) in enumerate(zip(pm, cols, crystals)):
        ax.errorbar(i, y, yerr=[[y - lo], [hi - y]], fmt="o", color=col, ecolor=col,
                    elinewidth=1.2, capsize=2.0, markersize=6.0, zorder=4)
        fig4b.append({"pair": "PIK3CA/mTOR", "crystal": cr, "y": y, "lo": lo, "hi": hi})
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["4L23\nPIK3CA", "4JPS\nPIK3CA", "5DXT\nPIK3CA", "4JSX\nmTOR"], fontsize=6.0)
    ax.set_ylabel(r"summary$_{\mathrm{min}}$")
    ax.set_ylim(0.12, 1.02)
    ax.set_xlim(-0.55, 3.45)
    ax.set_title("PIK3CA/mTOR receptor structures", fontsize=FS_AXIS, pad=3)
    PROVENANCE["plotted"]["fig4B"] = fig4b

    ax = axes[2]
    panel_label(ax, "C", x=-0.22, y=1.06)
    plotted_s = {}
    yy = np.arange(len(PRIMARY_PAIRS))
    for i, p in enumerate(PRIMARY_PAIRS):
        r = five_seed_range(D, p)
        ax.plot([r["min"], r["max"]], [i, i], color=C["vina"], lw=1.4, zorder=3)
        ax.plot(r["median"], i, "o", color=C["vina"], markersize=5.2, zorder=4)
        ax.plot(r["primary"], i, "D", color=C["desc"], markersize=4.4, zorder=5)
        plotted_s[p] = r
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(yy)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel(r"summary$_{\mathrm{min}}$ across five Vina seeds")
    ax.set_xlim(0.22, 0.82)
    ax.set_title("Five Vina seeds", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="D", color=C["desc"], ls="none", ms=5, label="primary seed"),
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5, label="median"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=6.0, frameon=False)
    PROVENANCE["plotted"]["fig4C"] = plotted_s

    fig.subplots_adjust(left=0.10, right=0.995, top=0.88, bottom=0.23)
    save_all(fig, "Fig5_computational_realization")
    plt.close(fig)


def fig5_mismatched(D: dict) -> None:
    fig, axes = plt.subplots(
        1, 3, figsize=(7.0, 4.55), sharey=True,
        gridspec_kw={"width_ratios": [1.00, 1.00, 1.18]},
    )

    ax = axes[0]
    panel_label(ax, "A", x=-0.28, y=1.04)
    recs = []
    for p in PRIMARY_PAIRS:
        r = wp_main(D, p)
        recs.append({"pair": p, "y": r["delta"], "lo": r["lo"], "hi": r["hi"], "excl": r["excl"]})
    forest_pairs(ax, PRIMARY_PAIRS, recs, "Main panels",
                 r"Δsummary$_{\mathrm{min}}$ (matched − mismatched)", (-0.22, 0.36))
    PROVENANCE["plotted"]["fig5A"] = recs

    ax = axes[1]
    panel_label(ax, "B", x=-0.28, y=1.04)
    recs_b = []
    for p in PRIMARY_PAIRS:
        if p == "EGFR/HER2":
            recs_b.append({"pair": p, "missing": True})
            continue
        r = wp_hold(D, p)
        recs_b.append({"pair": p, "y": r["delta"], "lo": r["lo"], "hi": r["hi"], "excl": r["excl"]})
    y_all = np.arange(len(PRIMARY_PAIRS))
    for i, rec in enumerate(recs_b):
        if rec.get("missing"):
            ax.text(-0.33, i, "n/a", color="#777777", fontsize=6.0, va="center")
            continue
        col = C["vina"] if rec["excl"] else "#888888"
        ax.plot([rec["lo"], rec["hi"]], [i, i], color=col, lw=1.45, zorder=3)
        ax.plot(rec["y"], i, "o", color=col, markersize=5.2, zorder=4)
    ax.axvline(0, color=C["ink"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y_all)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.0)
    ax.set_xlim(-0.36, 0.38)
    ax.set_title("Unused-pool holdout", fontsize=FS_AXIS, pad=3)
    ax.set_xlabel(r"Δsummary$_{\mathrm{min}}$ (matched − mismatched)")
    PROVENANCE["plotted"]["fig5B"] = recs_b

    ax = axes[2]
    panel_label(ax, "C", x=-0.28, y=1.04)
    recs_c = []
    for i, p in enumerate(PRIMARY_PAIRS):
        if p == "EGFR/HER2":
            recs_c.append({"pair": p, "missing": True})
            ax.text(0.17, i, "n/a", color="#777777", fontsize=6.0, va="center")
            continue
        main = primary_row(D, p)
        h = holdout_smin(D, p)
        recs_c.append({"pair": p, "main": main["smin"], "main_lo": main["lo"], "main_hi": main["hi"],
                       "hold": h["y"], "hold_lo": h["lo"], "hold_hi": h["hi"]})
        ax.errorbar(main["smin"], i + 0.14, xerr=[[main["smin"] - main["lo"]], [main["hi"] - main["smin"]]],
                    fmt="o", color=C["vina"], ecolor=C["vina"], elinewidth=1.1, capsize=1.8, markersize=5.0, zorder=4)
        ax.errorbar(h["y"], i - 0.14, xerr=[[h["y"] - h["lo"]], [h["hi"] - h["y"]]],
                    fmt="s", color=C["holdout"], ecolor=C["holdout"], elinewidth=1.1, capsize=1.8, markersize=4.8, zorder=4)
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(np.arange(len(PRIMARY_PAIRS)))
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.0)
    ax.set_xlabel(r"summary$_{\mathrm{min}}$")
    ax.set_xlim(0.12, 1.02)
    ax.set_title("Holdout vs main panel", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5.5, label="main"),
        Line2D([0], [0], marker="s", color=C["holdout"], ls="none", ms=5.0, label="holdout"),
    ], loc="upper right", fontsize=6.0, frameon=False)
    PROVENANCE["plotted"]["fig5C"] = recs_c

    for ax in axes[1:]:
        ax.tick_params(labelleft=True)
    fig.subplots_adjust(wspace=0.48, left=0.15, right=0.98, top=0.89, bottom=0.13)
    save_all(fig, "Fig4_mismatched_pocket")
    plt.close(fig)


def fig6_boundary(D: dict) -> None:
    from plot_jcim_si_composites_v1 import _pm48_e8

    fig = plt.figure(figsize=(7.0, 5.25))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.20, 1.20, 1.00], hspace=0.72, wspace=0.55)
    ax_a = fig.add_subplot(gs[:, :2])
    ax_b = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[1, 2])
    ax_d = fig.add_subplot(gs[2, 2])
    rules = ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"]
    rule_lab = ["θ=5.5", "θ=6.0", "θ=6.5", "strict"]

    ax = ax_a
    panel_label(ax, "A", x=-0.14, y=1.03)
    grid = {}
    matrix = []
    under = []
    for p in PRIMARY_PAIRS:
        recs = [theta_grid_record(D, p, rule) for rule in rules]
        vals = [r["value"] for r in recs]
        matrix.append(vals)
        under.append([r["under"] for r in recs])
        grid[p] = vals
    im = ax.imshow(matrix, cmap="RdBu", vmin=0.25, vmax=0.75, aspect="auto", interpolation="nearest")
    for i, vals in enumerate(matrix):
        for j, value in enumerate(vals):
            suffix = "†" if under[i][j] else ""
            ax.text(j, i, f"{value:.3f}{suffix}", ha="center", va="center", fontsize=5.9,
                    color="white" if value < 0.35 or value > 0.67 else C["ink"])
    ax.set_xticks(np.arange(len(rules)))
    ax.set_xticklabels(rule_lab, fontsize=6.4)
    ax.set_yticks(np.arange(len(PRIMARY_PAIRS)))
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.4)
    ax.set_title(r"Label rule and summary$_{\mathrm{min}}$", fontsize=FS_AXIS, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.045, pad=0.08)
    cbar.set_label(r"summary$_{\mathrm{min}}$", fontsize=FS_ANNO)
    cbar.ax.tick_params(labelsize=5.8, length=2)
    PROVENANCE["plotted"]["fig6A"] = grid

    ax = ax_b
    panel_label(ax, "B", x=-0.30, y=1.03)
    y48 = fnum(D["pm110"][("PM48", "vina")]["summary_min"])
    y110 = fnum(D["pm110"][("PM110", "vina")]["summary_min"])
    ax.plot([y48, y110], [0, 0], color="#B8B8B8", lw=1.0, zorder=2)
    ax.plot(y48, 0, "o", color=C["vina"], markersize=5.2, zorder=4, label="PM48")
    ax.plot(y110, 0, "s", color=C["holdout"], markersize=4.8, zorder=4, label="PM110")
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks([])
    ax.set_xlabel(r"Vina summary$_{\mathrm{min}}$")
    ax.set_xlim(0.45, 0.78)
    ax.set_title("PIK3CA/mTOR panel size", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", ncol=2, fontsize=5.7, frameon=False)
    PROVENANCE["plotted"]["fig6B"] = {"PM48": y48, "PM110": y110}

    ax = ax_c
    panel_label(ax, "C", x=-0.30, y=1.03)
    e8 = _pm48_e8()
    e16 = primary_row(D, "PIK3CA/mTOR")["smin"]
    ax.plot([e16, e8["summary_min"]], [0, 0], color="#B8B8B8", lw=1.0, zorder=2)
    ax.plot(e16, 0, "o", color=C["vina"], markersize=5.2, zorder=4, label="E=16")
    ax.plot(e8["summary_min"], 0, "s", color=C["gnina"], markersize=4.8, zorder=4, label="E=8")
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks([])
    ax.set_xlabel(r"Vina summary$_{\mathrm{min}}$")
    ax.set_xlim(0.45, 0.78)
    ax.set_title("Search exhaustiveness", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", ncol=2, fontsize=5.7, frameon=False)
    PROVENANCE["plotted"]["fig6C"] = {"E16": e16, "E8": e8["summary_min"], "e8_n": (e8["nD"], e8["nA"], e8["nB"])}

    ax = ax_d
    panel_label(ax, "D", x=-0.30, y=1.03)
    contract = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
    four = [next(r for r in D["native"] if r["pair"] == p) for p in contract]
    n_fail = sum(r["packaged_as_external_evaluation"] == "0" for r in four)
    n_pass = 4 - n_fail
    y_contract = np.arange(len(contract))
    ax.scatter(np.zeros(len(contract)), y_contract, marker="x", s=28, color=C["egfr"], linewidths=1.2, zorder=4)
    ax.set_yticks(y_contract)
    ax.set_yticklabels(["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB†", "PIK3CA/mTOR"], fontsize=5.7)
    ax.invert_yaxis()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["fail", "pass"], fontsize=6.0)
    ax.set_xlim(-0.25, 1.15)
    ax.set_title("BindingDB external gate", fontsize=FS_AXIS, pad=3)
    ax.text(0.0, -0.20, "† historical contract", transform=ax.transAxes,
            ha="left", fontsize=5.6, color="#555555")
    PROVENANCE["plotted"]["fig6D"] = {"n_fail": n_fail, "n_pass": n_pass, "contract": contract}

    fig.subplots_adjust(left=0.15, right=0.98, top=0.94, bottom=0.12)
    save_all(fig, "Fig6_evidence_boundary")
    plt.close(fig)


def toc_graphic() -> None:
    fig = plt.figure(figsize=(3.25, 1.75), dpi=300)
    ax = fig.add_axes([0.03, 0.08, 0.94, 0.84])
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 16)
    ax.axis("off")
    xs = [2.0, 6.2, 10.4, 14.6]
    names = ["Dual", "A-only", "B-only", "Neither"]
    cols = [(C["dual"], C["dual"]), (C["a_only"], "#DDDDDD"), ("#DDDDDD", C["b_only"]), ("#DDDDDD", "#DDDDDD")]
    for x, name, (c1, c2) in zip(xs, names, cols):
        ax.add_patch(Circle((x - 0.38, 10.6), 0.48, facecolor=c1, edgecolor=C["ink"], lw=0.45))
        ax.add_patch(Circle((x + 0.38, 10.6), 0.48, facecolor=c2, edgecolor=C["ink"], lw=0.45))
        ax.text(x, 8.7, name, ha="center", fontsize=6.5)
    ax.add_patch(FancyBboxPatch(
        (17.4, 8.55), 11.6, 3.9, boxstyle="round,pad=0.08,rounding_size=0.20",
        facecolor="#F4F7FA", edgecolor=C["vina"], lw=0.9,
    ))
    ax.text(23.2, 11.55, "pocket-matched", ha="center", fontsize=7, fontweight="bold")
    ax.text(23.2, 9.85, "directional evaluation", ha="center", fontsize=7)
    ax.text(15.0, 4.6, "Dual-vs-Neither  ≠  Dual-vs-selective", ha="center", fontsize=7.2, fontweight="bold")
    ax.text(15.0, 2.4, "Negative-class definition changes apparent performance",
            ha="center", fontsize=6.3, color="#555555")
    save_all(fig, "TOC_graphic", toc=True)
    plt.close(fig)


def fig_s4_forest(D: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 5.40))
    y = np.arange(len(PRIMARY_PAIRS))
    plotted = {}
    for i, p in enumerate(PRIMARY_PAIRS):
        r = primary_row(D, p)
        name, dval = best_desc(D, p)
        plotted[p] = {"vina": r, "desc_name": name, "desc": dval}
        if i % 2 == 0:
            ax.axhspan(i - 0.42, i + 0.42, color="#F4F7FA", zorder=0)
        ax.plot([r["lo"], r["hi"]], [i + 0.12, i + 0.12], color=C["vina"], lw=1.5, zorder=3)
        ax.plot(r["smin"], i + 0.12, "o", color=C["vina"], markersize=5.6, zorder=4)
        ax.plot(dval, i - 0.14, "s", color=C["desc"], markersize=4.8, zorder=4)
        ax.text(0.98, i, DESC_LABEL.get(name, name), transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=6.0, color=C["desc"])
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=7.0)
    ax.invert_yaxis()
    ax.set_xlabel(r"summary$_{\mathrm{min}}$ AUROC")
    ax.set_xlim(0.12, 1.02)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="Vina"),
        Line2D([0], [0], marker="s", color=C["desc"], ls="none", ms=5, label="best single descriptor"),
    ], loc="lower right", fontsize=6.4, frameon=False)
    PROVENANCE["plotted"]["figS4"] = plotted
    fig.subplots_adjust(left=0.16, right=0.86, top=0.94, bottom=0.12)
    save_all(fig, "FigS2_pocket_matched_forest")
    plt.close(fig)


def fig_s5_holdout(D: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.20))
    x = np.arange(len(HOLDOUT_PAIRS))
    fig_s = []
    for i, p in enumerate(HOLDOUT_PAIRS):
        m = primary_row(D, p)
        h = holdout_smin(D, p)
        fig_s.append({"pair": p, "main": m["smin"], "main_lo": m["lo"], "main_hi": m["hi"],
                      "hold": h["y"], "hold_lo": h["lo"], "hold_hi": h["hi"]})
        ax.errorbar(i - 0.12, m["smin"], yerr=[[m["smin"] - m["lo"]], [m["hi"] - m["smin"]]],
                    fmt="o", color=C["vina"], ecolor=C["vina"], elinewidth=1.2, capsize=2.0, markersize=5.8, zorder=4)
        ax.errorbar(i + 0.12, h["y"], yerr=[[h["y"] - h["lo"]], [h["hi"] - h["y"]]],
                    fmt="s", color=C["holdout"], ecolor=C["holdout"], elinewidth=1.2, capsize=2.0, markersize=5.4, zorder=4)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([PAIR_SHORT[p].replace("/", "/\n") for p in HOLDOUT_PAIRS], fontsize=6.0)
    ax.set_ylabel(r"Pocket-matched summary$_{\mathrm{min}}$")
    ax.set_ylim(0.10, 1.02)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="Main panel"),
        Line2D([0], [0], marker="s", color=C["holdout"], ls="none", ms=5.5, label="Unused-pool holdout"),
    ], loc="upper left", fontsize=6.4, frameon=False)
    ax.text(0.98, 0.03, "EGFR/HER2 has no holdout", transform=ax.transAxes, ha="right", fontsize=6.0, color="#666666")
    PROVENANCE["plotted"]["figS5"] = fig_s
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.16)
    save_all(fig, "FigS5_unused_pool_holdout")
    plt.close(fig)


def fig_s7_diagnostics(D: dict) -> None:
    census = _read(DATA / "jcim_novelty_v0/tables/theta6_pair_census_v1.csv")
    and_rows = _read(DATA / "jcim_novelty_v0/tables/and_filter_operating_point_v1.csv")
    ligand_rows = _read(DATA / "jcim_novelty_v0/tables/ligand_only_fullmap_auroc_v1.csv")
    n_pairs = len(census)
    n_dir = sum(int(r["directional_n10"]) for r in census)
    n_form = sum(int(r["formulation_n10"]) for r in census)
    n_dock_j0 = sum(int(r["docked_in_this_paper"]) for r in census)

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.45))
    ax = axes[0]
    panel_label(ax, "A", x=-0.18, y=1.06)
    counts = [n_pairs, n_dir, n_form, n_dock_j0, 8]
    labels = ["J0 pairs", "directional ≥10", "formulation ≥10", "docked in J0", "primary now"]
    yy = np.arange(5)
    ax.barh(yy, counts, color=[C["vina"], C["desc"], C["thick"], C["egfr"], C["dual"]], height=0.64, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels, fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("Pair count")
    ax.set_xlim(0, max(counts) + 8)
    ax.set_title("θ=6.0 census and primary set", fontsize=FS_AXIS, pad=3)
    for i, v in enumerate(counts):
        ax.text(v + 0.8, i, str(v), va="center", fontsize=6.5)

    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.06)
    for pair in UNIFIED_THRESHOLD_PAIRS:
        sub = [r for r in and_rows if r["pair"] == pair and r["score"] == "vina_worst"]
        rec = [fnum(r["recall_dual"]) for r in sub]
        prec = [fnum(r["precision_dual"]) for r in sub]
        ax.plot(rec, prec, marker="o", color=PAIR_COLOR[pair], label=pair, lw=1.1, markersize=4.0, zorder=3)
    ax.set_xlabel("Dual recall")
    ax.set_ylabel("Dual precision")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_title("AND-like filter (vina_worst)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=6.0, frameon=False)

    ax = axes[2]
    panel_label(ax, "C", x=-0.18, y=1.06)
    neither, directional = [], []
    for pair in UNIFIED_THRESHOLD_PAIRS:
        n = next(r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "D_vs_neither")
        s = next(r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "summary_min_ecfp4")
        neither.append(fnum(n["ecfp4_groupkfold_auroc"]))
        directional.append(fnum(s["ecfp4_groupkfold_auroc"]))
    x = np.arange(len(UNIFIED_THRESHOLD_PAIRS))
    ax.bar(x - 0.18, neither, 0.36, color=C["vina"], label="Dual vs neither", zorder=3)
    ax.bar(x + 0.18, directional, 0.36, color=C["egfr"], label="ECFP4 summary_min", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nmTOR"], fontsize=6.0)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_ylim(0.35, 1.05)
    ax.set_ylabel("GroupKFold AUROC")
    ax.set_title("Ligand-only full maps (EGFR/HER2, AChE/BChE, PIK3CA/mTOR)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=6.0, frameon=False)

    PROVENANCE["plotted"]["figS7"] = {
        "n_pairs": n_pairs, "n_dir": n_dir, "n_form": n_form, "n_dock_j0": n_dock_j0,
        "n_primary_now": 8, "neither": neither, "directional": directional,
    }
    fig.subplots_adjust(wspace=0.46, left=0.13, right=0.98, top=0.86, bottom=0.30)
    save_all(fig, "FigS1_posthoc_diagnostics")
    plt.close(fig)


def fig_s8_bindingdb(D: dict) -> None:
    contract = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
    four = [next(r for r in D["native"] if r["pair"] == p) for p in contract]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.40), gridspec_kw={"width_ratios": [1.2, 1.05]})
    ax = axes[0]
    panel_label(ax, "A", x=-0.16, y=1.06)
    stages = ["native_paired", "after_literature", "after_structure", "after_ecfp_lt_0.70"]
    labs = ["paired", "−literature", "−structure", "ECFP4<0.70"]
    x = np.arange(len(stages))
    plotted = {}
    cols = [C["egfr"], C["rtm"], C["gnina"], C["vina"]]
    for p, col in zip(contract, cols):
        r = next(row for row in four if row["pair"] == p)
        ys = [fnum(r[s]) for s in stages]
        ax.plot(x, ys, "-o", color=col, lw=1.1, markersize=4.2, label=p, zorder=3)
        plotted[p] = ys
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=6.2)
    ax.set_ylabel("Remaining InChIKeys")
    ax.set_title("BindingDB-native filter cascade", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2, fontsize=6.0, frameon=False)

    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.06)
    x = np.arange(len(contract))
    w = 0.18
    dual = [int(r["n_dual"]) for r in four]
    ao = [int(r["n_A_only"]) for r in four]
    bo = [int(r["n_B_only"]) for r in four]
    nei = [int(r["n_neither"]) for r in four]
    ax.bar(x - 1.5 * w, dual, w, color=C["dual"], label="dual", zorder=3)
    ax.bar(x - 0.5 * w, ao, w, color=C["a_only"], label="A-only", zorder=3)
    ax.bar(x + 0.5 * w, bo, w, color=C["b_only"], label="B-only", zorder=3)
    ax.bar(x + 1.5 * w, nei, w, color=C["neither"], label="neither", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"], fontsize=6.0)
    ax.set_ylabel("Count after ECFP4 filter")
    ax.set_title("Four-state remainder (not docked)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4, fontsize=6.0, frameon=False)
    n_fail = sum(r["packaged_as_external_evaluation"] == "0" for r in four)
    ax.text(0.98, 0.95, f"external gate pass = {4 - n_fail}/4", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.0, color="#555555")
    PROVENANCE["plotted"]["figS8"] = {
        "cascade": plotted, "n_dual": dual, "n_fail": n_fail,
        "after_ecfp": [fnum(r["after_ecfp_lt_0.70"]) for r in four],
    }
    fig.subplots_adjust(wspace=0.36, left=0.08, right=0.98, top=0.86, bottom=0.30)
    save_all(fig, "FigS7_bindingdb_native_slice")
    plt.close(fig)


def verify(D: dict) -> None:
    errors: list[str] = []

    p1 = PROVENANCE["plotted"]["fig1C"]
    if p1["n_pairs"] != 49:
        errors.append(f"fig1C n_pairs {p1['n_pairs']} != 49")
    if p1["n_thick"] != 4:
        errors.append(f"fig1C n_thick {p1['n_thick']} != 4")
    if p1["stages"]["primary rows after census"] != 8:
        errors.append("fig1C primary rows != 8")
    for pair, val in {"PIK3CA/MTOR": 80, "ACHE/BCHE": 78, "PIK3CA/PIK3CB": 56, "EGFR/HER2": 7, "HDAC1/HDAC6": 93}.items():
        _eq(errors, p1["j0_hardneg"][pair], val, 1e-6, f"fig1C J0 {pair}")
    for pair, val in {"F2/F10": 117, "JAK1/TYK2": 94, "JAK1/JAK2": 53, "PPARG/PPARA": 85, "PPARA/PPARD": 84}.items():
        _eq(errors, p1["dump_hardneg"][pair], val, 1e-6, f"fig1C dump {pair}")
    _eq(errors, p1["complete_case_min"], 0.145119, 1e-6, "fig1C coverage min")
    _eq(errors, p1["complete_case_max"], 0.340172, 1e-6, "fig1C coverage max")

    if "PIK3CA/PIK3CB" in PROVENANCE["plotted"]["fig2A"]["pairs"]:
        errors.append("fig2A must not include withdrawn PIK3CA/PIK3CB")
    if PROVENANCE["plotted"]["fig2A"]["pairs"] != list(PRIMARY_PAIRS):
        errors.append("fig2A pair order")
    for i, p in enumerate(PRIMARY_PAIRS):
        r = primary_row(D, p)
        _eq(errors, PROVENANCE["plotted"]["fig2A"]["DA"][i], r["da"], 5e-4, f"fig2A DA {p}")
        _eq(errors, PROVENANCE["plotted"]["fig2A"]["DB"][i], r["db"], 5e-4, f"fig2A DB {p}")
        rec = PROVENANCE["plotted"]["fig2B"][p]
        _eq(errors, rec["directional"]["y"], r["smin"], 5e-4, f"fig2B dir {p}")
        _eq(errors, rec["neither"]["y"], r["nei"], 5e-4, f"fig2B neither {p}")
        if rec["neither"]["n_neg"] != r["n_neg"]:
            errors.append(f"fig2B n_neg {p}")

    egfr = PROVENANCE["plotted"]["fig2C"]["EGFR/HER2"]
    jak = PROVENANCE["plotted"]["fig2C"]["JAK1/TYK2"]
    src_e = s34_row(D, "EGFR/HER2")
    src_j = s34_row(D, "JAK1/TYK2")
    _eq(errors, egfr["delta"], src_e["delta"], 5e-4, "fig2C EGFR vs CSV")
    _eq(errors, jak["delta"], src_j["delta"], 5e-4, "fig2C JAK1/TYK2 vs CSV")

    e_ecfp = ecfp_row(D, "EGFR/HER2", "D_vs_B")
    _eq(errors, PROVENANCE["plotted"]["fig3A"]["ecfp_db"][0], e_ecfp["ecfp"], 5e-4, "fig3A EGFR ECFP D/B vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig3A"]["vina_db"][0], e_ecfp["vina"], 5e-4, "fig3A EGFR Vina D/B vs CSV")
    ache_n = [primary_row(D, "AChE/BChE").get("n_dual", None)]
    # fig3C n must match scored AChE dual/A/B counts from the loaded tables
    t = D["theta6"]["AChE/BChE"]
    src_n = [int(t["n_dual"]), int(t["n_A_only"]), int(t["n_B_only"])]
    if PROVENANCE["plotted"]["fig3C"]["n"] != src_n:
        errors.append(f"fig3C n {PROVENANCE['plotted']['fig3C']['n']} != CSV {src_n}")

    g_egfr = gnina_indep(D, "EGFR/HER2")
    g_jak = gnina_indep(D, "JAK1/TYK2")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_smin"][0], g_egfr["smin"], 5e-4, "fig4A EGFR GNINA smin vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_neither"][0], g_egfr["nei"], 5e-4, "fig4A EGFR GNINA neither vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_smin"][1], g_jak["smin"], 5e-4, "fig4A JAK1/TYK2 GNINA smin vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_neither"][1], g_jak["nei"], 5e-4, "fig4A JAK1/TYK2 GNINA neither vs CSV")
    if "PIK3CA/PIK3CB" in {r["pair"] for r in PROVENANCE["plotted"]["fig4B"]}:
        errors.append("fig4B must not plot withdrawn PIK3CA/PIK3CB as a peer")
    b = {(r["pair"], r["crystal"]): r for r in PROVENANCE["plotted"]["fig4B"]}
    _eq(errors, b[("PIK3CA/mTOR", "4JPS")]["y"], fnum(D["jps"]["summary_min"]), 5e-4, "fig4B PM 4JPS vs CSV")
    _eq(errors, b[("PIK3CA/mTOR", "5DXT")]["y"], fnum(D["dxt"]["summary_min"]), 5e-4, "fig4B PM 5DXT vs CSV")
    _eq(errors, b[("PIK3CA/mTOR", "4JSX")]["y"], fnum(D["jsx"]["summary_min"]), 5e-4, "fig4B 4JSX vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig4C"]["EGFR/HER2"]["primary"], five_seed_range(D, "EGFR/HER2")["primary"], 5e-4, "fig4C EGFR primary vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig4C"]["F2/F10"]["primary"], five_seed_range(D, "F2/F10")["primary"], 5e-4, "fig4C F2 primary vs CSV")
    if PROVENANCE["plotted"]["fig4C"]["JAK1/TYK2"]["n"] != five_seed_range(D, "JAK1/TYK2")["n"]:
        errors.append("fig4C JAK1/TYK2 n_seeds vs CSV")

    a = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5A"]}
    _eq(errors, a["EGFR/HER2"]["y"], wp_main(D, "EGFR/HER2")["delta"], 5e-4, "fig5A EGFR delta vs CSV")
    _eq(errors, a["F2/F10"]["y"], wp_main(D, "F2/F10")["delta"], 5e-4, "fig5A F2 delta vs CSV")
    if any(r["pair"] == "PIK3CA/PIK3CB" for r in PROVENANCE["plotted"]["fig5A"]):
        errors.append("fig5A must not include withdrawn PIK3CA/PIK3CB")
    bhold = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5B"]}
    _eq(errors, bhold["PIK3CA/mTOR"]["y"], wp_hold(D, "PIK3CA/mTOR")["delta"], 5e-4, "fig5B PM delta vs CSV")
    _eq(errors, bhold["PPARG/PPARA"]["y"], wp_hold(D, "PPARG/PPARA")["delta"], 5e-4, "fig5B PPARG holdout delta vs CSV")
    c5 = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5C"]}
    _eq(errors, c5["PPARG/PPARA"]["hold"], wp_hold(D, "PPARG/PPARA")["smin"], 5e-4, "fig5C PPARG holdout smin vs CSV")
    _eq(errors, c5["JAK1/JAK2"]["hold"], wp_hold(D, "JAK1/JAK2")["smin"], 5e-4, "fig5C JAK1/JAK2 holdout smin vs CSV")

    _eq(errors, PROVENANCE["plotted"]["fig6A"]["PIK3CA/mTOR"][1], primary_row(D, "PIK3CA/mTOR")["smin"], 5e-4, "fig6A PM theta6 vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig6A"]["JAK1/TYK2"][1], primary_row(D, "JAK1/TYK2")["smin"], 5e-4, "fig6A JAK1/TYK2 theta6 vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig6B"]["PM48"], fnum(D["pm110"][("PM48", "vina")]["summary_min"]), 5e-4, "fig6B PM48 vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig6B"]["PM110"], fnum(D["pm110"][("PM110", "vina")]["summary_min"]), 5e-4, "fig6B PM110 vs CSV")
    _eq(errors, PROVENANCE["plotted"]["fig6C"]["E16"], primary_row(D, "PIK3CA/mTOR")["smin"], 5e-4, "fig6C E16 vs CSV")

    s4 = PROVENANCE["plotted"]["figS4"]
    _eq(errors, s4["EGFR/HER2"]["vina"]["smin"], primary_row(D, "EGFR/HER2")["smin"], 5e-4, "S4 EGFR vs CSV")
    _eq(errors, s4["PPARG/PPARA"]["vina"]["smin"], primary_row(D, "PPARG/PPARA")["smin"], 5e-4, "S4 PPARG vs CSV")
    if "PIK3CA/PIK3CB" in s4:
        errors.append("S4 must not include withdrawn PIK3CA/PIK3CB")
    s5 = {r["pair"]: r for r in PROVENANCE["plotted"]["figS5"]}
    _eq(errors, s5["AChE/BChE"]["hold"], wp_hold(D, "AChE/BChE")["smin"], 5e-4, "S5 AChE holdout vs CSV")
    if "PIK3CA/PIK3CB" in s5 or "EGFR/HER2" in s5:
        errors.append("S5 pair set")

    from PIL import Image
    for name, (w_in, h_in) in {
        "Fig2_negative_class_formulation.png": (7.0, None),
        "Fig4_mismatched_pocket.png": (7.0, None),
        "TOC_graphic.tif": (3.25, 1.75),
    }.items():
        im = Image.open(OUT / name)
        if im.mode != "RGB":
            errors.append(f"{name} mode {im.mode}")
        if h_in is not None:
            d = im.info.get("dpi", (300, 300))[0] or 300
            if abs(im.size[0] / d - w_in) > 0.08 or abs(im.size[1] / d - h_in) > 0.08:
                errors.append(f"{name} size {im.size} not {w_in}x{h_in}")
        im.close()

    if errors:
        raise SystemExit("VERIFICATION FAILED:\n" + "\n".join(errors))
    print("v3 verification OK")


def write_lock_and_captions() -> None:
    """Official lock and captions live in docs/; this plotter does not rewrite them."""
    return


def main() -> None:
    print("Official artwork writer is figures/jcim_article/scripts/update_figures_pr32.py")
    print("This module is a helper library (load/helpers). It does not write plotted_values.json.")
    raise SystemExit(0)


if __name__ == "__main__":
    main()
