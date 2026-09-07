#!/usr/bin/env python3
"""JCIM eight-row submission figures from frozen CSVs only.

No hand-typed AUROCs. No decorative data-unrelated arrows. No invented pairs.
Run: python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jcim_figure_style import (  # noqa: E402
    C,
    CENSUS_FIVE,
    FS_ANNO,
    FS_AXIS,
    HOLDOUT_PAIRS,
    ORIGINAL_THREE,
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
PROVENANCE: dict = {"source_files": {}, "plotted": {}}

S34_CONTRAST = "D_vs_B_or_neither_pocketA"
GNINA_INDEP_PAIRS = ["EGFR/HER2", "PIK3CA/mTOR", "JAK1/TYK2"]


def _read(path: Path) -> list[dict]:
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    PROVENANCE["source_files"].setdefault(str(path.relative_to(ROOT)).replace("\\", "/"), len(rows))
    return rows


def fnum(x) -> float:
    return float(x)


def _eq(errors: list[str], a, b, tol: float, msg: str) -> None:
    if abs(float(a) - float(b)) > tol:
        errors.append(f"{msg}: plotted {a} != source {b}")


def load() -> dict:
    j0 = _read(DATA / "jcim_j0j1_v0/tables/j0_strict_label_supply.csv")
    overlap = _read(DATA / "jcim_novelty_v0/tables/complete_case_usable_pchembl_overlap_v1.csv")
    theta = _read(DATA / "jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv")
    form = _read(DATA / "jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv")
    equal = _read(DATA / "jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv")
    ml = _read(DATA / "jcim_strengthen_t0t1_v0/tables/ligand_ml_baseline_scaffold_cv_v1.csv")
    incr = _read(DATA / "jcim_novelty_v0/tables/incremental_information_v1.csv")
    ache = _read(DATA / "jcim_bench_v0/tables/assembled_AChE_BChE.csv")
    gnina_ind = _read(DATA / "jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv")
    jps = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PM48_alt4JPS_v1.csv")[0]
    dxt = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PM48_alt5DXT_v1.csv")[0]
    jsx = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PM48_alt4JSX_v1.csv")[0]
    seeds = _read(DATA / "jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv")
    delta = _read(DATA / "jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv")
    hold_pm = _read(DATA / "jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv")
    pm110 = _read(DATA / "jcim_strengthen_t0t1_v0/tables/pm110_vs_pm48_pocket_matched_v1.csv")
    native = _read(DATA / "jcim_novelty_v0/tables/external_slice_summary_v1.csv")
    desc4 = _read(DATA / "jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv")

    five_t2 = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv")
    five_s34 = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/equal_score_negative_s34_v1.csv")
    five_ecfp = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv")
    five_grid = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/threshold_grid_v1.csv")
    five_ch = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv")
    five_seed_agg = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/fiveseed_summary_min_aggregate_v1.csv")
    five_wp = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/wrong_pocket_by_channel_v1.csv")
    five_xdb = _read(DATA / "jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_crossdb_v1/crossdb_strict_supply_v1.csv")

    theta6 = {r["pair"]: r for r in theta if r["label_rule"] == "theta_6.0"}
    form_by: dict = {}
    for r in form:
        form_by.setdefault(r["pair"], {})[r["contrast"]] = r
    tpsa = defaultdict(list)
    for r in ache:
        if r["cls"] in ("dual", "A_only", "B_only"):
            tpsa[r["cls"]].append(fnum(r["tpsa"]))
    return {
        "j0": j0,
        "overlap": {r["pair"]: r for r in overlap},
        "theta6": theta6,
        "theta_all": theta,
        "form_by": form_by,
        "equal": {(r["pair"], r["contrast"]): r for r in equal},
        "ml": {(r["pair"], r["contrast"]): r for r in ml},
        "incr": incr,
        "tpsa": tpsa,
        "gnina_ind": {(r["pair"], r["contrast"]): r for r in gnina_ind},
        "jps": jps,
        "dxt": dxt,
        "jsx": jsx,
        "seeds": seeds,
        "delta": {(r["pair"], r["set"]): r for r in delta},
        "hold_pm": {(r["pair"], r["variant"]): r for r in hold_pm},
        "pm110": {(r["panel"], r["arm"]): r for r in pm110},
        "native": native,
        "desc4": {r["pair"]: r for r in desc4},
        "five_t2": {r["pair"]: r for r in five_t2},
        "five_s34": {(r["pair"], r["contrast"]): r for r in five_s34},
        "five_ecfp": five_ecfp,
        "five_grid": five_grid,
        "five_ch": {(r["channel"], r["pair"]): r for r in five_ch},
        "five_seed_agg": {r["pair"]: r for r in five_seed_agg},
        "five_wp": {(r["channel"], r["pair"]): r for r in five_wp},
        "five_xdb": {(r["pair"], r["source"], r["rule"]): r for r in five_xdb},
    }


def primary_row(D: dict, pair: str) -> dict:
    """Table-2-comparable Vina θ=6.0 row for one primary pair."""
    if pair in ORIGINAL_THREE:
        t = D["theta6"][pair]
        n = D["form_by"][pair]["D_vs_neither_mean"]
        return {
            "da": fnum(t["auroc_D_vs_A"]),
            "db": fnum(t["auroc_D_vs_B"]),
            "smin": fnum(t["pocket_matched_summary_min"]),
            "lo": fnum(t["ci_lo"]),
            "hi": fnum(t["ci_hi"]),
            "nei": fnum(n["auroc"]),
            "nei_lo": fnum(n["ci_lo"]),
            "nei_hi": fnum(n["ci_hi"]),
            "n_neg": int(n["n_neg"]),
        }
    r = D["five_t2"][pair]
    return {
        "da": fnum(r["auroc_D_vs_A_pocketB"]),
        "db": fnum(r["auroc_D_vs_B_pocketA"]),
        "smin": fnum(r["summary_min"]),
        "lo": fnum(r["ci_lo"]),
        "hi": fnum(r["ci_hi"]),
        "nei": fnum(r["D_vs_neither_vina_mean"]),
        "nei_lo": fnum(r["D_vs_neither_ci_lo"]),
        "nei_hi": fnum(r["D_vs_neither_ci_hi"]),
        "n_neg": int(r["n_neither"]),
    }


def s34_row(D: dict, pair: str) -> dict:
    src = D["equal"] if pair in ORIGINAL_THREE else D["five_s34"]
    r = src[(pair, S34_CONTRAST)]
    return {
        "delta": fnum(r["delta_neither_minus_selective"]),
        "lo": fnum(r["delta_ci_lo"]),
        "hi": fnum(r["delta_ci_hi"]),
        "under": r["underpowered_neither"] == "1",
    }


def ecfp_row(D: dict, pair: str, contrast: str) -> dict:
    if pair in ORIGINAL_THREE:
        m = D["ml"][(pair, contrast)]
        base = next(r for r in D["incr"] if r["pair"] == pair and r["contrast"] == contrast and r["model"] == "ECFP4")
        plus = next(r for r in D["incr"] if r["pair"] == pair and r["contrast"] == contrast and r["model"] == "ECFP4+docking")
        return {
            "vina": fnum(m["auroc_dock_pocket_matched"]),
            "ecfp": fnum(m["auroc_ml"]),
            "ecfp_plus": fnum(plus["cv_auroc"]),
            "ecfp_base": fnum(base["cv_auroc"]),
        }
    base = next(r for r in D["five_ecfp"] if r["pair"] == pair and r["contrast"] == contrast and r["model"] == "ECFP4")
    plus = next(r for r in D["five_ecfp"] if r["pair"] == pair and r["contrast"] == contrast and r["model"] == "ECFP4+docking")
    return {
        "vina": fnum(base["rank_auroc_docking"]),
        "ecfp": fnum(base["cv_auroc"]),
        "ecfp_plus": fnum(plus["cv_auroc"]),
        "ecfp_base": fnum(base["cv_auroc"]),
    }


def best_desc(D: dict, pair: str) -> tuple[str, float]:
    if pair in ORIGINAL_THREE:
        r = D["desc4"][pair]
        return r["best_single_descriptor"], fnum(r["best_single_descriptor_summary_min"])
    r = D["five_t2"][pair]
    return r["best_single_descriptor"], fnum(r["best_single_descriptor_summary_min"])


def five_seed_range(D: dict, pair: str) -> dict:
    if pair in ORIGINAL_THREE:
        vals = [fnum(r["summary_min"]) for r in D["seeds"] if r["pair"] == pair]
        prim = primary_row(D, pair)["smin"]
        return {
            "primary": prim,
            "median": float(np.median(vals)),
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
            "n": len(vals),
        }
    r = D["five_seed_agg"][pair]
    return {
        "primary": fnum(r["primary_20260727"]),
        "median": fnum(r["median_of_per_seed_summary_min"]),
        "min": fnum(r["min_five_seed"]),
        "max": fnum(r["max_five_seed"]),
        "n": int(r["n_seeds"]),
        "crosses": r["crosses_0.5"] == "1",
    }


def wp_main(D: dict, pair: str) -> dict:
    if pair in ORIGINAL_THREE:
        r = D["delta"][(pair, "main_panel")]
        return {
            "delta": fnum(r["delta_matched_minus_wrong"]),
            "lo": fnum(r["delta_ci_lo"]),
            "hi": fnum(r["delta_ci_hi"]),
            "excl": r["ci_excludes_zero"] == "True",
        }
    r = D["five_wp"][("vina_20260727", pair)]
    return {
        "delta": fnum(r["delta_matched_minus_wrong"]),
        "lo": fnum(r["delta_ci_lo"]),
        "hi": fnum(r["delta_ci_hi"]),
        "excl": r["ci_excludes_zero"] == "True",
    }


def wp_hold(D: dict, pair: str) -> dict:
    if pair in ORIGINAL_THREE:
        r = D["delta"][(pair, "unused_pool_holdout")]
        return {
            "delta": fnum(r["delta_matched_minus_wrong"]),
            "lo": fnum(r["delta_ci_lo"]),
            "hi": fnum(r["delta_ci_hi"]),
            "excl": r["ci_excludes_zero"] == "True",
            "smin": fnum(r["matched_summary_min"]),
            "smin_lo": fnum(r["matched_ci_lo"]),
            "smin_hi": fnum(r["matched_ci_hi"]),
        }
    r = D["five_wp"][("holdout_vina_20260727", pair)]
    return {
        "delta": fnum(r["delta_matched_minus_wrong"]),
        "lo": fnum(r["delta_ci_lo"]),
        "hi": fnum(r["delta_ci_hi"]),
        "excl": r["ci_excludes_zero"] == "True",
        "smin": fnum(r["matched_summary_min"]),
        "smin_lo": fnum(r["matched_ci_lo"]),
        "smin_hi": fnum(r["matched_ci_hi"]),
    }


def holdout_smin(D: dict, pair: str) -> dict:
    if pair in ORIGINAL_THREE:
        r = D["hold_pm"][(pair, "pocket_matched_vina")]
        return {
            "y": fnum(r["summary_min"]),
            "lo": fnum(r["summary_min_ci_lo"]),
            "hi": fnum(r["summary_min_ci_hi"]),
        }
    r = D["five_ch"][("holdout_vina_20260727", pair)]
    return {"y": fnum(r["summary_min"]), "lo": fnum(r["ci_lo"]), "hi": fnum(r["ci_hi"])}


def theta_grid_smin(D: dict, pair: str, rule: str) -> float:
    if pair in ORIGINAL_THREE:
        r = next(row for row in D["theta_all"] if row["pair"] == pair and row["label_rule"] == rule)
        return fnum(r["pocket_matched_summary_min"])
    r = next(row for row in D["five_grid"] if row["pair"] == pair and row["label_rule"] == rule)
    return fnum(r["summary_min"])


def gnina_indep(D: dict, pair: str) -> dict:
    if pair == "JAK1/TYK2":
        r = D["five_ch"][("gnina_independent_jak1_tyk2", pair)]
        return {
            "smin": fnum(r["summary_min"]),
            "smin_lo": fnum(r["ci_lo"]),
            "smin_hi": fnum(r["ci_hi"]),
            "nei": fnum(r["D_vs_neither_mean"]),
            "nei_lo": fnum(r["D_vs_neither_ci_lo"]),
            "nei_hi": fnum(r["D_vs_neither_ci_hi"]),
        }
    return {
        "smin": fnum(D["gnina_ind"][(pair, "summary_min")]["auroc"]),
        "nei": fnum(D["gnina_ind"][(pair, "D_vs_neither_mean")]["auroc"]),
        "nei_lo": fnum(D["gnina_ind"][(pair, "D_vs_neither_mean")]["ci_lo"]),
        "nei_hi": fnum(D["gnina_ind"][(pair, "D_vs_neither_mean")]["ci_hi"]),
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
    fig = plt.figure(figsize=(7.0, 6.20))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.00, 1.38], hspace=0.34, wspace=0.16)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    panel_label(ax, "A", x=-0.02, y=1.04)
    ax.text(5.0, 9.55, "Four experimental states", ha="center", fontsize=FS_AXIS, fontweight="bold")
    for x, title, col in ((1.15, "Pocket A", C["a_only"]), (5.55, "Pocket B", C["b_only"])):
        ax.add_patch(FancyBboxPatch(
            (x, 5.85), 3.2, 3.15, boxstyle="round,pad=0.06,rounding_size=0.28",
            facecolor="#F4F7FA", edgecolor=col, lw=1.3, clip_on=False,
        ))
        ax.text(x + 1.6, 8.55, title, ha="center", fontsize=FS_AXIS, fontweight="bold")
    classes = [
        (1.7, "dual", C["dual"], C["dual"]),
        (3.9, "A-only", C["a_only"], "#DDDDDD"),
        (6.1, "B-only", "#DDDDDD", C["b_only"]),
        (8.3, "neither", "#DDDDDD", "#DDDDDD"),
    ]
    ax.text(5.0, 5.15, "four-state dual-target evaluation", ha="center", fontsize=FS_ANNO, color="#555555")
    for x, name, c1, c2 in classes:
        ax.add_patch(Circle((x - 0.28, 3.05), 0.36, facecolor=c1, edgecolor=C["ink"], lw=0.55, clip_on=False))
        ax.add_patch(Circle((x + 0.28, 3.05), 0.36, facecolor=c2, edgecolor=C["ink"], lw=0.55, clip_on=False))
        ax.text(x, 2.15, name, ha="center", fontsize=FS_ANNO)
    ax.text(5.0, 0.85, "A-only / B-only = selectivity hard negatives", ha="center", fontsize=6.5, color="#555555")

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
    y = np.arange(len(stages))
    vals = [s[0] for s in stages]
    cols = [C["other"], C["thick"], C["metal"], C["egfr"], C["vina"], C["dual"]]
    ax.barh(y, vals, color=cols, height=0.68, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([s[1] for s in stages], fontsize=6.2)
    ax.invert_yaxis()
    ax.set_xlabel("Pair count")
    ax.set_xlim(0, 56)
    ax.set_title("J0 scrape to eight-row primary set", fontsize=FS_AXIS, pad=4)
    for i, v in enumerate(vals):
        ax.text(v + 0.8, i, str(v), va="center", fontsize=6.4)
    ax.text(0.02, -0.18, "EGFR kept as supply-limited (B-only = 7), not a J0 thick pair.",
            transform=ax.transAxes, ha="left", fontsize=5.8, color="#555555")

    bx = fig.add_subplot(gs_c[0, 1])
    j0_keys = ["HDAC1/HDAC6", "PIK3CA/MTOR", "ACHE/BCHE", "PIK3CA/PIK3CB", "EGFR/HER2"]
    j0_lab = ["HDAC1/\nHDAC6", "PIK3CA/\nmTOR", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "EGFR/\nHER2"]
    j0_map = {r["pair"]: fnum(r["min_strict_hardneg"]) for r in rows}
    j0_vals = [j0_map[k] for k in j0_keys]
    j0_cols = [C["metal"], C["thick"], C["thick"], C["thick"], C["egfr"]]
    dump_vals = [fnum(D["five_xdb"][(p, "ChEMBL37_dump", "pChEMBL_STANDARD_OK")]["min_strict_hardneg"])
                 for p in CENSUS_FIVE]
    dump_lab = ["F2/\nF10", "JAK1/\nTYK2", "JAK1/\nJAK2", "PPARG/\nPPARA", "PPARA/\nPPARD"]
    x1 = np.arange(5)
    x2 = np.arange(5) + 6.2
    bx.bar(x1, j0_vals, color=j0_cols, width=0.78, zorder=3)
    bx.bar(x2, dump_vals, color=C["vina"], width=0.78, zorder=3)
    bx.axhline(50, color=C["ink"], ls="--", lw=0.8, zorder=2)
    bx.set_xticks(list(x1) + list(x2))
    bx.set_xticklabels(j0_lab + dump_lab, fontsize=5.4)
    bx.set_ylabel("min strict hard-neg.")
    bx.set_ylim(0, 140)
    bx.set_title("Hard-negative supply (two audits)", fontsize=FS_AXIS, pad=4)
    bx.text(2.0, 128, "J0 scrape", ha="center", fontsize=6.0, color="#555555")
    bx.text(8.2, 128, "later ChEMBL 37 dump", ha="center", fontsize=6.0, color="#555555")
    bx.text(0.02, 0.52, "gate ≥50", transform=bx.transAxes, ha="left", va="bottom",
            fontsize=6.0, color="#555555")
    bx.axvline(5.6, color="#DDDDDD", lw=0.8, zorder=1)

    fracs = [fnum(D["overlap"][p]["fraction_union_measured_both"]) for p in
             ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]]
    PROVENANCE["plotted"]["fig1C"] = {
        "n_pairs": n_pairs,
        "n_thick": n_thick,
        "stages": {lab: v for v, lab in stages},
        "j0_hardneg": {k: j0_map[k] for k in j0_keys},
        "dump_hardneg": {p: dump_vals[i] for i, p in enumerate(CENSUS_FIVE)},
        "complete_case_min": min(fracs),
        "complete_case_max": max(fracs),
    }
    fig.subplots_adjust(left=0.08, right=0.98, top=0.94, bottom=0.10)
    save_all(fig, "Fig1_four_state_and_supply")
    plt.close(fig)


def fig2_formulation(D: dict) -> None:
    fig = plt.figure(figsize=(7.0, 7.05))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.15, 1.12], hspace=0.38)
    y = np.arange(len(PRIMARY_PAIRS))
    h = 0.36

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A", x=-0.16, y=1.04)
    da = [primary_row(D, p)["da"] for p in PRIMARY_PAIRS]
    db = [primary_row(D, p)["db"] for p in PRIMARY_PAIRS]
    ax.barh(y + h / 2, da, h, color=C["vina"], label="D vs A-only (pocket B)", zorder=3)
    ax.barh(y - h / 2, db, h, color=C["a_only"], label="D vs B-only (pocket A)", zorder=3)
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("AUROC")
    ax.set_xlim(0, 1.0)
    ax.set_title("Directional pocket-matched arms", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="lower right", fontsize=6.2, frameon=False)
    ax.axhline(2.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig2A"] = {"DA": da, "DB": db, "pairs": list(PRIMARY_PAIRS)}

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B", x=-0.16, y=1.04)
    plotted = {}
    for i, p in enumerate(PRIMARY_PAIRS):
        r = primary_row(D, p)
        plotted[p] = r
        dir_err = np.array([[r["smin"] - r["lo"]], [r["hi"] - r["smin"]]])
        nei_err = np.array([[r["nei"] - r["nei_lo"]], [r["nei_hi"] - r["nei"]]])
        b1 = ax.barh(i + h / 2, r["smin"], h, xerr=dir_err, capsize=1.8,
                     color=C["vina"], ecolor=C["vina"], error_kw={"elinewidth": 0.85},
                     label="directional summary_min" if i == 0 else None, zorder=3)
        b2 = ax.barh(i - h / 2, r["nei"], h, xerr=nei_err, capsize=1.8,
                     color=C["desc"], ecolor=C["desc"], error_kw={"elinewidth": 0.85},
                     label="Dual vs neither (vina_mean)" if i == 0 else None, zorder=3)
        if p == "PIK3CA/mTOR":
            b2[0].set_hatch("///")
            b2[0].set_edgecolor(C["desc"])
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.05, 1.08)
    ax.set_title("Descriptive formulation contrast", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="lower right", fontsize=6.2, frameon=False)
    ax.text(0.98, 0.04, "PIK3CA/mTOR neither n=4, hatched", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.8, color="#666666")
    ax.axhline(2.5, color="#E6E6E6", lw=0.7, zorder=0)
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
    ax.axhline(2.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig2C"] = recs

    fig.subplots_adjust(left=0.18, right=0.97, top=0.96, bottom=0.07)
    save_all(fig, "Fig2_negative_class_formulation")
    plt.close(fig)


def fig3_chemistry(D: dict) -> None:
    fig = plt.figure(figsize=(7.0, 6.85))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 1.00], hspace=0.36, wspace=0.32)

    ax = fig.add_subplot(gs[0, :])
    panel_label(ax, "A", x=-0.08, y=1.04)
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
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.20, 1.05)
    ax.set_title("Vina rank AUROC versus ECFP4 scaffold GroupKFold", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5.5, label="Vina D/A"),
        Line2D([0], [0], marker="s", color=C["gnina"], ls="none", ms=5.0, label="Vina D/B"),
        Line2D([0], [0], marker="^", color=C["desc"], ls="none", ms=5.5, label="ECFP4 D/A"),
        Line2D([0], [0], marker="D", color=C["a_only"], ls="none", ms=5.0, label="ECFP4 D/B"),
    ], loc="lower right", ncol=2, fontsize=6.0, frameon=False)
    ax.axhline(2.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig3A"] = plotted

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B", x=-0.22, y=1.06)
    deltas, labels = [], []
    for p in PRIMARY_PAIRS:
        for contrast, tag in (("D_vs_A", "D/A"), ("D_vs_B", "D/B")):
            r = ecfp_row(D, p, contrast)
            dlt = r["ecfp_plus"] - r["ecfp_base"]
            deltas.append(dlt)
            labels.append(f"{PAIR_SHORT[p]} {tag}")
    yy = np.arange(len(deltas))
    bar_cols = [C["vina"] if abs(d) <= 0.025 else C["egfr"] for d in deltas]
    ax.barh(yy, deltas, color=bar_cols, height=0.72, zorder=3)
    ax.axvline(0, color=C["ink"], lw=0.8, zorder=2)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels, fontsize=5.3)
    ax.invert_yaxis()
    ax.set_xlabel("ΔAUROC (ECFP4+Vina − ECFP4)")
    ax.set_xlim(-0.040, 0.040)
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
    ax.set_title("AChE/BChE", fontsize=FS_AXIS, pad=3)
    PROVENANCE["plotted"]["fig3C"] = {
        "n": ns,
        "mean": [float(np.mean(d)) for d in data],
        "median": [float(np.median(d)) for d in data],
    }

    fig.subplots_adjust(left=0.16, right=0.98, top=0.94, bottom=0.08)
    save_all(fig, "Fig3_ligand_chemistry")
    plt.close(fig)


def fig4_realization(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.95), gridspec_kw={"width_ratios": [1.20, 1.05, 1.15]})
    w = 0.18
    x = np.arange(len(GNINA_INDEP_PAIRS))

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
    ax.bar(x - 1.5 * w, vina_smin, w, color=C["vina"], label="Vina summary_min", zorder=3)
    ax.bar(x - 0.5 * w, g_smin, w, color=C["gnina"], label="GNINA dock summary_min", zorder=3)
    ax.bar(x + 0.5 * w, vina_nei, w, color=C["desc"], label="Vina Dual vs neither", zorder=3)
    ax.bar(x + 1.5 * w, g_nei, w, color=C["a_only"], label="GNINA dock Dual vs neither", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/\nHER2", "PIK3CA/\nmTOR", "JAK1/\nTYK2"], fontsize=6.3)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.08)
    ax.set_title("Independent GNINA pose generation", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=5.4, frameon=False)
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
    ax.set_ylabel("summary_min")
    ax.set_ylim(0.12, 1.02)
    ax.set_title("PIK3CA/mTOR receptor swap", fontsize=FS_AXIS, pad=3)
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
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=6.1)
    ax.invert_yaxis()
    ax.set_xlabel("summary_min across 5 Vina seeds")
    ax.set_xlim(0.22, 0.82)
    ax.set_title("Five-seed range", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="D", color=C["desc"], ls="none", ms=5, label="primary seed"),
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5, label="median"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=5.6, frameon=False)
    ax.axhline(2.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig4C"] = plotted_s

    fig.subplots_adjust(wspace=0.42, left=0.08, right=0.98, top=0.88, bottom=0.30)
    save_all(fig, "Fig4_computational_realization")
    plt.close(fig)


def fig5_mismatched(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 5.55))

    ax = axes[0]
    panel_label(ax, "A", x=-0.28, y=1.04)
    recs = []
    for p in PRIMARY_PAIRS:
        r = wp_main(D, p)
        recs.append({"pair": p, "y": r["delta"], "lo": r["lo"], "hi": r["hi"], "excl": r["excl"]})
    forest_pairs(ax, PRIMARY_PAIRS, recs, "Main panels", "ΔAUROC (matched − mismatched)", (-0.22, 0.36))
    ax.axhline(2.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig5A"] = recs

    ax = axes[1]
    panel_label(ax, "B", x=-0.28, y=1.04)
    recs_b = []
    for p in HOLDOUT_PAIRS:
        r = wp_hold(D, p)
        recs_b.append({"pair": p, "y": r["delta"], "lo": r["lo"], "hi": r["hi"], "excl": r["excl"]})
    forest_pairs(ax, HOLDOUT_PAIRS, recs_b, "Unused-pool holdout Δ",
                 "ΔAUROC (matched − mismatched)", (-0.36, 0.38))
    ax.text(0.02, 0.02, "no EGFR/HER2 holdout", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=5.8, color="#666666")
    ax.axhline(1.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig5B"] = recs_b

    ax = axes[2]
    panel_label(ax, "C", x=-0.28, y=1.04)
    recs_c = []
    for i, p in enumerate(HOLDOUT_PAIRS):
        main = primary_row(D, p)
        h = holdout_smin(D, p)
        recs_c.append({"pair": p, "main": main["smin"], "main_lo": main["lo"], "main_hi": main["hi"],
                       "hold": h["y"], "hold_lo": h["lo"], "hold_hi": h["hi"]})
        ax.errorbar(main["smin"], i + 0.14, xerr=[[main["smin"] - main["lo"]], [main["hi"] - main["smin"]]],
                    fmt="o", color=C["vina"], ecolor=C["vina"], elinewidth=1.1, capsize=1.8, markersize=5.0, zorder=4)
        ax.errorbar(h["y"], i - 0.14, xerr=[[h["y"] - h["lo"]], [h["hi"] - h["y"]]],
                    fmt="s", color=C["holdout"], ecolor=C["holdout"], elinewidth=1.1, capsize=1.8, markersize=4.8, zorder=4)
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(np.arange(len(HOLDOUT_PAIRS)))
    ax.set_yticklabels([PAIR_SHORT[p] for p in HOLDOUT_PAIRS], fontsize=6.2)
    ax.invert_yaxis()
    ax.set_xlabel("summary_min")
    ax.set_xlim(0.12, 1.02)
    ax.set_title("Holdout vs main panel", fontsize=FS_AXIS, pad=3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5.5, label="main"),
        Line2D([0], [0], marker="s", color=C["holdout"], ls="none", ms=5.0, label="holdout"),
    ], loc="lower right", fontsize=5.8, frameon=False)
    ax.axhline(1.5, color="#E6E6E6", lw=0.7, zorder=0)
    PROVENANCE["plotted"]["fig5C"] = recs_c

    fig.subplots_adjust(wspace=0.55, left=0.16, right=0.98, top=0.90, bottom=0.10)
    save_all(fig, "Fig5_mismatched_pocket")
    plt.close(fig)


def fig6_boundary(D: dict) -> None:
    from plot_jcim_si_composites_v1 import _pm48_e8

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.15))
    rules = ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"]
    rule_lab = ["θ=5.5", "θ=6.0", "θ=6.5", "strict"]

    ax = axes[0, 0]
    panel_label(ax, "A", x=-0.18, y=1.04)
    x = np.arange(len(rules))
    grid = {}
    for p in PRIMARY_PAIRS:
        ys = [theta_grid_smin(D, p, rule) for rule in rules]
        ls = "-" if p in ORIGINAL_THREE else "--"
        ax.plot(x, ys, ls, color=PAIR_COLOR[p], lw=1.05, marker="o", markersize=3.6, label=PAIR_SHORT[p], zorder=3)
        grid[p] = ys
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(rule_lab, fontsize=6.5)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0.15, 0.85)
    ax.set_title("Label-threshold sensitivity", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2, fontsize=5.2, frameon=False)
    PROVENANCE["plotted"]["fig6A"] = grid

    ax = axes[0, 1]
    panel_label(ax, "B", x=-0.18, y=1.04)
    y48 = fnum(D["pm110"][("PM48", "vina")]["summary_min"])
    y110 = fnum(D["pm110"][("PM110", "vina")]["summary_min"])
    ax.bar([0], [y48], 0.55, color=C["vina"], zorder=3)
    ax.bar([1], [y110], 0.55, color=C["holdout"], zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["PM48\n(primary)", "PM110"], fontsize=6.5)
    ax.set_ylabel("Vina summary_min")
    ax.set_ylim(0, 1.0)
    ax.set_title("PIK3CA/mTOR panel size", fontsize=FS_AXIS, pad=3)
    ax.text(0, y48 + 0.035, f"{y48:.3f}", ha="center", fontsize=6.2)
    ax.text(1, y110 + 0.035, f"{y110:.3f}", ha="center", fontsize=6.2)
    PROVENANCE["plotted"]["fig6B"] = {"PM48": y48, "PM110": y110}

    ax = axes[1, 0]
    panel_label(ax, "C", x=-0.18, y=1.04)
    e8 = _pm48_e8()
    e16 = primary_row(D, "PIK3CA/mTOR")["smin"]
    ax.bar([0], [e16], 0.55, color=C["vina"], zorder=3)
    ax.bar([1], [e8["summary_min"]], 0.55, color=C["gnina"], zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["E=16\n(primary)", "E=8"], fontsize=6.5)
    ax.set_ylabel("PIK3CA/mTOR summary_min")
    ax.set_ylim(0, 1.0)
    ax.set_title("Search exhaustiveness", fontsize=FS_AXIS, pad=3)
    ax.text(0, e16 + 0.035, f"{e16:.3f}", ha="center", fontsize=6.2)
    ax.text(1, e8["summary_min"] + 0.035, f"{e8['summary_min']:.3f}", ha="center", fontsize=6.2)
    PROVENANCE["plotted"]["fig6C"] = {"E16": e16, "E8": e8["summary_min"], "e8_n": (e8["nD"], e8["nA"], e8["nB"])}

    ax = axes[1, 1]
    panel_label(ax, "D", x=-0.18, y=1.04)
    contract = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
    four = [next(r for r in D["native"] if r["pair"] == p) for p in contract]
    n_fail = sum(r["packaged_as_external_evaluation"] == "0" for r in four)
    n_pass = 4 - n_fail
    ax.bar(["gate fail", "gate pass"], [n_fail, n_pass], color=[C["egfr"], C["thick"]], width=0.55, zorder=3)
    ax.set_ylabel("Number of pairs")
    ax.set_ylim(0, 5)
    ax.set_title("BindingDB-native external gate", fontsize=FS_AXIS, pad=3)
    ax.text(0, n_fail + 0.12, str(n_fail), ha="center", fontsize=8, fontweight="bold")
    ax.text(1, n_pass + 0.12, str(n_pass), ha="center", fontsize=8, fontweight="bold")
    ax.text(0.5, 4.45, "original four-pair contract; not docked", ha="center", fontsize=5.8, color="#666666")
    PROVENANCE["plotted"]["fig6D"] = {"n_fail": n_fail, "n_pass": n_pass, "contract": contract}

    fig.subplots_adjust(wspace=0.38, hspace=0.62, left=0.10, right=0.98, top=0.94, bottom=0.14)
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
        ax.text(0.98, i, name, transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=5.6, color=C["desc"])
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([PAIR_SHORT[p] for p in PRIMARY_PAIRS], fontsize=7.0)
    ax.invert_yaxis()
    ax.set_xlabel("Pocket-matched summary_min AUROC (95% ligand bootstrap CI)")
    ax.set_xlim(0.12, 1.02)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="Vina"),
        Line2D([0], [0], marker="s", color=C["desc"], ls="none", ms=5, label="best single descriptor"),
    ], loc="lower right", fontsize=6.4, frameon=False)
    PROVENANCE["plotted"]["figS4"] = plotted
    fig.subplots_adjust(left=0.16, right=0.86, top=0.94, bottom=0.12)
    save_all(fig, "FigS_pocket_matched_forest")
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
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0.10, 1.02)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="Main panel"),
        Line2D([0], [0], marker="s", color=C["holdout"], ls="none", ms=5.5, label="Unused-pool holdout"),
    ], loc="upper left", fontsize=6.4, frameon=False)
    ax.text(0.98, 0.03, "EGFR/HER2 has no holdout", transform=ax.transAxes, ha="right", fontsize=6.0, color="#666666")
    PROVENANCE["plotted"]["figS5"] = fig_s
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.16)
    save_all(fig, "FigS_unused_pool_holdout")
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
    labels = ["audited\nJ0 pairs", "directional\nn≥10", "formulation\nn≥10", "J0-era\ndocked", "primary\nrows now"]
    ax.bar(range(5), counts, color=[C["vina"], C["desc"], C["thick"], C["egfr"], C["dual"]], width=0.68, zorder=3)
    ax.set_xticks(range(5))
    ax.set_xticklabels(labels, fontsize=5.6)
    ax.set_ylabel("Pair count")
    ax.set_ylim(0, max(counts) + 8)
    ax.set_title("θ=6.0 label census + later primary set", fontsize=FS_AXIS, pad=3)
    for i, v in enumerate(counts):
        ax.text(i, v + 1.0, str(v), ha="center", fontsize=6.5)

    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.06)
    for pair in ORIGINAL_THREE:
        sub = [r for r in and_rows if r["pair"] == pair and r["score"] == "vina_worst"]
        rec = [fnum(r["recall_dual"]) for r in sub]
        prec = [fnum(r["precision_dual"]) for r in sub]
        ax.plot(rec, prec, marker="o", color=PAIR_COLOR[pair], label=pair, lw=1.1, markersize=4.0, zorder=3)
    ax.set_xlabel("Dual recall")
    ax.set_ylabel("Dual precision")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_title("AND-like filter (vina_worst)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=5.4, frameon=False)

    ax = axes[2]
    panel_label(ax, "C", x=-0.18, y=1.06)
    neither, directional = [], []
    for pair in ORIGINAL_THREE:
        n = next(r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "D_vs_neither")
        s = next(r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "summary_min_ecfp4")
        neither.append(fnum(n["ecfp4_groupkfold_auroc"]))
        directional.append(fnum(s["ecfp4_groupkfold_auroc"]))
    x = np.arange(len(ORIGINAL_THREE))
    ax.bar(x - 0.18, neither, 0.36, color=C["vina"], label="Dual vs neither", zorder=3)
    ax.bar(x + 0.18, directional, 0.36, color=C["egfr"], label="ECFP4 summary_min", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nmTOR"], fontsize=6.0)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_ylim(0.35, 1.05)
    ax.set_ylabel("GroupKFold AUROC")
    ax.set_title("Ligand-only full maps (original three)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=5.4, frameon=False)

    PROVENANCE["plotted"]["figS7"] = {
        "n_pairs": n_pairs, "n_dir": n_dir, "n_form": n_form, "n_dock_j0": n_dock_j0,
        "n_primary_now": 8, "neither": neither, "directional": directional,
    }
    fig.subplots_adjust(wspace=0.42, left=0.08, right=0.98, top=0.86, bottom=0.30)
    save_all(fig, "FigS7_posthoc_diagnostics")
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
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2, fontsize=5.4, frameon=False)

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
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4, fontsize=5.4, frameon=False)
    n_fail = sum(r["packaged_as_external_evaluation"] == "0" for r in four)
    ax.text(0.98, 0.95, f"external gate pass = {4 - n_fail}/4", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.0, color="#555555")
    PROVENANCE["plotted"]["figS8"] = {
        "cascade": plotted, "n_dual": dual, "n_fail": n_fail,
        "after_ecfp": [fnum(r["after_ecfp_lt_0.70"]) for r in four],
    }
    fig.subplots_adjust(wspace=0.36, left=0.08, right=0.98, top=0.86, bottom=0.30)
    save_all(fig, "FigS_bindingdb_native_slice_v1")
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

    expected_dir = {
        "EGFR/HER2": (0.4297, 0.2818, 0.5775),
        "AChE/BChE": (0.6058, 0.4370, 0.7303),
        "PIK3CA/mTOR": (0.6921, 0.4702, 0.8133),
        "F2/F10": (0.3448, 0.2109, 0.4773),
        "JAK1/TYK2": (0.3649, 0.2306, 0.503),
        "JAK1/JAK2": (0.5884, 0.4444, 0.7246),
        "PPARG/PPARA": (0.6492, 0.5045, 0.7508),
        "PPARA/PPARD": (0.4463, 0.2958, 0.5841),
    }
    expected_nei = {
        "EGFR/HER2": (0.756, 12),
        "AChE/BChE": (0.6494, 15),
        "PIK3CA/mTOR": (0.5139, 4),
        "F2/F10": (0.5188, 12),
        "JAK1/TYK2": (0.7696, 14),
        "JAK1/JAK2": (0.7299, 14),
        "PPARG/PPARA": (0.6853, 14),
        "PPARA/PPARD": (0.5647, 14),
    }
    for p, (y, lo, hi) in expected_dir.items():
        rec = PROVENANCE["plotted"]["fig2B"][p]["directional"]
        _eq(errors, rec["y"], y, 5e-4, f"fig2B checksum dir {p}")
        _eq(errors, rec["lo"], lo, 5e-4, f"fig2B checksum lo {p}")
        _eq(errors, rec["hi"], hi, 5e-4, f"fig2B checksum hi {p}")
        _eq(errors, PROVENANCE["plotted"]["fig2B"][p]["neither"]["y"], expected_nei[p][0], 5e-4, f"fig2B nei {p}")
        if PROVENANCE["plotted"]["fig2B"][p]["neither"]["n_neg"] != expected_nei[p][1]:
            errors.append(f"fig2B checksum n_neg {p}")

    egfr = PROVENANCE["plotted"]["fig2C"]["EGFR/HER2"]
    jak = PROVENANCE["plotted"]["fig2C"]["JAK1/TYK2"]
    src_e = s34_row(D, "EGFR/HER2")
    src_j = s34_row(D, "JAK1/TYK2")
    _eq(errors, egfr["delta"], src_e["delta"], 5e-4, "fig2C EGFR vs CSV")
    _eq(errors, egfr["delta"], 0.3783, 5e-4, "fig2C EGFR 0.3783")
    _eq(errors, jak["delta"], src_j["delta"], 5e-4, "fig2C JAK1/TYK2 vs CSV")
    _eq(errors, jak["delta"], 0.4438, 5e-4, "fig2C JAK1/TYK2 0.4438")

    _eq(errors, PROVENANCE["plotted"]["fig3A"]["ecfp_db"][0], 0.8895, 5e-4, "fig3A EGFR ECFP D/B")
    _eq(errors, PROVENANCE["plotted"]["fig3A"]["vina_db"][0], 0.4297, 5e-4, "fig3A EGFR Vina D/B")
    if PROVENANCE["plotted"]["fig3C"]["n"] != [27, 25, 28]:
        errors.append(f"fig3C n {PROVENANCE['plotted']['fig3C']['n']}")

    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_smin"][0], 0.2199, 5e-4, "fig4A EGFR GNINA smin")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_neither"][0], 0.7825, 5e-4, "fig4A EGFR GNINA neither")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_smin"][2], 0.3172, 5e-4, "fig4A JAK1/TYK2 GNINA smin")
    _eq(errors, PROVENANCE["plotted"]["fig4A"]["gnina_neither"][2], 0.7048, 5e-4, "fig4A JAK1/TYK2 GNINA neither")
    if "PIK3CA/PIK3CB" in {r["pair"] for r in PROVENANCE["plotted"]["fig4B"]}:
        errors.append("fig4B must not plot withdrawn PIK3CA/PIK3CB as a peer")
    b = {(r["pair"], r["crystal"]): r for r in PROVENANCE["plotted"]["fig4B"]}
    _eq(errors, b[("PIK3CA/mTOR", "4JPS")]["y"], 0.4861, 5e-4, "fig4B PM 4JPS")
    _eq(errors, b[("PIK3CA/mTOR", "5DXT")]["y"], 0.5046, 5e-4, "fig4B PM 5DXT")
    _eq(errors, b[("PIK3CA/mTOR", "4JSX")]["y"], 0.6389, 5e-4, "fig4B 4JSX")
    _eq(errors, PROVENANCE["plotted"]["fig4C"]["EGFR/HER2"]["primary"], 0.4297, 5e-4, "fig4C EGFR primary")
    _eq(errors, PROVENANCE["plotted"]["fig4C"]["F2/F10"]["primary"], 0.3448, 5e-4, "fig4C F2 primary")
    if PROVENANCE["plotted"]["fig4C"]["JAK1/TYK2"]["n"] != 5:
        errors.append("fig4C JAK1/TYK2 n_seeds")
    if PROVENANCE["plotted"]["fig4C"]["F2/F10"].get("crosses"):
        errors.append("fig4C F2 five-seed range must not cross 0.5")

    a = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5A"]}
    _eq(errors, a["EGFR/HER2"]["y"], 0.1697, 5e-4, "fig5A EGFR delta")
    _eq(errors, a["F2/F10"]["y"], -0.0307, 5e-4, "fig5A F2 delta")
    if not a["EGFR/HER2"]["excl"] or a["PIK3CA/mTOR"]["excl"]:
        errors.append("fig5A CI exclude-zero pattern")
    if any(r["pair"] == "PIK3CA/PIK3CB" for r in PROVENANCE["plotted"]["fig5A"]):
        errors.append("fig5A must not include withdrawn PIK3CA/PIK3CB")
    bhold = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5B"]}
    _eq(errors, bhold["PIK3CA/mTOR"]["y"], -0.0225, 5e-4, "fig5B PM delta")
    _eq(errors, bhold["PPARG/PPARA"]["y"], 0.0061, 5e-4, "fig5B PPARG holdout delta")
    if any(r["excl"] for r in PROVENANCE["plotted"]["fig5B"]):
        errors.append("fig5B all CIs should include 0")
    c5 = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5C"]}
    _eq(errors, c5["PPARG/PPARA"]["hold"], 0.535, 5e-4, "fig5C PPARG holdout smin")
    _eq(errors, c5["JAK1/JAK2"]["hold"], 0.6194, 5e-4, "fig5C JAK1/JAK2 holdout smin")

    _eq(errors, PROVENANCE["plotted"]["fig6A"]["PIK3CA/mTOR"][1], 0.6921, 5e-4, "fig6A PM theta6")
    _eq(errors, PROVENANCE["plotted"]["fig6A"]["JAK1/TYK2"][1], 0.3649, 5e-4, "fig6A JAK1/TYK2 theta6")
    _eq(errors, PROVENANCE["plotted"]["fig6B"]["PM48"], 0.6921, 5e-4, "fig6B PM48")
    _eq(errors, PROVENANCE["plotted"]["fig6B"]["PM110"], 0.6483, 5e-4, "fig6B PM110")
    _eq(errors, PROVENANCE["plotted"]["fig6C"]["E16"], 0.6921, 5e-4, "fig6C E16")
    _eq(errors, PROVENANCE["plotted"]["fig6C"]["E8"], 0.6597, 5e-4, "fig6C E8")
    if PROVENANCE["plotted"]["fig6D"]["n_pass"] != 0 or PROVENANCE["plotted"]["fig6D"]["n_fail"] != 4:
        errors.append("fig6D BindingDB gate")

    s4 = PROVENANCE["plotted"]["figS4"]
    _eq(errors, s4["EGFR/HER2"]["vina"]["smin"], 0.4297, 5e-4, "S4 EGFR")
    _eq(errors, s4["PPARG/PPARA"]["vina"]["smin"], 0.6492, 5e-4, "S4 PPARG")
    if "PIK3CA/PIK3CB" in s4:
        errors.append("S4 must not include withdrawn PIK3CA/PIK3CB")
    s5 = {r["pair"]: r for r in PROVENANCE["plotted"]["figS5"]}
    _eq(errors, s5["AChE/BChE"]["hold"], 0.6175, 5e-4, "S5 AChE holdout")
    if "PIK3CA/PIK3CB" in s5 or "EGFR/HER2" in s5:
        errors.append("S5 pair set")
    s7 = PROVENANCE["plotted"]["figS7"]
    if s7["n_dir"] != 17 or s7["n_dock_j0"] != 4 or s7["n_primary_now"] != 8:
        errors.append(f"figS7 census {s7}")
    _eq(errors, s7["neither"][0], 0.9214, 5e-4, "figS7 EGFR Dual vs neither")
    s8 = PROVENANCE["plotted"]["figS8"]
    if s8["n_fail"] != 4:
        errors.append("figS8 gate")
    _eq(errors, s8["after_ecfp"][0], 216, 5e-4, "figS8 EGFR after ECFP")

    from PIL import Image
    for name, (w_in, h_in) in {
        "Fig2_negative_class_formulation.png": (7.0, None),
        "Fig5_mismatched_pocket.png": (7.0, None),
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
    lock = """# Eight-row panel lock (submission)

Branch: `cursor/chembl-exhaustive-pair-census-0b1a`  
Script: `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`  
Rule: every plotted number is read from the CSV in this table. No hand-typed AUROCs. No AI-drawn figures. No decorative arrows unrelated to the data. PIK3CA/PIK3CB is withdrawn and is not a primary-set row.

| Figure | Panel | Content | Unique source |
|---|---|---|---|
| 1 | A | Four ligand states (schematic) | none |
| 1 | B | Pocket-matched directional tasks (schematic; no arrows) | none |
| 1 | C | J0 49-pair scrape → thick gate → HDAC exclusion → historically docked 4 → PIK3CB withdrawal → census +5 → 8 primary; J0 hard-neg bars; later ChEMBL 37 dump hard-neg for the five new pairs | `j0_strict_label_supply.csv`; `five_pair_crossdb_v1/crossdb_strict_supply_v1.csv` ChEMBL37_dump; `complete_case_usable_pchembl_overlap_v1.csv` (J0-docked pairs only) |
| 2 | A | Directional D/A and D/B AUROC, eight primary pairs | original three: `unified_threshold_sensitivity_v2.csv` θ=6.0; five: `five_pair_stack_v1/table2_comparable_theta6_v1.csv` |
| 2 | B | Dual-vs-neither (`vina_mean`) vs directional `summary_min` | original three: formulation CSV; five: same table2 file |
| 2 | C | Fixed pocket-A score, negative-class ΔAUROC | original three: `formulation_equal_score_negative_v1.csv`; five: `equal_score_negative_s34_v1.csv` |
| 3 | A | ECFP4 GroupKFold vs Vina rank AUROC, both arms, eight pairs | original three: `ligand_ml_baseline_scaffold_cv_v1.csv`; five: `ecfp4_incremental_s20s24_v1.csv` |
| 3 | B | ECFP4 → ECFP4+docking ΔAUROC, 16 contrasts | original three: `incremental_information_v1.csv`; five: same ECFP4 file |
| 3 | C | AChE/BChE TPSA jitter + median/IQR | `assembled_AChE_BChE.csv` |
| 4 | A | Independent GNINA pose generation vs Vina: EGFR/HER2, PIK3CA/mTOR, JAK1/TYK2 only | `independent_dock_formulation_v1.csv`; `table2_comparable_by_channel_v1.csv` `gnina_independent_jak1_tyk2` |
| 4 | B | PIK3CA 4L23/4JPS/5DXT and mTOR 4JSX on the receptor-verified PIK3CA/mTOR pair only | alt receptor CSVs + unified_threshold |
| 4 | C | Five-seed `summary_min` range, eight pairs | original three: `multiseed_auroc_by_seed_v2.csv`; five: `fiveseed_summary_min_aggregate_v1.csv` |
| 5 | A | Main matched−mismatched Δ + 95% CI, eight pairs | original three: `wrong_pocket_paired_delta_bootstrap_v1.csv` `set=main_panel`; five: `wrong_pocket_by_channel_v1.csv` `vina_20260727` |
| 5 | B | Holdout matched−mismatched Δ + 95% CI (no EGFR; no withdrawn PIK3CB) | original two: same bootstrap CSV `unused_pool_holdout`; five: `wrong_pocket_by_channel_v1.csv` `holdout_vina_20260727` |
| 5 | C | Holdout vs main `summary_min` | original two: `holdout_pocket_matched_v1.csv`; five: `table2_comparable_by_channel_v1.csv` `holdout_vina_20260727` |
| 6 | A | θ-grid `summary_min`, eight pairs | original three: `unified_threshold_sensitivity_v2.csv`; five: `threshold_grid_v1.csv` |
| 6 | B | PM48 vs PM110 Vina | `pm110_vs_pm48_pocket_matched_v1.csv` |
| 6 | C | PM48 E=16 vs E=8 | E=16 from unified_threshold; E=8 from `scores_vina_E8_best.csv` |
| 6 | D | BindingDB-native gate on the original four-pair contract (0 pass). Five census pairs were not re-opened as BindingDB external. | `external_slice_summary_v1.csv` |
| S4 | — | Eight-row Vina forest + best single descriptor | same Table-2 sources + `descriptor_all_four_directional_v1.csv` |
| S5 | — | Unused-pool holdout vs main, seven pairs | same holdout sources as Fig 5C |
| S7 | — | J0 θ=6.0 census (docked=4 includes later-withdrawn PIK3CB) plus current primary n=8 | `theta6_pair_census_v1.csv` |
| S8 | — | BindingDB-native cascade; original four-pair contract; 0 pass | `external_slice_summary_v1.csv` |
| TOC | — | Four states and Dual-vs-neither ≠ Dual-vs-selective | schematic; no AUROCs; no arrows |

S1–S3, S9, S10 remain the original-set SI records and include the later-withdrawn PIK3CA/PIK3CB row where that is what those CSVs contain. They are not eight-row primary figures.
"""
    (ROOT / "docs" / "FIGURE_PANEL_LOCK_V3.md").write_text(lock, encoding="utf-8", newline="\n")

    text = """# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.
Regenerate: `python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py`

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four experimentally labeled ligand states: dual, A-only, B-only, and neither. A-only and B-only are selectivity hard negatives. (B) Primary tasks are pocket-matched directional AUROCs: Dual versus A-only scored in pocket B, Dual versus B-only scored in pocket A. `summary_min` is a descriptive worst-arm summary. (C) Left: J0 ChEMBL audit of 49 candidate pairs (`j0_strict_label_supply.csv`). Four pairs meet a thick hard-negative gate (min ≥50). HDAC1/HDAC6 is excluded as metal-dependent. EGFR/HER2 is retained as a supply-limited case (strict B-only = 7), giving four historically docked pairs. PIK3CA/PIK3CB was later withdrawn after a receptor-identity failure. A later ChEMBL 37 dump census added five ordinary pairs, leaving eight primary rows. This is collection completion, not a claim that eight pairs were in the J0 thick set. Right: min strict hard-negatives on the J0 scrape (left group) and on the later ChEMBL 37 dump for the five census pairs (right group; `crossdb_strict_supply_v1.csv`). Complete-case map coverage on the J0-docked pairs is 14.5%–34.0%. Cross-database counts for the original scrape are Figure S2.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Same frozen AutoDock Vina scores, unified θ = 6.0, eight primary rows. PIK3CA/PIK3CB is withdrawn and is not plotted. Horizontal gray rules separate the 2026-07-23 three from the five post-census pairs. (A) Directional Dual versus A-only (pocket B) and Dual versus B-only (pocket A). Original three: `unified_threshold_sensitivity_v2.csv`. Five: `table2_comparable_theta6_v1.csv`. (B) Descriptive comparison of directional `summary_min` with Dual versus neither using per-ligand `vina_mean`. These two columns differ in both negative class and score aggregation. PIK3CA/mTOR Dual versus neither is hatched (neither n = 4). (C) Pocket A score held fixed; only the negative class is replaced (B-only versus neither). EGFR/HER2 ΔAUROC = 0.378 [0.205, 0.547]; JAK1/TYK2 reproduces the gap (0.444 [0.263, 0.620]). Diamond, underpowered neither. Vertical dashed line, zero.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Scaffold GroupKFold ECFP4 logistic AUROC versus pocket-matched Vina rank AUROC on both directional arms for the eight primary pairs. EGFR/HER2 Dual versus B-only: ECFP4 0.8895 versus Vina 0.4297. Five-pair ECFP4 from `ecfp4_incremental_s20s24_v1.csv`. (B) Change in GroupKFold AUROC when the pocket-matched Vina score is added to ECFP4 (16 contrasts). (C) AChE/BChE TPSA by class: individual ligands (jittered) with median and IQR (`assembled_AChE_BChE.csv`). n = 27/25/28.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation (not CNN rescoring of Vina poses) on EGFR/HER2, PIK3CA/mTOR, and JAK1/TYK2 only. Independent search was not run on F2/F10, JAK1/JAK2, or the PPAR pairs. JAK1/TYK2 independent GNINA: `summary_min` 0.317 [0.183, 0.463], Dual versus neither 0.705. (B) Replacing PIK3CA 4L23 with 4JPS or 5DXT while holding mTOR frozen: PIK3CA/mTOR `summary_min` 0.692 → 0.486 / 0.505. 4JSX is an mTOR-pocket swap. The parallel swap on the withdrawn PIK3CA/PIK3CB pair is not plotted as a peer contrast. (C) Directional `summary_min` across five frozen Vina seeds; diamond, production seed 20260727. No post-census five-seed range crossed 0.5.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Δ = matched-pocket `summary_min` − mismatched-pocket `summary_min`, ligand bootstrap B = 2000. Matched uses Dual versus A-only in pocket B and Dual versus B-only in pocket A; mismatched swaps those score channels. This is a scoring-channel control, not redocking into a physically wrong site. Dark, CI excludes 0; gray, CI includes 0. (A) Main panels, eight primary pairs. EGFR/HER2 and AChE/BChE CIs exclude 0; the five post-census production-Vina CIs include 0. (B) Unused-pool holdout Δ. All seven CIs include 0. EGFR/HER2 has no holdout; withdrawn PIK3CA/PIK3CB is omitted. (C) Holdout versus main-panel `summary_min`. PPARG/PPARA holdout is 0.535 [0.350, 0.717]; JAK1/JAK2 stays same-direction (0.619 [0.420, 0.749]; drawn 20/20/18).

## Figure 6. Robustness checks and evidence boundary.

(A) Pocket-matched `summary_min` on the unified label-threshold grid for the eight primary pairs. Solid, 2026-07-23 three; dashed, post-census five. (B) PIK3CA/mTOR PM48 versus PM110 Vina. (C) PM48 exhaustiveness 16 versus 8, recomputed from `scores_vina_E8_best.csv` with the same pocket-matched definition. (D) BindingDB-native 202608 slice under the original four-pair contract: zero pairs meet the pre-frozen external gate; nothing was docked. The five census pairs were not re-opened as a BindingDB external set.

## Figure S1. Protocol and panel sensitivities.

Original-set protocol grid, GNINA CNN rescoring of Vina poses, PM48 versus PM110, and exhaustiveness. This SI record still includes the later-withdrawn PIK3CA/PIK3CB row where that is what the source CSVs contain. Independent GNINA pose generation is Figure 4A.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged original-scrape sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`. Five-pair ChEMBL 37 dump supply is Figure 1C (right group).

## Figure S3. Additional paired bootstrap differences.

Descriptor and scaffold-versus-random leakage checks on the original docked set. Matched-versus-mismatched main/holdout Δ CIs are Figure 5.

## Figure S4. Pocket-matched summary_min forest.

Vina CIs and the best single-descriptor reference on the eight primary rows. PIK3CA/PIK3CB is omitted.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` on the seven pairs that have a holdout. EGFR/HER2 has no holdout. Mismatched-pocket Δ CIs are Figure 5B.

## Figure S7. Post-hoc formulation and screening diagnostics.

θ = 6.0 J0 pair census (`docked_in_this_paper` = 4 includes the later-withdrawn PIK3CA/PIK3CB row), current primary n = 8, AND-like dual filter on the original three, and ligand-only full-map ECFP4 on the original three. Not docking upgrades and not a replacement for Table 2.

## Figure S8. BindingDB-native slice.

Filter cascade and remaining four-state counts after literature, structure, and ECFP4 < 0.70 on the original four-pair contract (`external_slice_summary_v1.csv`). Zero of four pairs meet the pre-frozen external gate; nothing was docked. Five census pairs were not re-opened as BindingDB external.

## Figure S9. Additional ligand-structure controls.

Prespecified descriptors, covariate-adjusted logistic AUROC, and matched-subset weak-arm tests on the original docked set. The Vina-only logistic AUROC is not the Table 2 rank AUROC.

## Figure S10. Matched versus mismatched point estimates.

Bar charts of matched versus mismatched `summary_min` on the original docked set. Paired Δ CIs are Figure 5.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs and no decorative arrows.
"""
    (OUT / "CAPTIONS.md").write_text(text, encoding="utf-8", newline="\n")


def main() -> None:
    apply_style()
    D = load()
    fig1_framework(D)
    fig2_formulation(D)
    fig3_chemistry(D)
    fig4_realization(D)
    fig5_mismatched(D)
    fig6_boundary(D)
    toc_graphic()
    fig_s4_forest(D)
    fig_s5_holdout(D)
    fig_s7_diagnostics(D)
    fig_s8_bindingdb(D)
    write_lock_and_captions()
    (OUT / "plotted_values.json").write_text(
        json.dumps(PROVENANCE, indent=2, default=str), encoding="utf-8"
    )
    verify(D)
    print("wrote", OUT)
    for p in sorted(OUT.glob("Fig*.png")) + sorted(OUT.glob("TOC*")):
        print(" ", p.name, p.stat().st_size)


if __name__ == "__main__":
    main()
