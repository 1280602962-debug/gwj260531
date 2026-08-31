#!/usr/bin/env python3
"""JCIM submission figures (six main + TOC) from frozen CSVs only.

No hand-typed AUROCs. No AI-drawn artwork.
Run: python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v2.py
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
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.lines import Line2D
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jcim_figure_style import (  # noqa: E402
    C,
    FS_ANNO,
    FS_AXIS,
    OUT,
    PAIR_ORDER,
    PAIR_SHORT,
    ROOT,
    apply_style,
    panel_label,
    save_all,
)

DATA = ROOT / "data"
PROVENANCE: dict = {"source_files": {}, "plotted": {}}
TICK = ["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"]
HOLD_PAIRS = ["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
HOLD_TICK = ["AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"]


def _read(path: Path) -> list[dict]:
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    PROVENANCE["source_files"].setdefault(str(path.relative_to(ROOT)).replace("\\", "/"), len(rows))
    return rows


def fnum(x) -> float:
    return float(x)


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
    pab_jps = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PAB_alt4JPS_v1.csv")[0]
    pab_dxt = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PAB_alt5DXT_v1.csv")[0]
    seeds = _read(DATA / "jcim_multiseed_v0/tables/multiseed_auroc_by_seed_v2.csv")
    delta = _read(DATA / "jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv")
    hold_m = _read(DATA / "jcim_holdout_v0/tables/holdout_matched_wrong_pocket_summary_v1.csv")
    pm110 = _read(DATA / "jcim_strengthen_t0t1_v0/tables/pm110_vs_pm48_pocket_matched_v1.csv")
    native = _read(DATA / "jcim_novelty_v0/tables/external_slice_summary_v1.csv")

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
        "pab_jps": pab_jps,
        "pab_dxt": pab_dxt,
        "seeds": seeds,
        "delta": {(r["pair"], r["set"]): r for r in delta},
        "hold_match": {(r["pair"], r["family"], r["aggregation"]): r for r in hold_m},
        "pm110": {(r["panel"], r["arm"]): r for r in pm110},
        "native": native,
    }


def fig1_framework(D: dict) -> None:
    fig = plt.figure(figsize=(7.0, 5.90))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.00, 1.32], hspace=0.32, wspace=0.14)

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
        (0.4, 5.55, "Dual vs A-only", "score in pocket B"),
        (0.4, 2.35, "Dual vs B-only", "score in pocket A"),
    ]
    for x, y, head, mid in boxes:
        ax.add_patch(FancyBboxPatch(
            (x, y), 9.2, 2.55, boxstyle="round,pad=0.05,rounding_size=0.22",
            facecolor="#F7F7F7", edgecolor="#CCCCCC", lw=0.8, clip_on=False,
        ))
        ax.text(x + 0.35, y + 1.75, head, ha="left", fontsize=FS_AXIS, fontweight="bold")
        ax.annotate("", xy=(x + 8.5, y + 1.05), xytext=(x + 4.2, y + 1.05),
                    arrowprops=dict(arrowstyle="-|>", color=C["vina"], lw=1.25))
        ax.text(x + 0.35, y + 0.95, mid, ha="left", va="center", fontsize=FS_ANNO, color=C["vina"])
    ax.text(5.0, 1.15, r"summary$_{\mathrm{min}}$ = min(AUROC$_{D/A}$, AUROC$_{D/B}$)",
            ha="center", fontsize=FS_ANNO)
    ax.text(5.0, 0.40, "Descriptive worst-arm summary; both arms are reported.",
            ha="center", fontsize=6.5, color="#555555")

    gs_c = gs[1, :].subgridspec(1, 2, width_ratios=[0.92, 1.28], wspace=0.24)
    ax = fig.add_subplot(gs_c[0, 0])
    panel_label(ax, "C", x=-0.08, y=1.08)
    rows = [r for r in D["j0"] if (r.get("min_strict_hardneg") or "").strip()]
    n_pairs = len(rows)
    n_thick = sum(1 for r in rows if fnum(r["min_strict_hardneg"]) >= 50)
    highlight = {
        "HDAC1/HDAC6": C["metal"],
        "PIK3CA/MTOR": C["thick"],
        "ACHE/BCHE": C["thick"],
        "PIK3CA/PIK3CB": C["thick"],
        "EGFR/HER2": C["egfr"],
    }
    plotted_hl = {}
    for r in rows:
        if r["pair"] in highlight:
            plotted_hl[r["pair"]] = fnum(r["min_strict_hardneg"])
    fracs = [fnum(D["overlap"][p]["fraction_union_measured_both"]) for p in PAIR_ORDER]

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    stages = [
        (7.15, f"{n_pairs} candidate pairs", "#A6CEE3"),
        (4.15, f"{n_thick} pairs with min hard-neg ≥50", C["thick"]),
        (1.15, "K = 4 docked  (HDAC metal-excluded;\nEGFR kept as supply-limited)", C["vina"]),
    ]
    for y, txt, col in stages:
        ax.add_patch(FancyBboxPatch(
            (0.35, y), 9.3, 2.15, boxstyle="round,pad=0.05,rounding_size=0.20",
            facecolor="#F7FBFF", edgecolor=col, lw=1.15, clip_on=False,
        ))
        ax.text(5.0, y + 1.05, txt, ha="center", va="center", fontsize=FS_ANNO)
    ax.annotate("", xy=(5.0, 6.30), xytext=(5.0, 7.15),
                arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=1.05))
    ax.annotate("", xy=(5.0, 3.30), xytext=(5.0, 4.15),
                arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=1.05))
    ax.set_title("Candidate pairs → directional hard-negative supply", loc="left", fontsize=FS_AXIS, pad=8)

    bx = fig.add_subplot(gs_c[0, 1])
    names = ["HDAC1/\nHDAC6", "PIK3CA/\nmTOR", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "EGFR/\nHER2"]
    keys = ["HDAC1/HDAC6", "PIK3CA/MTOR", "ACHE/BCHE", "PIK3CA/PIK3CB", "EGFR/HER2"]
    vals = [plotted_hl[k] for k in keys]
    cols = [highlight[k] for k in keys]
    xb = np.arange(5)
    bx.bar(xb, vals, color=cols, width=0.72, zorder=3)
    bx.axhline(50, color=C["ink"], ls="--", lw=0.8, zorder=2)
    bx.set_xticks(xb)
    bx.set_xticklabels(names, fontsize=6.0)
    bx.set_ylabel("min strict hard-neg.")
    bx.set_ylim(0, 115)
    bx.text(0.02, 0.52, "gate ≥50", transform=bx.transAxes, ha="left", va="bottom",
            fontsize=6.0, color="#555555")
    short = {"EGFR/HER2": "EGFR", "AChE/BChE": "AChE", "PIK3CA/PIK3CB": "PIK3CB", "PIK3CA/mTOR": "mTOR"}
    cov_txt = "; ".join(
        f"{short[p]} {100 * fnum(D['overlap'][p]['fraction_union_measured_both']):.1f}%"
        for p in PAIR_ORDER
    )
    fig.subplots_adjust(left=0.05, right=0.98, top=0.94, bottom=0.12)
    fig.text(
        0.52,
        0.025,
        f"Complete-case coverage: {cov_txt}. Cross-database counts: Figure S2.",
        ha="center",
        va="bottom",
        fontsize=6.0,
        color="#555555",
    )

    PROVENANCE["plotted"]["fig1C"] = {
        "n_pairs": n_pairs,
        "n_thick": n_thick,
        "highlighted": plotted_hl,
        "complete_case": {p: D["overlap"][p]["fraction_union_measured_both"] for p in PAIR_ORDER},
        "complete_case_min": min(fracs),
        "complete_case_max": max(fracs),
    }
    fig.subplots_adjust(left=0.05, right=0.98, top=0.94, bottom=0.12)
    save_all(fig, "Fig1_four_state_and_supply")
    plt.close(fig)


def fig2_formulation(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.65), gridspec_kw={"width_ratios": [1.05, 1.15, 1.15]})
    w = 0.36
    x = np.arange(len(PAIR_ORDER))

    ax = axes[0]
    panel_label(ax, "A", x=-0.22, y=1.06)
    da = [fnum(D["theta6"][p]["auroc_D_vs_A"]) for p in PAIR_ORDER]
    db = [fnum(D["theta6"][p]["auroc_D_vs_B"]) for p in PAIR_ORDER]
    ax.bar(x - w / 2, da, w, color=C["vina"], label="D vs A-only (pocket B)", zorder=3)
    ax.bar(x + w / 2, db, w, color=C["a_only"], label="D vs B-only (pocket A)", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.2)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.0)
    ax.set_title("Directional pocket-matched arms", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=5.8, frameon=False)
    PROVENANCE["plotted"]["fig2A"] = {"DA": da, "DB": db}

    ax = axes[1]
    panel_label(ax, "B", x=-0.22, y=1.06)
    dir_y, dir_lo, dir_hi, nei_y, nei_lo, nei_hi, n_neg = [], [], [], [], [], [], []
    plotted = {}
    for p in PAIR_ORDER:
        t = D["theta6"][p]
        r = D["form_by"][p]["D_vs_neither_mean"]
        dy, dlo, dhi = fnum(t["pocket_matched_summary_min"]), fnum(t["ci_lo"]), fnum(t["ci_hi"])
        ny, nlo, nhi = fnum(r["auroc"]), fnum(r["ci_lo"]), fnum(r["ci_hi"])
        nn = int(r["n_neg"])
        dir_y.append(dy); dir_lo.append(dlo); dir_hi.append(dhi)
        nei_y.append(ny); nei_lo.append(nlo); nei_hi.append(nhi); n_neg.append(nn)
        plotted[p] = {"directional": {"y": dy, "lo": dlo, "hi": dhi},
                      "neither": {"y": ny, "lo": nlo, "hi": nhi, "n_neg": nn}}
    dir_y = np.asarray(dir_y); nei_y = np.asarray(nei_y)
    dir_err = np.vstack([dir_y - dir_lo, np.asarray(dir_hi) - dir_y])
    nei_err = np.vstack([nei_y - nei_lo, np.asarray(nei_hi) - nei_y])
    ax.bar(x - w / 2, dir_y, w, yerr=dir_err, capsize=2.0, color=C["vina"], ecolor=C["vina"],
           error_kw={"elinewidth": 0.9}, label="directional summary_min", zorder=3)
    bars = ax.bar(x + w / 2, nei_y, w, yerr=nei_err, capsize=2.0, color=C["desc"], ecolor=C["desc"],
                  error_kw={"elinewidth": 0.9}, label="Dual vs neither (vina_mean)", zorder=3)
    bars[3].set_hatch("///")
    bars[3].set_edgecolor(C["desc"])
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.2)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.10, 1.12)
    ax.set_title("Descriptive formulation contrast", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=5.8, frameon=False)
    ax.text(0 - w / 2, dir_hi[0] + 0.025, f"{dir_y[0]:.3f}", ha="center", fontsize=5.5, color=C["vina"])
    ax.text(0 + w / 2, nei_hi[0] + 0.025, f"{nei_y[0]:.3f}", ha="center", fontsize=5.5, color=C["desc"])
    ax.text(3.0, 0.16, "neither n=4", ha="center", fontsize=5.8, color="#666666")
    PROVENANCE["plotted"]["fig2B"] = plotted

    ax = axes[2]
    panel_label(ax, "C", x=-0.22, y=1.06)
    contrast = "D_vs_B_or_neither_pocketA"
    ys, los, his, under = [], [], [], []
    for p in PAIR_ORDER:
        r = D["equal"][(p, contrast)]
        ys.append(fnum(r["delta_neither_minus_selective"]))
        los.append(fnum(r["delta_ci_lo"]))
        his.append(fnum(r["delta_ci_hi"]))
        under.append(r["underpowered_neither"] == "1")
    y = np.arange(len(PAIR_ORDER))
    for i, p in enumerate(PAIR_ORDER):
        col = C["egfr"] if p == "EGFR/HER2" else C["vina"]
        ax.plot([los[i], his[i]], [y[i], y[i]], color=col, lw=1.4, zorder=3)
        m = "D" if under[i] else "o"
        ax.plot(ys[i], y[i], m, color=col, markersize=6.0 if p == "EGFR/HER2" else 5.2, zorder=4)
    ax.axvline(0, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"], fontsize=6.5)
    ax.set_xlabel("ΔAUROC (neither − B-only)\nsame pocket A score")
    ax.set_xlim(-0.62, 0.72)
    ax.set_title("Fixed score; swap negative class", fontsize=FS_AXIS, pad=3)
    ax.invert_yaxis()
    ax.text(ys[0] + 0.02, 0.0, f"{ys[0]:.3f}", va="center", ha="left", fontsize=6.0, color=C["egfr"])
    PROVENANCE["plotted"]["fig2C"] = {
        p: {"delta": ys[i], "lo": los[i], "hi": his[i], "underpowered": under[i]}
        for i, p in enumerate(PAIR_ORDER)
    }

    fig.subplots_adjust(wspace=0.48, left=0.07, right=0.98, top=0.88, bottom=0.28)
    save_all(fig, "Fig2_negative_class_formulation")
    plt.close(fig)


def fig3_chemistry(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.70), gridspec_kw={"width_ratios": [1.2, 1.15, 0.95]})

    ax = axes[0]
    panel_label(ax, "A", x=-0.20, y=1.06)
    x = np.arange(len(PAIR_ORDER))
    w = 0.20
    vina_a = [fnum(D["ml"][(p, "D_vs_A")]["auroc_dock_pocket_matched"]) for p in PAIR_ORDER]
    vina_b = [fnum(D["ml"][(p, "D_vs_B")]["auroc_dock_pocket_matched"]) for p in PAIR_ORDER]
    ecfp_a = [fnum(D["ml"][(p, "D_vs_A")]["auroc_ml"]) for p in PAIR_ORDER]
    ecfp_b = [fnum(D["ml"][(p, "D_vs_B")]["auroc_ml"]) for p in PAIR_ORDER]
    ax.bar(x - 1.5 * w, vina_a, w, color=C["vina"], label="Vina D/A", zorder=3)
    ax.bar(x - 0.5 * w, vina_b, w, color="#56B4E9", label="Vina D/B", zorder=3)
    ax.bar(x + 0.5 * w, ecfp_a, w, color=C["desc"], label="ECFP4 D/A", zorder=3)
    ax.bar(x + 1.5 * w, ecfp_b, w, color="#D55E00", label="ECFP4 D/B", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.2)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2, fontsize=5.6, frameon=False)
    PROVENANCE["plotted"]["fig3A"] = {
        "vina_DA": vina_a, "vina_DB": vina_b, "ecfp_DA": ecfp_a, "ecfp_DB": ecfp_b,
    }

    ax = axes[1]
    panel_label(ax, "B", x=-0.20, y=1.06)
    deltas = []
    labels = []
    for p in PAIR_ORDER:
        for contrast in ("D_vs_A", "D_vs_B"):
            base = next(r for r in D["incr"] if r["pair"] == p and r["contrast"] == contrast and r["model"] == "ECFP4")
            plus = next(r for r in D["incr"] if r["pair"] == p and r["contrast"] == contrast and r["model"] == "ECFP4+docking")
            dlt = fnum(plus["cv_auroc"]) - fnum(base["cv_auroc"])
            deltas.append(dlt)
            labels.append(f"{p.split('/')[0]}\n{contrast[-1]}")
    y = np.arange(len(deltas))
    cols = [C["vina"] if abs(d) <= 0.0205 else C["egfr"] for d in deltas]
    ax.barh(y, deltas, color=cols, height=0.72, zorder=3)
    ax.axvline(0, color=C["ink"], lw=0.8, zorder=2)
    ylabels = []
    for p in PAIR_ORDER:
        for contrast in ("D/A", "D/B"):
            ylabels.append(f"{PAIR_SHORT.get(p, p)} {contrast}")
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=5.6)
    ax.set_xlabel("ΔAUROC (ECFP4+Vina − ECFP4)")
    ax.set_xlim(-0.035, 0.035)
    ax.invert_yaxis()
    PROVENANCE["plotted"]["fig3B_deltas"] = deltas
    PROVENANCE["plotted"]["fig3B_max_abs"] = float(max(abs(d) for d in deltas))

    ax = axes[2]
    panel_label(ax, "C", x=-0.24, y=1.06)
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
    means = [float(np.mean(d)) for d in data]
    medians = [float(np.median(d)) for d in data]
    PROVENANCE["plotted"]["fig3C"] = {"n": ns, "mean": means, "median": medians}

    fig.subplots_adjust(wspace=0.48, left=0.08, right=0.98, top=0.86, bottom=0.28)
    save_all(fig, "Fig3_ligand_chemistry")
    plt.close(fig)


def fig4_realization(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.75), gridspec_kw={"width_ratios": [1.15, 1.25, 0.95]})

    ax = axes[0]
    panel_label(ax, "A", x=-0.20, y=1.06)
    pairs = ["EGFR/HER2", "PIK3CA/mTOR"]
    x = np.arange(len(pairs))
    w = 0.18
    vina_smin = [fnum(D["theta6"][p]["pocket_matched_summary_min"]) for p in pairs]
    vina_nei = [fnum(D["form_by"][p]["D_vs_neither_mean"]["auroc"]) for p in pairs]
    g_smin = [fnum(D["gnina_ind"][(p, "summary_min")]["auroc"]) for p in pairs]
    g_nei = [fnum(D["gnina_ind"][(p, "D_vs_neither_mean")]["auroc"]) for p in pairs]
    ax.bar(x - 1.5 * w, vina_smin, w, color=C["vina"], label="Vina summary_min", zorder=3)
    ax.bar(x - 0.5 * w, g_smin, w, color=C["gnina"], label="GNINA dock summary_min", zorder=3)
    ax.bar(x + 0.5 * w, vina_nei, w, color=C["desc"], label="Vina Dual vs neither", zorder=3)
    ax.bar(x + 1.5 * w, g_nei, w, color="#D55E00", label="GNINA dock Dual vs neither", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/HER2", "PIK3CA/mTOR"], fontsize=6.5)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.08)
    ax.set_title("Independent GNINA pose generation", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=5.4, frameon=False)
    PROVENANCE["plotted"]["fig4A"] = {
        "vina_smin": vina_smin, "gnina_smin": g_smin, "vina_neither": vina_nei, "gnina_neither": g_nei,
    }

    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.06)
    crystals = ["4L23", "4JPS", "5DXT"]
    pm = [
        (fnum(D["theta6"]["PIK3CA/mTOR"]["pocket_matched_summary_min"]),
         fnum(D["theta6"]["PIK3CA/mTOR"]["ci_lo"]), fnum(D["theta6"]["PIK3CA/mTOR"]["ci_hi"])),
        (fnum(D["jps"]["summary_min"]), fnum(D["jps"]["summary_min_ci_lo"]), fnum(D["jps"]["summary_min_ci_hi"])),
        (fnum(D["dxt"]["summary_min"]), fnum(D["dxt"]["summary_min_ci_lo"]), fnum(D["dxt"]["summary_min_ci_hi"])),
    ]
    pab = [
        (fnum(D["theta6"]["PIK3CA/PIK3CB"]["pocket_matched_summary_min"]),
         fnum(D["theta6"]["PIK3CA/PIK3CB"]["ci_lo"]), fnum(D["theta6"]["PIK3CA/PIK3CB"]["ci_hi"])),
        (fnum(D["pab_jps"]["summary_min"]), fnum(D["pab_jps"]["summary_min_ci_lo"]), fnum(D["pab_jps"]["summary_min_ci_hi"])),
        (fnum(D["pab_dxt"]["summary_min"]), fnum(D["pab_dxt"]["summary_min_ci_lo"]), fnum(D["pab_dxt"]["summary_min_ci_hi"])),
    ]
    jsx = (fnum(D["jsx"]["summary_min"]), fnum(D["jsx"]["summary_min_ci_lo"]), fnum(D["jsx"]["summary_min_ci_hi"]))
    xp = np.arange(3)
    fig4b = []
    ax.plot(xp - 0.08, [pm[i][0] for i in range(3)], color=C["vina"], lw=0.9, zorder=2)
    ax.plot(xp + 0.08, [pab[i][0] for i in range(3)], color=C["a_only"], lw=0.9, zorder=2)
    for i, (y, lo, hi) in enumerate(pm):
        ax.errorbar(i - 0.08, y, yerr=[[y - lo], [hi - y]], fmt="o", color=C["vina"],
                    ecolor=C["vina"], elinewidth=1.2, capsize=2.0, markersize=6.0, zorder=4)
        fig4b.append({"pair": "PIK3CA/mTOR", "crystal": crystals[i], "y": y, "lo": lo, "hi": hi})
    for i, (y, lo, hi) in enumerate(pab):
        ax.errorbar(i + 0.08, y, yerr=[[y - lo], [hi - y]], fmt="s", color=C["a_only"],
                    ecolor=C["a_only"], elinewidth=1.2, capsize=2.0, markersize=5.5, zorder=4)
        fig4b.append({"pair": "PIK3CA/PIK3CB", "crystal": crystals[i], "y": y, "lo": lo, "hi": hi})
    ax.axvline(2.58, color="#CCCCCC", lw=0.75, ls=":", zorder=1)
    ax.errorbar(3.15, jsx[0], yerr=[[jsx[0] - jsx[1]], [jsx[2] - jsx[0]]], fmt="^",
                color=C["gnina"], ecolor=C["gnina"], elinewidth=1.2, capsize=2.0, markersize=6.0, zorder=4)
    fig4b.append({"pair": "PIK3CA/mTOR", "crystal": "4JSX", "y": jsx[0], "lo": jsx[1], "hi": jsx[2]})
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks([0, 1, 2, 3.15])
    ax.set_xticklabels(["4L23\nPIK3CA", "4JPS\nPIK3CA", "5DXT\nPIK3CA", "4JSX\nmTOR only"], fontsize=6.0)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0.12, 1.02)
    ax.set_xlim(-0.45, 3.55)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="PIK3CA/mTOR"),
        Line2D([0], [0], marker="s", color=C["a_only"], ls="none", ms=5.5, label="PIK3CA/PIK3CB"),
        Line2D([0], [0], marker="^", color=C["gnina"], ls="none", ms=6, label="mTOR swap only"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=5.6, frameon=False)
    PROVENANCE["plotted"]["fig4B"] = fig4b

    ax = axes[2]
    panel_label(ax, "C", x=-0.24, y=1.06)
    seed_vals = defaultdict(list)
    for r in D["seeds"]:
        seed_vals[r["pair"]].append(fnum(r["summary_min"]))
    y = np.arange(len(PAIR_ORDER))
    plotted_s = {}
    for i, p in enumerate(PAIR_ORDER):
        vals = np.asarray(seed_vals[p], float)
        med, lo, hi = float(np.median(vals)), float(np.min(vals)), float(np.max(vals))
        ax.plot([lo, hi], [i, i], color=C["vina"], lw=1.4, zorder=3)
        ax.plot(med, i, "o", color=C["vina"], markersize=5.5, zorder=4)
        prim = fnum(D["theta6"][p]["pocket_matched_summary_min"])
        ax.plot(prim, i, "D", color=C["desc"], markersize=4.5, zorder=5)
        plotted_s[p] = {"median": med, "min": lo, "max": hi, "primary": prim, "n": len(vals)}
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"], fontsize=6.2)
    ax.set_xlabel("summary_min across 5 Vina seeds")
    ax.set_xlim(0.25, 0.82)
    ax.invert_yaxis()
    ax.legend(handles=[
        Line2D([0], [0], marker="D", color=C["desc"], ls="none", ms=5, label="primary seed"),
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5, label="median"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=5.6, frameon=False)
    PROVENANCE["plotted"]["fig4C"] = plotted_s

    fig.subplots_adjust(wspace=0.42, left=0.08, right=0.98, top=0.86, bottom=0.30)
    save_all(fig, "Fig4_computational_realization")
    plt.close(fig)


def fig5_mismatched(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.55))

    def forest(ax, pairs, set_name, title):
        ys, los, his, excl = [], [], [], []
        for p in pairs:
            r = D["delta"][(p, set_name)]
            ys.append(fnum(r["delta_matched_minus_wrong"]))
            los.append(fnum(r["delta_ci_lo"]))
            his.append(fnum(r["delta_ci_hi"]))
            excl.append(r["ci_excludes_zero"] == "True")
        y = np.arange(len(pairs))
        for i in range(len(pairs)):
            col = C["vina"] if excl[i] else "#888888"
            ax.plot([los[i], his[i]], [i, i], color=col, lw=1.5, zorder=3)
            ax.plot(ys[i], i, "o", color=col, markersize=5.5, zorder=4)
        ax.axvline(0, color=C["ink"], ls="--", lw=0.85, zorder=1)
        ax.set_yticks(y)
        ax.set_yticklabels(pairs, fontsize=6.4)
        ax.set_xlabel("ΔAUROC (matched − mismatched)")
        ax.set_title(title, fontsize=FS_AXIS, pad=3)
        ax.invert_yaxis()
        return [{"pair": pairs[i], "delta": ys[i], "lo": los[i], "hi": his[i], "excludes_zero": excl[i]}
                for i in range(len(pairs))]

    ax = axes[0]
    panel_label(ax, "A", x=-0.22, y=1.06)
    PROVENANCE["plotted"]["fig5A"] = forest(ax, PAIR_ORDER, "main_panel", "Main panels")
    ax.set_xlim(-0.22, 0.36)

    ax = axes[1]
    panel_label(ax, "B", x=-0.22, y=1.06)
    PROVENANCE["plotted"]["fig5B"] = forest(ax, HOLD_PAIRS, "unused_pool_holdout", "Unused-pool holdout")
    ax.set_xlim(-0.36, 0.18)
    ax.text(0.98, 0.04, "no EGFR/HER2 holdout", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color="#666666")

    ax = axes[2]
    panel_label(ax, "C", x=-0.22, y=1.06)
    fams = ["unmatched", "potency_matched", "size_matched"]
    fam_lab = ["unmatched", "|Δp|≤0.5", "|ΔN|≤2"]
    x = np.arange(len(HOLD_PAIRS))
    w = 0.24
    plotted_c = {}
    for j, fam in enumerate(fams):
        vals = []
        for p in HOLD_PAIRS:
            rec = D["hold_match"][(p, fam, "pocket_matched")]
            vals.append(fnum(rec["gap_matched_minus_wrong"]))
        ax.bar(x + (j - 1) * w, vals, w, color=[C["vina"], C["desc"], C["a_only"]][j],
               label=fam_lab[j], zorder=3)
        plotted_c[fam] = vals
    ax.axhline(0, color=C["ink"], lw=0.8, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(HOLD_TICK, fontsize=6.2)
    ax.set_ylabel("Δ (matched − mismatched)")
    ax.set_title("Holdout matching (point Δ)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1, fontsize=5.6, frameon=False)
    ax.set_ylim(-0.28, 0.08)
    PROVENANCE["plotted"]["fig5C"] = plotted_c

    fig.subplots_adjust(wspace=0.48, left=0.14, right=0.98, top=0.86, bottom=0.28)
    save_all(fig, "Fig5_mismatched_pocket")
    plt.close(fig)


def fig6_boundary(D: dict) -> None:
    from plot_jcim_si_composites_v1 import _pm48_e8
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.85))
    rules = ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"]
    rule_lab = ["θ=5.5", "θ=6.0", "θ=6.5", "strict"]
    colors = {"EGFR/HER2": C["egfr"], "AChE/BChE": C["rtm"], "PIK3CA/PIK3CB": C["gnina"], "PIK3CA/mTOR": C["vina"]}

    ax = axes[0, 0]
    panel_label(ax, "A", x=-0.18, y=1.04)
    x = np.arange(len(rules))
    grid = {}
    for p in PAIR_ORDER:
        ys = []
        for rule in rules:
            r = next(row for row in D["theta_all"] if row["pair"] == p and row["label_rule"] == rule)
            ys.append(fnum(r["pocket_matched_summary_min"]))
        ax.plot(x, ys, "-o", color=colors[p], lw=1.1, markersize=4.2, label=p, zorder=3)
        grid[p] = ys
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(rule_lab, fontsize=6.5)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0.15, 0.85)
    ax.set_title("Label-threshold sensitivity", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=5.5, frameon=False)
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
    ax.text(0, y48 + 0.03, f"{y48:.3f}", ha="center", fontsize=6.2)
    ax.text(1, y110 + 0.03, f"{y110:.3f}", ha="center", fontsize=6.2)
    PROVENANCE["plotted"]["fig6B"] = {"PM48": y48, "PM110": y110}

    ax = axes[1, 0]
    panel_label(ax, "C", x=-0.18, y=1.04)
    e8 = _pm48_e8()
    e16 = fnum(D["theta6"]["PIK3CA/mTOR"]["pocket_matched_summary_min"])
    ax.bar([0], [e16], 0.55, color=C["vina"], zorder=3)
    ax.bar([1], [e8["summary_min"]], 0.55, color="#56B4E9", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["E=16\n(primary)", "E=8"], fontsize=6.5)
    ax.set_ylabel("PIK3CA/mTOR summary_min")
    ax.set_ylim(0, 1.0)
    ax.set_title("Search exhaustiveness", fontsize=FS_AXIS, pad=3)
    ax.text(0, e16 + 0.03, f"{e16:.3f}", ha="center", fontsize=6.2)
    ax.text(1, e8["summary_min"] + 0.03, f"{e8['summary_min']:.3f}", ha="center", fontsize=6.2)
    PROVENANCE["plotted"]["fig6C"] = {"E16": e16, "E8": e8["summary_min"], "e8_n": (e8["nD"], e8["nA"], e8["nB"])}

    ax = axes[1, 1]
    panel_label(ax, "D", x=-0.18, y=1.04)
    four = [r for r in D["native"] if r["pair"] in PAIR_ORDER]
    n_fail = sum(r["packaged_as_external_evaluation"] == "0" for r in four)
    n_pass = 4 - n_fail
    ax.bar(["gate fail", "gate pass"], [n_fail, n_pass], color=[C["egfr"], C["thick"]], width=0.55, zorder=3)
    ax.set_ylabel("Number of pairs")
    ax.set_ylim(0, 5)
    ax.set_title("BindingDB-native external gate", fontsize=FS_AXIS, pad=3)
    ax.text(0, n_fail + 0.12, str(n_fail), ha="center", fontsize=8, fontweight="bold")
    ax.text(1, n_pass + 0.12, str(n_pass), ha="center", fontsize=8, fontweight="bold")
    ax.text(0.5, 4.55, "not docked", ha="center", fontsize=6.0, color="#666666")
    PROVENANCE["plotted"]["fig6D"] = {"n_fail": n_fail, "n_pass": n_pass}

    fig.subplots_adjust(wspace=0.38, hspace=0.58, left=0.09, right=0.98, top=0.94, bottom=0.12)
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
    ax.annotate("", xy=(19.4, 10.6), xytext=(16.4, 10.6),
                arrowprops=dict(arrowstyle="-|>", color=C["vina"], lw=1.3))
    ax.text(22.6, 12.3, "pocket-matched", ha="center", fontsize=7, fontweight="bold")
    ax.text(22.6, 10.6, "directional evaluation", ha="center", fontsize=7)
    ax.text(15.0, 4.6, "Dual-vs-Neither  ≠  Dual-vs-selective", ha="center", fontsize=7.2, fontweight="bold")
    ax.text(15.0, 2.4, "Negative-class definition changes apparent performance",
            ha="center", fontsize=6.3, color="#555555")
    save_all(fig, "TOC_graphic", toc=True)
    plt.close(fig)


def fig_s7_diagnostics(D: dict | None = None) -> None:
    """Post-hoc formulation and screening diagnostics (not upgrades)."""
    census = _read(DATA / "jcim_novelty_v0/tables/theta6_pair_census_v1.csv")
    and_rows = _read(DATA / "jcim_novelty_v0/tables/and_filter_operating_point_v1.csv")
    ligand_rows = _read(DATA / "jcim_novelty_v0/tables/ligand_only_fullmap_auroc_v1.csv")
    n_pairs = len(census)
    n_dir = sum(int(r["directional_n10"]) for r in census)
    n_form = sum(int(r["formulation_n10"]) for r in census)
    n_dock = sum(int(r["docked_in_this_paper"]) for r in census)

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.35))
    ax = axes[0]
    panel_label(ax, "A", x=-0.18, y=1.06)
    counts = [n_pairs, n_dir, n_form, n_dock]
    labels = ["audited\npairs", "directional\nn≥10", "formulation\nn≥10", "docked\nhere"]
    ax.bar(range(4), counts, color=[C["vina"], C["desc"], C["thick"], C["egfr"]], width=0.62, zorder=3)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6.0)
    ax.set_ylabel("Pair count")
    ax.set_ylim(0, max(counts) + 8)
    ax.set_title("θ=6.0 label census", fontsize=FS_AXIS, pad=3)
    for i, v in enumerate(counts):
        ax.text(i, v + 1.0, str(v), ha="center", fontsize=6.5)

    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.06)
    for pair, col in zip(PAIR_ORDER, (C["egfr"], C["rtm"], C["gnina"], C["vina"])):
        sub = [r for r in and_rows if r["pair"] == pair and r["score"] == "vina_worst"]
        rec = [fnum(r["recall_dual"]) for r in sub]
        prec = [fnum(r["precision_dual"]) for r in sub]
        ax.plot(rec, prec, marker="o", color=col, label=pair, lw=1.1, markersize=4.0, zorder=3)
    ax.set_xlabel("Dual recall")
    ax.set_ylabel("Dual precision")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_title("AND-like filter (vina_worst)", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2, fontsize=5.4, frameon=False)

    ax = axes[2]
    panel_label(ax, "C", x=-0.18, y=1.06)
    neither, directional = [], []
    for pair in PAIR_ORDER:
        n = next(r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "D_vs_neither")
        s = next(r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "summary_min_ecfp4")
        neither.append(fnum(n["ecfp4_groupkfold_auroc"]))
        directional.append(fnum(s["ecfp4_groupkfold_auroc"]))
    x = np.arange(len(PAIR_ORDER))
    ax.bar(x - 0.18, neither, 0.36, color=C["vina"], label="Dual vs neither", zorder=3)
    ax.bar(x + 0.18, directional, 0.36, color=C["egfr"], label="ECFP4 summary_min", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.0)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.85, zorder=1)
    ax.set_ylim(0.35, 1.05)
    ax.set_ylabel("GroupKFold AUROC")
    ax.set_title("Ligand-only full maps", fontsize=FS_AXIS, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1, fontsize=5.4, frameon=False)

    PROVENANCE["plotted"]["figS7"] = {
        "n_pairs": n_pairs, "n_dir": n_dir, "n_form": n_form, "n_dock": n_dock,
        "neither": neither, "directional": directional,
    }
    fig.subplots_adjust(wspace=0.42, left=0.08, right=0.98, top=0.86, bottom=0.30)
    save_all(fig, "FigS7_posthoc_diagnostics")
    plt.close(fig)


def fig_s8_bindingdb() -> None:
    rows = _read(DATA / "jcim_novelty_v0/tables/external_slice_summary_v1.csv")
    four = [next(r for r in rows if r["pair"] == p) for p in PAIR_ORDER]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.40), gridspec_kw={"width_ratios": [1.2, 1.05]})

    ax = axes[0]
    panel_label(ax, "A", x=-0.16, y=1.06)
    stages = ["native_paired", "after_literature", "after_structure", "after_ecfp_lt_0.70"]
    labs = ["paired", "−literature", "−structure", "ECFP4<0.70"]
    x = np.arange(len(stages))
    plotted = {}
    for p, col in zip(PAIR_ORDER, (C["egfr"], C["rtm"], C["gnina"], C["vina"])):
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
    x = np.arange(len(PAIR_ORDER))
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
    ax.set_xticklabels(TICK, fontsize=6.0)
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


OBSOLETE_STEMS = [
    "Fig1_task_schematic",
    "Fig2_hardneg_supply",
    "Fig3_formulation_comparison",
    "Fig4_confounds",
    "Fig5_receptor_realization",
    "Fig6_wrong_pocket_paradox",
    "Fig7_confound_anatomy",
    "Fig8_diagnostic_workflow",
    "Fig3_pocket_matched_forest",
    "Fig5_holdout_and_crystal_swap",
    "FigS1_wrong_pocket_main_vs_holdout",
    "FigS2_protocol_sensitivity",
    "FigS3_confound_anatomy",
    "FigS4_holdout_mechanism_and_supply",
    "FigS_formulation_upgrades_v1",
]


def unlink_obsolete() -> None:
    for stem in OBSOLETE_STEMS:
        for ext in ("pdf", "png", "tif"):
            p = OUT / f"{stem}.{ext}"
            if p.exists():
                p.unlink()


def draw_remaining_si() -> None:
    import plot_jcim_article_figures_v1 as v1
    from plot_jcim_si_composites_v1 import draw_all, extend_load, verify_si

    sys.path.insert(0, str(ROOT / "data" / "jcim_novelty_v0" / "scripts"))
    from plot_detectable_effect_and_workflow_v1 import fig_s_detectable

    Dv1 = v1.load()
    v1.fig3_forest(Dv1)
    v1.figS_holdout(Dv1)
    extend_load(Dv1, v1._read, v1.PROVENANCE)
    draw_all(Dv1, v1.PROVENANCE)
    errors: list[str] = []
    verify_si(Dv1, v1.PROVENANCE, errors)
    if errors:
        raise SystemExit("SI VERIFICATION FAILED:\n" + "\n".join(errors))
    fig_s7_diagnostics()
    fig_s8_bindingdb()
    fig_s_detectable()
    print("SI verification OK")


def verify(D: dict) -> None:
    errors = []

    def eq(a, b, tol=5e-4, msg=""):
        if abs(float(a) - float(b)) > tol:
            errors.append(f"{msg}: plotted {a} != source {b}")

    p1 = PROVENANCE["plotted"]["fig1C"]
    if p1["n_pairs"] != 49:
        errors.append(f"fig1C n_pairs {p1['n_pairs']} != 49")
    if p1["n_thick"] != 4:
        errors.append(f"fig1C n_thick {p1['n_thick']} != 4")
    for pair, val in {"PIK3CA/MTOR": 80, "ACHE/BCHE": 78, "PIK3CA/PIK3CB": 56, "EGFR/HER2": 7, "HDAC1/HDAC6": 93}.items():
        eq(p1["highlighted"][pair], val, msg=f"fig1C {pair}")
    eq(p1["complete_case_min"], 0.145119, tol=1e-6, msg="fig1C coverage min")
    eq(p1["complete_case_max"], 0.340172, tol=1e-6, msg="fig1C coverage max")

    t6 = D["theta6"]
    for i, p in enumerate(PAIR_ORDER):
        eq(PROVENANCE["plotted"]["fig2A"]["DA"][i], t6[p]["auroc_D_vs_A"], msg=f"fig2A DA {p}")
        eq(PROVENANCE["plotted"]["fig2A"]["DB"][i], t6[p]["auroc_D_vs_B"], msg=f"fig2A DB {p}")
    expected_dir = {
        "EGFR/HER2": (0.4297, 0.2818, 0.5775),
        "AChE/BChE": (0.6058, 0.4370, 0.7303),
        "PIK3CA/PIK3CB": (0.5000, 0.3502, 0.6495),
        "PIK3CA/mTOR": (0.6921, 0.4702, 0.8133),
    }
    expected_nei = {
        "EGFR/HER2": (0.756, 0.5625, 0.9197, 12),
        "AChE/BChE": (0.6494, 0.484, 0.8123, 15),
        "PIK3CA/PIK3CB": (0.5592, 0.3728, 0.7456, 16),
        "PIK3CA/mTOR": (0.5139, 0.2222, 0.8056, 4),
    }
    for p in PAIR_ORDER:
        rec = PROVENANCE["plotted"]["fig2B"][p]
        eq(rec["directional"]["y"], t6[p]["pocket_matched_summary_min"], msg=f"fig2B dir {p}")
        eq(rec["directional"]["lo"], expected_dir[p][1], msg=f"fig2B dir lo {p}")
        eq(rec["directional"]["hi"], expected_dir[p][2], msg=f"fig2B dir hi {p}")
        nr = D["form_by"][p]["D_vs_neither_mean"]
        eq(rec["neither"]["y"], nr["auroc"], msg=f"fig2B neither {p}")
        eq(rec["neither"]["y"], expected_nei[p][0], msg=f"fig2B neither checksum {p}")
        if rec["neither"]["n_neg"] != expected_nei[p][3]:
            errors.append(f"fig2B n_neg {p}")
    egfr = PROVENANCE["plotted"]["fig2C"]["EGFR/HER2"]
    src = D["equal"][("EGFR/HER2", "D_vs_B_or_neither_pocketA")]
    eq(egfr["delta"], src["delta_neither_minus_selective"], msg="fig2C EGFR delta vs CSV")
    eq(egfr["delta"], 0.3783, msg="fig2C EGFR 0.3783")
    eq(egfr["lo"], 0.2050, msg="fig2C EGFR lo")
    eq(egfr["hi"], 0.5469, msg="fig2C EGFR hi")

    eq(PROVENANCE["plotted"]["fig3A"]["ecfp_DB"][0], 0.8895, msg="fig3A EGFR ECFP D/B")
    eq(PROVENANCE["plotted"]["fig3A"]["vina_DB"][0], 0.4297, msg="fig3A EGFR Vina D/B")
    if PROVENANCE["plotted"]["fig3B_max_abs"] > 0.0205:
        errors.append(f"fig3B max abs {PROVENANCE['plotted']['fig3B_max_abs']} > 0.020")
    if PROVENANCE["plotted"]["fig3C"]["n"] != [27, 25, 28]:
        errors.append(f"fig3C n {PROVENANCE['plotted']['fig3C']['n']}")

    eq(PROVENANCE["plotted"]["fig4A"]["gnina_smin"][0], 0.2199, msg="fig4A EGFR GNINA smin")
    eq(PROVENANCE["plotted"]["fig4A"]["gnina_neither"][0], 0.7825, msg="fig4A EGFR GNINA neither")
    b = {(r["pair"], r["crystal"]): r for r in PROVENANCE["plotted"]["fig4B"]}
    eq(b[("PIK3CA/mTOR", "4JPS")]["y"], 0.4861, msg="fig4B PM 4JPS")
    eq(b[("PIK3CA/mTOR", "5DXT")]["y"], 0.5046, msg="fig4B PM 5DXT")
    eq(b[("PIK3CA/PIK3CB", "4JPS")]["y"], 0.6905, msg="fig4B PAB 4JPS")
    eq(b[("PIK3CA/PIK3CB", "5DXT")]["y"], 0.6849, msg="fig4B PAB 5DXT")
    eq(b[("PIK3CA/mTOR", "4JSX")]["y"], 0.6389, msg="fig4B 4JSX")
    if ("PIK3CA/PIK3CB", "4JSX") in b:
        errors.append("fig4B must not put 4JSX on PIK3CA/PIK3CB")
    eq(PROVENANCE["plotted"]["fig4C"]["EGFR/HER2"]["primary"], 0.4297, msg="fig4C EGFR primary")
    if PROVENANCE["plotted"]["fig4C"]["EGFR/HER2"]["n"] != 5:
        errors.append("fig4C EGFR n_seeds")

    a = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5A"]}
    eq(a["EGFR/HER2"]["delta"], 0.1697, msg="fig5A EGFR delta")
    eq(a["AChE/BChE"]["delta"], 0.1614, msg="fig5A AChE delta")
    if not a["EGFR/HER2"]["excludes_zero"] or a["PIK3CA/mTOR"]["excludes_zero"]:
        errors.append("fig5A CI exclude-zero pattern")
    bhold = {r["pair"]: r for r in PROVENANCE["plotted"]["fig5B"]}
    eq(bhold["PIK3CA/mTOR"]["delta"], -0.0225, msg="fig5B PM delta")
    if any(r["excludes_zero"] for r in PROVENANCE["plotted"]["fig5B"]):
        errors.append("fig5B all CIs should include 0")
    eq(PROVENANCE["plotted"]["fig5C"]["unmatched"][2], -0.0225, msg="fig5C unmatched PM")

    eq(PROVENANCE["plotted"]["fig6A"]["PIK3CA/mTOR"][1], 0.6921, msg="fig6A PM theta6")
    eq(PROVENANCE["plotted"]["fig6B"]["PM48"], 0.6921, msg="fig6B PM48")
    eq(PROVENANCE["plotted"]["fig6B"]["PM110"], 0.6483, msg="fig6B PM110")
    eq(PROVENANCE["plotted"]["fig6C"]["E16"], 0.6921, msg="fig6C E16")
    eq(PROVENANCE["plotted"]["fig6C"]["E8"], 0.6597, msg="fig6C E8")
    if PROVENANCE["plotted"]["fig6D"]["n_pass"] != 0 or PROVENANCE["plotted"]["fig6D"]["n_fail"] != 4:
        errors.append("fig6D BindingDB gate")
    s7 = PROVENANCE["plotted"].get("figS7")
    if s7:
        if s7["n_dir"] != 17 or s7["n_dock"] != 4:
            errors.append(f"figS7 census {s7}")
        eq(s7["neither"][0], 0.9214, msg="figS7 EGFR Dual vs neither")
        eq(s7["directional"][0], 0.8013, msg="figS7 EGFR ECFP summary_min")
    s8 = PROVENANCE["plotted"].get("figS8")
    if s8:
        if s8["n_fail"] != 4:
            errors.append("figS8 gate")
        eq(s8["after_ecfp"][0], 216, msg="figS8 EGFR after ECFP")

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
    print("v2 verification OK")


def write_captions() -> None:
    text = """# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match `figures/jcim_article/`.
All numbers are read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V2.md`.
Regenerate: `python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v2.py`

## Figure 1. Four-state dual-target evaluation and data supply.

(A) Four experimentally labeled ligand states: dual, A-only, B-only, and neither. A-only and B-only are selectivity hard negatives. (B) Primary tasks are pocket-matched directional AUROCs: Dual versus A-only scored in pocket B, Dual versus B-only scored in pocket A. `summary_min` is a descriptive worst-arm summary. (C) J0 ChEMBL audit of 49 candidate pairs (`j0_strict_label_supply.csv`). Four pairs meet a thick hard-negative gate (min ≥50); HDAC1/HDAC6 is excluded as metal-dependent; EGFR/HER2 is retained as a supply-limited K=4 case (strict B-only = 7). Complete-case map coverage is 14.5%–34.0% (`complete_case_usable_pchembl_overlap_v1.csv`). Cross-database counts are Figure S2.

## Figure 2. Negative-class definition changes apparent dual-target evidence.

Same frozen AutoDock Vina scores, unified θ = 6.0. (A) Directional Dual versus A-only (pocket B) and Dual versus B-only (pocket A) from `unified_threshold_sensitivity_v2.csv`. (B) Descriptive comparison of directional `summary_min` (Table 2 CIs from the same unified-threshold file) with Dual versus neither using per-ligand `vina_mean` (`formulation_conventional_vs_directional_v1.csv`). These two columns differ in both negative class and score aggregation; the difference is not a paired test of one estimand. PIK3CA/mTOR Dual versus neither is hatched (neither n = 4). (C) Pocket A score held fixed; only the negative class is replaced (B-only versus neither). EGFR/HER2 ΔAUROC = 0.378 [0.205, 0.547] (`formulation_equal_score_negative_v1.csv`). Vertical dashed line in (C), zero.

## Figure 3. Ligand chemistry as a competing explanation.

(A) Scaffold GroupKFold ECFP4 logistic AUROC versus pocket-matched Vina rank AUROC on both directional arms (`ligand_ml_baseline_scaffold_cv_v1.csv`). EGFR/HER2 Dual versus B-only: ECFP4 0.8895 versus Vina 0.4297. (B) Change in GroupKFold AUROC when the pocket-matched Vina score is added to ECFP4 (`incremental_information_v1.csv`). The largest absolute change among the eight contrasts is ≤0.020. (C) AChE/BChE TPSA by class: individual ligands (jittered) with median and IQR (`assembled_AChE_BChE.csv`). n = 27/25/28.

## Figure 4. Computational realization.

(A) Independent GNINA 1.3.2 pose generation (not CNN rescoring of Vina poses) on EGFR/HER2 and PIK3CA/mTOR (`independent_dock_formulation_v1.csv`) versus the same-panel Vina values. (B) Replacing PIK3CA 4L23 with 4JPS or 5DXT while holding the second pocket frozen: PIK3CA/mTOR `summary_min` 0.692 → 0.486 / 0.505; PIK3CA/PIK3CB 0.500 → 0.691 / 0.685. 4JSX is an mTOR-pocket swap and is plotted with a distinct marker; it is not applied to PIK3CA/PIK3CB. CIs from the deposited swap tables and Table 2. (C) Directional `summary_min` across five frozen Vina seeds (`multiseed_auroc_by_seed_v2.csv`); diamond, production seed 20260727.

## Figure 5. Matched- versus mismatched-pocket scoring controls.

Δ = matched-pocket `summary_min` − mismatched-pocket `summary_min`, ligand bootstrap B = 2000 (`wrong_pocket_paired_delta_bootstrap_v1.csv`). Matched uses Dual versus A-only in pocket B and Dual versus B-only in pocket A; mismatched swaps those score channels. This is a scoring-channel control, not redocking into a physically wrong site. (A) Main panels. EGFR/HER2 and AChE/BChE CIs exclude 0; PIK3CA pairs include 0. (B) Unused-pool holdout. All three CIs include 0; point estimates are negative (mismatched ≥ matched). EGFR/HER2 has no holdout. (C) Holdout potency (|ΔpChEMBL| ≤ 0.5) and size (|Δheavy| ≤ 2) matching: point Δ only (`holdout_matched_wrong_pocket_summary_v1.csv`). Dark, CI excludes 0; gray, CI includes 0.

## Figure 6. Robustness checks and evidence boundary.

(A) Pocket-matched `summary_min` on the unified label-threshold grid (`unified_threshold_sensitivity_v2.csv`). (B) PIK3CA/mTOR PM48 versus PM110 Vina (`pm110_vs_pm48_pocket_matched_v1.csv`). (C) PM48 exhaustiveness 16 versus 8, recomputed from `scores_vina_E8_best.csv` with the same pocket-matched definition. (D) BindingDB-native 202608 slice: zero of four pairs meet the pre-frozen external gate; nothing was docked (`external_slice_summary_v1.csv`).

## Figure S1. Protocol and panel sensitivities.

(A) Unified label-threshold grid. Open markers, underpowered cells. (B) GNINA CNN rescoring of Vina poses (mode-1 versus best-of-9) versus primary Vina; this is not independent GNINA pose generation (Figure 4A / Table S32). (C) PM48 versus PM110 for Vina, RTMScore, and GNINA CNN best-of-9 rescoring. (D) Exhaustiveness 16 versus 8 and single-target enrichment on 4L23/4JT6.

## Figure S2. Equal-relation supply and holdout sampling shift.

Unchanged sources: `crossdb_strict_supply_v1.csv`; `holdout_vs_main_potency_size_v1.csv`.

## Figure S3. Additional paired bootstrap differences.

Descriptor and scaffold-versus-random leakage checks. Matched-versus-mismatched main/holdout Δ CIs are now Figure 5.

## Figure S4. Pocket-matched summary_min forest (former main figure).

Vina CIs from `unified_threshold_sensitivity_v2.csv`. GNINA in this figure is CNN rescoring of Vina poses (best-of-9), not independent pose generation.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched `summary_min` only; mismatched-pocket Δ CIs are Figure 5B.

## Figure S7. Post-hoc formulation and screening diagnostics.

θ = 6.0 pair census, AND-like dual filter, and full-map ligand-only ECFP4. Not docking upgrades and not a replacement for Table 2.

## Figure S8. BindingDB-native slice.

Filter cascade and remaining four-state counts after literature, structure, and ECFP4 < 0.70 (`external_slice_summary_v1.csv`). Zero of four pairs meet the pre-frozen external gate; nothing was docked.

## Figure S9. Additional ligand-structure controls.

Prespecified descriptors, covariate-adjusted logistic AUROC, and matched-subset weak-arm tests (former main Figure 7B–D). The Vina-only logistic AUROC is not the Table 2 rank AUROC.

## Figure S10. Matched versus mismatched point estimates.

Bar charts of matched versus mismatched `summary_min` on the main panel, unused-pool holdout, potency/size matching, and scoring-free contact counts. Paired Δ CIs are Figure 5. Contact count is exploratory and does not explain holdout reversal.

## TOC graphic (For Table of Contents Only).

Four experimental states, pocket-matched directional evaluation, and the qualitative statement that Dual-versus-neither is not Dual-versus-selective. No numerical AUROCs.
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
    draw_remaining_si()
    write_captions()
    unlink_obsolete()
    (OUT / "plotted_values.json").write_text(
        json.dumps(PROVENANCE, indent=2, default=str), encoding="utf-8"
    )
    verify(D)
    print("wrote", OUT)
    for p in sorted(OUT.glob("*")):
        print(" ", p.name, p.stat().st_size)


if __name__ == "__main__":
    main()
