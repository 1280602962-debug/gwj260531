#!/usr/bin/env python3
"""Draw DualFourClass-Bench JCIM main figures from frozen CSV sources only.

Every plotted number is read from a named source file. No hand-typed AUROCs.
Run: python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jcim_figure_style import (  # noqa: E402
    C,
    DESC_LABEL,
    FS_ANNO,
    FS_AXIS,
    OUT,
    PAIR_ORDER,
    ROOT,
    apply_style,
    panel_label,
    save_all,
)

DATA = ROOT / "data"
PROVENANCE: dict = {"source_files": {}, "plotted": {}}


def _read(path: Path) -> list[dict]:
    rows = list(csv.DictReader(path.open()))
    PROVENANCE["source_files"].setdefault(str(path.relative_to(ROOT)), len(rows))
    return rows


def fnum(x: str) -> float:
    return float(x)


def load() -> dict:
    j0 = _read(DATA / "jcim_j0j1_v0/tables/j0_strict_label_supply.csv")
    s12 = _read(DATA / "jcim_supply_crossdb_v0/tables/crossdb_strict_supply_v1.csv")
    theta = _read(DATA / "jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv")
    pm = _read(DATA / "jcim_bench_v0/tables/pocket_matched_directional_v1.csv")
    forest = _read(DATA / "jcim_bench_v0/tables/forest_summary_min_ci_v1.csv")
    gnina = _read(DATA / "jcim_bench_v0/tables/gnina_pocket_matched_mode01_vs_best9_k4_v1.csv")
    hold = _read(DATA / "jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv")
    jps = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PM48_alt4JPS_v1.csv")[0]
    dxt = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PM48_alt5DXT_v1.csv")[0]
    jsx = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PM48_alt4JSX_v1.csv")[0]
    pab_jps = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PAB_alt4JPS_v1.csv")[0]
    pab_dxt = _read(DATA / "jcim_structure_robust_v0/tables/pocket_matched_PAB_alt5DXT_v1.csv")[0]
    ache = _read(DATA / "jcim_bench_v0/tables/assembled_AChE_BChE.csv")
    form = _read(DATA / "jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv")

    theta6 = {r["pair"]: r for r in theta if r["label_rule"] == "theta_6.0"}
    form_by: dict = {}
    for r in form:
        form_by.setdefault(r["pair"], {})[r["contrast"]] = r
    pm_by = {(r["pair"], r["variant"]): r for r in pm}
    forest_by = {(r["pair"], r["arm"]): r for r in forest}
    gnina9 = {r["pair"]: r for r in gnina if r["channel"] == "best9"}
    hold_pm = {r["pair"]: r for r in hold if r["variant"] == "pocket_matched_vina"}
    hold_wp = {r["pair"]: r for r in hold if r["variant"] == "wrong_pocket_control_vina"}

    # Best single-descriptor reference = highest summary_min among the four
    # prespecified properties (descriptive, not a confirmatory competitor).
    best_desc = {}
    for pair in PAIR_ORDER:
        cands = []
        for arm in ("heavy", "mw", "clogp", "tpsa"):
            r = forest_by[(pair, arm)]
            cands.append((fnum(r["summary_min"]), arm, r))
        cands.sort(reverse=True)
        best_desc[pair] = cands[0]

    tpsa = defaultdict(list)
    for r in ache:
        if r["cls"] in ("dual", "A_only", "B_only"):
            tpsa[r["cls"]].append(fnum(r["tpsa"]))

    return {
        "j0": j0,
        "s12": s12,
        "theta6": theta6,
        "pm_by": pm_by,
        "forest_by": forest_by,
        "gnina9": gnina9,
        "hold_pm": hold_pm,
        "hold_wp": hold_wp,
        "jps": jps,
        "dxt": dxt,
        "jsx": jsx,
        "pab_jps": pab_jps,
        "pab_dxt": pab_dxt,
        "best_desc": best_desc,
        "tpsa": tpsa,
        "theta_all": theta,
        "form_by": form_by,
    }


def fig1_task() -> None:
    """Schematic — no experimental numbers."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.15))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis("off")

    ax = axes[0]
    panel_label(ax, "A", x=-0.02, y=1.02)
    # two pockets
    for x, title, col in ((1.3, "Pocket A", C["a_only"]), (5.7, "Pocket B", C["b_only"])):
        ax.add_patch(
            FancyBboxPatch(
                (x, 5.6), 3.0, 3.6,
                boxstyle="round,pad=0.08,rounding_size=0.35",
                facecolor="#F4F7FA", edgecolor=col, linewidth=1.4, clip_on=False,
            )
        )
        ax.text(x + 1.5, 8.7, title, ha="center", va="center", fontsize=FS_AXIS, fontweight="bold")

    # four ligand classes as two-site occupancy
    classes = [
        (1.6, "dual", C["dual"], C["dual"], True, True),
        (3.7, "A_only", C["a_only"], "#DDDDDD", True, False),
        (5.8, "B_only", "#DDDDDD", C["b_only"], False, True),
        (7.9, "neither", "#DDDDDD", "#DDDDDD", False, False),
    ]
    ax.text(5.0, 4.85, "Four ligand classes", ha="center", fontsize=FS_AXIS, fontweight="bold")
    ax.text(5.0, 4.35, "A_only / B_only = experimental hard negatives", ha="center", fontsize=FS_ANNO, color="#555555")
    for x, name, c1, c2, _, _ in classes:
        ax.add_patch(Circle((x - 0.28, 2.55), 0.38, facecolor=c1, edgecolor=C["ink"], lw=0.6, clip_on=False))
        ax.add_patch(Circle((x + 0.28, 2.55), 0.38, facecolor=c2, edgecolor=C["ink"], lw=0.6, clip_on=False))
        ax.text(x, 1.55, name, ha="center", va="top", fontsize=FS_ANNO)

    ax = axes[1]
    panel_label(ax, "B", x=-0.02, y=1.02)
    ax.text(5.0, 9.35, "Pocket-matched readout", ha="center", fontsize=FS_AXIS, fontweight="bold")

    boxes = [
        (0.5, 5.7, "dual vs A_only", "score in pocket B", "A_only is already potent at A"),
        (0.5, 2.55, "dual vs B_only", "score in pocket A", "B_only is already potent at B"),
    ]
    for x, y, head, mid, foot in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y), 9.0, 2.55,
                boxstyle="round,pad=0.06,rounding_size=0.25",
                facecolor="#F7F7F7", edgecolor="#CCCCCC", lw=0.8, clip_on=False,
            )
        )
        ax.text(x + 0.35, y + 1.85, head, ha="left", fontsize=FS_AXIS, fontweight="bold")
        ax.annotate(
            "",
            xy=(x + 8.4, y + 1.15),
            xytext=(x + 4.15, y + 1.15),
            arrowprops=dict(arrowstyle="-|>", color=C["vina"], lw=1.3),
        )
        ax.text(x + 0.35, y + 1.05, mid, ha="left", va="center", fontsize=FS_ANNO, color=C["vina"])
        ax.text(x + 0.35, y + 0.4, foot, ha="left", fontsize=6.5, color="#666666")

    ax.text(
        5.0, 1.55,
        r"summary_min  =  min(AUROC$_{\mathrm{D/A}}$, AUROC$_{\mathrm{D/B}}$)",
        ha="center", fontsize=FS_AXIS,
    )
    ax.text(
        5.0, 0.7,
        "Pooled mean of the two pockets can hide the weaker arm.",
        ha="center", fontsize=FS_ANNO, color="#555555",
    )

    fig.subplots_adjust(wspace=0.08, left=0.04, right=0.98, top=0.90, bottom=0.06)
    save_all(fig, "Fig1_task_schematic")
    plt.close(fig)


def fig2_supply(D: dict) -> None:
    rows = []
    for r in D["j0"]:
        raw = (r.get("min_strict_hardneg") or "").strip()
        if raw == "":
            continue
        rows.append(r)
    rows.sort(key=lambda r: fnum(r["min_strict_hardneg"]), reverse=True)

    def pretty(name: str) -> str:
        return (
            name.replace("ACHE/BCHE", "AChE/BChE")
            .replace("MTOR", "mTOR")
            .replace("BCL2L1_BclxL", "BCL2L1")
            .replace("RPS6KB1_p70S6K", "p70S6K")
        )

    highlight = {
        "PIK3CA/MTOR": ("thick", C["thick"]),
        "ACHE/BCHE": ("thick", C["thick"]),
        "PIK3CA/PIK3CB": ("thick", C["thick"]),
        "EGFR/HER2": ("supply-limited", C["egfr"]),
        "HDAC1/HDAC6": ("metal, excluded", C["metal"]),
    }

    fig, axes = plt.subplots(
        1, 2, figsize=(7.0, 6.15), gridspec_kw={"width_ratios": [1.55, 1.0]}
    )
    ax = axes[0]
    panel_label(ax, "A", x=-0.10, y=1.03)
    y = np.arange(len(rows))
    vals = np.array([fnum(r["min_strict_hardneg"]) for r in rows], dtype=float)
    colors = []
    for r in rows:
        if r["pair"] in highlight:
            colors.append(highlight[r["pair"]][1])
        elif r["metal_enzyme_risk"] == "True":
            colors.append(C["metal"])
        else:
            colors.append(C["other"])
    ax.barh(y, vals, height=0.78, color=colors, edgecolor="none", zorder=3)
    ax.axvline(50, color=C["ink"], ls="--", lw=0.9, zorder=2)
    ax.axvline(20, color="#888888", ls=":", lw=0.9, zorder=2)
    ax.text(50, 1.015, "≥50", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5, clip_on=False)
    ax.text(20, 1.015, "≥20", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5, color="#666666", clip_on=False)
    ax.set_yticks([])
    ax.set_xlabel("min(strict A_only, strict B_only)")
    ax.set_xlim(0, 118)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.invert_yaxis()
    # annotations only for highlighted pairs (to the right of each bar)
    for i, r in enumerate(rows):
        if r["pair"] not in highlight:
            continue
        v = fnum(r["min_strict_hardneg"])
        label = f"{pretty(r['pair'])}  {int(v)}"
        ax.text(min(v + 1.8, 100), i, label, va="center", ha="left", fontsize=6.5, color=C["ink"], clip_on=False)

    handles = [
        mpatches.Patch(color=C["thick"], label="Thick supply (≥50)"),
        mpatches.Patch(color=C["egfr"], label="EGFR/HER2 (supply-limited)"),
        mpatches.Patch(color=C["metal"], label="Metal enzyme"),
        mpatches.Patch(color=C["other"], label="Other audited pairs"),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        fontsize=6.5,
        bbox_to_anchor=(0.99, 0.0),
        borderaxespad=0.0,
        frameon=False,
    )

    PROVENANCE["plotted"]["fig2A_n_pairs"] = len(rows)
    PROVENANCE["plotted"]["fig2A_highlighted"] = {
        r["pair"]: fnum(r["min_strict_hardneg"]) for r in rows if r["pair"] in highlight
    }

    # Panel B: S12 count-level, equal_only
    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.02)
    # map S12 names
    want = [
        ("PIK3CA/MTOR", "PIK3CA/mTOR"),
        ("ACHE/BCHE", "AChE/BChE"),
        ("PIK3CA/PIK3CB", "PIK3CA/PIK3CB"),
        ("EGFR/HER2", "EGFR/HER2"),
    ]
    chem, bdb = [], []
    for src_name, _ in want:
        chem.append(
            fnum(
                next(
                    r["min_strict_hardneg"]
                    for r in D["s12"]
                    if r["pair"] == src_name and r["source"] == "ChEMBL_cache" and r["rule"] == "pChEMBL"
                )
            )
        )
        bdb.append(
            fnum(
                next(
                    r["min_strict_hardneg"]
                    for r in D["s12"]
                    if r["pair"] == src_name and r["source"] == "BindingDB" and r["rule"] == "equal_only"
                )
            )
        )
    x = np.arange(len(want))
    w = 0.36
    ax.bar(x - w / 2, chem, w, color=C["vina"], label="ChEMBL pChEMBL", zorder=3)
    ax.bar(x + w / 2, bdb, w, color=C["desc"], label="BindingDB equal_only", zorder=3)
    ax.axhline(50, color=C["ink"], ls="--", lw=0.9, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([p[1] for p in want], rotation=25, ha="right")
    ax.set_ylabel("min strict hard-negatives")
    ax.set_ylim(0, 110)
    ax.legend(loc="upper right", fontsize=6.5, frameon=False)
    ax.text(0.02, 0.02, "Count-level; no docking", transform=ax.transAxes, fontsize=6.5, color="#666666", va="bottom")
    PROVENANCE["plotted"]["fig2B_chembl"] = chem
    PROVENANCE["plotted"]["fig2B_bindingdb_equal_only"] = bdb

    fig.subplots_adjust(wspace=0.32, left=0.06, right=0.98, top=0.90, bottom=0.12)
    save_all(fig, "Fig2_hardneg_supply")
    plt.close(fig)


def fig3_forest(D: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.70))

    y_pair_gap = 1.15
    offsets = {"vina": 0.33, "rtm": 0.11, "gnina": -0.11, "desc": -0.33}
    markers = {"vina": "o", "rtm": "^", "gnina": "D", "desc": "s"}
    sizes = {"vina": 7.0, "rtm": 5.5, "gnina": 5.0, "desc": 5.5}
    lws = {"vina": 1.6, "rtm": 1.1, "gnina": 1.1, "desc": 1.1}

    plotted = {}
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    for i, pair in enumerate(PAIR_ORDER):
        y0 = -i * y_pair_gap
        if i % 2 == 0:
            ax.axhspan(y0 - 0.55, y0 + 0.55, color="#F4F7FA", zorder=0)
        t = D["theta6"][pair]
        vina = {
            "y": fnum(t["pocket_matched_summary_min"]),
            "lo": fnum(t["ci_lo"]),
            "hi": fnum(t["ci_hi"]),
            "da": fnum(t["auroc_D_vs_A"]),
            "db": fnum(t["auroc_D_vs_B"]),
        }
        rtm_r = D["pm_by"][(pair, "pocket_matched_rtm")]
        rtm = {
            "y": fnum(rtm_r["summary_min"]),
            "lo": fnum(rtm_r["summary_min_ci_lo"]),
            "hi": fnum(rtm_r["summary_min_ci_hi"]),
        }
        g = D["gnina9"][pair]
        gn = {
            "y": fnum(g["summary_min"]),
            "lo": fnum(g["ci_lo"]),
            "hi": fnum(g["ci_hi"]),
        }
        arm, dr = D["best_desc"][pair][1], D["best_desc"][pair][2]
        desc = {
            "arm": arm,
            "y": fnum(dr["summary_min"]),
            "lo": fnum(dr["ci_lo"]),
            "hi": fnum(dr["ci_hi"]),
        }
        plotted[pair] = {"vina": vina, "rtm": rtm, "gnina_best9": gn, "best_desc": desc}

        series = [
            ("vina", vina, C["vina"]),
            ("rtm", rtm, C["rtm"]),
            ("gnina", gn, C["gnina"]),
            ("desc", desc, C["desc"]),
        ]
        for key, d, col in series:
            yy = y0 + offsets[key]
            ax.errorbar(
                d["y"], yy,
                xerr=[[d["y"] - d["lo"]], [d["hi"] - d["y"]]],
                fmt=markers[key],
                color=col, ecolor=col,
                elinewidth=lws[key], capsize=2.2, capthick=0.8,
                markersize=sizes[key], markeredgecolor=col, markerfacecolor=col,
                zorder=4,
            )
        ax.text(
            1.02,
            y0,
            DESC_LABEL[arm],
            transform=trans,
            va="center",
            ha="left",
            fontsize=6.5,
            color=C["desc"],
            clip_on=False,
        )

    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_yticks([-i * y_pair_gap for i in range(len(PAIR_ORDER))])
    ax.set_yticklabels(PAIR_ORDER)
    ax.set_xlabel("Pocket-matched summary_min AUROC (95% ligand bootstrap CI)")
    ax.set_xlim(0.08, 0.95)
    ax.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_ylim(-3 * y_pair_gap - 0.70, 0.72)
    ax.text(1.02, 0.72, "best desc.", transform=trans, ha="left", va="bottom", fontsize=6.5, color=C["desc"], clip_on=False)

    handles = [
        Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=7, label="Vina (primary)"),
        Line2D([0], [0], marker="^", color=C["rtm"], ls="none", ms=6, label="RTMScore"),
        Line2D([0], [0], marker="D", color=C["gnina"], ls="none", ms=5.5, label="GNINA best-of-9"),
        Line2D([0], [0], marker="s", color=C["desc"], ls="none", ms=6, label="Best descriptor reference"),
    ]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.06),
        ncol=4,
        frameon=False,
        fontsize=6.5,
        borderaxespad=0.0,
        handletextpad=0.35,
        columnspacing=0.9,
    )

    PROVENANCE["plotted"]["figS_forest"] = plotted
    fig.subplots_adjust(left=0.20, right=0.84, top=0.82, bottom=0.14)
    save_all(fig, "FigS_pocket_matched_forest")
    plt.close(fig)


def fig3_formulation(D: dict) -> None:
    """Main Fig 3: Dual-versus-neither comparator vs directional summary_min.

    Directional CIs come from unified θ = 6.0 (summary_min rows in the
    formulation CSV have empty CIs). Dual-versus-neither CIs come from
    formulation_conventional_vs_directional_v1.csv. Dual versus all
    non-duals remains Table 3, not this figure.
    """
    fig, ax = plt.subplots(figsize=(7.0, 3.55))
    x = np.arange(len(PAIR_ORDER))
    width = 0.36
    plotted = {}
    dir_y, dir_lo, dir_hi = [], [], []
    nei_y, nei_lo, nei_hi, n_neg = [], [], [], []
    for p in PAIR_ORDER:
        t = D["theta6"][p]
        r = D["form_by"][p]["D_vs_neither_mean"]
        dy, dlo, dhi = fnum(t["pocket_matched_summary_min"]), fnum(t["ci_lo"]), fnum(t["ci_hi"])
        ny, nlo, nhi = fnum(r["auroc"]), fnum(r["ci_lo"]), fnum(r["ci_hi"])
        nn = int(r["n_neg"])
        dir_y.append(dy)
        dir_lo.append(dlo)
        dir_hi.append(dhi)
        nei_y.append(ny)
        nei_lo.append(nlo)
        nei_hi.append(nhi)
        n_neg.append(nn)
        plotted[p] = {
            "directional": {"y": dy, "lo": dlo, "hi": dhi},
            "neither": {"y": ny, "lo": nlo, "hi": nhi, "n_neg": nn},
        }
    dir_y = np.asarray(dir_y)
    nei_y = np.asarray(nei_y)
    dir_err = np.vstack([dir_y - np.asarray(dir_lo), np.asarray(dir_hi) - dir_y])
    nei_err = np.vstack([nei_y - np.asarray(nei_lo), np.asarray(nei_hi) - nei_y])
    ax.bar(
        x - width / 2, dir_y, width, yerr=dir_err, capsize=2.4,
        color=C["vina"], ecolor=C["vina"], error_kw={"elinewidth": 1.1},
        label="directional summary_min", zorder=3,
    )
    bars_n = ax.bar(
        x + width / 2, nei_y, width, yerr=nei_err, capsize=2.4,
        color=C["a_only"], ecolor=C["a_only"], error_kw={"elinewidth": 1.1},
        label="Dual vs neither comparator", zorder=3,
    )
    bars_n[3].set_hatch("///")
    bars_n[3].set_edgecolor(C["a_only"])
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"], fontsize=7)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.12, 1.05)
    ax.legend(loc="upper center", bbox_to_anchor=(0.64, 1.02), ncol=2, fontsize=6.5, frameon=False)
    ax.annotate(
        "neither n=4\nunderpowered",
        xy=(3.0 + width / 2, nei_y[3]), xytext=(3.08, 0.89),
        ha="center", va="bottom", fontsize=6.5, color=C["a_only"],
        arrowprops=dict(arrowstyle="-|>", color=C["a_only"], lw=0.8),
    )
    ax.annotate(
        "formulation gap",
        xy=(width / 2, 0.76), xytext=(-0.55, 0.92),
        fontsize=6.5, color=C["ink"],
        arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=0.8),
    )
    PROVENANCE["plotted"]["fig3"] = plotted
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.16)
    save_all(fig, "Fig3_formulation_comparison")
    plt.close(fig)


def fig4_confounds(D: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.55), gridspec_kw={"width_ratios": [1.15, 1.15, 0.95]})

    ax = axes[0]
    panel_label(ax, "A", x=-0.18, y=1.08)
    x = np.arange(len(PAIR_ORDER))
    w = 0.36
    da = [fnum(D["theta6"][p]["auroc_D_vs_A"]) for p in PAIR_ORDER]
    db = [fnum(D["theta6"][p]["auroc_D_vs_B"]) for p in PAIR_ORDER]
    ax.bar(x - w / 2, da, w, color=C["vina"], label="D vs A_only (pocket B)", zorder=3)
    ax.bar(x + w / 2, db, w, color=C["a_only"], label="D vs B_only (pocket A)", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"], fontsize=6.5)
    ax.set_ylabel("Pocket-matched AUROC")
    ax.set_ylim(0, 1.0)
    ax.legend(
        loc="lower center",
        fontsize=6.0,
        ncol=1,
        bbox_to_anchor=(0.5, 1.02),
        frameon=False,
        borderaxespad=0.0,
    )
    PROVENANCE["plotted"]["fig4A_DA"] = da
    PROVENANCE["plotted"]["fig4A_DB"] = db

    ax = axes[1]
    panel_label(ax, "B", x=-0.18, y=1.08)
    for i, pair in enumerate(PAIR_ORDER):
        t = D["theta6"][pair]
        vy, vlo, vhi = fnum(t["pocket_matched_summary_min"]), fnum(t["ci_lo"]), fnum(t["ci_hi"])
        _, arm, dr = D["best_desc"][pair]
        dy, dlo, dhi = fnum(dr["summary_min"]), fnum(dr["ci_lo"]), fnum(dr["ci_hi"])
        ax.errorbar(
            i - 0.14, vy, yerr=[[vy - vlo], [vhi - vy]],
            fmt="o", color=C["vina"], ecolor=C["vina"], elinewidth=1.3, capsize=2.2, markersize=6, zorder=4,
        )
        ax.errorbar(
            i + 0.14, dy, yerr=[[dy - dlo], [dhi - dy]],
            fmt="s", color=C["desc"], ecolor=C["desc"], elinewidth=1.3, capsize=2.2, markersize=5.5, zorder=4,
        )
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"], fontsize=6.5)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0.08, 1.0)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="Vina pocket-matched"),
            Line2D([0], [0], marker="s", color=C["desc"], ls="none", ms=6, label="Best descriptor"),
        ],
        loc="lower center",
        fontsize=6.0,
        ncol=1,
        bbox_to_anchor=(0.5, 1.02),
        frameon=False,
        borderaxespad=0.0,
    )

    ax = axes[2]
    panel_label(ax, "C", x=-0.22, y=1.08)
    data = [D["tpsa"]["dual"], D["tpsa"]["A_only"], D["tpsa"]["B_only"]]
    parts = ax.violinplot(data, positions=[1, 2, 3], widths=0.7, showmeans=False, showmedians=True, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor([C["dual"], C["a_only"], C["b_only"]][i])
        body.set_edgecolor(C["ink"])
        body.set_alpha(0.75)
        body.set_linewidth(0.6)
    parts["cmedians"].set_color(C["ink"])
    parts["cmedians"].set_linewidth(0.9)
    ns = [len(d) for d in data]
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([f"dual\nn={ns[0]}", f"A_only\nn={ns[1]}", f"B_only\nn={ns[2]}"])
    ax.set_ylabel(r"TPSA ($\mathrm{\AA}^2$)")
    ax.set_title("AChE/BChE", fontsize=FS_AXIS, pad=3)
    means = [float(np.mean(d)) for d in data]
    PROVENANCE["plotted"]["fig4C_tpsa_n"] = ns
    PROVENANCE["plotted"]["fig4C_tpsa_mean"] = means

    fig.subplots_adjust(wspace=0.42, left=0.08, right=0.98, top=0.78, bottom=0.20)
    save_all(fig, "Fig4_confounds")
    plt.close(fig)


def figS_holdout(D: dict) -> None:
    """Former Fig 5A: unused-pool holdout. Now SI; main Fig 5 is receptor-only."""
    fig, ax = plt.subplots(figsize=(7.0, 3.20))
    pairs = ["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
    x = np.arange(len(pairs))
    fig_s = []
    for i, pair in enumerate(pairs):
        t = D["theta6"][pair]
        h = D["hold_pm"][pair]
        my, mlo, mhi = fnum(t["pocket_matched_summary_min"]), fnum(t["ci_lo"]), fnum(t["ci_hi"])
        hy, hlo, hhi = fnum(h["summary_min"]), fnum(h["summary_min_ci_lo"]), fnum(h["summary_min_ci_hi"])
        fig_s.append(
            {"pair": pair, "main": my, "main_lo": mlo, "main_hi": mhi, "hold": hy, "hold_lo": hlo, "hold_hi": hhi}
        )
        ax.errorbar(
            i - 0.14, my, yerr=[[my - mlo], [mhi - my]],
            fmt="o", color=C["main"], ecolor=C["main"], elinewidth=1.4, capsize=2.4, markersize=6.5, zorder=4,
        )
        ax.errorbar(
            i + 0.14, hy, yerr=[[hy - hlo], [hhi - hy]],
            fmt="s", color=C["holdout"], ecolor=C["holdout"], elinewidth=1.4, capsize=2.4, markersize=6, zorder=4,
        )
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"], fontsize=7)
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0.12, 1.02)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=C["main"], ls="none", ms=6, label="Main panel"),
            Line2D([0], [0], marker="s", color=C["holdout"], ls="none", ms=6, label="Unused-pool holdout"),
        ],
        loc="upper left",
        fontsize=6.5,
        frameon=False,
    )
    ax.text(0.02, 0.03, "EGFR/HER2 has no holdout", transform=ax.transAxes, ha="left", fontsize=6.0, color="#666666")
    PROVENANCE["plotted"]["figS_holdout"] = fig_s
    fig.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.16)
    save_all(fig, "FigS_unused_pool_holdout")
    plt.close(fig)


def fig5_receptor(D: dict) -> None:
    """Main Fig 5: receptor-realization sensitivity, including the mTOR 4JSX swap on PIK3CA/mTOR."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.55), gridspec_kw={"width_ratios": [1.15, 1.0]})
    colors = {"4L23": C["main"], "4JPS": C["holdout"], "5DXT": C["a_only"], "4JSX": C["gnina"]}

    ax = axes[0]
    panel_label(ax, "A", x=-0.16, y=1.05)
    pm_points = [
        ("4L23", fnum(D["theta6"]["PIK3CA/mTOR"]["pocket_matched_summary_min"]),
         fnum(D["theta6"]["PIK3CA/mTOR"]["ci_lo"]), fnum(D["theta6"]["PIK3CA/mTOR"]["ci_hi"])),
        ("4JPS", fnum(D["jps"]["summary_min"]), fnum(D["jps"]["summary_min_ci_lo"]), fnum(D["jps"]["summary_min_ci_hi"])),
        ("5DXT", fnum(D["dxt"]["summary_min"]), fnum(D["dxt"]["summary_min_ci_lo"]), fnum(D["dxt"]["summary_min_ci_hi"])),
        ("4JSX", fnum(D["jsx"]["summary_min"]), fnum(D["jsx"]["summary_min_ci_lo"]), fnum(D["jsx"]["summary_min_ci_hi"])),
    ]
    fig5 = []
    for i, (crystal, y, lo, hi) in enumerate(pm_points):
        ax.errorbar(
            i, y, yerr=[[y - lo], [hi - y]],
            fmt="o", color=colors[crystal], ecolor=colors[crystal],
            elinewidth=1.5, capsize=2.4, markersize=6.5, zorder=4,
        )
        fig5.append({"pair": "PIK3CA/mTOR", "crystal": crystal, "y": y, "lo": lo, "hi": hi})
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["4L23\n(PIK3CA)", "4JPS\n(PIK3CA)", "5DXT\n(PIK3CA)", "4JSX\n(mTOR)"], fontsize=6.5)
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0.12, 1.02)
    ax.set_title("PIK3CA/mTOR", fontsize=FS_AXIS, pad=4)
    ax.set_xlim(-0.45, 3.45)

    ax = axes[1]
    panel_label(ax, "B", x=-0.16, y=1.05)
    pab_points = [
        ("4L23", fnum(D["theta6"]["PIK3CA/PIK3CB"]["pocket_matched_summary_min"]),
         fnum(D["theta6"]["PIK3CA/PIK3CB"]["ci_lo"]), fnum(D["theta6"]["PIK3CA/PIK3CB"]["ci_hi"])),
        ("4JPS", fnum(D["pab_jps"]["summary_min"]), fnum(D["pab_jps"]["summary_min_ci_lo"]), fnum(D["pab_jps"]["summary_min_ci_hi"])),
        ("5DXT", fnum(D["pab_dxt"]["summary_min"]), fnum(D["pab_dxt"]["summary_min_ci_lo"]), fnum(D["pab_dxt"]["summary_min_ci_hi"])),
    ]
    for i, (crystal, y, lo, hi) in enumerate(pab_points):
        ax.errorbar(
            i, y, yerr=[[y - lo], [hi - y]],
            fmt="o", color=colors[crystal], ecolor=colors[crystal],
            elinewidth=1.5, capsize=2.4, markersize=6.5, zorder=4,
        )
        fig5.append({"pair": "PIK3CA/PIK3CB", "crystal": crystal, "y": y, "lo": lo, "hi": hi})
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["4L23\n(PIK3CA)", "4JPS\n(PIK3CA)", "5DXT\n(PIK3CA)"], fontsize=6.5)
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0.12, 1.02)
    ax.set_title("PIK3CA/PIK3CB (B = 2WXF frozen)", fontsize=FS_AXIS, pad=4)
    ax.set_xlim(-0.45, 2.45)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=colors["4L23"], ls="none", ms=6, label="4L23 original"),
            Line2D([0], [0], marker="o", color=colors["4JPS"], ls="none", ms=6, label="4JPS"),
            Line2D([0], [0], marker="o", color=colors["5DXT"], ls="none", ms=6, label="5DXT"),
            Line2D([0], [0], marker="o", color=colors["4JSX"], ls="none", ms=6, label="4JSX (mTOR only)"),
        ],
        loc="upper left",
        fontsize=6.0,
        frameon=False,
    )
    PROVENANCE["plotted"]["fig5"] = fig5
    fig.subplots_adjust(wspace=0.34, left=0.09, right=0.98, top=0.86, bottom=0.18)
    save_all(fig, "Fig5_receptor_realization")
    plt.close(fig)


def toc_graphic() -> None:
    """ACS TOC: 3.25 x 1.75 in, no numerical results, not a crop of Fig 1."""
    fig = plt.figure(figsize=(3.25, 1.75), dpi=300)
    ax = fig.add_axes([0.02, 0.06, 0.96, 0.88])
    ax.set_xlim(0, 32.0)
    ax.set_ylim(0, 17.5)
    ax.axis("off")

    # LEFT: two pockets + four classes
    ax.add_patch(FancyBboxPatch((0.35, 7.3), 6.0, 9.2, boxstyle="round,pad=0.08,rounding_size=0.4",
                                facecolor="#F4F7FA", edgecolor=C["a_only"], lw=1.1))
    ax.add_patch(FancyBboxPatch((6.85, 7.3), 6.0, 9.2, boxstyle="round,pad=0.08,rounding_size=0.4",
                                facecolor="#F4F7FA", edgecolor=C["b_only"], lw=1.1))
    ax.text(3.35, 15.6, "Pocket A", ha="center", fontsize=8, fontweight="bold")
    ax.text(9.85, 15.6, "Pocket B", ha="center", fontsize=8, fontweight="bold")
    ax.text(6.6, 6.0, "dual vs hard-negatives", ha="center", fontsize=6.5, color="#444444")
    xs = [2.1, 5.1, 8.2, 11.1]
    names = ["dual", "A_only", "B_only", "neither"]
    cols = [(C["dual"], C["dual"]), (C["a_only"], "#DDDDDD"), ("#DDDDDD", C["b_only"]), ("#DDDDDD", "#DDDDDD")]
    for x, name, (c1, c2) in zip(xs, names, cols):
        ax.add_patch(Circle((x - 0.42, 3.35), 0.52, facecolor=c1, edgecolor=C["ink"], lw=0.5))
        ax.add_patch(Circle((x + 0.42, 3.35), 0.52, facecolor=c2, edgecolor=C["ink"], lw=0.5))
        ax.text(x, 1.55, name, ha="center", fontsize=6.5)

    # MIDDLE: opposite-pocket arrow
    ax.annotate("", xy=(21.2, 11.2), xytext=(14.5, 11.2),
                arrowprops=dict(arrowstyle="-|>", color=C["vina"], lw=1.5))
    ax.text(17.85, 12.55, "score opposite pocket", ha="center", fontsize=7, color=C["vina"])
    ax.text(17.85, 8.75, r"summary = min(arms)", ha="center", fontsize=7)

    # RIGHT: qualitative four dots, NO numerical AUROCs
    ax.plot([23.6, 30.4], [8.85, 8.85], ls="--", color=C["chance"], lw=0.9)
    ax.text(27.0, 9.55, "chance", ha="center", va="bottom", fontsize=7, color=C["chance"])
    ys = [7.45, 8.85, 9.55, 12.25]  # below / on / slightly above / modestly above — qualitative only
    cols_d = [C["egfr"], C["metal"], C["other"], C["thick"]]
    for y, col in zip(ys, cols_d):
        ax.scatter([25.6], [y], s=28, color=col, zorder=3, edgecolors=C["ink"], linewidths=0.3)
    ax.text(26.6, 15.45, "pair-dependent", ha="center", fontsize=8, fontweight="bold")
    ax.plot(26.85, 12.25, marker="*", color=C["desc"], markersize=6.5, linestyle="none")
    ax.text(26.6, 13.55, "size / crystal", ha="center", fontsize=6.5, color="#666666")

    save_all(fig, "TOC_graphic", toc=True)
    plt.close(fig)


def verify(D: dict) -> None:
    """Fail the run if any plotted value disagrees with its CSV source."""
    errors = []

    def eq(a, b, tol=1e-9, msg=""):
        if abs(float(a) - float(b)) > tol:
            errors.append(f"{msg}: plotted {a} != source {b}")

    t6 = D["theta6"]
    p3 = PROVENANCE["plotted"]["fig3"]
    expected_dir = {
        "EGFR/HER2": (0.4297, 0.284, 0.5759),
        "AChE/BChE": (0.6058, 0.4396, 0.74),
        "PIK3CA/PIK3CB": (0.5, 0.3468, 0.648),
        "PIK3CA/mTOR": (0.6921, 0.4638, 0.8015),
    }
    expected_neither = {
        "EGFR/HER2": (0.756, 0.5625, 0.9197, 12),
        "AChE/BChE": (0.6494, 0.484, 0.8123, 15),
        "PIK3CA/PIK3CB": (0.5592, 0.3728, 0.7456, 16),
        "PIK3CA/mTOR": (0.5139, 0.2222, 0.8056, 4),
    }
    for pair in PAIR_ORDER:
        eq(p3[pair]["directional"]["y"], t6[pair]["pocket_matched_summary_min"], msg=f"fig3 {pair} dir")
        eq(p3[pair]["directional"]["lo"], t6[pair]["ci_lo"], msg=f"fig3 {pair} dir lo")
        eq(p3[pair]["directional"]["hi"], t6[pair]["ci_hi"], msg=f"fig3 {pair} dir hi")
        ey, elo, ehi = expected_dir[pair]
        eq(p3[pair]["directional"]["y"], ey, msg=f"fig3 {pair} dir checksum")
        eq(p3[pair]["directional"]["lo"], elo, msg=f"fig3 {pair} dir lo checksum")
        eq(p3[pair]["directional"]["hi"], ehi, msg=f"fig3 {pair} dir hi checksum")
        nr = D["form_by"][pair]["D_vs_neither_mean"]
        eq(p3[pair]["neither"]["y"], nr["auroc"], msg=f"fig3 {pair} neither")
        eq(p3[pair]["neither"]["lo"], nr["ci_lo"], msg=f"fig3 {pair} neither lo")
        eq(p3[pair]["neither"]["hi"], nr["ci_hi"], msg=f"fig3 {pair} neither hi")
        ny, nlo, nhi, nn = expected_neither[pair]
        eq(p3[pair]["neither"]["y"], ny, msg=f"fig3 {pair} neither checksum")
        eq(p3[pair]["neither"]["lo"], nlo, msg=f"fig3 {pair} neither lo checksum")
        eq(p3[pair]["neither"]["hi"], nhi, msg=f"fig3 {pair} neither hi checksum")
        if p3[pair]["neither"]["n_neg"] != nn:
            errors.append(f"fig3 {pair} n_neg {p3[pair]['neither']['n_neg']} != {nn}")

    pS = PROVENANCE["plotted"]["figS_forest"]
    expected_vina = {
        "EGFR/HER2": (0.4297, 0.284, 0.5759, 0.6664, 0.4297),
        "AChE/BChE": (0.6058, 0.4396, 0.74, 0.6504, 0.6058),
        "PIK3CA/PIK3CB": (0.5, 0.3468, 0.648, 0.6905, 0.5),
        "PIK3CA/mTOR": (0.6921, 0.4638, 0.8015, 0.7143, 0.6921),
    }
    expected_desc_arm = {"EGFR/HER2": "clogp", "AChE/BChE": "tpsa", "PIK3CA/PIK3CB": "heavy", "PIK3CA/mTOR": "heavy"}
    expected_desc_y = {"EGFR/HER2": 0.4821, "AChE/BChE": 0.7333, "PIK3CA/PIK3CB": 0.6217, "PIK3CA/mTOR": 0.463}
    expected_gnina = {"EGFR/HER2": 0.2902, "AChE/BChE": 0.4127, "PIK3CA/PIK3CB": 0.5332, "PIK3CA/mTOR": 0.6548}
    expected_rtm = {"EGFR/HER2": 0.3527, "AChE/BChE": 0.5185, "PIK3CA/PIK3CB": 0.5421, "PIK3CA/mTOR": 0.6151}

    for pair in PAIR_ORDER:
        eq(pS[pair]["vina"]["y"], t6[pair]["pocket_matched_summary_min"], msg=f"figS_forest {pair} vina")
        eq(pS[pair]["vina"]["lo"], t6[pair]["ci_lo"], msg=f"figS_forest {pair} vina lo")
        eq(pS[pair]["vina"]["hi"], t6[pair]["ci_hi"], msg=f"figS_forest {pair} vina hi")
        ey, elo, ehi, eda, edb = expected_vina[pair]
        eq(pS[pair]["vina"]["y"], ey, msg=f"figS_forest {pair} vina checksum")
        eq(pS[pair]["vina"]["lo"], elo, msg=f"figS_forest {pair} vina lo checksum")
        eq(pS[pair]["vina"]["hi"], ehi, msg=f"figS_forest {pair} vina hi checksum")
        eq(pS[pair]["vina"]["da"], eda, msg=f"figS_forest {pair} DA checksum")
        eq(pS[pair]["vina"]["db"], edb, msg=f"figS_forest {pair} DB checksum")
        r = D["pm_by"][(pair, "pocket_matched_rtm")]
        eq(pS[pair]["rtm"]["y"], r["summary_min"], msg=f"figS_forest {pair} rtm")
        eq(pS[pair]["rtm"]["lo"], r["summary_min_ci_lo"], msg=f"figS_forest {pair} rtm lo")
        eq(pS[pair]["rtm"]["hi"], r["summary_min_ci_hi"], msg=f"figS_forest {pair} rtm hi")
        eq(pS[pair]["rtm"]["y"], expected_rtm[pair], msg=f"figS_forest {pair} rtm checksum")
        g = D["gnina9"][pair]
        eq(pS[pair]["gnina_best9"]["y"], g["summary_min"], msg=f"figS_forest {pair} gnina")
        eq(pS[pair]["gnina_best9"]["lo"], g["ci_lo"], msg=f"figS_forest {pair} gnina lo")
        eq(pS[pair]["gnina_best9"]["hi"], g["ci_hi"], msg=f"figS_forest {pair} gnina hi")
        eq(pS[pair]["gnina_best9"]["y"], expected_gnina[pair], msg=f"figS_forest {pair} gnina checksum")
        _, arm, dr = D["best_desc"][pair]
        eq(pS[pair]["best_desc"]["y"], dr["summary_min"], msg=f"figS_forest {pair} desc {arm}")
        if pS[pair]["best_desc"]["arm"] != expected_desc_arm[pair]:
            errors.append(f"figS_forest {pair} desc arm {pS[pair]['best_desc']['arm']} != {expected_desc_arm[pair]}")
        eq(pS[pair]["best_desc"]["y"], expected_desc_y[pair], msg=f"figS_forest {pair} desc checksum")

    for i, pair in enumerate(PAIR_ORDER):
        eq(PROVENANCE["plotted"]["fig4A_DA"][i], t6[pair]["auroc_D_vs_A"], msg=f"fig4A DA {pair}")
        eq(PROVENANCE["plotted"]["fig4A_DB"][i], t6[pair]["auroc_D_vs_B"], msg=f"fig4A DB {pair}")

    # Fig 2 highlights vs J0
    j0_map = {r["pair"]: r for r in D["j0"] if (r.get("min_strict_hardneg") or "").strip()}
    for pair, val in PROVENANCE["plotted"]["fig2A_highlighted"].items():
        eq(val, j0_map[pair]["min_strict_hardneg"], msg=f"fig2 {pair}")

    if PROVENANCE["plotted"]["fig2A_n_pairs"] != 49:
        errors.append(f"fig2A_n_pairs {PROVENANCE['plotted']['fig2A_n_pairs']} != 49")
    expected_hl = {"PIK3CA/MTOR": 80, "ACHE/BCHE": 78, "PIK3CA/PIK3CB": 56, "EGFR/HER2": 7, "HDAC1/HDAC6": 93}
    for pair, val in expected_hl.items():
        eq(PROVENANCE["plotted"]["fig2A_highlighted"][pair], val, msg=f"fig2 highlight checksum {pair}")

    # S12
    chem = PROVENANCE["plotted"]["fig2B_chembl"]
    bdb = PROVENANCE["plotted"]["fig2B_bindingdb_equal_only"]
    names = ["PIK3CA/MTOR", "ACHE/BCHE", "PIK3CA/PIK3CB", "EGFR/HER2"]
    for i, name in enumerate(names):
        c = next(r for r in D["s12"] if r["pair"] == name and r["source"] == "ChEMBL_cache" and r["rule"] == "pChEMBL")
        b = next(r for r in D["s12"] if r["pair"] == name and r["source"] == "BindingDB" and r["rule"] == "equal_only")
        eq(chem[i], c["min_strict_hardneg"], msg=f"fig2B chem {name}")
        eq(bdb[i], b["min_strict_hardneg"], msg=f"fig2B bdb {name}")

    eq(chem[0], 80, msg="fig2B chem PM checksum")
    eq(bdb[0], 76, msg="fig2B bdb PM checksum")
    eq(chem[1], 78, msg="fig2B chem AChE checksum")
    eq(bdb[1], 92, msg="fig2B bdb AChE checksum")
    eq(chem[2], 56, msg="fig2B chem PIK3CB checksum")
    eq(bdb[2], 58, msg="fig2B bdb PIK3CB checksum")
    eq(chem[3], 7, msg="fig2B chem EGFR checksum")
    eq(bdb[3], 31, msg="fig2B bdb EGFR checksum")

    # Fig S holdout (former Fig 5A)
    expected_hold = {
        "AChE/BChE": (0.6175, 0.4216, 0.7593),
        "PIK3CA/PIK3CB": (0.425, 0.2406, 0.6184),
        "PIK3CA/mTOR": (0.765, 0.6025, 0.8911),
    }
    for row in PROVENANCE["plotted"]["figS_holdout"]:
        pair = row["pair"]
        eq(row["main"], t6[pair]["pocket_matched_summary_min"], msg=f"figS_holdout {pair} main")
        eq(row["hold"], D["hold_pm"][pair]["summary_min"], msg=f"figS_holdout {pair} hold")
        eq(row["hold_lo"], D["hold_pm"][pair]["summary_min_ci_lo"], msg=f"figS_holdout {pair} hold lo")
        eq(row["hold_hi"], D["hold_pm"][pair]["summary_min_ci_hi"], msg=f"figS_holdout {pair} hold hi")
        ey, elo, ehi = expected_hold[pair]
        eq(row["hold"], ey, msg=f"figS_holdout {pair} hold checksum")
        eq(row["hold_lo"], elo, msg=f"figS_holdout {pair} hold lo checksum")
        eq(row["hold_hi"], ehi, msg=f"figS_holdout {pair} hold hi checksum")
    if "EGFR/HER2" in {row["pair"] for row in PROVENANCE["plotted"]["figS_holdout"]}:
        errors.append("figS_holdout must not include EGFR/HER2 (no holdout)")

    # Fig 5: receptor realization; 4JSX is the mTOR swap on PIK3CA/mTOR only
    b = {(row["pair"], row["crystal"]): row for row in PROVENANCE["plotted"]["fig5"]}
    eq(b[("PIK3CA/mTOR", "4L23")]["y"], t6["PIK3CA/mTOR"]["pocket_matched_summary_min"], msg="fig5 PM 4L23")
    eq(b[("PIK3CA/mTOR", "4JPS")]["y"], D["jps"]["summary_min"], msg="fig5 PM 4JPS")
    eq(b[("PIK3CA/mTOR", "5DXT")]["y"], D["dxt"]["summary_min"], msg="fig5 PM 5DXT")
    eq(b[("PIK3CA/mTOR", "4JSX")]["y"], D["jsx"]["summary_min"], msg="fig5 PM 4JSX")
    eq(b[("PIK3CA/mTOR", "4JSX")]["lo"], D["jsx"]["summary_min_ci_lo"], msg="fig5 PM 4JSX lo")
    eq(b[("PIK3CA/mTOR", "4JSX")]["hi"], D["jsx"]["summary_min_ci_hi"], msg="fig5 PM 4JSX hi")
    eq(b[("PIK3CA/PIK3CB", "4L23")]["y"], t6["PIK3CA/PIK3CB"]["pocket_matched_summary_min"], msg="fig5 PAB 4L23")
    eq(b[("PIK3CA/PIK3CB", "4JPS")]["y"], D["pab_jps"]["summary_min"], msg="fig5 PAB 4JPS")
    eq(b[("PIK3CA/PIK3CB", "5DXT")]["y"], D["pab_dxt"]["summary_min"], msg="fig5 PAB 5DXT")
    eq(b[("PIK3CA/mTOR", "4JPS")]["y"], 0.4861, msg="fig5 PM 4JPS checksum")
    eq(b[("PIK3CA/mTOR", "5DXT")]["y"], 0.5046, msg="fig5 PM 5DXT checksum")
    eq(b[("PIK3CA/mTOR", "4JSX")]["y"], 0.6389, msg="fig5 PM 4JSX checksum")
    eq(b[("PIK3CA/mTOR", "4JSX")]["lo"], 0.4178, msg="fig5 PM 4JSX lo checksum")
    eq(b[("PIK3CA/mTOR", "4JSX")]["hi"], 0.7758, msg="fig5 PM 4JSX hi checksum")
    eq(b[("PIK3CA/PIK3CB", "4JPS")]["y"], 0.6905, msg="fig5 PAB 4JPS checksum")
    eq(b[("PIK3CA/PIK3CB", "5DXT")]["y"], 0.6849, msg="fig5 PAB 5DXT checksum")
    eq(b[("PIK3CA/PIK3CB", "4JPS")]["lo"], 0.5156, msg="fig5 PAB 4JPS lo checksum")
    eq(b[("PIK3CA/PIK3CB", "5DXT")]["lo"], 0.5057, msg="fig5 PAB 5DXT lo checksum")
    if ("PIK3CA/PIK3CB", "4JSX") in b:
        errors.append("fig5 must not place the mTOR 4JSX swap on PIK3CA/PIK3CB")

    expected_wrong_main = {"EGFR/HER2": 0.26, "AChE/BChE": 0.4444, "PIK3CA/PIK3CB": 0.3489, "PIK3CA/mTOR": 0.6019}
    for i, pair in enumerate(PAIR_ORDER):
        eq(PROVENANCE["plotted"]["siA_matched"][i], t6[pair]["pocket_matched_summary_min"], msg=f"siA m {pair}")
        eq(PROVENANCE["plotted"]["siA_wrong"][i], D["pm_by"][(pair, "wrong_pocket_control_vina")]["summary_min"], msg=f"siA w {pair}")
        eq(PROVENANCE["plotted"]["siA_wrong"][i], expected_wrong_main[pair], msg=f"siA w checksum {pair}")
        if PROVENANCE["plotted"]["siA_matched"][i] <= PROVENANCE["plotted"]["siA_wrong"][i]:
            errors.append(f"siA {pair}: matched should exceed wrong-pocket")

    expected_wrong_hold = {"AChE/BChE": (0.6175, 0.6425), "PIK3CA/PIK3CB": (0.425, 0.52), "PIK3CA/mTOR": (0.765, 0.7875)}
    for i, pair in enumerate(["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]):
        eq(PROVENANCE["plotted"]["siB_matched"][i], expected_wrong_hold[pair][0], msg=f"siB m checksum {pair}")
        eq(PROVENANCE["plotted"]["siB_wrong"][i], expected_wrong_hold[pair][1], msg=f"siB w checksum {pair}")
        if PROVENANCE["plotted"]["siB_wrong"][i] < PROVENANCE["plotted"]["siB_matched"][i]:
            errors.append(f"siB {pair}: holdout wrong-pocket should be >= matched")

    from plot_jcim_si_composites_v1 import verify_si
    verify_si(D, PROVENANCE, errors)

    # TPSA means vs descriptor_by_class
    desc = list(csv.DictReader((DATA / "jcim_bench_v0/tables/ache_descriptor_by_class_v1.csv").open()))
    mean_map = {(r["class"], r["feature"]): fnum(r["mean"]) for r in desc}
    for cls, m in zip(["dual", "A_only", "B_only"], PROVENANCE["plotted"]["fig4C_tpsa_mean"]):
        eq(m, mean_map[(cls, "tpsa")], tol=1e-3, msg=f"fig4C mean {cls}")
    if PROVENANCE["plotted"]["fig4C_tpsa_n"] != [27, 25, 28]:
        errors.append(f"fig4C n {PROVENANCE['plotted']['fig4C_tpsa_n']} != [27, 25, 28]")
    eq(PROVENANCE["plotted"]["fig4C_tpsa_mean"][0], 75.2559, tol=1e-3, msg="fig4C dual mean checksum")
    eq(PROVENANCE["plotted"]["fig4C_tpsa_mean"][1], 53.7792, tol=1e-3, msg="fig4C A_only mean checksum")
    eq(PROVENANCE["plotted"]["fig4C_tpsa_mean"][2], 47.8796, tol=1e-3, msg="fig4C B_only mean checksum")

    # dpi / size / RGB checks
    from PIL import Image
    checks = {
        "Fig3_formulation_comparison.png": (7.0, None),
        "Fig6_wrong_pocket_paradox.png": (7.0, None),
        "Fig7_confound_anatomy.png": (7.0, None),
        "TOC_graphic.tif": (3.25, 1.75),
    }
    for name, (w_in, h_in) in checks.items():
        im = Image.open(OUT / name)
        if im.mode != "RGB":
            errors.append(f"{name} mode {im.mode} != RGB")
        dpi = im.info.get("dpi", (None, None))[0]
        if dpi is None:
            # matplotlib tiff may omit dpi; compute from pixels / figsize
            if name.startswith("TOC"):
                dpi = im.size[0] / 3.25
        if dpi is not None and dpi < 295:
            errors.append(f"{name} dpi {dpi} < 300")
        if h_in is not None:
            # TOC must be exact
            if abs(im.size[0] / 300 - w_in) > 0.05 or abs(im.size[1] / 300 - h_in) > 0.05:
                # allow using actual dpi
                d = dpi or 300
                if abs(im.size[0] / d - w_in) > 0.08 or abs(im.size[1] / d - h_in) > 0.08:
                    errors.append(f"{name} size {im.size} not {w_in}x{h_in} in at dpi={d}")

    if errors:
        raise SystemExit("VERIFICATION FAILED:\n" + "\n".join(errors))
    print("verification OK:", json.dumps({k: PROVENANCE["plotted"].get(k) for k in ("fig2A_n_pairs", "fig2A_highlighted")}, indent=2))


def write_captions() -> None:
    text = """# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match the files in `figures/jcim_article/`.
All numerical values are those plotted from the frozen CSVs (unrounded). Table 2 in the text may round to three decimals.

Regenerate (from `Dual_Target_Docking/`):
`python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py`
The script re-reads the CSVs, writes PDF/PNG/TIF, and fails if any plotted value disagrees with its source (`plotted_values.json`).

## Figure 1. Dual-target docking as dual-versus-selective discrimination.

(A) A strict dual-target benchmark distinguishes four experimentally labeled ligand states: dual-active (D), A-selective (A_only), B-selective (B_only), and neither. A_only and B_only are selectivity hard negatives: they are active on one target and can produce plausible docking scores, yet they lack activity on the other. Neither is curated as part of the four-state panel but is not used in the primary AUROCs. (B) The prespecified primary readout is two directional pairwise discriminations, not a four-class classifier and not a pooled score. Dual versus A_only is scored in pocket B; dual versus B_only is scored in pocket A. The pair-level summary is the weaker arm (summary_min), so a favorable score on one target cannot hide directional failure on the other. Pooled docking scores are retained only as a control.

## Figure 2. Public-data supply of strict hard negatives.

(A) Minimum of the two strict hard-negative counts (A_only, B_only) for every target pair in the J0 ChEMBL audit (`j0_strict_label_supply.csv`). Dashed line, thick-panel gate (≥50); dotted line, thin-panel gate (≥20). Highlighted: the three thick pairs used as K=4 main panels, EGFR/HER2 (7 B_only; supply-limited case), and HDAC1/HDAC6 (metal enzyme; excluded). (B) Count-level comparison of the same four pairs in ChEMBL pChEMBL versus BindingDB equal-relation measurements (Table S12). No docking.

## Figure 3. Benchmark formulation changes the apparent evidence for dual-target recognition.

Same frozen AutoDock Vina scores under two task formulations (unified θ = 6.0). Dark bars: directional pocket-matched `summary_min` with 95% ligand-bootstrap CIs from `unified_threshold_sensitivity_v2.csv` (Table 2). Orange bars: Dual-versus-neither comparator using pooled `vina_mean`, with CIs from `formulation_conventional_vs_directional_v1.csv` (Table 3). Dual versus neither is a nonselectivity-controlled comparator, not “the conventional dual-target benchmark.” EGFR/HER2 is the proof-of-principle gap (0.756 versus 0.430). AChE/BChE and PIK3CA/PIK3CB increments are small and overlapping. PIK3CA/mTOR Dual versus neither is hatched as underpowered (neither n = 4) and is not a reverse-overestimation result. Dual versus all non-duals is reported in Table 3, not here. Vertical dashed line, chance (0.5). The former pocket-matched forest (Vina / RTM / GNINA / descriptor) is Figure S4.

## Figure 4. Weak-arm asymmetry and physicochemical confounding.

(A) Directional Vina AUROCs at θ = 6.0: dual versus A_only (pocket B) and dual versus B_only (pocket A). (B) Vina pocket-matched summary_min versus the best single-descriptor reference, with 95% CIs. (C) TPSA on the AChE/BChE panel by class (individual ligands from `assembled_AChE_BChE.csv`; horizontal line, median). Dual ligands are more polar than either hard-negative class, matching the TPSA reference that exceeds Vina on this pair.

## Figure 5. Receptor realization can raise or lower apparent dual-target discrimination.

(A) PIK3CA/mTOR (PM48): pocket-matched summary_min after replacing one receptor at a time. 4L23 is the original PIK3CA structure (B = 4JT6 frozen); 4JPS and 5DXT replace PIK3CA; 4JSX replaces mTOR (A = 4L23 frozen). Point estimates: 0.692 (4L23), 0.486 (4JPS), 0.505 (5DXT), 0.639 (4JSX). (B) PIK3CA/PIK3CB: the same PIK3CA crystals with 2WXF held frozen. Point estimates: 0.500 (4L23), 0.691 (4JPS), 0.685 (5DXT). Error bars are 95% ligand-bootstrap CIs from the deposited swap tables and Table 2. 4JSX is an mTOR swap and is not applied to PIK3CA/PIK3CB. Receptor replacement is a realization effect, not a unidirectional collapse and not a robustness certificate. PIK3CA/PIK3CB uses the same 99-ligand set as Table 2 (PAB_034 timeout on original 4L23 and on both 4JPS and 5DXT). Unused-pool holdout is Figure S5.

## Figure 6. Wrong-pocket controls reveal an unresolved out-of-panel failure mode.

(A) Main K=4 panels: pocket-matched Vina summary_min versus the wrong-pocket control (`pocket_matched_directional_v1.csv`). Matched exceeds wrong-pocket on all four pairs. (B) Unused-pool holdout: the inequality reverses (wrong-pocket ≥ matched) on all three pairs with holdout supply (`holdout_pocket_matched_v1.csv`). EGFR/HER2 has no holdout. (C) The reversal remains after potency matching (|Δp| ≤ 0.5) or size matching (|Δheavy| ≤ 2) (`holdout_matched_wrong_pocket_summary_v1.csv`). Wrong-pocket remains ≥ matched on all nine cells; matching does not restore the main-panel inequality. (D) Scoring-free contact-count AUROC on pocket A (D vs A_only) and pocket B (D vs B_only) versus Vina wrong-pocket summary_min (`wrong_pocket_contact_v1_output.txt`; not a PLIF). B-arm contact is above chance; the magnitude does not reproduce Vina on PIK3CA/mTOR.

## Figure 7. Ligand-structure association and matched-subset tests.

(A) ECFP4 logistic regression under scaffold GroupKFold versus pocket-matched Vina on both directional arms (`ligand_ml_baseline_scaffold_cv_v1.csv`). Fingerprint AUROCs are chemotype–label association, not evidence of pocket physics. (B) Pocket-matched Vina versus all four prespecified physicochemical descriptors (heavy-atom count, MW, cLogP, TPSA) with 95% CIs. Descriptor CIs are from `forest_summary_min_ci_v1.csv`; Vina from θ = 6.0. Figure 4 reports only the best single-descriptor reference per pair; this panel shows all four. (C) Weak-arm (D vs B_only) logistic AUROC of Vina alone versus Vina plus heavy-atom count and TPSA, with the Vina odds ratio (`covariate_adjusted_v1.csv`). EGFR/HER2 score-only in that table is 0.5703 (the table’s logistic AUROC of feature `vina_A`), which is not the rank AUROC 0.4297 in Table 2. (D) D vs B_only after potency matching (|Δp| ≤ 0.5) or size matching (|Δheavy| ≤ 2) versus the unmatched full-panel contrast (`matched_subset_directional_v1.csv`). Error bars are the table’s single-contrast 95% CIs.

## Figure S1. Protocol knobs that do not change the ranking.

(A) Pocket-matched summary_min across the unified label-threshold grid (`unified_threshold_sensitivity_v2.csv`). Open markers are underpowered cells (EGFR/HER2 strict, n_B_only=7; PIK3CA/mTOR at θ=5.5, n_B_only=5, and strict, n_B_only=4). At the primary θ=6.0, PIK3CA/mTOR is the highest point estimate; AChE/BChE is flat at 0.6058 across the grid. The underpowered θ=5.5 PIK3CA/mTOR cell (0.5017) is not a ranking contradiction. (B) GNINA CNN mode01 versus best-of-9 versus the same-panel Vina reference (`gnina_pocket_matched_mode01_vs_best9_k4_v1.csv`). Best-of-9 versus mode01 moves summary_min by −0.04 to +0.08. EGFR/HER2 and AChE/BChE remain below chance on both GNINA channels. PIK3CA/PIK3CB GNINA best-of-9 is 0.533 versus Vina 0.500 (near chance). (C) PIK3CA/mTOR PM48 versus the PM110 expansion for Vina, RTMScore, and GNINA best-of-9 (`pm110_vs_pm48_pocket_matched_v1.csv`). (D) PM48 Vina at exhaustiveness 16 versus 8, computed from `scores_vina_E8_best.csv` (empty affinities skipped; ligands labeled neither were excluded) with the same pocket-matched definition, beside single-target enrichment AUROC and EF1% on 4L23 and 4JT6 (`single_target_enrichment_v1.csv`).

## Figure S2. Equal-relation supply and holdout sampling shift.

(A) Minimum strict hard-negative counts for the K=4 pairs in ChEMBL pChEMBL, BindingDB/PubChem `equal_only`, and BindingDB/PubChem `as_is` (`crossdb_strict_supply_v1.csv`). Count-level only; no docking. `as_is` lets EGFR/HER2 pass ≥50 because censored `>` values are treated as point estimates; `equal_only` does not. (B) Holdout minus main-panel mean pChEMBL for dual pA, A_only pA, and B_only pB (`holdout_vs_main_potency_size_v1.csv`). Sampling shift is real, especially on PIK3CA/mTOR, but does not reverse Figure 6C.

## Figure S3. Paired bootstrap differences that Figure 6 does not show.

All values are from `wrong_pocket_paired_delta_bootstrap_v1.csv` and `pocket_matched_vs_best_descriptor_delta_v1.csv` (B = 2000 ligand resamples, seed 20260729). Point Δ equals the rounded Table 2 / Figure 6 AUROCs subtracted at four decimals, not a separately rounded difference. Blue, 95% CI excludes 0; gray, CI includes 0. (A) Main K=4 panels: Δ = pocket-matched − wrong-pocket summary_min. Point Δ is positive on all four pairs (EGFR/HER2 0.1697, AChE/BChE 0.1614, PIK3CA/PIK3CB 0.1511, PIK3CA/mTOR 0.0902). Only EGFR/HER2 and AChE/BChE have CIs that exclude 0; PIK3CA/PIK3CB and PIK3CA/mTOR CIs include 0. (B) Unused-pool holdout: point Δ is negative on all three eligible pairs (wrong-pocket ≥ matched), and every CI includes 0. EGFR/HER2 has no holdout. This panel is the interval on the Figure 6B reversal, not a new docking experiment. (C) Pocket-matched Vina minus the best single-descriptor reference (EGFR/HER2 cLogP 0.4821; AChE/BChE TPSA 0.7333; PIK3CA/PIK3CB and PIK3CA/mTOR heavy-atom count). All four CIs include 0, including PIK3CA/mTOR +0.2291 [−0.0105, 0.4352]. This is not the pooled `vina_mean` gate (EGFR/HER2 0.2824). (D) ECFP4 logistic AUROC under scaffold GroupKFold versus random StratifiedKFold (`ligand_ml_scaffold_vs_random_v1.csv`). Mean (random − scaffold) across eight directional contrasts is 0.0112. Scaffold split remains the primary ML readout; this is a leakage check, not a search for a leakier split.

## Figure S4. Pocket-matched summary_min on the frozen K=4 set (former main Figure 3).

Vina (primary), RTMScore, GNINA CNN best-of-9, and the best single-descriptor reference (heavy-atom count, MW, cLogP, or TPSA) with 95% ligand-bootstrap CIs. Vina CIs are the θ = 6.0 values in `unified_threshold_sensitivity_v2.csv` (Table 2). Best descriptor (right column, from `forest_summary_min_ci_v1.csv`): EGFR/HER2 cLogP; AChE/BChE TPSA; PIK3CA/PIK3CB and PIK3CA/mTOR heavy-atom count. Vertical dashed line, chance (0.5). GNINA is a single CNN channel, not a three-engine competition. The main-text formulation comparison is Figure 3.

## Figure S5. Unused-pool holdout versus the main panel.

Pocket-matched summary_min on the main panel versus the unused-pool holdout (20/20/20; seed 20260731) for the three pairs with unused-pool supply. EGFR/HER2 has no holdout. PM110 is a same-family stability check shown in Figure S1C, not a third independent validation trajectory. This panel was formerly Figure 5A; the main-text Figure 5 is receptor realization only.

## TOC graphic (For Table of Contents Only).

DualFourClass-Bench asks whether docking can distinguish experimentally labeled dual-active ligands from single-target selective hard negatives in both pockets, rather than whether both docking scores are merely favorable. The graphic does not report numerical AUROCs and is not a reuse of Figure 1.
"""
    (OUT / "CAPTIONS.md").write_text(text)


def main() -> None:
    apply_style()
    D = load()
    fig1_task()
    fig2_supply(D)
    fig3_forest(D)
    fig3_formulation(D)
    fig4_confounds(D)
    figS_holdout(D)
    fig5_receptor(D)
    from plot_jcim_si_composites_v1 import draw_all, extend_load
    extend_load(D, _read, PROVENANCE)
    draw_all(D, PROVENANCE)
    toc_graphic()
    write_captions()
    (OUT / "plotted_values.json").write_text(json.dumps(PROVENANCE, indent=2, default=str))
    verify(D)
    obsolete = [
        "FigS1_wrong_pocket_main_vs_holdout",
        "FigS2_protocol_sensitivity",
        "FigS3_confound_anatomy",
        "FigS4_holdout_mechanism_and_supply",
        "Fig3_pocket_matched_forest",
        "Fig5_holdout_and_crystal_swap",
    ]
    for stem in obsolete:
        for ext in ("pdf", "png", "tif"):
            p = OUT / f"{stem}.{ext}"
            if p.exists():
                p.unlink()
    print("wrote", OUT)
    for p in sorted(OUT.glob("*")):
        print(" ", p.name, p.stat().st_size)


if __name__ == "__main__":
    main()
