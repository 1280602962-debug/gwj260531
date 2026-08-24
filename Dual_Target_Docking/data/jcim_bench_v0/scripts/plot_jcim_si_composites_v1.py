"""Extended JCIM figures from frozen CSVs (same style as Fig 1–5).

Main: Fig 6 wrong-pocket paradox · Fig 7 confound tests.
SI: S1 protocol knobs · S2 equal-relation supply + sampling shift.
Called from plot_jcim_article_figures_v1.py so one command regenerates everything.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from jcim_figure_style import (
    C,
    FS_AXIS,
    OUT,
    PAIR_ORDER,
    ROOT,
    panel_label,
    save_all,
)

DATA = ROOT / "data"

TICK = ["EGFR/\nHER2", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"]
HOLD_TICK = ["AChE/\nBChE", "PIK3CA/\nPIK3CB", "PIK3CA/\nmTOR"]
PAIR_COLOR = {
    "EGFR/HER2": C["egfr"],
    "AChE/BChE": C["rtm"],
    "PIK3CA/PIK3CB": C["gnina"],
    "PIK3CA/mTOR": C["vina"],
}
RULES = ["theta_5.5", "theta_6.0", "theta_6.5", "strict_6.5_5.5"]
RULE_LAB = ["θ=5.5", "θ=6.0", "θ=6.5", "strict\n6.5/5.5"]
S12_PAIRS = ["PIK3CA/MTOR", "ACHE/BCHE", "PIK3CA/PIK3CB", "EGFR/HER2"]
HOLD_PAIRS = ["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]


def legend_below(ax, ncol=2, fontsize=5.5, y=-0.22):
    """Place the legend under the x-axis so it cannot cover data or the panel title."""
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        fontsize=fontsize,
        frameon=False,
        borderaxespad=0.0,
        handlelength=1.2,
        columnspacing=0.75,
        handletextpad=0.35,
    )


def fnum(x) -> float:
    return float(x)


def _auroc(pos, neg) -> float:
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def extend_load(D: dict, _read, provenance: dict) -> dict:
    """Attach SI-only frozen tables onto the main-figure data dict."""
    theta = _read(DATA / "jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv")
    gnina = _read(DATA / "jcim_bench_v0/tables/gnina_pocket_matched_mode01_vs_best9_k4_v1.csv")
    pm110 = _read(DATA / "jcim_strengthen_t0t1_v0/tables/pm110_vs_pm48_pocket_matched_v1.csv")
    ml = _read(DATA / "jcim_strengthen_t0t1_v0/tables/ligand_ml_baseline_scaffold_cv_v1.csv")
    cov = _read(DATA / "jcim_strengthen_t0t1_v0/tables/covariate_adjusted_v1.csv")
    matched = _read(DATA / "jcim_strengthen_t0t1_v0/tables/matched_subset_directional_v1.csv")
    enrich = _read(DATA / "jcim_strengthen_t0t1_v0/tables/single_target_enrichment_v1.csv")
    hold_m = _read(DATA / "jcim_holdout_v0/tables/holdout_matched_wrong_pocket_summary_v1.csv")
    hold_ps = _read(DATA / "jcim_holdout_v0/tables/holdout_vs_main_potency_size_v1.csv")

    D["theta_all"] = theta
    D["gnina_mode"] = {r["pair"]: r for r in gnina if r["channel"] == "mode01"}
    D["gnina9"] = {r["pair"]: r for r in gnina if r["channel"] == "best9"}
    D["pm110"] = {(r["panel"], r["arm"]): r for r in pm110}
    D["ml_scaf"] = {(r["pair"], r["contrast"]): r for r in ml}
    D["cov"] = {(r["pair"], r["contrast"]): r for r in cov}
    D["matched"] = {(r["pair"], r["match_type"]): r for r in matched}
    D["enrich"] = {r["receptor"]: r for r in enrich}
    D["hold_match"] = {(r["pair"], r["family"], r["aggregation"]): r for r in hold_m}
    D["hold_ps"] = hold_ps
    D["e8"] = _pm48_e8()
    D["contact"] = _parse_contact()
    delta = _read(DATA / "jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv")
    desc_d = _read(DATA / "jcim_strengthen_t0t1_v0/tables/pocket_matched_vs_best_descriptor_delta_v1.csv")
    ml_both = _read(DATA / "jcim_strengthen_t0t1_v0/tables/ligand_ml_scaffold_vs_random_v1.csv")
    D["delta_wp"] = {(r["pair"], r["set"]): r for r in delta}
    D["delta_desc"] = {r["pair"]: r for r in desc_d}
    D["ml_both"] = {(r["pair"], r["contrast"]): r for r in ml_both}
    provenance.setdefault("source_files", {})
    provenance["source_files"][str((DATA / "pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv").relative_to(ROOT))] = len(
        list(csv.DictReader((DATA / "pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv").open()))
    )
    return D


def _pm48_e8() -> dict:
    """Pocket-matched Vina on PM48 at exhaustiveness=8, from frozen score files."""
    panel = {
        r["panel_id"]: r
        for r in csv.DictReader((DATA / "pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv").open())
    }
    e8 = {
        r["ligand"]: r
        for r in csv.DictReader((DATA / "pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv").open())
    }
    recs = []
    for lig, meta in panel.items():
        r = e8.get(lig)
        if not r:
            continue
        a, b = r.get("4L23_affinity_E8"), r.get("4JT6_affinity_E8")
        if not a or not b:
            continue
        recs.append({"cls": meta["class"], "A": -float(a), "B": -float(b)})
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = _auroc([r["B"] for r in D], [r["B"] for r in A])
    db = _auroc([r["A"] for r in D], [r["A"] for r in B])
    return {"da": da, "db": db, "summary_min": min(da, db), "nD": len(D), "nA": len(A), "nB": len(B)}


def _parse_contact() -> dict:
    path = DATA / "jcim_holdout_v0/analysis/wrong_pocket_contact_v1_output.txt"
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("HO"):
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            continue
        out[parts[1]] = {"A": float(parts[2]), "B": float(parts[3]), "min": min(float(parts[2]), float(parts[3]))}
    return out


def fig_s2_protocol(D: dict, P: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.70))

    # A: θ grid
    ax = axes[0, 0]
    panel_label(ax, "A", x=-0.18, y=1.04)
    plotted = {}
    x = np.arange(len(RULES))
    for pair in PAIR_ORDER:
        ys, lo, hi, under = [], [], [], []
        for rule in RULES:
            r = next(row for row in D["theta_all"] if row["pair"] == pair and row["label_rule"] == rule)
            ys.append(fnum(r["pocket_matched_summary_min"]))
            lo.append(fnum(r["ci_lo"]))
            hi.append(fnum(r["ci_hi"]))
            under.append(r["underpowered"] == "1")
        yerr = np.vstack([np.array(ys) - np.array(lo), np.array(hi) - np.array(ys)])
        ax.errorbar(
            x, ys, yerr=yerr, fmt="-o", color=PAIR_COLOR[pair], ecolor=PAIR_COLOR[pair],
            elinewidth=1.0, capsize=1.8, markersize=4.5, lw=1.1, zorder=4, label=pair,
        )
        for xi, yi, u in zip(x, ys, under):
            if u:
                ax.scatter([xi], [yi], s=28, facecolors="white", edgecolors=PAIR_COLOR[pair], linewidths=1.0, zorder=5)
        plotted[pair] = {"y": ys, "under": under}
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(RULE_LAB, fontsize=6.5)
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0.05, 0.95)
    ax.set_title("Label-threshold grid", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.24)
    ax.text(0.02, 0.04, "open marker = underpowered", transform=ax.transAxes, ha="left", fontsize=6.0, color="#666666")
    P["s2A"] = plotted

    # B: GNINA mode01 vs best9 vs Vina
    ax = axes[0, 1]
    panel_label(ax, "B", x=-0.18, y=1.04)
    x = np.arange(len(PAIR_ORDER))
    w = 0.24
    vina = [fnum(D["theta6"][p]["pocket_matched_summary_min"]) for p in PAIR_ORDER]
    m01 = [fnum(D["gnina_mode"][p]["summary_min"]) for p in PAIR_ORDER]
    b9 = [fnum(D["gnina9"][p]["summary_min"]) for p in PAIR_ORDER]
    ax.bar(x - w, vina, w, color=C["vina"], label="Vina (primary)", zorder=3)
    ax.bar(x, m01, w, color=C["gnina"], label="GNINA mode01", zorder=3)
    ax.bar(x + w, b9, w, color="#005F73", label="GNINA best-of-9", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.5)
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0, 1.0)
    ax.set_title("GNINA pose coverage", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=3, y=-0.24)
    P["s2B"] = {"vina": vina, "mode01": m01, "best9": b9}

    # C: PM48 vs PM110
    ax = axes[1, 0]
    panel_label(ax, "C", x=-0.18, y=1.04)
    arms = [("vina", "Vina", C["vina"]), ("rtm", "RTMScore", C["rtm"]), ("gnina_best9", "GNINA best-of-9", C["gnina"])]
    x = np.arange(len(arms))
    w = 0.32
    y48, y110 = [], []
    for arm, _, _ in arms:
        y48.append(fnum(D["pm110"][("PM48", arm)]["summary_min"]))
        y110.append(fnum(D["pm110"][("PM110", arm)]["summary_min"]))
    ax.bar(x - w / 2, y48, w, color=C["main"], label="PM48 (n=18/14/12)", zorder=3)
    ax.bar(x + w / 2, y110, w, color=C["holdout"], label="PM110 (n=30/30/30)", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([a[1] for a in arms], fontsize=7)
    ax.set_ylabel("Pocket-matched summary_min")
    ax.set_ylim(0, 1.0)
    ax.set_title("PIK3CA/mTOR panel expansion", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.24)
    P["s2C"] = {"PM48": y48, "PM110": y110}

    # D: E8 vs E16 + single-target enrichment
    ax = axes[1, 1]
    panel_label(ax, "D", x=-0.18, y=1.04)
    e16 = fnum(D["theta6"]["PIK3CA/mTOR"]["pocket_matched_summary_min"])
    e8 = D["e8"]["summary_min"]
    ax.bar([0], [e16], 0.55, color=C["vina"], zorder=3, label="E=16 (primary)")
    ax.bar([1], [e8], 0.55, color="#56B4E9", zorder=3, label="E=8")
    ax.bar([2.6], [fnum(D["enrich"]["4L23"]["auroc"])], 0.55, color=C["a_only"], zorder=3, label="4L23 single-target")
    ax.bar([3.6], [fnum(D["enrich"]["4JT6"]["auroc"])], 0.55, color=C["b_only"], zorder=3, label="4JT6 single-target")
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks([0, 1, 2.6, 3.6])
    ax.set_xticklabels(["E=16\n(main)", "E=8", "4L23\nPIK3CA", "4JT6\nmTOR"], fontsize=6.5)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.0)
    ax.set_title("Exhaustiveness + single-target EF", fontsize=FS_AXIS, pad=3)
    ax.text(0, e16 + 0.04, f"{e16:.3f}", ha="center", fontsize=6.0)
    ax.text(1, e8 + 0.04, f"{e8:.3f}", ha="center", fontsize=6.0)
    ax.text(2.6, fnum(D["enrich"]["4L23"]["auroc"]) + 0.04, f"EF1%={fnum(D['enrich']['4L23']['EF_1pct']):.2f}", ha="center", fontsize=5.5, color="#555555")
    ax.text(3.6, fnum(D["enrich"]["4JT6"]["auroc"]) + 0.04, f"EF1%={fnum(D['enrich']['4JT6']['EF_1pct']):.2f}", ha="center", fontsize=5.5, color="#555555")
    legend_below(ax, ncol=2, y=-0.24)
    P["s2D"] = {
        "e16": e16,
        "e8": e8,
        "e8_da": D["e8"]["da"],
        "e8_db": D["e8"]["db"],
        "e8_n": (D["e8"]["nD"], D["e8"]["nA"], D["e8"]["nB"]),
        "4L23": fnum(D["enrich"]["4L23"]["auroc"]),
        "4JT6": fnum(D["enrich"]["4JT6"]["auroc"]),
        "EF_4L23": fnum(D["enrich"]["4L23"]["EF_1pct"]),
        "EF_4JT6": fnum(D["enrich"]["4JT6"]["EF_1pct"]),
    }

    fig.subplots_adjust(wspace=0.38, hspace=0.62, left=0.09, right=0.98, top=0.94, bottom=0.14)
    save_all(fig, "FigS1_protocol_sensitivity")
    plt.close(fig)


def fig6_wrong_pocket(D: dict, P: dict) -> None:
    """Main-text Fig 6: matched vs wrong-pocket on main, holdout, matched holdout, contacts."""
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.70))
    w = 0.36

    ax = axes[0, 0]
    panel_label(ax, "A", x=-0.16, y=1.04)
    x = np.arange(len(PAIR_ORDER))
    matched, wrong = [], []
    for p in PAIR_ORDER:
        matched.append(fnum(D["theta6"][p]["pocket_matched_summary_min"]))
        wrong.append(fnum(D["pm_by"][(p, "wrong_pocket_control_vina")]["summary_min"]))
    ax.bar(x - w / 2, matched, w, color=C["vina"], label="Pocket-matched", zorder=3)
    ax.bar(x + w / 2, wrong, w, color="#999999", label="Wrong-pocket", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.5)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0, 1.05)
    ax.set_title("Main panels", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.24)
    P["fig6A_matched"] = matched
    P["fig6A_wrong"] = wrong
    P["siA_matched"] = matched
    P["siA_wrong"] = wrong

    ax = axes[0, 1]
    panel_label(ax, "B", x=-0.16, y=1.04)
    x = np.arange(len(HOLD_PAIRS))
    hm, hw = [], []
    for p in HOLD_PAIRS:
        hm.append(fnum(D["hold_pm"][p]["summary_min"]))
        hw.append(fnum(D["hold_wp"][p]["summary_min"]))
    ax.bar(x - w / 2, hm, w, color=C["vina"], label="Pocket-matched", zorder=3)
    ax.bar(x + w / 2, hw, w, color="#999999", label="Wrong-pocket", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(HOLD_TICK, fontsize=6.5)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0, 1.05)
    ax.set_title("Unused-pool holdout", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.24)
    ax.text(0.02, 0.03, "EGFR/HER2 has no holdout", transform=ax.transAxes, fontsize=5.5, color="#666666")
    P["fig6B_matched"] = hm
    P["fig6B_wrong"] = hw
    P["siB_matched"] = hm
    P["siB_wrong"] = hw

    ax = axes[1, 0]
    panel_label(ax, "C", x=-0.16, y=1.04)
    fams = ["unmatched", "potency_matched", "size_matched"]
    plotted = []
    x = np.arange(len(HOLD_PAIRS))
    for i, pair in enumerate(HOLD_PAIRS):
        for j, fam in enumerate(fams):
            m = D["hold_match"][(pair, fam, "pocket_matched")]
            wrow = D["hold_match"][(pair, fam, "wrong_pocket")]
            my, wy = fnum(m["summary_min"]), fnum(wrow["summary_min"])
            xx = i + (j - 1) * 0.22
            ax.plot([xx, xx], [my, wy], color="#CCCCCC", lw=0.8, zorder=2)
            ax.scatter([xx], [my], s=28, color=C["vina"], zorder=4, marker="o")
            ax.scatter([xx], [wy], s=28, color="#999999", zorder=4, marker="s")
            plotted.append({"pair": pair, "family": fam, "matched": my, "wrong": wy, "n_min": int(m["n_min"])})
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(HOLD_TICK, fontsize=6.5)
    ax.set_ylabel("summary_min")
    ax.set_ylim(0.15, 1.05)
    ax.set_title("Holdout matching does not reverse", fontsize=FS_AXIS, pad=3)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=5.5, label="Pocket-matched"),
            Line2D([0], [0], marker="s", color="#999999", ls="none", ms=5.5, label="Wrong-pocket"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=2,
        fontsize=5.5,
        frameon=False,
        borderaxespad=0.0,
    )
    ax.text(0.02, 0.03, "offsets: unmatched / potency / size", transform=ax.transAxes, fontsize=5.5, color="#666666")
    P["fig6C"] = plotted
    P["s4B"] = plotted

    ax = axes[1, 1]
    panel_label(ax, "D", x=-0.16, y=1.04)
    x = np.arange(len(HOLD_PAIRS))
    bw = 0.24
    cA, cB, vw = [], [], []
    for p in HOLD_PAIRS:
        cA.append(D["contact"][p]["A"])
        cB.append(D["contact"][p]["B"])
        vw.append(fnum(D["hold_wp"][p]["summary_min"]))
    ax.bar(x - bw, cA, bw, color=C["a_only"], label="contact pocket A", zorder=3)
    ax.bar(x, cB, bw, color=C["b_only"], label="contact pocket B", zorder=3)
    ax.bar(x + bw, vw, bw, color="#999999", label="Vina wrong-pocket", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(HOLD_TICK, fontsize=6.5)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.05)
    ax.set_title("Scoring-free contact count (not PLIF)", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=3, y=-0.28)
    P["fig6D"] = {"contact_A": cA, "contact_B": cB, "vina_wrong": vw}
    P["s4C"] = P["fig6D"]

    fig.subplots_adjust(wspace=0.38, hspace=0.64, left=0.10, right=0.98, top=0.94, bottom=0.15)
    save_all(fig, "Fig6_wrong_pocket_paradox")
    plt.close(fig)


def fig_s3_confounds(D: dict, P: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.80))

    # A: ECFP4 scaffold CV vs Vina (both arms)
    ax = axes[0, 0]
    panel_label(ax, "A", x=-0.18, y=1.04)
    x = np.arange(len(PAIR_ORDER))
    w = 0.18
    da_d, da_m, db_d, db_m = [], [], [], []
    for p in PAIR_ORDER:
        da = D["ml_scaf"][(p, "D_vs_A")]
        db = D["ml_scaf"][(p, "D_vs_B")]
        da_d.append(fnum(da["auroc_dock_pocket_matched"]))
        da_m.append(fnum(da["auroc_ml"]))
        db_d.append(fnum(db["auroc_dock_pocket_matched"]))
        db_m.append(fnum(db["auroc_ml"]))
    ax.bar(x - 1.5 * w, da_d, w, color=C["vina"], label="Vina D/A", zorder=3)
    ax.bar(x - 0.5 * w, da_m, w, color="#56B4E9", label="ECFP4 D/A", zorder=3)
    ax.bar(x + 0.5 * w, db_d, w, color=C["a_only"], label="Vina D/B", zorder=3)
    ax.bar(x + 1.5 * w, db_m, w, color=C["desc"], label="ECFP4 D/B", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.5)
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.15)
    ax.set_title("ECFP4 GroupKFold vs Vina", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.24)
    ax.text(0.98, 0.04, "chemotype–label association", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.0, color="#666666")
    P["s3A"] = {"da_dock": da_d, "da_ml": da_m, "db_dock": db_d, "db_ml": db_m}

    # B: all trivial descriptors vs pocket-matched Vina
    ax = axes[0, 1]
    panel_label(ax, "B", x=-0.18, y=1.04)
    desc_arms = [("heavy", C["thick"], "s"), ("mw", "#56B4E9", "D"), ("clogp", C["desc"], "^"), ("tpsa", C["a_only"], "v")]
    y_gap = 1.05
    plotted = {}
    for i, pair in enumerate(PAIR_ORDER):
        y0 = -i * y_gap
        t = D["theta6"][pair]
        vy, vlo, vhi = fnum(t["pocket_matched_summary_min"]), fnum(t["ci_lo"]), fnum(t["ci_hi"])
        ax.errorbar(vy, y0 + 0.32, xerr=[[vy - vlo], [vhi - vy]], fmt="o", color=C["vina"], ecolor=C["vina"], elinewidth=1.3, capsize=2.0, markersize=6, zorder=4)
        row = {"vina": vy}
        for k, (arm, col, mk) in enumerate(desc_arms):
            dr = D["forest_by"][(pair, arm)]
            dy, dlo, dhi = fnum(dr["summary_min"]), fnum(dr["ci_lo"]), fnum(dr["ci_hi"])
            yy = y0 + 0.12 - k * 0.14
            ax.errorbar(dy, yy, xerr=[[dy - dlo], [dhi - dy]], fmt=mk, color=col, ecolor=col, elinewidth=0.9, capsize=1.6, markersize=4.5, zorder=4)
            row[arm] = dy
        plotted[pair] = row
    ax.axvline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_yticks([-i * y_gap for i in range(4)])
    ax.set_yticklabels(PAIR_ORDER, fontsize=6.5)
    ax.set_xlabel("summary_min")
    ax.set_xlim(0.05, 0.95)
    ax.set_ylim(-3 * y_gap - 0.45, 0.55)
    ax.set_title("Vina vs all trivial descriptors", fontsize=FS_AXIS, pad=3)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=C["vina"], ls="none", ms=6, label="Vina"),
            Line2D([0], [0], marker="s", color=C["thick"], ls="none", ms=5, label="heavy"),
            Line2D([0], [0], marker="D", color="#56B4E9", ls="none", ms=5, label="MW"),
            Line2D([0], [0], marker="^", color=C["desc"], ls="none", ms=5, label="cLogP"),
            Line2D([0], [0], marker="v", color=C["a_only"], ls="none", ms=5, label="TPSA"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=5,
        fontsize=5.5,
        frameon=False,
        borderaxespad=0.0,
        columnspacing=0.7,
        handletextpad=0.25,
    )
    P["s3B"] = plotted

    # C: covariate-adjusted D vs B (weak arm)
    ax = axes[1, 0]
    panel_label(ax, "C", x=-0.18, y=1.04)
    x = np.arange(len(PAIR_ORDER))
    w = 0.32
    only, adj, ors = [], [], []
    for p in PAIR_ORDER:
        r = D["cov"][(p, "D_vs_B_only")]
        only.append(fnum(r["auroc_score_only"]))
        adj.append(fnum(r["auroc_score_plus_covariates"]))
        ors.append(fnum(r["or_score"]))
    ax.bar(x - w / 2, only, w, color=C["vina"], label="Vina only", zorder=3)
    ax.bar(x + w / 2, adj, w, color=C["desc"], label="Vina + heavy + TPSA", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.5)
    ax.set_ylabel("D vs B_only AUROC")
    ax.set_ylim(0, 1.18)
    ax.set_title("Covariate-adjusted weak arm", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.24)
    for i, o in enumerate(ors):
        ax.text(i, max(only[i], adj[i]) + 0.04, f"OR={o:.2f}", ha="center", fontsize=5.5, color="#555555")
    P["s3C"] = {"only": only, "adj": adj, "or": ors}

    # D: potency/size-matched D vs B
    ax = axes[1, 1]
    panel_label(ax, "D", x=-0.18, y=1.04)
    kinds = ["full_panel_pocket_matched", "potency_matched_D_vs_B", "size_matched_D_vs_B"]
    labs = ["full panel", "potency |Δp|≤0.5", "size |Δheavy|≤2"]
    cols = [C["vina"], C["desc"], C["rtm"]]
    x = np.arange(len(PAIR_ORDER))
    plotted = {k: [] for k in kinds}
    for i, pair in enumerate(PAIR_ORDER):
        for j, (k, col) in enumerate(zip(kinds, cols)):
            r = D["matched"][(pair, k)]
            y, lo, hi = fnum(r["auroc_contrast"]), fnum(r["ci_lo"]), fnum(r["ci_hi"])
            plotted[k].append(y)
            ax.errorbar(
                i + (j - 1) * 0.22, y, yerr=[[y - lo], [hi - y]],
                fmt="o", color=col, ecolor=col, elinewidth=1.1, capsize=2.0, markersize=5.5, zorder=4,
            )
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.5)
    ax.set_ylabel("D vs B_only AUROC")
    ax.set_ylim(0.10, 1.15)
    ax.set_title("Matched-subset weak arm", fontsize=FS_AXIS, pad=3)
    ax.legend(
        handles=[Line2D([0], [0], marker="o", color=c, ls="none", ms=5.5, label=l) for c, l in zip(cols, labs)],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=3,
        fontsize=5.0,
        frameon=False,
        borderaxespad=0.0,
        columnspacing=0.6,
        handletextpad=0.25,
    )
    P["s3D"] = plotted

    fig.subplots_adjust(wspace=0.40, hspace=0.68, left=0.11, right=0.98, top=0.94, bottom=0.16)
    save_all(fig, "Fig7_confound_anatomy")
    plt.close(fig)


def fig_s2_supply_sampling(D: dict, P: dict) -> None:
    """SI: leftover S4 panels that did not move into Fig 6."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.85))

    ax = axes[0]
    panel_label(ax, "A", x=-0.16, y=1.04)
    series = [
        ("ChEMBL_cache", "pChEMBL", C["vina"], "ChEMBL"),
        ("BindingDB", "equal_only", C["desc"], "BDB equal"),
        ("PubChem", "equal_only", C["rtm"], "PubChem equal"),
        ("BindingDB", "as_is", "#999999", "BDB as_is"),
        ("PubChem", "as_is", C["metal"], "PubChem as_is"),
    ]
    x = np.arange(4)
    bw = 0.15
    plotted = {s[3]: [] for s in series}
    for j, (src, rule, col, lab) in enumerate(series):
        vals = []
        for name in S12_PAIRS:
            r = next(row for row in D["s12"] if row["pair"] == name and row["source"] == src and row["rule"] == rule)
            vals.append(fnum(r["min_strict_hardneg"]))
        plotted[lab] = vals
        ax.bar(x + (j - 2) * bw, vals, bw, color=col, label=lab, zorder=3)
    ax.axhline(50, color=C["ink"], ls="--", lw=0.9, zorder=2)
    ax.axhline(20, color="#888888", ls=":", lw=0.9, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["PIK3CA/\nmTOR", "AChE/\nBChE", "PIK3CA/\nPIK3CB", "EGFR/\nHER2"], fontsize=6.5)
    ax.set_ylabel("min strict hard-negatives")
    ax.set_ylim(0, 170)
    ax.set_title("Count-level supply (no docking)", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=5, fontsize=5.0, y=-0.30)
    P["s4A"] = plotted

    ax = axes[1]
    panel_label(ax, "B", x=-0.16, y=1.04)
    feats = [("dual", "pA"), ("A_only", "pA"), ("B_only", "pB")]
    feat_lab = ["dual pA", "A_only pA", "B_only pB"]
    feat_col = [C["vina"], C["a_only"], C["b_only"]]
    x = np.arange(len(HOLD_PAIRS))
    bw = 0.24
    plotted = {lab: [] for lab in feat_lab}
    ps = {(r["pair"], r["cls"], r["feature"]): r for r in D["hold_ps"]}
    for j, ((cls, feat), lab, col) in enumerate(zip(feats, feat_lab, feat_col)):
        vals = [fnum(ps[(p, cls, feat)]["mean_delta_holdout_minus_main"]) for p in HOLD_PAIRS]
        plotted[lab] = vals
        ax.bar(x + (j - 1) * bw, vals, bw, color=col, label=lab, zorder=3)
    ax.axhline(0.0, color=C["ink"], lw=0.8, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(HOLD_TICK, fontsize=6.5)
    ax.set_ylabel("holdout − main (mean pChEMBL)")
    ax.set_ylim(-2.2, 1.2)
    ax.set_title("Holdout sampling shift (not a reversal)", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=3, y=-0.28)
    P["s4D"] = plotted

    fig.subplots_adjust(wspace=0.38, left=0.10, right=0.98, top=0.88, bottom=0.26)
    save_all(fig, "FigS2_equal_relation_and_sampling")
    plt.close(fig)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def _delta_forest(ax, rows, ylabels, xlabel):
    """Horizontal forest of paired Δ with 95% CIs. Gray if the CI includes 0."""
    plotted = []
    for i, (r, lab) in enumerate(zip(rows, ylabels)):
        x = fnum(r["delta"])
        lo, hi = fnum(r["lo"]), fnum(r["hi"])
        excl = bool(r["excl"])
        color = C["vina"] if excl else "#999999"
        ax.plot([lo, hi], [i, i], color=color, lw=1.7, zorder=3, solid_capstyle="round")
        ax.scatter([x], [i], s=38, color=color, zorder=4, marker="o", edgecolors="none")
        plotted.append({"label": lab, "x": x, "lo": lo, "hi": hi, "excl": excl})
    ax.axvline(0.0, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel(xlabel)
    ax.invert_yaxis()
    return plotted


def fig_s3_paired_deltas(D: dict, P: dict) -> None:
    """SI Fig S3: paired Δ forests (new). Does not reuse Fig 6 AUROC bars."""
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.70))

    ax = axes[0, 0]
    panel_label(ax, "A", x=-0.28, y=1.04)
    main_rows = []
    for p in PAIR_ORDER:
        r = D["delta_wp"][(p, "main_panel")]
        main_rows.append(
            {
                "delta": r["delta_matched_minus_wrong"],
                "lo": r["delta_ci_lo"],
                "hi": r["delta_ci_hi"],
                "excl": _truthy(r["ci_excludes_zero"]),
                "matched": r["matched_summary_min"],
                "wrong": r["wrong_summary_min"],
            }
        )
    P["s3newA"] = _delta_forest(
        ax,
        main_rows,
        ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"],
        "Δ summary_min (matched − wrong)",
    )
    for rec, src in zip(P["s3newA"], main_rows):
        rec["matched"] = fnum(src["matched"])
        rec["wrong"] = fnum(src["wrong"])
    ax.set_xlim(-0.20, 0.38)
    ax.set_title("Main panels", fontsize=FS_AXIS, pad=3)
    ax.text(0.02, 0.04, "blue: CI excludes 0", transform=ax.transAxes, fontsize=5.5, color=C["vina"])

    ax = axes[0, 1]
    panel_label(ax, "B", x=-0.28, y=1.04)
    hold_rows = []
    for p in HOLD_PAIRS:
        r = D["delta_wp"][(p, "unused_pool_holdout")]
        hold_rows.append(
            {
                "delta": r["delta_matched_minus_wrong"],
                "lo": r["delta_ci_lo"],
                "hi": r["delta_ci_hi"],
                "excl": _truthy(r["ci_excludes_zero"]),
                "matched": r["matched_summary_min"],
                "wrong": r["wrong_summary_min"],
            }
        )
    P["s3newB"] = _delta_forest(
        ax,
        hold_rows,
        ["AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"],
        "Δ summary_min (matched − wrong)",
    )
    for rec, src in zip(P["s3newB"], hold_rows):
        rec["matched"] = fnum(src["matched"])
        rec["wrong"] = fnum(src["wrong"])
    ax.set_xlim(-0.35, 0.18)
    ax.set_title("Unused-pool holdout", fontsize=FS_AXIS, pad=3)
    ax.text(0.02, 0.04, "EGFR/HER2 has no holdout", transform=ax.transAxes, fontsize=5.5, color="#666666")

    ax = axes[1, 0]
    panel_label(ax, "C", x=-0.28, y=1.04)
    desc_rows = []
    for p in PAIR_ORDER:
        r = D["delta_desc"][p]
        desc_rows.append(
            {
                "delta": r["delta_vina_minus_descriptor"],
                "lo": r["delta_ci_lo"],
                "hi": r["delta_ci_hi"],
                "excl": _truthy(r["ci_excludes_zero"]),
                "vina": r["vina_summary_min"],
                "desc": r["descriptor_summary_min"],
                "arm": r["best_descriptor"],
            }
        )
    P["s3newC"] = _delta_forest(
        ax,
        desc_rows,
        ["EGFR/HER2  (cLogP)", "AChE/BChE  (TPSA)", "PIK3CA/PIK3CB  (heavy)", "PIK3CA/mTOR  (heavy)"],
        "Δ summary_min (Vina − best descriptor)",
    )
    for rec, src in zip(P["s3newC"], desc_rows):
        rec["vina"] = fnum(src["vina"])
        rec["desc"] = fnum(src["desc"])
        rec["arm"] = src["arm"]
    ax.set_xlim(-0.40, 0.50)
    ax.set_title("Pocket-matched vs descriptor", fontsize=FS_AXIS, pad=3)

    ax = axes[1, 1]
    panel_label(ax, "D", x=-0.18, y=1.04)
    x = np.arange(len(PAIR_ORDER))
    w = 0.18
    sc_da, rd_da, sc_db, rd_db = [], [], [], []
    for p in PAIR_ORDER:
        da = D["ml_both"][(p, "D_vs_A")]
        db = D["ml_both"][(p, "D_vs_B")]
        sc_da.append(fnum(da["auroc_scaffold_GroupKFold"]))
        rd_da.append(fnum(da["auroc_random_StratifiedKFold"]))
        sc_db.append(fnum(db["auroc_scaffold_GroupKFold"]))
        rd_db.append(fnum(db["auroc_random_StratifiedKFold"]))
    ax.bar(x - 1.5 * w, sc_da, w, color=C["vina"], label="scaffold D/A", zorder=3)
    ax.bar(x - 0.5 * w, rd_da, w, color="#56B4E9", label="random D/A", zorder=3)
    ax.bar(x + 0.5 * w, sc_db, w, color=C["a_only"], label="scaffold D/B", zorder=3)
    ax.bar(x + 1.5 * w, rd_db, w, color=C["desc"], label="random D/B", zorder=3)
    ax.axhline(0.5, color=C["chance"], ls="--", lw=0.9, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK, fontsize=6.5)
    ax.set_ylabel("ECFP4 AUROC")
    ax.set_ylim(0, 1.15)
    ax.set_title("Scaffold vs random split (leakage check)", fontsize=FS_AXIS, pad=3)
    legend_below(ax, ncol=2, y=-0.26)
    dlt = [rd_da[i] - sc_da[i] for i in range(4)] + [rd_db[i] - sc_db[i] for i in range(4)]
    P["s3newD"] = {
        "sc_da": sc_da,
        "rd_da": rd_da,
        "sc_db": sc_db,
        "rd_db": rd_db,
        "delta_random_minus_scaffold": dlt,
        "mean_delta": float(np.mean(dlt)),
    }

    fig.subplots_adjust(wspace=0.42, hspace=0.58, left=0.16, right=0.98, top=0.94, bottom=0.14)
    save_all(fig, "FigS3_paired_delta_bootstrap")
    plt.close(fig)


def draw_all(D: dict, provenance: dict) -> None:
    P = provenance.setdefault("plotted", {})
    fig6_wrong_pocket(D, P)
    fig_s3_confounds(D, P)
    fig_s2_protocol(D, P)
    fig_s2_supply_sampling(D, P)
    fig_s3_paired_deltas(D, P)


def verify_si(D: dict, provenance: dict, errors: list) -> None:
    def eq(a, b, tol=1e-9, msg=""):
        if abs(float(a) - float(b)) > tol:
            errors.append(f"{msg}: plotted {a} != source {b}")

    P = provenance["plotted"]
    # S2A θ=6.0 matches Table 2
    for pair in PAIR_ORDER:
        eq(P["s2A"][pair]["y"][1], D["theta6"][pair]["pocket_matched_summary_min"], msg=f"s2A {pair} theta6")
    eq(P["s2A"]["EGFR/HER2"]["y"][3], 0.3242, msg="s2A EGFR strict checksum")
    eq(P["s2A"]["PIK3CA/mTOR"]["y"][0], 0.5017, msg="s2A PM theta5.5 checksum")
    eq(P["s2A"]["AChE/BChE"]["y"][0], 0.6058, msg="s2A AChE theta5.5")
    if P["s2A"]["EGFR/HER2"]["under"][0] or P["s2A"]["EGFR/HER2"]["under"][1]:
        errors.append("s2A EGFR theta_5.5/6.0 should not be underpowered")
    if not P["s2A"]["EGFR/HER2"]["under"][3]:
        errors.append("s2A EGFR strict should be underpowered")
    if not P["s2A"]["PIK3CA/mTOR"]["under"][0] or not P["s2A"]["PIK3CA/mTOR"]["under"][3]:
        errors.append("s2A PM theta_5.5 and strict should be underpowered")
    if any(P["s2A"]["AChE/BChE"]["under"]):
        errors.append("s2A AChE should never be underpowered")
    # θ=5.5 PM is underpowered and is NOT the highest point; do not caption it as a stable ranking
    if P["s2A"]["PIK3CA/mTOR"]["y"][0] >= P["s2A"]["AChE/BChE"]["y"][0]:
        errors.append("s2A theta_5.5: underpowered PM should sit below AChE")
    if P["s2A"]["PIK3CA/mTOR"]["y"][1] <= max(P["s2A"][p]["y"][1] for p in PAIR_ORDER if p != "PIK3CA/mTOR"):
        errors.append("s2A theta_6.0: PM should be the highest point estimate")

    eq(P["s2B"]["mode01"][0], 0.327, msg="s2B EGFR mode01")
    eq(P["s2B"]["best9"][0], 0.2902, msg="s2B EGFR best9")
    eq(P["s2B"]["best9"][2], 0.5332, msg="s2B PIK3CB best9")
    eq(P["s2B"]["vina"][3], 0.6921, msg="s2B PM vina")
    d_b9_m01 = [P["s2B"]["best9"][i] - P["s2B"]["mode01"][i] for i in range(4)]
    if min(d_b9_m01) < -0.041 or max(d_b9_m01) > 0.081:
        errors.append(f"s2B best9−mode01 outside −0.04..+0.08: {d_b9_m01}")

    eq(P["s2C"]["PM48"][0], 0.6921, msg="s2C PM48 vina")
    eq(P["s2C"]["PM110"][0], 0.6483, msg="s2C PM110 vina")
    eq(P["s2C"]["PM48"][2], 0.6548, msg="s2C PM48 gnina9")
    eq(P["s2C"]["PM110"][2], 0.6133, msg="s2C PM110 gnina9")

    eq(P["s2D"]["e16"], 0.6921, msg="s2D E16")
    eq(P["s2D"]["e8"], 0.6597, tol=5e-4, msg="s2D E8 checksum")
    eq(P["s2D"]["e8_da"], 0.7540, tol=5e-4, msg="s2D E8 DA")
    eq(P["s2D"]["e8_db"], 0.6597, tol=5e-4, msg="s2D E8 DB")
    if P["s2D"]["e8_n"] != (18, 14, 12):
        errors.append(f"s2D E8 n {P['s2D']['e8_n']} != (18,14,12)")
    eq(P["s2D"]["4L23"], 0.6027, msg="s2D 4L23 auroc")
    eq(P["s2D"]["4JT6"], 0.6288, msg="s2D 4JT6 auroc")
    eq(P["s2D"]["EF_4L23"], 2.0408, msg="s2D 4L23 EF")
    eq(P["s2D"]["EF_4JT6"], 2.0, msg="s2D 4JT6 EF")

    eq(P["s3A"]["db_ml"][0], 0.8527, msg="s3A EGFR ECFP D/B")
    eq(P["s3A"]["db_dock"][0], 0.4297, msg="s3A EGFR Vina D/B")
    eq(P["s3A"]["da_ml"][1], 0.9096, msg="s3A AChE ECFP D/A")
    eq(P["s3A"]["db_ml"][3], 0.8889, msg="s3A PM ECFP D/B")

    eq(P["s3B"]["EGFR/HER2"]["clogp"], 0.4821, msg="s3B EGFR clogp")
    eq(P["s3B"]["AChE/BChE"]["tpsa"], 0.7333, msg="s3B AChE tpsa")
    eq(P["s3B"]["PIK3CA/PIK3CB"]["heavy"], 0.6217, msg="s3B PIK3CB heavy")
    eq(P["s3B"]["PIK3CA/mTOR"]["heavy"], 0.463, msg="s3B PM heavy")
    eq(P["s3B"]["PIK3CA/mTOR"]["vina"], 0.6921, msg="s3B PM vina")

    eq(P["s3C"]["only"][0], 0.5703, msg="s3C EGFR score-only")
    if abs(P["s3C"]["only"][0] - 0.4297) < 1e-6:
        errors.append("s3C EGFR score-only must remain 0.5703, not Table 2 rank AUROC 0.4297")
    eq(P["s3C"]["adj"][1], 0.8069, msg="s3C AChE adj")
    eq(P["s3C"]["or"][1], 1.1755, msg="s3C AChE OR")
    eq(P["s3C"]["or"][3], 3.0807, msg="s3C PM OR")
    eq(P["s3C"]["only"][2], 0.5, msg="s3C PIK3CB score-only")

    eq(P["s3D"]["full_panel_pocket_matched"][0], 0.4297, msg="s3D EGFR full")
    eq(P["s3D"]["potency_matched_D_vs_B"][0], 0.4694, msg="s3D EGFR potency D/B")
    eq(P["s3D"]["size_matched_D_vs_B"][0], 0.5186, msg="s3D EGFR size D/B")
    eq(P["s3D"]["potency_matched_D_vs_B"][2], 0.4575, msg="s3D PIK3CB potency D/B")

    eq(P["s4A"]["ChEMBL"][3], 7, msg="s4A EGFR chembl")
    eq(P["s4A"]["BDB equal"][3], 31, msg="s4A EGFR bdb equal")
    eq(P["s4A"]["BDB as_is"][3], 85, msg="s4A EGFR bdb as_is")
    eq(P["s4A"]["PubChem as_is"][3], 88, msg="s4A EGFR pubchem as_is")
    eq(P["s4A"]["PubChem equal"][3], 30, msg="s4A EGFR pubchem equal")
    eq(P["s4A"]["BDB equal"][0], 76, msg="s4A PM bdb equal")

    # S4B: unmatched matches holdout table; potency AChE 0.5926 vs 0.642
    um = next(r for r in P["s4B"] if r["pair"] == "AChE/BChE" and r["family"] == "unmatched")
    eq(um["matched"], 0.6175, msg="s4B AChE unmatched matched")
    eq(um["wrong"], 0.6425, msg="s4B AChE unmatched wrong")
    pm = next(r for r in P["s4B"] if r["pair"] == "AChE/BChE" and r["family"] == "potency_matched")
    eq(pm["matched"], 0.5926, msg="s4B AChE potency matched")
    eq(pm["wrong"], 0.642, msg="s4B AChE potency wrong")
    pcb = next(r for r in P["s4B"] if r["pair"] == "PIK3CA/PIK3CB" and r["family"] == "potency_matched")
    eq(pcb["matched"], 0.3633, msg="s4B PIK3CB potency matched")
    eq(pcb["wrong"], 0.5625, msg="s4B PIK3CB potency wrong")
    for r in P["s4B"]:
        if r["wrong"] < r["matched"]:
            errors.append(f"s4B {r['pair']} {r['family']}: wrong should be >= matched")

    eq(P["s4C"]["contact_A"][0], 0.581, msg="s4C AChE contact A")
    eq(P["s4C"]["contact_B"][0], 0.706, msg="s4C AChE contact B")
    eq(P["s4C"]["contact_A"][2], 0.552, msg="s4C PM contact A")
    eq(P["s4C"]["contact_B"][2], 0.698, msg="s4C PM contact B")
    eq(P["s4C"]["contact_A"][1], 0.622, msg="s4C PIK3CB contact A")
    eq(P["s4C"]["contact_B"][1], 0.714, msg="s4C PIK3CB contact B")
    eq(P["s4C"]["vina_wrong"][0], 0.6425, msg="s4C AChE vina wrong")
    eq(P["s4C"]["vina_wrong"][1], 0.52, msg="s4C PIK3CB vina wrong")
    eq(P["s4C"]["vina_wrong"][2], 0.7875, msg="s4C PM vina wrong")
    if "fig6A_matched" not in P or "fig6D" not in P:
        errors.append("fig6 keys missing")
    else:
        for i, pair in enumerate(PAIR_ORDER):
            eq(P["fig6A_matched"][i], P["siA_matched"][i], msg=f"fig6A vs siA {pair}")
        eq(P["fig6D"]["contact_A"][0], P["s4C"]["contact_A"][0], msg="fig6D vs s4C")

    eq(P["s4D"]["dual pA"][2], -1.072, msg="s4D PM dual pA delta")
    eq(P["s4D"]["B_only pB"][2], -1.76, msg="s4D PM B_only pB delta")
    eq(P["s4D"]["A_only pA"][2], -1.259, msg="s4D PM A_only pA delta")

    if "s3newA" not in P or "s3newB" not in P or "s3newC" not in P or "s3newD" not in P:
        errors.append("Fig S3 paired-delta keys missing")
    else:
        expected_main_delta = {
            "EGFR/HER2": (0.1697, 0.06, 0.2803, True, 0.4297, 0.26),
            "AChE/BChE": (0.1614, 0.037, 0.269, True, 0.6058, 0.4444),
            "PIK3CA/PIK3CB": (0.1511, -0.0215, 0.3105, False, 0.5, 0.3489),
            "PIK3CA/mTOR": (0.0902, -0.1222, 0.2626, False, 0.6921, 0.6019),
        }
        for rec, pair in zip(P["s3newA"], PAIR_ORDER):
            src = D["delta_wp"][(pair, "main_panel")]
            eq(rec["x"], src["delta_matched_minus_wrong"], msg=f"s3A {pair} vs CSV")
            eq(rec["lo"], src["delta_ci_lo"], msg=f"s3A {pair} lo vs CSV")
            eq(rec["hi"], src["delta_ci_hi"], msg=f"s3A {pair} hi vs CSV")
            ex, elo, ehi, eex, em, ew = expected_main_delta[pair]
            eq(rec["x"], ex, msg=f"s3A {pair} delta checksum")
            eq(rec["lo"], elo, msg=f"s3A {pair} lo checksum")
            eq(rec["hi"], ehi, msg=f"s3A {pair} hi checksum")
            eq(rec["matched"], em, msg=f"s3A {pair} matched checksum")
            eq(rec["wrong"], ew, msg=f"s3A {pair} wrong checksum")
            eq(rec["x"], rec["matched"] - rec["wrong"], msg=f"s3A {pair} delta = matched−wrong")
            if rec["excl"] != eex:
                errors.append(f"s3A {pair} ci_excludes_zero {rec['excl']} != {eex}")
            eq(rec["matched"], D["theta6"][pair]["pocket_matched_summary_min"], msg=f"s3A {pair} vs Table 2")
            eq(rec["wrong"], D["pm_by"][(pair, "wrong_pocket_control_vina")]["summary_min"], msg=f"s3A {pair} vs Fig 6 wrong")

        expected_hold_delta = {
            "AChE/BChE": (-0.025, -0.1119, 0.0714, 0.6175, 0.6425),
            "PIK3CA/PIK3CB": (-0.095, -0.2814, 0.1143, 0.425, 0.52),
            "PIK3CA/mTOR": (-0.0225, -0.1165, 0.079, 0.765, 0.7875),
        }
        for rec, pair in zip(P["s3newB"], HOLD_PAIRS):
            src = D["delta_wp"][(pair, "unused_pool_holdout")]
            eq(rec["x"], src["delta_matched_minus_wrong"], msg=f"s3B {pair} vs CSV")
            ex, elo, ehi, em, ew = expected_hold_delta[pair]
            eq(rec["x"], ex, msg=f"s3B {pair} delta checksum")
            eq(rec["lo"], elo, msg=f"s3B {pair} lo checksum")
            eq(rec["hi"], ehi, msg=f"s3B {pair} hi checksum")
            eq(rec["matched"], em, msg=f"s3B {pair} matched checksum")
            eq(rec["wrong"], ew, msg=f"s3B {pair} wrong checksum")
            eq(rec["x"], rec["matched"] - rec["wrong"], msg=f"s3B {pair} delta = matched−wrong")
            if rec["excl"]:
                errors.append(f"s3B {pair}: holdout Δ CI should include 0")
            if rec["x"] >= 0:
                errors.append(f"s3B {pair}: holdout point Δ should be negative")

        expected_desc = {
            "EGFR/HER2": (-0.0524, -0.2, 0.1155, "clogp", 0.4821, 0.4297),
            "AChE/BChE": (-0.1275, -0.3039, 0.0493, "tpsa", 0.7333, 0.6058),
            "PIK3CA/PIK3CB": (-0.1217, -0.3197, 0.0891, "heavy", 0.6217, 0.5),
            "PIK3CA/mTOR": (0.2291, -0.0105, 0.4352, "heavy", 0.463, 0.6921),
        }
        for rec, pair in zip(P["s3newC"], PAIR_ORDER):
            src = D["delta_desc"][pair]
            eq(rec["x"], src["delta_vina_minus_descriptor"], msg=f"s3C {pair} vs CSV")
            ex, elo, ehi, earm, edesc, evina = expected_desc[pair]
            eq(rec["x"], ex, msg=f"s3C {pair} delta checksum")
            eq(rec["lo"], elo, msg=f"s3C {pair} lo checksum")
            eq(rec["hi"], ehi, msg=f"s3C {pair} hi checksum")
            eq(rec["desc"], edesc, msg=f"s3C {pair} descriptor checksum")
            eq(rec["vina"], evina, msg=f"s3C {pair} vina checksum")
            eq(rec["x"], rec["vina"] - rec["desc"], msg=f"s3C {pair} delta = vina−descriptor")
            if rec["arm"] != earm:
                errors.append(f"s3C {pair} arm {rec['arm']} != {earm}")
            if rec["excl"]:
                errors.append(f"s3C {pair}: descriptor Δ CI should include 0")

        eq(P["s3newD"]["sc_db"][0], 0.8527, msg="s3D EGFR scaffold D/B")
        eq(P["s3newD"]["rd_db"][0], 0.8884, msg="s3D EGFR random D/B")
        eq(P["s3newD"]["sc_da"][1], 0.9096, msg="s3D AChE scaffold D/A")
        eq(P["s3newD"]["rd_da"][1], 0.9096, msg="s3D AChE random D/A")
        eq(P["s3newD"]["mean_delta"], 0.0112375, tol=1e-6, msg="s3D mean random−scaffold")
        if abs(P["s3newD"]["mean_delta"]) > 0.03:
            errors.append("s3D mean random−scaffold leakage is not small")

    from PIL import Image

    for name in (
        "Fig6_wrong_pocket_paradox.png",
        "Fig7_confound_anatomy.png",
        "FigS1_protocol_sensitivity.png",
        "FigS2_equal_relation_and_sampling.png",
        "FigS3_paired_delta_bootstrap.png",
    ):
        im = Image.open(OUT / name)
        if im.mode != "RGB":
            errors.append(f"{name} mode {im.mode} != RGB")
        dpi = (im.info.get("dpi") or (300, 300))[0] or 300
        width_in = im.size[0] / dpi
        if abs(width_in - 7.0) > 0.08:
            errors.append(f"{name} width {width_in:.3f} in != 7.00 in (dpi={dpi}, px={im.size})")
        if im.size[1] / dpi > 9.167:
            errors.append(f"{name} height {im.size[1] / dpi:.3f} in exceeds ACS 9.167 in")
        im.close()
