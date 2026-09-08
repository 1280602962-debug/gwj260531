#!/usr/bin/env python3
"""JCAMD publication figures from archived experimental files only.

Every plotted number is read from a frozen CSV/JSON under data/ and asserted
against the locked protocol / funnel counts. Nothing is invented, interpolated,
or back-filled from manuscript prose.

Style: Springer large-journal column (84 mm / 174 mm), Arial-metric sans-serif
8 pt, Okabe–Ito colours, no grid, TrueType-embedded PDF, 600 dpi TIFF.
Captions live in figures/generated/CAPTIONS.md — never overlaid on the art.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "figures"))

from jmm_style import (  # noqa: E402
    AUDIT,
    FIG_OUT,
    FONT_SIZE_PT,
    GATE_FILL,
    KNOWN,
    MUTED,
    NEUTRAL,
    NLRP3_COLOR,
    PARETO,
    RANDOM_DECOY,
    THRESHOLD,
    TRUE_DECOY,
    URAT1_COLOR,
    URAT1_LIGHT,
    WARN,
    apply_panel_tags,
    apply_style,
    clean_axes,
    clean_hbar_axes,
    figsize_double,
    figsize_single,
    save_figure,
    set_axis_labels,
    ylim_headroom,
)

DATA = PROJECT_ROOT / "data"
DOCS = PROJECT_ROOT / "docs"
MANIFEST_PATH = FIG_OUT / "figure_manifest.json"
LOCK_PATH = FIG_OUT / "DATA_LOCK.json"

PROTOCOLS = ["P0", "P1", "P2", "P3", "P4", "P5"]
READOUT_SHORT = {
    "P0": "CNNscore",
    "P1": "Vina",
    "P2": "CNNaffinity",
    "P3": "minAff.",
    "P4": "RTM/Vina",
    "P5": "RTM/gnina",
}
AUDIT_DISPLAY = {
    "VECABRUTINIB": "Vecabrutinib",
    "ZELENIRSTAT": "Zelenirstat",
    "DEUCRICTIBANT": "Deucrictibant",
    "PRALICIGUAT": "Praliciguat",
    "GSK-3008348 FREE BASE": "GSK-3008348",
    "GSK-3008348": "GSK-3008348",
    "MLN-0415": "MLN-0415",
    "BI 653048": "BI 653048",
}
KNOWN_DISPLAY = {
    "LESINURAD": "lesinurad",
    "VERINURAD": "verinurad",
    "COLCHICINE": "colchicine",
}
# Textbook URAT1 drugs queried against the p>=6 active set (library + raw ChEMBL IDs).
URAT1_DRUG_QUERY = [
    ("lesinurad", ("CHEMBL2105720", "CHEMBL3301572"), "library"),
    ("benzbromarone", ("CHEMBL388590", "CHEMBL892"), "library"),
    ("dotinurad", ("CHEMBL4594446", "CHEMBL4594374"), "library"),
    ("probenecid", ("CHEMBL897",), "library"),
    ("verinurad", ("CHEMBL3707347", "CHEMBL3989871"), "library"),
    ("puliginurad", ("CHEMBL5314438",), "library"),
    ("SHR-4640", ("CHEMBL3746329",), "library"),
    ("isobavachin", ("CHEMBL5618175",), "raw_chembl"),
]


def _assert(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _close(a: float, b: float, tol: float, msg: str) -> None:
    if abs(float(a) - float(b)) > tol:
        raise AssertionError(f"{msg}: {a} vs {b} (tol={tol})")


def _display_name(raw: str) -> str:
    key = str(raw).strip()
    return AUDIT_DISPLAY.get(key, AUDIT_DISPLAY.get(key.upper(), KNOWN_DISPLAY.get(key.upper(), key.title())))


def load_and_lock() -> dict:
    """Load archived tables and assert frozen counts used in Results."""
    ef = pd.read_csv(DATA / "si" / "protocol_enrichment_ci" / "protocol_ef_ci.csv")
    redock = pd.read_csv(DATA / "redock_smoke" / "redock_results_lesinurad_9DKB.csv")
    pareto = pd.read_csv(DATA / "repurposing" / "p2" / "pareto_merged_scores.csv")
    nominated = pd.read_csv(DATA / "repurposing" / "p2" / "nominated_shortlist_diverse.csv")
    pose = pd.read_csv(DATA / "si" / "pose_qc" / "pose_qc_table.csv")
    pose_dual = pd.read_csv(DATA / "si" / "pose_qc" / "pose_qc_dual.csv")
    gates = pd.read_csv(DATA / "si" / "nomination_sensitivity" / "gate_counts.csv")
    screen = json.loads((DATA / "repurposing" / "screening" / "nlrp3_screening_summary_clinical_all.json").read_text())
    complete = json.loads((DATA / "si" / "complete_case_drop" / "summary.json").read_text())
    leak = json.loads((DATA / "si" / "decoy_leakage_audit" / "summary.json").read_text())
    nom_sum = json.loads((DATA / "repurposing" / "p2" / "candidate_nomination_summary.json").read_text())
    train = json.loads((DOCS / "MODEL_TRAINING_SUMMARY.json").read_text())
    data_sum = json.loads((DATA / "processed" / "data_summary.json").read_text())
    actives = pd.read_csv(DATA / "benchmarks" / "urat1_true_decoy" / "actives.csv")
    rand_tc = pd.read_csv(DATA / "si" / "decoy_leakage_audit" / "random_decoy_max_tc_to_active.csv")
    weak_tc = pd.read_csv(DATA / "si" / "decoy_leakage_audit" / "weak_active_max_tc_to_active.csv")
    paired_true = pd.read_csv(DATA / "si" / "protocol_paired_bootstrap" / "paired_bootstrap_vs_P2_true_decoy.csv")
    ml = pd.read_csv(DATA / "repurposing" / "screening" / "nlrp3_ml_scores_clinical_all.csv", low_memory=False)
    lib = pd.read_csv(DATA / "repurposing" / "repurposing_manifest.csv", low_memory=False)

    n_complete = len(pareto)
    n_valid = int(pareto["dock_score"].notna().sum())
    n_gate = int(((pareto["s_u_percentile"] >= 90) & (pareto["s_n_dock_percentile"] >= 90)).sum())
    n_pareto = int(pareto["pareto_front"].astype(bool).sum())
    missing = pareto.loc[pareto["dock_score"].isna(), "name"].tolist()

    _assert(n_complete == 1580, f"complete-case rows {n_complete}")
    _assert(n_valid == 1579, f"valid dual scores {n_valid}")
    _assert(n_gate == 51, f"dual-structure gate {n_gate}")
    _assert(n_pareto == 4, f"Pareto front {n_pareto}")
    _assert(missing == ["TAUROSELCHOLIC ACID"], f"empty-pose row {missing}")
    _assert(int(screen["n_scored"]) == 8319, "library n")
    _assert(int(screen["n_pred_active_ge_threshold"]) == 1588, "ML pool n")
    _assert(int(complete["n_dual_complete_case"]) == 1580, "complete-case json")
    _assert(int(nom_sum["n_dual_gate"]) == 51, "nomination gate")
    _assert(int(nom_sum["n_preferred_candidate"]) == 7, "preferred n")
    _assert(len(nominated) == 7, "shortlist rows")
    _assert(int(leak["random_decoy"]["n_scaffold_overlap_with_actives"]) == 0, "RandomDecoy scaffold leak")
    _assert(int(leak["random_decoy"]["n_tc_gt_0.5"]) == 0, "RandomDecoy TC>0.5")
    _assert(len(actives) == 469, "n actives")
    _assert(len(ml) == 8319, "ML score rows")

    p2_true = ef[(ef["protocol"] == "P2") & (ef["benchmark"] == "true")].iloc[0]
    p2_rand = ef[(ef["protocol"] == "P2") & (ef["benchmark"] == "random")].iloc[0]
    _close(p2_true["EF1pct"], 2.5867, 1e-4, "P2 True EF@1%")
    _close(p2_rand["EF1pct"], 0.2154, 1e-4, "P2 Random EF@1%")
    _close(p2_true["AUC"], 0.58, 1e-6, "P2 True AUC")
    _assert(p2_true["hits_at_1pct"] == "12/51", "P2 True hits")
    _assert(p2_rand["hits_at_1pct"] == "1/51", "P2 Random hits")

    p2_redock = redock[(redock["protocol_id"] == "P2") & (redock["exhaustiveness"] == 32)].iloc[0]
    _close(p2_redock["top1_rmsd_A"], 4.163, 5e-4, "P2 Top-1 RMSD")
    _close(p2_redock["best_ensemble_rmsd_A"], 0.994, 5e-4, "P2 best RMSD")

    gsk = pose[(pose["name"] == "GSK-3008348") & (pose["target"] == "URAT1")].iloc[0]
    les = pose[(pose["name"] == "LESINURAD") & (pose["target"] == "URAT1")].iloc[0]
    ver = pose[(pose["name"] == "VERINURAD") & (pose["target"] == "URAT1")].iloc[0]
    _close(gsk["acid_arg477"], 3.188923642861329, 1e-6, "GSK Arg477")
    _close(les["acid_arg477"], 14.200913104445076, 1e-6, "lesinurad Arg477")
    _close(ver["acid_arg477"], 2.859423718164201, 1e-6, "verinurad Arg477")
    _assert(bool(pose_dual["both_in_pocket"].all()), "preferred 7 both in pocket")

    les_row = pareto[pareto["name"].str.upper() == "LESINURAD"].iloc[0]
    ver_row = pareto[pareto["name"].str.upper() == "VERINURAD"].iloc[0]
    col_row = pareto[pareto["name"].str.upper() == "COLCHICINE"].iloc[0]
    _close(les_row["s_u_percentile"], 45.85443037974684, 1e-6, "lesinurad S_U")
    _close(ver_row["s_u_percentile"], 3.417721518987338, 1e-6, "verinurad S_U")
    _close(col_row["s_u_percentile"], 63.22784810126583, 1e-6, "colchicine S_U")

    tau90 = gates[gates["tau"] == 90].iloc[0]
    _assert(int(tau90["n_dual_gate"]) == 51, "gate_counts tau90 gate")
    _assert(int(tau90["n_preferred_mw_oral"]) == 7, "gate_counts tau90 preferred")

    sc = actives["scaffold"].value_counts()
    _assert(int(sc.iloc[0]) == 127, "top scaffold count")
    _assert(int(actives["scaffold"].nunique()) == 118, "n scaffolds")

    act_ids = set(actives["molecule_chembl_id"].astype(str))
    named_in_actives = {}
    for drug, cids, _src in URAT1_DRUG_QUERY:
        named_in_actives[drug] = any(c in act_ids for c in cids)
    _assert(named_in_actives["lesinurad"] is False, "lesinurad must be absent from actives")
    _assert(named_in_actives["benzbromarone"] is False, "benzbromarone absent")
    _assert(named_in_actives["dotinurad"] is False, "dotinurad absent")
    _assert(named_in_actives["verinurad"] is True, "verinurad present")
    _assert(named_in_actives["puliginurad"] is True, "puliginurad present")
    _assert(named_in_actives["SHR-4640"] is True, "SHR-4640 present")
    _assert(named_in_actives["isobavachin"] is True, "isobavachin present")

    lock = {
        "n_library": 8319,
        "n_ml_pool": 1588,
        "n_complete_case": 1580,
        "n_valid_dual": 1579,
        "n_dual_gate": 51,
        "n_preferred": 7,
        "n_pareto": 4,
        "empty_pose": "TAUROSELCHOLIC ACID",
        "pareto_names": pareto.loc[pareto["pareto_front"].astype(bool), "name"].tolist(),
        "preferred_names": nominated["name"].tolist(),
        "p2_true_ef1": float(p2_true["EF1pct"]),
        "p2_random_ef1": float(p2_rand["EF1pct"]),
        "p2_true_auc": float(p2_true["AUC"]),
        "p2_top1_rmsd": float(p2_redock["top1_rmsd_A"]),
        "gsk_acid_arg477": float(gsk["acid_arg477"]),
        "lesinurad_acid_arg477": float(les["acid_arg477"]),
        "n_actives": 469,
        "n_scaffolds": 118,
        "top_scaffold_n": 127,
        "named_in_actives": named_in_actives,
        "nlrp3_auroc": float(train["nlrp3"]["cv_metrics"]["auroc"]),
        "nlrp3_auprc": float(train["nlrp3"]["cv_metrics"]["auprc"]),
        "nlrp3_ef10": float(train["nlrp3"]["cv_metrics"]["ef_10pct"]),
        "font": "Liberation Sans (Arial-metric; Arial/Helvetica unavailable on this builder)",
    }
    _close(lock["nlrp3_auroc"], 0.8934017784193573, 1e-9, "NLRP3 AUROC")

    return {
        "ef": ef,
        "redock": redock,
        "pareto": pareto,
        "nominated": nominated,
        "pose": pose,
        "pose_dual": pose_dual,
        "gates": gates,
        "screen": screen,
        "leak": leak,
        "train": train,
        "data_sum": data_sum,
        "actives": actives,
        "rand_tc": rand_tc,
        "weak_tc": weak_tc,
        "paired_true": paired_true,
        "ml": ml,
        "lib": lib,
        "lock": lock,
        "named_in_actives": named_in_actives,
        "scaffold_counts": sc,
        "nom_sum": nom_sum,
        "complete": complete,
    }


def _ef_row(ef: pd.DataFrame, protocol: str, bench: str) -> pd.Series:
    hit = ef[(ef["protocol"] == protocol) & (ef["benchmark"] == bench)]
    _assert(len(hit) == 1, f"expected 1 row for {protocol}/{bench}")
    return hit.iloc[0]


def plot_fig01(d: dict) -> dict:
    """Protocol EF@1% (True vs Random) and lesinurad self-dock RMSD at exh=32."""
    ef, redock = d["ef"], d["redock"]
    fig, axes = plt.subplots(1, 2, figsize=figsize_double(88))
    x = np.arange(len(PROTOCOLS))
    width = 0.36

    ax = axes[0]
    true_y, rand_y = [], []
    true_err, rand_err = [], []
    for p in PROTOCOLS:
        t, r = _ef_row(ef, p, "true"), _ef_row(ef, p, "random")
        true_y.append(float(t["EF1pct"]))
        rand_y.append(float(r["EF1pct"]))
        true_err.append((float(t["EF1pct"]) - float(t["EF1pct_ci95_low"]), float(t["EF1pct_ci95_high"]) - float(t["EF1pct"])))
        rand_err.append((float(r["EF1pct"]) - float(r["EF1pct_ci95_low"]), float(r["EF1pct_ci95_high"]) - float(r["EF1pct"])))
    true_y, rand_y = np.array(true_y), np.array(rand_y)
    true_err, rand_err = np.array(true_err).T, np.array(rand_err).T
    # Highlight locked P2 behind the bars (not on top of data).
    p2_x = x[PROTOCOLS.index("P2")]
    ax.axvspan(p2_x - 0.52, p2_x + 0.52, color=GATE_FILL, alpha=0.10, zorder=0, lw=0)
    ax.bar(x - width / 2, true_y, width, color=TRUE_DECOY, edgecolor=NEUTRAL, linewidth=0.4, label="TrueDecoy", zorder=2)
    ax.bar(x + width / 2, rand_y, width, color=RANDOM_DECOY, edgecolor=NEUTRAL, linewidth=0.4, label="RandomDecoy", zorder=2)
    ax.errorbar(x - width / 2, true_y, yerr=true_err, fmt="none", ecolor=NEUTRAL, elinewidth=0.7, capsize=2.2, capthick=0.6, zorder=3)
    ax.errorbar(x + width / 2, rand_y, yerr=rand_err, fmt="none", ecolor=NEUTRAL, elinewidth=0.7, capsize=2.2, capthick=0.6, zorder=3)
    ax.axhline(1.0, color=NEUTRAL, linestyle=(0, (4, 3)), linewidth=0.7, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(["P0", "P1", "P2*", "P3", "P4", "P5"])
    clean_axes(ax)
    set_axis_labels(ax, "Protocol  (* locked Π*;  P0 = negative control)", "EF@1%")
    ax.set_ylim(0, 5.0)
    ax.set_xlim(-0.7, 5.7)
    ax.legend(loc="upper left", frameon=False, borderaxespad=0.2)
    ax.text(0.98, 0.04, "chance = 1", transform=ax.transAxes, ha="right", va="bottom", fontsize=FONT_SIZE_PT, color=NEUTRAL)

    ax = axes[1]
    sub = redock[redock["exhaustiveness"] == 32].set_index("protocol_id").loc[PROTOCOLS]
    top1 = sub["top1_rmsd_A"].to_numpy(dtype=float)
    best = sub["best_ensemble_rmsd_A"].to_numpy(dtype=float)
    ax.bar(x - width / 2, top1, width, color=TRUE_DECOY, edgecolor=NEUTRAL, linewidth=0.4, label="Top-1 RMSD")
    ax.bar(x + width / 2, best, width, color=URAT1_LIGHT, edgecolor=NEUTRAL, linewidth=0.4, label="Best-in-ensemble")
    ax.axhline(2.0, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(PROTOCOLS)
    clean_axes(ax)
    set_axis_labels(ax, "Protocol (exhaustiveness = 32)", "lesinurad self-dock RMSD (Å)")
    ax.set_ylim(0, 7.2)
    ax.set_xlim(-0.7, 5.7)
    ax.legend(loc="upper right", frameon=False, borderaxespad=0.2)
    ax.text(0.02, 0.96, "2 Å gate", transform=ax.transAxes, ha="left", va="top", fontsize=FONT_SIZE_PT, color=THRESHOLD)
    # Label only Top-1 poses that fail the 2 Å gate; passing values sit on the ensemble
    # bars and would cover them.
    for xi, val in zip(x - width / 2, top1):
        if val > 2.0:
            ax.text(xi, val + 0.12, f"{val:.2f}", ha="center", va="bottom", fontsize=FONT_SIZE_PT - 1, color=NEUTRAL, clip_on=False)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.22, wspace=0.38)
    apply_panel_tags(fig, axes, ("a", "b"))
    paths = save_figure(fig, "fig01_protocol_enrichment_selfdock", "main")
    return {"id": "fig01", "panel": "main", "description": "P0–P5 EF@1% True vs Random and lesinurad self-dock RMSD", **paths}


def plot_fig02(d: dict) -> dict:
    """Clinical-library funnel and dual-percentile scatter."""
    pareto, nominated, lock = d["pareto"], d["nominated"], d["lock"]
    fig, axes = plt.subplots(1, 2, figsize=figsize_double(96), gridspec_kw={"width_ratios": [0.92, 1.18]})

    ax = axes[0]
    stages = [
        "Clinical library",
        r"$q_N\geq0.5$ (NLRP3 ML)",
        "Dual-dock complete case",
        r"Gate $S_U$ and $S_{N,\mathrm{dock}}\geq90$",
        "Chemistry-filtered audit set",
    ]
    counts = [lock["n_library"], lock["n_ml_pool"], lock["n_complete_case"], lock["n_dual_gate"], lock["n_preferred"]]
    colors = [MUTED, NLRP3_COLOR, URAT1_COLOR, THRESHOLD, AUDIT]
    y = np.arange(len(stages))
    bars = ax.barh(y, counts, color=colors, edgecolor=NEUTRAL, linewidth=0.4, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(5, 2.0e4)
    for bar, n in zip(bars, counts):
        ax.text(n * 1.12, bar.get_y() + bar.get_height() / 2, f"{n:,}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
    clean_hbar_axes(ax)
    set_axis_labels(ax, "Number of compounds (log scale)", "")
    ax.tick_params(axis="y", length=0)

    ax = axes[1]
    valid = pareto[pareto["dock_score"].notna()].copy()
    _assert(len(valid) == 1579, "scatter n")
    pref_names = set(nominated["name"].astype(str).str.upper())
    pareto_mask = valid["pareto_front"].astype(bool)
    pref_mask = valid["name"].astype(str).str.upper().isin(pref_names)
    known_mask = valid["name"].astype(str).str.upper().isin(KNOWN_DISPLAY)
    gate = (valid["s_u_percentile"] >= 90) & (valid["s_n_dock_percentile"] >= 90)
    bg = valid[~pareto_mask & ~pref_mask & ~known_mask]
    ax.scatter(bg["s_u_percentile"], bg["s_n_dock_percentile"], s=7, alpha=0.28, color=MUTED, edgecolors="none", rasterized=True, zorder=1, label="Complete case (n = 1,579)")
    ax.add_patch(Rectangle((90, 90), 10, 10, fill=False, edgecolor=THRESHOLD, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2))
    ax.scatter(
        valid.loc[gate, "s_u_percentile"],
        valid.loc[gate, "s_n_dock_percentile"],
        s=14, alpha=0.85, color=THRESHOLD, edgecolors="none", zorder=3, label=r"Gate $\tau=90$ (n = 51)",
    )
    ax.scatter(
        valid.loc[pareto_mask, "s_u_percentile"],
        valid.loc[pareto_mask, "s_n_dock_percentile"],
        s=36, marker="D", color=PARETO, edgecolors=NEUTRAL, linewidth=0.4, zorder=4, label="Pareto front (n = 4, audit)",
    )
    ax.scatter(
        valid.loc[pref_mask, "s_u_percentile"],
        valid.loc[pref_mask, "s_n_dock_percentile"],
        s=42, marker="o", color=AUDIT, edgecolors=NEUTRAL, linewidth=0.5, zorder=5, label="Chemistry-filtered list (n = 7)",
    )
    ax.scatter(
        valid.loc[known_mask, "s_u_percentile"],
        valid.loc[known_mask, "s_n_dock_percentile"],
        s=44, marker="s", color=KNOWN, edgecolors="white", linewidth=0.5, zorder=6, label="Known gout drugs in pool",
    )
    ax.set_xlim(0, 103)
    ax.set_ylim(0, 103)
    clean_axes(ax)
    set_axis_labels(ax, r"$S_U$ percentile", r"$S_{N,\mathrm{dock}}$ percentile")
    # Legend below the scatter, clear of the 90/90 box and of the known-drug squares.
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.72, 0.045),
        ncol=2,
        frameon=False,
        fontsize=FONT_SIZE_PT,
        handletextpad=0.4,
        columnspacing=1.2,
        borderaxespad=0.0,
    )

    fig.subplots_adjust(left=0.24, right=0.98, top=0.90, bottom=0.22, wspace=0.42)
    apply_panel_tags(fig, axes, ("a", "b"))
    paths = save_figure(fig, "fig02_funnel_dual_percentiles", "main")
    return {"id": "fig02", "panel": "main", "description": "Clinical-library funnel and P2 dual-percentile scatter", **paths}


def plot_fig03(d: dict) -> dict:
    """Pose QC: carboxylate–Arg477, URAT1 COM, NLRP3 COM. No candidate language."""
    pose, pose_dual = d["pose"], d["pose_dual"]
    fig, axes = plt.subplots(1, 3, figsize=figsize_double(88))

    # (a) acid–Arg477 only where the archived column is populated.
    u = pose[pose["target"] == "URAT1"].copy()
    acid = u[u["acid_arg477"].notna()].copy()
    acid["label"] = acid["name"].map(_display_name)
    # Order: lesinurad (control, long), then verinurad, then GSK (audit-set acid).
    order = ["lesinurad", "verinurad", "GSK-3008348"]
    acid["label"] = pd.Categorical(acid["label"], categories=order, ordered=True)
    acid = acid.sort_values("label")
    colors = [KNOWN if lab != "GSK-3008348" else AUDIT for lab in acid["label"]]
    ax = axes[0]
    y = np.arange(len(acid))
    bars = ax.barh(y, acid["acid_arg477"].to_numpy(float), color=colors, edgecolor=NEUTRAL, linewidth=0.4, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(list(acid["label"]))
    ax.invert_yaxis()
    ax.set_xlim(0, 17.5)
    for bar, val in zip(bars, acid["acid_arg477"]):
        ax.text(val + 0.25, bar.get_y() + bar.get_height() / 2, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
    clean_hbar_axes(ax)
    set_axis_labels(ax, "Carboxylate O–Arg477 (Å)", "")
    ax.tick_params(axis="y", length=0)

    # (b) URAT1 COM-to-lesinurad for 7 + 2 controls.
    ax = axes[1]
    pref_order = ["VECABRUTINIB", "ZELENIRSTAT", "DEUCRICTIBANT", "PRALICIGUAT", "GSK-3008348", "MLN-0415", "BI 653048"]
    ctrl_order = ["LESINURAD", "VERINURAD"]
    rows = []
    for name in pref_order + ctrl_order:
        hit = u[u["name"].str.upper() == name]
        _assert(len(hit) == 1, f"URAT1 pose row {name}")
        rows.append((_display_name(name), float(hit.iloc[0]["com_to_ref"]), name in ctrl_order))
    y = np.arange(len(rows))
    vals = [r[1] for r in rows]
    cols = [KNOWN if r[2] else AUDIT for r in rows]
    bars = ax.barh(y, vals, color=cols, edgecolor=NEUTRAL, linewidth=0.4, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 6.6)
    ax.set_xticks([0, 2, 4, 6])
    ax.axvline(6.0, color=MUTED, linestyle=(0, (4, 3)), linewidth=0.7)
    for bar, val in zip(bars, vals):
        ax.text(val + 0.08, bar.get_y() + bar.get_height() / 2, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
    clean_hbar_axes(ax)
    set_axis_labels(ax, "URAT1 COM to lesinurad (Å)", "")
    ax.tick_params(axis="y", length=0)

    # (c) NLRP3 COM-to-RM5 for the 7 only (controls were not archived on 7ALV).
    ax = axes[2]
    n = pose[pose["target"] == "NLRP3"].copy()
    rows = []
    for name in pref_order:
        hit = n[n["name"].str.upper() == name]
        _assert(len(hit) == 1, f"NLRP3 pose row {name}")
        rows.append((_display_name(name), float(hit.iloc[0]["com_to_ref"])))
    y = np.arange(len(rows))
    vals = [r[1] for r in rows]
    bars = ax.barh(y, vals, color=NLRP3_COLOR, edgecolor=NEUTRAL, linewidth=0.4, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 4.0)
    for bar, val in zip(bars, vals):
        ax.text(val + 0.06, bar.get_y() + bar.get_height() / 2, f"{val:.2f}", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
    clean_hbar_axes(ax)
    set_axis_labels(ax, "NLRP3 COM to NP3-146 (Å)", "")
    ax.tick_params(axis="y", length=0)

    _assert(bool(pose_dual["both_in_pocket"].all()), "both_in_pocket")
    fig.subplots_adjust(left=0.16, right=0.97, top=0.90, bottom=0.16, wspace=0.70)
    apply_panel_tags(fig, axes, ("a", "b", "c"))
    paths = save_figure(fig, "fig03_pose_qc", "main")
    return {"id": "fig03", "panel": "main", "description": "Production-pose QC: Arg477 salt-bridge geometry and COM displacements", **paths}


def plot_fig04(d: dict) -> dict:
    """Active-set composition and decoy leakage — why the ranker is not an activity retriever."""
    fig, axes = plt.subplots(2, 2, figsize=figsize_double(148))

    # (a) named URAT1 drugs present vs absent in the 469-active set
    ax = axes[0, 0]
    drugs = [q[0] for q in URAT1_DRUG_QUERY]
    present = [d["named_in_actives"][name] for name in drugs]
    y = np.arange(len(drugs))
    ax.set_yticks(y)
    ax.set_yticklabels(drugs)
    ax.invert_yaxis()
    for yi, p in zip(y, present):
        ax.plot(0.18, yi, marker="o" if p else "x", markersize=7 if p else 8, color=THRESHOLD if p else WARN, markeredgecolor=NEUTRAL, markeredgewidth=0.4, linestyle="none")
        ax.text(0.32, yi, "in 469-set" if p else "absent", va="center", ha="left", fontsize=FONT_SIZE_PT, color=NEUTRAL, clip_on=False)
    ax.set_xlim(0, 1.35)
    ax.set_xticks([])
    clean_hbar_axes(ax)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0)
    ax.spines["bottom"].set_visible(False)
    ax.set_xlabel("")

    # (b) top Murcko scaffolds (labels are rank + count; SMILES in SI table)
    ax = axes[0, 1]
    sc = d["scaffold_counts"].head(5)
    y = np.arange(len(sc))
    bars = ax.barh(y, sc.values, color=TRUE_DECOY, edgecolor=NEUTRAL, linewidth=0.4, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([f"Scaffold {i+1}" for i in range(len(sc))])
    ax.invert_yaxis()
    ax.set_xlim(0, 160)
    for bar, n in zip(bars, sc.values):
        ax.text(n + 2, bar.get_y() + bar.get_height() / 2, f"{int(n)} / 469", va="center", ha="left", fontsize=FONT_SIZE_PT, clip_on=False)
    clean_hbar_axes(ax)
    set_axis_labels(ax, "Actives sharing Murcko scaffold", "")
    ax.tick_params(axis="y", length=0)

    # (c) RandomDecoy max TC to any active
    ax = axes[1, 0]
    tc = d["rand_tc"]["max_tc_active"].to_numpy(float)
    _assert(len(tc) == 4690, "random decoy TC rows")
    ax.hist(tc, bins=40, color=RANDOM_DECOY, edgecolor="white", linewidth=0.2, range=(0, 0.6))
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.8)
    clean_axes(ax)
    set_axis_labels(ax, "Max Tanimoto to any TrueDecoy active", "RandomDecoy count")
    ax.set_xlim(0, 0.6)
    ylim_headroom(ax, 0.12)
    ax.text(0.86, 0.94, "TC = 0.5\nn > 0.5: 0", transform=ax.transAxes, ha="left", va="top", fontsize=FONT_SIZE_PT, color=THRESHOLD)

    # (d) weak-active max TC (designed hard negatives; overlap is expected)
    ax = axes[1, 1]
    wtc = d["weak_tc"]["max_tc_active"].to_numpy(float)
    _assert(len(wtc) == 80, "weak-active TC rows")
    ax.hist(wtc, bins=20, color=WARN, edgecolor="white", linewidth=0.2, range=(0, 0.9))
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.8)
    clean_axes(ax)
    set_axis_labels(ax, "Max Tanimoto to any TrueDecoy active", "Weak-active count")
    ax.set_xlim(0, 0.9)
    ylim_headroom(ax, 0.12)
    n_gt = int((wtc > 0.5).sum())
    ax.text(0.70, 0.94, f"TC = 0.5\nn > 0.5: {n_gt}", transform=ax.transAxes, ha="left", va="top", fontsize=FONT_SIZE_PT, color=THRESHOLD)

    fig.subplots_adjust(left=0.16, right=0.96, top=0.92, bottom=0.10, hspace=0.48, wspace=0.42)
    apply_panel_tags(fig, axes, ("a", "b", "c", "d"))
    paths = save_figure(fig, "fig04_active_set_decoy_leakage", "main")
    return {"id": "fig04", "panel": "main", "description": "TrueDecoy active-set membership, scaffold bias, and decoy Tanimoto audit", **paths}


def plot_s01_gate(d: dict) -> dict:
    gates = d["gates"]
    fig, ax = plt.subplots(figsize=figsize_single(78))
    tau = gates["tau"].to_numpy(float)
    ax.plot(tau, gates["n_dual_gate"], marker="o", color=THRESHOLD, linewidth=1.0, markersize=4.5, label="Dual-structure gate")
    ax.plot(tau, gates["n_preferred_mw_oral"], marker="s", color=AUDIT, linewidth=1.0, markersize=4.5, label="Chemistry-filtered list")
    ax.axvline(90, color=MUTED, linestyle=(0, (4, 3)), linewidth=0.7)
    clean_axes(ax)
    set_axis_labels(ax, r"Percentile gate $\tau$", "Number of compounds")
    ax.set_xlim(68, 97)
    ax.set_ylim(0, float(gates["n_dual_gate"].max()) * 1.18)
    ax.legend(loc="upper right", frameon=False)
    ax.text(90, float(gates["n_dual_gate"].max()) * 1.05, "τ = 90", fontsize=FONT_SIZE_PT, color=MUTED, ha="center", va="bottom", clip_on=False)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.18)
    paths = save_figure(fig, "figS01_gate_sensitivity", "si")
    return {"id": "figS01", "panel": "si", "description": "Nomination-gate sensitivity (does not replace τ = 90)", **paths}


def plot_s02_auc(d: dict) -> dict:
    ef = d["ef"]
    fig, ax = plt.subplots(figsize=figsize_single(78))
    x = np.arange(len(PROTOCOLS))
    width = 0.36
    true_y, rand_y, true_err, rand_err = [], [], [], []
    for p in PROTOCOLS:
        t, r = _ef_row(ef, p, "true"), _ef_row(ef, p, "random")
        true_y.append(float(t["AUC"]))
        rand_y.append(float(r["AUC"]))
        true_err.append((float(t["AUC"]) - float(t["AUC_ci95_low"]), float(t["AUC_ci95_high"]) - float(t["AUC"])))
        rand_err.append((float(r["AUC"]) - float(r["AUC_ci95_low"]), float(r["AUC_ci95_high"]) - float(r["AUC"])))
    true_y, rand_y = np.array(true_y), np.array(rand_y)
    true_err, rand_err = np.array(true_err).T, np.array(rand_err).T
    p2_x = x[PROTOCOLS.index("P2")]
    ax.axvspan(p2_x - 0.52, p2_x + 0.52, color=GATE_FILL, alpha=0.10, zorder=0, lw=0)
    ax.bar(x - width / 2, true_y, width, color=TRUE_DECOY, edgecolor=NEUTRAL, linewidth=0.4, label="TrueDecoy")
    ax.bar(x + width / 2, rand_y, width, color=RANDOM_DECOY, edgecolor=NEUTRAL, linewidth=0.4, label="RandomDecoy")
    ax.errorbar(x - width / 2, true_y, yerr=true_err, fmt="none", ecolor=NEUTRAL, elinewidth=0.7, capsize=2.2, capthick=0.6)
    ax.errorbar(x + width / 2, rand_y, yerr=rand_err, fmt="none", ecolor=NEUTRAL, elinewidth=0.7, capsize=2.2, capthick=0.6)
    ax.axhline(0.5, color=NEUTRAL, linestyle=(0, (4, 3)), linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(PROTOCOLS)
    clean_axes(ax)
    set_axis_labels(ax, "Protocol", "ROC-AUC")
    ax.set_ylim(0.45, 0.72)
    ax.legend(loc="upper right", frameon=False)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.18)
    paths = save_figure(fig, "figS02_protocol_auc", "si")
    return {"id": "figS02", "panel": "si", "description": "P0–P5 ROC-AUC on TrueDecoy and RandomDecoy", **paths}


def plot_s03_nlrp3_cv(d: dict) -> dict:
    """NLRP3 scaffold-CV fold metrics from MODEL_TRAINING_SUMMARY.json (no OOF curve file)."""
    folds = d["train"]["nlrp3"]["fold_metrics"]
    overall = d["train"]["nlrp3"]["cv_metrics"]
    fig, ax = plt.subplots(figsize=figsize_single(78))
    x = np.arange(len(folds))
    width = 0.36
    auroc = [float(f["auroc"]) for f in folds]
    auprc = [float(f["auprc"]) for f in folds]
    ax.bar(x - width / 2, auroc, width, color=NLRP3_COLOR, edgecolor=NEUTRAL, linewidth=0.4, label="AUROC")
    ax.bar(x + width / 2, auprc, width, color="#E69F00", edgecolor=NEUTRAL, linewidth=0.4, label="AUPRC")
    ax.axhline(overall["auroc"], color=NLRP3_COLOR, linestyle=(0, (4, 3)), linewidth=0.7)
    ax.axhline(overall["auprc"], color="#E69F00", linestyle=(0, (4, 3)), linewidth=0.7)
    ax.axhline(0.5, color=NEUTRAL, linestyle=(0, (4, 3)), linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Fold {f['fold']+1}" for f in folds])
    clean_axes(ax)
    set_axis_labels(ax, "Scaffold-grouped CV fold", "OOF metric")
    ax.set_ylim(0.50, 1.05)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=2, frameon=False, borderaxespad=0.0)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.82, bottom=0.18)
    paths = save_figure(fig, "figS03_nlrp3_fold_metrics", "si")
    return {"id": "figS03", "panel": "si", "description": "NLRP3 assay-conditioned classifier scaffold-CV fold AUROC/AUPRC", **paths}


def plot_s04_known_drugs(d: dict) -> dict:
    """Known gout-related controls that entered the 1,580-row table."""
    pareto = d["pareto"]
    fig, ax = plt.subplots(figsize=figsize_single(78))
    names = ["LESINURAD", "VERINURAD", "COLCHICINE"]
    su, sn = [], []
    for name in names:
        row = pareto[pareto["name"].str.upper() == name].iloc[0]
        su.append(float(row["s_u_percentile"]))
        sn.append(float(row["s_n_dock_percentile"]))
    x = np.arange(len(names))
    width = 0.36
    ax.bar(x - width / 2, su, width, color=URAT1_COLOR, edgecolor=NEUTRAL, linewidth=0.4, label=r"$S_U$")
    ax.bar(x + width / 2, sn, width, color=NLRP3_COLOR, edgecolor=NEUTRAL, linewidth=0.4, label=r"$S_{N,\mathrm{dock}}$")
    ax.axhline(90, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([KNOWN_DISPLAY[n] for n in names])
    clean_axes(ax)
    set_axis_labels(ax, "Known gout-related drug in the complete case", "Percentile (n = 1,580)")
    ax.set_ylim(0, 105)
    ax.legend(loc="upper left", frameon=False)
    for xi, a, b in zip(x, su, sn):
        ax.text(xi - width / 2, a + 1.5, f"{a:.1f}", ha="center", va="bottom", fontsize=FONT_SIZE_PT - 1, clip_on=False)
        ax.text(xi + width / 2, b + 1.5, f"{b:.1f}", ha="center", va="bottom", fontsize=FONT_SIZE_PT - 1, clip_on=False)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.18)
    paths = save_figure(fig, "figS04_known_gout_percentiles", "si")
    return {"id": "figS04", "panel": "si", "description": "lesinurad, verinurad, colchicine P2 percentiles in the 1,580-row table", **paths}


def plot_s05_nlrp3_hist(d: dict) -> dict:
    ml = d["ml"]
    fig, ax = plt.subplots(figsize=figsize_single(76))
    p = ml["p_active_nlrp3"].to_numpy(float)
    n_ge = int((p >= 0.5).sum())
    _assert(n_ge == 1588, "hist n>=0.5")
    ax.hist(p, bins=40, color=NLRP3_COLOR, edgecolor="white", linewidth=0.2, range=(0, 1))
    ax.axvline(0.5, color=THRESHOLD, linestyle=(0, (4, 3)), linewidth=0.8)
    clean_axes(ax)
    set_axis_labels(ax, r"NLRP3 $q_N$ (predicted P(active))", "Clinical-library count")
    ax.set_xlim(0, 1)
    ylim_headroom(ax, 0.10)
    ax.text(0.54, 0.94, f"threshold 0.5\nn = {n_ge:,}", transform=ax.transAxes, ha="left", va="top", fontsize=FONT_SIZE_PT, color=THRESHOLD)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.18)
    paths = save_figure(fig, "figS05_nlrp3_qN_histogram", "si")
    return {"id": "figS05", "panel": "si", "description": "NLRP3 q_N distribution on the 8,319-compound clinical library", **paths}


def write_si_tables(d: dict) -> None:
    out = FIG_OUT / "tables"
    out.mkdir(parents=True, exist_ok=True)
    sc = d["scaffold_counts"].reset_index()
    sc.columns = ["murcko_smiles", "n_actives"]
    sc.insert(0, "rank", np.arange(1, len(sc) + 1))
    sc["fraction_of_469"] = sc["n_actives"] / 469
    sc.head(10).to_csv(out / "table_s_top_scaffolds.csv", index=False)

    rows = []
    act_ids = set(d["actives"]["molecule_chembl_id"].astype(str))
    for drug, cids, src in URAT1_DRUG_QUERY:
        rows.append({
            "name": drug,
            "chembl_ids_queried": ";".join(cids),
            "id_source": src,
            "in_p6_active_set": d["named_in_actives"][drug],
            "matching_chembl_id": next((c for c in cids if c in act_ids), ""),
        })
    pd.DataFrame(rows).to_csv(out / "table_s_named_urat1_in_actives.csv", index=False)

    # Audit-set scores exactly as plotted / tabulated.
    nom = d["nominated"][["name", "s_u_percentile", "s_n_dock_percentile", "s_n_ml_percentile", "mw", "max_phase"]].copy()
    nom["display_name"] = nom["name"].map(_display_name)
    nom.to_csv(out / "table_s_audit_set_percentiles.csv", index=False)


def write_captions(d: dict) -> None:
    lock = d["lock"]
    gsk = lock["gsk_acid_arg477"]
    les = lock["lesinurad_acid_arg477"]
    captions = f"""# Figure captions (JCAMD)

Lettering is 8 pt sans-serif (Liberation Sans; Arial/Helvetica metric-compatible substitute — Arial is not licensed on the builder). All numbers below are copied from archived files under `data/`; they are the same values asserted by `scripts/plot_jcamd_publication_figures.py`. Panel letters sit outside the data area. Captions are **not** drawn on the artwork.

## Figure 1. Locked URAT1 docking readout is a weak activity retriever and is not pose-accurate.

**a** Enrichment factor at the top 1% (EF@1%) for protocols P0–P5 on the pre-registered TrueDecoy and RandomDecoy benchmarks (9DKB). Bars are point estimates; whiskers are molecule-resampled bootstrap 95% percentile intervals (1,000 draws) from `data/si/protocol_enrichment_ci/protocol_ef_ci.csv`. The dashed line is chance (EF = 1). P2 (gnina CNNaffinity) is the production readout locked by the pre-registered rule (green band): TrueDecoy EF@1% = {lock['p2_true_ef1']:.2f} (12/51, hypergeometric *p* = 0.0016), RandomDecoy EF@1% = {lock['p2_random_ef1']:.2f} (1/51, 95% CI 0.00–1.04). P0 is the pre-registered negative-control readout (CNNscore), not a salvage protocol. **b** lesinurad self-docking on 9DKB at exhaustiveness = 32, `num_modes` = 9 (`data/redock_smoke/redock_results_lesinurad_9DKB.csv`). Top-1 heavy-atom RMSD versus the best RMSD in the nine-pose ensemble. The dashed line is the 2 Å Top-1 gate. P2 Top-1 RMSD = {lock['p2_top1_rmsd']:.2f} Å (fails the gate); every protocol’s ensemble minimum is ≤ 1.00 Å.

## Figure 2. Transfer of the locked P2 readout onto the clinical library is an audit, not a dual-node nomination.

**a** Compound counts at each frozen funnel step: ChEMBL clinical library *n* = {lock['n_library']:,}; NLRP3 ML shrink q_N ≥ 0.5, *n* = {lock['n_ml_pool']:,}; P2 dual-dock complete-case table *n* = {lock['n_complete_case']:,} (1,579 rows with valid dual scores; one empty-pose row, tauroselcholic acid, retained with percentile 0); dual-structure gate S_U ≥ 90 and S_N,dock ≥ 90, *n* = {lock['n_dual_gate']}; chemistry-filtered list (Veber + Ro5 HBD/HBA/logP, MW 200–550 Da, macrolide demotion), *n* = {lock['n_preferred']}. **b** S_U versus S_N,dock for the 1,579 valid dual-score rows. The dashed box is the τ = 90 gate. Red diamonds: raw docking Pareto front (*n* = 4; Idremcinal, Alemcinal, Cethromycin, Zamzetoclax — macrolide/erythromycin audit, not a follow-up list). Blue circles: chemistry-filtered list (*n* = 7). Black squares: gout-related drugs that are in the complete-case table (lesinurad, verinurad, colchicine). Names are omitted from the scatter so they do not cover points; identities are in Fig. 3, Fig. S4, and `figures/generated/tables/`.

## Figure 3. Production P2 poses sit in the crystal cavity but lesinurad has lost the Arg477 salt bridge.

**a** Shortest carboxylate-oxygen to Arg477 distance in the production first pose, reported only for molecules with an archived `acid_arg477` value (`data/si/pose_qc/pose_qc_table.csv`): lesinurad {les:.2f} Å, verinurad {d['pose'][(d['pose']['name']=='VERINURAD')&(d['pose']['target']=='URAT1')].iloc[0]['acid_arg477']:.2f} Å, GSK-3008348 {gsk:.2f} Å. The other six chemistry-filtered molecules are not carboxylic acids and have no value in that column. **b** URAT1 centre-of-mass displacement relative to co-crystal lesinurad for the seven chemistry-filtered molecules plus the two URAT1 controls. The dashed line is the 6 Å COM in-pocket cutoff used in pose QC. **c** NLRP3 COM displacement relative to co-crystal NP3-146/RM5 for the same seven molecules. lesinurad and verinurad were not archived on 7ALV. All seven have `both_in_pocket = True` and zero 2.2 Å clashes (`pose_qc_dual.csv`). These geometries do not constitute binding-mode or affinity evidence.

## Figure 4. The TrueDecoy active set is analog-biased and excludes textbook URAT1 drugs; RandomDecoy is not a near-neighbour leak.

**a** Membership of named URAT1-related compounds in the 469-compound p≥6 TrueDecoy active set, by ChEMBL identifier (library IDs plus the identifiers listed in `docs/DATA_FACT_CHECK.md`; isobavachin from `data/raw/URAT1_CHEMBL_cf12.csv`). lesinurad, benzbromarone, dotinurad and probenecid are absent; verinurad, puliginurad, SHR-4640 and isobavachin are present. **b** The five most frequent Murcko scaffolds among the 469 actives (118 unique scaffolds). Scaffold 1 accounts for 127/469 (27.1%); SMILES are in `tables/table_s_top_scaffolds.csv`. **c** Maximum Tanimoto similarity of each of 4,690 RandomDecoy molecules to any TrueDecoy active. None exceed 0.5; scaffold overlap with actives is 0 (`data/si/decoy_leakage_audit/`). **d** The same metric for 80 experimental weak actives (designed hard negatives): 14 have TC > 0.5 and none have TC > 0.85.

## Figure S1. Gate-threshold sensitivity (does not replace τ = 90).

Number of compounds passing the dual-structure gate and the subsequent chemistry filter as a function of percentile threshold τ (`data/si/nomination_sensitivity/gate_counts.csv`). The production analysis is locked at τ = 90 (51 → 7). Wider gates are monotonic sensitivity, not a second production shortlist.

## Figure S2. ROC-AUC for P0–P5 on TrueDecoy and RandomDecoy.

Point estimates and bootstrap 95% percentile intervals from the same file as Fig. 1a. The dashed line is AUC = 0.5. P2 TrueDecoy AUC = {lock['p2_true_auc']:.3f}. Absolute discrimination remains modest; AUC is not the locking metric.

## Figure S3. NLRP3 assay-conditioned classifier, scaffold-grouped CV.

Per-fold out-of-fold AUROC and AUPRC from `docs/MODEL_TRAINING_SUMMARY.json` (five folds). Dashed lines: pooled AUROC = {lock['nlrp3_auroc']:.3f}, AUPRC = {lock['nlrp3_auprc']:.3f}. EF@10% = {lock['nlrp3_ef10']:.2f} at prevalence ≈ 60%. No ROC/PR curve is drawn because the OOF prediction table is not in this archive.

## Figure S4. Known gout-related drugs that reached the dual-dock table.

P2 percentiles inside the archived 1,580-row denominator. None meet S_U ≥ 90 and S_N,dock ≥ 90. benzbromarone, dotinurad, allopurinol, febuxostat and probenecid have q_N = 0 and never entered the 1,588 docking pool (not plotted).

## Figure S5. NLRP3 q_N on the 8,319-compound clinical library.

Histogram of predicted P(active). The dashed line is the pre-registered shrink threshold 0.5 (*n* = 1,588).
"""
    (FIG_OUT / "CAPTIONS.md").write_text(captions)

    captions_cn = (
        "# 图注（中文，与英文 CAPTIONS.md 同一套数字）\n\n"
        "图面不写图注。字体 8 pt 无衬线。所有数字来自 `data/` 归档，由绘图脚本断言。\n\n"
        "**图 1.** a. P0–P5 在 TrueDecoy / RandomDecoy 上的 EF@1%（bootstrap 95% 区间）；"
        "虚线为随机（EF = 1）；P2 为按预注册规则锁定的生产读出，P0 为负对照。"
        "b. lesinurad@9DKB 自对接（exhaustiveness = 32）Top-1 与集合最低 RMSD；虚线为 2 Å。"
        "P2 Top-1 = 4.16 Å，不能当构象金标准。\n\n"
        "**图 2.** a. 冻结漏斗：8319 → 1588 → 1580 → 51 → 7。"
        "b. 1,579 个有效双靶分数的 S_U–S_N,dock 散点；虚线框 τ = 90；"
        "红钻为裸 Pareto 4（大环审计）；蓝圈为化学过滤名单 7；黑方为库内痛风相关对照。"
        "名单不作“双节点候选”标注。\n\n"
        "**图 3.** a. 仅对归档了 `acid_arg477` 的分子：lesinurad 14.20 Å，verinurad 2.86 Å，"
        "GSK-3008348 3.19 Å。b. URAT1 质心相对 lesinurad。c. NLRP3 质心相对 NP3-146。"
        "姿态在腔内 ≠ 结合模式或亲和力。\n\n"
        "**图 4.** a. 命名 URAT1 相关分子是否落入 469 个 p≥6 阳性："
        "lesinurad / 苯溴马隆 / dotinurad / 丙磺舒缺席；"
        "verinurad / puliginurad / SHR-4640 / isobavachin 在内。"
        "b. 阳性集 Top-5 Murcko 骨架（骨架 1 为 127/469）。"
        "c. RandomDecoy 相对阳性的最大 Tanimoto，无一 >0.5。"
        "d. 80 个实验弱活分子的同一指标。\n\n"
        "**图 S1–S5.** 门控敏感性、协议 AUC、NLRP3 折间 AUROC/AUPRC、库内对照百分位、"
        "临床库 q_N 直方图。不画 MD 轨迹，不画不存在的 OOF ROC 曲线。\n"
    )
    (FIG_OUT / "CAPTIONS_CN.md").write_text(captions_cn)


def main() -> None:
    apply_style()
    d = load_and_lock()
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    LOCK_PATH.write_text(json.dumps(d["lock"], indent=2))
    write_si_tables(d)

    entries = [
        plot_fig01(d),
        plot_fig02(d),
        plot_fig03(d),
        plot_fig04(d),
        plot_s01_gate(d),
        plot_s02_auc(d),
        plot_s03_nlrp3_cv(d),
        plot_s04_known_drugs(d),
        plot_s05_nlrp3_hist(d),
    ]
    write_captions(d)
    manifest = {
        "journal": "Journal of Computer-Aided Molecular Design",
        "style": {
            "font": "Liberation Sans 8 pt (Arial-metric substitute)",
            "columns_mm": {"single": 84, "double": 174},
            "grid": "none",
            "palette": "Okabe–Ito",
            "pdf_fonttype": 42,
            "tiff_dpi": "not written (vector PDF is the Springer line-art format)",
            "png_dpi": 300,
        },
        "data_lock": str(LOCK_PATH),
        "captions": str(FIG_OUT / "CAPTIONS.md"),
        "note": "Seven chemistry-filtered names are an audit of a failed ranker, not dual-node candidates. No MD trajectories. No fabricated OOF ROC.",
        "figures": entries,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {len(entries)} figure sets + captions to {FIG_OUT}")
    print(json.dumps({k: d["lock"][k] for k in ("n_library", "n_ml_pool", "n_complete_case", "n_preferred", "p2_true_ef1", "p2_random_ef1")}, indent=2))


if __name__ == "__main__":
    main()
