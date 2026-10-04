"""Read-only figures from this round's Stage 2 CSVs and frozen PRIMARY tables.

Does not recompute AUROC, bootstrap, chemistry models, or ranks.
Does not import jcim_stage2_lib or jcim_stage2_compute.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "jcim_stage2"
PRIMARY = ROOT / "results" / "formal_metrics"
OUT = ROOT / "results" / "jcim_stage2_figures"

ORDER = [
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
]
SHORT = {
    "EGFR/HER2": "EGFR/HER2",
    "JAK1/JAK2": "JAK1/JAK2",
    "JAK1/TYK2": "JAK1/TYK2",
    "PIK3CA/mTOR": "PIK3CA/mTOR",
    "AChE/BChE": "AChE/BChE",
    "F2/F10": "F2/F10",
    "PPARG/PPARA": "PPARG/PPARA",
    "PPARA/PPARD": "PPARA/PPARD",
}


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def fnum(value: str) -> float:
    return float(value)


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=160)
    plt.close(fig)


def stored_interval(ax, xs, points, lower, upper, *, color, marker, label):
    """Draw stored endpoints directly, without requiring CI containment of the point."""
    ax.vlines(xs, lower, upper, colors=color, linewidth=1.1)
    ax.plot(xs, lower, linestyle="none", marker="_", color=color, markersize=6)
    ax.plot(xs, upper, linestyle="none", marker="_", color=color, markersize=6)
    ax.plot(xs, points, linestyle="none", marker=marker, color=color, markersize=4.5, label=label)


def write_concordance(compare: list[dict]) -> None:
    states = [
        "CONCORDANT_IMPROVE",
        "CONCORDANT_WORSEN",
        "DISCORDANT",
        "NO_CHANGE",
    ]
    methods = ["M1", "M1b", "M2", "M3"]
    rows = []
    for method in methods:
        subset = [r for r in compare if r["method"] == method]
        counts = {state: sum(r["concordance"] == state for r in subset) for state in states}
        rows.append(
            {
                "method": method,
                "n_pairs": len(subset),
                "supporting_only": "YES" if method == "M1b" else "NO",
                "bootstrap_ci_added": "NO",
                **counts,
            }
        )
    path = OUT / "module_I_concordance_counts.csv"
    OUT.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "method",
                "n_pairs",
                "supporting_only",
                "bootstrap_ci_added",
                *states,
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def fig_primary(directional: list[dict], deltas: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.8))
    x = list(range(len(ORDER)))
    lookup = {(r["pair"], r["contrast"]): r for r in directional if r["method"] == "M0"}
    series = [
        ("dual_vs_A_only", "dual vs A_only (score_B)", "#0072B2", "o", -0.22),
        ("dual_vs_B_only", "dual vs B_only (score_A)", "#D55E00", "s", 0.0),
        ("summary_min", "summary_min", "#009E73", "^", 0.22),
    ]
    for contrast, label, color, marker, offset in series:
        rows = [lookup[pair, contrast] for pair in ORDER]
        stored_interval(
            axes[0], [i + offset for i in x],
            [fnum(r["auroc"]) for r in rows],
            [fnum(r["ci_lo"]) for r in rows], [fnum(r["ci_hi"]) for r in rows],
            color=color, marker=marker, label=label,
        )
    axes[0].axhline(0.5, color="0.6", lw=0.8)
    axes[0].set_ylim(0, 1)
    axes[0].set_title("PRIMARY M0: stored estimates and 95% CIs")
    axes[0].set_ylabel("AUROC")
    axes[0].legend(fontsize=7, loc="lower right")

    delta_names = ["M1-M0", "M1b-M1", "M2-M0", "M3-M0"]
    markers = ["o", "s", "^", "D"]
    colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
    offsets = [-0.27, -0.09, 0.09, 0.27]
    for name, marker, color, offset in zip(delta_names, markers, colors, offsets):
        lookup_delta = {r["pair"]: r for r in deltas if r["delta"] == name}
        rows = [lookup_delta[p] for p in ORDER]
        stored_interval(
            axes[1], [i + offset for i in x],
            [fnum(r["delta_summary_min"]) for r in rows],
            [fnum(r["ci_lo"]) for r in rows], [fnum(r["ci_hi"]) for r in rows],
            color=color, marker=marker, label=name,
        )
    axes[1].axhline(0, color="0.6", lw=0.8)
    axes[1].set_title("PRIMARY paired summary_min deltas and 95% CIs")
    axes[1].set_ylabel("delta summary_min")
    axes[1].legend(fontsize=7, ncol=2)
    for ax in axes:
        ax.set_xticks(list(x))
        ax.set_xticklabels([SHORT[p] for p in ORDER], rotation=35, ha="right")
    save(fig, "fig_primary_m0_and_deltas.png")


def fig_b(rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    width = 0.36
    x = list(range(len(ORDER)))
    for offset, arm, color in ((-width / 2, "D_vs_A", "#1f4e79"), (width / 2, "D_vs_B", "#b85c38")):
        ys, yerr_lo, yerr_hi = [], [], []
        for pair in ORDER:
            row = next(r for r in rows if r["pair"] == pair and r["direction"] == arm)
            y = fnum(row["delta_negative"])
            lo = fnum(row["delta_negative_ci_lo"])
            hi = fnum(row["delta_negative_ci_hi"])
            ys.append(y)
            yerr_lo.append(y - lo)
            yerr_hi.append(hi - y)
        ax.errorbar(
            [i + offset for i in x],
            ys,
            yerr=[yerr_lo, yerr_hi],
            fmt="o",
            color=color,
            label=f"{arm}: neither − single",
            capsize=3,
        )
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[p] for p in ORDER], rotation=35, ha="right")
    ax.set_ylabel("AUROC delta")
    ax.set_title("Module B negative-class contrast (paired bootstrap CI)")
    ax.legend(fontsize=8)
    save(fig, "fig_B_negative_class.png")


def fig_d_and_c(d_rows: list[dict], c_rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.0))
    x = list(range(len(ORDER)))
    series = [
        ("auroc_physchem", "physchem"),
        ("auroc_ecfp4", "ECFP4"),
        ("auroc_ecfp4_m0", "ECFP4+M0"),
        ("auroc_SUPPORTING_CHEMISTRY_MATCHED_M0", "supporting M0"),
    ]
    for col, arm in ((0, "D_vs_A"), (1, "D_vs_B")):
        ax = axes[0, col]
        for key, label in series:
            ys = [
                fnum(next(r for r in d_rows if r["pair"] == pair and r["arm"] == arm)[key])
                for pair in ORDER
            ]
            ax.plot(x, ys, marker="o", label=label)
        ax.set_ylim(0, 1)
        ax.axhline(0.5, color="0.75", lw=0.6)
        ax.set_title(f"D point estimates, {arm}")
        ax.set_ylabel("AUROC")
        ax.legend(fontsize=7, loc="lower left")
        ax.set_xticks(x)
        ax.set_xticklabels([SHORT[p] for p in ORDER], rotation=35, ha="right")
        cx = axes[1, col]
        ys, lo, hi = [], [], []
        for pair in ORDER:
            row = next(r for r in c_rows if r["pair"] == pair and r["arm"] == arm)
            y = fnum(row["delta_pocket"])
            ys.append(y)
            lo.append(y - fnum(row["delta_pocket_ci_lo"]))
            hi.append(fnum(row["delta_pocket_ci_hi"]) - y)
        cx.errorbar(x, ys, yerr=[lo, hi], fmt="o", capsize=3, color="#1f4e79")
        cx.axhline(0, color="0.6", lw=0.8)
        cx.set_title(f"C pocket delta CI, {arm}")
        cx.set_ylabel("correct − wrong")
        cx.set_xticks(x)
        cx.set_xticklabels([SHORT[p] for p in ORDER], rotation=35, ha="right")
    save(fig, "fig_D_chemistry_and_C_pocket.png")


def fig_i(desc: list[dict], compare: list[dict]) -> None:
    fig = plt.figure(figsize=(16.8, 14.2))
    grid = fig.add_gridspec(3, 4, width_ratios=[2.2, 1.0, 2.2, 1.0])
    point_axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in (0, 2)]
    table_axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in (1, 3)]
    composition_ax = fig.add_subplot(grid[2, :2])
    counts_ax = fig.add_subplot(grid[2, 2:])
    desc_by = {(r["pair"], r["method"], r["population_kind"]): r for r in desc}
    compare_by = {(r["pair"], r["method"]): r for r in compare}
    states = ["CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE"]
    state_colors = dict(zip(states, ["#0072B2", "#D55E00", "#CC79A7", "#777777"]))
    state_markers = dict(zip(states, ["o", "s", "D", "^"]))
    labels = dict(zip(states, ["Improve", "Worsen", "Discordant", "No change"]))
    methods = ["M1", "M1b", "M2", "M3"]
    all_dx = [fnum(r["delta_summary_min"]) for r in compare]
    all_dy = [fnum(r["delta_top10_single_target_fraction"]) for r in compare]
    for ax, table_ax, method in zip(point_axes, table_axes, methods):
        table_rows = []
        for pair_index, pair in enumerate(ORDER, start=1):
            row = compare_by[pair, method]
            baseline = desc_by[pair, "M0", f"pairwise_common_vs_{method}"]
            new = desc_by[pair, method, "pairwise_common_vs_M0"]
            state = row["concordance"]
            dx = fnum(row["delta_summary_min"])
            dy = fnum(row["delta_top10_single_target_fraction"])
            ax.scatter(dx, dy, color=state_colors[state], marker=state_markers[state], s=38, zorder=3)
            offset = (-6, 8) if pair_index % 2 == 0 else (6, 8)
            ax.annotate(
                str(pair_index),
                (dx, dy), xytext=offset, textcoords="offset points", fontsize=7,
                ha="right" if offset[0] < 0 else "left", va="center",
            )
            table_rows.append([str(pair_index), f"{baseline['n_dual_top']}→{new['n_dual_top']}/{new['k']}"])
        ax.axhline(0, color="0.65", lw=0.8)
        ax.axvline(0, color="0.65", lw=0.8)
        ax.set_xlim(min(all_dx) - 0.14, max(all_dx) + 0.14)
        ax.set_ylim(min(all_dy) - 0.16, max(all_dy) + 0.18)
        title = f"{method} vs M0" + (" (supporting only)" if method == "M1b" else "")
        ax.set_title(title + "\nPairwise common members; point labels: pair index", fontsize=10)
        source = "PRIMARY point sum" if method == "M1b" else "stored PRIMARY point"
        ax.set_xlabel(f"delta summary_min ({source})", fontsize=9)
        ax.set_ylabel("delta Top-k single-target fraction", fontsize=9)
        table_ax.set_axis_off()
        table = table_ax.table(cellText=table_rows, colLabels=["Pair", "dual M0→new/k"], cellLoc="center", colWidths=[0.2, 0.8], loc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.8)
    m0 = [
        r
        for r in desc
        if r["method"] == "M0" and r["population_kind"] == "method_specific" and r["table"] == "descriptive"
    ]
    by = {r["pair"]: r for r in m0}
    x = list(range(len(ORDER)))
    dual = [int(float(by[p]["n_dual_top"])) for p in ORDER]
    a_only = [int(float(by[p]["n_A_only_top"])) for p in ORDER]
    b_only = [int(float(by[p]["n_B_only_top"])) for p in ORDER]
    composition_ax.bar(x, dual, label="dual")
    composition_ax.bar(x, a_only, bottom=dual, label="A_only")
    bottom2 = [u + v for u, v in zip(dual, a_only)]
    composition_ax.bar(x, b_only, bottom=bottom2, label="B_only")
    for i, pair in enumerate(ORDER):
        composition_ax.text(i, int(float(by[pair]["k"])) + 0.15, f"k={by[pair]['k']}", ha="center", fontsize=7)
    composition_ax.set_title("M0 method-specific description: Top-k class counts")
    composition_ax.set_ylabel("ligands in Top-k")
    composition_ax.legend(fontsize=8)
    composition_ax.set_xticks(x)
    composition_ax.set_xticklabels([f"{i + 1}: {SHORT[p]}" for i, p in enumerate(ORDER)], rotation=35, ha="right")

    method_labels = ["M1", "M1b*", "M2", "M3"]
    method_keys = ["M1", "M1b", "M2", "M3"]
    bottoms = [0] * len(method_labels)
    for state in states:
        heights = []
        for key in method_keys:
            heights.append(sum(r["method"] == key and r["concordance"] == state for r in compare))
        counts_ax.bar(method_labels, heights, bottom=bottoms, color=state_colors[state], label=labels[state])
        bottoms = [b + h for b, h in zip(bottoms, heights)]
    counts_ax.set_ylim(0, 8.5)
    counts_ax.set_ylabel("pairs")
    counts_ax.set_title("Four-state counts (* M1b supporting only)")
    counts_ax.legend(fontsize=8)
    save(fig, "fig_I_topk_and_concordance.png")


def fig_e_f(seeds: list[dict], summary: list[dict], alt: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 10.0))
    x = list(range(len(ORDER)))
    seed_keys = [17, 29, 42, 71, 101]
    seed_colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
    seed_markers = ["o", "^", "s", "D", "v"]
    seed_offsets = [-0.18, -0.09, 0.0, 0.09, 0.18]
    seed_by = {(r["pair"], int(r["seed"])): r for r in seeds}
    panels = [
        (axes[0, 0], "auc_B", "E: auc_B = D_vs_A (score_B)"),
        (axes[0, 1], "auc_A", "E: auc_A = D_vs_B (score_A)"),
        (axes[1, 0], "summary_min", "E: summary_min = min(auc_B, auc_A)"),
    ]
    for ax, metric, title in panels:
        for seed, color, marker, offset in zip(seed_keys, seed_colors, seed_markers, seed_offsets):
            ys = [fnum(seed_by[pair, seed][metric]) for pair in ORDER]
            ax.scatter([i + offset for i in x], ys, marker=marker, s=24, color=color, label=f"seed {seed}", zorder=2)
        meds = [
            fnum(next(r for r in summary if r["pair"] == pair and r["metric"] == metric)["median"])
            for pair in ORDER
        ]
        ax.scatter(x, meds, marker="_", s=110, color="black", linewidths=1.6, zorder=3, label="stored median")
        iqr_labels = [
            f"IQR={fnum(next(r for r in summary if r['pair'] == pair and r['metric'] == metric)['iqr']):.3f}"
            for pair in ORDER
        ]
        ax.set_ylim(0, 1)
        ax.axhline(0.5, color="0.7", lw=0.8)
        ax.set_title(title + "\nFive stored seed points; IQR is not a CI", fontsize=11)
        ax.set_ylabel("AUROC")
        ax.legend(fontsize=7, ncol=2, loc="upper right")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{SHORT[p]}\n{lab}" for p, lab in zip(ORDER, iqr_labels)], rotation=35, ha="right", fontsize=8)

    alt_rows = sorted(alt, key=lambda r: (ORDER.index(r["pair"]), r["replaced_side"]))
    xf = list(range(len(alt_rows)))
    ax = axes[1, 1]
    for field, color, marker, offset, label in [
        ("delta_AUC_affected", "#0072B2", "o", -0.18, "affected-side delta and CI"),
        ("delta_summary_min", "#CC79A7", "D", 0.18, "summary_min delta and CI"),
    ]:
        stored_interval(
            ax, [i + offset for i in xf], [fnum(r[field]) for r in alt_rows],
            [fnum(r[field + "_ci_lo"]) for r in alt_rows],
            [fnum(r[field + "_ci_hi"]) for r in alt_rows],
            color=color, marker=marker, label=label,
        )
    ax.plot(xf, [fnum(r["delta_AUC_unreplaced"]) for r in alt_rows], "s", color="#E69F00", markersize=4.5, label="unreplaced-side delta")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_title("F: alternate receptor − primary\nSame members; stored paired 95% CIs", fontsize=11)
    ax.set_ylabel("AUROC / summary_min delta")
    ax.legend(fontsize=7)
    ax.set_xticks(xf)
    ax.set_xticklabels(
        [f"{SHORT[r['pair']]}\n{r['alt_id']}" for r in alt_rows], rotation=35, ha="right", fontsize=8
    )
    save(fig, "fig_E_iqr_and_F_alt.png")


def main() -> None:
    directional = read_csv(PRIMARY / "PRIMARY_DIRECTIONAL_METRICS.csv")
    deltas = read_csv(PRIMARY / "PRIMARY_METHOD_DELTA_METRICS.csv")
    b_rows = read_csv(RESULTS / "module_B_negative_class.csv")
    c_rows = read_csv(RESULTS / "module_C_wrong_pocket.csv")
    d_rows = read_csv(RESULTS / "module_D_chemistry_units.csv")
    desc = read_csv(RESULTS / "module_I_descriptive.csv")
    compare = read_csv(RESULTS / "module_I_method_compare.csv")
    seeds = read_csv(RESULTS / "module_E_five_seed.csv")
    summary = read_csv(RESULTS / "module_E_five_seed_summary.csv")
    alt = read_csv(RESULTS / "module_F_alt_receptor.csv")
    write_concordance(compare)
    fig_primary(directional, deltas)
    fig_b(b_rows)
    fig_d_and_c(d_rows, c_rows)
    fig_i(desc, compare)
    fig_e_f(seeds, summary, alt)
    print("FIGURES_WRITTEN", OUT)


if __name__ == "__main__":
    main()
