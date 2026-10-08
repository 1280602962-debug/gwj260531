"""Read-only manuscript figures from accepted PRIMARY and Stage 2 tables.

Does not import scientific compute modules. Does not recompute AUROC, bootstrap,
chemistry models, or ranks. Does not rewrite module_I_concordance_counts.csv.
"""

from __future__ import annotations

import csv
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "jcim_stage2"
PRIMARY = ROOT / "results" / "formal_metrics"
OUT = ROOT / "results" / "jcim_stage2_figures"
CONCORDANCE = OUT / "module_I_concordance_counts.csv"

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
ARMS = ["D_vs_A", "D_vs_B"]
METHODS = ["M0", "M1", "M1b", "M2", "M3"]
DELTAS = ["M1-M0", "M1b-M1", "M2-M0", "M3-M0"]
SEEDS = [17, 29, 42, 71, 101]
DESCRIPTORS = [
    ("MolWt", "MolWt (Da)"),
    ("MolLogP", "MolLogP (calculated)"),
    ("TPSA", "TPSA (A$^2$)"),
    ("NumHDonors", "H-bond donors"),
    ("NumHAcceptors", "H-bond acceptors"),
    ("NumRotatableBonds", "Rotatable bonds"),
    ("RingCount", "Ring count"),
    ("FormalCharge", "Formal charge"),
]
STATES = ["CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE"]
STATE_COLOR = {
    "CONCORDANT_IMPROVE": "#0072B2",
    "CONCORDANT_WORSEN": "#D55E00",
    "DISCORDANT": "#CC79A7",
    "NO_CHANGE": "#666666",
}
STATE_MARKER = {
    "CONCORDANT_IMPROVE": "o",
    "CONCORDANT_WORSEN": "s",
    "DISCORDANT": "D",
    "NO_CHANGE": "^",
}
SEED_COLOR = {
    17: "#0072B2",
    29: "#E69F00",
    42: "#009E73",
    71: "#CC79A7",
    101: "#56B4E9",
}
SEED_MARKER = {17: "o", 29: "^", 42: "s", 71: "D", 101: "v"}
C_BLUE, C_ORANGE, C_GREEN, C_PINK = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
PAGE_W = 7.0


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def num(value: str) -> float:
    return float(value)


def save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.svg")
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    plt.close(fig)


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def new_fig(height: float, rows: int = 1, cols: int = 1):
    return plt.subplots(rows, cols, figsize=(PAGE_W, height), layout="constrained")


def stored_interval(ax, xs, points, lower, upper, *, color, marker, label, hollow=False):
    ax.vlines(xs, lower, upper, colors=color, linewidth=0.8, zorder=2)
    ax.plot(xs, lower, linestyle="none", marker="_", color=color, markersize=5)
    ax.plot(xs, upper, linestyle="none", marker="_", color=color, markersize=5)
    ax.plot(
        xs,
        points,
        linestyle="none",
        marker=marker,
        color=color,
        markersize=4.2,
        markerfacecolor="none" if hollow else color,
        markeredgewidth=0.8,
        label=label,
        zorder=3,
    )


def arm_units():
    units = []
    for pair in ORDER:
        for arm in ARMS:
            units.append((pair, arm))
    return units


def unit_labels():
    labels = []
    for index, pair in enumerate(ORDER, start=1):
        labels.append(f"{index} {pair} D_vs_A")
        labels.append(f"{index} {pair} D_vs_B")
    return labels


def scheme() -> None:
    fig, ax = new_fig(4.6)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Scheme 1. Directional tasks and score workflows")

    def box(x, y, w, h, text, color):
        patch = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor=color, edgecolor="0.2", linewidth=0.6,
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5)

    box(0.3, 7.3, 4.4, 2.2, "Labels\ndual / A_only / B_only / neither\nprimary threshold 6.0", "#E6F1FA")
    box(5.3, 7.3, 4.4, 2.2, "Two directional tasks\nD_vs_A uses score_B\nD_vs_B uses score_A", "#FFF4E5")
    box(0.3, 4.2, 9.4, 2.4, "summary_min = min(AUROC of D_vs_A, AUROC of D_vs_B)\nnot AUROC of the per-ligand minimum score", "#F4F4F4")
    box(0.3, 0.4, 2.2, 3.2, "M0\nVina seed 42\nmode 1\nfixed pose", "#E5F5EF")
    box(2.7, 0.4, 2.2, 3.2, "M1 / M2\nrescore the\nsame M0 pose\nCNNscore / RTMScore", "#E5F5EF")
    box(5.1, 0.4, 2.2, 3.2, "M1b\nmax CNNscore\non saved poses\nrerank only", "#FFF6D8")
    box(7.5, 0.4, 2.2, 3.2, "M3\nindependent\nGNINA docking\nmode 1 CNNscore", "#FDE8E4")
    ax.annotate("", xy=(5.2, 8.4), xytext=(4.8, 8.4), arrowprops={"arrowstyle": "->", "lw": 0.7})
    save(fig, "Scheme01_Study_Design")


def fig01(directional, deltas) -> None:
    fig, axes = new_fig(4.8, 1, 2)
    x = np.arange(len(ORDER))
    lookup = {(r["pair"], r["contrast"]): r for r in directional if r["method"] == "M0"}
    series = [
        ("dual_vs_A_only", "D_vs_A (score_B)", C_BLUE, "o", -0.22),
        ("dual_vs_B_only", "D_vs_B (score_A)", C_ORANGE, "s", 0.0),
        ("summary_min", "summary_min", C_GREEN, "^", 0.22),
    ]
    for contrast, label, color, marker, offset in series:
        rows = [lookup[pair, contrast] for pair in ORDER]
        stored_interval(
            axes[0], x + offset, [num(r["auroc"]) for r in rows],
            [num(r["ci_lo"]) for r in rows], [num(r["ci_hi"]) for r in rows],
            color=color, marker=marker, label=label,
        )
    axes[0].axhline(0.5, color="0.6", lw=0.6)
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("AUROC")
    axes[0].set_title("a  M0 directional AUROC")
    axes[0].legend(loc="lower right", frameon=False)
    colors = [C_BLUE, "#E69F00", C_GREEN, C_PINK]
    markers = ["o", "s", "^", "D"]
    offsets = [-0.27, -0.09, 0.09, 0.27]
    for name, color, marker, offset in zip(DELTAS, colors, markers, offsets):
        rows = [next(r for r in deltas if r["pair"] == pair and r["delta"] == name) for pair in ORDER]
        stored_interval(
            axes[1], x + offset, [num(r["delta_summary_min"]) for r in rows],
            [num(r["ci_lo"]) for r in rows], [num(r["ci_hi"]) for r in rows],
            color=color, marker=marker, label=name,
        )
    axes[1].axhline(0, color="0.45", lw=0.6)
    axes[1].set_ylabel("delta summary_min")
    axes[1].set_title("b  Paired summary_min deltas")
    axes[1].legend(loc="lower left", frameon=False, ncol=2)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels([str(i) for i in range(1, 9)])
    save(fig, "Fig01_Primary_Directional_and_Deltas")


def fig02(rows) -> None:
    fig, axes = new_fig(6.2, 1, 2)
    units = arm_units()
    y = np.arange(len(units))[::-1]
    by = {(r["pair"], r["direction"]): r for r in rows}
    single = [num(by[u]["auc_single"]) for u in units]
    neither = [num(by[u]["auc_neither"]) for u in units]
    for yi, a, b in zip(y, single, neither):
        axes[0].plot([a, b], [yi, yi], color="0.75", lw=0.6, zorder=1)
    axes[0].scatter(single, y, marker="o", color=C_BLUE, s=16, label="single-target negative", zorder=2)
    axes[0].scatter(neither, y, marker="s", color=C_ORANGE, s=16, label="neither", zorder=2)
    axes[0].axvline(0.5, color="0.6", lw=0.6)
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("AUROC")
    axes[0].set_title("a  Same dual set, two negative classes")
    axes[0].legend(frameon=False, loc="lower right")
    points = [num(by[u]["delta_negative"]) for u in units]
    lo = [num(by[u]["delta_negative_ci_lo"]) for u in units]
    hi = [num(by[u]["delta_negative_ci_hi"]) for u in units]
    axes[1].hlines(y, lo, hi, color=C_BLUE, lw=0.8)
    axes[1].scatter(points, y, color=C_BLUE, s=16, zorder=2)
    axes[1].axvline(0, color="0.45", lw=0.6)
    axes[1].set_xlabel("neither − single")
    axes[1].set_title("b  Paired delta and 95% CI")
    labels = []
    for (pair, arm), row_key in zip(units, units):
        rec = by[row_key]
        index = ORDER.index(pair) + 1
        labels.append(f"{index} {pair} {arm}\nn={rec['n_dual']}/{rec['n_single']}/{rec['n_neither']}")
    for ax in axes:
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6)
    save(fig, "Fig02_Negative_Class_Contrast")


def fig03(rows) -> None:
    fig, axes = new_fig(6.4, 1, 2)
    units = arm_units()
    y = np.arange(len(units))[::-1]
    by = {(r["pair"], r["arm"]): r for r in rows}
    series = [
        ("auroc_physchem", "physchem", C_BLUE, "o", -0.18),
        ("auroc_ecfp4", "ECFP4", C_ORANGE, "s", -0.06),
        ("auroc_ecfp4_m0", "ECFP4+M0", C_GREEN, "^", 0.06),
        ("auroc_SUPPORTING_CHEMISTRY_MATCHED_M0", "matched M0", C_PINK, "D", 0.18),
    ]
    for key, label, color, marker, offset in series:
        axes[0].scatter(
            [num(by[u][key]) for u in units], y + offset,
            marker=marker, color=color, s=14, label=label, zorder=2,
        )
    axes[0].axvline(0.5, color="0.6", lw=0.6)
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("OOF AUROC")
    axes[0].set_title("a  Layer-3 models, no interval")
    axes[0].legend(frameon=False, loc="lower right", fontsize=6.5)
    deltas = [num(by[u]["delta_ecfp4_plus_m0"]) for u in units]
    axes[1].scatter(deltas, y, color=C_GREEN, s=16)
    axes[1].axvline(0, color="0.45", lw=0.6)
    axes[1].set_xlabel("ECFP4+M0 − ECFP4")
    axes[1].set_title("b  Point increment only")
    labels = [f"{ORDER.index(p)+1} {p} {a}" for p, a in units]
    for ax in axes:
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.5)
    save(fig, "Fig03_Chemistry_and_Docking_Increment")


def fig04(rows) -> None:
    fig, axes = new_fig(6.2, 1, 2)
    units = arm_units()
    y = np.arange(len(units))[::-1]
    by = {(r["pair"], r["arm"]): r for r in rows}
    correct = [num(by[u]["auc_correct"]) for u in units]
    wrong = [num(by[u]["auc_wrong"]) for u in units]
    for yi, a, b in zip(y, correct, wrong):
        axes[0].plot([a, b], [yi, yi], color="0.75", lw=0.6)
    axes[0].scatter(correct, y, marker="o", color=C_BLUE, s=16, label="correct pocket")
    axes[0].scatter(wrong, y, marker="s", color=C_ORANGE, s=16, label="opposite pocket")
    axes[0].axvline(0.5, color="0.6", lw=0.6)
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("AUROC")
    axes[0].set_title("a  Same members, two pockets")
    axes[0].legend(frameon=False, loc="lower right")
    points = [num(by[u]["delta_pocket"]) for u in units]
    lo = [num(by[u]["delta_pocket_ci_lo"]) for u in units]
    hi = [num(by[u]["delta_pocket_ci_hi"]) for u in units]
    axes[1].hlines(y, lo, hi, color=C_BLUE, lw=0.8)
    axes[1].scatter(points, y, color=C_BLUE, s=16)
    axes[1].axvline(0, color="0.45", lw=0.6)
    axes[1].set_xlabel("correct − wrong")
    axes[1].set_title("b  Paired delta and 95% CI")
    labels = [f"{ORDER.index(p)+1} {p} {a}" for p, a in units]
    for ax in axes:
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.5)
    save(fig, "Fig04_Pocket_Contrast")


def fig05(seeds, summary, alt) -> None:
    fig, axes = new_fig(4.8, 1, 2)
    x = np.arange(len(ORDER))
    seed_by = {(r["pair"], int(r["seed"])): r for r in seeds}
    for offset, seed in zip(np.linspace(-0.24, 0.24, 5), SEEDS):
        ys = [num(seed_by[pair, seed]["summary_min"]) for pair in ORDER]
        axes[0].scatter(x + offset, ys, s=16, color=SEED_COLOR[seed], marker=SEED_MARKER[seed], label=f"seed {seed}", zorder=2)
    med = [num(next(r for r in summary if r["pair"] == p and r["metric"] == "summary_min")["median"]) for p in ORDER]
    lo = [num(next(r for r in summary if r["pair"] == p and r["metric"] == "summary_min")["min"]) for p in ORDER]
    hi = [num(next(r for r in summary if r["pair"] == p and r["metric"] == "summary_min")["max"]) for p in ORDER]
    axes[0].vlines(x, lo, hi, color="0.35", lw=0.7, zorder=1)
    axes[0].scatter(x, med, marker="_", s=80, color="black", label="stored median", zorder=3)
    iqr = [f"IQR {num(next(r for r in summary if r['pair']==p and r['metric']=='summary_min')['iqr']):.3f}" for p in ORDER]
    axes[0].set_ylim(0, 1)
    axes[0].axhline(0.5, color="0.7", lw=0.5)
    axes[0].set_ylabel("summary_min AUROC")
    axes[0].set_title("a  Five stored seeds; line = stored min–max, not a CI")
    axes[0].legend(frameon=False, ncol=3, fontsize=6)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"{i}\n{lab}" for i, lab in enumerate(iqr, start=1)], fontsize=6)
    alt_rows = sorted(alt, key=lambda r: (ORDER.index(r["pair"]), r["replaced_side"]))
    xf = np.arange(len(alt_rows))
    for field, color, marker, offset, label in [
        ("delta_AUC_affected", C_BLUE, "o", -0.12, "affected side"),
        ("delta_summary_min", C_PINK, "D", 0.12, "summary_min"),
    ]:
        stored_interval(
            axes[1], xf + offset, [num(r[field]) for r in alt_rows],
            [num(r[field + "_ci_lo"]) for r in alt_rows],
            [num(r[field + "_ci_hi"]) for r in alt_rows],
            color=color, marker=marker, label=label,
        )
    axes[1].axhline(0, color="0.45", lw=0.6)
    axes[1].set_ylabel("alternate − primary")
    axes[1].set_title("b  Stored paired 95% CIs")
    axes[1].legend(frameon=False, fontsize=6.5)
    axes[1].set_xticks(xf)
    axes[1].set_xticklabels([f"{r['pair']}\n{r['alt_id']} {r['replaced_side']}" for r in alt_rows], fontsize=5.5, rotation=30, ha="right")
    save(fig, "Fig05_Seed_and_Receptor_Sensitivity")


def fig06(compare) -> None:
    fig, axes = new_fig(6.6, 2, 2)
    axes = axes.ravel()
    by = {(r["pair"], r["method"]): r for r in compare}
    xs, ys = [], []
    for row in compare:
        xs.append(num(row["delta_summary_min"]))
        ys.append(100 * num(row["delta_top10_single_target_fraction"]))
    pad_x = 0.08
    pad_y = 8
    for ax, method in zip(axes, ["M1", "M1b", "M2", "M3"]):
        for index, pair in enumerate(ORDER, start=1):
            row = by[pair, method]
            state = row["concordance"]
            x = num(row["delta_summary_min"])
            y = 100 * num(row["delta_top10_single_target_fraction"])
            ax.scatter(
                x, y, s=28, color=STATE_COLOR[state], marker=STATE_MARKER[state],
                facecolors="none" if method == "M1b" else STATE_COLOR[state],
                linewidths=0.9, zorder=3,
            )
            dx, dy = (5, 5) if index % 2 else (-7, -8)
            ax.annotate(str(index), (x, y), textcoords="offset points", xytext=(dx, dy), fontsize=6.5, ha="center")
        ax.axhline(0, color="0.55", lw=0.6)
        ax.axvline(0, color="0.55", lw=0.6)
        ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
        ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)
        title = f"{method} vs M0"
        if method == "M1b":
            title += " (supporting; open symbols)"
        ax.set_title(title)
        ax.set_xlabel("delta summary_min")
        ax.set_ylabel("Top-k single-target change (pp)")
    save(fig, "Fig06_Directional_and_Topk_Changes")


def descriptor_limits(d1):
    limits = {}
    for key, _label in DESCRIPTORS:
        vals = [num(r[key]) for r in d1]
        pad = 0.05 * (max(vals) - min(vals) + 1e-9)
        limits[key] = (min(vals) - pad, max(vals) + pad)
    return limits


def jitter(ligand_id: str, slot: int) -> float:
    raw = zlib.crc32(f"{ligand_id}:{slot}".encode()) % 1000
    return (raw / 999.0 - 0.5) * 0.28


def draw_boxes(ax, groups, positions):
    for pos, values in zip(positions, groups):
        if not values:
            continue
        arr = np.asarray(values, dtype=float)
        q1, med, q3 = np.percentile(arr, [25, 50, 75], method="linear")
        iqr = q3 - q1
        inside = arr[(arr >= q1 - 1.5 * iqr) & (arr <= q3 + 1.5 * iqr)]
        low = float(inside.min()) if len(inside) else float(q1)
        high = float(inside.max()) if len(inside) else float(q3)
        ax.plot([pos, pos], [low, q1], color="0.2", lw=0.6)
        ax.plot([pos, pos], [q3, high], color="0.2", lw=0.6)
        ax.plot([pos - 0.12, pos + 0.12], [low, low], color="0.2", lw=0.6)
        ax.plot([pos - 0.12, pos + 0.12], [high, high], color="0.2", lw=0.6)
        ax.add_patch(plt.Rectangle((pos - 0.16, q1), 0.32, q3 - q1, facecolor="#D6E6F5", edgecolor="0.15", lw=0.6))
        ax.plot([pos - 0.16, pos + 0.16], [med, med], color="0.1", lw=0.8)


def fig_descriptors(d1, arm: str, positive_class: str, stem: str, title: str) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(PAGE_W, 8.4), layout="constrained")
    limits = descriptor_limits(d1)
    subset = [r for r in d1 if r["arm"] == arm and r["class"] in {"dual", positive_class}]
    for ax, (key, label) in zip(axes.ravel(), DESCRIPTORS):
        for pair_index, pair in enumerate(ORDER):
            for slot, cls, color in ((0, "dual", C_BLUE), (1, positive_class, C_ORANGE)):
                rows = [r for r in subset if r["pair"] == pair and r["class"] == cls]
                vals = [num(r[key]) for r in rows]
                pos = pair_index + (-0.16 if slot == 0 else 0.16)
                draw_boxes(ax, [vals], [pos])
                ax.scatter(
                    [pos + jitter(r["canonical_ligand_id"], slot) for r in rows],
                    vals, s=6, color=color, alpha=0.75, linewidths=0, zorder=3,
                )
        ax.set_ylim(*limits[key])
        ax.set_title(label, fontsize=8)
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels([str(i) for i in range(1, 9)], fontsize=6.5)
    axes[0, 0].scatter([], [], s=12, color=C_BLUE, label="dual")
    axes[0, 0].scatter([], [], s=12, color=C_ORANGE, label=positive_class)
    axes[0, 0].legend(frameon=False, fontsize=6.5)
    fig.suptitle(title, fontsize=9)
    save(fig, stem)


def fig_s3(d1) -> None:
    fig, axes = new_fig(4.2, 1, 2)
    for ax, arm in zip(axes, ARMS):
        subset = [r for r in d1 if r["arm"] == arm]
        for pair_index, pair in enumerate(ORDER):
            rows = [r for r in subset if r["pair"] == pair]
            vals = [num(r["max_ecfp4_tanimoto_other_class"]) for r in rows]
            draw_boxes(ax, [vals], [pair_index])
            ax.scatter(
                [pair_index + jitter(r["canonical_ligand_id"], 0) for r in rows],
                vals, s=8, color=C_BLUE, alpha=0.7, linewidths=0,
            )
        ax.set_ylim(0, 1)
        ax.set_title(arm)
        ax.set_ylabel("max ECFP4 Tanimoto to other class")
        ax.set_xticks(range(8))
        ax.set_xticklabels([f"{i}" for i in range(1, 9)])
    save(fig, "FigS03_Cross_Class_Similarity")


def fig_s4(directional) -> None:
    fig, axes = new_fig(5.2, 1, 3)
    contrasts = [
        ("dual_vs_A_only", "D_vs_A"),
        ("dual_vs_B_only", "D_vs_B"),
        ("summary_min", "summary_min"),
    ]
    for ax, (contrast, title) in zip(axes, contrasts):
        grid = []
        for pair in ORDER:
            grid.append([
                num(next(r for r in directional if r["pair"] == pair and r["method"] == method and r["contrast"] == contrast)["auroc"])
                for method in METHODS
            ])
        image = np.asarray(grid)
        im = ax.imshow(image, vmin=0, vmax=1, cmap="cividis", aspect="auto")
        for i in range(image.shape[0]):
            for j in range(image.shape[1]):
                ax.text(j, i, f"{image[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color="white" if image[i, j] < 0.45 else "black")
        ax.set_xticks(range(5))
        ax.set_xticklabels(METHODS)
        ax.set_yticks(range(8))
        ax.set_yticklabels([f"{i} {p}" for i, p in enumerate(ORDER, start=1)], fontsize=6)
        ax.set_title(title)
    fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02, label="AUROC")
    save(fig, "FigS04_All_Method_AUROC")


def fig_s5(seeds, summary) -> None:
    fig, axes = new_fig(7.2, 3, 1)
    metrics = [("auc_B", "auc_B = D_vs_A"), ("auc_A", "auc_A = D_vs_B"), ("summary_min", "summary_min")]
    x = np.arange(len(ORDER))
    seed_by = {(r["pair"], int(r["seed"])): r for r in seeds}
    for ax, (metric, title) in zip(axes, metrics):
        for offset, seed in zip(np.linspace(-0.24, 0.24, 5), SEEDS):
            ys = [num(seed_by[pair, seed][metric]) for pair in ORDER]
            ax.scatter(x + offset, ys, s=14, color=SEED_COLOR[seed], marker=SEED_MARKER[seed], label=f"seed {seed}")
        rows = [next(r for r in summary if r["pair"] == p and r["metric"] == metric) for p in ORDER]
        ax.vlines(x, [num(r["min"]) for r in rows], [num(r["max"]) for r in rows], color="0.3", lw=0.7)
        ax.scatter(x, [num(r["median"]) for r in rows], marker="_", s=70, color="black", label="stored median")
        ax.set_ylim(0, 1)
        ax.axhline(0.5, color="0.7", lw=0.5)
        ax.set_title(title + "\nstored min–max; IQR is annotated, not drawn as a CI")
        ax.set_ylabel("AUROC")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{i}\nIQR {num(r['iqr']):.3f}" for i, r in enumerate(rows, start=1)], fontsize=6)
        ax.legend(frameon=False, ncol=3, fontsize=6)
    save(fig, "FigS05_Five_Seeds")


def fig_s6(alt) -> None:
    fig, ax = new_fig(4.4)
    rows = sorted(alt, key=lambda r: (ORDER.index(r["pair"]), r["replaced_side"]))
    y = np.arange(len(rows))[::-1]
    series = [
        ("delta_AUC_affected", "affected side", C_BLUE, "o", 0.18),
        ("delta_AUC_unreplaced", "unreplaced side", "#E69F00", "s", 0.0),
        ("delta_summary_min", "summary_min", C_PINK, "D", -0.18),
    ]
    for field, label, color, marker, offset in series:
        points = [num(r[field]) for r in rows]
        lo = [num(r[field + "_ci_lo"]) for r in rows]
        hi = [num(r[field + "_ci_hi"]) for r in rows]
        ax.hlines(y + offset, lo, hi, color=color, lw=0.7)
        ax.scatter(points, y + offset, marker=marker, color=color, s=14, label=label, zorder=2)
    ax.axvline(0, color="0.45", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['pair']} {r['alt_id']} side {r['replaced_side']}" for r in rows], fontsize=6.5)
    ax.set_xlabel("alternate − primary")
    ax.set_title("Stored paired 95% CIs; unreplaced zero is a construction check")
    ax.legend(frameon=False)
    save(fig, "FigS06_Alternative_Receptors")


def fig_s7(rows) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(PAGE_W, 8.2), layout="constrained")
    classes = ["dual", "A_only", "B_only", "neither"]
    targets = ["A", "B"]
    by = {(r["pair"], r["method"], r["class"], r["target"]): r for r in rows}
    panels = [(cls, target) for cls in classes for target in targets]
    for ax, (cls, target) in zip(axes.ravel(), panels):
        ax.set_axis_on()
        grid = np.zeros((8, 5))
        annot = []
        for i, pair in enumerate(ORDER):
            row_ann = []
            for j, method in enumerate(METHODS):
                rec = by[pair, method, cls, target]
                expected = num(rec["expected"])
                valid = num(rec["valid"])
                grid[i, j] = valid / expected if expected else np.nan
                row_ann.append(f"{int(valid)}/{int(expected)}")
            annot.append(row_ann)
        im = ax.imshow(grid, vmin=0, vmax=1, cmap="Blues", aspect="auto")
        for i in range(8):
            for j in range(5):
                ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=4.8, color="black")
        ax.set_title(f"{cls}, target {target}", fontsize=8)
        ax.set_xticks(range(5))
        ax.set_xticklabels(METHODS, fontsize=6)
        ax.set_yticks(range(8))
        ax.set_yticklabels([str(i) for i in range(1, 9)], fontsize=6)
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.01, label="valid/expected")
    save(fig, "FigS07_Missingness")


def fig_s8(desc, concordance) -> None:
    fig, axes = new_fig(4.6, 1, 2)
    m0 = [r for r in desc if r["method"] == "M0" and r["population_kind"] == "method_specific"]
    by = {r["pair"]: r for r in m0}
    x = np.arange(len(ORDER))
    dual = [int(float(by[p]["n_dual_top"])) for p in ORDER]
    a_only = [int(float(by[p]["n_A_only_top"])) for p in ORDER]
    b_only = [int(float(by[p]["n_B_only_top"])) for p in ORDER]
    axes[0].bar(x, dual, color=C_BLUE, label="dual")
    axes[0].bar(x, a_only, bottom=dual, color=C_ORANGE, label="A_only")
    axes[0].bar(x, b_only, bottom=[u + v for u, v in zip(dual, a_only)], color=C_GREEN, label="B_only")
    for i, pair in enumerate(ORDER):
        axes[0].text(i, int(float(by[pair]["k"])) + 0.15, f"k={by[pair]['k']}", ha="center", fontsize=6)
    axes[0].set_ylabel("ligands in Top-k")
    axes[0].set_title("M0 method-specific lists")
    axes[0].legend(frameon=False, fontsize=6.5)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([str(i) for i in range(1, 9)])
    methods = ["M1", "M1b", "M2", "M3"]
    counts = {r["method"]: r for r in concordance}
    bottoms = np.zeros(4)
    colors = [STATE_COLOR[s] for s in STATES]
    labels = ["Improve", "Worsen", "Discordant", "No change"]
    for state, color, label in zip(STATES, colors, labels):
        heights = [int(counts[m][state]) for m in methods]
        axes[1].bar(["M1", "M1b*", "M2", "M3"], heights, bottom=bottoms, color=color, label=label)
        bottoms = bottoms + np.asarray(heights)
    axes[1].set_ylim(0, 8.6)
    axes[1].set_ylabel("pairs")
    axes[1].set_title("Stored four-state counts")
    axes[1].legend(frameon=False, fontsize=6.5)
    save(fig, "FigS08_Topk_and_Concordance")


def main() -> None:
    style()
    directional = read_csv(PRIMARY / "PRIMARY_DIRECTIONAL_METRICS.csv")
    deltas = read_csv(PRIMARY / "PRIMARY_METHOD_DELTA_METRICS.csv")
    negative = read_csv(RESULTS / "module_B_negative_class.csv")
    pocket = read_csv(RESULTS / "module_C_wrong_pocket.csv")
    chemistry = read_csv(RESULTS / "module_D_chemistry_units.csv")
    d1 = read_csv(RESULTS / "module_D1_composition.csv")
    seeds = read_csv(RESULTS / "module_E_five_seed.csv")
    summary = read_csv(RESULTS / "module_E_five_seed_summary.csv")
    alt = read_csv(RESULTS / "module_F_alt_receptor.csv")
    missing = read_csv(RESULTS / "module_H_missingness.csv")
    desc = read_csv(RESULTS / "module_I_descriptive.csv")
    compare = read_csv(RESULTS / "module_I_method_compare.csv")
    concordance = read_csv(CONCORDANCE)
    scheme()
    fig01(directional, deltas)
    fig02(negative)
    fig03(chemistry)
    fig04(pocket)
    fig05(seeds, summary, alt)
    fig06(compare)
    fig_descriptors(d1, "D_vs_A", "A_only", "FigS01_Descriptors_D_vs_A", "D_vs_A descriptors on layer 1; dual and A_only")
    fig_descriptors(d1, "D_vs_B", "B_only", "FigS02_Descriptors_D_vs_B", "D_vs_B descriptors on layer 1; dual and B_only")
    fig_s3(d1)
    fig_s4(directional)
    fig_s5(seeds, summary)
    fig_s6(alt)
    fig_s7(missing)
    fig_s8(desc, concordance)
    print("FIGURES_WRITTEN", OUT)
    print("CONCORDANCE_CSV_UNCHANGED", CONCORDANCE)


if __name__ == "__main__":
    main()
