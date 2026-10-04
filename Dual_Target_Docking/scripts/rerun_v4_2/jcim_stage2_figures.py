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
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.6))
    x = range(len(ORDER))
    m0 = [r for r in directional if r["method"] == "M0" and r["contrast"] == "summary_min"]
    a = {
        r["pair"]: fnum(r["auroc"])
        for r in directional
        if r["method"] == "M0" and r["contrast"] == "dual_vs_A_only"
    }
    b = {
        r["pair"]: fnum(r["auroc"])
        for r in directional
        if r["method"] == "M0" and r["contrast"] == "dual_vs_B_only"
    }
    axes[0].plot(list(x), [a[p] for p in ORDER], "o-", label="dual vs A_only (score_B)")
    axes[0].plot(list(x), [b[p] for p in ORDER], "s-", label="dual vs B_only (score_A)")
    sm = {r["pair"]: fnum(r["auroc"]) for r in m0}
    axes[0].plot(list(x), [sm[p] for p in ORDER], "^-", label="summary_min")
    axes[0].axhline(0.5, color="0.6", lw=0.8)
    axes[0].set_ylim(0, 1)
    axes[0].set_title("PRIMARY M0, shared-dual universe")
    axes[0].set_ylabel("AUROC")
    axes[0].legend(fontsize=8)

    delta_names = ["M1-M0", "M1b-M1", "M2-M0", "M3-M0"]
    markers = ["o", "s", "^", "D"]
    for name, marker in zip(delta_names, markers):
        lookup = {r["pair"]: fnum(r["delta_summary_min"]) for r in deltas if r["delta"] == name}
        axes[1].plot(list(x), [lookup[p] for p in ORDER], marker=marker, linestyle="-", label=name)
    axes[1].axhline(0, color="0.6", lw=0.8)
    axes[1].set_title("PRIMARY paired summary_min deltas")
    axes[1].set_ylabel("delta summary_min")
    axes[1].legend(fontsize=8)
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
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8))
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
    axes[0].bar(x, dual, label="dual")
    axes[0].bar(x, a_only, bottom=dual, label="A_only")
    bottom2 = [u + v for u, v in zip(dual, a_only)]
    axes[0].bar(x, b_only, bottom=bottom2, label="B_only")
    for i, pair in enumerate(ORDER):
        axes[0].text(i, int(float(by[pair]["k"])) + 0.15, f"k={by[pair]['k']}", ha="center", fontsize=7)
    axes[0].set_title("Module I, M0 Top-k class counts")
    axes[0].set_ylabel("ligands in Top-k")
    axes[0].legend(fontsize=8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([SHORT[p] for p in ORDER], rotation=35, ha="right")

    states = ["CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE"]
    methods = ["M1", "M1b*", "M2", "M3"]
    method_keys = ["M1", "M1b", "M2", "M3"]
    bottoms = [0] * len(methods)
    for state in states:
        heights = []
        for key in method_keys:
            heights.append(sum(r["method"] == key and r["concordance"] == state for r in compare))
        axes[1].bar(methods, heights, bottom=bottoms, label=state)
        bottoms = [b + h for b, h in zip(bottoms, heights)]
    axes[1].set_ylim(0, 8.5)
    axes[1].set_ylabel("pairs")
    axes[1].set_title("Four-state counts (* M1b supporting only)")
    axes[1].legend(fontsize=7)
    save(fig, "fig_I_topk_and_concordance.png")


def fig_e_f(seeds: list[dict], summary: list[dict], alt: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8))
    x = list(range(len(ORDER)))
    # Five stored seed points. The written IQR stays a number under each pair.
    colors = {"auc_B": "#1f4e79", "summary_min": "#b85c38"}
    for metric, marker in (("auc_B", "o"), ("summary_min", "^")):
        xs, ys = [], []
        for i, pair in enumerate(ORDER):
            pts = [fnum(r[metric]) for r in seeds if r["pair"] == pair]
            xs.extend([i] * len(pts))
            ys.extend(pts)
        axes[0].scatter(xs, ys, marker=marker, s=22, color=colors[metric], label=metric, zorder=2)
        meds = [
            fnum(next(r for r in summary if r["pair"] == pair and r["metric"] == metric)["median"])
            for pair in ORDER
        ]
        axes[0].scatter(x, meds, marker="D", s=28, facecolors="none", edgecolors=colors[metric], linewidths=1.2, zorder=3, label=f"{metric} median")
    iqr_labels = []
    for pair in ORDER:
        row = next(r for r in summary if r["pair"] == pair and r["metric"] == "auc_B")
        iqr_labels.append(f"IQR {fnum(row['iqr']):.3f}")
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Module E stored seed points (IQR annotated, not a CI)")
    axes[0].set_ylabel("AUROC")
    axes[0].legend(fontsize=7)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(
        [f"{SHORT[p]}\n{lab}" for p, lab in zip(ORDER, iqr_labels)], rotation=35, ha="right", fontsize=7
    )

    alt_rows = sorted(alt, key=lambda r: (ORDER.index(r["pair"]), r["replaced_side"]))
    xf = list(range(len(alt_rows)))
    y = [fnum(r["delta_AUC_affected"]) for r in alt_rows]
    lo = [fnum(r["delta_AUC_affected"]) - fnum(r["delta_AUC_affected_ci_lo"]) for r in alt_rows]
    hi = [fnum(r["delta_AUC_affected_ci_hi"]) - fnum(r["delta_AUC_affected"]) for r in alt_rows]
    axes[1].errorbar(xf, y, yerr=[lo, hi], fmt="o", capsize=3, label="affected-side delta CI")
    axes[1].plot(xf, [fnum(r["delta_AUC_unreplaced"]) for r in alt_rows], "s", label="unreplaced delta")
    axes[1].axhline(0, color="0.6", lw=0.8)
    axes[1].set_title("Module F alt-receptor deltas")
    axes[1].set_ylabel("AUROC delta")
    axes[1].legend(fontsize=8)
    axes[1].set_xticks(xf)
    axes[1].set_xticklabels(
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
