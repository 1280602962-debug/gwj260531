#!/usr/bin/env python3
"""Publication figures from the computed tables. English labels; captions live in the Chinese report."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parents[1]
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def figure1():
    fig, ax = plt.subplots(figsize=(11.2, 4.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Two sources, not one funnel", loc="left", fontsize=12)

    def box(x, y, text, fc):
        ax.add_patch(plt.Rectangle((x, y), 2.1, 0.85, facecolor=fc, edgecolor="#333333", lw=0.6))
        ax.text(x + 1.05, y + 0.42, text, ha="center", va="center", fontsize=7.5)

    hist = [
        (0.3, 4.4, "Clinical library\n8319", "#d9e2ef"),
        (2.6, 4.4, "v1 p ≥ 0.5\n1588", "#d9e2ef"),
        (4.9, 4.4, "Acid equivalent\n303", "#d9e2ef"),
        (7.2, 4.4, "Docked\n156", "#d9e2ef"),
        (0.3, 2.6, "Dual structural\n≥2/3 seeds: 54", "#f4e1c1"),
        (2.6, 2.6, "∩ seed42-eligible\nchemistry: 37", "#f4e1c1"),
        (4.9, 2.6, "Frozen primary\n12", "#d8ead3"),
        (7.2, 2.6, "Frozen backup\n21", "#d8ead3"),
    ]
    for x, y, t, c in hist:
        box(x, y, t, c)
    ax.annotate("", xy=(2.55, 4.82), xytext=(2.45, 4.82), arrowprops=dict(arrowstyle="->", color="#444444"))
    ax.text(0.3, 5.45, "Historical v1 path", fontsize=9, fontweight="bold")
    ax.text(9.6, 4.6, "Excluded at freeze,\nonce: GSK-3008348\n+ 3 cephalosporins", fontsize=7.5, color="#333333")
    ax.text(9.6, 2.7, "Not renominated:\nAZD1981, sodelglitazar\n(seed 42 failed)", fontsize=7.5)
    ax.text(0.3, 1.55, "Stage 3B, separate v2-pool entrants", fontsize=9, fontweight="bold")
    box(0.3, 0.35, "New docked\n32", "#f7d6d0")
    box(2.6, 0.35, "Structural +\nchemistry 13", "#f7d6d0")
    box(4.9, 0.35, "Structural,\nchemistry fail 8", "#f7d6d0")
    box(7.2, 0.35, "Not structural\n11", "#f7d6d0")
    ax.text(9.6, 0.55, "Rule not applicable\n7 non-acids\n(not docked)", fontsize=7.5)
    save(fig, "figure1_study_flow")


def figure2():
    ref = pd.read_csv(OUT / "REFERENCE_DOCKING_CONSOLIDATED.csv")
    met = pd.read_csv(OUT / "AUROC_PRAUC_EF_ALL.csv")
    met = met[(met.population == "method_specific")]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    ax = axes[0]
    labels, top1, best = [], [], []
    for _, r in ref.iterrows():
        labels.append(f"{r.ligand[:4]}\n{int(r.seed)}")
        top1.append(r.top1_pose_rmsd_A)
        best.append(r.best_of_9_pose_rmsd_A)
    x = range(len(labels))
    ax.bar([i - 0.15 for i in x], top1, width=0.3, label="Selected pose", color="#4c78a8")
    ax.bar([i + 0.15 for i in x], best, width=0.3, label="Best of 9", color="#f2a07a")
    ax.axhline(2.0, color="#333333", ls="--", lw=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Heavy-atom RMSD (Å)")
    ax.set_title("Reference docking")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    order = ["P0", "P1", "P2", "P3", "P4", "P5"]
    true = met[met.benchmark == "true_decoy_benchmark"].set_index("protocol").loc[order]
    rand = met[met.benchmark == "random_decoy_benchmark"].set_index("protocol").loc[order]
    x = range(6)
    ax.bar([i - 0.15 for i in x], true.EF1, width=0.3, label="True-decoy set", color="#4c78a8")
    ax.bar([i + 0.15 for i in x], rand.EF1, width=0.3, label="Random-decoy set", color="#72b7b2")
    ax.axhline(1.0, color="#333333", ls="--", lw=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(order)
    ax.set_ylabel("EF at 1%  (k = floor(0.01 n))")
    ax.set_title("Scoring, pre-specified Top-K")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "figure2_docking_and_scores")
    true.reset_index()[["protocol", "EF1", "AUROC", "PR_AUC"]].to_csv(TAB / "figure2_true_decoy_metrics.csv", index=False)


def figure3():
    g = pd.read_csv(OUT / "CANDIDATE_GROUPS_FINAL.csv")
    labels = {
        "1_frozen_primary_12": "Frozen primary",
        "2_frozen_backup_21": "Frozen backup",
        "3_stage3b_new_13_not_historical_nomination": "Stage 3B new,\nnot renominated",
        "4_historical_structural_chemistry_clean_not_in_seed42_eligible": "Chemistry-clean,\nseed 42 failed",
        "5_structural_control_not_a_candidate": "Structural control",
        "5_structural_pass_chemistry_not_pass": "Structural pass,\nchemistry fail",
        "6_structural_not_pass": "Not structural",
        "7_rule_not_applicable": "Rule not applicable",
    }
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    ys = list(range(len(g)))
    colors = ["#2f6f4e" if i < 2 else "#b85c38" if i == 2 else "#8a8a8a" for i in range(len(g))]
    ax.barh(ys, g["n"], color=colors)
    ax.set_yticks(ys)
    ax.set_yticklabels([labels.get(x, x) for x in g["group"]], fontsize=8)
    ax.invert_yaxis()
    for y, n in zip(ys, g["n"]):
        ax.text(n + 1, y, str(int(n)), va="center", fontsize=8)
    ax.set_xlabel("Molecules")
    ax.set_title("Candidate strata. Green bars are the historical freeze only.")
    fig.tight_layout()
    save(fig, "figure3_candidate_strata")
    g.to_csv(TAB / "figure3_candidate_counts.csv", index=False)


def figure4():
    sc = pd.read_csv(OUT / "SCAFFOLD_COUNTS_BY_COHORT.csv")
    want = ["docked_156", "frozen_33", "historical_structural_21", "stage3b_new32", "stage3b_pass_13", "comparison_set"]
    sub = sc.set_index("cohort").loc[want]
    ctrl = pd.read_csv(OUT / "supplementary/SCAFFOLD_SAMPLE_SIZE_CONTROLS.csv")
    ctrl = ctrl[ctrl.parent == "clinical_8319"].groupby("n").n_murcko.mean()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    ax = axes[0]
    ax.bar(range(len(sub)), sub.n_murcko, color="#4c78a8")
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(["Docked 156", "Frozen 33", "Historical 21", "Stage 3B 32", "New pass 13", "Comparison 48"], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Unique Murcko scaffolds")
    ax.set_title("Scaffold counts")
    ax = axes[1]
    # new scaffolds relative to frozen 33 / docked 156
    labels = ["13 vs frozen 33", "13 vs docked 156", "32 vs docked 156"]
    new = [11, 11, 27]
    shared = [0, 0, 1]
    x = range(3)
    ax.bar(x, new, label="Murcko scaffolds absent from the reference", color="#4c78a8")
    ax.bar(x, shared, bottom=new, label="Shared Murcko scaffolds", color="#f2a07a")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Scaffolds")
    ax.set_title("New relative to an existing set")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "figure4_scaffold_diversity")
    sub.reset_index().to_csv(TAB / "figure4_scaffold_counts.csv", index=False)


def figure5():
    bias = pd.read_csv(OUT / "GATE_SCAFFOLD_BIAS.csv")
    order = [
        "pre_dock_acid_303", "docked_156", "urat1_a2_ge2", "nlrp3_structural_ge2",
        "dual_structural_ge2", "frozen_reportable_33",
    ]
    b = bias.set_index("cohort").loc[order]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    labels = ["Acid 303", "Docked 156", "URAT1 A2\n≥2 seeds", "NLRP3 gate\n≥2 seeds", "Dual gate\n≥2 seeds", "Frozen 33"]
    ax = axes[0]
    ax.plot(range(len(b)), b.fraction_n_carboxylate_ge1, marker="o", color="#4c78a8")
    ax.set_ylim(0.7, 1.05)
    ax.set_xticks(range(len(b)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Fraction with ≥1 carboxylate")
    ax.set_title("Historical path only")
    ax = axes[1]
    ax.plot(range(len(b)), b.median_mw, marker="o", color="#b85c38")
    ax.set_xticks(range(len(b)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Median molecular weight")
    ax.set_title("Size along the same path")
    fig.tight_layout()
    save(fig, "figure5_gate_bias")
    b.reset_index().to_csv(TAB / "figure5_gate_bias.csv", index=False)


def main():
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    print("figures", sorted(p.name for p in FIG.glob("figure*")))


if __name__ == "__main__":
    main()
