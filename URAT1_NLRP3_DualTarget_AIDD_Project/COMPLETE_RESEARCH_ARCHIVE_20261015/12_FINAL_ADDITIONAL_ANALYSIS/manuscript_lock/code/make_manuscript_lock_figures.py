#!/usr/bin/env python3
"""Redraw manuscript figures from deposited tables. Writes only under manuscript_lock/."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[4]
ANALYSIS = Path(__file__).resolve().parents[2]
LOCK = Path(__file__).resolve().parents[1]
FIG = LOCK / "figures"
TAB = LOCK / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

AD_THRESHOLD = 0.5672473867595819
PAL = {
    "primary": "#0072B2",
    "backup": "#E69F00",
    "new": "#009E73",
    "selected": "#0072B2",
    "best": "#D55E00",
    "hist": "#D6E3F0",
    "freeze": "#C8E6C9",
    "exclude": "#F4C7C3",
    "stage": "#F6E2C0",
    "neutral": "#E6E6E6",
    "ink": "#222222",
}


def save(fig, name: str) -> None:
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def style_ax(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def box(ax, x, y, w, h, text, fc, fs=6.6, ha="center"):
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.06",
            facecolor=fc,
            edgecolor=PAL["ink"],
            linewidth=0.6,
            mutation_aspect=0.6,
        )
    )
    tx = x + w / 2 if ha == "center" else x + 0.18
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, color=PAL["ink"], linespacing=1.2)


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.7,
            color="#444444",
        )
    )


def figure1() -> None:
    fig, ax = plt.subplots(figsize=(12.2, 9.0))
    ax.set_xlim(0, 13.15)
    ax.set_ylim(0, 11.35)
    ax.axis("off")
    ax.text(0.15, 11.05, "Historical path and a separate Stage 3B branch", fontsize=12, color=PAL["ink"], fontweight="bold")
    ax.text(
        0.15,
        10.68,
        "The 12 and 21 were fixed after structural evaluation of the historical 156, and before Stage 3B.",
        fontsize=8,
        color=PAL["ink"],
    )

    w, h, gap = 1.95, 1.05, 0.18
    xs = [0.15 + i * (w + gap) for i in range(6)]
    y1 = 9.3
    labels = [
        "Clinical library\n8319",
        "v1 output ≥ 0.5\n1588 members",
        "Acid equivalent\n303",
        "Docked\n156",
        "Dual structural\n≥2 of 3 seeds: 54",
        "and seed-42\neligible set: 37",
    ]
    for x, text in zip(xs, labels):
        box(ax, x, y1, w, h, text, PAL["hist"], fs=7.0)
    for i in range(5):
        arrow(ax, xs[i] + w + 0.01, y1 + h / 2, xs[i + 1] - 0.01, y1 + h / 2)

    ex_y, ex_h = 6.55, 2.15
    box(
        ax,
        0.15,
        ex_y,
        12.85,
        ex_h,
        "Removed once, 37 to 33\n"
        "REP_07907   GSK-3008348 FREE BASE    structural control\n"
        "REP_02721   CEFETRIZOLE    β-lactam\n"
        "REP_04581   CEFAZEDONE    β-lactam\n"
        "REP_04623   CEFOXAZOLE    β-lactam",
        PAL["exclude"],
        fs=7.6,
        ha="left",
    )
    src = xs[5] + w / 2
    arrow(ax, src, y1, src, ex_y + ex_h)

    box(ax, 1.3, 4.85, 4.6, 1.15, "Frozen primary 12\nnot a prospective validation set", PAL["freeze"], fs=7.4)
    box(ax, 7.15, 4.85, 4.6, 1.15, "Frozen backup 21\nnot a prospective validation set", PAL["freeze"], fs=7.4)
    arrow(ax, 3.6, ex_y, 3.6, 6.0)
    arrow(ax, 9.45, ex_y, 9.45, 6.0)

    ax.plot([0.15, 13.0], [4.5, 4.5], color="#888888", lw=0.6)
    ax.text(0.15, 4.15, "Stage 3B branch — not added to the frozen 12 or 21", fontsize=8.5, color=PAL["ink"], fontweight="bold")
    yb, hb, wb = 2.25, 1.4, 2.85
    stage_x = [0.15, 3.35, 6.55, 9.75]
    stage_labels = [
        ("Docked 32\n3 seeds × 2 targets", PAL["stage"]),
        ("Structural and\nchemistry pass\n13", "#B7E1CD"),
        ("Structural pass,\nchemistry not clean\n8", PAL["stage"]),
        ("Not structurally\nqualified\n11", PAL["neutral"]),
    ]
    for x, (text, fc) in zip(stage_x, stage_labels):
        box(ax, x, yb, wb, hb, text, fc, fs=7.2)
    for i in range(3):
        arrow(ax, stage_x[i] + wb + 0.02, yb + hb / 2, stage_x[i + 1] - 0.02, yb + hb / 2)
    box(ax, 0.15, 0.3, 4.15, 1.55, "Rule not applicable\n7 non-acids\nnot docked", PAL["neutral"], fs=7.2)
    ax.text(
        4.6,
        1.05,
        "13 + 8 + 11 = 32 docked molecules. No ID overlap with the frozen 33.",
        fontsize=7.6,
        color=PAL["ink"],
        va="center",
    )
    save(fig, "figure1_flow")


def figure2(ref: pd.DataFrame) -> None:
    order = ["lesinurad", "benzbromarone", "NP3-146"]
    titles = {
        "lesinurad": "Lesinurad, 9DKB",
        "benzbromarone": "Benzbromarone, 9DKA",
        "NP3-146": "NP3-146, 7ALV",
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.6), sharey=True)
    for ax, ligand in zip(axes, order):
        sub = ref[ref.ligand == ligand].sort_values("seed")
        x = np.arange(len(sub))
        ax.bar(x - 0.16, sub.top1_pose_rmsd_A, width=0.32, color=PAL["selected"], label="CNNscore-selected pose")
        ax.bar(x + 0.16, sub.best_of_9_pose_rmsd_A, width=0.32, color=PAL["best"], label="Best of 9")
        ax.axhline(2.0, color=PAL["ink"], ls="--", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"seed {int(s)}" for s in sub.seed], fontsize=8)
        ax.set_title(titles[ligand], fontsize=10)
        style_ax(ax)
        ax.set_ylim(0, 5.6)
    axes[0].set_ylabel("Heavy-atom RMSD in the receptor frame (Å)")
    axes[2].legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.suptitle("Production-protocol reference docking", fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "figure2_reference_rmsd")


def figure3() -> None:
    # Rows are predicted positive / negative; columns are labeled active / decoy.
    mats = {
        "A1, not the nomination rule": np.array([[102, 13], [126, 51]]),
        "A2, formal URAT1 entry rule": np.array([[217, 61], [11, 3]]),
    }
    stats = {
        "A1, not the nomination rule": "Sensitivity 0.447   Specificity 0.797   LR+ 2.20",
        "A2, formal URAT1 entry rule": "Sensitivity 0.952   Specificity 0.047   LR+ 0.999",
    }
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    for ax, (title, mat) in zip(axes, mats.items()):
        im = ax.imshow(mat, cmap="Blues", vmin=0, vmax=220)
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(["Labeled active\n228", "Decoy\n64"], fontsize=8)
        ax.set_yticklabels(["Gate pass", "Gate fail"], fontsize=8)
        for i in range(2):
            for j in range(2):
                value = int(mat[i, j])
                ax.text(j, i, str(value), ha="center", va="center", fontsize=12, color="white" if value >= 150 else "#111111")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(stats[title], fontsize=8)
        ax.tick_params(length=0)
    fig.suptitle("Retrospective URAT1 gate on 292 molecules (seed 42)", fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.02, 1, 0.90))
    save(fig, "figure3_a1_a2")


def morgan(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def figure4(coords: pd.DataFrame, tanimoto: pd.DataFrame, scaffolds: pd.DataFrame, var: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.0), gridspec_kw={"width_ratios": [1.25, 0.85, 1.15]})
    ax = axes[0]
    for group, color, label in (
        ("primary_12", PAL["primary"], "Frozen primary, 12"),
        ("backup_21", PAL["backup"], "Frozen backup, 21"),
        ("stage3b_13", PAL["new"], "Stage 3B pass, 13"),
    ):
        sub = coords[coords.group == group]
        ax.scatter(sub.pc1, sub.pc2, s=28, c=color, label=label, zorder=3, edgecolors="white", linewidths=0.3)
    ax.set_xlabel(f"PC1 ({100 * var[0]:.1f}% of centered-bit variance)")
    ax.set_ylabel(f"PC2 ({100 * var[1]:.1f}% of centered-bit variance)")
    ax.set_title("Display projection only", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="best")
    ax.text(
        0.02,
        0.02,
        "2D position is not a structural difference",
        transform=ax.transAxes,
        fontsize=7,
        color="#333333",
    )
    style_ax(ax)

    ax = axes[1]
    order = ["primary_12", "backup_21", "stage3b_13"]
    labels = ["Primary 12", "Backup 21", "New 13"]
    sc = scaffolds.set_index("group").loc[order]
    x = np.arange(3)
    ax.bar(x - 0.16, sc.n_molecules, width=0.32, color="#4C4C4C", label="Molecules")
    ax.bar(x + 0.16, sc.n_murcko, width=0.32, color="#56B4E9", label="Murcko scaffolds")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Count")
    ax.set_title("Murcko scaffolds", fontsize=10)
    ax.legend(frameon=False, fontsize=7)
    ax.set_ylim(0, 28)
    style_ax(ax)

    ax = axes[2]
    pieces, names, colors = [], [], []
    for key, name, color in (
        ("within_primary_12", "Within 12", PAL["primary"]),
        ("within_backup_21", "Within 21", PAL["backup"]),
        ("within_stage3b_13", "Within 13", PAL["new"]),
        ("new13_vs_frozen33", "13 vs 33", "#CC79A7"),
    ):
        vals = tanimoto.loc[tanimoto.contrast == key, "tanimoto"].to_numpy()
        pieces.append(vals)
        names.append(name)
        colors.append(color)
    bp = ax.boxplot(pieces, patch_artist=True, showfliers=True, medianprops=dict(color="#111111", lw=1.0), widths=0.55, flierprops=dict(marker="o", markersize=3, markerfacecolor="#555555", markeredgecolor="none", alpha=0.7))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    ax.set_xticklabels(names, fontsize=7.5)
    ax.set_ylabel("Morgan Tanimoto")
    ax.set_title("Fingerprint similarity", fontsize=10)
    ax.set_ylim(0, 1.05)
    style_ax(ax)
    fig.suptitle("Three groups, one fingerprint definition", fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "figure4_chemical_space")


def figure5(clinical: pd.DataFrame, marks: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    vals = clinical.max_tanimoto_train.to_numpy(dtype=float)
    ax.hist(vals, bins=40, color="#9ECAE1", edgecolor="white", linewidth=0.3)
    ax.axvline(AD_THRESHOLD, color="#D55E00", lw=1.2)
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.12)
    ax.text(
        0.62,
        ymax * 0.55,
        "0.5672473868\nstructural-similarity\nthreshold, not an\nactivity threshold\n\n5 of 8319 are\nat or above the line",
        color="#D55E00",
        fontsize=7.2,
        va="center",
    )
    ax.set_xlim(-0.02, 1.05)
    ax.set_xlabel("Maximum Morgan Tanimoto to the 513 training molecules")
    ax.set_ylabel("Clinical-library molecules")
    ax.set_title("8319-molecule distribution", fontsize=10)
    style_ax(ax)

    ax = axes[1]
    rng = np.random.default_rng(17)
    for i, (group, color, label) in enumerate(
        (
            ("primary_12", PAL["primary"], "Primary 12"),
            ("backup_21", PAL["backup"], "Backup 21"),
            ("stage3b_13", PAL["new"], "Stage 3B 13"),
        )
    ):
        sub = marks[marks.group == group]
        jitter = rng.uniform(-0.12, 0.12, size=len(sub))
        ax.scatter(np.full(len(sub), i) + jitter, sub.max_tanimoto_train, s=22, c=color, label=label, zorder=3, edgecolors="white", linewidths=0.3)
    ax.axhline(AD_THRESHOLD, color="#D55E00", lw=1.2)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Primary 12", "Backup 21", "Stage 3B 13"], fontsize=8)
    ax.set_ylabel("Maximum Morgan Tanimoto")
    ax.set_ylim(0, 0.75)
    ax.set_title("Candidates remain below the line", fontsize=10)
    ax.text(2.35, AD_THRESHOLD + 0.02, "0.567", color="#D55E00", fontsize=7.5, ha="right")
    style_ax(ax)
    fig.suptitle("v2 applicability domain is a structural-similarity rule", fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    save(fig, "figure5_applicability_domain")


def figure6(boot: pd.DataFrame) -> None:
    metrics = [
        ("delta_AUROC", "ΔAUROC"),
        ("delta_PR_AUC", "ΔPR-AUC"),
        ("delta_EF1", "ΔEF at 1%"),
        ("delta_EF5", "ΔEF at 5%"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), sharey=True)
    for ax, bench, title in (
        (axes[0], "true_decoy_benchmark", "True-decoy set"),
        (axes[1], "random_decoy_benchmark", "Random-decoy set"),
    ):
        sub = boot[boot.benchmark == bench]
        seeds = [17, 29, 42, 71, 101]
        for i, (col, label) in enumerate(metrics):
            for j, seed in enumerate(seeds):
                row = sub[sub.seed == seed].iloc[0]
                lo = row[f"{col}_ci95_low"]
                hi = row[f"{col}_ci95_high"]
                y = i + (j - 2) * 0.12
                ax.plot([lo, hi], [y, y], color="#0072B2", lw=1.0, alpha=0.85)
                ax.plot(row[col], y, "o", color="#D55E00", ms=3.5)
        ax.axvline(0, color=PAL["ink"], lw=0.8)
        ax.set_yticks(range(4))
        ax.set_yticklabels([m[1] for m in metrics], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("P5 − P2")
        style_ax(ax)
    fig.suptitle("Paired bootstrap, 2000 resamples; five RNG seeds are one analysis", fontsize=11, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    save(fig, "figure6_p2_p5_paired")


def load_groups() -> pd.DataFrame:
    primary = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/primary_candidates.csv")
    backup = pd.read_csv(ROOT / "data/campaigns/c5/04_shortlist_frozen/backup_candidates.csv")
    idcol = "ligand_id" if "ligand_id" in primary.columns else primary.columns[0]
    namecol = "name" if "name" in primary.columns else None
    rows = []
    for df, group in ((primary, "primary_12"), (backup, "backup_21")):
        for _, r in df.iterrows():
            rows.append({"ligand_id": r[idcol], "name": r[namecol] if namecol else "", "group": group})
    evidence = pd.read_csv(ANALYSIS / "CANDIDATE_EVIDENCE_MASTER.csv")
    new = evidence[evidence.group == "3_stage3b_new_13_not_historical_nomination"]
    for _, r in new.iterrows():
        rows.append({"ligand_id": r.ligand_id, "name": r["name"], "group": "stage3b_13"})
    out = pd.DataFrame(rows).drop_duplicates("ligand_id")
    if len(out) != 46:
        raise SystemExit(f"expected 46 molecules in the three groups, found {len(out)}")
    return out.sort_values(["group", "ligand_id"]).reset_index(drop=True)


def build_chemical_space(groups: pd.DataFrame):
    murcko = pd.read_csv(ANALYSIS / "MURCKO_SCAFFOLDS_ALL.csv", usecols=["repurposing_id", "standardized_smiles", "murcko_scaffold"])
    murcko = murcko.drop_duplicates("repurposing_id")
    df = groups.merge(murcko, left_on="ligand_id", right_on="repurposing_id", how="left")
    if df.standardized_smiles.isna().any():
        raise SystemExit("missing SMILES for a locked candidate")
    fps = []
    for smi in df.standardized_smiles:
        fp = morgan(smi)
        if fp is None:
            raise SystemExit("fingerprint failure")
        fps.append(fp)
    bits = np.zeros((len(fps), 2048), dtype=np.float64)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, bits[i])
    centered = bits - bits.mean(axis=0)
    # Column-centered SVD. Row order is group, then ligand_id.
    _, s, vt = np.linalg.svd(centered, full_matrices=False)
    scores = centered @ vt[:2].T
    var = (s ** 2) / (s ** 2).sum()
    coords = df[["ligand_id", "name", "group", "murcko_scaffold"]].copy()
    coords["pc1"] = scores[:, 0]
    coords["pc2"] = scores[:, 1]
    coords["fingerprint"] = "Morgan radius 2, 2048 bits, GetMorganFingerprintAsBitVect"
    coords["reduction"] = "PCA via SVD of the column-centered 46 x 2048 bit matrix; not variance-scaled"
    coords["row_order"] = "group then ligand_id"
    coords["pc1_variance_fraction"] = var[0]
    coords["pc2_variance_fraction"] = var[1]

    pair_rows = []
    index = {g: np.where(df.group.to_numpy() == g)[0] for g in df.group.unique()}
    contrasts = {
        "within_primary_12": (index["primary_12"], index["primary_12"], True),
        "within_backup_21": (index["backup_21"], index["backup_21"], True),
        "within_stage3b_13": (index["stage3b_13"], index["stage3b_13"], True),
        "new13_vs_frozen33": (
            index["stage3b_13"],
            np.concatenate([index["primary_12"], index["backup_21"]]),
            False,
        ),
    }
    for name, (ia, ib, within) in contrasts.items():
        for a in ia:
            for b in ib:
                if within and b <= a:
                    continue
                t = DataStructs.TanimotoSimilarity(fps[a], fps[b])
                pair_rows.append(
                    {
                        "contrast": name,
                        "id_a": df.ligand_id.iloc[a],
                        "id_b": df.ligand_id.iloc[b],
                        "tanimoto": t,
                        "same_murcko": df.murcko_scaffold.iloc[a] == df.murcko_scaffold.iloc[b],
                    }
                )
    tanimoto = pd.DataFrame(pair_rows)
    # Audit against the deposited pairwise table.
    deposited = pd.read_csv(ANALYSIS / "PAIRWISE_CANDIDATE_SIMILARITY.csv")
    dep = deposited.rename(columns={"tanimoto_ecfp4": "tanimoto_deposited"})
    merged = tanimoto.merge(dep[["id_a", "id_b", "tanimoto_deposited"]], on=["id_a", "id_b"], how="left")
    swapped = tanimoto.merge(
        dep[["id_a", "id_b", "tanimoto_deposited"]].rename(columns={"id_a": "id_b", "id_b": "id_a", "tanimoto_deposited": "tanimoto_swap"}),
        on=["id_a", "id_b"],
        how="left",
    )
    check = merged["tanimoto_deposited"].fillna(swapped["tanimoto_swap"])
    if check.isna().any():
        raise SystemExit(f"pairwise coverage gap: {int(check.isna().sum())}")
    max_abs = float(np.max(np.abs(check.to_numpy() - tanimoto.tanimoto.to_numpy())))
    if max_abs > 1e-9:
        raise SystemExit(f"Tanimoto mismatch vs deposited pairs: {max_abs}")

    sc_rows = []
    frozen_scaffolds = set(df.loc[df.group.isin(["primary_12", "backup_21"]), "murcko_scaffold"])
    for group, n_expected in (("primary_12", 12), ("backup_21", 21), ("stage3b_13", 13)):
        sub = df[df.group == group]
        scaffolds = set(sub.murcko_scaffold)
        sc_rows.append(
            {
                "group": group,
                "n_molecules": int(len(sub)),
                "n_murcko": int(len(scaffolds)),
                "n_shared_with_frozen_33": int(len(scaffolds & frozen_scaffolds)) if group == "stage3b_13" else int(len(scaffolds)),
            }
        )
        if len(sub) != n_expected:
            raise SystemExit(group)
    scaffolds = pd.DataFrame(sc_rows)
    if int(scaffolds.loc[scaffolds.group == "stage3b_13", "n_murcko"].iloc[0]) != 11:
        raise SystemExit("expected 11 Murcko scaffolds in the 13")
    if int(scaffolds.loc[scaffolds.group == "stage3b_13", "n_shared_with_frozen_33"].iloc[0]) != 0:
        raise SystemExit("expected 0 scaffolds shared with the frozen 33")
    return coords, tanimoto, scaffolds, var[:2]


def main() -> None:
    figure1()
    ref = pd.read_csv(ANALYSIS / "REFERENCE_DOCKING_CONSOLIDATED.csv")
    figure2(ref)
    figure3()

    groups = load_groups()
    coords, tanimoto, scaffolds, var = build_chemical_space(groups)
    figure4(coords, tanimoto, scaffolds, var)

    scores = pd.read_csv(
        ROOT / "results/nlrp3_model_v2/stage2/clinical_library_scores.csv",
        usecols=["repurposing_id", "name", "max_tanimoto_train", "in_applicability_domain"],
    )
    marks = groups.merge(
        scores.drop(columns=["name"]),
        left_on="ligand_id",
        right_on="repurposing_id",
        how="left",
    )
    if marks.max_tanimoto_train.isna().any():
        raise SystemExit("missing AD similarity for a locked candidate")
    figure5(scores, marks)

    boot = pd.read_csv(ANALYSIS / "P2_P5_PAIRED_BOOTSTRAP.csv")
    figure6(boot)

    # Stage 3B split and PF geometry, copied from deposited columns only.
    master = pd.read_csv(
        ROOT / "COMPLETE_RESEARCH_ARCHIVE_20261015/07_RANKING_AND_SELECTION/ALL_COMPOUNDS_MASTER_INDEX.csv",
        usecols=["repurposing_id", "name", "final_status", "unified_dual_pass_seeds", "in_stage3b_new32_docked", "in_frozen_primary_12", "in_frozen_backup_21"],
    )
    stage = master[master.in_stage3b_new32_docked].copy()
    stage.to_csv(TAB / "stage3b_32_status.csv", index=False)
    urat = pd.read_csv(ANALYSIS / "URAT1_ANCHOR_FEATURES.csv")
    nlrp = pd.read_csv(ANALYSIS / "NLRP3_INTERACTION_MATRIX.csv")
    pf_u = urat[urat.ligand_id == "REP_07580"][
        ["ligand_id", "seed", "acid_arg477_min_A", "centroid_A", "CNNscore", "urat1_pass"]
    ].copy()
    pf_n = (
        nlrp[nlrp.ligand_id == "REP_07580"]
        .groupby("seed", as_index=False)[["CNNscore", "overlap", "ifp_jaccard", "n_key_contacts", "structural_pass"]]
        .first()
        .rename(columns={"CNNscore": "nlrp3_CNNscore"})
    )
    pf = pf_u.merge(pf_n, on="seed")
    pf.to_csv(TAB / "pf03882845_docking_geometry.csv", index=False)

    summary_rows = []
    for contrast, g in tanimoto.groupby("contrast"):
        summary_rows.append(
            {
                "contrast": contrast,
                "n_pairs": int(len(g)),
                "median": float(g.tanimoto.median()),
                "max": float(g.tanimoto.max()),
                "min": float(g.tanimoto.min()),
            }
        )
    tanimoto_summary = pd.DataFrame(summary_rows)
    controls = pd.read_csv(ANALYSIS / "supplementary/SCAFFOLD_SAMPLE_SIZE_CONTROLS.csv")
    sample = (
        controls.groupby(["parent", "n"], as_index=False)
        .agg(mean_murcko=("n_murcko", "mean"), min_murcko=("n_murcko", "min"), max_murcko=("n_murcko", "max"), n_seeds=("seed", "nunique"))
    )

    coords.to_csv(TAB / "chemical_space_pca_coordinates.csv", index=False)
    tanimoto.to_csv(TAB / "chemical_space_tanimoto_pairs.csv", index=False)
    tanimoto_summary.to_csv(TAB / "chemical_space_tanimoto_summary.csv", index=False)
    scaffolds.to_csv(TAB / "chemical_space_murcko_counts.csv", index=False)
    sample.to_csv(TAB / "scaffold_sample_size_summary.csv", index=False)
    marks[["ligand_id", "name", "group", "max_tanimoto_train"]].to_csv(TAB / "candidate_max_tanimoto.csv", index=False)
    ref.to_csv(TAB / "reference_rmsd_source.csv", index=False)

    facts = {
        "ad_threshold": AD_THRESHOLD,
        "n_clinical_in_domain": int(scores.in_applicability_domain.sum()),
        "n_primary_in_domain": int((marks.group.eq("primary_12") & (marks.max_tanimoto_train >= AD_THRESHOLD)).sum()),
        "n_backup_in_domain": int((marks.group.eq("backup_21") & (marks.max_tanimoto_train >= AD_THRESHOLD)).sum()),
        "n_new_in_domain": int((marks.group.eq("stage3b_13") & (marks.max_tanimoto_train >= AD_THRESHOLD)).sum()),
        "pc_variance": [float(var[0]), float(var[1])],
        "murcko": scaffolds.to_dict(orient="records"),
        "tanimoto": tanimoto_summary.to_dict(orient="records"),
        "sample_size": sample.to_dict(orient="records"),
        "pf_urat1_cnnscore_mean": float(pf.CNNscore.mean()),
        "pf_arg": [float(x) for x in pf.acid_arg477_min_A],
        "pf_centroid": [float(x) for x in pf.centroid_A],
        "stage3b_status_counts": stage.final_status.value_counts().to_dict(),
        "stage3b_names": stage.loc[
            stage.final_status.eq("stage3b_new_structural_and_chemistry_pass"),
            ["repurposing_id", "name", "unified_dual_pass_seeds"],
        ]
        .sort_values("name")
        .to_dict(orient="records"),
        "frozen_overlap_with_stage3b": int(len(set(stage.repurposing_id) & set(groups.ligand_id[groups.group != "stage3b_13"]))),
    }
    (TAB / "locked_facts.json").write_text(json.dumps(facts, indent=2))
    print(json.dumps({k: facts[k] for k in ("pc_variance", "tanimoto", "murcko", "n_clinical_in_domain", "stage3b_status_counts", "pf_urat1_cnnscore_mean")}, indent=2))


if __name__ == "__main__":
    main()
