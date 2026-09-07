#!/usr/bin/env python3
"""Track B pocket-matched directional AUROC (DOCKING_PLAN_V1).

Primary: Dual vs A_only in pocket B; Dual vs B_only in pocket A;
summary_min = min(AUROC_D/A, AUROC_D/B). Ligands need both-end scores.
Also report Dual vs neither with vina_mean (formulation contrast).
Bootstrap: B=2000, seed 20260729. Does not replace Table 2.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
TABLES = ROOT / "tables"
OUT_TAB = LOCAL / "tables"
OUT_AN = LOCAL / "analysis"
N_BOOT = 2000
SEED = 20260729

PAIRS = [
    {
        "pair": "F2/F10",
        "system": "coagulation",
        "panel": TABLES / "track_b_panels" / "panel_F2_F10_v1.csv",
        "target_a": "4UDW",
        "target_b": "2JKH",
    },
    {
        "pair": "JAK1/TYK2",
        "system": "JAK",
        "panel": TABLES / "track_b_panels" / "panel_JAK1_TYK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "3LXP",
    },
    {
        "pair": "JAK1/JAK2",
        "system": "JAK",
        "panel": TABLES / "track_b_panels" / "panel_JAK1_JAK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "8BXH",
    },
    {
        "pair": "PPARG/PPARA",
        "system": "PPAR",
        "panel": TABLES / "track_b_panels" / "panel_PPARG_PPARA_v1.csv",
        "target_a": "9V8H",
        "target_b": "6LXA",
    },
    {
        "pair": "PPARA/PPARD",
        "system": "PPAR",
        "panel": TABLES / "track_b_panels" / "panel_PPARA_PPARD_v1.csv",
        "target_a": "6LXA",
        "target_b": "5U3Q",
    },
]


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def boot_auroc(pos, neg, n_boot=N_BOOT, seed=SEED):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        p = rng.choice(pos, size=len(pos), replace=True)
        n = rng.choice(neg, size=len(neg), replace=True)
        vals.append(auroc(p, n))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(auroc(pos, neg)), float(lo), float(hi)


def boot_summary_min(score_b_d, score_b_a, score_a_d, score_a_b, n_boot=N_BOOT, seed=SEED):
    """Resample each class arm and recompute min(D/A, D/B) each replicate."""
    sb_d = np.asarray(score_b_d, float)
    sb_a = np.asarray(score_b_a, float)
    sa_d = np.asarray(score_a_d, float)
    sa_b = np.asarray(score_a_b, float)
    pt = min(auroc(sb_d, sb_a), auroc(sa_d, sa_b))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        da = auroc(rng.choice(sb_d, size=len(sb_d), replace=True), rng.choice(sb_a, size=len(sb_a), replace=True))
        db = auroc(rng.choice(sa_d, size=len(sa_d), replace=True), rng.choice(sa_b, size=len(sa_b), replace=True))
        vals.append(min(da, db))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(pt), float(lo), float(hi)


def load_scores():
    rows = list(csv.DictReader((OUT_TAB / "scores_vina_mode1_v1.csv").open()))
    # pair -> ligand -> target -> score_S
    out = {}
    for r in rows:
        out.setdefault(r["pair"], {}).setdefault(r["ligand"], {})[r["target"]] = float(r["score_S"])
    return out


def assemble(spec, scores_by_pair, label_col: str):
    panel = list(csv.DictReader(spec["panel"].open()))
    sc = scores_by_pair.get(spec["pair"], {})
    recs = []
    missing = Counter()
    for row in panel:
        lig = row["panel_id"]
        cls = row[label_col]
        m = sc.get(lig, {})
        sa = m.get(spec["target_a"])
        sb = m.get(spec["target_b"])
        if sa is None or sb is None:
            missing[cls] += 1
            continue
        recs.append(
            {
                "ligand": lig,
                "cls": cls,
                "score_A": sa,
                "score_B": sb,
                "score_mean": (sa + sb) / 2.0,
            }
        )
    return recs, missing


def analyze_pair(spec, recs, label_rule: str):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    N = [r for r in recs if r["cls"] == "neither"]
    rows = []

    def add(contrast, score_name, pos, neg, note, seed_salt):
        h = int(hashlib.md5(f"{spec['pair']}|{label_rule}|{contrast}".encode()).hexdigest()[:8], 16)
        pt, lo, hi = boot_auroc(pos, neg, seed=SEED + (h % 99991) + seed_salt)
        rows.append(
            {
                "pair": spec["pair"],
                "system": spec["system"],
                "label_rule": label_rule,
                "engine": "vina_mode1",
                "formulation": "pocket_matched_directional",
                "contrast": contrast,
                "score": score_name,
                "n_pos": len(pos),
                "n_neg": len(neg),
                "n_scored_both_ends": len(recs),
                "auroc": "" if pt != pt else round(pt, 4),
                "ci_lo": "" if lo != lo else round(lo, 4),
                "ci_hi": "" if hi != hi else round(hi, 4),
                "note": note,
            }
        )
        return pt

    da = add(
        "D_vs_A_pocketB",
        "score_B",
        [r["score_B"] for r in D],
        [r["score_B"] for r in A],
        "dual vs A-only scored in pocket B",
        0,
    )
    db = add(
        "D_vs_B_pocketA",
        "score_A",
        [r["score_A"] for r in D],
        [r["score_A"] for r in B],
        "dual vs B-only scored in pocket A",
        1,
    )
    h = int(hashlib.md5(f"{spec['pair']}|{label_rule}|summary_min".encode()).hexdigest()[:8], 16)
    spt, slo, shi = boot_summary_min(
        [r["score_B"] for r in D],
        [r["score_B"] for r in A],
        [r["score_A"] for r in D],
        [r["score_A"] for r in B],
        seed=SEED + (h % 99991) + 2,
    )
    rows.append(
        {
            "pair": spec["pair"],
            "system": spec["system"],
            "label_rule": label_rule,
            "engine": "vina_mode1",
            "formulation": "pocket_matched_directional",
            "contrast": "summary_min",
            "score": "min(D/A,D/B)",
            "n_pos": len(D),
            "n_neg": min(len(A), len(B)),
            "n_scored_both_ends": len(recs),
            "auroc": round(spt, 4),
            "ci_lo": round(slo, 4),
            "ci_hi": round(shi, 4),
            "note": "worst-direction summary; CI from class-preserving bootstrap of min(DA,DB)",
        }
    )
    add(
        "D_vs_neither_vina_mean",
        "score_mean",
        [r["score_mean"] for r in D],
        [r["score_mean"] for r in N],
        "formulation contrast only; not primary",
        3,
    )
    # flip formulation tag for neither row
    rows[-1]["formulation"] = "dual_vs_neither"

    summary = {
        "pair": spec["pair"],
        "system": spec["system"],
        "label_rule": label_rule,
        "n_dual": len(D),
        "n_A_only": len(A),
        "n_B_only": len(B),
        "n_neither": len(N),
        "n_scored_both_ends": len(recs),
        "auroc_D_vs_A_pocketB": round(da, 4) if da == da else "",
        "auroc_D_vs_B_pocketA": round(db, 4) if db == db else "",
        "summary_min": round(spt, 4),
        "summary_min_ci_lo": round(slo, 4),
        "summary_min_ci_hi": round(shi, 4),
        "target_a": spec["target_a"],
        "target_b": spec["target_b"],
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    OUT_AN.mkdir(parents=True, exist_ok=True)
    scores = load_scores()
    all_rows = []
    summaries = []
    for label_col, label_rule in (("class", "strict_6.5_5.5_panel"), ("theta6_class", "theta6_0")):
        for spec in PAIRS:
            recs, missing = assemble(spec, scores, label_col)
            print(
                f"{spec['pair']} {label_rule}: scored={len(recs)} class={dict(Counter(r['cls'] for r in recs))} missing={dict(missing)}",
                flush=True,
            )
            rows, summary = analyze_pair(spec, recs, label_rule)
            all_rows.extend(rows)
            summaries.append(summary)

    write_csv(
        OUT_TAB / "track_b_directional_auroc_v1.csv",
        all_rows,
        [
            "pair",
            "system",
            "label_rule",
            "engine",
            "formulation",
            "contrast",
            "score",
            "n_pos",
            "n_neg",
            "n_scored_both_ends",
            "auroc",
            "ci_lo",
            "ci_hi",
            "note",
        ],
    )
    write_csv(
        OUT_TAB / "track_b_summary_min_v1.csv",
        summaries,
        [
            "pair",
            "system",
            "label_rule",
            "n_dual",
            "n_A_only",
            "n_B_only",
            "n_neither",
            "n_scored_both_ends",
            "auroc_D_vs_A_pocketB",
            "auroc_D_vs_B_pocketA",
            "summary_min",
            "summary_min_ci_lo",
            "summary_min_ci_hi",
            "target_a",
            "target_b",
        ],
    )

    # Primary readout table: strict panel labels (how panels were built)
    primary = [s for s in summaries if s["label_rule"] == "strict_6.5_5.5_panel"]
    lines = [
        "# Track B directional AUROC (local Vina pack)\n\n",
        "Engine: AutoDock Vina 1.2.7 mode-1; `S = −E`. Seed 20260727. Exhaustiveness 8.\n",
        "Estimand: pocket-matched directional AUROC; `summary_min = min(D/A, D/B)`.\n",
        "Bootstrap: B=2000, seed 20260729, class-preserving. **Does not replace Table 2.**\n",
        "Count **three systems** (coagulation, JAK, PPAR), not five pairs.\n\n",
        "## Primary (panel strict 6.5/5.5 class)\n\n",
        "| pair | system | n_scored | D/A (pocket B) | D/B (pocket A) | summary_min [95% CI] |\n",
        "|------|--------|---------:|---------------:|---------------:|---------------------:|\n",
    ]
    for s in primary:
        lines.append(
            f"| {s['pair']} | {s['system']} | {s['n_scored_both_ends']} | "
            f"{s['auroc_D_vs_A_pocketB']} | {s['auroc_D_vs_B_pocketA']} | "
            f"{s['summary_min']} [{s['summary_min_ci_lo']}, {s['summary_min_ci_hi']}] |\n"
        )
    lines.append("\n## Unified θ = 6.0 labels (same scores)\n\n")
    lines.append(
        "| pair | system | n_scored | D/A (pocket B) | D/B (pocket A) | summary_min [95% CI] |\n"
        "|------|--------|---------:|---------------:|---------------:|---------------------:|\n"
    )
    for s in summaries:
        if s["label_rule"] != "theta6_0":
            continue
        lines.append(
            f"| {s['pair']} | {s['system']} | {s['n_scored_both_ends']} | "
            f"{s['auroc_D_vs_A_pocketB']} | {s['auroc_D_vs_B_pocketA']} | "
            f"{s['summary_min']} [{s['summary_min_ci_lo']}, {s['summary_min_ci_hi']}] |\n"
        )
    lines.append(
        "\nArtifacts: `tables/track_b_summary_min_v1.csv`, `tables/track_b_directional_auroc_v1.csv`, "
        "`tables/scores_vina_mode1_v1.csv`, `tables/job_status.csv`, `tables/layer3_cognate_rmsd_v1.csv`.\n"
    )
    (OUT_AN / "TRACK_B_DIRECTIONAL_AUROC_V1.md").write_text("".join(lines))
    (OUT_AN / "track_b_summary_min_v1.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print("wrote", OUT_TAB / "track_b_summary_min_v1.csv")
    print("wrote", OUT_AN / "TRACK_B_DIRECTIONAL_AUROC_V1.md")
    for s in primary:
        print(
            f"PRIMARY {s['pair']}: summary_min={s['summary_min']} "
            f"[{s['summary_min_ci_lo']}, {s['summary_min_ci_hi']}] "
            f"DA={s['auroc_D_vs_A_pocketB']} DB={s['auroc_D_vs_B_pocketA']}"
        )


if __name__ == "__main__":
    main()
