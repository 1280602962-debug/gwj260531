#!/usr/bin/env python3
"""Analyze MCL1/Bcl-xL Vina scores under frozen θ=6.0 four-state labels."""
from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/mcl1_bclxl_panel_v0"
TAB = OUT / "tables"
ANALYSIS = OUT / "analysis"
ANALYSIS.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260827


def utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def auroc(y_true, scores) -> float:
    """Higher score = more positive. For Vina affinities, pass -affinity."""
    y = np.asarray(y_true, int)
    s = np.asarray(scores, float)
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # Mann-Whitney
    correct = 0.0
    for p in pos:
        correct += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return float(correct / (len(pos) * len(neg)))


def boot_ci(y, s, n=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, int)
    s = np.asarray(s, float)
    vals = []
    n_all = len(y)
    for _ in range(n):
        idx = rng.integers(0, n_all, n_all)
        vals.append(auroc(y[idx], s[idx]))
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    scores = list(csv.DictReader((TAB / "vina_scores_MBX_v1.csv").open()))
    ok = [r for r in scores if r["status"] in {"success", "exists"} and r["vina_mode1"] != ""]
    # pivot ligand -> pockets
    by_lig = {}
    meta = {}
    for r in ok:
        lid = r["panel_id"]
        meta[lid] = {
            "class": r["class"],
            "chembl": r["molecule_chembl_id"],
            "pA": float(r["pchembl_MCL1"]),
            "pB": float(r["pchembl_BCL2L1"]),
        }
        by_lig.setdefault(lid, {})[r["pdb"]] = -float(r["vina_mode1"])  # higher better
    complete = {
        lid: v for lid, v in by_lig.items() if "3WIY" in v and "3WIZ" in v
    }
    classes = Counter(meta[lid]["class"] for lid in complete)
    gate = json.loads((TAB / "lc6_gate_summary_v1.json").read_text())
    panel_role = gate.get("panel_role", "unknown")

    rows_form = []
    # Dual vs neither (mean)
    dual = [lid for lid in complete if meta[lid]["class"] == "dual"]
    neither = [lid for lid in complete if meta[lid]["class"] == "neither"]
    a_only = [lid for lid in complete if meta[lid]["class"] == "A_only"]
    b_only = [lid for lid in complete if meta[lid]["class"] == "B_only"]

    def mean_score(lid):
        return 0.5 * (complete[lid]["3WIY"] + complete[lid]["3WIZ"])

    def worst_score(lid):
        return min(complete[lid]["3WIY"], complete[lid]["3WIZ"])

    contrasts = [
        (
            "conventional_dual_vs_neither",
            "D_vs_neither_mean",
            "score_mean",
            dual,
            neither,
            mean_score,
            "dual vs neither using mean pocket score",
        ),
        (
            "dualfourclass_directional",
            "D_vs_A_pocketB",
            "score_B",
            dual,
            a_only,
            lambda lid: complete[lid]["3WIZ"],
            "dual vs A-only scored in pocket B (Bcl-xL/3WIZ)",
        ),
        (
            "dualfourclass_directional",
            "D_vs_B_pocketA",
            "score_A",
            dual,
            b_only,
            lambda lid: complete[lid]["3WIY"],
            "dual vs B-only scored in pocket A (MCL1/3WIY)",
        ),
        (
            "dualfourclass_directional",
            "summary_min",
            "min(D/A,D/B)",
            dual,
            None,  # special
            None,
            "worst-arm aggregation of the two directional AUROCs",
        ),
    ]

    arm_aurocs = {}
    for form, contrast, score_name, pos, neg, scorer, note in contrasts:
        if contrast == "summary_min":
            vals = [arm_aurocs.get("D_vs_A_pocketB"), arm_aurocs.get("D_vs_B_pocketA")]
            vals = [v for v in vals if v is not None and np.isfinite(v)]
            sm = float(min(vals)) if vals else float("nan")
            rows_form.append(
                {
                    "pair": "MCL1/Bcl-xL",
                    "engine": "vina_mode1",
                    "formulation": form,
                    "contrast": contrast,
                    "score": score_name,
                    "n_pos": len(dual),
                    "n_neg": "",
                    "auroc": round(sm, 4) if np.isfinite(sm) else "",
                    "ci_lo": "",
                    "ci_hi": "",
                    "note": note,
                    "panel_role": panel_role,
                }
            )
            continue
        y = [1] * len(pos) + [0] * len(neg)
        s = [scorer(lid) for lid in pos] + [scorer(lid) for lid in neg]
        auc = auroc(y, s)
        lo, hi = boot_ci(y, s)
        arm_aurocs[contrast] = auc
        rows_form.append(
            {
                "pair": "MCL1/Bcl-xL",
                "engine": "vina_mode1",
                "formulation": form,
                "contrast": contrast,
                "score": score_name,
                "n_pos": len(pos),
                "n_neg": len(neg),
                "auroc": round(auc, 4) if np.isfinite(auc) else "",
                "ci_lo": round(lo, 4) if np.isfinite(lo) else "",
                "ci_hi": round(hi, 4) if np.isfinite(hi) else "",
                "note": note,
                "panel_role": panel_role,
            }
        )

    # write formulation table
    keys = list(rows_form[0].keys())
    with (TAB / "formulation_auroc_MBX_v1.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows_form)

    # per-ligand wide table
    wide = []
    for lid, pockets in sorted(complete.items()):
        wide.append(
            {
                "panel_id": lid,
                "molecule_chembl_id": meta[lid]["chembl"],
                "class": meta[lid]["class"],
                "vina_3WIY": -pockets["3WIY"],
                "vina_3WIZ": -pockets["3WIZ"],
                "score_mean": mean_score(lid),
                "score_worst": worst_score(lid),
            }
        )
    with (TAB / "ligand_scores_wide_MBX_v1.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(wide[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(wide)

    sm = next(r for r in rows_form if r["contrast"] == "summary_min")
    dnn = next(r for r in rows_form if r["contrast"] == "D_vs_neither_mean")
    md = f"""# MCL1_BCLXL_DOCKING_VERDICT_V1

Updated: `{utc()}`

## Scope

Frozen ChEMBL θ=6.0 panel96 (24/24/24/24), receptors **3WIY / 3WIZ**, Vina mode-1.
LC6 pose-gold gate role: **`{panel_role}`**.

## Completeness

| class | complete ligands (both pockets) |
|-------|--------------------------------:|
| dual | {classes.get('dual', 0)} |
| A_only | {classes.get('A_only', 0)} |
| B_only | {classes.get('B_only', 0)} |
| neither | {classes.get('neither', 0)} |
| **total** | **{len(complete)}** / 96 |

Jobs scored ok: {len(ok)} / {len(scores)}.

## Formulation AUROCs (descriptive)

| contrast | n_pos | n_neg | AUROC | 95% CI |
|----------|------:|------:|------:|--------|
| Dual vs neither (mean) | {dnn['n_pos']} | {dnn['n_neg']} | {dnn['auroc']} | [{dnn['ci_lo']}, {dnn['ci_hi']}] |
| Dual vs A-only @ Bcl-xL | {arm_aurocs and rows_form[1]['n_pos']} | {rows_form[1]['n_neg']} | {rows_form[1]['auroc']} | [{rows_form[1]['ci_lo']}, {rows_form[1]['ci_hi']}] |
| Dual vs B-only @ MCL1 | {rows_form[2]['n_pos']} | {rows_form[2]['n_neg']} | {rows_form[2]['auroc']} | [{rows_form[2]['ci_lo']}, {rows_form[2]['ci_hi']}] |
| **summary_min** | {sm['n_pos']} | — | **{sm['auroc']}** | — |

## Interpretation rules

- If `panel_role=applicability_stress_test`: report honestly; do **not** claim target-general PPI screening performance.
- Homologous BCL-2 fold / BH3-groove domain shift — not a first non-kinase pair (AChE/BChE already is).
- Do not package as external validation.

## Files

- `tables/vina_scores_MBX_v1.csv`
- `tables/formulation_auroc_MBX_v1.csv`
- `tables/ligand_scores_wide_MBX_v1.csv`
- `analysis/MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1.md`
"""
    (ANALYSIS / "MCL1_BCLXL_DOCKING_VERDICT_V1.md").write_text(md)
    print(md)
    print(json.dumps(rows_form, indent=2))


if __name__ == "__main__":
    main()
