#!/usr/bin/env python3
"""AUROC endpoints for GNINA independent-docking mode-1 scores (θ=6.0 labels frozen)."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_independent_dock_v0"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 20260729
THETA = 6.0

SPEC = {
    "EGFR/HER2": {
        "scores": OUT / "tables/gnina_dock_scores_EGFR_HER2.csv",
        "panel": ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        "target_a": "3POZ",
        "target_b": "3RCD",
    },
    "PIK3CA/mTOR": {
        "scores": OUT / "tables/gnina_dock_scores_PIK3CA_mTOR.csv",
        "panel": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        "target_a": "4L23",
        "target_b": "4JT6",
    },
}

UNIFIED = ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
FORMULATION = ROOT / "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv"


def vina_refs(pair: str) -> tuple[float, float | None]:
    """Vina comparison numbers from committed Table 2 / Table 3 sources.

    Dual-versus-neither is filled only for EGFR/HER2. PIK3CA/mTOR neither n = 4
    is underpowered and is left blank, matching the deposited summary contract.
    """
    summary_min = None
    with UNIFIED.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["pair"] == pair and row["label_rule"] == "theta_6.0":
                summary_min = float(row["pocket_matched_summary_min"])
                break
    neither = None
    if pair == "EGFR/HER2":
        with FORMULATION.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if (
                    row["pair"] == pair
                    and row["formulation"] == "conventional_dual_vs_neither"
                    and row["contrast"] == "D_vs_neither_mean"
                ):
                    neither = float(row["auroc"])
                    break
    if summary_min is None:
        raise KeyError(f"no Table 2 summary_min for {pair}")
    return summary_min, neither


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


def add_contrast(rows, pair, formulation, contrast, score, pos, neg, n_pos, n_neg, note):
    h = hashlib.md5(f"{pair}|{formulation}|{contrast}".encode()).hexdigest()
    pt, lo, hi = boot_auroc(pos, neg, seed=SEED + (int(h[:8], 16) % 99991))
    rows.append(
        {
            "pair": pair,
            "engine": "gnina_dock_mode1",
            "formulation": formulation,
            "contrast": contrast,
            "score": score,
            "n_pos": n_pos,
            "n_neg": n_neg,
            "auroc": "" if pt != pt else round(pt, 4),
            "ci_lo": "" if lo != lo else round(lo, 4),
            "ci_hi": "" if hi != hi else round(hi, 4),
            "note": note,
        }
    )


def assemble_gnina(pair: str, cfg: dict) -> list[dict]:
    panel = pd.read_csv(cfg["panel"])
    cls_map = {}
    for _, r in panel.iterrows():
        lid = r.get("panel_id") or r.get("ligand")
        cls_map[lid] = r.get("class") or r.get("cls")

    scores = pd.read_csv(cfg["scores"])
    ok = scores[scores["status"].isin(["success", "exists"])].copy()
    pivot = ok.pivot_table(index="ligand", columns="target", values="gnina_mode1", aggfunc="first")
    recs = []
    for lig, row in pivot.iterrows():
        a = row.get(cfg["target_a"])
        b = row.get(cfg["target_b"])
        if pd.isna(a) or pd.isna(b):
            continue
        cls = cls_map.get(lig)
        if not cls:
            continue
        recs.append(
            {
                "ligand": lig,
                "cls": cls,
                "score_A": -float(a),
                "score_B": -float(b),
                "score_mean": -(float(a) + float(b)) / 2.0,
                "score_worst": min(-float(a), -float(b)),
            }
        )
    return recs


def formulation_table(pair: str, recs: list[dict]) -> list[dict]:
    rows = []
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    N = [r for r in recs if r["cls"] == "neither"]

    add_contrast(
        rows, pair, "dualfourclass_directional", "D_vs_A_pocketB", "score_B",
        [r["score_B"] for r in D], [r["score_B"] for r in A], len(D), len(A),
        "dual vs A-only scored in pocket B",
    )
    add_contrast(
        rows, pair, "dualfourclass_directional", "D_vs_B_pocketA", "score_A",
        [r["score_A"] for r in D], [r["score_A"] for r in B], len(D), len(B),
        "dual vs B-only scored in pocket A",
    )
    da = auroc([r["score_B"] for r in D], [r["score_B"] for r in A])
    db = auroc([r["score_A"] for r in D], [r["score_A"] for r in B])
    rows.append(
        {
            "pair": pair,
            "engine": "gnina_dock_mode1",
            "formulation": "dualfourclass_directional",
            "contrast": "summary_min",
            "score": "min(D/A,D/B)",
            "n_pos": len(D),
            "n_neg": min(len(A), len(B)),
            "auroc": round(min(da, db), 4),
            "ci_lo": "",
            "ci_hi": "",
            "note": "worst-arm aggregation",
        }
    )
    add_contrast(
        rows, pair, "conventional_dual_vs_neither", "D_vs_neither_mean", "score_mean",
        [r["score_mean"] for r in D], [r["score_mean"] for r in N], len(D), len(N),
        "dual vs neither using mean pocket score",
    )
    return rows


def enrichment_top10(pair: str, recs: list[dict]) -> list[dict]:
    rows = []
    n = len(recs)
    n_dual = sum(r["cls"] == "dual" for r in recs)
    if n_dual == 0:
        return rows
    random_hit = n_dual / n
    for score_name in ("score_mean", "score_worst"):
        ranked = sorted(recs, key=lambda r: r[score_name], reverse=True)
        k = min(10, n)
        top = ranked[:k]
        n_hit = sum(r["cls"] == "dual" for r in top)
        counts = Counter(r["cls"] for r in top)
        ef = (n_hit / k) / random_hit if random_hit > 0 else float("nan")
        rows.append(
            {
                "pair": pair,
                "engine": "gnina_dock_mode1",
                "score": score_name,
                "k": k,
                "n_dual_top": n_hit,
                "n_A_only_top": counts.get("A_only", 0),
                "n_B_only_top": counts.get("B_only", 0),
                "n_neither_top": counts.get("neither", 0),
                "EF": round(ef, 3),
            }
        )
    return rows


def verdict(gnina_sm: float, gnina_neither: float | None, vina_sm: float, vina_neither: float | None) -> str:
    """Label whether the Dual-versus-neither vs directional formulation gap remains.

    Compares formulation gaps (neither − summary_min), not GNINA vs Vina
    directional AUROC. EGFR/HER2 keeps `gap_remains` when both engines still
    show a large neither-minus-directional gap.
    """
    if vina_neither is not None and gnina_neither is not None:
        g_gap = gnina_neither - gnina_sm
        v_gap = vina_neither - vina_sm
        if g_gap >= 0.15 and v_gap >= 0.15:
            return "gap_remains"
        if g_gap > 0 and v_gap > 0:
            return "gap_smaller_same_sign"
        return "gap_gone_or_reversed"
    if abs(gnina_sm - vina_sm) <= 0.08:
        return "gap_smaller_same_sign"
    return "gap_remains"


def main() -> None:
    all_form = []
    all_enr = []
    summary_rows = []

    for pair, cfg in SPEC.items():
        if not cfg["scores"].exists():
            print(f"SKIP {pair}: missing {cfg['scores']}", flush=True)
            continue
        recs = assemble_gnina(pair, cfg)
        print(pair, Counter(r["cls"] for r in recs), "n=", len(recs), flush=True)
        form = formulation_table(pair, recs)
        enr = enrichment_top10(pair, recs)
        all_form.extend(form)
        all_enr.extend(enr)

        sm_row = next(r for r in form if r["contrast"] == "summary_min")
        dn_row = next((r for r in form if r["contrast"] == "D_vs_neither_mean"), None)
        sm = float(sm_row["auroc"])
        dn = float(dn_row["auroc"]) if dn_row and dn_row["auroc"] != "" else None
        vina_sm, vina_neither = vina_refs(pair)
        summary_rows.append(
            {
                "pair": pair,
                "engine": "gnina_dock_mode1",
                "gnina_summary_min": sm,
                "gnina_D_vs_neither_mean": dn if dn is not None else "",
                "vina_summary_min_ref": vina_sm,
                "vina_D_vs_neither_ref": vina_neither if vina_neither is not None else "",
                "delta_summary_min_vs_vina": round(sm - vina_sm, 4),
                "verdict": verdict(sm, dn, vina_sm, vina_neither),
            }
        )

    pd.DataFrame(all_form).to_csv(TAB / "independent_dock_formulation_v1.csv", index=False)
    pd.DataFrame(all_enr).to_csv(TAB / "independent_dock_enrichment_v1.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(TAB / "independent_dock_summary_v1.csv", index=False)

    meta = {
        "engine": "GNINA 1.3.2 docking search (mode-1 minimizedAffinity)",
        "n_boot": N_BOOT,
        "seed": SEED,
        "theta": THETA,
        "summary": summary_rows,
    }
    (OUT / "analysis").mkdir(exist_ok=True)
    (OUT / "analysis/independent_dock_run_meta_v1.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(summary_rows, indent=2), flush=True)


if __name__ == "__main__":
    main()
