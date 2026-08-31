#!/usr/bin/env python3
"""Holdout wrong-pocket diagnostic: potency/size matching + vs main-panel pChEMBL.

Zero new docking. Uses frozen holdout panels (pChEMBL) and vina_A/vina_B from
`holdout_ligand_scores_v1.csv`. Matching rule copied from
`jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py`:
  potency: |Δp| ≤ 0.5 on the shared-active end (D vs A_only: pA; D vs B_only: pB)
  size:    |Δheavy| ≤ 2

Question: is holdout wrong_pocket ≥ pocket_matched an artifact of unused-pool
sampling (potency/size shift vs the main panel), or does it survive matching?
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HOLDOUT = ROOT / "data" / "jcim_holdout_v0"
BENCH = ROOT / "data" / "jcim_bench_v0" / "tables"
TAB = HOLDOUT / "tables"
AN = HOLDOUT / "analysis"

N_BOOT = 2000
SEED = 20260729
DP = 0.5
DS = 2.0

PAIRS = {
    "HOAB": {
        "pair": "AChE/BChE",
        "pA": "pchembl_ACHE",
        "pB": "pchembl_BCHE",
        "assembled": "assembled_AChE_BChE.csv",
    },
    "HOAP": {
        "pair": "PIK3CA/PIK3CB",
        "pA": "pchembl_PIK3CA",
        "pB": "pchembl_PIK3CB",
        "assembled": "assembled_PIK3CA_PIK3CB.csv",
    },
    "HOPM": {
        "pair": "PIK3CA/mTOR",
        "pA": "pchembl_PIK3CA",
        "pB": "pchembl_MTOR",
        "assembled": "assembled_PIK3CA_mTOR.csv",
    },
}


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def summarize(xs: list[float]) -> dict:
    xs = [x for x in xs if x is not None and x == x]
    if not xs:
        return {"n": 0, "mean": None, "median": None, "iqr": None, "min": None, "max": None}
    s = sorted(xs)
    n = len(s)
    q1 = s[max(0, n // 4)]
    q3 = s[min(n - 1, (3 * n) // 4)]
    return {
        "n": n,
        "mean": round(sum(s) / n, 3),
        "median": round(statistics.median(s), 3),
        "iqr": round(q3 - q1, 3),
        "min": round(s[0], 3),
        "max": round(s[-1], 3),
    }


def nearest_match(duals, others, key_potency, dp, ds):
    kept_d, kept_o = [], []
    used = set()
    for d in duals:
        best, best_dist = None, 1e9
        for i, o in enumerate(others):
            if i in used:
                continue
            if d.get(key_potency) is None or o.get(key_potency) is None:
                continue
            dpot = abs(d[key_potency] - o[key_potency])
            dsz = abs(d["heavy"] - o["heavy"])
            if dpot <= dp and dsz <= ds:
                dist = dpot + 0.1 * dsz
                if dist < best_dist:
                    best_dist, best = dist, i
        if best is not None:
            used.add(best)
            kept_d.append(d)
            kept_o.append(others[best])
    return kept_d, kept_o


def boot_single(pos_scores, neg_scores, seed):
    if len(pos_scores) == 0 or len(neg_scores) == 0:
        return None, None
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(N_BOOT):
        p = rng.choice(pos_scores, size=len(pos_scores), replace=True)
        n = rng.choice(neg_scores, size=len(neg_scores), replace=True)
        v = auroc(p, n)
        if v == v:
            vals.append(v)
    if len(vals) < N_BOOT // 2:
        return None, None
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def load_holdout(prefix: str) -> list[dict]:
    meta = PAIRS[prefix]
    panel = {
        r["holdout_id"]: r
        for r in csv.DictReader((TAB / f"holdout_panel_{prefix}.csv").open())
    }
    scores = list(csv.DictReader((TAB / "holdout_ligand_scores_v1.csv").open()))
    recs = []
    for s in scores:
        if s["prefix"] != prefix:
            continue
        p = panel.get(s["ligand"])
        if p is None:
            continue
        recs.append(
            {
                "ligand": s["ligand"],
                "cls": s["cls"],
                "pA": fnum(p[meta["pA"]]),
                "pB": fnum(p[meta["pB"]]),
                "vina_A": fnum(s["vina_A"]),
                "vina_B": fnum(s["vina_B"]),
                "heavy": fnum(s["heavy"]),
                "source": "holdout",
            }
        )
    return [r for r in recs if r["vina_A"] is not None and r["vina_B"] is not None]


def load_main(prefix: str) -> list[dict]:
    meta = PAIRS[prefix]
    path = BENCH / meta["assembled"]
    recs = []
    for r in csv.DictReader(path.open()):
        recs.append(
            {
                "ligand": r["ligand"],
                "cls": r["cls"],
                "pA": fnum(r["pA"]),
                "pB": fnum(r["pB"]),
                "heavy": fnum(r["heavy"]),
                "source": "main_panel",
            }
        )
    return recs


def dist_rows(prefix: str, hold: list[dict], main: list[dict]) -> list[dict]:
    rows = []
    pair = PAIRS[prefix]["pair"]
    for cls in ("dual", "A_only", "B_only"):
        for key, lab in (("pA", "pA"), ("pB", "pB"), ("heavy", "heavy")):
            hs = [r[key] for r in hold if r["cls"] == cls]
            ms = [r[key] for r in main if r["cls"] == cls]
            sh, sm = summarize(hs), summarize(ms)
            delta = None
            if sh["mean"] is not None and sm["mean"] is not None:
                delta = round(sh["mean"] - sm["mean"], 3)
            rows.append(
                {
                    "prefix": prefix,
                    "pair": pair,
                    "cls": cls,
                    "feature": lab,
                    "holdout_n": sh["n"],
                    "holdout_mean": sh["mean"],
                    "holdout_median": sh["median"],
                    "holdout_iqr": sh["iqr"],
                    "holdout_min": sh["min"],
                    "holdout_max": sh["max"],
                    "main_n": sm["n"],
                    "main_mean": sm["mean"],
                    "main_median": sm["median"],
                    "main_iqr": sm["iqr"],
                    "mean_delta_holdout_minus_main": delta,
                }
            )
    return rows


def contrast_row(prefix, pair, match_type, aggregation, kd, ko, score_key, seed):
    pos = [r[score_key] for r in kd]
    neg = [r[score_key] for r in ko]
    point = auroc(pos, neg)
    lo, hi = boot_single(pos, neg, seed)
    n = min(len(kd), len(ko))
    return {
        "prefix": prefix,
        "pair": pair,
        "match_type": match_type,
        "aggregation": aggregation,
        "n_dual_matched": len(kd),
        "n_other_matched": len(ko),
        "n_contrast": n,
        "underpowered": int(n < 8),
        "auroc": round(point, 4) if point == point else "",
        "ci_lo": round(lo, 4) if lo is not None else "",
        "ci_hi": round(hi, 4) if hi is not None else "",
        "score_key": score_key,
    }


def matched_rows(prefix: str, recs: list[dict]) -> list[dict]:
    pair = PAIRS[prefix]["pair"]
    duals = [r for r in recs if r["cls"] == "dual"]
    aonly = [r for r in recs if r["cls"] == "A_only"]
    bonly = [r for r in recs if r["cls"] == "B_only"]
    specs = [
        ("potency_matched_D_vs_A", duals, aonly, "pA", DP, 999.0, "vina_B", "vina_A"),
        ("potency_matched_D_vs_B", duals, bonly, "pB", DP, 999.0, "vina_A", "vina_B"),
        ("size_matched_D_vs_A", duals, aonly, "pA", 999.0, DS, "vina_B", "vina_A"),
        ("size_matched_D_vs_B", duals, bonly, "pB", 999.0, DS, "vina_A", "vina_B"),
        ("unmatched_D_vs_A", duals, aonly, None, None, None, "vina_B", "vina_A"),
        ("unmatched_D_vs_B", duals, bonly, None, None, None, "vina_A", "vina_B"),
    ]
    out = []
    for i, (name, d, o, key, dp, ds, k_match, k_wrong) in enumerate(specs):
        if key is None:
            kd, ko = d, o
        else:
            kd, ko = nearest_match(d, o, key, dp, ds)
        out.append(
            contrast_row(
                prefix, pair, name, "pocket_matched", kd, ko, k_match, SEED + i * 17
            )
        )
        out.append(
            contrast_row(
                prefix,
                pair,
                name,
                "wrong_pocket",
                kd,
                ko,
                k_wrong,
                SEED + 1000 + i * 17,
            )
        )
    return out


def summary_min_rows(contrast_rows: list[dict]) -> list[dict]:
    """Pair the two arms of the same match family into summary_min."""
    families = {
        "potency_matched": ("potency_matched_D_vs_A", "potency_matched_D_vs_B"),
        "size_matched": ("size_matched_D_vs_A", "size_matched_D_vs_B"),
        "unmatched": ("unmatched_D_vs_A", "unmatched_D_vs_B"),
    }
    by = {(r["prefix"], r["match_type"], r["aggregation"]): r for r in contrast_rows}
    out = []
    for prefix, meta in PAIRS.items():
        for fam, (a_name, b_name) in families.items():
            for agg in ("pocket_matched", "wrong_pocket"):
                ra = by.get((prefix, a_name, agg))
                rb = by.get((prefix, b_name, agg))
                if ra is None or rb is None:
                    continue
                da = ra["auroc"]
                db = rb["auroc"]
                sm = ""
                if da != "" and db != "":
                    sm = round(min(float(da), float(db)), 4)
                nmin = min(ra["n_contrast"], rb["n_contrast"])
                out.append(
                    {
                        "prefix": prefix,
                        "pair": meta["pair"],
                        "family": fam,
                        "aggregation": agg,
                        "auroc_D_vs_A": da,
                        "auroc_D_vs_B": db,
                        "summary_min": sm,
                        "n_D_vs_A": ra["n_contrast"],
                        "n_D_vs_B": rb["n_contrast"],
                        "n_min": nmin,
                        "underpowered": int(nmin < 8),
                        "gap_matched_minus_wrong": "",
                    }
                )
        # fill gaps
        keyed = {(r["family"], r["aggregation"]): r for r in out if r["prefix"] == prefix}
        for fam in families:
            pm = keyed.get((fam, "pocket_matched"))
            wp = keyed.get((fam, "wrong_pocket"))
            if pm and wp and pm["summary_min"] != "" and wp["summary_min"] != "":
                gap = round(float(pm["summary_min"]) - float(wp["summary_min"]), 4)
                pm["gap_matched_minus_wrong"] = gap
                wp["gap_matched_minus_wrong"] = gap
    return out


def verdict_md(dist, contrast, summary) -> str:
    lines = [
        "# HOLDOUT_WRONG_POCKET_POTENCY_VERDICT_V1",
        "",
        "Zero new docking. Matching copied from T0 (`|Δp|≤0.5` potency; `|Δheavy|≤2` size).",
        "Pocket-matched: D/A uses vina_B; D/B uses vina_A. Wrong-pocket: the reverse.",
        "",
        "## Does unused-pool sampling shift potency vs the main panel?",
        "",
        "| pair | class | feature | holdout mean (n) | main mean (n) | Δ |",
        "|------|-------|---------|-----------------:|--------------:|--:|",
    ]
    for r in dist:
        if r["feature"] == "heavy" and r["cls"] == "neither":
            continue
        lines.append(
            f"| {r['pair']} | {r['cls']} | {r['feature']} | "
            f"{r['holdout_mean']} ({r['holdout_n']}) | "
            f"{r['main_mean']} ({r['main_n']}) | "
            f"{r['mean_delta_holdout_minus_main']} |"
        )
    lines += [
        "",
        "## Matched-subset directional AUROC (holdout)",
        "",
        "| pair | family | aggregation | D/A | D/B | summary_min | n_min | underpowered |",
        "|------|--------|-------------|----:|----:|------------:|------:|:------------:|",
    ]
    for r in summary:
        lines.append(
            f"| {r['pair']} | {r['family']} | {r['aggregation']} | "
            f"{r['auroc_D_vs_A']} | {r['auroc_D_vs_B']} | {r['summary_min']} | "
            f"{r['n_min']} | {'Y' if r['underpowered'] else 'N'} |"
        )
    lines += ["", "## Does matching flip wrong-pocket ≥ matched?", ""]
    flips = []
    for prefix, meta in PAIRS.items():
        rows = [r for r in summary if r["prefix"] == prefix]
        by = {(r["family"], r["aggregation"]): r for r in rows}
        for fam in ("unmatched", "potency_matched", "size_matched"):
            pm = by.get((fam, "pocket_matched"))
            wp = by.get((fam, "wrong_pocket"))
            if not pm or not wp or pm["summary_min"] == "" or wp["summary_min"] == "":
                continue
            gap = float(pm["summary_min"]) - float(wp["summary_min"])
            wp_ge = float(wp["summary_min"]) >= float(pm["summary_min"]) - 1e-9
            flips.append(
                f"- **{meta['pair']} / {fam}**: pocket-matched {pm['summary_min']} vs "
                f"wrong-pocket {wp['summary_min']} (gap matched−wrong = {gap:.3f}); "
                f"wrong ≥ matched: **{'yes' if wp_ge else 'no'}**"
                f"{'; underpowered (n_min<8)' if pm['underpowered'] or wp['underpowered'] else ''}."
            )
    lines.extend(flips)
    # headline
    pot_still = []
    pot_flipped = []
    for prefix, meta in PAIRS.items():
        rows = [r for r in summary if r["prefix"] == prefix]
        by = {(r["family"], r["aggregation"]): r for r in rows}
        pm = by.get(("potency_matched", "pocket_matched"))
        wp = by.get(("potency_matched", "wrong_pocket"))
        if not pm or not wp or pm["summary_min"] == "" or wp["summary_min"] == "":
            continue
        if float(wp["summary_min"]) >= float(pm["summary_min"]) - 1e-9:
            pot_still.append(meta["pair"])
        else:
            pot_flipped.append(meta["pair"])
    lines += ["", "### One-line verdict", ""]
    if pot_still and not pot_flipped:
        lines.append(
            "**Holdout wrong-pocket ≥ pocket-matched survives potency matching on all three pairs.** "
            "Unused-pool sampling is not a sufficient explanation of the paradox. "
            "This is a diagnostic, not a protocol change: primary holdout numbers stay unmatched."
        )
    elif pot_flipped and not pot_still:
        lines.append(
            "**Potency matching flips the holdout paradox (wrong-pocket drops below matched) on all three pairs.** "
            "The unmatched wrong-pocket ≥ matched pattern is consistent with a potency-distribution "
            "sampling difference rather than a pocket-specific inversion. Primary holdout numbers stay unmatched."
        )
    else:
        lines.append(
            "**Potency matching does not give a uniform flip across the three holdout pairs.** "
            "Report Table S13 in full; do not collapse the paradox into a single cause. "
            "Primary holdout numbers stay unmatched."
        )
    lines += [
        "",
        "## What this is not",
        "",
        "- Not a new docking run and not a replacement of Table S8 / HOLDOUT_VERDICT.",
        "- Matched subsets are small (n≈20 before matching); n_min<8 cells are underpowered.",
        "- Does not claim the main-panel pocket-matched > wrong-pocket gap is explained.",
        "",
        "```bash",
        "python3 Dual_Target_Docking/data/jcim_holdout_v0/scripts/wrong_pocket_potency_match_v1.py",
        "```",
        "",
    ]
    return "\n".join(lines) + "\n"


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)
    dist, contrast = [], []
    for prefix in PAIRS:
        hold = load_holdout(prefix)
        mainp = load_main(prefix)
        print(
            f"{prefix}: holdout {len(hold)}  main {len(mainp)}  "
            f"D/A/B={sum(r['cls']=='dual' for r in hold)}/"
            f"{sum(r['cls']=='A_only' for r in hold)}/"
            f"{sum(r['cls']=='B_only' for r in hold)}",
            flush=True,
        )
        dist.extend(dist_rows(prefix, hold, mainp))
        contrast.extend(matched_rows(prefix, hold))
    summary = summary_min_rows(contrast)
    write_csv(TAB / "holdout_vs_main_potency_size_v1.csv", dist)
    write_csv(TAB / "holdout_matched_wrong_pocket_v1.csv", contrast)
    write_csv(TAB / "holdout_matched_wrong_pocket_summary_v1.csv", summary)
    (AN / "HOLDOUT_WRONG_POCKET_POTENCY_VERDICT_V1.md").write_text(
        verdict_md(dist, contrast, summary)
    )
    print("wrote tables +", AN / "HOLDOUT_WRONG_POCKET_POTENCY_VERDICT_V1.md")
    for r in summary:
        print(
            f"{r['pair']:16s} {r['family']:16s} {r['aggregation']:15s} "
            f"min={r['summary_min']} n={r['n_min']} D/A={r['auroc_D_vs_A']} D/B={r['auroc_D_vs_B']}"
        )


if __name__ == "__main__":
    main()
