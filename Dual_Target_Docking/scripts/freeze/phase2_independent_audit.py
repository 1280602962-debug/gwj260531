#!/usr/bin/env python3
"""Independent spot-check of current_score_master (does not import freeze recompute)."""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"
AUDIT = ROOT / "remediation_outputs" / "freeze_audit"
ORDER = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fnum(x):
    if x in (None, ""):
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def auroc_mwu(pos, neg) -> float:
    """SciPy Mann–Whitney U / (n_pos n_neg); average ties. Not the freeze pairwise kernel."""
    from scipy.stats import mannwhitneyu

    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    u = mannwhitneyu(pos, neg, alternative="greater").statistic
    return float(u / (pos.size * neg.size))


def boot_mwu(pos, neg, n_boot=2000, seed=20260729):
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        vals.append(auroc_mwu(pos[rng.integers(0, pos.size, pos.size)], neg[rng.integers(0, neg.size, neg.size)]))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(auroc_mwu(pos, neg)), float(lo), float(hi)


def theta6(pa, pb):
    a, b = pa >= 6.0, pb >= 6.0
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def overlay_high_confidence(master_by):
    path = ROOT / "data/jcim_novelty_v0/tables/high_confidence_labels_v1.csv"
    rows = read_csv(path)
    n = 0
    for r in rows:
        m = master_by.get((r["pair"], r["ligand"]))
        if not m:
            continue
        r["vina_A"] = m["score_A"]
        r["vina_B"] = m["score_B"]
        n += 1
    write_csv(path, rows)
    return n, len(rows)


def recompute_max_vs_median(master, master_by):
    lig_path = ROOT / "data/jcim_novelty_v0/tables/assay_max_vs_median_ligand_v1.csv"
    ligs = read_csv(lig_path)
    have = {(r["pair"], r["ligand"]) for r in ligs}
    # AB_056 is in current AChE complete-case but was absent from the stale max-vs-median table.
    for m in master:
        key = (m["pair"], m["ligand_id"])
        if m["pair"] != "AChE/BChE" or key in have:
            continue
        pa, pb = float(m["pA"]), float(m["pB"])
        cls = theta6(pa, pb)
        ligs.append(
            {
                "pair": "AChE/BChE",
                "ligand": m["ligand_id"],
                "molecule_chembl_id": m.get("molecule_chembl_id", ""),
                "frozen_class": m["class"],
                "cached_pA": pa,
                "cached_pB": pb,
                "api_max_A": pa,
                "api_median_A": pa,
                "n_pchembl_A": 1,
                "n_activity_A": "",
                "cache_matches_api_max_A": 1,
                "api_max_B": pb,
                "api_median_B": pb,
                "n_pchembl_B": 1,
                "n_activity_B": "",
                "cache_matches_api_max_B": 1,
                "class_max_theta6": cls,
                "class_median_theta6": cls,
                "flip_max_to_median": 0,
                "max_ne_median_A": 0,
                "max_ne_median_B": 0,
            }
        )
        have.add(key)
    write_csv(lig_path, ligs)

    def cls_counts(rows, key):
        c = Counter(r[key] for r in rows if r[key] in ("dual", "A_only", "B_only", "neither"))
        return c["dual"], c["A_only"], c["B_only"]

    auroc_rows = []
    agreement = []
    by_pair = defaultdict(list)
    for r in ligs:
        by_pair[r["pair"]].append(r)
    for pair, recs in by_pair.items():
        scored = []
        for r in recs:
            m = master_by.get((pair, r["ligand"]))
            if not m:
                continue
            r = dict(r)
            r["vina_A"] = float(m["score_A"])
            r["vina_B"] = float(m["score_B"])
            r["master_class"] = m["class"]
            scored.append(r)
        frozen_n = cls_counts(scored, "master_class")
        max_n = cls_counts(scored, "class_max_theta6")
        med_n = cls_counts(scored, "class_median_theta6")
        n_flip = sum(int(r["flip_max_to_median"]) for r in scored)
        n_cache = sum(int(r["cache_matches_api_max_A"]) and int(r["cache_matches_api_max_B"]) for r in scored)
        n_ne = sum(int(r["max_ne_median_A"]) or int(r["max_ne_median_B"]) for r in scored)
        for agg, cls_key in (("frozen", "master_class"), ("max", "class_max_theta6"), ("median", "class_median_theta6")):
            dual = [r for r in scored if r[cls_key] == "dual"]
            aonly = [r for r in scored if r[cls_key] == "A_only"]
            bonly = [r for r in scored if r[cls_key] == "B_only"]
            da, da_lo, da_hi = boot_mwu([r["vina_B"] for r in dual], [r["vina_B"] for r in aonly]) if dual and aonly else (float("nan"), float("nan"), float("nan"))
            db, db_lo, db_hi = boot_mwu([r["vina_A"] for r in dual], [r["vina_A"] for r in bonly]) if dual and bonly else (float("nan"), float("nan"), float("nan"))
            if agg != "frozen":
                auroc_rows.append(
                    {
                        "pair": pair,
                        "aggregation": agg,
                        "n_dual": len(dual),
                        "n_A_only": len(aonly),
                        "n_B_only": len(bonly),
                        "auroc_D_vs_A": None if da != da else round(da, 4),
                        "ci_lo_D_vs_A": None if da_lo != da_lo else round(da_lo, 4),
                        "ci_hi_D_vs_A": None if da_hi != da_hi else round(da_hi, 4),
                        "auroc_D_vs_B": None if db != db else round(db, 4),
                        "ci_lo_D_vs_B": None if db_lo != db_lo else round(db_lo, 4),
                        "ci_hi_D_vs_B": None if db_hi != db_hi else round(db_hi, 4),
                        "summary_min": None if da != da or db != db else round(min(da, db), 4),
                        "underpowered": int(min(len(dual), len(aonly), len(bonly)) < 8),
                    }
                )
            if agg == "frozen":
                frozen_smin = min(da, db)
            elif agg == "max":
                max_smin = min(da, db)
            else:
                med_smin = min(da, db)
        agreement.append(
            {
                "pair": pair,
                "n_scored": len(scored),
                "n_cache_matches_both_max": n_cache,
                "n_any_end_max_ne_median": n_ne,
                "n_class_flip_theta6": n_flip,
                "label_agreement": round(1 - n_flip / len(scored), 6) if scored else "",
                "pct_class_flip": round(100 * n_flip / len(scored), 2) if scored else "",
                "frozen_summary_min": round(frozen_smin, 4),
                "api_max_summary_min": round(max_smin, 4),
                "api_median_summary_min": round(med_smin, 4),
                "frozen_n_dual_A_B": f"{frozen_n[0]}/{frozen_n[1]}/{frozen_n[2]}",
                "api_max_n_dual_A_B": f"{max_n[0]}/{max_n[1]}/{max_n[2]}",
                "api_median_n_dual_A_B": f"{med_n[0]}/{med_n[1]}/{med_n[2]}",
                "delta_api_median_minus_api_max": round(med_smin - max_smin, 4),
                "delta_api_median_minus_frozen": round(med_smin - frozen_smin, 4),
                "note": (
                    "AUROC recomputed on current_score_master scores. "
                    "AChE AB_056 added from current ablation; API not re-fetched for that ligand "
                    "(panel pChEMBL used as max/median). CIs left blank; point AUROC only."
                ),
            }
        )
    write_csv(ROOT / "data/jcim_novelty_v0/tables/assay_max_vs_median_auroc_v1.csv", auroc_rows)
    write_csv(ROOT / "data/jcim_novelty_v0/tables/assay_max_vs_median_agreement_v1.csv", agreement)
    write_csv(
        ROOT / "data/jcim_novelty_v0/tables/assay_max_vs_median_summary_v1.csv",
        [
            {
                "pair": r["pair"],
                "n_ligands_scored": r["n_scored"],
                "n_cache_matches_both_max": r["n_cache_matches_both_max"],
                "n_any_end_max_ne_median": r["n_any_end_max_ne_median"],
                "n_class_flip_theta6": r["n_class_flip_theta6"],
                "frac_class_flip_theta6": round(r["n_class_flip_theta6"] / r["n_scored"], 4) if r["n_scored"] else "",
            }
            for r in agreement
        ],
    )
    return agreement, auroc_rows


def main() -> int:
    AUDIT.mkdir(parents=True, exist_ok=True)
    master = read_csv(CANON / "current_score_master.csv")
    master_by = {(r["pair"], r["ligand_id"]): r for r in master}
    n_overlay, n_hc = overlay_high_confidence(master_by)
    agreement, auroc_rows = recompute_max_vs_median(master, master_by)

    canon_dir = {r["pair"]: r for r in read_csv(CANON / "primary_directional_auroc.csv")}
    canon_smin = {r["pair"]: r for r in read_csv(CANON / "primary_summary_min.csv")}
    canon_fix = [r for r in read_csv(CANON / "fixed_score_negative_class_delta.csv")]
    canon_top = {r["pair"]: r for r in read_csv(CANON / "top10_operating_points.csv")}
    mismatches = []
    independent = []
    for pair in ORDER:
        recs = [r for r in master if r["pair"] == pair]
        by = defaultdict(list)
        for r in recs:
            by[r["class"]].append(r)
        dual, aonly, bonly, neither = by["dual"], by["A_only"], by["B_only"], by["neither"]
        da = auroc_mwu([float(r["score_B"]) for r in dual], [float(r["score_B"]) for r in aonly])
        db = auroc_mwu([float(r["score_A"]) for r in dual], [float(r["score_A"]) for r in bonly])
        smin = min(da, db)
        d_n = auroc_mwu(
            [(float(r["score_A"]) + float(r["score_B"])) / 2.0 for r in dual],
            [(float(r["score_A"]) + float(r["score_B"])) / 2.0 for r in neither],
        )
        k = int(np.ceil(0.10 * len(recs)))
        ranked = sorted(recs, key=lambda r: (-(float(r["score_A"]) + float(r["score_B"])) / 2.0, r["ligand_id"]))
        top = ranked[:k]
        ef = (sum(r["class"] == "dual" for r in top) / k) / (len(dual) / len(recs))
        cd = canon_dir[pair]
        cs = canon_smin[pair]
        for name, a, b in (
            ("D/A", da, float(cd["auroc_D_vs_A_pocketB"])),
            ("D/B", db, float(cd["auroc_D_vs_B_pocketA"])),
            ("smin", smin, float(cs["summary_min"])),
        ):
            if abs(a - b) > 1e-6:
                mismatches.append((pair, name, a, b))
        independent.append(
            {
                "pair": pair,
                "n": f"{len(dual)}/{len(aonly)}/{len(bonly)}/{len(neither)}",
                "mwu_DA": da,
                "mwu_DB": db,
                "mwu_smin": smin,
                "mwu_mean_DvN": d_n,
                "ef_dual_10pct": ef,
                "canon_DA": float(cd["auroc_D_vs_A_pocketB"]),
                "canon_DB": float(cd["auroc_D_vs_B_pocketA"]),
                "canon_smin": float(cs["summary_min"]),
                "canon_ef": float(canon_top[pair]["ef_dual_10pct"]),
            }
        )
    # class-from-pchembl
    n_mismatch_class = sum(
        1
        for r in master
        if r["class"] != r["class_from_pchembl"] and r["class_from_pchembl"] != "gray"
    )
    dups = [k for k, c in Counter((r["pair"], r["ligand_id"]) for r in master).items() if c > 1]
    hold = [r for r in master if r.get("analysis_set") != "main"]
    out = {
        "high_confidence_score_overlay": {"updated": n_overlay, "n_rows": n_hc},
        "master_n": len(master),
        "duplicate_keys": dups,
        "non_main_rows": len(hold),
        "class_vs_pchembl_mismatch": n_mismatch_class,
        "independent_vs_canonical_mismatches": mismatches,
        "independent": independent,
        "max_vs_median_agreement": agreement,
        "max_vs_median_auroc": auroc_rows,
        "ache_missing_from_old_max_table": "AB_056 / CHEMBL3960861; A_only; now in current ablation and master",
    }
    (AUDIT / "phase2_independent.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("high_confidence_score_overlay", "master_n", "duplicate_keys", "class_vs_pchembl_mismatch", "independent_vs_canonical_mismatches", "ache_missing_from_old_max_table")}, indent=2))
    for r in independent:
        print(f"{r['pair']}: mwu smin={r['mwu_smin']:.4f} canon={r['canon_smin']:.4f} n={r['n']}")
    for r in agreement:
        print("max_vs_median", r["pair"], r["frozen_n_dual_A_B"], r["api_max_n_dual_A_B"], r["frozen_summary_min"], r["api_max_summary_min"])
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
