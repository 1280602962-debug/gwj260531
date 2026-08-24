#!/usr/bin/env python3
"""Claim-hardening analyses on frozen DualFourClass scores. No new docking.

1. summary_min vs arithmetic mean vs harmonic mean of the two directional AUROCs.
2. All four prespecified descriptor directional AUROCs (not best-of-4 only).
3. Docking attempted / successful / failed census on main panels and holdout.

Primary endpoint remains pocket-matched Vina summary_min. These tables are
sensitivity / transparency, not a second primary metric.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_novelty_v0"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SPEC = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        cls="class",
        ligand="ligand",
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        cls="class",
        ligand="ligand",
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        cls="class",
        ligand="ligand",
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        cls="class",
        ligand="ligand",
    ),
}

HOLDOUT = [
    dict(pair="AChE/BChE", prefix="HOAB", path="data/jcim_holdout_v0/tables/scores_vina_mode1_HOAB.csv"),
    dict(pair="PIK3CA/PIK3CB", prefix="HOAP", path="data/jcim_holdout_v0/tables/scores_vina_mode1_HOAP.csv"),
    dict(pair="PIK3CA/mTOR", prefix="HOPM", path="data/jcim_holdout_v0/tables/scores_vina_mode1_HOPM.csv"),
]

ORDER = ["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]
DIRECTIONAL = ROOT / "data" / "jcim_bench_v0" / "tables" / "pocket_matched_directional_v1.csv"
ASSEMBLED = ROOT / "data" / "jcim_bench_v0" / "tables" / "assembled_all_pairs_long.csv"
FORM = OUT / "tables" / "formulation_conventional_vs_directional_v1.csv"


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def load_csv(p: Path):
    with p.open() as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def harmonic_mean(a: float, b: float) -> float:
    if a <= 0 or b <= 0:
        return 0.0
    return 2.0 / (1.0 / a + 1.0 / b)


def aggregation_sensitivity() -> list[dict]:
    rows = [r for r in load_csv(DIRECTIONAL) if r["variant"] == "pocket_matched_vina"]
    by = {r["pair"]: r for r in rows}
    form = load_csv(FORM)
    nei = {}
    for r in form:
        if r["formulation"] == "conventional_dual_vs_neither" and r["contrast"] == "D_vs_neither_mean":
            nei[r["pair"]] = r
    out = []
    rank_store = {"min": [], "mean": [], "harmonic_mean": []}
    for pair in ORDER:
        r = by[pair]
        a = float(r["auroc_D_vs_A"])
        b = float(r["auroc_D_vs_B"])
        vals = {
            "min": min(a, b),
            "mean": (a + b) / 2.0,
            "harmonic_mean": harmonic_mean(a, b),
        }
        n = nei.get(pair, {})
        dual_neither = fnum(n.get("auroc"))
        rec = {
            "pair": pair,
            "n_dual": r["n_dual"],
            "n_A_only": r["n_A_only"],
            "n_B_only": r["n_B_only"],
            "auroc_D_vs_A": round(a, 4),
            "auroc_D_vs_B": round(b, 4),
            "summary_min": round(vals["min"], 4),
            "summary_mean": round(vals["mean"], 4),
            "summary_harmonic": round(vals["harmonic_mean"], 4),
            "dual_vs_neither_vina_mean": round(dual_neither, 4) if dual_neither is not None else "",
            "n_neither": n.get("n_neg", ""),
            "neither_underpowered": n.get("underpowered", ""),
            "gap_neither_minus_min": round(dual_neither - vals["min"], 4) if dual_neither is not None else "",
            "gap_neither_minus_mean": round(dual_neither - vals["mean"], 4) if dual_neither is not None else "",
            "gap_neither_minus_harmonic": round(dual_neither - vals["harmonic_mean"], 4) if dual_neither is not None else "",
            "note": "primary remains summary_min; mean/harmonic are sensitivity only; Dual-vs-neither uses a different negative set and is not a paired Δ",
        }
        out.append(rec)
        for k, v in vals.items():
            rank_store[k].append((v, pair))
    ranks = {k: {p: i + 1 for i, (_, p) in enumerate(sorted(vs, reverse=True))} for k, vs in rank_store.items()}
    for rec in out:
        rec["rank_min"] = ranks["min"][rec["pair"]]
        rec["rank_mean"] = ranks["mean"][rec["pair"]]
        rec["rank_harmonic"] = ranks["harmonic_mean"][rec["pair"]]
        rec["ranking_unchanged"] = int(
            rec["rank_min"] == rec["rank_mean"] == rec["rank_harmonic"]
        )
    return out


def descriptor_all_four() -> list[dict]:
    rows = load_csv(ASSEMBLED)
    desc = ["heavy", "mw", "clogp", "tpsa"]
    out = []
    for pair in ORDER:
        sub = [r for r in rows if r["pair"] == pair]
        dual = [r for r in sub if r["cls"] == "dual"]
        aonly = [r for r in sub if r["cls"] == "A_only"]
        bonly = [r for r in sub if r["cls"] == "B_only"]
        rec = {
            "pair": pair,
            "n_dual": len(dual),
            "n_A_only": len(aonly),
            "n_B_only": len(bonly),
        }
        mins = {}
        for d in desc:
            y_a = np.array([1] * len(dual) + [0] * len(aonly))
            s_a = np.array([float(r[d]) for r in dual + aonly])
            y_b = np.array([1] * len(dual) + [0] * len(bonly))
            s_b = np.array([float(r[d]) for r in dual + bonly])
            auc_a = float(roc_auc_score(y_a, s_a))
            auc_b = float(roc_auc_score(y_b, s_b))
            rec[f"{d}_D_vs_A"] = round(auc_a, 4)
            rec[f"{d}_D_vs_B"] = round(auc_b, 4)
            rec[f"{d}_summary_min"] = round(min(auc_a, auc_b), 4)
            mins[d] = min(auc_a, auc_b)
        best = max(mins, key=mins.get)
        rec["best_single_descriptor"] = best
        rec["best_single_descriptor_summary_min"] = round(mins[best], 4)
        rec["note"] = (
            "best_single_descriptor is the highest of four prespecified descriptors on this panel; "
            "descriptive reference only, not a confirmatory competitor"
        )
        out.append(rec)
    return out


def docking_census() -> list[dict]:
    out = []
    for pair, spec in SPEC.items():
        rows = load_csv(ROOT / spec["scores"])
        n_attempted = len(rows)
        ok_a = ok_b = both = fail_a = fail_b = fail_either = 0
        fail_ligands = []
        for r in rows:
            a = fnum(r.get(spec["vina_a"]))
            b = fnum(r.get(spec["vina_b"]))
            a_ok = a is not None
            b_ok = b is not None
            if a_ok:
                ok_a += 1
            else:
                fail_a += 1
            if b_ok:
                ok_b += 1
            else:
                fail_b += 1
            if a_ok and b_ok:
                both += 1
            else:
                fail_either += 1
                fail_ligands.append(
                    f"{r.get(spec['ligand'], '')}:{r.get(spec['cls'], '')}:A={'ok' if a_ok else 'fail'}:B={'ok' if b_ok else 'fail'}"
                )
        out.append(
            {
                "set": "main_panel",
                "pair": pair,
                "n_attempted": n_attempted,
                "n_success_pocket_A": ok_a,
                "n_success_pocket_B": ok_b,
                "n_success_both_ends": both,
                "n_fail_pocket_A": fail_a,
                "n_fail_pocket_B": fail_b,
                "n_fail_either_end": fail_either,
                "fail_rate_either": round(fail_either / n_attempted, 4) if n_attempted else "",
                "fail_rate_pocket_A": round(fail_a / n_attempted, 4) if n_attempted else "",
                "fail_rate_pocket_B": round(fail_b / n_attempted, 4) if n_attempted else "",
                "failed_ligands": ";".join(fail_ligands),
                "note": "AUROC tables condition on both-end scores; failures are engine/prep coverage, not silent missingness",
            }
        )
    for h in HOLDOUT:
        rows = load_csv(ROOT / h["path"])
        by_lig = {}
        for r in rows:
            by_lig.setdefault(r["ligand"], {})[r["receptor_key"]] = r
        n_attempted = len(by_lig)
        ok_a = ok_b = both = fail_a = fail_b = fail_either = 0
        fail_ligands = []
        for lig, d in sorted(by_lig.items()):
            ra, rb = d.get("A", {}), d.get("B", {})
            a_ok = ra.get("status") == "success" and fnum(ra.get("vina_mode1")) is not None
            b_ok = rb.get("status") == "success" and fnum(rb.get("vina_mode1")) is not None
            if a_ok:
                ok_a += 1
            else:
                fail_a += 1
            if b_ok:
                ok_b += 1
            else:
                fail_b += 1
            if a_ok and b_ok:
                both += 1
            else:
                fail_either += 1
                reason_a = (ra.get("reason") or "").replace("\n", " ")[:80]
                reason_b = (rb.get("reason") or "").replace("\n", " ")[:80]
                fail_ligands.append(f"{lig}:A={'ok' if a_ok else reason_a}:B={'ok' if b_ok else reason_b}")
        out.append(
            {
                "set": "holdout",
                "pair": h["pair"],
                "n_attempted": n_attempted,
                "n_success_pocket_A": ok_a,
                "n_success_pocket_B": ok_b,
                "n_success_both_ends": both,
                "n_fail_pocket_A": fail_a,
                "n_fail_pocket_B": fail_b,
                "n_fail_either_end": fail_either,
                "fail_rate_either": round(fail_either / n_attempted, 4) if n_attempted else "",
                "fail_rate_pocket_A": round(fail_a / n_attempted, 4) if n_attempted else "",
                "fail_rate_pocket_B": round(fail_b / n_attempted, 4) if n_attempted else "",
                "failed_ligands": ";".join(fail_ligands),
                "note": "HOAP_028 is a boron AutoDock atom-type B coverage failure on both ends, excluded from AUROC (59/60)",
            }
        )
    return out


def main():
    agg = aggregation_sensitivity()
    write_csv(TAB / "aggregation_min_mean_harmonic_v1.csv", agg)
    desc = descriptor_all_four()
    write_csv(TAB / "descriptor_all_four_directional_v1.csv", desc)
    census = docking_census()
    write_csv(TAB / "docking_failure_census_v1.csv", census)
    print("aggregation ranking unchanged on all pairs:", all(r["ranking_unchanged"] for r in agg))
    for r in agg:
        print(
            f"  {r['pair']}: min={r['summary_min']} mean={r['summary_mean']} "
            f"harm={r['summary_harmonic']} ranks={r['rank_min']}/{r['rank_mean']}/{r['rank_harmonic']}"
        )
    print("descriptor best-single:")
    for r in desc:
        print(f"  {r['pair']}: {r['best_single_descriptor']}={r['best_single_descriptor_summary_min']}")
    print("docking census:")
    for r in census:
        print(
            f"  {r['set']} {r['pair']}: attempted={r['n_attempted']} both={r['n_success_both_ends']} "
            f"fail_either={r['n_fail_either_end']} ligands={r['failed_ligands'][:120]}"
        )


if __name__ == "__main__":
    main()
