#!/usr/bin/env python3
"""Completeness audit, independent alternative score master, paired Δ statistics.

Official stats only on the common complete-case population with identical
labels, directional definition, and paired bootstrap indices.
Unreplaced direction is the internal negative control.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import auroc  # noqa: E402
from rerun_v4_2.alt_receptor_config import (  # noqa: E402
    ALTS,
    BOOT_B,
    BOOT_SEED,
    FORMAL_PAIRS,
    PRIMARY,
    RUN,
    TOP_FRACTION,
)

QA = RUN / "13_qa"
PRIMARY_MASTER = QA / "official_primary_seed42_score_master_formal4.csv"
ALT_JOBS = QA / "alt_production_seed42_master.csv"
OUT_MASTER = QA / "alt_score_master_independent.csv"
OUT_AUDIT = QA / "alt_production_completeness_audit.json"
OUT_STATS = QA / "alt_paired_delta_stats.json"
OUT_CSV = QA / "alt_paired_delta_stats.csv"


def _score(mode1, status) -> float | None:
    if status != "SUCCESS":
        return None
    if mode1 in ("", None):
        return None
    try:
        return -float(mode1)
    except (TypeError, ValueError):
        return None


def ef_at_k(y: np.ndarray, scores: np.ndarray, frac: float = TOP_FRACTION) -> float:
    n = y.size
    if n == 0 or int(y.sum()) == 0:
        return float("nan")
    k = max(1, int(math.ceil(frac * n)))
    order = np.argsort(-scores, kind="mergesort")
    hits = int(y[order[:k]].sum())
    base = float(y.mean())
    return (hits / k) / base if base > 0 else float("nan")


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2:
        return float("nan")
    ra = np.argsort(np.argsort(-a))
    rb = np.argsort(np.argsort(-b))
    return float(np.corrcoef(ra, rb)[0, 1])


def topk_overlap(scores_p: np.ndarray, scores_a: np.ndarray, eids: list[str], frac: float) -> dict:
    k = max(1, int(math.ceil(frac * len(eids))))
    op = [eids[i] for i in sorted(range(len(eids)), key=lambda i: (-scores_p[i], eids[i]))[:k]]
    oa = [eids[i] for i in sorted(range(len(eids)), key=lambda i: (-scores_a[i], eids[i]))[:k]]
    inter = set(op) & set(oa)
    return {"k": k, "overlap": len(inter), "jaccard": len(inter) / len(set(op) | set(oa))}


def digest(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_primary() -> dict[tuple[str, str], dict]:
    out = {}
    for r in csv.DictReader(PRIMARY_MASTER.open(encoding="utf-8-sig", newline="")):
        out[(r["pair"], r["global_ligand_entity_id"])] = r
    return out


def load_alt_jobs() -> dict[tuple[str, str, str], dict]:
    out = {}
    for r in csv.DictReader(ALT_JOBS.open(encoding="utf-8-sig", newline="")):
        out[(r["pair"], r["pdb_id"], r["global_ligand_entity_id"])] = r
    return out


def build_alt_master() -> tuple[list[dict], dict]:
    prim = load_primary()
    altj = load_alt_jobs()
    rows = []
    expected = 0
    accounted = 0
    for (pair, eid), pr in prim.items():
        for pdb, cfg in ALTS.items():
            if cfg["pair"] != pair:
                continue
            expected += 1
            job = altj.get((pair, pdb, eid))
            if job:
                accounted += 1
            status = (job or {}).get("status", "MISSING_JOB")
            sc = _score((job or {}).get("mode1_affinity"), status)
            rows.append({
                "pair": pair,
                "global_ligand_entity_id": eid,
                "canonical_ligand_id": pr["canonical_ligand_id"],
                "class": pr["class"],
                "activity_eligible": pr["activity_eligible"],
                "alt_pdb": pdb,
                "alt_id": cfg["id"],
                "side": cfg["side"],
                "design": cfg["design"],
                "kept_primary": cfg["kept_primary"],
                "replaced_primary": PRIMARY[pair][0] if cfg["side"] == "A" else PRIMARY[pair][1],
                "status_alt": status,
                "mode1_affinity_alt": (job or {}).get("mode1_affinity", ""),
                "score_alt": "" if sc is None else f"{sc:.6f}",
                "score_A_primary": pr["score_A"],
                "score_B_primary": pr["score_B"],
                "status_A_primary": pr["status_A"],
                "status_B_primary": pr["status_B"],
                "TIMEOUT_treated_as": "missing_not_imputed",
            })
    audit = {
        "jobs_expected_from_formal_mapping": expected,
        "jobs_accounted": accounted,
        "SUCCESS": sum(1 for r in rows if r["status_alt"] == "SUCCESS"),
        "TIMEOUT": sum(1 for r in rows if r["status_alt"] == "TIMEOUT"),
        "PERMANENT_FAIL": sum(1 for r in rows if r["status_alt"] == "PERMANENT_FAIL"),
        "MISSING_JOB": sum(1 for r in rows if r["status_alt"] == "MISSING_JOB"),
        "PRODUCTION_COMPLETE": accounted == expected and all(
            r["status_alt"] in {"SUCCESS", "TIMEOUT", "PERMANENT_FAIL"} for r in rows
        ),
        "AUROC_allowed": False,
    }
    audit["AUROC_allowed"] = bool(audit["PRODUCTION_COMPLETE"])
    return rows, audit


def common_population(prim_rows: list[dict], alt_rows: list[dict], pair: str, pdb: str) -> list[dict]:
    cfg = ALTS[pdb]
    alt_by = {(r["pair"], r["global_ligand_entity_id"]): r for r in alt_rows if r["alt_pdb"] == pdb}
    out = []
    for pr in prim_rows:
        if pr["pair"] != pair:
            continue
        if int(pr.get("in_official_population") or 0) != 1:
            continue
        ar = alt_by.get((pair, pr["global_ligand_entity_id"]))
        if not ar or ar["score_alt"] == "":
            continue
        sa_p = float(pr["score_A"])
        sb_p = float(pr["score_B"])
        salt = float(ar["score_alt"])
        if cfg["side"] == "A":
            sa_a, sb_a = salt, sb_p
        else:
            sa_a, sb_a = sa_p, salt
        out.append({
            "eid": pr["global_ligand_entity_id"],
            "class": pr["class"],
            "score_A_p": sa_p,
            "score_B_p": sb_p,
            "score_A_a": sa_a,
            "score_B_a": sb_a,
            "score_mean_p": 0.5 * (sa_p + sb_p),
            "score_mean_a": 0.5 * (sa_a + sb_a),
        })
    return sorted(out, key=lambda r: r["eid"])


def paired_stats(recs: list[dict], pair: str, pdb: str) -> dict:
    cls = np.array([r["class"] for r in recs])
    eids = [r["eid"] for r in recs]
    sap = np.array([r["score_A_p"] for r in recs])
    sbp = np.array([r["score_B_p"] for r in recs])
    saa = np.array([r["score_A_a"] for r in recs])
    sba = np.array([r["score_B_a"] for r in recs])
    d = np.flatnonzero(cls == "dual")
    a = np.flatnonzero(cls == "A_only")
    b = np.flatnonzero(cls == "B_only")
    da_p = auroc(sbp[d], sbp[a])
    db_p = auroc(sap[d], sap[b])
    da_a = auroc(sba[d], sba[a])
    db_a = auroc(saa[d], saa[b])
    sm_p = min(da_p, db_p)
    sm_a = min(da_a, db_a)
    rng = np.random.default_rng(BOOT_SEED)
    dda = np.empty(BOOT_B)
    ddb = np.empty(BOOT_B)
    dsm = np.empty(BOOT_B)
    for i in range(BOOT_B):
        di = rng.choice(d, size=d.size, replace=True) if d.size else d
        ai = rng.choice(a, size=a.size, replace=True) if a.size else a
        bi = rng.choice(b, size=b.size, replace=True) if b.size else b
        dda[i] = auroc(sba[di], sba[ai]) - auroc(sbp[di], sbp[ai])
        ddb[i] = auroc(saa[di], saa[bi]) - auroc(sap[di], sap[bi])
        smi_p = min(auroc(sbp[di], sbp[ai]), auroc(sap[di], sap[bi]))
        smi_a = min(auroc(sba[di], sba[ai]), auroc(saa[di], saa[bi]))
        dsm[i] = smi_a - smi_p
    y_dual = (cls == "dual").astype(float)
    side = ALTS[pdb]["side"]
    rec = {
        "pair": pair,
        "alt_pdb": pdb,
        "alt_id": ALTS[pdb]["id"],
        "design": ALTS[pdb]["design"],
        "replaced_side": side,
        "n_common": len(recs),
        "n_dual": int(d.size),
        "n_A_only": int(a.size),
        "n_B_only": int(b.size),
        "auroc_D_vs_A_primary": float(da_p),
        "auroc_D_vs_A_alt": float(da_a),
        "delta_auroc_D_vs_A": float(da_a - da_p),
        "delta_auroc_D_vs_A_ci_lo": float(np.percentile(dda, 2.5)),
        "delta_auroc_D_vs_A_ci_hi": float(np.percentile(dda, 97.5)),
        "auroc_D_vs_B_primary": float(db_p),
        "auroc_D_vs_B_alt": float(db_a),
        "delta_auroc_D_vs_B": float(db_a - db_p),
        "delta_auroc_D_vs_B_ci_lo": float(np.percentile(ddb, 2.5)),
        "delta_auroc_D_vs_B_ci_hi": float(np.percentile(ddb, 97.5)),
        "summary_min_primary": float(sm_p),
        "summary_min_alt": float(sm_a),
        "delta_summary_min": float(sm_a - sm_p),
        "delta_summary_min_ci_lo": float(np.percentile(dsm, 2.5)),
        "delta_summary_min_ci_hi": float(np.percentile(dsm, 97.5)),
        "affected_direction": "D_vs_B_pocketA" if side == "A" else "D_vs_A_pocketB",
        "negative_control_direction": "D_vs_A_pocketB" if side == "A" else "D_vs_B_pocketA",
        "ef10_dual_primary": ef_at_k(y_dual, 0.5 * (sap + sbp)),
        "ef10_dual_alt": ef_at_k(y_dual, 0.5 * (saa + sba)),
        "spearman_score_mean": spearman(0.5 * (sap + sbp), 0.5 * (saa + sba)),
        **{f"topk_{k}": v for k, v in topk_overlap(0.5 * (sap + sbp), 0.5 * (saa + sba), eids, TOP_FRACTION).items()},
        "bootstrap_B": BOOT_B,
        "bootstrap_seed": BOOT_SEED,
        "paired_identical_indices": True,
        "common_population_sha256": digest([r["eid"] for r in recs]),
    }
    rec["delta_affected"] = rec["delta_auroc_D_vs_B"] if side == "A" else rec["delta_auroc_D_vs_A"]
    rec["delta_control"] = rec["delta_auroc_D_vs_A"] if side == "A" else rec["delta_auroc_D_vs_B"]
    return rec


def main() -> int:
    alt_rows, audit = build_alt_master()
    fields = list(alt_rows[0].keys()) if alt_rows else []
    with OUT_MASTER.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(alt_rows)
    OUT_AUDIT.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["AUROC_allowed"]:
        print("STATS_BLOCKED_UNTIL_COMPLETENESS")
        return 1
    prim_rows = list(csv.DictReader(PRIMARY_MASTER.open(encoding="utf-8-sig", newline="")))
    stats = []
    for pdb, cfg in ALTS.items():
        recs = common_population(prim_rows, alt_rows, cfg["pair"], pdb)
        stats.append(paired_stats(recs, cfg["pair"], pdb))
    OUT_STATS.write_text(json.dumps({
        "official": True,
        "independent_of_primary_master_file": True,
        "promoted_to_primary": False,
        "rows": stats,
    }, indent=2) + "\n", encoding="utf-8")
    with OUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(stats[0].keys()), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(stats)
    for s in stats:
        print(
            f"{s['pair']} {s['alt_pdb']} n={s['n_common']} "
            f"dAff={s['delta_affected']:+.4f} dCtrl={s['delta_control']:+.4f} "
            f"dMin={s['delta_summary_min']:+.4f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
