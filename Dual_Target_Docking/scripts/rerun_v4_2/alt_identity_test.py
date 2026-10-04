#!/usr/bin/env python3
"""Build official A1/B1 seed42 master for the four formal pairs and replay
the receptor-sensitivity analysis path twice. Identity test only; no new science.

HISTORICAL_NON_FORMAL_ENTRY — documentary banner only.
This file is not a Stage 2 scientific input. Do not import it from
jcim_stage2_*. See 13_qa/HISTORICAL_NON_FORMAL_FILES.md.
Compute logic below is unchanged.
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
    BOOT_B,
    BOOT_SEED,
    FORMAL_PAIRS,
    PRIMARY,
    RUN,
    TOP_FRACTION,
)

QA = RUN / "13_qa"
PROTO = RUN / "00_protocol"
MAP = PROTO / "pair_ligand_mapping.csv"
PHASE7 = QA / "phase7_fiveseed_master.csv"
OUT_MASTER = QA / "official_primary_seed42_score_master_formal4.csv"
OUT_ID = QA / "alt_analysis_identity_test.json"


def _score(mode1) -> float | None:
    if mode1 in ("", None):
        return None
    try:
        return -float(mode1)
    except (TypeError, ValueError):
        return None


def load_mapping() -> list[dict]:
    rows = []
    for r in csv.DictReader(MAP.open(encoding="utf-8-sig", newline="")):
        if r["pair"] not in FORMAL_PAIRS:
            continue
        if r.get("parent_independent") != "1" or r.get("prepare_3d") != "1":
            continue
        if r.get("analysis_set") != "main":
            continue
        rows.append(r)
    return rows


def load_seed42() -> dict[tuple[str, str, str], dict]:
    out = {}
    for r in csv.DictReader(PHASE7.open(encoding="utf-8-sig", newline="")):
        if str(r.get("seed")) != "42":
            continue
        if r["pair"] not in FORMAL_PAIRS:
            continue
        out[(r["pair"], r["pdb_id"], r["global_ligand_entity_id"])] = r
    return out


def build_master() -> list[dict]:
    mapping = load_mapping()
    seed42 = load_seed42()
    rows = []
    for r in mapping:
        pair = r["pair"]
        eid = r["global_ligand_entity_id"]
        pdb_a, pdb_b = PRIMARY[pair]
        ja = seed42.get((pair, pdb_a, eid), {})
        jb = seed42.get((pair, pdb_b, eid), {})
        sa = _score(ja.get("mode1_affinity")) if ja.get("status") == "SUCCESS" else None
        sb = _score(jb.get("mode1_affinity")) if jb.get("status") == "SUCCESS" else None
        complete = sa is not None and sb is not None
        eligible = r.get("activity_eligible") == "1"
        rows.append({
            "pair": pair,
            "global_ligand_entity_id": eid,
            "canonical_ligand_id": r["canonical_ligand_id"],
            "class": r["class"],
            "activity_eligible": r.get("activity_eligible", ""),
            "analysis_set": r.get("analysis_set", ""),
            "pdb_A": pdb_a,
            "pdb_B": pdb_b,
            "status_A": ja.get("status", "MISSING_JOB"),
            "status_B": jb.get("status", "MISSING_JOB"),
            "mode1_affinity_A": ja.get("mode1_affinity", ""),
            "mode1_affinity_B": jb.get("mode1_affinity", ""),
            "score_A": "" if sa is None else f"{sa:.6f}",
            "score_B": "" if sb is None else f"{sb:.6f}",
            "score_mean": "" if not complete else f"{0.5 * (sa + sb):.6f}",
            "complete_case": int(complete),
            "in_official_population": int(complete and eligible),
            "D_vs_A_uses": "score_B",
            "D_vs_B_uses": "score_A",
            "TIMEOUT_treated_as": "missing_not_imputed",
        })
    return rows


def _f(x) -> float:
    return float(x)


def analysis_payload(master: list[dict], pair: str) -> dict:
    recs = [
        r for r in master
        if r["pair"] == pair and int(r["in_official_population"]) == 1
    ]
    recs = sorted(recs, key=lambda r: r["global_ligand_entity_id"])
    cls = np.array([r["class"] for r in recs])
    sa = np.array([_f(r["score_A"]) for r in recs], dtype=float)
    sb = np.array([_f(r["score_B"]) for r in recs], dtype=float)
    mean = 0.5 * (sa + sb)
    eids = [r["global_ligand_entity_id"] for r in recs]
    d = np.flatnonzero(cls == "dual")
    a = np.flatnonzero(cls == "A_only")
    b = np.flatnonzero(cls == "B_only")
    n = np.flatnonzero(cls == "neither")
    da_pos, da_neg = sb[d], sb[a]
    db_pos, db_neg = sa[d], sa[b]
    da = auroc(da_pos, da_neg)
    db = auroc(db_pos, db_neg)
    k = int(math.ceil(TOP_FRACTION * len(recs)))
    order = sorted(range(len(recs)), key=lambda i: (-mean[i], eids[i]))
    topk = [eids[i] for i in order[:k]]
    rng = np.random.default_rng(BOOT_SEED)
    boot = []
    for _ in range(BOOT_B):
        di = rng.choice(d, size=d.size, replace=True).tolist() if d.size else []
        ai = rng.choice(a, size=a.size, replace=True).tolist() if a.size else []
        bi = rng.choice(b, size=b.size, replace=True).tolist() if b.size else []
        boot.append({"dual": di, "A_only": ai, "B_only": bi})
    return {
        "pair": pair,
        "n_mapped_independent_prepare3d": sum(1 for r in master if r["pair"] == pair),
        "n_official_complete_case": len(recs),
        "n_dual": int(d.size),
        "n_A_only": int(a.size),
        "n_B_only": int(b.size),
        "n_neither": int(n.size),
        "membership_eids": eids,
        "class": cls.tolist(),
        "score_A": [round(float(x), 6) for x in sa],
        "score_B": [round(float(x), 6) for x in sb],
        "directional_mapping": {"D_vs_A": "score_B", "D_vs_B": "score_A"},
        "auroc_inputs": {
            "D_vs_A_pos_score_B": [round(float(x), 6) for x in da_pos],
            "D_vs_A_neg_score_B": [round(float(x), 6) for x in da_neg],
            "D_vs_B_pos_score_A": [round(float(x), 6) for x in db_pos],
            "D_vs_B_neg_score_A": [round(float(x), 6) for x in db_neg],
        },
        "auroc_D_vs_A_pocketB": None if math.isnan(da) else float(da),
        "auroc_D_vs_B_pocketA": None if math.isnan(db) else float(db),
        "summary_min": None if (math.isnan(da) or math.isnan(db)) else float(min(da, db)),
        "topk_k": k,
        "topk_sort_inputs": [{"eid": eids[i], "score_mean": round(float(mean[i]), 6)} for i in order],
        "topk_eids": topk,
        "bootstrap": {"B": BOOT_B, "seed": BOOT_SEED, "scheme": "class_stratified_shared_dual"},
        "bootstrap_indices": boot,
    }


def digest(obj) -> str:
    blob = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def main() -> int:
    master = build_master()
    fields = list(master[0].keys()) if master else []
    QA.mkdir(parents=True, exist_ok=True)
    with OUT_MASTER.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(master)

    pass1 = {p: analysis_payload(master, p) for p in FORMAL_PAIRS}
    pass2 = {p: analysis_payload(master, p) for p in FORMAL_PAIRS}
    cmp_keys = [
        "membership_eids", "class", "score_A", "score_B", "directional_mapping",
        "auroc_inputs", "topk_sort_inputs", "bootstrap_indices",
    ]
    pair_ok = {}
    for p in FORMAL_PAIRS:
        a = {k: pass1[p][k] for k in cmp_keys}
        b = {k: pass2[p][k] for k in cmp_keys}
        pair_ok[p] = digest(a) == digest(b)

    slim1 = {}
    for p, rec in pass1.items():
        slim = {k: rec[k] for k in rec if k != "bootstrap_indices"}
        slim["bootstrap_indices_sha256"] = digest(rec["bootstrap_indices"])
        slim1[p] = slim

    ok = all(pair_ok.values())
    report = {
        "IDENTITY_TEST_PASS": "YES" if ok else "NO",
        "official_master": str(OUT_MASTER),
        "n_master_rows": len(master),
        "bootstrap_B": BOOT_B,
        "bootstrap_seed": BOOT_SEED,
        "no_new_science": True,
        "pair_identity": pair_ok,
        "pair_summaries": slim1,
        "pass1_sha256": {p: digest({k: pass1[p][k] for k in cmp_keys}) for p in FORMAL_PAIRS},
        "pass2_sha256": {p: digest({k: pass2[p][k] for k in cmp_keys}) for p in FORMAL_PAIRS},
    }
    OUT_ID.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "IDENTITY_TEST_PASS": report["IDENTITY_TEST_PASS"],
        "n_master_rows": len(master),
        "pair_n": {p: slim1[p]["n_official_complete_case"] for p in FORMAL_PAIRS},
        "pair_identity": pair_ok,
    }, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
