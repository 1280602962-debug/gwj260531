#!/usr/bin/env python3
"""Audit candidate dual-target pairs against Dual-VSDS public-pair hard gates.

Reads/writes under Dual_Target_Docking/data/public_pair_selection/.
Re-run with --refetch to refresh ChEMBL molecule maps (slow).
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "public_pair_selection"

STANDARD_OK = {"IC50", "Ki", "Kd", "EC50", "IC50app", "Ki app", "Potency"}

PAIRS = [
    ("EGFR", "HER2"),
    ("MCL1", "BCL2L1_BclxL"),
    ("PIK3CA", "MTOR"),
    ("BRD4", "HDAC1"),
    ("BRD4", "HDAC6"),
    ("PARP1", "MET"),
    ("ACHE", "BCHE"),
    ("JAK2", "HDAC1"),
    ("CDK6", "BRD4"),
    ("AKT1", "RPS6KB1_p70S6K"),
    ("MCL1", "BCL2"),
    ("PIK3CA", "PIK3CB"),
]


def fetch_activities(chembl_id: str, max_pages: int = 200):
    best = {}
    n_act = 0
    offset = 0
    limit = 1000
    while True:
        url = (
            f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
            f"?target_chembl_id={chembl_id}"
            f"&pchembl_value__isnull=false"
            f"&limit={limit}&offset={offset}"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "dual-target-pair-audit/0.1"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = json.load(resp)
        acts = data.get("activities") or []
        if not acts:
            break
        for a in acts:
            st = a.get("standard_type") or ""
            if st not in STANDARD_OK:
                continue
            mol = a.get("molecule_chembl_id")
            pv = a.get("pchembl_value")
            if mol is None or pv is None:
                continue
            pv = float(pv)
            n_act += 1
            if mol not in best or pv > best[mol]:
                best[mol] = pv
        page_meta = data.get("page_meta") or {}
        total = page_meta.get("total_count") or 0
        offset += limit
        if offset >= total or len(acts) < limit:
            break
        if offset // limit >= max_pages:
            break
        time.sleep(0.05)
    return best, n_act


def classify_pair(map_a, map_b, act: float = 6.0):
    both = set(map_a) & set(map_b)
    dual = a_only = b_only = dual_weak = 0
    for m in both:
        pa, pb = map_a[m], map_b[m]
        a_act = pa >= act
        b_act = pb >= act
        if a_act and b_act:
            dual += 1
        elif a_act and not b_act:
            a_only += 1
        elif b_act and not a_act:
            b_only += 1
        else:
            dual_weak += 1
    return {
        "n_A": len(map_a),
        "n_B": len(map_b),
        "n_paired_measured_both": len(both),
        "dual": dual,
        "A_only": a_only,
        "B_only": b_only,
        "dual_weak": dual_weak,
        "hard_gate_nA_ge_200": len(map_a) >= 200,
        "hard_gate_nB_ge_200": len(map_b) >= 200,
        "hard_gate_has_dual": dual >= 10,
        "hard_gate_has_A_only": a_only >= 10,
        "hard_gate_has_B_only": b_only >= 10,
        "hard_gate_three_classes": dual >= 10 and a_only >= 10 and b_only >= 10,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refetch", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tid = {t["name"]: t for t in json.load(open(OUT / "chembl_target_ids.json"))}
    molmaps = {}
    for name, meta in tid.items():
        cache = OUT / f"mols_{name}.json"
        if cache.exists() and not args.refetch:
            molmaps[name] = {k: float(v) for k, v in json.load(open(cache)).items()}
            continue
        print(f"fetching {name} ...", flush=True)
        best, n_act = fetch_activities(meta["chembl"])
        molmaps[name] = best
        json.dump(best, open(cache, "w"))
        print(f"  {len(best)} mols / {n_act} rows", flush=True)

    rows = []
    for a, b in PAIRS:
        c = classify_pair(molmaps[a], molmaps[b])
        c.update(
            {
                "pair": f"{a}/{b}",
                "chembl_A": tid[a]["chembl"],
                "chembl_B": tid[b]["chembl"],
                "uniprot_A": tid[a]["uniprot"],
                "uniprot_B": tid[b]["uniprot"],
            }
        )
        rows.append(c)
        print(
            f"{c['pair']:28s} paired={c['n_paired_measured_both']:5d} "
            f"dual={c['dual']:4d} A={c['A_only']:4d} B={c['B_only']:4d} "
            f"3class={c['hard_gate_three_classes']}"
        )

    keys = [
        "pair",
        "uniprot_A",
        "uniprot_B",
        "chembl_A",
        "chembl_B",
        "n_A",
        "n_B",
        "n_paired_measured_both",
        "dual",
        "A_only",
        "B_only",
        "dual_weak",
        "hard_gate_nA_ge_200",
        "hard_gate_nB_ge_200",
        "hard_gate_has_dual",
        "hard_gate_has_A_only",
        "hard_gate_has_B_only",
        "hard_gate_three_classes",
    ]
    with open(OUT / "chembl_pair_fourclass.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in keys})


if __name__ == "__main__":
    main()
