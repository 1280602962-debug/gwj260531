#!/usr/bin/env python3
"""Leave-one-document-group Δ on pairs with a deposited ligand–document map.

Does not harvest ChEMBL. JAK1/TYK2 is reported unresolved if no map is present.
Missing classes after a drop are not filled as 0.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import auroc  # noqa: E402

CANON = ROOT / "results" / "canonical"
MASTER = CANON / "current_score_master.csv"
GROUPS = ROOT / "data/jcim_novelty_v0/tables/document_blocked_ligand_groups_v1.csv"
PAIRS = ("EGFR/HER2", "JAK1/TYK2")


def fnum(v):
    if v is None or v == "":
        return None
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def r4(x) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return ""
    return f"{float(x):.4f}"


def load_eligible(pair: str) -> list[dict]:
    out = []
    for r in read_csv(MASTER):
        if r["pair"] != pair or r.get("analysis_set") != "main":
            continue
        if r.get("complete_case") not in ("1", "True"):
            continue
        if str(r.get("activity_eligible", "1")) not in ("1", "True"):
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        out.append(rec)
    return out


def directional(recs):
    dual = [r for r in recs if r["cls"] == "dual"]
    a = [r for r in recs if r["cls"] == "A_only"]
    b = [r for r in recs if r["cls"] == "B_only"]
    nei = [r for r in recs if r["cls"] == "neither"]
    da = auroc([r["score_B"] for r in dual], [r["score_B"] for r in a]) if dual and a else float("nan")
    db = auroc([r["score_A"] for r in dual], [r["score_A"] for r in b]) if dual and b else float("nan")
    dn = auroc([r["score_A"] for r in dual], [r["score_A"] for r in nei]) if dual and nei else float("nan")
    delta = (dn - db) if math.isfinite(dn) and math.isfinite(db) else float("nan")
    return {
        "n_dual": len(dual),
        "n_A_only": len(a),
        "n_B_only": len(b),
        "n_neither": len(nei),
        "auroc_D_vs_A": da,
        "auroc_D_vs_B": db,
        "auroc_D_vs_neither_pocketA": dn,
        "delta_neither_minus_B": delta,
    }


def main() -> int:
    groups = defaultdict(dict)
    if GROUPS.is_file():
        for r in read_csv(GROUPS):
            groups[r["pair"]][r["ligand"]] = r["group_id"]
    rows = []
    for pair in PAIRS:
        recs = load_eligible(pair)
        gmap = groups.get(pair, {})
        if not gmap:
            full = directional(recs)
            rows.append(
                {
                    "pair": pair,
                    "dropped_group": "",
                    "n_dropped": 0,
                    "status": "unresolved_mapping_unavailable",
                    "n_dual": full["n_dual"],
                    "n_A_only": full["n_A_only"],
                    "n_B_only": full["n_B_only"],
                    "n_neither": full["n_neither"],
                    "auroc_D_vs_A": r4(full["auroc_D_vs_A"]),
                    "auroc_D_vs_B": r4(full["auroc_D_vs_B"]),
                    "delta_neither_minus_B": r4(full["delta_neither_minus_B"]),
                    "delta_change": "",
                    "note": "no deposited ligand-document map; ChEMBL sqlite not used",
                }
            )
            continue
        full = directional(recs)
        by_g = defaultdict(list)
        for r in recs:
            by_g[gmap.get(r["ligand_id"], r["ligand_id"])].append(r["ligand_id"])
        for gid, members in sorted(by_g.items()):
            keep = [r for r in recs if r["ligand_id"] not in set(members)]
            stats = directional(keep)
            missing = [k for k in ("n_dual", "n_A_only", "n_B_only") if stats[k] < 1]
            if missing:
                status = "not_estimable_missing_class"
                dlt = ""
                change = ""
            else:
                status = "ok"
                dlt = r4(stats["delta_neither_minus_B"])
                change = r4(stats["delta_neither_minus_B"] - full["delta_neither_minus_B"])
            rows.append(
                {
                    "pair": pair,
                    "dropped_group": gid,
                    "n_dropped": len(members),
                    "status": status,
                    "n_dual": stats["n_dual"],
                    "n_A_only": stats["n_A_only"],
                    "n_B_only": stats["n_B_only"],
                    "n_neither": stats["n_neither"],
                    "auroc_D_vs_A": r4(stats["auroc_D_vs_A"]) if status == "ok" else "",
                    "auroc_D_vs_B": r4(stats["auroc_D_vs_B"]) if status == "ok" else "",
                    "delta_neither_minus_B": dlt,
                    "delta_change": change,
                    "note": "leave-one-document-group on current eligible scores; missing class not filled as 0",
                }
            )
    out = CANON / "leave_one_document_delta.csv"
    write_csv(out, rows)
    print("wrote", out.relative_to(ROOT), "n=", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
