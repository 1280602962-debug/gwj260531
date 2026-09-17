#!/usr/bin/env python3
"""Apply explicit activity-record exclusions and rebuild ligand-level aggregates.

Reads the committed high-confidence assay audit and the paper-adjudicated
exclusion table. Does not fetch ChEMBL. Missing arms are unresolved, never
filled as inactive.

Writes:
  data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv
  data/processed/activity_adjudication/excluded_activity_rows_applied_v1.csv
  data/processed/activity_adjudication/ligand_status_v1.csv
"""
from __future__ import annotations

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import assign_fourclass  # noqa: E402

AUDIT = ROOT / "data/jcim_novelty_v0/tables/high_confidence_activity_audit_v1.csv"
EXCLUDE = ROOT / "data/jcim_novelty_v0/tables/primary_only_excluded_activity_rows_v1.csv"
OUT = ROOT / "data/processed/activity_adjudication"

TARGET = {
    "EGFR/HER2": ("CHEMBL203", "CHEMBL1824"),
    "PIK3CA/mTOR": ("CHEMBL4005", "CHEMBL2842"),
    "AChE/BChE": ("CHEMBL220", "CHEMBL1914"),
}
GENE_TO_TID = {
    "EGFR": "CHEMBL203",
    "HER2": "CHEMBL1824",
    "PIK3CA": "CHEMBL4005",
    "mTOR": "CHEMBL2842",
    "AChE": "CHEMBL220",
    "BChE": "CHEMBL1914",
}
PANEL_MAP = {
    "EGFR/HER2": "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
    "PIK3CA/mTOR": "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
    "AChE/BChE": "data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
}


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
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def r4(x):
    if x is None:
        return ""
    return f"{float(x):.4f}".rstrip("0").rstrip(".") if float(x) == int(float(x)) else f"{float(x):.4f}"


def load_panels():
    lig_by_cid = {}
    hist = {}
    for pair, rel in PANEL_MAP.items():
        rows = read_csv(ROOT / rel)
        id_col = "panel_id" if "panel_id" in rows[0] else "ligand"
        pa_col = [c for c in rows[0] if c.startswith("pchembl_") and "HER2" not in c and "BCHE" not in c and "MTOR" not in c][0]
        pb_col = [c for c in rows[0] if c.startswith("pchembl_") and c != pa_col][0]
        # explicit
        if pair == "EGFR/HER2":
            pa_col, pb_col = "pchembl_EGFR", "pchembl_HER2"
        elif pair == "PIK3CA/mTOR":
            pa_col, pb_col = "pchembl_PIK3CA", "pchembl_MTOR"
        else:
            pa_col, pb_col = "pchembl_ACHE", "pchembl_BCHE"
        for r in rows:
            cid = r["molecule_chembl_id"]
            lig = r[id_col]
            lig_by_cid.setdefault(cid, []).append((pair, lig))
            hist[(pair, lig)] = {
                "molecule_chembl_id": cid,
                "historical_class": r.get("class", ""),
                "historical_pA": fnum(r.get(pa_col)),
                "historical_pB": fnum(r.get(pb_col)),
            }
    return lig_by_cid, hist


def main() -> int:
    lig_by_cid, hist = load_panels()
    exclusions = read_csv(EXCLUDE)
    ex_keys = set()
    for r in exclusions:
        tid = GENE_TO_TID.get(r["target"])
        if not tid:
            continue
        ex_keys.add((r["ligand"], tid, r["document_chembl_id"], f"{fnum(r['pchembl_value']):.4f}"))

    applied = []
    values = defaultdict(lambda: {"A": [], "B": []})
    n_audit = 0
    n_keep = 0
    for row in read_csv(AUDIT):
        if str(row.get("keep")) not in {"1", "True", "true"}:
            continue
        cid = row["molecule_chembl_id"]
        tid = row["target_chembl_id"]
        pv = fnum(row.get("pchembl_value"))
        if pv is None:
            continue
        n_audit += 1
        for pair, lig in lig_by_cid.get(cid, []):
            ta, tb = TARGET[pair]
            end = "A" if tid == ta else ("B" if tid == tb else None)
            if end is None:
                continue
            key = (lig, tid, row["document_chembl_id"], f"{pv:.4f}")
            if key in ex_keys:
                applied.append(
                    {
                        "pair": pair,
                        "ligand": lig,
                        "molecule_chembl_id": cid,
                        "target_chembl_id": tid,
                        "document_chembl_id": row["document_chembl_id"],
                        "assay_chembl_id": row["assay_chembl_id"],
                        "pchembl_value": pv,
                        "standard_type": row.get("standard_type", ""),
                        "assay_type": row.get("assay_type", ""),
                        "reason": next(
                            e["reason"]
                            for e in exclusions
                            if e["ligand"] == lig
                            and GENE_TO_TID.get(e["target"]) == tid
                            and e["document_chembl_id"] == row["document_chembl_id"]
                            and abs(fnum(e["pchembl_value"]) - pv) < 1e-6
                        ),
                    }
                )
                continue
            values[(pair, lig)][end].append(pv)
            n_keep += 1

    status_rows = []
    agg_rows = []
    for (pair, lig), rec in sorted(hist.items()):
        va, vb = values[(pair, lig)]["A"], values[(pair, lig)]["B"]
        max_a = max(va) if va else None
        max_b = max(vb) if vb else None
        med_a = statistics.median(va) if va else None
        med_b = statistics.median(vb) if vb else None
        nA, nB = len(va), len(vb)
        if nA == 0 and nB == 0:
            # No high-confidence audit rows (e.g. AB_056). Keep panel pChEMBL
            # for primary labels; do not invent dump max/median provenance.
            max_a, max_b = rec["historical_pA"], rec["historical_pB"]
            med_a = med_b = None
            cls = assign_fourclass(max_a, max_b)
            eligible = 1
            st = "panel_pchembl_no_audit_rows"
        elif max_a is None or max_b is None:
            cls = None
            eligible = 0
            st = "unresolved_missing_arm"
        else:
            cls = assign_fourclass(max_a, max_b)
            eligible = 1
            st = "both_arms_present"
        hist_cls = rec["historical_class"]
        status_rows.append(
            {
                "pair": pair,
                "ligand": lig,
                "molecule_chembl_id": rec["molecule_chembl_id"],
                "historical_class": hist_cls,
                "historical_pA": rec["historical_pA"],
                "historical_pB": rec["historical_pB"],
                "n_act_A": len(va),
                "n_act_B": len(vb),
                "max_A": max_a if max_a is not None else "",
                "median_A": med_a if med_a is not None else "",
                "max_B": max_b if max_b is not None else "",
                "median_B": med_b if med_b is not None else "",
                "class_max": cls or "",
                "class_median": assign_fourclass(med_a, med_b) if med_a is not None and med_b is not None else "",
                "activity_eligible": eligible,
                "activity_status": st,
                "class_flip_vs_historical": int(bool(cls) and cls != hist_cls),
                "pA_changed": int(max_a is not None and rec["historical_pA"] is not None and abs(max_a - rec["historical_pA"]) > 1e-6),
                "pB_changed": int(max_b is not None and rec["historical_pB"] is not None and abs(max_b - rec["historical_pB"]) > 1e-6),
            }
        )
        agg_rows.append(status_rows[-1])

    write_csv(OUT / "excluded_activity_rows_applied_v1.csv", applied)
    write_csv(OUT / "ligand_activity_aggregate_v1.csv", agg_rows)
    write_csv(OUT / "ligand_status_v1.csv", status_rows)
    flips = [r for r in status_rows if r["class_flip_vs_historical"] or r["activity_status"] != "both_arms_present" or r["pA_changed"] or r["pB_changed"]]
    print(f"audit keep rows seen={n_audit} kept_for_agg={n_keep} exclusions_applied={len(applied)}")
    print(f"ligands={len(status_rows)} changed_or_unresolved={len(flips)}")
    for r in flips:
        print(
            f"  {r['pair']} {r['ligand']}: hist {r['historical_class']} {r['historical_pA']}/{r['historical_pB']} "
            f"-> {r['class_max']} {r['max_A']}/{r['max_B']} eligible={r['activity_eligible']} {r['activity_status']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
