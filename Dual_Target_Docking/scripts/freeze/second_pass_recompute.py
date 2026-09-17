#!/usr/bin/env python3
"""Second-pass residual freeze: AChE max-vs-median + real theta grid + census.

Zero-dock. Does not copy panel pChEMBL into dump-gated max/median.
Writes only the two canonical tables that this pass is allowed to refresh,
plus docking_failure_census_v1.csv.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "data" / "jcim_novelty_v0" / "scripts"))

from analysis.bootstrap_metrics import assign_fourclass, assign_strict  # noqa: E402
from analysis.compute_canonical_results import (  # noqa: E402
    CANON,
    compute_label_sensitivity,
    compute_max_median,
    load_main,
    read_csv,
    write_csv,
)
from claim_hardening_v1 import docking_census  # noqa: E402

OUT = ROOT / "remediation_outputs" / "freeze_audit"
DUMP = (
    ROOT
    / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/eight_pair_dump_gated_v1"
    / "max_vs_median_ligand_v1.csv"
)
MASTER = CANON / "current_score_master.csv"
AUDIT = ROOT / "data/jcim_novelty_v0/tables/high_confidence_activity_audit_v1.csv"
MOLS_A = ROOT / "data/public_pair_selection/mols_ACHE.json"
MOLS_B = ROOT / "data/public_pair_selection/mols_BCHE.json"


def _json_has(path: Path, cid: str) -> bool:
    if not path.is_file():
        return False
    data = json.loads(path.read_text(encoding="utf-8"))
    return cid in data


def diagnose() -> dict:
    cid = "CHEMBL3960861"
    lig = "AB_056"
    master = [r for r in read_csv(MASTER) if r["pair"] == "AChE/BChE" and r["analysis_set"] == "main"]
    row = next((r for r in master if r["ligand_id"] == lig), None)
    dump = read_csv(DUMP)
    dump_ab = [r for r in dump if r["pair"] == "AChE/BChE" and r["ligand"] == lig]
    dump_cid = [r for r in dump if r.get("molecule_chembl_id") == cid]
    dump_ache = [r for r in dump if r["pair"] == "AChE/BChE"]
    audit_hits = []
    if AUDIT.is_file():
        audit_hits = [r for r in read_csv(AUDIT) if r.get("molecule_chembl_id") == cid]
    cc = [r for r in master if r["complete_case"] in ("1", 1, "True")]
    dump_ids = {r["ligand"] for r in dump_ache}
    missing_from_dump = sorted(r["ligand_id"] for r in cc if r["ligand_id"] not in dump_ids)

    disagreements = []
    counts = {
        "construction_class": Counter(),
        "primary_class_theta6": Counter(),
        "theta6_from_pA_pB": Counter(),
        "strict_from_pA_pB": Counter(),
        "theta_5.5": Counter(),
        "theta_6.5": Counter(),
    }
    for r in cc:
        pa, pb = r.get("pA"), r.get("pB")
        t6 = assign_fourclass(pa, pb, 6.0)
        t55 = assign_fourclass(pa, pb, 5.5)
        t65 = assign_fourclass(pa, pb, 6.5)
        st = assign_strict(pa, pb)
        counts["construction_class"][r["construction_class"]] += 1
        counts["primary_class_theta6"][r["primary_class_theta6"]] += 1
        counts["theta6_from_pA_pB"][t6] += 1
        counts["strict_from_pA_pB"][st] += 1
        counts["theta_5.5"][t55] += 1
        counts["theta_6.5"][t65] += 1
        rec = {
            "ligand_id": r["ligand_id"],
            "chembl": r["molecule_chembl_id"],
            "pA": r["pA"],
            "pB": r["pB"],
            "construction_class": r["construction_class"],
            "primary_class_theta6": r["primary_class_theta6"],
            "theta6_from_pA_pB": t6,
            "strict_from_pA_pB": st,
            "theta_5.5": t55,
            "theta_6.5": t65,
        }
        if r["construction_class"] != r["primary_class_theta6"] or r["primary_class_theta6"] != t6 or t6 != st:
            disagreements.append(rec)

    return {
        "AB_056_master": row,
        "in_dump_gated_ligand_table": bool(dump_ab or dump_cid),
        "dump_ache_n": len(dump_ache),
        "current_complete_case_n": len(cc),
        "missing_from_dump": missing_from_dump,
        "in_high_confidence_activity_audit": len(audit_hits),
        "in_mols_ACHE_json": _json_has(MOLS_A, cid),
        "in_mols_BCHE_json": _json_has(MOLS_B, cid),
        "decision": "EXCLUDE_FROM_MAX_MEDIAN",
        "decision_reason": (
            "No ChEMBL37 dump-gated assay rows (n_act_A/n_act_B) for AB_056/"
            "CHEMBL3960861. Panel pChEMBL 6.94/4.3 is construction provenance only "
            "and is not copied into max/median aggregation."
        ),
        "class_counts": {k: dict(v) for k, v in counts.items()},
        "n_disagreement_ligands": len(disagreements),
        "disagreements": disagreements,
        "construction_vs_theta6": sum(
            1 for r in cc if r["construction_class"] != r["primary_class_theta6"]
        ),
        "primary_vs_recomputed_theta6": sum(
            1 for r in cc if r["primary_class_theta6"] != assign_fourclass(r.get("pA"), r.get("pB"), 6.0)
        ),
        "theta6_vs_strict": sum(
            1
            for r in cc
            if assign_fourclass(r.get("pA"), r.get("pB"), 6.0)
            != assign_strict(r.get("pA"), r.get("pB"))
        ),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    diag = diagnose()
    (OUT / "second_pass_diagnostic.json").write_text(
        json.dumps(diag, indent=2, default=str) + "\n", encoding="utf-8"
    )
    print("AB_056 in dump-gated table:", diag["in_dump_gated_ligand_table"])
    print("missing from dump:", diag["missing_from_dump"])
    print("construction_vs_theta6:", diag["construction_vs_theta6"])
    print("primary_vs_recomputed_theta6:", diag["primary_vs_recomputed_theta6"])
    print("theta6_vs_strict:", diag["theta6_vs_strict"])
    print("class_counts:", json.dumps(diag["class_counts"]))

    packs = load_main()
    labels = compute_label_sensitivity(packs)
    maxmed = compute_max_median(packs)
    write_csv(CANON / "label_aggregation_sensitivity.csv", labels)
    write_csv(CANON / "max_vs_median_sensitivity.csv", maxmed)
    ache_lab = [r for r in labels if r["pair"] == "AChE/BChE"]
    ache_mm = [r for r in maxmed if r["pair"] == "AChE/BChE"]
    print("AChE label sensitivity:")
    for r in ache_lab:
        print(
            f"  {r['label_rule']}: {r['n_dual']}/{r['n_A_only']}/{r['n_B_only']}/"
            f"{r['n_neither']} gray={r['n_gray']} D/A={r['auroc_D_vs_A']} "
            f"D/B={r['auroc_D_vs_B']} smin={r['summary_min']} "
            f"CI=[{r['ci_lo']},{r['ci_hi']}] status={r['status']}"
        )
    print("AChE max vs median:")
    for r in ache_mm:
        print(
            f"  {r['aggregation']}: {r['n_dual']}/{r['n_A_only']}/{r['n_B_only']} "
            f"D/A={r['auroc_D_vs_A']} D/B={r['auroc_D_vs_B']} "
            f"smin={r['summary_min']} CI=[{r['ci_lo']},{r['ci_hi']}]"
        )

    census = docking_census()
    census_path = ROOT / "data/jcim_novelty_v0/tables/docking_failure_census_v1.csv"
    write_csv(census_path, census)
    ache_c = next(r for r in census if r["set"] == "main_panel" and r["pair"] == "AChE/BChE")
    print("AChE census:", ache_c["n_success_both_ends"], ache_c["failed_ligands"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
