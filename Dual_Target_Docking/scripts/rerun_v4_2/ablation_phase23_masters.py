#!/usr/bin/env python3
"""Phase 2-3: official 8-pair M0 seed42 master + rescoring/independent manifests. No AUROC."""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    MAP,
    PAIRS,
    PHASE7,
    QA,
    RECEPTORS,
    SEC,
    GIND,
    VINA_JOBS,
)


def vina_jobdir(pair: str, pdb: str, eid: str) -> Path:
    return VINA_JOBS / f"{pair.replace('/', '_')}__{pdb}__{eid}__seed42"


def main() -> int:
    QA.mkdir(parents=True, exist_ok=True)
    SEC.mkdir(parents=True, exist_ok=True)
    GIND.mkdir(parents=True, exist_ok=True)

    mapping = [
        r for r in csv.DictReader(MAP.open(encoding="utf-8-sig", newline=""))
        if r["pair"] in PAIRS and r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"
    ]
    s42 = {}
    for r in csv.DictReader(PHASE7.open(encoding="utf-8-sig", newline="")):
        if str(r.get("seed")) != "42" or r["pair"] not in PAIRS:
            continue
        s42[(r["pair"], r["pdb_id"], r["global_ligand_entity_id"])] = r

    m0 = []
    for r in mapping:
        pair = r["pair"]
        eid = r["global_ligand_entity_id"]
        for side, pdb in zip(("A", "B"), RECEPTORS[pair]):
            j = s42.get((pair, pdb, eid), {})
            aff = j.get("mode1_affinity", "")
            status = j.get("status", "MISSING_JOB")
            score = ""
            if status == "SUCCESS" and aff not in ("", None):
                score = f"{-float(aff):.6f}"
            m0.append({
                "pair": pair,
                "target_side": side,
                "pdb_id": pdb,
                "global_ligand_entity_id": eid,
                "canonical_ligand_id": r["canonical_ligand_id"],
                "pair_ligand_id": r.get("panel_id", r["canonical_ligand_id"]),
                "class": r["class"],
                "activity_eligible": r.get("activity_eligible", ""),
                "analysis_set": r.get("analysis_set", ""),
                "parent_independent": r.get("parent_independent", ""),
                "vina_status": status,
                "vina_mode1_affinity": aff,
                "vina_score": score,
                "TIMEOUT_treated_as": "missing_not_imputed",
            })
    m0_path = QA / "official_primary_seed42_score_master_8pair.csv"
    with m0_path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(m0[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(m0)

    phys = {}
    for r in mapping:
        eid = r["global_ligand_entity_id"]
        for pdb in RECEPTORS[r["pair"]]:
            key = (pdb, eid)
            rec = phys.setdefault(key, {
                "pdb_id": pdb,
                "global_ligand_entity_id": eid,
                "pairs": set(),
                "canonical_ligand_ids": set(),
                "vina_status": "",
                "pose_path": "",
            })
            rec["pairs"].add(r["pair"])
            rec["canonical_ligand_ids"].add(r["canonical_ligand_id"])
            j = s42.get((r["pair"], pdb, eid), {})
            if j.get("status"):
                rec["vina_status"] = j["status"]
            pose = vina_jobdir(r["pair"], pdb, eid) / "out.pdbqt"
            if pose.is_file():
                rec["pose_path"] = str(pose)

    # prefer a pose path from any pair sharing the physical key
    for (pdb, eid), rec in phys.items():
        if rec["pose_path"]:
            continue
        for pair in rec["pairs"]:
            p = vina_jobdir(pair, pdb, eid) / "out.pdbqt"
            if p.is_file():
                rec["pose_path"] = str(p)
                break

    rescore_rows = []
    gnina_rows = []
    for (pdb, eid), rec in sorted(phys.items()):
        vina_ok = rec["vina_status"] == "SUCCESS" and bool(rec["pose_path"])
        rescore_rows.append({
            "pdb_id": pdb,
            "global_ligand_entity_id": eid,
            "pairs": ";".join(sorted(rec["pairs"])),
            "canonical_ligand_ids": ";".join(sorted(rec["canonical_ligand_ids"])),
            "vina_status": rec["vina_status"] or "MISSING_JOB",
            "pose_path": rec["pose_path"],
            "rescoring_eligible": int(vina_ok),
            "rescoring_status_pre": "ELIGIBLE" if vina_ok else "MISSING_NO_VINA_POSE",
        })
        gnina_rows.append({
            "pdb_id": pdb,
            "global_ligand_entity_id": eid,
            "pairs": ";".join(sorted(rec["pairs"])),
            "canonical_ligand_ids": ";".join(sorted(rec["canonical_ligand_ids"])),
            "seed": 42,
            "vina_status": rec["vina_status"] or "MISSING_JOB",
            "run_even_if_vina_timeout": 1,
        })

    with (SEC / "rescoring_job_manifest.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rescore_rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rescore_rows)
    with (GIND / "gnina_seed42_job_manifest.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(gnina_rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(gnina_rows)

    n_elig = sum(r["rescoring_eligible"] for r in rescore_rows)
    n_miss = sum(1 for r in rescore_rows if r["rescoring_status_pre"] == "MISSING_NO_VINA_POSE")
    shared = sum(1 for r in rescore_rows if ";" in r["pairs"])
    report = {
        "pair_expanded_rows": len(m0),
        "unique_physical_jobs": len(rescore_rows),
        "rescoring_eligible": n_elig,
        "MISSING_NO_VINA_POSE": n_miss,
        "shared_receptor_parent_jobs": shared,
        "gnina_independent_jobs": len(gnina_rows),
        "m0_path": str(m0_path),
        "AUROC": "NOT_COMPUTED",
    }
    (QA / "ablation_phase23_manifest_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if len(rescore_rows) != 1592 or len(m0) != 1614:
        print("MANIFEST_COUNT_HOLD")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
