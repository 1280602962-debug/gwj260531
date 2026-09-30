#!/usr/bin/env python3
"""M0/M2 lock, identity/reuse, activity eligibility. No AUROC."""
from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    GIND,
    MAP,
    PHASE7,
    PROTO,
    QA,
    RECEPTORS,
    REC,
    RTM_CKPT,
    RUN,
    SEC,
    UNIQUE_14,
)
from rerun_v4_2.pre_auroc_v2_lib import finite_float, split_models  # noqa: E402

M0 = QA / "official_primary_seed42_score_master_8pair.csv"
M2M = SEC / "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv"
M2L = SEC / "RTMSCORE_VINA_POSE_RESCORE_POSE_LONG.csv"
MANIFEST = SEC / "rescoring_job_manifest.csv"
EXPECTED_TIMEOUTS = {
    ("AB_001", "4EY7"), ("AB_001", "4BDS"),
    ("AB_053", "4EY7"), ("AB_053", "4BDS"),
    ("AB_054", "4EY7"), ("AB_054", "4BDS"),
}


def load(path: Path) -> list[dict]:
    return list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))


def m0_lock() -> dict:
    rows = load(M0)
    st = Counter(r["vina_status"] for r in rows)
    timeouts = {(r["canonical_ligand_id"], r["pdb_id"]) for r in rows if r["vina_status"] == "TIMEOUT"}
    score_ok = True
    bad = []
    for r in rows:
        if r["vina_status"] != "SUCCESS":
            if r.get("vina_score") not in ("", None):
                bad.append((r["pair"], r["pdb_id"], r["global_ligand_entity_id"], "timeout_has_score"))
            continue
        try:
            aff = float(r["vina_mode1_affinity"])
            sc = float(r["vina_score"])
            if abs(sc - (-aff)) > 1e-6:
                score_ok = False
                bad.append((r["pair"], r["pdb_id"], r["global_ligand_entity_id"], "score_not_neg_affinity"))
        except (TypeError, ValueError):
            score_ok = False
            bad.append((r["pair"], r["pdb_id"], r["global_ligand_entity_id"], "nonfinite"))
    return {
        "level": "pair-expanded",
        "rows": len(rows),
        "SUCCESS": st.get("SUCCESS", 0),
        "TIMEOUT": st.get("TIMEOUT", 0),
        "timeouts": sorted(list(timeouts)),
        "timeouts_exact_expected": timeouts == EXPECTED_TIMEOUTS,
        "score_is_neg_vina_mode1_affinity": score_ok and not bad,
        "imputation_column": all(r.get("TIMEOUT_treated_as") == "missing_not_imputed" for r in rows),
        "parent_independent_all_1": all(r.get("parent_independent") == "1" for r in rows),
        "bad_examples": bad[:10],
        "M0_LOCKED": (
            len(rows) == 1614 and st.get("SUCCESS") == 1608 and st.get("TIMEOUT") == 6
            and timeouts == EXPECTED_TIMEOUTS and score_ok and not bad
        ),
    }


def m2_lock() -> dict:
    master = load(M2M)
    long = load(M2L)
    man = load(MANIFEST)
    st = Counter(r["job_status"] for r in master)
    n_pose = Counter()
    for r in man:
        if r.get("rescoring_eligible") == "1" and r.get("pose_path"):
            n_pose[len(split_models(Path(r["pose_path"]).read_text(errors="replace")))] += 1
    expected_long = n_pose.get(9, 0) * 9 + n_pose.get(8, 0) * 8
    pockets = list((SEC / "rtmscore" / "pockets").glob("*_rtmscore_pocket10.pdb"))
    pocket_ok = {p.name.split("_")[0] for p in pockets}
    mode1_mismatch = []
    long_m1 = {
        (r["pdb_id"], r["global_ligand_entity_id"]): r["RTMScore"]
        for r in long if str(r.get("vina_mode")) == "1" or str(r.get("is_vina_mode1")) == "1"
    }
    for r in master:
        if r.get("job_status") != "SUCCESS":
            continue
        key = (r["pdb_id"], r["global_ligand_entity_id"])
        a, b = r.get("M2_RTMScore", ""), long_m1.get(key, "")
        try:
            if abs(float(a) - float(b)) > 1e-6:
                mode1_mismatch.append(key)
        except (TypeError, ValueError):
            mode1_mismatch.append(key)
    ckpt_ok = RTM_CKPT.name == "rtmscore_model1.pth" and RTM_CKPT.is_file()
    return {
        "level": "physical",
        "master_rows": len(master),
        "SUCCESS": st.get("SUCCESS", 0),
        "MISSING_NO_VINA_POSE": st.get("MISSING_NO_VINA_POSE", 0),
        "status_counts": dict(st),
        "long_rows": len(long),
        "expected_long_rows": expected_long,
        "pose_job_hist": dict(n_pose),
        "pockets_present": sorted(pocket_ok),
        "n_pockets": len(pockets),
        "all_14_pockets": set(UNIQUE_14) <= pocket_ok,
        "checkpoint": str(RTM_CKPT),
        "checkpoint_is_model1": ckpt_ok,
        "mode1_master_matches_long": len(mode1_mismatch) == 0,
        "mode1_mismatch_n": len(mode1_mismatch),
        "M2_LOCKED": (
            len(master) == 1592 and st.get("SUCCESS") == 1586 and len(long) == 14250
            and expected_long == 14250 and ckpt_ok and set(UNIQUE_14) <= pocket_ok
            and len(mode1_mismatch) == 0
        ),
    }


def identity() -> dict:
    mapping = load(MAP)
    in_pairs = [r for r in mapping if r["pair"] in RECEPTORS]
    indep = [r for r in in_pairs if r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"]
    all_parents = {r["global_ligand_entity_id"] for r in in_pairs}
    indep_parents = {r["global_ligand_entity_id"] for r in indep}
    ab040 = [r for r in in_pairs if r.get("canonical_ligand_id") == "AB_040" or r.get("panel_id") == "AB_040"]
    ab046 = [r for r in in_pairs if r.get("canonical_ligand_id") == "AB_046" or r.get("panel_id") == "AB_046"]
    m0 = load(M0)
    m0_ab040 = [r for r in m0 if r.get("canonical_ligand_id") == "AB_040" or r.get("pair_ligand_id") == "AB_040"]
    # shared physical jobs
    phys = defaultdict(set)
    for r in m0:
        phys[(r["pdb_id"], r["global_ligand_entity_id"])].add(r["pair"])
    shared = {k: sorted(v) for k, v in phys.items() if len(v) > 1}
    # reuse: same physical key must have same vina_score across pairs
    reuse_bad = []
    by_phys = defaultdict(list)
    for r in m0:
        by_phys[(r["pdb_id"], r["global_ligand_entity_id"])].append(r)
    for k, rs in by_phys.items():
        scores = {r.get("vina_score") for r in rs}
        statuses = {r.get("vina_status") for r in rs}
        if len(statuses) > 1 or len(scores) > 1:
            reuse_bad.append((k, list(statuses), list(scores)))
    return {
        "mapping_rows_in_8_pairs": len(in_pairs),
        "independent_pair_level_rows": len(indep),
        "independent_pair_level_parents": len({(r["pair"], r["global_ligand_entity_id"]) for r in indep}),
        "global_unique_parents_independent": len(indep_parents),
        "pair_expanded_m0": len(m0),
        "unique_physical_jobs": len(phys),
        "shared_physical_jobs": len(shared),
        "shared_6N7A": sum(1 for (pdb, _) in shared if pdb == "6N7A"),
        "shared_6LXA": sum(1 for (pdb, _) in shared if pdb == "6LXA"),
        "AB_040_mapping_rows": len(ab040),
        "AB_040_parent_independent_values": sorted({r.get("parent_independent") for r in ab040}),
        "AB_046_mapping_rows": len(ab046),
        "AB_040_in_m0_official": len(m0_ab040),
        "AB_040_not_independent_observation": len(m0_ab040) == 0 and all(r.get("parent_independent") != "1" for r in ab040),
        "m0_reuse_conflicts": len(reuse_bad),
        "reuse_conflict_examples": reuse_bad[:5],
        "expected": {
            "historical_records": 808,
            "independent_pair_level_parents": 807,
            "global_unique_parents": 785,
            "pair_expanded": 1614,
            "unique_physical": 1592,
        },
    }


def activity() -> dict:
    mapping = [r for r in load(MAP) if r["pair"] in RECEPTORS]
    m0 = load(M0)
    zero = [r for r in mapping if r.get("activity_eligible") == "0"]
    m0_zero = [r for r in m0 if r.get("activity_eligible") == "0"]
    # missing activity must not be class neither
    bad_neither = [
        r for r in m0
        if r.get("activity_eligible") == "0" and r.get("class") == "neither"
    ]
    watch = [r for r in mapping if r.get("canonical_ligand_id") in {"EH120_059", "AB_087"} or r.get("panel_id") in {"EH120_059", "AB_087"}]
    return {
        "mapping_activity_eligible_0": len(zero),
        "m0_activity_eligible_0": len(m0_zero),
        "eligible0_mapped_to_neither": len(bad_neither),
        "EH120_059_AB_087": [
            {k: r.get(k) for k in ("pair", "panel_id", "canonical_ligand_id", "class", "activity_eligible", "parent_independent")}
            for r in watch
        ],
        "fourclass_must_exclude_eligible0": True,
        "GATE_FAIL_if_eligible0_became_neither": len(bad_neither) > 0,
    }


def main() -> int:
    out = {
        "M0": m0_lock(),
        "M2": m2_lock(),
        "identity": identity(),
        "activity": activity(),
        "AUROC": "NOT_COMPUTED",
    }
    (QA / "PRE_AUROC_V2_LOCK_IDENTITY.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({
        "M0_LOCKED": out["M0"]["M0_LOCKED"],
        "M2_LOCKED": out["M2"]["M2_LOCKED"],
        "AB_040_ok": out["identity"]["AB_040_not_independent_observation"],
        "m0_reuse_conflicts": out["identity"]["m0_reuse_conflicts"],
        "eligible0_neither": out["activity"]["eligible0_mapped_to_neither"],
        "identity_counts": {
            "indep_parents": out["identity"]["independent_pair_level_parents"],
            "global_unique": out["identity"]["global_unique_parents_independent"],
            "pair_expanded": out["identity"]["pair_expanded_m0"],
            "physical": out["identity"]["unique_physical_jobs"],
            "shared": out["identity"]["shared_physical_jobs"],
        },
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
