#!/usr/bin/env python3
"""Read-only selection provenance for the seven frozen pairs.

Historical panels are read from commit a0ef2726. Current membership is read
from the V4.2 mapping. This script does not sample, relabel, or rewrite them.
"""
from __future__ import annotations

import csv
import io
import json
import random
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path("/home/gwj/q110_v42")
DTD = ROOT / "Dual_Target_Docking"
OUT = DTD / "results/current/provenance_audit_20261023"
FREEZE = "a0ef2726036e35f5bf2b78412f9a3ba43506e8e9"


def git_csv(path: str) -> list[dict]:
    text = subprocess.check_output(["git", "show", f"{FREEZE}:{path}"], cwd=ROOT, text=True)
    return list(csv.DictReader(io.StringIO(text)))


def write_csv(name: str, rows: list[dict]) -> None:
    with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def replay_ache_strict() -> set[str]:
    def load(name: str) -> dict[str, float]:
        text = subprocess.check_output(
            ["git", "show", f"{FREEZE}:Dual_Target_Docking/data/public_pair_selection/mols_{name}.json"],
            cwd=ROOT, text=True,
        )
        return {key: float(value) for key, value in json.loads(text).items()}

    ache_values, bche_values = load("ACHE"), load("BCHE")
    rows = []
    for chembl_id in sorted(set(ache_values) & set(bche_values)):
        pa, pb = ache_values[chembl_id], bche_values[chembl_id]
        if pa >= 6.5 and pb >= 6.5:
            label = "dual"
        elif pa >= 6.5 and pb <= 5.5:
            label = "A_only"
        elif pb >= 6.5 and pa <= 5.5:
            label = "B_only"
        elif pa <= 5.5 and pb <= 5.5:
            label = "neither"
        else:
            continue
        rows.append((chembl_id, label))
    rng = random.Random(20260729)
    picked = []
    caps = {}
    for label, need in (("dual", 28), ("A_only", 28), ("B_only", 28), ("neither", 16)):
        pool = [chembl_id for chembl_id, row_label in rows if row_label == label]
        rng.shuffle(pool)
        got = 0
        for chembl_id in pool:
            key = (label, chembl_id[:8])
            if caps.get(key, 0) >= 5:
                continue
            caps[key] = caps.get(key, 0) + 1
            picked.append(chembl_id)
            got += 1
            if got >= need:
                break
    return set(picked)


def main() -> None:
    mapping = list(csv.DictReader((DTD / "results/current/pair_ligand_mapping.csv").open()))
    cross = {
        (row["pair"], row["panel_id"]): row
        for row in csv.DictReader((DTD / "activity_label_source_audit_20261008/chembl37_dual_endpoint_crosswalk_805.csv").open())
    }
    current = [row for row in mapping if row["activity_eligible"] == "1" and row["parent_independent"] == "1"]
    by_pair = {}
    for row in current:
        by_pair.setdefault(row["pair"], []).append(row)

    provenance = []
    replay = []
    exceptions = []

    track = {
        "F2/F10": "Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panels/panel_F2_F10_v1.csv",
        "JAK1/JAK2": "Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_JAK2_v1.csv",
        "JAK1/TYK2": "Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_TYK2_v1.csv",
        "PPARG/PPARA": "Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARG_PPARA_v1.csv",
        "PPARA/PPARD": "Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARA_PPARD_v1.csv",
    }
    summary = {
        row["pair"]: row
        for row in git_csv("Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panel_summary_v1.csv")
    }
    for pair, path in track.items():
        hist = git_csv(path)
        hist_by = {row["molecule_chembl_id"]: row for row in hist}
        matched = 0
        class_mismatch = 0
        missing = []
        for row in by_pair[pair]:
            old = cross[(pair, row["panel_id"])]
            hist_row = hist_by.get(old["molecule_chembl_id"])
            if hist_row is None:
                missing.append(row["panel_id"])
                continue
            matched += 1
            if hist_row["theta6_class"] != row["class"]:
                class_mismatch += 1
        only_hist = sorted(set(hist_by) - {
            cross[(pair, row["panel_id"])]["molecule_chembl_id"] for row in by_pair[pair]
        })
        pool = summary[pair]
        provenance.append({
            "pair": pair,
            "historical_commit": FREEZE,
            "candidate_pool_source": "ChEMBL 37 dump census counts in track_b_panel_summary_v1.csv; molecule-level pool cache is gitignored and was not replayed",
            "activity_rule_at_sampling": "strict 6.5/5.5 construction; deposited theta6_class is the later primary label",
            "theta_final_label": "6.0",
            "random_seed": pool["seed"],
            "quota": pool["quota"],
            "pool_counts": f"theta dual/A/B/neither {pool['n_both_theta_dual']}/{pool['n_both_theta_A_only']}/{pool['n_both_theta_B_only']}/{pool['n_both_theta_neither']}; strict smallmol {pool['n_strict_dual_smallmol']}/{pool['n_strict_A_only_smallmol']}/{pool['n_strict_B_only_smallmol']}/{pool['n_strict_neither_smallmol']}",
            "chemical_filters_supported": "extract_track_b_panels_v1.py at d824eb6b records MOL, Small molecule, MW 150-750, heavy atoms 10-60 and no metal. It states there is no Bemis-Murcko cap.",
            "named_preference": "none recorded in the deposited panel",
            "historical_panel": path,
            "n_historical_panel": len(hist),
            "n_current_eligible": len(by_pair[pair]),
        })
        replay.append({
            "pair": pair,
            "audit_status": "DOCUMENTED_NOT_REPLAYED",
            "historical_commit": FREEZE,
            "sampler_program_found": "YES",
            "replay_executed": "NO",
            "n_historical": len(hist),
            "n_current_eligible": len(by_pair[pair]),
            "n_identity_match": matched,
            "n_identity_mismatch": len(only_hist) + len(missing),
            "n_theta6_class_mismatch": class_mismatch,
            "unresolved_reason": "extract_track_b_panels_v1.py at d824eb6b uses random.Random(20260729).shuffle by class. Its 2026-09-05 note says intersection ordering was sorted after the deposited panels and the script must not be rerun to replace them. The molecule-level candidate cache is not in git. Deposited identities and theta6 classes match the current 110.",
            "named_or_potency_preference": "NO",
            "further_repair": "NO",
        })
        if only_hist or missing or class_mismatch:
            exceptions.append({
                "pair": pair,
                "panel_id": ";".join(missing + only_hist),
                "kind": "UNEXPLAINED_PANEL_DIFFERENCE",
                "explanation": "Identity or theta6 class does not match the deposited Track B panel.",
                "affects_current_membership": "YES",
            })

    egfr = git_csv("Dual_Target_Docking/data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")
    egfr_by = {row["molecule_chembl_id"]: row for row in egfr}
    egfr_current = by_pair["EGFR/HER2"]
    matched = class_mismatch = 0
    for row in egfr_current:
        old = cross[("EGFR/HER2", row["panel_id"])]
        hist = egfr_by[old["molecule_chembl_id"]]
        matched += 1
        if hist["class"] != row["class"]:
            class_mismatch += 1
    excluded = next(row for row in egfr if row["panel_id"] == "EH120_059")
    role_counts = Counter(row["role_note"] for row in egfr)
    provenance.append({
        "pair": "EGFR/HER2",
        "historical_commit": FREEZE,
        "candidate_pool_source": "2026-07-23 ChEMBL REST panel recorded in FROZEN_PUBLIC_PAIRS.yaml; expanded panel deposited as panel_v0_120.csv",
        "activity_rule_at_sampling": "theta 6.0; panel40 retained and stage-1 expansion used hard-negative priority",
        "theta_final_label": "6.0",
        "random_seed": "20260727 in build_egfr_her2_panel120.py at d824eb6b",
        "quota": "deposited N=110 with classes 28/38/32/12 before the later unresolved-arm exclusion",
        "pool_counts": "yaml paired_both_measured 1751; dual 1182; A_only 207; B_only 46; dual_weak 316",
        "chemical_filters_supported": "Murcko maximum 5 per class and scaffold is stated in egfr_her2_panel120_v0/MANIFEST.md",
        "named_preference": "EH40_01 TAK-285 pose-gold anchor; EH40_02 lapatinib; one dual candidate recorded as replacing an imatinib label-noise candidate; 40 ligands retained from panel40; 70 labeled stage1 expand hardneg-priority",
        "historical_panel": "Dual_Target_Docking/data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        "n_historical_panel": len(egfr),
        "n_current_eligible": len(egfr_current),
    })
    replay.append({
        "pair": "EGFR/HER2",
        "audit_status": "DOCUMENTED_NOT_REPLAYED",
        "historical_commit": FREEZE,
        "sampler_program_found": "YES",
        "replay_executed": "NO",
        "n_historical": len(egfr),
        "n_current_eligible": len(egfr_current),
        "n_identity_match": matched,
        "n_identity_mismatch": 0,
        "n_theta6_class_mismatch": class_mismatch,
        "unresolved_reason": "build_egfr_her2_panel120.py at d824eb6b retains panel40 and fills quotas at theta 6.0, but molecule metadata can come from the ChEMBL API. This round did not call the API. The deposited 110 explain the current 109 after EH120_059 is marked activity_eligible=0.",
        "named_or_potency_preference": "YES",
        "further_repair": "NO",
    })
    exceptions.append({
        "pair": "EGFR/HER2",
        "panel_id": "EH120_059",
        "kind": "TECHNICAL_OR_IDENTITY_EXCLUSION",
        "explanation": f"{excluded['molecule_chembl_id']} / {excluded['inchi_key']} remains in the historical 110-member panel as A_only. LABEL_PROVENANCE.md marks it unresolved_missing_arm. Current mapping keeps activity_eligible=0, so the eligible count is 109. It was not removed because of a docking score.",
        "affects_current_membership": "NO",
    })
    exceptions.append({
        "pair": "EGFR/HER2",
        "panel_id": "EH40_01;EH40_02",
        "kind": "DOCUMENTED_NAMED_PREFERENCE",
        "explanation": "Role notes name TAK-285 as the pose-gold anchor and lapatinib as a classic dual inhibitor. " + "; ".join(f"{note}={count}" for note, count in role_counts.most_common()),
        "affects_current_membership": "NO",
    })

    ache = git_csv("Dual_Target_Docking/data/ache_bche_panel_v0/tables/panel_v0_strict.csv")
    if Counter(row["gray_excluded"] for row in ache) != Counter({"True": 100}):
        raise SystemExit("AChE gray_excluded column is no longer uniformly True")
    ache_by_panel = {row["panel_id"]: row for row in ache}
    ache_current_ids = {row["panel_id"] for row in by_pair["AChE/BChE"]}
    strict_versus_theta6 = 0
    for row in by_pair["AChE/BChE"]:
        if ache_by_panel[row["panel_id"]]["class"] != row["class"]:
            strict_versus_theta6 += 1
    provenance.append({
        "pair": "AChE/BChE",
        "historical_commit": FREEZE,
        "candidate_pool_source": "ache_bche_panel_v0 strict panel; earlier public mols_ACHE.json and mols_BCHE.json are candidate extracts, not the draw program",
        "activity_rule_at_sampling": "strict 6.5/5.5; gray zone excluded from the main panel by the manifest",
        "theta_final_label": "6.0",
        "random_seed": "20260729 in build_strict_panels.py at d824eb6b",
        "quota": "N=100; dual/A_only/B_only/neither 28/28/28/16",
        "pool_counts": "not independently reconstructed",
        "chemical_filters_supported": "ChEMBL canonical SMILES; salts desalted to the largest fragment at preparation, stated in MANIFEST.md",
        "named_preference": "AB_040 and AB_046 are the salt and free-base records of one parent; no separate named seed was found",
        "historical_panel": "Dual_Target_Docking/data/ache_bche_panel_v0/tables/panel_v0_strict.csv",
        "n_historical_panel": len(ache),
        "n_current_eligible": len(by_pair["AChE/BChE"]),
    })
    replay_ids = replay_ache_strict()
    deposited_ids = {row["molecule_chembl_id"] for row in ache}
    replay.append({
        "pair": "AChE/BChE",
        "audit_status": "DOCUMENTED_NOT_REPLAYED",
        "historical_commit": FREEZE,
        "sampler_program_found": "YES",
        "replay_executed": "YES",
        "n_historical": len(ache),
        "n_current_eligible": len(by_pair["AChE/BChE"]),
        "n_identity_match": len(ache_current_ids & set(ache_by_panel)),
        "n_identity_mismatch": 0,
        "n_theta6_class_mismatch": strict_versus_theta6,
        "unresolved_reason": (
            "Read-only replay of build_strict_panels.py logic with random.Random(20260729), strict 6.5/5.5, "
            "and the chembl-id eight-character cap returned 100 molecules. Intersection with the deposited panel is "
            f"{len(replay_ids & deposited_ids)}. Replay-only {sorted(replay_ids - deposited_ids)}; "
            f"deposited-only {sorted(deposited_ids - replay_ids)}. The script comment says SMILES were absent and the cap is not Bemis-Murcko. "
            "The seed was not changed and the draw was not repeated. Current 98 eligible observations remain inside the deposited 100 "
            "after AB_040 is aliased to AB_046 and AB_087 is activity_eligible=0. Strict class equals the current theta=6.0 class on those 98. "
            "gray_excluded=True is written on every emitted row by the script."
        ),
        "named_or_potency_preference": "NO",
        "further_repair": "NO",
    })
    exceptions.append({
        "pair": "AChE/BChE",
        "panel_id": "AB_040",
        "kind": "PARENT_ALIAS",
        "explanation": "Historical AB_040 is CHEMBL222904, the oxalate salt of the same structure as AB_046 CHEMBL254016. Current mapping sets AB_040 parent_independent=0 and canonical_ligand_id=AB_046. AB_046 is the eligible observation. The fold file still contains an AB_040 row; FORMAL_AUTHORITY_PATHS.yaml says it is not an independent observation.",
        "affects_current_membership": "NO",
    })
    exceptions.append({
        "pair": "AChE/BChE",
        "panel_id": "AB_056",
        "kind": "REPLAY_DIFFERS_FROM_DEPOSITED_PANEL",
        "explanation": "The read-only replay selects CHEMBL173607 where the deposited panel has AB_056 / CHEMBL3960861. Both are strict A_only. AB_056 remains in the current eligible set. The script's diversity key is the first eight characters of the ChEMBL ID, while its comment says a real Murcko cap needed SMILES. This does not remove AB_056 from the frozen evaluation set.",
        "affects_current_membership": "NO",
    })
    exceptions.append({
        "pair": "AChE/BChE",
        "panel_id": "AB_087",
        "kind": "TECHNICAL_OR_IDENTITY_EXCLUSION",
        "explanation": f"Historical {ache_by_panel['AB_087']['molecule_chembl_id']} is strict-class {ache_by_panel['AB_087']['class']}. LABEL_PROVENANCE.md marks AB_087 unresolved_missing_arm. Current activity_eligible=0. The eligible count is 98.",
        "affects_current_membership": "NO",
    })

    if sum(int(row["n_current_eligible"]) for row in replay) != 757:
        raise SystemExit("seven-pair eligible total is not 757")
    write_csv("SEVEN_PAIR_SELECTION_PROVENANCE.csv", provenance)
    write_csv("SEVEN_PAIR_EXACT_REPLAY_RESULTS.csv", replay)
    write_csv("SEVEN_PAIR_SELECTION_EXCEPTIONS.csv", exceptions)
    print("seven-pair rows", len(replay), "exceptions", len(exceptions), "theta6_or_strict_differences", {row["pair"]: row["n_theta6_class_mismatch"] for row in replay})


if __name__ == "__main__":
    main()
