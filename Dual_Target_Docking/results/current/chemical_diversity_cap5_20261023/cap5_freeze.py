#!/usr/bin/env python3
"""Freeze the cap-5 scaffold subset before any sensitivity statistic.

Membership selection reads identity and canonical SMILES only.
It does not read docking scores, AUROC, model probabilities, or potency.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, __version__ as rdkit_version

ROOT = Path("/home/gwj/q110_v42")
DTD = ROOT / "Dual_Target_Docking"
OUT = DTD / "results/current/chemical_diversity_cap5_20261023"
SCRIPTS = DTD / "scripts/rerun_v4_2"
import sys

sys.path.insert(0, str(SCRIPTS))
from jcim_stage2_compute import _murcko  # noqa: E402

PREFIX = "BM5_SENSITIVITY_20261023_V1"
SCRIPT_VERSION = "cap5_freeze.py/1"
NULL_AUROC = 0.5
NULL_DELTA = 0.0
POINT_BAND = 0.02


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def connectivity_key(smiles: str) -> str:
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return ""
    return Chem.MolToInchiKey(mol).split("-")[0]


def selection_hash(pair: str, class_label: str, entity: str, panel_id: str) -> str:
    text = PREFIX + "|" + pair + "|" + class_label + "|" + entity + "|" + panel_id
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    head = __import__("subprocess").check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    pop_path = DTD / "results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
    map_path = DTD / "results/current/pair_ligand_mapping.csv"
    fold_path = DTD / "results/canonical/model_fold_assignments.csv"
    primary_path = DTD / "results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv"
    delta_path = DTD / "results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv"
    q110_path = DTD / "extension_jcim_validation/pik3ca_mtor_q110_v42/PIK3CA_MTOR_Q110_MEMBERSHIP.csv"
    authority_path = DTD / "CURRENT_RESULT_AUTHORITIES.json"
    population = read_csv(pop_path)
    mapping = read_csv(map_path)
    folds = read_csv(fold_path)
    eligible = [row for row in population if row.get("activity_eligible") == "1"]
    classes = {"dual": 0, "A_only": 0, "B_only": 0, "neither": 0}
    for row in eligible:
        classes[row["class"]] = classes.get(row["class"], 0) + 1
    q110 = [row for row in eligible if row["pair"] == "PIK3CA/mTOR"]
    q110_classes = {name: sum(row["class"] == name for row in q110) for name in classes}
    acceptance = []

    def add(check, ok, observed, expected):
        acceptance.append({
            "check": check,
            "result": "PASS" if ok else "FAIL",
            "observed": observed,
            "expected": expected,
        })

    add("population_rows", len(population) == 869, str(len(population)), "869")
    add("eligible_rows", len(eligible) == 867, str(len(eligible)), "867")
    add("class_dual", classes.get("dual") == 248, str(classes.get("dual")), "248")
    add("class_A_only", classes.get("A_only") == 256, str(classes.get("A_only")), "256")
    add("class_B_only", classes.get("B_only") == 252, str(classes.get("B_only")), "252")
    add("class_neither", classes.get("neither") == 111, str(classes.get("neither")), "111")
    add("q110_n", len(q110) == 110, str(len(q110)), "110")
    add("q110_classes", q110_classes == {"dual": 32, "A_only": 32, "B_only": 32, "neither": 14}, str(q110_classes), "32/32/32/14")
    add("unique_eligible_keys", len({(r["pair"], r["panel_id"]) for r in eligible}) == 867, "867", "867")
    add("ab040_not_eligible_observation", not any(r["panel_id"] == "AB_040" and r["activity_eligible"] == "1" and r.get("parent_independent", "1") == "1" for r in population), "checked", "alias excluded")
    authority = json.loads(authority_path.read_text(encoding="utf-8"))
    add("authority_population", authority.get("population", "").endswith("PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"), authority.get("population", ""), "results/current population")
    add("authority_q110", authority.get("pair_override", {}).get("PIK3CA/mTOR") == "extension_jcim_validation/pik3ca_mtor_q110_v42", str(authority.get("pair_override")), "Q110 override")
    if any(row["result"] == "FAIL" for row in acceptance):
        write_csv(OUT / "INPUT_ACCEPTANCE.csv", acceptance, ["check", "result", "observed", "expected"])
        raise SystemExit("INPUT_BLOCKED")

    smiles = {
        (row["pair"], row["canonical_ligand_id"]): row.get("canonical_smiles") or ""
        for row in mapping
        if str(row.get("parent_independent", "1")) not in {"0", "false", "False"}
    }
    import numpy
    import sklearn
    manifest = {
        "head": head,
        "remote_checked": "196639e9c5c407c3656f93eee0d14fb26be1b784",
        "script_version": SCRIPT_VERSION,
        "rdkit": rdkit_version,
        "numpy": numpy.__version__,
        "sklearn": sklearn.__version__,
        "python": platform.python_version(),
        "files": {
            "population": {"path": str(pop_path.relative_to(ROOT)), "sha256": sha256_file(pop_path)},
            "mapping": {"path": str(map_path.relative_to(ROOT)), "sha256": sha256_file(map_path)},
            "folds": {"path": str(fold_path.relative_to(ROOT)), "sha256": sha256_file(fold_path)},
            "primary": {"path": str(primary_path.relative_to(ROOT)), "sha256": sha256_file(primary_path)},
            "primary_deltas": {"path": str(delta_path.relative_to(ROOT)), "sha256": sha256_file(delta_path)},
            "q110_membership": {"path": str(q110_path.relative_to(ROOT)), "sha256": sha256_file(q110_path)},
            "stage2_B": {"path": "Dual_Target_Docking/results/jcim_stage2/module_B_negative_class.csv", "sha256": sha256_file(DTD / "results/jcim_stage2/module_B_negative_class.csv")},
            "stage2_C": {"path": "Dual_Target_Docking/results/jcim_stage2/module_C_wrong_pocket.csv", "sha256": sha256_file(DTD / "results/jcim_stage2/module_C_wrong_pocket.csv")},
            "stage2_D_units": {"path": "Dual_Target_Docking/results/jcim_stage2/module_D_chemistry_units.csv", "sha256": sha256_file(DTD / "results/jcim_stage2/module_D_chemistry_units.csv")},
            "stage2_D_oof": {"path": "Dual_Target_Docking/results/jcim_stage2/module_D_chemistry_oof_ligands.csv", "sha256": sha256_file(DTD / "results/jcim_stage2/module_D_chemistry_oof_ligands.csv")},
            "stage2_I_desc": {"path": "Dual_Target_Docking/results/jcim_stage2/module_I_descriptive.csv", "sha256": sha256_file(DTD / "results/jcim_stage2/module_I_descriptive.csv")},
            "stage2_I_lig": {"path": "Dual_Target_Docking/results/jcim_stage2/module_I_ligand_ranks.csv", "sha256": sha256_file(DTD / "results/jcim_stage2/module_I_ligand_ranks.csv")},
        },
        "murcko": "jcim_stage2_compute._murcko = MolToSmiles(GetScaffoldForMol)",
        "descriptive_thresholds_frozen_before_sensitivity": {
            "point_band": POINT_BAND,
            "auroc_null": NULL_AUROC,
            "delta_null": NULL_DELTA,
            "point_change_class": ["STABLE", "ATTENUATED", "AMPLIFIED", "REVERSED", "NEAR_NULL", "NOT_ESTIMABLE"],
            "role": "descriptive labels only; not a significance or equivalence test",
        },
    }
    (OUT / "INPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_csv(OUT / "INPUT_ACCEPTANCE.csv", acceptance, ["check", "result", "observed", "expected"])

    assigned = []
    blocking = []
    for row in eligible:
        key = (row["pair"], row["canonical_ligand_id"])
        smi = smiles.get(key, "")
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            blocking.append((row["pair"], row["panel_id"], "SMILES_PARSE_FAILURE"))
            scaffold, rings = "", -1
        else:
            scaffold = _murcko(mol)
            rings = mol.GetRingInfo().NumRings()
            if scaffold == "" and rings != 0:
                blocking.append((row["pair"], row["panel_id"], "EMPTY_SCAFFOLD_WITH_RINGS"))
            if scaffold != "" and rings == 0:
                blocking.append((row["pair"], row["panel_id"], "NONEMPTY_SCAFFOLD_WITHOUT_RINGS"))
        assigned.append({
            "pair": row["pair"],
            "panel_id": row["panel_id"],
            "global_ligand_entity_id": row["global_ligand_entity_id"],
            "canonical_ligand_id": row["canonical_ligand_id"],
            "class": row["class"],
            "canonical_smiles": smi,
            "murcko_scaffold": scaffold,
            "scaffold_is_acyclic": "YES" if scaffold == "" and rings == 0 else "NO",
            "ring_count": rings,
        })
    if blocking:
        write_csv(OUT / "SCAFFOLD_ASSIGNMENTS_867.csv", assigned, list(assigned[0]))
        raise SystemExit("SCAFFOLD_PARSE_BLOCKED " + str(blocking[:10]))
    write_csv(OUT / "SCAFFOLD_ASSIGNMENTS_867.csv", assigned, list(assigned[0]))
    by_key = {(row["pair"], row["canonical_ligand_id"]): row for row in assigned}

    concordance = []
    for fold in folds:
        ligand = fold["ligand_id"]
        obs = by_key.get((fold["pair"], ligand))
        fold_s = fold.get("scaffold") or ""
        if obs is None:
            mapping_row = next((row for row in mapping if row["pair"] == fold["pair"] and row["panel_id"] == ligand), None)
            kind = "FOLD_LIGAND_NOT_IN_ELIGIBLE_867"
            if mapping_row and str(mapping_row.get("parent_independent", "1")) in {"0", "false", "False"}:
                kind = "ALIAS_FOLD_ROW_NOT_AN_OBSERVATION"
            concordance.append({
                "pair": fold["pair"], "arm": fold["arm"], "ligand_id": ligand,
                "class": fold.get("class", ""),
                "fold_scaffold": fold_s,
                "computed_scaffold": "",
                "relation": kind,
                "blocks_freeze": "NO",
                "note": "Fold coverage is not the full four-class set. Alias rows are not eligible observations.",
            })
            continue
        computed = obs["murcko_scaffold"]
        fold_token = fold_s.strip().lower()
        fold_means_empty = fold_token in {"", "acyclic"}
        computed_empty = computed == ""
        if computed == fold_s:
            relation, blocks = "EXACT_STRING_MATCH", "NO"
        elif computed_empty and fold_means_empty and obs["ring_count"] == 0:
            relation, blocks = "REPRESENTATION_DIFFERENCE", "NO"
        else:
            left = connectivity_key(computed)
            right = "" if fold_means_empty else connectivity_key(fold_s)
            if left and left == right:
                relation, blocks = "REPRESENTATION_DIFFERENCE", "NO"
            elif computed_empty != (right == "" and fold_means_empty):
                relation, blocks = "ACYCLIC_BOUNDARY_CONFLICT", "YES"
            else:
                relation, blocks = "CONNECTIVITY_DIFFERENCE", "YES"
        concordance.append({
            "pair": fold["pair"], "arm": fold["arm"], "ligand_id": ligand,
            "class": obs["class"],
            "fold_scaffold": fold_s,
            "computed_scaffold": computed,
            "relation": relation,
            "blocks_freeze": blocks,
            "note": (
                "Historical fold token 'acyclic' matches an empty _murcko scaffold; ring count is 0."
                if fold_token == "acyclic" and relation == "REPRESENTATION_DIFFERENCE"
                else "Computed with jcim_stage2_compute._murcko from current mapping canonical_smiles. String differences with the same InChIKey connectivity block are representation differences."
            ),
        })
    write_csv(OUT / "SCAFFOLD_CONCORDANCE_AUDIT.csv", concordance, list(concordance[0]))
    conflicts = [row for row in concordance if row["blocks_freeze"] == "YES"]
    if conflicts:
        raise SystemExit(f"SCAFFOLD_GROUPING_CONFLICT {len(conflicts)}")

    groups = defaultdict(list)
    for row in assigned:
        if row["scaffold_is_acyclic"] == "YES":
            continue
        groups[(row["pair"], row["class"], row["murcko_scaffold"])].append(row)
    ranks = {}
    for key, members in groups.items():
        ordered = sorted(members, key=lambda row: (
            selection_hash(row["pair"], row["class"], row["global_ligand_entity_id"], row["panel_id"]),
            row["panel_id"],
        ))
        for index, row in enumerate(ordered, start=1):
            ranks[(row["pair"], row["panel_id"])] = (index, len(ordered))

    membership = []
    for row in assigned:
        if row["scaffold_is_acyclic"] == "YES":
            rank, count, retained, reason = "", 0, "YES", "RETAINED"
            count = sum(
                1 for item in assigned
                if item["pair"] == row["pair"] and item["class"] == row["class"] and item["scaffold_is_acyclic"] == "YES"
            )
        else:
            rank, count = ranks[(row["pair"], row["panel_id"])]
            retained = "YES" if count <= 5 or rank <= 5 else "NO"
            reason = "RETAINED" if retained == "YES" else "EXCESS_SAME_CLASS_SCAFFOLD"
        membership.append({
            "pair": row["pair"],
            "panel_id": row["panel_id"],
            "global_ligand_entity_id": row["global_ligand_entity_id"],
            "class": row["class"],
            "canonical_smiles": row["canonical_smiles"],
            "murcko_scaffold": row["murcko_scaffold"],
            "scaffold_is_acyclic": row["scaffold_is_acyclic"],
            "original_class_scaffold_count": count,
            "hash_rank": rank,
            "retained_cap5": retained,
            "exclusion_reason": reason,
        })
    fields = list(membership[0])
    write_csv(OUT / "CAP5_MEMBERSHIP_867.csv", membership, fields)
    membership_sha = sha256_file(OUT / "CAP5_MEMBERSHIP_867.csv")

    excluded = [row for row in membership if row["retained_cap5"] == "NO"]
    write_csv(OUT / "CAP5_EXCLUDED_LIGANDS.csv", excluded or [{key: "" for key in fields}], fields)

    counts = []
    diversity = []
    for pair in sorted({row["pair"] for row in membership}):
        for class_label in ("dual", "A_only", "B_only", "neither"):
            rows = [row for row in membership if row["pair"] == pair and row["class"] == class_label]
            acyclic = [row for row in rows if row["scaffold_is_acyclic"] == "YES"]
            cyclic = [row for row in rows if row["scaffold_is_acyclic"] == "NO"]
            scaffolds = defaultdict(int)
            for row in cyclic:
                scaffolds[row["murcko_scaffold"]] += 1
            over = {name: n for name, n in scaffolds.items() if n > 5}
            retained_rows = [row for row in rows if row["retained_cap5"] == "YES"]
            retained_scaffolds = {row["murcko_scaffold"] for row in retained_rows if row["scaffold_is_acyclic"] == "NO"}
            counts.append({
                "pair": pair,
                "class": class_label,
                "n_original": len(rows),
                "n_acyclic": len(acyclic),
                "n_nonempty_scaffolds": len(scaffolds),
                "max_scaffold_repeat": max(scaffolds.values()) if scaffolds else 0,
                "n_scaffolds_over_5": len(over),
                "n_excluded": sum(row["retained_cap5"] == "NO" for row in rows),
                "n_cap5": len(retained_rows),
                "n_cap5_scaffolds": len(retained_scaffolds),
                "retention_fraction": f"{(len(retained_rows) / len(rows)):.10f}" if rows else "",
            })
            diversity.append(counts[-1])
    write_csv(OUT / "CAP5_CLASS_COUNTS.csv", counts, list(counts[0]))
    write_csv(OUT / "CAP5_SCAFFOLD_DIVERSITY_SUMMARY.csv", diversity, list(diversity[0]))

    freeze = {
        "rule": PREFIX,
        "murcko": "jcim_stage2_compute._murcko",
        "hash_fields": ["BM5_SENSITIVITY_20261023_V1", "pair", "class", "global_ligand_entity_id", "panel_id"],
        "hash_encoding": "utf-8",
        "hash_separator": "|",
        "hash_order": "ascending hex digest, then panel_id",
        "cap": 5,
        "cap_scope": "pair x class x nonempty Murcko scaffold",
        "acyclic_policy": "not subject to the nonempty-scaffold cap",
        "scores_read": False,
        "membership_sha256": membership_sha,
        "n_eligible": len(membership),
        "n_retained": sum(row["retained_cap5"] == "YES" for row in membership),
        "n_excluded": len(excluded),
        "n_acyclic": sum(row["scaffold_is_acyclic"] == "YES" for row in membership),
        "script_version": SCRIPT_VERSION,
        "descriptive_thresholds": manifest["descriptive_thresholds_frozen_before_sensitivity"],
        "sensitivity_statistics": "not_run_in_this_script",
    }
    (OUT / "CAP5_SELECTION_FREEZE.json").write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"retained": freeze["n_retained"], "excluded": freeze["n_excluded"], "acyclic": freeze["n_acyclic"], "conflicts": len(conflicts)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
