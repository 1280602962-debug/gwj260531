#!/usr/bin/env python3
"""Stage 2 prep audit. Does not train models or write scientific results.

Population acceptance recounts members from raw CSVs in this file.
It does not call jcim_stage2_lib.build_layer_members.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
from collections import Counter, defaultdict
from pathlib import Path

from jcim_stage2_lib import (
    ALT_PAIRS,
    ALT_PDBS,
    BLOCKED_MODULES,
    EXECUTABLE_MODULES,
    FORBIDDEN_IMPORT_SUBSTRINGS,
    PATHS,
    four_sided_membership_identical,
    m1b_link_membership,
    read_csv,
    resolve_project_root,
    synthetic_self_check,
)

ARMS = (("D_vs_A", {"dual", "A_only"}), ("D_vs_B", {"dual", "B_only"}))
PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _m0_ok_independent(row: dict, arm: str) -> bool:
    if row["class"] == "dual":
        return row.get("M0_A_valid") == "1" and row.get("M0_B_valid") == "1"
    if arm == "D_vs_A":
        return row.get("M0_B_valid") == "1"
    return row.get("M0_A_valid") == "1"


def independent_layer_table(pop, alias_ids, fold_keys) -> list[dict]:
    rows = []
    for pair in PAIRS:
        for arm, classes in ARMS:
            labeled = [
                r
                for r in pop
                if r["pair"] == pair
                and r.get("activity_eligible") == "1"
                and r.get("canonical_ligand_id") not in alias_ids
                and r.get("class") in classes
            ]
            m0 = [r for r in labeled if _m0_ok_independent(r, arm)]
            un_l = [r["canonical_ligand_id"] for r in labeled if (pair, arm, r["canonical_ligand_id"]) not in fold_keys]
            un_m = [r["canonical_ligand_id"] for r in m0 if (pair, arm, r["canonical_ligand_id"]) not in fold_keys]
            chem = [r for r in m0 if (pair, arm, r["canonical_ligand_id"]) in fold_keys]
            miss = [
                r["canonical_ligand_id"]
                for r in labeled
                if not _m0_ok_independent(r, arm)
            ]
            rows.append(
                {
                    "pair": pair,
                    "arm": arm,
                    "n_labeled_population": len(labeled),
                    "n_m0_directional_population": len(m0),
                    "n_m0_missing": len(miss),
                    "n_fold_unassigned_in_labeled_population": len(un_l),
                    "n_fold_unassigned_after_m0_filter": len(un_m),
                    "n_final_chemistry_oof": len(chem),
                    "labeled_ids": tuple(sorted(r["canonical_ligand_id"] + ":" + r["class"] for r in labeled)),
                    "m0_ids": tuple(sorted(r["canonical_ligand_id"] + ":" + r["class"] for r in m0)),
                    "chem_ids": tuple(sorted(r["canonical_ligand_id"] + ":" + r["class"] for r in chem)),
                    "fold_unassigned_in_labeled": ";".join(un_l),
                    "fold_unassigned_after_m0": ";".join(un_m),
                    "m0_missing_ids": ";".join(miss),
                }
            )
    return rows


def forbidden_imports(script: Path) -> list[str]:
    tree = ast.parse(script.read_text(), filename=str(script))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [((node.module or "") + "." + a.name) for a in node.names]
            names.append(node.module or "")
        else:
            continue
        for name in names:
            for bad in FORBIDDEN_IMPORT_SUBSTRINGS:
                if bad in name:
                    found.append(f"{script.name}:{name}")
    return found


def smiles_ok(mapping_rows) -> tuple[int, int]:
    try:
        from rdkit import Chem
    except ImportError:
        return -1, -1
    ok = 0
    for row in mapping_rows:
        mol = Chem.MolFromSmiles(row.get("canonical_smiles") or "")
        if mol is not None:
            ok += 1
    return ok, len(mapping_rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    root = resolve_project_root(args.project_root)
    script_dir = Path(__file__).resolve().parent
    checks = []

    def add(name, ok, detail):
        checks.append({"check": name, "ok": "PASS" if ok else "FAIL", "detail": detail})

    pop = read_csv(root / PATHS["population"])
    folds = read_csv(root / PATHS["folds"])
    alias = read_csv(root / PATHS["alias"])
    mapping = read_csv(root / PATHS["mapping"])
    m0_master = read_csv(root / PATHS["m0_master"])
    fiveseed = read_csv(root / PATHS["fiveseed"])
    alt = read_csv(root / PATHS["alt_independent"])
    primary = read_csv(root / PATHS["primary_metrics"])
    delta = read_csv(root / PATHS["primary_delta"])
    alias_ids = {r["alias_id"] for r in alias}
    fold_keys = {(r["pair"], r["arm"], r["ligand_id"]) for r in folds}

    layers = independent_layer_table(pop, alias_ids, fold_keys)
    n1 = sum(r["n_labeled_population"] for r in layers)
    n2 = sum(r["n_m0_directional_population"] for r in layers)
    n3 = sum(r["n_final_chemistry_oof"] for r in layers)
    un_l = sum(r["n_fold_unassigned_in_labeled_population"] for r in layers)
    un_m = sum(r["n_fold_unassigned_after_m0_filter"] for r in layers)
    un_l_ids = sorted({x.split(":")[0] for r in layers for x in r["fold_unassigned_in_labeled"].split(";") if x})
    un_m_ids = sorted({x for r in layers for x in r["fold_unassigned_after_m0"].split(";") if x})
    add("layer_counts_938_934_928", n1 == 938 and n2 == 934 and n3 == 928, f"{n1}/{n2}/{n3}")
    add("fold_unassigned_labeled_10_7", un_l == 10 and len(un_l_ids) == 7, f"{un_l}/{len(un_l_ids)}:{','.join(un_l_ids)}")
    add("fold_unassigned_after_m0_6_4", un_m == 6 and len(un_m_ids) == 4, f"{un_m}/{len(un_m_ids)}:{','.join(un_m_ids)}")
    add("alias_ab040_only", alias_ids == {"AB_040"}, str(sorted(alias_ids)))
    add("ab040_not_in_pop", all(r["canonical_ligand_id"] != "AB_040" for r in pop), "")

    published = root / "analysis_plan_jcim" / "CHEMISTRY_POPULATION_LAYERS.csv"
    pub_rows = { (r["pair"], r["arm"]): r for r in read_csv(published) } if published.exists() else {}
    layer_match = True
    if not pub_rows:
        layer_match = False
    for row in layers:
        pub = pub_rows.get((row["pair"], row["arm"]))
        if not pub:
            layer_match = False
            break
        for key in (
            "n_labeled_population",
            "n_m0_directional_population",
            "n_fold_unassigned_in_labeled_population",
            "n_fold_unassigned_after_m0_filter",
            "n_final_chemistry_oof",
        ):
            if str(row[key]) != str(pub[key]):
                layer_match = False
    add("published_layer_table_matches_independent", layer_match, str(published))

    ache_a = [r for r in pop if r["pair"] == "AChE/BChE" and r.get("activity_eligible") == "1" and r["class"] == "A_only"]
    ache_a_m0 = [r for r in ache_a if r.get("M0_B_valid") == "1"]
    ache_a_both = [r for r in ache_a_m0 if r.get("M0_A_valid") == "1" and r.get("M0_B_valid") == "1"]
    extra = len(ache_a_m0) - len(ache_a_both)
    add("wrongpocket_ache_extra_drop_0", extra == 0 and len(ache_a) == 27 and len(ache_a_m0) == 25, f"lab={len(ache_a)} m0={len(ache_a_m0)} both={len(ache_a_both)} extra={extra}")

    add("primary_metrics_120", len(primary) == 120, str(len(primary)))
    add("primary_delta_32", len(delta) == 32, str(len(delta)))
    add("population_807", len(pop) == 807, str(len(pop)))

    m0_lookup = {(r["pair"], r["canonical_ligand_id"], r["target_side"]): r for r in m0_master}
    m0_mismatch = 0
    for row in pop:
        for side, score_key, valid_key in (("A", "M0_score_A", "M0_A_valid"), ("B", "M0_score_B", "M0_B_valid")):
            src = m0_lookup.get((row["pair"], row["canonical_ligand_id"], side))
            if src is None:
                m0_mismatch += 1
                continue
            if src.get("vina_status") == "SUCCESS" and row.get(valid_key) == "1":
                if abs(float(src["vina_score"]) - float(row[score_key])) > 1e-9:
                    m0_mismatch += 1
            elif src.get("vina_status") == "TIMEOUT" and row.get(valid_key) == "1":
                m0_mismatch += 1
    add("seed42_master_matches_formal_m0", m0_mismatch == 0, f"mismatches={m0_mismatch}")

    fs42 = [r for r in fiveseed if r["seed"] == "42"]
    add("fiveseed_seed42_1614", len(fs42) == 1614, str(len(fs42)))
    fs_mismatch = 0
    pdb_side = {(r["pair"], r["pdb_id"]): r["target_side"] for r in m0_master}
    for row in fs42:
        side = pdb_side.get((row["pair"], row["pdb_id"]))
        src = m0_lookup.get((row["pair"], row["canonical_ligand_id"], side)) if side else None
        if src is None:
            fs_mismatch += 1
            continue
        if row["status"] != src["vina_status"]:
            fs_mismatch += 1
            continue
        if row["status"] == "SUCCESS":
            expected = -float(row["mode1_affinity"])
            if abs(expected - float(src["vina_score"])) > 1e-6:
                fs_mismatch += 1
    add("fiveseed_seed42_matches_m0_master", fs_mismatch == 0, f"mismatches={fs_mismatch}")

    ab046 = [r for r in fiveseed if r["canonical_ligand_id"] == "AB_046" and r["seed"] == "17" and r["pdb_id"] == "4BDS"]
    add("ab046_seed17_4bds_timeout", len(ab046) == 1 and ab046[0]["status"] == "TIMEOUT", str(ab046[0] if ab046 else None))

    alt_mismatch = 0
    for row in alt:
        if row["pair"] not in ALT_PAIRS:
            alt_mismatch += 1
        for side, key in (("A", "score_A_primary"), ("B", "score_B_primary")):
            src = m0_lookup.get((row["pair"], row["canonical_ligand_id"], side))
            if src is None:
                alt_mismatch += 1
                continue
            if src.get("vina_status") == "SUCCESS" and row.get(key) not in ("", None):
                if abs(float(src["vina_score"]) - float(row[key])) > 1e-6:
                    alt_mismatch += 1
    add("alt_primary_scores_match_formal_m0", alt_mismatch == 0, f"mismatches={alt_mismatch}")
    add("alt_pdbs_frozen", set(r["alt_pdb"] for r in alt) <= set(ALT_PDBS), str(sorted(set(r["alt_pdb"] for r in alt))))

    by_pair = defaultdict(list)
    for row in pop:
        by_pair[row["pair"]].append(row)
    ident = all(four_sided_membership_identical(rows, "M1b", "M1", "M0") for rows in by_pair.values())
    add("m1b_m0_membership_equals_m1_m0", ident, "")
    three = {pair: m1b_link_membership(rows) for pair, rows in by_pair.items()}
    three_fail = [f"{pair}:{';'.join(v['diffs'])}" for pair, v in three.items() if not v["identical"]]
    add("m1b_three_pairwise_membership_identical", not three_fail, ",".join(three_fail) if three_fail else "all_pairs_identical")

    fold_dups = [k for k, n in Counter((r["pair"], r["arm"], r["ligand_id"]) for r in folds).items() if n > 1]
    add("fold_keys_unique", not fold_dups, str(fold_dups[:3]))
    scaffold_fold = defaultdict(set)
    for row in folds:
        scaffold_fold[(row["pair"], row["arm"], row.get("scaffold") or "")].add(row["fold_id"])
    leak = [k for k, v in scaffold_fold.items() if k[2] and len(v) > 1]
    add("scaffold_not_split_across_folds", not leak, str(leak[:3]))

    parse_ok, parse_n = smiles_ok(mapping)
    if parse_ok < 0:
        add("smiles_parse_808", False, "DEPENDENCY_MISSING/NOT_VERIFIED")
    else:
        add("smiles_parse_808", parse_ok == parse_n == 808, f"{parse_ok}/{parse_n}")

    for name in ("jcim_stage2_lib.py", "jcim_stage2_audit.py", "jcim_stage2_compute.py"):
        hits = forbidden_imports(script_dir / name)
        add(f"no_historical_imports_{name}", not hits, ",".join(hits))

    issues = synthetic_self_check()
    add("synthetic_unit_checks", not issues, ",".join(issues))
    add("blocked_modules_gjk", BLOCKED_MODULES == ("G", "J", "K"), str(BLOCKED_MODULES))
    add("executable_modules_bcdefhi", EXECUTABLE_MODULES == ("B", "C", "D", "E", "F", "H", "I"), str(EXECUTABLE_MODULES))

    import subprocess

    repo = root.parent
    base = "d365c90183473c566162552fff1483006519195a"
    authority = [
        PATHS["primary_metrics"],
        PATHS["primary_delta"],
        PATHS["population"],
        PATHS["m0_master"],
        PATHS["fiveseed"],
        PATHS["alt_independent"],
        PATHS["folds"],
        PATHS["mapping"],
        Path("reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml"),
        Path("reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml"),
        Path("reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/FORMAL_AUTHORITY_PATHS.yaml"),
        Path("reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml"),
    ]
    present = all((root / rel).exists() for rel in authority)
    add("authority_files_present", present, "presence_only_not_an_edit_check")
    dirty = []
    for rel in authority:
        repo_rel = Path("Dual_Target_Docking") / rel
        diff = subprocess.check_output(["git", "diff", "--name-only", base, "--", str(repo_rel)], cwd=repo, text=True).strip()
        work = subprocess.check_output(["git", "diff", "--name-only", "--", str(repo_rel)], cwd=repo, text=True).strip()
        if diff or work:
            dirty.append(str(rel))
    add("authority_files_unchanged_vs_d365c901", not dirty, ",".join(dirty))

    valid_finite = 0
    for row in pop:
        for method in ("M0", "M1", "M1b", "M2", "M3"):
            for side in ("A", "B"):
                flag = row.get(f"{method}_{side}_valid") == "1"
                try:
                    val = row.get(f"{method}_score_{side}")
                    is_f = val not in (None, "") and float(val) == float(val) and abs(float(val)) != float("inf")
                except (TypeError, ValueError):
                    is_f = False
                if flag != is_f:
                    valid_finite += 1
    add("valid_flag_matches_finite_score", valid_finite == 0, str(valid_finite))
    pop_cls = {(r["pair"], r["canonical_ligand_id"]): r["class"] for r in pop}
    fold_cls_mismatch = sum(1 for r in folds if pop_cls.get((r["pair"], r["ligand_id"]), r["class"]) != r["class"])
    add("fold_class_matches_population", fold_cls_mismatch == 0, str(fold_cls_mismatch))
    indep_map = [
        r for r in mapping
        if r["canonical_ligand_id"] not in alias_ids
        and str(r.get("parent_independent", "1")) not in {"0", "false", "False"}
    ]
    map_counts = Counter((r["pair"], r["canonical_ligand_id"]) for r in indep_map)
    map_dups = [k for k, n in map_counts.items() if n > 1]
    add("mapping_pair_id_unique_among_independent", not map_dups, str(map_dups[:5]))

    failed = [c for c in checks if c["ok"] != "PASS"]
    report_path = root / "analysis_plan_jcim" / "STAGE2_PREP_ACCEPTANCE.md"
    if args.write_report:
        lines = [
            "# STAGE2_PREP_ACCEPTANCE",
            "",
            "This report is read-only prep/acceptance only. It does not prove module implementations.",
            "Implementation evidence is `jcim_stage2_impl_check.py` on artificial data.",
            "",
            "SECOND_STAGE_COMPUTE_EXECUTED=NO",
            "",
            "| check | result | detail |",
            "|---|---|---|",
        ]
        for c in checks:
            lines.append(f"| {c['check']} | {c['ok']} | {c['detail']} |")
        lines.extend(
            [
                "",
                f"Independent layer totals: labeled={n1} m0={n2} chemistry_oof={n3}",
                f"Labeled fold-unassigned records/ligands: {un_l}/{len(un_l_ids)} ({', '.join(un_l_ids)})",
                f"After-M0 fold-unassigned records/ligands: {un_m}/{len(un_m_ids)} ({', '.join(un_m_ids)})",
                "",
                "Population members were recounted in this audit file from the raw CSVs.",
                "The audit does not call `jcim_stage2_lib.build_layer_members`.",
                "",
                "Blocked modules remain G, J, K.",
            ]
        )
        report_path.write_text("\n".join(lines) + "\n")
        print(f"WROTE {report_path}")
    for c in checks:
        print(f"{c['ok']:4} {c['check']} {c['detail']}")
    print(f"FAILED {len(failed)}/{len(checks)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
