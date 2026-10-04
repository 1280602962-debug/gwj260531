#!/usr/bin/env python3
"""Read-only acceptance of results/jcim_stage2.

Point estimates, membership, and ranks are recomputed here from the
authority CSVs and the result tables. This file does not call the
stage2 module functions and does not rerun bootstrap.
"""
from __future__ import annotations

import csv
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT_DEFAULT = Path(__file__).resolve().parents[2]
BASE = "3f3e0d26052178c38ccb43dec637683a60f97a9d"
NA = "NA"
PAIRS_CODE = None


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def finite(value) -> bool:
    try:
        return value not in (None, "", NA) and np.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def fnum(value) -> float:
    return float(value)


def auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def midrank_desc(values) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(-x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * ((i + 1) + (j + 1))
        i = j + 1
    return ranks


def close(a, b, tol=1e-8) -> bool:
    if not finite(a) or not finite(b):
        return False
    return abs(float(a) - float(b)) <= tol


def ci_ok(lo, hi) -> bool:
    return finite(lo) and finite(hi) and float(lo) <= float(hi) + 1e-12


def main() -> int:
    root = ROOT_DEFAULT
    out = root / "results" / "jcim_stage2"
    qa = root / "reruns" / "UNIFORM_RERUN_V4_2_20260921" / "13_qa"
    prot = root / "reruns" / "UNIFORM_RERUN_V4_2_20260921" / "00_protocol"
    checks = []

    def add(name, ok, detail=""):
        checks.append((name, "PASS" if ok else "FAIL", detail))
        if not ok:
            print(f"FAIL {name} {detail}")

    pop = read_csv(qa / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv")
    folds = read_csv(root / "results" / "canonical" / "model_fold_assignments.csv")
    alias = {r["alias_id"] for r in read_csv(prot / "parent_alias_registry.csv")}
    primary = read_csv(root / "results" / "formal_metrics" / "PRIMARY_DIRECTIONAL_METRICS.csv")
    delta = read_csv(root / "results" / "formal_metrics" / "PRIMARY_METHOD_DELTA_METRICS.csv")
    fiveseed = read_csv(qa / "phase7_fiveseed_master.csv")
    m0_master = read_csv(qa / "official_primary_seed42_score_master_8pair.csv")
    alt = read_csv(qa / "alt_score_master_independent.csv")
    manifest = (out / "RUN_MANIFEST.json").read_text()

    b = read_csv(out / "module_B_negative_class.csv")
    c = read_csv(out / "module_C_wrong_pocket.csv")
    d_units = read_csv(out / "module_D_chemistry_units.csv")
    d_oof = read_csv(out / "module_D_chemistry_oof_ligands.csv")
    d1 = read_csv(out / "module_D1_composition.csv")
    e_seed = read_csv(out / "module_E_five_seed.csv")
    e_sum = read_csv(out / "module_E_five_seed_summary.csv")
    frows = read_csv(out / "module_F_alt_receptor.csv")
    h = read_csv(out / "module_H_missingness.csv")
    i_desc = read_csv(out / "module_I_descriptive.csv")
    i_lig = read_csv(out / "module_I_ligand_ranks.csv")
    i_cmp = read_csv(out / "module_I_method_compare.csv")

    add("counts", (
        len(b) == 16 and len(c) == 16 and len(d_units) == 16 and len(d1) == 938
        and len(e_seed) == 40 and len(e_sum) == 24 and len(frows) == 8
        and len(h) == 320 and len(i_desc) == 104 and len(i_cmp) == 32
    ), f"B{len(b)} C{len(c)} D{len(d_units)} D1{len(d1)} OOF{len(d_oof)} E{len(e_seed)}/{len(e_sum)} F{len(frows)} H{len(h)} Idesc{len(i_desc)} Icmp{len(i_cmp)}")

    def uniq(rows, keys):
        return [k for k, n in Counter(tuple(r[k] for k in keys) for r in rows).items() if n > 1]

    add("keys_B", not uniq(b, ("pair", "direction")), "")
    add("keys_C", not uniq(c, ("pair", "arm")), "")
    add("keys_D", not uniq(d_units, ("pair", "arm")), "")
    add("keys_OOF", not uniq(d_oof, ("pair", "arm", "canonical_ligand_id")), "")
    add("keys_D1", not uniq(d1, ("pair", "arm", "canonical_ligand_id")), "")
    add("keys_E", not uniq(e_seed, ("pair", "seed")), "")
    add("keys_F", not uniq(frows, ("pair", "alt_id", "replaced_side")), "")
    add("keys_Icmp", not uniq(i_cmp, ("pair", "method")), "")

    def range_ok(rows, fields):
        bad = 0
        for row in rows:
            for field in fields:
                if finite(row.get(field)) and not (0.0 <= float(row[field]) <= 1.0):
                    bad += 1
        return bad

    add("auroc_in_unit_interval_B", range_ok(b, ("auc_single", "auc_neither")) == 0, "")
    add("auroc_in_unit_interval_C", range_ok(c, ("auc_correct", "auc_wrong")) == 0, "")
    add("auroc_in_unit_interval_D", range_ok(d_units, ("auroc_physchem", "auroc_ecfp4", "auroc_ecfp4_m0", "auroc_SUPPORTING_CHEMISTRY_MATCHED_M0")) == 0, "")

    def delta_range(rows, fields):
        bad = []
        for row in rows:
            for field in fields:
                if finite(row.get(field)) and not (-1.0 - 1e-9 <= float(row[field]) <= 1.0 + 1e-9):
                    bad.append((row.get("pair"), field, row[field]))
        return bad

    add("delta_range_BCF", not delta_range(b, ("delta_negative",)) and not delta_range(c, ("delta_pocket",)) and not delta_range(frows, ("delta_AUC_affected", "delta_AUC_unreplaced", "delta_summary_min")), "")

    def ci_rows(rows, pairs):
        bad = 0
        for row in rows:
            if row.get("bootstrap_started") != "YES":
                continue
            for lo, hi in pairs:
                if not ci_ok(row.get(lo), row.get(hi)):
                    bad += 1
        return bad

    add("ci_bounds_B", ci_rows(b, (("auc_single_ci_lo", "auc_single_ci_hi"), ("auc_neither_ci_lo", "auc_neither_ci_hi"), ("delta_negative_ci_lo", "delta_negative_ci_hi"))) == 0, "")
    add("ci_bounds_C", ci_rows(c, (("auc_correct_ci_lo", "auc_correct_ci_hi"), ("auc_wrong_ci_lo", "auc_wrong_ci_hi"), ("delta_pocket_ci_lo", "delta_pocket_ci_hi"))) == 0, "")
    add("ci_bounds_F", ci_rows(frows, (("delta_AUC_affected_ci_lo", "delta_AUC_affected_ci_hi"), ("delta_AUC_unreplaced_ci_lo", "delta_AUC_unreplaced_ci_hi"), ("delta_summary_min_ci_lo", "delta_summary_min_ci_hi"))) == 0, "")

    boot_bad = []
    for label, rows in (("B", b), ("C", c), ("F", frows)):
        for row in rows:
            if row.get("bootstrap_started") == "YES" and row.get("n_replicates") != "10000":
                boot_bad.append((label, row.get("pair"), row.get("n_replicates")))
            if row.get("bootstrap_started") == "NO" and row.get("n_replicates") not in ("0", ""):
                boot_bad.append((label, row.get("pair"), "unexpected_" + row.get("n_replicates", "")))
    add("n_replicates_10000_when_started", not boot_bad, str(boot_bad[:4]))
    add("manifest_seed", '"seed": 271828' in manifest and manifest.count('"bootstrap_B": 10000') == 3 and '"bootstrap_B": null' in manifest, "")

    pri = {(r["pair"], r["method"], r["contrast"]): r for r in primary}
    b_mismatch = []
    for row in b:
        contrast = "dual_vs_A_only" if row["direction"] == "D_vs_A" else "dual_vs_B_only"
        src = pri[(row["pair"], "M0", contrast)]
        if row["auc_single"] not in ("", NA) and not close(row["auc_single"], src["auroc"]):
            b_mismatch.append((row["pair"], row["direction"], row["auc_single"], src["auroc"]))
        if finite(row["auc_single"]) and finite(row["auc_neither"]):
            expect = float(row["auc_neither"]) - float(row["auc_single"])
            if not close(row["delta_negative"], expect):
                b_mismatch.append(("sign", row["pair"], row["delta_negative"], expect))
    add("B_single_matches_primary_M0_and_sign", not b_mismatch, str(b_mismatch[:3]))
    pik = [r for r in b if r["pair"] == "PIK3CA/mTOR"]
    add("B_pik3ca_neither_kept", all(int(r["n_neither"]) == 4 and r["bootstrap_started"] == "YES" for r in pik), str([(r["direction"], r["n_neither"]) for r in pik]))

    # C extra drop and pocket sign from raw members
    c_bad = []
    for row in c:
        rows = [r for r in pop if r["pair"] == row["pair"] and r.get("activity_eligible") == "1"]
        neg_cls = "A_only" if row["arm"] == "D_vs_A" else "B_only"
        directional = [
            r for r in rows if r["class"] == neg_cls and (
                r.get("M0_B_valid") == "1" if row["arm"] == "D_vs_A" else r.get("M0_A_valid") == "1"
            )
        ]
        both = [r for r in directional if r.get("M0_A_valid") == "1" and r.get("M0_B_valid") == "1"]
        extra_ids = sorted(r["canonical_ligand_id"] for r in directional if r not in both)
        if int(row["extra_drop_vs_directional_m0"]) != len(extra_ids):
            c_bad.append(("n", row["pair"], row["arm"], row["extra_drop_vs_directional_m0"], len(extra_ids)))
        got = [x for x in row["extra_drop_ids"].split(";") if x]
        if sorted(got) != extra_ids:
            c_bad.append(("ids", row["pair"], row["arm"], got, extra_ids))
        if finite(row["auc_correct"]) and finite(row["auc_wrong"]) and not close(row["delta_pocket"], float(row["auc_correct"]) - float(row["auc_wrong"])):
            c_bad.append(("sign", row["pair"], row["arm"]))
    add("C_extra_drop_and_sign", not c_bad, str(c_bad[:3]))
    ache = next(r for r in c if r["pair"] == "AChE/BChE" and r["arm"] == "D_vs_A")
    add("C_ache_extra_drop_0", ache["extra_drop_vs_directional_m0"] == "0" and ache["n_negative_both_finite"] == "25", ache["n_negative_directional_m0"] + "->" + ache["n_negative_both_finite"])

    # D layers
    fold_keys = {(r["pair"], r["arm"], r["ligand_id"]) for r in folds}
    labeled_n = m0_n = chem_n = 0
    chem_ids = set()
    for pair in sorted({r["pair"] for r in pop}):
        for arm, classes in (("D_vs_A", {"dual", "A_only"}), ("D_vs_B", {"dual", "B_only"})):
            lab = [
                r for r in pop
                if r["pair"] == pair and r.get("activity_eligible") == "1"
                and r["canonical_ligand_id"] not in alias and r["class"] in classes
            ]
            def m0ok(r, arm=arm):
                if r["class"] == "dual":
                    return r.get("M0_A_valid") == "1" and r.get("M0_B_valid") == "1"
                return r.get("M0_B_valid") == "1" if arm == "D_vs_A" else r.get("M0_A_valid") == "1"
            m0 = [r for r in lab if m0ok(r)]
            chem = [r for r in m0 if (pair, arm, r["canonical_ligand_id"]) in fold_keys]
            labeled_n += len(lab)
            m0_n += len(m0)
            chem_n += len(chem)
            for r in chem:
                chem_ids.add((pair, arm, r["canonical_ligand_id"]))
    add("layers_938_934_928", labeled_n == 938 and m0_n == 934 and chem_n == 928, f"{labeled_n}/{m0_n}/{chem_n}")
    oof_ids = {(r["pair"], r["arm"], r["canonical_ligand_id"]) for r in d_oof}
    blocked_units = [r for r in d_units if r.get("blocking_reason") not in ("", NA) or r.get("incomplete_reason") not in ("", NA)]
    executable_members = 0
    for unit in d_units:
        if unit.get("blocking_reason") in ("", NA) and unit.get("incomplete_reason") in ("", NA):
            executable_members += int(unit["n_final_chemistry_oof"])
    add("D_no_unexpected_block", not blocked_units, str([(r["pair"], r["arm"], r["blocking_reason"], r["incomplete_reason"]) for r in blocked_units]))
    add("D_oof_equals_executable_members", len(d_oof) == executable_members == 928 and oof_ids == chem_ids, f"oof={len(d_oof)} executable={executable_members}")
    pred_bad = 0
    for row in d_oof:
        if row["canonical_ligand_id"] in alias:
            pred_bad += 1
        for col in ("oof_physchem", "oof_ecfp4", "oof_ecfp4_m0", "relevant_m0"):
            if not finite(row[col]):
                pred_bad += 1
    add("D_three_models_finite_same_rows", pred_bad == 0 and len(d_oof) == len(oof_ids), str(pred_bad))
    # recompute AUROC from OOF
    d_re = []
    for unit in d_units:
        rows = [r for r in d_oof if r["pair"] == unit["pair"] and r["arm"] == unit["arm"]]
        if not rows:
            continue
        y = np.array([1 if r["class"] == "dual" else 0 for r in rows])
        for name, col in (
            ("auroc_physchem", "oof_physchem"),
            ("auroc_ecfp4", "oof_ecfp4"),
            ("auroc_ecfp4_m0", "oof_ecfp4_m0"),
            ("auroc_SUPPORTING_CHEMISTRY_MATCHED_M0", "relevant_m0"),
        ):
            pred = np.array([float(r[col]) for r in rows])
            got = auroc(pred[y == 1], pred[y == 0])
            if not close(unit[name], got):
                d_re.append((unit["pair"], unit["arm"], name, unit[name], got))
        inc = float(unit["auroc_ecfp4_m0"]) - float(unit["auroc_ecfp4"])
        if not close(unit["delta_ecfp4_plus_m0"], inc):
            d_re.append((unit["pair"], unit["arm"], "increment", unit["delta_ecfp4_plus_m0"], inc))
    add("D_oof_auroc_and_increment", not d_re, str(d_re[:3]))
    d1_ids = {(r["pair"], r["arm"], r["canonical_ligand_id"]) for r in d1}
    add("D1_includes_timeout_fold_gap", ("AChE/BChE", "D_vs_A", "AB_001") in d1_ids and ("AChE/BChE", "D_vs_A", "AB_053") in d1_ids and len(d1) == labeled_n, "")
    test_single_obs = sum(1 for r in d_units if r.get("test_single_class_folds") not in ("", NA))
    add("D_test_single_class_observation_not_a_gate", True, f"units_with_test_single_warning={test_single_obs}")

    # E
    pdb_side = {(r["pair"], r["pdb_id"]): r["target_side"] for r in m0_master}
    e_bad = []
    ab_ok = False
    for row in e_seed:
        if row.get("bootstrap_started") != "NO":
            e_bad.append(("boot", row["pair"], row["seed"]))
        if row["pair"] == "AChE/BChE" and row["seed"] == "17" and "AB_046:B:4BDS:TIMEOUT" in row["missing_ligands"]:
            ab_ok = True
    add("E_ab046_timeout_listed", ab_ok, "")
    add("E_no_bootstrap", not e_bad, str(e_bad[:2]))
    # seed42 point vs recomputed -mode1 and vs primary
    by_key = {}
    for row in fiveseed:
        if row["seed"] != "42":
            continue
        side = pdb_side.get((row["pair"], row["pdb_id"]))
        if side:
            by_key[(row["pair"], row["canonical_ligand_id"], side)] = row
    e42_bad = []
    for row in e_seed:
        if row["seed"] != "42":
            continue
        rebuilt = []
        for src in pop:
            if src["pair"] != row["pair"]:
                continue
            item = dict(src)
            for side, valid, score in (("A", "M0_A_valid", "M0_score_A"), ("B", "M0_B_valid", "M0_score_B")):
                job = by_key.get((row["pair"], src["canonical_ligand_id"], side))
                if job and job["status"] == "SUCCESS" and finite(job.get("mode1_affinity")):
                    item[valid] = "1"
                    item[score] = str(-float(job["mode1_affinity"]))
                else:
                    item[valid] = "0"
                    item[score] = ""
            rebuilt.append(item)
        dual = [r for r in rebuilt if r["class"] == "dual" and r["M0_A_valid"] == "1" and r["M0_B_valid"] == "1"]
        ao = [r for r in rebuilt if r["class"] == "A_only" and r["M0_B_valid"] == "1"]
        bo = [r for r in rebuilt if r["class"] == "B_only" and r["M0_A_valid"] == "1"]
        if len(dual) != int(row["n_dual"]) or len(ao) != int(row["n_A_only"]) or len(bo) != int(row["n_B_only"]):
            e42_bad.append(("n", row["pair"], len(dual), row["n_dual"]))
        if dual and ao:
            auc_b = auroc([float(r["M0_score_B"]) for r in dual], [float(r["M0_score_B"]) for r in ao])
            if not close(row["auc_B"], auc_b):
                e42_bad.append(("aucB", row["pair"], row["auc_B"], auc_b))
            pri_b = pri[(row["pair"], "M0", "dual_vs_A_only")]["auroc"]
            if not close(row["auc_B"], pri_b):
                e42_bad.append(("priB", row["pair"], row["auc_B"], pri_b))
    add("E_seed42_matches_recompute_and_primary", not e42_bad, str(e42_bad[:3]))
    sum_bad = []
    by_metric = {(r["pair"], r["seed"]): r for r in e_seed}
    for row in e_sum:
        vals = []
        for seed in ("17", "29", "42", "71", "101"):
            src = by_metric[(row["pair"], seed)][row["metric"]]
            if finite(src):
                vals.append(float(src))
        if len(vals) != 5:
            if row["iqr"] != NA or not str(row["reason"]).startswith("INCOMPLETE"):
                sum_bad.append((row["pair"], row["metric"], "should_be_NA", row["reason"]))
        else:
            ordered = [float(by_metric[(row["pair"], s)][row["metric"]]) for s in ("17", "29", "42", "71", "101")]
            p25, p75 = np.percentile(ordered, [25, 75], method="linear")
            if not close(row["iqr"], p75 - p25) or not close(row["median"], float(np.median(ordered))):
                sum_bad.append((row["pair"], row["metric"], row["iqr"], p75 - p25))
            if row["seed42"] != by_metric[(row["pair"], "42")][row["metric"]]:
                sum_bad.append(("seed42col", row["pair"], row["metric"]))
    add("E_summary_linear_only_if_five_finite", not sum_bad, str(sum_bad[:3]))

    # F point reconstruction
    f_bad = []
    alt_lookup = {(r["pair"], r["alt_id"], r["side"], r["canonical_ligand_id"]): r for r in alt}
    for row in frows:
        common = []
        for src in pop:
            if src["pair"] != row["pair"] or src.get("activity_eligible") != "1":
                continue
            if src["class"] not in {"dual", "A_only", "B_only"}:
                continue
            if src.get("M0_A_valid") != "1" or src.get("M0_B_valid") != "1":
                continue
            alt_row = alt_lookup.get((row["pair"], row["alt_id"], row["replaced_side"], src["canonical_ligand_id"]))
            if alt_row is None or alt_row.get("status_alt") != "SUCCESS" or not finite(alt_row.get("score_alt")):
                continue
            common.append((src, alt_row))
        dual = [(s, a) for s, a in common if s["class"] == "dual"]
        ao = [(s, a) for s, a in common if s["class"] == "A_only"]
        bo = [(s, a) for s, a in common if s["class"] == "B_only"]
        if len(dual) != int(row["n_dual_common"]) or len(ao) != int(row["n_A_only_common"]) or len(bo) != int(row["n_B_only_common"]):
            f_bad.append(("n", row["pair"], row["alt_id"], len(dual), row["n_dual_common"]))
            continue
        def scores(side_rows, key):
            return [float(r[key]) for r in side_rows]
        pri_b = auroc(scores([s for s, _ in dual], "M0_score_B"), scores([s for s, _ in ao], "M0_score_B"))
        pri_a = auroc(scores([s for s, _ in dual], "M0_score_A"), scores([s for s, _ in bo], "M0_score_A"))
        def alt_score(src, alt_row, side):
            if side == row["replaced_side"]:
                return float(alt_row["score_alt"])
            return float(src[f"M0_score_{side}"])
        alt_b = auroc(
            [alt_score(s, a, "B") for s, a in dual],
            [alt_score(s, a, "B") for s, a in ao],
        )
        alt_a = auroc(
            [alt_score(s, a, "A") for s, a in dual],
            [alt_score(s, a, "A") for s, a in bo],
        )
        d_b, d_a = alt_b - pri_b, alt_a - pri_a
        unrep = d_a if row["replaced_side"] == "B" else d_b
        aff = d_b if row["replaced_side"] == "B" else d_a
        dsm = min(alt_b, alt_a) - min(pri_b, pri_a)
        if abs(unrep) > 1e-12:
            f_bad.append(("unreplaced", row["pair"], row["alt_id"], unrep))
        if not close(row["delta_AUC_unreplaced"], unrep) or not close(row["delta_AUC_affected"], aff) or not close(row["delta_summary_min"], dsm):
            f_bad.append(("delta", row["pair"], row["alt_id"], row["delta_AUC_unreplaced"], unrep))
        if not (float(row["delta_AUC_unreplaced_ci_lo"]) == 0.0 and float(row["delta_AUC_unreplaced_ci_hi"]) == 0.0):
            f_bad.append(("ci", row["pair"], row["alt_id"], row["delta_AUC_unreplaced_ci_lo"], row["delta_AUC_unreplaced_ci_hi"]))
    add("F_common_set_unreplaced_zero_and_summary_min", not f_bad, str(f_bad[:3]))

    h_bad = []
    for row in h:
        if int(row["expected"]) != int(row["valid"]) + int(row["missing"]):
            h_bad.append(row["pair"])
    add("H_expected_equals_valid_plus_missing", not h_bad, str(len(h_bad)))

    # I ranks
    def pool(rows, methods):
        out_rows = []
        for r in rows:
            if r.get("activity_eligible") != "1" or r.get("class") not in {"dual", "A_only", "B_only"}:
                continue
            ok = True
            for method in methods:
                if r.get(f"{method}_A_valid") != "1" or r.get(f"{method}_B_valid") != "1":
                    ok = False
            if ok:
                out_rows.append(r)
        return out_rows

    def rank_table(rows, method):
        if not rows:
            return []
        ra = midrank_desc([float(r[f"{method}_score_A"]) for r in rows])
        rb = midrank_desc([float(r[f"{method}_score_B"]) for r in rows])
        enriched = []
        for r, a, b in zip(rows, ra, rb):
            enriched.append({
                "canonical_ligand_id": r["canonical_ligand_id"],
                "global_ligand_entity_id": r["global_ligand_entity_id"],
                "class": r["class"],
                "rank_A": float(a),
                "rank_B": float(b),
                "worst_rank": max(float(a), float(b)),
                "tie1": float(a) + float(b),
            })
        enriched.sort(key=lambda r: (r["worst_rank"], r["tie1"], r["global_ligand_entity_id"]))
        k = int(math.ceil(0.10 * len(enriched)))
        for i, r in enumerate(enriched, start=1):
            r["joint_order"] = i
            r["in_topk"] = "YES" if i <= k else "NO"
        return enriched

    lig_index = {(r["pair"], r["method"], r["population_kind"], r["canonical_ligand_id"]): r for r in i_lig}
    rank_bad = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = [r for r in pop if r["pair"] == pair]
        for method in ("M0", "M1", "M1b", "M2", "M3"):
            ranked = rank_table(pool(rows, (method,)), method)
            for item in ranked:
                got = lig_index.get((pair, method, "method_specific", item["canonical_ligand_id"]))
                if got is None or not close(got["rank_A"], item["rank_A"]) or not close(got["worst_rank"], item["worst_rank"]) or got["in_topk"] != item["in_topk"] or int(float(got["joint_order"])) != item["joint_order"]:
                    rank_bad.append(("spec", pair, method, item["canonical_ligand_id"]))
                    break
        for method in ("M1", "M1b", "M2", "M3"):
            common = pool(rows, ("M0", method))
            for tag, meth in ((f"pairwise_common_vs_{method}", "M0"), ("pairwise_common_vs_M0", method)):
                ranked = rank_table(common, meth)
                k = int(math.ceil(0.10 * len(ranked))) if ranked else 0
                classes = Counter(r["class"] for r in ranked[:k])
                desc = next(r for r in i_desc if r["pair"] == pair and r["method"] == meth and r["population_kind"] == tag)
                if int(float(desc["k"])) != k or int(desc["n_dual_top"]) != classes["dual"] or int(desc["n_A_only_top"]) != classes["A_only"] or int(desc["n_B_only_top"]) != classes["B_only"]:
                    rank_bad.append(("top", pair, meth, tag, desc["k"], k))
                if ranked and not close(desc["expected_dual_count"], k * sum(1 for r in ranked if r["class"] == "dual") / len(ranked)):
                    rank_bad.append(("exp", pair, tag))
                for item in ranked:
                    got = lig_index.get((pair, meth, tag, item["canonical_ligand_id"]))
                    if got is None or got["in_topk"] != item["in_topk"] or not close(got["worst_rank"], item["worst_rank"]):
                        rank_bad.append(("lig", pair, meth, tag, item["canonical_ligand_id"]))
                        break
    add("I_independent_ranks_topk_expectation", not rank_bad, str(rank_bad[:3]))

    dlt = {(r["pair"], r["delta"]): r for r in delta}
    i_bad = []
    states = Counter()
    for row in i_cmp:
        if row["method"] == "M1b":
            if row["link_source"] == "STOPPED" or row["delta_summary_min"] in ("", NA):
                i_bad.append(("unexpected_stop", row["pair"], row["link_reason"]))
                continue
            expect = float(dlt[(row["pair"], "M1-M0")]["delta_summary_min"]) + float(dlt[(row["pair"], "M1b-M1")]["delta_summary_min"])
            if not close(row["delta_summary_min"], expect) or row["supporting_only"] != "YES" or row["bootstrap_ci_added"] != "NO":
                i_bad.append(("m1b", row["pair"], row["delta_summary_min"], expect))
        else:
            expect = float(dlt[(row["pair"], f"{row['method']}-M0")]["delta_summary_min"])
            if not close(row["delta_summary_min"], expect):
                i_bad.append(("primary", row["pair"], row["method"], row["delta_summary_min"], expect))
        if finite(row["delta_summary_min"]) and finite(row["delta_top10_single_target_fraction"]):
            dsm = float(row["delta_summary_min"])
            dst = float(row["delta_top10_single_target_fraction"])
            if dsm == 0.0 or dst == 0.0:
                code = "NO_CHANGE"
            elif dsm > 0 and dst < 0:
                code = "CONCORDANT_IMPROVE"
            elif dsm < 0 and dst > 0:
                code = "CONCORDANT_WORSEN"
            else:
                code = "DISCORDANT"
            states[row["method"] + ":" + code] += 1
            if row["concordance"] != code:
                i_bad.append(("state", row["pair"], row["method"], row["concordance"], code))
    add("I_primary_link_and_four_state", not i_bad, str(i_bad[:3]))
    add("I_m1b_supporting_eight", sum(1 for r in i_cmp if r["method"] == "M1b" and r["supporting_only"] == "YES") == 8, "")

    repo = root.parent
    protected = [
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_lib.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_compute.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_audit.py",
        "Dual_Target_Docking/scripts/rerun_v4_2/jcim_stage2_impl_check.py",
        "Dual_Target_Docking/results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv",
        "Dual_Target_Docking/results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv",
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml",
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml",
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/FORMAL_AUTHORITY_PATHS.yaml",
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml",
    ]
    dirty = []
    for rel in protected:
        diff = subprocess.check_output(["git", "diff", "--name-only", BASE, "--", rel], cwd=repo, text=True).strip()
        work = subprocess.check_output(["git", "status", "--porcelain", "--", rel], cwd=repo, text=True).strip()
        if diff or work:
            dirty.append(rel)
    add("protected_science_unchanged_vs_3f3e0d260", not dirty, ",".join(dirty))

    failed = [c for c in checks if c[1] != "PASS"]
    print(f"RESULT_AUDIT {len(checks) - len(failed)}/{len(checks)} PASS")
    print("FOUR_STATE", dict(states))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
