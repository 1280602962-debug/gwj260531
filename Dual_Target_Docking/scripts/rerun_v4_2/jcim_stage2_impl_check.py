#!/usr/bin/env python3
"""Artificial-data checks of the real Stage 2 module functions and writers.

Does not read official project members or scores. Does not write results/jcim_stage2.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

from jcim_stage2_compute import (
    module_b,
    module_c,
    module_d,
    module_d1,
    module_e,
    module_f,
    module_h,
    module_i,
    run_stage2_modules,
)
from jcim_stage2_lib import (
    ALT_PAIRS,
    BOOT_SEED,
    NA,
    SCHEMA_B,
    SCHEMA_C,
    SCHEMA_D1,
    SCHEMA_D_OOF,
    SCHEMA_D_UNIT,
    SCHEMA_E_SEED,
    SCHEMA_E_SUM,
    SCHEMA_F,
    SCHEMA_H,
    SCHEMA_I_CMP,
    SCHEMA_I_DESC,
    SCHEMA_I_LIG,
    SEEDS,
    concordance_code,
    paired_negative_class_bootstrap,
    paired_pocket_bootstrap,
    read_csv,
    write_csv,
)


def _row(pair, lig, cls, **scores):
    rec = {
        "pair": pair,
        "canonical_ligand_id": lig,
        "global_ligand_entity_id": lig,
        "class": cls,
        "activity_eligible": "1",
    }
    for method in ("M0", "M1", "M1b", "M2", "M3"):
        rec[f"{method}_A_valid"] = "0"
        rec[f"{method}_B_valid"] = "0"
        rec[f"{method}_score_A"] = ""
        rec[f"{method}_score_B"] = ""
        rec[f"{method}_A_missing"] = ""
        rec[f"{method}_B_missing"] = ""
    rec.update(scores)
    return rec


def _valid(rec, method, side, value):
    rec[f"{method}_{side}_valid"] = "1"
    rec[f"{method}_score_{side}"] = str(value)
    return rec


def check(name, ok, detail, rows):
    rows.append({"check": name, "ok": "PASS" if ok else "FAIL", "detail": detail})
    if not ok:
        print(f"FAIL {name} {detail}")


def test_b(rows):
    pair = "SYN/B"
    dual = [_valid(_valid(_row(pair, f"D{i}", "dual"), "M0", "A", 3), "M0", "B", 5 + i) for i in range(4)]
    single = [_valid(_valid(_row(pair, f"A{i}", "A_only"), "M0", "A", 1), "M0", "B", 4 + i) for i in range(4)]
    neither = [_valid(_valid(_row(pair, f"N{i}", "neither"), "M0", "A", 1), "M0", "B", 0) for i in range(4)]
    pop = dual + single + neither
    rec = next(r for r in module_b(pop, n_boot=80, seed=BOOT_SEED) if r["direction"] == "D_vs_A")
    check("B_positive_delta_sign", float(rec["delta_negative"]) > 0, rec["delta_negative"], rows)
    single_easy = [_valid(_valid(_row(pair, f"A{i}", "A_only"), "M0", "A", 1), "M0", "B", 0) for i in range(4)]
    neither_hard = [_valid(_valid(_row(pair, f"N{i}", "neither"), "M0", "A", 1), "M0", "B", 4 + i) for i in range(4)]
    rec_n = next(r for r in module_b(dual + single_easy + neither_hard, n_boot=80, seed=BOOT_SEED) if r["direction"] == "D_vs_A")
    check("B_negative_delta_sign", float(rec_n["delta_negative"]) < 0, rec_n["delta_negative"], rows)
    # zero: same scores
    same = [_valid(_valid(_row(pair, f"X{i}", "A_only"), "M0", "A", 1), "M0", "B", 1) for i in range(4)]
    same_n = [_valid(_valid(_row(pair, f"Y{i}", "neither"), "M0", "A", 1), "M0", "B", 1) for i in range(4)]
    rec0 = next(r for r in module_b(dual[:1] + [_valid(_valid(_row(pair, "D0b", "dual"), "M0", "A", 3), "M0", "B", 1)] + same + same_n, n_boot=40, seed=BOOT_SEED) if r["direction"] == "D_vs_A")
    # ties in auroc helper
    from jcim_stage2_lib import auroc
    check("B_tie_auroc", abs(auroc([1, 1], [1, 1]) - 0.5) < 1e-12, "", rows)
    empty = module_b(dual + single, n_boot=20, seed=BOOT_SEED)
    rec_e = next(r for r in empty if r["direction"] == "D_vs_A")
    check("B_empty_neither_na", rec_e["delta_negative"] == NA and rec_e["reason"] == "EMPTY_NEITHER_CLASS", rec_e["reason"], rows)
    check("B_empty_keeps_n_single", int(rec_e["n_single"]) == 4, rec_e["n_single"], rows)
    record = []
    packed = paired_negative_class_bootstrap(dual, single, neither, "M0_score_B", n_boot=6, seed=1, record=record)
    shared = all(item[0] == "dual" and item[1] is not None for item in record)
    distinct_neg = any(item[3] != item[5] for item in record)
    check("B_shared_dual_index_recorded", shared and len(record) == 6, str(len(record)), rows)
    check("B_neg_classes_drawn_separately", distinct_neg, "", rows)
    check("B_delta_ci_from_replicate_deltas", "delta_negative_ci_lo" in packed and packed["n_replicates"] == 6, "", rows)
    check("B_ci_matches_delta_reps", abs(packed["delta_negative_ci_lo"] - float(np.percentile(packed["delta_reps"], 2.5, method="linear"))) < 1e-12, "", rows)


def test_c(rows):
    pair = "SYN/C"
    dual = []
    negs = []
    for i in range(4):
        dual.append(_valid(_valid(_row(pair, f"D{i}", "dual"), "M0", "A", 5), "M0", "B", 5))
        negs.append(_valid(_valid(_row(pair, f"A{i}", "A_only"), "M0", "A", 1), "M0", "B", 1))
    rec = next(r for r in module_c(dual + negs, n_boot=40, seed=BOOT_SEED) if r["arm"] == "D_vs_A")
    check("C_identical_pockets_zero_delta", abs(float(rec["delta_pocket"])) < 1e-12, rec["delta_pocket"], rows)
    record = []
    packed = paired_pocket_bootstrap(dual, negs, "M0_score_B", "M0_score_A", n_boot=8, seed=2, record=record)
    check("C_explicit_shared_indices", len(record) == 8 and all(len(x) == 2 for x in record), str(len(record)), rows)
    # swap: make B informative, A anti
    dual2, negs2 = [], []
    for i in range(6):
        dual2.append(_valid(_valid(_row(pair, f"E{i}", "dual"), "M0", "A", 0), "M0", "B", 9))
        negs2.append(_valid(_valid(_row(pair, f"F{i}", "A_only"), "M0", "A", 9), "M0", "B", 0))
    rec2 = next(r for r in module_c(dual2 + negs2, n_boot=80, seed=BOOT_SEED) if r["arm"] == "D_vs_A")
    check("C_correct_better_positive", float(rec2["delta_pocket"]) > 0, rec2["delta_pocket"], rows)
    swapped = paired_pocket_bootstrap(dual2, negs2, "M0_score_A", "M0_score_B", n_boot=80, seed=BOOT_SEED)
    check("C_swap_flips_point", abs(float(swapped["delta_pocket"]) + float(rec2["delta_pocket"])) < 1e-12, str(swapped["delta_pocket"]), rows)
    check(
        "C_swap_flips_ci",
        abs(float(swapped["delta_pocket_ci_lo"]) + float(rec2["delta_pocket_ci_hi"])) < 1e-9
        and abs(float(swapped["delta_pocket_ci_hi"]) + float(rec2["delta_pocket_ci_lo"])) < 1e-9,
        f"{swapped['delta_pocket_ci_lo']},{rec2['delta_pocket_ci_hi']}",
        rows,
    )
    extra_neg = _valid(_row(pair, "DROP1", "A_only"), "M0", "B", 1)
    rec3 = next(r for r in module_c(dual + negs + [extra_neg], n_boot=10, seed=1) if r["arm"] == "D_vs_A")
    check("C_extra_drop_from_members", int(rec3["extra_drop_vs_directional_m0"]) == 1 and "DROP1" in rec3["extra_drop_ids"], rec3["extra_drop_ids"], rows)


def test_d(rows):
    try:
        from rdkit import Chem  # noqa: F401
        from sklearn.linear_model import LogisticRegression  # noqa: F401
    except ImportError as exc:
        check("D_dependency", False, f"DEPENDENCY_MISSING:{exc}", rows)
        return
    smiles_a = ["CCO", "CCN", "CCC", "CCCO", "CCCN", "CCCC", "c1ccccc1", "c1ccncc1", "C1CCCCC1", "CC(C)O"]
    pair_x, pair_y = "SYN/DX", "SYN/DY"
    pop, mapping, folds = [], [], []
    # 5 dual + 5 A_only, folds 1-5; fold 5 test is all dual (last dual only? need one class in test)
    # dual D0-D4 folds 1-5, A_only A0-A4 folds 1-4 and A4 also fold 4 so fold 5 test = D4 only
    for i in range(5):
        lig = f"DXD{i}"
        rec = _valid(_valid(_row(pair_x, lig, "dual"), "M0", "A", 2 + i), "M0", "B", 3 + i)
        pop.append(rec)
        mapping.append({"pair": pair_x, "canonical_ligand_id": lig, "canonical_smiles": smiles_a[i], "parent_independent": "1"})
        folds.append({"pair": pair_x, "arm": "D_vs_A", "ligand_id": lig, "fold_id": str(i + 1), "class": "dual", "scaffold": f"s{i}"})
    for i in range(5):
        lig = f"DXA{i}"
        rec = _valid(_valid(_row(pair_x, lig, "A_only"), "M0", "A", 0.1), "M0", "B", 0.2)
        pop.append(rec)
        mapping.append({"pair": pair_x, "canonical_ligand_id": lig, "canonical_smiles": smiles_a[5 + i], "parent_independent": "1"})
        folds.append({"pair": pair_x, "arm": "D_vs_A", "ligand_id": lig, "fold_id": str(min(i + 1, 4)), "class": "A_only", "scaffold": f"t{i}"})
    # same ligand id different pair/SMILES
    pop.append(_valid(_valid(_row(pair_y, "DXD0", "dual"), "M0", "A", 1), "M0", "B", 1))
    pop.append(_valid(_valid(_row(pair_y, "DYA0", "A_only"), "M0", "A", 0), "M0", "B", 0))
    mapping.append({"pair": pair_y, "canonical_ligand_id": "DXD0", "canonical_smiles": "NCCCCN", "parent_independent": "1"})
    mapping.append({"pair": pair_y, "canonical_ligand_id": "DYA0", "canonical_smiles": "CCCCCCCC", "parent_independent": "1"})
    folds.append({"pair": pair_y, "arm": "D_vs_A", "ligand_id": "DXD0", "fold_id": "1", "class": "dual", "scaffold": "u"})
    folds.append({"pair": pair_y, "arm": "D_vs_A", "ligand_id": "DYA0", "fold_id": "1", "class": "A_only", "scaffold": "v"})
    tables = {"population": pop, "mapping": mapping, "alias": [], "folds": folds, "scaler_probe": []}
    units, oof = module_d(tables)
    unit = next(r for r in units if r["pair"] == pair_x and r["arm"] == "D_vs_A")
    check("D_test_single_class_warning", "5" in str(unit["test_single_class_folds"]) and unit["blocking_reason"] in ("", NA), unit["test_single_class_folds"], rows)
    oof_x = [r for r in oof if r["pair"] == pair_x]
    check("D_all_oof_predictions", len(oof_x) == 10 and all(r["oof_physchem"] not in ("", NA) for r in oof_x), str(len(oof_x)), rows)
    fold5 = [r for r in oof_x if int(r["fold_id"]) == 5]
    check("D_single_class_test_still_predicted", len(fold5) == 1 and fold5[0]["oof_ecfp4"] not in ("", NA), str(fold5), rows)
    ids_p = {r["canonical_ligand_id"] for r in oof_x}
    check("D_three_models_same_members", ids_p == {f"DXD{i}" for i in range(5)} | {f"DXA{i}" for i in range(5)}, str(ids_p), rows)
    check("D_supporting_same_n", unit["auroc_SUPPORTING_CHEMISTRY_MATCHED_M0"] not in ("", NA), str(unit["auroc_SUPPORTING_CHEMISTRY_MATCHED_M0"]), rows)
    d1 = module_d1(tables)
    dxd0 = [r for r in d1 if r["pair"] == pair_x and r["canonical_ligand_id"] == "DXD0"][0]
    dyd0 = [r for r in d1 if r["pair"] == pair_y and r["canonical_ligand_id"] == "DXD0"][0]
    check("D_pair_key_no_cross_pair_smiles", float(dxd0["MolWt"]) != float(dyd0["MolWt"]), f"{dxd0['MolWt']} vs {dyd0['MolWt']}", rows)
    check("D1_acyclic_murcko_empty", dxd0["murcko"] == "", dxd0["murcko"], rows)
    benz = next(r for r in d1 if r["canonical_ligand_id"] == "DXA1")
    check("D1_ring_murcko_nonempty", benz["murcko"] != "", benz["murcko"], rows)
    # train single class
    block_pop, block_map, block_fold = [], [], []
    for i in range(4):
        lig = f"BD{i}"
        block_pop.append(_valid(_valid(_row("SYN/BL", lig, "dual"), "M0", "A", 2), "M0", "B", 2))
        block_map.append({"pair": "SYN/BL", "canonical_ligand_id": lig, "canonical_smiles": smiles_a[i], "parent_independent": "1"})
        block_fold.append({"pair": "SYN/BL", "arm": "D_vs_A", "ligand_id": lig, "fold_id": str(i + 2), "class": "dual", "scaffold": f"b{i}"})
    for i in range(3):
        lig = f"BA{i}"
        block_pop.append(_valid(_valid(_row("SYN/BL", lig, "A_only"), "M0", "A", 0), "M0", "B", 0))
        block_map.append({"pair": "SYN/BL", "canonical_ligand_id": lig, "canonical_smiles": smiles_a[5 + i], "parent_independent": "1"})
        block_fold.append({"pair": "SYN/BL", "arm": "D_vs_A", "ligand_id": lig, "fold_id": "1", "class": "A_only", "scaffold": f"c{i}"})
    b_units, b_oof = module_d({"population": block_pop, "mapping": block_map, "alias": [], "folds": block_fold})
    bu = next(r for r in b_units if r["pair"] == "SYN/BL")
    check("D_train_single_class_blocks", bu["blocking_reason"] == "TRAIN_SINGLE_CLASS" and not b_oof, bu["blocking_reason"], rows)
    probe = tables["scaler_probe"]
    check("D_scaler_probe_from_module_d", bool(probe), str(len(probe)), rows)
    scaler_ok = True
    detail = ""
    for item in probe:
        if item["pair"] != pair_x:
            continue
        train = np.asarray(item["train_m0"], dtype=float)
        test = np.asarray(item["test_m0"], dtype=float)
        if train.size == 0:
            scaler_ok = False
            detail = "empty_train"
            break
        expected_mean = float(train.mean())
        if abs(item["scaler_mean"] - expected_mean) > 1e-12:
            scaler_ok = False
            detail = f"mean {item['scaler_mean']} != {expected_mean}"
            break
        if test.size and np.isclose(test, expected_mean).all() is False:
            if abs(float(test.mean()) - expected_mean) > 1e-12 and abs(item["scaler_mean"] - float(test.mean())) < 1e-12:
                scaler_ok = False
                detail = "scaler_used_test_mean"
                break
        if set(item["test_ids"]) & set(item["train_ids"]):
            scaler_ok = False
            detail = "train_test_id_overlap"
            break
    check("D_scaler_uses_module_d_train_m0", scaler_ok, detail or f"n_probe={len(probe)}", rows)


def test_e(rows):
    pair = "PIK3CA/mTOR"
    pop = [
        _valid(_valid(_row(pair, "SYN1", "dual"), "M0", "A", 1), "M0", "B", 1),
        _valid(_valid(_row(pair, "SYN2", "A_only"), "M0", "A", 0), "M0", "B", 0),
        _valid(_valid(_row(pair, "SYN3", "B_only"), "M0", "A", 0), "M0", "B", 0),
    ]
    master = []
    m0_master = [
        {"pair": pair, "pdb_id": "PDBA", "target_side": "A", "canonical_ligand_id": "SYN1"},
        {"pair": pair, "pdb_id": "PDBB", "target_side": "B", "canonical_ligand_id": "SYN1"},
    ]
    for seed in SEEDS:
        for lig, cls, pdb, side in (
            ("SYN1", "dual", "PDBA", "A"),
            ("SYN1", "dual", "PDBB", "B"),
            ("SYN2", "A_only", "PDBB", "B"),
            ("SYN3", "B_only", "PDBA", "A"),
        ):
            timeout = seed == 17 and lig == "SYN1" and pdb == "PDBB"
            master.append({
                "job_id": f"{lig}_{pdb}_{seed}",
                "pair": pair, "pdb_id": pdb, "canonical_ligand_id": lig, "seed": str(seed),
                "status": "TIMEOUT" if timeout else "SUCCESS",
                "mode1_affinity": "" if timeout else "-2.0",
                "primary_score": "1" if seed == 42 else "0",
            })
            m0_master.append({"pair": pair, "pdb_id": pdb, "target_side": side, "canonical_ligand_id": lig})
    seed_rows, sums = module_e({"population": pop, "fiveseed": master, "m0_master": m0_master})
    s17 = next(r for r in seed_rows if int(r["seed"]) == 17)
    check("E_timeout_listed_not_imputed", "SYN1:B:PDBB:TIMEOUT" in s17["missing_ligands"] and s17["n_dual"] == 0, s17["missing_ligands"], rows)
    check("E_no_primary_score_field_as_affinity", all("primary_score" not in (r.get("missing_jobs") or "") for r in seed_rows), "", rows)
    sm = next(r for r in sums if r["metric"] == "summary_min")
    check("E_incomplete_seeds_no_iqr", sm["reason"].startswith("INCOMPLETE_SEEDS") and sm["iqr"] == NA, sm["reason"], rows)
    # five finite
    master2 = []
    for seed in SEEDS:
        for lig, pdb in (("SYN1", "PDBA"), ("SYN1", "PDBB"), ("SYN2", "PDBB"), ("SYN3", "PDBA")):
            master2.append({
                "job_id": f"{lig}_{pdb}_{seed}", "pair": pair, "pdb_id": pdb,
                "canonical_ligand_id": lig, "seed": str(seed), "status": "SUCCESS",
                "mode1_affinity": str(-1.0 - seed / 100), "primary_score": "0",
            })
    _sr, sums2 = module_e({"population": pop, "fiveseed": master2, "m0_master": m0_master})
    sm2 = next(r for r in sums2 if r["metric"] == "auc_B")
    check("E_five_finite_has_iqr", sm2["iqr"] != NA and sm2["n_finite_seeds"] == 5, str(sm2), rows)


def test_f(rows):
    pair = ALT_PAIRS[0]
    duals = []
    aos = []
    bos = []
    for i in range(4):
        duals.append(_valid(_valid(_row(pair, f"FD{i}", "dual"), "M0", "A", 2), "M0", "B", 2))
        aos.append(_valid(_valid(_row(pair, f"FA{i}", "A_only"), "M0", "A", 2), "M0", "B", 2))
        bos.append(_valid(_valid(_row(pair, f"FB{i}", "B_only"), "M0", "A", 2), "M0", "B", 2))
    alt = []
    for row in duals + aos + bos:
        status = "TIMEOUT" if row["canonical_ligand_id"] == "FD3" else "SUCCESS"
        if status != "SUCCESS":
            score = ""
        elif row["class"] == "dual":
            score = "9"
        else:
            score = "1"
        alt.append({
            "pair": pair, "alt_id": "ALT1", "side": "B", "canonical_ligand_id": row["canonical_ligand_id"],
            "status_alt": status, "score_alt": score, "alt_pdb": "4L2Y",
        })
    tables = {"population": duals + aos + bos, "alt_independent": alt}
    recs = module_f(tables, n_boot=40, seed=BOOT_SEED)
    rec = recs[0]
    check("F_dropped_timeout_dual", rec["n_dual_common"] == 3 and "FD3" not in rec["common_dual_ids"], rec["common_dual_ids"], rows)
    check("F_unreplaced_zero_when_A_unchanged", abs(float(rec["delta_AUC_unreplaced"])) < 1e-12, rec["delta_AUC_unreplaced"], rows)
    check("F_affected_can_change", float(rec["delta_AUC_affected"]) != 0.0, rec["delta_AUC_affected"], rows)


def test_i(rows):
    pair = "SYN/I"
    pop = []
    for i in range(6):
        r = _row(pair, f"ID{i}", "dual")
        m0ab = (10, 10) if i == 0 else (3, 3)
        m1ab = (10, 10) if i <= 1 else (1, 1)
        m2ab = (1, 1)
        _valid(_valid(r, "M0", "A", m0ab[0]), "M0", "B", m0ab[1])
        _valid(_valid(r, "M1", "A", m1ab[0]), "M1", "B", m1ab[1])
        _valid(_valid(r, "M1b", "A", m1ab[0]), "M1b", "B", m1ab[1])
        _valid(_valid(r, "M2", "A", m2ab[0]), "M2", "B", m2ab[1])
        _valid(_valid(r, "M3", "A", m0ab[0]), "M3", "B", m0ab[1])
        pop.append(r)
    for i in range(6):
        r = _row(pair, f"IA{i}", "A_only")
        m0ab = (10, 10) if i == 0 else (1, 1)
        m1ab = (1, 1)
        m2ab = (10, 10) if i <= 1 else (1, 1)
        _valid(_valid(r, "M0", "A", m0ab[0]), "M0", "B", m0ab[1])
        _valid(_valid(r, "M1", "A", m1ab[0]), "M1", "B", m1ab[1])
        _valid(_valid(r, "M1b", "A", m1ab[0]), "M1b", "B", m1ab[1])
        _valid(_valid(r, "M2", "A", m2ab[0]), "M2", "B", m2ab[1])
        _valid(_valid(r, "M3", "A", m0ab[0]), "M3", "B", m0ab[1])
        pop.append(r)
    for i in range(6):
        r = _row(pair, f"IB{i}", "B_only")
        _valid(_valid(r, "M0", "A", 1), "M0", "B", 1)
        _valid(_valid(r, "M1", "A", 1), "M1", "B", 1)
        _valid(_valid(r, "M1b", "A", 1), "M1b", "B", 1)
        _valid(_valid(r, "M2", "A", 1), "M2", "B", 1)
        if i < 5:
            _valid(_valid(r, "M3", "A", 1), "M3", "B", 1)
        pop.append(r)
    primary = [
        {"pair": pair, "delta": "M1-M0", "delta_summary_min": "0.2"},
        {"pair": pair, "delta": "M1b-M1", "delta_summary_min": "0.05"},
        {"pair": pair, "delta": "M2-M0", "delta_summary_min": "-0.2"},
        {"pair": pair, "delta": "M3-M0", "delta_summary_min": "0.0"},
    ]
    desc, ligs, cmps = module_i(pop, primary)
    m1 = next(r for r in cmps if r["method"] == "M1")
    m1b = next(r for r in cmps if r["method"] == "M1b")
    m2 = next(r for r in cmps if r["method"] == "M2")
    m3 = next(r for r in cmps if r["method"] == "M3")
    check("I_m1_reads_primary", abs(float(m1["delta_summary_min"]) - 0.2) < 1e-12, m1["delta_summary_min"], rows)
    check("I_m1b_sum_when_membership_ok", abs(float(m1b["delta_summary_min"]) - 0.25) < 1e-12, m1b["delta_summary_min"], rows)
    check("I_four_state_improve", m1["concordance"] == "CONCORDANT_IMPROVE", m1["concordance"], rows)
    check("I_four_state_worsen", m2["concordance"] == "CONCORDANT_WORSEN", m2["concordance"], rows)
    check("I_four_state_no_change", m3["concordance"] == "NO_CHANGE", m3["concordance"], rows)
    check("I_concordance_same_sign_discordant", concordance_code(0.1, 0.2) == "DISCORDANT", "", rows)
    # membership stop
    pop2 = [dict(r) for r in pop]
    for r in pop2:
        if r["canonical_ligand_id"] == "IB5":
            r["M1b_A_valid"] = "0"
            r["M1b_score_A"] = ""
    _d, _l, cmps2 = module_i(pop2, primary)
    m1b2 = next(r for r in cmps2 if r["method"] == "M1b")
    check("I_m1b_stops_if_membership_differs", m1b2["link_source"] == "STOPPED" and m1b2["delta_summary_min"] == NA, m1b2["link_reason"], rows)
    # Third pairwise population differs: extra dual has M1b+M1 but no M0.
    pop3 = [dict(r) for r in pop]
    extra = _row(pair, "ID_EXTRA", "dual")
    _valid(_valid(extra, "M1", "A", 4), "M1", "B", 4)
    _valid(_valid(extra, "M1b", "A", 4), "M1b", "B", 4)
    pop3.append(extra)
    _d, _l, cmps3 = module_i(pop3, primary)
    m1b3 = next(r for r in cmps3 if r["method"] == "M1b")
    check(
        "I_m1b_stops_if_m1b_m1_third_group_differs",
        m1b3["link_source"] == "STOPPED" and m1b3["delta_summary_min"] == NA and "dual" in m1b3["link_reason"] and "ID_EXTRA" in m1b3["link_reason"],
        m1b3["link_reason"],
        rows,
    )
    check("I_ligand_ranks_present", len(ligs) > 0 and all(r.get("canonical_ligand_id") and r.get("rank_A") != "" for r in ligs), str(len(ligs)), rows)
    check("I_descriptive_counts", any(r["n_A_only_top"] != "" for r in desc), "", rows)


def test_writer_and_path(rows):
    tmp = Path(tempfile.mkdtemp(prefix="jcim_stage2_impl_"))
    # first-row NA then later deltas
    na = {k: NA for k in SCHEMA_B}
    na.update({"module": "B", "pair": "P0", "direction": "D_vs_A", "reason": "EMPTY_NEITHER_CLASS", "bootstrap_started": "NO"})
    ok = {k: NA for k in SCHEMA_B}
    ok.update({"module": "B", "pair": "P1", "direction": "D_vs_A", "delta_negative": 0.2, "delta_negative_ci_lo": 0.1, "delta_negative_ci_hi": 0.3, "bootstrap_started": "YES"})
    write_csv(tmp / "b.csv", [na, ok], SCHEMA_B)
    back = read_csv(tmp / "b.csv")
    check("CSV_first_na_keeps_delta_fields", back[1]["delta_negative"] == "0.2" and "delta_negative_ci_lo" in back[0], str(back[0].keys()), rows)
    write_csv(tmp / "empty.csv", [], SCHEMA_H)
    empty = read_csv(tmp / "empty.csv")
    check("CSV_empty_has_header", empty == [], "header_only_ok", rows)
    # full path
    pair = ALT_PAIRS[0]
    pop = []
    mapping, folds, fiveseed, m0_master, alt, primary = [], [], [], [], [], []
    smiles = ["CCO", "CCN", "CCC", "CCCO", "c1ccccc1", "C1CCCCC1", "CC(C)O", "CCCC", "CCCN", "CCCl"]
    for i in range(4):
        r = _row(pair, f"PD{i}", "dual")
        for m in ("M0", "M1", "M1b", "M2", "M3"):
            _valid(_valid(r, m, "A", 4 + i), m, "B", 4 + i)
        pop.append(r)
        mapping.append({"pair": pair, "canonical_ligand_id": f"PD{i}", "canonical_smiles": smiles[i], "parent_independent": "1"})
        folds.append({"pair": pair, "arm": "D_vs_A", "ligand_id": f"PD{i}", "fold_id": str((i % 4) + 1), "class": "dual", "scaffold": f"g{i}"})
        folds.append({"pair": pair, "arm": "D_vs_B", "ligand_id": f"PD{i}", "fold_id": str((i % 4) + 1), "class": "dual", "scaffold": f"g{i}"})
    for i in range(4):
        r = _row(pair, f"PA{i}", "A_only")
        for m in ("M0", "M1", "M1b", "M2", "M3"):
            _valid(_valid(r, m, "A", 0.2), m, "B", 0.2)
        pop.append(r)
        mapping.append({"pair": pair, "canonical_ligand_id": f"PA{i}", "canonical_smiles": smiles[4 + i], "parent_independent": "1"})
        folds.append({"pair": pair, "arm": "D_vs_A", "ligand_id": f"PA{i}", "fold_id": str((i % 4) + 1), "class": "A_only", "scaffold": f"h{i}"})
    for i in range(4):
        r = _row(pair, f"PB{i}", "B_only")
        for m in ("M0", "M1", "M1b", "M2", "M3"):
            _valid(_valid(r, m, "A", 0.2), m, "B", 0.2)
        pop.append(r)
        mapping.append({"pair": pair, "canonical_ligand_id": f"PB{i}", "canonical_smiles": smiles[6 + (i % 4)], "parent_independent": "1"})
        folds.append({"pair": pair, "arm": "D_vs_B", "ligand_id": f"PB{i}", "fold_id": str((i % 4) + 1), "class": "B_only", "scaffold": f"i{i}"})
    n = _row(pair, "PN0", "neither")
    _valid(_valid(n, "M0", "A", 0), "M0", "B", 0)
    pop.append(n)
    for seed in SEEDS:
        for lig, pdb, side in (("PD0", "PA", "A"), ("PD0", "PB", "B"), ("PA0", "PB", "B"), ("PB0", "PA", "A")):
            fiveseed.append({"job_id": f"{lig}{pdb}{seed}", "pair": pair, "pdb_id": pdb, "canonical_ligand_id": lig if lig != "PA0" else "PA0", "seed": str(seed), "status": "SUCCESS", "mode1_affinity": "-1.5", "primary_score": "0"})
            m0_master.append({"pair": pair, "pdb_id": pdb, "target_side": side, "canonical_ligand_id": lig})
    # complete fiveseed for all pop ligands/sides
    for seed in SEEDS:
        for r in pop:
            for side, pdb in (("A", "PA"), ("B", "PB")):
                if r["class"] == "neither":
                    continue
                fiveseed.append({"job_id": f"{r['canonical_ligand_id']}{pdb}{seed}x", "pair": pair, "pdb_id": pdb, "canonical_ligand_id": r["canonical_ligand_id"], "seed": str(seed), "status": "SUCCESS", "mode1_affinity": "-1.2", "primary_score": "0"})
                m0_master.append({"pair": pair, "pdb_id": pdb, "target_side": side, "canonical_ligand_id": r["canonical_ligand_id"]})
    for r in pop:
        if r["class"] == "neither":
            continue
        alt.append({"pair": pair, "alt_id": "ALTZ", "side": "B", "canonical_ligand_id": r["canonical_ligand_id"], "status_alt": "SUCCESS", "score_alt": r["M0_score_B"], "alt_pdb": "4L2Y"})
    primary.extend([
        {"pair": pair, "delta": "M1-M0", "delta_summary_min": "0.1"},
        {"pair": pair, "delta": "M1b-M1", "delta_summary_min": "0.0"},
        {"pair": pair, "delta": "M2-M0", "delta_summary_min": "-0.1"},
        {"pair": pair, "delta": "M3-M0", "delta_summary_min": "0.05"},
    ])
    tables = {
        "population": pop, "mapping": mapping, "alias": [], "folds": folds,
        "fiveseed": fiveseed, "m0_master": m0_master, "alt_independent": alt, "primary_delta": primary,
    }
    run_stage2_modules(tables, tmp / "out", n_boot=30, seed=BOOT_SEED)
    expected = [
        ("module_B_negative_class.csv", SCHEMA_B),
        ("module_C_wrong_pocket.csv", SCHEMA_C),
        ("module_D_chemistry_units.csv", SCHEMA_D_UNIT),
        ("module_D_chemistry_oof_ligands.csv", SCHEMA_D_OOF),
        ("module_D1_composition.csv", SCHEMA_D1),
        ("module_E_five_seed.csv", SCHEMA_E_SEED),
        ("module_E_five_seed_summary.csv", SCHEMA_E_SUM),
        ("module_F_alt_receptor.csv", SCHEMA_F),
        ("module_H_missingness.csv", SCHEMA_H),
        ("module_I_descriptive.csv", SCHEMA_I_DESC),
        ("module_I_ligand_ranks.csv", SCHEMA_I_LIG),
        ("module_I_method_compare.csv", SCHEMA_I_CMP),
    ]
    for name, schema in expected:
        path = tmp / "out" / name
        data = read_csv(path)
        check(f"PATH_reread_{name}", path.exists() and list(data[0].keys()) == list(schema) if data else list(open(path).readline().strip().split(",")) == list(schema), name, rows)
    man = json.loads((tmp / "out" / "RUN_MANIFEST.json").read_text())
    check("PATH_manifest_e_no_false_bootstrap", man["module_statistics"]["E"]["bootstrap_B"] is None, str(man["module_statistics"]["E"]), rows)
    check("PATH_no_real_outdir", True, str(tmp), rows)
    print(f"IMPL_TMP {tmp}")


def main() -> int:
    rows = []
    test_b(rows)
    test_c(rows)
    test_d(rows)
    test_e(rows)
    test_f(rows)
    test_i(rows)
    test_writer_and_path(rows)
    failed = [r for r in rows if r["ok"] != "PASS"]
    print(f"IMPL_CHECKS {len(rows) - len(failed)}/{len(rows)} PASS")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
