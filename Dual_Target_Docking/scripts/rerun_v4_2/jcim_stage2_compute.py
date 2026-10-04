#!/usr/bin/env python3
"""Stage 2 compute entry for modules B, C, D, E, F, H, I.

Does nothing unless --execute-compute is passed.
Forbidden: modules G, J, K; PRIMARY overwrite; historical alt stats imports.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from jcim_stage2_lib import (
    ALT_PAIRS,
    BLOCKED_MODULES,
    BOOT_B,
    BOOT_SEED,
    ECFP_NBITS,
    ECFP_RADIUS,
    EXECUTABLE_MODULES,
    LR_C,
    LR_MAX_ITER,
    LR_RANDOM_STATE,
    LR_SOLVER,
    METHODS,
    NA,
    PHYSCHEM_NAMES,
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
    alias_ids,
    auroc,
    build_layer_members,
    concordance_code,
    empty_record,
    finite,
    five_seed_iqr,
    fold_keyset,
    four_sided_membership_sets,
    m1b_link_membership,
    joint_rank_rows,
    load_stage2_tables,
    mapping_smiles,
    member_keys,
    paired_alt_bootstrap,
    paired_negative_class_bootstrap,
    paired_pocket_bootstrap,
    pairwise_ranking_pool,
    ranking_pool,
    resolve_project_root,
    scores,
    shared_dual_universe,
    topk_composition,
    write_csv,
)

OUTPUT_REL = Path("results/jcim_stage2")


def _rows_for_pair(pop, pair):
    return [r for r in pop if r["pair"] == pair]


def module_b(pop, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> list[dict]:
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        uni = shared_dual_universe(rows, "M0")
        dual = uni["dual_shared"]
        for direction, single, neither, key in (
            ("D_vs_A", uni["A_only_for_AUC_B"],
             [r for r in rows if r.get("activity_eligible") == "1" and r["class"] == "neither" and r.get("M0_B_valid") == "1"],
             "M0_score_B"),
            ("D_vs_B", uni["B_only_for_AUC_A"],
             [r for r in rows if r.get("activity_eligible") == "1" and r["class"] == "neither" and r.get("M0_A_valid") == "1"],
             "M0_score_A"),
        ):
            rec = empty_record(
                SCHEMA_B, module="B", pair=pair, direction=direction,
                n_dual=len(dual), n_single=len(single), n_neither=len(neither),
                bootstrap_started="NO", n_replicates=0,
            )
            packed = paired_negative_class_bootstrap(dual, single, neither, key, n_boot=n_boot, seed=seed)
            rec["reason"] = packed.get("reason", "")
            rec["bootstrap_started"] = packed.get("bootstrap_started", "NO")
            rec["n_replicates"] = packed.get("n_replicates", 0)
            for field in (
                "auc_single", "auc_neither", "delta_negative",
                "auc_single_ci_lo", "auc_single_ci_hi",
                "auc_neither_ci_lo", "auc_neither_ci_hi",
                "delta_negative_ci_lo", "delta_negative_ci_hi",
            ):
                if field in packed:
                    rec[field] = packed[field]
            if not dual:
                rec["reason"] = "INSUFFICIENT_n_dual"
            out.append(rec)
    return out


def module_c(pop, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> list[dict]:
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        dual = [
            r for r in rows
            if r.get("activity_eligible") == "1" and r["class"] == "dual"
            and r.get("M0_A_valid") == "1" and r.get("M0_B_valid") == "1"
        ]
        for arm, neg_cls, correct, wrong in (
            ("D_vs_A", "A_only", "M0_score_B", "M0_score_A"),
            ("D_vs_B", "B_only", "M0_score_A", "M0_score_B"),
        ):
            directional = [
                r for r in rows
                if r.get("activity_eligible") == "1" and r["class"] == neg_cls
                and (r.get("M0_B_valid") == "1" if arm == "D_vs_A" else r.get("M0_A_valid") == "1")
            ]
            both = [r for r in directional if r.get("M0_A_valid") == "1" and r.get("M0_B_valid") == "1"]
            extra = [r["canonical_ligand_id"] for r in directional if r not in both]
            rec = empty_record(
                SCHEMA_C, module="C", pair=pair, arm=arm,
                n_dual=len(dual), n_negative_directional_m0=len(directional),
                n_negative_both_finite=len(both),
                extra_drop_vs_directional_m0=len(extra),
                extra_drop_ids=";".join(extra),
                bootstrap_started="NO", n_replicates=0,
            )
            if not dual or not both:
                rec["reason"] = "INSUFFICIENT"
                out.append(rec)
                continue
            packed = paired_pocket_bootstrap(dual, both, correct, wrong, n_boot=n_boot, seed=seed)
            rec.update({k: packed[k] for k in packed if k != "delta_reps"})
            out.append(rec)
    return out


def _physchem(mol):
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    return {
        "MolWt": float(Descriptors.MolWt(mol)),
        "MolLogP": float(Crippen.MolLogP(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "NumHDonors": float(Lipinski.NumHDonors(mol)),
        "NumHAcceptors": float(Lipinski.NumHAcceptors(mol)),
        "NumRotatableBonds": float(Lipinski.NumRotatableBonds(mol)),
        "RingCount": float(rdMolDescriptors.CalcNumRings(mol)),
        "FormalCharge": float(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
    }


def _murcko(mol):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    core = MurckoScaffold.GetScaffoldForMol(mol)
    if core is None or core.GetNumAtoms() == 0:
        return ""
    return Chem.MolToSmiles(core)


def _ecfp(mol):
    from rdkit import DataStructs
    from rdkit.Chem import AllChem

    return AllChem.GetMorganFingerprintAsBitVect(
        mol, ECFP_RADIUS, nBits=ECFP_NBITS, useChirality=False, useFeatures=False
    )


def _ecfp_arr(mol):
    import numpy as np

    fp = _ecfp(mol)
    arr = np.zeros((ECFP_NBITS,), dtype=float)
    on = fp.GetOnBits() if hasattr(fp, "GetOnBits") else []
    for idx in on:
        arr[int(idx)] = 1.0
    return arr


def module_d1(tables) -> list[dict]:
    from rdkit import Chem, DataStructs

    pop = tables["population"]
    aliases = alias_ids(tables["alias"])
    smiles = mapping_smiles(tables["mapping"], aliases)
    folds = fold_keyset(tables["folds"])
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        for arm in ("D_vs_A", "D_vs_B"):
            labeled = build_layer_members(pop, aliases, folds, pair, arm)["labeled"]
            mols = {}
            for row in labeled:
                smi = smiles.get((pair, row["canonical_ligand_id"]), "")
                mol = Chem.MolFromSmiles(smi) if smi else None
                mols[row["canonical_ligand_id"]] = mol
            fps = {lig: _ecfp(mol) for lig, mol in mols.items() if mol is not None}
            for row in labeled:
                rec = empty_record(SCHEMA_D1, pair=pair, arm=arm,
                                   canonical_ligand_id=row["canonical_ligand_id"])
                rec["class"] = row["class"]
                mol = mols.get(row["canonical_ligand_id"])
                if mol is None:
                    rec["max_ecfp4_tanimoto_other_class"] = NA
                    out.append(rec)
                    continue
                rec.update(_physchem(mol))
                rec["murcko"] = _murcko(mol)
                others = [fps[o["canonical_ligand_id"]] for o in labeled
                          if o["class"] != row["class"] and o["canonical_ligand_id"] in fps]
                if others and row["canonical_ligand_id"] in fps:
                    rec["max_ecfp4_tanimoto_other_class"] = max(
                        DataStructs.TanimotoSimilarity(fps[row["canonical_ligand_id"]], other) for other in others
                    )
                else:
                    rec["max_ecfp4_tanimoto_other_class"] = NA
                out.append(rec)
    return out


def module_d(tables) -> tuple[list[dict], list[dict]]:
    from rdkit import Chem
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    pop = tables["population"]
    aliases = alias_ids(tables["alias"])
    fold_keys = fold_keyset(tables["folds"])
    fold_of = {(r["pair"], r["arm"], r["ligand_id"]): r["fold_id"] for r in tables["folds"]}
    smiles = mapping_smiles(tables["mapping"], aliases)
    units, oof_rows = [], []
    for pair in sorted({r["pair"] for r in pop}):
        for arm in ("D_vs_A", "D_vs_B"):
            layers = build_layer_members(pop, aliases, fold_keys, pair, arm)
            chem = layers["chemistry_oof"]
            un_l = [r["canonical_ligand_id"] for r in layers["labeled"] if (pair, arm, r["canonical_ligand_id"]) not in fold_keys]
            un_m = [r["canonical_ligand_id"] for r in layers["m0"] if (pair, arm, r["canonical_ligand_id"]) not in fold_keys]
            rec = empty_record(
                SCHEMA_D_UNIT, module="D", pair=pair, arm=arm,
                n_labeled_population=len(layers["labeled"]),
                n_m0_directional_population=len(layers["m0"]),
                n_final_chemistry_oof=len(chem),
                fold_unassigned_in_labeled=";".join(un_l),
                fold_unassigned_after_m0=";".join(un_m),
                blocking_reason="", evidence_status="", incomplete_reason="",
            )
            if not chem:
                rec["blocking_reason"] = "EMPTY_CHEMISTRY_OOF"
                units.append(rec)
                continue
            y, fold_ids, phys, fps, m0s, ids, classes = [], [], [], [], [], [], []
            parse_fail = False
            for row in chem:
                mol = Chem.MolFromSmiles(smiles.get((pair, row["canonical_ligand_id"]), ""))
                if mol is None:
                    parse_fail = True
                    break
                y.append(1 if row["class"] == "dual" else 0)
                fold_ids.append(int(fold_of[(pair, arm, row["canonical_ligand_id"])]))
                phys.append([_physchem(mol)[k] for k in PHYSCHEM_NAMES])
                fps.append(_ecfp_arr(mol))
                key = "M0_score_B" if arm == "D_vs_A" else "M0_score_A"
                m0s.append(float(row[key]))
                ids.append(row["canonical_ligand_id"])
                classes.append(row["class"])
            if parse_fail:
                rec["blocking_reason"] = "CHEMISTRY_SMILES_PARSE_FAILURE"
                units.append(rec)
                continue
            y = np.asarray(y)
            fold_arr = np.asarray(fold_ids)
            phys = np.asarray(phys, dtype=float)
            fps = np.asarray(fps, dtype=float)
            m0s = np.asarray(m0s, dtype=float).reshape(-1, 1)
            uniq = sorted(set(fold_arr.tolist()))
            test_warn = []
            blocked = ""
            for fid in uniq:
                tr = fold_arr != fid
                te = fold_arr == fid
                if len(set(y[tr].tolist())) < 2:
                    blocked = "TRAIN_SINGLE_CLASS"
                if len(set(y[te].tolist())) < 2:
                    test_warn.append(str(fid))
            rec["test_single_class_folds"] = ";".join(test_warn)
            if test_warn:
                rec["evidence_status"] = "TEST_SINGLE_CLASS_WARNING"
            if blocked:
                rec["blocking_reason"] = blocked
                units.append(rec)
                continue

            def fit_lr(Xt, yt):
                model = LogisticRegression(
                    penalty="l2", C=LR_C, solver=LR_SOLVER, class_weight=None,
                    max_iter=LR_MAX_ITER, random_state=LR_RANDOM_STATE,
                )
                model.fit(Xt, yt)
                return model

            def oof_physchem():
                pred = np.full(len(y), np.nan)
                for fid in uniq:
                    tr, te = fold_arr != fid, fold_arr == fid
                    scaler = StandardScaler()
                    Xt = scaler.fit_transform(phys[tr])
                    Xe = scaler.transform(phys[te])
                    pred[te] = fit_lr(Xt, y[tr]).predict_proba(Xe)[:, 1]
                return pred

            def oof_ecfp():
                pred = np.full(len(y), np.nan)
                for fid in uniq:
                    tr, te = fold_arr != fid, fold_arr == fid
                    pred[te] = fit_lr(fps[tr], y[tr]).predict_proba(fps[te])[:, 1]
                return pred

            def oof_ecfp_m0():
                pred = np.full(len(y), np.nan)
                probe = tables.get("scaler_probe")
                for fid in uniq:
                    tr, te = fold_arr != fid, fold_arr == fid
                    scaler = StandardScaler()
                    m0_tr = scaler.fit_transform(m0s[tr])
                    m0_te = scaler.transform(m0s[te])
                    if probe is not None:
                        probe.append({
                            "pair": pair,
                            "arm": arm,
                            "fold_id": int(fid),
                            "train_ids": [ids[i] for i, flag in enumerate(tr) if flag],
                            "train_m0": [float(m0s[i, 0]) for i, flag in enumerate(tr) if flag],
                            "test_ids": [ids[i] for i, flag in enumerate(te) if flag],
                            "test_m0": [float(m0s[i, 0]) for i, flag in enumerate(te) if flag],
                            "scaler_mean": float(scaler.mean_[0]),
                            "scaler_scale": float(scaler.scale_[0]),
                        })
                    Xt = np.hstack([fps[tr], m0_tr])
                    Xe = np.hstack([fps[te], m0_te])
                    pred[te] = fit_lr(Xt, y[tr]).predict_proba(Xe)[:, 1]
                return pred

            preds = {
                "physchem": oof_physchem(),
                "ecfp4": oof_ecfp(),
                "ecfp4_m0": oof_ecfp_m0(),
            }
            support = m0s.ravel()
            incomplete = []
            for name, pred in preds.items():
                if pred.size != len(chem) or not np.isfinite(pred).all():
                    incomplete.append(name)
            if incomplete:
                rec["incomplete_reason"] = "INCOMPLETE_OOF:" + ",".join(incomplete)
                units.append(rec)
                continue
            if len(set(y.tolist())) < 2:
                rec["blocking_reason"] = "OOF_SINGLE_CLASS"
                units.append(rec)
                continue
            rec["auroc_physchem"] = auroc(preds["physchem"][y == 1], preds["physchem"][y == 0])
            rec["auroc_ecfp4"] = auroc(preds["ecfp4"][y == 1], preds["ecfp4"][y == 0])
            rec["auroc_ecfp4_m0"] = auroc(preds["ecfp4_m0"][y == 1], preds["ecfp4_m0"][y == 0])
            rec["auroc_SUPPORTING_CHEMISTRY_MATCHED_M0"] = auroc(support[y == 1], support[y == 0])
            rec["delta_ecfp4_plus_m0"] = float(rec["auroc_ecfp4_m0"]) - float(rec["auroc_ecfp4"])
            units.append(rec)
            for i, lig in enumerate(ids):
                oof_rows.append({
                    "pair": pair, "arm": arm, "canonical_ligand_id": lig, "class": classes[i],
                    "fold_id": int(fold_arr[i]),
                    "oof_physchem": float(preds["physchem"][i]),
                    "oof_ecfp4": float(preds["ecfp4"][i]),
                    "oof_ecfp4_m0": float(preds["ecfp4_m0"][i]),
                    "relevant_m0": float(support[i]),
                })
    return units, oof_rows


def module_e(tables) -> tuple[list[dict], list[dict]]:
    pop = tables["population"]
    master = tables["fiveseed"]
    m0 = tables["m0_master"]
    pdb_side = {(r["pair"], r["pdb_id"]): r["target_side"] for r in m0}
    seed_rows = []
    by_pair_metric = defaultdict(lambda: {seed: NA for seed in SEEDS})
    for seed in SEEDS:
        by_key = {}
        for row in master:
            if int(row["seed"]) != int(seed):
                continue
            side = pdb_side.get((row["pair"], row["pdb_id"]))
            if side is None:
                continue
            by_key[(row["pair"], row["canonical_ligand_id"], side)] = row
        for pair in sorted({r["pair"] for r in pop}):
            rebuilt, missing_jobs, missing_ligs = [], [], []
            for row in _rows_for_pair(pop, pair):
                item = dict(row)
                for side, valid, score in (("A", "M0_A_valid", "M0_score_A"), ("B", "M0_B_valid", "M0_score_B")):
                    src = by_key.get((pair, row["canonical_ligand_id"], side))
                    if src and src["status"] == "SUCCESS" and finite(src.get("mode1_affinity")):
                        item[valid] = "1"
                        item[score] = str(-float(src["mode1_affinity"]))
                    else:
                        item[valid] = "0"
                        item[score] = ""
                        if src:
                            missing_jobs.append(src.get("job_id") or f"{pair}:{row['canonical_ligand_id']}:{side}:{seed}")
                            missing_ligs.append(f"{row['canonical_ligand_id']}:{side}:{src.get('pdb_id','')}:{src.get('status','')}")
                        else:
                            missing_jobs.append(f"NO_ROW:{pair}:{row['canonical_ligand_id']}:{side}:{seed}")
                            missing_ligs.append(f"{row['canonical_ligand_id']}:{side}:NO_ROW")
                rebuilt.append(item)
            uni = shared_dual_universe(rebuilt, "M0")
            rec = empty_record(
                SCHEMA_E_SEED, module="E", pair=pair, seed=seed,
                n_dual=uni["n_dual"], n_A_only=uni["n_A_only"], n_B_only=uni["n_B_only"],
                missing_jobs=";".join(missing_jobs), missing_ligands=";".join(missing_ligs),
                reason="", bootstrap_started="NO",
            )
            if uni["dual_shared"] and uni["A_only_for_AUC_B"]:
                rec["auc_B"] = auroc(scores(uni["dual_shared"], "M0_score_B"), scores(uni["A_only_for_AUC_B"], "M0_score_B"))
            else:
                rec["reason"] = "INSUFFICIENT_AUC_B"
            if uni["dual_shared"] and uni["B_only_for_AUC_A"]:
                rec["auc_A"] = auroc(scores(uni["dual_shared"], "M0_score_A"), scores(uni["B_only_for_AUC_A"], "M0_score_A"))
            else:
                rec["reason"] = (rec["reason"] + ";INSUFFICIENT_AUC_A").strip(";")
            if finite(rec["auc_B"]) and finite(rec["auc_A"]):
                rec["summary_min"] = min(float(rec["auc_B"]), float(rec["auc_A"]))
            seed_rows.append(rec)
            for metric in ("auc_B", "auc_A", "summary_min"):
                by_pair_metric[(pair, metric)][seed] = rec[metric]
    summaries = []
    for pair in sorted({r["pair"] for r in pop}):
        for metric in ("auc_B", "auc_A", "summary_min"):
            vals = by_pair_metric[(pair, metric)]
            rec = empty_record(SCHEMA_E_SUM, module="E", pair=pair, metric=metric, seed42=vals[42])
            finite_vals = [float(vals[s]) for s in SEEDS if finite(vals[s])]
            rec["n_finite_seeds"] = len(finite_vals)
            if len(finite_vals) != 5:
                rec["reason"] = f"INCOMPLETE_SEEDS:{len(finite_vals)}/5"
            else:
                ordered = [float(vals[s]) for s in SEEDS]
                med, mn, mx, iqr = five_seed_iqr(ordered)
                rec.update({"median": med, "min": mn, "max": mx, "iqr": iqr, "reason": ""})
            summaries.append(rec)
    return seed_rows, summaries


def module_f(tables, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> list[dict]:
    pop = tables["population"]
    alt = tables["alt_independent"]
    by_alt = defaultdict(list)
    for row in alt:
        by_alt[(row["pair"], row["alt_id"], row["side"])].append(row)
    out = []
    for pair in ALT_PAIRS:
        rows = _rows_for_pair(pop, pair)
        if not rows:
            continue
        for (p, alt_id, side), items in sorted(by_alt.items()):
            if p != pair:
                continue
            lookup = {r["canonical_ligand_id"]: r for r in items}
            common, exclusions = [], []
            for row in rows:
                if row.get("activity_eligible") != "1" or row.get("class") not in {"dual", "A_only", "B_only", "neither"}:
                    continue
                if row.get("M0_A_valid") != "1" or row.get("M0_B_valid") != "1":
                    exclusions.append(f"{row['canonical_ligand_id']}:PRIMARY_NOT_BOTH_FINITE")
                    continue
                src = lookup.get(row["canonical_ligand_id"])
                if src is None:
                    exclusions.append(f"{row['canonical_ligand_id']}:NO_ALT_ROW")
                    continue
                if src.get("status_alt") != "SUCCESS" or not finite(src.get("score_alt")):
                    exclusions.append(f"{row['canonical_ligand_id']}:ALT_NOT_VALID:{src.get('status_alt','')}")
                    continue
                common.append(row)
            pri_rows, alt_rows = [], []
            for row in common:
                src = lookup[row["canonical_ligand_id"]]
                pri = dict(row)
                alt_r = dict(row)
                if side == "A":
                    alt_r["M0_score_A"] = src["score_alt"]
                    alt_r["M0_A_valid"] = "1"
                    alt_r["M0_score_B"] = row["M0_score_B"]
                    alt_r["M0_B_valid"] = "1"
                else:
                    alt_r["M0_score_B"] = src["score_alt"]
                    alt_r["M0_B_valid"] = "1"
                    alt_r["M0_score_A"] = row["M0_score_A"]
                    alt_r["M0_A_valid"] = "1"
                pri_rows.append(pri)
                alt_rows.append(alt_r)
            uni = shared_dual_universe(pri_rows, "M0")
            rec = empty_record(
                SCHEMA_F, module="F", pair=pair, alt_id=alt_id, replaced_side=side,
                n_dual_common=uni["n_dual"], n_A_only_common=uni["n_A_only"], n_B_only_common=uni["n_B_only"],
                common_dual_ids=";".join(sorted(r["canonical_ligand_id"] for r in uni["dual_shared"])),
                common_A_only_ids=";".join(sorted(r["canonical_ligand_id"] for r in uni["A_only_for_AUC_B"])),
                common_B_only_ids=";".join(sorted(r["canonical_ligand_id"] for r in uni["B_only_for_AUC_A"])),
                exclusion_reasons=";".join(exclusions),
                bootstrap_started="NO", n_replicates=0,
            )
            if not (uni["dual_shared"] and uni["A_only_for_AUC_B"] and uni["B_only_for_AUC_A"]):
                rec["reason"] = "INSUFFICIENT_COMMON_SIDES"
                out.append(rec)
                continue
            packed = paired_alt_bootstrap(pri_rows, alt_rows, side, n_boot=n_boot, seed=seed)
            rec.update({k: packed[k] for k in packed})
            out.append(rec)
    return out


def module_h(pop) -> list[dict]:
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        for method in METHODS:
            for cls in ("dual", "A_only", "B_only", "neither"):
                subset = [r for r in rows if r.get("activity_eligible") == "1" and r["class"] == cls]
                for side in ("A", "B"):
                    valid = sum(1 for r in subset if r.get(f"{method}_{side}_valid") == "1")
                    reasons = sorted({r.get(f"{method}_{side}_missing") or "" for r in subset if r.get(f"{method}_{side}_valid") != "1" and r.get(f"{method}_{side}_missing")})
                    out.append({
                        "module": "H", "pair": pair, "method": method, "class": cls, "target": side,
                        "expected": len(subset), "valid": valid, "missing": len(subset) - valid,
                        "missing_reason": ";".join(reasons),
                    })
    return out


def _desc_record(pair, method, kind, ranked) -> dict:
    top = topk_composition(ranked)
    rec = empty_record(SCHEMA_I_DESC, module="I", table="descriptive", pair=pair, method=method, population_kind=kind)
    rec.update(top)
    return rec


def _lig_records(pair, method, kind, ranked) -> list[dict]:
    k = int(topk_composition(ranked)["k"] or 0)
    rows = []
    for i, row in enumerate(ranked):
        rows.append({
            "module": "I", "table": "ligand_ranks", "pair": pair, "method": method,
            "population_kind": kind, "canonical_ligand_id": row["canonical_ligand_id"],
            "global_ligand_entity_id": row["global_ligand_entity_id"], "class": row["class"],
            "rank_A": row["rank_A"], "rank_B": row["rank_B"], "worst_rank": row["worst_rank"],
            "joint_order": row["joint_order"], "in_topk": "YES" if i < k else "NO",
        })
    return rows


def module_i(pop, primary_delta) -> tuple[list[dict], list[dict], list[dict]]:
    delta_lookup = {(r["pair"], r["delta"]): r for r in (primary_delta or [])}
    desc, ligs, cmps = [], [], []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        for method in METHODS:
            pool = ranking_pool(rows, method)
            ranked = joint_rank_rows(pool, method)
            desc.append(_desc_record(pair, method, "method_specific", ranked))
            ligs.extend(_lig_records(pair, method, "method_specific", ranked))
        link_mem = m1b_link_membership(rows)
        for method in ("M1", "M1b", "M2", "M3"):
            common = pairwise_ranking_pool(rows, "M0", method)
            ranked_m0 = joint_rank_rows(common, "M0")
            ranked_mj = joint_rank_rows(common, method)
            desc.append(_desc_record(pair, "M0", f"pairwise_common_vs_{method}", ranked_m0))
            desc.append(_desc_record(pair, method, f"pairwise_common_vs_M0", ranked_mj))
            ligs.extend(_lig_records(pair, "M0", f"pairwise_common_vs_{method}", ranked_m0))
            ligs.extend(_lig_records(pair, method, "pairwise_common_vs_M0", ranked_mj))
            top_m0 = topk_composition(ranked_m0)
            top_mj = topk_composition(ranked_mj)
            rec = empty_record(SCHEMA_I_CMP, module="I", table="method_compare", pair=pair, method=method,
                               n_rank_common=len(common), supporting_only="NO", bootstrap_ci_added="NO")
            dir_m0 = four_sided_membership_sets(rows, method, "M0")
            rank_keys = member_keys(common)
            dir_keys = dir_m0["dual"] | dir_m0["A_only"] | dir_m0["B_only"]
            rec["n_directional_common_dual"] = len(dir_m0["dual"])
            only_rank = sorted(f"{a}:{b}:{c}" for a, b, c in (rank_keys - dir_keys))
            only_dir = sorted(f"{a}:{b}:{c}" for a, b, c in (dir_keys - rank_keys))
            rec["rank_vs_directional_member_diff"] = (
                "" if not (only_rank or only_dir)
                else "rank_only=" + ",".join(only_rank) + ";dir_only=" + ",".join(only_dir)
            )
            if method == "M1b":
                rec["supporting_only"] = "YES"
                rec["bootstrap_ci_added"] = "NO"
                if not link_mem["identical"]:
                    rec["link_source"] = "STOPPED"
                    rec["link_reason"] = "M1b_LINK_MEMBERSHIP_NE:" + ";".join(link_mem["diffs"])
                    rec["delta_summary_min"] = NA
                    rec["delta_top10_single_target_fraction"] = NA
                    rec["concordance"] = NA
                    cmps.append(rec)
                    continue
                d1 = delta_lookup.get((pair, "M1-M0"))
                d2 = delta_lookup.get((pair, "M1b-M1"))
                if not d1 or not d2 or not finite(d1.get("delta_summary_min")) or not finite(d2.get("delta_summary_min")):
                    rec["link_source"] = "STOPPED"
                    rec["link_reason"] = "MISSING_PRIMARY_PAIRED_POINT"
                    cmps.append(rec)
                    continue
                rec["delta_summary_min"] = float(d1["delta_summary_min"]) + float(d2["delta_summary_min"])
                rec["link_source"] = "PRIMARY_POINT_SUM_M1b-M1_PLUS_M1-M0"
                rec["link_reason"] = ""
            else:
                src = delta_lookup.get((pair, f"{method}-M0"))
                if not src or not finite(src.get("delta_summary_min")):
                    rec["link_source"] = "STOPPED"
                    rec["link_reason"] = "MISSING_PRIMARY_PAIRED_POINT"
                    cmps.append(rec)
                    continue
                rec["delta_summary_min"] = float(src["delta_summary_min"])
                rec["link_source"] = f"PRIMARY_{method}-M0"
                rec["link_reason"] = ""
            if finite(top_mj["top10_single_target_fraction"]) and finite(top_m0["top10_single_target_fraction"]):
                rec["delta_top10_single_target_fraction"] = (
                    float(top_mj["top10_single_target_fraction"]) - float(top_m0["top10_single_target_fraction"])
                )
                rec["concordance"] = concordance_code(
                    float(rec["delta_summary_min"]), float(rec["delta_top10_single_target_fraction"])
                )
            cmps.append(rec)
    return desc, ligs, cmps


def run_stage2_modules(tables, outdir: Path, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    pop = tables["population"]
    write_csv(outdir / "module_B_negative_class.csv", module_b(pop, n_boot=n_boot, seed=seed), SCHEMA_B)
    write_csv(outdir / "module_C_wrong_pocket.csv", module_c(pop, n_boot=n_boot, seed=seed), SCHEMA_C)
    d_units, d_oof = module_d(tables)
    write_csv(outdir / "module_D_chemistry_units.csv", d_units, SCHEMA_D_UNIT)
    write_csv(outdir / "module_D_chemistry_oof_ligands.csv", d_oof, SCHEMA_D_OOF)
    write_csv(outdir / "module_D1_composition.csv", module_d1(tables), SCHEMA_D1)
    e_seed, e_sum = module_e(tables)
    write_csv(outdir / "module_E_five_seed.csv", e_seed, SCHEMA_E_SEED)
    write_csv(outdir / "module_E_five_seed_summary.csv", e_sum, SCHEMA_E_SUM)
    write_csv(outdir / "module_F_alt_receptor.csv", module_f(tables, n_boot=n_boot, seed=seed), SCHEMA_F)
    write_csv(outdir / "module_H_missingness.csv", module_h(pop), SCHEMA_H)
    i_desc, i_lig, i_cmp = module_i(pop, tables.get("primary_delta") or [])
    write_csv(outdir / "module_I_descriptive.csv", i_desc, SCHEMA_I_DESC)
    write_csv(outdir / "module_I_ligand_ranks.csv", i_lig, SCHEMA_I_LIG)
    write_csv(outdir / "module_I_method_compare.csv", i_cmp, SCHEMA_I_CMP)
    manifest = {
        "executable_modules": EXECUTABLE_MODULES,
        "blocked_modules": BLOCKED_MODULES,
        "module_statistics": {
            "B": {"bootstrap_B": n_boot, "seed": seed},
            "C": {"bootstrap_B": n_boot, "seed": seed},
            "D": {"bootstrap_B": None, "note": "no chemistry bootstrap"},
            "E": {"bootstrap_B": None, "note": "five-seed point estimates only; IQR is not bootstrap"},
            "F": {"bootstrap_B": n_boot, "seed": seed},
            "H": {"bootstrap_B": None, "note": "counts only"},
            "I": {"bootstrap_B": None, "note": "ranking point estimates; directional deltas read from PRIMARY"},
        },
        "primary_overwritten": False,
    }
    (outdir / "RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def execute(project_root: Path) -> None:
    outdir = project_root / OUTPUT_REL
    if outdir.exists():
        raise SystemExit(f"refuse: output dir already exists: {outdir}")
    tables = load_stage2_tables(project_root)
    missing = [k for k in ("population", "mapping", "alias", "folds", "fiveseed", "alt_independent", "m0_master", "primary_delta") if not tables.get(k)]
    if missing:
        raise SystemExit(f"missing tables: {missing}")
    run_stage2_modules(tables, outdir)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--execute-compute", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    root = resolve_project_root(args.project_root)
    if args.dry_run and not args.execute_compute:
        print("DRY_RUN_OK")
        print("SECOND_STAGE_COMPUTE_EXECUTED=NO")
        print(f"would_write {root / OUTPUT_REL}")
        print("blocked", ",".join(BLOCKED_MODULES))
        return 0
    if not args.execute_compute:
        print("REFUSED: pass --execute-compute only after a separate compute authorization")
        print("SECOND_STAGE_COMPUTE_EXECUTED=NO")
        return 2
    execute(root)
    print("SECOND_STAGE_COMPUTE_EXECUTED=YES")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
