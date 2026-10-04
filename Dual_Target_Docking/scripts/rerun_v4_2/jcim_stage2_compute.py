#!/usr/bin/env python3
"""Stage 2 compute entry for modules B, C, D, E, F, H, I.

This script does nothing unless --execute-compute is passed.
This preparation commit must not invoke that flag.

Forbidden: modules G, J, K; PRIMARY overwrite; historical alt stats imports.
"""
from __future__ import annotations

import argparse
import csv
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
    PATHS,
    PHYSCHEM_NAMES,
    SEEDS,
    auroc,
    alias_ids,
    build_layer_members,
    concordance_code,
    fold_keyset,
    four_sided_membership_identical,
    joint_rank_rows,
    load_stage2_tables,
    percentile_ci,
    resolve_project_root,
    scores,
    shared_dual_universe,
    topk_composition,
)

OUTPUT_REL = Path("results/jcim_stage2")


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _rows_for_pair(pop, pair):
    return [r for r in pop if r["pair"] == pair]


def _bootstrap_two_auc(dual, neg, dual_key, neg_key, n_boot=BOOT_B, seed=BOOT_SEED):
    pt = auroc(scores(dual, dual_key), scores(neg, neg_key))
    rng = np_rng(seed)
    reps = []
    for _ in range(n_boot):
        d_s = [dual[i] for i in rng.integers(0, len(dual), size=len(dual))]
        n_s = [neg[i] for i in rng.integers(0, len(neg), size=len(neg))]
        reps.append(auroc(scores(d_s, dual_key), scores(n_s, neg_key)))
    lo, hi = percentile_ci(reps)
    return pt, lo, hi


def np_rng(seed):
    import numpy as np

    return np.random.default_rng(seed)


def module_b(pop) -> list[dict]:
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        uni = shared_dual_universe(rows, "M0")
        neither_a = [
            r for r in rows if r.get("activity_eligible") == "1" and r["class"] == "neither" and r.get("M0_A_valid") == "1"
        ]
        neither_b = [
            r for r in rows if r.get("activity_eligible") == "1" and r["class"] == "neither" and r.get("M0_B_valid") == "1"
        ]
        dual = uni["dual_shared"]
        tasks = (
            ("dual_vs_A_only", dual, uni["A_only_for_AUC_B"], "M0_score_B"),
            ("dual_vs_neither_B", dual, neither_b, "M0_score_B"),
            ("dual_vs_B_only", dual, uni["B_only_for_AUC_A"], "M0_score_A"),
            ("dual_vs_neither_A", dual, neither_a, "M0_score_A"),
        )
        packed = {}
        for name, pos, neg, key in tasks:
            if not pos or not neg:
                packed[name] = {"value": "", "ci_lo": "", "ci_hi": "", "n_neg": len(neg), "reason": "INSUFFICIENT"}
                continue
            pt, lo, hi = _bootstrap_two_auc(pos, neg, key, key)
            packed[name] = {"value": pt, "ci_lo": lo, "ci_hi": hi, "n_neg": len(neg), "reason": ""}
        for direction, single, neither, _key in (
            ("B", "dual_vs_A_only", "dual_vs_neither_B", "M0_score_B"),
            ("A", "dual_vs_B_only", "dual_vs_neither_A", "M0_score_A"),
        ):
            s, n = packed[single], packed[neither]
            delta = ""
            if s["value"] != "" and n["value"] != "":
                delta = float(s["value"]) - float(n["value"])
            out.append(
                {
                    "module": "B",
                    "pair": pair,
                    "direction": direction,
                    "n_dual": uni["n_dual"],
                    "n_single": s["n_neg"],
                    "n_neither": n["n_neg"],
                    "auc_single": s["value"],
                    "auc_neither": n["value"],
                    "delta_negative": delta,
                    "auc_single_ci_lo": s["ci_lo"],
                    "auc_single_ci_hi": s["ci_hi"],
                    "auc_neither_ci_lo": n["ci_lo"],
                    "auc_neither_ci_hi": n["ci_hi"],
                }
            )
    return out


def module_c(pop) -> list[dict]:
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        dual = [
            r
            for r in rows
            if r.get("activity_eligible") == "1"
            and r["class"] == "dual"
            and r.get("M0_A_valid") == "1"
            and r.get("M0_B_valid") == "1"
        ]
        for arm, neg_cls, correct, wrong in (
            ("D_vs_A", "A_only", "M0_score_B", "M0_score_A"),
            ("D_vs_B", "B_only", "M0_score_A", "M0_score_B"),
        ):
            negs = [
                r
                for r in rows
                if r.get("activity_eligible") == "1"
                and r["class"] == neg_cls
                and r.get("M0_A_valid") == "1"
                and r.get("M0_B_valid") == "1"
            ]
            rec = {
                "module": "C",
                "pair": pair,
                "arm": arm,
                "n_dual": len(dual),
                "n_negative_both_finite": len(negs),
                "auc_correct": "",
                "auc_wrong": "",
                "delta_pocket": "",
                "extra_drop_vs_directional_m0": "",
            }
            directional_negs = [
                r
                for r in rows
                if r.get("activity_eligible") == "1"
                and r["class"] == neg_cls
                and (r.get("M0_B_valid") == "1" if arm == "D_vs_A" else r.get("M0_A_valid") == "1")
            ]
            rec["extra_drop_vs_directional_m0"] = len(directional_negs) - len(negs)
            if dual and negs:
                pt_c, lo_c, hi_c = _bootstrap_two_auc(dual, negs, correct, correct)
                rng = np_rng(BOOT_SEED)
                wrong_reps = []
                for _ in range(BOOT_B):
                    d_s = [dual[i] for i in rng.integers(0, len(dual), size=len(dual))]
                    n_s = [negs[i] for i in rng.integers(0, len(negs), size=len(negs))]
                    wrong_reps.append(auroc(scores(d_s, wrong), scores(n_s, wrong)))
                pt_w = auroc(scores(dual, wrong), scores(negs, wrong))
                lo_w, hi_w = percentile_ci(wrong_reps)
                rec.update(
                    {
                        "auc_correct": pt_c,
                        "auc_correct_ci_lo": lo_c,
                        "auc_correct_ci_hi": hi_c,
                        "auc_wrong": pt_w,
                        "auc_wrong_ci_lo": lo_w,
                        "auc_wrong_ci_hi": hi_w,
                        "delta_pocket": pt_c - pt_w,
                    }
                )
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


def _ecfp(mol):
    from rdkit.Chem import AllChem
    import numpy as np

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, ECFP_RADIUS, nBits=ECFP_NBITS, useChirality=False, useFeatures=False)
    arr = np.zeros((ECFP_NBITS,), dtype=float)
    for idx in fp.GetOnBits():
        arr[idx] = 1.0
    return arr


def module_d(tables) -> list[dict]:
    from rdkit import Chem
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    pop = tables["population"]
    aliases = alias_ids(tables["alias"])
    fold_keys = fold_keyset(tables["folds"])
    fold_of = {(r["pair"], r["arm"], r["ligand_id"]): r["fold_id"] for r in tables["folds"]}
    smiles = {r["canonical_ligand_id"]: r.get("canonical_smiles") or "" for r in tables["mapping"]}
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        for arm in ("D_vs_A", "D_vs_B"):
            layers = build_layer_members(pop, aliases, fold_keys, pair, arm)
            chem = layers["chemistry_oof"]
            rec = {
                "module": "D",
                "pair": pair,
                "arm": arm,
                "n_labeled_population": len(layers["labeled"]),
                "n_m0_directional_population": len(layers["m0"]),
                "n_final_chemistry_oof": len(chem),
                "blocking_reason": "",
            }
            if not chem:
                rec["blocking_reason"] = "EMPTY_CHEMISTRY_OOF"
                out.append(rec)
                continue
            y, folds, phys, fps, m0s, ids = [], [], [], [], [], []
            parse_fail = False
            for row in chem:
                mol = Chem.MolFromSmiles(smiles.get(row["canonical_ligand_id"], ""))
                if mol is None:
                    parse_fail = True
                    break
                y.append(1 if row["class"] == "dual" else 0)
                folds.append(int(fold_of[(pair, arm, row["canonical_ligand_id"])]))
                phys.append([_physchem(mol)[k] for k in PHYSCHEM_NAMES])
                fps.append(_ecfp(mol))
                key = "M0_score_B" if arm == "D_vs_A" else "M0_score_A"
                m0s.append(float(row[key]))
                ids.append(row["canonical_ligand_id"])
            if parse_fail:
                rec["blocking_reason"] = "CHEMISTRY_SMILES_PARSE_FAILURE"
                out.append(rec)
                continue
            y = np.asarray(y)
            folds = np.asarray(folds)
            phys = np.asarray(phys, dtype=float)
            fps = np.asarray(fps, dtype=float)
            m0s = np.asarray(m0s, dtype=float).reshape(-1, 1)
            fold_ids = sorted(set(folds.tolist()))
            blocked = ""
            for fid in fold_ids:
                train = folds != fid
                if len(set(y[train].tolist())) < 2:
                    blocked = "TRAIN_SINGLE_CLASS"
            if blocked:
                rec["blocking_reason"] = blocked
                out.append(rec)
                continue

            def oof_proba(X, scale=False):
                pred = np.full(len(y), np.nan)
                for fid in fold_ids:
                    tr, te = folds != fid, folds == fid
                    if len(set(y[te].tolist())) < 2:
                        continue
                    Xt, Xe = X[tr], X[te]
                    if scale:
                        scaler = StandardScaler()
                        Xt = scaler.fit_transform(Xt)
                        Xe = scaler.transform(Xe)
                    model = LogisticRegression(
                        penalty="l2",
                        C=LR_C,
                        solver=LR_SOLVER,
                        class_weight=None,
                        max_iter=LR_MAX_ITER,
                        random_state=LR_RANDOM_STATE,
                    )
                    model.fit(Xt, y[tr])
                    pred[te] = model.predict_proba(Xe)[:, 1]
                return pred

            models = {
                "physchem": oof_proba(phys, scale=True),
                "ecfp4": oof_proba(fps, scale=False),
                "ecfp4_m0": oof_proba(np.hstack([fps, m0s]), scale=False),
                "supporting_m0": m0s.ravel(),
            }
            if len(set(y.tolist())) < 2:
                rec["blocking_reason"] = "OOF_SINGLE_CLASS"
                out.append(rec)
                continue
            for name, pred in models.items():
                mask = np.isfinite(pred)
                if mask.sum() == 0 or len(set(y[mask].tolist())) < 2:
                    rec[f"auroc_{name}"] = ""
                    rec[f"{name}_reason"] = "OOF_SINGLE_CLASS" if name != "supporting_m0" else "INSUFFICIENT"
                    continue
                rec[f"auroc_{name}"] = auroc(pred[mask][y[mask] == 1], pred[mask][y[mask] == 0])
            rec["delta_ecfp4_plus_m0"] = (
                float(rec["auroc_ecfp4_m0"]) - float(rec["auroc_ecfp4"])
                if rec.get("auroc_ecfp4_m0") != "" and rec.get("auroc_ecfp4") != ""
                else ""
            )
            out.append(rec)
    return out


def module_e(tables) -> list[dict]:
    pop = tables["population"]
    master = tables["fiveseed"]
    m0 = tables["m0_master"]
    pdb_side = {(r["pair"], r["pdb_id"]): r["target_side"] for r in m0}
    out = []
    for seed in SEEDS:
        by_key = {}
        for row in master:
            if int(row["seed"]) != seed:
                continue
            side = pdb_side.get((row["pair"], row["pdb_id"]))
            if side is None:
                continue
            by_key[(row["pair"], row["canonical_ligand_id"], side)] = row
        for pair in sorted({r["pair"] for r in pop}):
            rebuilt = []
            for row in _rows_for_pair(pop, pair):
                item = dict(row)
                for side, valid, score in (("A", "M0_A_valid", "M0_score_A"), ("B", "M0_B_valid", "M0_score_B")):
                    src = by_key.get((pair, row["canonical_ligand_id"], side))
                    if src and src["status"] == "SUCCESS":
                        item[valid] = "1"
                        item[score] = str(-float(src["mode1_affinity"]))
                    else:
                        item[valid] = "0"
                        item[score] = ""
                rebuilt.append(item)
            uni = shared_dual_universe(rebuilt, "M0")
            rec = {
                "module": "E",
                "pair": pair,
                "seed": seed,
                "n_dual": uni["n_dual"],
                "n_A_only": uni["n_A_only"],
                "n_B_only": uni["n_B_only"],
                "auc_B": "",
                "auc_A": "",
                "summary_min": "",
            }
            if uni["dual_shared"] and uni["A_only_for_AUC_B"]:
                rec["auc_B"] = auroc(scores(uni["dual_shared"], "M0_score_B"), scores(uni["A_only_for_AUC_B"], "M0_score_B"))
            if uni["dual_shared"] and uni["B_only_for_AUC_A"]:
                rec["auc_A"] = auroc(scores(uni["dual_shared"], "M0_score_A"), scores(uni["B_only_for_AUC_A"], "M0_score_A"))
            if rec["auc_B"] != "" and rec["auc_A"] != "":
                rec["summary_min"] = min(float(rec["auc_B"]), float(rec["auc_A"]))
            out.append(rec)
    return out


def module_f(tables) -> list[dict]:
    pop = tables["population"]
    alt = tables["alt_independent"]
    out = []
    by_alt = defaultdict(list)
    for row in alt:
        by_alt[(row["pair"], row["alt_id"], row["side"])].append(row)
    for pair in ALT_PAIRS:
        rows = _rows_for_pair(pop, pair)
        for (p, alt_id, side), items in sorted(by_alt.items()):
            if p != pair:
                continue
            lookup = {(r["canonical_ligand_id"]): r for r in items}
            rebuilt = []
            for row in rows:
                item = dict(row)
                src = lookup.get(row["canonical_ligand_id"])
                if src and src.get("status_alt") == "SUCCESS":
                    if side == "A":
                        item["M0_score_A"] = src["score_alt"]
                        item["M0_A_valid"] = "1"
                    else:
                        item["M0_score_B"] = src["score_alt"]
                        item["M0_B_valid"] = "1"
                elif src:
                    if side == "A":
                        item["M0_A_valid"] = "0"
                        item["M0_score_A"] = ""
                    else:
                        item["M0_B_valid"] = "0"
                        item["M0_score_B"] = ""
                rebuilt.append(item)
            uni_alt = shared_dual_universe(rebuilt, "M0")
            uni_pri = shared_dual_universe(rows, "M0")
            rec = {
                "module": "F",
                "pair": pair,
                "alt_id": alt_id,
                "replaced_side": side,
                "n_dual_alt": uni_alt["n_dual"],
                "n_dual_primary": uni_pri["n_dual"],
            }
            def pack(uni, tag):
                auc_b = auc_a = ""
                if uni["dual_shared"] and uni["A_only_for_AUC_B"]:
                    auc_b = auroc(scores(uni["dual_shared"], "M0_score_B"), scores(uni["A_only_for_AUC_B"], "M0_score_B"))
                if uni["dual_shared"] and uni["B_only_for_AUC_A"]:
                    auc_a = auroc(scores(uni["dual_shared"], "M0_score_A"), scores(uni["B_only_for_AUC_A"], "M0_score_A"))
                rec[f"auc_B_{tag}"] = auc_b
                rec[f"auc_A_{tag}"] = auc_a
                rec[f"summary_min_{tag}"] = min(float(auc_b), float(auc_a)) if auc_b != "" and auc_a != "" else ""

            pack(uni_alt, "alt")
            pack(uni_pri, "primary")
            if rec["auc_B_alt"] != "" and rec["auc_B_primary"] != "":
                rec["delta_AUC_B"] = float(rec["auc_B_alt"]) - float(rec["auc_B_primary"])
            if rec["auc_A_alt"] != "" and rec["auc_A_primary"] != "":
                rec["delta_AUC_A"] = float(rec["auc_A_alt"]) - float(rec["auc_A_primary"])
            rec["delta_AUC_affected"] = rec.get("delta_AUC_B") if side == "B" else rec.get("delta_AUC_A")
            rec["delta_AUC_unreplaced"] = rec.get("delta_AUC_A") if side == "B" else rec.get("delta_AUC_B")
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
                    reasons = [r.get(f"{method}_{side}_missing") or "" for r in subset if r.get(f"{method}_{side}_valid") != "1"]
                    out.append(
                        {
                            "module": "H",
                            "pair": pair,
                            "method": method,
                            "class": cls,
                            "target": side,
                            "expected": len(subset),
                            "valid": valid,
                            "missing": len(subset) - valid,
                            "missing_reason": ";".join(sorted({x for x in reasons if x})),
                        }
                    )
    return out


def module_i(pop) -> list[dict]:
    out = []
    for pair in sorted({r["pair"] for r in pop}):
        rows = _rows_for_pair(pop, pair)
        by_method = {}
        for method in METHODS:
            ranked = joint_rank_rows(rows, method)
            by_method[method] = (ranked, topk_composition(ranked) if ranked else None)
        ident = four_sided_membership_identical(rows, "M1b", "M1", "M0")
        uni = {}
        for method in METHODS:
            u = shared_dual_universe(rows, method)
            sm = ""
            if u["dual_shared"] and u["A_only_for_AUC_B"] and u["B_only_for_AUC_A"]:
                auc_b = auroc(scores(u["dual_shared"], f"{method}_score_B"), scores(u["A_only_for_AUC_B"], f"{method}_score_B"))
                auc_a = auroc(scores(u["dual_shared"], f"{method}_score_A"), scores(u["B_only_for_AUC_A"], f"{method}_score_A"))
                sm = min(auc_b, auc_a)
            uni[method] = sm
        for method in METHODS:
            ranked, top = by_method[method]
            rec = {
                "module": "I",
                "pair": pair,
                "method": method,
                "n_ranked": len(ranked),
                "k": top["k"] if top else "",
                "top10_dual_fraction": top["top10_dual_fraction"] if top else "",
                "top10_single_target_fraction": top["top10_single_target_fraction"] if top else "",
                "dual_recall": top["dual_recall"] if top else "",
                "expected_dual_count": top["expected_dual_count"] if top else "",
                "expected_dual_fraction": top["expected_dual_fraction"] if top else "",
                "summary_min_point": uni[method],
                "m1b_m0_membership_equals_m1_m0": "YES" if ident else "NO",
            }
            if method != "M0" and top and by_method["M0"][1] and uni[method] != "" and uni["M0"] != "":
                d_sm = float(uni[method]) - float(uni["M0"])
                d_stf = float(top["top10_single_target_fraction"]) - float(by_method["M0"][1]["top10_single_target_fraction"])
                rec["delta_summary_min"] = d_sm
                rec["delta_top10_single_target_fraction"] = d_stf
                rec["concordance"] = concordance_code(d_sm, d_stf)
                if method == "M1b":
                    rec["supporting_only"] = "YES"
                    rec["bootstrap_ci_added"] = "NO"
            out.append(rec)
    return out


def execute(project_root: Path) -> None:
    outdir = project_root / OUTPUT_REL
    if outdir.exists():
        raise SystemExit(f"refuse: output dir already exists: {outdir}")
    tables = load_stage2_tables(project_root)
    missing = [k for k in ("population", "mapping", "alias", "folds", "fiveseed", "alt_independent", "m0_master") if not tables[k]]
    if missing:
        raise SystemExit(f"missing tables: {missing}")
    pop = tables["population"]
    _write_csv(outdir / "module_B_negative_class.csv", module_b(pop))
    _write_csv(outdir / "module_C_wrong_pocket.csv", module_c(pop))
    _write_csv(outdir / "module_D_chemistry.csv", module_d(tables))
    _write_csv(outdir / "module_E_five_seed.csv", module_e(tables))
    _write_csv(outdir / "module_F_alt_receptor.csv", module_f(tables))
    _write_csv(outdir / "module_H_missingness.csv", module_h(pop))
    _write_csv(outdir / "module_I_joint_ranking.csv", module_i(pop))
    (outdir / "RUN_MANIFEST.json").write_text(
        json.dumps(
            {
                "executable_modules": EXECUTABLE_MODULES,
                "blocked_modules": BLOCKED_MODULES,
                "bootstrap": {"B": BOOT_B, "seed": BOOT_SEED},
                "primary_overwritten": False,
            },
            indent=2,
        )
        + "\n"
    )


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
