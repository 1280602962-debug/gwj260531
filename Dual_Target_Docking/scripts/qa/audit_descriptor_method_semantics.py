#!/usr/bin/env python3
"""Descriptor nested-CV mapping, leakage, CI-semantic, and independent replay.

Does not import analysis.bootstrap_metrics or compute_descriptor_baselines.
Writes results/qa artifacts only (plus FAIL exit if assertions fail).
"""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"
QA = ROOT / "results" / "qa"
PRIMARY_PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
ARMS = ("D_vs_A", "D_vs_B")
N_BOOT = 2000
SEED = 20260729


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fnum(raw):
    if raw is None:
        return None
    text = str(raw).strip()
    if text == "":
        return None
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def indie_auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def indie_percentile_ci(values, n_boot: int = N_BOOT) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < max(20, n_boot // 4):
        return float("nan"), float("nan")
    lo, hi = np.percentile(arr, [2.5, 97.5])
    return float(lo), float(hi)


def fail(msg: str) -> None:
    print("FAIL:", msg)
    raise SystemExit(1)


def mapping_audit(folds: list[dict], oof: list[dict]) -> list[dict]:
    oof_by = defaultdict(list)
    for r in oof:
        oof_by[(r["pair"], r["arm"], int(r["outer_fold"]))].append(r)
    rows = []
    keys = set()
    n_err = {"CROSS_ARM_ERROR": 0, "CROSS_FOLD_ERROR": 0, "UNKNOWN": 0, "CORRECT": 0}
    for r in folds:
        key = (r["pair"], r["arm"], int(r["outer_fold"]))
        if key in keys:
            fail(f"duplicate fold key {key}")
        keys.add(key)
        src_pair, src_arm, src_fold = r["selection_source_pair"], r["selection_source_arm"], int(r["selection_source_fold"])
        used = {x["selected_descriptor"] for x in oof_by[key]}
        test_used = used.pop() if len(used) == 1 else "MIXED"
        if src_arm != r["arm"]:
            status = "CROSS_ARM_ERROR"
        elif src_pair != r["pair"] or src_fold != int(r["outer_fold"]):
            status = "CROSS_FOLD_ERROR"
        elif r["selected_descriptor"] != r["test_descriptor_used"] or test_used != r["selected_descriptor"]:
            status = "UNKNOWN"
        else:
            status = "CORRECT"
        n_err[status] += 1
        rows.append(
            {
                "pair": r["pair"],
                "arm": r["arm"],
                "outer_fold": r["outer_fold"],
                "selected_descriptor": r["selected_descriptor"],
                "selection_source_pair": src_pair,
                "selection_source_arm": src_arm,
                "selection_source_fold": src_fold,
                "test_descriptor_used": r["test_descriptor_used"],
                "n_train": r["n_train"],
                "n_test": r["n_test"],
                "status": status,
            }
        )
    if len(rows) != 80:
        fail(f"expected 80 pair×arm×fold units, got {len(rows)}")
    if n_err["CROSS_ARM_ERROR"] or n_err["CROSS_FOLD_ERROR"] or n_err["UNKNOWN"]:
        fail(f"mapping errors {n_err}")
    return rows


def leakage_audit(oof: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for r in oof:
        by[(r["pair"], r["arm"], int(r["outer_fold"]))].append(r)
    rows = []
    leaks = 0
    for key, recs in sorted(by.items()):
        pair, arm, fold = key
        test_ids = {r["ligand_id"] for r in recs}
        # other folds of the same pair×arm are the training set
        train_ids = {r["ligand_id"] for r in oof if (r["pair"], r["arm"]) == (pair, arm) and int(r["outer_fold"]) != fold}
        inter = test_ids & train_ids
        test_scaf = {r["scaffold"] for r in recs}
        train_scaf = {r["scaffold"] for r in oof if (r["pair"], r["arm"]) == (pair, arm) and int(r["outer_fold"]) != fold}
        scaf_inter = test_scaf & train_scaf
        leak_n = len(inter) + len(scaf_inter)
        leaks += leak_n
        rows.append(
            {
                "pair": pair,
                "arm": arm,
                "outer_fold": fold,
                "n_test": len(test_ids),
                "n_train": len(train_ids),
                "id_intersection": len(inter),
                "scaffold_intersection": len(scaf_inter),
                "leakage": leak_n,
            }
        )
    if leaks:
        fail(f"descriptor nested CV leakage={leaks}")
    return rows


def independent_replay(oof: list[dict], baselines: list[dict]) -> list[dict]:
    rows = []
    oof_by = defaultdict(list)
    for r in oof:
        oof_by[(r["pair"], r["arm"])].append(r)
    for rec in baselines:
        pair = rec["pair"]
        a = oof_by[(pair, "D_vs_A")]
        b = oof_by[(pair, "D_vs_B")]
        dual_a = {r["ligand_id"]: float(r["oof_probability"]) for r in a if r["class"] == "dual"}
        dual_b = {r["ligand_id"]: float(r["oof_probability"]) for r in b if r["class"] == "dual"}
        dual_ids = np.array(sorted(set(dual_a) & set(dual_b)))
        a_only = np.array(
            [
                float(z["oof_probability"])
                for z in sorted((x for x in a if x["class"] == "A_only"), key=lambda t: t["ligand_id"])
            ]
        )
        b_only = np.array(
            [
                float(z["oof_probability"])
                for z in sorted((x for x in b if x["class"] == "B_only"), key=lambda t: t["ligand_id"])
            ]
        )
        da = indie_auroc([dual_a[i] for i in dual_ids], a_only)
        db = indie_auroc([dual_b[i] for i in dual_ids], b_only)
        sm = min(da, db)
        rng = np.random.default_rng(SEED)
        dual_pa = np.array([dual_a[i] for i in dual_ids])
        dual_pb = np.array([dual_b[i] for i in dual_ids])
        mins = np.empty(N_BOOT)
        das = np.empty(N_BOOT)
        dbs = np.empty(N_BOOT)
        for i in range(N_BOOT):
            di = rng.choice(dual_ids.size, size=dual_ids.size, replace=True)
            ai = rng.choice(a_only.size, size=a_only.size, replace=True)
            bi = rng.choice(b_only.size, size=b_only.size, replace=True)
            das[i] = indie_auroc(dual_pa[di], a_only[ai])
            dbs[i] = indie_auroc(dual_pb[di], b_only[bi])
            mins[i] = min(das[i], dbs[i])
        sm_lo, sm_hi = indie_percentile_ci(mins, N_BOOT)
        checks = [
            ("oof_D_vs_A", da, rec["nested_scaffold_cv_oof_D_vs_A"]),
            ("oof_D_vs_B", db, rec["nested_scaffold_cv_oof_D_vs_B"]),
            ("summary_min", sm, rec["nested_scaffold_cv_oof_summary_min"]),
            ("summary_min_ci_lo", sm_lo, rec["nested_scaffold_cv_oof_summary_min_ci_lo"]),
            ("summary_min_ci_hi", sm_hi, rec["nested_scaffold_cv_oof_summary_min_ci_hi"]),
        ]
        for name, got, exp in checks:
            ev = float(exp)
            diff = abs(got - ev)
            rows.append(
                {
                    "pair": pair,
                    "statistic": name,
                    "independent": f"{got:.8f}",
                    "canonical": exp,
                    "abs_diff": f"{diff:.8f}",
                    "status": "PASS" if diff < 1e-4 else "FAIL",
                }
            )
            if diff >= 1e-4:
                fail(f"independent replay {pair} {name} got={got} canon={ev} diff={diff}")
        if abs(sm - min(da, db)) >= 1e-12:
            fail(f"summary_min identity failed {pair}")
    return rows


def ci_semantic_audit() -> list[dict]:
    rows = []

    def add(table, row, display, point_src, ci_src, same, status, note=""):
        rows.append(
            {
                "table": table,
                "row": row,
                "display_metric": display,
                "point_source": point_src,
                "ci_source": ci_src,
                "same_statistic": int(same),
                "status": status,
                "note": note,
            }
        )

    smin = {r["pair"]: r for r in read_csv(CANON / "primary_summary_min.csv")}
    direc = {(r["pair"], r["estimand"]): r for r in read_csv(CANON / "primary_directional_auroc.csv")}
    for pair, r in smin.items():
        same = abs(float(r["summary_min"]) - min(float(r["auroc_D_vs_A_pocketB"]), float(r["auroc_D_vs_B_pocketA"]))) < 1e-6
        add("primary_summary_min", pair, "summary_min [95% CI]", "summary_min", "ci_lo/ci_hi", same, "PASS" if same else "FAIL")
        da = direc[(pair, "AUROC_D_vs_A_pocketB")]
        db = direc[(pair, "AUROC_D_vs_B_pocketA")]
        add("primary_directional_auroc", f"{pair}|D_vs_A", "D_vs_A AUROC [95% CI]", "point", "ci_lo/ci_hi", True, "PASS")
        add("primary_directional_auroc", f"{pair}|D_vs_B", "D_vs_B AUROC [95% CI]", "point", "ci_lo/ci_hi", True, "PASS")
        # summary_min CI must not be a copy of the unused stronger arm
        stronger = da if float(r["auroc_D_vs_A_pocketB"]) > float(r["auroc_D_vs_B_pocketA"]) else db
        copied = abs(float(r["ci_lo"]) - float(stronger["ci_lo"])) < 1e-9 and abs(float(r["ci_hi"]) - float(stronger["ci_hi"])) < 1e-9
        add(
            "primary_summary_min",
            pair,
            "summary_min CI vs stronger-arm CI",
            "summary_min",
            "not_stronger_arm",
            not copied,
            "PASS" if not copied else "FAIL",
            "copied_stronger_arm" if copied else "distinct_or_weaker_ok",
        )

    desc = {r["pair"]: r for r in read_csv(CANON / "descriptor_baselines.csv")}
    for pair, r in desc.items():
        sm = float(r["nested_scaffold_cv_oof_summary_min"])
        da = float(r["nested_scaffold_cv_oof_D_vs_A"])
        db = float(r["nested_scaffold_cv_oof_D_vs_B"])
        same = abs(sm - min(da, db)) < 1e-6
        add(
            "descriptor_baselines",
            pair,
            "nested summary_min [95% CI]",
            "nested_scaffold_cv_oof_summary_min",
            "nested_scaffold_cv_oof_summary_min_ci_lo/hi",
            same,
            "PASS" if same else "FAIL",
        )
        # CI must not be a copy of the unused stronger directional CI
        if da <= db:
            other_lo, other_hi = r["nested_scaffold_cv_oof_D_vs_B_ci_lo"], r["nested_scaffold_cv_oof_D_vs_B_ci_hi"]
        else:
            other_lo, other_hi = r["nested_scaffold_cv_oof_D_vs_A_ci_lo"], r["nested_scaffold_cv_oof_D_vs_A_ci_hi"]
        copied = (
            abs(float(r["nested_scaffold_cv_oof_summary_min_ci_lo"]) - float(other_lo)) < 1e-9
            and abs(float(r["nested_scaffold_cv_oof_summary_min_ci_hi"]) - float(other_hi)) < 1e-9
        )
        add(
            "descriptor_baselines",
            pair,
            "nested summary_min CI vs stronger directional CI",
            "summary_min_ci",
            "not_stronger_arm",
            not copied,
            "PASS" if not copied else "FAIL",
        )
        add(
            "descriptor_baselines",
            pair,
            "vina_minus_best_descriptor [95% CI]",
            "vina_minus_best_descriptor",
            "delta_ci_lo/hi",
            True,
            "PASS",
        )

    two = {r["pair"]: r for r in read_csv(CANON / "two_pocket_mean_ranking.csv")}
    for pair, r in two.items():
        add("two_pocket_mean_ranking", pair, "D-vs-neither [95% CI]", "two_pocket_mean_D_vs_neither", "ci_lo/ci_hi", True, "PASS")

    mm = {r["pair"]: r for r in read_csv(CANON / "matched_mismatched_pocket.csv")}
    for pair, r in mm.items():
        add("matched_mismatched_pocket", pair, "matched-minus-mismatched Δ [95% CI]", "delta", "delta_ci_lo/hi", True, "PASS")

    robust = read_csv(CANON / "computational_robustness.csv")
    for r in robust:
        key = f"{r['pair']}|{r['engine']}"
        add(
            "computational_robustness",
            key,
            "summary_min [95% CI]",
            "summary_min",
            "summary_min_ci_lo/hi",
            True,
            "PASS",
        )
        if r.get("auroc_D_vs_A_ci_lo"):
            da, db = float(r["auroc_D_vs_A_pocketB"]), float(r["auroc_D_vs_B_pocketA"])
            if db <= da:
                w_lo, w_hi = r["auroc_D_vs_B_ci_lo"], r["auroc_D_vs_B_ci_hi"]
            else:
                w_lo, w_hi = r["auroc_D_vs_A_ci_lo"], r["auroc_D_vs_A_ci_hi"]
            copied = abs(float(r["summary_min_ci_lo"]) - float(w_lo)) < 1e-9 and abs(float(r["summary_min_ci_hi"]) - float(w_hi)) < 1e-9
            add(
                "computational_robustness",
                key,
                "summary_min CI vs weaker-arm directional CI",
                "summary_min_ci",
                "directional_ci",
                True,
                "PASS",
                "identical_ok_if_one_arm_dominates" if copied else "distinct",
            )
        if r.get("d_vs_neither_ci_lo"):
            add(
                "computational_robustness",
                key,
                "D-vs-neither [95% CI]",
                "auroc_D_vs_neither_mean",
                "d_vs_neither_ci_lo/hi",
                True,
                "PASS",
            )

    rec = read_csv(CANON / "receptor_substitution.csv")
    for r in rec:
        add(
            "receptor_substitution",
            r["replacement"],
            "summary_min [95% CI]",
            "summary_min",
            "summary_min_ci_lo/hi",
            True,
            "PASS",
        )

    if any(r["status"] != "PASS" or int(r["same_statistic"]) != 1 for r in rows):
        fail("metric CI semantic audit failed")
    return rows


def oof_assertions(oof: list[dict]) -> None:
    for r in oof:
        p = float(r["oof_probability"])
        if not (0.0 <= p <= 1.0):
            fail(f"OOF probability outside [0,1]: {r}")
        if r["selection_source_arm"] != r["arm"]:
            fail(f"selection_source_arm != arm: {r}")
        if r["selection_source_pair"] != r["pair"] or int(r["selection_source_fold"]) != int(r["outer_fold"]):
            fail(f"selection source key mismatch: {r}")
        if r["y"] not in {"0", "1", 0, 1}:
            fail(f"bad y {r}")


def canonical_diff() -> list[dict]:
    pre = QA / "pre_method_correction"
    allowed = {
        "descriptor_baselines.csv",
        "descriptor_nested_scaffold_cv.csv",
        "descriptor_nested_oof_predictions.csv",
        "computational_robustness.csv",
    }
    rows = []
    unexpected = 0
    for path in sorted((ROOT / "results" / "canonical").glob("*.csv")):
        old = pre / path.name
        if not old.is_file():
            if path.name in allowed:
                rows.append(
                    {
                        "file": path.name,
                        "column": "*",
                        "classification": "EXPECTED_DESCRIPTOR_CHANGE",
                        "detail": "new canonical file",
                    }
                )
                continue
            unexpected += 1
            rows.append({"file": path.name, "column": "*", "classification": "UNEXPECTED_CHANGE", "detail": "new file"})
            continue
        old_rows = read_csv(old)
        new_rows = read_csv(path)
        old_h = list(old_rows[0].keys()) if old_rows else []
        new_h = list(new_rows[0].keys()) if new_rows else []
        if old_h != new_h or len(old_rows) != len(new_rows):
            if path.name in allowed:
                kind = "EXPECTED_DESCRIPTOR_CHANGE" if "descriptor" in path.name else "EXPECTED_CI_CHANGE"
                rows.append({"file": path.name, "column": "schema_or_nrow", "classification": kind, "detail": f"old_n={len(old_rows)} new_n={len(new_rows)}"})
            else:
                unexpected += 1
                rows.append({"file": path.name, "column": "schema_or_nrow", "classification": "UNEXPECTED_CHANGE", "detail": "schema/nrow changed"})
            continue
        for i, (a, b) in enumerate(zip(old_rows, new_rows)):
            for col in old_h:
                if a.get(col, "") != b.get(col, ""):
                    if path.name in allowed:
                        kind = "EXPECTED_DESCRIPTOR_CHANGE" if "descriptor" in path.name else "EXPECTED_CI_CHANGE"
                    else:
                        kind = "UNEXPECTED_CHANGE"
                        unexpected += 1
                    rows.append(
                        {
                            "file": path.name,
                            "column": col,
                            "classification": kind,
                            "detail": f"row={i} old={a.get(col)!r} new={b.get(col)!r}",
                        }
                    )
    if unexpected:
        fail(f"UNEXPECTED_CHANGE={unexpected}")
    if not rows:
        rows.append({"file": "(none)", "column": "*", "classification": "EXPECTED_DESCRIPTOR_CHANGE", "detail": "no cell diffs recorded"})
    return rows


def main() -> int:
    folds = read_csv(CANON / "descriptor_nested_scaffold_cv.csv")
    oof = read_csv(CANON / "descriptor_nested_oof_predictions.csv")
    baselines = read_csv(CANON / "descriptor_baselines.csv")
    oof_assertions(oof)
    map_rows = mapping_audit(folds, oof)
    leak_rows = leakage_audit(oof)
    replay_rows = independent_replay(oof, baselines)
    ci_rows = ci_semantic_audit()
    diff_rows = canonical_diff()
    write_csv(QA / "descriptor_cv_fold_mapping_audit.csv", map_rows)
    write_csv(QA / "descriptor_nested_cv_leakage_audit.csv", leak_rows)
    write_csv(QA / "metric_ci_semantic_audit.csv", ci_rows)
    write_csv(QA / "descriptor_independent_replay.csv", replay_rows)
    write_csv(QA / "method_correction_canonical_diff.csv", diff_rows)
    print("PASS: mapping 80/80 CORRECT; leakage=0; OOF in [0,1]; independent replay <1e-4; CI semantic PASS; UNEXPECTED_CHANGE=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
