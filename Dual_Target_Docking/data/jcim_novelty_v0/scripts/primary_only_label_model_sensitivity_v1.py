#!/usr/bin/env python3
"""Adjudicated primary-source label and model sensitivity analysis.

This is deliberately *not* called a complete primary-only rebuild.  The
paper-level review table covers the decisive records inspected so far, whereas
the assay audit covers only priority ligands.  The analysis therefore removes
records that have been explicitly adjudicated as review-derived, incorrectly
mapped, or cellular/pathway surrogates and leaves unreviewed records unchanged.

Three label scenarios are reported side by side:
  1. frozen_cached: labels used by the frozen main analysis;
  2. api_max_all_sources: current high-confidence ChEMBL maximum labels;
  3. adjudicated_primary_only: (2) after explicit non-primary exclusions.

The primary Vina readout and the secondary RTMscore readout are recomputed
without redocking.  Missing target arms created by exclusion are unresolved;
they are never silently treated as inactive.
"""
from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_novelty_v0"
TAB = OUT / "tables"
ANA = OUT / "analysis"
THETA = 6.0
N_BOOT = 2000
SEED = 20260828

EXCLUDE_DECISIONS = {
    "cellular_surrogate_only",
    "exclude",
    "exclude_from_biochemical_max",
    "exclude_from_primary_only",
    "secondary_only",
}

SCORE_SPEC = {
    "EGFR/HER2": {
        "path": "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        "vina_A": "3POZ_affinity", "vina_B": "3RCD_affinity",
        "rtmscore_A": "rtm_3POZ", "rtmscore_B": "rtm_3RCD",
    },
    "AChE/BChE": {
        "path": "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "vina_A": "vina_ACHE", "vina_B": "vina_BCHE",
        "rtmscore_A": "rtm_ACHE", "rtmscore_B": "rtm_BCHE",
    },
    "PIK3CA/PIK3CB": {
        "path": "data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        "vina_A": "vina_PIK3CA", "vina_B": "vina_PIK3CB",
        "rtmscore_A": "rtm_PIK3CA", "rtmscore_B": "rtm_PIK3CB",
    },
    "PIK3CA/mTOR": {
        "path": "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        "vina_A": "4L23_affinity", "vina_B": "4JT6_affinity",
        "rtmscore_A": "rtm_4L23", "rtmscore_B": "rtm_4JT6",
    },
}


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value):
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def classify(pa, pb):
    if pa is None or pb is None:
        return "unresolved_missing_arm"
    if pa >= THETA and pb >= THETA:
        return "dual"
    if pa >= THETA:
        return "A_only"
    if pb >= THETA:
        return "B_only"
    return "neither"


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    delta = p[:, None] - n[None, :]
    return float(((delta > 0).sum() + 0.5 * (delta == 0).sum()) / delta.size)


def bootstrap_metrics(records, seed):
    d = [r for r in records if r["class"] == "dual"]
    a = [r for r in records if r["class"] == "A_only"]
    b = [r for r in records if r["class"] == "B_only"]
    auc_da = auroc([r["score_B"] for r in d], [r["score_B"] for r in a])
    auc_db = auroc([r["score_A"] for r in d], [r["score_A"] for r in b])
    if min(len(d), len(a), len(b)) == 0:
        return auc_da, (float("nan"), float("nan")), auc_db, (float("nan"), float("nan")), float("nan"), (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    da_boot, db_boot, sm_boot = [], [], []
    for _ in range(N_BOOT):
        dd = [d[i] for i in rng.integers(0, len(d), len(d))]
        aa = [a[i] for i in rng.integers(0, len(a), len(a))]
        bb = [b[i] for i in rng.integers(0, len(b), len(b))]
        x = auroc([r["score_B"] for r in dd], [r["score_B"] for r in aa])
        y = auroc([r["score_A"] for r in dd], [r["score_A"] for r in bb])
        da_boot.append(x)
        db_boot.append(y)
        sm_boot.append(min(x, y))
    ci = lambda values: tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))
    return auc_da, ci(da_boot), auc_db, ci(db_boot), min(auc_da, auc_db), ci(sm_boot)


def fmt(value):
    return "" if value is None or value != value else round(float(value), 4)


def load_scores():
    scores = {}
    for pair, spec in SCORE_SPEC.items():
        for row in read_csv(ROOT / spec["path"]):
            ligand = row.get("ligand") or row.get("panel_id")
            rec = {}
            for model in ("vina", "rtmscore"):
                a, b = fnum(row.get(spec[f"{model}_A"])), fnum(row.get(spec[f"{model}_B"]))
                if model == "vina":
                    a = None if a is None else -a
                    b = None if b is None else -b
                rec[model] = (a, b)
            scores[(pair, ligand)] = rec
    return scores


def main():
    labels = read_csv(TAB / "high_confidence_labels_v1.csv")
    audit = read_csv(TAB / "assay_context_audit.csv")
    review = read_csv(TAB / "decisive_source_human_review_v1.csv")
    scores = load_scores()

    excluded_keys = defaultdict(set)
    reviewed_keys = set()
    for row in review:
        key = (row["pair"], row["ligand"], row["target"], row["document_chembl_id"])
        reviewed_keys.add(key)
        if row["decision"] in EXCLUDE_DECISIONS:
            excluded_keys[(row["pair"], row["ligand"], row["target"])].add(row["document_chembl_id"])

    audit_by_arm = defaultdict(list)
    excluded_rows = []
    for row in audit:
        arm = (row["pair"], row["ligand"], row["target_name"])
        is_source_excluded = row["document_chembl_id"] in excluded_keys.get(arm, set())
        is_machine_excluded = row.get("human_include_exclude") == "exclude"
        # Source-level include decisions override the earlier coarse assay-type
        # machine flag (e.g., source-verified Ellman assays stored as type A).
        exact_key = (*arm, row["document_chembl_id"])
        source_included = any(
            r["pair"] == arm[0] and r["ligand"] == arm[1] and r["target"] == arm[2]
            and r["document_chembl_id"] == row["document_chembl_id"] and r["decision"] == "include"
            for r in review
        )
        drop = is_source_excluded or (is_machine_excluded and not source_included)
        if drop:
            excluded_rows.append({
                "pair": arm[0], "ligand": arm[1], "target": arm[2],
                "document_chembl_id": row["document_chembl_id"],
                "pchembl_value": row["pchembl_value"],
                "reason": "source_adjudication" if is_source_excluded else "machine_human_exclude",
            })
        else:
            audit_by_arm[arm].append(row)

    label_rows = []
    class_by_scenario = {}
    for row in labels:
        pair, ligand = row["pair"], row["ligand"]
        base_a, base_b = fnum(row["high_conf_max_A"]), fnum(row["high_conf_max_B"])
        adjusted = []
        arm_status = []
        for target, base in ((pair.split("/")[0], base_a), (pair.split("/")[1], base_b)):
            arm = (pair, ligand, target)
            if arm in excluded_keys:
                values = [fnum(x["pchembl_value"]) for x in audit_by_arm.get(arm, [])]
                values = [x for x in values if x is not None]
                adjusted.append(max(values) if values else None)
                arm_status.append("adjudicated_recomputed" if values else "unresolved_after_exclusion")
            else:
                adjusted.append(base)
                arm_status.append("unchanged_not_adjudicated")
        adj_a, adj_b = adjusted
        api_class = classify(base_a, base_b)
        adj_class = classify(adj_a, adj_b)
        frozen_class = row["frozen_class"]
        class_by_scenario[(pair, ligand, "frozen_cached")] = frozen_class
        class_by_scenario[(pair, ligand, "api_max_all_sources")] = api_class
        class_by_scenario[(pair, ligand, "adjudicated_primary_only")] = adj_class
        label_rows.append({
            "pair": pair, "ligand": ligand, "molecule_chembl_id": row["molecule_chembl_id"],
            "frozen_class": frozen_class, "api_max_class": api_class,
            "adjudicated_primary_only_class": adj_class,
            "api_max_A": base_a, "api_max_B": base_b,
            "primary_only_max_A": adj_a, "primary_only_max_B": adj_b,
            "arm_A_status": arm_status[0], "arm_B_status": arm_status[1],
            "flip_vs_api_max": int(adj_class != api_class),
            "flip_vs_frozen": int(adj_class != frozen_class),
            "model_eligible_primary_only": int(adj_class in {"dual", "A_only", "B_only", "neither"}),
        })

    scenarios = ("frozen_cached", "api_max_all_sources", "adjudicated_primary_only")
    model_rows = []
    for pair in SCORE_SPEC:
        for model in ("vina", "rtmscore"):
            for scenario in scenarios:
                records = []
                for row in labels:
                    key = (pair, row["ligand"])
                    ab = scores.get(key, {}).get(model)
                    if row["pair"] != pair or not ab or None in ab:
                        continue
                    cls = class_by_scenario[(pair, row["ligand"], scenario)]
                    records.append({"class": cls, "score_A": ab[0], "score_B": ab[1]})
                digest = hashlib.md5(f"{pair}|{model}|{scenario}".encode()).hexdigest()
                seed = SEED + int(digest[:8], 16) % 100003
                da, da_ci, db, db_ci, sm, sm_ci = bootstrap_metrics(records, seed)
                counts = {c: sum(r["class"] == c for r in records) for c in ("dual", "A_only", "B_only", "neither", "unresolved_missing_arm")}
                model_rows.append({
                    "pair": pair, "model": model, "label_scenario": scenario,
                    "n_scored_total": len(records), "n_dual": counts["dual"],
                    "n_A_only": counts["A_only"], "n_B_only": counts["B_only"],
                    "n_neither": counts["neither"], "n_unresolved": counts["unresolved_missing_arm"],
                    "auroc_D_vs_A_pocketB": fmt(da), "ci_D_vs_A_lo": fmt(da_ci[0]), "ci_D_vs_A_hi": fmt(da_ci[1]),
                    "auroc_D_vs_B_pocketA": fmt(db), "ci_D_vs_B_lo": fmt(db_ci[0]), "ci_D_vs_B_hi": fmt(db_ci[1]),
                    "summary_min": fmt(sm), "summary_min_ci_lo": fmt(sm_ci[0]), "summary_min_ci_hi": fmt(sm_ci[1]),
                    "underpowered": int(min(counts["dual"], counts["A_only"], counts["B_only"]) < 8),
                })

    # Add point deltas against the matching API-max and frozen scenarios.
    index = {(r["pair"], r["model"], r["label_scenario"]): r for r in model_rows}
    for row in model_rows:
        if row["label_scenario"] != "adjudicated_primary_only" or row["summary_min"] == "":
            row["delta_summary_min_vs_api_max"] = ""
            row["delta_summary_min_vs_frozen"] = ""
            continue
        api = index[(row["pair"], row["model"], "api_max_all_sources")]
        frozen = index[(row["pair"], row["model"], "frozen_cached")]
        row["delta_summary_min_vs_api_max"] = fmt(float(row["summary_min"]) - float(api["summary_min"]))
        row["delta_summary_min_vs_frozen"] = fmt(float(row["summary_min"]) - float(frozen["summary_min"]))

    write_csv(TAB / "primary_only_label_sensitivity_v1.csv", label_rows)
    write_csv(TAB / "primary_only_excluded_activity_rows_v1.csv", excluded_rows)
    write_csv(TAB / "primary_only_model_sensitivity_v1.csv", model_rows)

    flips_api = [r for r in label_rows if r["flip_vs_api_max"]]
    unresolved = [r for r in label_rows if r["adjudicated_primary_only_class"] == "unresolved_missing_arm"]
    primary_models = [r for r in model_rows if r["label_scenario"] == "adjudicated_primary_only"]
    lines = [
        "# Adjudicated primary-only label and model sensitivity v1", "",
        "## Scope", "",
        "This is a decision-targeted primary-source sensitivity analysis, not a claim that every ChEMBL record has been paper-level verified. Explicitly adjudicated review-derived, target-mapping-error, and cellular/pathway-surrogate rows are removed. Unreviewed high-confidence rows remain unchanged. A missing arm after removal is `unresolved_missing_arm`, never inactive.", "",
        f"- Threshold: pActivity >= {THETA:.1f}.",
        f"- Bootstrap: ligand resampling, B={N_BOOT}, seed base={SEED}.",
        f"- Explicitly excluded activity rows: {len(excluded_rows)}.",
        f"- Class changes versus API-max: {len(flips_api)}.",
        f"- Unresolved ligands after exclusion: {len(unresolved)}.", "",
        "## Changed or unresolved labels", "",
        "| Pair | Ligand | API-max | Adjudicated primary-only | A max | B max |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for r in flips_api:
        lines.append(f"| {r['pair']} | {r['ligand']} | {r['api_max_class']} | {r['adjudicated_primary_only_class']} | {r['primary_only_max_A']} | {r['primary_only_max_B']} |")
    if not flips_api:
        lines.append("| — | none | — | — | — | — |")
    lines += ["", "## Model sensitivity", "", "| Pair | Model | D/A | D/B | summary_min [95% CI] | delta vs API-max | delta vs frozen | unresolved |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in primary_models:
        lines.append(
            f"| {r['pair']} | {r['model']} | {r['auroc_D_vs_A_pocketB']} | {r['auroc_D_vs_B_pocketA']} | "
            f"{r['summary_min']} [{r['summary_min_ci_lo']}, {r['summary_min_ci_hi']}] | "
            f"{r.get('delta_summary_min_vs_api_max', '')} | {r.get('delta_summary_min_vs_frozen', '')} | {r['n_unresolved']} |"
        )
    lines += [
        "", "## Interpretation guardrail", "",
        "The `adjudicated_primary_only` scenario is suitable as a transparent SI sensitivity table. It cannot yet support the sentence 'all labels were reconstructed exclusively from primary papers'. That stronger statement requires paper-level tier assignment for every record that can determine a ligand's maximum on either target.", "",
    ]
    ANA.mkdir(parents=True, exist_ok=True)
    (ANA / "PRIMARY_ONLY_LABEL_MODEL_SENSITIVITY_V1.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"excluded rows={len(excluded_rows)}; flips vs API-max={len(flips_api)}; unresolved={len(unresolved)}")
    for row in primary_models:
        print(row["pair"], row["model"], "summary_min", row["summary_min"], "delta_api", row.get("delta_summary_min_vs_api_max"))


if __name__ == "__main__":
    main()
