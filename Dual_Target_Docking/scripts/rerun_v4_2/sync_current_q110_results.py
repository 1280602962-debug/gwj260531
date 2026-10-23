"""Promote accepted Q110 rows into current eight-pair result tables.

No docking, model fitting, AUROC or bootstrap is run here. Baseline rows are
read from an immutable Git commit, so repeated runs cannot accumulate rows.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parent
PAIR = "PIK3CA/mTOR"
BASELINE = "dc028b3350e804cbf8f40acea08261ed4f89c4fa"
QREL = Path("extension_jcim_validation/pik3ca_mtor_q110_v42")
Q = ROOT / QREL
RERUN = Path("reruns/UNIFORM_RERUN_V4_2_20260921")
STAGE = Path("results/jcim_stage2")
CURRENT = Path("results/current")
MODULES = {
    "module_B_negative_class.csv": "module_B.csv",
    "module_C_wrong_pocket.csv": "module_C.csv",
    "module_D1_composition.csv": "module_D1.csv",
    "module_D_chemistry_oof_ligands.csv": "module_D_oof.csv",
    "module_D_chemistry_units.csv": "module_D_units.csv",
    "module_E_five_seed.csv": "module_E_seed.csv",
    "module_E_five_seed_summary.csv": "module_E_summary.csv",
    "module_H_missingness.csv": "module_H.csv",
    "module_I_descriptive.csv": "module_I_desc.csv",
    "module_I_ligand_ranks.csv": "module_I_lig.csv",
    "module_I_method_compare.csv": "module_I_cmp.csv",
}


def parse(text):
    reader = csv.DictReader(io.StringIO(text))
    return list(reader.fieldnames), list(reader)


def read(rel):
    return parse((ROOT / rel).read_text(encoding="utf-8"))


def write(rel, fields, rows):
    path = ROOT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def old(rel, ref):
    text = subprocess.check_output(
        ["git", "show", f"{ref}:Dual_Target_Docking/{Path(rel).as_posix()}"],
        cwd=REPO, text=True, encoding="utf-8",
    )
    return parse(text)


def replace_pair(base, replacement):
    result = []
    inserted = False
    for row in base:
        if row["pair"] == PAIR:
            if not inserted:
                result.extend(replacement)
                inserted = True
        else:
            result.append(row)
    assert inserted
    assert [r for r in result if r["pair"] != PAIR] == [r for r in base if r["pair"] != PAIR]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", default=BASELINE)
    args = parser.parse_args()
    ref = args.baseline_ref
    checks = []

    for filename, qfile, keys in (
        ("PRIMARY_DIRECTIONAL_METRICS.csv", "Q110_PRIMARY_METRICS.csv", ("method", "contrast")),
        ("PRIMARY_METHOD_DELTA_METRICS.csv", "Q110_PRIMARY_DELTAS.csv", ("delta",)),
    ):
        rel = Path("results/formal_metrics") / filename
        fields, base = old(rel, ref)
        templates = {tuple(r[k] for k in keys): r for r in base if r["pair"] == PAIR}
        qrows = read(QREL / qfile)[1]
        replacement = []
        for qr in qrows:
            row = dict(templates[tuple(qr[k] for k in keys)])
            row.update({k: v for k, v in qr.items() if k in fields})
            if "contrast" in row:
                row.update(n_A_only="32", n_B_only="32")
            replacement.append(row)
        write(rel, fields, replace_pair(base, replacement))
        checks.append({"table": str(rel), "rows": len(base), "other_seven_pairs_unchanged": True})

    for name, qname in MODULES.items():
        rel = STAGE / name
        fields, base = old(rel, ref)
        qfields, qrows = read(QREL / "q110_stage2" / qname)
        assert set(fields) == set(qfields), name
        merged = replace_pair(base, qrows)
        write(rel, fields, merged)
        checks.append({"table": str(rel), "rows": len(merged), "other_seven_pairs_unchanged": True})

    rel = STAGE / "module_F_alt_receptor.csv"
    fields, base = old(rel, ref)
    alt = [r for r in base if r["pair"] != PAIR]
    assert len(alt) == 6
    write(rel, fields, alt)

    oldpop = RERUN / "13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
    fields, base = old(oldpop, ref)
    qrows = read(QREL / "Q110_POPULATION.csv")[1]
    pop = replace_pair(base, qrows)
    write(CURRENT / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv", fields, pop)
    # This old derived population is superseded; raw docking files stay intact.
    (ROOT / oldpop).unlink(missing_ok=True)

    fields, base = old(RERUN / "00_protocol/pair_ligand_mapping.csv", ref)
    membership = read(QREL / "PIK3CA_MTOR_Q110_MEMBERSHIP.csv")[1]
    # Identity and SMILES come from the frozen membership. Fields not provided
    # there stay blank rather than being guessed from unrelated old IDs.
    mapping_q = []
    for m, p in zip(membership, qrows):
        assert m["set_id"] == p["canonical_ligand_id"]
        row = {k: "" for k in fields}
        row.update(run_id="Q110_V4_2", pair=PAIR, panel_id=p["panel_id"],
                   global_ligand_entity_id=p["global_ligand_entity_id"],
                   canonical_ligand_id=p["canonical_ligand_id"], parent_independent="1",
                   parent_role="canonical", analysis_set=p["analysis_set"],
                   **{"class": p["class"]}, activity_eligible="1", prepare_3d="1",
                   canonical_smiles=m["canonical_smiles"])
        mapping_q.append(row)
    mapping = replace_pair(base, mapping_q)
    write(CURRENT / "pair_ligand_mapping.csv", fields, mapping)

    foldrel = Path("results/canonical/model_fold_assignments.csv")
    fields, base = old(foldrel, ref)
    fold_q = read(QREL / "Q110_FOLD_ASSIGNMENTS.csv")[1]
    write(foldrel, fields, replace_pair(base, fold_q))

    layerrel = Path("analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv")
    fields, base = old(layerrel, ref)
    layerq = []
    for arm in ("D_vs_A", "D_vs_B"):
        row = {k: "" for k in fields}
        row.update(pair=PAIR, arm=arm, n_labeled_population="64",
                   n_m0_directional_population="64", n_m0_missing="0",
                   n_fold_unassigned_in_labeled_population="0",
                   n_fold_unassigned_after_m0_filter="0", n_final_chemistry_oof="64")
        layerq.append(row)
    layers = replace_pair(base, layerq)
    write(layerrel, fields, layers)

    cmp = read(STAGE / "module_I_method_compare.csv")[1]
    states = ("CONCORDANT_IMPROVE", "CONCORDANT_WORSEN", "DISCORDANT", "NO_CHANGE")
    counts = []
    for method in ("M1", "M1b", "M2", "M3"):
        count = Counter(r["concordance"] for r in cmp if r["method"] == method)
        counts.append({"method": method, "n_pairs": "8", "supporting_only": "YES" if method == "M1b" else "NO", "bootstrap_ci_added": "NO", **{s: str(count[s]) for s in states}})
    write(Path("results/jcim_stage2_figures/module_I_concordance_counts.csv"),
          ["method", "n_pairs", "supporting_only", "bootstrap_ci_added", *states], counts)

    eligible = [r for r in pop if r["activity_eligible"] == "1" and r["class"]]
    totals = Counter(r["class"] for r in eligible)
    assert len(pop) == 869 and len(eligible) == 867 and len(mapping) == 870
    assert totals == {"dual": 248, "A_only": 256, "B_only": 252, "neither": 111}
    for m in ("M0", "M1", "M1b", "M2", "M3"):
        for side in ("A", "B"):
            assert all(r[f"{m}_{side}_valid"] == "1" for r in qrows)
    report = {
        "state": "current_results_synchronized", "source_commit": BASELINE,
        "population_rows": len(pop), "eligible_observations": len(eligible),
        "classes": dict(totals), "mapping_rows": len(mapping),
        "layer_counts": {k: sum(int(r[k]) for r in layers) for k in
                         ("n_labeled_population", "n_m0_directional_population", "n_final_chemistry_oof")},
        "alt_receptor_pairs": sorted({r["pair"] for r in alt}),
        "alt_receptor_conditions": len(alt), "table_checks": checks,
        "scientific_recomputation": False,
    }
    (ROOT / CURRENT / "RESULT_SYNC_REPORT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    manifest = {
        "state": "current_primary", "source_commit": BASELINE,
        "pair_override": {PAIR: str(QREL)}, "other_seven_pairs": "frozen V4.2",
        "population": str(CURRENT / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"),
        "mapping": str(CURRENT / "pair_ligand_mapping.csv"),
        "folds": str(foldrel), "primary_metrics": "results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv",
        "primary_deltas": "results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv",
        "stage2": str(STAGE), "historical_pm48": f"git:{BASELINE}",
        "PIK3CA_mTOR_alt_receptor": "not_computed_for_Q110; excluded_from_current_results",
    }
    (ROOT / "CURRENT_RESULT_AUTHORITIES.json").write_text(json.dumps(manifest, indent=2)+"\n")
    (ROOT / STAGE / "RUN_MANIFEST.json").write_text(json.dumps({
        "state": "synchronized_from_accepted_results", "source_commit": BASELINE,
        "pair_override": {PAIR: str(QREL / "q110_stage2")},
        "other_seven_pairs": "unchanged", "scientific_recomputation": False,
        "PIK3CA_mTOR_module_F": "not_computed_for_Q110; excluded",
        "blocked_modules": ["G", "J", "K"],
    }, indent=2)+"\n")
    print(json.dumps({k: v for k, v in report.items() if k != "table_checks"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
