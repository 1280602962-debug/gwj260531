#!/usr/bin/env python3
"""Read current files and check the eight-pair scientific chain.

Numeric replay only. Manuscript prose is not checked.
Does not use SHA/hash/checksum/git-head as a PASS/FAIL gate.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import (  # noqa: E402
    CANONICAL_CSV_NAMES,
    EGFR_PRODUCTION_SEED,
    EGFR_UNIFORM_VINA_CSV,
    FIVE_SEED_EGFR_CORRECTED_BOX_SCORES,
    current_egfr_score_source,
    egfr_uniform_ready,
    is_primary_row,
    parse_finite,
)
from analysis.bootstrap_metrics import auroc  # noqa: E402

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
TRACK_B = {"JAK1/JAK2", "JAK1/TYK2", "F2/F10", "PPARG/PPARA", "PPARA/PPARD"}
FALLBACK = ("AB_001", "AB_053", "AB_054", "AB_056", "AB_097")


def r3(x) -> str:
    d = Decimal(str(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return f"{d:.3f}"


def fail(msg: str) -> None:
    print("FAIL:", msg)
    raise SystemExit(1)


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        fail(f"missing {path.relative_to(ROOT)}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def check_pack_tables() -> None:
    canon = ROOT / "results/canonical"
    pack = ROOT / "submission_pack/tables/canonical"
    names = {p.name for p in canon.glob("*.csv")}
    if names != set(CANONICAL_CSV_NAMES):
        fail(f"canonical CSV set {sorted(names)} != {list(CANONICAL_CSV_NAMES)}")
    for name in sorted(names):
        left = canon / name
        right = pack / name
        if not right.is_file():
            fail(f"submission pack missing {name}")
        if left.read_bytes() != right.read_bytes():
            fail(f"submission pack {name} differs from results/canonical")
    print(f"PASS: table pack {len(names)}/{len(names)} canonical CSV identical")


def check_leftovers() -> None:
    if (ROOT / "results/freeze_rebuild").exists():
        fail("results/freeze_rebuild/ must not exist")
    if (ROOT / "figures/jcim_article/plotted_values.json").exists():
        fail("figures/jcim_article/plotted_values.json must not exist")
    leftover_0112 = list((ROOT / "results").rglob("*0.0112*")) if (ROOT / "results").is_dir() else []
    if leftover_0112:
        fail(f"leftover 0.0112 current-source candidate {leftover_0112}")


def check_writing_index_paths() -> None:
    """Existence only. Not a SHA/hash gate."""
    required = [
        ROOT / "figures/jcim_article/plotted_values_postfix.json",
        ROOT / "docs/PR39_SCIENTIFIC_DATA_AUDIT.md",
        ROOT / "docs/WRITING_INDEX_FREEZE.md",
        ROOT / "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
        ROOT / "data/jcim_chembl_universe_v0/tables/track_b_local_run_v1.yaml",
        ROOT / "results/canonical/ecfp4_incremental_information.csv",
        ROOT / "results/canonical/computational_robustness.csv",
        ROOT / "results/canonical/primary_summary_min.csv",
        ROOT / "results/canonical/primary_directional_auroc.csv",
        ROOT / "results/canonical/five_seed_summary_min.csv",
        ROOT / "results/canonical/five_seed_fixed_membership_sensitivity.csv",
        ROOT / "results/canonical/current_score_master.csv",
        ROOT / "data/provenance/receptor_input_registry.csv",
    ]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.is_file()]
    if missing:
        fail(f"writing-index/lock path missing: {missing}")
    index = (ROOT / "docs/WRITING_INDEX_FREEZE.md").read_text(encoding="utf-8")
    lock = (ROOT / "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md").read_text(encoding="utf-8")
    if "HISTORICAL_LEFTOVER_FILES.md" in index or "HISTORICAL_LEFTOVER_FILES.md" in lock:
        fail("HISTORICAL_LEFTOVER_FILES.md is still cited")
    if "plotted_values.json" in index and "plotted_values_postfix.json" not in index:
        fail("WRITING_INDEX_FREEZE still treats plotted_values.json as authority")


def provenance_summary(master: list[dict]) -> list[dict]:
    rows = []
    for pair in PAIRS:
        recs = [r for r in master if r["pair"] == pair and r.get("analysis_set") == "main"]
        kinds = sorted({r.get("activity_source_kind") or "" for r in recs})
        for kind in kinds:
            sub = [r for r in recs if (r.get("activity_source_kind") or "") == kind]
            n_complete = sum(1 for r in sub if r.get("complete_case") in ("1", "True"))
            n_primary = sum(1 for r in sub if is_primary_row(r))
            rows.append(
                {
                    "pair": pair,
                    "activity_source_kind": kind,
                    "n_main": len(sub),
                    "n_complete_case": n_complete,
                    "n_primary_used": n_primary,
                }
            )
            print(
                f"PROVENANCE {pair:14} {kind:28} n_main={len(sub):3} "
                f"n_complete={n_complete:3} n_primary={n_primary:3}"
            )
    return rows


def check_five_seed_comparable() -> None:
    seeds = read_csv(ROOT / "results/canonical/five_seed_summary_min.csv")
    smin = {r["pair"]: r for r in read_csv(ROOT / "results/canonical/primary_summary_min.csv")}
    fixed = read_csv(ROOT / "results/canonical/five_seed_fixed_membership_sensitivity.csv")
    comparable_vals = []
    egfr_primary = None
    for r in seeds:
        comparable = str(r.get("comparable_to_current_primary", "")).strip()
        if comparable != "1":
            fail(f"five-seed {r['pair']} seed {r['seed']} not comparable ({comparable})")
        if str(r.get("same_protocol_as_primary", "")).strip() != "1":
            fail(f"five-seed {r['pair']} seed {r['seed']} same_protocol_as_primary={r.get('same_protocol_as_primary')}")
        membership = str(r.get("same_membership_as_primary", "")).strip()
        if r["pair"] == "AChE/BChE":
            if membership != "0":
                fail(f"AChE five-seed must have same_membership_as_primary=0, got {membership}")
        elif membership != "1":
            fail(f"{r['pair']} seed {r['seed']} same_membership_as_primary={membership}")
        if r.get("realization_status") != "current_or_compatible":
            fail(f"five-seed {r['pair']} seed {r['seed']} status={r.get('realization_status')}")
        comparable_vals.append(float(r["summary_min"]))
        if r["pair"] == "EGFR/HER2" and str(r["seed"]) == "20260727":
            egfr_primary = r
            if abs(float(r["summary_min"]) - float(smin["EGFR/HER2"]["summary_min"])) > 1e-4:
                fail(f"EGFR five-seed production summary_min {r['summary_min']} vs Table 2 {smin['EGFR/HER2']['summary_min']}")
            if (int(r["n_dual"]), int(r["n_A_only"]), int(r["n_B_only"])) != (
                int(smin["EGFR/HER2"]["n_dual"]),
                int(smin["EGFR/HER2"]["n_A_only"]),
                int(smin["EGFR/HER2"]["n_B_only"]),
            ):
                fail(f"EGFR five-seed production class counts {r['n_dual']}/{r['n_A_only']}/{r['n_B_only']}")
    if egfr_primary is None:
        fail("EGFR five-seed production seed 20260727 missing")
    ache_fixed = [r for r in fixed if r["pair"] == "AChE/BChE"]
    if not ache_fixed:
        fail("AChE fixed-membership sensitivity missing")
    n_inter = {int(r["n_intersection"]) for r in ache_fixed}
    if n_inter != {88}:
        fail(f"AChE five-seed intersection expected n=88, got {n_inter}")
    if any(r.get("qualitative_change") != "no" for r in ache_fixed):
        fail("AChE fixed-membership qualitative_change is not no")
    plotted = json.loads((ROOT / "figures/jcim_article/plotted_values_postfix.json").read_text(encoding="utf-8"))
    nested = plotted.get("nested_plotted") or {}
    fig4c = nested.get("fig4C") or {}
    egfr = fig4c.get("EGFR/HER2") or {}
    if int(egfr.get("comparable_to_current_primary", 0)) != 1:
        fail(f"Figure 5C EGFR not marked comparable {egfr}")
    if egfr.get("primary") in (None, "", "NA"):
        fail("Figure 5C missing EGFR five-seed primary overlay")
    print(
        f"PASS: EGFR five-seed corrected-box comparable; AChE membership documented n=88; range "
        f"{min(comparable_vals):.4f}–{max(comparable_vals):.4f}"
    )


def check_ab056(master: list[dict]) -> None:
    recs = [r for r in master if is_primary_row(r) and r["pair"] == "AChE/BChE"]
    dual = [r for r in recs if r["primary_class_theta6"] == "dual"]
    a_only = [r for r in recs if r["primary_class_theta6"] == "A_only"]
    b_only = [r for r in recs if r["primary_class_theta6"] == "B_only"]
    da = auroc([float(r["score_B"]) for r in dual], [float(r["score_B"]) for r in a_only])
    db = auroc([float(r["score_A"]) for r in dual], [float(r["score_A"]) for r in b_only])
    sm = min(da, db)
    a_ex = [r for r in a_only if r["ligand_id"] != "AB_056"]
    da_ex = auroc([float(r["score_B"]) for r in dual], [float(r["score_B"]) for r in a_ex])
    db_ex = auroc([float(r["score_A"]) for r in dual], [float(r["score_A"]) for r in b_only])
    sm_ex = min(da_ex, db_ex)
    print(
        f"QA AChE AB_056 sensitivity current D_vs_A={da:.4f} D_vs_B={db:.4f} smin={sm:.4f}; "
        f"exclude AB_056 D_vs_A={da_ex:.4f} D_vs_B={db_ex:.4f} smin={sm_ex:.4f}"
    )
    if abs(sm_ex - sm) > 1e-4:
        fail(f"excluding AB_056 changed summary_min {sm} -> {sm_ex}")
    if abs(da_ex - 0.6504) > 5e-4 or abs(db_ex - 0.6058) > 5e-4:
        fail(f"exclude-AB_056 directional {da_ex:.4f}/{db_ex:.4f}")


def check_folds() -> None:
    folds = read_csv(ROOT / "results/canonical/model_fold_assignments.csv")
    oof = read_csv(ROOT / "results/canonical/ecfp4_oof_predictions.csv")
    print(f"QA ECFP fold map row count={len(folds)}")
    scaf = defaultdict(set)
    for r in folds:
        scaf[(r["pair"], r["arm"], r["scaffold"])].add(r["fold_id"])
    split = {k: v for k, v in scaf.items() if len(v) != 1}
    if split:
        fail(f"scaffold split across folds {list(split)[:3]}")
    by_fold = {(r["pair"], r["arm"], r["ligand_id"]): r["fold_id"] for r in folds}
    models_by = defaultdict(set)
    for r in oof:
        if str(r.get("scaled", "0")) not in ("0", "False"):
            continue
        models_by[(r["pair"], r["arm"], r["ligand_id"])].add(r["model"])
        if (r["pair"], r["arm"], r["ligand_id"]) not in by_fold:
            fail(f"OOF ligand missing fold {r['pair']} {r['arm']} {r['ligand_id']}")
    needed = {"ECFP4", "docking", "ECFP4+docking"}
    for key, models in models_by.items():
        if models != needed:
            fail(f"{key} models {models} do not share the contrast")
    print("PASS: scaffolds stay in one fold; three models share fold_id per contrast")


def main() -> int:
    master = read_csv(ROOT / "results/canonical/current_score_master.csv")
    pairs = {r["pair"] for r in master}
    if pairs != set(PAIRS):
        fail(f"master pairs {sorted(pairs)}")
    counts = Counter((r["pair"], r["ligand_id"], r.get("analysis_set")) for r in master)
    dups = [k for k, n in counts.items() if n > 1]
    if dups:
        fail(f"duplicate pair/ligand keys {dups[:5]}")

    panel = {
        r["panel_id"]: r
        for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")
    }
    egfr = [r for r in master if r["pair"] == "EGFR/HER2" and r.get("analysis_set") == "main"]
    if {r["ligand_id"] for r in egfr} != set(panel):
        fail(f"EGFR master ligands != panel_v0_120 ({len(egfr)} vs {len(panel)})")
    if egfr_uniform_ready():
        want = {}
        for r in read_csv(EGFR_UNIFORM_VINA_CSV):
            if str(r.get("seed")) != str(EGFR_PRODUCTION_SEED):
                continue
            if r.get("status") not in {"ok", "success"}:
                continue
            energy = parse_finite(r.get("vina_mode1"))
            if energy is None:
                continue
            rec = want.setdefault(r["ligand"], {})
            if r.get("pdb") == "3POZ":
                rec["A"] = -energy
            elif r.get("pdb") == "3RCD":
                rec["B"] = -energy
        for row in egfr:
            src = want.get(row["ligand_id"], {})
            if "A" in src and "B" in src:
                if abs(float(row["score_A"]) - src["A"]) > 1e-6 or abs(float(row["score_B"]) - src["B"]) > 1e-6:
                    fail(f"EGFR {row['ligand_id']} score mismatch vs uniform Vina")
            elif row.get("complete_case") in ("1", "True", 1):
                fail(f"EGFR {row['ligand_id']} complete_case without both uniform scores")
            if row.get("score_source") != current_egfr_score_source():
                fail(f"EGFR {row['ligand_id']} score_source={row.get('score_source')}")
    else:
        ab = {
            r["ligand"]: r
            for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")
        }
        if len(egfr) != len(ab):
            fail(f"EGFR n={len(egfr)} ablation={len(ab)}")
        for row in egfr:
            old = ab.get(row["ligand_id"])
            if old is None:
                fail(f"EGFR {row['ligand_id']} missing from ablation scores")
            sa = -float(old["3POZ_affinity"])
            sb = -float(old["3RCD_affinity"])
            if abs(float(row["score_A"]) - sa) > 1e-6 or abs(float(row["score_B"]) - sb) > 1e-6:
                fail(f"EGFR {row['ligand_id']} score mismatch")

    for pair, lig in (("EGFR/HER2", "EH120_059"), ("AChE/BChE", "AB_087")):
        recs = [r for r in master if r["pair"] == pair and r["ligand_id"] == lig]
        if not recs or recs[0].get("activity_status") != "unresolved_missing_arm":
            fail(f"{pair} {lig} not unresolved_missing_arm")
        if any(is_primary_row(r) for r in recs):
            fail(f"{pair} {lig} entered primary")

    primary_fb = []
    for lig in FALLBACK:
        recs = [r for r in master if r["pair"] == "AChE/BChE" and r["ligand_id"] == lig and r.get("analysis_set") == "main"]
        if not recs or recs[0].get("activity_source_kind") != "panel_pchembl_no_audit_rows":
            fail(f"AChE {lig} provenance {recs}")
        if any(is_primary_row(r) for r in recs):
            primary_fb.append(lig)
    if primary_fb != ["AB_056"]:
        fail(f"AChE fallback in primary {primary_fb}")
    for r in master:
        if r.get("analysis_set") == "main" and r["pair"] in TRACK_B:
            if r.get("activity_source_kind") != "chembl37_dump_panel":
                fail(f"{r['pair']} {r['ligand_id']} Track-B kind={r.get('activity_source_kind')}")

    provenance_summary(master)
    check_ab056(master)

    direc = {(r["pair"], r["estimand"]): r for r in read_csv(ROOT / "results/canonical/primary_directional_auroc.csv")}
    smin = {r["pair"]: r for r in read_csv(ROOT / "results/canonical/primary_summary_min.csv")}
    packs: dict[str, dict[str, list]] = {}
    for row in master:
        if not is_primary_row(row):
            continue
        cls = row.get("primary_class_theta6") or ""
        if cls not in {"dual", "A_only", "B_only"}:
            continue
        packs.setdefault(row["pair"], {}).setdefault(cls, []).append(row)
    for pair in PAIRS:
        dual = packs[pair]["dual"]
        a_only = packs[pair]["A_only"]
        b_only = packs[pair]["B_only"]
        da = auroc([float(r["score_B"]) for r in dual], [float(r["score_B"]) for r in a_only])
        db = auroc([float(r["score_A"]) for r in dual], [float(r["score_A"]) for r in b_only])
        rec_da = direc[(pair, "AUROC_D_vs_A_pocketB")]
        rec_db = direc[(pair, "AUROC_D_vs_B_pocketA")]
        if abs(da - float(rec_da["point"])) > 1e-4 or abs(db - float(rec_db["point"])) > 1e-4:
            fail(f"{pair} directional AUROC replay mismatch")
        sm = min(da, db)
        if abs(sm - float(smin[pair]["summary_min"])) > 1e-4:
            fail(f"{pair} summary_min replay mismatch")
        table_n = (int(smin[pair]["n_dual"]), int(smin[pair]["n_A_only"]), int(smin[pair]["n_B_only"]))
        master_n = (len(dual), len(a_only), len(b_only))
        if table_n != master_n:
            fail(f"{pair} class counts table {table_n} vs master {master_n}")

    inc = read_csv(ROOT / "results/canonical/ecfp4_incremental_information.csv")
    max_abs = max(abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])) for r in inc)
    loc = max(inc, key=lambda r: abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])))
    print(f"QA ECFP max |Δ AUROC|={max_abs:.4f} at {loc['pair']} {loc['contrast']}")
    plotted = json.loads((ROOT / "figures/jcim_article/plotted_values_postfix.json").read_text(encoding="utf-8"))
    recs = plotted.get("records") or []
    plotted_max = None
    for rec in recs:
        if rec.get("metric") == "fig3B_max_abs":
            plotted_max = abs(float(rec["raw_value"]))
            break
    if plotted_max is None:
        nested = plotted.get("nested_plotted") or {}
        plotted_max = abs(float(nested["fig3B_max_abs"]))
    if abs(plotted_max - max_abs) > 1e-6:
        fail(f"Figure 3B plotted max |Δ|={plotted_max} != canonical {max_abs}")

    check_folds()
    check_five_seed_comparable()
    check_pack_tables()
    check_leftovers()
    check_writing_index_paths()

    pack_en = (ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    if r3(smin["EGFR/HER2"]["summary_min"]) not in pack_en:
        fail("submission_pack EN manuscript missing EGFR summary_min")

    det = read_csv(ROOT / "results/canonical/detectable_effect_simulation.csv")
    det_pairs = {r["pair"] for r in det}
    if det_pairs != set(PAIRS):
        fail(f"detectable-effect pairs {sorted(det_pairs)}")
    if any(r.get("n_mc") not in {"1000", 1000} for r in det):
        fail("detectable-effect n_mc is not 1000")
    n_a = len(packs["AChE/BChE"]["A_only"])
    n_b = len(packs["AChE/BChE"]["B_only"])
    ache = next(
        r
        for r in det
        if r["pair"] == "AChE/BChE" and r["contrast"] == "summary_min" and r["true_auroc"] == "0.50"
    )
    if ache["n_neg"] != f"{n_a}/{n_b}":
        fail(f"AChE detectable-effect n_A/n_B={ache['n_neg']} vs master {n_a}/{n_b}")
    if any(r.get("bootstrap") != "class_stratified_shared_dual" for r in det):
        fail("detectable-effect bootstrap is not class_stratified_shared_dual")

    scores = EGFR_UNIFORM_VINA_CSV if egfr_uniform_ready() else FIVE_SEED_EGFR_CORRECTED_BOX_SCORES
    if not scores.is_file():
        fail(f"EGFR five-seed score file missing: {scores}")

    print("PASS: current score master, Table 2 replay, provenance, five-seed, plotted ECFP, pack tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
