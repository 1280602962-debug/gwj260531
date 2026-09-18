#!/usr/bin/env python3
"""Independent checks of a freeze-rebuild directory.

Truth is the freeze score master plus deposited score/pose files, not
results/canonical and not file existence. Old QA PASS is not acceptance.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import (  # noqa: E402
    COGNATE_RMSD_SOURCE,
    D_VS_A_SCORE,
    D_VS_B_SCORE,
    FORBIDDEN_CURRENT_TOKENS,
    N_BOOT,
    N_MC_DETECTABLE,
    PRIMARY_PAIRS,
    SEED,
    TOP_FRACTION,
    is_primary_row,
    parse_finite,
)
from analysis.bootstrap_metrics import auroc, summary_min_stratified  # noqa: E402

FREEZE = ROOT / "results" / "freeze_rebuild"
REPORT = FREEZE / "VERIFICATION_REPORT.md"
FAILS: list[str] = []
NOTES: list[str] = []


def fail(msg: str) -> None:
    FAILS.append(msg)


def r3(x) -> str:
    d = Decimal(str(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return f"{d:.3f}"


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        fail(f"missing {path}")
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def load_packs(master_rows):
    packs = {p: [] for p in PRIMARY_PAIRS}
    for r in master_rows:
        if not is_primary_row(r):
            continue
        rec = dict(r)
        rec["score_A"] = float(r["score_A"])
        rec["score_B"] = float(r["score_B"])
        rec["cls"] = r["primary_class_theta6"]
        rec["score_mean"] = (rec["score_A"] + rec["score_B"]) / 2
        packs[r["pair"]].append(rec)
    return packs


def check_parse_finite() -> None:
    assert parse_finite("") is None
    assert parse_finite("nan") is None
    assert parse_finite("inf") is None
    assert parse_finite("-inf") is None
    assert parse_finite("not-a-number") is None
    assert parse_finite("1.5") == 1.5


def check_master(rows):
    pairs = {r["pair"] for r in rows}
    if pairs != set(PRIMARY_PAIRS):
        fail(f"master pairs {sorted(pairs)}")
    keys = Counter((r["pair"], r["ligand_id"], r["analysis_set"]) for r in rows)
    dups = [k for k, n in keys.items() if n > 1]
    if dups:
        fail(f"duplicate master keys {dups[:5]}")
    egfr = [r for r in rows if r["pair"] == "EGFR/HER2" and r["analysis_set"] == "main"]
    if len(egfr) != 110:
        fail(f"EGFR main n={len(egfr)}")
    ablation = {r["ligand"]: r for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")}
    for row in egfr:
        old = ablation.get(row["ligand_id"])
        if old is None:
            fail(f"EGFR {row['ligand_id']} missing ablation")
            continue
        sa = -float(old["3POZ_affinity"])
        sb = -float(old["3RCD_affinity"])
        if abs(float(row["score_A"]) - sa) > 1e-6 or abs(float(row["score_B"]) - sb) > 1e-6:
            fail(f"EGFR {row['ligand_id']} not corrected-box ablation")
    inelig = [r for r in rows if r["analysis_set"] == "main" and r.get("activity_eligible") not in ("1", "True")]
    names = {(r["pair"], r["ligand_id"]) for r in inelig}
    if ("EGFR/HER2", "EH120_059") not in names:
        fail("EH120_059 not flagged ineligible")
    if ("AChE/BChE", "AB_087") not in names:
        fail("AB_087 not flagged ineligible")
    for pair, lig in names:
        recs = [r for r in rows if r["pair"] == pair and r["ligand_id"] == lig]
        if any(is_primary_row(r) for r in recs):
            fail(f"{pair} {lig} still in primary sample")


def check_table2(packs, smin_rows, direc_rows):
    smin = {r["pair"]: r for r in smin_rows}
    direc = {(r["pair"], r["estimand"]): r for r in direc_rows}
    expected_n = {
        "EGFR/HER2": (28, 37, 32),
        "AChE/BChE": (27, 26, 28),
    }
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        dual = [r for r in recs if r["cls"] == "dual"]
        a_only = [r for r in recs if r["cls"] == "A_only"]
        b_only = [r for r in recs if r["cls"] == "B_only"]
        n = (len(dual), len(a_only), len(b_only))
        if pair in expected_n and n != expected_n[pair]:
            fail(f"{pair} n_scored {n} expected {expected_n[pair]}")
        da = auroc([r[D_VS_A_SCORE] for r in dual], [r[D_VS_A_SCORE] for r in a_only])
        db = auroc([r[D_VS_B_SCORE] for r in dual], [r[D_VS_B_SCORE] for r in b_only])
        rec_da = direc.get((pair, "AUROC_D_vs_A_pocketB"))
        rec_db = direc.get((pair, "AUROC_D_vs_B_pocketA"))
        if rec_da is None or abs(da - float(rec_da["point"])) > 1e-4:
            fail(f"{pair} D-vs-A replay {da} vs table {rec_da}")
        if rec_db is None or abs(db - float(rec_db["point"])) > 1e-4:
            fail(f"{pair} D-vs-B replay {db} vs table {rec_db}")
        sm = min(da, db)
        if abs(sm - float(smin[pair]["summary_min"])) > 1e-4:
            fail(f"{pair} summary_min replay {sm} vs {smin[pair]['summary_min']}")
        weaker = "D_vs_A_pocketB" if da <= db else "D_vs_B_pocketA"
        if smin[pair]["weaker_arm"] != weaker:
            fail(f"{pair} weaker_arm {smin[pair]['weaker_arm']} vs {weaker}")
        if smin[pair].get("bootstrap") != "class_stratified_shared_dual":
            fail(f"{pair} bootstrap {smin[pair].get('bootstrap')}")
        if str(smin[pair].get("B")) != str(N_BOOT) or str(smin[pair].get("seed")) != str(SEED):
            fail(f"{pair} B/seed {smin[pair].get('B')} {smin[pair].get('seed')}")
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        cls = np.array([r["cls"] for r in recs])
        stats = summary_min_stratified(sa, sb, cls, n_boot=N_BOOT, seed=SEED)
        if abs(stats["summary_min_ci_lo"] - float(smin[pair]["ci_lo"])) > 1e-4:
            fail(f"{pair} summary_min CI lo {stats['summary_min_ci_lo']} vs {smin[pair]['ci_lo']}")
        if abs(stats["summary_min_ci_hi"] - float(smin[pair]["ci_hi"])) > 1e-4:
            fail(f"{pair} summary_min CI hi {stats['summary_min_ci_hi']} vs {smin[pair]['ci_hi']}")


def check_ranking(packs, rank_rows, two_rows):
    rank = {r["pair"]: r for r in rank_rows}
    two = {r["pair"]: r for r in two_rows}
    for pair, recs in packs.items():
        n = len(recs)
        k = max(1, math.ceil(TOP_FRACTION * n))
        if int(rank[pair]["top_k"]) != k:
            fail(f"{pair} top_k {rank[pair]['top_k']} vs {k}")
        ranked = sorted(recs, key=lambda r: (-r["score_mean"], r["ligand_id"]))
        top = ranked[:k]
        n_dual = sum(1 for r in recs if r["cls"] == "dual")
        top_dual = sum(1 for r in top if r["cls"] == "dual")
        if int(rank[pair]["top_dual"]) != top_dual:
            fail(f"{pair} top_dual {rank[pair]['top_dual']} vs {top_dual}")
        ef = (top_dual / k) / (n_dual / n)
        if abs(ef - float(rank[pair]["ef_dual_10pct"])) > 1e-4:
            fail(f"{pair} EF {ef} vs {rank[pair]['ef_dual_10pct']}")
        if "ci_lo" not in two[pair] or two[pair]["ci_lo"] in ("", None):
            fail(f"{pair} two-pocket mean missing CI")
        dual_m = [r["score_mean"] for r in recs if r["cls"] == "dual"]
        nei_m = [r["score_mean"] for r in recs if r["cls"] == "neither"]
        pt = auroc(dual_m, nei_m)
        if abs(pt - float(two[pair]["two_pocket_mean_D_vs_neither"])) > 1e-4:
            fail(f"{pair} two-pocket mean {pt} vs {two[pair]['two_pocket_mean_D_vs_neither']}")


def check_max_median(maxmed):
    ache = [r for r in maxmed if r["pair"] == "AChE/BChE" and r["aggregation"] == "max"]
    if not ache or int(ache[0]["n_ligands"]) != 94:
        fail(f"AChE max/median n={ache}")
    if (int(ache[0]["n_dual"]), int(ache[0]["n_A_only"]), int(ache[0]["n_B_only"])) != (27, 25, 28):
        fail(f"AChE max class counts {ache[0]}")
    egfr = [r for r in maxmed if r["pair"] == "EGFR/HER2" and r["aggregation"] == "max"]
    if not egfr or int(egfr[0]["n_ligands"]) != 109:
        fail(f"EGFR max/median n={egfr}")


def check_models(inc, folds):
    vals = [abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])) for r in inc]
    if not vals:
        fail("empty ECFP incremental table")
        return 0.0
    max_abs = max(vals)
    loc = inc[vals.index(max_abs)]
    NOTES.append(f"ECFP max |Δ|={max_abs:.4f} at {loc['pair']} {loc['contrast']}")
    if any(int(r["fold_id"]) < 0 for r in folds if r["arm"] == "D_vs_A"):
        fail("negative fold_id in primary ECFP folds")
    by = Counter((r["pair"], r["arm"], r["ligand_id"]) for r in folds)
    if any(n != 1 for n in by.values()):
        fail("duplicate fold assignments")
    return max_abs


def check_manuscript(max_abs: float):
    zh = (ROOT / "docs/MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    en = (ROOT / "docs/MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    token = r3(max_abs)
    if f"最大绝对变化为 {token}" not in zh:
        fail(f"ZH Discussion missing incremental {token}")
    if f"at most {token}" not in en:
        fail(f"EN missing incremental {token}")
    if "最大绝对变化为 0.023" in zh or "at most 0.023 across" in en:
        if token != "0.023":
            fail("stale ECFP incremental 0.023 remains")
    for token_bad in FORBIDDEN_CURRENT_TOKENS:
        if token_bad in zh or token_bad in en:
            fail(f"pre-fix token {token_bad} in manuscript")


def check_robustness(robust, rec_sub, cluster, rmsd, det):
    engines = {(r["pair"], r["engine"]) for r in robust}
    for key in (
        ("EGFR/HER2", "gnina_dock_mode1"),
        ("PIK3CA/mTOR", "gnina_dock_mode1"),
        ("JAK1/TYK2", "gnina_dock_mode1"),
        ("PIK3CA/mTOR", "vina_alt_4JPS"),
        ("PIK3CA/mTOR", "vina_alt_5DXT"),
        ("PIK3CA/mTOR", "vina_alt_4JSX"),
    ):
        if key not in engines:
            fail(f"missing robustness row {key}")
    if len(rec_sub) != 3:
        fail(f"receptor_substitution n={len(rec_sub)}")
    jak_doc = [r for r in cluster if r["pair"] == "JAK1/TYK2" and r["estimator"] == "document_cluster"]
    if not jak_doc or jak_doc[0]["status"] != "unresolved_mapping_unavailable":
        fail(f"JAK document cluster {jak_doc}")
    if any(r.get("pdb") == "3POZ" and abs(float(r.get("calcrrms_best_A") or r.get("best_A") or 0) - 0.760) < 1e-6 for r in rmsd):
        fail("PR #35 reconstructed EGFR 0.760 Å used as current RMSD")
    egfr_rmsd = next((r for r in rmsd if r.get("pdb") == "3POZ" or r.get("protein") == "EGFR"), None)
    if egfr_rmsd is None:
        fail("EGFR RMSD row missing")
    else:
        best = float(egfr_rmsd.get("calcrrms_best_A") or egfr_rmsd.get("best_A"))
        if abs(best - 1.019) > 1e-3:
            fail(f"EGFR current RMSD best={best}")
    if {r["pair"] for r in det} != set(PRIMARY_PAIRS):
        fail("detectable-effect pair set")
    if any(str(r.get("n_mc")) != str(N_MC_DETECTABLE) for r in det):
        fail("detectable-effect n_mc")
    ache = next(r for r in det if r["pair"] == "AChE/BChE" and r["contrast"] == "summary_min" and r["true_auroc"] == "0.50")
    if ache["n_neg"] != "26/28":
        fail(f"detectable-effect AChE n_neg={ache['n_neg']}")
    source_n = len(read_csv(COGNATE_RMSD_SOURCE))
    if source_n != 14 or len(rmsd) != 14:
        fail(f"RMSD n source={source_n} freeze={len(rmsd)}")


def write_report():
    lines = ["# Freeze verification", ""]
    lines.append(f"freeze_dir: `{FREEZE}`")
    lines.append(f"result: {'FAIL' if FAILS else 'PASS'}")
    lines.append("")
    if FAILS:
        lines.append("## Failures")
        lines.extend(f"- {m}" for m in FAILS)
        lines.append("")
    lines.append("## Notes")
    lines.extend(f"- {m}" for m in NOTES)
    lines.append("")
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


def main() -> int:
    if not FREEZE.is_dir():
        print("FAIL: freeze directory missing")
        return 1
    check_parse_finite()
    master = read_csv(FREEZE / "current_score_master.csv")
    check_master(master)
    packs = load_packs(master)
    check_table2(packs, read_csv(FREEZE / "primary_summary_min.csv"), read_csv(FREEZE / "primary_directional_auroc.csv"))
    check_ranking(packs, read_csv(FREEZE / "top10_operating_points.csv"), read_csv(FREEZE / "two_pocket_mean_ranking.csv"))
    check_max_median(read_csv(FREEZE / "max_vs_median_sensitivity.csv"))
    max_abs = check_models(read_csv(FREEZE / "ecfp4_incremental_information.csv"), read_csv(FREEZE / "model_fold_assignments.csv"))
    check_robustness(
        read_csv(FREEZE / "computational_robustness.csv"),
        read_csv(FREEZE / "receptor_substitution.csv"),
        read_csv(FREEZE / "cluster_bootstrap_sensitivity.csv"),
        read_csv(FREEZE / "cognate_rmsd.csv"),
        read_csv(FREEZE / "detectable_effect_simulation.csv"),
    )
    check_manuscript(max_abs)
    write_report()
    return 1 if FAILS else 0


if __name__ == "__main__":
    raise SystemExit(main())
