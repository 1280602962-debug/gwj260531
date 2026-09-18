#!/usr/bin/env python3
"""Independent checks of a freeze-rebuild directory.

Truth is the freeze score master plus deposited score/pose files, not
results/canonical and not file existence. Old QA PASS is not acceptance.
Manuscript prose is not checked.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.analysis_config import (  # noqa: E402
    CANONICAL_CSV_NAMES,
    COGNATE_RMSD_SOURCE,
    D_VS_A_SCORE,
    D_VS_B_SCORE,
    N_BOOT,
    N_MC_DETECTABLE,
    PRIMARY_PAIRS,
    SEED,
    TOP_FRACTION,
    is_primary_row,
    parse_finite,
)
from analysis.bootstrap_metrics import auroc, summary_min_stratified  # noqa: E402

FAILS: list[str] = []
NOTES: list[str] = []
FREEZE = ROOT / "results" / "freeze_rebuild"


def fail(msg: str) -> None:
    FAILS.append(msg)


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
        rec["score_worst"] = min(rec["score_A"], rec["score_B"])
        packs[r["pair"]].append(rec)
    return packs


def check_parse_finite() -> None:
    assert parse_finite("") is None
    assert parse_finite("nan") is None
    assert parse_finite("inf") is None
    assert parse_finite("-inf") is None
    assert parse_finite("not-a-number") is None
    assert parse_finite("1.5") == 1.5


def check_expected_files(freeze: Path) -> None:
    have = {p.name for p in freeze.glob("*.csv")}
    missing = sorted(set(CANONICAL_CSV_NAMES) - have)
    extra = sorted(have - set(CANONICAL_CSV_NAMES))
    if missing:
        fail(f"rebuild missing CSV {missing}")
    if extra:
        fail(f"rebuild unexpected CSV {extra}")


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
        sa = np.array([r["score_A"] for r in recs])
        sb = np.array([r["score_B"] for r in recs])
        cls = np.array([r["cls"] for r in recs])
        stats = summary_min_stratified(sa, sb, cls, n_boot=N_BOOT, seed=SEED)
        if abs(stats["auroc_D_vs_A_ci_lo"] - float(rec_da["ci_lo"])) > 1e-4 or abs(stats["auroc_D_vs_A_ci_hi"] - float(rec_da["ci_hi"])) > 1e-4:
            fail(f"{pair} D-vs-A CI replay {stats['auroc_D_vs_A_ci_lo']},{stats['auroc_D_vs_A_ci_hi']} vs {rec_da['ci_lo']},{rec_da['ci_hi']}")
        if abs(stats["auroc_D_vs_B_ci_lo"] - float(rec_db["ci_lo"])) > 1e-4 or abs(stats["auroc_D_vs_B_ci_hi"] - float(rec_db["ci_hi"])) > 1e-4:
            fail(f"{pair} D-vs-B CI replay {stats['auroc_D_vs_B_ci_lo']},{stats['auroc_D_vs_B_ci_hi']} vs {rec_db['ci_lo']},{rec_db['ci_hi']}")
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
        if abs(stats["summary_min_ci_lo"] - float(smin[pair]["ci_lo"])) > 1e-4:
            fail(f"{pair} summary_min CI lo {stats['summary_min_ci_lo']} vs {smin[pair]['ci_lo']}")
        if abs(stats["summary_min_ci_hi"] - float(smin[pair]["ci_hi"])) > 1e-4:
            fail(f"{pair} summary_min CI hi {stats['summary_min_ci_hi']} vs {smin[pair]['ci_hi']}")


def check_fixed_score(packs, fixed_rows):
    by = {(r["pair"], r["contrast"]): r for r in fixed_rows}
    egfr = by[("EGFR/HER2", "D_vs_B_or_neither_pocketA")]
    if abs(float(egfr["auroc_dual_vs_selective"]) - 0.3237) > 1e-3:
        fail(f"EGFR pocket A D-vs-B {egfr['auroc_dual_vs_selective']}")
    if abs(float(egfr["auroc_dual_vs_neither"]) - 0.7857) > 1e-3:
        fail(f"EGFR pocket A D-vs-neither {egfr['auroc_dual_vs_neither']}")
    if abs(float(egfr["delta_neither_minus_selective"]) - 0.4621) > 1e-3:
        fail(f"EGFR pocket A delta {egfr['delta_neither_minus_selective']}")
    jak = by[("JAK1/TYK2", "D_vs_B_or_neither_pocketA")]
    if abs(float(jak["delta_neither_minus_selective"]) - 0.4438) > 1e-3:
        fail(f"JAK1/TYK2 pocket A delta {jak['delta_neither_minus_selective']}")
    for pair, recs in packs.items():
        dual = [r for r in recs if r["cls"] == "dual"]
        neither = [r for r in recs if r["cls"] == "neither"]
        for contrast, key, sel_cls in (
            ("D_vs_B_or_neither_pocketA", "score_A", "B_only"),
            ("D_vs_A_or_neither_pocketB", "score_B", "A_only"),
        ):
            selective = [r for r in recs if r["cls"] == sel_cls]
            auc_s = auroc([r[key] for r in dual], [r[key] for r in selective])
            auc_n = auroc([r[key] for r in dual], [r[key] for r in neither])
            rec = by[(pair, contrast)]
            if abs(auc_s - float(rec["auroc_dual_vs_selective"])) > 1e-4:
                fail(f"{pair} {contrast} selective replay {auc_s} vs {rec['auroc_dual_vs_selective']}")
            if abs(auc_n - float(rec["auroc_dual_vs_neither"])) > 1e-4:
                fail(f"{pair} {contrast} neither replay {auc_n} vs {rec['auroc_dual_vs_neither']}")
            if abs((auc_n - auc_s) - float(rec["delta_neither_minus_selective"])) > 1e-4:
                fail(f"{pair} {contrast} delta replay")


def check_ranking(packs, rank_rows, two_rows, and_rows):
    rank = {r["pair"]: r for r in rank_rows}
    two = {r["pair"]: r for r in two_rows}
    and_by = {r["pair"]: r for r in and_rows}
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
        dual_m = [r["score_mean"] for r in recs if r["cls"] == "dual"]
        nei_m = [r["score_mean"] for r in recs if r["cls"] == "neither"]
        pt = auroc(dual_m, nei_m)
        if abs(pt - float(two[pair]["two_pocket_mean_D_vs_neither"])) > 1e-4:
            fail(f"{pair} two-pocket mean {pt} vs {two[pair]['two_pocket_mean_D_vs_neither']}")
        dual_w = [r["score_worst"] for r in recs if r["cls"] == "dual"]
        thresh = float(np.median(np.asarray(dual_w, dtype=float)))
        usable = [r for r in recs if r["cls"] in {"dual", "A_only", "B_only"}]
        kept = [r for r in usable if r["score_worst"] >= thresh]
        if abs(thresh - float(and_by[pair]["and_threshold"])) > 1e-4:
            fail(f"{pair} AND threshold {thresh} vs {and_by[pair]['and_threshold']}")
        if int(and_by[pair]["and_n_pass"]) != len(kept):
            fail(f"{pair} AND n_pass {and_by[pair]['and_n_pass']} vs {len(kept)}")


def check_max_median(maxmed):
    ache = [r for r in maxmed if r["pair"] == "AChE/BChE" and r["aggregation"] == "max"]
    if not ache or int(ache[0]["n_ligands"]) != 94:
        fail(f"AChE max/median n={ache}")
    if (int(ache[0]["n_dual"]), int(ache[0]["n_A_only"]), int(ache[0]["n_B_only"])) != (27, 25, 28):
        fail(f"AChE max class counts {ache[0]}")
    egfr = [r for r in maxmed if r["pair"] == "EGFR/HER2" and r["aggregation"] == "max"]
    if not egfr or int(egfr[0]["n_ligands"]) != 109:
        fail(f"EGFR max/median n={egfr}")


def check_models(packs, inc, folds, oof):
    vals = [abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])) for r in inc]
    if not vals:
        fail("empty ECFP incremental table")
        return 0.0
    max_abs = max(vals)
    loc = inc[vals.index(max_abs)]
    NOTES.append(f"ECFP max |Δ|={max_abs:.4f} at {loc['pair']} {loc['contrast']}")
    fold_keys = Counter((r["pair"], r["arm"], r["ligand_id"]) for r in folds)
    if any(n != 1 for n in fold_keys.values()):
        fail("duplicate fold assignments")
    by_fold = {(r["pair"], r["arm"], r["ligand_id"]): r for r in folds}
    scaf_fold = defaultdict(set)
    for r in folds:
        scaf_fold[(r["pair"], r["arm"], r["scaffold"])].add(r["fold_id"])
    for key, fset in scaf_fold.items():
        if len(fset) != 1:
            fail(f"scaffold split across folds {key} {fset}")
    oof_by = defaultdict(list)
    models_by = defaultdict(set)
    ligands_by = defaultdict(set)
    for r in oof:
        if str(r.get("scaled", "0")) not in ("0", "False"):
            continue
        key = (r["pair"], r["arm"], r["model"])
        oof_by[key].append(r)
        models_by[(r["pair"], r["arm"])].add(r["model"])
        ligands_by[(r["pair"], r["arm"], r["model"])].add(r["ligand_id"])
    needed = {"ECFP4", "docking", "ECFP4+docking"}
    inc_by = {(r["pair"], r["contrast"]): r for r in inc}
    for pair in PRIMARY_PAIRS:
        recs = packs[pair]
        for arm, pos_cls, neg_cls in (("D_vs_A", "dual", "A_only"), ("D_vs_B", "dual", "B_only")):
            kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            if models_by[(pair, arm)] != needed:
                fail(f"{pair} {arm} models {models_by[(pair, arm)]}")
            ids = {r["ligand_id"] for r in kept}
            for model in needed:
                if ligands_by[(pair, arm, model)] != ids:
                    fail(f"{pair} {arm} {model} ligand set mismatch")
            rec = inc_by.get((pair, arm))
            if rec is None:
                fail(f"{pair} {arm} missing incremental row")
                continue
            if int(rec["n_pos"]) != sum(1 for r in kept if r["cls"] == pos_cls):
                fail(f"{pair} {arm} n_pos {rec['n_pos']} vs master")
            if int(rec["n_neg"]) != sum(1 for r in kept if r["cls"] == neg_cls):
                fail(f"{pair} {arm} n_neg {rec['n_neg']} vs master")
            replay = {}
            for model in needed:
                rows = oof_by[(pair, arm, model)]
                if len(rows) != len(kept):
                    fail(f"{pair} {arm} {model} OOF n={len(rows)} vs {len(kept)}")
                    continue
                y = np.array([int(r["y"]) for r in rows], dtype=int)
                p = np.array([float(r["oof_prob"]) for r in rows], dtype=float)
                replay[model] = float(roc_auc_score(y, p))
            if abs(replay["ECFP4"] - float(rec["cv_auroc_ECFP4"])) > 1e-4:
                fail(f"{pair} {arm} ECFP4 OOF AUROC {replay['ECFP4']} vs {rec['cv_auroc_ECFP4']}")
            if abs(replay["ECFP4+docking"] - float(rec["cv_auroc_ECFP4_docking"])) > 1e-4:
                fail(f"{pair} {arm} ECFP4+dock OOF AUROC {replay['ECFP4+docking']} vs {rec['cv_auroc_ECFP4_docking']}")
            delta = replay["ECFP4+docking"] - replay["ECFP4"]
            if abs(delta - float(rec["delta_ECFP4_plus_docking_minus_ECFP4"])) > 1e-4:
                fail(f"{pair} {arm} delta replay {delta} vs {rec['delta_ECFP4_plus_docking_minus_ECFP4']}")
            for r in rows:
                fk = by_fold.get((pair, arm, r["ligand_id"]))
                if fk is None or int(fk["fold_id"]) < 0:
                    fail(f"{pair} {arm} {r['ligand_id']} missing/negative fold")
    replay_max = max(abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])) for r in inc)
    replay_loc = max(inc, key=lambda r: abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])))
    if abs(replay_max - max_abs) > 1e-12:
        fail(f"max |Δ| inconsistency {replay_max} vs {max_abs}")
    NOTES.append(f"ECFP replay max |Δ|={replay_max:.4f} at {replay_loc['pair']} {replay_loc['contrast']}")
    return max_abs


def check_pocket(packs, matched, hold):
    hold_by = {r["pair"]: r for r in hold}
    mat_by = {r["pair"]: r for r in matched}
    ache_excl = False
    egfr_excl = True
    hold_excl = []
    for pair in PRIMARY_PAIRS:
        rec = mat_by[pair]
        lo, hi = float(rec["delta_ci_lo"]), float(rec["delta_ci_hi"])
        excl = not (lo <= 0 <= hi)
        if pair == "AChE/BChE":
            ache_excl = excl
        if pair == "EGFR/HER2":
            egfr_excl = excl
        if pair == "EGFR/HER2":
            if pair in hold_by:
                fail("EGFR has a holdout row")
            continue
        h = hold_by.get(pair)
        if h is None:
            fail(f"{pair} holdout missing")
            continue
        hlo, hhi = float(h["mm_ci_lo"]), float(h["mm_ci_hi"])
        hold_excl.append(not (hlo <= 0 <= hhi))
    if not ache_excl:
        fail("AChE/BChE main-panel matched-mismatched CI includes 0")
    if egfr_excl:
        fail("EGFR/HER2 main-panel matched-mismatched CI excludes 0")
    if any(hold_excl):
        fail(f"holdout interval excludes 0: {hold_excl}")


def check_robustness(robust, rec_sub, cluster, rmsd, det):
    engines = {(r["pair"], r["engine"]) for r in robust}
    for key in (
        ("EGFR/HER2", "gnina_dock_mode1"),
        ("PIK3CA/mTOR", "gnina_dock_mode1"),
        ("JAK1/TYK2", "gnina_dock_mode1"),
        ("PIK3CA/mTOR", "vina_alt_4JPS"),
        ("PIK3CA/mTOR", "vina_alt_5DXT"),
        ("PIK3CA/mTOR", "vina_alt_4JSX"),
        ("PPARG/PPARA", "rtm_best9"),
        ("PPARG/PPARA", "gnina_cnn_affinity"),
    ):
        if key not in engines:
            fail(f"missing robustness row {key}")
    rtm = next((r for r in robust if r["pair"] == "PPARG/PPARA" and r["engine"] == "rtm_best9"), None)
    cnn = next((r for r in robust if r["pair"] == "PPARG/PPARA" and r["engine"] == "gnina_cnn_affinity"), None)
    if rtm is None or abs(float(rtm["summary_min"]) - 0.3691) > 1e-3:
        fail(f"PPARG RTM summary_min {rtm}")
    if cnn is None or abs(float(cnn["summary_min"]) - 0.5000) > 1e-3:
        fail(f"PPARG CNN summary_min {cnn}")
    if "not independent" not in (rtm.get("note") or "").lower():
        fail("RTM note does not mark same-pose rescoring")
    if len(rec_sub) != 3:
        fail(f"receptor_substitution n={len(rec_sub)}")
    primary_ref = {r["replacement"]: r["primary_summary_min"] for r in rec_sub}
    if len(set(primary_ref.values())) != 1:
        fail(f"receptor substitution primary reference not unique {primary_ref}")
    jak_doc = [r for r in cluster if r["pair"] == "JAK1/TYK2" and r["estimator"] == "document_cluster"]
    if not jak_doc or jak_doc[0]["status"] != "unresolved_mapping_unavailable":
        fail(f"JAK document cluster {jak_doc}")
    if any(r.get("pdb") == "3POZ" and abs(float(r.get("calcrrms_best_A") or r.get("best_A") or 0) - 0.760) < 1e-6 for r in rmsd):
        fail("PR #35 reconstructed EGFR 0.760 Å used as current RMSD")
    by_prot = {r.get("protein"): r for r in rmsd}
    egfr = by_prot.get("EGFR")
    her2 = by_prot.get("HER2")
    if egfr is None or abs(float(egfr["calcrrms_top1_A"]) - 1.019) > 1e-3:
        fail(f"EGFR top-1 RMSD {egfr}")
    if her2 is None or abs(float(her2["calcrrms_top1_A"]) - 1.947) > 1e-3:
        fail(f"HER2 top-1 RMSD {her2}")
    for prot in ("JAK2", "mTOR", "BChE", "PPARG", "PPARA"):
        rec = by_prot.get(prot)
        if rec is None:
            fail(f"{prot} RMSD missing")
            continue
        top1 = float(rec["calcrrms_top1_A"])
        best = float(rec["calcrrms_best_A"])
        if top1 < 2.0:
            fail(f"{prot} top-1 {top1} unexpectedly < 2 Å")
        if best >= 2.0:
            fail(f"{prot} lowest saved pose {best} not < 2 Å")
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


def write_report(freeze: Path):
    if freeze.resolve() == (ROOT / "results" / "canonical").resolve():
        report = Path("/tmp/canonical_verification_report.md")
    else:
        report = freeze / "VERIFICATION_REPORT.md"
    lines = ["# Freeze verification", ""]
    lines.append(f"freeze_dir: `{freeze}`")
    lines.append(f"result: {'FAIL' if FAILS else 'PASS'}")
    lines.append("")
    if FAILS:
        lines.append("## Failures")
        lines.extend(f"- {m}" for m in FAILS)
        lines.append("")
    lines.append("## Notes")
    lines.extend(f"- {m}" for m in NOTES)
    lines.append("")
    report.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


def main() -> int:
    global FREEZE
    parser = argparse.ArgumentParser(description="Independently verify a freeze rebuild directory.")
    parser.add_argument("--outdir", default="/tmp/dual_target_freeze_rebuild")
    args = parser.parse_args()
    FREEZE = Path(args.outdir)
    if not FREEZE.is_absolute():
        FREEZE = (ROOT / FREEZE).resolve()
    if not FREEZE.is_dir():
        print("FAIL: freeze directory missing", FREEZE)
        return 1
    check_parse_finite()
    check_expected_files(FREEZE)
    master = read_csv(FREEZE / "current_score_master.csv")
    check_master(master)
    packs = load_packs(master)
    check_table2(packs, read_csv(FREEZE / "primary_summary_min.csv"), read_csv(FREEZE / "primary_directional_auroc.csv"))
    check_fixed_score(packs, read_csv(FREEZE / "fixed_score_negative_class_delta.csv"))
    check_ranking(
        packs,
        read_csv(FREEZE / "top10_operating_points.csv"),
        read_csv(FREEZE / "two_pocket_mean_ranking.csv"),
        read_csv(FREEZE / "and_filter_operating_points.csv"),
    )
    check_max_median(read_csv(FREEZE / "max_vs_median_sensitivity.csv"))
    check_models(
        packs,
        read_csv(FREEZE / "ecfp4_incremental_information.csv"),
        read_csv(FREEZE / "model_fold_assignments.csv"),
        read_csv(FREEZE / "ecfp4_oof_predictions.csv"),
    )
    check_pocket(packs, read_csv(FREEZE / "matched_mismatched_pocket.csv"), read_csv(FREEZE / "holdout_metrics.csv"))
    check_robustness(
        read_csv(FREEZE / "computational_robustness.csv"),
        read_csv(FREEZE / "receptor_substitution.csv"),
        read_csv(FREEZE / "cluster_bootstrap_sensitivity.csv"),
        read_csv(FREEZE / "cognate_rmsd.csv"),
        read_csv(FREEZE / "detectable_effect_simulation.csv"),
    )
    write_report(FREEZE)
    return 1 if FAILS else 0


if __name__ == "__main__":
    raise SystemExit(main())
