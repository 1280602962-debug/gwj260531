#!/usr/bin/env python3
"""Expose directional CIs and compare bootstrap schemes on deposited scores.

Replays the locked Table 2 draws without overwriting their source tables.
The additional class-stratified analysis was introduced during manuscript
revision; it is not a prespecified replacement for the primary analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import platform
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data/jcim_novelty_v0/tables"
UNIV = ROOT / "data/jcim_chembl_universe_v0"
LOCAL = UNIV / "local_track_b_v0"
N_BOOT = 2000
SEED = 20260729
ORDER = ("EGFR/HER2", "JAK1/JAK2", "JAK1/TYK2", "PIK3CA/mTOR",
         "AChE/BChE", "F2/F10", "PPARG/PPARA", "PPARA/PPARD")
CLASSES = ("dual", "A_only", "B_only")


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def offset(*parts):
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:16], 16) % 99991


def auc(pos, neg):
    if not len(pos) or not len(neg):
        return float("nan")
    d = np.asarray(pos)[:, None] - np.asarray(neg)[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size)


def metrics(recs):
    group = {c: [r for r in recs if r["cls"] == c] for c in CLASSES}
    da = auc([r["vina_B"] for r in group["dual"]], [r["vina_B"] for r in group["A_only"]])
    db = auc([r["vina_A"] for r in group["dual"]], [r["vina_A"] for r in group["B_only"]])
    return np.array([da, db, min(da, db)])


def bootstrap(recs, seed, stratified=False):
    use = [r for r in recs if r["cls"] in CLASSES]
    rng = np.random.default_rng(seed)
    idx = np.arange(len(use))
    class_idx = [np.array([i for i, r in enumerate(use) if r["cls"] == c]) for c in CLASSES]
    assert all(len(a) for a in class_idx)
    draws = []
    for _ in range(N_BOOT):
        if stratified:
            # One dual draw is shared by both score channels and both arms.
            ii = np.concatenate([rng.choice(a, len(a), replace=True) for a in class_idx])
        else:
            ii = rng.choice(idx, len(idx), replace=True)
        sub = [use[int(i)] for i in ii]
        value = metrics(sub)
        if np.isfinite(value).all():
            draws.append(value)
    assert len(draws) >= N_BOOT // 2
    return np.percentile(draws, [2.5, 97.5], axis=0), len(draws)


def side(lo, hi):
    return "above_0.5" if lo > 0.5 else "below_0.5" if hi < 0.5 else "includes_0.5"


def build():
    m3 = module("review_t0", ROOT / "data/jcim_strengthen_t0t1_v0/scripts/build_t0_strengthen_v1.py")
    sys.path.insert(0, str(UNIV / "scripts"))
    import analyze_five_pair_stack_v1 as m5
    import analyze_five_pair_local_channels_v1 as mc
    packs = {}
    sources = [Path(__file__), Path(m3.__file__), Path(m5.__file__), Path(mc.__file__)]
    for pair, recs in m3.assemble_records().items():
        packs[pair] = [{**r, "cls": m3.assign_fourclass(r["pA"], r["pB"], 6.0)}
                       for r in recs if r["pA"] is not None and r["pB"] is not None]
        sources.append(ROOT / m3.PAIR_SPEC[pair]["scores"])
        if m3.PAIR_SPEC[pair]["panel"]:
            sources.append(ROOT / m3.PAIR_SPEC[pair]["panel"])
    scores = m5.load_scores()
    sources.append(LOCAL / "tables/scores_vina_mode1_v1.csv")
    for spec in m5.PAIRS:
        packs[spec["pair"]], _ = m5.assemble(spec, scores, {})
        sources.append(spec["panel"])

    lock3_path = ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
    lock5_path = LOCAL / "tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv"
    sources += [lock3_path, lock5_path]
    locked = {r["pair"]: r for r in rows(lock3_path) if r["label_rule"] == "theta_6.0"}
    locked.update({r["pair"]: r for r in rows(lock5_path)})
    primary, sensitivity, membership, checks = [], [], [], []
    for pair in ORDER:
        recs = packs[pair]
        assert len({r["ligand"] for r in recs}) == len(recs), pair
        assert all(np.isfinite([r["vina_A"], r["vina_B"], r["pA"], r["pB"]]).all() for r in recs), pair
        use = [r for r in recs if r["cls"] in CLASSES]
        counts = Counter(r["cls"] for r in recs)
        lock = locked[pair]
        for cls, col in zip(CLASSES, ("n_dual", "n_A_only", "n_B_only")):
            assert counts[cls] == int(lock[col]), (pair, cls, counts[cls], lock[col])
        point = metrics(use)
        expected = [float(lock.get("auroc_D_vs_A") or lock["auroc_D_vs_A_pocketB"]),
                    float(lock.get("auroc_D_vs_B") or lock["auroc_D_vs_B_pocketA"]),
                    float(lock.get("pocket_matched_summary_min") or lock["summary_min"])]
        assert np.allclose(point, expected, atol=0.000051, rtol=0), (pair, point, expected)
        seed = SEED + offset(pair, "theta_6.0")
        ci, n = bootstrap(use, seed)
        assert np.allclose(ci[:, 2], [float(lock["ci_lo"]), float(lock["ci_hi"])], atol=0.000051, rtol=0), (pair, ci)
        primary.append(dict(pair=pair, n_dual=counts["dual"], n_A_only=counts["A_only"],
                            n_B_only=counts["B_only"], n_neither=counts["neither"],
                            auroc_D_vs_A=point[0], D_vs_A_ci_lo=ci[0, 0], D_vs_A_ci_hi=ci[1, 0],
                            auroc_D_vs_B=point[1], D_vs_B_ci_lo=ci[0, 1], D_vs_B_ci_hi=ci[1, 1],
                            summary_min=point[2], summary_min_ci_lo=ci[0, 2], summary_min_ci_hi=ci[1, 2],
                            bootstrap="pooled_non_stratified", n_boot=N_BOOT, n_valid=n, seed=seed))
        strat_seed = SEED + offset(pair, "theta_6.0", "class_stratified_review")
        sci, sn = bootstrap(use, strat_seed, stratified=True)
        sensitivity.append(dict(pair=pair, summary_min=point[2],
                                pooled_ci_lo=ci[0, 2], pooled_ci_hi=ci[1, 2],
                                stratified_ci_lo=sci[0, 2], stratified_ci_hi=sci[1, 2],
                                pooled_relation_to_0_5=side(*ci[:, 2]),
                                stratified_relation_to_0_5=side(*sci[:, 2]),
                                n_boot=N_BOOT, n_valid=sn, seed=strat_seed,
                                D_vs_A_ci_lo=sci[0, 0], D_vs_A_ci_hi=sci[1, 0],
                                D_vs_B_ci_lo=sci[0, 1], D_vs_B_ci_hi=sci[1, 1]))
        membership += [dict(pair=pair, ligand=r["ligand"], cls=r["cls"], pA=r["pA"], pB=r["pB"],
                            score_A=r["vina_A"], score_B=r["vina_B"]) for r in recs]
        checks.append(dict(check="canonical_point_count_and_summary_CI_replay", pair=pair, status="PASS"))

    # Verify JAK1/TYK2 GNINA interval provenance, retaining the original joint CI.
    gpath = LOCAL / "tables/scores_gnina_independent_jak1_tyk2_v1.csv"
    gmap, _ = mc.load_status_scores(gpath, "pair", "score_S", ok=("success",))
    spec = next(s for s in m5.PAIRS if s["pair"] == "JAK1/TYK2")
    grec, _ = mc.recs_from_panel(spec, gmap[spec["pair"]])
    gseed = SEED + offset(spec["pair"], "gnina_independent_jak1_tyk2", "theta_6.0")
    gci, gn = bootstrap(grec, gseed)
    gpoint = metrics(grec)
    gtable = LOCAL / "tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv"
    grow = next(r for r in rows(gtable) if r["pair"] == spec["pair"] and r["channel"] == "gnina_independent_jak1_tyk2")
    assert np.allclose(gci[:, 2], [float(grow["ci_lo"]), float(grow["ci_hi"])], atol=0.000051, rtol=0)
    gnina = [dict(pair=spec["pair"], metric=metric, point=gpoint[i], ci_lo=gci[0, i], ci_hi=gci[1, i],
                  bootstrap="pooled_non_stratified_shared_dual", seed=gseed, n_valid=gn)
             for i, metric in enumerate(("D_vs_A_pocketB", "D_vs_B_pocketA", "summary_min"))]
    sources += [gpath, gtable]
    checks.append(dict(check="GNINA_JAK1_TYK2_summary_CI_is_joint_minimum", status="PASS"))

    # Identical retrospective ranking rules for the original and standard-supply examples.
    operations = []
    for pair in ("EGFR/HER2", "JAK1/TYK2"):
        recs = packs[pair]
        top = sorted(recs, key=lambda r: (-r["vina_mean"], r["ligand"]))[:10]
        threshold = float(np.median([r["vina_worst"] for r in recs if r["cls"] == "dual"]))
        use = [r for r in recs if r["cls"] in CLASSES]
        retained = [r for r in use if r["vina_worst"] >= threshold]
        ct = Counter(r["cls"] for r in top)
        cr = Counter(r["cls"] for r in retained)
        operations.append(dict(pair=pair, n_ranked=len(recs), top_k=10, top_dual=ct["dual"],
                               top_A_only=ct["A_only"], top_B_only=ct["B_only"], top_neither=ct["neither"],
                               top_ligand_ids=";".join(r["ligand"] for r in top), threshold=threshold,
                               n_filter_input=len(use), retained_dual=cr["dual"], retained_A_only=cr["A_only"],
                               retained_B_only=cr["B_only"], dual_recall=cr["dual"]/sum(r["cls"] == "dual" for r in use),
                               dual_precision=cr["dual"]/len(retained),
                               filter_rule="score_worst >= median_dual_score_worst; neither excluded",
                               tie_rule="descending score_mean then ascending ligand ID"))

    def csv_text(data):
        f = io.StringIO(newline="")
        w = csv.DictWriter(f, fieldnames=list(data[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(data)
        return f.getvalue()

    import rdkit
    provenance = dict(n_boot=N_BOOT, base_seed=SEED,
                      environment=dict(python=platform.python_version(), numpy=np.__version__, rdkit=rdkit.__version__),
                      input_sha256={str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(set(sources))}, checks=checks,
                      interpretation="Pointwise percentile intervals; no multiplicity correction. Stratified bootstrap is a revision sensitivity analysis.")
    out = {"primary_directional_intervals_review_v1.csv": csv_text(primary),
           "summary_min_stratified_sensitivity_review_v1.csv": csv_text(sensitivity),
           "review_scored_membership_v1.csv": csv_text(membership),
           "gnina_jak_interval_provenance_review_v1.csv": csv_text(gnina),
           "operating_point_examples_review_v1.csv": csv_text(operations),
           "review_statistics_provenance_v1.json": json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="recompute in memory and compare without writing")
    args = parser.parse_args()
    out = build()
    for name, content in out.items():
        path = TAB / name
        if args.check:
            assert path.read_text(encoding="utf-8") == content, f"replay differs: {path}"
        else:
            path.write_text(content, encoding="utf-8", newline="")
        print(("checked " if args.check else "wrote ") + name)
    print("PASS: eight canonical points, counts and pooled summary intervals; GNINA joint CI provenance")


if __name__ == "__main__":
    main()
