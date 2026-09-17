#!/usr/bin/env python3
"""Build the unique current per-ligand score master (zero-dock).

EGFR/HER2: corrected-box ablation scores (not review_scored_membership).
AChE/BChE: corrected panel + scores.
PIK3CA/mTOR: PM48 rdkit production scores.
Track B five pairs: scores_vina_mode1_v1.csv.

Writes:
  data/processed/current_score_master.csv
  results/canonical/current_score_master.csv
  results/canonical/score_master_migration_diff.md
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.bootstrap_metrics import assign_fourclass  # noqa: E402

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
RECEPTORS = {
    "EGFR/HER2": ("3POZ", "3RCD"),
    "JAK1/JAK2": ("6N7A", "8BXH"),
    "JAK1/TYK2": ("6N7A", "3LXP"),
    "PIK3CA/mTOR": ("4L23", "4JT6"),
    "AChE/BChE": ("4EY7", "4BDS"),
    "F2/F10": ("4UDW", "2JKH"),
    "PPARG/PPARA": ("9V8H", "6LXA"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}
FIELDS = [
    "pair",
    "ligand_id",
    "molecule_chembl_id",
    "smiles",
    "construction_class",
    "primary_class_theta6",
    "pA",
    "pB",
    "score_A",
    "score_B",
    "affinity_A",
    "affinity_B",
    "complete_case",
    "analysis_set",
    "receptor_A",
    "receptor_B",
    "score_source",
    "postfix_status",
    "activity_eligible",
    "activity_status",
    "historical_pA",
    "historical_pB",
    "n_act_A",
    "n_act_B",
]


def fnum(v):
    if v is None or v == "":
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    return x


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row_of(**kwargs) -> dict:
    rec = {k: "" for k in FIELDS}
    rec.update(kwargs)
    recA = RECEPTORS[rec["pair"]]
    rec["receptor_A"], rec["receptor_B"] = recA
    pa, pb = fnum(rec.get("pA")), fnum(rec.get("pB"))
    rec["pA"] = "" if pa is None else pa
    rec["pB"] = "" if pb is None else pb
    rec["primary_class_theta6"] = assign_fourclass(pa, pb) or ""
    sa, sb = fnum(rec.get("score_A")), fnum(rec.get("score_B"))
    rec["complete_case"] = int(sa is not None and sb is not None)
    rec["activity_eligible"] = rec.get("activity_eligible") if rec.get("activity_eligible") not in ("", None) else 1
    rec["activity_status"] = rec.get("activity_status") or "panel_default"
    rec["historical_pA"] = rec.get("historical_pA") if rec.get("historical_pA") not in ("", None) else rec["pA"]
    rec["historical_pB"] = rec.get("historical_pB") if rec.get("historical_pB") not in ("", None) else rec["pB"]
    rec["n_act_A"] = rec.get("n_act_A") if rec.get("n_act_A") not in ("", None) else ""
    rec["n_act_B"] = rec.get("n_act_B") if rec.get("n_act_B") not in ("", None) else ""
    if sa is not None:
        rec["score_A"] = sa
        rec["affinity_A"] = rec.get("affinity_A") if rec.get("affinity_A") not in ("", None) else -sa
    if sb is not None:
        rec["score_B"] = sb
        rec["affinity_B"] = rec.get("affinity_B") if rec.get("affinity_B") not in ("", None) else -sb
    return rec


def load_egfr() -> list[dict]:
    scores = {r["ligand"]: r for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")}
    panel = {r["panel_id"]: r for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")}
    out = []
    for lig, s in scores.items():
        p = panel.get(lig, {})
        sa = fnum(s.get("vina_3POZ_hb"))
        sb = fnum(s.get("vina_3RCD_hb"))
        out.append(
            row_of(
                pair="EGFR/HER2",
                ligand_id=lig,
                molecule_chembl_id=s.get("molecule_chembl_id") or p.get("molecule_chembl_id", ""),
                smiles=p.get("smiles", ""),
                construction_class=s.get("class") or p.get("class", ""),
                pA=p.get("pchembl_EGFR"),
                pB=p.get("pchembl_HER2"),
                score_A=sa,
                score_B=sb,
                affinity_A=fnum(s.get("3POZ_affinity")),
                affinity_B=fnum(s.get("3RCD_affinity")),
                analysis_set="main",
                score_source="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
                postfix_status="corrected_box",
            )
        )
    return out


def load_ache() -> list[dict]:
    scores = {r["ligand"]: r for r in read_csv(ROOT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv")}
    panel = {r["panel_id"]: r for r in read_csv(ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict.csv")}
    out = []
    for lig, p in panel.items():
        s = scores.get(lig, {})
        sa = fnum(s.get("vina_ACHE_hb"))
        sb = fnum(s.get("vina_BCHE_hb"))
        out.append(
            row_of(
                pair="AChE/BChE",
                ligand_id=lig,
                molecule_chembl_id=p.get("molecule_chembl_id", ""),
                smiles=p.get("smiles") or s.get("smiles", ""),
                construction_class=p.get("class", ""),
                pA=p.get("pchembl_ACHE"),
                pB=p.get("pchembl_BCHE"),
                score_A=sa,
                score_B=sb,
                affinity_A=fnum(s.get("vina_ACHE")),
                affinity_B=fnum(s.get("vina_BCHE")),
                analysis_set="main",
                score_source="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
                postfix_status="corrected_panel",
            )
        )
    return out


def load_pm48() -> list[dict]:
    scores = {r["ligand"]: r for r in read_csv(ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv")}
    panel = {r["panel_id"]: r for r in read_csv(ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv")}
    out = []
    for lig, s in scores.items():
        p = panel.get(lig, {})
        sa = fnum(s.get("vina_4L23_hb"))
        sb = fnum(s.get("vina_4JT6_hb"))
        out.append(
            row_of(
                pair="PIK3CA/mTOR",
                ligand_id=lig,
                molecule_chembl_id=s.get("molecule_chembl_id") or p.get("molecule_chembl_id", ""),
                smiles=p.get("smiles", ""),
                construction_class=s.get("class") or p.get("class", ""),
                pA=p.get("pchembl_PIK3CA"),
                pB=p.get("pchembl_MTOR"),
                score_A=sa,
                score_B=sb,
                affinity_A=fnum(s.get("4L23_affinity")),
                affinity_B=fnum(s.get("4JT6_affinity")),
                analysis_set="main",
                score_source="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
                postfix_status="production",
            )
        )
    return out


TRACK_B_PANELS = {
    "JAK1/JAK2": "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_JAK2_v1.csv",
    "JAK1/TYK2": "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_TYK2_v1.csv",
    "F2/F10": "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_F2_F10_v1.csv",
    "PPARG/PPARA": "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARG_PPARA_v1.csv",
    "PPARA/PPARD": "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARA_PPARD_v1.csv",
}
TRACK_B_HOLDOUT = {
    "JAK1/JAK2": "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOJ1J2_v1.csv",
    "JAK1/TYK2": "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOJ1TYK2_v1.csv",
    "F2/F10": "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOF2F10_v1.csv",
    "PPARG/PPARA": "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOPGPA_v1.csv",
    "PPARA/PPARD": "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOPAPD_v1.csv",
}


def pivot_track_b(path: Path, pair: str) -> dict[str, dict]:
    pdb_a, pdb_b = RECEPTORS[pair]
    wide: dict[str, dict] = defaultdict(dict)
    for r in read_csv(path):
        if r.get("pair") != pair:
            continue
        if r.get("status") and r["status"] not in {"success", ""}:
            continue
        lig = r["ligand"]
        tgt = r["target"]
        sc = fnum(r.get("score_S"))
        en = fnum(r.get("mode1_energy"))
        if tgt == pdb_a:
            wide[lig]["score_A"] = sc
            wide[lig]["affinity_A"] = en
        elif tgt == pdb_b:
            wide[lig]["score_B"] = sc
            wide[lig]["affinity_B"] = en
    return wide


def load_track_b() -> list[dict]:
    scores_path = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv"
    out = []
    for pair, panel_rel in TRACK_B_PANELS.items():
        panel = {r["panel_id"]: r for r in read_csv(ROOT / panel_rel)}
        wide = pivot_track_b(scores_path, pair)
        for lig, p in panel.items():
            s = wide.get(lig, {})
            out.append(
                row_of(
                    pair=pair,
                    ligand_id=lig,
                    molecule_chembl_id=p.get("molecule_chembl_id", ""),
                    smiles=p.get("canonical_smiles", ""),
                    construction_class=p.get("class") or p.get("theta6_class", ""),
                    pA=p.get("pchembl_A"),
                    pB=p.get("pchembl_B"),
                    score_A=s.get("score_A"),
                    score_B=s.get("score_B"),
                    affinity_A=s.get("affinity_A"),
                    affinity_B=s.get("affinity_B"),
                    analysis_set="main",
                    score_source="data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv",
                    postfix_status="production",
                )
            )
    return out


def load_holdout_ab_pm() -> list[dict]:
    rows = read_csv(ROOT / "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv")
    out = []
    src = "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv"
    for r in rows:
        pair = r["pair"]
        if pair not in RECEPTORS:
            continue
        sa, sb = fnum(r.get("vina_A")), fnum(r.get("vina_B"))
        status = "corrected_panel" if pair == "AChE/BChE" else "production"
        out.append(
            row_of(
                pair=pair,
                ligand_id=r["ligand"],
                molecule_chembl_id=r.get("chembl", ""),
                smiles=r.get("smiles", ""),
                construction_class=r.get("cls", ""),
                pA="",
                pB="",
                score_A=sa,
                score_B=sb,
                affinity_A=fnum(r.get("vina_A_raw")),
                affinity_B=fnum(r.get("vina_B_raw")),
                analysis_set="holdout",
                score_source=src,
                postfix_status=status,
            )
        )
    return out


def load_holdout_track_b() -> list[dict]:
    scores = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/allpairs_stack/holdout/tables/holdout_scores_vina_mode1_v1.csv"
    out = []
    for pair, panel_rel in TRACK_B_HOLDOUT.items():
        panel_path = ROOT / panel_rel
        if not panel_path.is_file() or not scores.is_file():
            continue
        panel = {r["holdout_id"]: r for r in read_csv(panel_path)}
        wide = pivot_track_b(scores, pair)
        for lig, p in panel.items():
            s = wide.get(lig, {})
            out.append(
                row_of(
                    pair=pair,
                    ligand_id=lig,
                    molecule_chembl_id=p.get("molecule_chembl_id", ""),
                    smiles=p.get("canonical_smiles", ""),
                    construction_class=p.get("class", ""),
                    pA=p.get("pchembl_A"),
                    pB=p.get("pchembl_B"),
                    score_A=s.get("score_A"),
                    score_B=s.get("score_B"),
                    affinity_A=s.get("affinity_A"),
                    affinity_B=s.get("affinity_B"),
                    analysis_set="holdout",
                    score_source=str(scores.relative_to(ROOT)),
                    postfix_status="production",
                )
            )
    return out


def compare_membership(master: list[dict]) -> str:
    lines = ["# Score master migration diff", "", "One-time report. Authoritative file: `results/canonical/current_score_master.csv`.", ""]
    review = read_csv(ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv")
    by = {(r["pair"], r["ligand_id"]): r for r in master if r["analysis_set"] == "main"}
    mismatches = []
    missing = []
    for r in review:
        key = (r["pair"], r["ligand"])
        m = by.get(key)
        if m is None:
            missing.append(key)
            continue
        if not m["complete_case"]:
            continue
        old_a, old_b = fnum(r["score_A"]), fnum(r["score_B"])
        new_a, new_b = fnum(m["score_A"]), fnum(m["score_B"])
        if old_a is None or new_a is None:
            continue
        if abs(old_a - new_a) > 1e-6 or abs(old_b - new_b) > 1e-6:
            mismatches.append((key, old_a, old_b, new_a, new_b))
    by_pair = Counter(k[0] for k, *_ in mismatches)
    lines.append("## review_scored_membership_v1.csv vs master")
    lines.append(f"- ligands in review: {len(review)}")
    lines.append(f"- missing from master main: {len(missing)}")
    lines.append(f"- score mismatches (abs > 1e-6): {len(mismatches)}")
    for pair, n in sorted(by_pair.items()):
        lines.append(f"  - {pair}: {n}")
    if by_pair:
        lines.append("- **EGFR/HER2 review scores are pre-fix** if EGFR is in the mismatch list. Do not sort or bootstrap from that file.")
    lines.append("")
    # EGFR vs legacy
    legacy = ROOT / "data/_legacy_archive/egfr_her2_scores_pre_fix/ablation_ligand_scores.csv"
    if legacy.is_file():
        old = {r["ligand"]: r for r in read_csv(legacy)}
        n_match_legacy = 0
        n_match_review = 0
        for r in review:
            if r["pair"] != "EGFR/HER2":
                continue
            o = old.get(r["ligand"])
            if o and abs(fnum(r["score_A"]) - fnum(o["vina_3POZ_hb"])) < 1e-6:
                n_match_legacy += 1
        lines.append(f"- EGFR review rows matching pre-fix ablation vina_*_hb: {n_match_legacy}")
        _ = n_match_review
    lines.append("")
    lines.append("## Files that can be regenerated from the master")
    lines.append("- eight-pair ranking / EF tables")
    lines.append("- directional AUROC / summary_min / scheme-B CIs")
    lines.append("- fixed-score ΔAUROC")
    lines.append("- matched-minus-mismatched")
    lines.append("- descriptor and ECFP4 models (SMILES come with the master)")
    lines.append("")
    lines.append("## Files that still contain old EGFR scores (do not use as current analysis input)")
    lines.append("- `data/jcim_novelty_v0/tables/review_scored_membership_v1.csv` (EGFR pre-fix scores)")
    lines.append("- `data/_legacy_archive/egfr_her2_scores_pre_fix/ablation_ligand_scores.csv`")
    lines.append("- `data/egfr_her2_panel120_v0/tables/scores_vina.csv` (stale sibling; production is ablation_ligand_scores.csv)")
    lines.append("")
    main = [r for r in master if r["analysis_set"] == "main"]
    lines.append("## Master complete-case counts (main, theta=6 class)")
    for pair in PRIMARY_PAIRS:
        recs = [
            r
            for r in main
            if r["pair"] == pair
            and r["complete_case"] in (1, "1")
            and str(r.get("activity_eligible", "1")) in ("1", "True")
        ]
        c = Counter(r["primary_class_theta6"] for r in recs)
        n_score = sum(1 for r in main if r["pair"] == pair and r["complete_case"] in (1, "1"))
        lines.append(
            f"- {pair}: scored={n_score} activity-eligible={len(recs)} "
            f"D={c['dual']} A={c['A_only']} B={c['B_only']} N={c['neither']}"
        )
    return "\n".join(lines) + "\n"


def apply_adjudication(rows: list[dict]) -> list[dict]:
    path = ROOT / "data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv"
    if not path.is_file():
        raise SystemExit("missing activity adjudication table; run scripts/analysis/adjudicate_activity_records.py first")
    by = {(r["pair"], r["ligand"]): r for r in read_csv(path)}
    out = []
    for rec in rows:
        adj = by.get((rec["pair"], rec["ligand_id"]))
        if adj is None or rec.get("analysis_set") != "main":
            out.append(rec)
            continue
        rec = dict(rec)
        rec["historical_pA"] = adj.get("historical_pA", rec["pA"])
        rec["historical_pB"] = adj.get("historical_pB", rec["pB"])
        rec["n_act_A"] = adj.get("n_act_A", "")
        rec["n_act_B"] = adj.get("n_act_B", "")
        rec["activity_eligible"] = int(adj.get("activity_eligible") or 0)
        rec["activity_status"] = adj.get("activity_status") or rec.get("activity_status") or ""
        pa, pb = fnum(adj.get("max_A")), fnum(adj.get("max_B"))
        rec["pA"] = "" if pa is None else pa
        rec["pB"] = "" if pb is None else pb
        rec["primary_class_theta6"] = assign_fourclass(pa, pb) or ""
        out.append(rec)
    return out


def main() -> int:
    rows = []
    rows.extend(load_egfr())
    rows.extend(load_ache())
    rows.extend(load_pm48())
    rows.extend(load_track_b())
    rows.extend(load_holdout_ab_pm())
    rows.extend(load_holdout_track_b())
    rows = apply_adjudication(rows)
    # uniqueness
    seen = set()
    for r in rows:
        key = (r["pair"], r["ligand_id"], r["analysis_set"])
        if key in seen:
            raise SystemExit(f"duplicate master key {key}")
        seen.add(key)
    rows.sort(key=lambda r: (PRIMARY_PAIRS.index(r["pair"]) if r["pair"] in PRIMARY_PAIRS else 99, r["analysis_set"], r["ligand_id"]))
    dests = [
        ROOT / "data/processed/current_score_master.csv",
        ROOT / "results/canonical/current_score_master.csv",
    ]
    for dest in dests:
        write_csv(dest, rows)
        print("wrote", dest.relative_to(ROOT), len(rows))
    report = compare_membership(rows)
    rep = ROOT / "results/canonical/score_master_migration_diff.md"
    rep.parent.mkdir(parents=True, exist_ok=True)
    rep.write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
