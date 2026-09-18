#!/usr/bin/env python3
"""Read current files and check the eight-pair scientific chain.

Not a release-engineering or checksum framework. Exit 0 only if the
current score master, canonical Table 2 points, and manuscripts agree.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
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
FORBIDDEN = (
    "submission_pack_pre_v4_archive",
    "data/_legacy_archive",
    "pik3ca_pik3cb",
    "mcl1_bclxl",
    "remediation_outputs",
)
PREFIX = ("0.430", "0.808", "0.378")


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


def main() -> int:
    master = read_csv(ROOT / "results/canonical/current_score_master.csv")
    pairs = {r["pair"] for r in master}
    if pairs != set(PAIRS):
        fail(f"master pairs {sorted(pairs)}")
    counts = Counter((r["pair"], r["ligand_id"]) for r in master)
    dups = [k for k, n in counts.items() if n > 1]
    if dups:
        fail(f"duplicate pair/ligand keys {dups[:5]}")

    ab = {
        r["ligand"]: r
        for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")
    }
    egfr = [r for r in master if r["pair"] == "EGFR/HER2"]
    if len(egfr) != 110:
        fail(f"EGFR n={len(egfr)}")
    for row in egfr:
        old = ab.get(row["ligand_id"])
        if old is None:
            fail(f"EGFR {row['ligand_id']} missing from ablation scores")
        sa = -float(old["3POZ_affinity"])
        sb = -float(old["3RCD_affinity"])
        if abs(float(row["score_A"]) - sa) > 1e-6 or abs(float(row["score_B"]) - sb) > 1e-6:
            fail(f"EGFR {row['ligand_id']} score mismatch")

    ache = [
        r
        for r in master
        if r["pair"] == "AChE/BChE"
        and r.get("analysis_set") == "main"
        and r.get("complete_case") == "1"
        and str(r.get("activity_eligible", "1")) in ("1", "True")
    ]
    n = Counter(r["primary_class_theta6"] for r in ache)
    if (n["dual"], n["A_only"], n["B_only"]) != (27, 26, 28):
        fail(f"AChE n_scored {dict(n)}")

    direc = {(r["pair"], r["estimand"]): r for r in read_csv(ROOT / "results/canonical/primary_directional_auroc.csv")}
    smin = {r["pair"]: r for r in read_csv(ROOT / "results/canonical/primary_summary_min.csv")}
    packs: dict[str, dict[str, list]] = {}
    for row in master:
        if row.get("analysis_set") != "main" or row.get("complete_case") != "1":
            continue
        if str(row.get("activity_eligible", "1")) not in ("1", "True"):
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

    en = (ROOT / "docs/MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    zh = (ROOT / "docs/MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    for token in PREFIX:
        if token in en or token in zh:
            fail(f"pre-fix token {token} still in manuscript")
    for pair, rec in smin.items():
        token = r3(rec["summary_min"])
        if token not in en or token not in zh:
            fail(f"{pair} summary_min {token} missing from manuscripts")
        line_en = next((ln for ln in en.splitlines() if ln.startswith(f"| {pair} |") and "[" in ln and ln.count("|") >= 5), "")
        if not line_en:
            fail(f"EN Table 2 missing {pair}")

    inc = read_csv(ROOT / "results/canonical/ecfp4_incremental_information.csv")
    max_abs = max(abs(float(r["delta_ECFP4_plus_docking_minus_ECFP4"])) for r in inc)
    inc_token = r3(max_abs)
    if inc_token not in en or inc_token not in zh:
        fail(f"ECFP incremental |Δ| {inc_token} missing from manuscripts")
    if "最大绝对变化为 0.023" in zh or "at most 0.023 across" in en:
        if inc_token != "0.023":
            fail("stale ECFP incremental 0.023 remains in manuscript")

    pack_en = (ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    if r3(smin["EGFR/HER2"]["summary_min"]) not in pack_en:
        fail("submission_pack EN manuscript missing EGFR summary_min")

    det = read_csv(ROOT / "results/canonical/detectable_effect_simulation.csv")
    det_pairs = {r["pair"] for r in det}
    if det_pairs != set(PAIRS):
        fail(f"detectable-effect pairs {sorted(det_pairs)}")
    if any(r.get("n_mc") not in {"1000", 1000} for r in det):
        fail("detectable-effect n_mc is not 1000")
    ache = next(
        r
        for r in det
        if r["pair"] == "AChE/BChE" and r["contrast"] == "summary_min" and r["true_auroc"] == "0.50"
    )
    if ache["n_neg"] != "26/28":
        fail(f"AChE detectable-effect n_A/n_B={ache['n_neg']} (expected 26/28)")
    if any(r.get("bootstrap") != "class_stratified_shared_dual" for r in det):
        fail("detectable-effect bootstrap is not class_stratified_shared_dual")

    scan = list((ROOT / "scripts/analysis").glob("*.py"))
    scan += [
        ROOT / "docs/MANUSCRIPT_JCIM_EN.md",
        ROOT / "docs/MANUSCRIPT_JCIM_ZH.md",
        ROOT / "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
        ROOT / "docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
        ROOT / "figures/jcim_article/scripts/update_figures_pr32.py",
        ROOT / "figures/jcim_article/scripts/inject_s2_boxes_from_json.py",
    ]
    for path in scan:
        text = path.read_text(encoding="utf-8")
        for token in FORBIDDEN:
            if token in text and not (path.name == "build_current_score_master.py" and token == "data/_legacy_archive"):
                fail(f"{path.relative_to(ROOT)} still mentions {token}")

    print("PASS: current score master, Table 2 replay, manuscripts, and current analysis paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
