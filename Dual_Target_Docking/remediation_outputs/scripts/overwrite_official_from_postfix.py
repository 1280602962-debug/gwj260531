#!/usr/bin/env python3
"""Adopt post-fix EGFR/AChE docking tables as official and patch manuscripts.

EGFR ligand-only ECFP4/physchem rows stay at the pre-fix values.
Figure scripts must not keep hardcoded pre-fix AUROCs.
"""
from __future__ import annotations

import csv
import shutil
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
STG = ROOT / "remediation_outputs/canonical_tables/jcim_novelty_v0"
OFF = ROOT / "data/jcim_novelty_v0/tables"
SUB = ROOT / "submission_pack/tables"
SNAP = ROOT / "figures/jcim_article/input_snapshot/data/jcim_novelty_v0/tables"


def replace_pair_rows(dest: Path, src: Path, pairs: set[str], pair_key: str = "pair") -> None:
    if not dest.exists() or not src.exists():
        print("skip missing", dest, src)
        return
    src_rows = [r for r in csv.DictReader(src.open()) if r[pair_key] in pairs]
    dest_rows = list(csv.DictReader(dest.open()))
    fields = dest_rows[0].keys()
    kept = [r for r in dest_rows if r[pair_key] not in pairs]
    # keep dest column order; fill from src
    new_pair_rows = []
    for r in src_rows:
        new_pair_rows.append({k: r.get(k, "") for k in fields})
    out = []
    inserted = set()
    for r in dest_rows:
        if r[pair_key] in pairs:
            if r[pair_key] not in inserted:
                out.extend([x for x in new_pair_rows if x[pair_key] == r[pair_key]])
                inserted.add(r[pair_key])
        else:
            out.append(r)
    for r in new_pair_rows:
        if r[pair_key] not in inserted:
            out.append(r)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fields))
        w.writeheader()
        w.writerows(out)
    print("updated", dest)


def patch_incremental_official() -> None:
    dest = OFF / "incremental_information_v1.csv"
    src = STG / "incremental_information_v1.csv"
    src_rows = list(csv.DictReader(src.open()))
    dest_rows = list(csv.DictReader(dest.open()))
    src_by = {(r["pair"], r["contrast"], r["model"]): r for r in src_rows}
    ligand_only = {"ECFP4", "physchem", "ECFP4+physchem"}
    out = []
    for r in dest_rows:
        key = (r["pair"], r["contrast"], r["model"])
        s = src_by.get(key)
        if r["pair"] == "AChE/BChE" and s:
            out.append(s)
        elif r["pair"] == "EGFR/HER2" and s:
            # keep official ligand-only cv_auroc; update docking rank + docking-using models' rank column
            merged = dict(r)
            merged["rank_auroc_docking"] = s["rank_auroc_docking"]
            if r["model"] not in ligand_only:
                # docking-containing logistic used new scores; keep official ECFP4-only
                if "ECFP4" in r["model"]:
                    merged["cv_auroc"] = r["cv_auroc"]
                    merged["note"] = (
                        r["note"]
                        + "; ligand-only ECFP4 frozen pre-fix; rank_auroc_docking post-fix Vina"
                    )
                else:
                    merged["cv_auroc"] = s["cv_auroc"]
            out.append(merged)
        else:
            out.append(r)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dest_rows[0].keys()))
        w.writeheader()
        w.writerows(out)
    print("updated incremental with frozen EGFR ligand-only", dest)


def patch_unified() -> None:
    path = ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv"
    rows = list(csv.DictReader(path.open()))
    for r in rows:
        if r["pair"] == "EGFR/HER2" and r["label_rule"] == "theta_6.0":
            r["pocket_matched_summary_min"] = "0.3237"
            r["auroc_D_vs_A"] = "0.6607"
            r["auroc_D_vs_B"] = "0.3237"
            r["ci_lo"] = "0.1953"
            r["ci_hi"] = "0.4710"
        if r["pair"] == "AChE/BChE":
            r["n_A_only"] = "26"
            r["pocket_matched_summary_min"] = "0.6058"
            r["auroc_D_vs_A"] = "0.6524"
            r["auroc_D_vs_B"] = "0.6058"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("updated", path)


def copy_if(src: Path, dest: Path) -> None:
    if src.exists() and dest.parent.exists():
        shutil.copy2(src, dest)
        print("copied", dest)


def main() -> int:
    pairs = {"EGFR/HER2", "AChE/BChE"}
    for name in (
        "formulation_equal_score_negative_v1.csv",
        "formulation_conventional_vs_directional_v1.csv",
        "formulation_summary_v1.csv",
    ):
        replace_pair_rows(OFF / name, STG / name, pairs)
        if (SUB / name).exists():
            replace_pair_rows(SUB / name, STG / name, pairs)
        if (SNAP / name).exists():
            replace_pair_rows(SNAP / name, STG / name, pairs)
    patch_incremental_official()
    patch_unified()

    gnina_src = ROOT / "remediation_outputs/phase_gnina_independent/independent_dock_formulation_EGFR_HER2.csv"
    gnina_off = ROOT / "data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv"
    if gnina_src.exists() and gnina_off.exists():
        src = {r["contrast"]: r for r in csv.DictReader(gnina_src.open())}
        rows = list(csv.DictReader(gnina_off.open()))
        for r in rows:
            if r["pair"] != "EGFR/HER2":
                continue
            m = src.get(r["contrast"])
            if not m:
                continue
            r["auroc"] = f"{float(m['auroc']):.4f}"
            if m.get("ci_lo"):
                r["ci_lo"] = f"{float(m['ci_lo']):.4f}"
                r["ci_hi"] = f"{float(m['ci_hi']):.4f}"
            r["n_pos"] = m["n_pos"]
            r["n_neg"] = m["n_neg"]
        with gnina_off.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("updated", gnina_off)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
