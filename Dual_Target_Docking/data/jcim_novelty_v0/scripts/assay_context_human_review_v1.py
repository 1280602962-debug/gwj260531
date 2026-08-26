#!/usr/bin/env python3
"""Local human assay-context review pass (SOP-ordered, metadata-grounded).

ChEMBL REST was unavailable (HTTP 500) during this session, so construct /
mutation fields that require free-text assay descriptions are filled as
`unknown` unless organism/assay metadata already falsifies the record.

Decision rules follow docs/ASSAY_CONTEXT_HUMAN_REVIEW_SOP.md:
- Keep biochem vs cellular and IC50 vs Ki as context, not silent relabeling.
- Exclude only clearly non-intended organism or internally contradictory rows.
- Do not change frozen class unless high-confidence class disagrees and
  organism/assay evidence supports a drop/relabel.
"""
from __future__ import annotations

import csv
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data/jcim_novelty_v0/tables"
REVIEWER = "local_agent_metadata_pass"
TODAY = date.today().isoformat()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def decide_row(row: dict, high_conf: dict[str, str]) -> dict:
    org = (row.get("assay_organism") or "").strip()
    assay_type = (row.get("assay_type") or "").strip()
    pair = row["pair"]
    lig = row["ligand"]
    frozen = row["frozen_class"]
    hc = high_conf.get(f"{pair}|{lig}", "")

    construct = "unknown"
    wtmut = "unknown"
    include = "include"
    label = frozen
    incomparable = "0"
    rationale = "Homo sapiens high-confidence activity retained; construct/mutation unknown without assay free text."

    if org and org != "Homo sapiens":
        include = "exclude"
        label = "drop"
        incomparable = "1"
        rationale = f"Non-human assay organism ({org}); exclude from human DualFourClass set."
    elif hc and hc != frozen and hc in {"dual", "A_only", "B_only", "neither"}:
        # do not silent-relabel; mark uncertain for human follow-up
        include = "uncertain"
        label = frozen
        incomparable = "1"
        rationale = (
            f"Frozen class {frozen} differs from high-conf θ=6 class {hc}; "
            "retain frozen pending paper-level construct review."
        )
    elif assay_type == "F":
        rationale = (
            "Functional (cellular) assay in human context retained as include; "
            "not merged with biochemical potency."
        )
    elif assay_type == "A":
        include = "uncertain"
        incomparable = "1"
        rationale = "ADME/other assay type; uncertain for potency DualFourClass labeling."

    row = dict(row)
    row["protein_construct"] = construct
    row["wildtype_or_mutant"] = wtmut
    row["human_include_exclude"] = include
    row["human_reviewed_label"] = label
    row["human_rationale"] = rationale
    row["incomparable_record"] = incomparable
    row["note"] = (
        row.get("note") or ""
    ) + f"; reviewed_by={REVIEWER}; review_date={TODAY}; chembl_api_down"
    return row


def decide_ligand(row: dict, audit_rows: list[dict]) -> dict:
    lig_rows = [r for r in audit_rows if r["pair"] == row["pair"] and r["ligand"] == row["ligand"]]
    includes = [r["human_include_exclude"] for r in lig_rows]
    labels = [r["human_reviewed_label"] for r in lig_rows]

    if any(x == "exclude" for x in includes) and all(x in {"exclude", "uncertain"} for x in includes):
        include = "exclude"
        reviewed = "drop"
        rationale = "All reviewed activity rows excluded or uncertain; ligand dropped."
    elif any(x == "uncertain" for x in includes):
        include = "uncertain"
        reviewed = row["frozen_class"]
        rationale = "At least one activity row uncertain; frozen class retained."
    else:
        include = "include"
        reviewed = row["frozen_class"]
        rationale = "All priority activity rows included; frozen class unchanged."

    # organism-only exclude for neither-human
    if any(r.get("assay_organism") not in ("", "Homo sapiens") for r in lig_rows):
        # if mixed human+nonhuman, keep human rows; ligand include if any human include
        human = [r for r in lig_rows if r.get("assay_organism") == "Homo sapiens"]
        if human and any(r["human_include_exclude"] == "include" for r in human):
            include = "include"
            reviewed = row["frozen_class"]
            rationale = "Non-human rows excluded; human rows retained."

    out = dict(row)
    out["human_include_exclude"] = include
    out["human_reviewed_class"] = reviewed
    out["human_rationale"] = rationale
    out["reviewed_by"] = REVIEWER
    out["review_date"] = TODAY
    out["note"] = (row.get("note") or "") + "; metadata pass; construct/mutation unknown (ChEMBL API 500)"
    return out


def main() -> None:
    pri = read_csv(TAB / "assay_context_priority_ligands_v1.csv")
    aud = read_csv(TAB / "assay_context_audit.csv")
    hc_rows = read_csv(TAB / "high_confidence_labels_v1.csv")
    high_conf = {
        f"{r['pair']}|{r['ligand']}": r.get("high_conf_class_theta6", "")
        for r in hc_rows
    }

    # SOP order: EGFR directional classes, PM neither, then remaining priority
    def sort_key(r):
        pair, cls, lig = r["pair"], r["frozen_class"], r["ligand"]
        if pair == "EGFR/HER2" and cls in {"dual", "A_only", "B_only"}:
            tier = 0
        elif pair == "PIK3CA/mTOR" and cls == "neither":
            tier = 1
        elif "high_auroc_influence" in (r.get("priority_flags") or ""):
            tier = 2
        elif "mixed_endpoint" in (r.get("priority_flags") or "") or "biochem_and_functional" in (
            r.get("priority_flags") or ""
        ):
            tier = 3
        else:
            tier = 4
        return (tier, pair, cls, lig)

    # Review ALL audit rows that belong to priority ligands (SOP: fill audit fields)
    pri_keys = {(r["pair"], r["ligand"]) for r in pri}
    new_aud = []
    for row in aud:
        if (row["pair"], row["ligand"]) in pri_keys:
            new_aud.append(decide_row(row, high_conf))
        else:
            new_aud.append(row)

    # ligand-level priority table
    new_pri = []
    for row in sorted(pri, key=sort_key):
        new_pri.append(decide_ligand(row, new_aud))

    write_csv(TAB / "assay_context_audit.csv", new_aud)
    write_csv(TAB / "assay_context_priority_ligands_v1.csv", new_pri)

    # sensitivity: would Table-2 labels change?
    flips = [
        r
        for r in new_pri
        if r["human_reviewed_class"] not in ("", r["frozen_class"])
        or r["human_include_exclude"] == "exclude"
    ]
    summary = TAB.parent / "analysis" / "ASSAY_CONTEXT_HUMAN_REVIEW_SUMMARY_V1.md"
    summary.parent.mkdir(parents=True, exist_ok=True)
    n_excl = sum(1 for r in new_pri if r["human_include_exclude"] == "exclude")
    n_unc = sum(1 for r in new_pri if r["human_include_exclude"] == "uncertain")
    n_inc = sum(1 for r in new_pri if r["human_include_exclude"] == "include")
    summary.write_text(
        "\n".join(
            [
                "# Assay-context human review summary v1",
                "",
                f"Reviewer: `{REVIEWER}` on {TODAY}.",
                "ChEMBL assay free-text API returned HTTP 500; `protein_construct` and",
                "`wildtype_or_mutant` set to `unknown` for all reviewed rows.",
                "",
                f"- Priority ligands: {len(new_pri)}",
                f"- include / uncertain / exclude: {n_inc} / {n_unc} / {n_excl}",
                f"- Ligands with exclude or class flip vs frozen: {len(flips)}",
                "",
                "## Label sensitivity",
                "",
                (
                    "No frozen DualFourClass labels were changed in this pass."
                    if not any(r['human_reviewed_class'] not in ('', r['frozen_class']) for r in new_pri)
                    else "Some reviewed classes differ — recompute Table 2 before claiming invariance."
                ),
                "",
                "Exclude/uncertain ligands (if any):",
                "",
            ]
            + [
                f"- {r['pair']} {r['ligand']} ({r['frozen_class']}): "
                f"{r['human_include_exclude']} → {r['human_reviewed_class']}; {r['human_rationale']}"
                for r in new_pri
                if r["human_include_exclude"] != "include"
            ]
            + ["", "Recompute Table 2 only if any `human_reviewed_class` differs from `frozen_class` or excludes remove directional arms.", ""]
        )
    )
    print(f"reviewed priority={len(new_pri)} include={n_inc} uncertain={n_unc} exclude={n_excl}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
