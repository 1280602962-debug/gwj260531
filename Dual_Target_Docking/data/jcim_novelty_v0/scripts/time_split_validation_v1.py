#!/usr/bin/env python3
"""Frozen time-split protocol and scored-ligand evaluation.

Cutoffs 2015, 2018, and 2020 were specified before AUROC was computed.
A ligand's year is the earliest retained high-confidence document year.
Late-document ligands never enter threshold, receptor, or metric selection.
If too few pairs meet the sample gate, the result is not packaged as
external validation.
"""
from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
CACHE = ROOT / "data" / "jcim_novelty_v0" / "cache" / "document_year_v1"
TAB.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

BASE = "https://www.ebi.ac.uk/chembl/api/data"
SEED = 20260729
N_BOOT = 2000
PRIMARY_CUTOFF = 2018
CUTOFFS = (2015, 2018, 2020)
MIN_CLASS_AUROC = 10
MIN_CLASS_PREFERRED = 15
MIN_EVALUABLE_PAIRS = 2
BATCH = 20

TARGETS = {
    "EGFR/HER2": ("CHEMBL203", "CHEMBL1824"),
    "AChE/BChE": ("CHEMBL220", "CHEMBL1914"),
    "PIK3CA/PIK3CB": ("CHEMBL4005", "CHEMBL3145"),
    "PIK3CA/mTOR": ("CHEMBL4005", "CHEMBL2842"),
}

CONTRASTS = (
    ("D_vs_A", "dual", "A_only", "vina_B"),
    ("D_vs_B", "dual", "B_only", "vina_A"),
)


def stable_offset(*parts, modulus=99991):
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value):
    try:
        if value in ("", None):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def auroc(pos, neg) -> float:
    if not pos or not neg:
        return float("nan")
    p = np.asarray(pos, float)
    n = np.asarray(neg, float)
    delta = p[:, None] - n[None, :]
    return float(((delta > 0).sum() + 0.5 * (delta == 0).sum()) / (len(p) * len(n)))


def request_json(endpoint: str, params: dict, cache_key: str) -> dict:
    cache = CACHE / f"{cache_key}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    query = urllib.parse.urlencode(params)
    url = f"{BASE}/{endpoint}.json?{query}"
    last_error = None
    for attempt in range(6):
        try:
            req = urllib.request.Request(
                url,
                headers={"Accept": "application/json", "User-Agent": "DualFourClass-Bench/1.0"},
            )
            with urllib.request.urlopen(req, timeout=90) as response:
                result = json.loads(response.read().decode("utf-8"))
            cache.write_text(json.dumps(result), encoding="utf-8")
            return result
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"ChEMBL request failed: {url}: {last_error}")


def fetch_years(document_ids: list[str]) -> dict[str, int | None]:
    lookup_path = TAB / "document_year_lookup_v1.csv"
    known = {}
    if lookup_path.exists():
        for row in read_csv(lookup_path):
            year = fnum(row.get("year"))
            known[row["document_chembl_id"]] = int(year) if year is not None else None
    missing = [doc for doc in document_ids if doc not in known]
    for start in range(0, len(missing), BATCH):
        batch = missing[start : start + BATCH]
        payload = request_json(
            "document",
            {"document_chembl_id__in": ",".join(batch), "limit": max(len(batch), 20)},
            "docs_" + hashlib.sha256(",".join(batch).encode()).hexdigest()[:16],
        )
        hits = {row["document_chembl_id"]: row for row in payload.get("documents") or []}
        for doc in batch:
            year = hits.get(doc, {}).get("year")
            try:
                known[doc] = int(year) if year not in (None, "") else None
            except (TypeError, ValueError):
                known[doc] = None
    rows = [
        {
            "document_chembl_id": doc,
            "year": "" if year is None else year,
            "source": "ChEMBL document.year API; cached lookup committed for reproducibility",
        }
        for doc, year in sorted(known.items())
    ]
    write_csv(lookup_path, rows)
    return known


def boot_ci(pos, neg, seed):
    if len(pos) < 2 or len(neg) < 2:
        return None, None
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(N_BOOT):
        p = [pos[i] for i in rng.integers(0, len(pos), len(pos))]
        n = [neg[i] for i in rng.integers(0, len(neg), len(neg))]
        value = auroc(p, n)
        if value == value:
            values.append(value)
    if len(values) < N_BOOT // 2:
        return None, None
    lo, hi = np.percentile(values, [2.5, 97.5])
    return float(lo), float(hi)


def gate(counts: dict[str, int]) -> str:
    dual, a_only, b_only = counts["dual"], counts["A_only"], counts["B_only"]
    if min(dual, a_only, b_only) >= MIN_CLASS_PREFERRED:
        return "evaluable"
    if min(dual, a_only, b_only) >= MIN_CLASS_AUROC:
        return "underpowered_report"
    if min(dual, a_only, b_only) > 0:
        return "descriptive_only"
    return "unevaluable"


def main() -> None:
    labels = read_csv(TAB / "high_confidence_labels_v1.csv")
    activities = [
        row
        for row in read_csv(TAB / "high_confidence_activity_audit_v1.csv")
        if row.get("keep") in ("1", "True", "true")
    ]
    by_mol_target = defaultdict(list)
    for row in activities:
        by_mol_target[(row["molecule_chembl_id"], row["target_chembl_id"])].append(row)

    document_ids = sorted({row["document_chembl_id"] for row in activities if row.get("document_chembl_id")})
    years = fetch_years(document_ids)

    ligand_rows = []
    recs_by_pair = defaultdict(list)
    for row in labels:
        pair = row["pair"]
        target_a, target_b = TARGETS[pair]
        mol = row["molecule_chembl_id"]
        docs = []
        for target in (target_a, target_b):
            for rec in by_mol_target[(mol, target)]:
                doc = rec.get("document_chembl_id")
                if doc:
                    docs.append(doc)
        doc_years = [years.get(doc) for doc in docs if years.get(doc) is not None]
        first_year = min(doc_years) if doc_years else None
        rec = {
            "pair": pair,
            "ligand": row["ligand"],
            "molecule_chembl_id": mol,
            "cls": row["frozen_class"],
            "first_year": first_year,
            "n_documents": len(set(docs)),
            "vina_A": fnum(row["vina_A"]),
            "vina_B": fnum(row["vina_B"]),
        }
        recs_by_pair[pair].append(rec)
        ligand_rows.append(
            {
                "pair": pair,
                "ligand": rec["ligand"],
                "molecule_chembl_id": mol,
                "frozen_class": rec["cls"],
                "first_document_year": "" if first_year is None else first_year,
                "n_dated_documents": len(doc_years),
                "n_documents": rec["n_documents"],
                "note": "year = min retained high-confidence document.year; missing year excluded from dated splits",
            }
        )
    write_csv(TAB / "time_split_ligand_years_v1.csv", ligand_rows)

    count_rows = []
    auroc_rows = []
    for cutoff in CUTOFFS:
        for pair, recs in recs_by_pair.items():
            dated = [rec for rec in recs if rec["first_year"] is not None]
            train = [rec for rec in dated if rec["first_year"] < cutoff]
            test = [rec for rec in dated if rec["first_year"] >= cutoff]
            for split_name, split in (("train_before", train), ("test_on_or_after", test)):
                counts = {
                    cls: sum(rec["cls"] == cls for rec in split)
                    for cls in ("dual", "A_only", "B_only", "neither")
                }
                count_rows.append(
                    {
                        "cutoff_year": cutoff,
                        "primary_cutoff": int(cutoff == PRIMARY_CUTOFF),
                        "pair": pair,
                        "split": split_name,
                        "n_dual": counts["dual"],
                        "n_A_only": counts["A_only"],
                        "n_B_only": counts["B_only"],
                        "n_neither": counts["neither"],
                        "n_total": len(split),
                        "gate": gate(counts) if split_name.startswith("test") else "train",
                        "note": "cutoffs frozen before AUROC; late ligands unused for method selection",
                    }
                )
            test_gate = gate(
                {
                    cls: sum(rec["cls"] == cls for rec in test)
                    for cls in ("dual", "A_only", "B_only", "neither")
                }
            )
            for contrast, pos_cls, neg_cls, score_key in CONTRASTS:
                pos = [rec[score_key] for rec in test if rec["cls"] == pos_cls and rec[score_key] is not None]
                neg = [rec[score_key] for rec in test if rec["cls"] == neg_cls and rec[score_key] is not None]
                report_auroc = test_gate in {"evaluable", "underpowered_report"}
                value = auroc(pos, neg) if report_auroc else float("nan")
                lo = hi = None
                if report_auroc:
                    lo, hi = boot_ci(pos, neg, SEED + stable_offset(pair, cutoff, contrast))
                auroc_rows.append(
                    {
                        "cutoff_year": cutoff,
                        "primary_cutoff": int(cutoff == PRIMARY_CUTOFF),
                        "pair": pair,
                        "contrast": contrast,
                        "n_pos_test": len(pos),
                        "n_neg_test": len(neg),
                        "gate": test_gate,
                        "auroc": round(value, 4) if value == value else "",
                        "ci_lo": round(lo, 4) if lo is not None else "",
                        "ci_hi": round(hi, 4) if hi is not None else "",
                        "packaged_as_external_validation": 0,
                        "note": (
                            "AUROC omitted because a directional class has n<10"
                            if not report_auroc
                            else "scored-panel time split; not a new docking campaign"
                        ),
                    }
                )

    primary_test = [
        row
        for row in count_rows
        if row["cutoff_year"] == PRIMARY_CUTOFF and row["split"] == "test_on_or_after"
    ]
    n_evaluable_pairs = len(
        {
            row["pair"]
            for row in primary_test
            if row["gate"] in {"evaluable", "underpowered_report"}
        }
    )
    package = int(n_evaluable_pairs >= MIN_EVALUABLE_PAIRS)
    for row in auroc_rows:
        if row["cutoff_year"] == PRIMARY_CUTOFF and row["gate"] in {
            "evaluable",
            "underpowered_report",
        }:
            row["packaged_as_external_validation"] = package

    write_csv(TAB / "time_split_class_counts_v1.csv", count_rows)
    write_csv(TAB / "time_split_auroc_v1.csv", auroc_rows)

    lines = [
        "# Time-split protocol freeze and result",
        "",
        "Frozen before seeing AUROC:",
        "",
        f"- Primary cutoff year: **{PRIMARY_CUTOFF}** (train first_year < {PRIMARY_CUTOFF}; test first_year ≥ {PRIMARY_CUTOFF}).",
        "- Sensitivity cutoffs: 2015 and 2020. These were not chosen after looking at class counts that favor AUROC.",
        "- Ligand year = minimum `document.year` among retained high-confidence records.",
        "- Late ligands are not used to choose thresholds, receptors, or endpoints.",
        "- Minimum AUROC gate: dual, A-only, and B-only each ≥10 in the test split; ≥15 preferred.",
        "- Below 10: descriptive counts only. Cutoff is not moved to recover a class.",
        f"- External-validation package requires ≥{MIN_EVALUABLE_PAIRS} pairs passing the AUROC gate at the primary cutoff.",
        "",
        f"Primary cutoff evaluable/underpowered pairs: **{n_evaluable_pairs}**.",
        (
            "Packaged as external validation: **yes**."
            if package
            else "Packaged as external validation: **no**. Keep the internal formulation-audit claim."
        ),
        "",
        "| cutoff | pair | test dual/A/B/neither | gate | D_vs_A | D_vs_B |",
        "|-------:|------|----------------------:|------|--------|--------|",
    ]
    auroc_index = {(r["cutoff_year"], r["pair"], r["contrast"]): r for r in auroc_rows}
    for row in count_rows:
        if row["split"] != "test_on_or_after":
            continue
        da = auroc_index[(row["cutoff_year"], row["pair"], "D_vs_A")]
        db = auroc_index[(row["cutoff_year"], row["pair"], "D_vs_B")]
        lines.append(
            f"| {row['cutoff_year']} | {row['pair']} | "
            f"{row['n_dual']}/{row['n_A_only']}/{row['n_B_only']}/{row['n_neither']} | "
            f"{row['gate']} | {da['auroc']} | {db['auroc']} |"
        )
    lines += [
        "",
        "If this split is insufficient, BindingDB is the next option (`docs/BINDINGDB_EXTERNAL_SOP.md`).",
        "Do not dock new ligands until the independent set and evaluation contract are frozen.",
        "",
    ]
    (ANALYSIS / "TIME_SPLIT_VERDICT.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
