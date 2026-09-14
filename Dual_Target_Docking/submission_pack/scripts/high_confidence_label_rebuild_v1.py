#!/usr/bin/env python3
"""Post-hoc high-confidence ChEMBL label robustness view.

The primary panels remain frozen. This script queries the current ChEMBL API
for panel compounds and retains exact quantitative records on Homo sapiens
SINGLE PROTEIN targets from assays with confidence_score >= 8. Results are a
dated robustness analysis, not a replacement for the frozen 2026-07-23 labels.
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
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from assay_aggregation_max_vs_median_v1 import assemble_jobs, auroc, boot_auroc, classify


ROOT = Path(__file__).resolve().parents[3]
BASE = "https://www.ebi.ac.uk/chembl/api/data"
CACHE = ROOT / "data" / "jcim_novelty_v0" / "cache" / "high_confidence_v1"
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
CACHE.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

ALLOWED_TYPES = {"IC50", "KI", "KD", "EC50", "POTENCY"}
BATCH = 40


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def request_json(endpoint: str, params: dict, cache_key: str) -> dict:
    cache = CACHE / f"{cache_key}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    query = urllib.parse.urlencode(params)
    url = f"{BASE}/{endpoint}.json?{query}"
    last_error = None
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "DualFourClass-Bench/1.0"})
            with urllib.request.urlopen(req, timeout=90) as response:
                result = json.loads(response.read().decode("utf-8"))
            cache.write_text(json.dumps(result), encoding="utf-8")
            return result
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"ChEMBL request failed: {url}: {last_error}")


def fetch_pages(endpoint: str, params: dict, item_key: str, prefix: str) -> list[dict]:
    rows = []
    offset = 0
    while True:
        page_params = {**params, "limit": 1000, "offset": offset}
        digest = hashlib.sha256(json.dumps(page_params, sort_keys=True).encode()).hexdigest()[:20]
        payload = request_json(endpoint, page_params, f"{prefix}_{digest}")
        items = payload.get(item_key) or []
        rows.extend(items)
        if not (payload.get("page_meta") or {}).get("next") or not items:
            break
        offset += 1000
    return rows


def chunks(values: list[str], size: int = BATCH):
    for start in range(0, len(values), size):
        yield values[start:start + size]


def main() -> None:
    ligands, jobs = assemble_jobs()
    target_to_molecules: dict[str, set[str]] = defaultdict(set)
    for molecule, target in jobs:
        target_to_molecules[target].add(molecule)

    target_meta = {}
    for target in sorted(target_to_molecules):
        payload = request_json("target", {"target_chembl_id": target, "limit": 10}, f"target_{target}")
        hits = payload.get("targets") or []
        if len(hits) != 1:
            raise RuntimeError(f"expected one target record for {target}, got {len(hits)}")
        target_meta[target] = hits[0]

    activities = []
    for target, molecules in sorted(target_to_molecules.items()):
        mols = sorted(molecules)
        for batch_index, batch in enumerate(chunks(mols), 1):
            activities.extend(
                fetch_pages(
                    "activity",
                    {"target_chembl_id": target, "molecule_chembl_id__in": ",".join(batch)},
                    "activities",
                    f"activity_{target}_{batch_index}",
                )
            )
        print(f"{target}: {len(mols)} panel molecules")

    assay_ids = sorted({a.get("assay_chembl_id") for a in activities if a.get("assay_chembl_id")})
    assays = []
    for batch_index, batch in enumerate(chunks(assay_ids), 1):
        assays.extend(
            fetch_pages(
                "assay",
                {"assay_chembl_id__in": ",".join(batch)},
                "assays",
                f"assay_{batch_index}",
            )
        )
    assay_meta = {a["assay_chembl_id"]: a for a in assays}

    retained: dict[tuple[str, str], list[dict]] = defaultdict(list)
    audit_rows = []
    for activity in activities:
        molecule = activity.get("molecule_chembl_id")
        target = activity.get("target_chembl_id")
        assay = assay_meta.get(activity.get("assay_chembl_id"), {})
        target_row = target_meta.get(target, {})
        reasons = []
        try:
            value = float(activity.get("pchembl_value"))
        except (TypeError, ValueError):
            value = None
            reasons.append("missing_pchembl")
        if target_row.get("target_type") != "SINGLE PROTEIN":
            reasons.append("target_not_single_protein")
        if target_row.get("organism") != "Homo sapiens":
            reasons.append("target_not_human")
        if int(assay.get("confidence_score") or 0) < 8:
            reasons.append("confidence_lt_8")
        if str(activity.get("standard_relation") or "").strip() != "=":
            reasons.append("relation_not_equal")
        if str(activity.get("standard_type") or "").upper() not in ALLOWED_TYPES:
            reasons.append("endpoint_not_allowed")
        if activity.get("data_validity_comment") not in (None, ""):
            reasons.append("data_validity_flag")
        if activity.get("potential_duplicate") in (True, 1, "1"):
            reasons.append("potential_duplicate")
        keep = not reasons
        compact = {
            "molecule_chembl_id": molecule,
            "target_chembl_id": target,
            "assay_chembl_id": activity.get("assay_chembl_id"),
            "document_chembl_id": activity.get("document_chembl_id"),
            "pchembl_value": value,
            "standard_type": activity.get("standard_type"),
            "standard_relation": activity.get("standard_relation"),
            "assay_type": assay.get("assay_type"),
            "assay_organism": assay.get("assay_organism"),
            "confidence_score": assay.get("confidence_score"),
            "keep": int(keep),
            "exclusion_reasons": ";".join(reasons),
        }
        audit_rows.append(compact)
        if keep:
            retained[(molecule, target)].append(compact)

    ligand_rows = []
    for ligand in ligands:
        values_a = [r["pchembl_value"] for r in retained[(ligand["molecule_chembl_id"], ligand["target_A"])]]
        values_b = [r["pchembl_value"] for r in retained[(ligand["molecule_chembl_id"], ligand["target_B"])]]
        pa = max(values_a) if values_a else None
        pb = max(values_b) if values_b else None
        cls = classify(pa, pb)
        ligand_rows.append(
            {
                "pair": ligand["pair"],
                "ligand": ligand["ligand"],
                "molecule_chembl_id": ligand["molecule_chembl_id"],
                "frozen_class": ligand["frozen_class"],
                "high_conf_class_theta6": cls,
                "high_conf_max_A": pa,
                "high_conf_max_B": pb,
                "n_high_conf_A": len(values_a),
                "n_high_conf_B": len(values_b),
                "complete_high_conf": int(pa is not None and pb is not None),
                "class_matches_frozen": int(cls == ligand["frozen_class"]) if cls else "",
                "vina_A": ligand["vina_A"],
                "vina_B": ligand["vina_B"],
            }
        )

    summary = []
    for pair in dict.fromkeys(r["pair"] for r in ligand_rows):
        all_pair = [r for r in ligand_rows if r["pair"] == pair]
        complete = [r for r in all_pair if r["complete_high_conf"]]
        D = [r for r in complete if r["high_conf_class_theta6"] == "dual"]
        A = [r for r in complete if r["high_conf_class_theta6"] == "A_only"]
        B = [r for r in complete if r["high_conf_class_theta6"] == "B_only"]
        N = [r for r in complete if r["high_conf_class_theta6"] == "neither"]
        da = boot_auroc([r["vina_B"] for r in D], [r["vina_B"] for r in A])
        db = boot_auroc([r["vina_A"] for r in D], [r["vina_A"] for r in B])
        summary.append(
            {
                "pair": pair,
                "n_frozen_scored": len(all_pair),
                "n_complete_high_conf": len(complete),
                "coverage_fraction": round(len(complete) / len(all_pair), 4),
                "n_class_matches_frozen": sum(r["class_matches_frozen"] == 1 for r in complete),
                "n_dual": len(D),
                "n_A_only": len(A),
                "n_B_only": len(B),
                "n_neither": len(N),
                "auroc_D_vs_A": None if da[0] != da[0] else round(da[0], 4),
                "auroc_D_vs_B": None if db[0] != db[0] else round(db[0], 4),
                "summary_min": None if da[0] != da[0] or db[0] != db[0] else round(min(da[0], db[0]), 4),
                "underpowered": int(min(len(D), len(A), len(B)) < 8),
            }
        )

    retrieved = datetime.now(timezone.utc).isoformat()
    write_csv(TAB / "high_confidence_activity_audit_v1.csv", audit_rows)
    write_csv(TAB / "high_confidence_labels_v1.csv", ligand_rows)
    write_csv(TAB / "high_confidence_summary_v1.csv", summary)
    run_metadata = {
        "retrieved_utc": retrieved,
        "filters": {
            "target_type": "SINGLE PROTEIN",
            "target_organism": "Homo sapiens",
            "confidence_score_gte": 8,
            "standard_relation": "=",
            "standard_types": sorted(ALLOWED_TYPES),
            "exclude_data_validity_flags": True,
            "exclude_potential_duplicates": True,
        },
        "post_hoc_robustness_not_primary": True,
        "n_activity_rows": len(audit_rows),
        "n_assays": len(assay_meta),
    }
    metadata_text = json.dumps(run_metadata, indent=2) + "\n"
    (CACHE / "RUN_METADATA.json").write_text(metadata_text, encoding="utf-8")
    (OUT / "analysis" / "high_confidence_run_meta_v1.json").write_text(
        metadata_text, encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
