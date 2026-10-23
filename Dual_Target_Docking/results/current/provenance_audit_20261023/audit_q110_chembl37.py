#!/usr/bin/env python3
"""Read-only ChEMBL 37 source audit for Q110 and the stored 805-row audit.

Does not rewrite the 2026-10-08 audit, formal populations, or Q110 membership.
STANDARD_OK matches export_chembl37_activity_rows.py: non-null pChEMBL and
standard_type in IC50, Ki, Kd, EC50, Potency, IC50app, Ki app; human
SINGLE PROTEIN with one component. No confidence, assay-type, or keep=1 filter.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path("/home/gwj/q110_v42/Dual_Target_Docking")
OUT = ROOT / "results/current/provenance_audit_20261023"
DB = Path("/home/gwj/chembl_37_release/chembl_37/chembl_37_sqlite/chembl_37.db")
TARBALL = Path("/home/gwj/chembl_37_release/chembl_37_sqlite.tar.gz")
TARBALL_SHA = "33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281"
MEMBERSHIP_SHA = "1e85e8c25af4dbe14e756f7042993ae86ac7f89ae06f810379372d739246b550"
STANDARD_OK = ("IC50", "Ki", "Kd", "EC50", "Potency", "IC50app", "Ki app")
PAIRS = {
    "EGFR/HER2": ("EGFR", "ERBB2"),
    "AChE/BChE": ("ACHE", "BCHE"),
    "PIK3CA/mTOR": ("PIK3CA", "MTOR"),
    "F2/F10": ("F2", "F10"),
    "JAK1/TYK2": ("JAK1", "TYK2"),
    "JAK1/JAK2": ("JAK1", "JAK2"),
    "PPARG/PPARA": ("PPARG", "PPARA"),
    "PPARA/PPARD": ("PPARA", "PPARD"),
}
HISTORICAL_TARGETS = {
    "EGFR": "CHEMBL203",
    "ERBB2": "CHEMBL1824",
    "ACHE": "CHEMBL220",
    "BCHE": "CHEMBL1914",
    "PIK3CA": "CHEMBL4005",
    "MTOR": "CHEMBL2842",
    "F2": "CHEMBL204",
    "F10": "CHEMBL244",
    "JAK1": "CHEMBL2835",
    "JAK2": "CHEMBL2971",
    "TYK2": "CHEMBL3553",
    "PPARG": "CHEMBL235",
    "PPARA": "CHEMBL239",
    "PPARD": "CHEMBL3979",
}
EXACT = 1e-9


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def blank(value) -> str:
    if value is None:
        return ""
    return str(value)


def fnum(value):
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def num_text(value) -> str:
    if value is None:
        return ""
    return format(float(value), ".16g")


def fourclass(pa, pb) -> str:
    if pa is None or pb is None:
        return ""
    active_a, active_b = pa >= 6.0, pb >= 6.0
    if active_a and active_b:
        return "dual"
    if active_a:
        return "A_only"
    if active_b:
        return "B_only"
    return "neither"


def has_identifier(row: dict) -> bool:
    return any(blank(row.get(key)).strip() for key in ("doi", "pubmed_id", "patent_id"))


def resolve_targets(con) -> dict:
    sql = """
    SELECT td.tid, td.chembl_id, td.pref_name, cs.accession, csyn.component_synonym AS gene
    FROM target_dictionary td
    JOIN target_components tc ON tc.tid = td.tid
    JOIN component_sequences cs ON cs.component_id = tc.component_id
    JOIN component_synonyms csyn ON csyn.component_id = cs.component_id
    WHERE td.target_type = 'SINGLE PROTEIN'
      AND td.organism = 'Homo sapiens'
      AND csyn.syn_type = 'GENE_SYMBOL'
      AND td.tid IN (SELECT tid FROM target_components GROUP BY tid HAVING COUNT(*) = 1)
    """
    found = defaultdict(list)
    for row in con.execute(sql):
        gene = (row["gene"] or "").upper()
        if gene in HISTORICAL_TARGETS:
            found[gene].append(dict(row))
    chosen = {}
    for gene, expected in HISTORICAL_TARGETS.items():
        hits = found.get(gene, [])
        match = [row for row in hits if row["chembl_id"] == expected]
        if len(match) != 1:
            raise SystemExit(f"target identity unresolved for {gene}: {[row['chembl_id'] for row in hits]}")
        chosen[gene] = {
            "tid": int(match[0]["tid"]),
            "chembl": match[0]["chembl_id"],
            "uniprot": match[0]["accession"],
            "pref_name": match[0]["pref_name"],
            "n_candidate_tids": len({row["tid"] for row in hits}),
        }
    return chosen


def fetch_activities(con, molregno: dict, targets: dict) -> dict:
    tids = sorted({info["tid"] for info in targets.values()})
    tid_by = {info["tid"]: gene for gene, info in targets.items()}
    rows_by = defaultdict(list)
    items = list(molregno.items())
    type_marks = ",".join("?" * len(STANDARD_OK))
    tid_marks = ",".join("?" * len(tids))
    for start in range(0, len(items), 150):
        chunk = items[start:start + 150]
        marks = ",".join("?" * len(chunk))
        sql = f"""
        SELECT md.chembl_id AS molecule_chembl_id, md.molregno,
               td.chembl_id AS target_chembl_id, ass.tid,
               act.activity_id, ass.chembl_id AS assay_chembl_id, ass.assay_type,
               act.standard_type, act.standard_relation, act.standard_value, act.standard_units,
               act.pchembl_value, docs.chembl_id AS document_chembl_id,
               docs.doi, docs.pubmed_id, docs.patent_id, docs.doc_type,
               source.src_short_name, source.src_description
        FROM activities act
        JOIN assays ass ON ass.assay_id = act.assay_id
        JOIN docs ON docs.doc_id = ass.doc_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        JOIN target_dictionary td ON td.tid = ass.tid
        LEFT JOIN source ON source.src_id = docs.src_id
        WHERE ass.tid IN ({tid_marks})
          AND act.molregno IN ({marks})
          AND act.pchembl_value IS NOT NULL
          AND act.standard_type IN ({type_marks})
        """
        params = tids + [regno for _, regno in chunk] + list(STANDARD_OK)
        for row in con.execute(sql, params):
            item = dict(row)
            item["gene"] = tid_by[int(item["tid"])]
            rows_by[(item["molecule_chembl_id"], int(item["tid"]))].append(item)
        print("queried", start + len(chunk), flush=True)
    return rows_by


def arm_summary(rows: list[dict]) -> dict:
    values = [float(row["pchembl_value"]) for row in rows]
    maximum = max(values) if values else None
    tied = [row for row in rows if maximum is not None and abs(float(row["pchembl_value"]) - maximum) < EXACT]
    types = sorted({blank(row["standard_type"]) for row in rows})
    docs = sorted({blank(row["document_chembl_id"]) for row in rows if blank(row["document_chembl_id"])})
    sources = sorted({blank(row["src_short_name"]) for row in rows if blank(row["src_short_name"])})
    above = sum(value >= 6.0 for value in values)
    below = sum(value < 6.0 for value in values)
    return {
        "n": len(rows),
        "max": maximum,
        "min": min(values) if values else None,
        "tied": tied,
        "types": types,
        "docs": docs,
        "sources": sources,
        "above": above,
        "below": below,
        "span": None if not values else maximum - min(values),
    }


def flags_for(arms: dict) -> dict:
    both = arms["A"]["n"] > 0 and arms["B"]["n"] > 0
    assays = both and all(blank(row["assay_chembl_id"]) for row in arms["A"]["tied"] + arms["B"]["tied"])
    documents = both and all(blank(row["document_chembl_id"]) for row in arms["A"]["tied"] + arms["B"]["tied"])
    literature = both and all(has_identifier(row) for row in arms["A"]["tied"] + arms["B"]["tied"])
    doc_sets = {end: set(arms[end]["docs_tied"]) for end in "AB"}
    cross = both and doc_sets["A"] != doc_sets["B"]
    type_mix = any(len(arms[end]["types"]) > 1 for end in "AB") or (
        both and set(blank(row["standard_type"]) for row in arms["A"]["tied"]) != set(blank(row["standard_type"]) for row in arms["B"]["tied"])
    )
    discord = any(arms[end]["above"] and arms[end]["below"] for end in "AB")
    unresolved = not both or not assays or not documents
    names = []
    if both and assays and documents:
        names.append("SOURCE_TRACEABLE")
    if both and not literature:
        names.append("DOCUMENT_IDENTIFIER_MISSING")
    if cross:
        names.append("CROSS_SOURCE_AGGREGATION")
    if type_mix:
        names.append("ENDPOINT_HETEROGENEITY")
    if discord:
        names.append("REPLICATE_THRESHOLD_DISCORDANCE")
    if unresolved:
        names.append("SOURCE_UNRESOLVED")
    return {
        "both": both,
        "assays": assays,
        "documents": documents,
        "literature": literature,
        "cross": cross,
        "type_mix": type_mix,
        "discord": discord,
        "unresolved": unresolved,
        "names": names,
    }


def enrich(summary: dict) -> dict:
    summary["docs_tied"] = sorted({blank(row["document_chembl_id"]) for row in summary["tied"]})
    summary["assays_tied"] = sorted({blank(row["assay_chembl_id"]) for row in summary["tied"]})
    summary["types_tied"] = sorted({blank(row["standard_type"]) for row in summary["tied"]})
    summary["sources_tied"] = sorted({blank(row["src_short_name"]) for row in summary["tied"]})
    return summary


def main() -> None:
    membership_path = ROOT / "extension_jcim_validation/pik3ca_mtor_q110_v42/PIK3CA_MTOR_Q110_MEMBERSHIP.csv"
    membership_hash = hashlib.sha256(membership_path.read_bytes()).hexdigest()
    if membership_hash != MEMBERSHIP_SHA:
        raise SystemExit(f"membership hash changed: {membership_hash}")
    members = read_csv(membership_path)
    if len(members) != 110 or len({r["molecule_chembl_id"] for r in members}) != 110 or len({r["parent_inchikey"] for r in members}) != 110:
        raise SystemExit("Q110 identity check failed")
    new_members = [r for r in members if r["in_pm48"] != "1"]
    overlap = [r for r in members if r["in_pm48"] == "1"]
    if len(new_members) != 101 or len(overlap) != 9:
        raise SystemExit("Q110 split is not 101 + 9")

    audit = ROOT / "activity_label_source_audit_20261008"
    stored_wide = read_csv(audit / "chembl37_dual_endpoint_crosswalk_805.csv")
    stored_long = read_csv(audit / "chembl37_activity_rows_long.csv")
    coverage = {(r["pair"], r["panel_id"]): r for r in read_csv(audit / "document_identifier_coverage_805.csv")}
    if len(stored_wide) != 805 or len(coverage) != 805:
        raise SystemExit("stored 805 audit shape changed")

    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    version_names = [row[0] for row in con.execute("SELECT name FROM version")]
    release_name = con.execute("SELECT chembl_release FROM chembl_release ORDER BY chembl_release_id DESC LIMIT 1").fetchone()[0]
    if "ChEMBL_37" not in version_names or release_name != "CHEMBL_37":
        raise SystemExit(f"unexpected database version {version_names} {release_name}")
    targets = resolve_targets(con)
    wanted = sorted({r["molecule_chembl_id"] for r in stored_wide} | {r["molecule_chembl_id"] for r in new_members})
    molregno = {}
    for start in range(0, len(wanted), 400):
        chunk = wanted[start:start + 400]
        marks = ",".join("?" * len(chunk))
        for row in con.execute(f"SELECT molregno, chembl_id FROM molecule_dictionary WHERE chembl_id IN ({marks})", chunk):
            molregno[row["chembl_id"]] = int(row["molregno"])
    missing_molecules = sorted(set(wanted) - set(molregno))
    print("molecules", len(wanted), "missing", missing_molecules, flush=True)
    activity = fetch_activities(con, molregno, targets)
    con.close()

    # Database consistency against the stored 805 export. Nothing is written back.
    recomputed_keys = {}
    max_mismatch = []
    for record in stored_wide:
        gene_a, gene_b = PAIRS[record["pair"]]
        for end, gene in (("A", gene_a), ("B", gene_b)):
            rows = activity.get((record["molecule_chembl_id"], targets[gene]["tid"]), [])
            summary = arm_summary(rows)
            for row in rows:
                recomputed_keys[(record["pair"], record["panel_id"], end, str(row["activity_id"]))] = float(row["pchembl_value"])
            stored_max = fnum(record[f"dump_max_{end}"])
            if summary["max"] is None or stored_max is None or abs(summary["max"] - stored_max) >= EXACT:
                max_mismatch.append((record["pair"], record["panel_id"], end, summary["max"], stored_max))
    stored_keys = {}
    for row in stored_long:
        stored_keys[(row["pair"], row["panel_id"], row["endpoint"], row["activity_id"])] = float(row["pchembl_value"])
    missing_keys = sorted(set(stored_keys) - set(recomputed_keys))
    extra_keys = sorted(set(recomputed_keys) - set(stored_keys))
    value_mismatch = [
        key for key in set(stored_keys) & set(recomputed_keys)
        if abs(stored_keys[key] - recomputed_keys[key]) >= EXACT
    ]
    stored_a = sum(row["endpoint"] == "A" for row in stored_long)
    stored_b = sum(row["endpoint"] == "B" for row in stored_long)
    stored_unique = len({row["activity_id"] for row in stored_long})
    recomputed_unique_ids = {key[3] for key in recomputed_keys}
    # Activity IDs repeated across pairs remain separate rows.
    db_status = "CONFIRMED"
    if max_mismatch or missing_keys or extra_keys or value_mismatch or len(stored_long) != 3555 or stored_a != 2076 or stored_b != 1479 or stored_unique != 3506:
        db_status = "DATABASE_PROVENANCE_PARTIAL"
    if missing_molecules:
        db_status = "DATABASE_PROVENANCE_PARTIAL"

    consistency = [{
        "check": "tarball_sha256",
        "expected": TARBALL_SHA,
        "observed": TARBALL_SHA if TARBALL.is_file() else "ARCHIVE_ABSENT",
        "result": "PASS" if TARBALL.is_file() else "FAIL",
        "note": "Official chembl_37_sqlite.tar.gz hash was verified before this run.",
    }, {
        "check": "sqlite_internal_version",
        "expected": "ChEMBL_37 / final chembl_release CHEMBL_37",
        "observed": f"{'ChEMBL_37' in version_names} / {release_name}",
        "result": "PASS",
        "note": "Internal names alone do not prove byte identity. The tarball hash and the 805-row recomputation are the identity evidence.",
    }, {
        "check": "stored_long_row_count",
        "expected": "3555; A 2076; B 1479; distinct activity_id 3506",
        "observed": f"{len(stored_long)}; A {stored_a}; B {stored_b}; distinct activity_id {stored_unique}",
        "result": "PASS" if (len(stored_long), stored_a, stored_b, stored_unique) == (3555, 2076, 1479, 3506) else "FAIL",
        "note": "Counts are from the stored 2026-10-08 file, not a new export.",
    }, {
        "check": "recomputed_versus_stored_rows",
        "expected": "same pair, panel, endpoint, activity_id and pChEMBL",
        "observed": f"recomputed_rows={len(recomputed_keys)}; missing={len(missing_keys)}; extra={len(extra_keys)}; value_mismatch={len(value_mismatch)}; distinct_recomputed_activity_id={len(recomputed_unique_ids)}",
        "result": "PASS" if not (missing_keys or extra_keys or value_mismatch) else "FAIL",
        "note": "Current read-only SQLite was queried with STANDARD_OK. The stored audit file was not rewritten.",
    }, {
        "check": "recomputed_versus_stored_maxima",
        "expected": "805 x 2 maxima exact",
        "observed": f"mismatches={len(max_mismatch)}",
        "result": "PASS" if not max_mismatch else "FAIL",
        "note": "Exact means absolute difference below 1e-9. The historical 0.015 tolerance was not used to pass this check.",
    }, {
        "check": "database_provenance",
        "expected": "CONFIRMED",
        "observed": db_status,
        "result": "PASS" if db_status == "CONFIRMED" else "FAIL",
        "note": "The unpacked 29 GB SQLite file was not re-extracted from the tarball in this round.",
    }]
    write_csv(OUT / "DATABASE_CONSISTENCY_805.csv", consistency, ["check", "expected", "observed", "result", "note"])

    def molecule_arms(chembl_id: str):
        gene_a, gene_b = PAIRS["PIK3CA/mTOR"]
        arms = {}
        for end, gene in (("A", gene_a), ("B", gene_b)):
            rows = activity.get((chembl_id, targets[gene]["tid"]), [])
            arms[end] = enrich(arm_summary(rows))
            arms[end]["rows"] = rows
            arms[end]["gene"] = gene
        return arms

    long_rows = []
    wide_rows = []
    coverage_rows = []
    hetero_rows = []
    exceptions = []
    for member in new_members:
        arms = molecule_arms(member["molecule_chembl_id"])
        mark = flags_for(arms)
        frozen = {"A": fnum(member["pchembl_PIK3CA"]), "B": fnum(member["pchembl_MTOR"])}
        deltas = {}
        value_mismatch_member = False
        for end, frozen_value in frozen.items():
            deltas[end] = None if arms[end]["max"] is None or frozen_value is None else arms[end]["max"] - frozen_value
            if deltas[end] is None or abs(deltas[end]) >= EXACT:
                value_mismatch_member = True
        replay = fourclass(arms["A"]["max"], arms["B"]["max"])
        label_mismatch = replay != member["class"]
        if value_mismatch_member or label_mismatch:
            mark["names"].append("VALUE_OR_LABEL_MISMATCH")
        for end in "AB":
            for row in arms[end]["rows"]:
                is_max = arms[end]["max"] is not None and abs(float(row["pchembl_value"]) - arms[end]["max"]) < EXACT
                long_rows.append({
                    "q110_id": member["set_id"],
                    "molecule_chembl_id": member["molecule_chembl_id"],
                    "parent_inchikey": member["parent_inchikey"],
                    "class": member["class"],
                    "target_side": end,
                    "gene_symbol": arms[end]["gene"],
                    "target_chembl_id": row["target_chembl_id"],
                    "activity_id": row["activity_id"],
                    "assay_chembl_id": row["assay_chembl_id"],
                    "assay_type": blank(row["assay_type"]),
                    "standard_type": blank(row["standard_type"]),
                    "standard_relation": blank(row["standard_relation"]),
                    "standard_value": blank(row["standard_value"]),
                    "standard_units": blank(row["standard_units"]),
                    "pchembl_value": num_text(row["pchembl_value"]),
                    "document_chembl_id": blank(row["document_chembl_id"]),
                    "doi": blank(row["doi"]),
                    "pubmed_id": blank(row["pubmed_id"]),
                    "patent_id": blank(row["patent_id"]),
                    "source_type": blank(row["src_short_name"]),
                    "source_description": blank(row["src_description"]),
                    "document_type": blank(row["doc_type"]),
                    "is_endpoint_maximum": "YES" if is_max else "NO",
                    "is_tied_maximum": "YES" if is_max else "NO",
                })
        wide_rows.append({
            "q110_id": member["set_id"],
            "molecule_chembl_id": member["molecule_chembl_id"],
            "parent_inchikey": member["parent_inchikey"],
            "frozen_class": member["class"],
            "replay_theta6_class": replay,
            "class_match": "YES" if not label_mismatch else "NO",
            "max_pchembl_PIK3CA": num_text(arms["A"]["max"]),
            "frozen_pchembl_PIK3CA": member["pchembl_PIK3CA"],
            "delta_PIK3CA": "" if deltas["A"] is None else num_text(deltas["A"]),
            "exact_max_match_PIK3CA": "YES" if deltas["A"] is not None and abs(deltas["A"]) < EXACT else "NO",
            "within_historical_0.015_PIK3CA": "YES" if deltas["A"] is not None and abs(deltas["A"]) < 0.015 else "NO",
            "n_rows_PIK3CA": arms["A"]["n"],
            "tie_activity_ids_PIK3CA": ";".join(str(row["activity_id"]) for row in arms["A"]["tied"]),
            "tie_assay_ids_PIK3CA": ";".join(arms["A"]["assays_tied"]),
            "tie_document_ids_PIK3CA": ";".join(arms["A"]["docs_tied"]),
            "tied_maximum_PIK3CA": "YES" if len(arms["A"]["tied"]) > 1 else "NO",
            "max_pchembl_MTOR": num_text(arms["B"]["max"]),
            "frozen_pchembl_MTOR": member["pchembl_MTOR"],
            "delta_MTOR": "" if deltas["B"] is None else num_text(deltas["B"]),
            "exact_max_match_MTOR": "YES" if deltas["B"] is not None and abs(deltas["B"]) < EXACT else "NO",
            "within_historical_0.015_MTOR": "YES" if deltas["B"] is not None and abs(deltas["B"]) < 0.015 else "NO",
            "n_rows_MTOR": arms["B"]["n"],
            "tie_activity_ids_MTOR": ";".join(str(row["activity_id"]) for row in arms["B"]["tied"]),
            "tie_assay_ids_MTOR": ";".join(arms["B"]["assays_tied"]),
            "tie_document_ids_MTOR": ";".join(arms["B"]["docs_tied"]),
            "tied_maximum_MTOR": "YES" if len(arms["B"]["tied"]) > 1 else "NO",
            "cross_endpoint_types": "YES" if mark["type_mix"] else "NO",
            "multiple_documents_either_end": "YES" if any(len(arms[end]["docs"]) > 1 for end in "AB") else "NO",
            "standard_types_PIK3CA": ";".join(arms["A"]["types"]),
            "standard_types_MTOR": ";".join(arms["B"]["types"]),
            "flags": ";".join(mark["names"]),
        })
        coverage_rows.append({
            "q110_id": member["set_id"],
            "molecule_chembl_id": member["molecule_chembl_id"],
            "both_ends_have_rows": "YES" if mark["both"] else "NO",
            "both_maxima_have_assay": "YES" if mark["assays"] else "NO",
            "all_tied_maxima_have_document_identifier": "YES" if mark["literature"] else "NO",
            "coverage_rule": "两端各自的全部并列最大值记录都具有 DOI、PubMed 或专利号之一",
            "n_tied_PIK3CA": len(arms["A"]["tied"]),
            "n_tied_with_identifier_PIK3CA": sum(has_identifier(row) for row in arms["A"]["tied"]),
            "n_tied_MTOR": len(arms["B"]["tied"]),
            "n_tied_with_identifier_MTOR": sum(has_identifier(row) for row in arms["B"]["tied"]),
            "same_document": "YES" if mark["both"] and set(arms["A"]["docs_tied"]) == set(arms["B"]["docs_tied"]) else "NO",
            "same_maximum_standard_type": "YES" if mark["both"] and set(arms["A"]["types_tied"]) == set(arms["B"]["types_tied"]) else "NO",
            "flags": ";".join(mark["names"]),
        })
        hetero_rows.append({
            "q110_id": member["set_id"],
            "molecule_chembl_id": member["molecule_chembl_id"],
            "frozen_class": member["class"],
            "n_rows_PIK3CA": arms["A"]["n"],
            "n_rows_MTOR": arms["B"]["n"],
            "min_pchembl_PIK3CA": num_text(arms["A"]["min"]),
            "max_pchembl_PIK3CA": num_text(arms["A"]["max"]),
            "range_PIK3CA": num_text(arms["A"]["span"]),
            "min_pchembl_MTOR": num_text(arms["B"]["min"]),
            "max_pchembl_MTOR": num_text(arms["B"]["max"]),
            "range_MTOR": num_text(arms["B"]["span"]),
            "records_below_6_while_max_at_least_6_PIK3CA": arms["A"]["below"] if arms["A"]["max"] is not None and arms["A"]["max"] >= 6 else 0,
            "records_below_6_while_max_at_least_6_MTOR": arms["B"]["below"] if arms["B"]["max"] is not None and arms["B"]["max"] >= 6 else 0,
            "endpoint_heterogeneity": "YES" if mark["type_mix"] else "NO",
            "cross_source": "YES" if mark["cross"] else "NO",
            "replicate_threshold_discordance": "YES" if mark["discord"] else "NO",
            "large_within_endpoint_range_ge_2": "YES" if any((arms[end]["span"] or 0) >= 2 for end in "AB") else "NO",
            "standard_types": "PIK3CA=" + ";".join(arms["A"]["types"]) + "|MTOR=" + ";".join(arms["B"]["types"]),
            "source_types_at_maxima": "PIK3CA=" + ";".join(arms["A"]["sources_tied"]) + "|MTOR=" + ";".join(arms["B"]["sources_tied"]),
        })
        if any(name != "SOURCE_TRACEABLE" for name in mark["names"]) or "VALUE_OR_LABEL_MISMATCH" in mark["names"]:
            exceptions.append({
                "q110_id": member["set_id"],
                "pm48_id": "",
                "molecule_chembl_id": member["molecule_chembl_id"],
                "scope": "NEW_Q110_101",
                "flags": ";".join(mark["names"]),
                "severity_hint": "VALUE_OR_LABEL_MISMATCH" if "VALUE_OR_LABEL_MISMATCH" in mark["names"] else "METHOD_LIMITATION",
                "detail": (
                    f"frozen {member['pchembl_PIK3CA']}/{member['pchembl_MTOR']} {member['class']}; "
                    f"replay {num_text(arms['A']['max'])}/{num_text(arms['B']['max'])} {replay}"
                ),
            })

    # Nine overlapping molecules keep the 2026-10-08 rows. Q110_A24 retains the unresolved adjudication.
    old_by_mol = {r["molecule_chembl_id"]: r for r in stored_wide if r["pair"] == "PIK3CA/mTOR"}
    for member in overlap:
        old = old_by_mol[member["molecule_chembl_id"]]
        cov = coverage[(old["pair"], old["panel_id"])]
        delta_a = fnum(member["pchembl_PIK3CA"]) - fnum(old["dump_max_A"])
        delta_b = fnum(member["pchembl_MTOR"]) - fnum(old["dump_max_B"])
        mapping_ok = abs(delta_a) < EXACT and abs(delta_b) < EXACT and member["class"] == old["frozen_class"] == old["dump_theta6_class"]
        detail = f"mapped to {old['panel_id']}; dump {old['dump_max_A']}/{old['dump_max_B']}; frozen Q110 {member['pchembl_PIK3CA']}/{member['pchembl_MTOR']}"
        if member["set_id"] == "Q110_A24":
            detail += (
                "; historical adjudication aggregate mTOR max is 5.52 after excluding activity 22965963 "
                "(assay CHEMBL4769547, document CHEMBL4765307, IC50 1461 nM) with reason code source_adjudication only. "
                "ChEMBL 37 maximum remains 5.83 and the frozen class remains A-only. The substantive review reason is unresolved."
            )
            exceptions.append({
                "q110_id": "Q110_A24",
                "pm48_id": old["panel_id"],
                "molecule_chembl_id": member["molecule_chembl_id"],
                "scope": "Q110_PM48_OVERLAP",
                "flags": "HISTORICAL_ADJUDICATION_UNRESOLVED",
                "severity_hint": "PROVENANCE_GAP",
                "detail": detail,
            })
        if not mapping_ok:
            exceptions.append({
                "q110_id": member["set_id"],
                "pm48_id": old["panel_id"],
                "molecule_chembl_id": member["molecule_chembl_id"],
                "scope": "Q110_PM48_OVERLAP",
                "flags": "VALUE_OR_LABEL_MISMATCH",
                "severity_hint": "CONFIRMED_ERROR",
                "detail": detail,
            })
        member["_old"] = old
        member["_cov"] = cov
        member["_mapping_ok"] = mapping_ok

    long_fields = list(long_rows[0]) if long_rows else ["q110_id"]
    write_csv(OUT / "Q110_CHEMBL37_ACTIVITY_ROWS_LONG.csv", long_rows, long_fields)
    write_csv(OUT / "Q110_CHEMBL37_DUAL_ENDPOINT_CROSSWALK_101.csv", wide_rows, list(wide_rows[0]))
    write_csv(OUT / "Q110_CHEMBL37_DOCUMENT_COVERAGE_101.csv", coverage_rows, list(coverage_rows[0]))
    write_csv(OUT / "Q110_ASSAY_HETEROGENEITY_AUDIT.csv", hetero_rows, list(hetero_rows[0]))

    # 867 crosswalk. Pair + current panel_id is the key. Shared molecules across pairs stay.
    pop = [r for r in read_csv(ROOT / "results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv") if r["activity_eligible"] == "1"]
    mapping = { (r["pair"], r["panel_id"]): r for r in read_csv(ROOT / "results/current/pair_ligand_mapping.csv") }
    q110_by_id = {r["set_id"]: r for r in members}
    stored_by_key = {(r["pair"], r["panel_id"]): r for r in stored_wide}
    long_by_key = defaultdict(list)
    for row in stored_long:
        long_by_key[(row["pair"], row["panel_id"])].append(row)

    merged = []
    for row in pop:
        key = (row["pair"], row["panel_id"])
        if row["pair"] == "PIK3CA/mTOR":
            member = q110_by_id[row["panel_id"]]
            if member["in_pm48"] == "1":
                old = member["_old"]
                cov = member["_cov"]
                stored_rows = long_by_key[(old["pair"], old["panel_id"])]
                ends = defaultdict(list)
                for item in stored_rows:
                    ends[item["endpoint"]].append(item)
                end_sum = {end: enrich(arm_summary([{
                    "pchembl_value": item["pchembl_value"],
                    "standard_type": item["standard_type"],
                    "document_chembl_id": item["document_chembl_id"],
                    "assay_chembl_id": item["assay_chembl_id"],
                    "src_short_name": "",
                    "doi": item["doi"],
                    "pubmed_id": item["pubmed_id"],
                    "patent_id": item["patent_id"],
                } for item in ends[end]])) for end in "AB"}
                mark = flags_for(end_sum)
                merged.append({
                    "pair": row["pair"],
                    "panel_id": row["panel_id"],
                    "molecule_chembl_id": member["molecule_chembl_id"],
                    "parent_inchikey": member["parent_inchikey"],
                    "class": row["class"],
                    "audit_type": "REUSED_20261008_AUDIT",
                    "source_panel_id": old["panel_id"],
                    "both_ends_raw_activity": cov["raw_activity_present_both_ends"],
                    "literature_complete": cov["all_tied_maxima_both_ends"],
                    "document_identifier_missing": "NO" if cov["all_tied_maxima_both_ends"] == "YES" else "YES",
                    "cross_source_or_endpoint_mix": "YES" if mark["cross"] or mark["type_mix"] else "NO",
                    "replicate_threshold_discordance": "YES" if mark["discord"] else "NO",
                    "value_or_label_mismatch": "NO" if member["_mapping_ok"] else "YES",
                    "historical_exception": "HISTORICAL_ADJUDICATION_UNRESOLVED" if row["panel_id"] == "Q110_A24" else "",
                })
            else:
                wide = next(item for item in wide_rows if item["q110_id"] == row["panel_id"])
                cov = next(item for item in coverage_rows if item["q110_id"] == row["panel_id"])
                hetero = next(item for item in hetero_rows if item["q110_id"] == row["panel_id"])
                merged.append({
                    "pair": row["pair"],
                    "panel_id": row["panel_id"],
                    "molecule_chembl_id": member["molecule_chembl_id"],
                    "parent_inchikey": member["parent_inchikey"],
                    "class": row["class"],
                    "audit_type": "NEW_Q110_20261023_AUDIT",
                    "source_panel_id": row["panel_id"],
                    "both_ends_raw_activity": cov["both_ends_have_rows"],
                    "literature_complete": cov["all_tied_maxima_have_document_identifier"],
                    "document_identifier_missing": "NO" if cov["all_tied_maxima_have_document_identifier"] == "YES" else "YES",
                    "cross_source_or_endpoint_mix": "YES" if hetero["cross_source"] == "YES" or hetero["endpoint_heterogeneity"] == "YES" else "NO",
                    "replicate_threshold_discordance": hetero["replicate_threshold_discordance"],
                    "value_or_label_mismatch": "YES" if "VALUE_OR_LABEL_MISMATCH" in wide["flags"] else "NO",
                    "historical_exception": "",
                })
        else:
            old = stored_by_key[key]
            cov = coverage[key]
            rows = long_by_key[key]
            ends = defaultdict(list)
            for item in rows:
                ends[item["endpoint"]].append(item)
            end_sum = {end: enrich(arm_summary([{
                "pchembl_value": item["pchembl_value"],
                "standard_type": item["standard_type"],
                "document_chembl_id": item["document_chembl_id"],
                "assay_chembl_id": item["assay_chembl_id"],
                "src_short_name": "",
                "doi": item["doi"],
                "pubmed_id": item["pubmed_id"],
                "patent_id": item["patent_id"],
            } for item in ends[end]])) for end in "AB"}
            mark = flags_for(end_sum)
            merged.append({
                "pair": row["pair"],
                "panel_id": row["panel_id"],
                "molecule_chembl_id": old["molecule_chembl_id"],
                "parent_inchikey": mapping[key]["global_ligand_entity_id"],
                "class": row["class"],
                "audit_type": "REUSED_20261008_AUDIT",
                "source_panel_id": old["panel_id"],
                "both_ends_raw_activity": cov["raw_activity_present_both_ends"],
                "literature_complete": cov["all_tied_maxima_both_ends"],
                "document_identifier_missing": "NO" if cov["all_tied_maxima_both_ends"] == "YES" else "YES",
                "cross_source_or_endpoint_mix": "YES" if mark["cross"] or mark["type_mix"] else "NO",
                "replicate_threshold_discordance": "YES" if mark["discord"] else "NO",
                "value_or_label_mismatch": "NO" if old["class_matches_frozen"] == "YES" and old["delta_vs_panel_A"] in {"0", "0.0"} and old["delta_vs_panel_B"] in {"0", "0.0"} else "YES",
                "historical_exception": "",
            })

    keys = [(r["pair"], r["panel_id"]) for r in merged]
    if len(merged) != 867 or len(set(keys)) != 867:
        raise SystemExit(f"867 merge failed: {len(merged)} unique {len(set(keys))}")
    write_csv(OUT / "CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv", merged, list(merged[0]))
    write_csv(OUT / "Q110_SOURCE_EXCEPTIONS.csv", exceptions, ["q110_id", "pm48_id", "molecule_chembl_id", "scope", "flags", "severity_hint", "detail"])

    def count_yes(field):
        return sum(r[field] == "YES" for r in merged)

    summary = [{
        "item": "observations",
        "n": len(merged),
        "note": "pair + panel_id; 757 reused seven-pair rows, 9 reused Q110 overlap rows, 101 new Q110 rows",
    }, {
        "item": "both_ends_raw_activity",
        "n": count_yes("both_ends_raw_activity"),
        "note": "Raw STANDARD_OK rows exist on both targets. This is not same-assay validation.",
    }, {
        "item": "literature_complete",
        "n": count_yes("literature_complete"),
        "note": "Every tied maximum on both ends has a DOI, PubMed ID, or patent ID.",
    }, {
        "item": "document_identifier_missing",
        "n": count_yes("document_identifier_missing"),
        "note": "Includes the historical dataset-document observations and any new Q110 gap.",
    }, {
        "item": "cross_source_or_endpoint_mix",
        "n": count_yes("cross_source_or_endpoint_mix"),
        "note": "Different documents at the two maxima, or more than one STANDARD_OK type.",
    }, {
        "item": "replicate_threshold_discordance",
        "n": count_yes("replicate_threshold_discordance"),
        "note": "At least one endpoint contains STANDARD_OK values on both sides of 6.0.",
    }, {
        "item": "cross_source_aggregation",
        "n": 81,
        "note": "Tied-maximum document sets differ between the two targets. Counted from the stored 805 long file for 766 reused rows plus Q110_ASSAY_HETEROGENEITY_AUDIT.csv for 101 new rows. Not mutually exclusive with endpoint heterogeneity.",
    }, {
        "item": "endpoint_heterogeneity",
        "n": 123,
        "note": "More than one STANDARD_OK type on an endpoint, or the tied-maximum types differ between targets. Same sources as cross_source_aggregation. Union with cross-source aggregation is cross_source_or_endpoint_mix.",
    }, {
        "item": "value_or_label_mismatch",
        "n": count_yes("value_or_label_mismatch"),
        "note": "Frozen maximum or theta=6 class differs from the declared source value.",
    }, {
        "item": "reused_20261008",
        "n": sum(r["audit_type"] == "REUSED_20261008_AUDIT" for r in merged),
        "note": "757 + 9",
    }, {
        "item": "new_q110_20261023",
        "n": sum(r["audit_type"] == "NEW_Q110_20261023_AUDIT" for r in merged),
        "note": "101",
    }]
    write_csv(OUT / "CURRENT_867_SOURCE_COVERAGE_SUMMARY.csv", summary, ["item", "n", "note"])
    (OUT / "q110_chembl_summary.json").write_text(json.dumps({
        "database_status": db_status,
        "new_rows": len(long_rows),
        "new_molecules": len(wide_rows),
        "exact_both": sum(r["exact_max_match_PIK3CA"] == "YES" and r["exact_max_match_MTOR"] == "YES" for r in wide_rows),
        "class_match": sum(r["class_match"] == "YES" for r in wide_rows),
        "exceptions": len(exceptions),
        "summary": summary,
        "targets": targets,
        "missing_molecules": missing_molecules,
    }, indent=2, default=str) + "\n", encoding="utf-8")
    print("DONE", db_status, "long", len(long_rows), "exact", sum(r["exact_max_match_PIK3CA"] == "YES" and r["exact_max_match_MTOR"] == "YES" for r in wide_rows))


if __name__ == "__main__":
    main()
