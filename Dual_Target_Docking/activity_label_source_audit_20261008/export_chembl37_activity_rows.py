#!/usr/bin/env python3
"""Export ChEMBL 37 activity rows for the frozen 805 records.

Reads the restored SQLite and the existing mapping. Does not write formal inputs.
The row filter is the historical STANDARD_OK predicate used by harvest_activities:
non-null pChEMBL and standard_type in IC50, Ki, Kd, EC50, Potency, IC50app, Ki app.
No confidence, assay-type, or keep=1 filter is added.
"""
from __future__ import annotations

import csv
import io
import os
import sqlite3
import subprocess
from collections import defaultdict
from pathlib import Path

DB = Path(os.environ.get("CHEMBL37_DB", "/home/gwj/chembl_37_release/chembl_37/chembl_37_sqlite/chembl_37.db"))
REPO = Path("/tmp/v4_2_final_audited_release_20260930")
OUT = Path("/mnt/d/CADD paper exercise/dual target docking/activity_label_source_audit_20261008")
MAPPING = REPO / "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/pair_ligand_mapping.csv"
CROSSWALK = OUT / "activity_label_source_crosswalk_805.csv"
STANDARD_OK = ("IC50", "Ki", "Kd", "EC50", "Potency", "IC50app", "Ki app")
TOL = 0.015
PAIRS = {
    "EGFR/HER2": ("EGFR", "ERBB2", "2026-07-23 ChEMBL REST panel; dump rows use STANDARD_OK"),
    "AChE/BChE": ("ACHE", "BCHE", "2026-07-23 ChEMBL REST panel; dump rows use STANDARD_OK"),
    "PIK3CA/mTOR": ("PIK3CA", "MTOR", "2026-07-23 ChEMBL REST panel; dump rows use STANDARD_OK"),
    "F2/F10": ("F2", "F10", "ChEMBL 37 dump extract; MAX of STANDARD_OK pChEMBL"),
    "JAK1/TYK2": ("JAK1", "TYK2", "ChEMBL 37 dump extract; MAX of STANDARD_OK pChEMBL"),
    "JAK1/JAK2": ("JAK1", "JAK2", "ChEMBL 37 dump extract; MAX of STANDARD_OK pChEMBL"),
    "PPARG/PPARA": ("PPARG", "PPARA", "ChEMBL 37 dump extract; MAX of STANDARD_OK pChEMBL"),
    "PPARA/PPARD": ("PPARA", "PPARD", "ChEMBL 37 dump extract; MAX of STANDARD_OK pChEMBL"),
}


def fnum(value):
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def fourclass(pa, pb):
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


def blank(value):
    if value is None:
        return ""
    return str(value)


def literature_flag(doi, pubmed, patent):
    if doi or pubmed or patent:
        return "present"
    return "missing"


def resolve_targets(con, genes):
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
    found = {}
    for row in con.execute(sql):
        gene = (row["gene"] or "").upper()
        if gene in genes and gene not in found:
            found[gene] = {
                "tid": int(row["tid"]),
                "chembl": row["chembl_id"],
                "uniprot": row["accession"],
                "pref_name": row["pref_name"],
            }
    missing = sorted(genes - set(found))
    if missing:
        raise SystemExit(f"unresolved genes: {missing}")
    return found


def main() -> None:
    if not DB.exists():
        raise SystemExit(f"database not ready: {DB}")
    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    release_tables = [row[0] for row in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name IN ('version','chembl_release')")]
    release = ""
    if "version" in release_tables:
        cols = [row[1] for row in con.execute("PRAGMA table_info(version)")]
        print("version columns", cols)
        for row in con.execute("SELECT * FROM version"):
            print("VERSION", dict(row))
            text = " ".join(str(value) for value in dict(row).values())
            if "chembl_37" in text.lower() or text.strip() in {"37", "CHEMBL_37"}:
                release = text
    if not release:
        raise SystemExit("internal version table did not identify chembl_37; export stopped")

    genes = {gene for pair in PAIRS.values() for gene in pair[:2]}
    targets = resolve_targets(con, genes)
    for gene, info in sorted(targets.items()):
        print(f"TARGET {gene} {info['chembl']} tid={info['tid']} {info['uniprot']}")

    index = list(csv.DictReader(CROSSWALK.open()))
    if len(index) != 805:
        raise SystemExit(f"expected 805 index rows, got {len(index)}")
    aggregate_text = subprocess.check_output(
        ["git", "show", "ab8a709ec:Dual_Target_Docking/data/processed/activity_adjudication/ligand_activity_aggregate_v1.csv"],
        cwd=REPO, text=True,
    )
    aggregate = {
        (row["pair"], row["ligand"]): row
        for row in csv.DictReader(io.StringIO(aggregate_text))
    }
    frozen = {
        (row["pair"], row["panel_id"]): row["class"]
        for row in csv.DictReader(MAPPING.open())
        if row["activity_eligible"] == "1" and row["parent_independent"] == "1"
    }

    mol_ids = sorted({row["molecule_chembl_id"] for row in index})
    molregno = {}
    for start in range(0, len(mol_ids), 400):
        chunk = mol_ids[start:start + 400]
        marks = ",".join("?" * len(chunk))
        for row in con.execute(f"SELECT molregno, chembl_id FROM molecule_dictionary WHERE chembl_id IN ({marks})", chunk):
            molregno[row["chembl_id"]] = int(row["molregno"])
    missing_mols = sorted(set(mol_ids) - set(molregno))
    print("molecules", len(mol_ids), "resolved", len(molregno), "missing", len(missing_mols))

    doc_cols = {row[1] for row in con.execute("PRAGMA table_info(docs)")}
    act_cols = {row[1] for row in con.execute("PRAGMA table_info(activities)")}
    for required in ("activity_id", "pchembl_value", "standard_type", "standard_relation", "standard_value", "standard_units"):
        if required not in act_cols:
            raise SystemExit(f"activities.{required} absent")
    doi_col = "doi" if "doi" in doc_cols else None
    pubmed_col = "pubmed_id" if "pubmed_id" in doc_cols else None
    patent_col = "patent_id" if "patent_id" in doc_cols else None

    tids = sorted({info["tid"] for info in targets.values()})
    tid_marks = ",".join("?" * len(tids))
    type_marks = ",".join("?" * len(STANDARD_OK))
    rows_by_mol_tid = defaultdict(list)
    mol_list = list(molregno.items())
    for start in range(0, len(mol_list), 200):
        chunk = mol_list[start:start + 200]
        id_by_regno = {regno: chembl for chembl, regno in chunk}
        marks = ",".join("?" * len(chunk))
        sql = f"""
        SELECT md.chembl_id AS molecule_chembl_id, md.molregno,
               td.chembl_id AS target_chembl_id, ass.tid,
               act.activity_id, ass.chembl_id AS assay_chembl_id, ass.assay_type,
               act.standard_type, act.standard_relation, act.standard_value, act.standard_units,
               act.pchembl_value, docs.chembl_id AS document_chembl_id,
               {('docs.' + doi_col) if doi_col else "NULL"} AS doi,
               {('docs.' + pubmed_col) if pubmed_col else "NULL"} AS pubmed_id,
               {('docs.' + patent_col) if patent_col else "NULL"} AS patent_id
        FROM activities act
        JOIN assays ass ON ass.assay_id = act.assay_id
        JOIN docs ON docs.doc_id = ass.doc_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        JOIN target_dictionary td ON td.tid = ass.tid
        WHERE ass.tid IN ({tid_marks})
          AND act.molregno IN ({marks})
          AND act.pchembl_value IS NOT NULL
          AND act.standard_type IN ({type_marks})
        """
        params = tids + [regno for _, regno in chunk] + list(STANDARD_OK)
        for row in con.execute(sql, params):
            rows_by_mol_tid[(row["molecule_chembl_id"], int(row["tid"]))].append(row)
        print("queried", start + len(chunk), "molecules", flush=True)

    long_path = OUT / "chembl37_activity_rows_long.csv"
    wide_path = OUT / "chembl37_dual_endpoint_crosswalk_805.csv"
    long_fields = [
        "pair", "panel_id", "molecule_chembl_id", "molregno", "endpoint", "gene",
        "target_chembl_id", "activity_id", "assay_chembl_id", "assay_type",
        "standard_type", "standard_relation", "standard_value", "standard_units",
        "pchembl_value", "document_chembl_id", "doi", "pubmed_id", "patent_id",
        "literature_status", "is_tied_maximum", "extraction_rule",
    ]
    wide_fields = [
        "pair", "panel_id", "molecule_chembl_id", "frozen_class", "extraction_rule",
        "gene_A", "target_chembl_A", "gene_B", "target_chembl_B",
        "n_rows_A", "n_rows_B", "dump_max_A", "dump_max_B",
        "tie_activity_ids_A", "tie_activity_ids_B",
        "tie_literature_A", "tie_literature_B",
        "historical_panel_pA", "historical_panel_pB",
        "delta_vs_panel_A", "delta_vs_panel_B",
        "historical_aggregate_max_A", "historical_aggregate_max_B",
        "delta_vs_aggregate_A", "delta_vs_aggregate_B",
        "dump_theta6_class", "class_matches_frozen",
        "source_conclusion", "pause_reason",
    ]
    paused = []
    n_both = n_lit = 0
    with long_path.open("w", encoding="utf-8", newline="") as long_handle, wide_path.open("w", encoding="utf-8", newline="") as wide_handle:
        long_writer = csv.DictWriter(long_handle, fieldnames=long_fields, lineterminator="\n")
        wide_writer = csv.DictWriter(wide_handle, fieldnames=wide_fields, lineterminator="\n")
        long_writer.writeheader()
        wide_writer.writeheader()
        for record in index:
            pair = record["pair"]
            gene_a, gene_b, rule = PAIRS[record["pair"]]
            target_a, target_b = targets[gene_a], targets[gene_b]
            mol = record["molecule_chembl_id"]
            arms = {
                "A": rows_by_mol_tid.get((mol, target_a["tid"]), []),
                "B": rows_by_mol_tid.get((mol, target_b["tid"]), []),
            }
            maxima = {}
            tie_ids = {}
            tie_lit = {}
            for end, gene, target in (("A", gene_a, target_a), ("B", gene_b, target_b)):
                values = [float(row["pchembl_value"]) for row in arms[end]]
                maxima[end] = max(values) if values else None
                tied = []
                if maxima[end] is not None:
                    tied = [row for row in arms[end] if abs(float(row["pchembl_value"]) - maxima[end]) < 1e-9]
                tie_ids[end] = ";".join(str(row["activity_id"]) for row in tied)
                tie_lit[end] = "missing" if not tied else (
                    "present" if all(literature_flag(row["doi"], row["pubmed_id"], row["patent_id"]) == "present" for row in tied) else "missing"
                )
                for row in arms[end]:
                    is_tie = maxima[end] is not None and abs(float(row["pchembl_value"]) - maxima[end]) < 1e-9
                    long_writer.writerow({
                        "pair": pair,
                        "panel_id": record["panel_id"],
                        "molecule_chembl_id": mol,
                        "molregno": row["molregno"],
                        "endpoint": end,
                        "gene": gene,
                        "target_chembl_id": row["target_chembl_id"],
                        "activity_id": row["activity_id"],
                        "assay_chembl_id": row["assay_chembl_id"],
                        "assay_type": blank(row["assay_type"]),
                        "standard_type": blank(row["standard_type"]),
                        "standard_relation": blank(row["standard_relation"]),
                        "standard_value": blank(row["standard_value"]),
                        "standard_units": blank(row["standard_units"]),
                        "pchembl_value": row["pchembl_value"],
                        "document_chembl_id": blank(row["document_chembl_id"]),
                        "doi": blank(row["doi"]),
                        "pubmed_id": blank(row["pubmed_id"]),
                        "patent_id": blank(row["patent_id"]),
                        "literature_status": literature_flag(row["doi"], row["pubmed_id"], row["patent_id"]),
                        "is_tied_maximum": "YES" if is_tie else "NO",
                        "extraction_rule": rule,
                    })
            panel_a, panel_b = fnum(record["aggregated_pA"]), fnum(record["aggregated_pB"])
            dump_class = fourclass(maxima["A"], maxima["B"])
            frozen_class = frozen.get((pair, record["panel_id"]), "")
            reasons = []
            if not arms["A"] or not arms["B"]:
                reasons.append("missing_endpoint_rows")
            else:
                n_both += 1
            if dump_class != frozen_class:
                reasons.append(f"class {dump_class or 'missing'} != frozen {frozen_class}")
            delta_panel = {}
            for end, stored in (("A", panel_a), ("B", panel_b)):
                if maxima[end] is None or stored is None:
                    delta_panel[end] = ""
                    reasons.append(f"no_panel_comparison_{end}")
                else:
                    delta_panel[end] = maxima[end] - stored
                    if abs(delta_panel[end]) >= TOL:
                        reasons.append(f"panel_delta_{end}={delta_panel[end]:.4f}")
            if arms["A"] and arms["B"] and tie_lit["A"] == "present" and tie_lit["B"] == "present":
                n_lit += 1
            agg = aggregate.get((pair, record["panel_id"]))
            agg_a = fnum(agg["max_A"]) if agg else None
            agg_b = fnum(agg["max_B"]) if agg else None
            delta_agg = {}
            for end, stored in (("A", agg_a), ("B", agg_b)):
                if stored is None or maxima[end] is None:
                    delta_agg[end] = ""
                else:
                    delta_agg[end] = maxima[end] - stored
                    if abs(delta_agg[end]) >= TOL:
                        reasons.append(f"audit_aggregate_delta_{end}={delta_agg[end]:.4f}")
            conclusion = "PAUSED" if reasons else "OK"
            if conclusion == "PAUSED":
                paused.append((pair, record["panel_id"], ";".join(reasons)))
            wide_writer.writerow({
                "pair": pair,
                "panel_id": record["panel_id"],
                "molecule_chembl_id": mol,
                "frozen_class": frozen_class,
                "extraction_rule": rule,
                "gene_A": gene_a,
                "target_chembl_A": target_a["chembl"],
                "gene_B": gene_b,
                "target_chembl_B": target_b["chembl"],
                "n_rows_A": len(arms["A"]),
                "n_rows_B": len(arms["B"]),
                "dump_max_A": "" if maxima["A"] is None else f"{maxima['A']:.10g}",
                "dump_max_B": "" if maxima["B"] is None else f"{maxima['B']:.10g}",
                "tie_activity_ids_A": tie_ids["A"],
                "tie_activity_ids_B": tie_ids["B"],
                "tie_literature_A": tie_lit["A"],
                "tie_literature_B": tie_lit["B"],
                "historical_panel_pA": record["aggregated_pA"],
                "historical_panel_pB": record["aggregated_pB"],
                "delta_vs_panel_A": "" if delta_panel["A"] == "" else f"{delta_panel['A']:.10g}",
                "delta_vs_panel_B": "" if delta_panel["B"] == "" else f"{delta_panel['B']:.10g}",
                "historical_aggregate_max_A": "" if agg_a is None else f"{agg_a:.10g}",
                "historical_aggregate_max_B": "" if agg_b is None else f"{agg_b:.10g}",
                "delta_vs_aggregate_A": "" if delta_agg["A"] == "" else f"{delta_agg['A']:.10g}",
                "delta_vs_aggregate_B": "" if delta_agg["B"] == "" else f"{delta_agg['B']:.10g}",
                "dump_theta6_class": dump_class,
                "class_matches_frozen": "YES" if dump_class == frozen_class else "NO",
                "source_conclusion": conclusion,
                "pause_reason": ";".join(reasons),
            })
    print("BOTH_ENDPOINTS", n_both, "TIE_LITERATURE_BOTH", n_lit, "PAUSED", len(paused))
    for item in paused:
        print("PAUSE", *item)
    print("WROTE", long_path)
    print("WROTE", wide_path)


if __name__ == "__main__":
    main()
