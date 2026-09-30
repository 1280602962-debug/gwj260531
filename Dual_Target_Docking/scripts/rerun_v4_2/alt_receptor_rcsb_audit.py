#!/usr/bin/env python3
"""RCSB availability audit for V4.2 alternative-receptor design.

Structure/identity only. Does not read docking scores, RMSD QC, or AUROC.
"""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

QA = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa")
OUT = QA / "alt_receptor_rcsb_raw.json"

TARGETS = [
    {"target": "EGFR", "uniprot": "P00533", "primary_pdb": "3POZ", "domain": "kinase", "family": "kinase_ATP"},
    {"target": "HER2", "uniprot": "P04626", "primary_pdb": "3RCD", "domain": "kinase", "family": "kinase_ATP"},
    {"target": "JAK1", "uniprot": "P23458", "primary_pdb": "6N7A", "domain": "JH1_kinase", "family": "kinase_ATP"},
    {"target": "JAK2", "uniprot": "O60674", "primary_pdb": "8BXH", "domain": "JH1_kinase", "family": "kinase_ATP"},
    {"target": "TYK2", "uniprot": "P29597", "primary_pdb": "3LXP", "domain": "JH1_kinase", "family": "kinase_ATP"},
    {"target": "PIK3CA", "uniprot": "P42336", "primary_pdb": "4L23", "domain": "p110alpha_catalytic", "family": "kinase_ATP"},
    {"target": "mTOR", "uniprot": "P42345", "primary_pdb": "4JT6", "domain": "FRB_kinase", "family": "kinase_ATP"},
    {"target": "AChE", "uniprot": "P22303", "primary_pdb": "4EY7", "domain": "catalytic", "family": "AChE_BChE"},
    {"target": "BChE", "uniprot": "P06276", "primary_pdb": "4BDS", "domain": "catalytic", "family": "AChE_BChE"},
    {"target": "F2", "uniprot": "P00734", "primary_pdb": "4UDW", "domain": "serine_protease", "family": "F2_F10"},
    {"target": "F10", "uniprot": "P00742", "primary_pdb": "2JKH", "domain": "serine_protease", "family": "F2_F10"},
    {"target": "PPARG", "uniprot": "P37231", "primary_pdb": "9V8H", "domain": "LBD", "family": "PPAR_NR"},
    {"target": "PPARA", "uniprot": "Q07869", "primary_pdb": "6LXA", "domain": "LBD", "family": "PPAR_NR"},
    {"target": "PPARD", "uniprot": "Q03181", "primary_pdb": "5U3Q", "domain": "LBD", "family": "PPAR_NR"},
]


def post_json(url: str, payload: dict) -> dict:
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode())


def search_uniprot(accession: str) -> list[str]:
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": accession,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator": "exact_match",
                        "value": "Homo sapiens",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.structure_determination_methodology",
                        "operator": "exact_match",
                        "value": "experimental",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 10000},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }
    out = post_json("https://search.rcsb.org/rcsbsearch/v2/query", query)
    return [r["identifier"] for r in out.get("result_set", [])]


GQL = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct { title }
    rcsb_accession_info { deposit_date initial_release_date }
    exptl { method }
    rcsb_entry_info { resolution_combined experimental_method polymer_entity_count deposited_polymer_entity_instance_count }
    rcsb_primary_citation { pdbx_database_id_DOI pdbx_database_id_PubMed journal_abbrev title }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_strand_id rcsb_entity_polymer_type rcsb_sample_sequence_length }
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      rcsb_polymer_entity_feature { type name feature_positions { beg_seq_id end_seq_id } }
      rcsb_entity_source_organism { ncbi_scientific_name }
      uniprots { rcsb_id }
      rcsb_polymer_entity_align {
        aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
        reference_database_accession
      }
    }
    nonpolymer_entities {
      rcsb_id
      nonpolymer_comp { chem_comp { id name type } }
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}
"""


def gql_entries(ids: list[str]) -> list[dict]:
    rows = []
    for i in range(0, len(ids), 25):
        chunk = ids[i:i + 25]
        data = post_json("https://data.rcsb.org/graphql", {"query": GQL, "variables": {"ids": chunk}})
        rows.extend(data.get("data", {}).get("entries") or [])
        time.sleep(0.15)
    return rows


def main() -> None:
    QA.mkdir(parents=True, exist_ok=True)
    blob = {"targets": []}
    for spec in TARGETS:
        print("search", spec["target"], spec["uniprot"], flush=True)
        ids = search_uniprot(spec["uniprot"])
        print(" ", spec["target"], "n_entries", len(ids), flush=True)
        entries = gql_entries(ids)
        blob["targets"].append({**spec, "n_entries": len(ids), "pdb_ids": ids, "entries": entries})
        time.sleep(0.2)
    OUT.write_text(json.dumps(blob), encoding="utf-8")
    print("wrote", OUT, "bytes", OUT.stat().st_size)


if __name__ == "__main__":
    main()
