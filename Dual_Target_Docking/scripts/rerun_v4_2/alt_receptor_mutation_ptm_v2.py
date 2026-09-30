#!/usr/bin/env python3
"""Independent mutation / PTM / validation extract. No docking or AUROC."""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import gemmi

CIFDIR = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/alt_cif")
QA = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa")
PDBS = [
    "3POZ", "1XKK", "3RCD", "3PP0", "6N7A", "4EI4", "8BXH", "8BM2", "3LXP",
    "4L23", "4L2Y", "4JPS", "5XGH", "5DXT", "4JT6", "4JT5", "4JSX", "4JSV",
    "4EY7", "4EY6", "4BDS", "1P0M", "4UDW", "3SHC", "2JKH", "2Y5F",
    "9V8H", "7AWC", "6LXA", "6KAX", "5U3Q", "5U46",
]


def as_s(x) -> str:
    if x is None:
        return ""
    s = gemmi.cif.as_string(x) if not isinstance(x, str) else x
    return "" if s in {".", "?", ""} else s


def table(block, prefix, cols):
    return block.find(prefix, cols)


def extract(pdb: str) -> dict:
    path = CIFDIR / f"{pdb}.cif"
    b = gemmi.cif.read(str(path)).sole_block()
    rec = {"pdb": pdb}

    # entity mutation field
    muts = []
    et = table(b, "_entity.", ["id", "type", "pdbx_description", "pdbx_mutation"])
    if et:
        for r in et:
            muts.append({
                "entity": as_s(r[0]), "type": as_s(r[1]),
                "desc": as_s(r[2]), "pdbx_mutation": as_s(r[3]) or "NONE_IN_entity.pdbx_mutation",
            })
    rec["entity"] = muts

    # struct_ref_seq_dif full
    difs = []
    dt = table(b, "_struct_ref_seq_dif.", [
        "align_id", "pdbx_pdb_strand_id", "mon_id", "db_mon_id",
        "db_align_seq_id", "pdbx_auth_seq_num", "details",
        "pdbx_seq_db_name", "pdbx_seq_db_accession_code", "pdbx_seq_db_seq_num",
    ])
    if dt:
        for r in dt:
            difs.append({
                "align": as_s(r[0]), "strand": as_s(r[1]), "pdb_aa": as_s(r[2]),
                "db_aa": as_s(r[3]), "db_align": as_s(r[4]), "auth": as_s(r[5]),
                "details": as_s(r[6]), "db": as_s(r[7]), "acc": as_s(r[8]), "db_seq": as_s(r[9]),
            })
    rec["n_struct_ref_seq_dif"] = len(difs)
    rec["struct_ref_seq_dif"] = difs

    # modified residues
    mods = []
    mt = table(b, "_pdbx_struct_mod_residue.", [
        "label_comp_id", "auth_comp_id", "auth_asym_id", "auth_seq_id",
        "parent_comp_id", "details",
    ])
    if mt:
        for r in mt:
            mods.append({
                "label": as_s(r[0]), "auth_comp": as_s(r[1]), "chain": as_s(r[2]),
                "seq": as_s(r[3]), "parent": as_s(r[4]), "details": as_s(r[5]),
            })
    rec["mod_residues"] = mods

    # struct_site (sometimes binding site)
    sites = []
    st = table(b, "_struct_site.", ["id", "details"])
    if st:
        for r in st:
            sites.append({"id": as_s(r[0]), "details": as_s(r[1])})
    rec["struct_sites"] = sites[:20]

    # software / deposition
    rec["deposition_date"] = as_s(b.find_value("_pdbx_database_status.recvd_initial_deposition_date"))
    rec["revision"] = as_s(b.find_value("_pdbx_audit_revision_history.revision_date"))

    # entity_poly mutation / one letter
    poly = []
    pt = table(b, "_entity_poly.", ["entity_id", "type", "pdbx_strand_id", "pdbx_seq_one_letter_code_can"])
    if pt:
        for r in pt:
            seq = as_s(r[3])
            poly.append({
                "entity": as_s(r[0]), "type": as_s(r[1]), "strands": as_s(r[2]),
                "seq_len": len(seq.replace("\n", "").replace(" ", "")) if seq else 0,
            })
    rec["entity_poly"] = poly

    # keyword
    rec["keywords"] = as_s(b.find_value("_struct_keywords.pdbx_keywords"))
    rec["text"] = as_s(b.find_value("_struct_keywords.text"))
    return rec


def fetch_rcsb(pdb: str) -> dict:
    q = {
        "query": """query($id: String!) {
          entry(entry_id: $id) {
            rcsb_id
            struct { title }
            rcsb_entry_info { resolution_combined experimental_method }
            rcsb_primary_citation { pdbx_database_id_DOI pdbx_database_id_PubMed title journal_abbrev }
            polymer_entities {
              rcsb_id
              entity_poly { pdbx_seq_one_letter_code_can type }
              rcsb_polymer_entity { pdbx_description pdbx_mutation rcsb_source_organism_members }
              rcsb_polymer_entity_container_identifiers { uniprot_ids auth_asym_ids }
              rcsb_polymer_entity_feature { type name description }
            }
            nonpolymer_entities {
              nonpolymer_comp { chem_comp { id name } }
            }
          }
        }""",
        "variables": {"id": pdb},
    }
    req = urllib.request.Request(
        "https://data.rcsb.org/graphql",
        data=json.dumps(q).encode(),
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode())
    except Exception as exc:
        return {"error": str(exc), "pdb": pdb}


def fetch_validation(pdb: str) -> dict:
    urls = [
        f"https://files.rcsb.org/pub/pdb/validation_reports/{pdb[1:3].lower()}/{pdb.lower()}/{pdb.lower()}_validation.xml.gz",
        f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb.lower()}_validation.xml",
        f"https://files.wwpdb.org/pub/pdb/validation_reports/{pdb[1:3].lower()}/{pdb.lower()}/{pdb.lower()}_validation.xml.gz",
    ]
    dest = CIFDIR / f"{pdb}_validation.xml"
    rec = {"pdb": pdb, "validation_fetched": "NO", "urls_tried": urls}
    for url in urls:
        try:
            urllib.request.urlretrieve(url, dest)
            raw = dest.read_bytes()
            if raw[:2] == b"\x1f\x8b":
                import gzip
                text = gzip.decompress(raw).decode("utf-8", "replace")
                dest.write_text(text, encoding="utf-8")
            else:
                text = raw.decode("utf-8", "replace")
            rec["validation_fetched"] = "YES"
            rec["url"] = url
            rec["bytes"] = len(text)
            for key in ("clashscore", "percent-rama-outliers", "percent-rota-outliers", "DCC_R", "DCC_Rfree", "EDS_R", "EDS_Rfree"):
                marker = f'{key}="'
                i = text.find(marker)
                if i >= 0:
                    j = text.find('"', i + len(marker))
                    rec[key] = text[i + len(marker):j]
            break
        except Exception as exc:
            rec.setdefault("errors", []).append(f"{url}: {exc}")
    return rec


def main():
    out = {}
    for pdb in PDBS:
        print("cif", pdb, flush=True)
        out[pdb] = {"cif": extract(pdb), "rcsb": fetch_rcsb(pdb), "val": fetch_validation(pdb)}
        mut = [e.get("pdbx_mutation") for e in out[pdb]["cif"]["entity"]]
        mods = out[pdb]["cif"]["mod_residues"]
        print("  mut", mut, "mods", mods, "val", out[pdb]["val"].get("validation_fetched"), flush=True)
    (QA / "alt_receptor_mutation_ptm_v2.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("wrote", QA / "alt_receptor_mutation_ptm_v2.json")


if __name__ == "__main__":
    main()
