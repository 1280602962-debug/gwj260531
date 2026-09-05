#!/usr/bin/env python3
"""Collect RCSB/PDBe evidence for Tier-1 proposed receptors.

Does not dock. Does not set pass_fail by itself — writes a structured
evidence JSON/CSV that Layer-2 signoff can cite. Cache is gitignored.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache" / "site_verify_rcsb"
TABLES = ROOT / "tables"
CACHE.mkdir(parents=True, exist_ok=True)

DATA = "https://data.rcsb.org/rest/v1/core"
GRAPH = "https://data.rcsb.org/graphql"

# UniProt domain windows used only as a reading aid (inclusive residue numbers).
DOMAIN_WINDOWS = {
    "P42336": {"label": "PIK3CA p110", "windows": [("p110", 1, 1068)]},
    "P42345": {
        "label": "mTOR",
        "windows": [("FRB", 2015, 2114), ("kinase_ATP", 2181, 2516)],
    },
    "P22303": {"label": "AChE", "windows": [("catalytic", 1, 614)]},
    "P06276": {"label": "BChE", "windows": [("catalytic", 1, 602)]},
    "P00734": {
        "label": "prothrombin",
        "windows": [("light", 328, 363), ("heavy_catalytic", 364, 622)],
    },
    "P00742": {
        "label": "factor X",
        "windows": [("light", 41, 179), ("heavy_catalytic", 235, 488)],
    },
    "P23458": {
        "label": "JAK1",
        "windows": [("JH2", 583, 855), ("JH1", 866, 1154)],
    },
    "P29597": {
        "label": "TYK2",
        "windows": [("JH2", 589, 875), ("JH1", 888, 1176)],
    },
    "O60674": {
        "label": "JAK2",
        "windows": [("JH2", 543, 839), ("JH1", 849, 1124)],
    },
    "P37231": {"label": "PPARG", "windows": [("LBD", 204, 505)]},
    "Q07869": {"label": "PPARA", "windows": [("LBD", 171, 468)]},
    "Q03181": {"label": "PPARD", "windows": [("LBD", 140, 441)]},
    "P43235": {"label": "CTSK", "windows": [("papain", 115, 329)]},
    "P25774": {"label": "CTSS", "windows": [("papain", 115, 331)]},
}

TARGETS = [
    ("PIK3CA/mTOR", "A", "PIK3CA", "P42336", "4L23", "X6K", "class I PI3K ATP site (p110alpha)"),
    ("PIK3CA/mTOR", "B", "MTOR", "P42345", "4JT6", "X6K", "mTOR kinase ATP site (NOT the FKBP-rapamycin/FRB site)"),
    ("AChE/BChE", "A", "ACHE", "P22303", "4EY7", "E20", "catalytic gorge / CAS-PAS"),
    ("AChE/BChE", "B", "BCHE", "P06276", "4BDS", "THA", "catalytic gorge"),
    ("F2/F10", "A", "F2", "P00734", "4UDW", "N6L", "thrombin S1 pocket (catalytic domain)"),
    ("F2/F10", "B", "F10", "P00742", "2JKH", "BI7", "factor Xa S1/S4 pockets (catalytic domain)"),
    ("JAK1/TYK2", "A", "JAK1", "P23458", "6N7A", "KEV", "JH1 kinase ATP site"),
    ("JAK1/TYK2", "B", "TYK2", "P29597", "3LXP", "IZA", "JH1 kinase ATP site (declare if JH2 pseudokinase)"),
    ("JAK1/JAK2", "A", "JAK1", "P23458", "6N7A", "KEV", "JH1 kinase ATP site"),
    ("JAK1/JAK2", "B", "JAK2", "O60674", "8BXH", "C87", "JH1 kinase ATP site"),
    ("PPARG/PPARA", "A", "PPARG", "P37231", "9V8H", "BRL", "ligand-binding domain"),
    ("PPARG/PPARA", "B", "PPARA", "Q07869", "6LXA", "EPA", "ligand-binding domain"),
    ("PPARA/PPARD", "A", "PPARA", "Q07869", "6LXA", "EPA", "ligand-binding domain"),
    ("PPARA/PPARD", "B", "PPARD", "Q03181", "5U3Q", "7UJ", "ligand-binding domain"),
    ("CTSK/CTSS", "A", "CTSK", "P43235", "4X6H", "3XT", "S2 pocket, papain-fold active site"),
    ("CTSK/CTSS", "B", "CTSS", "P25774", "9GJ2", "KH0", "S2 pocket, papain-fold active site"),
]


def get_json(url: str, payload: bytes | None = None) -> dict:
    req = urllib.request.Request(
        url,
        data=payload,
        headers={"User-Agent": "DualFourClass-site-verify/1.0", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def cached_get(name: str, url: str, payload: bytes | None = None) -> dict:
    path = CACHE / name
    if path.exists() and path.stat().st_size > 20:
        return json.loads(path.read_text())
    last = None
    for i in range(4):
        try:
            data = get_json(url, payload)
            path.write_text(json.dumps(data))
            time.sleep(0.15)
            return data
        except urllib.error.HTTPError as e:
            last = e
            if e.code == 404:
                path.write_text(json.dumps({"error": 404, "url": url}))
                return {"error": 404, "url": url}
            time.sleep(1.5 * (i + 1))
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def graphql_entry(pdb: str) -> dict:
    q = {
        "query": """
        query ($id: String!) {
          entry(entry_id: $id) {
            rcsb_id
            struct { title }
            rcsb_entry_info { resolution_combined experimental_method }
            polymer_entities {
              rcsb_id
              entity_poly { pdbx_seq_one_letter_code_can type }
              rcsb_polymer_entity { pdbx_description }
              rcsb_polymer_entity_container_identifiers { uniprot_ids auth_asym_ids }
              rcsb_entity_source_organism { ncbi_scientific_name }
              rcsb_polymer_entity_align {
                reference_database_accession
                aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
              }
            }
            nonpolymer_entities {
              rcsb_id
              rcsb_nonpolymer_entity { pdbx_description }
              nonpolymer_comp {
                chem_comp { id name formula_weight type }
              }
            }
          }
        }
        """,
        "variables": {"id": pdb},
    }
    return cached_get(f"gql_{pdb}.json", GRAPH, json.dumps(q).encode())


def ligand_binding(pdb: str, ccd: str) -> dict:
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb.lower()}"
    return cached_get(f"pdbe_lig_{pdb}.json", url)


def pdbe_uniprot(pdb: str) -> dict:
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb.lower()}"
    return cached_get(f"pdbe_up_{pdb}.json", url)


def pdbe_molecules(pdb: str) -> dict:
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb.lower()}"
    return cached_get(f"pdbe_mol_{pdb}.json", url)


def rcsb_entry(pdb: str) -> dict:
    return cached_get(f"rcsb_entry_{pdb}.json", f"{DATA}/entry/{pdb}")


def rcsb_chemcomp(ccd: str) -> dict:
    return cached_get(f"ccd_{ccd}.json", f"{DATA}/chemcomp/{ccd}")


def overlap_windows(ref_start: int, ref_end: int, accession: str) -> list[str]:
    wins = DOMAIN_WINDOWS.get(accession, {}).get("windows", [])
    hits = []
    for name, a, b in wins:
        lo, hi = max(ref_start, a), min(ref_end, b)
        if hi >= lo:
            hits.append(f"{name}:{lo}-{hi}(ov={hi-lo+1})")
    return hits


def summarise(pair, end, gene, acc, pdb, ccd, intended) -> dict:
    entry = rcsb_entry(pdb)
    title = ((entry.get("struct") or {}).get("title")) or ""
    reso_list = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    reso = reso_list[0] if reso_list else None
    method = (entry.get("exptl") or [{}])[0].get("method") if entry.get("exptl") else None
    if not method:
        method = (entry.get("rcsb_entry_info") or {}).get("experimental_method")

    mols = pdbe_molecules(pdb)
    mol_list = (mols.get(pdb.lower()) or []) if isinstance(mols, dict) else []
    polymers = []
    ligands = []
    for m in mol_list:
        mt = (m.get("molecule_type") or "").lower()
        rec = {
            "entity": m.get("entity_id"),
            "uniprots": [
                x.get("accession")
                for x in (m.get("cross_references") or [])
                if str(x.get("accession", "")).startswith(("P", "Q", "O", "A"))
            ] or ([m.get("uniprot_accession")] if m.get("uniprot_accession") else []),
            "orgs": [m.get("source_organism") or m.get("organism_scientific_name") or ""],
            "length": m.get("length") or (len(m.get("sequence") or "")),
            "desc": (m.get("molecule_name") or [m.get("synonym") or ""])[0]
            if isinstance(m.get("molecule_name"), list)
            else (m.get("molecule_name") or m.get("chem_comp_name") or ""),
            "chains": m.get("in_chains") or [],
            "type": m.get("molecule_type"),
        }
        # PDBe sometimes puts accession at top-level
        if m.get("accession"):
            rec["uniprots"] = list({*(rec["uniprots"] or []), m["accession"]})
        if "polypeptide" in mt:
            polymers.append(rec)
        elif mt in {"bound", "ligand", "carbohydrate polymer"} or m.get("chem_comp_ids"):
            for cid in m.get("chem_comp_ids") or []:
                ligands.append({"id": cid, "name": rec["desc"], "mw": m.get("weight"), "type": mt})

    up = pdbe_uniprot(pdb)
    up_map = ((up.get(pdb.lower()) or {}).get("UniProt")) or {}
    aligns = []
    all_uniprots = set()
    for uacc, block in up_map.items():
        all_uniprots.add(uacc)
        for mapping in block.get("mappings") or []:
            start = mapping.get("uniprot_start") or mapping.get("unp_start")
            stop = mapping.get("uniprot_end") or mapping.get("unp_end")
            if start is None:
                continue
            rec_al = {
                "accession": uacc,
                "ref": f"{start}-{stop}",
                "len": int(stop) - int(start) + 1,
                "chain": mapping.get("chain_id"),
                "windows": overlap_windows(int(start), int(stop), acc) if uacc == acc else [],
                "identity": mapping.get("identity"),
            }
            aligns.append(rec_al)

    matching_aligns = [a for a in aligns if a["accession"] == acc]
    window_hits = sorted({w for a in matching_aligns for w in a["windows"]})
    longest = max((a["len"] for a in matching_aligns), default=0)
    matching_orgs = []
    for p in polymers:
        # attach uniprot from mapping via chain
        for a in aligns:
            if a["chain"] in (p["chains"] or []):
                p["uniprots"] = list({*(p["uniprots"] or []), a["accession"]})
        if acc in (p["uniprots"] or []):
            matching_orgs.extend(p["orgs"])

    # organism from PDBe molecule source
    human = False
    if matching_orgs:
        human = any("sapiens" in (o or "").lower() for o in matching_orgs)
    # fallback: PDBe uniprot name
    if acc in up_map:
        name = (up_map[acc].get("name") or "") + " " + (up_map[acc].get("identifier") or "")
        if "HUMAN" in name.upper():
            human = True

    ccd_info = rcsb_chemcomp(ccd) if ccd else {}
    cognate = {
        "id": ccd,
        "name": ((ccd_info.get("chem_comp") or {}).get("name")) or ccd,
        "mw": (ccd_info.get("chem_comp") or {}).get("formula_weight"),
        "type": (ccd_info.get("chem_comp") or {}).get("type"),
        "formula": (ccd_info.get("chem_comp") or {}).get("formula"),
    }
    pdbe_lig = ligand_binding(pdb, ccd)
    lig_rows = (pdbe_lig.get(pdb.lower()) or []) if isinstance(pdbe_lig, dict) else []
    cognate_in_entry = any(
        (row.get("chem_comp_id") or row.get("chem_comp_unid") or "").upper() == ccd.upper()
        for row in lig_rows
    ) or any(x.get("id") == ccd for x in ligands)

    partners = []
    for uacc, block in up_map.items():
        if uacc == acc:
            continue
        partners.append(f"{block.get('name') or uacc}|{uacc}")

    return {
        "pair": pair,
        "end": end,
        "gene": gene,
        "uniprot": acc,
        "pdb": pdb,
        "ccd": ccd,
        "intended_site": intended,
        "title": title,
        "resolution": reso,
        "method": method,
        "n_polymers": len(polymers),
        "matching_entities": matching_aligns,
        "longest_matching_aa": longest,
        "label_acc_present": acc in all_uniprots,
        "human_on_matching": human,
        "domain_window_hits": ";".join(window_hits),
        "uniprot_span": ";".join(f"{a['chain']}:{a['ref']}" for a in matching_aligns),
        "all_uniprots": sorted(all_uniprots),
        "partners": partners,
        "cognate": cognate,
        "cognate_in_entry": cognate_in_entry,
        "all_ligands": sorted({x.get("chem_comp_id") for x in lig_rows if x.get("chem_comp_id")}),
        "rcsb_3d": f"https://www.rcsb.org/3d-view/{pdb}",
        "rcsb_entry": f"https://www.rcsb.org/structure/{pdb}",
        "rcsb_ligand": f"https://www.rcsb.org/ligand/{ccd}",
    }


def main() -> None:
    rows = []
    bundle = []
    seen_pdb = {}
    for args in TARGETS:
        pdb = args[4]
        if pdb in seen_pdb:
            rec = dict(seen_pdb[pdb])
            rec["pair"], rec["end"], rec["gene"] = args[0], args[1], args[2]
            rec["intended_site"] = args[6]
        else:
            rec = summarise(*args)
            seen_pdb[pdb] = rec
        rows.append(rec)
        bundle.append(rec)
        print(
            f"{rec['pair']:16} {rec['end']} {rec['pdb']} {rec['ccd']:4} "
            f"acc={rec['label_acc_present']} human={rec['human_on_matching']} "
            f"len={rec['longest_matching_aa']} win={rec['domain_window_hits'] or '-'} "
            f"lig={rec['cognate']['name'][:40] if rec['cognate'] else 'MISSING'}"
        )

    (CACHE / "evidence_bundle_v1.json").write_text(json.dumps(bundle, indent=2))
    out_csv = TABLES / "site_verification_evidence_v1.csv"
    fields = [
        "pair", "end", "gene", "uniprot", "pdb", "ccd", "intended_site",
        "title", "resolution", "method", "label_acc_present", "human_on_matching",
        "longest_matching_aa", "uniprot_span", "domain_window_hits", "all_uniprots", "partners",
        "cognate_name", "cognate_mw", "cognate_in_entry", "all_ligands", "rcsb_3d", "rcsb_ligand",
    ]
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rec in rows:
            cog = rec.get("cognate") or {}
            w.writerow({
                "pair": rec["pair"],
                "end": rec["end"],
                "gene": rec["gene"],
                "uniprot": rec["uniprot"],
                "pdb": rec["pdb"],
                "ccd": rec["ccd"],
                "intended_site": rec["intended_site"],
                "title": rec["title"],
                "resolution": rec["resolution"],
                "method": rec["method"],
                "label_acc_present": int(rec["label_acc_present"]),
                "human_on_matching": int(rec["human_on_matching"]),
                "longest_matching_aa": rec["longest_matching_aa"],
                "uniprot_span": rec.get("uniprot_span", ""),
                "domain_window_hits": rec["domain_window_hits"],
                "all_uniprots": ";".join(rec["all_uniprots"]),
                "partners": " || ".join(rec["partners"]),
                "cognate_name": cog.get("name", ""),
                "cognate_mw": cog.get("mw", ""),
                "cognate_in_entry": int(rec.get("cognate_in_entry") or 0),
                "all_ligands": ";".join(x for x in rec["all_ligands"] if x),
                "rcsb_3d": rec["rcsb_3d"],
                "rcsb_ligand": rec["rcsb_ligand"],
            })
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
