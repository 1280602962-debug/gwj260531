#!/usr/bin/env python3
"""Audit the protein identity of every receptor actually used for docking.

For each frozen receptor PDB: RCSB entry title, resolution, polymer entity
description, SIFTS UniProt accession, source organism, and the cognate ligand
identity/MW. Compares the resolved accession against the accession the
DualFourClass label set was harvested from.

Catches the 3T8M-class failure mode (right family, wrong protein/species)
on structures that were actually used, not only on rejected candidates.

Does not dock. Does not rewrite Table 2.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.request
from pathlib import Path

TABLES = Path(__file__).resolve().parents[1] / "tables"
DATA = "https://data.rcsb.org/rest/v1/core"

# (pair, end label, expected gene, expected human UniProt used for ChEMBL labels,
#  pdb, cognate comp_id, where the receptor file lives)
RECEPTORS = [
    ("PIK3CA/mTOR", "A", "PIK3CA", "P42336", "4L23", "X6K", "pik3ca_mtor_panel48_rdkit_v0"),
    ("PIK3CA/mTOR", "B", "MTOR", "P42345", "4JT6", "X6K", "pik3ca_mtor_panel48_rdkit_v0"),
    ("AChE/BChE", "A", "ACHE", "P22303", "4EY7", "E20", "ache_bche_panel_v0"),
    ("AChE/BChE", "B", "BCHE", "P06276", "4BDS", "THA", "ache_bche_panel_v0"),
    ("PIK3CA/PIK3CB", "A", "PIK3CA", "P42336", "4L23", "X6K", "pik3ca_pik3cb_panel_v0"),
    ("PIK3CA/PIK3CB", "B", "PIK3CB", "P42338", "2WXF", "039", "pik3ca_pik3cb_panel_v0"),
    ("EGFR/HER2", "A", "EGFR", "P00533", "3POZ", "03P", "egfr_her2_panel120_v0"),
    ("EGFR/HER2", "B", "ERBB2", "P04626", "3RCD", "03P", "egfr_her2_panel120_v0"),
    # Receptor-realization alternates (Table S30)
    ("PIK3CA/mTOR alt", "A", "PIK3CA", "P42336", "4JPS", "", "jcim_structure_robust_v0"),
    ("PIK3CA/mTOR alt", "A", "PIK3CA", "P42336", "5DXT", "", "jcim_structure_robust_v0"),
    ("PIK3CA/mTOR alt", "B", "MTOR", "P42345", "4JSX", "", "jcim_structure_robust_v0"),
    # Rejected candidates, kept for the record
    ("rejected cand", "B", "PIK3CB", "P42338", "2Y3A", "GD9", "not_used"),
    ("rejected cand", "B", "PIK3CB", "P42338", "4BFR", "J82", "not_used"),
    ("rejected cand", "B", "PIK3CB", "P42338", "3T8M", "3T8", "not_used"),
    ("stress test", "A", "MCL1", "Q07820", "3WIY", "", "mcl1_bclxl_panel_v0"),
    ("stress test", "B", "BCL2L1", "Q07817", "3WIZ", "", "mcl1_bclxl_panel_v0"),
]


def getj(url: str, tries: int = 4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "DualFourClass-receptor-audit", "Accept": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception:
            if attempt == tries - 1:
                return None
            time.sleep(1.0 * (attempt + 1))
    return None


def audit_pdb(pdb: str) -> dict:
    entry = getj(f"{DATA}/entry/{pdb}")
    if entry is None:
        return {"pdb": pdb, "error": "entry_fetch_failed"}
    info = entry.get("rcsb_entry_info") or {}
    reso = (info.get("resolution_combined") or [None])[0]
    title = (entry.get("struct") or {}).get("title") or ""
    ids = entry.get("rcsb_entry_container_identifiers") or {}
    ents = []
    for pid in ids.get("polymer_entity_ids") or []:
        pe = getj(f"{DATA}/polymer_entity/{pdb}/{pid}")
        if pe is None:
            continue
        desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or ""
        refs = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "reference_sequence_identifiers"
        ) or []
        orgs = sorted({(s.get("scientific_name") or "").strip() for s in (pe.get("rcsb_entity_source_organism") or [])})
        for r in refs:
            ents.append(
                {
                    "entity": pid,
                    "desc": desc,
                    "accession": r.get("database_accession"),
                    "coverage": r.get("entity_sequence_coverage"),
                    "organism": "; ".join(o for o in orgs if o),
                }
            )
        if not refs:
            ents.append({"entity": pid, "desc": desc, "accession": "", "coverage": "", "organism": "; ".join(orgs)})
    ligs = []
    for nid in ids.get("non_polymer_entity_ids") or []:
        ne = getj(f"{DATA}/nonpolymer_entity/{pdb}/{nid}")
        if ne is None:
            continue
        cid = (ne.get("pdbx_entity_nonpoly") or {}).get("comp_id")
        if not cid:
            continue
        cc = getj(f"{DATA}/chemcomp/{cid}") or {}
        c = cc.get("chem_comp") or {}
        ligs.append({"comp_id": cid, "mw": c.get("formula_weight"), "name": c.get("name"), "type": c.get("type")})
    return {"pdb": pdb, "reso": reso, "title": title, "entities": ents, "ligands": ligs}


def main() -> int:
    cache: dict[str, dict] = {}
    rows = []
    for pair, end, gene, expect_acc, pdb, cognate, pack in RECEPTORS:
        if pdb not in cache:
            print(f"fetching {pdb} ...", flush=True)
            cache[pdb] = audit_pdb(pdb)
        rec = cache[pdb]
        if rec.get("error"):
            rows.append({"pair": pair, "end": end, "expected_gene": gene, "pdb": pdb, "verdict": "FETCH_FAILED"})
            continue
        ents = rec["entities"]
        accs = [e["accession"] for e in ents if e["accession"]]
        orgs = sorted({e["organism"] for e in ents if e["organism"]})
        primary = max(ents, key=lambda e: (e["coverage"] or 0)) if ents else {}
        cog = next((l for l in rec["ligands"] if l["comp_id"] == cognate), None)
        acc_match = expect_acc in accs
        human = any("homo sapiens" in (o or "").lower() for o in orgs)
        if acc_match and human:
            verdict = "OK_expected_human_protein"
        elif acc_match and not human:
            verdict = "SPECIES_MISMATCH_accession_ok"
        elif not acc_match and human:
            verdict = "ACCESSION_MISMATCH_human"
        else:
            verdict = "WRONG_PROTEIN_AND_SPECIES"
        rows.append(
            {
                "pair": pair,
                "end": end,
                "expected_gene": gene,
                "expected_uniprot_for_labels": expect_acc,
                "pdb": pdb,
                "resolution": rec["reso"],
                "entry_title": rec["title"][:150],
                "resolved_accessions": ";".join(accs),
                "primary_entity_desc": (primary.get("desc") or "")[:100],
                "primary_coverage": primary.get("coverage"),
                "source_organism": "; ".join(orgs),
                "accession_matches_labels": int(acc_match),
                "is_human": int(human),
                "cognate_comp_id": cognate,
                "cognate_mw": (cog or {}).get("mw", ""),
                "cognate_name": ((cog or {}).get("name") or "")[:70],
                "cognate_is_small_molecule": int(
                    bool(cog) and (cog.get("type") == "non-polymer") and 100.0 <= float(cog.get("mw") or 0) <= 900.0
                ),
                "pack": pack,
                "verdict": verdict,
            }
        )
    fields = list(rows[0].keys())
    out = TABLES / "receptor_identity_audit_v1.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}\n")
    for r in rows:
        print(
            f"{r['pair']:18s} {r['end']} exp={r['expected_gene']:8s} {r['pdb']:5s} "
            f"acc={r.get('resolved_accessions',''):20s} org={r.get('source_organism','')[:22]:22s} {r['verdict']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
