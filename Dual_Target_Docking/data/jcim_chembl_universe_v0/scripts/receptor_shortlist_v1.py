#!/usr/bin/env python3
"""Prospective receptor shortlist for the Tier-1 dockable roster.

Runs the identity audit BEFORE any docking, which is the ordering the 2WXF
failure showed is mandatory: accession + organism + sequence coverage first,
cognate ligand drug-likeness second, cognate redocking only afterwards.

For each target: rank candidate holo PDB entries by resolution, then verify
  - SIFTS UniProt accession == the accession the ChEMBL labels came from
  - source organism is human
  - entity sequence coverage
  - engineered mutation count
  - a drug-like non-polymer cognate exists (MW 150-750, not additive/cofactor)
  - the matched entity is a protein, not a presented peptide/epitope

Domain-only constructs (kinase JH1, thrombin catalytic domain, NR LBD) are
legitimate, so entity length is reported but only peptides are rejected.
Metadata alone cannot tell an ATP-site holo from an allosteric/interface holo:
`intended_site` is declared per target and must be verified by a human before
cognate redocking. mTOR is the worked example - every mTOR entry at <=2.5 A is
an FKBP-rapamycin/FRB complex, not the ATP pocket.

Also reports, per pair, the drug-like chemical components co-crystallised on
both ends (same-ligand pose-gold candidates).

Does not dock. Does not change Table 2 or K = 4.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.request
from pathlib import Path

TABLES = Path(__file__).resolve().parents[1] / "tables"
CACHE = Path(__file__).resolve().parents[1] / "cache" / "rcsb_shortlist_v1"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA = "https://data.rcsb.org/rest/v1/core"

RESO_STAGES = (2.5, 3.6)
TOP_N = 12
MIN_ENTITY_LENGTH = 80  # peptide/epitope cut; domain constructs are legitimate
UNIPROT = "https://rest.uniprot.org/uniprotkb"
MW_MIN, MW_MAX = 150.0, 750.0

# Tier-1 roster: conventional soluble pockets that passed supply + H3 + ligand identity.
INTENDED_SITE = {
    "PIK3CA": "class I PI3K ATP site (p110alpha)",
    "MTOR": "mTOR kinase ATP site (NOT the FKBP-rapamycin/FRB site)",
    "ACHE": "catalytic gorge / CAS-PAS",
    "BCHE": "catalytic gorge",
    "F2": "thrombin S1 pocket (catalytic domain)",
    "F10": "factor Xa S1/S4 pockets (catalytic domain)",
    "JAK1": "JH1 kinase ATP site",
    "TYK2": "JH1 kinase ATP site (declare if JH2 pseudokinase)",
    "JAK2": "JH1 kinase ATP site",
    "PPARG": "ligand-binding domain",
    "PPARA": "ligand-binding domain",
    "PPARD": "ligand-binding domain",
    "CTSK": "S2 pocket, papain-fold active site",
    "CTSS": "S2 pocket, papain-fold active site",
}

TIER1_PAIRS = [
    ("PIK3CA/mTOR", ("PIK3CA", "P42336"), ("MTOR", "P42345"), "already_docked_verified"),
    ("AChE/BChE", ("ACHE", "P22303"), ("BCHE", "P06276"), "already_docked_verified"),
    ("F2/F10", ("F2", "P00734"), ("F10", "P00742"), "new"),
    ("JAK1/TYK2", ("JAK1", "P23458"), ("TYK2", "P29597"), "new"),
    ("JAK1/JAK2", ("JAK1", "P23458"), ("JAK2", "O60674"), "new"),
    ("PPARG/PPARA", ("PPARG", "P37231"), ("PPARA", "Q07869"), "new"),
    ("PPARA/PPARD", ("PPARA", "Q07869"), ("PPARD", "Q03181"), "new"),
    ("CTSK/CTSS", ("CTSK", "P43235"), ("CTSS", "P25774"), "new"),
]

ADDITIVE = {
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "PE8", "PE4",
    "P6G", "P33", "12P", "2PE", "DMS", "EOH", "BME", "DTT", "TRS", "EPE", "MES", "ACT",
    "ACY", "ACE", "CIT", "FLC", "TLA", "TAR", "MLT", "MLA", "FMT", "NO3", "CL", "NA",
    "K", "MG", "CA", "ZN", "MN", "FE", "CU", "NI", "CO", "CD", "IOD", "BR", "F", "NH4",
    "NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "SIA", "NDG", "FUL", "BOG", "LMT", "DGA",
    "OLA", "OLC", "PLM", "STE", "MYR", "DAO", "PEE", "PC1", "PTY", "CLR", "CHL", "Y01",
    "UMQ", "C8E", "D10", "HEX", "OCT", "D12", "UNL", "UNK", "UNX", "IMD", "B3P", "BU3",
    "SCN", "HP6", "MPD", "SIN", "BEN", "AZI", "CAC", "MRD", "PGO", "IPA",
}
COFACTOR = {
    "ADP", "ATP", "AMP", "ANP", "ACP", "GNP", "GTP", "GDP", "GSP", "NAD", "NAI", "NAP",
    "NDP", "NAX", "FAD", "FMN", "SAH", "SAM", "COA", "ACO", "HEM", "HEC", "PLP", "TPP",
    "UDP", "UMP", "CMP", "TTP", "DTP", "5GP", "BGC",
}


def post(payload: dict) -> dict:
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        SEARCH, data=data,
        headers={"Content-Type": "application/json", "User-Agent": "DualFourClass-shortlist"},
        method="POST",
    )
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                body = resp.read()
            return json.loads(body) if body else {}
        except Exception:
            if attempt == 3:
                return {}
            time.sleep(1.5 * (attempt + 1))
    return {}


def getj(url: str) -> dict | None:
    key = url.split("/core/")[-1].replace("/", "_")
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    for attempt in range(4):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "DualFourClass-shortlist", "Accept": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                rec = json.load(resp)
            path.write_text(json.dumps(rec))
            return rec
        except Exception:
            if attempt == 3:
                return None
            time.sleep(1.0 * (attempt + 1))
    return None


def uniprot_length(acc: str) -> int | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"uniprot_len_{acc}.json"
    if path.exists():
        return json.loads(path.read_text()).get("length")
    try:
        req = urllib.request.Request(
            f"{UNIPROT}/{acc}.json", headers={"User-Agent": "DualFourClass-shortlist", "Accept": "application/json"}
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            rec = json.load(resp)
        n = (rec.get("sequence") or {}).get("length")
        path.write_text(json.dumps({"length": n}))
        return n
    except Exception:
        return None


def candidate_entries(acc: str, reso_max: float) -> list[str]:
    nodes = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
            "operator": "exact_match", "value": "UniProt"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal", "value": reso_max}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.deposited_nonpolymer_entity_instance_count",
            "operator": "greater", "value": 0}},
    ]
    rec = post({
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": TOP_N},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    })
    return [h["identifier"] for h in rec.get("result_set") or []]


def chemcomp(cid: str) -> dict:
    rec = getj(f"{DATA}/chemcomp/{cid}") or {}
    c = rec.get("chem_comp") or {}
    return {"comp_id": cid, "mw": c.get("formula_weight"), "name": c.get("name"), "type": c.get("type")}


def is_druglike(lig: dict) -> bool:
    cid = (lig.get("comp_id") or "").upper()
    if cid in ADDITIVE or cid in COFACTOR:
        return False
    if (lig.get("type") or "") != "non-polymer":
        return False
    mw = lig.get("mw")
    return mw is not None and MW_MIN <= float(mw) <= MW_MAX


def audit_entry(pdb: str, expect_acc: str, ref_len: int | None) -> dict:
    entry = getj(f"{DATA}/entry/{pdb}")
    if entry is None:
        return {"pdb": pdb, "verdict": "FETCH_FAILED"}
    info = entry.get("rcsb_entry_info") or {}
    ids = entry.get("rcsb_entry_container_identifiers") or {}
    reso = (info.get("resolution_combined") or [None])[0]
    pids = ids.get("polymer_entity_ids") or []
    accs, orgs, cov, muts, descs = [], set(), None, None, []
    ent_len = None
    for pid in pids:
        pe = getj(f"{DATA}/polymer_entity/{pdb}/{pid}")
        if pe is None:
            continue
        ent = pe.get("rcsb_polymer_entity") or {}
        descs.append((ent.get("pdbx_description") or "")[:60])
        for s in pe.get("rcsb_entity_source_organism") or []:
            if s.get("scientific_name"):
                orgs.add(s["scientific_name"])
        for r in (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "reference_sequence_identifiers"
        ) or []:
            a = r.get("database_accession")
            accs.append(a)
            if a == expect_acc:
                poly = pe.get("entity_poly") or {}
                this_len = poly.get("rcsb_sample_sequence_length")
                if this_len is None:
                    seq = (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
                    this_len = len(seq) or None
                # zymogen-derived proteases split catalytic and light chains across
                # entities that share one accession; keep the longest (catalytic) one
                if this_len is not None and (ent_len is None or this_len > ent_len):
                    ent_len = this_len
                    cov = r.get("entity_sequence_coverage")
                    muts = ent.get("rcsb_mutation_count", ent.get("mutation_count"))
                elif ent_len is None:
                    cov = r.get("entity_sequence_coverage")
                    muts = ent.get("rcsb_mutation_count", ent.get("mutation_count"))
    ligs = []
    for nid in ids.get("non_polymer_entity_ids") or []:
        ne = getj(f"{DATA}/nonpolymer_entity/{pdb}/{nid}")
        if ne is None:
            continue
        cid = (ne.get("pdbx_entity_nonpoly") or {}).get("comp_id")
        if cid:
            ligs.append(chemcomp(cid))
    drug = [l for l in ligs if is_druglike(l)]
    drug.sort(key=lambda l: -(float(l["mw"] or 0)))
    acc_ok = expect_acc in accs
    human = any("homo sapiens" in o.lower() for o in orgs)
    ratio = (ent_len / ref_len) if (ent_len and ref_len) else None
    if not acc_ok:
        verdict = "REJECT_wrong_accession"
    elif not human:
        verdict = "REJECT_non_human"
    elif ent_len is not None and ent_len < MIN_ENTITY_LENGTH:
        verdict = "REJECT_peptide_entity_not_protein"
    elif not drug:
        verdict = "REJECT_no_druglike_cognate"
    elif len(pids) > 1:
        verdict = "CAUTION_complex_site_check_required"
    elif muts not in (0, None):
        verdict = "CAUTION_engineered_mutations"
    else:
        verdict = "OK_candidate"
    return {
        "pdb": pdb,
        "resolution": reso,
        "n_polymer_entities": len(pids),
        "resolved_accessions": ";".join(a for a in accs if a),
        "accession_matches_labels": int(acc_ok),
        "source_organism": "; ".join(sorted(orgs)),
        "is_human": int(human),
        "entity_sequence_coverage": cov,
        "entity_length": ent_len if ent_len else "",
        "uniprot_length": ref_len if ref_len else "",
        "entity_length_ratio": round(ratio, 3) if ratio is not None else "",
        "mutation_count": muts if muts is not None else "",
        "entity_descriptions": " | ".join(descs),
        "site_check_required": 1,
        "site_verified_by_human": 0,
        "cognate_comp_id": drug[0]["comp_id"] if drug else "",
        "cognate_mw": drug[0]["mw"] if drug else "",
        "cognate_name": (drug[0]["name"] or "")[:60] if drug else "",
        "n_druglike_ligands": len(drug),
        "all_druglike_comp_ids": ";".join(l["comp_id"] for l in drug),
        "verdict": verdict,
    }


def main() -> int:
    per_target: dict[str, list[dict]] = {}
    targets: list[tuple[str, str]] = []
    for _, a, b, _ in TIER1_PAIRS:
        for t in (a, b):
            if t not in targets:
                targets.append(t)

    order = {
        "OK_candidate": 0,
        "CAUTION_engineered_mutations": 1,
        "CAUTION_complex_site_check_required": 2,
    }
    rows = []
    for gene, acc in targets:
        ref_len = uniprot_length(acc)
        audits: list[dict] = []
        seen: set[str] = set()
        for stage, reso_max in enumerate(RESO_STAGES, 1):
            print(f"[{gene} {acc}] searching \u2264{reso_max} \u00c5 ...", flush=True)
            for pdb in candidate_entries(acc, reso_max):
                if pdb in seen:
                    continue
                seen.add(pdb)
                rec = audit_entry(pdb, acc, ref_len)
                rec.update({
                    "gene": gene,
                    "uniprot": acc,
                    "search_stage_reso_max": reso_max,
                    "intended_site": INTENDED_SITE.get(gene, ""),
                })
                audits.append(rec)
            if any(order.get(r.get("verdict", ""), 9) <= 2 for r in audits):
                break
        audits.sort(key=lambda r: (order.get(r.get("verdict", ""), 9), float(r.get("resolution") or 99)))
        for i, rec in enumerate(audits, 1):
            rec["rank"] = i
            rows.append(rec)
        per_target[gene] = audits
        ok = [r for r in audits if r["verdict"].startswith(("OK_", "CAUTION_"))]
        print(f"    candidates={len(audits)} OK={len(ok)} best={ok[0]['pdb'] if ok else '-'}", flush=True)

    fields = [
        "gene", "uniprot", "rank", "pdb", "resolution", "verdict",
        "accession_matches_labels", "is_human", "entity_sequence_coverage",
        "entity_length", "uniprot_length", "entity_length_ratio",
        "mutation_count", "n_polymer_entities", "intended_site",
        "site_check_required", "site_verified_by_human",
        "search_stage_reso_max", "resolved_accessions", "source_organism",
        "cognate_comp_id", "cognate_mw", "cognate_name", "n_druglike_ligands",
        "all_druglike_comp_ids", "entity_descriptions",
    ]
    with (TABLES / "tier1_receptor_shortlist_v1.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    pair_rows = []
    for pair, (ga, aa), (gb, ab), status in TIER1_PAIRS:
        oa = [r for r in per_target[ga] if r["verdict"].startswith(("OK_", "CAUTION_"))]
        ob = [r for r in per_target[gb] if r["verdict"].startswith(("OK_", "CAUTION_"))]
        la = {c for r in per_target[ga] for c in (r["all_druglike_comp_ids"] or "").split(";") if c}
        lb = {c for r in per_target[gb] for c in (r["all_druglike_comp_ids"] or "").split(";") if c}
        shared = sorted(la & lb)
        pair_rows.append({
            "pair": pair,
            "status": status,
            "gene_A": ga, "uniprot_A": aa,
            "gene_B": gb, "uniprot_B": ab,
            "n_ok_receptors_A": len(oa),
            "n_ok_receptors_B": len(ob),
            "top_receptor_A": oa[0]["pdb"] if oa else "",
            "top_reso_A": oa[0]["resolution"] if oa else "",
            "top_cognate_A": oa[0]["cognate_comp_id"] if oa else "",
            "top_receptor_B": ob[0]["pdb"] if ob else "",
            "top_reso_B": ob[0]["resolution"] if ob else "",
            "top_cognate_B": ob[0]["cognate_comp_id"] if ob else "",
            "verdict_A": oa[0]["verdict"] if oa else "NONE",
            "verdict_B": ob[0]["verdict"] if ob else "NONE",
            "shared_druglike_comp_ids": ";".join(shared[:15]),
            "n_shared_druglike": len(shared),
            "ready_to_dock": int(bool(oa) and bool(ob)),
        })
    pfields = list(pair_rows[0].keys())
    with (TABLES / "tier1_pair_receptor_plan_v1.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pfields)
        w.writeheader()
        w.writerows(pair_rows)

    print("\npair                 A receptor        B receptor        shared drug-like  ready")
    for r in pair_rows:
        print(
            f"{r['pair']:20s} {r['top_receptor_A']:5s} {str(r['top_reso_A']):>5} {r['top_cognate_A']:6s} "
            f"{r['top_receptor_B']:5s} {str(r['top_reso_B']):>5} {r['top_cognate_B']:6s} "
            f"n={r['n_shared_druglike']:<4} ready={r['ready_to_dock']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
