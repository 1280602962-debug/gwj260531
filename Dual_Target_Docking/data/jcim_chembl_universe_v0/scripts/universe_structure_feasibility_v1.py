#!/usr/bin/env python3
"""PDB/holo feasibility screen for every ChEMBL 37 strict-thick pair.

Does not dock. Does not change Table 2 / K = 4.

H3 matches Dual_Target_Docking/docs/PUBLIC_TARGET_PAIR_SELECTION_REPORT.md:
RCSB experimental entries, resolution ≤ 3.5 Å, ≥1 non-polymer ligand instance.
Suggested pass: ≥5 holo entries on each end.

Ligand overlap uses RCSB chem_comp facets, excluding common crystallization
additives (same spirit as the original report's warning that raw holo counts
include solvent noise).
"""
from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

TABLES = Path(__file__).resolve().parents[1] / "tables"
CACHE = Path(__file__).resolve().parents[1] / "cache" / "rcsb_holo_v1"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

RESO_MAX = 3.5
H3_MIN = 5
SLEEP = 0.12
UA = "DualFourClass-universe-structure-v1"

SOLVENT = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE",
    "PE8", "PE4", "P6G", "P33", "P15", "P4C", "12P", "2PE", "DMS", "DMSO", "EOH",
    "BME", "DTT", "TRS", "EPE", "HEP", "MES", "ACT", "ACY", "ACE", "CIT", "FLC",
    "TLA", "TAR", "MLT", "MLA", "FMT", "NO3", "CL", "NA", "K", "MG", "CA", "ZN",
    "MN", "FE", "CU", "NI", "CO", "CD", "IOD", "BR", "F", "NH4", "NH3", "NAG",
    "BMA", "MAN", "FUC", "GAL", "GLC", "SIA", "NDG", "FUL", "BOG", "LMT", "DGA",
    "OLA", "OLC", "PLM", "STE", "MYR", "DAO", "PEE", "PC1", "PTY", "CLR", "CHL",
    "UMQ", "C8E", "D10", "HEX", "OCT", "D12", "UNL", "UNK", "UNX",
}

# Family-level docking class. Used only after numeric H3.
POCKET = {
    "PIK3CA": "kinase_atp",
    "PIK3CB": "kinase_atp",
    "PIK3CG": "kinase_atp",
    "MTOR": "kinase_atp",
    "JAK1": "kinase_atp",
    "JAK2": "kinase_atp",
    "JAK3": "kinase_atp",
    "TYK2": "kinase_atp",
    "MAPK1": "kinase_atp",
    "EGFR": "kinase_atp",
    "ACHE": "hydrolase_gorge",
    "BCHE": "hydrolase_gorge",
    "MAOA": "flavin_oxidase",
    "MAOB": "flavin_oxidase",
    "F2": "serine_protease",
    "F10": "serine_protease",
    "PRSS1": "serine_protease",
    "CTSK": "cysteine_protease",
    "CTSS": "cysteine_protease",
    "CNR1": "gpcr",
    "CNR2": "gpcr",
    "HCRTR1": "gpcr",
    "HCRTR2": "gpcr",
    "ADORA1": "gpcr",
    "ADORA2A": "gpcr",
    "ADORA3": "gpcr",
    "OPRM1": "gpcr",
    "OPRD1": "gpcr",
    "OPRK1": "gpcr",
    "HTR6": "gpcr",
    "HTR7": "gpcr",
    "S1PR1": "gpcr",
    "S1PR3": "gpcr",
    "NPSR1": "gpcr",
    "TSHR": "gpcr",
    "SLC6A2": "slc_transporter",
    "SLC6A3": "slc_transporter",
    "SLC6A4": "slc_transporter",
    "PPARA": "nuclear_receptor_lbd",
    "PPARD": "nuclear_receptor_lbd",
    "PPARG": "nuclear_receptor_lbd",
    "THRB": "nuclear_receptor_lbd",
    "BRD4": "bromodomain",
    "CREBBP": "hat_and_bromodomain",
    "L3MBTL1": "chromatin_reader",
    "HDAC1": "zn_hydrolase",
    "HDAC2": "zn_hydrolase",
    "HDAC3": "zn_hydrolase",
    "HDAC4": "zn_hydrolase",
    "HDAC6": "zn_hydrolase",
    "HDAC8": "zn_hydrolase",
    "CA1": "zn_lyase",
    "CA2": "zn_lyase",
    "CA4": "zn_lyase",
    "CA9": "zn_lyase",
    "CYP2C9": "heme_p450",
    "CYP2C19": "heme_p450",
    "CYP2D6": "heme_p450",
    "CYP3A4": "heme_p450",
    "ALDH1A1": "qhts_counterscreen",
    "MAPT": "qhts_counterscreen",
    "LMNA": "qhts_counterscreen",
    "SMN1": "qhts_counterscreen",
    "POLB": "qhts_counterscreen",
    "GAA": "qhts_counterscreen",
    "HTT": "qhts_counterscreen",
    "RAB9A": "qhts_counterscreen",
    "TDP1": "qhts_counterscreen",
    "HSD17B10": "qhts_counterscreen",
    "TP53": "qhts_counterscreen",
    "BLM": "qhts_counterscreen",
    "NPC1": "qhts_counterscreen",
}

CONVENTIONAL = {
    "kinase_atp",
    "hydrolase_gorge",
    "flavin_oxidase",
    "serine_protease",
    "cysteine_protease",
    "nuclear_receptor_lbd",
}


def post_search(payload: dict) -> dict:
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        SEARCH,
        data=data,
        headers={"Content-Type": "application/json", "User-Agent": UA},
        method="POST",
    )
    last = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                body = resp.read()
            if not body:
                return {}
            return json.loads(body)
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code == 204:
                return {}
            time.sleep(1.5 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"RCSB search failed: {last}")


def uniprot_nodes(accession: str) -> list[dict]:
    return [
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
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                "operator": "exact_match",
                "value": "UniProt",
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal",
                "value": RESO_MAX,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.deposited_nonpolymer_entity_instance_count",
                "operator": "greater",
                "value": 0,
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
    ]


def holo_query(accession: str) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"holo_{accession}.json"
    if path.exists():
        return json.loads(path.read_text())
    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": uniprot_nodes(accession)},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 25},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }
    rec = post_search(payload)
    path.write_text(json.dumps(rec))
    time.sleep(SLEEP)
    return rec


def ligand_query(accession: str) -> dict:
    path = CACHE / f"lig_{accession}.json"
    if path.exists():
        return json.loads(path.read_text())
    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": uniprot_nodes(accession)},
        "return_type": "non_polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 0},
            "results_content_type": ["experimental"],
            "facets": [
                {
                    "name": "comp_ids",
                    "aggregation_type": "terms",
                    "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                    "max_num_intervals": 3000,
                    "min_interval_population": 1,
                }
            ],
        },
    }
    rec = post_search(payload)
    path.write_text(json.dumps(rec))
    time.sleep(SLEEP)
    return rec


def parse_holo(rec: dict) -> tuple[int, list[str]]:
    n = int(rec.get("total_count") or 0)
    ids = [h.get("identifier", "") for h in rec.get("result_set") or []]
    ids = [i for i in ids if i]
    return n, ids


def parse_ligands(rec: dict) -> set[str]:
    out = set()
    for facet in rec.get("facets") or []:
        if facet.get("name") != "comp_ids":
            continue
        for bucket in facet.get("buckets") or facet.get("terms") or []:
            cid = str(bucket.get("label") or "").upper()
            if cid and cid not in SOLVENT:
                out.add(cid)
    return out


def pocket(gene: str) -> str:
    return POCKET.get(gene, "unclassified")


def pair_verdict(row: dict) -> tuple[str, str]:
    bucket = row["supply_bucket"]
    gene_a, gene_b = row["gene_A"], row["gene_B"]
    pa, pb = pocket(gene_a), pocket(gene_b)
    h3 = int(row["h3_pass"])
    n_ol = int(row["n_shared_non_solvent_ligands"])
    if bucket == "qhts_or_common_counter_screen":
        return "exclude_qhts_counterscreen", "Not a designed dual-target pair; PubChem-style shared libraries."
    if bucket == "cyp_adme_panel":
        return "exclude_cyp_adme", "ADME panel, not a dual-target design problem."
    if bucket == "metal_enzyme":
        return "exclude_metal", "Zn-HDAC or carbonic anhydrase; original H4 excludes as sole primary."
    if not h3:
        return "fail_H3_holo", "Fewer than 5 holo entries (reso ≤ 3.5 Å, non-polymer ligand) on at least one end."
    # PPAR is same NR1C family despite ChEMBL subclass strings.
    if {gene_a, gene_b} <= {"PPARA", "PPARD", "PPARG"}:
        return (
            "include_candidate_same_family_nr",
            "Nuclear-receptor LBD, conventional small-molecule pocket, H3 pass. Same PPAR family (isoform duals exist). Not a cross-class dual.",
        )
    if pa == "bromodomain" or pb == "bromodomain" or pa == "hat_and_bromodomain" or pb == "hat_and_bromodomain":
        return (
            "include_candidate_epigenetic_caveat",
            "H3 pass. CREBBP is HAT+bromodomain; BRD4 is bromodomain. Dual literature exists, but pockets are not kinase-ATP equivalent. Needs domain/pocket freeze before docking.",
        )
    if pa == "gpcr" and pb == "gpcr":
        extra = " Same-ligand co-crystal chem_comp overlap is an upper bound (membrane lipids may remain)."
        if n_ol:
            extra = f" Shared non-solvent chem_comp n={n_ol}."
        return (
            "include_candidate_gpcr_homolog",
            "H3 pass. Class A GPCR 7TM pocket is dockable in principle (membrane protein, lipid/detergent holos)." + extra,
        )
    if pa == "slc_transporter" or pb == "slc_transporter":
        return (
            "include_candidate_transporter_homolog",
            "H3 pass. SLC6 structures exist but are membrane proteins with detergent/lipid ligands; harder than soluble kinases.",
        )
    if pa in CONVENTIONAL and pb in CONVENTIONAL:
        fam = "same-family homolog/isoform" if row["same_class"] == "1" or (
            pa == pb
        ) else "cross-family"
        gold = " Same-end ligand chem_comp overlap detected." if n_ol else " No shared non-solvent chem_comp in this coarse screen (does not prove absence of dual co-crystal)."
        return (
            "include_candidate_conventional",
            f"H3 pass. Conventional small-molecule pockets ({pa}/{pb}); {fam}." + gold,
        )
    return "review", f"H3 pass but pocket class {pa}/{pb} needs manual review."


def main() -> int:
    thick = list(csv.DictReader((TABLES / "universe_pairs_strict_thick_annotated_v1.csv").open()))
    uniprots: dict[str, str] = {}
    for r in thick:
        uniprots[r["uniprot_A"]] = r["gene_A"]
        uniprots[r["uniprot_B"]] = r["gene_B"]
    print(f"thick pairs={len(thick)} unique_uniprots={len(uniprots)}", flush=True)

    holo_map = {}
    lig_map = {}
    for i, (acc, gene) in enumerate(sorted(uniprots.items(), key=lambda x: x[1]), 1):
        print(f"[{i}/{len(uniprots)}] {gene} {acc}", flush=True)
        hrec = holo_query(acc)
        n, examples = parse_holo(hrec)
        holo_map[acc] = {"n": n, "examples": examples}
        lrec = ligand_query(acc)
        lig_map[acc] = parse_ligands(lrec)
        print(f"    holo={n} example={examples[:8]} ligands={len(lig_map[acc])}", flush=True)

    tgt_rows = []
    for acc, gene in sorted(uniprots.items(), key=lambda x: x[1]):
        h = holo_map[acc]
        ligs = sorted(lig_map[acc])
        tgt_rows.append(
            {
                "gene": gene,
                "uniprot": acc,
                "pocket_class": pocket(gene),
                "n_holo_reso_le_3.5": h["n"],
                "h3_end_pass": int(h["n"] >= H3_MIN),
                "example_pdbs": ";".join(h["examples"][:8]),
                "n_non_solvent_chem_comp": len(ligs),
                "example_chem_comp": ";".join(ligs[:12]),
            }
        )
    tgt_fields = list(tgt_rows[0].keys())
    with (TABLES / "universe_target_holo_v1.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=tgt_fields)
        w.writeheader()
        w.writerows(tgt_rows)

    pair_rows = []
    for r in thick:
        ua, ub = r["uniprot_A"], r["uniprot_B"]
        ha, hb = holo_map[ua], holo_map[ub]
        shared = sorted(lig_map[ua] & lig_map[ub])
        rec = dict(r)
        rec.update(
            {
                "pocket_A": pocket(r["gene_A"]),
                "pocket_B": pocket(r["gene_B"]),
                "n_holo_A": ha["n"],
                "n_holo_B": hb["n"],
                "example_pdbs_A": ";".join(ha["examples"][:8]),
                "example_pdbs_B": ";".join(hb["examples"][:8]),
                "h3_A": int(ha["n"] >= H3_MIN),
                "h3_B": int(hb["n"] >= H3_MIN),
                "h3_pass": int(ha["n"] >= H3_MIN and hb["n"] >= H3_MIN),
                "n_shared_non_solvent_ligands": len(shared),
                "shared_chem_comp": ";".join(shared[:20]),
            }
        )
        decision, note = pair_verdict(rec)
        rec["structure_decision"] = decision
        rec["structure_note"] = note
        rec["include_for_fresh_roster"] = int(decision.startswith("include_candidate"))
        pair_rows.append(rec)

    pair_fields = [
        "gene_A",
        "gene_B",
        "uniprot_A",
        "uniprot_B",
        "supply_bucket",
        "same_class",
        "metal_either",
        "min_strict_hardneg",
        "n_both_measured",
        "pocket_A",
        "pocket_B",
        "n_holo_A",
        "n_holo_B",
        "h3_A",
        "h3_B",
        "h3_pass",
        "example_pdbs_A",
        "example_pdbs_B",
        "n_shared_non_solvent_ligands",
        "shared_chem_comp",
        "structure_decision",
        "include_for_fresh_roster",
        "structure_note",
    ]
    pair_rows.sort(key=lambda x: (-x["include_for_fresh_roster"], -int(x["min_strict_hardneg"])))
    with (TABLES / "universe_pairs_structure_feasibility_v1.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pair_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(pair_rows)

    counts = defaultdict(int)
    for rec in pair_rows:
        counts[rec["structure_decision"]] += 1
    print("DECISIONS", dict(counts), flush=True)
    print("include_for_fresh_roster", sum(r["include_for_fresh_roster"] for r in pair_rows), flush=True)
    print("H3 pass among 86", sum(r["h3_pass"] for r in pair_rows), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
