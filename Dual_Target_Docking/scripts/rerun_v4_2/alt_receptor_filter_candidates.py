#!/usr/bin/env python3
"""Filter RCSB audit to experimental holo alternatives. No docking/AUROC."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

RAW = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/alt_receptor_rcsb_raw.json")
OUT = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/alt_receptor_filtered.json")

# UniProt domain ranges used only to require same functional domain coverage.
# Sources: UniProt domain annotations (kinase, catalytic, NR LBD).
DOMAIN = {
    "EGFR": (712, 979),
    "HER2": (720, 987),
    "JAK1": (868, 1154),
    "JAK2": (839, 1124),
    "TYK2": (889, 1178),
    "PIK3CA": (696, 1068),
    "mTOR": (2182, 2516),
    "AChE": (32, 574),
    "BChE": (29, 602),
    "F2": (364, 622),
    "F10": (235, 467),
    "PPARG": (207, 477),
    "PPARA": (230, 468),
    "PPARD": (171, 441),
}

# Common crystallization / buffer / ion components — not binding-site ligands.
ADDITIVE = {
    "HOH", "DOD", "NA", "K", "MG", "CA", "ZN", "CL", "BR", "IOD", "SO4", "PO4",
    "GOL", "EDO", "PEG", "PGE", "PG4", "PE4", "1PE", "P6G", "PE5", "PE8",
    "DMS", "DMSO", "BME", "DTT", "ACT", "ACY", "ACE", "FMT", "EOH", "IPA",
    "MPD", "MRD", "BOG", "NAG", "MAN", "BMA", "FUC", "GAL", "GLC", "SIA",
    "NH2", "NH4", "CO3", "CAC", "CIT", "TLA", "TAR", "MLA", "MAL", "SUC",
    "TRS", "EPE", "HEP", "MES", "IMD", "OH", "OXY", "NO3", "SCN", "AZI",
    "CD", "NI", "CO", "MN", "FE", "CU", "HG", "PB", "SR", "BA",
}

# Extra partners that change background vs V4.2 primary extract policy.
PRIMARY_PARTNER = {
    "EGFR": "none_kinase_only",
    "HER2": "none_kinase_only",
    "JAK1": "none_JH1_only",
    "JAK2": "none_JH1_only",
    "TYK2": "none_JH1_only",
    "PIK3CA": "p85_niSH2_present_in_crystal_REMOVED_in_V42_extract",
    "mTOR": "mLST8_present_in_crystal_REMOVED_in_V42_extract",
    "AChE": "homodimer_AB_retained",
    "BChE": "monomer",
    "F2": "hirudin_I_RETAINED_plus_light_L",
    "F10": "light_chain_L_RETAINED",
    "PPARG": "PG08NL_peptide_B_RETAINED",
    "PPARA": "none_LBD_only",
    "PPARD": "none_LBD_copyB",
}


def coverage(aligns, uniprot, lo, hi) -> float:
    if not aligns:
        return 0.0
    covered = set()
    for a in aligns:
        if (a.get("reference_database_accession") or "").split("_")[0] != uniprot:
            # accession may be bare UniProt
            if a.get("reference_database_accession") != uniprot:
                continue
        for reg in a.get("aligned_regions") or []:
            beg = int(reg["ref_beg_seq_id"])
            n = int(reg["length"])
            for i in range(beg, beg + n):
                if lo <= i <= hi:
                    covered.add(i)
    return len(covered) / float(hi - lo + 1)


def mutations(entity) -> list[str]:
    out = []
    for feat in entity.get("rcsb_polymer_entity_feature") or []:
        name = (feat.get("type") or "") + " " + (feat.get("name") or "")
        if "mutat" in name.lower() or "variant" in name.lower() or "engineered" in name.lower():
            for pos in feat.get("feature_positions") or []:
                out.append(f"{feat.get('name')}:{pos.get('beg_seq_id')}-{pos.get('end_seq_id')}")
            if not feat.get("feature_positions"):
                out.append(feat.get("name") or feat.get("type") or "")
    return out


def ligands(entry) -> list[dict]:
    rows = []
    for ne in entry.get("nonpolymer_entities") or []:
        cc = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        cid = cc.get("id") or ""
        if not cid or cid in ADDITIVE:
            continue
        rows.append({
            "ccd": cid,
            "name": cc.get("name") or "",
            "type": cc.get("type") or "",
            "auth_asym_ids": (ne.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [],
        })
    return rows


def main() -> None:
    blob = json.loads(RAW.read_text())
    out = []
    for tgt in blob["targets"]:
        lo, hi = DOMAIN[tgt["target"]]
        for e in tgt["entries"] or []:
            if not e:
                continue
            pdb = e.get("rcsb_id")
            method = ",".join((x or {}).get("method") or "" for x in (e.get("exptl") or []))
            reso = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or []
            reso_v = reso[0] if reso else None
            # matching polymer entities
            entities = []
            for pe in e.get("polymer_entities") or []:
                ups = [u.get("rcsb_id") for u in (pe.get("uniprots") or [])]
                if tgt["uniprot"] not in ups:
                    continue
                cov = coverage(pe.get("rcsb_polymer_entity_align") or [], tgt["uniprot"], lo, hi)
                entities.append({
                    "entity": pe.get("rcsb_id"),
                    "desc": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
                    "auth_asym_ids": (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [],
                    "length": (pe.get("entity_poly") or {}).get("rcsb_sample_sequence_length"),
                    "coverage": round(cov, 3),
                    "mutations": mutations(pe),
                    "uniprots": ups,
                })
            if not entities:
                continue
            best_cov = max(x["coverage"] for x in entities)
            ligs = ligands(e)
            cit = e.get("rcsb_primary_citation") or {}
            partners = []
            for pe in e.get("polymer_entities") or []:
                ups = [u.get("rcsb_id") for u in (pe.get("uniprots") or [])]
                if tgt["uniprot"] in ups:
                    continue
                partners.append({
                    "desc": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
                    "uniprots": ups,
                    "auth_asym_ids": (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [],
                })
            rec = {
                "target": tgt["target"],
                "uniprot": tgt["uniprot"],
                "primary_pdb": tgt["primary_pdb"],
                "family": tgt["family"],
                "domain": tgt["domain"],
                "pdb_id": pdb,
                "is_primary": pdb == tgt["primary_pdb"],
                "title": (e.get("struct") or {}).get("title"),
                "method": method,
                "resolution": reso_v,
                "deposit_date": (e.get("rcsb_accession_info") or {}).get("deposit_date"),
                "doi": cit.get("pdbx_database_id_DOI"),
                "pmid": cit.get("pdbx_database_id_PubMed"),
                "citation_title": cit.get("title"),
                "domain_coverage": best_cov,
                "entities": entities,
                "ligands": ligs,
                "n_nonadditive_ligands": len(ligs),
                "partners": partners,
                "primary_partner_background": PRIMARY_PARTNER[tgt["target"]],
            }
            rec["pass_experimental"] = "X-RAY" in method.upper() or "ELECTRON" in method.upper() or "NEUTRON" in method.upper()
            rec["pass_domain"] = best_cov >= 0.70
            rec["pass_holo"] = len(ligs) >= 1
            rec["pass_hard"] = rec["pass_experimental"] and rec["pass_domain"] and rec["pass_holo"]
            out.append(rec)
    OUT.write_text(json.dumps(out), encoding="utf-8")
    # summary
    by = defaultdict(lambda: {"n": 0, "hard": 0, "holo": 0})
    for r in out:
        by[r["target"]]["n"] += 1
        by[r["target"]]["holo"] += int(r["pass_holo"])
        by[r["target"]]["hard"] += int(r["pass_hard"] and not r["is_primary"])
    print("target n_human_exp holo hard_alt")
    for t, c in by.items():
        print(t, c["n"], c["holo"], c["hard"])
    print("wrote", OUT, "n", len(out))


if __name__ == "__main__":
    main()
