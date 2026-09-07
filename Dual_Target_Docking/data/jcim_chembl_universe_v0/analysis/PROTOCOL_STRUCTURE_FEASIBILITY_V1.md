# Protocol — PDB/holo screen of all 86 strict-thick pairs

**Status:** executed after the supply census. Does not dock. Does not change K = 4.  
**Script:** `scripts/universe_structure_feasibility_v1.py`

## Frame

Every pair in `universe_pairs_strict_thick_annotated_v1.csv` (N = 86) is assigned a structure decision. Pairs are not silently dropped.

## H3 (identical to the 2026-07-23 public-pair report)

- Source: RCSB Search API v2, `return_type=entry`, `results_content_type=experimental`
- UniProt accession exact match on `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers`
- `rcsb_entry_info.resolution_combined ≤ 3.5`
- `deposited_nonpolymer_entity_instance_count > 0`
- Pass: ≥ 5 such entries **on each end**

Human UniProt PDB cross-references (UniProt REST) are used only to explain H3 failures (e.g. PIK3CB P42338 has zero xrefs). They do not replace the RCSB count.

## Ligand overlap (diagnostic)

RCSB `non_polymer_entity` facet on `nonpolymer_comp_id`, minus a solvent/ion/detergent stoplist. This is an **upper bound**. It does not prove a dual-end same-ligand co-crystal (the original report already warned that raw holo counts include additives). PI-103 / X6K on 4JT6 (3.60 Å) is outside this cut.

## H4 pocket class

Assigned from gene family, not from docking:

- conventional: kinase ATP, hydrolase gorge, flavin oxidase, serine/cysteine protease, nuclear-receptor LBD
- gpcr / slc_transporter / bromodomain / HAT: flagged, not auto-excluded
- metal HDAC/CA and CYP: excluded as in the original H4 / ADME rules
- qHTS hub proteins: excluded as not designed dual-target pairs

## Forbidden

- Do not treat the 19 H3-pass pairs as a new DualFourClass-Bench.
- Do not dock them in this paper.
- Do not replace K = 4 or Table 2.
- Do not count mouse 2Y3A as a human PIK3CB H3 pass without writing the ortholog exception.
