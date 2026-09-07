#!/usr/bin/env python3
"""Independent geometric re-verification of the 14 Tier-1 receptor PDBs.

Does not trust RCSB titles, UniProt domain tables, or any cached CSV in this
repo. Re-fetches identity from RCSB REST and computes real 3D numbers from the
downloaded mmCIF coordinates:

  - entity -> UniProt residue-span alignment (rcsb_polymer_entity_align)
  - primary-citation DOI/PMID
  - depositor-authored `_struct_site_gen` binding-site residues, where present
  - a generic <=4.2 A heavy-atom contact search, where `_struct_site` is absent
  - literal Cys25(Sgamma)-ligand distance for the two covalent-capable
    cathepsin cognates (3XT/I37 on 4X6H; KH0 on 9GJ2)

Produces the numbers reported in analysis/GEOMETRIC_POCKET_VERIFICATION_V1.md.
Requires `gemmi` (coordinate parsing). Network access to files.rcsb.org and
data.rcsb.org is required; nothing is committed to the repo by this script.

Usage:
    python3 verify_geometric_pockets_v1.py --cache-dir /tmp/verify_cache
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

import gemmi

DATA = "https://data.rcsb.org/rest/v1/core"
FILES = "https://files.rcsb.org/download"

# (pdb, expected_uniprot, expected_gene, cognate_ccd, ligand_chain_or_None)
TARGETS = [
    ("4UDW", "P00734", "F2", "N6L", "H"),
    ("2JKH", "P00742", "F10", "BI7", "A"),
    ("6N7A", "P23458", "JAK1", "KEV", None),
    ("3LXP", "P29597", "TYK2", "IZA", "A"),
    ("8BXH", "O60674", "JAK2", "C87", "A"),
    ("9V8H", "P37231", "PPARG", "BRL", "A"),
    ("6LXA", "Q07869", "PPARA", "EPA", "A"),
    ("5U3Q", "Q03181", "PPARD", "7UJ", "A"),
    ("4X6H", "P43235", "CTSK", "3XT", "A"),
    ("9GJ2", "P25774", "CTSS", "KH0", None),
    ("4L23", "P42336", "PIK3CA", "X6K", "A"),
    ("4JT6", "P42345", "MTOR", "X6K", "A"),
    ("4EY7", "P22303", "ACHE", "E20", "A"),
    ("4BDS", "P06276", "BCHE", "THA", "A"),
]


def getj(url: str, cache: Path) -> dict:
    if cache.exists():
        return json.loads(cache.read_text())
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geom-verify/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.load(r)
            cache.write_text(json.dumps(d))
            return d
        except Exception:
            if attempt == 2:
                return {}
            time.sleep(1.0 * (attempt + 1))
    return {}


def fetch_cif(pdb: str, cache_dir: Path) -> Path:
    path = cache_dir / f"{pdb.lower()}.cif"
    if not path.exists():
        urllib.request.urlretrieve(f"{FILES}/{pdb.lower()}.cif", path)
    return path


def identity_check(pdb: str, exp_acc: str, cache_dir: Path) -> None:
    entry = getj(f"{DATA}/entry/{pdb}", cache_dir / f"entry_{pdb}.json")
    cit = entry.get("rcsb_primary_citation") or {}
    print(f"  title: {(entry.get('struct') or {}).get('title', '')}")
    print(f"  DOI={cit.get('pdbx_database_id_DOI')} PMID={cit.get('pdbx_database_id_PubMed')} "
          f"journal={cit.get('rcsb_journal_abbrev')}")
    pids = (entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []
    for pid in pids:
        pe = getj(f"{DATA}/polymer_entity/{pdb}/{pid}", cache_dir / f"pe_{pdb}_{pid}.json")
        ent = pe.get("rcsb_polymer_entity") or {}
        orgs = [s.get("scientific_name") for s in pe.get("rcsb_entity_source_organism") or []]
        for al in pe.get("rcsb_polymer_entity_align") or []:
            if al.get("reference_database_accession") == exp_acc:
                for a in al.get("aligned_regions") or []:
                    lo, ln = a.get("entity_beg_seq_id"), a.get("length")
                    ulo = a.get("ref_beg_seq_id")
                    print(f"    entity {pid} ({ent.get('pdbx_description', '')[:40]!r}, "
                          f"organism={orgs}): entity {lo}-{lo + ln - 1} -> UniProt {ulo}-{ulo + ln - 1}")


def struct_site_residues(cif_path: Path) -> dict[str, list[str]]:
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block()
    gen = block.find(
        ["_struct_site_gen.site_id", "_struct_site_gen.label_comp_id",
         "_struct_site_gen.auth_asym_id", "_struct_site_gen.auth_seq_id"]
    )
    out: dict[str, list[str]] = {}
    for sid, comp, asym, seq in gen:
        out.setdefault(sid, []).append(f"{comp}{seq}({asym})")
    return out


def contact_search(cif_path: Path, lig_name: str, lig_chain: str | None, cutoff: float = 4.2) -> str:
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    lig_atoms = []
    for chain in model:
        if lig_chain and chain.name != lig_chain:
            continue
        for res in chain:
            if res.name == lig_name:
                lig_atoms.extend(a for a in res if not a.element.is_hydrogen)
    if not lig_atoms:
        return f"ligand {lig_name} not found"
    hits = set()
    for chain in model:
        for res in chain:
            if res.name == lig_name or res.is_water():
                continue
            for atom in res:
                if atom.element.is_hydrogen:
                    continue
                if any(atom.pos.dist(la.pos) <= cutoff for la in lig_atoms):
                    hits.add((chain.name, res.name, res.seqid.num))
                    break
    return ", ".join(f"{n}{s}({c})" for c, n, s in sorted(hits, key=lambda x: (x[0], x[2])))


def covalent_distance(cif_path: Path, chain_id: str, cys_seq: int, lig_name: str) -> float | None:
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    cys, lig = None, []
    for chain in model:
        if chain.name != chain_id:
            continue
        for res in chain:
            if res.seqid.num == cys_seq and res.name == "CYS":
                cys = res
            if res.name == lig_name:
                lig.extend(a for a in res if not a.element.is_hydrogen)
    if cys is None or not lig:
        return None
    sg = next((a for a in cys if a.name == "SG"), None)
    if sg is None:
        return None
    return min(sg.pos.dist(a.pos) for a in lig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default="/tmp/verify_cache")
    args = ap.parse_args()
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    for pdb, acc, gene, ccd, lig_chain in TARGETS:
        print("=" * 90)
        print(f"{pdb}  expected {gene} / {acc}")
        identity_check(pdb, acc, cache_dir)
        cif_path = fetch_cif(pdb, cache_dir)
        sites = struct_site_residues(cif_path)
        matched = [(sid, res) for sid, res in sites.items() if any(r.startswith(ccd) for r in res)]
        if matched:
            for sid, res in matched:
                print(f"  depositor struct_site {sid}: {', '.join(res)}")
        else:
            print(f"  contact search (<= 4.2 A): {contact_search(cif_path, ccd, lig_chain)}")

    # Covalent-state distances for the two cathepsin cognates
    print("=" * 90)
    print("Covalent-state distances (Cys25 Sgamma to ligand heavy atoms):")
    x6h = fetch_cif("4X6H", cache_dir)
    print("  4X6H Cys25 vs 3XT:", covalent_distance(x6h, "A", 25, "3XT"), "A")
    print("  4X6H Cys25 vs I37:", covalent_distance(x6h, "A", 25, "I37"), "A")
    gj2 = fetch_cif("9GJ2", cache_dir)
    print("  9GJ2 Cys25(A) vs KH0:", covalent_distance(gj2, "A", 25, "KH0"), "A")
    print("  9GJ2 Cys25(B) vs KH0:", covalent_distance(gj2, "B", 25, "KH0"), "A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
