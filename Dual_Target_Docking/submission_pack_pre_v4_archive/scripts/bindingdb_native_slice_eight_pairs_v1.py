#!/usr/bin/env python3
"""BindingDB-native independence remainder for all eight primary pairs.

Same cascade as bindingdb_native_slice_v1.py: native BindingDB/Patent
curation, human wild-type single-chain exact IC50/Ki/Kd, θ = 6.0 four-state
labels, drop shared literature, drop shared InChIKey/ChEMBL ID, drop ECFP4
Tanimoto ≥ 0.70, then n ≥ 20 and ≥ 3 sources per directional class.

Does not dock. Does not package an external evaluation set.
PIK3CA/PIK3CB and MCL1/Bcl-xL are not primary pairs in this remainder.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import sqlite3
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import bindingdb_native_slice_v1 as b  # noqa: E402

TAB = b.TAB
CACHE = b.CACHE
ANALYSIS = b.ANALYSIS
FIG = b.FIG
CHEMBL_DB = Path("/tmp/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
STANDARD_OK = ("IC50", "Ki", "Kd", "EC50", "Potency", "IC50app", "Ki app")

RELEASE = "202609"
ARCHIVES = (
    "BindingDB_BindingDB_Articles_202609_tsv.zip",
    "BindingDB_Patents_202609_tsv.zip",
    "BindingDB_Assays_202609_tsv.zip",
    "BindingDB_rsid_eaids_202609_tsv.zip",
)
EXPECTED_MD5 = {
    "BindingDB_BindingDB_Articles_202609_tsv.zip": "8c41a1fcf8b828d99070f4bca6bd7a86",
    "BindingDB_Patents_202609_tsv.zip": "cc4ef414211e64303108a6dfbfd02906",
    "BindingDB_Assays_202609_tsv.zip": "b80b1ea6064aa25fb44387b10f205ce8",
    "BindingDB_rsid_eaids_202609_tsv.zip": "f2a9567064ad4d893dd3706d5edaf2c8",
}

TARGETS = {
    "EGFR": "P00533",
    "HER2": "P04626",
    "ACHE": "P22303",
    "BCHE": "P06276",
    "PIK3CA": "P42336",
    "MTOR": "P42345",
    "F2": "P00734",
    "F10": "P00742",
    "JAK1": "P23458",
    "TYK2": "P29597",
    "JAK2": "O60674",
    "PPARG": "P37231",
    "PPARA": "Q07869",
    "PPARD": "Q03181",
}
UNIPROT_TO_TARGET = {v: k for k, v in TARGETS.items()}
PAIRS = (
    ("EGFR/HER2", "EGFR", "HER2"),
    ("AChE/BChE", "ACHE", "BCHE"),
    ("PIK3CA/mTOR", "PIK3CA", "MTOR"),
    ("F2/F10", "F2", "F10"),
    ("JAK1/TYK2", "JAK1", "TYK2"),
    ("JAK1/JAK2", "JAK1", "JAK2"),
    ("PPARG/PPARA", "PPARG", "PPARA"),
    ("PPARA/PPARD", "PPARA", "PPARD"),
)
PAIR_ROLES = {
    "EGFR/HER2": "thin_or_primary",
    "AChE/BChE": "primary_external",
    "PIK3CA/mTOR": "primary_external",
    "F2/F10": "primary_external",
    "JAK1/TYK2": "primary_external",
    "JAK1/JAK2": "primary_external",
    "PPARG/PPARA": "primary_external",
    "PPARA/PPARD": "primary_external",
}
CHEMBL_TARGET = {
    "EGFR": "CHEMBL203",
    "HER2": "CHEMBL1824",
    "ACHE": "CHEMBL220",
    "BCHE": "CHEMBL1914",
    "PIK3CA": "CHEMBL4005",
    "MTOR": "CHEMBL2842",
    "F2": "CHEMBL204",
    "F10": "CHEMBL244",
    "JAK1": "CHEMBL2835",
    "TYK2": "CHEMBL3553",
    "JAK2": "CHEMBL2971",
    "PPARG": "CHEMBL235",
    "PPARA": "CHEMBL239",
    "PPARD": "CHEMBL3979",
}
SQLITE_MAP_TARGETS = ("F2", "F10", "JAK1", "TYK2", "JAK2", "PPARG", "PPARA", "PPARD")
ORIGINAL_PAIRS = {"EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"}
PANEL_SMILES = {
    "EGFR/HER2": [
        "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
    ],
    "AChE/BChE": [
        "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv",
    ],
    "PIK3CA/mTOR": [
        "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        "data/pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv",
        "data/jcim_holdout_v0/tables/holdout_panel_HOPM.csv",
    ],
    "F2/F10": [
        "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_F2_F10_v1.csv",
        "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOF2F10_v1.csv",
    ],
    "JAK1/TYK2": [
        "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_TYK2_v1.csv",
        "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOJ1TYK2_v1.csv",
    ],
    "JAK1/JAK2": [
        "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_JAK2_v1.csv",
        "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOJ1J2_v1.csv",
    ],
    "PPARG/PPARA": [
        "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARG_PPARA_v1.csv",
        "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOPGPA_v1.csv",
    ],
    "PPARA/PPARD": [
        "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARA_PPARD_v1.csv",
        "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOPAPD_v1.csv",
    ],
}


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lock_archives() -> list[dict]:
    rows = []
    for name in ARCHIVES:
        dest = CACHE / name
        if not dest.exists() or dest.stat().st_size < 1000:
            raise FileNotFoundError(f"missing BindingDB archive {dest}")
        got = md5_file(dest)
        expected = EXPECTED_MD5[name]
        if got != expected:
            raise RuntimeError(f"md5 mismatch {name}: {got} != {expected}")
        rows.append(
            {
                "file": name,
                "url": f"{b.BDB_BASE}/{name}",
                "release": RELEASE,
                "nbytes": dest.stat().st_size,
                "md5": got,
                "sha256": b.sha256_file(dest),
            }
        )
    b.write_csv(TAB / "bindingdb_archive_lock_v1.csv", rows)
    return rows


def stream_native_records() -> list[dict]:
    extract_path = CACHE / "native_uniprot_extract_eight_pairs_v1.csv.gz"
    if extract_path.exists() and extract_path.stat().st_size > 1000:
        with gzip.open(extract_path, "rt", encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle))
        print(f"reuse native extract {len(records)} from {extract_path}", flush=True)
        return records
    records: list[dict] = []
    wanted = set(TARGETS.values())
    keep_names = (
        "BindingDB Reactant_set_id",
        "Ligand SMILES",
        "Ligand InChI Key",
        "BindingDB MonomerID",
        "Target Name",
        "Target Source Organism According to Curator or DataSource",
        "Ki (nM)",
        "IC50 (nM)",
        "Kd (nM)",
        "EC50 (nM)",
        "Curation/DataSource",
        "Article DOI",
        "PMID",
        "Patent Number",
        "ChEMBL ID of Ligand",
        "Number of Protein Chains in Target (>1 implies a multichain complex)",
        "UniProt (SwissProt) Primary ID of Target Chain 1",
    )
    with gzip.open(extract_path, "wt", encoding="utf-8", newline="") as out:
        writer = None
        for zip_name, source_tag in (
            (ARCHIVES[0], "articles"),
            (ARCHIVES[1], "patents"),
        ):
            zpath = CACHE / zip_name
            print(f"stream {zip_name}", flush=True)
            with zipfile.ZipFile(zpath) as zf:
                inner = zf.namelist()[0]
                with zf.open(inner) as raw:
                    handle = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
                    reader = csv.reader(handle, delimiter="\t")
                    header = next(reader)
                    index = b.colmap(header)
                    if writer is None:
                        writer = csv.DictWriter(out, fieldnames=["source_file", *keep_names], lineterminator="\n")
                        writer.writeheader()
                    n_keep = 0
                    for i, row in enumerate(reader, start=1):
                        uniprot = b.get(row, index, "UniProt (SwissProt) Primary ID of Target Chain 1").strip()
                        if uniprot not in wanted:
                            continue
                        rec = {name: b.get(row, index, name).strip() for name in keep_names}
                        rec["source_file"] = source_tag
                        writer.writerow(rec)
                        records.append(rec)
                        n_keep += 1
                        if i % 500_000 == 0:
                            print(f"  {source_tag} scanned {i:,} kept {n_keep:,}", flush=True)
                    print(f"  {source_tag} kept {n_keep:,}", flush=True)
    print(f"native extract rows {len(records)} -> {extract_path}", flush=True)
    return records


def load_rsid_assays() -> dict[str, str]:
    mapping = {}
    zpath = CACHE / ARCHIVES[3]
    with zipfile.ZipFile(zpath) as zf:
        with zf.open(zf.namelist()[0]) as raw:
            handle = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                mapping[row["REACTANT_SET_ID"].strip()] = (row.get("ENTRYID_ASSAYID") or "").strip()
    return mapping


def development_molecules(pair: str) -> tuple[set[str], dict[str, object], set[str]]:
    keys: set[str] = set()
    fps: dict[str, object] = {}
    chembl_ids: set[str] = set()
    for rel in PANEL_SMILES.get(pair, []):
        path = ROOT / rel
        if not path.exists():
            raise FileNotFoundError(path)
        for row in b.read_csv(path):
            smiles = row.get("smiles") or row.get("canonical_smiles") or ""
            chembl = row.get("molecule_chembl_id") or ""
            if chembl:
                chembl_ids.add(chembl)
            key = (row.get("inchi_key") or "").strip()
            if not key:
                key, smiles = b.largest_fragment_inchikey(smiles)
            else:
                key, smiles = b.largest_fragment_inchikey(smiles, key)
            if key:
                keys.add(key)
                fp = b.fingerprint(smiles)
                if fp is not None:
                    fps[key] = fp
    return keys, fps, chembl_ids


def harvest_sqlite_maps() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    cache_path = CACHE / "chembl_sqlite_maps_eight_pairs_v1.csv.gz"
    maps: dict[str, dict[str, str]] = {name: {} for name in SQLITE_MAP_TARGETS}
    inchikeys: dict[str, dict[str, str]] = {name: {} for name in SQLITE_MAP_TARGETS}
    if cache_path.exists() and cache_path.stat().st_size > 1000:
        with gzip.open(cache_path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                maps[row["target"]][row["chembl_id"]] = row["smiles"]
                if row.get("inchikey"):
                    inchikeys[row["target"]][row["chembl_id"]] = row["inchikey"]
        print("reuse sqlite maps", {k: len(v) for k, v in maps.items()}, flush=True)
        return maps, inchikeys
    if not CHEMBL_DB.exists():
        raise FileNotFoundError(CHEMBL_DB)
    con = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    rows_out = []
    placeholders = ",".join("?" * len(STANDARD_OK))
    sql = f"""
        SELECT md.chembl_id, cs.canonical_smiles, cs.standard_inchi_key
        FROM activities act
        JOIN assays ass ON ass.assay_id = act.assay_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        JOIN compound_structures cs ON cs.molregno = md.molregno
        WHERE ass.tid = ?
          AND act.pchembl_value IS NOT NULL
          AND act.standard_type IN ({placeholders})
        GROUP BY md.chembl_id
    """
    for name in SQLITE_MAP_TARGETS:
        chembl_id = CHEMBL_TARGET[name]
        tid = con.execute(
            "SELECT tid FROM target_dictionary WHERE chembl_id = ? AND target_type = 'SINGLE PROTEIN' AND organism = 'Homo sapiens'",
            [chembl_id],
        ).fetchone()
        if tid is None:
            raise RuntimeError(f"no ChEMBL tid for {name} {chembl_id}")
        n = 0
        for cid, smiles, inchikey in con.execute(sql, [tid[0], *STANDARD_OK]):
            if not cid:
                continue
            maps[name][cid] = smiles or ""
            if inchikey:
                inchikeys[name][cid] = inchikey
            rows_out.append({"target": name, "chembl_id": cid, "smiles": smiles or "", "inchikey": inchikey or ""})
            n += 1
        print(f"sqlite map {name} {chembl_id} tid={tid[0]} n={n}", flush=True)
    con.close()
    with gzip.open(cache_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["target", "chembl_id", "smiles", "inchikey"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)
    return maps, inchikeys


def map_chembl_ids(target: str, sqlite_maps: dict[str, dict[str, str]]) -> set[str]:
    if target in sqlite_maps:
        return set(sqlite_maps[target])
    return b.map_chembl_ids(target)


def map_inchikeys(
    target: str,
    smiles_cache: dict[str, str],
    sqlite_maps: dict[str, dict[str, str]],
    sqlite_inchikeys: dict[str, dict[str, str]],
    inchikey_cache: dict[str, str],
) -> set[str]:
    if target in sqlite_maps:
        keys = set()
        for cid, smiles in sqlite_maps[target].items():
            smiles_cache.setdefault(cid, smiles)
            key = sqlite_inchikeys.get(target, {}).get(cid) or inchikey_cache.get(cid)
            if not key:
                key, _ = b.largest_fragment_inchikey(smiles)
            if key:
                inchikey_cache[cid] = key
                keys.add(key)
        return keys
    return b.map_inchikeys(target, smiles_cache)


def pair_chembl_ids(pair: str) -> list[str]:
    ids = []
    for rel in PANEL_SMILES[pair]:
        for row in b.read_csv(ROOT / rel):
            cid = row.get("molecule_chembl_id") or ""
            if cid:
                ids.append(cid)
    return sorted(set(ids))


def harvest_pair_documents() -> dict[str, tuple[set[str], set[str], set[str]]]:
    cache_path = TAB / "chembl_document_pubmed_doi_eight_pairs_v1.csv"
    by_pair: dict[str, tuple[set[str], set[str], set[str]]] = {}
    if cache_path.exists() and cache_path.stat().st_size > 100:
        grouped: dict[str, list[dict]] = defaultdict(list)
        for row in b.read_csv(cache_path):
            grouped[row["pair"]].append(row)
        for pair in (p[0] for p in PAIRS):
            pmids, dois, patents = set(), set(), set()
            for row in grouped.get(pair, []):
                if row.get("pubmed_id"):
                    pmid = row["pubmed_id"]
                    pmids.add(str(int(pmid)) if str(pmid).isdigit() else str(pmid))
                if row.get("doi"):
                    dois.add(b.normalize_doi(row["doi"]))
                if row.get("patent_id"):
                    patents.add(b.normalize_patent(row["patent_id"]))
            by_pair[pair] = (pmids, dois, patents)
        print("reuse pair documents", {k: tuple(len(x) for x in v) for k, v in by_pair.items()}, flush=True)
        return by_pair
    if not CHEMBL_DB.exists():
        raise FileNotFoundError(CHEMBL_DB)
    con = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    sql = """
        SELECT DISTINCT d.pubmed_id, d.doi, d.patent_id
        FROM activities act
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        JOIN docs d ON d.doc_id = act.doc_id
        WHERE md.chembl_id IN ({placeholders})
    """
    out_rows = []
    for pair, _a, _b in PAIRS:
        ids = pair_chembl_ids(pair)
        pmids, dois, patents = set(), set(), set()
        if ids:
            placeholders = ",".join("?" * len(ids))
            for pubmed_id, doi, patent_id in con.execute(sql.format(placeholders=placeholders), ids):
                pmid = str(pubmed_id).strip() if pubmed_id else ""
                doi_n = b.normalize_doi(doi or "")
                pat = b.normalize_patent(patent_id or "")
                if pmid and pmid not in {"None", "0"}:
                    pmids.add(str(int(pmid)) if pmid.isdigit() else pmid)
                if doi_n:
                    dois.add(doi_n)
                if pat:
                    patents.add(pat)
                out_rows.append(
                    {
                        "pair": pair,
                        "pubmed_id": pmid if pmid not in {"None", "0"} else "",
                        "doi": doi_n,
                        "patent_id": patent_id or "",
                    }
                )
        by_pair[pair] = (pmids, dois, patents)
        print(f"docs {pair} panel_ligands={len(ids)} pmid={len(pmids)} doi={len(dois)} patent={len(patents)}", flush=True)
    con.close()
    b.write_csv(cache_path, out_rows)
    return by_pair


def build_pair_ligands(records: list[dict], assays: dict[str, str]) -> dict[str, list[dict]]:
    b.TARGETS = TARGETS
    b.UNIPROT_TO_TARGET = UNIPROT_TO_TARGET
    b.PAIRS = PAIRS
    return b.build_pair_ligands(records, assays)


def apply_independence(
    pair_ligands,
    smiles_cache,
    sqlite_maps,
    sqlite_inchikeys,
    original_lit,
    pair_lit,
    pubmed_status,
):
    inchikey_cache: dict[str, str] = {}
    flow_rows = []
    panel_rows = []
    summary_rows = []
    orig_pmids, orig_dois, orig_patents = original_lit
    for pair, name_a, name_b in PAIRS:
        ligands = pair_ligands[pair]
        dev_keys, dev_fps, _dev_cids = development_molecules(pair)
        map_keys = map_inchikeys(
            name_a, smiles_cache, sqlite_maps, sqlite_inchikeys, inchikey_cache
        ) | map_inchikeys(name_b, smiles_cache, sqlite_maps, sqlite_inchikeys, inchikey_cache)
        map_cids = map_chembl_ids(name_a, sqlite_maps) | map_chembl_ids(name_b, sqlite_maps)
        structure_keys = set(dev_keys) | set(map_keys)
        dev_fp_list = list(dev_fps.values())
        if pair in ORIGINAL_PAIRS:
            chembl_pmids, chembl_dois, chembl_patents = orig_pmids, orig_dois, orig_patents
        else:
            chembl_pmids, chembl_dois, chembl_patents = pair_lit[pair]

        def count_classes(rows):
            c = Counter(r["class"] for r in rows if r["class"])
            return c["dual"], c["A_only"], c["B_only"], c["neither"]

        both = [r for r in ligands if r["class"]]
        layers = [("native_paired_theta6", both)]
        no_lit = []
        for row in both:
            pmids = set((row["pmids"] or "").split(";")) - {""}
            dois = set((row["dois"] or "").split(";")) - {""}
            patents = set((row["patents"] or "").split(";")) - {""}
            hit = bool(pmids & chembl_pmids) or bool(dois & chembl_dois) or bool(patents & chembl_patents)
            if not hit:
                no_lit.append(row)
        layers.append(("drop_shared_literature", no_lit))
        no_struct = []
        for row in no_lit:
            chembl_hit = any(cid in map_cids for cid in (row.get("chembl_ids") or "").split(";") if cid)
            if row["inchikey"] in structure_keys or chembl_hit:
                continue
            no_struct.append(row)
        layers.append(("drop_shared_structure", no_struct))
        chem_primary = []
        chem_mid = []
        for row in no_struct:
            tan = b.max_tanimoto(b.fingerprint(row["smiles"]), dev_fp_list)
            row = dict(row)
            row["max_ecfp4_to_development"] = "" if tan < 0 else round(tan, 4)
            if tan < 0:
                row["chem_stratum"] = "no_fp"
                chem_primary.append(row)
            elif tan < b.NEAR_TANIMOTO:
                row["chem_stratum"] = "<0.50"
                chem_primary.append(row)
            elif tan < b.PRIMARY_TANIMOTO:
                row["chem_stratum"] = "0.50-0.70"
                chem_mid.append(row)
            else:
                row["chem_stratum"] = ">=0.70"
        layers.append(("drop_neighbors_lt_0.70", chem_primary))
        layers.append(("lt_0.70_plus_0.50_0.70_sensitivity", chem_primary + chem_mid))

        for layer_name, rows in layers:
            d, a, bb, n = count_classes(rows)
            flow_rows.append(
                {
                    "pair": pair,
                    "layer": layer_name,
                    "n_ligands": len(rows),
                    "n_dual": d,
                    "n_A_only": a,
                    "n_B_only": bb,
                    "n_neither": n,
                }
            )

        primary = chem_primary
        stats = b.class_source_stats(primary)
        role, note = b.gate_status(stats, pair)
        summary_rows.append(
            {
                "pair": pair,
                "role": PAIR_ROLES[pair],
                "pubmed_status": pubmed_status,
                "n_dev_inchikeys": len(dev_keys),
                "n_map_inchikeys": len(map_keys),
                "native_paired": len(both),
                "after_literature": len(no_lit),
                "after_structure": len(no_struct),
                "after_ecfp_lt_0.70": len(primary),
                "n_dual": stats["dual"]["n"],
                "n_A_only": stats["A_only"]["n"],
                "n_B_only": stats["B_only"]["n"],
                "n_neither": stats["neither"]["n"],
                "n_sources_dual": stats["dual"]["n_sources"],
                "n_sources_A_only": stats["A_only"]["n_sources"],
                "n_sources_B_only": stats["B_only"]["n_sources"],
                "top_doc_frac_dual": stats["dual"]["top_source_fraction"],
                "top_doc_frac_A_only": stats["A_only"]["top_source_fraction"],
                "top_doc_frac_B_only": stats["B_only"]["top_source_fraction"],
                "gate": role,
                "packaged_as_external_evaluation": 0,
                "note": note + "; no docking in this session",
            }
        )
        for row in primary:
            if not row["class"]:
                continue
            panel_rows.append(
                {
                    "pair": pair,
                    "inchikey": row["inchikey"],
                    "class": row["class"],
                    "strict_class": row["strict_class"],
                    "pA": row["pA"],
                    "pB": row["pB"],
                    "smiles": row["smiles"],
                    "scaffold": row["scaffold"],
                    "source_key": row["source_key"],
                    "pmids": row["pmids"],
                    "dois": row["dois"],
                    "patents": row["patents"],
                    "monomerids": row["monomerids"],
                    "chembl_ids": row["chembl_ids"],
                    "max_ecfp4_to_development": row.get("max_ecfp4_to_development", ""),
                    "chem_stratum": row.get("chem_stratum", ""),
                    "gate": role,
                }
            )
        b.write_csv(
            TAB / f"external_panel_{pair.replace('/', '_').replace('-', '')}_v1.csv",
            [r for r in panel_rows if r["pair"] == pair],
        )
        print(
            pair,
            "paired",
            len(both),
            "lit",
            len(no_lit),
            "struct",
            len(no_struct),
            "ecfp",
            len(primary),
            "D/A/B",
            stats["dual"]["n"],
            stats["A_only"]["n"],
            stats["B_only"]["n"],
            role,
            flush=True,
        )

    b.write_csv(TAB / "external_candidate_flow.csv", flow_rows)
    b.write_csv(TAB / "external_slice_summary_v1.csv", summary_rows)
    b.write_csv(TAB / "external_slice_ligands_v1.csv", panel_rows)
    for row in summary_rows:
        path = TAB / f"external_panel_{row['pair'].replace('/', '_').replace('-', '')}_v1.csv"
        if path.exists() and path.stat().st_size:
            (TAB / f"external_panel_{row['pair'].replace('/', '_').replace('-', '')}_v1.sha256").write_text(
                b.sha256_file(path) + "\n", encoding="utf-8"
            )
    return summary_rows, flow_rows, panel_rows


def plot_flow(flow_rows: list[dict], summary_rows: list[dict]) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    layers = [
        "native_paired_theta6",
        "drop_shared_literature",
        "drop_shared_structure",
        "drop_neighbors_lt_0.70",
    ]
    pairs = [p[0] for p in PAIRS]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.2))
    x = np.arange(len(pairs))
    width = 0.18
    colors = ["#4C78A8", "#F58518", "#54A24B", "#E45756"]
    for i, layer in enumerate(layers):
        mins = []
        for pair in pairs:
            rec = next((r for r in flow_rows if r["pair"] == pair and r["layer"] == layer), None)
            mins.append(0 if rec is None else min(int(rec["n_dual"]), int(rec["n_A_only"]), int(rec["n_B_only"])))
        axes[0].bar(x + (i - 1.5) * width, mins, width, label=layer.replace("_", "\n"), color=colors[i])
    axes[0].axhline(20, color="0.3", ls="--", lw=0.8, label="primary n=20")
    axes[0].axhline(15, color="0.5", ls=":", lw=0.8, label="thin n=15")
    axes[0].set_xticks(x, [p.replace("/", "/\n") for p in pairs], fontsize=6)
    axes[0].set_ylabel("min(dual, A-only, B-only)")
    axes[0].set_title("A  BindingDB-native independence flow")
    axes[0].legend(fontsize=5.5, frameon=False, loc="upper right")

    labels = [r["pair"].replace("/", "/\n") for r in summary_rows]
    dual = [int(r["n_dual"]) for r in summary_rows]
    a_only = [int(r["n_A_only"]) for r in summary_rows]
    b_only = [int(r["n_B_only"]) for r in summary_rows]
    x2 = np.arange(len(labels))
    axes[1].bar(x2 - 0.2, dual, 0.2, label="dual", color="#4C78A8")
    axes[1].bar(x2, a_only, 0.2, label="A-only", color="#F58518")
    axes[1].bar(x2 + 0.2, b_only, 0.2, label="B-only", color="#54A24B")
    axes[1].axhline(20, color="0.3", ls="--", lw=0.8)
    axes[1].set_xticks(x2, labels, fontsize=6)
    axes[1].set_title("B  After literature+structure+ECFP<0.70")
    axes[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "FigS_bindingdb_native_slice_v1.png", dpi=160)
    fig.savefig(FIG / "FigS_bindingdb_native_slice_v1.pdf")
    plt.close(fig)


def write_verdict(summary_rows, flow_rows, pubmed_status, lock_rows) -> None:
    primary = [r for r in summary_rows if r["gate"] == "primary_external"]
    thin = [r for r in summary_rows if r["gate"] == "thin_replication"]
    lines = [
        "# BindingDB-native external-slice freeze (eight primary pairs)",
        "",
        "Same independence cascade as `protocol/external_slice_contract.yaml`.",
        "This session did **not** dock and does **not** package external evaluation.",
        f"BindingDB archive `{RELEASE}` md5-verified. Document lookup: `{pubmed_status}`.",
        "",
        f"Primary-gate pairs: **{len(primary)}**. Thin EGFR-style replications: **{len(thin)}**.",
        "Packaged as external evaluation: **no**.",
        "",
        "| pair | dual/A/B/neither after ECFP<0.70 | gate | sources dual/A/B |",
        "|------|--------------------------------:|------|------------------|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['pair']} | {row['n_dual']}/{row['n_A_only']}/{row['n_B_only']}/{row['n_neither']} | {row['gate']} | {row['n_sources_dual']}/{row['n_sources_A_only']}/{row['n_sources_B_only']} |"
        )
    lines += [
        "",
        "No BindingDB AUROC was computed. Archive sha256 values are in",
        f"`tables/bindingdb_archive_lock_v1.csv` ({len(lock_rows)} files, release {RELEASE}).",
    ]
    (ANALYSIS / "EXTERNAL_SLICE_EIGHT_PAIRS_V1.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    frozen = TAB / "external_slice_summary_202608_contract_v1.csv"
    current = TAB / "external_slice_summary_v1.csv"
    if current.exists() and not frozen.exists():
        frozen.write_text(current.read_text(encoding="utf-8"), encoding="utf-8")
    print("lock BindingDB archives", flush=True)
    lock_rows = lock_archives()
    print("harvest sqlite maps", flush=True)
    sqlite_maps, sqlite_inchikeys = harvest_sqlite_maps()
    print("stream native UniProt rows", flush=True)
    records = stream_native_records()
    assays = load_rsid_assays()
    print("aggregate ligands", flush=True)
    pair_ligands = build_pair_ligands(records, assays)
    smiles_cache = b.load_smiles_cache()
    print("original-pair ChEMBL documents", flush=True)
    orig_pmids, orig_dois, orig_patents, pubmed_status = b.collect_dev_documents()
    print(
        f"original lit pmids={len(orig_pmids)} dois={len(orig_dois)} patents={len(orig_patents)} status={pubmed_status}",
        flush=True,
    )
    print("eight-pair panel/holdout documents from ChEMBL 37 sqlite", flush=True)
    pair_lit = harvest_pair_documents()
    print("independence layers", flush=True)
    summary_rows, flow_rows, _panel = apply_independence(
        pair_ligands,
        smiles_cache,
        sqlite_maps,
        sqlite_inchikeys,
        (orig_pmids, orig_dois, orig_patents),
        pair_lit,
        pubmed_status,
    )
    plot_flow(flow_rows, summary_rows)
    write_verdict(summary_rows, flow_rows, pubmed_status, lock_rows)
    print("wrote eight-pair BindingDB-native slice tables", flush=True)
    for row in summary_rows:
        print(row["pair"], row["n_dual"], row["n_A_only"], row["n_B_only"], row["n_neither"], row["gate"], flush=True)


if __name__ == "__main__":
    main()
