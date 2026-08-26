#!/usr/bin/env python3
"""BindingDB-native independent-slice supply audit (Phase 1; no docking).

Rules were frozen in protocol/external_slice_contract.yaml and
docs/JCIM_NO_WETLAB_DEEP_PLAN_V2.md before independent class counts were
used to change filters. This script does not dock, does not compute AUROC,
and does not package a slice as external evaluation.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, inchi
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
PROTOCOL = ROOT / "data" / "jcim_novelty_v0" / "protocol"
CACHE = ROOT / "data" / "jcim_novelty_v0" / "cache" / "bindingdb_native_v1"
FIG = ROOT / "figures" / "jcim_article"
for path in (TAB, ANALYSIS, PROTOCOL, CACHE, FIG):
    path.mkdir(parents=True, exist_ok=True)

RELEASE = "202608"
BDB_BASE = "https://www.bindingdb.org/rwd/bind/downloads"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_POLY = "https://data.rcsb.org/rest/v1/core/polymer_entity"
UA = "DualFourClass-Bench/1.0 (BindingDB-native slice; academic; no docking)"
THETA = 6.0
SEED = 20260729
ECFP_BITS = 2048
ECFP_RADIUS = 2
PRIMARY_TANIMOTO = 0.70
NEAR_TANIMOTO = 0.50

ARCHIVES = (
    "BindingDB_BindingDB_Articles_202608_tsv.zip",
    "BindingDB_Patents_202608_tsv.zip",
    "BindingDB_Assays_202608_tsv.zip",
    "BindingDB_rsid_eaids_202608_tsv.zip",
)

TARGETS = {
    "EGFR": "P00533",
    "HER2": "P04626",
    "ACHE": "P22303",
    "BCHE": "P06276",
    "PIK3CA": "P42336",
    "PIK3CB": "P42338",
    "MTOR": "P42345",
    "MCL1": "Q07820",
    "BCL2L1": "Q07817",
}
UNIPROT_TO_TARGET = {v: k for k, v in TARGETS.items()}
PAIRS = (
    ("EGFR/HER2", "EGFR", "HER2"),
    ("AChE/BChE", "ACHE", "BCHE"),
    ("PIK3CA/PIK3CB", "PIK3CA", "PIK3CB"),
    ("PIK3CA/mTOR", "PIK3CA", "MTOR"),
    ("MCL1/Bcl-xL", "MCL1", "BCL2L1"),
)
PAIR_ROLES = {
    "EGFR/HER2": "thin_or_primary",
    "AChE/BChE": "primary_external",
    "PIK3CA/PIK3CB": "backup",
    "PIK3CA/mTOR": "primary_external",
    "MCL1/Bcl-xL": "ppi_bh3_extension",
}
ENDPOINTS = (("Ki (nM)", "Ki"), ("IC50 (nM)", "IC50"), ("Kd (nM)", "Kd"))
KEEP_CURATION = ("bindingdb", "patent")
DROP_CURATION = ("chembl", "pubchem", "pdsp", "csar")
HUMAN = {"homo sapiens", "human"}
MUT_WORDS = ("mutant", "mutation", "mutated", "fusion", "chimera", "truncate")
MUT_RE = re.compile(r"\[[^\]]*[A-Z]\d+[A-Z]|del\d|ins[A-Z]", re.I)
NUM_RE = re.compile(r"([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)")
PMID_RE = re.compile(r"\d{4,9}")

PANEL_SMILES = {
    "EGFR/HER2": [
        "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
    ],
    "AChE/BChE": [
        "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        "data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv",
    ],
    "PIK3CA/PIK3CB": [
        "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "data/jcim_holdout_v0/tables/holdout_panel_HOAP.csv",
    ],
    "PIK3CA/mTOR": [
        "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        "data/pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv",
        "data/jcim_holdout_v0/tables/holdout_panel_HOPM.csv",
    ],
    "MCL1/Bcl-xL": [],
}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def request_json(url: str, cache_key: str, timeout: int = 30, attempts: int = 3):
    cache = CACHE / "http" / f"{cache_key}.json"
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists() and cache.stat().st_size > 2:
        return json.loads(cache.read_text(encoding="utf-8"))
    last_error = None
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            cache.write_text(json.dumps(payload), encoding="utf-8")
            return payload
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            time.sleep(1.2 * (attempt + 1))
    print(f"skip {cache_key}: {last_error}", flush=True)
    return None


def download_archives() -> list[dict]:
    lock_rows = []
    for name in ARCHIVES:
        dest = CACHE / name
        url = f"{BDB_BASE}/{name}"
        md5_url = url.replace(".zip", ".md5")
        if not dest.exists() or dest.stat().st_size < 100:
            print(f"download {name}", flush=True)
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=600) as response:
                dest.write_bytes(response.read())
        md5_path = CACHE / name.replace(".zip", ".md5")
        if not md5_path.exists():
            try:
                req = urllib.request.Request(md5_url, headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=60) as response:
                    md5_path.write_bytes(response.read())
            except Exception:
                pass
        expected_md5 = md5_path.read_text(encoding="utf-8").strip().split()[0] if md5_path.exists() else ""
        import hashlib as _hl

        got_md5 = _hl.md5(dest.read_bytes()).hexdigest()
        if expected_md5 and got_md5 != expected_md5:
            raise RuntimeError(f"md5 mismatch {name}: {got_md5} != {expected_md5}")
        lock_rows.append(
            {
                "file": name,
                "url": url,
                "release": RELEASE,
                "nbytes": dest.stat().st_size,
                "md5": got_md5,
                "sha256": sha256_file(dest),
            }
        )
    write_csv(TAB / "bindingdb_archive_lock_v1.csv", lock_rows)
    return lock_rows


def parse_affinity(raw: str):
    s = (raw or "").strip().replace(",", "")
    if not s:
        return "", None
    qual = "="
    rest = s
    if s[:2] in (">=", "<=", "~="):
        qual, rest = s[:2], s[2:].strip()
    elif s[0] in "<>~":
        qual, rest = s[0], s[1:].strip()
    match = NUM_RE.search(rest)
    if not match:
        return qual, None
    return qual, float(match.group(1))


def nm_to_p(nm: float | None) -> float | None:
    if nm is None or nm <= 0:
        return None
    return 9.0 - math.log10(nm)


def is_human(organism: str) -> bool:
    return (organism or "").strip().lower() in HUMAN


def is_mutant(target_name: str) -> bool:
    name = target_name or ""
    lower = name.lower()
    if "wild-type" in lower or "wild type" in lower:
        return False
    if any(word in lower for word in MUT_WORDS):
        return True
    return bool(MUT_RE.search(name))


def native_curation(source: str) -> bool:
    text = (source or "").lower()
    if any(token in text for token in DROP_CURATION):
        return False
    return any(token in text for token in KEEP_CURATION)


def largest_fragment_inchikey(smiles: str, fallback: str = "") -> tuple[str, str]:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return fallback, smiles or ""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return fallback, smiles or ""
    best = max(frags, key=lambda fragment: fragment.GetNumHeavyAtoms())
    canon = Chem.MolToSmiles(best) or (smiles or "")
    try:
        key = inchi.MolToInchiKey(best) or fallback
    except Exception:
        key = fallback
    return key, canon


def fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    best = max(frags, key=lambda fragment: fragment.GetNumHeavyAtoms())
    return AllChem.GetMorganFingerprintAsBitVect(best, ECFP_RADIUS, nBits=ECFP_BITS)


def murcko(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or ""
    except Exception:
        return ""


def colmap(header: list[str]) -> dict[str, int]:
    return {name: idx for idx, name in enumerate(header)}


def get(row: list[str], index: dict[str, int], name: str) -> str:
    idx = index.get(name)
    if idx is None or idx >= len(row):
        return ""
    return row[idx]


def stream_native_records() -> list[dict]:
    extract_path = CACHE / "native_uniprot_extract_v1.csv.gz"
    if extract_path.exists() and extract_path.stat().st_size > 1000:
        records = []
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
            ("BindingDB_BindingDB_Articles_202608_tsv.zip", "articles"),
            ("BindingDB_Patents_202608_tsv.zip", "patents"),
        ):
            zpath = CACHE / zip_name
            with zipfile.ZipFile(zpath) as zf:
                inner = zf.namelist()[0]
                with zf.open(inner) as raw:
                    handle = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
                    reader = csv.reader(handle, delimiter="\t")
                    header = next(reader)
                    index = colmap(header)
                    if writer is None:
                        writer = csv.DictWriter(
                            out,
                            fieldnames=["source_file", *keep_names],
                            lineterminator="\n",
                        )
                        writer.writeheader()
                    for row in reader:
                        uniprot = get(row, index, "UniProt (SwissProt) Primary ID of Target Chain 1").strip()
                        if uniprot not in wanted:
                            continue
                        rec = {name: get(row, index, name).strip() for name in keep_names}
                        rec["source_file"] = source_tag
                        writer.writerow(rec)
                        records.append(rec)
    print(f"native extract rows {len(records)} -> {extract_path}", flush=True)
    return records


def load_rsid_assays() -> dict[str, str]:
    mapping = {}
    zpath = CACHE / "BindingDB_rsid_eaids_202608_tsv.zip"
    with zipfile.ZipFile(zpath) as zf:
        with zf.open(zf.namelist()[0]) as raw:
            handle = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                mapping[row["REACTANT_SET_ID"].strip()] = (row.get("ENTRYID_ASSAYID") or "").strip()
    return mapping


def four_state(p_a: float | None, p_b: float | None, theta: float = THETA) -> str:
    if p_a is None or p_b is None:
        return ""
    a_on, b_on = p_a >= theta, p_b >= theta
    if a_on and b_on:
        return "dual"
    if a_on and not b_on:
        return "A_only"
    if b_on and not a_on:
        return "B_only"
    return "neither"


def strict_state(p_a: float | None, p_b: float | None) -> str:
    if p_a is None or p_b is None:
        return ""
    if p_a >= 6.5 and p_b >= 6.5:
        return "dual"
    if p_a >= 6.5 and p_b <= 5.5:
        return "A_only"
    if p_b >= 6.5 and p_a <= 5.5:
        return "B_only"
    if p_a <= 5.5 and p_b <= 5.5:
        return "neither"
    return "gray"


def aggregate_target(values: dict[str, list[float]]) -> tuple[float | None, str, str]:
    endpoint_medians = []
    discordant = []
    used = []
    for endpoint, nums in values.items():
        if not nums:
            continue
        lo, hi = min(nums), max(nums)
        if lo < THETA <= hi:
            discordant.append(endpoint)
            continue
        med = sorted(nums)[len(nums) // 2] if len(nums) % 2 else 0.5 * (
            sorted(nums)[len(nums) // 2 - 1] + sorted(nums)[len(nums) // 2]
        )
        endpoint_medians.append(med)
        used.append(f"{endpoint}:{med:.4f}")
    if discordant and not endpoint_medians:
        return None, "discordant", ";".join(discordant)
    if not endpoint_medians:
        return None, "no_exact", ""
    endpoint_medians.sort()
    n = len(endpoint_medians)
    p = endpoint_medians[n // 2] if n % 2 else 0.5 * (endpoint_medians[n // 2 - 1] + endpoint_medians[n // 2])
    return p, "ok", ",".join(used)


def development_molecules(pair: str) -> tuple[set[str], dict[str, object], set[str]]:
    keys: set[str] = set()
    fps: dict[str, object] = {}
    chembl_ids: set[str] = set()
    for rel in PANEL_SMILES.get(pair, []):
        path = ROOT / rel
        if not path.exists():
            continue
        for row in read_csv(path):
            smiles = row.get("smiles") or ""
            chembl = row.get("molecule_chembl_id") or ""
            if chembl:
                chembl_ids.add(chembl)
            key = (row.get("inchi_key") or "").strip()
            if not key:
                key, smiles = largest_fragment_inchikey(smiles)
            else:
                key, smiles = largest_fragment_inchikey(smiles, key)
            if key:
                keys.add(key)
                fp = fingerprint(smiles)
                if fp is not None:
                    fps[key] = fp
    return keys, fps, chembl_ids


MOL_JSON = {
    "EGFR": "mols_EGFR.json",
    "HER2": "mols_HER2.json",
    "ACHE": "mols_ACHE.json",
    "BCHE": "mols_BCHE.json",
    "PIK3CA": "mols_PIK3CA.json",
    "PIK3CB": "mols_PIK3CB.json",
    "MTOR": "mols_MTOR.json",
    "MCL1": "mols_MCL1.json",
    "BCL2L1": "mols_BCL2L1_BclxL.json",
}


def map_chembl_ids(target: str) -> set[str]:
    path = ROOT / "data" / "public_pair_selection" / MOL_JSON[target]
    if not path.exists():
        return set()
    return set(json.loads(path.read_text(encoding="utf-8")))


def map_inchikeys(target: str, smiles_by_chembl: dict[str, str]) -> set[str]:
    path = ROOT / "data" / "public_pair_selection" / MOL_JSON[target]
    if not path.exists():
        return set()
    mols = json.loads(path.read_text(encoding="utf-8"))
    keys = set()
    for chembl_id in mols:
        smiles = smiles_by_chembl.get(chembl_id, "")
        if not smiles:
            continue
        key, _ = largest_fragment_inchikey(smiles)
        if key:
            keys.add(key)
    return keys


def load_smiles_cache() -> dict[str, str]:
    path = TAB / "chembl_smiles_cache_theta6_v1.csv"
    if not path.exists():
        return {}
    return {row["chembl_id"]: row["smiles"] for row in read_csv(path) if row.get("chembl_id")}


def fetch_smiles(chembl_ids: list[str], cache: dict[str, str]) -> dict[str, str]:
    missing = [cid for cid in chembl_ids if cid and cid not in cache]
    for i in range(0, len(missing), 20):
        batch = missing[i : i + 20]
        query = urllib.parse.urlencode(
            {"molecule_chembl_id__in": ",".join(batch), "limit": max(len(batch), 5)}
        )
        payload = request_json(f"{CHEMBL}/molecule.json?{query}", f"mol_{hashlib.md5(','.join(batch).encode()).hexdigest()[:12]}")
        hits = {row["molecule_chembl_id"]: row for row in (payload or {}).get("molecules") or []}
        for cid in batch:
            rec = hits.get(cid) or {}
            struct = rec.get("molecule_structures") or {}
            smiles = struct.get("canonical_smiles") or ""
            if smiles:
                cache[cid] = smiles
        time.sleep(0.15)
    return cache


def collect_dev_documents() -> tuple[set[str], set[str], set[str], str]:
    docs = set()
    for name in (
        "high_confidence_activity_audit_v1.csv",
        "document_blocked_ligand_groups_v1.csv",
        "document_year_lookup_v1.csv",
    ):
        path = TAB / name
        if not path.exists():
            continue
        for row in read_csv(path):
            if row.get("document_chembl_id"):
                docs.add(row["document_chembl_id"])
            raw = row.get("documents") or ""
            for item in raw.split(";"):
                if item.startswith("CHEMBL"):
                    docs.add(item.strip())
    known_pmid = {}
    lookup = TAB / "document_pubmed_lookup_v1.csv"
    if lookup.exists():
        for row in read_csv(lookup):
            if row.get("document_chembl_id"):
                known_pmid[row["document_chembl_id"]] = row.get("pubmed_id") or ""
    docs = sorted(docs)
    pmids: set[str] = set()
    dois: set[str] = set()
    patents: set[str] = set()
    fetched_rows = []
    have = {}
    http_cache = CACHE / "http"
    if http_cache.exists():
        for path in http_cache.glob("doc*.json"):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except Exception:
                continue
            for rec in payload.get("documents") or []:
                if rec.get("document_chembl_id"):
                    have[rec["document_chembl_id"]] = rec
    missing = [doc for doc in docs if doc not in have]
    status = "ok"
    print(f"document cache hits {len(docs) and len([d for d in docs if d in have])}; missing {len(missing)}/{len(docs)}", flush=True)
    if missing:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        def one_doc(doc: str):
            return doc, request_json(
                f"{CHEMBL}/document.json?document_chembl_id={doc}",
                f"doc1_{doc}",
                timeout=8,
                attempts=1,
            )

        resolved_missing = 0
        with ThreadPoolExecutor(max_workers=6) as pool:
            futures = [pool.submit(one_doc, doc) for doc in missing]
            for i, fut in enumerate(as_completed(futures), start=1):
                doc, payload = fut.result()
                recs = (payload or {}).get("documents") or []
                if recs:
                    have[doc] = recs[0]
                    resolved_missing += 1
                if i % 40 == 0 or i == len(futures):
                    print(f"  missing fetch {i}/{len(futures)} extra_resolved={resolved_missing}", flush=True)
        if resolved_missing < len(missing):
            status = "chembl_document_api_partial"
    for doc in docs:
        rec = have.get(doc) or {}
        pmid = rec.get("pubmed_id") or known_pmid.get(doc) or ""
        doi = (rec.get("doi") or "").strip().lower()
        patent = (rec.get("patent_id") or "").strip().upper()
        fetched_rows.append(
            {
                "document_chembl_id": doc,
                "pubmed_id": pmid or "",
                "doi": doi,
                "patent_id": patent,
                "year": rec.get("year") or "",
            }
        )
        if pmid:
            pmids.add(str(int(pmid)) if str(pmid).isdigit() else str(pmid))
        if doi:
            dois.add(doi.replace("https://doi.org/", "").strip())
        if patent:
            patents.add(re.sub(r"[^A-Z0-9]", "", patent.upper()))
    if fetched_rows:
        write_csv(TAB / "chembl_document_pubmed_doi_v1.csv", fetched_rows)
    n_resolved = sum(1 for r in fetched_rows if r["pubmed_id"] or r["doi"] or r["patent_id"])
    if n_resolved < len(docs):
        status = "chembl_document_api_partial"
    print(f"resolved {n_resolved}/{len(docs)} development documents", flush=True)
    return pmids, dois, patents, status


def normalize_patent(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", (value or "").upper())


def normalize_doi(value: str) -> str:
    text = (value or "").strip().lower()
    text = text.replace("https://doi.org/", "").replace("http://doi.org/", "")
    return text


def source_key(pmid: str, doi: str, patent: str) -> str:
    if pmid and PMID_RE.fullmatch(str(pmid).strip()):
        return f"pmid:{int(pmid)}"
    if doi:
        return f"doi:{normalize_doi(doi)}"
    if patent:
        return f"pat:{normalize_patent(patent)}"
    return ""


def max_tanimoto(fp, library: list) -> float:
    if fp is None or not library:
        return -1.0
    return float(max(DataStructs.BulkTanimotoSimilarity(fp, library)))


def class_source_stats(ligands: list[dict]) -> dict[str, dict]:
    out = {}
    for cls in ("dual", "A_only", "B_only", "neither"):
        members = [row for row in ligands if row["class"] == cls]
        sources = Counter()
        for row in members:
            key = row.get("source_key") or f"unknown:{row['inchikey'][:8]}"
            sources[key] += 1
        n = len(members)
        top_frac = (max(sources.values()) / n) if n else 0.0
        out[cls] = {
            "n": n,
            "n_sources": len(sources),
            "top_source_fraction": round(top_frac, 4),
        }
    return out


def gate_status(stats: dict[str, dict], pair: str) -> tuple[str, str]:
    dual, a_only, b_only = stats["dual"], stats["A_only"], stats["B_only"]
    n_ok = dual["n"] >= 20 and a_only["n"] >= 20 and b_only["n"] >= 20
    src_ok = all(stats[c]["n_sources"] >= 3 for c in ("dual", "A_only", "B_only"))
    frac_ok = all(stats[c]["top_source_fraction"] <= 0.50 for c in ("dual", "A_only", "B_only") if stats[c]["n"])
    thin = dual["n"] >= 15 and a_only["n"] >= 15 and b_only["n"] >= 15
    if n_ok and src_ok and frac_ok:
        return "primary_external", "passes A3 primary gate"
    if pair == "EGFR/HER2" and thin:
        return "thin_replication", "15-19 directional n; not a cross-pair primary"
    if thin and not n_ok:
        return "below_primary", "directional n 15-19 but not a designated thin pair"
    return "insufficient", "does not meet primary or thin gate"


def build_pair_ligands(records: list[dict], assays: dict[str, str]) -> dict[str, list[dict]]:
    by_target = defaultdict(list)
    provenance_rows = []
    for rec in records:
        uniprot = rec["UniProt (SwissProt) Primary ID of Target Chain 1"]
        target = UNIPROT_TO_TARGET[uniprot]
        curation = rec["Curation/DataSource"]
        reasons = []
        if not native_curation(curation):
            reasons.append("not_native_curation")
        if not is_human(rec["Target Source Organism According to Curator or DataSource"]):
            reasons.append("not_human")
        if is_mutant(rec["Target Name"]):
            reasons.append("mutant_or_fusion")
        try:
            n_chains = int(float(rec["Number of Protein Chains in Target (>1 implies a multichain complex)"] or "1"))
        except ValueError:
            n_chains = 1
        if n_chains != 1:
            reasons.append("multi_chain")
        key, smiles = largest_fragment_inchikey(rec["Ligand SMILES"], rec["Ligand InChI Key"])
        exact_values = defaultdict(list)
        relations = []
        for col, endpoint in ENDPOINTS:
            qual, nm = parse_affinity(rec.get(col, ""))
            if nm is None:
                continue
            relations.append(f"{endpoint}:{qual}")
            if qual == "=":
                pval = nm_to_p(nm)
                if pval is not None:
                    exact_values[endpoint].append(pval)
            else:
                reasons.append(f"non_equal_{endpoint}")
        if not exact_values:
            reasons.append("no_exact_IC50_Ki_Kd")
        keep_primary = (
            native_curation(curation)
            and is_human(rec["Target Source Organism According to Curator or DataSource"])
            and not is_mutant(rec["Target Name"])
            and n_chains == 1
            and bool(exact_values)
            and bool(key)
        )
        row = {
            **rec,
            "target": target,
            "inchikey": key,
            "canon_smiles": smiles,
            "keep_primary": int(keep_primary),
            "exclusion": ";".join(reasons),
            "entry_assay": assays.get(rec["BindingDB Reactant_set_id"], ""),
            "exact_endpoints": ",".join(sorted(exact_values)),
        }
        provenance_rows.append(
            {
                "rsid": rec["BindingDB Reactant_set_id"],
                "source_file": rec["source_file"],
                "target": target,
                "uniprot": uniprot,
                "monomerid": rec["BindingDB MonomerID"],
                "inchikey": key,
                "chembl_ligand": rec["ChEMBL ID of Ligand"],
                "curation": curation,
                "organism": rec["Target Source Organism According to Curator or DataSource"],
                "target_name": rec["Target Name"],
                "pmid": rec["PMID"],
                "doi": rec["Article DOI"],
                "patent": rec["Patent Number"],
                "entry_assay": row["entry_assay"],
                "exact_endpoints": row["exact_endpoints"],
                "keep_primary": int(keep_primary),
                "exclusion": ";".join(reasons),
            }
        )
        if keep_primary:
            by_target[target].append(
                {
                    "inchikey": key,
                    "smiles": smiles,
                    "monomerid": rec["BindingDB MonomerID"],
                    "chembl_ligand": rec["ChEMBL ID of Ligand"],
                    "pmid": "".join(PMID_RE.findall(rec["PMID"] or "")[:1]),
                    "doi": normalize_doi(rec["Article DOI"]),
                    "patent": rec["Patent Number"],
                    "source_key": source_key(
                        "".join(PMID_RE.findall(rec["PMID"] or "")[:1]),
                        rec["Article DOI"],
                        rec["Patent Number"],
                    ),
                    "exact_values": exact_values,
                    "curation": curation,
                    "rsid": rec["BindingDB Reactant_set_id"],
                }
            )
    # persist a compact provenance (may be large; gzip)
    prov_path = TAB / "bindingdb_native_record_provenance_v1.csv.gz"
    with gzip.open(prov_path, "wt", encoding="utf-8", newline="") as handle:
        if provenance_rows:
            writer = csv.DictWriter(handle, fieldnames=list(provenance_rows[0].keys()), lineterminator="\n")
            writer.writeheader()
            writer.writerows(provenance_rows)
    print(f"provenance rows {len(provenance_rows)} kept_primary {sum(r['keep_primary'] for r in provenance_rows)}", flush=True)

    pair_ligands = {}
    for pair, name_a, name_b in PAIRS:
        merged = defaultdict(lambda: {"A": defaultdict(list), "B": defaultdict(list), "meta": {}})
        for side, tname in (("A", name_a), ("B", name_b)):
            for rec in by_target.get(tname, []):
                item = merged[rec["inchikey"]]
                for endpoint, vals in rec["exact_values"].items():
                    item[side][endpoint].extend(vals)
                meta = item["meta"]
                meta.setdefault("smiles", rec["smiles"])
                meta.setdefault("monomerids", set()).add(rec["monomerid"])
                meta.setdefault("pmids", set())
                meta.setdefault("dois", set())
                meta.setdefault("patents", set())
                meta.setdefault("source_keys", set())
                meta.setdefault("chembl", set())
                if rec["pmid"]:
                    meta["pmids"].add(rec["pmid"])
                if rec["doi"]:
                    meta["dois"].add(rec["doi"])
                if rec["patent"]:
                    meta["patents"].add(normalize_patent(rec["patent"]))
                if rec["source_key"]:
                    meta["source_keys"].add(rec["source_key"])
                if rec["chembl_ligand"]:
                    meta["chembl"].add(rec["chembl_ligand"])
        ligands = []
        for key, payload in merged.items():
            p_a, status_a, detail_a = aggregate_target(payload["A"])
            p_b, status_b, detail_b = aggregate_target(payload["B"])
            cls = four_state(p_a, p_b)
            meta = payload["meta"]
            source_keys = sorted(meta.get("source_keys") or [])
            ligands.append(
                {
                    "pair": pair,
                    "inchikey": key,
                    "smiles": meta.get("smiles", ""),
                    "class": cls,
                    "strict_class": strict_state(p_a, p_b),
                    "pA": "" if p_a is None else round(p_a, 4),
                    "pB": "" if p_b is None else round(p_b, 4),
                    "status_A": status_a,
                    "status_B": status_b,
                    "detail_A": detail_a,
                    "detail_B": detail_b,
                    "n_sources": len(source_keys),
                    "source_key": source_keys[0] if source_keys else "",
                    "pmids": ";".join(sorted(meta.get("pmids") or [])),
                    "dois": ";".join(sorted(meta.get("dois") or [])),
                    "patents": ";".join(sorted(meta.get("patents") or [])),
                    "monomerids": ";".join(sorted(meta.get("monomerids") or [])),
                    "chembl_ids": ";".join(sorted(meta.get("chembl") or [])),
                    "scaffold": murcko(meta.get("smiles", "")),
                }
            )
        pair_ligands[pair] = ligands
    return pair_ligands


def apply_independence(pair_ligands, smiles_cache, chembl_pmids, chembl_dois, chembl_patents, pubmed_status):
    flow_rows = []
    panel_rows = []
    summary_rows = []
    for pair, name_a, name_b in PAIRS:
        ligands = pair_ligands[pair]
        dev_keys, dev_fps, _ = development_molecules(pair)
        map_keys = map_inchikeys(name_a, smiles_cache) | map_inchikeys(name_b, smiles_cache)
        map_cids = map_chembl_ids(name_a) | map_chembl_ids(name_b)
        structure_keys = set(dev_keys) | set(map_keys)
        if pair == "MCL1/Bcl-xL" and not dev_fps:
            mols = map_cids
            for cid in mols:
                fp = fingerprint(smiles_cache.get(cid, ""))
                if fp is not None:
                    dev_fps[cid] = fp
        dev_fp_list = list(dev_fps.values())

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
            tan = max_tanimoto(fingerprint(row["smiles"]), dev_fp_list)
            row = dict(row)
            row["max_ecfp4_to_development"] = "" if tan < 0 else round(tan, 4)
            if tan < 0:
                row["chem_stratum"] = "no_fp"
                chem_primary.append(row)
            elif tan < NEAR_TANIMOTO:
                row["chem_stratum"] = "<0.50"
                chem_primary.append(row)
            elif tan < PRIMARY_TANIMOTO:
                row["chem_stratum"] = "0.50-0.70"
                chem_mid.append(row)
            else:
                row["chem_stratum"] = ">=0.70"
        layers.append(("drop_neighbors_lt_0.70", chem_primary))
        layers.append(("lt_0.70_plus_0.50_0.70_sensitivity", chem_primary + chem_mid))

        for layer_name, rows in layers:
            d, a, b, n = count_classes(rows)
            flow_rows.append(
                {
                    "pair": pair,
                    "layer": layer_name,
                    "n_ligands": len(rows),
                    "n_dual": d,
                    "n_A_only": a,
                    "n_B_only": b,
                    "n_neither": n,
                }
            )

        primary = chem_primary
        stats = class_source_stats(primary)
        role, note = gate_status(stats, pair)
        n_primary_pairs_placeholder = 0
        packaged = 0
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
                "packaged_as_external_evaluation": packaged,
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
        write_csv(
            TAB / f"external_panel_{pair.replace('/', '_').replace('-', '')}_v1.csv",
            [r for r in panel_rows if r["pair"] == pair],
        )

    write_csv(TAB / "external_candidate_flow.csv", flow_rows)
    write_csv(TAB / "external_slice_summary_v1.csv", summary_rows)
    write_csv(TAB / "external_slice_ligands_v1.csv", panel_rows)
    for row in summary_rows:
        path = TAB / f"external_panel_{row['pair'].replace('/', '_').replace('-', '')}_v1.csv"
        if path.exists() and path.stat().st_size:
            (TAB / f"external_panel_{row['pair'].replace('/', '_').replace('-', '')}_v1.sha256").write_text(
                sha256_file(path) + "\n", encoding="utf-8"
            )
    return summary_rows, flow_rows, panel_rows


def freeze_mcl1_panel(smiles_cache: dict[str, str]) -> list[dict]:
    mols_a = json.loads((ROOT / "data/public_pair_selection/mols_MCL1.json").read_text(encoding="utf-8"))
    mols_b = json.loads((ROOT / "data/public_pair_selection/mols_BCL2L1_BclxL.json").read_text(encoding="utf-8"))
    both = sorted(set(mols_a) & set(mols_b))
    fetch_smiles(both, smiles_cache)
    by_class = defaultdict(list)
    for cid in both:
        p_a, p_b = float(mols_a[cid]), float(mols_b[cid])
        cls = four_state(p_a, p_b)
        if not cls:
            continue
        smiles = smiles_cache.get(cid, "")
        key, canon = largest_fragment_inchikey(smiles)
        by_class[cls].append(
            {
                "panel_id": "",
                "molecule_chembl_id": cid,
                "class": cls,
                "pchembl_MCL1": round(p_a, 4),
                "pchembl_BCL2L1": round(p_b, 4),
                "smiles": canon,
                "inchi_key": key,
                "murcko_scaffold": murcko(canon),
                "seed": SEED,
                "docked": 0,
                "note": "ChEMBL map freeze; B-only is exhaustive at n=24; no unused-pool holdout",
            }
        )
    census = {cls: len(rows) for cls, rows in by_class.items()}
    selected = []
    quotas = {"dual": 24, "A_only": 24, "B_only": 24, "neither": 24}
    for cls, quota in quotas.items():
        rows = sorted(by_class.get(cls, []), key=lambda r: r["molecule_chembl_id"])
        rng = hashlib.sha256(f"{SEED}|MCL1|{cls}".encode()).digest()
        # deterministic shuffle
        order = sorted(range(len(rows)), key=lambda i: hashlib.sha256(f"{SEED}|{cls}|{rows[i]['molecule_chembl_id']}".encode()).hexdigest())
        rows = [rows[i] for i in order]
        picked = []
        scaffold_count = Counter()
        cap = 2 if cls != "B_only" else 99
        for row in rows:
            scaf = row["murcko_scaffold"] or row["molecule_chembl_id"]
            if scaffold_count[scaf] >= cap:
                continue
            scaffold_count[scaf] += 1
            picked.append(row)
            if len(picked) >= min(quota, len(rows)):
                break
        if cls == "B_only":
            picked = rows[:]  # exhaustive
        for i, row in enumerate(picked, start=1):
            row = dict(row)
            row["panel_id"] = f"MX{cls[0].upper()}{i:03d}"
            selected.append(row)
    write_csv(TAB / "mcl1_bclxl_chembl_panel96_v1.csv", selected)
    (TAB / "mcl1_bclxl_chembl_panel96_v1.sha256").write_text(sha256_file(TAB / "mcl1_bclxl_chembl_panel96_v1.csv") + "\n", encoding="utf-8")
    summary = [
        {
            "pair": "MCL1/Bcl-xL",
            "theta": THETA,
            "map_paired": len(both),
            "map_dual": census.get("dual", 0),
            "map_A_only": census.get("A_only", 0),
            "map_B_only": census.get("B_only", 0),
            "map_neither": census.get("neither", 0),
            "panel_dual": sum(r["class"] == "dual" for r in selected),
            "panel_A_only": sum(r["class"] == "A_only" for r in selected),
            "panel_B_only": sum(r["class"] == "B_only" for r in selected),
            "panel_neither": sum(r["class"] == "neither" for r in selected),
            "b_only_exhaustive": int(census.get("B_only", 0) <= 24),
            "docked": 0,
            "pose_gold_gate": "not_run_no_vina",
            "domain_role": "PPI_BH3_groove; homologous BCL-2 fold; not a first non-kinase pair",
        }
    ]
    write_csv(TAB / "mcl1_bclxl_panel_freeze_v1.csv", summary)
    return summary


def rcsb_receptor_rows() -> list[dict]:
    specs = [
        ("MCL1", "3WIY", "primary", "LC6", "compound 10 Tanaka 2013"),
        ("BCL2L1", "3WIZ", "primary", "LC6", "compound 10 Tanaka 2013"),
        ("MCL1", "6UDV", "alternate", "", "preselected on resolution/holo WT, not AUROC"),
        ("BCL2L1", "3SP7", "alternate", "", "preselected on resolution/holo WT, not AUROC"),
    ]
    rows = []
    for target, pdb, role, ligand, note in specs:
        entry = request_json(f"{RCSB_ENTRY}/{pdb}", f"rcsb_entry_{pdb}") or {}
        info = entry.get("rcsb_entry_info") or {}
        resolution = ""
        res_list = info.get("resolution_combined") or []
        if res_list:
            resolution = res_list[0]
        method = ""
        if entry.get("exptl"):
            method = entry["exptl"][0].get("method") or ""
        polymer = request_json(f"{RCSB_POLY}/{pdb}/1", f"rcsb_poly_{pdb}_1") or {}
        entity = polymer.get("entity_poly") or {}
        rcsb_entity = polymer.get("rcsb_polymer_entity") or {}
        container = polymer.get("rcsb_polymer_entity_container_identifiers") or {}
        chains = container.get("auth_asym_ids") or container.get("asym_ids") or []
        mutations = rcsb_entity.get("mutation_count") or entity.get("rcsb_mutation_count") or 0
        uniprot = ""
        uniprots = container.get("uniprot_ids") or []
        if uniprots:
            uniprot = uniprots[0]
        rows.append(
            {
                "target": target,
                "pdb_id": pdb,
                "role": role,
                "method": method,
                "resolution_A": resolution,
                "selected_entity": 1,
                "auth_asym_ids": ",".join(chains),
                "primary_chain": chains[0] if chains else "",
                "uniprot": uniprot,
                "mutation_count": mutations,
                "deposited_polymer_residues": info.get("deposited_polymer_monomer_count") or "",
                "deposited_nonpolymer": info.get("deposited_nonpolymer_entity_instance_count") or "",
                "cognate_ligand": ligand,
                "pose_gold_gate_run": 0,
                "note": note + "; chain recorded from RCSB entity 1, not assumed to be A",
            }
        )
        time.sleep(0.1)
    write_csv(TAB / "mcl1_bclxl_receptor_freeze_v1.csv", rows)
    return rows


def literature_comparator() -> list[dict]:
    rows = [
        {
            "resource": "DualFourClass (this work)",
            "doi": "",
            "negative_definition": "experimental A-only and B-only selectives; neither retained descriptively",
            "paired_measurements": "yes, complete-case both targets",
            "chemical_controls": "descriptors, ECFP4, caliper, AND filter, document-blocked CV",
            "receptor_sensitivity": "yes, two PIK3CA crystals and one mTOR crystal",
            "externality": "ChEMBL internal; BindingDB-native slice frozen but not docked",
            "target_scope": "K=4 docking; MCL1/Bcl-xL panel frozen, not docked",
        },
        {
            "resource": "Zhou, Li, Hou 2013 dual-kinase docking",
            "doi": "10.1021/ci400065e",
            "negative_definition": "inactives / non-duals rather than experimental selectives as primary hard negatives",
            "paired_measurements": "dual-kinase actives vs others; not four-state selectives",
            "chemical_controls": "limited relative to later decoy-bias literature",
            "receptor_sensitivity": "not a receptor-realization audit of the DualFourClass type",
            "externality": "ChEMBL/literature dual kinase set of that era",
            "target_scope": "kinase dual inhibitors",
        },
        {
            "resource": "DUD-E",
            "doi": "10.1021/jm300687e",
            "negative_definition": "property-matched putative decoys, not experimental selectives",
            "paired_measurements": "no dual-target paired labels",
            "chemical_controls": "property matching is the design",
            "receptor_sensitivity": "single crystal per target in the original set",
            "externality": "standard VS benchmark",
            "target_scope": "102 targets; single-target VS",
        },
        {
            "resource": "LIT-PCBA",
            "doi": "10.1021/acs.jcim.0c00155",
            "negative_definition": "experimental assay inactives from dose-response PubChem sets",
            "paired_measurements": "no dual-target four-state protocol",
            "chemical_controls": "assay-derived, bias-aware single-target sets",
            "receptor_sensitivity": "not the primary design",
            "externality": "PubChem bioassay origin",
            "target_scope": "15 single targets",
        },
        {
            "resource": "CASF-2016",
            "doi": "10.1021/acs.jcim.8b00545",
            "negative_definition": "scoring / docking / screening power on crystal complexes, not four-state duals",
            "paired_measurements": "protein-ligand complexes",
            "chemical_controls": "core set quality filters",
            "receptor_sensitivity": "scoring power on deposited poses",
            "externality": "PDBbind core set",
            "target_scope": "pose/score benchmark, not dual-target recognition",
        },
        {
            "resource": "DOCKSTRING",
            "doi": "10.1021/acs.jcim.1c01334",
            "negative_definition": "ExCAPE-DB activity labels; docking scores as dataset features",
            "paired_measurements": "many targets, not a four-state dual protocol",
            "chemical_controls": "standardized docking pipeline at scale",
            "receptor_sensitivity": "one prepared receptor per target",
            "externality": "fixed public dataset release",
            "target_scope": "58 targets; target diversity, not dual-selective hard negatives",
        },
    ]
    write_csv(TAB / "benchmark_literature_comparator_v1.csv", rows)
    return rows


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
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.0))
    x = np.arange(len(pairs))
    width = 0.18
    colors = ["#4C78A8", "#F58518", "#54A24B", "#E45756"]
    for i, layer in enumerate(layers):
        mins = []
        for pair in pairs:
            rec = next((r for r in flow_rows if r["pair"] == pair and r["layer"] == layer), None)
            if rec is None:
                mins.append(0)
            else:
                mins.append(min(int(rec["n_dual"]), int(rec["n_A_only"]), int(rec["n_B_only"])))
        axes[0].bar(x + (i - 1.5) * width, mins, width, label=layer.replace("_", "\n"), color=colors[i])
    axes[0].axhline(20, color="0.3", ls="--", lw=0.8, label="primary n=20")
    axes[0].axhline(15, color="0.5", ls=":", lw=0.8, label="thin n=15")
    axes[0].set_xticks(x, [p.replace("/", "/\n") for p in pairs], fontsize=7)
    axes[0].set_ylabel("min(dual, A-only, B-only)")
    axes[0].set_title("A  BindingDB-native independence flow")
    axes[0].legend(fontsize=6, frameon=False, loc="upper right")

    labels = [r["pair"].replace("/", "/\n") for r in summary_rows]
    dual = [int(r["n_dual"]) for r in summary_rows]
    a_only = [int(r["n_A_only"]) for r in summary_rows]
    b_only = [int(r["n_B_only"]) for r in summary_rows]
    x2 = np.arange(len(labels))
    axes[1].bar(x2 - 0.2, dual, 0.2, label="dual", color="#4C78A8")
    axes[1].bar(x2, a_only, 0.2, label="A-only", color="#F58518")
    axes[1].bar(x2 + 0.2, b_only, 0.2, label="B-only", color="#54A24B")
    axes[1].axhline(20, color="0.3", ls="--", lw=0.8)
    axes[1].set_xticks(x2, labels, fontsize=7)
    axes[1].set_title("B  After literature+structure+ECFP<0.70")
    axes[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "FigS_bindingdb_native_slice_v1.png", dpi=160)
    fig.savefig(FIG / "FigS_bindingdb_native_slice_v1.pdf")
    plt.close(fig)


def write_verdict(summary_rows, flow_rows, mcl1, receptors, pubmed_status, lock_rows):
    primary = [r for r in summary_rows if r["gate"] == "primary_external" and r["pair"] != "MCL1/Bcl-xL"]
    thin = [r for r in summary_rows if r["gate"] == "thin_replication"]
    lines = [
        "# BindingDB-native external-slice freeze",
        "",
        "Protocol frozen in `protocol/external_slice_contract.yaml` before docking.",
        "This session did **not** dock and does **not** package external evaluation.",
        f"BindingDB archive `{RELEASE}` md5-verified. ChEMBL document lookup: `{pubmed_status}`.",
        "",
        f"Primary-gate pairs (directional n≥20, ≥3 sources/class, top-document ≤50%): **{len(primary)}**.",
        f"Thin EGFR-style replications: **{len(thin)}**.",
        f"Packaged as external evaluation: **no** (no docking; gates require ≥2 primary pairs plus later pose deposit).",
        "",
        "| pair | dual/A/B/neither after ECFP<0.70 | gate | sources dual/A/B |",
        "|------|--------------------------------:|------|------------------|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['pair']} | {row['n_dual']}/{row['n_A_only']}/{row['n_B_only']}/{row['n_neither']} | {row['gate']} | {row['n_sources_dual']}/{row['n_sources_A_only']}/{row['n_sources_B_only']} |"
        )
    if len(primary) < 2:
        lines += [
            "",
            "Stop rule: fewer than two primary BindingDB-native pairs. Keep the manuscript as a",
            "four-pair formulation audit. Do not call the remaining ligands a database-external set.",
        ]
    else:
        lines += [
            "",
            "Supply is enough to *consider* a later local docking of the frozen panels.",
            "Until those ligands are docked under the frozen receptors/boxes, they remain a supply freeze.",
        ]
    m = mcl1[0]
    lines += [
        "",
        "## MCL1/Bcl-xL",
        "",
        f"ChEMBL map at θ=6.0: dual/A/B/neither {m['map_dual']}/{m['map_A_only']}/{m['map_B_only']}/{m['map_neither']}.",
        f"Frozen panel: {m['panel_dual']}/{m['panel_A_only']}/{m['panel_B_only']}/{m['panel_neither']}.",
        "B-only is exhaustive on the cached map. No same-library holdout. LC6 pose-gold gate was not run (no Vina).",
        "Do not call this pair a disparate-fold pair. It is a PPI/BH3 groove domain shift.",
        "",
        "## Primary receptors",
        "",
    ]
    for rec in receptors:
        lines.append(
            f"- {rec['target']} {rec['pdb_id']} ({rec['role']}): chain `{rec['primary_chain']}` entity 1, {rec['resolution_A']} Å, mutations={rec['mutation_count']}."
        )
    lines += [
        "",
        "No BindingDB AUROC was computed. Do not inspect these counts and then change θ, boxes, or receptors.",
        f"Archive sha256 values are in `tables/bindingdb_archive_lock_v1.csv` ({len(lock_rows)} files).",
    ]
    (ANALYSIS / "EXTERNAL_SLICE_FREEZE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    print("lock BindingDB archives", flush=True)
    lock_rows = download_archives()
    print("stream native UniProt rows", flush=True)
    records = stream_native_records()
    assays = load_rsid_assays()
    print("aggregate ligands", flush=True)
    pair_ligands = build_pair_ligands(records, assays)
    smiles_cache = load_smiles_cache()
    mcl1_ids = json.loads((ROOT / "data/public_pair_selection/mols_MCL1.json").read_text(encoding="utf-8"))
    bcl_ids = json.loads((ROOT / "data/public_pair_selection/mols_BCL2L1_BclxL.json").read_text(encoding="utf-8"))
    paired = sorted(set(mcl1_ids) & set(bcl_ids))
    print(f"fetch {len(paired)} paired MCL1/Bcl-xL SMILES", flush=True)
    fetch_smiles(paired, smiles_cache)
    write_csv(
        TAB / "chembl_smiles_cache_mcl1_v1.csv",
        [{"chembl_id": cid, "smiles": smiles_cache[cid]} for cid in paired if cid in smiles_cache],
    )
    print("ChEMBL document PMID/DOI", flush=True)
    chembl_pmids, chembl_dois, chembl_patents, pubmed_status = collect_dev_documents()
    print(f"dev literature pmids={len(chembl_pmids)} dois={len(chembl_dois)} patents={len(chembl_patents)} status={pubmed_status}", flush=True)
    print("independence layers", flush=True)
    summary_rows, flow_rows, _panel = apply_independence(
        pair_ligands, smiles_cache, chembl_pmids, chembl_dois, chembl_patents, pubmed_status
    )
    print("MCL1 panel freeze", flush=True)
    mcl1 = freeze_mcl1_panel(smiles_cache)
    print("receptor metadata", flush=True)
    receptors = rcsb_receptor_rows()
    literature_comparator()
    plot_flow(flow_rows, summary_rows)
    write_verdict(summary_rows, flow_rows, mcl1, receptors, pubmed_status, lock_rows)
    print("wrote BindingDB-native slice tables", flush=True)
    for row in summary_rows:
        print(
            row["pair"],
            row["n_dual"],
            row["n_A_only"],
            row["n_B_only"],
            row["n_neither"],
            row["gate"],
            flush=True,
        )


if __name__ == "__main__":
    main()
