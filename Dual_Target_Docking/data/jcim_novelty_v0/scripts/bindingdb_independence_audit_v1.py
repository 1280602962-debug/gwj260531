#!/usr/bin/env python3
"""BindingDB independence audit after the frozen time-split failed.

No docking. Rules were frozen in docs/BINDINGDB_EXTERNAL_SOP.md before
independent class counts were inspected: equal-relation IC50/Ki/Kd/EC50,
θ = 6.0 four-state labels, drop panel InChIKeys, drop panel PMIDs.
A ligand is database-external only if UniChem maps no ChEMBL ID in the
usable-pChEMBL map and no retained BindingDB PMID matches a panel document.
Even if class counts pass, this file does not compute AUROC and is not
packaged as external validation without docking.

ChEMBL molecule.json batch lookups are avoided (HTTP 500 on this snapshot).
Structure overlap uses UniChem InChIKey → ChEMBL ID against local maps.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
PUBLIC = ROOT / "data" / "public_pair_selection"
CACHE = ROOT / "data" / "jcim_novelty_v0" / "cache" / "bindingdb_independence_v1"
TAB.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNICHEM = "https://www.ebi.ac.uk/unichem/rest/inchikey"
BDB = "https://bindingdb.org/rest/getLigandsByUniprots"
UA = "DualFourClass-Bench/1.0 (independence audit; no docking)"
KEEP_TYPES = {"IC50", "KI", "KD", "EC50"}
THETA = 6.0
CUTOFF_NM = 1_000_000
MIN_CLASS = 10
MIN_PAIRS = 2
NUM_RE = re.compile(r"([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)")

TARGETS = {
    "EGFR": "P00533",
    "HER2": "P04626",
    "ACHE": "P22303",
    "BCHE": "P06276",
    "PIK3CA": "P42336",
    "PIK3CB": "P42338",
    "MTOR": "P42345",
}
PAIRS = (
    ("EGFR/HER2", "EGFR", "HER2"),
    ("AChE/BChE", "ACHE", "BCHE"),
    ("PIK3CA/PIK3CB", "PIK3CA", "PIK3CB"),
    ("PIK3CA/mTOR", "PIK3CA", "MTOR"),
)
CHEMBL_FILES = {
    "EGFR": "mols_EGFR.json",
    "HER2": "mols_HER2.json",
    "ACHE": "mols_ACHE.json",
    "BCHE": "mols_BCHE.json",
    "PIK3CA": "mols_PIK3CA.json",
    "PIK3CB": "mols_PIK3CB.json",
    "MTOR": "mols_MTOR.json",
}
PANEL_SMILES = {
    "EGFR/HER2": "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
    "AChE/BChE": "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
    "PIK3CA/PIK3CB": "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
    "PIK3CA/mTOR": "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
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


def request_json(url: str, cache_key: str, timeout: int = 180, required: bool = True, attempts: int = 5):
    cache = CACHE / f"{cache_key}.json"
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
            time.sleep(1.5 * (attempt + 1))
    if required:
        raise RuntimeError(f"request failed {url}: {last_error}")
    print(f"skip {cache_key}: {last_error}", flush=True)
    return None


def inchikey_from_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return ""
    try:
        return inchi.MolToInchiKey(mol) or ""
    except Exception:
        return ""


def parse_affinity(raw: str):
    s = (raw or "").strip().replace(",", "")
    if not s:
        return "", None
    qual, rest = "", s
    if s[0] in "<>~":
        qual, rest = s[0], s[1:].strip()
    elif s[:2] in (">=", "<=", "~="):
        qual, rest = s[0], s[2:].strip()
    match = NUM_RE.search(rest)
    if not match:
        return qual, None
    return qual, float(match.group(1))


def nm_to_p(nm: float | None) -> float | None:
    if nm is None or nm <= 0:
        return None
    return 9.0 - math.log10(nm)


def classify(p_a: float, p_b: float) -> str:
    a, b = p_a >= THETA, p_b >= THETA
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def chunks(values: list[str], size: int):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def chembl_map_ids() -> dict[str, set[str]]:
    return {
        name: set(json.loads((PUBLIC / filename).read_text(encoding="utf-8")))
        for name, filename in CHEMBL_FILES.items()
    }


def unichem_chembl_ids(key: str) -> tuple[set[str], str]:
    """Return ChEMBL IDs mapped to an InChIKey, plus a status string."""
    if not key:
        return set(), "no_inchikey"
    cache_key = f"unichem_{hashlib.sha256(key.encode()).hexdigest()[:16]}"
    payload = request_json(f"{UNICHEM}/{key}", cache_key, timeout=60, required=False)
    if payload is None:
        return set(), "unichem_fail"
    ids = set()
    rows = payload if isinstance(payload, list) else []
    for row in rows:
        if str(row.get("src_id")) != "1":
            continue
        cid = str(row.get("src_compound_id") or "").strip()
        if not cid:
            continue
        ids.add(cid if cid.startswith("CHEMBL") else f"CHEMBL{cid}")
    return ids, "ok"


def document_pubmed(doc_ids: list[str]) -> tuple[dict[str, str], str]:
    lookup_path = TAB / "document_pubmed_lookup_v1.csv"
    known = {}
    if lookup_path.exists():
        for row in read_csv(lookup_path):
            known[row["document_chembl_id"]] = row.get("pubmed_id") or ""
    missing = [doc for doc in doc_ids if doc not in known]
    status = "ok"
    if missing:
        probe = missing[0]
        payload = request_json(
            f"{CHEMBL}/document.json?document_chembl_id={probe}",
            f"docpub1_{probe}",
            timeout=20,
            required=False,
            attempts=2,
        )
        if payload is None:
            status = "chembl_document_api_unavailable"
            missing = []
    for batch in chunks(missing, 5):
        query = urllib.parse.urlencode(
            {"document_chembl_id__in": ",".join(batch), "limit": max(len(batch), 5)}
        )
        digest = hashlib.sha256(",".join(batch).encode()).hexdigest()[:16]
        payload = request_json(
            f"{CHEMBL}/document.json?{query}",
            f"docpub_{digest}",
            timeout=30,
            required=False,
            attempts=2,
        )
        if payload is None:
            status = "chembl_document_api_unavailable"
            break
        hits = {row["document_chembl_id"]: row for row in payload.get("documents") or []}
        for doc in batch:
            pmid = hits.get(doc, {}).get("pubmed_id")
            known[doc] = "" if pmid in (None, "") else str(pmid)
    write_csv(
        lookup_path,
        [
            {"document_chembl_id": doc, "pubmed_id": pmid}
            for doc, pmid in sorted(known.items())
        ],
    )
    return known, status


def fetch_bindingdb(uniprot: str) -> list[dict]:
    url = f"{BDB}?uniprot={uniprot}&cutoff={CUTOFF_NM}&response=application/json"
    payload = request_json(url, f"bdb_{uniprot}_{CUTOFF_NM}", timeout=300)
    block = payload.get("getLindsByUniprotsResponse") or payload.get("getLigandsByUniprotsResponse") or {}
    affinities = block.get("affinities")
    if affinities is None:
        return []
    if isinstance(affinities, dict):
        return [affinities]
    return affinities


def aggregate_target(rows: list[dict]):
    pmax = defaultdict(dict)
    smiles = {}
    pmids = defaultdict(set)
    dois = defaultdict(set)
    types = defaultdict(set)
    for row in rows:
        atype = str(row.get("affinity_type") or "").strip().upper()
        if atype not in KEEP_TYPES:
            continue
        qual, nm = parse_affinity(str(row.get("affinity") or ""))
        if qual not in {"", "="}:
            continue
        p = nm_to_p(nm)
        if p is None:
            continue
        mid = str(row.get("monomerid") or "").strip()
        if not mid:
            continue
        prev = pmax[mid].get(atype)
        if prev is None or p > prev:
            pmax[mid][atype] = p
        mixed_prev = pmax[mid].get("MIXED")
        if mixed_prev is None or p > mixed_prev:
            pmax[mid]["MIXED"] = p
        if row.get("smile") and mid not in smiles:
            smiles[mid] = row["smile"]
        if row.get("pmid"):
            pmids[mid].add(str(row["pmid"]))
        if row.get("doi"):
            dois[mid].add(str(row["doi"]).lower())
        types[mid].add(atype)
    return pmax, smiles, pmids, dois, types


def panel_inchikeys() -> dict[str, set[str]]:
    out = defaultdict(set)
    for pair, relative in PANEL_SMILES.items():
        for row in read_csv(ROOT / relative):
            key = inchikey_from_smiles(row.get("smiles") or "")
            if key:
                out[pair].add(key)
    return out


def gate(counts: dict[str, int]) -> str:
    if min(counts["dual"], counts["A_only"], counts["B_only"]) >= MIN_CLASS:
        return "supply_enough_to_dock"
    if min(counts["dual"], counts["A_only"], counts["B_only"]) > 0:
        return "descriptive_only"
    return "unevaluable"


def probe_unichem() -> str:
    payload = request_json(
        f"{UNICHEM}/BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "unichem_probe_aspirin",
        timeout=12,
        required=False,
        attempts=1,
    )
    return "ok" if payload is not None else "unichem_unavailable"


def main() -> None:
    map_ids = chembl_map_ids()
    print("loaded ChEMBL maps", {k: len(v) for k, v in map_ids.items()}, flush=True)
    unichem_status = "unichem_skipped_bulk"
    print(
        "UniChem skipped for bulk BindingDB ligands (per-key REST is too slow for this audit). "
        "ChEMBL-map overlap is therefore not claimed.",
        flush=True,
    )

    panel_keys = panel_inchikeys()
    bdb_by_target = {}
    for name, uniprot in TARGETS.items():
        print(f"BindingDB {name} {uniprot}", flush=True)
        bdb_by_target[name] = aggregate_target(fetch_bindingdb(uniprot))
        print(f"  monomers={len(bdb_by_target[name][0])}", flush=True)

    activities = [
        row
        for row in read_csv(TAB / "high_confidence_activity_audit_v1.csv")
        if row.get("keep") in ("1", "True", "true")
    ]
    panel_docs = sorted({row["document_chembl_id"] for row in activities if row.get("document_chembl_id")})
    pubmed, pubmed_status = document_pubmed(panel_docs)
    panel_pmids = {pmid for pmid in pubmed.values() if pmid}
    print(f"panel documents={len(panel_docs)} pubmed_mapped={len(panel_pmids)} status={pubmed_status}", flush=True)

    unichem_cache: dict[str, tuple[set[str], str]] = {}
    ligand_rows = []
    summary_rows = []
    for pair, name_a, name_b in PAIRS:
        pmax_a, smiles_a, pmids_a, dois_a, types_a = bdb_by_target[name_a]
        pmax_b, smiles_b, pmids_b, dois_b, types_b = bdb_by_target[name_b]
        chembl_ids = map_ids[name_a] | map_ids[name_b]
        pair_panel = panel_keys[pair]
        both = set(pmax_a) & set(pmax_b)
        buckets = defaultdict(lambda: defaultdict(int))
        n_pmid_overlap = 0
        n_struct_panel = 0
        n_struct_map = 0
        n_parse_fail = 0
        n_unichem_fail = 0
        print(f"{pair}: both-end BindingDB monomers={len(both)}", flush=True)
        for mid in sorted(both):
            p_a = pmax_a[mid].get("MIXED")
            p_b = pmax_b[mid].get("MIXED")
            if p_a is None or p_b is None:
                continue
            cls = classify(p_a, p_b)
            buckets["all_equal_mixed"][cls] += 1
            smiles = smiles_a.get(mid) or smiles_b.get(mid) or ""
            key = inchikey_from_smiles(smiles)
            if not key:
                n_parse_fail += 1
            if unichem_status != "ok":
                chembl_hits, uni_status = set(), unichem_status
            elif key not in unichem_cache:
                unichem_cache[key] = unichem_chembl_ids(key)
                chembl_hits, uni_status = unichem_cache[key]
            else:
                chembl_hits, uni_status = unichem_cache[key]
            if uni_status == "unichem_fail":
                n_unichem_fail += 1
            in_panel = bool(key and key in pair_panel)
            in_map = bool(chembl_hits & chembl_ids)
            lig_pmids = pmids_a.get(mid, set()) | pmids_b.get(mid, set())
            pmid_hit = bool(lig_pmids & panel_pmids)
            if in_panel:
                n_struct_panel += 1
            if in_map:
                n_struct_map += 1
            if pmid_hit:
                n_pmid_overlap += 1
            if not in_panel:
                buckets["not_panel_structure"][cls] += 1
            if uni_status == "ok" and not in_map:
                buckets["not_chembl_map_structure"][cls] += 1
            if pubmed_status == "ok" and not pmid_hit:
                buckets["not_panel_pmid"][cls] += 1
            independent = (
                bool(key)
                and uni_status == "ok"
                and (not in_map)
                and (pubmed_status == "ok")
                and (not pmid_hit)
            )
            structure_indep = bool(key) and uni_status == "ok" and (not in_map)
            panel_struct_indep = bool(key) and (not in_panel)
            if panel_struct_indep:
                buckets["not_panel_structure_class"][cls] += 1
            if structure_indep:
                buckets["structure_independent"][cls] += 1
            if independent:
                buckets["structure_and_literature_independent"][cls] += 1
                ic50_a = pmax_a[mid].get("IC50")
                ic50_b = pmax_b[mid].get("IC50")
                if ic50_a is not None and ic50_b is not None:
                    buckets["ic50_only_independent"][classify(ic50_a, ic50_b)] += 1
            ligand_rows.append(
                {
                    "pair": pair,
                    "bindingdb_monomerid": mid,
                    "inchikey": key,
                    "frozen_class_theta6_mixed": cls,
                    "pA_mixed": round(p_a, 4),
                    "pB_mixed": round(p_b, 4),
                    "in_scored_panel_structure": int(in_panel),
                    "in_chembl_map_structure": int(in_map),
                    "panel_pmid_overlap": int(pmid_hit),
                    "unichem_status": uni_status,
                    "n_pmids": len(lig_pmids),
                    "independent_structure_and_literature": int(independent),
                    "affinity_types": ";".join(sorted((types_a.get(mid) or set()) | (types_b.get(mid) or set()))),
                }
            )
        n_both = sum(buckets["all_equal_mixed"].values())
        indep = buckets["structure_and_literature_independent"]
        struct_only = buckets["structure_independent"]
        panel_only = buckets["not_panel_structure_class"]
        mixed = buckets["all_equal_mixed"]
        summary_rows.append(
            {
                "pair": pair,
                "n_bindingdb_both_equal_mixed": n_both,
                "mixed_dual": mixed["dual"],
                "mixed_A_only": mixed["A_only"],
                "mixed_B_only": mixed["B_only"],
                "mixed_neither": mixed["neither"],
                "n_inchikey_parse_fail": n_parse_fail,
                "n_unichem_fail": n_unichem_fail,
                "unichem_status": unichem_status,
                "n_structure_in_scored_panel": n_struct_panel,
                "n_structure_in_chembl_map": n_struct_map,
                "n_pmid_overlap_panel_docs": n_pmid_overlap,
                "frac_pmid_overlap": round(n_pmid_overlap / n_both, 4) if n_both else "",
                "frac_structure_in_chembl_map": round(n_struct_map / n_both, 4) if n_both else "",
                "pubmed_lookup_status": pubmed_status,
                "panel_indep_dual": panel_only["dual"],
                "panel_indep_A_only": panel_only["A_only"],
                "panel_indep_B_only": panel_only["B_only"],
                "panel_indep_neither": panel_only["neither"],
                "gate_panel_structure": gate(panel_only),
                "struct_indep_dual": struct_only["dual"],
                "struct_indep_A_only": struct_only["A_only"],
                "struct_indep_B_only": struct_only["B_only"],
                "struct_indep_neither": struct_only["neither"],
                "indep_dual": indep["dual"],
                "indep_A_only": indep["A_only"],
                "indep_B_only": indep["B_only"],
                "indep_neither": indep["neither"],
                "ic50_only_indep_dual": buckets["ic50_only_independent"]["dual"],
                "ic50_only_indep_A_only": buckets["ic50_only_independent"]["A_only"],
                "ic50_only_indep_B_only": buckets["ic50_only_independent"]["B_only"],
                "gate_independent": gate(indep),
                "packaged_as_external_validation": 0,
                "note": (
                    "no docking; literature+ChEMBL-map independence requires UniChem and ChEMBL document APIs; "
                    f"unichem={unichem_status}; pubmed={pubmed_status}"
                ),
            }
        )

    n_pass = sum(row["gate_independent"] == "supply_enough_to_dock" for row in summary_rows)
    write_csv(TAB / "bindingdb_independence_ligands_v1.csv", ligand_rows)
    write_csv(TAB / "bindingdb_independence_summary_v1.csv", summary_rows)

    lines = [
        "# BindingDB independence audit",
        "",
        "Protocol frozen before inspecting independent class counts.",
        "This is **not** external validation: no new docking was run.",
        f"Panel PMID lookup: `{pubmed_status}`. UniChem: `{unichem_status}`.",
        "",
        f"Pairs with database+literature independent dual/A-only/B-only each ≥{MIN_CLASS}: **{n_pass}**.",
        "Packaged as external validation: **no**.",
        "",
        "| pair | BDB both dual/A/B/neither | panel overlap n | not-in-panel dual/A/B/neither | panel-structure gate | lit+map gate |",
        "|------|---------------------------|----------------:|-------------------------------|----------------------|--------------|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['pair']} | "
            f"{row['mixed_dual']}/{row['mixed_A_only']}/{row['mixed_B_only']}/{row['mixed_neither']} | "
            f"{row['n_structure_in_scored_panel']} | "
            f"{row['panel_indep_dual']}/{row['panel_indep_A_only']}/{row['panel_indep_B_only']}/{row['panel_indep_neither']} | "
            f"{row['gate_panel_structure']} | {row['gate_independent']} |"
        )
    if unichem_status != "ok" or pubmed_status != "ok":
        lines.append(
            "\nChEMBL-map and/or PMID overlap could not be completed in this session. "
            "Panel-structure overlap is local and is **not** database-external independence. "
            "Do not package as external validation. Keep the internal formulation-audit claim."
        )
    elif n_pass < MIN_PAIRS:
        lines.append(
            "\nStop condition: fewer than two pairs have an independent directional set with n≥10. "
            "Keep the internal formulation-audit claim. Do not dock a leftover slice to chase AUROC."
        )
    else:
        lines.append(
            "\nSupply would permit a docking campaign on the independent slice. "
            "That campaign was not run here and is a local step."
        )
    (ANALYSIS / "BINDINGDB_INDEPENDENCE_VERDICT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
