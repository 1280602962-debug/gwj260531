#!/usr/bin/env python3
"""Build the first-pass, paper-level assay source audit queue.

The queue contains every source document tied for the maximum pChEMBL value
within each priority ligand-target group.  It is intentionally narrower than
the set of every historical supporting document, but it is the minimum source
set that can directly change the frozen max-aggregation labels.
"""

from __future__ import annotations

import csv
import html
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
AUDIT = TABLES / "assay_context_audit.csv"
DOC_MAP = TABLES / "chembl_document_pubmed_doi_v1.csv"
OUT_CSV = TABLES / "assay_source_reading_list_v1.csv"
OUT_FULL_CSV = TABLES / "assay_source_full_inventory_v1.csv"
OUT_MD = ROOT.parents[1] / "docs" / "ASSAY_SOURCE_READING_LIST_V1.md"


# Entries absent from the repository's frozen ChEMBL document mapping but
# independently resolved against the publisher/PubMed record.
OVERRIDES = {
    "CHEMBL4665812": {
        "title": "Targeting Her2-insYVMA with Covalent Inhibitors-A Focused Compound Screening and Structure-Based Design Approach",
        "pubmed_id": "32931277",
        "doi": "10.1021/acs.jmedchem.0c00870",
        "source_type": "primary_article",
    },
    "CHEMBL4680246": {
        "title": "Design, synthesis and biological evaluation of novel O-carbamoyl ferulamide derivatives as multi-target-directed ligands for the treatment of Alzheimer's disease",
        "doi": "10.1016/j.ejmech.2020.112265",
        "source_type": "primary_article",
    },
    "CHEMBL4765307": {
        "title": "Structure-Based Drug Design and Synthesis of PI3Kalpha-Selective Inhibitor (PF-06843195)",
        "pubmed_id": "33356246",
        "doi": "10.1021/acs.jmedchem.0c01652",
        "source_type": "primary_article",
    },
    "CHEMBL5131445": {
        "title": "Optimization of a novel piperazinone series as potent selective peripheral covalent BTK inhibitors",
        "pubmed_id": "35041943",
        "doi": "10.1016/j.bmcl.2022.128549",
        "source_type": "primary_article",
    },
    "CHEMBL5500428": {
        "title": "Identification of 6-Anilino Imidazo[4,5-c]pyridin-2-ones as Selective DNA-Dependent Protein Kinase Inhibitors and Their Application as Radiosensitizers",
        "pubmed_id": "39007759",
        "doi": "10.1021/acs.jmedchem.4c01120",
        "source_type": "primary_article",
    },
    "CHEMBL5620391": {
        "title": "Discovery of N-(2-chloro-5-(3-(pyridin-4-yl)-1H-pyrazolo[3,4-b]pyridin-5-yl)pyridin-3-yl)-4-fluorobenzenesulfonamide (FD274) as a highly potent PI3K/mTOR dual inhibitor for the treatment of acute myeloid leukemia",
        "pubmed_id": "37329712",
        "doi": "10.1016/j.ejmech.2023.115543",
        "source_type": "primary_article",
    },
    "CHEMBL1909046": {
        "title": "DrugMatrix protein kinase inhibition dataset",
        "source_type": "deposited_dataset",
    },
    "CHEMBL1201862": {
        "title": "ChEMBL kinase profiling dataset/source",
        "source_type": "deposited_dataset",
    },
    "CHEMBL5446079": {
        "title": "Affinity Biochemical Literature for EUbOPEN Chemogenomic Library",
        "source_type": "deposited_dataset",
    },
    "CHEMBL5465560": {
        "title": "Selectivity Literature for EUbOPEN Chemogenomic Library",
        "source_type": "deposited_dataset",
    },
    "CHEMBL5252533": {
        "title": "Dual Kinase-Bromodomain Inhibitors in Anticancer Drug Discovery: A Structural and Pharmacological Perspective",
        "source_type": "review_article",
    },
    "CHEMBL5214883": {
        "source_type": "review_article",
    },
    "CHEMBL5728076": {
        "title": "ErbB receptor inhibitors",
        "source_type": "patent",
    },
    "CHEMBL5728361": {
        "title": "(S)-2-(1-cyclopropylethyl)-5-(4-methyl-2-((6-(2-oxopyrrolidin-1-yl)pyridin-2-yl)amino)thiazol-5-yl)-7-(methylsulfonyl)isoindolin-1-one as a phosphatidylinositol 3-kinase inhibitor",
        "source_type": "patent",
    },
    "CHEMBL5729445": {
        "title": "Tetrahydroquinoline derivatives and a process of preparation thereof",
        "source_type": "patent",
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def fetch_pubmed(pmids: list[str]) -> dict[str, dict[str, str]]:
    if not pmids:
        return {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids), "retmode": "json"}
    )
    try:
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = json.load(response)
    except Exception:
        return {}
    result = payload.get("result", {})
    return {
        pmid: {
            "title": html.unescape(result.get(pmid, {}).get("title", "")).strip(),
            "journal": result.get(pmid, {}).get("fulljournalname", "").strip(),
            "publication_date": result.get(pmid, {}).get("pubdate", "").strip(),
        }
        for pmid in pmids
        if pmid in result
    }


def fetch_pubmed_types(pmids: list[str]) -> dict[str, str]:
    if not pmids:
        return {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    )
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            root = ET.parse(response).getroot()
    except Exception:
        return {}
    result = {}
    for article in root.findall(".//PubmedArticle"):
        pmid_node = article.find(".//MedlineCitation/PMID")
        if pmid_node is None or not pmid_node.text:
            continue
        types = sorted(
            {
                node.text.strip()
                for node in article.findall(".//PublicationTypeList/PublicationType")
                if node.text and node.text.strip()
            }
        )
        result[pmid_node.text.strip()] = ";".join(types)
    return result


def fetch_crossref(doi: str) -> dict[str, str]:
    if not doi:
        return {}
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(url, headers={"User-Agent": "JCIM-assay-audit/1.0 (mailto:repository-audit@example.invalid)"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            item = json.load(response).get("message", {})
    except Exception:
        return {}
    titles = item.get("title") or []
    journals = item.get("container-title") or []
    return {
        "title": html.unescape(titles[0]).strip() if titles else "",
        "journal": html.unescape(journals[0]).strip() if journals else "",
    }


def fetch_chembl_document(doc_id: str) -> dict[str, str]:
    url = f"https://www.ebi.ac.uk/chembl/api/data/document/{doc_id}.json"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            item = json.load(response)
    except Exception:
        return {}
    fields = (
        "title",
        "journal",
        "journal_full_title",
        "year",
        "pubmed_id",
        "doi",
        "patent_id",
        "doc_type",
        "src_id",
    )
    return {field: str(item[field]).strip() for field in fields if item.get(field) is not None}


def infer_source_type(meta: dict[str, str]) -> str:
    if meta.get("source_type"):
        return meta["source_type"]
    if meta.get("patent_id"):
        return "patent"
    title = meta.get("title", "").lower()
    if "review" in meta.get("publication_types", "").lower():
        return "review_article"
    if "review" in title or "perspective" in title:
        return "review_article"
    if meta.get("doi") or meta.get("pubmed_id"):
        return "article_unclassified"
    return "unresolved_chembl_source"


def links(meta: dict[str, str], doc_id: str) -> tuple[str, str]:
    chembl_url = f"https://www.ebi.ac.uk/chembl/explore/document/{doc_id}"
    if meta.get("doi"):
        direct = "https://doi.org/" + meta["doi"]
    elif meta.get("pubmed_id"):
        direct = "https://pubmed.ncbi.nlm.nih.gov/" + meta["pubmed_id"] + "/"
    elif meta.get("patent_id"):
        direct = "https://patents.google.com/patent/" + meta["patent_id"].replace("-", "")
    else:
        direct = ""
    return direct, chembl_url


def main() -> None:
    audit = read_csv(AUDIT)
    maxima: dict[tuple[str, str, str], float] = {}
    for row in audit:
        key = (row["pair"], row["ligand"], row["target_chembl_id"])
        value = float(row["pchembl_value"])
        maxima[key] = max(maxima.get(key, value), value)

    decisive = []
    for row in audit:
        key = (row["pair"], row["ligand"], row["target_chembl_id"])
        if abs(float(row["pchembl_value"]) - maxima[key]) < 1e-12:
            decisive.append(row)

    doc_meta = {row["document_chembl_id"]: row for row in read_csv(DOC_MAP)}
    if OUT_CSV.exists():
        for cached in read_csv(OUT_CSV):
            meta = doc_meta.setdefault(
                cached["document_chembl_id"],
                {"document_chembl_id": cached["document_chembl_id"]},
            )
            for field in ("title", "journal", "year", "pubmed_id", "doi", "patent_id"):
                if not meta.get(field) and cached.get(field):
                    meta[field] = cached[field]
            if cached.get("source_type") not in {"", "unresolved_chembl_source", "article_unclassified"}:
                meta.setdefault("source_type", cached["source_type"])
    for doc_id, override in OVERRIDES.items():
        doc_meta.setdefault(doc_id, {"document_chembl_id": doc_id})
        doc_meta[doc_id].update({k: v for k, v in override.items() if v})

    pmids = sorted({m.get("pubmed_id", "") for m in doc_meta.values()} - {""})
    pubmed_meta: dict[str, dict[str, str]] = {}
    pubmed_types: dict[str, str] = {}
    for start in range(0, len(pmids), 100):
        batch = pmids[start : start + 100]
        pubmed_meta.update(fetch_pubmed(batch))
        pubmed_types.update(fetch_pubmed_types(batch))
        time.sleep(0.35)

    for meta in doc_meta.values():
        pm = pubmed_meta.get(meta.get("pubmed_id", ""), {})
        for field in ("title", "journal", "publication_date"):
            if not meta.get(field) and pm.get(field):
                meta[field] = pm[field]
        if pubmed_types.get(meta.get("pubmed_id", "")):
            meta["publication_types"] = pubmed_types[meta["pubmed_id"]]

    decisive_doc_ids = {row["document_chembl_id"] for row in decisive}
    unresolved_ids = [
        doc_id
        for doc_id in sorted(decisive_doc_ids)
        if not doc_meta.setdefault(doc_id, {"document_chembl_id": doc_id}).get("title")
    ]
    with ThreadPoolExecutor(max_workers=8) as pool:
        live_records = dict(zip(unresolved_ids, pool.map(fetch_chembl_document, unresolved_ids)))
    for doc_id, live in live_records.items():
        meta = doc_meta[doc_id]
        for key, value in live.items():
            if not meta.get(key):
                meta[key] = value

    by_doc: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in decisive:
        by_doc[row["document_chembl_id"]].append(row)

    records = []
    for doc_id in sorted(by_doc):
        rows = by_doc[doc_id]
        meta = dict(doc_meta.get(doc_id, {"document_chembl_id": doc_id}))
        pm = pubmed_meta.get(meta.get("pubmed_id", ""), {})
        for field in ("title", "journal", "publication_date"):
            if not meta.get(field) and pm.get(field):
                meta[field] = pm[field]
        if not meta.get("title") and meta.get("doi"):
            meta.update({k: v for k, v in fetch_crossref(meta["doi"]).items() if v})
            time.sleep(0.12)
        source_type = infer_source_type(meta)
        direct_url, chembl_url = links(meta, doc_id)
        pairs = sorted({row["pair"] for row in rows})
        ligands = sorted({row["ligand"] for row in rows})
        targets = sorted({row["target_name"] for row in rows})
        assay_ids = sorted({row["assay_chembl_id"] for row in rows})
        if source_type == "primary_article":
            action = "read_full_text_and_supporting_information"
        elif source_type in {"review_article", "possible_review_check_primary"}:
            action = "trace_each_value_to_cited_primary_source"
        elif source_type == "patent":
            action = "audit_patent_examples_and_assay_protocol"
        elif source_type == "deposited_dataset":
            action = "audit_dataset_assay_description_and_provenance"
        elif source_type == "article_unclassified":
            action = "classify_then_read_full_text_and_supporting_information"
        else:
            action = "resolve_bibliography_from_chembl_then_audit"
        records.append(
            {
                "document_chembl_id": doc_id,
                "source_type": source_type,
                "title": meta.get("title", ""),
                "journal": meta.get("journal", ""),
                "year": meta.get("year", "") or meta.get("publication_date", ""),
                "pubmed_id": meta.get("pubmed_id", ""),
                "doi": meta.get("doi", ""),
                "patent_id": meta.get("patent_id", ""),
                "publication_types": meta.get("publication_types", ""),
                "direct_url": direct_url,
                "chembl_url": chembl_url,
                "pairs": ";".join(pairs),
                "ligands": ";".join(ligands),
                "targets": ";".join(targets),
                "assay_chembl_ids": ";".join(assay_ids),
                "n_decisive_rows": str(len(rows)),
                "required_action": action,
            }
        )

    fields = list(records[0])
    with OUT_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    by_all_doc: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in audit:
        by_all_doc[row["document_chembl_id"]].append(row)
    full_records = []
    for doc_id in sorted(by_all_doc):
        rows = by_all_doc[doc_id]
        meta = dict(doc_meta.get(doc_id, {"document_chembl_id": doc_id}))
        pm = pubmed_meta.get(meta.get("pubmed_id", ""), {})
        for field in ("title", "journal", "publication_date"):
            if not meta.get(field) and pm.get(field):
                meta[field] = pm[field]
        if pubmed_types.get(meta.get("pubmed_id", "")):
            meta["publication_types"] = pubmed_types[meta["pubmed_id"]]
        source_type = infer_source_type(meta)
        direct_url, chembl_url = links(meta, doc_id)
        mandatory = doc_id in by_doc
        full_records.append(
            {
                "document_chembl_id": doc_id,
                "audit_tier": "mandatory_first_pass" if mandatory else "conditional_follow_up",
                "source_type": source_type,
                "title": meta.get("title", ""),
                "journal": meta.get("journal", ""),
                "year": meta.get("year", "") or meta.get("publication_date", ""),
                "pubmed_id": meta.get("pubmed_id", ""),
                "doi": meta.get("doi", ""),
                "patent_id": meta.get("patent_id", ""),
                "publication_types": meta.get("publication_types", ""),
                "direct_url": direct_url,
                "chembl_url": chembl_url,
                "pairs": ";".join(sorted({row["pair"] for row in rows})),
                "ligands": ";".join(sorted({row["ligand"] for row in rows})),
                "targets": ";".join(sorted({row["target_name"] for row in rows})),
                "assay_chembl_ids": ";".join(sorted({row["assay_chembl_id"] for row in rows})),
                "n_activity_rows": str(len(rows)),
                "review_trigger": (
                    "always_review_now"
                    if mandatory
                    else "review_if_decisive_source_is_excluded_or_for_all-measurement_sensitivity"
                ),
            }
        )
    with OUT_FULL_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(full_records[0]))
        writer.writeheader()
        writer.writerows(full_records)

    counts: dict[str, int] = defaultdict(int)
    for record in records:
        counts[record["source_type"]] += 1
    md = [
        "# Assay source reading list v1",
        "",
        "This is the minimum first-pass paper/source audit set for the 186 priority ligands. It includes every source tied for the maximum pChEMBL value in each ligand-target group; a rejected decisive source triggers review of the next-highest source and is therefore not the final upper bound.",
        "",
        f"- Priority assay rows: {len(audit)}",
        f"- Decisive max-tied rows: {len(decisive)}",
        f"- Unique decisive source documents: {len(records)}",
        f"- Complete historical source inventory: {len(full_records)} (see `{OUT_FULL_CSV.name}`)",
        "- Source-type counts: " + ", ".join(f"{k}={counts[k]}" for k in sorted(counts)),
        "",
        "Do not treat reviews or deposited datasets as primary experimental articles. For a review, trace the exact compound/value to the cited original paper. For a dataset, preserve the dataset identifier and assay metadata and do not invent a paper citation.",
        "",
        "| # | ChEMBL document | Type | Citation/title | Affected scope | Required action |",
        "|---:|---|---|---|---|---|",
    ]
    for index, record in enumerate(records, 1):
        title = record["title"] or "Bibliography unresolved in local freeze"
        if record["direct_url"]:
            citation = f"[{title}]({record['direct_url']})"
        else:
            citation = f"[{title}]({record['chembl_url']})"
        doc_link = f"[{record['document_chembl_id']}]({record['chembl_url']})"
        scope = f"{record['pairs']} / {record['targets']} / {record['ligands']}"
        md.append(
            f"| {index} | {doc_link} | {record['source_type']} | {citation} | {scope} | {record['required_action']} |"
        )
    md.extend(
        [
            "",
            "## Completion rule",
            "",
            "For every affected ligand-target value, record assay endpoint, biochemical/cellular format, species, protein construct/domain boundaries, wild type or mutation, cofactors/substrate, incubation time, readout, replicate/statistical information, exact compound identity, relation qualifier, unit conversion, and whether the reported value is directly comparable to the frozen threshold. Record a page/table/figure/SI locator for every decision.",
            "",
            "## Important limitation",
            "",
            "This list is the decision-changing first pass, not all historical evidence. The complete audit table contains additional non-maximal documents. Review them when a decisive source is excluded, when a maximum is inconsistent with the paper, or when robust aggregation/sensitivity analysis requires all measurements.",
        ]
    )
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"wrote {OUT_CSV} ({len(records)} sources)")
    print(f"wrote {OUT_FULL_CSV} ({len(full_records)} sources)")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
