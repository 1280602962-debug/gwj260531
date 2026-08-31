#!/usr/bin/env python3
"""Inventory and identity-check downloaded source PDFs against the audit queue."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from difflib import SequenceMatcher
from pathlib import Path

from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
SOURCE_LIST = TABLES / "assay_source_reading_list_v1.csv"
OUTPUT = TABLES / "downloaded_source_pdf_inventory_v1.csv"

CHEMBL_RE = re.compile(r"CHEMBL\d+", re.I)
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.I)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_pdf(path: Path) -> dict[str, str]:
    result = {
        "pdf_valid": "0",
        "encrypted": "",
        "n_pages": "",
        "text_chars": "0",
        "metadata_title": "",
        "extracted_dois": "",
        "text_extract_status": "not_attempted",
        "error": "",
        "text": "",
    }
    try:
        reader = PdfReader(str(path), strict=False)
        result["pdf_valid"] = "1"
        result["encrypted"] = "1" if reader.is_encrypted else "0"
        result["n_pages"] = str(len(reader.pages))
        metadata = reader.metadata or {}
        result["metadata_title"] = str(metadata.get("/Title", "") or "").strip()
        page_text = []
        failed_pages = 0
        for page in reader.pages:
            try:
                page_text.append(page.extract_text() or "")
            except Exception:
                page_text.append("")
                failed_pages += 1
        text = "\n".join(page_text)
        result["text"] = text
        result["text_chars"] = str(len(text))
        result["extracted_dois"] = ";".join(sorted({d.rstrip(".,;)") for d in DOI_RE.findall(text)}))
        if failed_pages:
            result["text_extract_status"] = f"partial:{failed_pages}_pages_failed"
        elif len(text.strip()) < 500:
            result["text_extract_status"] = "image_or_low_text_pdf"
        else:
            result["text_extract_status"] = "ok"
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["text_extract_status"] = "failed"
    return result


def identity_status(expected: dict[str, str], extracted: dict[str, str]) -> tuple[str, str]:
    text_norm = norm(extracted["text"][:100000])
    expected_doi = expected.get("doi", "").lower().strip()
    found_dois = extracted["extracted_dois"].lower()
    patent = norm(expected.get("patent_id", ""))
    title = expected.get("title", "")
    title_similarity = 0.0
    if title:
        candidate = extracted["metadata_title"] or extracted["text"][:3000]
        title_similarity = SequenceMatcher(None, norm(title), norm(candidate)).ratio()
    if extracted["pdf_valid"] != "1":
        return "invalid_or_non_pdf", f"title_similarity={title_similarity:.3f}"
    if expected_doi and expected_doi in found_dois:
        return "exact_doi_match", f"title_similarity={title_similarity:.3f}"
    if patent and patent in text_norm:
        return "exact_patent_match", f"title_similarity={title_similarity:.3f}"
    if title_similarity >= 0.60:
        return "probable_title_match", f"title_similarity={title_similarity:.3f}"
    if expected.get("source_type") == "deposited_dataset":
        return "dataset_attachment_manual_check", f"title_similarity={title_similarity:.3f}"
    return "manual_identity_check", f"title_similarity={title_similarity:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_dir", type=Path)
    args = parser.parse_args()

    expected_rows = read_csv(SOURCE_LIST)
    expected = {row["document_chembl_id"].upper(): row for row in expected_rows}
    files_by_id: dict[str, Path] = {}
    for path in sorted(args.source_dir.iterdir()):
        if not path.is_file():
            continue
        match = CHEMBL_RE.search(path.name)
        if match:
            files_by_id[match.group(0).upper()] = path

    output_rows = []
    hashes: dict[str, list[str]] = {}
    for doc_id, row in expected.items():
        path = files_by_id.get(doc_id)
        if path is None:
            output_rows.append(
                {
                    "document_chembl_id": doc_id,
                    "source_type": row["source_type"],
                    "download_status": "missing",
                    "filename": "",
                    "bytes": "",
                    "sha256": "",
                    "duplicate_sha256_group": "",
                    "pdf_valid": "",
                    "encrypted": "",
                    "n_pages": "",
                    "text_chars": "",
                    "text_extract_status": "",
                    "expected_title": row["title"],
                    "metadata_title": "",
                    "expected_doi": row["doi"],
                    "extracted_dois": "",
                    "identity_status": "missing",
                    "identity_evidence": "",
                    "assay_terms_present": "",
                    "construct_terms_present": "",
                    "pairs": row["pairs"],
                    "ligands": row["ligands"],
                    "targets": row["targets"],
                    "error": "",
                }
            )
            continue
        digest = sha256(path)
        hashes.setdefault(digest, []).append(doc_id)
        extracted = extract_pdf(path)
        status, evidence = identity_status(row, extracted)
        lowered = extracted["text"].lower()
        assay_terms = [term for term in ("ic50", "ki", "kd", "ec50", "assay", "inhibition") if term in lowered]
        construct_terms = [
            term
            for term in ("wild type", "wild-type", "mutant", "mutation", "recombinant", "domain", "construct", "atp")
            if term in lowered
        ]
        output_rows.append(
            {
                "document_chembl_id": doc_id,
                "source_type": row["source_type"],
                "download_status": "downloaded",
                "filename": path.name,
                "bytes": str(path.stat().st_size),
                "sha256": digest,
                "duplicate_sha256_group": "",
                "pdf_valid": extracted["pdf_valid"],
                "encrypted": extracted["encrypted"],
                "n_pages": extracted["n_pages"],
                "text_chars": extracted["text_chars"],
                "text_extract_status": extracted["text_extract_status"],
                "expected_title": row["title"],
                "metadata_title": extracted["metadata_title"],
                "expected_doi": row["doi"],
                "extracted_dois": extracted["extracted_dois"],
                "identity_status": status,
                "identity_evidence": evidence,
                "assay_terms_present": ";".join(assay_terms),
                "construct_terms_present": ";".join(construct_terms),
                "pairs": row["pairs"],
                "ligands": row["ligands"],
                "targets": row["targets"],
                "error": extracted["error"],
            }
        )

    duplicates = {key: value for key, value in hashes.items() if len(value) > 1}
    for row in output_rows:
        if row["sha256"] in duplicates:
            row["duplicate_sha256_group"] = ";".join(duplicates[row["sha256"]])

    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"wrote {OUTPUT}")
    print(f"expected={len(expected_rows)} downloaded={sum(r['download_status'] == 'downloaded' for r in output_rows)}")
    for key in ("exact_doi_match", "exact_patent_match", "probable_title_match", "manual_identity_check", "invalid_or_non_pdf", "missing"):
        print(f"{key}={sum(r['identity_status'] == key for r in output_rows)}")
    print(f"duplicate_hash_groups={len(duplicates)}")


if __name__ == "__main__":
    main()
