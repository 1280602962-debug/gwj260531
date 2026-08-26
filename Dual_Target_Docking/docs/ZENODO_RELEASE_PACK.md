# Zenodo / GitHub release pack (do this last)

Do **not** mint a DOI while document-blocked, assay-context, or time-split tables may still change.

## Cloud-prepared now

- SHA-256 manifest: `data/jcim_novelty_v0/tables/REVISION_CHECKSUM_MANIFEST_v1.csv`
- CI: `.github/workflows/revision-validate.yml`
- ChEMBL current-view date: 2026-08-26 (high-confidence rebuild)
- Supply-audit freeze date: 2026-07-23
- Environment pin: `requirements-analysis.txt` and `data/jcim_strengthen_t0t1_v0/ENV_PIN.md`

## Local / maintainer actions

1. Finish human assay-context decisions if any labels change, then freeze tables.
2. `git tag -a dualfourclass-jcim-v1.0 -m "Frozen DualFourClass formulation-audit release"`
3. GitHub Release from that tag (not from a moving branch).
4. Upload to Zenodo: panels, labels, assay provenance CSVs, receptors, boxes, scores, cognate artifacts that exist, checksums, `requirements-analysis.txt`.
5. Replace Data Availability URLs with the tag and the Zenodo DOI.

Temporary files such as `__pycache__` are gitignored and must not go into the archive.
