# Revision items: cloud-complete vs local

Wet-lab prospective testing is out of scope. Everything else is listed below. Cloud work uses only frozen scores and public ChEMBL metadata; it does not redock.

| Item | Cloud now | Local required | Status in this revision |
|------|-----------|----------------|-------------------------|
| 1. Document-blocked CV + document-cluster bootstrap | Yes. Same folds for ECFP4, physicochemical, and docking logistic models. | None, unless a pair is marked not stably estimable and a later human wants a different grouping (do not regroup to chase AUROC). | Tables S39–S40 |
| 2. Assay-context audit | Machine extraction, risk flags, leave-one-out influence, empty include/exclude columns. | Read source papers for construct, WT/mutant, and include/exclude. Recompute Table 2 if labels change. | Machine extract done; **local metadata review filled include/exclude** (`ASSAY_CONTEXT_HUMAN_REVIEW_SUMMARY_V1.md`). Construct/mutation remain `unknown` until ChEMBL free-text is readable. No frozen-label flips → Table 2 not recomputed. |
| 3. Freeze time-split protocol | Yes. Cutoffs 2015/2018/2020 frozen in code before AUROC. | None for the protocol. | `TIME_SPLIT_VERDICT.md` |
| 4. Evaluate time-split on already-scored ligands | Yes, using ChEMBL `document.year`. | Docking any new late-era ligand not already in the panels. | Table S41; package as external validation only if ≥2 pairs pass the n≥10 gate |
| 5. BindingDB external set | Native 202608 article/patent TSV rebuild: literature, structure, ECFP4 < 0.70. 0 primary-gate pairs; 0 thin replications. Remaining n are upper bounds (`chembl_document_api_partial`). | Docking only if ≥2 primary pairs had passed (they did not). MCL1 LC6 pose-gold gate if a later session docks that frozen panel. | Tables S48–S49 frozen; **not packaged**; **not docked**. Table S43 kept as REST historical count. |
| 6. EGFR/HER2 cognate poses | Inventory. Crystal/receptor files present; nine-mode cognate PDBQT absent. | Recover original Vina 9-mode files or re-redock as reconstructed QC. | **Reconstructed QC deposited** under `data/egfr_her2_panel40_v0/cognate_qc/` (not historical production). 3RCD passes top-1; 3POZ fails top-1/top-3 but recovers best-of-9 (0.68 Å). |
| 7. Checksum, numeric audit, gitignore, CI | Yes. Native-slice downloader is **not** in CI (BindingDB zips ~180 MB). | GitHub Release + Zenodo DOI (maintainer). | Manifest + workflow; no DOI yet |
| 8. New target pairs | θ = 6.0 census done (Table S44: 17 unique pairs with Dual/A/B n≥10). MCL1/Bcl-xL panel/receptors frozen (Tables S50–S51). | Docking extra pairs only after a pose-gold gate. Do not dock the full MCL1 panel before LC6. | Census + MCL1 freeze; **not docked** |
| 9. Receptor ensemble | Not first priority. | New docking. | Not started |
| 10. Prospective experiment | Impossible in this project. | Wet lab. | Out of scope |
| 11. Property-caliper / AND filter / ligand-only full maps | Yes. Frozen scores + cached maps. | None. | Tables S45–S47; Figure S7 |

Recommended next local session: run the MCL1/Bcl-xL LC6 pose-gold gate on 3WIY/3WIZ before any panel docking; do not regroup document-blocked folds to chase AUROC. Mint Zenodo from a frozen tag before JCIM. Do not re-open BindingDB θ or chemical gates after seeing Tables S48–S49.
