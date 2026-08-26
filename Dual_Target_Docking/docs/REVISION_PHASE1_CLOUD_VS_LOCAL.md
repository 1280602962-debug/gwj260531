# Revision items: cloud-complete vs local

Wet-lab prospective testing is out of scope. Everything else is listed below. Cloud work uses only frozen scores and public ChEMBL metadata; it does not redock.

| Item | Cloud now | Local required | Status in this revision |
|------|-----------|----------------|-------------------------|
| 1. Document-blocked CV + document-cluster bootstrap | Yes. Same folds for ECFP4, physicochemical, and docking logistic models. | None, unless a pair is marked not stably estimable and a later human wants a different grouping (do not regroup to chase AUROC). | Tables S39–S40 |
| 2. Assay-context audit | Machine extraction, risk flags, leave-one-out influence, empty include/exclude columns. | Read source papers for construct, WT/mutant, and include/exclude. Recompute Table 2 if labels change. | Machine extract done; **local metadata review filled include/exclude** (`ASSAY_CONTEXT_HUMAN_REVIEW_SUMMARY_V1.md`). Construct/mutation remain `unknown` until ChEMBL free-text is readable. No frozen-label flips → Table 2 not recomputed. |
| 3. Freeze time-split protocol | Yes. Cutoffs 2015/2018/2020 frozen in code before AUROC. | None for the protocol. | `TIME_SPLIT_VERDICT.md` |
| 4. Evaluate time-split on already-scored ligands | Yes, using ChEMBL `document.year`. | Docking any new late-era ligand not already in the panels. | Table S41; package as external validation only if ≥2 pairs pass the n≥10 gate |
| 5. BindingDB external set | Protocol only. | Accession lock, Ki/Kd/IC50 split, structure and literature dedup, then docking of true new ligands. | SOP; start only if time-split fails |
| 6. EGFR/HER2 cognate poses | Inventory. Crystal/receptor files present; nine-mode cognate PDBQT absent. | Recover original Vina 9-mode files or re-redock as reconstructed QC. | **Reconstructed QC deposited** under `data/egfr_her2_panel40_v0/cognate_qc/` (not historical production). 3RCD passes top-1; 3POZ fails top-1/top-3 but recovers best-of-9 (0.68 Å). |
| 7. Checksum, numeric audit, gitignore, CI | Yes. | GitHub Release + Zenodo DOI (maintainer). | Manifest + workflow; no DOI yet |
| 8. New target pairs | Data-supply audit can be done later from cached maps. | Docking, receptors, boxes. | Not first priority |
| 9. Receptor ensemble | Not first priority. | New docking. | Not started |
| 10. Prospective experiment | Impossible in this project. | Wet lab. | Out of scope |

Recommended next local session: when ChEMBL assay free-text is available, fill `protein_construct` / `wildtype_or_mutant` for uncertain rows; optionally deposit historical EGFR/HER2 cognate poses if recovered from backup (replace reconstructed QC only with SI note).
