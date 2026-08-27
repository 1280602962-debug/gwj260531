# JCIM Submission Remediation Plan

This plan separates repairs that can be made with the frozen repository from claims that require new data. Completed docking calculations do not by themselves make the study submission-ready.

## P0: Scientific Blockers

1. **High-confidence activity view — partially completed.** The 2026-08-26 current-ChEMBL robustness view retains human SINGLE PROTEIN targets, confidence≥8, exact quantitative relations, allowed endpoints, no validity flag, and no potential-duplicate flag. It preserves 352/352 scored labels. A metadata include/exclude pass covers 186 priority ligands (179/7/0); construct/mutation remain unknown. It must remain a dated post-hoc sensitivity rather than silently replacing the frozen analysis.
2. **Complete-case, provenance, and document-blocked analyses — completed on frozen scores.** Table S37 reports usable-pChEMBL coverage and source-document concentration. Tables S39–S40 report document-blocked CV: EGFR/HER2 Dual versus B-only remains 0.430; PIK3CA/mTOR Dual versus B-only cannot be stably estimated. This does not recover labels for unmeasured structures.
3. **Put the core contrast on comparable scales.** Separate negative-class choice from score aggregation. Recompute Dual-versus-neither with pocket-specific, mean, and worst-pocket scores and report the corresponding directional arms. Descriptive bootstrap intervals for differences may be reported, but different negative sets prevent a paired significance interpretation.
4. **Add a genuinely external test or lower the claim.** The unused-pool resample is not external validation. A pre-frozen 2018 literature-year split on the scored panels yielded zero pairs with dual/A-only/B-only each n≥10 and is not packaged as external validation. BindingDB REST counts (Table S43) remain historical supply. A BindingDB-native 202608 archive rebuild (Tables S48–S49) applied literature, structure, and ECFP4 < 0.70 filters and yielded **zero pairs** meeting the pre-frozen primary external gate; remaining n are upper bounds; the slice was not docked and is not packaged as external evaluation. The manuscript claim stays an internal four-pair formulation audit.

## P1: High-Value Strengthening

1. **Completed descriptively.** Table S38 reports class-wise MW, heavy atoms, cLogP, TPSA, formal charge, rotatable bonds, scaffold count, singleton-scaffold fraction, and nearest-dual Tanimoto; Tables S35 and S37 report measurement count and source-document concentration. These are confounding diagnostics, not causal adjustments.
2. **Partially completed.** Top-1/top-3/all-deposited RMSD was topology-checked for 4EY7, 4BDS, 2WXF, and reconstructed EGFR/HER2 QC. 4BDS and reconstructed 3POZ fail at top-1; 3POZ also fails top-3 but recovers best-of-9 (0.760 Å). 3RCD reconstructed QC passes top-1 (1.855 Å). Best-of-nine remains a search-coverage rather than score-ranking check.
3. **Completed for the main directional panels.** Table S27 reports failed ligands by class and property, arm-available AUROC, and deterministic rank-extreme bounds. Failures are concentrated among large/flexible ligands and are not missing at random. The stress test preserves current pair-level conclusions but does not extend applicability to unsupported chemistry.
4. Define fixed train/test files and a machine-readable evaluation contract if DualFourClass-Bench is released as a reusable benchmark. **Partially completed:** `DUALFOURCLASS_EVALUATION_CONTRACT_v1.json` freezes the estimands, pairs, seeds, and table paths. Tagged Zenodo files remain a maintainer action.
5. Add at least one non-kinase-family target pair beyond AChE/BChE before making target-general statements. A θ = 6.0 census found 17 unique pairs with Dual/A-only/B-only n ≥ 10, including PPI pairs. MCL1/Bcl-xL was docked only as an LC6-gate-fail applicability stress-test and is **not** a fifth Table 2 pair.

## P2: Submission and Reproducibility

1. Publish a versioned Zenodo archive containing the exact panels, assay provenance, receptor files, boxes, scores, retained poses required by the paper, environment lock, and checksums.
2. **Partially completed.** The zero-docking analysis was rerun in a separate clean copy, cross-platform UTF-8 and headless plotting failures were repaired, and core statistical outputs were reproduced byte-for-byte across two independent process runs. `REVISION_CHECKSUM_MANIFEST_v1.csv` and `.github/workflows/revision-validate.yml` now pin manuscript-facing tables. A GitHub Release and Zenodo DOI remain maintainer actions from a tagged snapshot. See `docs/JCIM_SUBMISSION_FIT_ASSESSMENT.md` for the 2023–2026 JCIM bar judgment.
3. Replace mutable repository links in Data Availability with a release tag and DOI.
4. **Completed for the current zero-docking reanalysis.** `requirements-analysis.txt` pins NumPy 2.5.2, pandas 3.0.5, Matplotlib 3.11.1, SciPy 1.18.1, scikit-learn 1.9.0, RDKit 2026.3.5, Meeko 0.7.1, and gemmi 0.7.5. Canonical secondary tables were regenerated in that environment. Bootstrap sub-seeds now use stable SHA-256 offsets rather than Python's process-randomized `hash()`. The historical docking engine and preparation environment remains separately documented in `ENV_PIN.md` and was not rerun.

## Claims Permitted before P0 Completion

- The four-state formulation is a useful evaluation design.
- On the frozen EGFR/HER2 case, Dual-versus-neither gave a more favorable descriptive result than directional selective-hard-negative evaluation under both Vina and one GNINA protocol.
- Performance depended on ligand composition and receptor realization in these panels.

## Claims Not Permitted before P0 Completion

- Target-general reliability, robustness, or systematic overestimation.
- Assay-harmonized experimental ground truth.
- External validation or prospective utility. The 2018 time split is a negative result on sample supply, not external validation.
- A comprehensive or representative dual-target docking benchmark.
