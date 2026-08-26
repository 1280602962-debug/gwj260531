# JCIM Submission Remediation Plan

This plan separates repairs that can be made with the frozen repository from claims that require new data. Completed docking calculations do not by themselves make the study submission-ready.

## P0: Scientific Blockers

1. **High-confidence activity view — partially completed.** The 2026-08-26 current-ChEMBL robustness view retains human SINGLE PROTEIN targets, confidence≥8, exact quantitative relations, allowed endpoints, no validity flag, and no potential-duplicate flag. It preserves 352/352 scored labels and records assay/document identifiers plus exclusion reasons. Remaining work is construct/mutation-context review and, if required, document-level assay harmonization. It must remain a dated post-hoc sensitivity rather than silently replacing the frozen analysis.
2. **Complete-case and provenance audit — completed descriptively.** Table S37 reports usable-pChEMBL A-only-measured, B-only-measured, and both-target map coverage, plus class-wise source-document concentration from the dated high-confidence records. Only 14.5%–34.0% of the usable-value union is measured at both targets; the PIK3CA/mTOR neither class is a four-ligand, single-document sample. This audit cannot recover labels for unmeasured structures or substitute for document-blocked/external validation.
3. **Put the core contrast on comparable scales.** Separate negative-class choice from score aggregation. Recompute Dual-versus-neither with pocket-specific, mean, and worst-pocket scores and report the corresponding directional arms. Descriptive bootstrap intervals for differences may be reported, but different negative sets prevent a paired significance interpretation.
4. **Add a genuinely external test or lower the claim.** The unused-pool resample is not external validation. A time-split ChEMBL harvest, database-external structure-matched set, or blinded experimental panel is required for an external-validation claim.

## P1: High-Value Strengthening

1. **Completed descriptively.** Table S38 reports class-wise MW, heavy atoms, cLogP, TPSA, formal charge, rotatable bonds, scaffold count, singleton-scaffold fraction, and nearest-dual Tanimoto; Tables S35 and S37 report measurement count and source-document concentration. These are confounding diagnostics, not causal adjustments.
2. **Partially completed.** Top-1/top-3/all-deposited RMSD was topology-checked and recomputed for 4EY7, 4BDS, and 2WXF; deposited symmetry-aware results already give the corresponding ranks for 4L23 and 4JT6. 4BDS and 4JT6 fail at top-1 but pass within top-3. Historical EGFR/HER2 summaries remain, but their cognate pose artifacts are absent from the repository, so top-3 cannot be independently recomputed. Best-of-nine remains a search-coverage rather than score-ranking check.
3. **Completed for the main directional panels.** Table S27 reports failed ligands by class and property, arm-available AUROC, and deterministic rank-extreme bounds. Failures are concentrated among large/flexible ligands and are not missing at random. The stress test preserves current pair-level conclusions but does not extend applicability to unsupported chemistry.
4. Define fixed train/test files and a machine-readable evaluation contract if DualFourClass-Bench is released as a reusable benchmark.
5. Add at least one non-kinase-family target pair beyond AChE/BChE before making target-general statements.

## P2: Submission and Reproducibility

1. Publish a versioned Zenodo archive containing the exact panels, assay provenance, receptor files, boxes, scores, retained poses required by the paper, environment lock, and checksums.
2. **Partially completed.** The zero-docking analysis was rerun in a separate clean copy, cross-platform UTF-8 and headless plotting failures were repaired, and core statistical outputs were reproduced byte-for-byte across two independent process runs. A repository-wide frozen checksum manifest and CI job are still required for release-grade verification.
3. Replace mutable repository links in Data Availability with a release tag and DOI.
4. **Completed for the current zero-docking reanalysis.** `requirements-analysis.txt` pins NumPy 2.5.2, pandas 3.0.5, Matplotlib 3.11.1, SciPy 1.18.1, scikit-learn 1.9.0, RDKit 2026.3.5, Meeko 0.7.1, and gemmi 0.7.5. Canonical secondary tables were regenerated in that environment. Bootstrap sub-seeds now use stable SHA-256 offsets rather than Python's process-randomized `hash()`. The historical docking engine and preparation environment remains separately documented in `ENV_PIN.md` and was not rerun.

## Claims Permitted before P0 Completion

- The four-state formulation is a useful evaluation design.
- On the frozen EGFR/HER2 case, Dual-versus-neither gave a more favorable descriptive result than directional selective-hard-negative evaluation under both Vina and one GNINA protocol.
- Performance depended on ligand composition and receptor realization in these panels.

## Claims Not Permitted before P0 Completion

- Target-general reliability, robustness, or systematic overestimation.
- Assay-harmonized experimental ground truth.
- External validation or prospective utility.
- A comprehensive or representative dual-target docking benchmark.
