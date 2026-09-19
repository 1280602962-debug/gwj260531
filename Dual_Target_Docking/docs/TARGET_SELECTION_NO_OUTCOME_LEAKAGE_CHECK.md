# Target-selection no-outcome-leakage check

Date: 2026-09-18  
Authority for the eight pairs: `scripts/analysis/analysis_config.py` `PRIMARY_PAIRS`.  
Audit table: `data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv`.

## Columns actually used for inclusion/exclusion

`pair_eligibility_audit_s14_v1.csv` columns: `pair`, `last_gate_reached`, `included_excluded`, `reason`, `evidence_used`.

Token scan of that file: **no** `AUROC`, `summary_min`, `GNINA`, `matched`, `mismatched`, `ECFP`, or `ranking` strings.

Inclusion reasons are activity supply, human holo availability, site interpretability, and protocol compatibility (noncovalent rigid-receptor Vina). EGFR/HER2 is a documented **supply-limited exception** at θ=6.0, not a performance-based rescue.

Excluded pairs fail G3/G4/G5 for covalent cysteine protease, dual-domain ambiguity, antitarget relationship, GPCR/membrane/SLC6 construct issues, or ligand-identity QC (`pair_ligand_identity_qc_v1.csv`). None of those reasons cite docking AUROC.

Track-B yaml header: “Do not edit after seeing AUROCs.” That is a lock statement, not evidence that AUROCs were used to pick the five pairs.

## Downstream files that contain AUROC / summary_min / GNINA / ECFP

These are **outcomes** on the fixed eight pairs, not selection inputs:

- `results/canonical/primary_summary_min.csv`
- `results/canonical/computational_robustness.csv`
- `results/canonical/matched_minus_mismatched.csv`
- `results/canonical/ecfp4_incremental_information.csv`
- `results/canonical/top10_operating_points.csv`

## Broken historical citations (not outcome leakage)

`evidence_used` cites `TIER1_DOCKING_ROSTER_V1.md` and `FEASIBLE_PAIR_LADDER_V1.md`, which are **git history only / historical / unavailable**. Current footnotes must use `pair_eligibility_audit_s14_v1.csv`. That is a provenance gap for the original audit table’s footnotes, not evidence that selection used docking metrics.

## Verdict

Pair inclusion/exclusion as recorded in the current eligibility CSV is **not outcome-driven**. Runtime membership is the hard-coded `PRIMARY_PAIRS` tuple.
