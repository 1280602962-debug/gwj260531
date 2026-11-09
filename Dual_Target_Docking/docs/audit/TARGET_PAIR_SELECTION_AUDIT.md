# TARGET PAIR SELECTION AUDIT

Date: 2026-09-19  
Evidence A: `data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv`  
Evidence B: `scripts/analysis/analysis_config.py` `PRIMARY_PAIRS`  
Evidence C: `data/jcim_chembl_universe_v0/tables/universe_pairs_structure_feasibility_v1.csv`  
Prior memo `docs/TARGET_SELECTION_NO_OUTCOME_LEAKAGE_CHECK.md` is not treated as truth.

## Runtime roster

Hard-coded eight pairs:

EGFR/HER2, JAK1/JAK2, JAK1/TYK2, PIK3CA/mTOR, AChE/BChE, F2/F10, PPARG/PPARA, PPARA/PPARD.

This tuple, not the eligibility CSV, is what analysis code actually loops.

## Inclusion reasons in the eligibility CSV

Included pairs reach G5 protocol compatibility, except EGFR/HER2 (`supply-limited exception`).

Reasons cite: conventional noncovalent pockets; human holo; cognate-defined site; Vina protocol feasibility; activity supply.

Token scan of the eligibility CSV: **no** `AUROC`, `summary_min`, `GNINA`, `matched`, `mismatched`, `ECFP`, or `ranking`.

## Exclusions

CTSK/CTSS, CREBBP/BRD4, F2/PRSS1, GPCR/SLC6 pairs: ligand-identity / membrane / covalent / domain-ambiguity gates.

OPRM1/OPRK1 and JAK3/TYK2: G3 after drug-like filter reduced strict selective counts.

None of those reasons cite docking performance.

## EGFR/HER2 exception

CSV reason (design-period language): failed strict 6.5/5.5 minimum selective-class count of 7, but had human holos, a cognate site, and sufficient θ=6.0 dual / A-only / B-only for directional evaluation.

This is a **supply** exception. No AUROC is in the inclusion file.  
Cannot prove from this tree that the exception was written before any docking was seen; the recorded reason does not use outcomes.

## PIK3CA/PIK3CB is not a current primary pair

`universe_pairs_structure_feasibility_v1.csv`:

> FAIL H3: human PIK3CB (P42338) has zero PDB xrefs. Only mouse ortholog holos exist (2Y3A, 4BFR).

That is a **human-holo structural** exclusion, not a performance exclusion.  
Historical PIK3CA/PIK3CB figures are not typeset. WRITING_INDEX forbids presenting it as current.

## Outcome-driven add/drop

No evidence in the eligibility CSV or `PRIMARY_PAIRS` that a pair was removed because AUROC was poor, or added because AUROC was good.

Historical TIER1 / FEASIBLE markdowns cited in `evidence_used` are git-history only. That is a **footnote provenance gap**, not outcome leakage.

Track-B yaml “Do not edit after seeing AUROCs” is a lock statement, not proof that AUROCs were used to pick pairs.

## Verdict

**NO_OUTCOME_LEAKAGE**

Scope: recorded inclusion/exclusion tables and runtime `PRIMARY_PAIRS`.  
Limitation: original design-meeting notes are NOT_RECOVERABLE; cannot reconstruct informal discussion.
