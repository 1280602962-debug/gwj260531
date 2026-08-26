# JCIM 2023–2026 submission-bar assessment (2026-08-26)

This is an internal editorial judgment, not a claim in the manuscript.
It uses the frozen DualFourClass results after the local assay-context
metadata pass and reconstructed EGFR/HER2 cognate QC.

## Verdict

**As a JCIM Research Article that claims a general dual-target docking
benchmark: no.**

**As a tightly scoped evaluation Article whose claim is a four-pair
formulation/failure-mode audit: borderline, and still below the median
empirical scale of recent JCIM docking papers.** Honest claim language
makes desk-rejection less likely than overclaiming, but it does not by
itself meet the last-three-year bar.

Do not submit until a Zenodo DOI exists. Even then, the remaining
scientific gaps below are larger than the remaining writing gaps.

## What recent JCIM docking papers typically bring

Representative 2023–2026 JCIM docking/evaluation papers are larger or
more method-forward than this case panel:

- Schaller et al., JCIM 2024: kinase cross-docking on 589 structures /
  423 ATP-competitive ligands, with automated pipelines.
- Challenge and comparison papers (Polaris/ASAP 2025 and related JCIM
  docking evaluations) use community test sets, pose-success rates, and
  released code/data.
- Method papers still usually report CASF, DUD-E, DEKOIS, LIT-PCBA, or
  similarly large established sets, plus a DOI.

The journal has published negative and formulation-critical work, so a
small, honest audit is not automatically out of scope. Scale, independent
testing, and deposit remain the usual difference between “publishable
evaluation” and “internal methods note.”

## What this manuscript actually contributes

The durable contribution is not a performance ranking. It is that
dual-target docking evidence changes with the negative class and with
receptor realization:

- EGFR/HER2 Dual versus neither 0.756 versus directional `summary_min`
  0.430; independent GNINA pose generation 0.783 versus 0.220.
- ECFP4/property baselines often match or beat docking; adding docking
  to ECFP4 changed CV AUROC by at most 0.020.
- The same PIK3CA receptor swap lowered PIK3CA/mTOR `summary_min` and
  raised PIK3CA/PIK3CB.
- Document-blocked CV left the EGFR/HER2 weak arm at 0.430 and showed
  that PIK3CA/mTOR Dual versus B-only cannot be stably estimated.

Those findings are JCIM-relevant. They are also empirically thin.

## What the new local results change — and what they do not

Completed locally on 2026-08-26:

1. Assay-context metadata review of 186 priority ligands: 179 include /
   7 uncertain / 0 exclude; **zero frozen-label flips**. Construct and
   mutation remain `unknown`. This is **not** assay harmonization and
   does **not** recompute Table 2.
2. Reconstructed EGFR/HER2 cognate QC (not historical production):
   topology-checked 3POZ top-1 9.505 Å / top-3 6.227 Å / best-of-9
   0.760 Å (ranking fails; search coverage passes); 3RCD top-1 1.855 Å
   (passes). This closes a reproducibility hole and strengthens the
   search-coverage versus ranking distinction already made for 4BDS/4JT6.

These items lower a specific reviewer complaint (“missing cognate
artifacts / empty include-exclude columns”). They do not move the paper
across the JCIM bar.

## Remaining gaps versus the 2023–2026 bar

| Gap | Status after local deposit | Effect on JCIM fit |
|-----|----------------------------|--------------------|
| K = 4; three kinase ATP sites; shared PIK3CA | Unchanged | Looks like a case series, not a benchmark |
| All four primary `summary_min` CIs include 0.5 | Unchanged | No pair has a statistically resolved directional effect |
| No docked database-external or time-split test | 2018 split unevaluable; BindingDB not-in-panel supply exists on all four pairs but UniChem/PMID independence was not finished and nothing was docked | Reviewers can reject for missing independent test |
| Assay not paper-harmonized | Metadata pass only | Labels remain operational, not experimental ground truth |
| No Zenodo DOI | Still a moving branch | Data Availability is below current JCIM practice |
| No prospective experiment | Out of scope | Acceptable only if claims stay computational |

## Permitted versus fatal claims at submission

Permitted:

- Four-pair formulation audit of docking-based dual-target recognition.
- Selective hard negatives and confounder/receptor checks are necessary
  evaluation controls in these panels.
- EGFR/HER2 shows a descriptive Dual-versus-neither versus directional
  gap under Vina and one GNINA protocol.

Fatal if stated:

- A general, representative, or comprehensive dual-target docking
  benchmark.
- Target-general reliability or systematic overestimation.
- Assay-harmonized ground truth.
- External validation, including shopping the 2015 AChE/BChE time split.
- That reconstructed 3POZ QC is the original production pose set.

## Practical recommendation

Do **not** send this as a methods/benchmark Research Article claiming
broad docking performance.

A realistic path is a short, tightly titled evaluation Article (current
title is already the right scope) after:

1. minting a Zenodo DOI from a frozen tag;
2. either docking a true BindingDB-independent slice on ≥2 pairs or
   stating in the cover letter that independent docking was impossible
   because supply/literature overlap failed the pre-frozen gate;
3. keeping every CI-includes-0.5 sentence.

Even then, a 2024–2026 JCIM reviewer can still reject for N and for the
unresolved independent-test gap. That risk is scientific, not editorial
wording.

## If not JCIM: journals where this manuscript is already in range

“Already sufficient” here means: the current four-pair formulation audit,
with the present controls and honest limitations, matches what that
journal actually publishes. It does **not** mean guaranteed acceptance.
Mint a Zenodo DOI from a frozen tag before any of these submissions.
Do not widen the claim to a general benchmark.

### First choice (best fit, already sufficient)

1. **Journal of Computer-Aided Molecular Design (JCAMD).**
   This is the closest intellectual home. JCAMD regularly publishes
   docking validation, failure-mode, and reproducibility papers on
   smaller, protocol-defined panels. A 2026 JCAMD review on docking
   failure modes (conditional modeling, ligand-state definition,
   out-of-distribution tests) is the same problem this manuscript
   operationalizes. Keep the current title and limitations. Expect
   questions about K = 4, not a desk-reject for “not a method paper.”

2. **ACS Omega.**
   ACS sister journal to JCIM, with a soundness-and-transparency bar
   rather than a large-benchmark bar. Recent CADD papers there include
   few-target docking validations weaker than this audit. The present
   statistics, SI, and checksums are already above the median ACS Omega
   docking article, provided the cover letter says “four-pair audit”
   and not “JCIM-style benchmark.”

### Also in range (submit as-is after DOI)

3. **Molecular Informatics.** Cheminformatics evaluation and dual-target
   case studies are normal content. Shorter than the JCIM draft is
   fine; do not add new pairs just to look larger.
4. **Journal of Molecular Modeling.** Docking case panels of this size
   are routine. Lower prestige than JCAMD/ACS Omega; highest
   claim-to-venue match if the goal is to place the record, not to
   maximize IF.

### Stretch, not “already fully sufficient”

5. **Journal of Cheminformatics.** Right culture (open evaluation, FAIR,
   negative controls), but many 2024–2025 papers are software, large
   benchmarks, or reusable toolkits. K = 4 without a docked external
   set can still be rejected as too small. Only try after Zenodo and
   only if DualFourClass is framed as a reusable evaluation protocol
   with frozen files, not as a four-target performance ranking.

### Soundness floor (acceptance-likely, lower citation venue)

6. **Scientific Reports**, **PLOS ONE**, or **RSC Advances.**
   These journals ask whether the study is technically sound, not
   whether it is a field-defining benchmark. The current manuscript
   already meets that bar if claims stay computational. Use these if
   the priority is to publish the negative/formulation record quickly.

### Do not send as currently written

- **JCIM / JCTC / Chemical Science / Briefings in Bioinformatics** as a
  general docking or methods paper: empirical scale is below the recent
  median.
- MDPI *Molecules* / *IJMS* / *Pharmaceuticals* as a first choice: they
  would likely accept, but they are a floor, not a match for the
  paper’s actual contribution.

### Cover-letter sentence that matches the viable venues

This is a four-pair, ChEMBL-constrained formulation audit of
docking-based dual-target recognition. It reports failure modes of
negative-class choice, ligand confounding, and receptor realization; it
does not claim a general dual-target docking benchmark or prospective
utility.

## Can this be tried at JCIM first?

Yes, as a **sequential ACS strategy**, not as the venue where the paper is already fully sufficient.

Submit only as an **Article** with the current title. Do not use Application Note (that type is for software, databases, and web servers, ≤5000 words, with the tool name in the title). Do not use Letter (≤~3500 words; this is not a preliminary finding).

The only good reason to try JCIM before JCAMD/ACS Omega is ACS manuscript transfer: a JCIM reject can often move to ACS Omega with the referee reports. That makes one JCIM attempt rational if the authors accept that **rejection is the modal outcome**.

Do not submit without a Zenodo DOI. The expected path is: DOI → JCIM Article → likely reject or heavy revision → transfer to ACS Omega or resubmit to JCAMD. Direct ACS Omega or JCAMD remains the higher-probability first submission.
