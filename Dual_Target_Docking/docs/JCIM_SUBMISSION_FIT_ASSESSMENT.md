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

Zero-docking upgrades added 2026-08-26 (still no wet lab, no new docking):

1. θ = 6.0 four-state census of 49 unique J0 pairs: 17 pairs have
   Dual/A-only/B-only n ≥ 10 (Table S44). Label supply exists beyond K = 4;
   docking evaluation does not.
2. Property-caliper matching and AND-filter operating points on the frozen
   scores (Tables S45–S46). EGFR/HER2 AND precision at the median Dual
   `vina_worst` is 0.298.
3. Ligand-only ECFP4 on the **full ChEMBL maps** of the four pairs (Table S47).
   Dual versus neither remains easier than Dual versus selectives at map
   scale. This is a chemical-label result, not a docking scale-up.

These items make the formulation-audit claim harder to dismiss as an
n ≈ 28 sampling artifact. They do **not** create a general docking
benchmark, dock BindingDB, or exclude 0.5 from Table 2 CIs. JCIM remains
a borderline evaluation Article after Zenodo, not a methods/benchmark lock.

## Remaining gaps versus the 2023–2026 bar

| Gap | Status after local deposit | Effect on JCIM fit |
|-----|----------------------------|--------------------|
| K = 4; three kinase ATP sites; shared PIK3CA | Unchanged | Looks like a case series, not a benchmark |
| All four primary `summary_min` CIs include 0.5 | Unchanged | No pair has a statistically resolved directional effect |
| No docked database-external or time-split test | 2018 split unevaluable; BindingDB-native 202608 slice frozen after literature/structure/ECFP4 < 0.70; **0 pairs pass the primary gate**; remaining n are upper bounds; **not docked**; Table S43 REST counts are historical supply only | Reviewers can still reject for missing *docked* independent test; the independence hole is now a documented supply stop, not an unfinished PMID/UniChem task |
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
2. stating in the cover letter that a BindingDB-native independent slice was frozen, that zero pairs met the pre-frozen primary external gate, and that the slice was therefore not docked;
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
   This is the closest intellectual home. See the 2023–2026 comparator
   note below. Keep the current title and limitations. Expect
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

## JCAMD 2023–2026 comparators (what actually published)

There is **no 2023–2026 JCAMD original article** that is a four-pair
experimental dual-target docking formulation audit (four ligand states,
document-blocked CV, hard-negative vs easy-decoy, receptor realization,
assay-context, frozen files). Dual-target *application* papers exist;
formulation *audits* of docking-based dual-target recognition do not.
Springer search of JCAMD 2023–2026 for “dual target” returns 114 hits,
almost all applied CADD (QSAR + dock + MD + named candidate), not
four-state docking audits.

### Closest JCAMD neighbors (same questions, not the same paper)

1. **Kittelson, Martins, Santos, Celante & Gomes, *J Comput Aided Mol Des*
   40, 137 (2026).**
   Review: *Reproducibility, validation, and failure modes across
   classical and AI-driven molecular docking*
   (https://doi.org/10.1007/s10822-026-00849-8).
   Strongest intellectual neighbor. Treats docking as conditional
   modeling whose interpretability depends on ligand-state definition,
   decoys, OOD tests, and FAIR reporting. DualFourClass is an empirical
   operationalization of one of those failure modes (negative-class /
   Dual vs B-only), not a review. Cite it in the JCAMD cover letter.

2. **Vázquez, García, Llinares, Luque & Herrero, *J Comput Aided Mol Des*
   38, 18 (2024).**
   *On the relevance of query definition in the performance of 3D
   ligand-based virtual screening* (https://doi.org/10.1007/s10822-024-00561-5).
   Same scientific sentence as DualFourClass: **how the query/class is
   defined changes the performance number.** They use DUD-E+ and
   ligand-based 3D VS, not dual-target docking. Larger and cleaner as a
   public-set evaluation; DualFourClass is narrower but unique on
   four-state dual-target docking.

3. **Ugurlu & He, *J Comput Aided Mol Des* 40, 45 (2026).**
   *Prodrug-ML: prodrug-likeness prediction via machine learning on
   sampled negative decoys* (https://doi.org/10.1007/s10822-025-00725-x).
   Same negative-class lesson (easy decoys inflate scores). ML/ADMET,
   not docking, not dual-target. DualFourClass is the docking analogue
   of this argument.

4. **Bozkır, İbişoğlu, Güler & Bozkır, *J Comput Aided Mol Des* 40, 54
   (2026).**
   *Computational prioritization of multi-target inhibitors: explainable
   QSAR and docking-based discovery of dual AChE/BACE1 chemotypes*
   (https://doi.org/10.1007/s10822-025-00757-3).
   Closest **applied dual-target** neighbor on JCAMD: ChEMBL QSAR,
   scaffold nested CV, bootstrap CIs, docking of named chemotypes,
   ADMET. Dual recognition is QSAR + two independent docking scores,
   not Dual vs B-only. DualFourClass is stronger as docking evidence
   and weaker as a product (no named lead).

5. **Andola & Doble, *J Comput Aided Mol Des* 40, 77 (2026).**
   *Design of novel PI3Kα and PI3Kγ inhibitors* using pharmacophore,
   protein–ligand contacts, and ML (https://doi.org/10.1007/s10822-025-00734-w).
   Same kinase family as PIK3CA/PIK3CB, but inhibitor *design*, not a
   four-state docking audit. Do not cite as a DualFourClass analogue.

### Head-to-head (honest)

| Axis | Typical JCAMD dual-target CADD (e.g. Bozkır 2026) | Kittelson 2026 review | Vázquez 2024 query-definition | DualFourClass |
|---|---|---|---|---|
| Type | Applied pipeline | Review | Ligand-based VS evaluation | Experimental docking audit |
| Dual-target claim | Dual QSAR + two dock scores | General docking | Not dual-target | Four ligand states, Dual vs B-only |
| Negative-class control | Usually implicit / SMOTE | Discussed | Query definition is the paper | Easy decoy vs hard B-only, measured |
| Scale | One pair + named candidates | Broad literature | DUD-E+ | K = 4 pairs |
| Outcome | Named chemotypes | Recommendations | Definition changes AUC | Dual vs B-only CIs include 0.5 |
| Fit to JCAMD | High (common type) | High | High | High if titled as an audit |

DualFourClass is **not** a better version of the typical JCAMD dual-target
hit-finding paper. It is a better version of the *evaluation* papers
JCAMD also publishes (Vázquez 2024; Kittelson 2026). Against those
evaluation papers it is original on dual-target docking and smaller in N.
Against Bozkır-style applied dual-target papers it is more rigorous and
more negative, which is an advantage at JCAMD only if the title stays
an audit.

Do not tell JCAMD editors this is the first dual-target docking paper.
Tell them it is, to our knowledge, the first **four-pair experimental
audit of docking-based dual-target recognition under a frozen four-state
label protocol**. That claim is supportable from the 2023–2026 JCAMD
record; a claim of being the first dual-target CADD paper is not.
