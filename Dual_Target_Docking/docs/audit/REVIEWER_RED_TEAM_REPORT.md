# REVIEWER RED-TEAM REPORT

Date: 2026-09-19  
Standard: JCIM computational / cheminformatics / reproducibility review.

Legend: CURRENTLY ANSWERED | ANSWERABLE_AFTER_DOCUMENTATION | REQUIRES_ANALYSIS | TRUE_LIMITATION | POTENTIAL_FATAL_ISSUE

## Reviewer 1 — computational chemistry / docking

1. How do you know the selected receptor conformation is representative?  
   **TRUE_LIMITATION.** One primary holo per target (plus three PM swaps). No ensemble / MD.

2. Were ligand protonation and tautomer states standardized?  
   **TRUE_LIMITATION** for 7/8 pairs (PDBQT not deposited). EGFR uniform: largest fragment, ETKDGv3, MMFF, Meeko — **CURRENTLY ANSWERED** as a campaign, not as a tautomer enumeration.

3. Why is PM run at E=16 while other pairs use E=8?  
   **CURRENTLY ANSWERED.** 4JT6 failed the cognate-coverage gate at E=8. Sensitivity E=8 exists and is not Table 2.

4. Were waters / metals / peptides treated consistently?  
   **TRUE_LIMITATION.** 9V8H retains peptide (confirmed in PDBQT). Waters absent from PDBQT. Removal protocol NOT_RECOVERABLE.

5. Is 4EY7 docked as a dimer? Does that change the box?  
   **ANSWERABLE_AFTER_DOCUMENTATION.** Dimer AB is in the PDBQT; box is the E20 cognate AABB. No monomer control.

6. TYK2: JH1 vs JH2?  
   **CURRENTLY ANSWERED** in yaml (JH1 / 3LXP / IZA). Sequence-level re-annotation was not repeated in this audit.

7. Cognate top-1 RMSD > 2 Å on several slots — is the site wrong?  
   **CURRENTLY ANSWERED** if worded as coverage: best-of-saved < 2 Å. **OVERSTATED** if used as pose-ranking validation.

8. Are GNINA results independent pose generation?  
   **CURRENTLY ANSWERED** for three pairs. RTM/CNN are same-pose. Must not be collapsed.

9. Was the EGFR box hydrogen-inclusive?  
   **CURRENTLY ANSWERED** for current files: corrected JSON declares heavy-atom AABB. Legacy files exist but match coordinates here.

10. Can the original docking be reconstructed from deposited inputs?  
    **TRUE_LIMITATION** except EGFR uniform. See REPRODUCIBILITY_BOUNDARY.

11. Why not re-prepare receptors under one Meeko protocol now?  
    **CURRENTLY ANSWERED** as a freeze rule: re-preparation would create a new experiment, not recover the old one.

12. Is rigid-receptor docking an appropriate primary engine for dual-target selectivity?  
    **TRUE_LIMITATION.** That is the study object, not a hidden assumption. Must stay in the claim scope.

13. Did timeout ligands get imputed?  
    **CURRENTLY ANSWERED.** EGFR EH40_31 skipped at 600 s; n_B_only 32→31. Not zero-filled.

14. Same ligand, two pockets — complete cases only?  
    **CURRENTLY ANSWERED** for primary. AChE five-seed available-case is not complete vs primary n=95.

15. Could crystal cognates bias the box away from the dual-ligand ensemble?  
    **TRUE_LIMITATION.** One cognate-defined box per slot. No ligand-ensemble box.

## Reviewer 2 — cheminformatics / statistics

1. Was the pair roster influenced by observed docking performance?  
   **CURRENTLY ANSWERED** at the file level: **NO_OUTCOME_LEAKAGE** in the eligibility CSV. Informal discussion NOT_RECOVERABLE.

2. Are eight pairs independent observations?  
   **TRUE_LIMITATION.** Shared JAK1 and PPARA. Must not be meta-analyzed as n=8.

3. Was descriptor selection performed on the evaluation set?  
   **CURRENTLY ANSWERED.** Full-panel best = descriptive. Nested scaffold-CV = predictive control. If SI/main text slips, that is P1.

4. Feature-scaling leakage in ECFP?  
   **CURRENTLY ANSWERED.** Primary ECFP is unscaled; scaler arm is sensitivity.

5. Scaffold split vs GroupKFold?  
   **CURRENTLY ANSWERED.** GroupKFold on Murcko groups; not a single hold-out scaffold split.

6. AUROC tie handling?  
   **CURRENTLY ANSWERED.** Mann–Whitney mid-ranks (0.5 on ties). Independent replay matches.

7. Is the point estimate the bootstrap mean?  
   **CURRENTLY ANSWERED.** No. Original sample; percentile CI.

8. Overlapping CIs treated as a test?  
   **ANSWERABLE_AFTER_DOCUMENTATION.** Text is mostly careful. Must not convert overlapping CIs into “no difference” claims.

9. Multiple comparisons across 8×2 directions?  
   **TRUE_LIMITATION.** Pointwise CIs. No family-wise correction. Do not write “significant after correction.”

10. Are seed effects confounded by changing membership?  
    **CURRENTLY ANSWERED** if Figure 5C/caption use the two flags. **P1** while plotted_values still uses fuzzy `comparable`.

11. Is the negative-class effect chemistry rather than docking?  
    **REQUIRES_ANALYSIS** for a full decomposition; ECFP already shows chemistry carries class signal. Fixed-score Δ holds the score channel fixed, so it is not only a chemistry-composition story, but chemistry still differs by class.

12. Max-pChEMBL bias?  
    **CURRENTLY ANSWERED** as a named sensitivity. Not re-aggregated from raw dump in this audit.

13. Train/test leakage in ECFP?  
    **CURRENTLY ANSWERED** at the fold-assignment file level. Same ligands across model types share folds. No row duplication found in the OOF table keys.

14. Detectable-effect vs post-hoc power?  
    **CURRENTLY ANSWERED** in SI/main if wording is kept. Calling it observed power would be P1.

15. Small n / PM neither n=4?  
    **TRUE_LIMITATION.** Already annotated. Must stay visible.

## Reviewer 3 — data provenance / reproducibility

1. Can original docking be reconstructed?  
   **TRUE_LIMITATION.** See boundary. Not a fatal issue if the paper states the boundary.

2. Are Track-B activities assay-adjudicated?  
   **CURRENTLY ANSWERED** if three kinds are kept. Claiming uniformity would be P1.

3. Receptor protonation / pH / loops?  
   **TRUE_LIMITATION / NOT_RECOVERABLE.**

4. Which EGFR scores are current?  
   **CURRENTLY ANSWERED** in code/tables: uniform RDKit/Meeko seed 20260727. Historical 0.3237 must not reappear as current.

5. Why do root docs still say 0.3237?  
   **P1.** Freeze-check memos were not archived. Reviewer can find them.

6. Why does SI S6 prose contradict SI S6 table?  
   **P1 / POTENTIAL_FATAL_ISSUE for writing quality.** Table is correct; prose is pre-uniform leftover. A reviewer will treat this as the authors not reading their own SI.

7. Why do manuscript cluster CIs disagree with SI S4?  
   **P1.** Same substitution event; main text not updated.

8. Two score-master paths?  
   **ANSWERABLE_AFTER_DOCUMENTATION.** Currently identical. Still a duplicated authority.

9. JAK1/TYK2 document cluster?  
   **CURRENTLY ANSWERED** as `unresolved_mapping_unavailable`. Must not copy an old CI.

10. External validation failed?  
    **CURRENTLY ANSWERED** only as 0/8 eligibility. “Failed validation” would be false.

11. Personal absolute paths in ENV_PIN / recovered PIK3CA script?  
    **ANSWERABLE_AFTER_DOCUMENTATION.** Historical snapshots.

12. Is analysis RDKit the docking RDKit?  
    **CURRENTLY ANSWERED** if kept split. 2026.03.5 vs 2026.3.1.

13. Who writes canonical files?  
    **CURRENTLY ANSWERED.** Named writers; promote path exists.

14. Hidden sample exclusions?  
    **CURRENTLY ANSWERED** for EH120_059, AB_087, EH40_31, AChE five-seed n=88. Membership registry now lists IDs.

15. Can a third party clean-rebuild Table 2 without docking?  
    **CURRENTLY ANSWERED** in principle from frozen scores + pins. This audit did the arithmetic independently. A full temp-dir rebuild was not re-run here.

## Fatal vs not

Nothing found that falsifies the deposited Table 2 numbers.  
The closest **writing-level fatal** risk is SI S6 prose vs table, plus stale cluster CIs in the main text. Those are fixable in a remediation writing pass without new experiments.
