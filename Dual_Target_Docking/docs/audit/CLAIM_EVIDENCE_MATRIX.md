# CLAIM–EVIDENCE MATRIX

Date: 2026-09-19  
Evidence grade: A frozen input / B executable independent replay / C protocol / D log / E generated table / F markdown / G manuscript.

| # | Claim | Supporting analyses | Contradictory / limiting evidence | Scope | Strength | Allowed wording | Class |
|---|--------|---------------------|-----------------------------------|-------|----------|-----------------|-------|
| 1 | Docking can/cannot distinguish dual from single-target-active ligands | Table 2 directional AUROCs; independent replay MATCH | Pair-specific: PPARG/PPARA above 0.5; EGFR and F2/F10 below 0.5; 5 pairs cross 0.5 | 8 pairs, θ=6.0, Vina mode-1 | Pair-specific rank AUROC only | “On this panel, directional AUROC ranged …” | SUPPORTED_WITH_LIMIT |
| 2 | Negative-class definition changes apparent performance | Fixed-score Δ Table S4; Fig 2A; EGFR pocket A 0.446 [0.259, 0.632]; JAK1/TYK2 0.444 [0.261, 0.630] | Other pairs’ Δ CIs include 0; PM neither n=4 | Same score channel, different negatives | Strong on two directions; not universal | “Changing the control class altered AUROC on some pairs/directions” | SUPPORTED |
| 3 | Chemistry already separates the experimental classes | ECFP4 OOF; descriptor screen | Docking-only logistic OOF is not rank AUROC; must not be compared as if it were Table 2 | Same ligands, scaffold GroupKFold | Strong that fingerprints carry class signal | “Receptor-free ECFP4 carried substantial class information on several directions” | SUPPORTED |
| 4 | Docking adds little incremental information beyond ECFP4 | ECFP4+docking − ECFP4; max \|Δ\| = 0.0234 at PPARA/PPARD D vs B | Incremental is OOF logistic, not rank AUROC; 0.0112 is superseded | 16 arms | Descriptive incremental | “Under this model, adding the matched docking score did not produce a consistent AUROC gain” | SUPPORTED_WITH_LIMIT |
| 5 | Best descriptor is a predictive baseline | Nested CV exists | Full-panel best is selected on the evaluation set | If manuscript treats full-panel best as predictive | Optimism | Must call it descriptive; use nested OOF as the predictive number | DESCRIPTIVE_ONLY for full-panel best; SUPPORTED_WITH_LIMIT for nested |
| 6 | Matched pocket is better than mismatched | MM Δ; AChE 0.177 [0.053, 0.291]; EGFR 0.107 [0.006, 0.220] | Other 6 main CIs include 0; all 7 holdouts include 0; homologous pockets are imperfect controls | Main vs holdout | Pair-specific | “AChE and EGFR main-panel intervals excluded 0” | SUPPORTED_WITH_LIMIT |
| 7 | “Only AChE excludes 0” | — | Contradicted by current EGFR MM and by SI Table S6 cells | Current uniform EGFR | False | Do not use | OVERSTATED / STALE |
| 8 | PIK3CA/mTOR is more reproducible | Five-seed range on n=48; E=16 | Neither n=4; receptor swap moves summary_min; top-1 RMSD on 4JT6 is 7.118 Å | One pair | Weak as a general claim | “PM used E=16 because 4JT6 failed the E=8 coverage gate” is supported; “more reproducible” is not a primary claim | DESCRIPTIVE_ONLY |
| 9 | Receptor sensitivity does not replace primary | 4JPS/5DXT/4JSX keep the other pocket from the master | Membership unchanged (n=48); scores change | PM48 only | Correct as sensitivity | “Sensitivity, not a second Table 2” | SUPPORTED |
| 10 | Results are robust | Five-seed; label; cluster | AChE membership changes; EGFR GNINA n=89; cluster document JAK1/TYK2 unresolved | Named sensitivities only | Partial | Name the sensitivity; do not say “robust in general” | SUPPORTED_WITH_LIMIT |
| 11 | Generalizable / validated externally | Eligibility CSV 0/8 packaged | No external docking was run | BindingDB/PubChem filter | Eligibility only | “0/8 pairs met the external-eligibility gate; no external docking was performed” | SUPPORTED as eligibility; UNSUPPORTED as validation |
| 12 | Eight pairs are independent replicates | — | Shared JAK1 and PPARA receptors; shared chemical series possible | Multiplicity | False if worded as n=8 independent | “Eight pair-specific evaluations; shared targets are not independent” | OVERSTATED if independent |
| 13 | GNINA confirms the design effect | Independent GNINA on 3 pairs | Same-pose CNN/RTM are not independent; EGFR GNINA n=89 | 3 pairs | Design-difference recurrence only | “Independent pose generation on three pairs; not an eight-pair engine ranking” | SUPPORTED_WITH_LIMIT |
| 14 | Cognate redocking validates the site | CalcRMS table; EGFR 1.019 Å; HER2 1.947 Å | Several top-1 ≥ 2 Å with best-of-saved < 2 Å | 14 slots | Search coverage, not ranking | “Lowest saved RMSD is coverage, not top-1 ranking quality” | SUPPORTED_WITH_LIMIT |
| 15 | All ligands/receptors were prepared uniformly | EGFR uniform campaign | Receptor prep NOT_RECOVERABLE; 7/8 ligand PDBQT absent; LigPrep is historical | Methods | False if absolute | State deposited PDBQT + recovered scripts + EGFR CASE 2 | UNSUPPORTED if uniform |
| 16 | Detectable-effect is observed power | Simulation CSV; n_mc=1000 | Manuscript/SI already say it is not observed power | Class sizes from master | Correct if named | “detectable-effect simulation” | SUPPORTED |

## Writing errors that break claim–evidence alignment

- Main-text cluster CIs for EGFR (G) do not match E (`cluster_bootstrap_sensitivity.csv`).
- SI S6 prose (G) does not match E / the SI table in the same section.

These are P1. They do not change the independent Table 2 replay.
