# Score master migration diff

One-time report. Authoritative file: `results/canonical/current_score_master.csv`.

## review_scored_membership_v1.csv vs master
- ligands in review: 797
- missing from master main: 0
- score mismatches (abs > 1e-6): 0


## Files that can be regenerated from the master
- eight-pair ranking / EF tables
- directional AUROC / summary_min / scheme-B CIs
- fixed-score ΔAUROC
- matched-minus-mismatched
- descriptor and ECFP4 models (SMILES come with the master)

## Files that still contain old EGFR scores (do not use as current analysis input)
- `data/jcim_novelty_v0/tables/review_scored_membership_v1.csv` (EGFR pre-fix scores)
- `data/_legacy_archive/egfr_her2_scores_pre_fix/ablation_ligand_scores.csv`
- `data/egfr_her2_panel120_v0/tables/scores_vina.csv` (stale sibling; production is ablation_ligand_scores.csv)

## Master complete-case counts (main, theta=6 class)
- EGFR/HER2: scored=110 activity-eligible=109 D=28 A=37 B=32 N=12
- JAK1/JAK2: scored=110 activity-eligible=110 D=32 A=32 B=32 N=14
- JAK1/TYK2: scored=109 activity-eligible=109 D=31 A=32 B=32 N=14
- PIK3CA/mTOR: scored=48 activity-eligible=48 D=18 A=14 B=12 N=4
- AChE/BChE: scored=96 activity-eligible=95 D=27 A=26 B=28 N=14
- F2/F10: scored=107 activity-eligible=107 D=31 A=32 B=32 N=12
- PPARG/PPARA: scored=109 activity-eligible=109 D=32 A=31 B=32 N=14
- PPARA/PPARD: scored=110 activity-eligible=110 D=32 A=32 B=32 N=14
