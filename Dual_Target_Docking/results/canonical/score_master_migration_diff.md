# Score master migration diff

One-time report. Authoritative file: `results/canonical/current_score_master.csv`.

## review_scored_membership_v1.csv vs master
- ligands in review: 799
- missing from master main: 0
- score mismatches (abs > 1e-6): 110
  - EGFR/HER2: 110
- **EGFR/HER2 review scores are pre-fix** if EGFR is in the mismatch list. Do not sort or bootstrap from that file.

- EGFR review rows matching pre-fix ablation vina_*_hb: 110

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
- EGFR/HER2: n=110 D=28 A=38 B=32 N=12
- JAK1/JAK2: n=110 D=32 A=32 B=32 N=14
- JAK1/TYK2: n=109 D=31 A=32 B=32 N=14
- PIK3CA/mTOR: n=48 D=18 A=14 B=12 N=4
- AChE/BChE: n=96 D=27 A=26 B=28 N=15
- F2/F10: n=107 D=31 A=32 B=32 N=12
- PPARG/PPARA: n=109 D=32 A=31 B=32 N=14
- PPARA/PPARD: n=110 D=32 A=32 B=32 N=14
