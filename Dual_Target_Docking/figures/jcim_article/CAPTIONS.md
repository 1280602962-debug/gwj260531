# Figure captions (manuscript; not printed on the image)

JCIM: captions are self-contained; panel letters match the files in `figures/jcim_article/`.
All numerical values are those plotted from the frozen CSVs (unrounded). Table 2 in the text may round to three decimals.

Regenerate (from `Dual_Target_Docking/`):
`python3 data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py`
The script re-reads the CSVs, writes PDF/PNG/TIF, and fails if any plotted value disagrees with its source (`plotted_values.json`).

## Figure 1. DualFourClass task and pocket-matched readout.

(A) Dual-target evaluation requires four ligand classes measured on both targets: dual-actives and experimental single-end hard negatives (A_only, B_only), plus neither. (B) The primary metric scores dual versus A_only in pocket B and dual versus B_only in pocket A; summary_min is the smaller of the two AUROCs so that a weak arm cannot be hidden by pooling.

## Figure 2. Public-data supply of strict hard negatives.

(A) Minimum of the two strict hard-negative counts (A_only, B_only) for every target pair in the J0 ChEMBL audit (`j0_strict_label_supply.csv`). Dashed line, thick-panel gate (≥50); dotted line, thin-panel gate (≥20). Highlighted: the three thick pairs used as K=4 main panels, EGFR/HER2 (7 B_only; supply-limited case), and HDAC1/HDAC6 (metal enzyme; excluded). (B) Count-level comparison of the same four pairs in ChEMBL pChEMBL versus BindingDB equal-relation measurements (Table S12). No docking.

## Figure 3. Pocket-matched summary_min on the frozen K=4 set.

Vina (primary), RTMScore, GNINA CNN best-of-9, and the strongest trivial descriptor (heavy-atom count, MW, cLogP, or TPSA) with 95% ligand-bootstrap CIs. Vina CIs are the θ = 6.0 values in `unified_threshold_sensitivity_v2.csv` (Table 2). Best descriptor (right column, from `forest_summary_min_ci_v1.csv`): EGFR/HER2 cLogP; AChE/BChE TPSA; PIK3CA/PIK3CB and PIK3CA/mTOR heavy-atom count. Vertical dashed line, chance (0.5). GNINA is a single CNN channel, not a three-engine competition.

## Figure 4. Weak-arm asymmetry and physicochemical confounding.

(A) Directional Vina AUROCs at θ = 6.0: dual versus A_only (pocket B) and dual versus B_only (pocket A). (B) Vina pocket-matched summary_min versus the strongest trivial descriptor, with 95% CIs. (C) TPSA on the AChE/BChE panel by class (individual ligands from `assembled_AChE_BChE.csv`; horizontal line, median). Dual ligands are more polar than either hard-negative class, matching the TPSA baseline that exceeds Vina on this pair.

## Figure 5. Ligand-set holdout versus receptor swap.

(A) Pocket-matched summary_min on the main panel versus the unused-pool holdout (20/20/20; seed 20260731) for the three pairs with unused-pool supply. EGFR/HER2 has no holdout. (B) PM48 crystal swap: replacing PIK3CA 4L23 with 4JPS or 5DXT (mTOR held at 4JT6) or replacing mTOR 4JT6 with 4JSX (PIK3CA held at 4L23). Error bars are 95% ligand-bootstrap CIs.

## Figure S1. Wrong-pocket control on the main panels and the holdout.

(A) Main K=4 panels: pocket-matched Vina summary_min versus the wrong-pocket control (`pocket_matched_directional_v1.csv`). Matched exceeds wrong-pocket on all four pairs. (B) Unused-pool holdout: the inequality reverses (wrong-pocket ≥ matched) on all three pairs (`holdout_pocket_matched_v1.csv`).

## TOC graphic (For Table of Contents Only).

DualFourClass-Bench evaluates whether docking ranks dual-target ligands above experimental single-end hard negatives in both pockets. The graphic does not report numerical AUROCs and is not a reuse of Figure 1.
