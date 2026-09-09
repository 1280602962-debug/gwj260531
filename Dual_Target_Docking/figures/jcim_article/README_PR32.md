# PR #32 figure update

Pinned numerical-data snapshot: `abb61a20a04eb6a085ad526876624eadb518c4cc`.
Artwork is generated from that snapshot; `plotted_values.json` records both
`data_snapshot_commit` and `artwork_git_head` (git HEAD when the generator ran).
Input SHA-256 values are in the same JSON file.

Manuscript captions: `MANUSCRIPT_FIGURE_CAPTIONS.md`.
Repository notes: `FIGURE_ASSET_NOTES.md`.

Run from any directory:

```text
python path/to/figures/jcim_article/scripts/update_figures_pr32.py --source-root Dual_Target_Docking
python path/to/figures/jcim_article/scripts/audit_figures_pr32.py
```

The generator uses the established local v3 style (Arial when available,
Okabe–Ito colors, 7-inch maximum width, 300-dpi RGB PNG/TIFF, embedded-font PDF).

Target pairs follow a protein-system display order on Figures 2–6 and Tables 1–3.
Historical CSV splits remain routing keys only and are not drawn.

## Content and chart choices

| Figure | Current content | Presentation reason |
|---|---|---|
| 1 | Four states, two directional tasks, whole-database census | State matrix and count boxes preserve exact counts over several orders of magnitude without a misleading linear funnel. |
| 2 | Fixed-score ΔAUROC (hero), directional dumbbells, Dual–neither | Forest and paired points show that the same scores change interpretation with the experimental-state comparison. |
| 3 | Sixteen ECFP4/Vina comparisons, incremental AUROC (±0.03), illustrative AChE/BChE TPSA | Point comparisons and individual-compound distributions preserve the paired design. |
| 4 | Three independent GNINA systems, PIK3CA/mTOR receptor changes, five Vina seeds | 2-up / 1-full-width layout; point plots emphasize the task gap, not engine ranking. |
| 5 | Main and holdout matched−mismatched Δ on the same rows; holdout vs main summary | Color encodes ligand set (main vs holdout), not interval exclusion of zero. |
| 6 | Label-rule heatmap, cluster resampling of fixed-score Δ, BindingDB compound and source gates | Evidence boundary is label robustness, source dependence, and external supply. |
| S1 | PIK3CA/mTOR panel size and exhaustiveness | Protocol sensitivity demoted from the main evidence-boundary figure. |
| S4 | Vina intervals and best descriptor points | Descriptive reference. |
| S5 | Seven holdout/main comparisons | Matches Table S7. |
| S6 | Detectability simulation in the three available systems | Full true-AUROC grid 0.50–0.75. |
| S7 | EGFR/HER2 Top-10 and AND-filter composition | Exploratory; remains SI. |
| S8 | Filtered compound and source counts | Same matrices as Figure 6C/D. |
| S11 | EGFR/HER2 and JAK1/TYK2 cluster intervals | SI copy of Figure 6B. |
| S12 | Cognate RMSD for 14 primary receptors (top-1 and best-of-9) | Protocol QC, not a main-text structural claim. |

S11 is an additional **figure** identifier and does not renumber SI tables.
Tables S1–S13 remain unchanged. Historical original-set figure files are not kept
in this directory.

Figure 1 split deliverables:

- `../editable/Figure1_AB_editable_v2.pptx`: editable conceptual A–B schematic.
- `Fig1_C_chEMBL_supply.png/.tif/.pdf`: Python-rendered census panel from the pinned CSV.
- `Fig1_four_state_and_supply.*`: retained combined reference composite.

## Interpretation details

- EGFR/HER2 GNINA Dual-versus-neither has 11 scored negative ligands, versus
  12 for Vina. PIK3CA/mTOR GNINA Dual-versus-A-only has 13 negatives versus
  14 for Vina. These are available-score comparisons, not common-ligand
  performance tests between algorithms.
- Figure 4A shows point estimates and does not substitute weak-arm intervals
  for summary intervals.
- Simulation Figure S6 reports two-sided CI exclusion of 0.5, including at
  true AUROC 0.50. It is not a probability of demonstrating superiority.
- BindingDB counts describe the supplied snapshot and eligibility rules,
  not successful independent docking validation.
