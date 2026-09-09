# PR #32 figure update

Source: PR #32, commit `abb61a20a04eb6a085ad526876624eadb518c4cc`.
The source manuscripts and every table read by the generator are retained in
`input_snapshot/`, with SHA-256 values in `plotted_values.json`.

Run from any directory:

```text
python path/to/figures/jcim_article/scripts/update_figures_pr32.py
```

The generator uses the established local v3 style (Arial when available,
Okabe–Ito colors, 7-inch maximum width, 300-dpi RGB PNG/TIFF, embedded-font PDF).
The original checkout's older data directory is not used for these outputs.

## Content and chart choices

| Figure | Current content | Presentation reason |
|---|---|---|
| 1 | Four states, two tasks, whole-database census | State matrix and count boxes preserve exact counts over several orders of magnitude without a misleading linear funnel. The conceptual A–B schematic is supplied separately as an editable PowerPoint; the C census panel is supplied as a data-rendered figure. |
| 2 | Eight directional AUROCs, Table 3 comparison, both fixed-score negative-class contrasts | Paired points and confidence intervals show magnitude and uncertainty; no zero-based bars for AUROC. |
| 3 | Sixteen ECFP4/Vina comparisons, incremental AUROC, AChE/BChE TPSA | Point comparisons and individual-compound distributions preserve the paired design and sample spread. |
| 4 | Three independent GNINA systems, PIK3CA/mTOR receptor changes, five Vina seeds | Separate panels distinguish engine point estimates, bootstrap intervals, and seed ranges. |
| 5 | Matched-minus-mismatched summary differences and seven internal holdouts | Forest plots distinguish paired differences from absolute AUROC. |
| 6 | Label sensitivity, panel size, exhaustiveness, eight-pair BindingDB count gate | Categorical heatmaps show missing class coverage; color saturation marks the gate, not a performance measure. |
| S4 | Vina intervals and best descriptor points | Descriptive reference, without inventing descriptor confidence intervals. |
| S5 | Seven holdout/main comparisons | Matches Table S7. |
| S6 | Detectability simulation in the three available systems | Full true-AUROC grid 0.50–0.75; other five pairs were not simulated in this source. |
| S7 | EGFR/HER2 Top-10 and AND-filter composition | Separate count charts preserve different denominators (110 and 98). |
| S8 | Filtered compound and source counts for all eight pairs | Separate threshold-relative matrices show both external admission conditions. |
| S11 | EGFR/HER2 and JAK1/TYK2 fixed-score difference intervals | Ligand, scaffold and document intervals shown together, exposing JAK1/TYK2 document dependence. |

S11 is an additional **figure** identifier and does not renumber SI tables.
Tables S1–S13 remain unchanged. Figure numbers 1–6 and the already cited S4/S6
are preserved. No manuscript text is overwritten by this figure-only update.

`archive_before_pr32/` contains the previous local outputs. Original-set
S1–S3/S9/S10 are historical records and are excluded from the current top-level
figure set. They must not be mixed into the eight-pair submission figures.

Figure 1 split deliverables:

- `../editable/Figure1_AB_editable_v2.pptx`: editable conceptual A–B schematic.
- `Fig1_C_chEMBL_supply.png/.tif/.pdf`: Python-rendered census panel from the pinned CSV.
- `Fig1_four_state_and_supply.*`: retained combined reference composite; it is not required when the submission system accepts the two split assets separately.

## Interpretation details

- EGFR/HER2 GNINA Dual-versus-neither has 11 scored negative ligands, versus
  12 for Vina. PIK3CA/mTOR GNINA Dual-versus-A-only has 13 negatives versus
  14 for Vina. These are available-score comparisons, not common-ligand
  performance tests between algorithms.
- The original three-pair GNINA table supplies intervals for individual arms,
  but not a bootstrapped minimum-of-two interval. Figure 4A shows point
  estimates and does not substitute weak-arm intervals for summary intervals.
- Simulation Figure S6 reports two-sided CI exclusion of 0.5, including at
  true AUROC 0.50. It is not a probability of demonstrating superiority.
- BindingDB counts describe the supplied snapshot and eligibility rules,
  not successful independent docking validation.
