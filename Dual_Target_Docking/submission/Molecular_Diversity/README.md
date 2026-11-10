# Molecular Diversity figures and tables

This is the current presentation package prepared for Molecular Diversity on 10 November 2026 and refined from presentation commit `41d48a0a6cf7d112b805b299bf2618d74654ab63`. Scientific estimates, analysis membership and confidence intervals come from the accepted branch at `ce81b011bab449c6cf5e2b9664626ed808be5857`. The main text has eight figures and two tables, meeting the journal's ten-item display limit. Previous Figure 8 is now Figure S10; previous Figure 9 is now Figure 8. Table numbering is unchanged.

## Files to use

| Content | Location |
| --- | --- |
| Main figures 1–8 | `Main_Figures/Fig1`–`Fig8` |
| Supplementary figures S1–S10 and Scheme S1 | `Supplementary_Figures/` |
| All figures in one PDF | [Figures_Review.pdf](Figures_Review.pdf) |
| Paginated table proof | [Tables_Review.pdf](Tables_Review.pdf) |
| Editable main tables | [Tables/Main_Tables.xlsx](Tables/Main_Tables.xlsx) |
| Editable supplementary tables | [Tables/Supplementary_Tables.xlsx](Tables/Supplementary_Tables.xlsx) |
| Native HTML tables for copying into a manuscript editor | [Tables/Editable_Tables.html](Tables/Editable_Tables.html) |
| English figure captions | [Figure_Captions_EN.md](Figure_Captions_EN.md) |
| Chinese figure captions | [Figure_Captions_CN.md](Figure_Captions_CN.md) |
| Figure-specific data and full-precision numbered-table CSVs | `Source_Data/` |
| Numeric, source-cell and export verification | `Checks/` |
| Figure-by-figure chart-selection rationale | [Visualization_Review_CN.md](Visualization_Review_CN.md) |
| Previous and revised encodings for four figures | [Visualization_Comparison.pdf](Visualization_Comparison.pdf) |
| Numbering changes | [Figure_Renumbering.csv](Figure_Renumbering.csv) |

EPS is the preferred submission vector export; the package also includes embedded-font PDF, editable-text SVG and RGB 600-dpi PNG. Figures are 174 mm wide, with lettering of at least 8 pt. Insert the full captions in the manuscript; the individual artwork files contain panel labels and plot annotations. The tables remain editable. Their PDF is a layout proof rather than a replacement for native table cells.

## Presentation changes

Direction labels separate the target pair from the single-target negative class. The two tasks remain dual versus A-only using target B and dual versus B-only using target A. Method-specific estimates, paired method differences and joint-ranking populations retain their separate denominators. Minimum AUROC denotes the smaller directional AUROC; its bootstrap value is recalculated within each resample.

Figure 5 now separates the model matrix from signed, labelled score increments. Figure 7a shows changes from seed 42 with absolute baseline values in row labels; the full absolute results remain in Figure S7. Figure S2 uses unavailable-score counts for color and retains every valid/expected annotation. Figure S6 has separate class-specific distribution summaries. All four source-data exports are byte-identical to the preceding presentation version. These changes improve specific reading tasks; they do not establish a universally optimal chart type.

Table 2 retains all 32 matched comparisons in one page. Counts and recall are shown as S1 → comparison scheme; recall is displayed as a percentage. S4 versus S1 is marked as a supporting comparison. Supplementary tables combine estimates with their own 95% confidence intervals, separate receptor identity from box geometry, and keep each target-pair group on one continuation page. Display rounding does not replace the full-precision source CSVs.

The previous submission exports remain as provenance. Use this directory for the revised presentation. No new docking, model fitting, bootstrap analysis or selection of study members was performed.

## Journal and reference comparison

The official [Molecular Diversity submission guidelines](https://link.springer.com/journal/11030/submission-guidelines) were checked on 10 November 2026. The journal-specific double-column range is 105–177.8 mm; 174 mm was chosen to also fit the publisher's general large-journal artwork guidance. The package follows the requirements for font embedding, RGB color, lowercase panel letters, external captions, sequential numbering and editable tables. Separate supplemental figure and table numbering is retained; uploaded supplemental file naming can use the journal's Online Resource convention.

The user's three supplied reference PDFs were read, and their corresponding published figure images were inspected. Their panels and table hierarchy informed the presentation choices; no published image was reused.

| Reference | Relevant example | Application here |
| --- | --- | --- |
| Yu et al., GRU-based de novo design and in-silico prioritization of EZH2 inhibitors, DOI [10.1007/s11030-026-11481-x](https://doi.org/10.1007/s11030-026-11481-x) | Figs. 5 and 8; Tables 2–3 | Related comparisons grouped into panels with consistent axes; full numerical detail in tables |
| Jiang et al., Artificial intelligence-driven rational design and optimization of a potent terpenoid-derived PCSK9 inhibitor, DOI [10.1007/s11030-026-11564-9](https://doi.org/10.1007/s11030-026-11564-9) | Fig. 9; Tables 4–5 | Model comparisons grouped consistently; circular plots were not adopted because long labels and interval comparison are clearer on Cartesian axes |
| Bhambri et al., Structure-guided design and DFT-based optimization of natural product-derived dual inhibitors targeting CDK-1 and PARP-1, DOI [10.1007/s11030-025-11456-4](https://doi.org/10.1007/s11030-025-11456-4) | Fig. 5; Tables 3 and 5 | Matched target panels, explicit target identities and units; no unrelated molecular structures were added to the benchmark figures |

## Reproduction

Run the following from a checkout containing the accepted source data:

```text
python Dual_Target_Docking/scripts/rerun_v4_2/refine_molecular_diversity_visuals.py
python Dual_Target_Docking/scripts/rerun_v4_2/build_molecular_diversity_tables.py
python Dual_Target_Docking/scripts/rerun_v4_2/audit_molecular_diversity_artwork.py
python Dual_Target_Docking/scripts/rerun_v4_2/finalize_molecular_diversity_artwork.py
```

Dependencies: Python, Matplotlib, NumPy, PyMuPDF, Pillow and openpyxl. The existing figure builder supplies the accepted joins and figure-specific source definitions. New workbook authoring uses openpyxl because the artifact-tool runtime was unavailable in this session.
