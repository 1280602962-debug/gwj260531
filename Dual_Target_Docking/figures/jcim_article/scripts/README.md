# JCIM figure scripts

This directory is the organized, reproducible entry point for the figures in the
parent directory. All scripts read the frozen CSV/JSON inputs under `data/` and
write the generated files to `figures/jcim_article/` through the shared style
module.

## Current generators

- `plot_jcim_article_figures_v3.py`: Figures 1–6, S4, S5, S7, S8, and the TOC graphic.
- `plot_jcim_si_composites_v1.py`: Figures S1–S3, S9, and S10.
- `plot_detectable_effect_and_workflow_v1.py`: Figure S6 (the optional workflow
  schematic is retained in the source script but is not part of the current
  figure set unless explicitly regenerated).
- `jcim_figure_style.py`: shared fonts, colors, output paths, and `save_all`.

The canonical SI filenames are `FigS1_...` through `FigS10_...` where present.
The three-file PDF/PNG/TIF triplets are intentional output formats, not duplicates.
Legacy source copies remain under `data/*/scripts/` so existing commands and
provenance links continue to work; use the copies here for figure-folder
organization.
