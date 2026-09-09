# Figure asset notes (repository, not manuscript captions)

JCIM panel letters match `figures/jcim_article/`. Every plotted number is read from the frozen CSVs named in `docs/FIGURE_PANEL_LOCK_V3.md`.

## Provenance

- **Pinned numerical-data snapshot:** `abb61a20a04eb6a085ad526876624eadb518c4cc`
- Input files and SHA-256 values: `plotted_values.json` (`data_snapshot_commit`, `artwork_git_head`)
- Figure file hashes: `FIGURE_SHA256.json`
- Manuscript captions: `MANUSCRIPT_FIGURE_CAPTIONS.md`

Regenerate:

```text
python3 figures/jcim_article/scripts/update_figures_pr32.py
python3 figures/jcim_article/scripts/audit_figures_pr32.py
```

The snapshot commit pins tables and assembled manuscripts. Artwork may be newer than that commit; `artwork_git_head` records the git HEAD when the generator ran.

## Figure 1 split assets

- `../editable/Figure1_AB_editable_v2.pptx`: editable conceptual A–B schematic.
- `Fig1_C_chEMBL_supply.png/.tif/.pdf`: Python-rendered census panel.
- `Fig1_four_state_and_supply.*`: combined reference composite.

## Historical SI figure IDs

S2, S3, S9, and S10 remain SI **table** records. No artwork for those original-set panels is shipped. Table numbers S1–S13 are unchanged.

S8 and S11 repeat Figure 6C/D and 6B for SI readers who want a standalone page.
