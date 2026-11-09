# JCIM figure scripts

Current generator:

```bash
python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root Dual_Target_Docking
```

That script imports `plot_jcim_article_figures_v3.py` for shared loaders and writes Figures 1–5, S1–S5, and the TOC graphic into `figures/jcim_article/`.

- `jcim_figure_style.py`: fonts, colors, `save_all`
- `inject_s2_boxes_from_json.py`: retained helper scanned by `scripts/qa/check_current_chain.py`

Plotted values: `../plotted_values_postfix.json`.
Numbering lock: `docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md`.
