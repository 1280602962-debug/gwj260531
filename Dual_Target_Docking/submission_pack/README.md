# JCIM submission pack

This directory is the publication-facing slice of DualFourClass-Bench.
It is not a second copy of the docking pose workspaces.

## Use these files for submission

| Path | Role |
|---|---|
| `manuscript/MANUSCRIPT_JCIM_EN.md` | English manuscript (assemble from section drafts; do not hand-edit the assembled file as the source of truth) |
| `manuscript/MANUSCRIPT_JCIM_ZH.md` | Chinese working manuscript |
| `manuscript/SUPPORTING_INFORMATION_JCIM_EN_V1.md` | English SI Tables S1–S9 |
| `manuscript/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` | Chinese SI |
| `manuscript/STATISTICAL_LOCK_V1.md` | Table 2 / Table 3 estimand lock |
| `manuscript/FIGURE_PANEL_LOCK_V3.md` | Figure-to-CSV map |
| `manuscript/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md` | Five-round numeric audit |
| `tables/` | Frozen CSVs cited by Tables 1–3 and S1–S9 |
| `figures/` | Regenerated main and SI figures from `figures/jcim_article/scripts/update_figures_pr32.py` |
| `scripts/` | Assemble, validate, freeze, figure audit, pack |

## Do not submit as primary evidence

- `plot_jcim_article_figures_v1.py` / `v2.py` (withdrawn-pair leftovers)
- `plot_jcim_si_composites_v1.py` S1–S3 / S9 / S10 (historical original-set artwork; removed from this repository; may still tick PIK3CA/PIK3CB)
- `data/pik3ca_pik3cb_panel_v0/` (withdrawn pair archive)
- `external_slice_summary_202608_contract_v1.csv` (legacy three-pair BindingDB snapshot)
- Pose workspaces and multi-GB score dumps (indexed in the repository, not copied here)
- A minted Zenodo DOI (not created)

## Rebuild / submission freeze

Do not replot figures or recompute scientific tables. SI markdown is hand-edited; the freeze runner only reassembles the main manuscripts.

```bash
python3 scripts/audit/freeze_submission_v1.py
```
