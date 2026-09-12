# JCIM submission pack

This directory is the publication-facing slice of DualFourClass-Bench.
It is not a second copy of the docking pose workspaces.

## Use these files for submission

| Path | Role |
|---|---|
| `manuscript/MANUSCRIPT_JCIM_EN.md` | English manuscript (assemble from section drafts; do not hand-edit the assembled file as the source of truth) |
| `manuscript/MANUSCRIPT_JCIM_ZH.md` | Chinese working manuscript |
| `manuscript/SUPPORTING_INFORMATION_JCIM_EN_V1.md` | English SI Tables S1–S14 |
| `manuscript/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md` | Chinese SI |
| `manuscript/STATISTICAL_LOCK_V1.md` | Table 2 / Table 3 estimand lock |
| `manuscript/FIGURE_PANEL_LOCK_V3.md` | Figure-to-CSV map |
| `manuscript/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md` | Five-round numeric audit |
| `tables/` | Frozen CSVs cited by Tables 1–3 and S2–S11 |
| `figures/` | Regenerated main and SI figures from `figures/jcim_article/scripts/update_figures_pr32.py` |
| `scripts/` | Assemble, validate, bootstrap-lock reader, figure v3, audit, pack |

## Do not submit as primary evidence

- `plot_jcim_article_figures_v1.py` / `v2.py` (withdrawn-pair leftovers)
- `plot_jcim_si_composites_v1.py` S1–S3 / S9 / S10 (historical original-set artwork; removed from this repository; may still tick PIK3CA/PIK3CB)
- `data/pik3ca_pik3cb_panel_v0/` (withdrawn pair archive)
- `external_slice_summary_202608_contract_v1.csv` (legacy three-pair BindingDB snapshot)
- Pose workspaces and multi-GB score dumps (indexed in the repository, not copied here)
- A minted Zenodo DOI (not created)

## Rebuild

```bash
python3 docs/assemble_manuscript_en.py
python3 docs/assemble_manuscript_zh.py
python3 scripts/primary/bootstrap_primary.py
python3 data/jcim_novelty_v0/scripts/validate_revision_v1.py
python3 data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py --check
python3 scripts/audit/audit_submission_five_rounds_v1.py
python3 figures/jcim_article/scripts/update_figures_pr32.py --source-root .
python3 figures/jcim_article/scripts/audit_figures_pr32.py
python3 scripts/audit/pack_submission_v1.py
```
