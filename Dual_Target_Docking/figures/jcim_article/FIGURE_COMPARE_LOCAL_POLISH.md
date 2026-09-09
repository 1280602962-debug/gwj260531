# Local polish vs eight-row submission figures

Source of the polish set: Windows `Dual_Target_Docking/figures/jcim_article` (2026-09-08 redraw).
Full copy: `figures/jcim_article/local_polish_20260908/`.
Scripts only: `figures/jcim_article/scripts/`.

**Keep the existing PR-32 / `submission_pack` figures as the official eight-row set.**
Do not replace Fig 1–6, S4, S5, S7, S8, or the TOC with the local redraw.

## Why the PR-32 figures stay

The five-round audit (`docs/SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md`) locked three facts that the local redraw reverses:

| Lock | Official PR-32 figure | Local polish |
|---|---|---|
| J0 scrape pair count | 48; S7 bar “audited J0 pairs = 48” | 49; S7 “J0 pairs = 49” |
| Withdrawn PIK3CA/PIK3CB is not a J0 supply bar | Fig 1C right: HDAC, PIK3CA/mTOR, AChE/BChE, EGFR only | Fig 1C right adds a PIK3CA/PIK3CB bar (~56) |
| Fig 6D is the eight-pair external gate | 8 fail / 0 pass; “not docked” | only the four historical contract pairs, with a PIK3CB† row |

Local `plot_jcim_article_figures_v3.py` even asserts `n_pairs != 49` as an error. The branch script asserts `!= 48` and rejects a PIK3CB J0 bar. Those are not style differences.

S7 local “docked in J0 = 4” also reopens the historically docked / withdrawn-PIK3CB count that the audit moved to narrative. Official S7 uses “J0-era docked rows = 3” plus “primary now = 8”.

## Where the local redraw is better

Typography and layout are tighter (shorter canvases, forest/dot plots instead of some bars, 2×2 state matrix in Fig 1A). That is useful as a style reference, not as a replacement until the local script is re-locked to 48 / no PIK3CB bar / eight-pair Fig 6D.

S1–S3, S6, S9, S10 remain original-set SI. The local numbered files are the better-looking copies of those archive panels (S1 still marks PIK3CA/PIK3CB†). They were copied into `archive_original_set/` and must not be used as eight-row SI.

Numbered SI names (`FigS4_…`, `FigS5_…`, `FigS8_…`) are clearer than `FigS_pocket_matched_forest` etc. Rename only after a regenerate from the locked v3 script.

## Keep map

| Role | Keep |
|---|---|
| Main text Fig 1–6, TOC | PR-32 / `submission_pack` |
| Eight-row SI S4, S5, S7, S8 | PR-32 / `submission_pack` |
| Original-set SI S1–S3, S6, S9, S10 | local polish, stored under `archive_original_set/` |
| Figure scripts that regenerate the submission set | `data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py` (branch lock) |
| Local scripts | `figures/jcim_article/scripts/` for organization only; do not run them as the submission generator |

## If the local style should become official later

Re-run the local style module against the frozen CSVs **after** changing the Fig 1C J0 list and the `n_pairs == 48` check, then pass `scripts/audit/audit_submission_five_rounds_v1.py` and `plot_jcim_article_figures_v3.py` verification before replacing `submission_pack/figures/`.
