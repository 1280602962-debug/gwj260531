# Submission integration policy

This is one GitHub repository (`1280602962-debug/gwj260531`) with three DualFourClass lines:

| Line | Ref | Role |
|---|---|---|
| Freeze snapshot | `cursor/pik3ca-mtor-structure-freeze-0b1a@236fd60c` | Historical. Compressed manuscript, literature positioning, lock inventory. **Do not edit.** |
| ChatGPT snapshot | `chatgpt/jcim-manuscript-review-20260828@d6cc6c12` | Historical. Five-seed scores plus an uncorrected Dual-versus-neither estimand in the article. **Do not edit. Do not merge.** |
| Submission integration | `cursor/jcim-final-integration-0b1a` | **Only branch that continues.** Freeze data + freeze prose + corrected multi-seed v2. |

A full `git merge` of the ChatGPT branch is forbidden: it would restore the long audit-style Chinese Results, session language, MCL1-in-main-text placement, and the wrong Dual-versus-neither numbers in the abstract.

## What was absorbed from the ChatGPT snapshot

- Frozen five-seed Vina score table (`multiseed_scores_long_v1.csv`) — already present on freeze; reused, not recopied from a merge.
- The idea of a formal SI Table S54 for seed sensitivity.

## What was not absorbed

- ChatGPT abstract five-seed medians
- ChatGPT Results/Discussion gap conclusions written from `mean(AUC_A, AUC_B)`
- ChatGPT validator pins for the v1 S54 numbers
- ChatGPT checksums
- Any overwrite of freeze TITLE/ABSTRACT/INTRODUCTION/METHODS/RESULTS/DISCUSSION/CONCLUSIONS

## Working tree vs release

The Git working tree keeps research history. The journal/Zenodo snapshot is built by:

```bash
python3 scripts/package_jcim_release_v1.py
```

which copies only `KEEP_FOR_RELEASE=1` rows from `data/manuscript_lock/ARTICLE_ASSET_INDEX_v1.csv` into `release/JCIM_submission_v1/` (gitignored).
