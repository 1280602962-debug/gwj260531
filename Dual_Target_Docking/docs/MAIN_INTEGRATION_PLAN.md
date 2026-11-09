# Integration onto main (not executed)

Freeze branch: `cursor/scientific-freeze-c7cc`  
Baseline: PR #38 `cursor/methods-sentence-audit-c7cc` @ `17435413410290960da3437de640509705f00b51`  
This freeze PR: #39 (draft).  
**Do not merge, close, or delete anything without an explicit confirmation.**

## What would go to main

Only `Dual_Target_Docking/`. `JNK1_Selectivity_Project/` is untouched on this branch.

Relative to `origin/main` (`f2b50243`), Dual_Target_Docking is the JCIM evaluation stack that never landed (PRs #23–#38). The freeze is a scientific lock on top of that stack, not a new docking campaign.

Suggested merge shape, when confirmed:

1. Keep main’s JNK1 history as-is.
2. Bring in Dual_Target_Docking from this freeze (or from PR #38 + this freeze) as one directory add/update.
3. Do not merge PR #35 as a parent. RMSD unification is already represented by `all14_cognate_rmsd_calcrrms_v1.csv`; reconstructed EGFR 0.760 Å stays out.
4. Leave historical Dual_Target PRs open until the user closes them.

## Other-topic check

| Path | This freeze |
|---|---|
| `JNK1_Selectivity_Project/` | no edits |
| `.github/workflows/revision-validate.yml` | unchanged; still runs Dual_Target `check_current_chain.py` |
| repo `README.md` | unchanged |

`git diff origin/main -- JNK1_Selectivity_Project` from this branch should be the existing main-vs-stack difference, not freeze work.

## Old Dual_Target PRs

| PR | Recommendation |
|---|---|
| #38 | Baseline of this freeze. Keep open. After freeze is accepted, #38+ freeze is the science tree. |
| #37 language polish | Superseded for science; optional language-only follow-up after freeze numbers. |
| #36 manuscript consolidation | Same. |
| #35 RMSD | Do not merge. Item decisions in `docs/PR35_ITEM_DECISIONS.md`. |
| #32–#34, #23–#30 | Historical stack. Do not merge piecemeal. |
| #26 census | Open; not this freeze. |

## Reviewable diff

Freeze-only: `git diff 17435413..cursor/scientific-freeze-c7cc -- Dual_Target_Docking`  
Main-facing Dual_Target tree: `git diff origin/main...cursor/scientific-freeze-c7cc -- Dual_Target_Docking`

No push to main. No `--delete` of remote branches.
