# Project identity lock (2026-09-05)

This file is the destination identity of the **same JCIM article**.
It supersedes “frozen K=4 + optional Track B chapter”, the intermediate
“Three-Pair” manuscript title edit of 2026-09-05, and any “prespecified
eight/nine-pair freeze from 2026-07-23” wording.

Read this before answering later questions about what to analyse,
what belongs in Table 2, or what the title may eventually say.

## Locked destination set

Historically docked = **9 pairs** (original 4 + 5 post-census additions).
Main Table 2 / Table 3 destination = **8 rows**.
Not 8 because CTSK/CTSS was pulled in. Not 9 rows in the main table
because PIK3CA/PIK3CB is a documented receptor-identity failure.

| # | Pair | Role in destination main table | How it entered |
|---|------|--------------------------------|----------------|
| 1 | PIK3CA/mTOR | ordinary thick | 2026-07-23 J0 freeze |
| 2 | AChE/BChE | ordinary thick | 2026-07-23 J0 freeze |
| 3 | EGFR/HER2 | supply-limited, **still a main-table row** | 2026-07-23 J0 freeze (same role as the original paper) |
| 4 | F2/F10 | ordinary noncovalent | ChEMBL 37 census, then docked |
| 5 | JAK1/TYK2 | ordinary noncovalent | ChEMBL 37 census, then docked |
| 6 | JAK1/JAK2 | ordinary noncovalent | ChEMBL 37 census, then docked |
| 7 | PPARG/PPARA | ordinary noncovalent | ChEMBL 37 census, then docked |
| 8 | PPARA/PPARD | ordinary noncovalent | ChEMBL 37 census, then docked |

| Pair | Status | Must not be done |
|------|--------|------------------|
| PIK3CA/PIK3CB | SI-only receptor-identity failure (2WXF = murine PIK3CD / O35904, not human PIK3CB / P42338). Existing scores stay in the repo. | Re-dock; 2Y3A production swap; put back as a peer Table 2 row; call it a human p110β isoform control |
| CTSK/CTSS | Census G1–G5 on labels/holos; reversible-covalent crystals; out of unified AutoDock Vina | Ordinary Vina to “make 8 census systems” |
| CREBBP/BRD4, GPCRs, SLC6, F2/PRSS1, JAK3/TYK2 | Failed a declared gate | Silent restock after seeing AUROCs |

Arithmetic that must stay visible:

- Original freeze = 4.
- New ordinary noncovalent = 5.
- 4 + 5 = **9 historically docked**.
- Withdraw PIK3CA/PIK3CB from the main table → **8 main-table rows**.
- Census dockable-system slot 8 was CTSK/CTSS. That slot is **not** a substitute for the withdrawn PIK3CB row.

EGFR stays **in** the main table as the supply-limited case the original paper already treated that way. It is not moved to a footnote to manufacture a cleaner “ordinary-only” table.

## Locked timeline (honesty constraint)

Experiments on the five new pairs may copy the **already written** K=4
analysis stack. The timeline must **not** be rewritten as:

> “We froze these eight (or nine) pairs on 2026-07-23 before scores were seen.”

That would replace a collection-incompleteness gap with a HARKing /
selection-bias gap.

Required Methods 2.0 sentence when the manuscript is eventually updated
(only after the analysis stack is actually complete; do not change the
title or Table 2 identity again before that):

> The original primary set was the 2026-07-23 four-pair freeze. After a
> post-hoc ChEMBL 37 universe census, five additional ordinary
> noncovalent pairs that passed the pre-declared G1–G5 gates were added
> to the **same** analysis stack. This is collection completion, not a
> pre-registered eight-pair freeze.

Do not call the five pairs “Track B” in the destination article.
“Track B” remains a repository folder name only
(`local_track_b_v0/`).

## What “same analysis stack” means

Copy **rules**, not a fantasy all-pairs menu.

All 8 destination main-table pairs must have every item that Methods
2.4–2.6 already reported for **every then-primary pair**.

Do **not** spread items the original paper only ran on a subset:

| Item | Original rule | Destination rule |
|------|---------------|------------------|
| Independent GNINA search | EGFR/HER2 and PIK3CA/mTOR only (formulation-gap pairs) | Add **JAK1/TYK2 only**. Do not run on JAK1/JAK2, F2/F10, or either PPAR pair |
| unused-pool holdout | 20/20/20 after excluding main (+ PM110 on PM); EGFR ineligible | Same numeric gate. JAK1/JAK2 leftover strict B-only small-mol = 21; eligible but thin (margin = 1). EGFR remains ineligible. PIK3CA/PIK3CB is not carried forward |
| Crystal swap | PIK3CA 4JPS/5DXT and mTOR 4JSX only | Do not add swaps to the five new pairs unless a new written protocol names the crystals first |
| PM110 / E=16 | PIK3CA/mTOR only | Do not apply |

## Current manuscript is not the destination

The 2026-09-05 edit that retitled the article “Three-Pair” and removed
the PIK3CA/PIK3CB row from Tables 1–3 is an **intermediate** withdrawal
of a failed receptor, not the destination identity. Do not treat
“Three-Pair Formulation Audit” as the final title. Do not further
rewrite Abstract / Table 2 / title until the five new pairs have the
all-pairs stack (zero-dock + five-seed / RTM / GNINA-CNN). Then expand
the main table to the 8 rows above and retitle as a multi-pair
four-state formulation audit (not “Bench”, not a pre-frozen eight-pair
study).

## Forbidden recurrences

- Answering a later question as if the live project were still only the
  original four pairs.
- Opening a separate Track B chapter instead of the same article.
- Averaging eight or nine pairs into one master AUROC.
- Pulling CTSK/CTSS or 2Y3A into ordinary Vina to hit a round number.
- Changing LigPrep or the dock seed to 42.
- Hard-docking the failed BindingDB external slice.
