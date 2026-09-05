# Five-pair dump-gated stack

Destination: same JCIM article, 8-row main table after withdrawing
PIK3CA/PIK3CB (`PROJECT_IDENTITY_LOCK_V1.md`). These five pairs were added
after the ChEMBL 37 census. This run does **not** restock Table 2.

sqlite: `/tmp/chembl/chembl_37/chembl_37_sqlite/chembl_37.db`

Tarball SHA-256:
`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`
(matches `jcim_chembl_universe_v0/README.md`).

Document / activity join uses the same STANDARD_OK + non-null pChEMBL
endpoints as panel harvest. That is the dump analogue of
`TIME_SPLIT_PROTOCOL_FREEZE.md`. It is **not** the later K=4 API
high-confidence audit.

## Checks that had to pass before IDs were frozen

| check | result |
|---|---|
| dump max vs panel pChEMBL (tol 0.015) | **0 / 545** mismatches |
| leftover D/A/B vs frozen `track_b_panel_summary_v1.csv` | **MATCH** on all five pairs |
| holdout overlap with *same-pair* main panel | **0** (by construction) |
| undated scored panel ligands | **0** |

If leftover had not matched the frozen summary, holdout IDs would not have
been drawn.

## max vs median (θ = 6.0, same frozen Vina)

One ligand flips class: `PPARA/PPARD` `PAPD_002` / CHEMBL121
(max A = 8.40, median A = 5.46, 18 PPARA activities → A_only → neither).
Directional AUROC is unchanged on four pairs; PPARA/PPARD D-vs-A moves
0.6465 → 0.6361. Do not replace max-pChEMBL labels with median.

## Literature-year split

Ligand year = **earliest** `docs.year` among retained STANDARD_OK records.
Train: first year < cutoff. Test: first year ≥ cutoff. Cutoffs 2015 / 2018 /
2020 frozen before AUROC. Test AUROC only if dual, A-only, and B-only each
n ≥ 10. Train is counts only. Do not shop cutoffs.

2018 test (primary):

| pair | test D/A/B/N | AUROC reportable | D/A | D/B | min [95% CI] |
|---|---:|---:|---:|---:|---|
| F2/F10 | 4/4/0/1 | no | — | — | counts only |
| JAK1/TYK2 | 16/21/27/5 | yes | 0.6220 | 0.3681 | 0.3681 [0.2000, 0.5452] |
| JAK1/JAK2 | 21/21/24/6 | yes | 0.6531 | 0.7937 | 0.6531 [0.4755, 0.8100] |
| PPARG/PPARA | 0/0/2/2 | no | — | — | counts only |
| PPARA/PPARD | 6/1/2/4 | no | — | — | counts only |

The protocol’s *necessary* pair-count gate (≥2 pairs with n ≥ 10 at 2018)
is met. That is **not** sufficient to package this as manuscript external
validation:

- the two passing pairs share JAK1 and are not independent
- JAK1/TYK2 test `summary_min` is anti-directional (0.368; CI [0.20, 0.55]
  includes 0.5)
- JAK1/JAK2 test min CI also includes 0.5
- this is still the same scored ChEMBL panel, not a new harvest
- identity lock: do not restock Table 2 or retitle the article

2020 fails the n ≥ 10 gate on every pair (JAK1/TYK2 test dual = 9). Do not
move the cutoff to recover 2015.

## Document-cluster bootstrap and document-blocked CV

Document groups = ligands sharing a STANDARD_OK document, union-find
connected. Point AUROCs match the zero-dock stack. Cluster CIs are wide
(JAK1/TYK2 D-vs-A [0.16, 0.80]). Document-blocked GroupKFold: ligand-only
ECFP4 still beats docking; several docking CV values are far below 0.5
(JAK1/TYK2 0.10 / 0.15). Those are diagnostics, not a new primary table.

## Leftover holdout IDs (seed 20260731, Murcko cap 3)

| pair | leftover D/A/B | drawn | note |
|---|---:|---:|---|
| F2/F10 | 312/76/245 | 20/20/20 | |
| JAK1/TYK2 | 1874/59/80 | 20/20/20 | |
| JAK1/JAK2 | 5953/76/21 | **20/20/18** | leftover B = 21 is eligible; Murcko cap blocked 2 |
| PPARG/PPARA | 408/50/59 | 20/20/20 | |
| PPARA/PPARD | 187/50/68 | 20/20/20 | |

JAK1/JAK2 leftover B-only remains eligible (21 ≥ 20). The *drawn* panel is
20/20/18 because one leftover scaffold already hit cap 3. Do **not** relax
the cap to invent two more ligands. Do **not** mark the pair EGFR-style
ineligible. If this holdout is docked later, report 18 B-only, not 20.

Cross-pair chemistry (expected; leftover is per-pair):

- 13 holdout IDs appear in a *different* pair’s main panel (JAK or PPAR)
- 3 molecule IDs appear in two holdout panels
  (`CHEMBL5831731`, `CHEMBL5901583` JAK B-only both pairs;
  `CHEMBL148774` PPARG/PPARA B-only and PPARA/PPARD A-only)

See `holdout_cross_pair_overlap_v1.csv`. IDs are frozen. Do not re-draw
after seeing scores.

## What this is not

- Not Table 2 and not a title change
- Not BindingDB docking
- Not a reason to re-dock PIK3CA/PIK3CB or run ordinary Vina on CTSK/CTSS
