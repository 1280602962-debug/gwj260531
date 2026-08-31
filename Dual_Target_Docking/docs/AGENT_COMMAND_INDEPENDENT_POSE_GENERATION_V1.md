# Local follow-up: independent pose-generation protocol

**Status (2026-08-26):** local GNINA 1.3.2 docking-search results are deposited in
`data/jcim_independent_dock_v0/`. EGFR/HER2 verdict: **gap remains**
(Dual versus neither 0.783 vs directional `summary_min` 0.220).
PIK3CA/mTOR: `summary_min` 0.633 vs Vina 0.692 (smaller, same sign).
This file remains the protocol record; do not rewrite it as an engine bake-off.

Zero-cloud-docking. This pack is for a **local** machine with GNINA 1.3.x
(or smina) and the frozen DualFourClass receptors/ligands. It is **not** a
multi-engine bake-off and must not be written up as “we benchmarked docking
engines.”

## Scientific question

Does the **principal formulation effect** persist when pose generation is
changed?

Primary frozen finding to protect:

- EGFR/HER2 Dual versus neither AUROC **0.756** vs directional `summary_min` **0.430**
  (same Vina scores; Table 3).

Secondary directional case:

- PIK3CA/mTOR PM48 `summary_min` **0.692** (Table 2).

GNINA CNN / RTMScore already in the paper **rescore frozen Vina poses**. They
are independent **scoring channels**, not independent **pose-generation
engines**. This pack fills that gap.

## Engine choice (one is enough)

Preferred (open, reproducible, already used in this project for rescoring):

```text
GNINA 1.3.2  docking mode
  --cnn_scoring rescore   # still a docking search; do NOT reuse Vina PDBQT poses
```

If GNINA docking mode is unavailable, use **smina** with the same frozen
receptor PDBQT, boxes, exhaustiveness, seed, and `num_modes = 9`.

Do **not** add Glide / GOLD / AutoDock4 unless a later revision needs a
second independent engine.

## Frozen inputs (do not rebuild panels)

| Pair | Ligands | Receptor A | Receptor B | Exhaustiveness | Box JSON |
|------|---------|------------|------------|----------------:|----------|
| EGFR/HER2 | `data/egfr_her2_panel120_v0/` panel + ligand PDBQT | 3POZ | 3RCD | 8 | existing `boxes/*.json` |
| PIK3CA/mTOR PM48 | `data/pik3ca_mtor_panel48_rdkit_v0/` | 4L23 | 4JT6 | 16 | existing boxes |

Reuse the **same** Meeko ligand PDBQT and receptor PDBQT as the Vina primary
run. Seed `20260727`. `n_modes = 9`. `energy_range = 3`.

Skip AChE/BChE and PIK3CA/PIK3CB unless EGFR/HER2 + PM48 already answer the
formulation question.

## Endpoints (fixed before scores are seen)

On each engine’s own mode-1 scores, with unified θ = 6.0 labels:

1. Dual versus neither (`mean` of the two pocket scores).
2. Pocket-matched directional AUROCs and `summary_min`.
3. Mixed-library Top-10 composition on EGFR/HER2 (optional, same definition as Table S25).

**Do not** compare “which engine has higher AUROC” as a winner.

## How to write the result (claim freeze)

| Outcome | Allowed sentence |
|---------|------------------|
| Gap remains (neither ≫ directional on EGFR/HER2) | The formulation effect is not Vina-specific under this independent pose-generation protocol. |
| Gap smaller, same sign | Magnitude is engine-dependent; the qualitative formulation effect is partly reproducible. |
| Gap gone or reversed | Apparent dual-target discrimination is scoring/pose-generation dependent. |

Never: “GNINA/smina is better than Vina for dual-target docking.”

## Suggested local command skeleton

```bash
# Example: GNINA docking search (paths are local).
gnina --receptor 3POZ_receptor.pdbqt \
      --ligand EH120_001.pdbqt \
      --center_x ... --center_y ... --center_z ... \
      --size_x ... --size_y ... --size_z ... \
      --exhaustiveness 8 --num_modes 9 --energy_range 3 \
      --seed 20260727 --cpu 8 \
      --out EH120_001_3POZ_gnina_dock.pdbqt
```

After both pockets finish, compute AUROCs with the same `auroc` / bootstrap
helpers as `data/jcim_novelty_v0/scripts/benchmark_formulation_v1.py`
(`N_BOOT = 2000`, seed `20260729`).

Deposit CSVs under:

```text
data/jcim_independent_dock_v0/tables/
```

and point `MASTER_RESULTS_TABLE.csv` at those files only after the numbers
exist. Do not invent placeholders.
