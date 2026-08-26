# Independent pose-generation result (local GNINA docking)

Engine: GNINA 1.3.2 docking search (`minimizedAffinity` mode-1), seed 20260727.
Not a Vina-pose rescore.

## EGFR/HER2 (primary)

| Formulation | Vina (frozen) | GNINA independent dock |
|-------------|---------------|------------------------|
| Dual vs neither (mean) | 0.756 | **0.783** |
| Directional `summary_min` | 0.430 | **0.220** |

**Verdict: gap remains** (neither ≫ directional).

Allowed sentence:
> The formulation effect is not Vina-specific under this independent pose-generation protocol.

## PIK3CA/mTOR (secondary)

| Metric | Vina | GNINA dock |
|--------|------|------------|
| `summary_min` | 0.692 | 0.633 |

## Notes

- Failed docks: `EH120_109` (both pockets), `PM48_19` (both pockets). Analyses use complete ligands only (n=109 and n=47).
- Tables: `tables/independent_dock_*.csv`, `tables/gnina_dock_scores_*.csv`.
