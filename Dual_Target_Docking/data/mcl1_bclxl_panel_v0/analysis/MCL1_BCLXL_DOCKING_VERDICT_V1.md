# MCL1_BCLXL_DOCKING_VERDICT_V1

Updated: `2026-08-27T04:05:41Z`

## Scope

Frozen ChEMBL θ=6.0 panel96 (24/24/24/24), receptors **3WIY / 3WIZ**, Vina mode-1.
LC6 pose-gold gate role: **`applicability_stress_test`**.

## Completeness

| class | complete ligands (both pockets) |
|-------|--------------------------------:|
| dual | 23 |
| A_only | 24 |
| B_only | 22 |
| neither | 24 |
| **total** | **93** / 96 |

Jobs scored ok: 186 / 192.

## Formulation AUROCs (descriptive)

| contrast | n_pos | n_neg | AUROC | 95% CI |
|----------|------:|------:|------:|--------|
| Dual vs neither (mean) | 23 | 24 | 0.6277 | [0.4618, 0.7857] |
| Dual vs A-only @ Bcl-xL | 23 | 24 | 0.7935 | [0.6549, 0.9146] |
| Dual vs B-only @ MCL1 | 23 | 22 | 0.6087 | [0.4385, 0.7757] |
| **summary_min** | 23 | — | **0.6087** | — |

## Interpretation rules

- If `panel_role=applicability_stress_test`: report honestly; do **not** claim target-general PPI screening performance.
- Homologous BCL-2 fold / BH3-groove domain shift — not a first non-kinase pair (AChE/BChE already is).
- Do not package as external validation.

## Files

- `tables/vina_scores_MBX_v1.csv`
- `tables/formulation_auroc_MBX_v1.csv`
- `tables/ligand_scores_wide_MBX_v1.csv`
- `analysis/MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1.md`
