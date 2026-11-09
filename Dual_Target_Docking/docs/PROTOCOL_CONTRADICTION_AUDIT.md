# Protocol contradiction audit

Date: 2026-09-18  
Rule: if two versions exist, exactly one is `current primary`; others are `sensitivity` or `historical/superseded`.

## EGFR/HER2

| Topic | Version A | Version B | Classification after this audit |
|-------|-----------|-----------|----------------------------------|
| Ligand prep | Historical yaml: mixed LigPrep/RDKit | CASE 2 uniform RDKit rebuild `egfr_her2_uniform_rdkit_v1` | Deposited production scores: **not_recoverable**. Rebuild: not yet primary. |
| Box | `3POZ_box.json` / `3RCD_box.json` | `*_box_corrected.json` | Corrected = **current primary**. Uncorrected = superseded. |
| Five-seed scores | `multiseed_scores_long_v1.csv` EGFR rows | `scores_vina_mode1_EGFR_corrected_box_fiveseed.csv` | Corrected-box file **overrides** in `compute_canonical_results._load_five_seed_scores`. |
| Production scores | `ablation_ligand_scores.csv` | uniform rebuild tables (not yet written) | Ablation file is still `EGFR_SCORE_SOURCE`. Rebuild must not look like current Table 2. |
| Box JSON status | was `corrected_not_yet_canonical` | now `current_primary_corrected_heavy_atom` | Stale status field fixed this round. |

## PIK3CA/mTOR

| Topic | Version A | Version B | Classification |
|-------|-----------|-----------|----------------|
| Ligand prep | Predecessor LigPrep campaign (local log; `prep_delta_vs_ligprep.csv`) | `ligand_prep: rdkit_etkdg_meeko` + `ablation_ligand_scores.csv` | RDKit = **current primary**. LigPrep = superseded. Protocol note neutralized this round (no longer states LigPrep as fact for current receptors/boxes/E). |
| Exhaustiveness | E=16 primary | E=8 `scores_vina_E8_best.csv` | E=16 = primary. E=8 = sensitivity (`protocol_sensitivity.csv`). |
| Panel size | PM48 | PM110 | PM48 = primary. PM110 = sensitivity. |
| Receptors | 4L23 / 4JT6 | 4JPS, 5DXT, 4JSX | Primary vs receptor-substitution sensitivity. |

## AChE/BChE

| Topic | Version A | Version B | Classification |
|-------|-----------|-----------|----------------|
| Labels | protocol `label_rule: strict_6.5_5.5` | primary θ=6.0 adjudication | Construction quota vs evaluation threshold. **Not two current class definitions.** Yaml now states `label_rule_scope: panel_construction_quota_only`. |
| Receptors | 4EY7 / 4BDS | 5DYW, 6QAA, 6ZWI | Alt cognate-fail receptors = sensitivity only. |
| Five-seed n | Table 2 n=95 | five-seed production seed n=94; intersection n=88 | Per-seed complete-case is the current five-seed realization. Fixed-membership is QA-only (`results/qa/five_seed_fixed_membership_sensitivity.csv`). |

## Track-B five pairs

| Topic | Version A | Version B | Classification |
|-------|-----------|-----------|----------------|
| Activity | `chembl37_dump_panel` | Track-A assay adjudication | Dump panel is current for the five pairs. Must not be called assay-adjudicated. |
| GNINA | Independent search JAK1/TYK2 only | CNN best-of-9 rescore of Vina poses | Independent = Fig 5A subset. CNN = Table S7 sensitivity. |
| RTMScore | best-of-9 on Vina poses | n/a | Table S7 sensitivity, not independent docking. |
| README poses | previously claimed poses not in slim tree | `poses/` has 9742 PDBQT | README corrected this round. |

## Exhaustiveness / seed / n_modes / score convention

Shared current primary: Vina 1.2.7, `score_S = −mode1`, n_modes=9, production seed 20260727, five-seed list 20260727/20260811–14, D vs A uses pocket B, D vs B uses pocket A, `summary_min` = min of those arms, bootstrap B=2000 seed 20260729.

Exception: PIK3CA/mTOR primary **E=16**; all other primary pairs **E=8**. That is documented protocol sensitivity, not a hidden second primary.
