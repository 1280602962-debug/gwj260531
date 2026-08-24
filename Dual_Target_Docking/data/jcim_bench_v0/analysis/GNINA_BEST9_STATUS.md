# GNINA rescore — STATUS: DONE (best-of-9 CNN minimize/rescore)

**Updated:** after fair all-mode GNINA rescore

## Binary / protocol
`/mnt/d/CADD paper exercise/gnina/bin/gnina` (v1.3.2), CPU `--no_gpu`

- Input: **all** Vina `mode_01`…`mode_09` PDBQT → SDF (Open Babel)
- `gnina --cnn_scoring rescore --minimize --seed 20260727 --cpu 1`
- Per ligand–target: take **max CNNscore** over up to 9 modes (ties → first max)
- Pocket-matched arm still uses `min(score_A, score_B)` on the per-end best-of-9 CNN
- mode_01-only tables retained as `scores_gnina_*_mode01_backup.csv`

## Packs
- **AChE/BChE**: panel=`/home/gwj/repos/gwj260531/Dual_Target_Docking/data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv` best9=96 mode01=96 long=1702
- **PIK3CA/PIK3CB**: panel=`/home/gwj/repos/gwj260531/Dual_Target_Docking/data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv` best9=100 mode01=100 long=1787
- **PIK3CA/mTOR**: panel=`/home/gwj/repos/gwj260531/Dual_Target_Docking/data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv` best9=48 mode01=48 long=864
- **EGFR/HER2**: panel=`/home/gwj/repos/gwj260531/Dual_Target_Docking/data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv` best9=110 mode01=110 long=1978

## mode_01 vs best-of-9 directional AUROC

| pair | summary_min mode01 | summary_min best9 | Δ | frac mode01 wins |
|------|-------------------:|------------------:|--:|-----------------:|
| AChE/BChE | 0.3719 | 0.3585 | -0.0133 | 0.194 |
| PIK3CA/PIK3CB | 0.5064 | 0.4337 | -0.0727 | 0.271 |
| PIK3CA/mTOR | 0.5635 | 0.5952 | +0.0317 | 0.292 |
| EGFR/HER2 | 0.2634 | 0.2645 | +0.0011 | 0.286 |

## Artifacts
- `data/jcim_bench_v0/tables/gnina_mode01_vs_best9_ligand.csv`
- `data/jcim_bench_v0/tables/gnina_mode01_vs_best9_auroc.csv`
- Per-pack: `tables/scores_gnina_long.csv`, `tables/scores_gnina_best.csv`

## Claim update
RTMScore and GNINA now share the same pose coverage (best-of-9 over the same Vina modes).
Three-engine contrast is pose-symmetric; still do **not** claim a universal docking decision rule.
