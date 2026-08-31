# POSE_UPLOAD_CHECKLIST — Zenodo 姿态包

> 目标：审稿可复现 top pose + 分数；路径以 `/mnt/d/CADD paper exercise/dual target docking/results/` 为准。

## 必须上传（每面板 × 每受体）

| 面板 | 受体 | top1 pose | 分数表 | 盒子/受体 |
|------|------|-----------|--------|-----------|
| PM48 RDKit | 4L23, 4JT6 | `pik3ca_mtor_panel48_rdkit_v0/poses/{target}/{ligand}/mode_01.pdbqt` | `tables/scores_vina_long.csv`, `ablation_ligand_scores.csv` | `receptors/`, `boxes/` |
| PM48 E8 对照 | 4L23, 4JT6 | `logs/vina_E8/{target}_{lig}_out.pdbqt` | `tables/scores_vina_E8_best.csv` | 同上 |
| PM110 RDKit | 4L23, 4JT6 | `pik3ca_mtor_panel110_rdkit_v0/poses/...` | `tables/scores_vina_long.csv` | 同上 |
| AChE/BChE | 4EY7, 4BDS | `ache_bche_panel_v0/poses/` | `tables/scores_vina_long.csv` | `receptors/`, `boxes/` |
| PIK3CA/PIK3CB | 4L23 (reused), 2WXF | `pik3ca_pik3cb_panel_v0/poses/` | 同上 | 同上 |
| EGFR/HER2 | 3POZ, 3RCD | `egfr_her2_panel120_v0/poses/` | 同上 | 同上 |

## GNINA 重打分姿态

- 输入：Vina `mode_01.pdbqt` → obabel SDF → GNINA minimized SDF
- 表：`tables/scores_gnina_best.csv`, `scores_gnina_long.csv`

## 可选 SI

- E8 vs E16 exhaustiveness 全 9 modes（体积大时可只留 best + log）
- 单靶 enrichment decoy/actives pose 子集（top10 + random10）

## 不要上传

- 未收敛/失败 ligand 的空 pose
- LigPrep 主姿态（仅作 PM48 prep sensitivity delta）
