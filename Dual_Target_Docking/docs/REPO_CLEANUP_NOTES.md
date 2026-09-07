# Repo cleanup notes（Dual_Target_Docking）

> 日期：2026-07-30  
> 目的：去掉前期探索/废弃路线，只保留 DualFourClass JCIM 评测文所需材料。

## 已从 git 删除

### 数据 / 运行
- `data/cross_pair_bootstrap_v0/` — 早期 bootstrap，已被 `jcim_bench_v0` 取代  
- `data/diag_egfr_her2/` — 40 分子诊断面板入口（分数已并入 panel packs）  
- `data/dual_target_structures/` — 早期共晶目录草稿  
- `data/jcim_feasibility_v0/` — 可行性试算，结论已进 j0j1 / bench  
- `data/literature/` — 文献 CSV 大目录（NMI/调研用）  
- `data/plan_v2_redteam_v0/` — 红队试算  
- `data/track_a_starter_v0/` — Track A 起稿包（已被 strengthen / bench 取代）  
- `data/schema/` — 旧模板  
- `research_runs/` — 方向重定临时 scoop  

### 体积瘦身（保留协议与分数表）
- `egfr_her2_panel40_v0/poses/`, `logs/`  
- `pik3ca_mtor_panel48_v0/poses/`, `logs/`, `ligands_pdbqt/`, `ligands_sdf/`  

姿态仍在本地 results 目录；Zenodo 时再上传 top1。

### 文档（废弃主张 / 过时 agent 命令 / 调研）
- 全部 `NMI_*`、乘客 moiety 方案  
- EGFR 扩面 / Stage M / J0 旧 agent 命令（已完成）  
- 文献大目录 md、早期 critique、SOP 流水账、plan v1/v2 全文等  

## 明确保留

- K=4 面板分数 + 受体/盒子 + cognate QC 记录  
- `jcim_bench_v0` + `jcim_strengthen_t0t1_v0` + `jcim_j0j1_v0`  
- `stage_m_v0`（测量审计与 Track B=Weak）  
- `public_pair_selection` + `protocols`  
- 现行 JCIM 路线 / P0 指南 / CLAIM_CEILING  

## 未动

- Dual_Target_Docking/ 不导入、不依赖仓库根目录下的其他课题。根目录若仍有其他项目目录，与本评测文无关。
