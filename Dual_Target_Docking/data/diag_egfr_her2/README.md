# EGFR/HER2 对接小面板 panel_v0_40（40 分子）

> 用途：第一枪协议测试（3POZ / 3RCD + Vina）+ whole / moiety / passenger 诊断。  
> **不是**全库虚拟筛选集。  
> 详细操作见 [`../../docs/EGFR_HER2_DIAGNOSTIC_DEMO.md`](../../docs/EGFR_HER2_DIAGNOSTIC_DEMO.md)。

## 组成

| class | N | 定义 | 角色 |
|-------|---|------|------|
| **dual** | 15 | 两端 pChEMBL ≥ 6 | 正例 |
| **A_only** | 12 | EGFR≥6 且 HER2 **测得** &lt;6 | 硬负（主攻击） |
| **B_only** | 10 | HER2≥6 且 EGFR **测得** &lt;6 | 硬负 |
| **neither** | 3 | 两端都测且都 &lt;6 | 底噪 |
| **合计** | **40** | | |

来源：`../public_pair_selection/egfr_her2_fourclass_chembl_ids.csv` 分层抽样（按 `min_pchembl` 跨档取样）+ 强制纳入锚点。

## 必含锚点

| ChEMBL ID | 名称 | class | EGFR / HER2 pChEMBL | 说明 |
|-----------|------|-------|---------------------|------|
| CHEMBL1614725 | **TAK-285** | dual | 9.00 / 8.52 | 双端共晶 3POZ/3RCD；姿态 QC + moiety 金标准 |
| CHEMBL554 | **Lapatinib** | dual | 10.22 / 8.80 | 经典双 EGFR/HER2 TKI |
| CHEMBL483321 | CP-724714 | dual | 8.19 / 8.00 | 文献双 TKI 之一 |

## 对接怎么用这张表

1. 读 `panel_v0_40.csv`（已含 canonical SMILES）  
2. 每个分子 → **EGFR(3POZ)** 与 **HER2(3RCD)** 各对接一次（整分子）  
3. 先只对 **TAK-285** 算自对接 RMSD；过关再跑全面板  
4. 有把握的分子再补 moiety 标注；切不清的先标 `ambiguous` 不进主结论  

任务数：40 × 2 = **80 次对接**（另加 TAK-285 QC）。

## 文件

| 文件 | 内容 |
|------|------|
| `panel_v0_40.csv` | 40 分子名单 + SMILES + 活性 + 角色注释 |

字段：`panel_id, molecule_chembl_id, class, pchembl_EGFR, pchembl_HER2, min_pchembl, pref_name, smiles, inchi_key, max_phase, role_note`

## 已知局限（用前必读）

1. **标签是 ChEMBL 操作定义**，不是「每个都是教科书 EGFR/HER2 双 TKI」。抽样里可能混入交叉标注噪声（已剔除明显离谱的 imatinib `CHEMBL941`）。  
2. 命名药里若出现非典型 EGFR/HER2 化学型（如部分激酶工具药），对接前可再人工 QC，标 `docking_eligible=false`。  
3. 本面板 **尚未** 完成 moiety 原子划分；先跑 whole-mol 协议，再对锚点+可标注子集做乘客分析。  
4. 同源 ATP 口袋：本面板适合 **跑通流程**；乘客假说是否成立还需异质对复验（见 `docs/CRITIQUE_AND_NEXT_STEPS.md`）。

## 一页清单

- [x] 40 分子四类面板已导出  
- [x] 含 TAK-285 / lapatinib  
- [x] SMILES 齐全  
- [ ] TAK-285 自对接 RMSD  
- [ ] 全面板 80 次对接  
- [ ] moiety 标注（至少锚点）  
- [ ] whole vs moiety 诊断表  

Related freeze pack: `../egfr_her2_panel40_v0/` (uploaded 2026-07-27).
