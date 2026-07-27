# PIK3CA/mTOR：小分子准备与对接数量冻结

> 对照 EGFR/HER2 `panel40_v0`（15 dual / 12 A_only / 10 B_only / 3 neither）。  
> **不要一上来对接 ChEMBL 全部 2713 配对分子。**

---

## 1. 到底对接多少？

**现在先做档 1：已冻结为 `pik3ca_mtor_panel48_v0`（N=48）。**  
清单：[`../data/pik3ca_mtor_panel48_v0/tables/panel_v0_48.csv`](../data/pik3ca_mtor_panel48_v0/tables/panel_v0_48.csv) · 说明：[`../data/pik3ca_mtor_panel48_v0/README.md`](../data/pik3ca_mtor_panel48_v0/README.md)

| 档 | 名称 | 分子数 N | 对接任务 | 用途 |
|----|------|----------|----------|------|
| **1（已冻结）** | 协议迁移面板 | **48** | 96（两端） | 与 EGFR/HER2 同量级比失败模式；先冻协议 |
| **2（主表）** | 外推主面板 | **120–200** | 240–400 | JCIM 级主结果；四类更稳 |
| **3（后期）** | 大规模基准 | 500–2713 | 千级 | 全量/分层；协议锁死后才做 |

### 档 1 组成（已落地）

| 类 | 定义 | n | 必含 |
|----|------|---|------|
| dual | 两端 pChEMBL ≥ 6 | 18 | **PI-103 = PM48_01** |
| A_only | PIK3CA≥6 且 mTOR **测得** &lt;6 | 14 | Alpelisib, Taselisib 等 |
| B_only | mTOR≥6 且 PIK3CA **测得** &lt;6 | 12 | AZD-8055, Ku-0063794, WYE-132, OSI-027 |
| neither | 两端都测且都 &lt;6 | 4 | 底噪 |

### 明确不做

- 不对接全部 2002 dual / 2713 paired  
- 不把 **PROTAC / 降解剂** 混进主面板（另表）  
- 不把「只测过一端」的分子标成 A_only/B_only  
- 不为凑数加入未测 decoy（档 3 再说）

---

## 2. 小分子怎么准备（与 EGFR/HER2 对齐）

目标：**同一套配体准备，两端同一构象输入**，只换受体/盒子。

### 2.1 选分子（准备之前）

1. 从 `mols_PIK3CA.json` / `mols_MTOR.json` 做配对四类表（阈值 pChEMBL ≥ 6）。  
2. 强制纳入：**PI-103**（ChEMBL 查 ID；PDB chem_comp **X6K**）。  
3. 优先经典 ATP 双抑制剂；排除明显 PROTAC（大 MW、长 linker、E3 binder 特征）。  
4. 记录：`panel_id`, ChEMBL ID, SMILES, InChIKey, class, pChEMBL_A/B, scaffold, role_note。  
5. `architecture`：有文献再标 merged/linked；否则写 `unknown`，**不要猜**。

### 2.2 LigPrep（与 panel40 同一工具链）

| 项 | 冻结建议 |
|----|----------|
| 工具 | Schrödinger **LigPrep**（与 EGFR/HER2 一致） |
| 输入 | 面板 CSV 的 canonical SMILES（盐剥离后母核） |
| pH | 7.0 ± 0.5（与蛋白准备一致） |
| 立体 | 指定立体保留；未指定则生成合理对映/非对映（**每分子最终只保留 1 个代表构象进对接**，多构象另做敏感性） |
| 互变异构 | LigPrep 默认；记录选用态 |
| 输出 | `.maegz` → RDKit/Meeko → **每个 panel_id 一个 PDBQT** |
| 禁止 | 档 1 中途改 LigPrep 参数；敏感度另开版本号 |

可选对照：Open Babel / RDKit ETKDG 只作「准备器敏感性」，不作主面板。

### 2.3 对接输入规则

- 每个分子：**整分子** 一份 3D → 对 **4L23** 与 **4JT6** 各跑一次（任务数 = 2N）。  
- 全局 seed：`20260727`（与协议冻结一致）；exhaustiveness 用 v0.1 已选定值（默认先跟 panel40 的 E=8）。  
- n_modes=9，energy_range=3。  
- **不要**在敏感度里重新 LigPrep；只改 exhaustiveness/seed。

### 2.4 姿态 QC（开面板打分前）

- **PI-103** 双端 cognate/self-dock：重原子 RMSD **< 2 Å**（Vina top1 或 RTM-best，与 EGFR 侧记录方式一致）。  
- 不过关：先改盒子/质子化/蛋白准备，**不动分子面板名单**。

---

## 3. 和「第一次 40 个」的关系

| | EGFR/HER2 panel40 | PIK3CA/mTOR 现在 |
|--|-------------------|------------------|
| 角色 | 诊断 demo | **同尺寸协议迁移**，不是立刻全规模 |
| N | 40 | **40–60（推荐 48）** |
| 金标准配体 | TAK-285 | PI-103 |
| 结构 | 3POZ / 3RCD | 4L23 / 4JT6 |
| 通过后再扩 | — | → 档 2（120–200） |

**一句话：** 第二对先再做一张「约 40–50 分子」的四类面板验证尺子能否外推；主文规模表放到协议锁死之后的 120–200，而不是现在就上千。
