# 本地 Agent 操作单：PI-103 QC No-Go 之后（升 E=16 → Phase B′ → 全面板）

> **角色：** 本地对接 agent（有 `/mnt/d/CADD paper exercise/...` 与 Vina 环境）  
> **状态：** Phase A 完成；Phase B @E=8 = **No-Go**（4JT6 best9≈5.0 Å）；E=16 诊断已显示 4JT6 best9≈0.445 Å  
> **决策（已冻结）：** 不自动全面板；走 **轻量 sanity → 书面 E=16 → Phase B′ → Go 才 48×2**  
> **权威短决策：** [`PIK3CA_MTOR_COGNATE_QC_NEXT_STEP.md`](PIK3CA_MTOR_COGNATE_QC_NEXT_STEP.md)  
> **原 SOP：** [`PIK3CA_MTOR_PANEL48_LOCAL_DOCKING_SOP.md`](PIK3CA_MTOR_PANEL48_LOCAL_DOCKING_SOP.md)

---

## 0. 硬约束（每次动手前默念）

1. **禁止**在当前 E=8 No-Go 下启动 48×2。  
2. **禁止**修改 EGFR/HER2 panel40 的 `exhaustiveness=8`。  
3. **禁止**重新 LigPrep / 改面板分子名单（除非 Phase A 映射发现错误）。  
4. **禁止**为刷分把全面板升到 E=32；E=32 仅允许 `PM48_01@4JT6` 诊断。  
5. 升 E 的唯一合法理由：**cognate 采样 QC**（已满足：E=8 失败、E=16 成功）。  
6. 工作根固定：

```text
ROOT=/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0
```

7. Vina 二进制：与 EGFR 敏感性实验相同路径（本机曾用 `/home/gwj/miniconda3/bin/vina`）；若不同，写入 `protocol/paths.json` 并全程只用一个。

---

## 1. 目录与已有产物（先盘点，勿覆盖历史）

```bash
cd "$ROOT"
ls -la protocol/ receptors/ boxes/ ligands_pdbqt/ analysis/cognate_redock_v0/
```

**必须已存在（Phase A/B）：**

| 路径 | 用途 |
|------|------|
| `ligands_pdbqt/PM48_01.pdbqt` … `PM48_48.pdbqt` | 对接输入 |
| `receptors/4L23_receptor.pdbqt` | PIK3CA |
| `receptors/4JT6_receptor.pdbqt` | mTOR |
| `boxes/4L23_box.json` / `boxes/4JT6_box.json` | 盒子 |
| `tables/4L23_cocrystal_X6K.pdb` / `tables/4JT6_cocrystal_X6K.pdb` | RMSD 参考 |
| `analysis/cognate_redock_v0/COGNATE_QC_VERDICT.md` | E=8 No-Go 裁决 |
| `protocol/protocol.yaml` | 当前仍可能写着 E=8 |

**新建（本轮）：**

```text
analysis/cognate_redock_v0/
  SANITY_4JT6.md
  tables/pm48_01_rmsd_E16.csv
  COGNATE_QC_VERDICT_E16.md
protocol/
  EXHAUSTIVENESS_FREEZE_v0_1.md
poses/cognate_E16/          # 与 E=8 poses 分开，勿覆盖
logs/cognate_E16/
```

E=8 的 poses/logs **只读归档**，不要删、不要覆盖。

---

## 2. Step 0 — 4JT6 轻量 sanity（必须先做，≤30–60 min）

### 2.1 检查清单 → 写 `analysis/cognate_redock_v0/SANITY_4JT6.md`

对每一项写 **PASS / FAIL / 备注**：

| ID | 检查 | 如何做 | 失败时 |
|----|------|--------|--------|
| S1 | 受体是 mTOR 激酶 ATP 位点，不是 FRB/FKBP | 读 prepared PDB 标题/链；确认有激酶域残基，无 rapamycin 故事 | 换回正确 prepared；重做 Phase B@E=8 |
| S2 | 盒子来自 **同坐标系** X6K | `4JT6_cocrystal_X6K.pdb` 坐标应落在 receptor 口袋内；box center 应靠近 X6K 质心 | 从正确复合物重提 X6K + 重建 box |
| S3 | PM48_01 = PI-103 | 与面板 CSV / InChIKey `TUVCWJQQGGETHL-UHFFFAOYSA-N` 一致；与 X6K 重原子数匹配 | 修正配体映射，重导出 PDBQT |
| S4 | 盒子尺寸合理 | 读 `boxes/4JT6_box.json`：padding 5 Å，min edge≥20 | 重建 box（定义不变） |
| S5 | E=16 诊断可复现线索 | 确认已有 E=16@seed20260727 的 log/pose，或本轮将重跑 | — |

### 2.2 分支

```text
若 S1–S4 任一项 FAIL
  → 修正后：仅用 E=8、seed=20260727 重跑 PM48_01×两端 Phase B
  → 两端 best9<2 → Go → 全面板仍用 E=8
  → 若仍 4JT6 失败 → 再进入 Step 1（升 E=16）

若 S1–S4 全 PASS（预期）
  → 进入 Step 1（书面 E=16）
```

**不要**在 sanity 阶段重做 Protein Prep Wizard 全流程，除非明确发现错结构。

---

## 3. Step 1 — 书面冻结 exhaustiveness=16

### 3.1 更新 `protocol/protocol.yaml`

把对接相关字段改为（保留 Maestro 原路径记录）：

```yaml
freeze_id: pik3ca_mtor_panel48_v0
engine: AutoDock_Vina
vina_version: "1.2.7"
seed_policy: fixed_global
seed_fixed_global: 20260727
n_modes: 9
energy_range: 3
exhaustiveness: 16
exhaustiveness_v0_1: 16
exhaustiveness_note: >
  Pair-specific. EGFR/HER2 panel40 remains exhaustiveness=8.
  Raised because PI-103@4JT6 E=8 best_of_9 failed (~5.0 A) while
  E=16 recovered best_of_9 ~0.445 A (seed 20260727).
cognate_qc_gate: "both ends best_of_9 heavy-atom RMSD < 2.0 A"
full_panel_blocked_until: "COGNATE_QC_VERDICT_E16.md says Go"
```

### 3.2 写 `protocol/EXHAUSTIVENESS_FREEZE_v0_1.md`

必须含：

- 一句结论：`exhaustiveness_v0_1 = 16`（仅 DTPAIR-01）  
- E=8 / E=16 数字表（从已有 VERDICT 抄）  
- 明确：`egfr_her2_panel40_v0` **不改**  
- 签字等价：agent 运行时间戳 + hostname

### 3.3 更新 `protocol/SEED_POLICY.md`（若无则新建）

```text
seed_fixed_global = 20260727
noise_seeds_for_cognate_only = [7, 42]
full_panel_seed = 20260727 only
```

---

## 4. Step 2 — Phase B′：E=16 cognate 确认

### 4.1 作业表（只跑这些，先别全面板）

| job_id | target | ligand | exhaustiveness | seed | 输出目录 |
|--------|--------|--------|----------------|------|----------|
| B16_A | 4L23 | PM48_01 | 16 | 20260727 | `poses/cognate_E16/4L23/PM48_01_seed20260727/` |
| B16_B | 4JT6 | PM48_01 | 16 | 20260727 | `poses/cognate_E16/4JT6/PM48_01_seed20260727/` |
| B16_C | 4JT6 | PM48_01 | 16 | 7 | `poses/cognate_E16/4JT6/PM48_01_seed7/` |
| B16_D | 4JT6 | PM48_01 | 16 | 42 | `poses/cognate_E16/4JT6/PM48_01_seed42/` |

共 **4** 个 Vina 作业。`n_modes=9`，`energy_range=3`，`cpu=1`。

### 4.2 conf 文件模板

对每个 job 写 `logs/cognate_E16/confs/<job_id>.txt`：

```text
receptor = <ROOT>/receptors/<TARGET>_receptor.pdbqt
ligand = <ROOT>/ligands_pdbqt/PM48_01.pdbqt
center_x = <from boxes/<TARGET>_box.json>
center_y = ...
center_z = ...
size_x = ...
size_y = ...
size_z = ...
exhaustiveness = 16
num_modes = 9
energy_range = 3
cpu = 1
seed = <SEED>
out = <pose_dir>/PM48_01_all_modes.pdbqt
```

运行：

```bash
vina --config logs/cognate_E16/confs/<job_id>.txt \
  > logs/cognate_E16/<job_id>.log 2>&1
```

### 4.3 跑后检查（每个 job）

- [ ] log 含 `seed` / `random seed` 等于设定值  
- [ ] `PM48_01_all_modes.pdbqt` 存在  
- [ ] 拆成 `mode_01.pdbqt` … `mode_09.pdbqt`（可用 EGFR 敏感性脚本同款 `split_models`）  
- [ ] 恰好 9 个 mode（若不足，记 FAIL 并重跑该 job）

### 4.4 RMSD（定义冻结，与 EGFR 一致）

参考：`tables/<TARGET>_cocrystal_X6K.pdb`  
查询：各 `mode_XX.pdbqt`（转 SDF/MOL 后算）  
规则：

- **重原子 only**  
- 模板约束下 **min CalcRMS**（对称校正）  
- 不做蛋白叠合  

写出：`analysis/cognate_redock_v0/tables/pm48_01_rmsd_E16.csv`

```text
target,seed,exhaustiveness,rmsd_mode1,rmsd_best_of_9,best_of_9_mode,pass_mode1_lt2,pass_best_of_9_lt2
4L23,20260727,16,...
4JT6,20260727,16,...
4JT6,7,16,...
4JT6,42,16,...
```

### 4.5 写 `COGNATE_QC_VERDICT_E16.md`（必须含一句 Go/No-Go）

**Go 条件（同时满足）：**

1. `4L23` @ seed `20260727`：`best_of_9 < 2.0`  
2. `4JT6` @ seed `20260727`：`best_of_9 < 2.0`  
3. 噪声：`4JT6` 在 seeds `{20260727,7,42}` 中 **≥2/3** 的 `best_of_9 < 2.0`（含主 seed）

**仍算 Go 但必须注明：**

- `rmsd_mode1 ≥ 2` 但 `best_of_9 < 2`（预期 4JT6 mode1 仍可能 ~7 Å）  
- → 全面板 **必须** 输出 9 modes；后续 **必须** 计划 RTM best-of-9

**No-Go：**

- 主 seed 任一端 best9≥2 → 不要开全面板；按下文「失败分支」  
- 仅主 seed 过、7 与 42 都≥2 → 视为不稳：可加跑 `PM48_01@4JT6` E=32 诊断；仍不稳则考虑 4JT5 对照，**不改面板**

---

## 5. Step 3 — 仅 Go 之后：全面板 48×2 @ E=16

### 5.1 启动前门禁（全部打勾才跑）

- [ ] `COGNATE_QC_VERDICT_E16.md` 首行或显式写有 **`Verdict: Go`**  
- [ ] `protocol.yaml` 中 `exhaustiveness: 16`  
- [ ] E=8 No-Go 文件仍保留（历史）  
- [ ] 48 个 `ligands_pdbqt/PM48_XX.pdbqt` 齐全  

### 5.2 作业矩阵

```text
ligands: PM48_01 .. PM48_48
targets: 4L23, 4JT6
exhaustiveness: 16
seed: 20260727          # 全面板只用这一个
n_modes: 9
energy_range: 3
total_jobs: 96
```

输出布局（与 panel40 对齐，但标 E16）：

```text
poses/4L23/<PM48_XX>/mode_01.pdbqt ... mode_09.pdbqt
poses/4L23/<PM48_XX>/<PM48_XX>_all_modes.pdbqt
poses/4JT6/<PM48_XX>/...
logs/vina/4L23_<PM48_XX>.log
logs/vina/4JT6_<PM48_XX>.log
logs/vina_confs/...
```

可用并行（多作业），但 **每作业 cpu=1**。建议写 `scripts/run_panel48_vina_E16.py`（仿 EGFR `run_exhaustiveness_sensitivity.py`），并产出：

- `tables/job_status.csv`  
- `tables/scores_vina_long.csv`（每 mode 一行）  
- `tables/scores_vina.csv`（每配体汇总）

### 5.3 `scores_vina.csv` 最低列

```text
panel_id,class,vina_4L23_mode1,vina_4JT6_mode1,vina_mean,vina_min,vina_delta,
best_mode_4L23,best_mode_4JT6
```

（Vina 分保持原生：越负越好。）

### 5.4 跑完验收

- [ ] 96/96 jobs `OK`  
- [ ] 每 log 可 grep 到 seed `20260727`  
- [ ] 每配体每靶 9 modes  
- [ ] `MANIFEST.md` 更新：`exhaustiveness=16`，cognate Go 链接  

### 5.5 全面板后（同一 agent 可续，或停等指令）

1. RTMScore best-of-9（若环境可用）→ `rtm_mean` / `rtm_min` / `rtm_min_z`  
2. `tables/ablation_metrics.csv`：AUROC Dual vs A∪B；Top10 dual / hardneg  
3. `tables/hardneg_cases.md`：至少 alpelisib/taselisib、AZD-8055/Ku-0063794  

**若 RTM 暂不可用：** 先交付 Vina 表 + cognate QC 包，在 MANIFEST 标明 RTM pending。

---

## 6. 失败分支（Agent 决策树）

```text
Step0 FAIL (硬错误)
  → 修正结构/盒子 → Phase B @ E=8 → 过则全面板 E=8；不过再升 E=16

Step2 主 seed 4JT6 best9≥2 @ E=16
  → 复查 box/X6K；重跑一次 B16_B
  → 仍失败：E=32 仅 PM48_01@4JT6
  → 仍失败：停全面板；报告并建议 4JT5 对照 QC

Step2 主 seed 过但 7&42 都失败
  → E=32 噪声三 seed；写不稳说明
  → 若 E=32 稳定：可书面冻结 E=32（需用户确认）；默认先报告等指令

Step3 中途作业失败
  → 只重跑失败 job；不要改 E/seed
```

---

## 7. 交付清单（本轮结束时必须有）

**无论 Go / No-Go：**

- [ ] `SANITY_4JT6.md`  
- [ ] `EXHAUSTIVENESS_FREEZE_v0_1.md`  
- [ ] `pm48_01_rmsd_E16.csv`  
- [ ] `COGNATE_QC_VERDICT_E16.md`（含 **Verdict: Go|No-Go**）  
- [ ] 4 个 cognate E16 logs + poses  

**仅 Go 且完成全面板时另加：**

- [ ] 96 vina logs/poses  
- [ ] `scores_vina.csv` / `scores_vina_long.csv` / `job_status.csv`  
- [ ] 更新 `MANIFEST.md` / `protocol.yaml`  

可选：把 `analysis/cognate_redock_v0/` 与 `protocol/` 小文件同步回 GitHub 仓库对应路径（大体量 poses 可只留本地 + MANIFEST 写绝对路径）。

---

## 8. 给本地 Agent 的可粘贴任务块

```text
工作根：
/mnt/d/CADD paper exercise/dual target docking/results/pik3ca_mtor_panel48_v0

任务：执行 PIK3CA/mTOR panel48 在 PI-103 cognate E=8 No-Go 之后的升 E 路径。
严格按仓库 Dual_Target_Docking/docs/PIK3CA_MTOR_AGENT_OPS_E16.md。

顺序：
1) 写 SANITY_4JT6.md（轻量检查，不大改蛋白）
2) 若无硬错误：冻结 protocol exhaustiveness=16（勿改 EGFR E=8）
3) 只跑 PM48_01：4L23@20260727、4JT6@20260727/7/42，全部 E=16,n_modes=9
4) 按重原子 min-CalcRMS 写 pm48_01_rmsd_E16.csv
5) 写 COGNATE_QC_VERDICT_E16.md；仅当 Verdict=Go 才启动 48×2@E=16@seed20260727
6) 禁止在 Go 前开全面板；禁止覆盖 E=8 历史 poses

完成后回报：Verdict、两端 best9 数字、噪声结果、是否已启动/完成全面板。
```

---

## 9. 一句话

**先 sanity → 书面 E=16 → 只对 PI-103 做 E=16 双端+噪声确认 → Go 才 96 作业全面板；全程固定 seed `20260727`，全面板不要混 E。**
