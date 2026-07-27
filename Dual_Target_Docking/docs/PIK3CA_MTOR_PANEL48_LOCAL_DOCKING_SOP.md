# 本地操作单：`pik3ca_mtor_panel48_v0` 对接（先重对接 QC，再全面板）

> 给本地 agent / 本机执行。  
> 对齐 EGFR/HER2 `panel40_v0` + 已冻结协议：`seed=20260727`，`exhaustiveness_v0_1=8`，`n_modes=9`。  
> **硬门槛：PI-103 双端 cognate/self-dock 未过关前，禁止开全面板 48×2。**

---

## 0. 已给定的输入（不要另找文件）

| 角色 | 路径 |
|------|------|
| 配体 LigPrep 输出 | `D:\CADD paper exercise\dual target docking\Maestro doc\pik3ca_mtor_panel48_v0_ligprep\pik3ca_mtor_panel48_v0_ligprep-out.maegz` |
| PIK3CA 蛋白 | `D:\CADD paper exercise\dual target docking\Maestro doc\4L23_PIK3CA_prepared.pdb` |
| mTOR 蛋白 | `D:\CADD paper exercise\dual target docking\Maestro doc\4JT6_mTOR_prepared.pdb` |
| 面板名单 | 仓库 `Dual_Target_Docking/data/pik3ca_mtor_panel48_v0/tables/panel_v0_48.csv` |
| 姿态金标准 | **PM48_01 = PI-103**（CHEMBL573339；PDB chem_comp **X6K**） |

建议工作根目录：

```text
D:\CADD paper exercise\dual target docking\results\pik3ca_mtor_panel48_v0\
```

建议 Linux/WSL 挂载写法（若与 EGFR 一致）：

```text
/mnt/d/CADD paper exercise/dual target docking/...
```

---

## 1. 协议冻结值（直接照抄，不要改）

来自 EGFR/HER2 exhaustiveness 敏感性裁决（`exhaustiveness_v0_1=8`）：

| 项 | 值 |
|----|-----|
| 引擎 | AutoDock Vina **1.2.7** |
| seed | **`20260727`**（`fixed_global`；每作业同一 seed） |
| exhaustiveness | **8** |
| n_modes | **9** |
| energy_range | **3** |
| cpu_per_job | 1（可并行多作业，但单作业 1 CPU） |
| 盒子定义 | 共晶配体 **X6K** 的 AABB + **5 Å** padding；任一边 **min edge = 20 Å** |
| RMSD 门槛（重对接） | 重原子；**best_of_9 < 2.0 Å** 两端都过才算 QC 通过 |
| 主报告指标 | 同时记录 `rmsd_mode1` 与 `rmsd_best_of_9`（及 RTM-best，若可用） |

**不要**在本轮改 LigPrep、不要为“好看分数”升 exhaustiveness。

---

## 2. 目录骨架（先建）

在 `results/pik3ca_mtor_panel48_v0/` 下：

```text
protocol/
  protocol.yaml
  SEED_POLICY.md
receptors/
boxes/
ligands_sdf/
ligands_pdbqt/
poses/
  4L23/
  4JT6/
logs/
  vina/
  vina_confs/
  cognate_qc/
tables/
analysis/
  cognate_redock_v0/
```

把两份 prepared PDB 复制进 `receptors/`，并保留原 Maestro 路径记录在 `protocol.yaml`。

---

## 3. Phase A — 输入检查与格式转换（全面板前）

### A1. 核对 maegz 与面板 ID

1. 用 RDKit `MaeMolSupplier`（或 Maestro 导出）打开  
   `pik3ca_mtor_panel48_v0_ligprep-out.maegz`
2. 每个结构标题 / `s_m_title` / 自定义属性必须能映射到 **`PM48_01`…`PM48_48`**
3. 写出：
   - `tables/ligand_input_manifest.csv`  
     列：`panel_id,chembl_id,pref_name,n_confs_in_maegz,chosen_conf_index,sdf_path,pdbqt_path,notes`
4. **每个 panel_id 只保留 1 个对接输入构象**（与 EGFR panel40 相同；多构象不要偷偷全开）
5. 确认 **PM48_01** 存在且可识别为 PI-103

若 maegz 里某 ID 缺失 / 一对多无法消歧 → **停**，先修映射，不对接。

### A2. 导出 SDF + PDBQT

对每个 `PM48_XX`：

1. `ligands_sdf/PM48_XX.sdf`
2. Meeko（或与 EGFR 相同的转换链）→ `ligands_pdbqt/PM48_XX.pdbqt`

受体：

```bash
# 示例（路径按本机调整）
mk_prepare_receptor.py --read_pdb receptors/4L23_PIK3CA_prepared.pdb -o receptors/4L23_receptor.pdbqt
mk_prepare_receptor.py --read_pdb receptors/4JT6_mTOR_prepared.pdb -o receptors/4JT6_receptor.pdbqt
```

### A3. 从共晶配体建盒子（必须用 X6K，不要用手点）

1. 从原始/prepared 复合物中提取 **X6K**（PI-103）坐标：
   - `tables/4L23_cocrystal_X6K.pdb`
   - `tables/4JT6_cocrystal_X6K.pdb`
2. 盒子 = 配体重原子 AABB + 5 Å；边长不足 20 则扩到 20
3. 写出：
   - `boxes/4L23_box.json`
   - `boxes/4JT6_box.json`
   - `boxes/all_boxes.json`

JSON 字段与 EGFR 一致：`center_x/y/z`, `size_x/y/z`, `n_ligand_atoms`。

**检查：** prepared PDB 若已去掉配体，必须从**去配体前的同一坐标系复合物**取 X6K；禁止把配体从别的结构硬贴过来。

---

## 4. Phase B — PI-103 双端重对接 QC（先做，只做 PM48_01）

### B1. 只跑 2 个作业

| job | receptor | ligand | seed | E | n_modes |
|-----|----------|--------|------|---|---------|
| 1 | 4L23 | PM48_01 | 20260727 | 8 | 9 |
| 2 | 4JT6 | PM48_01 | 20260727 | 8 | 9 |

输出：

```text
poses/4L23/PM48_01/mode_01.pdbqt … mode_09.pdbqt + PM48_01_all_modes.pdbqt
poses/4JT6/PM48_01/...
logs/vina/4L23_PM48_01.log
logs/vina/4JT6_PM48_01.log
```

每个 log 必须能 grep 到 `random seed: 20260727`（或 Vina 打印的等价字段）。

### B2. 算 RMSD（定义冻结）

对齐 EGFR `rmsd_definition.md`：

- 参考：`tables/4L23_cocrystal_X6K.pdb` / `tables/4JT6_cocrystal_X6K.pdb`
- 原子：**重原子 only**
- 对称：模板约束下 **min CalcRMS**（不要裸 GetBestRMS 糊弄）
- 不做蛋白重叠合；坐标系必须与 docking 受体一致

写出：`analysis/cognate_redock_v0/tables/pm48_01_rmsd.csv`

建议列：

```text
target,seed,exhaustiveness,rmsd_mode1,rmsd_best_of_9,best_of_9_mode,pass_mode1_lt2,pass_best_of_9_lt2,rmsd_rtm_best,rtm_best_mode
```

### B3. 裁决（Go / No-Go）

| 结果 | 动作 |
|------|------|
| **两端** `rmsd_best_of_9 < 2.0` | **Go** → 进入 Phase C 全面板 |
| 任一端 best_of_9 ≥ 2.0，但 9 个 mode 里有接近失败 | 先查盒子/质子化/配体映射/是否用错构象；**可**在同一 seed 下试 E=16 仅 PM48_01 诊断，**不**据此自动全面板升 E |
| mode1 失败但 best_of_9 成功 | **仍算采样 QC 通过**（与 EGFR/TAK-285 在 3POZ 上相同现象）；全面板照开，但必须在报告里写明，并计划 RTM best-of-9 |
| 两端 mode1 与 best_of_9 都失败 | **No-Go**：停全面板；修蛋白准备/盒子/配体后再重做 Phase B |

写裁决文件：`analysis/cognate_redock_v0/COGNATE_QC_VERDICT.md`（必须含 Go/No-Go 一句结论）。

可选：若 RTMScore 环境可用，对 PM48_01 的 9 个 mode 做 RTM，并填 `rmsd_rtm_best`。

---

## 5. Phase C — 全面板对接（仅 Go 之后）

### C1. 作业矩阵

- 配体：`PM48_01` … `PM48_48`（48）
- 靶：`4L23`, `4JT6`（2）
- **总作业 = 96**
- 每作业：`seed=20260727`, `exhaustiveness=8`, `n_modes=9`, `energy_range=3`

### C2. 输出规范（与 panel40 对齐）

```text
poses/<4L23|4JT6>/<PM48_XX>/mode_01.pdbqt … mode_09.pdbqt
poses/<4L23|4JT6>/<PM48_XX>/<PM48_XX>_all_modes.pdbqt
logs/vina/<target>_<PM48_XX>.log
logs/vina_confs/<target>_<PM48_XX>.txt
```

跑完检查：

- 96 个 log 都有 seed `20260727`
- 96 个配体目录都有 9 个 mode
- `tables/job_status.csv` + `tables/scores_vina.csv` + `tables/scores_vina_long.csv`

### C3. Vina 汇总列（最低集）

`scores_vina.csv` 建议：

```text
panel_id,class,vina_4L23_mode1,vina_4JT6_mode1,vina_mean,vina_min,vina_delta
```

（符号保持 Vina 原生：越负越好；融合前勿随便取负。）

---

## 6. Phase D — RTMScore 重打分 + 决策臂（全面板后）

与 EGFR panel40 相同主臂：

1. 每靶从 prepared 蛋白切 **10 Å pocket** → `receptors/4L23_pocket_10.0.pdb` / `4JT6_pocket_10.0.pdb`
2. 每配体 9 个 mode → SDF，跑 RTMScore
3. 取每端 **max RTM** 的 mode 为 `best_rtm_mode`
4. 双端聚合：`rtm_mean`, `rtm_min`, `rtm_min_z`
5. 主指标（面板内）：
   - AUROC dual vs (A_only ∪ B_only)
   - Top10 中 hardneg 个数 / dual 个数
6. 输出：`tables/ablation_metrics.csv`, `tables/ablation_ranks.csv`, `tables/hardneg_cases.md`

**不要**为打掉某个硬负去改 clash 阈值或事后改面板。

---

## 7. `protocol.yaml` 最小字段（跑前写好）

```yaml
freeze_id: pik3ca_mtor_panel48_v0
engine: AutoDock_Vina
vina_version: "1.2.7"
seed_policy: fixed_global
seed_fixed_global: 20260727
exhaustiveness: 8
n_modes: 9
energy_range: 3
targets:
  PIK3CA:
    pdb: 4L23
    prepared_protein: "D:/CADD paper exercise/dual target docking/Maestro doc/4L23_PIK3CA_prepared.pdb"
    ligand_cocrystal_resname: X6K
    pose_gold_panel_id: PM48_01
  MTOR:
    pdb: 4JT6
    prepared_protein: "D:/CADD paper exercise/dual target docking/Maestro doc/4JT6_mTOR_prepared.pdb"
    ligand_cocrystal_resname: X6K
    pose_gold_panel_id: PM48_01
ligand_prep:
  tool: Schrodinger_LigPrep
  maegz: "D:/CADD paper exercise/dual target docking/Maestro doc/pik3ca_mtor_panel48_v0_ligprep/pik3ca_mtor_panel48_v0_ligprep-out.maegz"
box_definition: "AABB(X6K) + 5A padding; min edge 20A"
cognate_qc_gate: "both ends best_of_9 heavy-atom RMSD < 2.0 A"
```

---

## 8. 执行顺序（检查清单）

- [ ] 复制蛋白到 `receptors/`；记录 Maestro 原路径  
- [ ] 解析 maegz → 48× SDF/PDBQT；`ligand_input_manifest.csv` 齐全  
- [ ] 提取 X6K；写出 `4L23/4JT6` box JSON  
- [ ] **只跑 PM48_01 × 2 端** cognate dock  
- [ ] 算 RMSD；写 `COGNATE_QC_VERDICT.md`  
- [ ] **Go 才**跑 96 作业全面板  
- [ ] 汇总 Vina 表  
- [ ] RTM + ablation（环境可用时）  
- [ ] 打包 freeze 回传仓库（tables/poses 元数据/protocol；大体量 poses 可本地保留并在 MANIFEST 写路径）

---

## 9. 与 EGFR 第一枪的对照（避免做错）

| | EGFR/HER2 panel40 | 本任务 PIK3CA/mTOR panel48 |
|--|-------------------|----------------------------|
| 金标准配体 | EH40_01 TAK-285 / 03P | **PM48_01 PI-103 / X6K** |
| 结构 | 3POZ / 3RCD | **4L23 / 4JT6** |
| seed | as-run 乱 seed → 正向已改固定 | **一开始就用 20260727** |
| exhaustiveness | 敏感性后定为 8 | **直接用 8** |
| 先重对接？ | 已做过 | **必须先做，作为 Go 门** |

---

## 10. 本地 agent 一句话任务

> 使用给定 Maestro 蛋白与 `pik3ca_mtor_panel48_v0_ligprep-out.maegz`，按 Vina 1.2.7 / seed `20260727` / E=8 / n_modes=9，先完成 **PM48_01(PI-103) 对 4L23 与 4JT6 的 cognate 重对接 RMSD 评估**；仅当两端 `best_of_9 < 2 Å` 时再跑完全 48×2 面板，并按 panel40 方式汇总 Vina（及可用时的 RTM）结果。
