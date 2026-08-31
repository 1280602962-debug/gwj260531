# COGNATE QC VERDICT — PM48_01 (PI-103 / X6K)

**结论：No-Go（按协议 E=8 硬门槛）**

PIK3CA（4L23）在 E=8 已过关；mTOR（4JT6）在协议 exhaustiveness=8 下 `best_of_9 ≥ 2.0 Å`。  
**禁止开全面板 48×2。**

同 seed 下 E=16 诊断显示 4JT6 可以采到近晶构象（mode3，0.445 Å）→ 更像 **mTOR 端 E=8 采样不足**，不像盒子坐标系错误。按 SOP：**不据此自动全面板升 E**；需你确认下一步（修协议升 E / 再诊断 / 或接受风险后开面板）。

---

## 协议参数（主 QC）

- engine: AutoDock Vina 1.2.7
- seed: `20260727`
- exhaustiveness: **8**
- n_modes: 9 / energy_range: 3
- box: AABB(X6K)+5 Å，min edge 20 Å
- RMSD: 重原子；meeko `REMARK SMILES IDX` 映射；模板自同构 min CalcRMS；**同坐标系不叠合**

---

## 主结果（E=8）

| target | rmsd_mode1 | rmsd_best_of_9 | best_mode | mode1&lt;2 | best9&lt;2 |
|--------|------------|----------------|-----------|-----------|-----------|
| 4L23 | **0.624** | **0.624** | 1 | True | **True** |
| 4JT6 | 7.118 | **5.003** | 2 | False | **False** |

表：`analysis/cognate_redock_v0/tables/pm48_01_rmsd.csv`

---

## 诊断（仅 PM48_01，E=16，同 seed）

| target | rmsd_mode1 | rmsd_best_of_9 | best_mode | best9&lt;2 |
|--------|------------|----------------|-----------|-----------|
| 4L23 | 0.624 | 0.624 | 1 | True |
| 4JT6 | 7.118 | **0.445** | 3 | **True** |

表：`analysis/cognate_redock_v0/tables/pm48_01_rmsd_E16_diag.csv`

判读：
- 4L23：E=8 已 mode1 近晶，协议友好。
- 4JT6：E=8 未找到近晶；E=16 找到（best_of_9 通过，mode1 仍失败）——与 EGFR/TAK-285「采样到了但打分未排第一」同类，但本任务 **Go 门要求的是 E=8 两端 best_of_9**。

---

## 建议下一步（需你选）

1. **保持 No-Go**：检查 4JT6 质子化/受体准备/盒子 padding，再重跑 Phase B（E=8）。
2. **协议例外（需书面确认）**：将本项目 `exhaustiveness` 升为 16 后重做 Phase B；通过后再 Phase C（**不要**在未确认时自动全面板）。
3. 不做全面板，先只把 Phase A/B 产物归档。

---

## Phase A 完成摘要

- 48/48 配体 SDF+PDBQT（maegz 映射 `s_canvas_panel_id`，每 ID 1 conf）
- 受体：`receptors/4L23_receptor.pdbqt`、`4JT6_receptor.pdbqt`
- 盒子：`boxes/all_boxes.json`（X6K）
- 协议：`protocol/protocol.yaml`（seed=20260727，E=8）
