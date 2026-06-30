# Teacher Gate 质控对接数据集说明

> 对应 `config/docking_ensemble.yaml` → `teacher_gate`  
> 对接方式：**SP → XP**；三态 grid：**9DKB / 9B1K / 9B1L**

---

## 1. 三项质控与数据对应

| Gate | 检验内容 | 通过标准 | 数据集 | 文件 |
|------|----------|----------|--------|------|
| **Gate 1** | lesinurad redock | RMSD ≤ 2.0 Å | **1 条** | `teacher_gate_qc_panel_b_direction.csv` 中 lesinurad |
| **Gate 2** | 四药构象方向 | **4/4** 满足 $S_\pi>0$ | **4 条** | 同上（四药 SMILES） |
| **Gate 3** | 负样本 vs 活性集 | median($\pi_{in}+\pi_{occ}$)_D < median($\pi_{in}+\pi_{occ}$)_A | **822 + 8000 条** | `distill_subset_a.csv` + `distill_subset_d.csv` |

清单总表：`data/distill/teacher_gate_qc_manifest.csv`

---

## 2. 分档执行（推荐顺序）

### 档 A — 网格调试（Gate 1 + 2，必须先过）

| 项目 | 值 |
|------|-----|
| 化合物数 | **5**（lesinurad 用于 redock + 四药方向；lesinurad 只对接一次） |
| 唯一 SMILES | **4**（四药） |
| Grid | 3（9DKB, 9B1K, 9B1L）；Gate 1 redock 仅 **9DKB** |
| 对接任务 | 4 × 3 = **12**（SP→XP） + redock 分析 |
| 输入文件 | `teacher_gate_qc_panel_b_direction.csv` |

**注意**：`verinurad` 在 `distill_manifest.csv` 中**缺失**，但在 `distill_subset_e.csv` 中有。Gate 2 必须手动加入该 SMILES。

### 档 B — 统计 Gate 3（Grid 通过后）

| 集合 | 角色 | 条数 | 文件 |
|------|------|------|------|
| **活性集 A** | 参照（有 pActivity 的 URAT1 训练活性） | **822** | `data/distill/distill_subset_a.csv` |
| **负样本集 D** | 无活性多样性阴性 | **8000** | `data/distill/distill_subset_d.csv` |

对接任务：$(822+8000) \times 3 =$ **26,466**（SP→XP）

可先对 A、D 各 **随机抽 200 条**（seed=42）做预检，再跑全量。

### 档 C — 通过后扩展

全量 `distill_manifest.csv`（8973 条）→ Teacher 标签生成。

---

## 3. 四药 SMILES（Gate 2）

| 药 | SMILES | 共晶 PDB（参考，非对接 grid） |
|----|--------|------------------------------|
| lesinurad | `O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12` | 9DKB |
| benzbromarone | `CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1` | 9DKA |
| verinurad | `CC(C)(Sc1ccncc1-c1ccc(C#N)c2ccccc12)C(=O)O` | 9JDY |
| dotinurad | `O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21` | 9JE1 |

方向检验用**通用三态 grid**（9DKB/9B1K/9B1L），不用药物专属 PDB。

---

## 4. 分数计算（三项 Gate 共用）

对每个分子、每个 grid 取 XP 最终 pose 的 **GlideScore** → $u_{in}, u_{occ}, u_{out}$：

$$
\pi_s = \frac{e^{-u_s}}{\sum_{s'} e^{-u_{s'}}}, \quad s \in \{in, occ, out\}
$$

$$
S_\pi = \pi_{in} + \pi_{occ} - \pi_{out}
$$

| Gate | 判定 |
|------|------|
| Gate 1 | lesinurad @ 9DKB：pose 与共晶 **heavy-atom RMSD ≤ 2.0 Å**（NTD 对齐后） |
| Gate 2 | 四药各自 $S_\pi > 0$（**4/4**） |
| Gate 3 | median($\pi_{in}+\pi_{occ}$) over **D** < median($\pi_{in}+\pi_{occ}$) over **A** |

---

## 5. Schrödinger 输入建议

```text
Phase 0 (Gate 1+2):
  input: teacher_gate_qc_panel_b_direction.csv  (4 SMILES)
  grids: grid_inward_9DKB, grid_occluded_9B1K, grid_outward_9B1L
  mode:  SP → XP each grid
  redock: lesinurad vs 9DKB co-crystal only

Phase 1 (Gate 3):
  input_A: distill_subset_a.csv      (822)
  input_D: distill_subset_d.csv      (8000)
  mode:  SP → XP × 3 grids
  output: per-mol pi_in, pi_occ, pi_out → aggregate medians
```

---

## 6. 失败处理

| 失败项 | 动作 |
|--------|------|
| Gate 1 redock | 检查 9DKB Prep、grid 中心（lesinurad）、质子化；勿换 9JDZ |
| Gate 2 方向 <4/4 | 检查 occ/out grid 是否用 urate 位点；四药是否用 XP 分 |
| Gate 3 中位数倒挂 | 检查 π softmax 是否用三态同批归一化；D 是否误标为有活性 |

**任一 Gate 失败 → 停用 Teacher M-CPDL 标签，回退 v2 手工 $S_{trap}$。**
