# Glide / AF3 Handoff（L7b triage 输出）

**性质**：本包是 chemotype triage 输出，用于放宽阈值的共价对接预筛，**不是活性命中**，待实验验证。

**输入来源**：`prefilter_L0L7/L7b/L7b_dock_ready_*.csv`（Step B 收紧 watch 后）

---

## 文件布局

| 文件 | 用途 |
|------|------|
| `handoff_glide_af3/L7b_{sim_yl,sim_56d,novel,pan}.smi` | Glide/其他对接：`smiles\\tid` |
| `handoff_glide_af3/af3_manifest_{track}.csv` | AF3 分仓清单 |
| `handoff_glide_af3/af3_manifest_ALL.csv` | 合并清单（仅作索引；**禁止按合并总分榜排序选型**） |
| `L7b/L7b_summary.json` | 各仓数量、ok/soft_watch、保留率、ID 交集 |
| `watch_reason_counts.json` | watch 原因审计 |

---

## Glide（共价）

1. **靶点**：JNK2；共价残基 **Cys116**；模板建议 **8ELC**。
2. **阈值**：相对常规 covalent docking **放宽**（本阶段目标是 enrichment / triage，不是精排 IC50）。
3. **分仓分别跑**：`sim_yl` / `sim_56d` / `novel` / `pan` 各自独立作业与榜单。
4. **禁止**：合并四仓做总分榜再截断；pan 仅作对照。
5. **Warhead**：丙烯酰胺；清单中 `warhead_atom_hint` = acrylamide β-碳（末端 =CH2）。

---

## AF3

1. **主门控**：`mPAE`（与此前 AF3 vs 薛定谔对比结论一致）。
2. **勿用** mPAE 细排预测 IC50。
3. **56d 仓**：提高进入 AF3 的比例（相对 sim_yl）；novel 仓按 triage 名单全量或二次抽样均可，但需保留 track 标签。
4. 清单列：`id, smiles, warhead_atom_hint, target=JNK2, cys=116, template=8ELC, track=`。

---

## 配额说明（L7b）

- `sim_yl`：优先 `ok`，不足用 `soft_watch`；目标 ≥3000。
- `sim_56d`：目标 ≥1200。
- `novel`：仅 Step A 启发式 `keep=yes`；`unsure` 另存 `novel_unsure_hold.csv`，未进对接包。
- `pan`：≤300，作对照。

**hard_bad**（氯乙酰胺、乙烯砜、马来酰亚胺等）已剔除；**soft_watch**（芳香胺衍生丙烯酰胺等）可保留但降权。

---

## 明确不做（除非另行授权）

- 本交接包 **不自动启动** Glide 或 AF3 计算。
- 不把任何分子表述为“活性命中”。
