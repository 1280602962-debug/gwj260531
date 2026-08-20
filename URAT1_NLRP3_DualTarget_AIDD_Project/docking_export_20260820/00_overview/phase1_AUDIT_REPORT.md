# 对接包预检报告（2026-07-27）

## 结论

**受体 / 配体格式 / 分片 / Vina+gnina 冒烟均通过。**  
已修复 2 个会影响公平对比或后续跑崩的脚本问题（见下）。上传集群后只需改 `config.sh` 里的**软件路径**（及 GPU 开关），不必再改对接逻辑。

---

## 已修复（本次检查中改掉）

| 问题 | 原状 | 现况 |
|------|------|------|
| gnina 搜索盒 | `--autobox_ligand` → 约 15–18 Å，**小于** Vina 的 20 Å | 与 Vina **同一** `center/size` |
| RTMScore | 一次合并全部 pose；`\| tee` 吞掉超时退出码 | **按 500 分子分块**；`pipefail`；写出 `scores.csv` |

---

## 检查通过项

### 受体
- `9DKB_receptor.pdb`：7883 atoms，**无 UNL、无水**，含 H（Maestro 准备）
- `9DKB_receptor.pdbqt`：刚性受体，4688 ATOM，无 BRANCH/TORSDOF
- 与 redock 烟雾试验原子数一致（7883）

### 搜索盒
- center `(99.980, 102.958, 105.657)`，size `20 / 20.01 / 20`
- 晶体配体 COM 与盒子中心距离 **&lt; 0.001 Å**，晶体原子 100% 在盒内
- RTM 口袋 109 残基全部落在受体上；口袋 COM 距盒子中心 ~1.1 Å（正常）

### 配体
- SDF **9839**，PDBQT **9838**（`mol_09300` Meeko 电荷失败，仅跳过 Vina）
- 抽样 44 个：构象 / ROOT+TORSDOF 均正常
- `mol_index_map.csv`：9839 行、无重复；`source_file_index` 与 pool SMILES **200/200 一致**
- 角色：active 469 + true_decoy 4687 + random_decoy 4683 = 9839
- LigPrep 坐标在原点附近属正常；对接引擎会放入盒子（log 中 `initial pose not within box` 可忽略）

### 分片
- 197 shards，覆盖 `[0,9839)`，无空隙；SLURM `--array=0-196` 正确

### 冒烟
- **Vina** `mol_00010`：9 modes，top affinity −7.827 kcal/mol
- **gnina** `mol_00000`（CPU `--no_gpu`）：9 poses，含 `minimizedAffinity` / `CNNscore` / `CNNaffinity`

---

## 已知可接受项（不必改文件）

1. **`mol_09300`**：无 PDBQT → Vina 自动 skip；gnina 仍有 SDF 可跑  
2. LigPrep 丢弃 10 个输入分子（原 9849→9839）  
3. gnina 输出 SDF 被 RDKit 报 “tagged as 2D”：属性仍完整，不影响打分采集  
4. `inputs/ligands_all.sdf`、`lesinurad_crystal_UNL.pdb`：对接脚本不读，可不上传

---

## 集群上只需改这些（`config.sh`）

```bash
export VINA_BIN=vina                    # 或绝对路径
export GNINA_BIN=/path/to/gnina
export GNINA_NO_GPU=0                   # 装好 cuDNN + module load cuda 后
export OBABEL_BIN=obabel                # P4 需要
export RTMSCORE_PY=/path/to/rtmscore.py
export RTMSCORE_MODEL=/path/to/rtmscore_model1.pth
export RTMSCORE_ENV_ACTIVATE='source ... && conda activate rtmscore'  # 按实际改

# 超时（一般保持默认即可）
# VINA_TIMEOUT_SEC=600
# GNINA_TIMEOUT_SEC=1800
```

cuDNN 安装后：

```bash
module load nvidia/cuda/12.2   # 按站点模块名
export LD_LIBRARY_PATH=$(python -c "import os,nvidia.cudnn; print(os.path.join(os.path.dirname(nvidia.cudnn.__file__),'lib'))"):$LD_LIBRARY_PATH
```

---

## 上传内容（最小完整集）

```text
server_dock_maestro_prep/
  config.sh                 # 从 config.example.sh 复制后改路径
  scripts/  slurm/  shards/
  inputs/   # 受体 pdb+pdbqt、口袋、benchmark csv、crystal_ref（不必传 ligands_all.sdf）
  work/ligands_sdf/
  work/ligands_pdbqt/
  work/mol_index_map.csv
```
