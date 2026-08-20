# Maestro 预准备 → 服务器对接包

> 预检结论见 **`AUDIT_REPORT.md`**（2026-07-27）。脚本已按预检结果修好；集群上主要改 `config.sh` 路径/GPU。

本目录把 **LigPrep 输出 (.maegz)** 与 **9DKB_prepared.pdb** 转成 Vina/gnina/RTMScore 可直接跑的格式，并带 **单分子超时跳过**。


## 已准备好的内容

| 路径 | 说明 |
|------|------|
| `inputs/ligands_all.sdf` | LigPrep maegz → SDF（9839 个构象） |
| `inputs/9DKB_receptor.pdb` | 去掉 UNL 配体与水后的蛋白 |
| `inputs/9DKB_receptor.pdbqt` | Vina 受体 |
| `inputs/lesinurad_crystal_ref.sdf` | 搜索盒 / autobox 参考 |
| `inputs/vina_box.txt` | center/size |
| `inputs/9DKB_pocket_10.0.pdb` | RTMScore 口袋 |
| `work/ligands_sdf/mol_XXXXX.sdf` | 单分子 SDF（gnina），9839 |
| `work/ligands_pdbqt/mol_XXXXX.pdbqt` | 单分子 PDBQT（vina），9838（1 个电荷异常跳过） |
| `work/mol_index_map.csv` | mol_id → 原始 pool 行（1-based Source_File_Index） |

**不再在服务器上从 SMILES 嵌 3D**——直接用 LigPrep 构象。LigPrep 丢掉 10 个输入；打分 join 用 `mol_index_map.csv`，不要假设 `mol_i` = CSV 第 i 行。

## 超时跳过（重要）

在 `config.sh` 中修改（单位：秒）：

```bash
VINA_TIMEOUT_SEC=600    # 默认 10 分钟/分子；超时则跳过并记入 logs
GNINA_TIMEOUT_SEC=1800  # 默认 30 分钟/分子
RTM_TIMEOUT_SEC=3600    # RTMScore 整批
```

实现：GNU `timeout -k 15 <秒>`；超时分子写入 `work/logs/{vina,gnina}/timeouts.txt`，删除不完整输出，继续下一个。

更严示例（服务器调试用）：

```bash
export VINA_TIMEOUT_SEC=300
export GNINA_TIMEOUT_SEC=900
```

## 上机步骤

```bash
cp config.example.sh config.sh
# 编辑 VINA_BIN / GNINA_BIN / RTMSCORE_* / 超时秒数

bash scripts/00_check_env.sh

# 若 work/ligands_* 尚未生成（本机已生成可跳过）:
# python3 scripts/split_ligands_meeko.py --sdf inputs/ligands_all.sdf --outdir work --write-pdbqt --nproc 8

python3 scripts/make_shards.py --ligands-dir work/ligands_sdf --outdir shards --shard-size 50
# 按打印结果改 slurm --array=0-(n-1)

# SLURM
sbatch slurm/vina_array.slurm
sbatch slurm/gnina_array.slurm

# 或无调度器
bash scripts/02_run_vina_local.sh
bash scripts/03_run_gnina_local.sh

# 打分 + 选 Π*
bash scripts/04_rtmscore_rescore.sh vina
bash scripts/04_rtmscore_rescore.sh gnina
bash scripts/05_join_and_select_pi.sh
```

冒烟（片 0）：

```bash
bash scripts/run_vina_shard.sh shards/shard_0000.txt
bash scripts/run_gnina_shard.sh shards/shard_0000.txt
```

## 参数（与协议筛选一致）

- `exhaustiveness=32`，`num_modes=9`，`seed=42`
- 盒子：(99.980, 102.958, 105.657)，≈20 Å（**Vina 与 gnina 使用同一显式盒子**）
- gnina：`--cnn_scoring rescore`；无 GPU 时 `GNINA_NO_GPU=1`；有 cuDNN 后改 `0`
- 已知跳过：`mol_09300` 无 PDBQT（仅 Vina）


## 重新从 Maestro 转换（仅需时）

在 Windows（本机已有 Schrödinger）:

```bat
cd /d "D:\CADD paper exercise\NLRP3_URAT1\URAT1 docking methond test"
D:\Schrodinger2025\utilities\structconvert.exe ^
  unique_docking_pool_ligprep\unique_docking_pool_ligprep-out.maegz ^
  server_dock_maestro_prep\inputs\ligands_all.sdf
```

然后本地再跑 `split_ligands_meeko.py`。
