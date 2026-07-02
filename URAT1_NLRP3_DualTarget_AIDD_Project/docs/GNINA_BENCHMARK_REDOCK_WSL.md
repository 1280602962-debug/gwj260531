# GNINA Benchmark 对接测试（WSL · 你的本地 gnina 环境）

适用于已安装在：

```text
/mnt/d/CADD paper exercise/gnina/
```

论文项目 GitHub 仓库：

```text
https://github.com/1280602962-debug/gwj260531
```

子目录：`URAT1_NLRP3_DualTarget_AIDD_Project/`

---

## 一、复制给 WSL Agent 的完整指令

```text
你在 WSL 中工作。GNINA 已安装在：
  GNINA_ROOT="/mnt/d/CADD paper exercise/gnina"

论文项目从 GitHub 获取：
  REPO_URL="https://github.com/1280602962-debug/gwj260531.git"
  WORK_DIR="${WORK_DIR:-$HOME/gwj260531}"   # 或 /mnt/d/CADD paper exercise/gwj260531
  git clone --depth 1 "$REPO_URL" "$WORK_DIR" 2>/dev/null || (cd "$WORK_DIR" && git pull)
  PROJECT="$WORK_DIR/URAT1_NLRP3_DualTarget_AIDD_Project"

任务：对 9DKB（URAT1）和 7ALV（NLRP3）做 **benchmark 单分子 GNINA 对接测试**（先不要跑 1588 全库）。每次对接前必须：
  source "$GNINA_ROOT/activate.sh"

【结构文件】
- 9DKB：$PROJECT/data/structures/pdb/9DKB.cif（若无则从 RCSB 下载）
- 7ALV：$PROJECT/data/structures/pdb/7ALV.pdb

若 prepare_docking.py 需要 PDB 而非 CIF，请用 Open Babel 或 gemmi 将 9DKB.cif 转为 9DKB.pdb：
  obabel $PROJECT/data/structures/pdb/9DKB.cif -O $GNINA_ROOT/output/benchmark/9DKB.pdb

【参考配体定义搜索盒 — autobox】
- 9DKB：共晶配体残基名 A1AIL（lesinurad）
- 7ALV：共晶配体残基名 RM5（MCC950-class analog）

【Benchmark 分子 — 先用你方 prepare_docking.py + run_docking.sh】

靶点 9DKB（URAT1）对接这 4 个：
1. lesinurad   SMILES: O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12
2. benzbromarone SMILES: CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1
3. dotinurad   SMILES: O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21
4. EGCG        SMILES: O=C(O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1

靶点 7ALV（NLRP3）对接这 2 个：
5. MCC950      SMILES: CC(C)(O)c1coc(S(=O)(=O)NC(=O)Nc2c3c(cc4c2CCC4)CCC3)c1
6. EGCG        （同上 SMILES）

【推荐目录结构】
$GNINA_ROOT/output/benchmark/
  9dkb/prepare/          # prepare_docking.py 输出
  9dkb/dock_lesinurad/
  9dkb/dock_benzbromarone/
  ...
  7alv/prepare/
  7alv/dock_mcc950/
  ...

【每个 benchmark 的执行模板 — 9DKB 以 lesinurad 为例】

mkdir -p "$GNINA_ROOT/output/benchmark/9dkb/dock_lesinurad"
source "$GNINA_ROOT/activate.sh"

# Step A: 准备受体 + 参考配体（从 9DKB 复合物）
python "$GNINA_ROOT/scripts/prepare_docking.py" \
  --receptor "$GNINA_ROOT/output/benchmark/9DKB.pdb" \
  --ligand-resname A1AIL \
  --ligand-smiles "O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12" \
  --out-dir "$GNINA_ROOT/output/benchmark/9dkb/prepare_lesinurad"

# Step B: 对接（Glide-XP–like：exhaustiveness=32，只保留 1 个最佳 pose，CNN rescore）
gnina -r "$GNINA_ROOT/output/benchmark/9dkb/prepare_lesinurad/receptor.pdb" \
  -l "$GNINA_ROOT/output/benchmark/9dkb/prepare_lesinurad/query_ligand.sdf" \
  --autobox_ligand "$GNINA_ROOT/output/benchmark/9dkb/prepare_lesinurad/ref_ligand.sdf" \
  -o "$GNINA_ROOT/output/benchmark/9dkb/dock_lesinurad/docked.sdf" \
  --exhaustiveness 32 --num_modes 1 --cnn_scoring rescore --no_gpu \
  --log "$GNINA_ROOT/output/benchmark/9dkb/dock_lesinurad/dock.log"

7ALV 模板相同，换：
  --receptor 7ALV.pdb
  --ligand-resname RM5
  MCC950 / EGCG 的 SMILES

【汇总表 — 必须输出 CSV】
路径：$PROJECT/results/gnina_benchmark/benchmark_redock_summary.csv

列：target,pdb_id,compound,affinity_kcal_mol,cnn_pose_score,cnn_affinity,pose_file,log_file,status,notes

【成功标准】
- 6 个对接均生成 docked.sdf 且 mode 1 亲和力为负值
- lesinurad @ 9DKB：最好能与晶体 pose 比 RMSD（若 prepare 脚本保留了 ref）
- 汇总表 + 每个 log 前 30 行贴到汇报

【不要执行】
- 1588 全库批量

【若 prepare_docking.py 参数与上文不一致】
先运行 python .../prepare_docking.py --help，按实际参数改，但保持「autobox_ligand = 共晶配体」原则不变。
```

---

## 二、更短版 Agent 指令（推荐 · 一键脚本 · Glide-XP–like）

```text
你在 WSL 中工作。建议在 tmux 中运行（CPU 约 1.5–3 小时）。

GNINA_ROOT="/mnt/d/CADD paper exercise/gnina"
REPO_URL="https://github.com/1280602962-debug/gwj260531.git"
WORK_DIR="$HOME/gwj260531"
git clone --depth 1 "$REPO_URL" "$WORK_DIR" 2>/dev/null || (cd "$WORK_DIR" && git pull)
PROJECT="$WORK_DIR/URAT1_NLRP3_DualTarget_AIDD_Project"

任务：对 9DKB（URAT1）和 7ALV（NLRP3）做 6 个 benchmark GNINA 对接测试，不要跑 1588 全库。

对接列表：
- 9DKB：lesinurad, benzbromarone, dotinurad, EGCG（autobox 残基 A1AIL）
- 7ALV：MCC950, EGCG（autobox 残基 RM5）

参数（Glide-XP–like，只保留最佳构象）：
  exhaustiveness=32
  num_modes=1
  cnn_scoring=rescore
  --no_gpu（无 NVIDIA GPU 时）

执行：
  export GNINA_ROOT
  export PROJECT_ROOT="$PROJECT"
  export EXHAUST=32 NUM_MODES=1 CNN_SCORING=rescore
  cd "$PROJECT"
  bash scripts/benchmark_gnina_redock.sh

说明：分数不与历史 Glide XP 混用；log 中 mode 1 的 affinity (kcal/mol) 为 dock_score。

完成后汇报：
- $PROJECT/results/gnina_benchmark/benchmark_redock_summary.csv
- 每个 dock.log 前 30 行
- lesinurad @ 9DKB：若可行，报告重对接 RMSD（相对晶体 pose）
```

---

## 三、Benchmark 一览

| # | 靶点 | PDB | 化合物 | autobox 残基 | 论文角色 |
|---|------|-----|--------|--------------|----------|
| 1 | URAT1 | 9DKB | lesinurad | A1AIL | 重对接 / 阳性 |
| 2 | URAT1 | 9DKB | benzbromarone | A1AIL | MD 基准 |
| 3 | URAT1 | 9DKB | dotinurad | A1AIL | MD 基准 |
| 4 | URAT1 | 9DKB | EGCG | A1AIL | Pareto lead |
| 5 | NLRP3 | 7ALV | MCC950 | RM5 | MD 基准 |
| 6 | NLRP3 | 7ALV | EGCG | RM5 | 探索性 MD |

可选加测（非必须）：verinurad @ 9DKB；colchicine @ 7ALV（预期 NLRP3 对接一般）。

---

## 四、与论文项目脚本的衔接

Benchmark 通过后，全库对接可二选一：

1. **继续用你的** `run_docking.sh` 批处理（推荐与本地环境一致）  
2. **或用项目内** `scripts/run_gnina_batch.py`（需把 `GNINA_ROOT/bin/gnina` 写入 `config/docking_open_source.yaml` 的 `gnina.binary`）

将 benchmark 汇总 CSV 复制到：

```text
$PROJECT/results/gnina_benchmark/
```

供 `merge_docking_pareto.py` 调试前验证分数列格式。

---

## 五、预期耗时（CPU --no_gpu，exhaustiveness=32）

单分子约 15–30 分钟 × 6 ≈ **1.5–3 小时**。建议 `tmux` 挂机。

---

## 六、GNINA 参数与 Glide XP 对照

| 目的 | GNINA 设置 | 说明 |
|------|------------|------|
| 搜索深度（≈ Glide SP 量级） | `--exhaustiveness 32` | 非 XP 经验势，但采样更充分 |
| 只保留最佳 pose | `--num_modes 1` | 虚筛 / Pareto 用 mode 1 affinity |
| CNN 精修打分 | `--cnn_scoring rescore` | Vina 采样 + CNN 重打分 |
| CPU | `--no_gpu` | 有 GPU 时可去掉以加速 |

**勿将** GNINA affinity 与历史 `r_glide_XP_GScore` 混在同一分析中。

---

## 七、常见问题

| 问题 | 处理 |
|------|------|
| 9DKB 是 CIF | `obabel -icif -opdb` 转 PDB |
| libcudnn | 已用 activate.sh 解决 |
| prepare 脚本参数名不同 | 先 `--help`，保持 autobox 逻辑 |
| EGCG 失败 | 检查 SMILES 手性；尝试简化质子化 |
