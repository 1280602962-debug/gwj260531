# URAT1 对接能力矩阵（本机执行）

目的：在现有 GNINA 生产协议下，拆开搜索、姿势选择和受体构象依赖。不重开 `rank_track`，不改 12/21，不跑 MD，不对 156 临床库重对接。

旧 W1 闸门（苯溴马隆 Top-1 ≤ 2 Å 否则停）已取消。9DK9 apo 不进入主矩阵。

## D1（先跑）：9DK 同系列 3×3×3

27 格 = 苯溴马隆 / lesinurad / TD-3 × 9DKA / 9DKB / 9DKC × seeds 42/43/44。

已复用 7 格，本机新跑 **20** 格：

```bash
# 只看作业清单，不调用 gnina
python scripts/run_urat1_capability_dock.py --dry-run

# 有 GPU
python scripts/run_urat1_capability_dock.py --gnina /path/to/gnina --cpu 8

# 仅 CPU
python scripts/run_urat1_capability_dock.py --no-gpu --cpu 8

# 已有 SDF 后只重算表
python scripts/score_urat1_capability_dock.py
```

`gnina` 查找顺序：`--gnina`、环境变量 `GNINA`、`tools/gnina`、PATH。

引擎锁：`exhaustiveness=32`，`num_modes=9`，`cnn_scoring=rescore`。盒子与已准备受体来自 `config/docking_c5_w1.yaml`。输出：

- `data/campaigns/c5/06_capability_dock/d1/d1_job_inventory.csv`
- `data/campaigns/c5/06_capability_dock/d1/d1_capability_summary.csv`
- 逐格 SDF / `.log` / `_stdout.txt`

分类规则（预登记，2.0 Å，`pose_rmsd`，不对配体先叠合）：

| 条件 | `capability_class` |
|---|---|
| Top-1 ≤ 2 Å | `selection_ok` |
| best-of-9 ≤ 2 Å 且 Top-1 > 2 Å | `search_ok_selection_fail` |
| best-of-9 > 2 Å | `search_incomplete`（不得称选择失败） |

交叉对接会先把沉积参考配体经链 A Cα Kabsch 变到该受体坐标系，再算 RMSD。不要用这张表挑“最佳筛库受体”。

## D2（可选，另表）：原生自对接

dotinurad@9B1G、verinurad@9IRY、lesinurad@9B1H，各 3 种子，共 9 个新作业。不要与 D1 拼成一张公平 RMSD 联赛。TD-3@9DKC 已在 D1 对角线，这里不重跑。

```bash
python scripts/run_urat1_capability_d2.py --prepare-only --fetch
python scripts/run_urat1_capability_d2.py --gnina /path/to/gnina
```

dotinurad 使用已准备的中性酚 PDBQT；看过 RMSD 后不要改成酚负离子。

## 不要做

156 库重对、换评分函数、MM/GBSA 补排序、未授权 MD、把 9B1 与 9DK 混成一个 benchmark。
