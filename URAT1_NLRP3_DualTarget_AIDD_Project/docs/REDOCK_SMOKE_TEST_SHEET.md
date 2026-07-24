# URAT1 重对接烟雾测试表（本机快速版）

> 目的：先确认 **Vina / gnina / RTMScore** 能否在 9DKB 上采样到 near-native，以及各打分读出是否把正确姿排前。  
> 这不是 TrueDecoy 富集考试；**只有 lesinurad@9DKB 是严格自对接**。其余分子是对照/扩展。  
> 填写模板：`data/redock_smoke/redock_results_template.csv`（复制到 `results/redock_smoke/redock_results_filled.csv`）

---

## 0. 本地怎么跑（5–30 分钟量级）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project

# 一次性：受体 + 依赖（若尚未准备）
python3 scripts/prepare_receptor_vina.py --target urat1_9dkb
# gnina（可选）：bash scripts/setup_gnina_wsl_cpu.sh

# 快速烟雾（exhaustiveness=8，多构象）
bash scripts/run_redock_smoke_local.sh

# 正式门控（与 Methods 一致：exhaustiveness=32）
EXHAUST=32 bash scripts/run_redock_smoke_local.sh
```

产出目录：`results/redock_smoke/`  
分数汇总：`results/redock_smoke/scores_summary.csv`  
空白成绩单：复制 `data/redock_smoke/redock_results_template.csv` → `results/redock_smoke/redock_results_filled.csv` 填写。

**本机时间粗估（单分子 lesinurad）：**

| 引擎 | exhaustiveness | 大致耗时（笔记本 4–8 核） |
|------|----------------|--------------------------|
| Vina | 8 | 约 1–5 min |
| Vina | 32 | 约 5–20 min |
| gnina CPU | 8 | 约 5–20 min |
| gnina CPU | 32 | 约 20–60+ min |

先用 `EXHAUST=8` 验证流水线能通；门控结论以 `EXHAUST=32` 为准。

---

## 1. 测试分子（优先顺序）

| 优先级 | 分子 | PDB / 参考 | SMILES（务必用此式） | 角色 |
|--------|------|------------|----------------------|------|
| **必做** | lesinurad | **9DKB / A1AIL** | `O=C(O)CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12` | **唯一严格自对接门控** |
| 建议 | benzbromarone | 9DKB 盒（非共晶） | `CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1` | 强 URAT1 药；看分数能否给出合理姿态 |
| 建议 | dotinurad | 9DKB 盒（非共晶） | `O=C(c1cc(Cl)c(O)c(Cl)c1)N1CS(=O)(=O)c2ccccc21` | 同上 |
| 可选 | MCC950 | **7ALV / RM5** | `CC(C)(O)c1coc(S(=O)(=O)NC(=O)Nc2c3c(cc4c2CCC4)CCC3)c1` | NLRP3 侧烟雾；类似物口袋，非 MCC950 自身共晶 |

> 注意：仓库旧配置里曾出现错误 lesinurad SMILES；以本表与 `data/redock_smoke/redock_pool.csv` 为准。

---

## 2. 对接打分方式 ↔ 要填的格子

每个引擎先产出 **构象集合**（建议 `num_modes=9`），再按读出填表。

| 协议 ID | 构象引擎 | 排序读出 | 方向 | 本机先测？ | 对应脚本 |
|---------|----------|----------|------|------------|----------|
| **P1** | Vina | Vina affinity (kcal/mol) | 越低越好 | ✅ | `run_vina_batch.py` |
| **P2** | gnina | CNNaffinity | 越高越好 | ✅ | `run_gnina_batch.py` |
| **P3** | gnina | gnina affinity (kcal/mol) | 越低越好 | ✅ | 同上 log |
| **P0** | gnina | CNNscore | 越高越好 | ✅ 负对照 | 同上 log |
| **P4** | Vina 构象集 | RTMScore | 越高越好 | ⏳ 装好 RTMScore 再填 | 对 Vina `*_out.pdbqt` 重打分 |
| **P5** | gnina 构象集 | RTMScore | 越高越好 | ⏳ | 对 gnina `*_out.sdf` 重打分 |

**每个协议填 3 个几何指标（lesinurad@9DKB）：**

| 指标 | 含义 | 通过线 |
|------|------|--------|
| **Top-1 RMSD (Å)** | 该读出排名第 1 的姿 vs 晶体 | ≤ **2.0** → 可宣称 pose-accurate |
| **Best-in-ensemble RMSD (Å)** | 9 个姿里相对晶体最低的 RMSD | ≤ **2.0** → 采样够；若 Top-1 失败则是**排序问题** |
| **RTMScore-selected RMSD (Å)** | RTMScore 最高分姿 vs 晶体 | ≤ 2.0 作结构用姿候选 |

判读口诀：

- Best ≤2 且 Top-1 ≤2 → 搜索+打分都可用  
- Best ≤2 但 Top-1 ≫2 → **能对接，但该读出不宜作 pose 证据**；仍可进 TrueDecoy 富集比较  
- Best ≫2 → 搜索盒/质子化/配体准备有问题，先修流水线再比协议  

---

## 3. 成绩单（打印或复制到 CSV）

### 3.1 主表：lesinurad @ 9DKB（必填）

设置：`exhaustiveness=____`（烟雾建议 8；门控填 32）；`num_modes=9`；日期：________

| 协议 | Top-1 分数 | Top-1 RMSD (Å) | Best RMSD (Å) | RTMScore 选姿 RMSD (Å) | 采样 OK? | 打分 OK? | 备注 |
|------|------------|----------------|---------------|------------------------|----------|----------|------|
| P1 Vina affinity | | | | — | ☐ | ☐ | |
| P2 CNNaffinity | | | | — | ☐ | ☐ | |
| P3 gnina affinity | | | | — | ☐ | ☐ | |
| P0 CNNscore（负对照） | | | | — | ☐ | ☐ | |
| P4 Vina→RTMScore | | | | | ☐ | ☐ | |
| P5 gnina→RTMScore | | | | | ☐ | ☐ | |

### 3.2 扩展表：同盒对照（可选，无晶体 RMSD 则只记分数/目视）

| 分子 | 引擎 | 读出 | Top-1 分数 | 目视口袋？ | 备注 |
|------|------|------|------------|------------|------|
| benzbromarone | Vina / gnina | affinity / CNNaff | | ☐ | 非自对接 |
| dotinurad | Vina / gnina | affinity / CNNaff | | ☐ | 非自对接 |
| MCC950 @ 7ALV | gnina | affinity / CNNaff | | ☐ | 类似物模板 |

---

## 4. RMSD 怎么快速量（本机）

1. 从 9DKB 导出共晶配体 **A1AIL**（PyMOL / ChimeraX / gemmi）。  
2. 打开对接输出：Vina `results/redock_smoke/vina/poses/LESINURAD_out.pdbqt` 或 gnina `.../gnina/poses/LESINURAD_out.sdf`。  
3. 蛋白叠合后，对配体重原子算 RMSD（勿用全原子含 H）。  
4. 对 mode 1…9 各算一次 → 填 Top-1 与 Best。  

若暂时只验证“能不能跑通”，可先只填 **分数列**，几何列用 PyMOL 补。

---

## 5. 与后续 TrueDecoy 的关系

| 阶段 | 问题 | 本表回答 |
|------|------|----------|
| 重对接烟雾 | 引擎能否采到近晶体姿？读出是否排对？ | ✅ 本表 |
| 协议筛选 | 哪种读出在 TrueDecoy 上 EF@1% 最好？ | 下一阶段全库 |
| 生产 Π* | 临床库怎么排序？ | 锁定后再跑 1588 |

**即使 Top-1 RMSD 未过 2 Å，只要 Best-in-ensemble 过关，仍可继续做富集比较**；正文需诚实报告三类 RMSD，并限制 pose 主张。

---

## 6. 相关文件

| 文件 | 用途 |
|------|------|
| `data/redock_smoke/redock_pool.csv` | 正确 SMILES 测试池 |
| `data/redock_smoke/redock_results_template.csv` | 可填成绩单 |
| `scripts/run_redock_smoke_local.sh` | 本机一键 Vina(+gnina) |
| `results/redock_smoke/` | 本地运行产出（gitignore） |
| `config/docking_open_source.yaml` | 搜索盒与默认参数 |
| `docs/GNINA_BENCHMARK_REDOCK_WSL.md` | 你本机已有 gnina 目录时的备用流程 |
