# C1 本地战役：可测 URAT1–NLRP3 双靶候选（预注册执行书）

> 你已接受 [`PROJECT_REPLAN_MOLECULAR_DIVERSITY.md`](PROJECT_REPLAN_MOLECULAR_DIVERSITY.md)，目标改为**可实验验证的双靶候选**。  
> 本文是该规划的**本地执行层**：哪些作业必须在有 gnina / GPU / MD 的机器上跑、按什么顺序、什么时候停。切面分析见第 12–16 节。  
> 冻结的 P2 表（`data/repurposing/p2/`，1,580 行）**只作旧战役阴性基线**，C1 不覆盖、不事后抬分。

科学锁：[`config/campaign_c1.yaml`](../config/campaign_c1.yaml)（**不要**传给 `run_gnina_batch.py`）。  
gnina 引擎：[`config/docking_c1.yaml`](../config/docking_c1.yaml)。

---

## 0. 先把两条轨道写死（否则会做成旧 7 个名字的翻案）

对接作为**活性检索器**，和对接作为**羧酸姿态发生器**，过关标准不同。C1 两条轨道都预注册；**Rank 轨不过关就关闭，Acid 轨仍可出假说**。

| 轨道 | 问的问题 | URAT1 证据 | Vecabrutinib | 过关后能否当 Mol Divers 主结果 |
|------|----------|------------|--------------|--------------------------------|
| **Rank** | 羧酸根准备 + 多构象之后，GNINA 能不能在 RandomDecoy 上排出活性？ | CNNscore 选姿后的 CNNaffinity / CNN_VS 百分位 | 仅当 Rank 过关才可能留下 | 可以写 “prioritized testable candidates”，仍须写未做湿实验 |
| **Acid** | 不要求打分器会检索活性；只要求酸根能靠近 Arg477，且 NLRP3 侧有细胞模型或文献 | 羧酸（或等价酸根）+ Arg477 ≤ 4 Å + 口袋内质控 | **默认排除**（不是酸，入选只靠旧 CNN 百分位） | 可以写 “acid-pose hypotheses”，**不能**写成对接排出了活性 |

规划书第 4.7 节把 Vecabrutinib 放在 Tier 1，第 4.5 节又要求酸根–Arg477，这两句不能同时当真。C1 的裁决是：

- **不预承诺**旧 7 个名字进入 C1 短名单；
- Rank 不过关 → Vecabrutinib 不得当 URAT1 候选；
- GSK-3008348 只在 Acid 轨（或 Rank 也过关时）按羧酸规则重新判定，不因为旧 \(S_U=97.5\) 自动入选。

旧漏斗 8319→1588→1580→51→7 整段保留在文中，角色是：**失败排序器迁到临床库会吐出什么**。C1 短名单若与这 7 个重叠，必须是新规则独立通过，不是“还是这几个”。

---

## 1. 为什么这些必须本地

当前云端没有 `gnina` / `vina` / GROMACS。更关键的是科学缺口都在**新的三维作业**，不是再画一张旧表。

| 缺口 | 归档现状 | 不重跑的后果 |
|------|----------|----------------|
| 配体准备 | `prepare_ligands_vina.py` 只 `AddHs` + Meeko，**不在 pH 7.4 枚举羧酸根** | lesinurad 生产姿 Arg477 = 14.20 Å；只改 GSK/lesinurad 是事后修补 |
| `num_modes` | 基准 9，生产 1 | 规划书 4.1 的“基准=生产”不成立 |
| 选姿 vs 排序 | 生产只留 1 个构象；基准 SDF 虽有 9 构象，但是**旧准备** | 不能把旧 SDF 的 CNN_VS 当成羧酸根战役已过关 |
| NLRP3 自对接 | 只有 MCC950 类似物对照 | 7ALV 结构协议没有金标准 |
| 多种子 | 仅 seed 42 | 规划书的稳定性门无法做 |
| MD | 无生产轨迹 | 不能编 RMSD / MM-GBSA |

基准 gnina SDF 里已经有 `CNNscore` / `CNNaffinity` / `CNN_VS` 字段——那只适合在云端做**旧准备**的指标补表，用来和 C1 对比，**不能代替**本地羧酸根重对接。

---

## 2. 机器怎么选

| 阶段 | 分子规模（约） | WSL2 CPU（`--no_gpu`） | 本地 GPU |
|------|----------------|-------------------------|----------|
| L0 环境 + 烟雾 | 几个对照 | 适合 | 适合 |
| L1 羧酸根准备 | ~1 万 SMILES（CPU，RDKit；**不需要 gnina**） | 适合；也可在云端写 PDBQT 再拷到本机 | 适合 |
| L2 自对接 + 强制回收 | ~15 个配体 × 9DKB（+ NP3-146@7ALV） | **必须先在这里过关** | 更快 |
| L3 全诱饵 URAT1 | **9,849** 独特 SMILES（冻结协议完整案例 9,839）× `exh=32` × `num_modes=9` | 不建议作为主路径（量级远大于当年 1,588×2） | **主路径** |
| L3 三种子全库 | L3 ×3 | 不要 | 仅当 seed 42 已过关再加 |
| L4 NLRP3 位点基准 | 视策展规模，通常几百～几千 | 小集合可以 | 大集合用 GPU |
| L5 临床库双靶 | 过关后再跑；池子约 1,588×2，或 Acid 轨下羧酸子集 | GPU 更稳；CPU 可挂机但不要与 L3 抢 | 主路径 |
| L7 MD | 6 体系 × 3 条重复，URAT1 要膜 | 不适合 | 必须 |

结论：**L2 不过关，不要开 L3。L3 用 GPU。不要在 WSL CPU 上把“全诱饵 + 三种子 + 临床库”一次排进去。**

---

## 3. 目录（全部新建，禁止写入 `data/repurposing/p2/`）

```
data/campaigns/c1/
  00_preregistration/     # 本文件副本、campaign_c1.yaml、时间戳
  01_ligand_prep/         # 羧酸根 PDBQT + 准备失败清单
  02_selfdock/            # lesinurad@9DKB，NP3-146@7ALV
  03_forced_recovery/     # 教科书 URAT1 药，不受 q_N 挡住
  04_decoy_dock_9dkb/     # TrueDecoy + RandomDecoy，C1 读出
  05_metrics/             # EF/BEDROC/LogAUC/bootstrap；pass_fail.json
  06_nlrp3_benchmark/     # 位点相关集合（若 L3 已过关或并行策展）
  07_clinical_dock/       # 仅 pass_fail.json 允许后才创建
  08_nomination/          # Acid / Rank 两份短名单
  09_md/                  # 仅短名单冻结后
```

分析结果脚本可以仍在云端跑（读你推上来的 CSV）。**gnina 作业本身必须本地。**

---

## 4. 阶段与停止规则

信息流必须单向：诱饵集与自对接 → `pass_fail.json` → 才允许打开临床库。看过 C1 临床名单之后，**不得**改 L3 过关线、不得把 P0 升为生产读出、不得把 lesinurad 加进阳性后用旧分重算 EF。

### L0 环境（本地，短）

```text
cd URAT1_NLRP3_DualTarget_AIDD_Project
git checkout cursor/urat1-nlrp3-dualtarget-aidd-e43d
git pull

which tools/gnina || bash scripts/setup_gnina_wsl_cpu.sh
tools/gnina --version    # 锁定 1.3.1
nvidia-smi || echo "CPU_ONLY"

python3 -c "import rdkit, meeko, pandas, yaml; print('ok')"
ls data/structures/prepared/9DKB_receptor.pdbqt data/structures/prepared/7ALV_receptor.pdbqt
ls data/benchmarks/urat1_true_decoy/unique_docking_pool.csv
```

受体 PDBQT **沿用**现有 9DKB / 7ALV（搜索盒与冻结战役相同，见 `config/campaign_c1.yaml`）。C1 改的是配体微状态与 `num_modes`，不是换口袋。

### L1 羧酸根配体准备（CPU；不需要 gnina，云端或本地均可）

对以下集合用**同一套** pH 7.4 羧酸根规则（Dimorphite-DL 或等价物，再 Meeko）：

1. `unique_docking_pool.csv`（9,849 SMILES；冻结协议完整案例是 9,839，差的 10 个当时无 gnina 姿态）  
2. 强制回收名单（下节）  
3. 临床库 8,319（可先做 SMILES 枚举，PDBQT 按后续池子再写，避免一次准备全部）

记录：每个分子是否含羧酸、枚举出几个微状态、准备失败原因。  
**禁止：** 只给 GSK 和 lesinurad 手改质子化，诱饵集沿用旧 PDBQT。

**L1 出口：** `01_ligand_prep/prep_summary.json`。lesinurad / GSK-3008348 必须出现**脱质子羧酸根**。若枚举失败，停，先修准备脚本，不要对接。

### L2 自对接 + 强制回收（本地，小，先于全库）

同一 gnina 设置：`exh=32`，`num_modes=9`，`cnn_scoring=rescore`，种子 42/43/44。

| 作业 | 受体 | 通过标准（预注册） |
|------|------|-------------------|
| lesinurad 羧酸根 | 9DKB | CNNscore 选出构象 RMSD ≤ 2.0 Å **且** 羧酸 O–Arg477 ≤ 4.0 Å |
| NP3-146/RM5 | 7ALV | CNNscore 选出构象 RMSD ≤ 2.0 Å（相对共晶） |
| 强制回收表 | 9DKB | 见 L3 的回收条款；L2 先出姿态，确认酸药不再丢盐桥 |

MCC950@7ALV 已有，**仍不得**代替 NP3-146。

**L2 失败怎么走**

- lesinurad 羧酸根后 Arg477 仍 > 4 Å：自由对接的准备漏洞不是主因。下一步是**预注册的 Arg477 约束对接**（仍先只跑自对接+回收，不跑全诱饵池）。约束对接过关后，文章必须改口为“晶体酸根姿匹配”，不能再说“对接排出活性”。
- NP3-146 自对接失败：NLRP3 结构臂降为探索性；双靶候选的 NLRP3 侧只能靠细胞模型 + 口袋占位，主张再降一档。

**L2 不过关，禁止 L3 全诱饵。**

### L3 URAT1 全诱饵（本地 GPU，战役核心）

输入：L1 准备好的独特结构（池文件 9,849 SMILES）。标签用 `data/benchmarks/protocol_selection/mol_protocol_scores.csv` 的 `mol_id` / `role` 对齐，不要另编一套 ID。  
读出（看分之前锁定，见 yaml）：

| ID | 定义 | 角色 |
|----|------|------|
| C1_P2* | 9 构象里 CNNscore 最高者的 CNNaffinity | Rank 轨主排序 |
| C1_VS | 同一构象的 CNN_VS | 并列主排序（规划书未测过的乘积） |
| C1_P2max | 9 构象 CNNaffinity 最大 | 敏感性，不单独锁协议 |
| C1_P0 | 9 构象 CNNscore 最大 | 姿态/合理性比较，**不是**把旧负对照翻案 |

先跑 **seed=42**。三种子全库只在 seed 42 已过关且有 GPU 时加。

**Rank 轨过关（必须同时满足，写进 `05_metrics/pass_fail.json`）：**

1. RandomDecoy 上 C1_P2* **或** C1_VS：EF@1% 点估计 > 1，且前 1% 超几何 *p* < 0.05（切片 \(n_f=\lfloor 0.01N\rfloor\)，与冻结表同一算法）；  
2. 强制回收：lesinurad、verinurad、benzbromarone **至少一个**进入该读出的前 10%；  
3. L2 的 lesinurad 自对接已通过。

TrueDecoy EF、BEDROC(α=80.5 与 160.9)、adjusted LogAUC、PR-AUC、骨架配对 bootstrap **全部要报**，但**不得**用来在 RandomDecoy 失败后改选读出。  
P5/RTM 本轮不作为 C1 生产读出（覆盖与 Random 失败问题未消失，除非你另开预注册）。

**L3 失败：** Rank 轨关闭。允许：(a) 停在方法学阴性 + 旧审计 7；(b) 仅 Acid 轨，用 L2 已通过的羧酸姿规则筛临床酸，**明确不是活性检索**。不允许：改 EF 切片、只报 TrueDecoy、把 C1_P0 升成生产排序、用旧 1,580 分给 lesinurad 补标签。

### L4 NLRP3 位点基准（本地，可与 L3 后半并行）

- 自对接：L2 已做 NP3-146。  
- 另建“直接/位点相关”阳性（磺酰脲口袋配体，与 IL-1 细胞全集分开）。阴性不足则性质匹配诱饵，并在文中写**推定阴性**。  
- 若直接阳性太少：NLRP3 结构验证保持探索性，提名时 NLRP3 臂以 \(q_N\) + 口袋占位为主，不写 “NLRP3 对接已验证”。

### L5 临床库（仅 `pass_fail.json` 允许后）

**Rank 过关时**

- 配体：L1 羧酸根规则重准备的对接池。默认仍从 \(q_N\ge 0.5\) 的 1,588 出发（NLRP3 生物学门），但 **URAT1 强制回收药即使 \(q_N=0\) 也要对接**。  
- 两靶：`num_modes=9`，与 L3 同一 gnina。  
- 排序分量：过关的那个读出（C1_P2* 或 C1_VS，L3 预先指定并列，**取先写在 pass_fail 里通过的主读出**；两个都过则主报 C1_VS，C1_P2* 作敏感性——这条也预先写死，避免看名单后挑）。  
  预注册裁决：**两个都过 → 主读出 = C1_VS**；只有一个过 → 用那一个。

**仅 Acid 轨（Rank 失败）时**

- 从临床库取羧酸（或四唑/酰磺酰胺等价物）∩（\(q_N\ge 0.5\) 或文献 NLRP3/IL-1 证据）；  
- 9DKB 对接只服务 Arg477 / 口袋质控，**不用百分位排名**；  
- 7ALV 要求口袋内、无严重冲突。

### L6 偏倚、指纹、药化（本地出姿态后，分析可回云端）

规划书模块 D/E：MW/电荷相关、配体效率、残差排名、相互作用指纹。  
残差排名在 Rank 轨只作**敏感性**，不能单独把分子抬进短名单（与规划书“不能仅因校正后上升而重定义规则”一致）。

### L7 MD（本地 GPU，短名单冻结之后）

体系数随 **C1 实际短名单** 定，不再预绑 Vecabrutinib：

- 对照永远要：lesinurad 晶体羧酸根@9DKB（膜）；NP3-146 共晶@7ALV（水）。对照失败则该靶候选轨迹一律不解释。  
- 每个 C1 主候选：两靶都跑。  
- 3 条独立重复；有资源每条 100 ns，否则 3×50 ns 并写明探索性质。  
- 只报 RMSD、口袋保留、关键接触占据；**不做 MM-GBSA 重排**，不比谁更亲 lesinurad。  
- MD 不是入选门控，是冻结候选的压力测试。

旧 `MD_RUN_PLAN.md` 里 Vecabrutinib×两靶 **仅当它通过 C1 Rank 轨** 才执行；Acid 轨不要为了“凑双节点”给非酸配体开 URAT1 膜体系。

### L8 湿实验（实验室，非计算）

每个 C1 主候选仍要：**URAT1 摄取抑制** + **NLRP3/IL-1**。计算文章只能写 testable / prioritized，不能写 identified dual inhibitors。

---

## 5. 强制回收表（预注册，不受 \(q_N\) 挡住）

URAT1，9DKB，C1 准备 + C1 读出：

lesinurad，benzbromarone，dotinurad，verinurad，probenecid，puliginurad，SHR-4640。

过关条款只用前三个经典酸药的“至少一个进前 10%”，其余报出但不事后扩成过关条件。

---

## 6. 提名规则（看 C1 临床名单之前锁定）

### Rank 轨（L3 已通过）

1. NLRP3：\(q_N\ge 0.5\) 或文献炎症证据；域外只标记不删。  
2. 两靶 C1 主读出均进入池内前 10%。  
3. 两靶 CNNscore 选姿：无 2.2 Å 冲突、口袋内。  
4. URAT1：若分子含羧酸，Arg477 ≤ 4 Å；**非酸分子不因百分位自动获得与酸药同等的 URAT1 主张**（须在表中单列“非酸、仅排序证据”）。  
5. 偏倚：原始分前 10% **且**（配体效率或残差）前 15%——敏感性；主短名单以 2+3+药化为主。  
6. 药化：无 PAINS/Brenk、Veber、Ro5 氢键/logP、MW 200–550、降级大环。  
7. 稳定性：L2 三种子自对接已过；临床主候选至少再跑 2 个种子或两种羧酸微状态，不跌出前 20%。

主短名单目标：**2 个主候选 + 至多 3 个备选**。不保证与旧 7 个重合。

### Acid 轨（Rank 关闭时的唯一候选路径）

羧酸（或等价物）+ Arg477 ≤ 4 Å + 两靶口袋内 + NLRP3 生物学门 + 同一套药化。  
**没有百分位 AND。** 产品名称必须是 acid-pose dual-node hypotheses。

---

## 7. 规划书第九节补算包：本地 vs 非本地

| # | 条目 | 哪里做 |
|---|------|--------|
| 1 | 基准与生产参数统一 | **本地** L3+L5（同一 yaml） |
| 2 | URAT1 多姿态重算 | **本地** L3；旧 SDF 仅作对照 |
| 3 | NP3-146 自对接 | **本地** L2 |
| 4 | BEDROC/LogAUC/EF/PR-AUC | 旧分可在云端补 SI；**C1 分必须等本地 CSV** |
| 5 | 骨架 bootstrap | 同上 |
| 6–7 | 偏倚 / 残差排名 | 有 C1 分之后云端可算 |
| 8 | 7 个旧名字多种子 | **不要**把旧 7 个当必须作业；改为 C1 短名单的稳定性 |
| 9 | 相互作用指纹 | 本地出 SDF 后，分析可云端 |
| 10 | 按新门控冻结候选 | **本地 L5 之后**；允许空名单 |
| 建议 1 | NLRP3 位点基准 | 名单可云端策展，对接 **本地** L4 |
| 建议 2 | 第二姿态方法 | **本地**，仅 Rank 失败或 L2 失败后的换路线，不要与 L3 同时开四条栈 |
| 建议 3 | MD | **本地** L7 |
| 建议 4 | 可购性/暴露 | 云端/公开数据库即可 |

路线 B/C（KarmaDock、rDock 等）**不要与 L3 并行**。它们是 L3 停损后的换栈，每换一次都要重新预注册，否则又是看结果挑引擎。

---

## 8. 给本地 Agent 的最短任务书（可整段粘贴）

```text
你在本机工作（需要 gnina 1.3.1）。不要覆盖 data/repurposing/p2/。
战役文件：docs/LOCAL_C1_CANDIDATE_CAMPAIGN.md、config/campaign_c1.yaml（科学锁）、
config/docking_c1.yaml（传给 run_gnina_batch.py 的 --config）。

只做 L0→L2，不要开 unique_docking_pool 全诱饵，除非我书面确认 L2 已通过。

L0：安装/确认 tools/gnina、9DKB 与 7ALV 受体 PDBQT、Python 依赖。
L1：实现 pH 7.4 羧酸根枚举（Dimorphite-DL 或等价）+ Meeko；先准备
    lesinurad、GSK-3008348、NP3-146、强制回收名单。
    禁止只改这两个酸而沿用旧诱饵 PDBQT。
    不要把 campaign_c1.yaml 传给 run_gnina_batch.py。
L2：exh=32，num_modes=9，cnn_scoring=rescore，种子 42/43/44。
    lesinurad@9DKB：CNNscore 选姿 RMSD≤2 Å 且 acid–Arg477≤4 Å。
    NP3-146@7ALV：CNNscore 选姿 RMSD≤2 Å。
    读出从 SDF 解析，不要只用 batch CSV 的 max CNNaffinity。
把姿态、距离、RMSD 写成 data/campaigns/c1/02_selfdock/ 下的 CSV/JSON。
失败就停并说明；成功后再问我是否开 L3（GPU 全诱饵）。
```

L3 开跑前再单独授权。命令形态（GPU 示例，路径以本机为准）：

```bash
# 不要使用 run_funnel_p2.sh（旧 num_modes=1）。不要把 campaign_c1.yaml 传给 --config。
python3 scripts/run_gnina_batch.py \
  --config config/docking_c1.yaml \
  --target urat1_9dkb \
  --manifest data/campaigns/c1/01_ligand_prep/unique_pool/ligand_manifest.csv \
  --output-dir data/campaigns/c1/04_decoy_dock_9dkb \
  --jobs 8
# 然后必须从 poses/*.sdf 解析 C1_P2star / C1_VS；batch CSV 只是 C1_P2max。
```

`config/docking_c1.yaml` 已按 `run_gnina_batch.py` 的 schema 写好（`num_modes: 9`，`--seed 42`）。开 L3 前仍须：(1) 羧酸根准备脚本；(2) SDF 多构象解析；(3) 10 分子烟雾确认 SDF 里有 `CNNscore`/`CNNaffinity`。

---

## 9. 云端可以并行、但不算战役过关

这些**不能**替代 L2/L3：

- 从旧基准 SDF 汇总 CNN_VS / CNNscore-选姿的 EF（标签：frozen-prep，不是 C1）；  
- BEDROC/LogAUC 补在冻结 P0–P5 上；  
- 1,580 行 MW–分数偏倚图；  
- NLRP3 位点阳性的 ChEMBL 策展表；  
- 按 C1 规则改 Methods 草稿（Results 数字空着，等 `pass_fail.json`）。

---

## 10. 禁止清单（战役级）

- 覆盖或“更新” `pareto_merged_scores.csv` / 旧 7 个优选表当作 C1 产品；  
- L3 失败后改超几何切片、只展示 TrueDecoy、启用 P0/P5 为生产排序；  
- 把 lesinurad 加进阳性后用**旧**对接分重算 EF；  
- 未写 `pass_fail.json` 就对接临床库并为 Vecabrutinib 开 MD；  
- 用一条 MD 或一次 MM-GBSA 证明双靶活性；  
- 与 L3 同时换 KarmaDock/DiffDock 再看谁能救出旧名单。

---

## 11. 战役结束后文章怎么用这些数字

- **Rank 过关：** 主投 *Molecular Diversity* 可以按规划书第七章写候选发现，但必须把冻结 P2 战役作为“参数不统一 / 中性酸准备”的阴性对照，而不是删掉。C1 短名单是主表。  
- **仅 Acid 过关：** 主结果是羧酸姿态假说 + 方法学失败；期刊更接近 JCAMD，或 Mol Divers 但标题不能写 identifies dual-node candidates from docking ranks。  
- **两轨都失败：** 停。产品就是已有 Fig. 1–4 的方法学阴性文。不要再加引擎。

---

## 12. 现有脚本不能直接当 C1（本地开跑前必须先改）

这是规划里最容易踩空的一层：本机有 gnina 也不等于把 `run_funnel_p2.sh` 再跑一遍就是 C1。

| 现有资产 | 实际行为 | C1 需要 | 谁来改 |
|----------|----------|---------|--------|
| `scripts/prepare_ligands_vina.py` | `AddHs` + Meeko，**不枚举 pH 7.4 羧酸根** | Dimorphite-DL（或等价）→ 脱质子羧酸 → 再 Meeko；记录微状态数与失败原因 | **云端可写脚本并在小分子上测**；全池 PDBQT 云端或本地 CPU 均可 |
| `scripts/run_gnina_batch.py` | 读生产 yaml；`score_mode=cnnaff` 时取 **max CNNaffinity**；不传 `--seed`（除非 extra_args） | 引擎配置已有 [`config/docking_c1.yaml`](../config/docking_c1.yaml)；**读出仍必须从 SDF 解析** | 云端写解析器；本地只跑 gnina |
| `config/campaign_c1.yaml` | 科学锁（过关线、轨道、名单） | **不要**当作 `--config` 传给 batch | 已拆分：科学锁 vs 引擎 yaml |
| `config/docking_production_p2.yaml` | `num_modes: 1` | C1 用 `docking_c1.yaml`（`num_modes: 9`） | 已提供；禁止改生产 yaml |
| `scripts/run_funnel_p2.sh` | 旧漏斗：中性准备 + 单构象 | **禁止用于 C1** | — |
| 冻结协议 P2 | `collect_dock_scores.py` 对 9 构象取 **CNNaffinity 最大**（= C1_P2max） | Rank 主读出是 **CNNscore 最高构象的 CNNaffinity**（C1_P2star）以及该构象的 CNN_VS | 规划书 4.1 与冻结 P2 **不是同一读出**；不能把旧 EF 当作 C1 已过关 |
| `unique_docking_pool.csv` | 只有 SMILES / role，无 `mol_id` | 与 `data/benchmarks/protocol_selection/mol_protocol_scores.csv` 按 SMILES 对齐，沿用 `mol_XXXXX` | 云端写 manifest |
| 7ALV 自对接 | 只有 MCC950 类对照 | 必须从 `7ALV.pdb` 抽 RM5/NP3-146 作参考坐标 | 本地 L2；参考 mol 可在云端先抽 |

**读出对照（看分之前锁定）：**

| ID | 定义 | 冻结战役里有没有 |
|----|------|------------------|
| C1_P2star | 9 构象中 CNNscore 最高者的 CNNaffinity | **没有。** 冻结 P2 = max CNNaffinity |
| C1_VS | 上述同一构象的 `CNNscore × CNNaffinity`（SDF 字段 `CNN_VS`） | 基准 SDF 里有字段，从未当协议 |
| C1_P2max | max CNNaffinity | 就是冻结 P2 的定义（旧准备、旧质子化） |
| C1_P0 | max CNNscore | 冻结负对照 P0 |

因此：即使本地把 `num_modes` 改成 9 再跑现有 batch，得到的仍是 **C1_P2max**，不是规划书要的选姿后再排序。Rank 轨过关必须以 SDF 解析的 C1_P2star / C1_VS 为准。

---

## 13. 三切面：必须本地 / 云端可做但不算过关 / 不要做

判据不是“这台机器方不方便”，而是**缺了 gnina/GPU/GROMACS 是否根本得不到战役数字**。

### A. 必须在有 gnina 的机器上执行

没有 gnina 的云端作业**不能**代替这些。WSL CPU 只适合 L0/L2 烟雾，不适合 L3。

| 作业 | 为什么必须本地 | 停止规则 |
|------|----------------|----------|
| L2 lesinurad@9DKB 自对接 | 羧酸根 + 多构象 + Arg477 距离是新三维实验 | 不过关 → 禁止 L3；可转预注册约束对接（仍本地、仍先自对接） |
| L2 NP3-146@7ALV 自对接 | 仓库里没有该金标准 | 失败 → NLRP3 结构臂降为探索性 |
| L3 全诱饵 URAT1 | ~1 万分子 × exh=32 × 9 构象；要 GPU | Rank 轨的唯一过关实验 |
| L4 NLRP3 位点对接 | 策展表可在云端，姿态必须 gnina | 阳性太少则结构臂保持探索性 |
| L5 临床库双靶 | 必须等 `pass_fail.json`；同一套 gnina | 未过关就对接 = 看名单改规则 |
| L7 MD | 需要 GROMACS（或等价）+ URAT1 膜体系 | 短名单冻结后才开；对照失败则该靶轨迹不解释 |

**算力信封（量级，仓库里没有逐分子计时，按 GPU gnina CNN-rescore 经验）：**

- L2：约 15 配体 × 3 种子，单卡数小时内。
- L3 seed 42：~9,849 配体。若 30–90 秒/分子，约 **80–250 GPU·小时**；三种子再 ×3。不要在 WSL `--no_gpu` 上排。
- L5：Rank 过关后约 1,588×2（加强制回收）。比 L3 小，仍应用 GPU。
- L7：6 体系 × 3 重复 × 50–100 ns 膜蛋白，按本机 GROMACS 吞吐量另排，**不要与 L3 抢同一张卡**。

### B. 云端现在就可以做（不需要 gnina；**不算 C1 过关**）

| 作业 | 产出 | 诚实标签 |
|------|------|----------|
| 羧酸根准备脚本 + 强制回收 7 个药的 PDBQT 烟雾 | `prepare_ligands_c1.py`；lesinurad/GSK 必须出现脱质子羧酸 | 脚本就绪 ≠ L2 通过 |
| 从 SDF 解析 C1_P2star / C1_VS / C1_P2max / C1_P0 | 解析器；可先在冻结基准 SDF 上跑 | 结果必须标 **frozen-prep**，不是羧酸根战役 |
| 冻结 P0–P5 上补 BEDROC(α=80.5/160.9)、LogAUC | SI 表 | 旧读出，不重锁 Π\* |
| 1,580 行 MW/电荷–分数偏倚 | 规划书模块 D 的**旧漏斗**图 | 不能用来给旧 7 个翻案 |
| NLRP3 位点相关阳性策展（ChEMBL 磺酰脲口袋 vs IL-1 细胞全集） | L4 输入表 | 对接仍本地 |
| Methods 按 C1 规则改草稿 | 数字空着 | 等 `pass_fail.json` |
| 从 `7ALV.pdb` 抽出 RM5 参考坐标（RDKit/gemmi） | L2 参考 mol | 自对接仍本地 |

冻结基准 SDF（`docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/poses/gnina_sdf/`，9,839 个）已经带 `CNNscore`/`CNNaffinity`/`CNN_VS`。那只证明**解析器**和“旧准备下 CNN_VS 会不会过 RandomDecoy”，**代替不了**羧酸根重对接。

### C. 不要做（看起来像进度，其实是翻案或并行换栈）

- 覆盖 `data/repurposing/p2/`。
- 只改 GSK/lesinurad 质子化，诱饵沿用旧 PDBQT。
- 把 lesinurad 加进阳性后用**旧**分重算 EF。
- 未写 `pass_fail.json` 就对接临床库 / 为 Vecabrutinib 开 URAT1 膜 MD。
- 与 L3 同时开 KarmaDock / DiffDock / rDock 看谁能救出旧 7 个名字。
- 用 `run_funnel_p2.sh` 或 `campaign_c1.yaml` 当 gnina `--config`。
- 把规划书 4.7 的 Vecabrutinib Tier 1 当成 C1 已提名。

---

## 14. 规划书模块如何落到本地（含 4.5 vs 4.7）

[`PROJECT_REPLAN_MOLECULAR_DIVERSITY.md`](PROJECT_REPLAN_MOLECULAR_DIVERSITY.md) 是产品目标；本节是**执行裁决**，避免模块之间互相打架。

| 规划书 | C1 落地 | 必须本地？ |
|--------|---------|------------|
| 4.1 基准=生产、`num_modes` 9、三种子 | `docking_c1.yaml`；L2 三种子；L3 先 seed 42 | gnina 作业本地；yaml 已在仓库 |
| 4.1 CNNscore 选姿再读 CNNaffinity | **C1_P2star**，从 SDF 解析；不是冻结 P2 | 解析可云端；对接本地 |
| 4.2 URAT1 lesinurad 自对接 + 双诱饵 | L2 + L3 | 本地 |
| 4.2 NLRP3 NP3-146 自对接；MCC950 不能代替 | L2 | 本地 |
| 4.3 BEDROC/LogAUC 等 | 旧分可云端 SI；C1 分等本地 CSV | 指标脚本云端 |
| 4.4 尺寸偏倚 | 有 C1 分之后 | 分析云端 |
| 4.5 酸根–Arg477 为 URAT1 相互作用证据 | **Acid 轨硬门**；Rank 轨对含酸分子同样要求，非酸单列 | 距离计算在本地 SDF 上做，或把 SDF 推回云端 |
| 4.7 Vecabrutinib = Tier 1 | **作废为预承诺。** Rank 不过关不得当 URAT1 候选；Acid 轨默认排除非酸 | 见 §0 |
| 4.7 GSK = Tier 2 | 只按羧酸规则重新判定，不因旧 \(S_U=97.5\) 自动入选 | 本地 L2/L5 |
| 第九节“7 个候选多种子” | 改为 **C1 短名单** 的稳定性，不是旧 7 个必做 | 本地，短名单冻结后 |
| 路线 B/C 换引擎 | L3 停损后才预注册，禁止与 L3 并行 | 若启动，仍本地 |

第十节“两周/两至三周”的日历不要当战役时钟。真正的门是：**L2 通过 → L3 `pass_fail.json` → 才允许 L5**。L2 失败就不要为了赶投稿去筛临床库。

---

## 15. 推荐工作切面（先脚本、后三维、再名单）

单向信息流，和规划书第十节兼容，但按**机器**切开。

```
云端（现在，无 gnina）          本地 GPU（你授权后）
─────────────────────────      ─────────────────────────
羧酸根准备脚本 + 小分子测试  →  L0 锁定 gnina 1.3.1
SDF 解析器（C1_P2star/VS）  →  L1 全池 PDBQT（若云端未写完）
frozen-prep CNN_VS SI（可选） →  L2 自对接 + 强制回收   ← 停在这里等书面确认
NLRP3 位点阳性策展表         →  L3 全诱饵 seed 42
Methods 空数字草稿            →  05_metrics/pass_fail.json
偏倚图（旧 1580 行，SI）     →  仅通过后：L5 临床库
                              →  短名单冻结后：L7 MD
```

**本地 Agent 第一次只接 L0–L2。** L3 是一次独立授权。不要把“规划书第九节必须完成 1–10”理解成本地一次排进全诱饵 + 临床库 + MD。

**Vecabrutinib / GSK：** 写入 C1 短名单之前，它们只是冻结 P2 审计名单上的两个名字。MD 清单见 [`MD_RUN_PLAN.md`](MD_RUN_PLAN.md)——Vecabrutinib 的 URAT1 膜体系仅当 Rank 轨把它留下来才跑。

---

## 16. 本地机器清单（开 L0 前）

```text
□ NVIDIA GPU + 驱动（L3）；nvidia-smi 可见
□ gnina 1.3.1（tools/gnina 或 PATH）；与冻结战役同一 CNN ensemble
□ 受体：data/structures/prepared/9DKB_receptor.pdbqt
         data/structures/prepared/7ALV_receptor.pdbqt
□ Python：rdkit, meeko, pandas, yaml；羧酸根枚举用 dimorphite-dl（或已审查的等价物）
□ 磁盘：L3 SDF ~数千×多构象，预留数十 GB；不要写进 data/repurposing/p2/
□ 分支：cursor/urat1-nlrp3-dualtarget-aidd-e43d，已 pull
□ 配置：config/docking_c1.yaml 作为 --config；config/campaign_c1.yaml 只作科学锁
□ 对照 SMILES：lesinurad / 强制回收表 / NP3-146(RM5) 已能从结构文件或文献结构拿到
```

缺 GPU 时：只做 L0–L2 CPU 烟雾。不要开始 L3。
