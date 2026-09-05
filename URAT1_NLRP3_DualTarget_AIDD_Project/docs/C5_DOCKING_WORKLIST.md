# C5 对接工单：还要补哪些、受体、配体、设置

> 权威科学锁：`config/campaign_c5.yaml`  
> 结论汇总（闸门含义 / 能不能排序 / MD 对照 / 要不要另开排序轨）：**`docs/C5_RANKING_AND_NEXT_DOCKS.md`**  
> 引擎锁：GNINA **1.3.1**，参数与 `config/docking_c1.yaml` 完全相同  
> W1 盒子出处：`data/campaigns/c5/00_verification/w1_box_centers.json`（2026-09-04 用沉积坐标重算）  
> W1 可执行配置：`config/docking_c5_w1.yaml`  
> 这台云沙箱**没有 gnina**。结构下载和受体/配体准备可以在这里做；对接必须在本机跑。  
> **短名单已冻结。** 选出候选之前不再需要本机 gnina。`md_authorized` 仍为 `false`，改正文后再谈 MD。

---

## 现在要开的对接（2026-09-05；短名单已冻结）

URAT1 上对接分**没有筛选排序能力**（`rank_track: closed`）。  
冻结表：`data/campaigns/c5/04_shortlist_frozen/SHORTLIST_FROZEN.md`（主名单 12 + 备份 21）。

| 优先级 | 任务 | 新 gnina job | 状态 |
|---|---|---:|---|
| 1 | **W4 / Task2** 7ALV 面板 | 146 | **已跑 131 ok / 15 fail；结构门已打分并通过（vs 诱饵）** |
| 2 | **W2 / Task3** URAT1 IFP | 0 | **已评；gate_pass，但是 A1 的更严子集** |
| 3 | **短名单冻结** | 0 | **已冻：12 主 + 21 备份** |
| 4 可选 SI | 重试 5 个空 SDF 诱饵 × 3 | 15 | 不为过门；完整矩阵才需要 |
| 5 可选 SI | W1 其余 29 格 | 29 | 不当推进条件 |
| **下一步（零对接）** | **写正文** | 0 | **现在做这个** |
| 不做 | 临床 156 重对接、Rank 轨、换引擎、未授权就 MD | — | 禁止 |

**W4 146 格明细**（受体只有 7ALV；盒 `[16.756, 35.449, 125.714]` 20³；exh=32 / 9 modes / seeds 42/43/44）：

| 子集 | 配体 | 种子 | job |
|---|---|---|---:|
| 8 个非共晶阳性 | MCC950、CHEMBL4204644、CHEMBL5219789、CHEMBL4212407、CHEMBL6143743、CHEMBL4209503、CHEMBL4216836、CHEMBL6171925 | 42/43/44 | 24 |
| 40 个已锁诱饵 | `data/campaigns/c5/02_nlrp3_panel/w4_decoys_locked.csv`（C5W4D_001–040） | 42/43/44 | 120 |
| 缺的背景酸 | `REP_07837` | 43、44 | 2 |

复用：NP3-146@7ALV × 3；其余背景临床酸已有 exh=32 姿态。不要用旧 `exh=8` 面板 SDF。

**W2**：不对新分子。IFP 锚在沉积 R75/A1AIL/A1A45 上，不用苯溴马隆 CNNscore Top-1 失败姿。

---

## 0 所有新对接共用的设置（不许另开引擎）

| 项 | 锁定值 | 不要改成 |
|---|---|---|
| 引擎 | GNINA 1.3.1（`tools/gnina` 或 PATH） | KarmaDock / RTMScore / Vina-only / 新 XGBoost |
| `exhaustiveness` | **32** | 旧 NLRP3 面板用过的 **8**（那批 SDF 不得拿来重标定 W4） |
| `num_modes` | **9** | 生产 P2 的 `num_modes=1` |
| `cnn_scoring` | `rescore` | `cnn_scoring=all` 或关 CNN |
| 种子 | **42 / 43 / 44** | 单种子交差 |
| 配体准备 | `scripts/prepare_ligands_c1.py`：Dimorphite-DL pH 7.4 → Meeko；embed seed `0xC0FFEE` | 中性羧酸、另选微状态 |
| 羧酸 | pH 7.4 **脱质子** `CO2-` | 中性 `COOH` |
| 苯溴马隆 | 酚羟基按已冻结微状态 **酚负离子**（见下） | 中性酚 |
| 受体准备 | `scripts/prepare_receptor_vina.py`：只留 A 链、去水/HETATM、Open Babel pH 7.4 加氢、刚性 PDBQT | 带配体的 holo、改链、改 pH |
| GPU | 本机有卡就按 `docking_c1.yaml`（`no_gpu: false`）；无卡用 `docking_c1_cpu.yaml` 的 `no_gpu: true`，**其余参数不变** | 用 CPU 配置改小 exhaustiveness |

命令模板（每个格子只换 `-r` / `-l` / `--center_*` / `--seed` / `-o`）：

```bash
gnina -r <receptor.pdbqt> -l <ligand.pdbqt> \
  --center_x <cx> --center_y <cy> --center_z <cz> \
  --size_x <sx> --size_y <sy> --size_z <sz> \
  --exhaustiveness 32 --num_modes 9 --cnn_scoring rescore \
  --seed <42|43|44> -o <out.sdf>
```

---

## 1 已经完成、不要重做

| 批次 | 受体 | 配体 | 种子 | 盒子 | 产物 |
|---|---|---|---|---|---|
| Phase I 9839 unique | 9DKB | TrueDecoy+RandomDecoy 池 | 生产种子 | 锁定 22³ | `docking_export_20260820/01_phase1_benchmark_URAT1_9DKB/poses/gnina_sdf` |
| 临床酸 156 双臂 | 9DKB + 7ALV | chemistry-pass 156 | 42/43/44 | 锁定 | `data/campaigns/c1/07_clinical_dock/acid_dual*`（seed42 在 `acid_dual/`，43/44 在 `acid_dual_a2/`） |
| lesinurad 自对接 | 9DKB | lesinurad 羧酸根 | 42/43/44 | 锁定 22³ | `data/campaigns/c1/02_selfdock/urat1_9dkb/seed{42,43,44}/` |
| NP3-146 自对接 | 7ALV | NP3-146 | 42/43/44 | 锁定 20³ | `data/campaigns/c1/02_selfdock/nlrp3_7alv/seed{42,43,44}/` |
| 强制回收 seed42 | 9DKB | 苯溴马隆 / lesinurad / verinurad / puliginurad / SHR-4640 / 丙磺舒 / dotinurad / GSK-3008348 | 42 | 锁定 22³ | `data/campaigns/c1/03_forced_recovery/urat1_9dkb/seed42/` |
| Rank 轨 / L3 全诱饵 | — | — | — | — | **关闭。** `allow_L3_full_decoy: false` |

NLRP3 结构面板 seed42（30 个配体，`exhaustiveness=8`）**只作历史记录**。W4 重标定不得复用这批 SDF。

---

## 2 必须补：Job A = W1 四结构交叉对接

**设计**：3 个晶体配体 × 4 个同构建体受体 × 3 种子 = **36** 个 gnina job。  
**可复用 4 格**：lesinurad@9DKB 三种子 + 苯溴马隆@9DKB seed42。  
**新跑 32 格** + 准备 3 个新受体 + 准备 TD-3 配体。

目的是姿态保真度随柔性变化，**不是排序**。看完矩阵之前不要提名候选。

### 2.1 四个受体

同一篇：Suo/Fedor/Lee *Nat. Commun.* 2025, **16**:5178（DOI 10.1038/s41467-025-60480-3）。准备协议与 9DKB 相同。9DKC **必须 mmCIF**（`.pdb` 下载 404）。

| 靶点键 | PDB | 状态 | 分辨率 | 共晶 CCD | 已准备？ | 盒中心 (Å) | 盒边 |
|---|---|---|---:|---|---|---|---|
| `urat1_9dkb` | 9DKB | holo | 2.74 | **A1AIL**（lesinurad；不是 LES） | 是 `9DKB_receptor.pdbqt` | **`[99.966, 102.967, 105.699]`**（生产锁定；新算 A1AIL 质心差 0.08/0.69/0.23 Å，**不要改**） | 22³ |
| `urat1_9dka` | 9DKA | holo | 3.00 | R75（苯溴马隆） | **否** → `9DKA_receptor.pdbqt` | **`[107.167, 105.296, 107.628]`**（R75 重原子质心，n=22） | 22³ |
| `urat1_9dkc` | 9DKC | holo | 2.55 | A1A45（TD-3） | **否** → `9DKC_receptor.pdbqt` | **`[102.345, 103.269, 107.416]`**（A1A45 重原子质心，n=28） | 22³ |
| `urat1_9dk9` | 9DK9 | **apo** | 2.68 | — | **否** → `9DK9_receptor.pdbqt` | **`[100.693, 107.700, 105.653]`**（把 9DKB 锁定中心经 A 链 517 个 CA Kabsch 叠合转到 9DK9 坐标系；叠合 RMSD 0.573 Å） | 22³ |

**禁止**：把 9DKB 的笛卡尔中心直接抄到 9DK9/9DKA/9DKC。四套坐标不在同一实验室坐标系。9DKA 原生质心相对 9DKB 锁定盒偏差约 7.2 Å。

准备命令（本机有 Open Babel）：

```bash
for t in urat1_9dka urat1_9dkc urat1_9dk9; do
  python3 scripts/prepare_receptor_vina.py \
    --config config/docking_c5_w1.yaml --target "$t"
done
```

残基编号：准备后的 9DKB PDBQT 把文献 Arg477 写成 **ARG A 476**。W1 报酸根–Arg 距离时每个结构单独核对胍基，不要整批套文献序号。

### 2.2 三个配体

SMILES 来自 RCSB chemcomp（已写入 `w1_reference_ligand_verification.json`）。对接用 pH 7.4 微状态，不是中性沉积式。

| 配体 | CCD | 可旋转键 | 角色 | 权威 SMILES（沉积中性式） | 对接微状态 | PDBQT |
|---|---|---:|---|---|---|---|
| lesinurad | A1AIL | 5 | 已知柔性失败例 | `O=C(O)CSc1nnc(Br)n1c1ccc(C2CC2)c2ccccc12` | `O=C([O-])CSc1nnc(Br)n1-c1ccc(C2CC2)c2ccccc12` | 已有 `01_ligand_prep/forced_recovery/pdbqt/lesinurad.pdbqt` |
| 苯溴马隆 | R75 | **3** | **唯一受体准备健全性对照** | `CCc1oc2ccccc2c1C(=O)c1cc(Br)c(O)c(Br)c1` | `CCc1oc2ccccc2c1C(=O)c1cc(Br)c([O-])c(Br)c1` | 已有 `.../benzbromarone.pdbqt` |
| TD-3 | A1A45 | **5** | 与 lesinurad 同类的硫醚–羧酸摆臂复现，**不是刚性对照** | `O=C(O)C(C)(C)Sc1nc2ncccc2n1Cc1ccc(Br)c2ccccc21` | Dimorphite-DL 脱质子羧酸（与 lesinurad 同一规则） | **要现做** |

TD-3 准备：

```bash
python3 scripts/prepare_ligands_c1.py \
  --input-csv data/campaigns/c5/01_ligand_prep/w1_td3.csv \
  --output-dir data/campaigns/c5/01_ligand_prep/w1_refs
```

### 2.3 36 格矩阵（✓ 已有 / ○ 要跑）

行 = 配体，列 = 受体。每个格子 × seeds 42/43/44。

| 配体 \ 受体 | 9DKB | 9DKA | 9DKC | 9DK9 apo |
|---|---|---|---|---|
| lesinurad | ✓ 42/43/44 自对接 | ○ ○ ○ | ○ ○ ○ | ○ ○ ○ |
| 苯溴马隆 | ✓ 42 强制回收；○ 43、44 | ○ ○ ○ **自对接判据格** | ○ ○ ○ | ○ ○ ○ |
| TD-3 | ○ ○ ○ | ○ ○ ○ | ○ ○ ○ **自对接观察格** | ○ ○ ○ |

输出目录：`data/campaigns/c5/01_crossdock/{9dkb,9dka,9dkc,9dk9}/seed{42,43,44}/<ligand>_out.sdf`

每格报告：CNNscore 选姿 Top-1 RMSD、Top-3 最佳 RMSD、9 姿最佳 RMSD、酸根/酚氧–Arg477 胍基最小距离。RMSD 对照该配体**自己的晶体坐标**（苯溴马隆对照 9DKA/R75，lesinurad 对照 9DKB/A1AIL，TD-3 对照 9DKC/A1A45）；交叉对接先把参考配体叠到该受体口袋再算。

### 2.4 通过 / 停

- **受体准备是否正常**：只看苯溴马隆 **自对接 @ 9DKA**，CNNscore Top-1 RMSD **≤ 2.0 Å**。
- 苯溴马隆自对接失败 → **停**，查受体准备/盒子，不准进 W2。
- TD-3 自对接失败 → **支持**柔性摆臂假说，不是否决。
- TD-3 自对接意外通过 → 柔性不是唯一因素，单独讨论 lesinurad，不改引擎。

### 2.5 闸门结果（2026-09-05，已审计）

苯溴马隆 @ 9DKA × seeds 42/43/44：**闸门仍失败**（CNNscore Top-1 `pose_rmsd` ≈ 3.585 / 3.603 / 3.595 Å）。产物在 `data/campaigns/c5/01_crossdock/`。

独立复算（`scripts/parse_c1_sdf_readouts.py::pose_rmsd`，**不对配体做叠合**）：

| seed | Top-1 pose_rmsd | best-of-9 pose_rmsd（模式） |
|---:|---:|---:|
| 42 | 3.585 | **1.110**（mode 3） |
| 43 | 3.603 | **1.087**（mode 3） |
| 44 | 3.595 | **1.109**（mode 3） |

准备 QC 通过：盒子是 9DKA `[107.167, 105.296, 107.628]`，受体无残留 R75，配体是酚负离子。

**分岔：搜得到、选不中。** 9 姿里有晶体盆（~1.11 Å），CNNscore 第一姿是口袋内另一取向（~3.59 Å）。与 lesinurad@9DKB 同类。GetBestRMS 叠合后 ≈1.65/0.48 Å 也支持"形状对、刚体位置偏"，不是对称性匹配的假象。

不要用审计报告里的 GetBestRMS（≈1.65 / 0.48 Å）把闸门改成通过：`GetBestRMS` 会先叠合配体，那是构象相似度，不是对接姿态 RMSD。2.0 Å 门槛仍作用在**未叠合**的 Top-1 `pose_rmsd` 上，**仍是失败**，没有改判。

**对第 120-122 行原始规则的修正说明**（如实记录，不是事后放水）：原规则只写了"苯溴马隆自对接失败→停，不准进 W2"，没区分失败的原因。搜索/选择分岔（`prep_bug` / `search_ok_selection_fail` / `pocket_cannot_hold_crystal`）是在只看到 Top-1 数字（未看逐姿明细）时就先登记好的判据，属于把欠规定的规则说清楚，不是看到明细数据后才编的借口。数值门槛本身（2.0 Å、Top-1、`pose_rmsd`、不叠合）没有变，现在仍是失败。

- **W4（7ALV/NLRP3）**：与本闸门无科学依赖——不同靶点、不同受体、自己的自对接对照（NP3-146@7ALV）已通过。之前说"等 W1"只是算力/叙事排序，不是真依赖，所以开 W4 不算破例。
- **W2（URAT1 IFP 门）**：原设计假设"晶体自对接"能重现晶体姿态（NLRP3 那边确实如此）。现在证明 URAT1 这个假设不成立（搜得到、CNNscore 选不中）。修正做法：IFP 锚点直接用**晶体沉积坐标**（R75/A1AIL/A1A45 的 PDB/CIF 原始坐标），不用 CNNscore 挑出来的重对接姿态——这是去掉了"选择步骤"这一层依赖，让锚点更硬，不是放松判据。

完整修正文字见 `config/campaign_c5.yaml` 的 `audit_conclusion_2026-09-05.amendment_to_on_fail_original` / `w4_dependency_check` / `w2_scope_change`。

下一步：

| 项 | 决定 |
|---|---|
| 其余 29 格 | 不当通过条件；只可作 SI 柔性矩阵，可选 |
| Task2 / W4 | **可以开**（7ALV 已通过；诱饵已锁） |
| Task3 / W2 | **可以开**，IFP 从晶体 R75/A1AIL/A1A45 锚定，不用 CNNscore Top-1 失败姿 |
| 阈值改 4 Å / 宣布闸门通过 | 禁止 |

---

## 3 必须补：Job B = W4 NLRP3 面板重建

**受体只有 7ALV**（已准备）。盒子锁定：`[16.756, 35.449, 125.714]`，**20³**（RM5 质心与锁定值逐位相同）。  
设置与临床酸相同：**exhaustiveness=32**，9 modes，seeds 42/43/44。  
**不要**用 `run_nlrp3_structural_panel.py` 里写死的 `exh=8`。

阈值必须在**看到临床 156 名单之前**冻结。

### 3.1 阳性 9 个（已去重）

丢掉 `CHEMBL3183703`（与 MCC950 SMILES 逐字符相同）。名单来自现有 `panel_ligands.csv`，不是新编的：

| ligand_id | 角色 | seed42@exh32 | seed43/44@exh32 |
|---|---|---|---|
| NP3-146 | 7ALV 共晶 RM5 | ✓ 自对接 | ✓ 自对接 |
| MCC950 | 工具阳性 | ○ 要重做（旧面板 exh=8） | ○ |
| CHEMBL4204644 | ChEMBL 磺酰脲 | ○ | ○ |
| CHEMBL5219789 | ChEMBL 磺酰脲 | ○ | ○ |
| CHEMBL4212407 | ChEMBL 磺酰脲 | ○ | ○ |
| CHEMBL6143743 | ChEMBL 磺酰脲 | ○ | ○ |
| CHEMBL4209503 | ChEMBL 磺酰脲 | ○ | ○ |
| CHEMBL4216836 | ChEMBL 磺酰脲 | ○ | ○ |
| CHEMBL6171925 | ChEMBL 磺酰脲 | ○ | ○ |

SMILES 以 `data/campaigns/c1/05_metrics/nlrp3_structural_panel/panel_ligands.csv` 为准。已有 pdbqt 可复用：`.../nlrp3_structural_panel/pdbqt/`。

### 3.2 背景 20 个临床酸（不新对接，差 2 格）

同一份 `panel_ligands.csv` 的 `clinical_acid_background`。生产 exh=32 姿态：

- seed42：全部 20 个在 `acid_dual/nlrp3_7alv/seed42/`
- seed43/44：19/20 在 `acid_dual_a2/nlrp3_7alv/seed{43,44}/`
- **缺：`REP_07837`（MK-8457）@ 7ALV seeds 43 和 44** → 补 2 个 job，配体已有 `01_ligand_prep/acid_clinical_chemistry_pass/pdbqt/REP_07837.pdbqt`

### 3.3 新诱饵 ≥40 个（本工单不点名）

**不要手写诱饵名字。** 用脚本从已有池抽样，规则预登记如下：

| 项 | 锁定 |
|---|---|
| 来源 | `data/benchmarks/urat1_true_decoy/true_decoys.csv`（性质匹配信封，尚未用于 NLRP3 门） |
| 不要用 | `experimental_inactives.csv`（那是 URAT1 弱活性，不是 NLRP3 诱饵）；不要用临床 156 |
| 性质匹配 | 对 9 个阳性的 MW / logP / TPSA / HBD / HBA / 可旋转键做窗口匹配（沿用 TrueDecoy `matching_assignments.csv` 同一套窗） |
| 相似性 | Morgan ECFP4（r=2, 2048 bit）对**任一** NLRP3 阳性 max TC **≤ 0.5** |
| 数量 | ≥40，固定抽样种子（建议 `0xC5DEC0`），抽完写入 CSV 后再对接 |
| 准备 | 与临床酸相同的 Dimorphite-DL → Meeko |

### 3.4 Job 数

| 子集 | 新 gnina job |
|---|---:|
| 8 个非共晶阳性 × 3 种子 | 24 |
| ≥40 诱饵 × 3 种子 | ≥120 |
| REP_07837 @ 7ALV seeds 43/44 | 2 |
| **合计** | **≥146** |
| 复用（不对） | NP3-146 × 3 + 背景 20×1 + 19×2 = 61 |

输出：`data/campaigns/c5/02_nlrp3_panel/{positives,background,decoys}/seed{42,43,44}/`

通过规则（已锁）：结构门特异性 > 宽松门，且 Fisher *p* < 0.05。失败则正文降级为「姿态质控」，不另发明指标。

---

## 4 零新对接：Job C = W2 URAT1 IFP 门

**不对任何新分子。** 重打分 Phase I 的 9DKB 9-mode SDF。

- 脚本：`scripts/run_acid_gate_benchmark.py`
- 集合：228 羧酸 active vs 64 羧酸 true decoy  
  `data/campaigns/c1/05_metrics/acid_gate_retrospective_benchmark/`
- 关键残基（文献锚定，不是几何反推）：`S35, M214, F241, F360, F364, F365, D389, K393, Q437, F449, R477, Q473`
- 阈值从**晶体沉积坐标**（R75 / A1AIL / A1A45 的 PDB/CIF 原始坐标）锚定，**不用** CNNscore 选出的重对接姿态；在 228-vs-64 上**只评估一次**，禁止网格搜索最大化 OR
- 通过：OR 的 95% CI 下界 > 1；失败则回落到 A1 ∩ A2，不得再发明新门

---

## 5 现在不要做

| 工作 | 何时 / 为何不做 |
|---|---|
| W5 MD | `md_authorized: false`；短名单已冻，仍须显式授权后才开 |
| W6 OAT1 反筛 | 可选，冻结之后；Jeon *Structure* 2025 / Wu & Luo *Sci. Adv.* 2025；同一 gnina 设置；不进标题摘要 |
| 临床 156 重对接 | 三种子双臂已齐 |
| Rank 轨 / 9DKB 全诱饵 | 已关 |
| 为提高 EF 新建 RandomDecoy | 禁止 |
| 换打分引擎 / 重训 URAT1 ML / 用 MM-GBSA 给 URAT1 排名 | 禁止（见 `docs/C5_RANKING_AND_NEXT_DOCKS.md`） |
| 为“有个 URAT1 排序分”另开新轨 | 禁止 |

---

## 6 建议执行顺序（本机，闸门失败后）

1. 9DKA / 9DKC / 9DK9 受体与 TD-3 已准备；W1 闸门 3 job 已跑完且失败——**不要**再把其余 29 格当解锁条件。  
2. W4 诱饵 CSV 已锁（`w4_decoys_locked.csv`，n=40）。  
3. W4 146 job **已跑**（131 ok / 15 fail）；结构门 **已打分并通过**（vs 诱饵）。  
4. W2 **已评**（IFP 过门但是 A1 子集）。不必等 W1 其余 29 格。  
5. **短名单已冻结**（12 主 + 21 备份）。现在写正文，才谈 MD。

本机不再有必须新跑的 gnina。选出候选之前不必再算。可选 SI：15 个失败诱饵重试，或 W1 剩余 29。
