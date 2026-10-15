# 第三阶段 B 本地执行报告

范围锁定为 194 次 GNINA，没有扩大集合，没有重做已有 156 个分子，没有改阈值，没有用分子动力学结果改规则。本报告之后停止，不进入文稿整合。

## 任务完成情况

194 次全部完成。状态只有 `completed`：32 个新分子 × 2 个靶点 × 3 个种子 = 192，另加 REP_07837 在 NLRP3 7ALV 上的种子 43 和 44。没有对接失败，没有超时，没有少于 9 个姿势的输出。`num_modes 9` 在这批任务里都达到了；若某次少于 9 个姿势，仍会记为对接成功，本批没有这种情形。

单次耗时 67.4–280.8 秒，中位数 146.3 秒。并行保持 2 个任务 × 4 线程。软件与盒子与严格复现一致：GNINA v1.3.2 master:f23dd2b，`--no_gpu --exhaustiveness 32 --cnn_scoring rescore`，URAT1 盒子中心 (99.966, 102.967, 105.699)、边长 22 Å，NLRP3 盒子中心 (16.756, 35.449, 125.714)、边长 20 Å。受体沿用已有 PDBQT，没有重新准备。

## 配体准备

32/32 准备成功，0 个失败。Dimorphite-DL pH 7.4 后，32 个都按已有规则选了去质子羧酸微物种（`prefer_deprotonated_carboxylate`）。微物种数：1 个的 18 个分子，2 个的 11 个，3 个的 1 个，4 个的 2 个。形式电荷：−1 有 24 个，−2 有 3 个（OLSALAZINE、MEBROFENIN、ARCLOFENIN），0 有 5 个。电荷为 0 的五个是两性离子，羧酸仍为 `[O-]`，同时含 `[NH+]`：RAZUPENEM、ISOPROPICILLIN、CEDIROGANT、EVATEGREL、TEBIPENEM。没有因为对接结果改选微物种。

6 个分子的选定微物种仍有未指定立体化学，保持原样：RAZUPENEM（1）、ISOPROPICILLIN（1）、BENFURODIL HEMISUCCINATE（1）、FENBENICILLIN（1）、CEDIROGANT（2）、EVATEGREL（1）。

## 32 个新分子的三种子评价

URAT1 用 A2：羧酸氧到 Arg477 ≤ 7.7027 Å 且重原子质心到晶体配体 ≤ 6.0 Å 的姿势里，取 CNNscore 最高者。NLRP3 用完整结构门：最高 CNNscore 姿势需 CNNscore ≥ 0.5、质心 ≤ 6.0 Å、口袋重叠 ≥ 0.50、IFP Jaccard ≥ 0.50、7 个关键残基中至少 5 个接触、且无 < 2.2 Å 冲突。同一种子上两个靶点都通过，才算该种子双通过。结构合格要求至少 2 个双通过种子。化学过滤沿用第三阶段 A 已冻结的 β-内酰胺、PAINS、Brenk 标记。

96 个 URAT1 对接和 96 个 NLRP3 对接全部完成。URAT1 A2 通过 88/96 个种子，8 个种子是对接完成但没有几何合格姿势。NLRP3 结构门通过 65/96 个种子，31 个种子是对接完成但结构门未过。双通过种子 62/96。

分子层面：

- 结构合格且化学标记干净：13
- 结构合格但化学标记不通过：8
- 对接完成但未达到 2 个双通过种子：11
- 对接失败：0
- 超时：0

这 11 个不是“没有 URAT1 活性”。其中多数 URAT1 A2 已通过，未入选是因为同一种子上 NLRP3 结构门未过，或双通过种子不足 2 个。

### 结构合格且化学标记干净（13）

| ID | 名称 | 双通过种子 |
| --- | --- | --- |
| REP_06153 | ARONIXIL | 3 |
| REP_03766 | BENFURODIL HEMISUCCINATE | 3 |
| REP_04062 | BROXITALAMIC ACID | 3 |
| REP_06887 | CEDIROGANT | 3 |
| REP_07890 | GSK-065 | 3 |
| REP_07657 | GSK-625433 | 3 |
| REP_01548 | INDEGLITAZAR | 3 |
| REP_01934 | MEBROFENIN | 3 |
| REP_03072 | METANIXIN | 3 |
| REP_05371 | PIRINIXIC ACID | 3 |
| REP_01012 | RECOFLAVONE | 3 |
| REP_00834 | TETOMILAST | 3 |
| REP_00025 | FEBUXOSTAT | 2 |

CEDIROGANT 和 BENFURODIL HEMISUCCINATE 的选定微物种仍有未指定立体化学，记录在案，没有补指派。

### 结构合格但化学标记不通过（8）

沿用冻结标记，没有重调阈值。β-内酰胺：RAZUPENEM（2 个双通过种子）、FENBENICILLIN（2）、TEBIPENEM（3）。`alert_clean_pains_brenk` 为假：SUXIBUZONE、OLSALAZINE、N6022、ENOXAMAST、SUCCISULFONE，均为 3 个双通过种子。

### 对接完成但结构未合格（11）

| ID | 名称 | 种子 42 | 种子 43 | 种子 44 | 说明 |
| --- | --- | --- | --- | --- | --- |
| REP_02627 | ARCLOFENIN | U过 / N未过 | 双通过 | U过 / N未过 | 只有 1 个双通过种子 |
| REP_00694 | RITOBEGRON | U无几何 / N过 | 双通过 | U无几何 / N过 | 只有 1 个双通过种子 |
| REP_04189 | MALEYLSULFATHIAZOLE | U过 / N未过 | U过 / N未过 | U无几何 / N过 | 没有同一种子双通过 |
| REP_01829 | SUCCINYLSULFATHIAZOLE | U过 / N未过 | U过 / N未过 | U过 / N未过 | NLRP3 三次均为 CNNscore < 0.5 |
| REP_05301 | CINALUKAST | U过 / N未过 | U过 / N未过 | U过 / N未过 | NLRP3 三次均为 CNNscore < 0.5 |
| REP_06737 | CUDETAXESTAT | U过 / N未过 | U过 / N未过 | U过 / N未过 | NLRP3 三次 CNNscore 约 0.18 |
| REP_03673 | GALTIFENIN | U过 / N未过 | U过 / N未过 | U过 / N未过 | NLRP3 为 CNNscore 或关键接触不足 |
| REP_07314 | EVATEGREL | U过 / N未过 | U过 / N未过 | U过 / N未过 | NLRP3 的 IFP 与关键接触不足 |
| REP_03178 | ISOPROPICILLIN | U过 / N未过 | U过 / N未过 | U过 / N未过 | NLRP3 关键接触不足；β-内酰胺标记仍保留，但未进入化学提名 |
| REP_03265 | DIODONE | U无几何 / N未过 | U过 / N未过 | U无几何 / N未过 | NLRP3 的 IFP 与关键接触不足 |
| REP_05525 | GSK-2636771 | 双未过 | 双未过 | 双未过 | URAT1 三次都没有几何合格姿势；NLRP3 三次 CNNscore 略低于 0.5 |

URAT1 没有几何合格姿势的 8 个种子是：RITOBEGRON 42/44，DIODONE 42/44，MALEYLSULFATHIAZOLE 44，GSK-2636771 42/43/44。这些种子都写出了 9 个姿势。

## REP_07837（MK-8457）补算

历史日志在 NLRP3 种子 43 和 44 停在“Using random seed”，没有分数表，也没有 SDF。那是未完成的计算，不是几何不合格。这两次已补算，各 9 个姿势，262 秒左右，输出在 `gap_REP_07837/rerun/`，与历史未完成日志分开存放。URAT1 三个种子仍用已存姿势，没有重对接。

统一规则下：URAT1 A2 三个种子都通过；NLRP3 三个种子都未通过，原因都是最高 CNNscore 姿势约 0.27–0.28，低于 0.5。重叠、IFP 和关键接触在这三次里是够的。双通过种子为 0，结构不合格。种子 43 和 44 的 NLRP3 姿势来源标记为 `stage3b_gap_fill`，不能当成历史原件。

## 与已有候选的对照

已存 156 个分子的 SDF 用同一套 A2 和 NLRP3 结构门重读，没有重对接。156 × 3 个种子的姿势文件都在，没有缺文件。

冻结的 12 个主候选和 21 个备用候选，共 33 个，在这套统一规则下全部仍是结构合格且化学标记通过。主候选里 PF-03882845 等 10 个为 3 个双通过种子，ADMILPARANT 和 TONAPOFYLLINE 为 2 个。没有用分子动力学轨迹改变任何阈值。

另有 21 个已存分子达到结构门，但不在冻结化学表里，因此记为 `structural_pass_chemistry_not_in_frozen_table`，没有自动提名。新的 13 个化学通过分子与这 33 个已有候选没有 ID 重叠。

### PF-03882845（REP_07580）

三个种子都是同一种子双通过，结构合格，化学标记通过。姿势全部来自已存 SDF，没有重对接，也没有用其分子动力学结果改规则。

| 种子 | URAT1 A2 | 羧酸–Arg477 | 质心 | NLRP3 CNNscore | 重叠 | IFP | 关键接触 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 42 | 通过，CNNscore 0.643 | 3.24 Å | 1.33 Å | 0.624 | 0.867 | 0.85 | 7/7 |
| 43 | 通过，CNNscore 0.644 | 3.23 Å | 1.33 Å | 0.620 | 0.867 | 0.85 | 7/7 |
| 44 | 通过，CNNscore 0.643 | 3.24 Å | 1.33 Å | 0.629 | 0.867 | 0.85 | 7/7 |

三次 NLRP3 都没有冲突。

## 规则不适用的 7 个非羧酸

未对接。酰磺胺或四唑，羧酸–Arg477 规则不适用。这不是“没有 URAT1 活性”。

REP_00615 TASISULAM SODIUM，REP_00700 TASISULAM，REP_01949 REPINOTAN，REP_03009 QUINOTOLAST，REP_04219 GALOSEMIDE，REP_07443 SULOFENUR，REP_07457 DU 125530。

## 写出位置

表格在 `results/nlrp3_model_v2/stage3b_local/analysis/`：`task_execution.csv`、`new32_per_seed.csv`、`new32_summary.csv`、`historical_per_seed_rescored.csv`、`historical_summary_rescored.csv`、`PF03882845_REP_07580_per_seed.csv`、`excluded_noncarboxylate.csv`、`unified_structural_candidates.csv`、`counts.json`。工作区里新增内容只在 `stage3b_local/` 下。历史姿势、1588 库、1525 v2 库、冻结 12/21 名单、模型文件和 PF-03882845 分子动力学轨迹都没有改写。
