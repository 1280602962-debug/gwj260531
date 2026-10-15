# 七次参照重对接审计（停在本批）

日期：2026-10-15  
工作树：`/home/gwj/work/stage3b-refdock`  
分支：`cursor/stage3b-refdock-20261015`  
起点提交：`95ee46a0e6c9b7bd126bcd2942b2c304db19f4aa`（`cursor/nlrp3-model-calibration-audit-20261014`）  
原检出 `/home/gwj/work/c1-campaign` 仍停在 `cursor/urat1-nlrp3-dualtarget-aidd-e43d`，历史对接文件未改写。

本批之后没有启动 32 个新分子或 192 次正式对接。种子 44 没有重跑。

## 软件与输入

GNINA：`/mnt/d/CADD paper exercise/gnina/bin/gnina`，`v1.3.2 master:f23dd2b`，Built Jul 8 2025。与历史 stdout 中的构建字符串相同。运行时带上该安装自带的 `conda_env/lib`，并使用 `--no_gpu`。

四个配体 PDBQT 与两份受体 PDBQT，提交内文件与历史命令所用路径逐字节相同，没有重新制备配体。

| 任务 | 配体 PDBQT | 受体 | 盒子 | 种子 |
|---|---|---|---|---|
| Lesinurad × 9DKB | `forced_recovery/pdbqt/lesinurad.pdbqt` | `9DKB_receptor.pdbqt` | 99.966, 102.967, 105.699 / 22³ | 42 |
| NP3-146 × 7ALV | `np3_146/pdbqt/NP3-146.pdbqt` | `7ALV_receptor.pdbqt` | 16.756, 35.449, 125.714 / 20³ | 42 |
| REP_00015、REP_00300 | `acid_clinical_chemistry_pass/pdbqt/` | 同上 | 同上 | 42 或 43 |

参数：`exhaustiveness=32`，`num_modes=9`，`cnn_scoring=rescore`，`cpu=4`。

URAT1 后处理使用 `c1_acid_pose_selection.evaluate_urat1_acid_sdf(..., rule="a2")`：先要求羧酸氧–Arg477 ≤ 7.7027 Å 且质心 ≤ 6.0 Å，再在合格姿势中取最高 CNNscore。NLRP3 使用 `c1_nlrp3_pose_metrics.evaluate_nlrp3_structural`：最高 CNNscore 姿势上同时检查 CNNscore ≥ 0.5、质心 ≤ 6.0 Å、重叠 ≥ 0.50、IFP ≥ 0.50、关键残基 ≥ 5/7、无 < 2.2 Å 碰撞。

## 已有酸性库覆盖

提交内计数与 `docking_protocol_audit.json` 一致：9DKB 与 7ALV 的种子 42 各 156 个 SDF；9DKB 种子 43、44 各 156 个；7ALV 种子 43、44 各 155 个。缺的是 REP_07837 在 7ALV 上的种子 43 和 44。

REP_07837 的 9DKB 种子 42/43/44 和 7ALV 种子 42 的 SDF 都在。7ALV 种子 43、44 的 SDF 不在提交里。旧工作目录里有两份未纳入提交的日志，命令已写出，但内容停在 `Using random seed`，没有打分表，也没有对应 SDF。这是计算未完成，不是几何不通过。本阶段没有补算。

阶段 3A 把已入库队列记为 118 个分子，其中羧酸 103 个；磁盘上的临床 SDF 是 156 个化学通过分子。这两套数字不是同一次清点。本报告的文件计数按 156/155，没有把缺文件记成几何失败。

## 七次对照

姿势按重原子空间匹配，不按 SDF 行号。七个任务的文件顺序也恰好一致。SDF 正文与历史文件的差别只有每条记录的 OpenBabel 时间戳（例如 `OpenBabel08262610173D` 与 `OpenBabel10152610253D`）。坐标和分值相同。

严格数值标准：匹配姿势的 CNNscore 绝对差 ≤ 0.0001。七次的最大差都是 0，最大重原子 RMSD 都是 0。

| 任务 | 组别 | 姿势 | 最大 CNNscore 差 | 最大 RMSD (Å) | 历史门 | 新门 | 结论 |
|---|---|---:|---:|---:|---|---|---|
| Lesinurad × 9DKB × 42 | 晶体对照 | 9 | 0 | 0 | A2 过，模式 6，4/9 几何合格 | 相同 | 分数与 A2 结论都一致 |
| NP3-146 × 7ALV × 42 | 晶体对照 | 9 | 0 | 0 | 结构门过，模式 1 | 相同 | 分数与结构门都一致 |
| REP_00015 × 9DKB × 42 | 临床生产 | 9 | 0 | 0 | A2 过，模式 3，5/9 | 相同 | 一致 |
| REP_00015 × 7ALV × 42 | 临床生产 | 9 | 0 | 0 | 结构门过，模式 1 | 相同 | 一致 |
| REP_00300 × 9DKB × 42 | 临床生产 | 9 | 0 | 0 | A2 无合格几何，0/9 | 相同 | 一致，结论是不通过 |
| REP_00015 × 9DKB × 43 | 临床生产 | 9 | 0 | 0 | A2 过，模式 1，7/9 | 相同 | 一致 |
| REP_00300 × 9DKB × 43 | 临床生产 | 9 | 0 | 0 | A2 无合格几何，0/9 | 相同 | 一致，结论是不通过 |

### 晶体对照（不与临床生产混写）

Lesinurad：A2 选中模式 6，CNNscore 0.670793176，羧酸氧–Arg477 3.93 Å，质心 2.91 Å。九个姿势里 4 个过几何。历史 `lesinurad_metrics.json` 记录的是最高 CNNscore（模式 1）相对晶体失败：酸根距离 14.11 Å，RMSD 4.30 Å。那是另一条选择规则。同一套九个姿势在 A2 下新旧都通过。数值复现成立，不能据此改写当年那份 CNNscore-first 记录。

NP3-146：最高 CNNscore 为模式 1，0.971881032；质心 0.36 Å，口袋重叠 1.00，IFP 0.842，关键接触 6/7，无严重碰撞。结构门通过。历史 `NP3-146_metrics.json` 写有 `rmsd_error: no RMSD match`，那是当时 RMSD 程序没配上原子，不是这次结构门的结论。

### 临床生产

REP_00015 在种子 42 的两个靶点都通过各自的门，种子 43 的 9DKB 也通过。种子 44 没有重跑。REP_00300 在种子 42 和 43 的 9DKB 上都是九个姿势里没有一个同时满足酸根距离和质心条件，A2 不选出姿势。这与当初把它当作几何不通过参照的记录一致，而且姿势文件是完整的。

## 可否并入后续新分子

这七次使用的是同一 GNINA 构建、同一受体、同一配体 PDBQT 和同一盒子。匹配姿势的 CNNscore 差为 0，A2 与 NLRP3 完整结构门的通过与否、选中模式都与历史 SDF 相同。按预定的 0.0001 标准，这七次可以视为严格数值复现，筛选结论也可比。

仍不能外推的部分：

- 种子 44 没有复现。
- REP_07837 在 7ALV 种子 43、44 的姿势仍然缺失，需单独补算，不能并进 32 个新分子。
- 早期约 20 Å 盒子、Vina、RTMScore 和最高 CNNaffinity 汇总不在本次可比范围内。
- 临床任务的原始命令行来自旧工作目录中未提交的 stdout；分数以提交内 SDF 为准，stdout 与 SDF 一致。

逐姿势表：`pose_matches.csv`。任务表：`task_comparison.csv`。完整字段：`comparison_summary.json`。每次的命令、stdout、`gnina.log`、输入 PDBQT 和 `poses_out.sdf` 在对应子目录。
