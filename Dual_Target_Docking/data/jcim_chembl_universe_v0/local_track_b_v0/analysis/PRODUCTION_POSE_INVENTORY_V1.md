# 新增五对主 Vina 姿态清单

基线提交：`1ea95aa5`。不重对接。

## 仓库里有什么

| 文件 | 状态 |
|---|---|
| `tables/scores_vina_mode1_v1.csv` | 在 git 中，1094 条成功分数 |
| `tables/job_status.csv` | 在 git 中，1100 条任务（6 条超时 skip） |
| `poses/` | **未纳入 git**（`.gitignore` 排除） |
| `logs/vina/` | **未纳入 git** |

## 本地磁盘核验（2026-09-12）

回放脚本：`scripts/replay_track_b_vina_mode1_v1.py`  
结果表：`tables/scores_vina_mode1_pose_replay_v1.csv`

| 项目 | 数字 |
|---|---:|
| 任务总数 | 1100 |
| 成功并写入分数 | 1094 |
| 本地 `mode_01.pdbqt` 找到 | 1094 |
| `REMARK VINA RESULT` 与 CSV 一致 | **1094 / 1094** |
| 姿态缺失 | 0（仅限本机当前磁盘） |
| 超时 skip | 6 |

6 条超时（不是 `TORSDOF ≥ 25`）：

- `F2F10_080` @ 2JKH、4UDW
- `F2F10_105`、`F2F10_106` @ 4UDW
- `J1TYK2_092` @ 3LXP
- `PGPA_030` @ 9V8H

## 含义

本机可以按保存姿态回放主分数，说明这些数字来自已跑过的 Vina 输出，不是后写的空表。  
**公开 git 提交里仍然没有这 1094 个姿态文件。** 审稿人或第三方只克隆仓库时，不能从 `REMARK VINA RESULT` 重建分数。姿态应另作只读归档（例如 Zenodo），不能用一次新对接冒充旧运行。

共晶姿态不在此限制内：8 个 `cognate_qc/*_out_E8.pdbqt` 已在 git 中。
