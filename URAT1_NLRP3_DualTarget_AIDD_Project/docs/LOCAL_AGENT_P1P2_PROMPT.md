# 本地 Agent 离线完成 P1/P2 指南

> **⚠️ 已归档（2026-07）**：本文档描述 **OAT 迁移 / P1P2 旧路线**，**不再作为论文主流程**。  
> 当前工作流：[`WORKFLOW_CURRENT.md`](WORKFLOW_CURRENT.md) | 归档说明：[`LEGACY_ARCHIVE.md`](LEGACY_ARCHIVE.md)

**适用场景**：本地 agent 无法访问 GitHub，项目已通过 U 盘 / ZIP / `git pull`（由你手动执行）同步到本机。

---

## 一、先决条件（人工完成，不需 agent 连网）

1. 把整份项目拷到本地，目录名保持 `URAT1_NLRP3_DualTarget_AIDD_Project/`
2. 确认以下文件存在：
   - `data/processed/urat1_curated.csv`（没有则先跑数据准备）
   - `data/auxiliary/oat_combined_transfer.csv`（没有则跑 `00b_prepare_auxiliary_data.py`）
3. Python 3.10+，`pip install -r requirements.txt`

---

## 二、直接复制给本地 Agent 的命令（推荐）

把下面整段粘贴给本地 Cursor / Agent：

```
你在离线环境工作，不要尝试 git clone / push / 访问 GitHub。

项目根目录：URAT1_NLRP3_DualTarget_AIDD_Project/

任务：完成 P1 + P2

【P1】OAT 迁移训练
- 脚本：scripts/02_train_asymmetric_models.py
- 已支持 --oat-transfer / --no-oat-transfer
- OAT 数据：data/auxiliary/oat_combined_transfer.csv（73 条，与 URAT1 重叠 13 SMILES）
- 方法：每折 CV 先在 OAT 上 XGBoost 预训练，再在 URAT1 train fold 上 sequential fine-tune
- NLRP3 训练逻辑不变

【P2】消融对比
- 脚本：scripts/run_oat_transfer_ablation.py
- 自动跑 baseline（无 OAT）vs transfer（有 OAT），各做一次 benchmark
- 输出：results/training/oat_transfer_ablation.json

执行顺序（在项目根目录的 scripts/ 下或从项目根用 python3 scripts/...）：

cd URAT1_NLRP3_DualTarget_AIDD_Project
pip install -r requirements.txt

# 若缺 processed 数据：
python3 scripts/00_prepare_data.py

# 若缺 auxiliary OAT 数据（且你有 data/raw/auxiliary/OAT1_ALL*.csv）：
python3 scripts/00b_prepare_auxiliary_data.py --copy-raw

# P2 一键消融（含 P1 两次训练 + benchmark）：
python3 scripts/run_oat_transfer_ablation.py

# 或分步：
python3 scripts/02_train_asymmetric_models.py --no-oat-transfer --output results/training/ablation_no_oat
python3 scripts/02_train_asymmetric_models.py --oat-transfer --output results/training/ablation_oat_transfer
python3 scripts/07_benchmark_backtest.py --model-dir results/training/ablation_no_oat --output results/benchmark_backtest/ablation_no_oat
python3 scripts/07_benchmark_backtest.py --model-dir results/training/ablation_oat_transfer --output results/benchmark_backtest/ablation_oat_transfer

完成后汇报：
1. oat_transfer_ablation.json 中 baseline vs transfer 的 Spearman、R²、benchmark x/4
2. URAT1 是否仍为 URAT1_NO_GO
3. 若有报错，贴完整 traceback
```

---

## 三、若本地代码较旧（还没有 P1/P2 实现）

把下面「实现规格」段交给 agent，让它**只改本地文件、不连 GitHub**：

```
在 scripts/02_train_asymmetric_models.py 实现 OAT 迁移：

1. 读取 data/auxiliary/oat_combined_transfer.csv
2. 新增 CLI：--oat-transfer（默认开）/ --no-oat-transfer
3. train_urat1_cv() 每折：
   - scaler 在 URAT1 fit fold 上 fit
   - OAT 预训练排除当前 test fold 的 SMILES（防泄漏）
   - model.fit(x_fit, y, xgb_model=pretrain_booster) 做 sequential fine-tune
4. fit_urat1_final() 同样支持 OAT 预训练
5. training_report.json 记录 transfer_learning 字段

新建 scripts/run_oat_transfer_ablation.py：
- 跑 --no-oat-transfer 和 --oat-transfer 各一次
- 各跑 07_benchmark_backtest.py
- 写 results/training/oat_transfer_ablation.json

参考 config/model_hierarchy.yaml 中 urat1_arm.transfer_learning.method=sequential_finetune
事实数据以 docs/DATA_FACT_CHECK.md 为准，不要夸大 OAT 样本量（合并约 73 条）。
```

---

## 四、预期输出文件

| 文件 | 说明 |
|------|------|
| `results/training/ablation_no_oat/training_report.json` | 无迁移 CV + 筛选判定 |
| `results/training/ablation_oat_transfer/training_report.json` | 有迁移 |
| `results/training/oat_transfer_ablation.json` | **P2 对比摘要** |
| `results/benchmark_backtest/ablation_*/benchmark_backtest_report.json` | benchmark 2/4 等 |

---

## 五、常见问题

**Q: `data/processed/` 为空？**  
先运行 `python3 scripts/00_prepare_data.py`（需要 `data/raw/` 下 ChEMBL 导出）。

**Q: 找不到 `oat_combined_transfer.csv`？**  
运行 `python3 scripts/00b_prepare_auxiliary_data.py`，或从已处理机器拷贝 `data/auxiliary/`。

**Q: agent 想 push 到 GitHub？**  
告诉他跳过；你手动 `git pull` / 拷贝 `results/` 回云仓库即可。
