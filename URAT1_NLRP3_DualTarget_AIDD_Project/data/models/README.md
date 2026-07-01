# 训练模型（已提交 Git）

由 `scripts/02_train_asymmetric_models.py --no-oat-transfer` 生成。

| 文件 | 说明 |
|------|------|
| `nlrp3_model.joblib` | NLRP3 assay-conditioned 分类器（筛选用） |
| `urat1_model.joblib` | URAT1 回归（benchmark / SI only） |
| `training_report.json` | CV 指标与门槛检查 |

运行时脚本默认从 `results/training/` 读取；本目录为 **Git 备份副本**，便于克隆后直接筛选。
