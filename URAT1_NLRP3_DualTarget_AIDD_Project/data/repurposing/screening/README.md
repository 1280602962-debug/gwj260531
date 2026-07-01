# NLRP3 ML 筛选输出（已提交 Git）

由 `scripts/screen_repurposing_library.py` 生成，与 `results/repurposing/` 内容同步。

| 文件 | 说明 |
|------|------|
| `nlrp3_ml_scores_clinical_all.csv` | 全临床库 **8319** 条 NLRP3 ML 分数 |
| `docking_pool_p05.csv` | **P(active)≥0.5 对接池，1588 条**（主流程 Maestro 输入） |
| `nlrp3_screening_summary_clinical_all.json` | 全库统计 + 对照药排名 |
| `nlrp3_top_for_dual_docking_clinical_all.csv` | Top 5% / top-N 短名单（可选对照） |
| `nlrp3_ml_scores_phase_ge3.csv` | III 期+上市子集（1283）分数 |
| `docking_pool_p05_phase_ge3.csv` | 子集 P≥0.5，**247 条**（SI 敏感性） |
| `nlrp3_screening_summary_phase_ge3.json` | 子集统计 |

## 重新生成

```bash
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto

python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel phase_ge3 --export-p05-pool --skip-tanimoto
```
