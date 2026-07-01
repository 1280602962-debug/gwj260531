# ChEMBL 临床药物重定位库

| 文件 | 说明 |
|------|------|
| `repurposing_manifest.csv` | 主库 **8319** 条（ChEMBL Phase + ATC 导出合并） |

## 生成对接池

```bash
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto
```

输出在 `results/repurposing/`（gitignore）：

- `nlrp3_ml_scores_clinical_all.csv`
- `docking_pool_p05.csv`（P≥0.5，约 1588 条）

## 勿与以下混淆

- `data/distill/distill_manifest.csv` — 8973，**仅 URAT1 回顾**
