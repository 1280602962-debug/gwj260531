# ChEMBL 临床药物重定位库

| 文件 | 说明 |
|------|------|
| `repurposing_manifest.csv` | 主库 **8319** 条（ChEMBL Phase + ATC 导出合并） |
| **`screening/docking_pool_p05.csv`** | **P≥0.5 对接池，1588 条**（已提交 Git） |
| `screening/` | NLRP3 ML 全库分数与 SI 子集（见 `screening/README.md`） |
| **`p2/`** | **gnina P2 生产漏斗归档**：双靶完整案例 1,580、Pareto 审计、化学提名 |

## 生成对接池

```bash
python3 scripts/screen_repurposing_library.py \
  --input data/repurposing/repurposing_manifest.csv \
  --panel clinical_all --export-p05-pool --skip-tanimoto
```

输出默认写入 `results/repurposing/`；请同步复制到 `data/repurposing/screening/` 后提交 Git。

## 勿与以下混淆

- `data/distill/distill_manifest.csv` — 8973，**仅 URAT1 回顾**
