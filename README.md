# gwj260531

## 项目目录

| 项目 | 说明 |
|------|------|
| [**JNK1_Selectivity_Project/**](JNK1_Selectivity_Project/README.md) | JNK1/2/3 亚型选择性抑制剂 CADD/AIDD 完整流程 |
| [**URAT1_NLRP3_DualTarget_AIDD_Project/**](URAT1_NLRP3_DualTarget_AIDD_Project/README.md) | 痛风 URAT1–NLRP3 双节点重定位（**V2**：TrueDecoy 协议筛选 + 不对称漏斗） |

### JNK1 快速开始

```bash
cd JNK1_Selectivity_Project
pip install -r requirements.txt
python3 scripts/07_compare_models.py --skip-prepare --skip-similarity --skip-chemprop
python3 scripts/06_virtual_screening.py --library data/libraries/your_library.csv --output results/screening_v2
```

### URAT1/NLRP3 双节点项目（当前 V2）

```bash
cd URAT1_NLRP3_DualTarget_AIDD_Project
# 论文与写作入口（勿打开已过时的 MANUSCRIPT_DRAFT_CN.md / Glide XP 稿）
cat docs/MANUSCRIPT_OUTLINE_V2.md
cat docs/INTRO_DRAFT_CN.md
cat docs/METHODS_DRAFT_CN.md
```

> 注意：仓库默认分支 `main` 仍以 JNK 项目为主说明；**URAT1 最新文稿在 PR 分支** `cursor/urat1-nlrp3-dualtarget-aidd-e43d`（或合并后的对应目录）。
