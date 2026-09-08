# 双靶重定位药物库指南（当前路线）

> **主流程**：[`WORKFLOW.md`](WORKFLOW.md)  
> **不是**：8973 distill；不是 Enamine 百万库

---

## 1. 主库文件

| 文件 | n | 说明 |
|------|---|------|
| `data/repurposing/repurposing_manifest.csv` | 8319 | ChEMBL 临床阶段 + ATC 合并 |
| `library_panel=primary_atc_phase` | 449 | 疾病 ATC ∩ III期/上市（SI） |

---

## 2. 筛选流程（当前标准）

```
8319 manifest
  → NLRP3 ML（screen_repurposing_library.py）
  → P(active) ≥ 0.5  →  docking_pool_p05.csv（n≈1588）
  → URAT1 9DKB + NLRP3 7ALV（gnina P2）
  → Pareto
```

**为何用 P≥0.5**：与训练 binary active 阈值一致；规模对标重定位文献（~1500–2700）；算力可接受。

**SI 敏感性**：`max_phase≥3` 子集（247）；Top 5%（416）。

---

## 3. 构建 manifest（从 ChEMBL Excel）

见 [`LOCAL_AGENT_REPURPOSING_LIBRARY_PROMPT.md`](LOCAL_AGENT_REPURPOSING_LIBRARY_PROMPT.md) 与 `scripts/build_repurposing_library.py`。

---

## 4. 与 8973 分工

| 库 | 用途 |
|----|------|
| **repurposing_manifest** | NLRP3 ML + 双靶对接 + Pareto |
| **distill_manifest 8973** | 仅 URAT1 9DKB 回顾（A vs D） |

---

## 5. 脚本

| 脚本 | 作用 |
|------|------|
| `screen_repurposing_library.py` | NLRP3 ML + 导出对接池 |
| `merge_docking_pareto.py` | 对接分 + Pareto（对接完成后） |
| `build_repurposing_library.py` | 从 ChEMBL Excel 建 manifest |
