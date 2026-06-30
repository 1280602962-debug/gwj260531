# Benchmark 选择标准与合理性说明

详见数据文件：[`data/benchmarks/literature_benchmarks.csv`](../data/benchmarks/literature_benchmarks.csv)

---

## 六条选择标准

1. **领域锚定**：临床批准药、临床候选或共晶工具化合物  
2. **骨架/机制多样**：URAT1 四药化学型与结合模式不同  
3. **独立性分层**：Tier1a（训练集外）vs Tier1b（训练集内 sanity）  
4. **标签可靠性**：文献原刊 IC50 优先于 ChEMBL 中位值  
5. **测定一致**：NLRP3 用 IL-1β 细胞读出；URAT1 用摄取/转运  
6. **阴性机制特异**：allopurinol（XO）、colchicine（间接 NLRP3）

---

## 分层考试

| 考试 | 分子 | Pass |
|------|------|------|
| 结构 URAT1 | 四药 | redock ≤2.0–2.5 Å；π_in > π_out ≥3/4 |
| ML 外推 URAT1 | lesinurad, dotinurad, benzbromarone | L2 Top-500 ≥2/3 |
| Sanity | verinurad, MCC950, GDC-2394 | 必过（不证明外推） |
| 阴性 | allopurinol, colchicine | 排名后 20% |

---

## 已知数据冲突（须在论文 Limitations 说明）

| 分子 | 冲突 | 处理 |
|------|------|------|
| lesinurad | ChEMBL ~5.1 vs Burns 3.53 µM vs Dai EM 39 µM | 分 assay 报告；结构考试优先 |
| dotinurad | ML 预测偏低 | 作 scaffold-novel 硬测试，非 benchmark 选错 |
