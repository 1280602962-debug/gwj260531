# JNK1 项目汇报 — 演讲备注（备用幻灯片）

> 与 `JNK1_Project_Presentation_v2.9.pptx` 内演讲者备注同步。重新生成：`python3 scripts/build_project_pptx.py`

## 第 1 页｜【备用】ML 外部 Decoy 验证

### 幻灯片要点

- Recall@6.0：99.3%
- Decoy FPR：95.3%；Specificity：4.7%
- ROC-AUC：0.876；EF1%：9.20
- 结论：F1 是高召回粗筛，特异性由对接承担

### 演讲备注

回应 Q1：9/9 benchmark 循环验证质疑。

---

## 第 2 页｜【备用】9-Compound Benchmark 全表

### 幻灯片要点

- SP600125 / CC-930 / E1 / TCS JNK 6O / JNK-IN-8
- CC-90001 / Q63 / AS602801 / CC-401
- 酶学 IC50 见 literature_benchmarks.csv
- 对接方向见 benchmark_deltas_51c1.csv

### 演讲备注

被问标定细节时展开。

---

## 第 3 页｜【备用】选择性探索遗留数字（§5.4）

### 幻灯片要点

- pass_selectivity：233
- Tier 1′：57；Tier 2：92；Tier 3：1191
- 均未用于 MD 短名单或采购排序
- 保留作失败探索记录

### 演讲备注

回应「233 个选择性 hit 去哪了」。

---

## 第 4 页｜【备用】§1.3 文献抑制剂证据分级

### 幻灯片要点

- A 级：氨基吡唑、YL5084、CC-930、JNK-IN-8
- B 级：E1、CC-90001、TCS JNK 6O
- C 级：SP600125（pan 对照）
- 详见报告 §1.3 与 §11 参考文献

### 演讲备注

领域背景深问时使用。

---
