# Novel 人工抽查打分清单（Step A）

**面板文件：** `novel_review_panel.csv`（50 条，seed=42）  
**来源：** `L7/L7_dock_ready_novel.csv`（3500，全部 `reactivity=ok`）  
**用途：** chemotype triage 质量抽查——**不是活性命中判定**。填完 `keep` 列后通知继续 Step B。

## 分层构成

| stratum | n | 含义 |
|---------|---|------|
| `erg_top15` | 15 | `erg_max` 最高 |
| `erg_median15` | 15 | `erg_max` 靠近全仓中位（≈0.791）随机 |
| `tc_core_boundary10` | 10 | `max_tc_core` 最靠近 0.22 下界 |
| `hinge0_erg_only10` | 10 | `hinge_hits=0`，仅靠 ErG 进 Novel |

## 填表说明

在 `novel_review_panel.csv` 填写：

- `keep`：`yes` / `no` / `unsure`
- `review_note`：可选短注（如 “染料样” / “铰链像” / “SA 高”）

## 六条打分规则（逐条过）

对每个分子，按顺序自问：

1. **激酶 ATP/铰链样？**  
   是否有可作铰链结合的杂环氮、偏平面的芳香/杂芳体系？若完全是脂肪链/糖肽样/无药效团骨架 → 倾向 `no`。

2. **丙烯酰胺单位点且几何上像能伸到 Cys？**  
   弹头应是单一丙烯酰胺；连接子是否短到离谱或埋在分子中间导致 Cys116 不可及（凭 2D 粗判即可）→ 明显不合理标 `no`/`unsure`。

3. **明显 PAINS / 染料 / 金属螯合 / 多肽样？**  
   多酚、偶氮、儿茶酚、长肽样、强螯合结构 → `no`（即使 `pains_flag=0` 也可人工判）。

4. **与 YL5084/56d 只是“看着不同但无药效团”？**  
   `max_tc_core` 低且 ErG 高，但仍无激酶样药效团 → 假 Novel，倾向 `no`。真正 scaffold-hop 应保留某种铰链/疏水药效团线索。

5. **合成难度**  
   `sa > 5` → 在 note 标 `risky`；不自动否决，但可降为 `unsure` 或 `no`（若同时有其他问题）。

6. **总体 keep？**  
   综合 1–5：值得进入后续松 Glide / AF3 配额 → `yes`；明显垃圾 → `no`；吃不准 → `unsure`。

## 建议节奏

- 先扫 `hinge0_erg_only10`（假阳性风险最高）与 `erg_top15`（是否真优）。  
- `tc_core_boundary10` 看边界是否“Sim 漏网”或合理 Novel。  
- 不必追求 50 个全 `yes`；如实 `no`/`unsure` 才能收紧 Novel。

## 完成后

把填好的 `novel_review_panel.csv` 放回同目录（或告知路径），授权执行 **Step B（收紧 watch + 重建 L7b）**。  
在此之前：**不启动 Glide / AF3**。
