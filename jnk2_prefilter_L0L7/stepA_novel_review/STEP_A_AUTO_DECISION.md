# Step A 自动 keep 判定报告（启发式）

**方法：** 按 `novel_review_checklist.md` 六条规则，用 RDKit SMARTS/物化启发式对 50 条面板自动打分。  
**输出：** `novel_review_panel_filled.csv`（含 `keep` / `review_note`）；审计版 `novel_review_panel_auto_scored.csv`（含 r1–r5 布尔列）。

> 注：原 `novel_review_panel.csv` 若被 Excel 占用可能未覆盖，请以 `novel_review_panel_filled.csv` 为准。

## 面板 50 条结果（seed=42）

| keep | n | 比例 |
|------|---|------|
| **yes** | 31 | 62% |
| **unsure** | 16 | 32% |
| **no** | 3 | 6% |

### 分层

| stratum | yes | unsure | no |
|---------|-----|--------|-----|
| erg_top15 | 9 | 6 | 0 |
| erg_median15 | 10 | 5 | 0 |
| tc_core_boundary10 | 7 | 3 | 0 |
| hinge0_erg_only10 | 5 | 2 | **3** |

**解读：** 假阳性风险最高的 `hinge0_erg_only10` 中 3/10 被判 `no`（无铰链药效团 / PAINS / 儿茶酚）；高 ErG 顶层并非全过，说明 ErG 单独不足以保 Novel。

### 自动否决的 3 条（no）

| id | 原因 |
|----|------|
| CHEMBL517121_ACR0 | 无激酶药效团 + 假 Novel（低 core-Tc 高 ErG） |
| CHEMBL188782_ACR0 | 无铰链样 + PAINS |
| CHEMBL1309819_ACR0 | 儿茶酚样 bad motif |

## 外推到全仓 3500 Novel（L7）

若将同一启发式用于全部 `L7_dock_ready_novel.csv`：

| keep | n |
|------|---|
| yes | **2,246** (64.2%) |
| unsure | 940 (26.9%) |
| no | 314 (9.0%) |

## 是否进入 Step B？

**建议：可以进入 Step B**，但 Novel 仓采用**保守策略**：

1. **Novel L7b** 仅保留启发式 **`keep=yes`**（约 **2,246**），**不纳入 unsure**；不足 3,500 配额则如实报告（符合 ARS「禁止回填 discard」）。
2. **sim_yl / sim_56d / pan** 的 watch 收紧与 Novel 判定独立，按 Step B 原计划执行。
3. **unsure（940）** 不进入对接主名单；可另存 `novel_unsure_hold.csv` 供你后续人工翻案。
4. 本判定为 **chemotype triage**，**不是活性命中**；对接后仍以 AF3 mPAE 门控。

若你同意，下一步执行 Step B + Step C（仍不跑 Glide/AF3，除非另授权）。

## 启发式规则摘要（可复现）

- **R1 激酶/铰链样：** `hinge_hits≥1` 或 hinge SMARTS 命中，或 ≥2 芳香环且含芳香杂原子  
- **R2 丙烯酰胺几何：** 单 acrylamide；酰胺 N 到芳香体系最短路径 ≤8；RB≤12  
- **R3 排除：** PAINS catalog、偶氮、儿茶酚、肽样多酰胺、≥3 酚 OH  
- **R4 假 Novel：** `max_tc_core<0.15` 且 `erg≥0.78` 且无 R1 → no  
- **R5 SA>5：** 标 risky → 降为 unsure（无其他问题时）  
- **R6 综合：** 见 `novel_review_panel_auto_scored.csv`
