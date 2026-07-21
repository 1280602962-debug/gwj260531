# 中文先行大纲（V2 配套写作稿，2026-07-21）

> 配合 `MANUSCRIPT_OUTLINE_V2.md`。标 ✅ 者数据已就位、现在即可动笔；标 ⏳ 者需等 URAT1 协议筛选（P0–P5 + RTMScore）跑完。
> 当前对接进度：URAT1 Vina ~70%，gnina 待跑，RTMScore 旁路进行中 → **协议 Π\* 尚未锁定**。

---

## 0. 现在能写 / 需等的分工

| 板块 | 状态 | 说明 |
|------|------|------|
| NLRP3 模型性能与缩库 | ✅ | 数字齐全 |
| 不对称设计论证（URAT1 ML 不足 vs NLRP3 ML 可用） | ✅ | benchmark 回收对比已有 |
| TrueDecoy/RandomDecoy 方法框架（Methods 写法） | ✅ | 集已建好，规则已定 |
| 期刊定位、claim 边界、Differentiation | ✅ | 见 V2 + paper_spine 分析 |
| URAT1 协议筛选结果表 / Π\* 选定 | ⏳ | 等对接跑完 |
| redock 双指标最终数值 | ⏳ | lesinurad top-1/best/RTMScore |
| Pareto/提名用 Π\* 终版重算 | ⏳ | 协议定后 |
| MD 构象讨论 | ⏳ | 视算力 |

---

## 1. 标题（暂定）

痛风 URAT1–NLRP3 双节点临床药物重定位：基于 TrueDecoy/RandomDecoy 的 URAT1 对接协议筛选与不对称 NLRP3 机器学习漏斗

英文见 V2 §3.3。

---

## 2. 摘要（五句，先写 NLRP3 与框架句，URAT1 数字留空）

1. 痛风涉及 URAT1（尿酸重吸收）与 NLRP3（炎症）双节点，临床库缺乏尊重两靶数据不对称、且对接协议经系统筛选的可复现流程。
2. 我们先在 URAT1 inward-open 结构（9DKB）上，用性质匹配 TrueDecoy 与随机 RandomDecoy 两套基准，比较开源对接与重打分协议（AutoDock Vina、gnina、RTMScore），并诚实报告自对接排序表现。【⏳ 数字待填】
3. 随后以 NLRP3 assay-conditioned 分类模型（骨架 CV **AUROC 0.893 / AUPRC 0.914**）缩小 8,319 临床化合物至 1,588，再在双靶上并行对接并用 Pareto 取非支配前沿。
4. 漏斗回收已知对照、将 PAINS 型盲筛命中（EGCG）主动降级，并提名机制更清晰的候选（canagliflozin）作为可检验假说。
5. 本工作输出可证伪的计算假说与方法学负结果，非双靶抑制剂发现或临床推荐。

---

## 3. 引言（四段）

1. **疾病双轴**：URAT1 降尿酸 vs NLRP3 抗炎；临床常分轴用药，双节点干预有理论吸引力。
2. **计算难点**：转运体 cryo-EM 虽新（9DKB），但对接富集 ≠ 结合位点证明；RandomDecoy 易虚高；活性标签未必都作用于所选口袋 → 需要 TrueDecoy 与诚实 redock。
3. **相关工作与区分**：Gu/Hou VSDS-VD 双诱饵范式（借其评估逻辑，工具裁剪为开源）；湿法双靶抗痛风（不同范式，互补）；PLK1/NLRP3 不对称 VS（**必须写清区别，防"换靶增量"批评**）。
4. **本文贡献**：C1 URAT1 对接协议筛选；C2 不对称双节点漏斗；C3 Pareto≠提名的成药性审计；C4 canagliflozin 假说。claim 边界明确。

---

## 4. NLRP3 回顾（本次重点，✅ 可全写）

### 4.1 数据与模型
- 任务：assay-conditioned 二分类（多 assay 集成，见 `nlrp3_screening_summary`）。
- 交叉验证（骨架 GroupKFold, 5 折）：**AUROC 0.893，AUPRC 0.914，EF@10%≈1.57**，判定 screening-ready。
- benchmark 回收：**MCC950 P=1.0；GDC-2394 P=0.917 → 2/2 通过**（均在训练集，作 sanity，不宣称外推）。

### 4.2 缩库行为（8319 → 1588，阈值 0.5）
- 对照药百分位（NLRP3 ML）：
  - verinurad 90.4、colchicine 90.4（高：结构相似/间接调节，符合预期，也是"NLRP3 ML 不宜单独定夺"的警示）
  - lesinurad 81.8（中）
  - benzbromarone / dotinurad / allopurinol / febuxostat ≈ 35.4（低）
- 解读：分类器能把"炎症相关"从纯降尿酸药中分层；colchicine 高分正是**为何还需对接与 Pareto 双轴、避免单一 ML 假阳性**的论据。

### 4.3 NLRP3 对接的定位（回答"用什么对接"）
- 7ALV 配体为 **MCC950 类似物 NP3-146（非 MCC950 本体）**，属药效团模板 → NLRP3 pose 可信度天然弱于 URAT1。
- 因此 NLRP3 对接是 **双轴中的结构佐证**，非主证据；主轴是 ML。
- 协议：**直接沿用 URAT1 选定的 Π\***（同引擎/同打分）@ 7ALV，不单独做完整 True/Random 选拔（如审稿要求可补 NLRP3 专属诱饵集）。
- S_N 取 max(NLRP3 ML 百分位, 7ALV 对接百分位)，理由：Spearman(ML, 7ALV)≈ −0.04，两者近正交，取 max 避免单轴埋没。

---

## 5. URAT1 对接协议筛选（⏳ 等跑完，先写 Methods 与占位）

### 5.1 为什么 URAT1 要做协议筛选而 NLRP3 不做
- URAT1 benchmark 回收仅 **2/4**（lesinurad、dotinurad 失败），ML 不可主排序 → 对接必须靠得住 → 必须先选协议。
- NLRP3 ML 已 2/2、AUROC 0.89 → 有可靠主轴，对接仅佐证。

### 5.2 基准集（✅ 已建）
- `true_decoy_benchmark.csv`：469 actives（pActivity≥6）+ 6073 性质匹配 TrueDecoy
- `random_decoy_benchmark.csv`：同 actives + 6073 RandomDecoy
- 受体 9DKB，统一搜索盒。

### 5.3 候选协议（✅ 已锁定设计，结果 ⏳）
P1 Vina affinity｜P2 gnina CNNaffinity｜P3 gnina affinity(kcal)｜P4 Vina+RTMScore｜P5 gnina+RTMScore｜P0 gnina CNNscore(负对照)｜P6–P8 Glide±RTMScore(有许可)。

### 5.4 选优规则（✅ 预先锁定）
主判 TrueDecoy EF@1%（并列 EF@5%、AUC）；否决 RandomDecoy 明显变差；平局四药回收百分位。→ 输出 Π\*（可含分层：快引擎粗排 + RTMScore 精排）。

### 5.5 redock 门控（✅ 写法定，数值 ⏳）
- 报告三值：Top-1 RMSD、Best-in-ensemble RMSD、RTMScore 选姿 RMSD。
- 已知现象（诚实写）：Vina/gnina 原生 **top-1 未过 2 Å**（≈4.5–5 Å），但 ensemble 内存在近晶体姿（≈0.84 Å），说明**采样可达、原生排序不可靠** → 富集结论限于"排序能力"，pose 叙事改用 RTMScore/晶体姿。

---

## 6. 不对称漏斗（应用 Π\*，⏳ 终值待协议定）

1. 8319 →（NLRP3 ML P≥0.5）→ 1588 ✅
2. Π\* @ 9DKB + Π\* @ 7ALV ⏳
3. Pareto(S_U, S_N) → 短名单 ⏳（现有 Glide 开发跑：1451→Pareto 6，仅作占位）
4. 模块 A–F：过滤/适用域/稳健性 → EGCG 降级、canagliflozin 提名 ✅逻辑，⏳终表

---

## 7. 构象讨论 / MD（⏳ 视算力，非发现证明）
- 体系：benz/dot@9DKB；MCC950@7ALV；canagliflozin@9DKB+7ALV。
- 初始姿用晶体或 RTMScore struct_pose，不用失败 top-1。
- 时长 50–100 ns（JCAMD 可接受，不冲 JMM）。

---

## 8. 讨论
1. TrueDecoy 相对旧 A vs D 如何改变协议选择。
2. redock 排序失败对 pose 叙事的限制；富集仍指导排序。
3. 两靶差异：为何 URAT1 对接主导、NLRP3 ML 主导——同一套对接协议、不同证据权重。
4. EGCG vs canagliflozin：审计如何改提名。
5. 机制边界与证伪路径（URAT1 摄取实验 / MSU–IL-1β）。
6. 局限：诱饵池规模、开源 vs Glide、单态口袋、7ALV 非 MCC950 共晶、无湿实验。

---

## 9. 结论
先按双诱饵框架选定 URAT1 对接排序协议，再嵌入不对称 NLRP3 漏斗，产出可审计、可证伪的双节点重定位假说；非双靶药发现。

---

## 10. 图表
- Fig1 总流程；Fig2 True vs Random 各协议 EF/AUC ⏳；Fig3 redock RMSD 三值 ⏳；Fig4 漏斗与对照 ✅/⏳；Fig5 Pareto 与 Module F（EGCG 降级）；Fig6 MD ⏳。
- 表1 协议定义与选优结果 ⏳；表2 Pareto 6；表3 提名 Top；表S redock 数值。

---

## 11. 现在就能落笔的写作顺序（建议）
1. §4 NLRP3 全节（数字齐）
2. §5.1–5.4 URAT1 协议筛选的 Methods 与 rationale（除结果表）
3. §3 引言（含 Differentiation）
4. §5.5 redock 写法 + §8 讨论骨架
5. 摘要框架句
→ 待 URAT1 跑完，回填 §5 结果表、Π\*、§6 终值、Fig2/3。
