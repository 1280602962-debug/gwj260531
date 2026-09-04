# 路线 C2：以"转运循环阻断"重建课题

> 本文是**课题重建方案 + 预登记草案**。机器可读锁见 [`config/campaign_c2.yaml`](../config/campaign_c2.yaml)。
>
> 本文**不覆盖**任何冻结数据：`data/repurposing/p2/` 保持原样；C1 的 `pass_fail.json`、A1/A2 产物保持原样。
> C1 Rank 轨仍为 closed；C1 Acid 短名单（PF-04620110 / Admilparant / PSI-697 / PF-03882845 / Lanifibranor）在本路线中**降级为输入池中的先验假设**，不是产品。

---

## 0. 对上一轮分析的裁定

上一轮结论中**成立、可写进论文**的部分：

| 结论 | 状态 | 依据 |
|---|---|---|
| 羧酸根重对接不能作为活性排序器 | **成立**，且本文给出更强的机制性理由（见 §2） | A2 回顾 OR = 0.970、Fisher p = 1.0、sens 0.952、spec 0.047（`acid_gate_benchmark_summary.json`） |
| RandomDecoy 构造在 Gu 框架内自洽 | **成立** | 6 万包络剩余随机抽 4690；True∩Random = 0；骨架重叠 0；无 TC > 0.5 |
| 用 RandomDecoy 单独否决"临床库可排序" | **过强，需修正** | P2: True EF@1% = 2.587 (12/51, p = 0.0016) vs Random EF@1% = 0.215 (1/51, p = 0.992)；P0 两侧均 ≈1.94 |
| 单构象打分排不出 IC50 | **成立，但上一轮给的理由不够硬** | 上一轮只说"任务误指定"；§2 把它升级为可证伪的结构论证 |
| 换 KarmaDock/CarsiDock/DiffDock 无用 | **成立** | LIT-PCBA 2026 中位 EF@1%：AutoDock-GNINA 2.14、DiffDock 0.84、ML 重排 4.49 |

需要**撤回或改写**的部分：

1. 上一轮建议的"双态差分打分"（inward vs outward 打分之差）**不新颖**。Costanzi & Vilar (*J. Comput. Chem.* 2011, 33:561) 已用 β2AR active/inactive 对接分数之差做 agonist/blocker 判别（ROC + LDA）；Männel & Kolb (*Mol. Pharmacol.* 2019) 的 "selective reranking" 同理。更致命的是 **HOLIgraph**（*J. Cheminform.* 2025, 17:1020）已经在 SLC 上做完了：把 222 个配体对接到 OATP1B1 的 inward / outward cryo-EM 结构，用相互作用图 GNN 判别 inhibitor / non-inhibitor，平衡准确率 > 90%。若按上一轮写法投出去，会被判为 HOLIgraph 的换靶复现。
2. 上一轮把 conformal FDR 选择列为主要创新之一，**在本课题的样本量下不成立**。可用于校准的"机制标注"化合物只有个位数（§5），n ≈ 8 时 BH-FDR 没有实际功效。conformal 只保留在生成阶段（N 大）作为可选项，不作为卖点。

因此本文换掉的不是打分函数，是**判别对象**。

---

## 1. 一句话路线

> **不再问"谁与口袋结合更紧"，而问"谁与转运循环不相容"。**
> 抑制剂的定义性特征不是亲和力，而是**占据尿酸位点后无法沿 outward → occluded → inward 的异构化路径被带过去**。
> 据此建立一个**对比型**（而非绝对型）可观测量，用已有 cryo-EM 药物阶梯 + 被转运阴离子做正/负校准，再用它在临床库上提名候选，最后用本地 MD 做机制确证。

---

## 2. 为什么必须换判别对象（本课题的核心论证）

2024–2025 的四组 cryo-EM 把 URAT1 的机制钉死了：

- **Dai & Lee**, *Cell Res.* 2024, 34:776–787（9B1F–9B1O）：捕获 outward → inward 转变；尿酸结合在富 Phe 口袋并"engage with key gating residues to **drive the transport cycle**"；反向阴离子**吡嗪酸占据三个不同的功能性位点**；三个药物都"compete with substrates and **halt the transport cycle**"，其中 verinurad 与 dotinurad 进一步"**hijack gating residues** to achieve high potency"。
- **Guo/Chen**, *Nat. Commun.* 2025, 16:1512（9IRW/9IRX/9IRY）：在 NTD 对齐下，**benzbromarone 与 outward-facing 态的 F449 侧链冲突，verinurad 与 R447 冲突**；结论是抑制机制为"blocking not only the binding of urate but also **the structural isomerization** of hURAT1"。
- **Cell Discovery** 2025, 11（9JE1 等）：四个药物均在 inward-open；尿酸在 inward-open (3.3 Å) / outward-open (4.1 Å) / occluded (4.7 Å) 三态。
- **Nat. Commun.** 2025, 16 (60480-3 / 61226-x)：抑制剂选择性结合 inward-open 态；不同药物的构象与相互作用网络差异决定特异性。

由此得到两条对本课题致命的推论：

### 2.1 "阴离子/羧酸识别"是被转运的签名，不是抑制的签名

URAT1 是**尿酸 / 单羧酸反向转运体**：它的生理职责就是结合并释放阴离子（尿酸进、乳酸/烟酸/吡嗪酸出）。因此任何"能被富 Phe 腔与电正性残基（Arg477 一侧）识别的阴离子"这一特征，被**底物、反向阴离子、抑制剂三类同时共享**。

- 吡嗪酸是**抗尿酸排泄剂**（升高尿酸重吸收），却在 URAT1 上占据三个位点。
- 尿酸本身与抑制剂占同一中心位点。

所以以"羧酸根 + Arg 距离 + 在腔内"为门的 A2，其判别力上界是被结构决定的。实测 **OR = 0.970、p = 1.0、sens 0.952、spec 0.047** 正是这个上界的体现——不是门开得太松，而是**这个特征在物理上不含判别信息**。这条论证把上一轮"A2 没新意"升级成"A2 不可能有判别力"，是可发表的负结果。

> **附带的药理危险**：以酸性/阴离子相似性为主的筛选，其假阳性的一部分恰恰是**吡嗪酸样反向阴离子**，即潜在的**抗尿酸排泄**分子——方向与治疗目标相反。任何只报对接分数的流程都不检查这一项。

### 2.2 单构象打分与"抑制"之间不存在单调关系

设 \(S_{\text{in}}(L)\) 为配体 \(L\) 在 inward-open 态的任意占据型打分。上述结构事实意味着：

\[
S_{\text{in}}(\text{urate}) \sim S_{\text{in}}(\text{pyrazinoate}) \sim S_{\text{in}}(\text{inhibitor})
\]

三类分子在同一态、同一位点都"打得好"。区分它们的量是**跨态的**：

\[
\Delta(L) \;=\; \big[\text{与 outward / occluded 态的不相容度}\big] \;-\; \big[\text{与 inward 态的不相容度}\big]
\]

底物与反向阴离子 \(\Delta \approx 0\)（三态皆可容纳，故被带过去）；抑制剂 \(\Delta \gg 0\)（占位后卡死异构化）。

这同时解释了 P2 的失败模式：\(S_{\text{in}}\) 中与分子大小、疏水体积、CNN 训练偏倚相关的成分是**共模误差**，在 TrueDecoy（物化匹配）上被部分抵消所以 EF = 2.59，在包络随机集上暴露所以 EF = 0.215。而 \(\Delta\) 是**同一分子、同一姿态在两个受体态之间的差**，共模成分被结构性地抵消。**共模抑制**是本路线相对"再换一个打分函数"的物理理由，不是修辞。

---

## 3. 可观测量：三层，从确定性到增强采样

对每个分子，先在 inward-open 态得到姿态（这一步只要求**召回**，不做排序声称），然后计算下列量。

### T1 — 跨态姿态移植的立体不相容度（确定性，无 MD，可跑全库）

1. 以 **NTD（TM1–6）Cα** 将 outward-open / occluded 结构叠合到 inward-open 结构；
2. 将 inward 姿态的配体坐标**原样移植**进被叠合的目标态；
3. 计算重原子立体重叠：

\[
\Delta_{\text{state}} \;=\; \frac{1}{N_{\text{heavy}}}\sum_{i \in L}\sum_{j \in P} \max\!\big(0,\; r_i + r_j - d_{ij}\big)
\]

4. **以 CTD（TM7–12）叠合重复一次**；只有在两种叠合下都判为不相容，才计为 outward-incompatible（保守规则）。

设计约束（针对 P2 的失败教训，必须预登记）：
- \(\Delta_{\text{state}}\) 按重原子数归一，且必须报告它与 MW / 重原子数 / logP 的相关；若 \(|\rho| > 0.4\)，改用对 apo 参考归一的 z-score 或做偏相关校正后再判读。
- 主指标是**对比量** \(\Delta_{\text{out}} - \Delta_{\text{in}}\)，不是 \(\Delta_{\text{out}}\) 本身。

**T1 的内建真值**：必须复现 Guo/Chen 已发表的两处冲突——benzbromarone↔F449、verinurad↔R447。复现不了就说明叠合/移植实现有错，不许继续。

### T2 — 态特异复合物稳定性（常规膜 MD）

在 inward 与（分辨率允许范围内的）outward/occluded 态分别跑膜 MD，读出：
- 配体是否被排出、位移 RMSD；
- 门控残基接触占据率（Arg477、Phe365、Ser35、F449、R447、5-Phe 笼）；
- 胞外/胞内门距离分布是否被配体钉住。

力场按已发表 URAT1 MD 配方：**Amber14SB（蛋白+脂）/ GAFF2（配体）/ POPC / TIP3P，PME，vdW cutoff 9 Å，2 fs，Parrinello–Rahman**，生产段 **3 × 150 ns**（*Nat. Commun.* 2025, 16:61226-x 用的就是这一套）。

### T3 — 异构化惩罚（增强采样，本路线的方法学产品）

沿 **NTD–CTD 相对刚体转动 / 门距离** 定义 CV，比较 apo、+urate、+抑制剂：
- 首选低成本路线：多副本慢生长 SMD + Jarzynski 平均，报**相对功** \(W_{\text{iso}}\)；
- 若收敛允许，再做窗口伞形采样得 PMF，报**相对势垒移动** \(\Delta\Delta G^{\ddagger}_{\text{iso}}\)。

**只报相对量，不报绝对结合自由能。**不做 MM-GBSA 重排。

> **分辨率红线**：Cell Discovery 的 outward-open (4.1 Å) 与 occluded (4.7 Å) 且 **ECD 缺失**；Dai & Lee 的 occluded 9B1K 为 3.3 Å。因此：MD 与 CV **优先用 9B1K/9B1L**；任何**侧链级**结论不得建立在 4.1/4.7 Å 模型上；T1 在低分辨率态上只允许做骨架/叶片级判读，并须同时报告在 9B1K/9B1L 上的结果。

---

## 4. 先例与新颖性边界（必须在论文里自己写清楚）

| 先例 | 它做了什么 | 本路线的差异 |
|---|---|---|
| Costanzi & Vilar 2011；Männel & Kolb 2019 | GPCR active/inactive **对接分数之差**判别 agonist/antagonist | 我们用的不是两个静态分数之差，而是**同一姿态跨态移植的立体不相容**与**异构化路径惩罚** |
| **HOLIgraph** 2025（OATP1B1） | 对接到 inward/outward 双态，用相互作用图 GNN 判别 **inhibitor vs non-inhibitor**，BA > 90% | 判别对象不同：我们判别 **inhibitor vs 被转运配体（底物/反向阴离子）**。双态在 HOLIgraph 中是**特征来源**，在本文中是**机制对比的两端** |
| NIS *JCIM* 2025 | 多离子占据态对接 + ML | 同上：仍是 inhibitor/non-inhibitor，且无路径量 |
| MFS 外排泵抑制剂建模 | 用"稳定某构象"作定性叙述 | 我们把"锁定"变成**预登记的定量判据**，并用 cryo-EM 药物阶梯校准 |
| Love 等 2011；Yang *RSC Adv.* 2023；Du 2024（URAT1） | 对接出姿态 + 药效团/IFD + MM-GBSA → HEK | 我们不把分数当亲和力，也不再做一遍 e-pharmacophore 漏斗 |

**可主张的新颖点（三条，缺一条则降级为应用文）**
1. 在抗转运体上，把"抑制"重新形式化为**转运循环不相容**，并给出**共模抑制**的对比型可观测量；
2. 用**被转运阴离子（尿酸、吡嗪酸、烟酸、乳酸）作为机制匹配负控制**——这是任何 inhibitor/non-inhibitor 数据集都不具备的负例本体；
3. 证明**阴离子/羧酸门在结构上不可判别**（用自己的 OR ≈ 0.97 作实证），并用同一批姿态做 T1 的**同场对照**。

---

## 5. 校准集（正/负控制）

### 5.1 主阶梯（同一实验室、同一 ¹⁴C-尿酸摄取测定，均有人源 cryo-EM 共结构）

| 化合物 | IC50 | 态 | PDB |
|---|---|---|---|
| dotinurad | ~8 nM | inward-open | 9JE1 |
| verinurad | 40 nM | inward-open | 9IRY / Cell Discovery 系列 |
| benzbromarone | ~200 nM | inward-open | 9IRX / Cell Discovery 系列 |
| lesinurad | ~12 µM | inward-open | 9DKB（本课题生产态）/ 9JDZ |

来源：*Cell Discovery* 2025, 11（8 nM / 40 nM / ~200 nM / ~12 µM）。跨度约 3.2 log，**同测定**，是本课题唯一可用的自洽阶梯。

### 5.2 扩展阶梯（跨测定，须标注测定差异）

Taniguchi 等 *JPET* 2019, 371:162（MDCKII-URAT1）：dotinurad **0.0372 µM**、benzbromarone **0.190 µM**、lesinurad **30.0 µM**、probenecid **165 µM**。probenecid 仅从此源进入，且本课题 A2 回顾中 probenecid 一直不过门——它是最有价值的困难正例。

### 5.3 机制匹配负控制（本路线的关键资产）

| 分子 | 身份 | 为什么是负例 |
|---|---|---|
| urate | 底物 | 结合同一中心位点，三态均可容纳（9B1K/9B1L 及 Cell Discovery 三态） |
| pyrazinoate | 反向阴离子，**抗尿酸排泄** | 占据三个功能位点却不是抑制剂（Dai & Lee 2024） |
| nicotinate | 反向阴离子 | 同上类别 |
| lactate | 生理单羧酸反向底物 | 同上类别 |

**判据**：任何声称有效的可观测量，必须把 5.1/5.2 的抑制剂与 5.3 的被转运配体分开。做不到 → 该可观测量作废。这是本路线最硬的一道闸，也是它与所有 inhibitor/non-inhibitor 工作的分界。

### 5.4 有功效的大样本对照（用已有姿态，几乎零新增算力）

T1 足够便宜，因此直接跑在**已有的 phase-1 gnina 姿态**上：469 个 p≥6 活性 + 80 个实验弱活（p<5）+ TrueDecoy / RandomDecoy 诱饵。产出与 A2 **同场可比**的 OR / AUC：

> A2 酸门：OR = 0.970，p = 1.0 ⟶ T1 状态门：OR = ?，p = ?

这是本课题的**决胜实验**，且现在就能做（见 S1）。必须按 Murcko 骨架分层报告（Top-1 骨架占 127/469 = 27%），主指标为**骨架簇内**统计量。

---

## 6. 实施阶段

> 顺序不可打乱：**对照先行**。任何候选体系在对照未过关前不得开跑。

### S0 — 结构清单核验（本地，无算力）
- 逐条核验并落盘：9DKB、9B1K（occluded, 3.3 Å, urate）、9B1L（outward, urate）、9IRW/9IRX/9IRY、9JE1、9JDZ；NLRP3：7ALV、8EJ4、7PZC、8ETR。
- **本仓库已有过 9JDZ 误标事故**，因此每个 accession 必须记录标题、态、配体、分辨率，不得凭记忆使用。
- 产物：`data/campaigns/c2/00_preregistration/structure_manifest.csv`

### S1 — T1 实现与决胜回顾（关键路径，先做这一步）
1. 实现 NTD/CTD 双叠合 + 姿态移植 + \(\Delta_{\text{state}}\)；
2. **真值检查**：复现 benzbromarone↔F449、verinurad↔R447；
3. 5.3 负控制检查：urate/pyrazinoate/nicotinate/lactate 必须 \(\Delta \approx 0\)；
4. 5.4 大样本回顾：469 + 80 + 诱饵，骨架分层 OR/AUC，与 A2 同表并列；
5. 尺寸共模检查：\(\Delta_{\text{state}}\) vs MW/重原子/logP 相关系数。
- **门**：2、3 任一失败 → 停止，不进入 S2。4 若 OR 的 95% CI 覆盖 1 → 路线判负，按 §7.2 出稿。

### S2 — 生成阶段（召回优先，不做排序声称）
- 池：临床/再定位池（**不预先按酸性过滤**——§2.1 已说明酸性不是抑制签名）；
- 对 inward 态（9DKB + 9JE1）多种子对接，保留占据尿酸位点者；
- 记录召回而非精度；C1 A2 的 40 个 eligible 与竞争短名单作为**先验子集**一并带入，接受同样裁决。

### S3 — T1 状态门筛选
- 对 S2 存活者算 \(\Delta_{\text{out}}\)、\(\Delta_{\text{occ}}\)，按 S1 校准出的阈值分三类：
  **inhibitor-like**（inward 占位 + 跨态不相容）/ **substrate-like**（三态皆容纳，**标记为抗尿酸排泄风险**）/ **non-binder**。
- 产物：候选短名单（预期 20–60），**substrate-like 必须单独列表并显式警示**，不得混入候选。

### S4 — MD 确证（本地，对照先行）
必跑对照（按顺序）：
1. apo 9DKB 膜体系
2. urate @ inward（底物对照）
3. lesinurad 晶体姿 @ 9DKB（弱抑制剂对照）
4. dotinurad @ 9JE1（强抑制剂对照）
5. pyrazinoate @ inward（反向阴离子对照）
6. NP3-146 @ 7ALV（NLRP3 结构对照）

对照通过后再跑候选，每体系 3 × 150 ns。
- **预写死判读**：体系 3/4 的门控残基占据未复现晶体接触 → URAT1 侧全部不解释；体系 5 若与抑制剂无法区分 → T2 判负，退回只用 T1。

### S5 — T3 异构化惩罚（方法学产品）
- 先在 apo / +urate / +dotinurad / +lesinurad 上建立 CV 与 \(W_{\text{iso}}\)；
- **验收**：\(W_{\text{iso}}\) 必须把 urate 与 dotinurad 分开，且 dotinurad > lesinurad（与 3 log 的 IC50 差同向）；
- 通过后才算候选。**n = 4 无法支撑相关性主张**，因此主张限于"分离 + 序一致"，并显式声明样本量限制。

---

## 7. 预登记通过 / 失败标准

### 7.1 通过（可提名候选）
- T1 复现两处已发表冲突；
- T1 把被转运阴离子与抑制剂分开；
- S1 大样本回顾中 T1 状态门的骨架分层 OR 显著 > 1，且显著优于 A2 酸门（OR ≈ 0.97）；
- \(\Delta_{\text{state}}\) 与分子尺寸无强相关（或校正后仍显著）；
- S4 对照体系全部过关。

### 7.2 失败（仍然出稿，方向不同）
若 S1 判负，本课题的产品变为**方法学负结果 + 结构性解释**：
> 在 SLC22 反向转运体上，(i) 单构象占据型打分与抑制无单调关系；(ii) 阴离子/羧酸几何门在结构上不可判别（OR ≈ 0.97）；(iii) 即便把双态不相容形式化，也不足以判别——并给出三态分辨率与姿态误差的定量归因。
配合已有的 TrueDecoy/RandomDecoy 分解、lesinurad 自对接失败（4.30/4.31/4.88 Å）、A2 OR ≈ 0.97、临床库负迁移，这是一篇完整的方法学论文。**两种结局都可发表**，这是本路线相对"再筛一轮"的最大优势。

---

## 8. NLRP3 对称臂：同一个物理主张

NLRP3 侧的机制与 URAT1 **同构**，这让"双节点"第一次有了统一的物理论点，而不是两个拼在一起的单靶故事。

- 7ALV：NP3-146 与 ADP 共结晶，抑制剂充当 **intramolecular glue**，把 NBD/HD1/WHD/HD2 四个亚域锁在**非活性闭合态**；
- 活化需要 NBD 相对其余亚域**转动约 85–90°**（ATP 态活性寡聚体 8EJ4）；MCC950 稳定闭合 ADP 态（7PZC 闭合十聚体）。

因此 NLRP3 抑制同样是**构象锁定**而非单纯占位，可用同一套 T1/T3 逻辑：以闭合态（7ALV）姿态移植进活性态（8EJ4）NACHT，度量不相容；对照为 NP3-146（自对接已达 **0.82 / 0.67 / 0.68 Å**，是本课题最可靠的结构对照）。

论文主标题因此可以写成：**构象锁定型虚拟筛选（conformational-lock virtual screening）**——一个跨两个完全不同折叠的方法主张，而不是"痛风双靶筛选"。

---

## 9. 与已有资产的关系

| 资产 | 在 C2 中的角色 |
|---|---|
| 冻结 P2 表（`data/repurposing/p2/`） | 负迁移基线，**只读** |
| TrueDecoy / RandomDecoy | 生成阶段对照 + 诱饵本体分解（方法学章节） |
| phase-1 gnina 姿态 SDF | **S1 决胜回顾的输入**，零新增对接 |
| C1 `pass_fail.json`、A1/A2 | 只读历史；A2 的 OR ≈ 0.97 升级为 §2.1 的实证 |
| C1 A2 短名单（40 eligible + 5 竞争） | 输入池中的先验子集，接受 S3 同等裁决 |
| NLRP3 ML（AUROC 0.893） | 生成阶段先验，不作排序 |
| URAT1 配体 ML（OOF R² 0.53，命名药 2/4） | 已判失败，不复活 |

## 10. 禁止事项（沿用并扩展）

- 不覆盖 `data/repurposing/p2/`，不重锁 Π\*；
- 不用诱饵标签训练任何用于协议选择的模型；
- 不把 \(\Delta_{\text{state}}\)、\(W_{\text{iso}}\) 称作亲和力或 IC50；
- 不做 MM-GBSA 重排，不比"谁更像 lesinurad"；
- 不在 4.1 / 4.7 Å 模型上做侧链级结论；
- 不把 substrate-like 分子写进候选（须作为抗尿酸排泄风险单列）；
- 不声称双靶抑制剂——产品是**待实验验证的双节点构象锁定假说**；
- 不在 n = 4 的阶梯上声称相关性。

## 11. 风险与缓解

| 风险 | 缓解 |
|---|---|
| outward/occluded 分辨率低、ECD 缺失 | 优先 9B1K/9B1L；低分辨率态只做叶片级判读；双叠合保守规则 |
| 姿态误差传播进 T1（生产姿 Arg477 曾达 14 Å） | T1 输入姿态须先过 inward 占位检查；对 5.1 的 4 个药物用**晶体姿**而非对接姿建立阈值 |
| \(\Delta_{\text{state}}\) 仍被尺寸驱动 | 归一 + 相关性检查 + 偏相关校正，预登记阈值 \(|\rho| \le 0.4\) |
| T3 收敛困难 | 分层：SMD/Jarzynski 相对功优先，PMF 为可选；不收敛则只报 T1/T2 |
| 校准样本量小 | 主张限于"分离 + 序一致"；大样本功效由 S1 的 469 + 80 + 诱饵提供 |
| 与 HOLIgraph 重叠被质疑 | §4 表格自陈差异；判别对象为 inhibitor vs 被转运配体，且含路径量 |

## 12. 参考文献锚点

- Dai Y, Lee CH. *Cell Res.* 2024, 34:776–787（9B1F–9B1O；三态、吡嗪酸三位点、gating hijack）
- Guo W 等. *Nat. Commun.* 2025, 16:1512（9IRW/9IRX/9IRY；F449 / R447 跨态冲突；阻断 structural isomerization）
- *Cell Discovery* 2025, 11（9JE1 等；四药 inward-open + 尿酸三态；IC50 8 nM / 40 nM / ~200 nM / ~12 µM）
- *Nat. Commun.* 2025, 16:60480-3 与 16:61226-x（inward-open 选择性；Amber14SB/GAFF2/POPC/TIP3P，3 × 150 ns）
- Taniguchi T 等. *JPET* 2019, 371:162（dotinurad 0.0372 / benzbromarone 0.190 / lesinurad 30.0 / probenecid 165 µM）
- Cai 等. *Sci. Rep.* 2017, 7:706（verinurad；依赖 Phe365 与 Ser35）
- Dekker C 等（7ALV；intramolecular glue，锁定非活性闭合态）
- **先例声明**：Costanzi & Vilar *J. Comput. Chem.* 2011, 33:561；Männel & Kolb *Mol. Pharmacol.* 2019；HOLIgraph *J. Cheminform.* 2025, 17:1020；NIS *JCIM* 2025
- 评价方法学：Wallach & Heifets *JCIM* 2018（AVE）；Chen 等 *PLoS One* 2019（DUD-E 偏倚）；Bayes 富集 arXiv:2403.10478；LIT-PCBA 泄漏 arXiv:2507.21404；LIT-PCBA 基准 arXiv:2605.01681；Gu 等 *Nat. Mach. Intell.* 2025
