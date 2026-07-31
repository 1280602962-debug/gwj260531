# JCIM 补充实验方案 v2 —— 外部 holdout + 结构稳健性 + 剩余厚度

> 日期：2026-07-31
> 前提：环境已验证可实际对接（AutoDock Vina 1.2.7 python 绑定 + meeko 0.7.1 + RDKit 均已装好，ChEMBL/RCSB 网络可达）；`CLAIM_CEILING.md` 与 `PAIR_ROLES_APPROVED_JCIM.yaml` 的冻结协议不变。
> 关系：本文件是对 `JCIM_STRENGTHENING_PLAN_V1.md`（Wave 0-4 总规划）第 3 步（文献 holdout）与 Wave 2.4/2.5（结构稳健）的**可执行落地版**，并回应「和近三年 JCIM 比还差什么」的复核结论。

---

## 0. 先回答：上一轮说法哪些对、哪些要改

| 上轮说法 | 复核结论 | 依据 |
|----------|----------|------|
| 没有外部 holdout 是最大缺口 | **成立** | `B_GROUP_VERDICT.md`："C 组（decoy 受体/换构象/holdout）未开" |
| 结构稳健性（换 PDB）没做 | **成立** | 仅有失败的候选记录（BChE: 6ZWI/6QAA/5DYW 未过 QC；PIK3CB: 2Y3A/4BFR 未过 QC），没有"换了结构、结论仍同向"的正向证据 |
| 规模感偏案例集、PM 的 Δ CI 仍含 0 | **成立，但需要更准确的措辞** | `PRIMARY_METRIC_V2.md`/Results 3.3：PM 点估计超基线，Δ 95% CI 仍含 0；这是**功效不足**，不是"结果不成立" |
| 建议再扩 EGFR 或找第 5 对来加规模 | **需要修正：本轮不做** | `CLAIM_CEILING.md` 明令禁止"expand_EGFR_same_protocol_for_significance"；供给审计（J0/`strict_label_supply.csv`）已证明非金属对里只剩这 3 对够格，冒然找第 5 对风险是新的"薄面板" |
| 失败学分型未做成主文级 | **成立，但范围要修正** | 目前 `pik3ca_mtor_panel48_v0/analysis/failure_typology_v0/` 只对 **PM** 做了完整分型（T1/T2/T5 + 案例），EGFR 只有诊断级材料，AChE/BChE、PIK3CA/PIK3CB **完全没有**失败分型 |
| 标签层做不到 median/confidence 敏感性 | **成立且不可解** | `T0_SKIPS.md`：本地缓存每个 ChEMBL ID 只存了一个最大 pChEMBL 浮点数，没有逐条 assay 字段；这不是"没做"，是"数据源不支持"，只能写 Limitations，不能作为可补实验 |
| 单靶 sanity 只覆盖 PM 两端 | **成立** | `SINGLE_TARGET_ENRICHMENT_SANITY.md` 只有 4L23/4JT6 两个受体；AChE/BChE/PIK3CB/EGFR/HER2 五个受体没有单靶富集对照 |
| scaffold 指纹基线暴露的问题 | **需要补充一句**：这不是"实验缺口"，而是**已完成且已经诚实处理**的分析——`ML_BASELINE_LEAKAGE_CHECK.md` 已用 GroupKFold 修过泄漏，仍报告 0.78–0.91，且已写入 Limitations（singleton 支架、控制有限）。这条不需要再补实验，只需要在正文措辞里不要淡化 |

**结论：上一轮判断方向正确，需要两处修正**——(a) 不建议再找第 5 对/扩 EGFR；(b) 失败分型和单靶 sanity 的缺口范围比原说法更具体（分别缺 3 对里的 2 对）。真正该补的实验，按性价比排序是：

1. **ChEMBL 未用池 holdout**（零新增数据源风险，本方案已给出具体做法）
2. **结构稳健性**（PM 优先，AChE/BChE 或 PIK3CB 二选一）
3. **失败分型扩到 AChE/BChE 与 PIK3CA/PIK3CB**
4. **单靶 sanity 扩到非 PM 口袋**（至少 AChE 或 EGFR 一端）

---

## 1. 「文献双靶 holdout」怎么整理 —— 两条可选路径

### 1.1 为什么不直接去文献里现凑一张双靶药物表

用户原话是"文献双靶 holdout"，但直接去 SAR 论文/临床管线里现摘化合物名单，在本仓库的工作原则下有两个硬伤：

1. **无法只用真实数据、不编造** —— 论文里报告的 IC50/Ki 往往跟 ChEMBL 里的活性口径（pChEMBL 聚合规则）不一致，逐个手抄容易引入未审计的数字。
2. **多数"著名"双靶化合物已经在你们的冻结面板里** —— 例如 PI-103、Torin1、Omipalisib、WYE-132 已经是 PM48 里的条目（见 `failure_typology_v0`）；再把它们当"独立外部 holdout"是**信息重复**，不是真正的外推检验。

因此本方案给出两条路径，**推荐路径 A 作为主 holdout**，路径 B 作为可选、範圍更小的补充说明；两者互不冲突。

### 1.2 路径 A（推荐，已验证可行）：ChEMBL 冻结面板之外的「未用严格池」

**核心思想**：panel 建造时（`build_strict_panels.py` / `pm48_rdkit_prep_dock.py` 等）用严格标签规则（pChEMBL ≥6.5 / ≤5.5）在 ChEMBL 缓存里筛出的候选，永远比实际抽进面板的多。没被抽中的那部分，从来没有被用来定规则、选受体、调协议——它们是**真正意义上的"面板构建时未见过"的数据**，可以当作诚实的、有边界说明的 holdout。

**已验证的具体产出**（本轮已跑通，见 §4）：

| 靶对 | 严格池总量 (dual/A/B) | 已用于冻结面板 | 剩余未用池 (dual/A/B) | 抽取 holdout (dual/A/B) |
|------|------------------------|----------------|------------------------|--------------------------|
| PIK3CA/mTOR | 1552 / 80 / 81 | 74（含 PM48⊂PM110） | 1523 / 57 / 59 | **20 / 20 / 20** |
| AChE/BChE | 687 / 189 / 78 | 84 | 659 / 161 / 50 | **20 / 20 / 20** |
| PIK3CA/PIK3CB | 602 / 56 / 67 | 84 | 574 / 28 / 39 | **20 / 20 / 20** |

**整理步骤（脚本化、可复现，已在仓库落地为 `data/jcim_holdout_v0/scripts/build_holdout_candidate_pool_v1.py`）：**

1. 复用与冻结面板**完全相同**的标签规则（pChEMBL ≥6.5 / ≤5.5），从 `mols_*.json` 重算全量严格池 —— 不发明新规则。
2. 用 ChEMBL ID 精确排除已进入对应冻结面板 CSV 的全部条目（PM 用 panel110 的超集排除，覆盖 PM48）。
3. 在剩余池里，用**新种子**（`HOLDOUT_SEED=20260731`，区别于面板建造用过的所有种子）按 dual/A_only/B_only 各配额抽样，Murcko 支架封顶 3 个/类，避免同一化学型刷穿 holdout。
4. 只对被抽中的条目现查 ChEMBL API 拿 SMILES（不批量重下整库）。
5. **冻结顺序是规则→抽样→建列表→（此后才允许）对接→算分**；抽样清单写死后不得因为分数不好看而回去改配额或换种子。

**诚实的措辞边界（写入正文/SI 时必须遵守）：**
- 这是 **"未参与面板构建与协议调优的 ChEMBL 池"holdout**，不是跨数据库/跨机构的独立数据源；标题建议用 *"post-hoc unused-pool holdout"* 或中文"面板外冻结验证集"，不要写成"独立文献验证集"，避免过度宣称。
- 它仍然共享同一 ChEMBL 抓取批次（2026-07-23 锁定）和同一潜在的策展偏差；这一点必须写进 Limitations，不能让读者误以为这是完全异质的外部检验。
- 它可以回答的问题：**"分数规则和受体/协议是否只在建面板时凑出来的效果，换一批同规则、同来源但未见过的配体是否还同向？"**——这是该 holdout 的准确定位，不能拔高。

### 1.3 路径 B（可选、范围要收窄）：命名化学型交叉核对

如果仍想保留"文献"色彩，可选路径 B 作为路径 A 的**补充旁注**，而不是替代：

- 只在 PI3K/mTOR 双抑制剂这一个化学空间上做（该领域文献命名化合物多，如 NVP-BEZ235/dactolisib、GDC-0980、apitolisib(GSK2126458)、gedatolisib、voxtalisib(XL765) 等），检查这些**是否落在路径 A 抽出的未用池里**。
- 若命中：在 SI 用一两行点名"holdout 中包含文献报道的 X、Y 双抑制剂化学型"，增加可读性，但**不单独算它们的 AUROC**（样本太小、会被读成挑好例子）。
- 若像 PI-103/Torin1/Omipalisib 这类已经被主面板用掉，直接跳过，不去别处找替身——否则又变成"现凑名单"。

**结论：路径 A 是唯一进入主分析的 holdout；路径 B 至多作为 SI 里一句旁注。**

---

## 2. Holdout 对接与评估方案（防泄漏协议）

### 2.1 对接协议 —— 完全复用冻结协议，零调参

| 靶对 | 受体（复用现有 PDBQT/盒子） | Exhaustiveness | 配体制备 |
|------|------------------------------|-----------------|----------|
| PIK3CA/mTOR | 4L23 / 4JT6（`pik3ca_mtor_panel48_rdkit_v0/receptors`） | **16**（同主面板） | RDKit ETKDGv3 seed=20260727 + MMFF200 + meeko 默认，同 Methods 2.5 |
| AChE/BChE | 4EY7 / 4BDS（`ache_bche_panel_v0/receptors`） | 8 | 同上 |
| PIK3CA/PIK3CB | 4L23 / 2WXF（`pik3ca_pik3cb_panel_v0/receptors`） | 8 | 同上 |

- 打分读出：mode 1（Top-1）的 `REMARK VINA RESULT` 分数，与主表 `vina_A`/`vina_B` 定义一致。
- **不**重新做 cognate QC、不换盒子、不改 exhaustiveness —— 这些都是本方案要检验的"已冻结协议"，如果连协议都跟着 holdout 重新调，holdout 就失去意义。

### 2.2 评估规则 —— 与主表同一把尺子，不能二次调参

1. 主指标：**口袋匹配方向 AUROC**（dual vs A_only 用对端口袋分数；dual vs B_only 用另一对端口袋分数），summary_min 取两臂较小值 —— 与 Methods 2.6 定义逐字相同。
2. 置信区间：配体层 bootstrap，与主表同方法。
3. **不允许**：因为 holdout 上某个基线/某个臂"更好看"就换汇总方式（如临时改用池化或 worst-pocket）；holdout 只能用来验证已经写进 Methods 的规则,不能反过来筛选规则。
4. 报告口径：holdout summary_min 与主表 summary_min **并列展示**（森林图或表格并排一行),不合并进主表的 bootstrap 池,以保持"训练/验证-独立评价"的边界清晰。
5. 若 holdout 上某靶对结论方向翻转（如主表方向信号消失或反转）：**如实报告**，写入 Discussion，不重新采样holdout去"救"结果。

### 2.3 与平凡基线的对照

Holdout 同时计算 heavy_atoms / MW / cLogP / TPSA 四个平凡基线的 summary_min（与主表方法一致），保持"对接是否仍跑赢最强基线"的检验在 holdout 上同样成立/不成立。

---

## 3. 结构稳健性对接方案（Wave 2.4 落地版）

### 3.1 目标与优先级

回应"结论是口袋物理特异性、还是这一套晶体结构的偶然产物"。优先做 **PIK3CA/mTOR**（当前唯一过描述符基线的靶对，最值得证伪/证实），AChE/BChE 或 PIK3CB 二选一作为第二证据。

### 3.2 候选受体盘点（本轮已用 RCSB 检索 API 实查，非纸面猜测）

| 靶标 | 现用 | 已试且未过 cognate QC 的候选（**不要重复**） | 本轮实查到的新候选（未做 cognate QC，需下一步验证） |
|------|------|-------------------------------------------------|----------------------------------------------------|
| PIK3CA | 4L23 | 3T8M（看似过 RMSD 但确认为 **p110γ 嵌合体**，见 §3.3.1，**排除**） | **4JPS**（2.2 Å，真实 PIK3CA α，含配体）、**5DXT**（2.25 Å，p110α + GDC-0326）、5UBR（2.4 Å）、6GVF（2.5 Å）、5UK8（2.5 Å）——均经 RCSB `polymer_entity` 描述核实为 *Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform*，非嵌合体 |
| mTOR | 4JT6（E16 才过 QC） | — | **4JSX**（mTORΔN-mLST8-**Torin2** 复合物，3.5 Å）——与 4JT6 同一 ΔN-mLST8 截短构建体家族（mTOR 全长不结晶，这是领域标准做法），但结合的是不同配体（Torin2 vs PI-103）、不同晶型，可作为"同构建体家族、不同晶体"的稳健性检验 |
| BChE | 4BDS | 6ZWI（best9≈2.3–2.5Å）、6QAA/5DYW（PDBQT 解析失败） | 本轮未查到新候选，需下一步用 RCSB 全文检索 + 手动核对 entity 描述（避免重蹈 3T8M 覆辙） |
| PIK3CB | 2WXF | 2Y3A（≈3.85Å）、4BFR（解析失败）、**3T8M（本轮新排除，见下）** | 本轮未查到新候选；下一步需专门检索"PIK3CB/p110β + 小分子 + 分辨率 ≤3 Å + entity 描述含 beta"，并逐条核对是否为嵌合体 |

#### 3.2.1 本轮教训案例：3T8M 看似免费可用，实为嵌合体

`data/pik3ca_pik3cb_panel_v0/cognate_qc/3T8M_cognate_out_{E8,E16}.pdbqt` 早就在仓库里，本轮补算 RMSD 得到 **mode1 = best_of_9 = 0.461 Å**（远优于门槛），一度以为是"零成本"的第二个 PIK3CB 结构。但用 RCSB `polymer_entity` API 核对多肽实体描述后发现其为 **PI3K-γ (PIK3CG) catalytic subunit** 而非 PIK3CB，条目标题「Rational Design of PI3K-alpha Inhibitors that Exhibit Selectivity Over the PI3K-beta Isoform」对应的是**以 p110γ 为骨架、嫁接 α/β 特异性口袋残基的嵌合体**——这是 PI3K 结构生物学领域绕过 α/β 全长蛋白难结晶问题的常规策略,但嵌合体不能代表真实 PIK3CB 蛋白环境。**该案例已从结构稳健性候选中排除**，仅作为方法论附注保留（见 §3.3 步骤 1）。

### 3.3 执行步骤

1. **零成本步骤已做完**：3T8M 的 cognate RMSD 已补算并排除，见 §3.2.1（教训：cognate 门槛之外必须核对 `polymer_entity` 描述，排除嵌合体）。
2. 对仍缺候选的靶标（PIK3CA 第二结构、PIK3CB 第二结构、mTOR 第二结构、AChE 或 BChE 第二结构），按现有协议检索 PDB：优先选**分辨率 ≤2.5 Å、含类药小分子配体、多肽实体描述明确为目标基因本身（非嵌合体/非同源支架）、非同一课题组同批次沉积**的条目，与已用受体构象/配体化学型有区别（避免"换了 PDB 号但其实是同一构象"）。
3. 每个新候选先做 cognate redock QC（best_of_9 <2 Å 门槛，与 Methods 2.4 相同协议、相同 seed），**只有过 QC 的结构才进入对接**，未过的和已失败的一样记录在案，不重复尝试。
4. 通过 QC 的替代受体，对该靶对**冻结面板的全部配体子集**（建议取 dual+A_only+B_only 全量，PM48/AChE-BChE-100/PIK3CB-100 规模都可承受）重新对接，仅替换该受体，另一端不变。
5. 用相同口袋匹配公式重算 summary_min，与主表并列报告（"原受体 vs 替代受体"两行）。

### 3.4 判据

- 若替代结构下 summary_min 方向不变（PM 仍最高、其余仍不过基线）→ 写入正文作为"结构稳健"的正向证据。
- 若方向翻转 → 如实报告为"受体依赖"，这本身也是一个可发表的诊断结论（不是失败），但要在 Discussion 里明确降级相应 claim 的强度。

---

## 4. 本轮已落地的实操产出，及交给本地 Agent 的执行交接

云端本轮**只做到"零/低风险的数据与流水线准备"**，实际大批量对接按你的要求交给本地 Agent 执行，避免云端算力/时长限制影响真实结果。已完成部分：

1. **环境可行性已验证**：云端曾临时装过 `vina==1.2.7`（与冻结协议版本一致）、`meeko==0.7.1`、`gemmi`、`pandas`，确认 ChEMBL REST API 与 RCSB 均可直接访问，且用 Vina Python 绑定（`from vina import Vina`）+ meeko `mk_prepare_ligand.py` 跑通了完整"RDKit 建构象 → meeko 转 PDBQT → Vina 对接 → 解析 mode1 分数"链路，分数量级合理（AChE/BChE 端 mode1 约 −8 至 −13 kcal/mol，与主面板同量级）。**本地 Agent 只需自行安装同版本工具链即可复用全部脚本，无需改动路径**（脚本用 `Path(__file__).resolve().parents[3]` 定位仓库根，具有可移植性）。
2. **Holdout 候选池构建已完整跑完**（零对接，产出为最终版，不需要本地重新生成）：`data/jcim_holdout_v0/scripts/build_holdout_candidate_pool_v1.py` 产出：
   - `data/jcim_holdout_v0/tables/strict_pool_full_{HOPM,HOAB,HOAP}.csv` —— 三对靶标的严格标签全量池（含是否已用于冻结面板的标记）。
   - `data/jcim_holdout_v0/tables/holdout_panel_{HOPM,HOAB,HOAP}.csv` —— 冻结种子（20260731）抽出的 20/20/20（dual/A_only/B_only）holdout 配体，含 SMILES 与 Murcko 支架，**这份清单是最终清单，本地对接时不要重新抽样**。
3. **Holdout 对接脚本已就绪并做过小规模可行性验证**（`data/jcim_holdout_v0/scripts/dock_holdout_v1.py`），云端仅跑了 AChE/BChE 端的一小部分（54/120 个受体-配体组合，验证流程无误后已按你的要求停止），**未完成的部分需要本地 Agent 接着跑完**，脚本本身支持断点续跑（按已记录的 (受体, 配体) 组合去重，重跑同一命令即可继续，不会重复计算已有结果）。
4. **结构稳健性候选已完成一轮真实检索与一次重要证伪**：见 §3.2/§3.2.1（3T8M 看似可用、实为 PI3Kγ 嵌合体，已排除；PIK3CA 与 mTOR 各找到经身份核实的真实候选）。

### 4.1 交给本地 Agent 的具体操作步骤

```bash
cd Dual_Target_Docking

# 1. 环境（与冻结协议版本一致；rdkit/pandas/scipy 若已装可跳过）
pip install vina==1.2.7 meeko==0.7.1 gemmi rdkit pandas scipy

# 2. 直接对接（候选池与 holdout 清单已在仓库里，不要重新跑 build_holdout_candidate_pool_v1.py，
#    除非你要重新冻结抽样种子——那样会破坏"抽样规则先冻结、后看分数"的防泄漏约定）
python3 data/jcim_holdout_v0/scripts/dock_holdout_v1.py --prefix HOAB --workers 4   # 断点续跑，会跳过已完成的 54 个组合
python3 data/jcim_holdout_v0/scripts/dock_holdout_v1.py --prefix HOAP --workers 4
python3 data/jcim_holdout_v0/scripts/dock_holdout_v1.py --prefix HOPM --workers 4   # E=16，最慢，建议放最后
```

- `--workers` 按本地 CPU 核数调整（云端 4 核时单job E=8 约 15–100s，E=16 更慢）。
- 三个 CSV（`scores_vina_mode1_{HOAB,HOAP,HOPM}.csv`）跑完后即可按 §2.2 的口袋匹配公式重算 summary_min + bootstrap CI，与主表并列。
- 若要做 §3 的结构稳健性对接，需要先对新候选（4JPS/5DXT 等 PIK3CA；4JSX mTOR；BChE/PIK3CB 待补）跑一遍 cognate redock QC（复用 `ache_bche_panel_v0/scripts/freeze_receptors_cognate_qc_v2.py` 的逻辑，把里面写死的本地 conda 路径换成本地环境的 `vina`/`mk_prepare_receptor.py` 路径），只有过 <2 Å 门槛的才能进入下一步换受体重对接。

### 4.2 产出去向（本地跑完后回填到稿件）

三对 holdout 全部跑完后：新增 `HOLDOUT_VERDICT.md`（口袋匹配 summary_min + 95% CI + 与主表并列的森林图行）、Methods 补一段"Holdout evaluation"、Results 补一节"3.9 面板外冻结验证集"、SI 增 Table S7（holdout 配体清单 + 标签 + 分数）。结构稳健性对接跑完后同理补 Results/SI 段落（§3.4 判据）。

---

## 5. 剩余厚度补项（非本轮对接对象，按优先级列出留作后续）

| 优先级 | 内容 | 对接量级 | 说明 |
|--------|------|----------|------|
| P0 | 结构稳健性对接（§3） | PM 全量 ×1 替代受体端 ≈ 全面板×1端；AChE/BChE 或 PIK3CB 二选一同理 | 需要先完成候选检索+QC，方能定对接量 |
| P1 | 失败分型扩到 AChE/BChE、PIK3CA/PIK3CB | 0 新对接（复用已有分数，做姿态级解读） | 复刻 PM48 failure typology 的分析框架，不需新增 Vina 运行，只需导出 RTM-best 姿态做几何/化学型分析 |
| P1 | 单靶 sanity 扩到非 PM 口袋 | 每受体 ~50 活性 + 150 decoy ×1 受体 | 建议先补 AChE（供给最厚），EGFR 次之（供给受限案例本身也需要"对接没坏"的旁证） |
| P2 | 路径 B 文献化学型旁注 | 0 新对接 | 纯文本核对，见 §1.3 |

**明确不做**（与 `CLAIM_CEILING.md`、`JCIM_STRENGTHENING_PLAN_V1.md` §5 一致）：再找第 5 对靶标；扩大 EGFR 面板赌显著；把 holdout 结果用来反调标签规则或协议参数；把 warning flags/化学型旗标折进主分数。

---

## 6. 文件索引

| 路径 | 内容 |
|------|------|
| `data/jcim_holdout_v0/scripts/build_holdout_candidate_pool_v1.py` | Holdout 候选池构建（路径 A，零对接） |
| `data/jcim_holdout_v0/scripts/dock_holdout_v1.py` | Holdout 对接（复用冻结协议，可断点续跑） |
| `data/jcim_holdout_v0/tables/holdout_panel_{HOPM,HOAB,HOAP}.csv` | 三对靶标的冻结 holdout 配体清单 |
| `data/jcim_holdout_v0/tables/scores_vina_mode1_{HOPM,HOAB,HOAP}.csv` | Holdout Vina mode-1 分数（滚动写入） |
| `docs/JCIM_STRENGTHENING_PLAN_V1.md` | 总规划（Wave 0-4），本方案是其第 3 步的可执行细化 |
| `data/jcim_strengthen_t0t1_v0/analysis/B_GROUP_VERDICT.md` | 记录"C 组未开"的原始判断依据 |
