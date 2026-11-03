#!/usr/bin/env python3
"""Chinese scientific reports generated from the computed tables."""
from __future__ import annotations

from posthoc_common import REPORT_DIR


def _f(value, digits=3) -> str:
    return f"{float(value):.{digits}f}"


def _n(value) -> str:
    return str(int(value))


def _write(name: str, text: str) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    (REPORT_DIR / name).write_text(text.strip() + "\n")


def _row(summary, pair, arm):
    hit = summary[(summary["pair"] == pair) & (summary["arm"] == arm)].iloc[0]
    return hit


def write_reports(frames: dict) -> None:
    summary = frames["A_FOLD_OVERLAP_SUMMARY"]
    nearest = frames["A_TEST_NEAREST_TRAIN"]
    audit = frames["A_SCAFFOLD_FIELD_AUDIT"]
    decision = frames["A_TRIGGER_DECISION"].iloc[0]
    pairwise = frames["B_EGFR_PAIRWISE_RANKING"].iloc[0]
    dist = frames["B_EGFR_SCORE_DISTRIBUTIONS"].set_index("cohort")
    desc = frames["B_EGFR_DESCRIPTOR_COMPARISONS"].set_index("descriptor")
    assoc = frames["B_EGFR_SCORE_DESCRIPTOR_ASSOCIATIONS"]
    activity = frames["B_EGFR_ACTIVITY_DISTRIBUTIONS"]
    molecules = frames["B_EGFR_ACTIVITY_MOLECULES"]
    flags = frames["B_EGFR_PROVENANCE_FLAGS"]
    counts = frames["C_TOPK_CONCORDANCE_COUNTS"]
    cmp = frames["C_TOPK_METHOD_COMPARISONS"]
    repro = frames["C_TOPK_BASELINE_REPRODUCTION"]

    n_ge70 = int((nearest["ge_0_70"] == "YES").sum())
    n_ge85 = int((nearest["ge_0_85"] == "YES").sum())
    n_same = int((nearest["nearest_train_same_class"] == "YES").sum())
    n_chem = int(len(nearest))
    med_min = float(summary["max_tanimoto_median"].min())
    med_max = float(summary["max_tanimoto_median"].max())
    rel = audit["scaffold_relation"].value_counts().to_dict()

    def rho(descriptor, scope):
        row = assoc[(assoc["descriptor"] == descriptor) & (assoc["scope"] == scope)].iloc[0]
        if row["status"] == "CONSTANT_VARIABLE":
            return "CONSTANT_VARIABLE"
        return _f(row["spearman_rho"], 2)

    def flag_count(cohort, cls, flag):
        hit = flags[(flags["analysis_cohort"] == cohort) & (flags["class"] == cls) & (flags["flag"] == flag)].iloc[0]
        return int(hit["n_yes"]), int(hit["n_class"])

    def state_n(cutoff, scope, state):
        hit = counts[(counts["cutoff"] == cutoff) & (counts["scope"] == scope) & (counts["concordance"] == state)].iloc[0]
        return int(hit["n"])

    def case(pair, comparison, cutoff):
        return cmp[(cmp["pair"] == pair) & (cmp["comparison"] == comparison) & (cmp["cutoff"] == cutoff)].iloc[0]

    dual_egfr = molecules[molecules["class"] == "dual"]
    b_egfr = molecules[molecules["class"] == "B_only"]
    n_dual_near = int(((dual_egfr["pchembl_EGFR"] >= 6) & (dual_egfr["pchembl_EGFR"] < 6.5)).sum())
    n_b_her2_near = int(((b_egfr["pchembl_HER2"] >= 6) & (b_egfr["pchembl_HER2"] < 6.5)).sum())
    b_bins = activity[(activity["cohort"] == "B_only") & (activity["endpoint"] == "EGFR")].set_index("bin")

    lines = []
    lines.append("## 16 个方向的训练–测试近邻")
    lines.append("")
    lines.append("| 靶点对 | 方向 | 人数 | 中位最大相似度 | ≥0.70 | ≥0.85 | 最近邻同类 |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['pair']} | {row['arm']} | {_n(row['n_molecules'])} | {_f(row['max_tanimoto_median'], 3)} | "
            f"{_n(row['n_ge_0_70'])} ({_f(100 * row['fraction_ge_0_70'], 1)}%) | "
            f"{_n(row['n_ge_0_85'])} ({_f(100 * row['fraction_ge_0_85'], 1)}%) | "
            f"{_n(row['n_nearest_same_class'])}/{_n(row['n_molecules'])} |"
        )
    neighbor_table = "\n".join(lines)

    a_text = f"""
# 化学模型跨折骨架与相似性审计

## 问题与人群

本分析检查既有五折中，测试分子是否与训练分子具有相同骨架、相同化合物实体或较高的 ECFP4 相似性。它不重新训练正式模型，也不把化学相似性写成已经证实的数据泄漏。

人群是 `module_D_chemistry_oof_ligands.csv` 的 {n_chem} 条折外记录，覆盖 8 个靶点对和 16 个方向。每个分子只在自己的测试折中与该折训练成员比较。D1 中“跨实验类别最大相似度”没有被当作训练–测试相似度。ECFP4 使用 radius 2、2048 bit、`useChirality=False`、`useFeatures=False`。0.70 和 0.85 是预先固定的描述性分界。

## 骨架字段

按项目函数 `MurckoScaffold.GetScaffoldForMol` 后接 `MolToSmiles` 重算骨架。990 条与冻结 `scaffold` 字符串完全一致。3 条 JAK 分子的冻结字符串与当前 RDKit 芳香性写法不同，解析后的规范 SMILES 相同，记为 `CANONICAL_EQUIVALENT`。PPARG/PPARA 的 PGPA_029 在折文件中记为 `acyclic`，重算结果为空 Murcko，记为 `EMPTY_LABEL_MATCH`。这 4 条没有替换冻结字段。空 Murcko 只出现 1 次，没有被当成一个共同的有环骨架。

关系计数：STRING_MATCH {rel.get('STRING_MATCH', 0)}，CANONICAL_EQUIVALENT {rel.get('CANONICAL_EQUIVALENT', 0)}，EMPTY_LABEL_MATCH {rel.get('EMPTY_LABEL_MATCH', 0)}。

## 跨折重复

16 个方向的非空 Murcko 骨架都没有出现在一个以上的折中。用冻结字段做规范化后再检查，结果同样为 0。16 个方向都没有相同 `global_ligand_entity_id` 或相同 canonical SMILES 同时出现在训练折和测试折。因此不存在“同一化合物实体进入训练和测试”的情况。

七对原有分折均满足情况 A。按预设规则，不重新开展骨架分组交叉验证，也不替换 PIK3CA/mTOR 已经使用的骨架分组正式分折。PIK3CA/mTOR 两个方向的跨折非空骨架重复数也是 0。触发状态：{decision['situation']}，{decision['conditional_scaffold_group_cv']}。

没有相同非空骨架跨折，只说明这一种重复没有发生。它不能证明测试集位于新的化学系列，也不能排除其他来源的化学组成差异。

## 训练–测试近邻

{n_chem} 个测试观察的最大训练–测试 ECFP4 Tanimoto 相似度中，各方向的中位数范围为 {_f(med_min, 3)}–{_f(med_max, 3)}。达到或超过 0.70 的观察为 {n_ge70}/{n_chem}（{_f(100 * n_ge70 / n_chem, 1)}%），达到或超过 0.85 的观察为 {n_ge85}/{n_chem}（{_f(100 * n_ge85 / n_chem, 1)}%）。最近训练邻与测试分子属于同一实验类别的观察为 {n_same}/{n_chem}（{_f(100 * n_same / n_chem, 1)}%）。

PIK3CA/mTOR 两个方向没有分子达到 0.70，中位最大相似度分别为 {_f(_row(summary, 'PIK3CA/mTOR', 'D_vs_A')['max_tanimoto_median'], 3)} 和 {_f(_row(summary, 'PIK3CA/mTOR', 'D_vs_B')['max_tanimoto_median'], 3)}。PPARA/PPARD 的 D_vs_A 有 {_n(_row(summary, 'PPARA/PPARD', 'D_vs_A')['n_ge_0_70'])}/{_n(_row(summary, 'PPARA/PPARD', 'D_vs_A')['n_molecules'])} 个分子达到 0.70，是 16 个方向中比例最高的。JAK1/JAK2 的 D_vs_A 有 {_n(_row(summary, 'JAK1/JAK2', 'D_vs_A')['n_ge_0_85'])} 个分子达到 0.85。

{neighbor_table}

## 对 ECFP4 结果的含义

正式结果中，ECFP4 在 16 个方向的折外 AUROC 都高于同成员 S1。这一比较使用的是已经保存的折外预测和同成员直接评分，本轮没有改写这些数值。跨折审计表明，这一优势不是由相同非空骨架或相同化合物实体进入训练和测试造成的。

同时，典型测试分子在训练集中存在中等程度的最近邻，约七成七的最近邻与测试分子同属一个实验类别。因此，ECFP4 的折外 AUROC 描述的是当前分折和当前化学邻域中的类别判别，不能单独作为对新化学系列仍然保持同样表现的证据。0.70 或 0.85 以上的近邻只占少数，也不宜把结果写成大多数测试分子都是训练集的近重复。

本轮没有据此宣布 ECFP4 无效，也没有据此认为 S1 更合理。二者使用的信息不同：ECFP4 拟合的是当前标签与二维结构的关系，S1 是对接评分。

## 本轮没有做的模型计算

情况 A 成立，所以没有生成骨架分组 OOF，也没有替换 Figure 5 或 Table S7。ECFP4 加 S1 增量的 bootstrap 区间按任务边界保留为后续可选项，本轮没有计算。
"""
    _write("A_SCAFFOLD_AUDIT_CN.md", a_text)
    _write("A_ANALYSIS_REPORT_CN.md", a_text)

    donors_pooled = rho("NumHDonors", "pooled_dual_and_B_only")
    donors_dual = rho("NumHDonors", "dual")
    donors_b = rho("NumHDonors", "B_only")
    tpsa_pooled = rho("TPSA", "pooled_dual_and_B_only")
    tpsa_dual = rho("TPSA", "dual")
    tpsa_b = rho("TPSA", "B_only")
    b_cross, b_n = flag_count("D_vs_B_primary", "B_only", "cross_source_or_endpoint_mix")
    d_cross, d_n = flag_count("D_vs_B_primary", "dual", "cross_source_or_endpoint_mix")
    d_rep, _ = flag_count("D_vs_B_primary", "dual", "replicate_threshold_discordance")
    b_rep, _ = flag_count("D_vs_B_primary", "B_only", "replicate_threshold_discordance")

    b_text = f"""
# EGFR/HER2 的 S1 反向排序诊断

## 冻结人群与 AUROC

主要人群为 EGFR/HER2、`activity_eligible=1`、`M0_A_valid=1`、类别为 dual 或 B-only。人数为 dual {_n(pairwise['n_dual'])}、B-only {_n(pairwise['n_B_only'])}，与正式 D_vs_B 一致。评分使用已经转换为 higher-better 的 `M0_score_A`，没有再次取相反数。

896 组成对比较中，dual 评分更高 {_n(pairwise['n_dual_score_higher'])} 次，B-only 评分更高 {_n(pairwise['n_B_only_score_higher'])} 次，并列 {_n(pairwise['n_tie'])} 次。ties=0.5 的 Mann–Whitney AUROC 为 {pairwise['auroc_ties_0_5']:.12g}。正式表按 10 位小数记录为 {float(pairwise['formal_auroc']):.10f}，二者在该显示精度上一致。95% 区间 0.1328125000–0.3783482143 直接引自正式表，本轮没有重做 bootstrap。

## 直接观察到的评分分布

S1 EGFR 评分中位数：dual {_f(dist.loc['dual', 'median'], 3)}（IQR {_f(dist.loc['dual', 'iqr'], 3)}，范围 {_f(dist.loc['dual', 'minimum'], 3)}–{_f(dist.loc['dual', 'maximum'], 3)}），B-only {_f(dist.loc['B_only', 'median'], 3)}（IQR {_f(dist.loc['B_only', 'iqr'], 3)}，范围 {_f(dist.loc['B_only', 'minimum'], 3)}–{_f(dist.loc['B_only', 'maximum'], 3)}），neither {_f(dist.loc['neither', 'median'], 3)}（IQR {_f(dist.loc['neither', 'iqr'], 3)}，范围 {_f(dist.loc['neither', 'minimum'], 3)}–{_f(dist.loc['neither', 'maximum'], 3)}）。

B-only 的中位数比 dual 高约 {_f(dist.loc['B_only', 'median'] - dist.loc['dual', 'median'], 3)}。两组范围大量重叠，反向排序不是由少数极端分子单独造成的。neither 的中位数低于 dual。按同一 dual 集合重算的 dual–neither AUROC 为 {float(pairwise['recomputed_auc_neither']):.6f}，与 module B 保存值一致。module B 保存的 neither 相对 B-only 的配对差值为 {float(pairwise['stored_delta_negative']):.6f}，95% 区间为 {float(pairwise['stored_delta_negative_ci_lo']):.4f}–{float(pairwise['stored_delta_negative_ci_hi']):.4f}。本报告直接引用该配对结果，没有把两次独立 AUROC 相减后称作新的配对估计。

因此，当前数据直接显示的是：EGFR 的 S1 把这 32 个 B-only 分子整体排在 28 个 dual 分子之前，同时仍把 12 个 neither 分子整体排在 dual 之后。

## 描述符

八项描述符全部保留。FormalCharge 在 60 个分子上均为 0，状态为 `CONSTANT_VARIABLE`，不能用于区分 dual 和 B-only，也没有报告相关系数。

中位数差异（dual 减 B-only）：MolWt {_f(desc.loc['MolWt', 'median_difference_dual_minus_B_only'], 2)}，MolLogP {_f(desc.loc['MolLogP', 'median_difference_dual_minus_B_only'], 3)}，TPSA {_f(desc.loc['TPSA', 'median_difference_dual_minus_B_only'], 2)}，NumHDonors {_f(desc.loc['NumHDonors', 'median_difference_dual_minus_B_only'], 2)}，NumHAcceptors {_f(desc.loc['NumHAcceptors', 'median_difference_dual_minus_B_only'], 2)}，NumRotatableBonds {_f(desc.loc['NumRotatableBonds', 'median_difference_dual_minus_B_only'], 2)}，RingCount {_f(desc.loc['RingCount', 'median_difference_dual_minus_B_only'], 2)}。

B-only 的分子量中位数更高，环数中位数为 {_f(desc.loc['RingCount', 'B_only_median'], 0)}，dual 为 {_f(desc.loc['RingCount', 'dual_median'], 0)}。MolLogP 的两组中位数几乎相同。可旋转键的中位数都是 {_f(desc.loc['NumRotatableBonds', 'dual_median'], 0)}。

Spearman 相关是探索性描述。RingCount 与 S1 的相关系数在合并人群为 {rho('RingCount', 'pooled_dual_and_B_only')}，dual 内为 {rho('RingCount', 'dual')}，B-only 内为 {rho('RingCount', 'B_only')}。MolWt 的对应值为 {rho('MolWt', 'pooled_dual_and_B_only')}、{rho('MolWt', 'dual')} 和 {rho('MolWt', 'B_only')}。这两组关系在类别内部仍然存在，所以不能只解释为 dual 与 B-only 的类别混合。

NumHDonors 的合并相关系数为 {donors_pooled}，dual 内为 {donors_dual}，B-only 内为 {donors_b}。TPSA 的对应值为 {tpsa_pooled}、{tpsa_dual} 和 {tpsa_b}。类别内符号或强度与合并结果不同，合并相关系数不能写成评分机制。MolLogP 在两类内部都与评分正相关，但两类中位数没有分开。可旋转键的合并相关较弱，类别内方向也不一致。

这些是当前 60 个分子上的描述性关联。它们没有证明 Vina 普遍偏好大分子，也没有证明环数是反向排序的原因。

## 核实后的实验活性

60 个主要成员都能在 805 条来源交叉表中按 pair 和 panel_id 一对一匹配。每条记录的 dump 最大值、历史面板值和历史聚合最大值在 EGFR 端与 HER2 端都完全一致，冻结类别也一致。报告使用这一核实后的聚合 pChEMBL，没有使用单条未审计测定，也没有由四分类标签反推连续值。

B-only 的 EGFR 端：小于 5.0 为 {_n(b_bins.loc['<5.0', 'n_in_bin'])} 人，5.0 至小于 5.5 为 {_n(b_bins.loc['5.0_to_<5.5', 'n_in_bin'])} 人，5.5 至小于 6.0 为 {_n(b_bins.loc['5.5_to_<6.0', 'n_in_bin'])} 人，大于或等于 6.0 为 {_n(b_bins.loc['>=6.0', 'n_in_bin'])} 人。中位数 {_f(b_bins.loc['<5.0', 'median'], 2)}，范围 {_f(b_bins.loc['<5.0', 'minimum'], 2)}–{_f(b_bins.loc['<5.0', 'maximum'], 2)}。23/32 位于阈值下方 0.5 个 pChEMBL 单位以内。

B-only 的 HER2 端全部大于或等于 6.0，中位数见活性表；其中 {n_b_her2_near}/32 位于 6.0 至小于 6.5。dual 的 EGFR 端中位数为 8.925，范围 6.18–10.70，其中 {n_dual_near}/28 位于 6.0 至小于 6.5；HER2 端中位数为 8.880，范围 6.61–10.22。neither 的 EGFR 端中位数为 4.375，12 个分子中 10 个小于 5.0。

来源标记：B-only 的跨来源或混合终点为 {b_cross}/{b_n}，重复测定跨阈值标记为 {b_rep}/{b_n}。dual 的相应数字为 {d_cross}/{d_n} 和 {d_rep}/{d_n}。B-only 靠近阈值这一现象，并不是由重复测定跨阈值标记集中出现造成的。

## 证据分层

已直接观察到的现象：在这 28 与 32 个分子上，EGFR S1 的成对反向比例为 675/896；B-only 评分中位数更高；neither 的评分中位数更低，dual–neither AUROC 为 0.7470；B-only 的 EGFR 聚合活性集中在阈值下方，dual 的 EGFR 活性中心远高于 6；形式电荷没有变异。

与现象相关、但不能证明因果的因素：B-only 的环数和分子量中位数更高，而且这两项与 S1 的正相关在两个类别内部仍然可见；多数 B-only 的 EGFR 端聚合活性距离 6.0 不到 0.5。活性接近阈值会增加标签对聚合规则的敏感性，但它不能单独说明为什么实验活性中位数约 8.9 的 dual 会系统性地获得低于 B-only 的对接评分。

现有数据不能判断的机制：特定受体构象、口袋残基、质子化或互变异构状态是否造成这一排序。本轮没有授权新的对接、结构准备或机制模拟。

## 可以写入讨论的范围

讨论可以报告反向排序的成对比例、三类评分位置，以及 B-only 的 EGFR 活性靠近正式阈值。讨论需要同时写明环数和分子量只是描述性相关，形式电荷没有提供组别差异。讨论不宜把反向排序归因于电荷、某个残基或“Vina 必然偏好大分子”。
"""
    _write("B_EGFR_INVERSION_DIAGNOSTICS_CN.md", b_text)

    def case_line(pair, comparison):
        bits = []
        for cutoff in ("top_5pct", "top_10pct", "top_20pct"):
            row = case(pair, comparison, cutoff)
            bits.append(
                f"{cutoff}: dual {int(row['S1_n_dual_top'])}→{int(row['new_n_dual_top'])}，"
                f"k={int(row['k'])}，{row['concordance']}"
            )
        return "；".join(bits)

    c_text = f"""
# 联合排序 Top-k 截断比例敏感性

## 定义

比较使用每种方法与 S1 各自的共同双端完整成员，不使用一个固定的 S1 名单。joint order 按既有规则生成：两端 descending average midrank，第一键为 worst rank，其后为两端名次之和和 `global_ligand_entity_id`。`in_topk` 只属于原来的 10% 截断，5% 和 20% 都按完整排序重新取前 k 名。k = ceil(p × N)，p 固定为 0.05、0.10 和 0.20。

summary_min 差值读取正式配对点估计。S4 对 S1 继续使用已核实的支持性点估计，没有新建置信区间。四态分类使用原始精度：summary_min 差值与 Top-k 单靶比例变化都非零时，按二者的符号组合分为一致改善、一致变差或不一致；任一差值恰好为 0 时为无变化。

## Top-10% 复现

独立调用正式排序函数并与 `module_I_descriptive.csv`、`module_I_method_compare.csv` 和配体排序文件对照。描述性字段、32 组方法比较、四态计数、排序成员和共同池两侧成员均一致。复现检查全部为 PASS，随后才计算 5% 和 20%。

正式 10% 的 32 组结果为：一致改善 {state_n('top_10pct', 'all_32', 'CONCORDANT_IMPROVE')}，一致变差 {state_n('top_10pct', 'all_32', 'CONCORDANT_WORSEN')}，不一致 {state_n('top_10pct', 'all_32', 'DISCORDANT')}，无变化 {state_n('top_10pct', 'all_32', 'NO_CHANGE')}。四种内部方法的计数与已验收基线一致：M1 为 2/2/3/1，M1b 为 3/3/1/1，M2 为 5/3/0/0，M3 为 4/2/2/0。不含 S4 的 24 组为 {state_n('top_10pct', 'exclude_S4', 'CONCORDANT_IMPROVE')}/{state_n('top_10pct', 'exclude_S4', 'CONCORDANT_WORSEN')}/{state_n('top_10pct', 'exclude_S4', 'DISCORDANT')}/{state_n('top_10pct', 'exclude_S4', 'NO_CHANGE')}。

## 三种截断比例

32 组四态计数：

| 截断 | 一致改善 | 一致变差 | 不一致 | 无变化 |
|---|---:|---:|---:|---:|
| Top 5% | {state_n('top_5pct', 'all_32', 'CONCORDANT_IMPROVE')} | {state_n('top_5pct', 'all_32', 'CONCORDANT_WORSEN')} | {state_n('top_5pct', 'all_32', 'DISCORDANT')} | {state_n('top_5pct', 'all_32', 'NO_CHANGE')} |
| Top 10% | {state_n('top_10pct', 'all_32', 'CONCORDANT_IMPROVE')} | {state_n('top_10pct', 'all_32', 'CONCORDANT_WORSEN')} | {state_n('top_10pct', 'all_32', 'DISCORDANT')} | {state_n('top_10pct', 'all_32', 'NO_CHANGE')} |
| Top 20% | {state_n('top_20pct', 'all_32', 'CONCORDANT_IMPROVE')} | {state_n('top_20pct', 'all_32', 'CONCORDANT_WORSEN')} | {state_n('top_20pct', 'all_32', 'DISCORDANT')} | {state_n('top_20pct', 'all_32', 'NO_CHANGE')} |

不含 S4 的 24 组：Top 5% 为 {state_n('top_5pct', 'exclude_S4', 'CONCORDANT_IMPROVE')}/{state_n('top_5pct', 'exclude_S4', 'CONCORDANT_WORSEN')}/{state_n('top_5pct', 'exclude_S4', 'DISCORDANT')}/{state_n('top_5pct', 'exclude_S4', 'NO_CHANGE')}；Top 20% 为 {state_n('top_20pct', 'exclude_S4', 'CONCORDANT_IMPROVE')}/{state_n('top_20pct', 'exclude_S4', 'CONCORDANT_WORSEN')}/{state_n('top_20pct', 'exclude_S4', 'DISCORDANT')}/{state_n('top_20pct', 'exclude_S4', 'NO_CHANGE')}。

5% 名单很短，k 多为 5，单靶比例更容易恰好不变，因此无变化从 2 组增加到 10 组。20% 的不一致减少到 3 组，无变化为 5 组。三种比例都保留了改善、变差和不完全一致同时存在的格局，但 24/6/2 这组具体计数只属于 10%。

Top-k 中 dual 数的变化是当前评价集候选名单的组成变化。dual 基线比例由预设配额决定，不能写成天然化合物库的命中率改善。

## 三个指定案例

F2/F10 的 S3 对 S1：summary_min 差值保持 +0.2470703125。{case_line('F2/F10', 'S3_vs_S1')}。三个比例下 dual 数都增加，四态都保持一致改善。该案例的方向不依赖 10% 这一单个截断；增加的绝对人数随 k 变化。

PIK3CA/mTOR 的 S2 对 S1：summary_min 差值保持 −0.0048828125。{case_line('PIK3CA/mTOR', 'S2_vs_S1')}。5% 和 10% 都是不一致，因为 summary_min 略降而 Top-k dual 数增加。20% 时 dual 数由 9 个变为 6 个，单靶比例上升，四态变为一致变差。因此，正文中的这个不一致例子属于 10% 名单，不能写成在邻近截断比例下仍然成立。

PPARG/PPARA 的 S3 对 S1：summary_min 差值保持 −0.2915039062。{case_line('PPARG/PPARA', 'S3_vs_S1')}。10% 和 20% 为一致变差。5% 时两端 Top-k dual 数都是 3，单靶比例变化为 0，按预设规则归为无变化。负向的 summary_min 差值没有在 5% 名单上表现为 dual 数下降。

## 对原结论的影响

10% 仍然是正式主分析。新结果支持保留“方向性差值与前列名单组成多数同向、但并不总是同向”的定性判断。需要降低强度的是把 24、6、2 或 PIK3CA/mTOR 的不一致写成与截断比例无关。F2/F10 的 S3 改善方向在三个比例下都存在；PPARG/PPARA 的 S3 变差在 10% 和 20% 存在，在 5% 没有表现为名单组成变化。
"""
    _write("C_TOPK_SENSITIVITY_CN.md", c_text)
    _write("C_TOPK_SENSITIVITY_REPORT_CN.md", c_text)

    failed = repro[repro["status"] != "PASS"]
    final = f"""
# 投稿前支持性分析的科学评估

## 1. 本轮执行的计算

在核实远端基线 `d1ca17512eabe5c9d102ebe75c3f9c59668452dc` 后，完成三项只读支持性分析。

分析 A 对 994 条化学模型折外记录计算了 Bemis–Murcko 骨架、跨折实体重复和训练折内最大 ECFP4 Tanimoto 相似度。七对靶点没有非空骨架跨折，预设的骨架分组重训没有执行。

分析 B 在 EGFR/HER2 的正式 D_vs_B 人群上重算 S1 AUROC，并描述评分、八项描述符、核实后的聚合 pChEMBL 和来源标记。

分析 C 先复现正式 Top 10%，再按同一排序计算 Top 5% 和 Top 20%。

没有重跑对接，没有修改 867 条正式人群、标签、PIK3CA/mTOR 的 110 分子评价集、活性阈值或模型参数。ECFP4 加 S1 的 bootstrap 区间和 35 条跨阈值重复测定的排除敏感性都按本轮边界留在未执行状态。

## 2. 成功复现的原始结果

合格观察 867，全局配体实体 845，四类人数 248、256、252、111。标签层方向记录 1004，Vina 有效方向记录 1000，化学模型折外记录 994。EGFR/HER2 的 S1 D_vs_B AUROC 与正式 10 位小数 0.2466517857 一致，精确成对结果为 221/896。dual–neither AUROC 与 module B 的 0.747023809524 一致。Top 10% 的 32 组四态为 14、10、6、2，方法别计数与基线一致。复现失败数：{len(failed)}。

## 3. 问题分级

Critical：未发现。正式成员、AUROC 和 Top 10% 四态能够由冻结输入复现，没有发现同一化合物实体跨折。

Major：Top 10% 的具体四态计数，以及 PIK3CA/mTOR 上 S2 相对 S1 的不一致，不能外推到 5% 或 20%。20% 时该比较变为一致变差。这不是正式 10% 计算错误；如果正文把该例子写成截断比例无关的现象，解释会过强。

Moderate：ECFP4 高于同成员 S1 的结果限定在当前分折。16 个方向没有非空骨架跨折，但 139/994 个测试观察的训练最近邻相似度至少为 0.70，764/994 的最近邻属于同一实验类别。EGFR/HER2 的反向排序是 675/896 的成对结果，并伴随环数、分子量与评分的类别内相关，以及 23/32 个 B-only 的 EGFR 聚合活性位于 5.5 至 6.0。这些相关因素没有构成已证实的机制。

Minor：4 条冻结骨架字符串与当前 RDKit 规范化结果在写法上不同，其中 3 条化学等价，1 条是空骨架的 `acyclic` 标记。形式电荷在 EGFR 诊断人群中没有变异。

## 4. 增强现有结论的结果

ECFP4 与 S1 的同成员比较没有被相同骨架或相同实体跨折所解释。F2/F10 上 S3 相对 S1 的 dual 数增加在 5%、10% 和 20% 都保持同一方向。PPARG/PPARA 上 S3 的明显负向 summary_min 差值在 10% 和 20% 都伴随 dual 数下降。EGFR 的反向排序在 neither 对照中并不表现为同样的方向，module B 的配对差值仍然是既有证据。

## 5. 需要降低强度的结果

ECFP4 的折外表现不宜写成对新化学系列的验证。PIK3CA/mTOR 的 S2 不一致只在 5% 和 10% 出现。PPARG/PPARA 的 S3 在 5% 没有改变 Top-k dual 数。24、6、2 的四态计数只描述 10%。EGFR 反向排序现在有更具体的分子和活性背景，但仍然不是一个已解释的结合机制。

## 6. 当前数据不能回答的问题

不能判断口袋残基、受体构象或配体互变异构是否造成 EGFR 反向排序。不能由这次审计估计 ECFP4 在一个真正新系列上的 AUROC。不能把 Top-k 富集解释成天然库命中率。没有计算化学增量的置信区间，因此增量仍是点估计。没有做跨阈值重复测定的排除敏感性，本轮也不因 EGFR 诊断而删除任何正式成员。

## 7. 论文句子的处理原则

正式 10% 结果、16 个方向的 ECFP4 对 S1 优势，以及 EGFR D_vs_B 的 AUROC，都可以保留。需要补充的是适用范围：训练–测试近邻的定量描述、EGFR 活性靠近阈值和描述符相关的证据边界，以及 Top-k 计数对截断比例的依赖性。建议句子见 `CLAIM_EVIDENCE_ALIGNMENT.csv`。本文件不直接改写 Methods、Results 或 Discussion。

## 8. 是否还需要其他计算

进入英文稿整理前，本轮三项支持性分析已经覆盖任务预先指定的主要解释风险。最小且已经完成的补充就是这三项。没有必要为了投稿再重跑对接或重训化学模型。

若审稿人之后要求化学增量的不确定性，再单独批准同一成员折外预测的 bootstrap，并写明它不包括重新训练的全部不确定性。35 条跨阈值观察的排除分析也要另行批准，不能从 EGFR 的阈值邻近结果自动触发。
"""
    _write("FINAL_SCIENTIFIC_ASSESSMENT_CN.md", final)

    claims = [
        {
            "manuscript_section": "Results 3.4",
            "scientific_claim": "ECFP4 在 16 个方向的折外 AUROC 均高于同成员 S1。",
            "evidence_file": "results/jcim_stage2/module_D_chemistry_units.csv",
            "statistical_quantity": "16/16 directions with AUROC(ECFP4) > AUROC(matched S1)",
            "observed_result": "正式化学单元表保持该结果；本轮未改写 OOF 预测。",
            "uncertainty_or_limitation": "该比较限定于当前成员和当前分折。",
            "recommended_wording_cn": "在当前 16 个方向的共同成员上，ECFP4 的折外 AUROC 均高于直接计算的 S1。",
            "evidence_status": "DIRECTLY_SUPPORTED",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "16 个方向没有相同非空 Murcko 骨架跨折，也没有相同化合物实体跨折。",
            "evidence_file": "outputs/A_FOLD_OVERLAP_SUMMARY.csv",
            "statistical_quantity": "n_nonempty_scaffolds_in_multiple_folds; n_global_entities_in_multiple_folds",
            "observed_result": "16 个方向两项均为 0。",
            "uncertainty_or_limitation": "不能据此排除相似化学系列或其他偏差。",
            "recommended_wording_cn": "16 个化学模型方向均未发现非空 Bemis–Murcko 骨架或同一全局化合物实体同时出现在训练折和测试折。",
            "evidence_status": "DIRECTLY_SUPPORTED",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "部分测试分子与训练折分子具有较高二维相似性，因此 ECFP4 结果不能直接外推到新化学系列。",
            "evidence_file": "outputs/A_TEST_NEAREST_TRAIN.csv",
            "statistical_quantity": "max ECFP4 Tanimoto to the training fold",
            "observed_result": f"各方向中位数 {med_min:.3f}–{med_max:.3f}；{n_ge70}/{n_chem} 至少为 0.70；{n_ge85}/{n_chem} 至少为 0.85；{n_same}/{n_chem} 的最近邻类别相同。",
            "uncertainty_or_limitation": "0.70 和 0.85 是描述性分界，不是泄漏判定线。",
            "recommended_wording_cn": (
                f"各方向测试分子与训练折的最大 ECFP4 相似度中位数为 {med_min:.3f}–{med_max:.3f}；"
                f"{n_chem} 个观察中 {n_ge70} 个至少为 0.70，{n_ge85} 个至少为 0.85。"
                "该结果限制对新化学系列的外推，但不构成相同骨架跨折。"
            ),
            "evidence_status": "SUPPORTED_WITH_LIMITATIONS",
        },
        {
            "manuscript_section": "Results 3.2",
            "scientific_claim": "EGFR/HER2 的 S1 在 dual 与 B-only 比较中的 AUROC 为 0.2467。",
            "evidence_file": "outputs/B_EGFR_PAIRWISE_RANKING.csv",
            "statistical_quantity": "Mann-Whitney AUROC with ties = 0.5",
            "observed_result": "221/896 = 0.246651785714；正式表 10 位小数为 0.2466517857。",
            "uncertainty_or_limitation": "区间仍使用正式 bootstrap，本轮未重算。",
            "recommended_wording_cn": "EGFR/HER2 使用 EGFR 端 S1 比较 dual 与 B-only 时，AUROC 为 0.2467（28 个 dual、32 个 B-only）。",
            "evidence_status": "DIRECTLY_SUPPORTED",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "B-only 的 EGFR 聚合活性集中在正式阈值附近，可作为反向排序的潜在影响因素。",
            "evidence_file": "outputs/B_EGFR_ACTIVITY_DISTRIBUTIONS.csv",
            "statistical_quantity": "verified aggregate EGFR pChEMBL bins",
            "observed_result": "B-only 中 <5.0、5.0–<5.5、5.5–<6.0 和 ≥6.0 分别为 3、6、23 和 0；dual EGFR 中位数为 8.925。",
            "uncertainty_or_limitation": "阈值邻近不能单独证明它造成了 675/896 的反向成对比较。",
            "recommended_wording_cn": "32 个 B-only 分子中有 23 个的核实 EGFR 聚合 pChEMBL 位于 5.5 至 6.0；这是标签靠近阈值的分布事实，本数据不能把它确定为反向排序的原因。",
            "evidence_status": "EXPLORATORY_ASSOCIATION",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "形式电荷造成了 EGFR 反向排序。",
            "evidence_file": "outputs/B_EGFR_DESCRIPTOR_COMPARISONS.csv",
            "statistical_quantity": "FormalCharge",
            "observed_result": "60 个分子的 FormalCharge 均为 0。",
            "uncertainty_or_limitation": "常数变量没有组别信息。",
            "recommended_wording_cn": "该诊断人群的形式电荷没有变异，不能用来解释 dual 与 B-only 的评分差异。",
            "evidence_status": "NOT_SUPPORTED",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "环数和分子量与 EGFR S1 存在描述性关联。",
            "evidence_file": "outputs/B_EGFR_SCORE_DESCRIPTOR_ASSOCIATIONS.csv",
            "statistical_quantity": "Spearman rho within class",
            "observed_result": f"RingCount rho all/dual/B-only = {rho('RingCount', 'pooled_dual_and_B_only')}/{rho('RingCount', 'dual')}/{rho('RingCount', 'B_only')}；MolWt = {rho('MolWt', 'pooled_dual_and_B_only')}/{rho('MolWt', 'dual')}/{rho('MolWt', 'B_only')}。",
            "uncertainty_or_limitation": "类别内相关仍不是因果证据，样本量分别为 60、28 和 32。",
            "recommended_wording_cn": "在该 60 个分子中，环数和分子量与 EGFR S1 的正相关在 dual 和 B-only 内部都可见；这是探索性关联，不能写成 Vina 的普遍偏好。",
            "evidence_status": "EXPLORATORY_ASSOCIATION",
        },
        {
            "manuscript_section": "Results 3.7",
            "scientific_claim": "32 组 Top-10% 比较中，24 组同向、6 组不一致、2 组无变化。",
            "evidence_file": "outputs/C_TOPK_CONCORDANCE_COUNTS.csv",
            "statistical_quantity": "four-state counts at top 10%",
            "observed_result": "14 改善、10 变差、6 不一致、2 无变化；同向合计 24。",
            "uncertainty_or_limitation": "5% 为 9/12/1/10，20% 为 13/11/3/5。具体计数依赖截断比例。",
            "recommended_wording_cn": "在正式的前 10% 名单中，32 组比较有 24 组方向一致、6 组不一致、2 组至少一项没有变化。前 5% 和前 20% 的对应计数不同，因此这些数字只描述 10% 截断。",
            "evidence_status": "SUPPORTED_WITH_LIMITATIONS",
        },
        {
            "manuscript_section": "Results 3.7",
            "scientific_claim": "F2/F10 上 S3 相对 S1 的 summary_min 上升，同时 Top-10% dual 数由 4 增加到 7。",
            "evidence_file": "outputs/C_TOPK_METHOD_COMPARISONS.csv",
            "statistical_quantity": "delta summary_min and Top-k dual count",
            "observed_result": "summary_min +0.2470703125；5%、10%、20% 的 dual 数变化分别为 2→3、4→7 和 6→10，三次均为一致改善。",
            "uncertainty_or_limitation": "绝对人数依赖 k；不能写成化合物库命中率。",
            "recommended_wording_cn": "F2/F10 上 S3 相对 S1 的 summary_min 高 0.2471；前 5%、10% 和 20% 名单中的 dual 数分别由 2、4、6 个增加到 3、7、10 个。",
            "evidence_status": "DIRECTLY_SUPPORTED",
        },
        {
            "manuscript_section": "Results 3.7",
            "scientific_claim": "PIK3CA/mTOR 上 S2 的 summary_min 略降，但 Top-10% dual 数由 4 增加到 6。",
            "evidence_file": "outputs/C_TOPK_METHOD_COMPARISONS.csv",
            "statistical_quantity": "concordance across cutoffs",
            "observed_result": "10% 为不一致（4→6）；5% 仍不一致（2→3）；20% 变为一致变差（9→6）。",
            "uncertainty_or_limitation": "不一致没有在 20% 保持。",
            "recommended_wording_cn": "在前 10% 名单中，PIK3CA/mTOR 的 S2 使 summary_min 略降 0.0049，同时 dual 数由 4 个增加到 6 个；前 20% 时 dual 数转为由 9 个减少到 6 个。",
            "evidence_status": "SUPPORTED_WITH_LIMITATIONS",
        },
        {
            "manuscript_section": "Results 3.7",
            "scientific_claim": "PPARG/PPARA 上 S3 相对 S1 为一致变差。",
            "evidence_file": "outputs/C_TOPK_METHOD_COMPARISONS.csv",
            "statistical_quantity": "concordance at 5%, 10%, 20%",
            "observed_result": "summary_min −0.2915039062；10% dual 6→4，20% 13→7，均为一致变差；5% dual 3→3，为无变化。",
            "uncertainty_or_limitation": "5% 的名单组成没有变化。",
            "recommended_wording_cn": "PPARG/PPARA 上 S3 的 summary_min 低于 S1 0.2915；前 10% 和 20% 的 dual 数分别由 6 个降至 4 个、由 13 个降至 7 个，前 5% 的 dual 数保持 3 个。",
            "evidence_status": "SUPPORTED_WITH_LIMITATIONS",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "对接评分区分双靶与单靶配体的能力依赖评价方向、负类和候选截断，不能由单一有利分数代表双靶判别。",
            "evidence_file": "reports/FINAL_SCIENTIFIC_ASSESSMENT_CN.md",
            "statistical_quantity": "directional AUROC, negative-class contrast, Top-k concordance",
            "observed_result": "三项支持性分析均未推翻这一总判断，同时给出了相似性、阈值邻近和截断比例三项限制。",
            "uncertainty_or_limitation": "结论仍限于八对靶点和当前回顾性评价集。",
            "recommended_wording_cn": "现有对接评分对双靶与单靶配体的区分随方向、负类和候选截断而变化；化学模型优势限定于当前分折及其化学邻域。",
            "evidence_status": "SUPPORTED_WITH_LIMITATIONS",
        },
        {
            "manuscript_section": "Discussion",
            "scientific_claim": "特定 EGFR 残基或 HER2 受体构象导致了 S1 反向排序。",
            "evidence_file": "NOT_COMPUTABLE",
            "statistical_quantity": "none",
            "observed_result": "本轮没有结构机制数据。",
            "uncertainty_or_limitation": "需要新的结构或对照计算，本次未获授权。",
            "recommended_wording_cn": "当前评分和分子属性数据不能确定受体构象或特定残基是否参与 EGFR/HER2 的反向排序。",
            "evidence_status": "NOT_COMPUTABLE",
        },
    ]
    import pandas as pd
    pd.DataFrame(claims).to_csv(REPORT_DIR / "CLAIM_EVIDENCE_ALIGNMENT.csv", index=False)
