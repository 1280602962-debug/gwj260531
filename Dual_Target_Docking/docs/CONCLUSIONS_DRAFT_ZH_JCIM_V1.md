# Conclusions（中文工作稿 · JCIM Articles）

> 投稿以英文为准：[`CONCLUSIONS_SECTION_JCIM_EN_V1.md`](CONCLUSIONS_SECTION_JCIM_EN_V1.md)。  
> 两段：做了什么 + 得到了什么；意味着什么 + 未来评价标准。不复制 Results，不写 PDB/引擎数字，不写 validated / robust / 通用决策规则。  
> 主张边界：[`CLAIM_CEILING.md`](../data/jcim_bench_v0/CLAIM_CEILING.md)。术语分层见文末。

---

## 5. Conclusions

本研究建立 DualFourClass-Bench，作为对接双靶识别的实验锚定评价环境，显式检验对接能否将 dual-active 配体与 A-selective、B-selective 硬负样本区分开来。在四对冻结靶标上，口袋匹配对接表现出有限且高度依赖靶对的判别，summary_min AUROC 介于 0.430 与 0.692 之间。PIK3CA/mTOR 给出最高点估计，并在未参与建面的配体池中保持正向方向信号；但主面板估计的不确定度及其对受体结构的敏感性，排除将其解释为可迁移的双靶决策规则。

更广泛的分析表明，表观双靶对接性能由任务定义、配体化学组成和受体实现方式共同决定。EGFR/HER2 上，Dual-versus-neither comparator（AUROC 0.756）会支持对接双靶识别，而方向性弱臂仍为 0.430；该 formulation gap 依赖靶对，不是四对统一反转，也是描述性对照而非配对显著性检验。若干靶对上，物化或化学型参考达到或超过对接判别；未使用配体池还暴露出错口袋对照不低于口袋匹配的未解决反转，效价或尺寸匹配未能消除，尽管相应配对置信区间仍包含零。受体替换可按靶对升高或降低成对判别，这是 realization effect，不是稳健性证书。这些发现反对把两个口袋上的有利分数当作双靶活性的充分证据。双靶虚拟筛选应纳入实验定义的单靶硬负样本、配体层混淆对照、面板外配体评价以及受体敏感性分析。因此，DualFourClass-Bench 的主要贡献不是产生一个普适的 docking winner，而是提供系统协议，用以界定 docking-based dual-target recognition 的证据与可靠性边界。

---

## 术语分层（全文只在 Conclusion 用 grounded）

| 章节 | 推荐用语 |
|------|----------|
| Introduction | experimentally defined dual-target recognition task |
| Methods | experimentally derived activity labels |
| Results | experimentally defined hard negatives |
| Conclusion | experimentally grounded evaluation setting（仅一次） |

禁止：validated；robust performance；generalizable dual-target docking strategy；docking is ineffective；docking can identify dual-target ligands。
