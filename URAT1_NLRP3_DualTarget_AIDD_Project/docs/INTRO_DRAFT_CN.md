# 引言初稿（中文；二次去 AI 腔润色）

> 相对上一版，进一步压缩套话过渡、对称排比、“不是……而是……”句式和段末复述；保留可核实事实与引用。

痛风属于慢性炎症性关节炎。关节及周围组织中沉积的单钠尿酸盐（MSU）晶体，可通过 NLRP3 炎症小体驱动 IL-1β 释放，诱发局部炎症（Dalbeth et al., *Lancet* 2021, https://doi.org/10.1016/S0140-6736(21)00569-9；Leask et al., *Nat Rev Rheumatol* 2024, https://doi.org/10.1038/s41584-024-01137-1）。高尿酸血症是其必要前提，多数情况下来自尿酸排泄不足，生成过多相对少见。尿酸由肝脏嘌呤代谢产生，经肾脏与肠道排出；在肾近端小管，URAT1（*SLC22A12*）负责约 90% 的尿酸重吸收，因此是降尿酸药的重要靶点（Dai & Lee, *Cell Res* 2024, https://doi.org/10.1038/s41422-024-01023-1；Lin et al., *Signal Transduct Target Ther* 2024, https://doi.org/10.1038/s41392-024-01916-y）。MSU 被巨噬细胞摄取后，NLRP3 组装并激活 caspase-1，将 pro-IL-1β 加工为活性 IL-1β，这是急性发作的关键环节（Martinon et al., *Nature* 2006, https://doi.org/10.1038/nature04516；Zhao et al., *Front Immunol* 2023, https://doi.org/10.3389/fimmu.2023.1137822）。临床治疗因而同时面对两件事：把血尿酸降下来，以及把晶体诱发的炎症压住。

降尿酸仍是长期管理的基础。黄嘌呤氧化酶抑制剂（别嘌醇、非布司他）减少尿酸生成；促尿酸排泄药则抑制肾小管重吸收，其中 URAT1 抑制剂包括苯溴马隆、lesinurad、dotinurad 和仍在研究中的 verinurad。近年冷冻电镜结构显示，这些抑制剂可占据 URAT1 的底物结合区并稳定内向构象，从而阻断转运（Fedor/Suo et al., *Nat Commun* 2025, https://doi.org/10.1038/s41467-025-60480-3；Wu et al., *Cell Discov* 2025, https://doi.org/10.1038/s41421-025-00779-z）。急性炎症期常用秋水仙碱、非甾体抗炎药、糖皮质激素，重症也可使用 IL-1 相关生物制剂。问题在于现有药物并不好用尽。别嘌醇在 HLA-B\*5801 阳性患者中有严重超敏风险，东亚人群筛查需求更高；非布司他的心血管安全性在 CARES 与 FAST 之间并不完全同向，高危患者仍需权衡（Borghi et al., *Cardiol J* 2023, https://journals.viamedica.pl/cardiology_journal/article/view/98254）；苯溴马隆的肝毒性限制了使用范围，lesinurad 的肾脏安全性也使其应用受限。实践中还有不少患者长期达不到目标血尿酸，依从性差会进一步放大差距。抗炎药起效快，但多停留在对症层面；针对 NLRP3 的直接抑制剂如 MCC950 虽机制清楚，却因安全性未能进入常规治疗。降尿酸与抗炎如果简单叠加，又容易带来额外的相互作用和不良反应。对合并多种慢病的患者来说，这一矛盾尤其突出。

临床上，痛风常与其他慢病并存。心血管疾病、代谢综合征、脂肪肝、慢性肾脏病和 2 型糖尿病都很常见（Cleveland Clinic J Med 2024, https://www.ccjm.org/content/91/7/392；Front Cardiovasc Med 2023, https://doi.org/10.3389/fcvm.2023.1190069）。在这种背景下，全新活性分子从发现到可用往往耗时很长。已经进入临床阶段的化合物通常已有人体安全性和药代信息，若能从中找出与 URAT1、NLRP3 相关的新用途线索，转化成本通常低于从头做药，也更贴近真实用药环境。我们因此以临床阶段药物库做计算筛选，看其中有没有值得拿去实验验证的双节点候选。

要做这件事，先得承认两个靶点并不对称。URAT1 是多构象膜转运体（Dai & Lee, *Cell Res* 2024, https://doi.org/10.1038/s41422-024-01023-1），公开活性数据噪声不小，我们已有的回归模型对已知尿酸药回收不稳，排序不能主要靠机器学习，得靠对接，但哪种对接和打分更适合 URAT1，本身就需要先比较。NLRP3 的小分子抑制位点在 NACHT。全长蛋白体积大、常以寡聚体出现，不太适合直接拿来布对接网格，所以我们用约 2.8 Å 的 7ALV NACHT 结构；其上结合的是 MCC950 类似物 NP3-146，口袋信息可用于磺酰脲类相关分析，并作为后续分子动力学起点（Dekker et al., *J Mol Biol* 2021, https://doi.org/10.1016/j.jmb.2021.167189）。NLRP3 一侧我们更信任分类模型来缩库，对接和动力学用来检查结构是否说得通。URAT1 一侧则先按 Gu 等提出的 TrueDecoy / RandomDecoy 思路（Gu et al., *Nat Mach Intell* 2025, https://doi.org/10.1038/s42256-025-00993-0）建本靶点诱饵集，比较 Vina、gnina 和 RTMScore 等开源方案，再把表现更好的协议接到临床库筛选里。

整条流程希望把协议选择、缩库、双靶对接和成药性检查串起来，最后给出可以拿去验证的假说，而不是一份“已证实双靶药”的名单。

---

## 引用清单

1. Dalbeth N, et al. *Lancet* 2021. doi:10.1016/S0140-6736(21)00569-9
2. Leask MP, et al. *Nat Rev Rheumatol* 2024. doi:10.1038/s41584-024-01137-1
3. Dai Y, Lee CH. *Cell Res* 2024. doi:10.1038/s41422-024-01023-1
4. Lin X, et al. *Signal Transduct Target Ther* 2024. doi:10.1038/s41392-024-01916-y
5. Martinon F, et al. *Nature* 2006. doi:10.1038/nature04516
6. Zhao J, et al. *Front Immunol* 2023. doi:10.3389/fimmu.2023.1137822
7. Fedor JG/Suo Y, et al. *Nat Commun* 2025. doi:10.1038/s41467-025-60480-3
8. Wu C, et al. *Cell Discov* 2025. doi:10.1038/s41421-025-00779-z
9. Borghi C, et al. *Cardiol J* 2023 update
10. Dekker A, et al. *J Mol Biol* 2021. doi:10.1016/j.jmb.2021.167189
11. Gu S, et al. *Nat Mach Intell* 2025. doi:10.1038/s42256-025-00993-0
12. Cleveland Clinic J Med 2024;91(7):392
13. Front Cardiovasc Med 2023. doi:10.3389/fcvm.2023.1190069
