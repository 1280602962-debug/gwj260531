# 引言初稿（中文；原创表述，基于本项目首次五段内容）

> **更正说明**：此前一版过多沿用他文句式，接近换皮，已废弃。  
> 本版只保留本项目首次引言的信息骨架，全部重新措辞；不套用他文原句与“天然产物双靶”路线。

痛风是成人常见的晶体性关节炎。血尿酸长期偏高时，单钠尿酸盐（MSU）可在关节及其周围沉积；巨噬细胞识别晶体后，NLRP3 炎症小体被激活，caspase-1 将 pro-IL-1β 加工为活性 IL-1β，从而引发红肿热痛的急性发作（Dalbeth et al., *Lancet* 2021, https://doi.org/10.1016/S0140-6736(21)00569-9；Leask et al., *Nat Rev Rheumatol* 2024, https://doi.org/10.1038/s41584-024-01137-1；Martinon et al., *Nature* 2006, https://doi.org/10.1038/nature04516）。尿酸来自肝脏嘌呤代谢。多数患者高尿酸血症主要因肾脏或肠道排泄不足，其中近端小管 URAT1（*SLC22A12*）承担约九成尿酸重吸收，因此是降尿酸药的关键靶点（Dai & Lee, *Cell Res* 2024, https://doi.org/10.1038/s41422-024-01023-1；Lin et al., *Signal Transduct Target Ther* 2024, https://doi.org/10.1038/s41392-024-01916-y）。与此相应，临床管理通常要同时处理两件事：降低血尿酸，以及抑制晶体驱动的炎症（Zhao et al., *Front Immunol* 2023, https://doi.org/10.3389/fimmu.2023.1137822）。

降尿酸药物大致包括减少生成的黄嘌呤氧化酶抑制剂（别嘌醇、非布司他）和促进排泄的药物。后者中，URAT1 抑制剂如苯溴马隆、lesinurad、dotinurad 及处于临床研究中的 verinurad，可占据底物结合区并稳定内向构象，阻断尿酸重吸收（Fedor/Suo et al., *Nat Commun* 2025, https://doi.org/10.1038/s41467-025-60480-3；Wu et al., *Cell Discov* 2025, https://doi.org/10.1038/s41421-025-00779-z）。急性期抗炎则常用秋水仙碱、非甾体抗炎药、糖皮质激素，必要时使用 IL-1 相关生物制剂。现有方案仍有明显短板。别嘌醇在 HLA-B\*5801 阳性人群中有严重超敏风险；非布司他的心血管安全性在 CARES 与 FAST 之间结论并不完全一致，高危患者仍需权衡（Borghi et al., *Cardiol J* 2023, https://journals.viamedica.pl/cardiology_journal/article/view/98254）；苯溴马隆受肝毒性限制，lesinurad 也因肾脏安全性影响使用。长期达标率并不理想，依从性差会进一步拉低效果。抗炎药起效快，但多偏对症；MCC950 等 NLRP3 直接抑制剂虽机制清楚，却因安全性未能进入常规治疗。降尿酸与抗炎若简单联用，还可能增加相互作用与不良反应。

痛风还常与心血管病、代谢综合征、脂肪肝、慢性肾脏病和糖尿病等并存（Cleveland Clinic J Med 2024, https://www.ccjm.org/content/91/7/392；Front Cardiovasc Med 2023, https://doi.org/10.3389/fcvm.2023.1190069）。在这种合并症背景下，全新分子从发现到可用周期长；已进入临床阶段的化合物通常已有人体安全性与药代信息，更适合作为重定位筛选对象。我们因此以临床阶段药物库开展计算筛选，评估其中是否存在同时关联 URAT1 与 NLRP3、并值得实验验证的候选，而不是事先指定某一分子。

本项目要建立的是一条可重复核查的 URAT1–NLRP3 双节点计算流程，输出可检验假说，而不是已经证实的双靶抑制剂。两靶条件并不对称。URAT1 是多构象膜转运体，公开活性噪声较大，机器学习回归对已知尿酸药回收不稳，排序应主要依靠对接，但需先比较哪套搜索与打分更合适（Dai & Lee, *Cell Res* 2024, https://doi.org/10.1038/s41422-024-01023-1）。我们按 Gu 等提出的 TrueDecoy / RandomDecoy 思路（Gu et al., *Nat Mach Intell* 2025, https://doi.org/10.1038/s42256-025-00993-0）建立本靶点诱饵集，比较 Vina、gnina 和 RTMScore 等开源方案，再把较好的协议接到后续筛选。NLRP3 小分子位点在 NACHT。全长蛋白体积大、常呈寡聚体，不便直接布网格，故采用约 2.8 Å 的 7ALV 结构；其上为 MCC950 类类似物 NP3-146，可作为对接与分子动力学模板（Dekker et al., *J Mol Biol* 2021, https://doi.org/10.1016/j.jmb.2021.167189）。NLRP3 一侧以分类模型缩库，对接与动力学用于检查结构是否合理。

整条流程将协议选择、缩库、双靶对接和成药性审查串起来，形成可拿去验证的假说列表。后续仍需实验，不能等同于双靶药发现或临床用药建议。

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
