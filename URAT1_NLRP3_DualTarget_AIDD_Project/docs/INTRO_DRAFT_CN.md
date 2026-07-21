# 引言初稿（中文；原创表述）

> **叙事口径**：本研究按**事先规定的评估框架**设计（双诱饵选协议 → 不对称漏斗 → 审计提名），不是“先试某引擎失败再改路线”的补救日记。  
> **下一部分 Methods**：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)  
> **勿打开**：[`MANUSCRIPT_DRAFT_CN.md`](MANUSCRIPT_DRAFT_CN.md)（Glide XP 旧稿，已过时）

痛风是成人常见的晶体性关节炎。血尿酸长期偏高时，单钠尿酸盐（MSU）可在关节及其周围沉积；巨噬细胞识别晶体后，NLRP3 炎症小体被激活，caspase-1 将 pro-IL-1β 加工为活性 IL-1β，从而引发红肿热痛的急性发作（Dalbeth et al., *Lancet* 2021, https://doi.org/10.1016/S0140-6736(21)00569-9；Leask et al., *Nat Rev Rheumatol* 2024, https://doi.org/10.1038/s41584-024-01137-1；Martinon et al., *Nature* 2006, https://doi.org/10.1038/nature04516）。尿酸来自肝脏嘌呤代谢。多数患者高尿酸血症主要因肾脏或肠道排泄不足，其中近端小管 URAT1（*SLC22A12*）承担约九成尿酸重吸收，因此是降尿酸药的关键靶点（Dai & Lee, *Cell Res* 2024, https://doi.org/10.1038/s41422-024-01023-1；Lin et al., *Signal Transduct Target Ther* 2024, https://doi.org/10.1038/s41392-024-01916-y）。与此相应，临床管理通常要同时处理两件事：降低血尿酸，以及抑制晶体驱动的炎症（Zhao et al., *Front Immunol* 2023, https://doi.org/10.3389/fimmu.2023.1137822）。

降尿酸药物大致包括减少生成的黄嘌呤氧化酶抑制剂（别嘌醇、非布司他）和促进排泄的药物。后者中，URAT1 抑制剂如苯溴马隆、lesinurad、dotinurad 及处于临床研究中的 verinurad，可占据底物结合区并稳定内向构象，阻断尿酸重吸收（Fedor/Suo et al., *Nat Commun* 2025, https://doi.org/10.1038/s41467-025-60480-3；Wu et al., *Cell Discov* 2025, https://doi.org/10.1038/s41421-025-00779-z）。急性期抗炎则常用秋水仙碱、非甾体抗炎药、糖皮质激素，必要时使用 IL-1 相关生物制剂。现有方案仍有明显短板。别嘌醇在 HLA-B\*5801 阳性人群中有严重超敏风险；非布司他的心血管安全性在 CARES 与 FAST 之间结论并不完全一致，高危患者仍需权衡（Borghi et al., *Cardiol J* 2023, https://journals.viamedica.pl/cardiology_journal/article/view/98254）；苯溴马隆受肝毒性限制，lesinurad 也因肾脏安全性影响使用。长期达标率并不理想，依从性差会进一步拉低效果。抗炎药起效快，但多偏对症；MCC950 等 NLRP3 直接抑制剂虽机制清楚，却因安全性未能进入常规治疗。降尿酸与抗炎若简单联用，还可能增加相互作用与不良反应。

痛风还常与心血管病、代谢综合征、脂肪肝、慢性肾脏病和糖尿病等并存（Cleveland Clinic J Med 2024, https://www.ccjm.org/content/91/7/392；Front Cardiovasc Med 2023, https://doi.org/10.3389/fcvm.2023.1190069）。在这种合并症背景下，全新分子从发现到可用周期长；已进入临床阶段的化合物通常已有人体安全性与药代信息，更适合作为重定位筛选对象。我们因此以临床阶段药物库为化学空间，评估其中是否存在同时关联 URAT1 与 NLRP3、并值得实验跟进的候选，而不是事先指定某一分子作为终点。

两靶在数据结构与机制标签上并不对称，因此不宜套用“双靶统一机器学习排序”的对称流程。URAT1 是多构象膜转运体，公开活性噪声较大；回归模型若对已知尿酸药基准回收不足，则预先规定**不得**以其分数主排临床库，而改以结构对接提供 URAT1 轴证据（Dai & Lee, *Cell Res* 2024, https://doi.org/10.1038/s41422-024-01023-1）。对接排序本身亦不能默认某一引擎或某一诱饵设定：随机诱饵往往夸大富集，性质匹配诱饵更能暴露协议弱点。本研究因此把 Gu 等提出的 TrueDecoy / RandomDecoy 双诱饵评估（Gu et al., *Nat Mach Intell* 2025, https://doi.org/10.1038/s42256-025-00993-0）写入方法设计的**第一阶段**：在 URAT1 inward-open 结构（9DKB）上，按事先锁定的富集与对照回收规则，于开源搜索–打分组合（AutoDock Vina、gnina、RTMScore）中选定生产用排序协议 Π\*，再将该协议用于后续漏斗。NLRP3 一侧则以分类模型负责缩库；结构模板采用约 2.8 Å 的 NACHT 结构 7ALV（共晶配体为 MCC950 类类似物 NP3-146），对接与分子动力学仅作结构佐证，而非单独替代机器学习主轴（Dekker et al., *J Mol Biol* 2021, https://doi.org/10.1016/j.jmb.2021.167189）。

据此，本文报告一条**先验规定、可逐步复核**的双节点计算流程：双诱饵协议筛选 → NLRP3 机器学习缩库 → 双靶对接百分位与 Pareto 整合 → 成药性与结构警报审计后的假说提名。流程输出的是可证伪的计算假说与方法学对照结果，不等于已证实的双靶抑制剂，也不构成临床用药建议。

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
