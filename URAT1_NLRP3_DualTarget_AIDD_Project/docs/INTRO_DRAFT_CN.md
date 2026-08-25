# 引言初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)。Results：[`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)。Discussion：[`DISCUSSION_DRAFT_CN.md`](DISCUSSION_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用"作者，年份"，DOI 见文末清单。

痛风是成人最常见的晶体性关节炎之一，其发生发展与持续性高尿酸血症及单钠尿酸（monosodium urate，MSU）晶体沉积密切相关。多数患者的高尿酸血症并非源于尿酸生成过多，而主要归因于肾脏尿酸排泄不足（Dalbeth 等，2021）。当血尿酸水平长期升高并超过尿酸盐溶解度阈值时，MSU 晶体可在关节及其周围组织中形成和沉积；晶体形成后被巨噬细胞等先天免疫细胞识别，诱导 NLRP3 炎症小体活化并启动炎症级联反应，最终导致急性痛风性关节炎发作（Martinon 等，2006；Dalbeth 等，2021；Leask 等，2024）。可见痛风的病理过程并非单一环节所致，而是呈现由尿酸稳态失衡、MSU 晶体形成向晶体驱动的炎症反应连续发展的病理链条。介导肾近端小管尿酸重吸收的关键转运体 URAT1 位于该链条上游的尿酸稳态调控环节，影响高尿酸血症程度及后续晶体形成风险；NLRP3 则位于下游晶体诱导炎症环节，介导 MSU 晶体所触发的炎症反应。二者分别处于痛风病理过程的尿酸稳态调控与晶体驱动炎症两个阶段，构成具有互补性的干预节点，本文以下统称为"双节点"（Dalbeth 等，2019；Liu 等，2023）。

现行痛风药物治疗正是围绕这两个节点分别展开。降尿酸一端，黄嘌呤氧化酶抑制剂减少尿酸生成，URAT1 抑制剂促进尿酸排泄，包括苯溴马隆、lesinurad、dotinurad 以及仍处临床试验阶段的 verinurad（FitzGerald 等，2020；Dai 与 Lee，2024）；抗炎一端，急性发作以秋水仙碱、非甾体抗炎药或糖皮质激素为一线，白细胞介素-1 抑制剂仅在上述药物无效、不耐受或禁忌时使用（FitzGerald 等，2020）。然而两端在治疗时程内的作用范围并不重叠：降尿酸治疗可降低尿酸负荷，却对已沉积 MSU 晶体所诱导的急性炎症缺乏直接、快速的消解作用；抗炎治疗能控制发作期炎症，却不降低尿酸负荷，也不清除已沉积的晶体（Dalbeth 等，2021）。二者还各自受到安全性与适用范围的限制：别嘌醇存在 HLA-B\*5801 相关严重皮肤不良反应风险（Hung 等，2005），非布司他在高心血管风险人群中的安全性仍需谨慎评估（White 等，2018；Mackenzie 等，2020；Borghi 等，2024），苯溴马隆与 lesinurad 则分别受肝脏和肾脏安全性问题制约（Lee 等，2008；Tausche 等，2017）。

正因两端所覆盖的病理环节不同，指南在启动降尿酸治疗时建议同时给予抗炎预防、疗程至少 3–6 个月（FitzGerald 等，2020）。这意味着在相当长的时间窗内，患者需并行接受两类作用机制不同的药物，这一共给药情形为同时覆盖两个病理环节的双节点分子提供了潜在的应用场景。痛风患者常合并心血管疾病、代谢性疾病与慢性肾脏病（Dalbeth 等，2021；Du 等，2024），本身多处于多重用药状态，联合方案会进一步叠加药物相互作用风险——例如秋水仙碱经 CYP3A4 与 P-糖蛋白代谢转运，与他汀类、大环内酯类等常用药物存在明确相互作用（FitzGerald 等，2020）。若单一分子能够同时降低尿酸负荷并抑制晶体驱动的炎症，则理论上可在部分治疗阶段以更少的药物条目覆盖两个病理环节，从而具有简化药理干预、减少联合用药相互作用暴露的潜在价值。需要说明的是，这一价值属于治疗方案层面的潜在优势，并不意味着双节点分子在疗效上优于现行联合治疗。

这样的分子在化学层面并非不可企及。近期已有研究通过天然产物衍生与结构优化，获得可同时影响尿酸盐转运与 NLRP3 炎症通路的双作用分子，为该策略提供了实验层面的可行性先例（Zhang 等，2025）。另一方面，已报道的 URAT1 抑制剂呈现较为多样的化学骨架，其中包括若干并非为该靶点设计的临床药物（Enomoto 等，2002）；近期冷冻电镜结构进一步显示，骨架差异显著的多类抑制剂可结合于同一中央配体结合腔并稳定 inward-open 构象（Wu 等，2025；Suo 等，2025）。这提示在既有药物化学空间中检索 URAT1 结合线索具有现实可能性，也为以结构为基础的筛选提供了实验依据。

在实现路径上，从头开发全新双靶分子面临周期长、淘汰率高的固有困难，NLRP3 一端的转化经验尤其说明了这一点。尽管 NLRP3 抑制具有明确的抗炎药理基础（Coll 等，2015），专门针对该靶点开发的小分子仍在安全性与疗效验证上受挫——NLRP3 抑制剂 GDC-2394 即在首次人体研究中因严重肝损伤事件终止开发（Tang 等，2023），迄今尚无 NLRP3 小分子抑制剂成为痛风的常规治疗药物。相较之下，在已进入临床研究或已上市的小分子中检索双节点线索构成一条风险结构不同的补充路径：其中相当一部分分子已积累人体药代动力学、暴露或安全性信息，可在一定程度上降低成药性与制剂层面的不确定性，并为后续实验验证提供更明确的开发背景（Pushpakom 等，2019）。需要说明的是，已有人体暴露数据并不等同于消除靶点相关毒性风险，GDC-2394 的经历正说明这一点；重定位所降低的主要是药代与成药性维度的不确定性，具体安全性结论仍需在目标适应证与相应暴露水平下重新评估。

据此，本研究聚焦于一个具体的转化问题：在已进入临床研究或已上市的小分子中，是否存在能够同时获得 URAT1 结构结合支持、以及 NLRP3 配体基活性证据与结构证据支持的候选分子？由于 URAT1 与 NLRP3 两侧的公开活性数据在数据规模、测定条件与化学空间覆盖方面并不对等，两侧无法采用同一套筛选策略。URAT1 侧公开活性数据在化学空间覆盖上相对有限，与本研究所用临床重定位库之间存在明显的化学空间差异（该差异在结果部分予以定量刻画），单纯依赖数据驱动模型存在较高的外推风险，因此本研究以近期冷冻电镜结构为基础的分子对接作为该侧的主要筛选依据；NLRP3 侧公开活性数据规模较大，故先以配体基优先级模型对临床库进行缩减与排序，再对缩减后的分子进行结构对接以获得独立的结构层面证据，其中机器学习输出仅用于候选缩库与优先级排序，不作为活性确证。

考虑到分子对接评分与单一预测模型均存在固有偏差，本研究在生产性筛选前先构建 URAT1 活性物–诱饵基准集，在性质匹配诱饵与随机诱饵两种条件下比较多套对接及重打分协议的富集能力，据此锁定正式筛选所用协议。随后将经 NLRP3 模型缩减后的临床阶段及已上市小分子分别对 URAT1 与 NLRP3 进行结构对接，并将各证据来源统一转换为池内百分位，整合为 URAT1 结构证据与 NLRP3 综合证据两条可比较的评价轴；同时纳入已知 URAT1 抑制剂与抗炎药物作为内部对照，用以审计两条证据轴的实际区分能力。在此基础上通过 Pareto 非支配分析考察双轴证据分布，并进一步结合结构警报、类药性、分子量与骨架去冗余等药物化学标准提出优先候选。需要强调的是，对接协议比较若达不到多样化库可用门槛，则后续百分位名单不能写成双节点候选鉴定；活性与作用机制只能由与排名脱钩的实验回答。

---

## 引用清单（定稿时改为期刊格式）

1. Dalbeth N, Choi HK, Joosten LAB, et al. Gout. *Lancet.* 2021;397(10287):1843-1855. doi:10.1016/S0140-6736(21)00569-9
2. Martinon F, Pétrilli V, Mayor A, Tardivel A, Tschopp J. Gout-associated uric acid crystals activate the NALP3 inflammasome. *Nature.* 2006;440(7081):237-241. doi:10.1038/nature04516
3. Leask MP, et al. The pathogenesis of gout: molecular insights from genetic, epigenomic and transcriptomic studies. *Nat Rev Rheumatol.* 2024. doi:10.1038/s41584-024-01137-1
4. Dalbeth N, Choi HK, Joosten LAB, et al. Gout. *Nat Rev Dis Primers.* 2019;5:69. doi:10.1038/s41572-019-0115-y
4a. Liu Y-r, Wang J-q, Li J. Role of NLRP3 in the pathogenesis and treatment of gout arthritis. *Front Immunol.* 2023;14:1137822. doi:10.3389/fimmu.2023.1137822
5. FitzGerald JD, Dalbeth N, Mikuls T, et al. 2020 American College of Rheumatology guideline for the management of gout. *Arthritis Rheumatol.* 2020;72(6):879-895. doi:10.1002/art.41247
6. Dai Y, Lee CH. Transport mechanism and structural pharmacology of human urate transporter URAT1. *Cell Res.* 2024. doi:10.1038/s41422-024-01023-1
7. Hung SI, Chung WH, Liou LB, et al. HLA-B\*5801 allele as a genetic marker for severe cutaneous adverse reactions caused by allopurinol. *Proc Natl Acad Sci USA.* 2005;102(11):4134-4139. doi:10.1073/pnas.0409500102
8. White WB, Saag KG, Becker MA, et al. Cardiovascular safety of febuxostat or allopurinol in patients with gout (CARES). *N Engl J Med.* 2018;378(13):1200-1210. doi:10.1056/NEJMoa1710895
9. Mackenzie IS, Ford I, Nuki G, et al. Long-term cardiovascular safety of febuxostat compared with allopurinol in patients with gout (FAST). *Lancet.* 2020;396(10264):1745-1757. doi:10.1016/S0140-6736(20)32234-0
10. Borghi C, Domienik-Karłowicz J, Tykarski A, et al. Expert consensus for the diagnosis and treatment of patients with hyperuricemia and high cardiovascular risk: 2023 update. *Cardiol J.* 2024;31(1):1-14. doi:10.5603/cj.98254
11. Lee MH, Graham GG, Williams KM, Day RO. A benefit-risk assessment of benzbromarone in the treatment of gout. *Drug Saf.* 2008;31(8):643-665. doi:10.2165/00002018-200831080-00002
12. Tausche AK, Alten R, Dalbeth N, et al. Lesinurad monotherapy in gout patients intolerant to a xanthine oxidase inhibitor. *Rheumatology (Oxford).* 2017;56(12):2170-2178. doi:10.1093/rheumatology/kex305
13. Du L, Zong Y, Li H, et al. Hyperuricemia and its related diseases: mechanisms and advances in therapy. *Signal Transduct Target Ther.* 2024;9:212. doi:10.1038/s41392-024-01916-y
14. Zhang Z, Shi X, Wu T, et al. Discovery of multi-target anti-gout agents from *Eurycoma longifolia* Jack. *Nat Commun.* 2025;16:7430. doi:10.1038/s41467-025-62645-6
15. Enomoto A, Kimura H, Chairoungdua A, et al. Molecular identification of a renal urate–anion exchanger that regulates blood urate levels. *Nature.* 2002;417(6887):447-452. doi:10.1038/nature742
16. Wu C, Zhang C, Jin S, et al. Molecular mechanisms of urate transport by the native human URAT1 and its inhibition by anti-gout drugs. *Cell Discov.* 2025;11:33. doi:10.1038/s41421-025-00779-z
17. Suo Y, Fedor JG, Zhang H, et al. *Nat Commun.* 2025;16:5178. doi:10.1038/s41467-025-60480-3 （Suo 与 Fedor 同等贡献；期刊作 Suo et al.）
18. Coll RC, Robertson AAB, Chae JJ, et al. A small-molecule inhibitor of the NLRP3 inflammasome for the treatment of inflammatory diseases. *Nat Med.* 2015;21(3):248-255. doi:10.1038/nm.3806
19. Tang F, Kunder R, Chu T, et al. First-in-human phase 1 trial of NLRP3 inhibitor GDC-2394. *Clin Transl Sci.* 2023;16(9):1653-1666. doi:10.1111/cts.13576
20. Pushpakom S, Iorio F, Eyers PA, et al. Drug repurposing: progress, challenges and recommendations. *Nat Rev Drug Discov.* 2019;18(1):41-58. doi:10.1038/nrd.2018.168
