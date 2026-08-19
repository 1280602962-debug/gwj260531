# 引言初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用“作者，年份”，DOI 见文末清单。

痛风是成人最常见的晶体性关节炎之一。当血尿酸长期升高，单钠尿酸盐可在关节及其周围沉积；巨噬细胞识别晶体后激活 NLRP3 炎症小体，经 caspase-1 将前体白细胞介素-1β 加工为活性形式，从而触发急性炎症发作（Dalbeth 等，2021；Leask 等，2024；Martinon 等，2006）。尿酸主要来自肝脏嘌呤代谢。多数患者的高尿酸血症并非生成过多，而是肾脏或肠道排泄不足；近端小管尿酸盐转运体 URAT1（由 *SLC22A12* 编码）承担约九成尿酸重吸收，因而是降尿酸治疗的关键靶点（Dai 与 Lee，2024；Lin 等，2024）。临床管理往往需要同时应对两个相互关联却机制不同的节点：降低血尿酸，以及抑制晶体驱动的炎症（Liu 等，2023）。

现有药物分降尿酸和控制发作两条路径。降尿酸一端，黄嘌呤氧化酶抑制剂减少尿酸生成；促排泄一端则用 URAT1 抑制剂，包括苯溴马隆、lesinurad、dotinurad，以及仍在临床试验中的 verinurad（FitzGerald 等，2020；Dai 与 Lee，2024）。冷冻电镜显示，这四种抑制剂结合在尿酸盐所在的中央口袋，并把转运体稳定在内向开放构象（Wu 等，2025）；苯溴马隆与 lesinurad 的内向开放结合在另一项结构研究中也得到证实（Fedor 等，2025）。急性发作的一线药物是秋水仙碱、非甾体抗炎药或糖皮质激素；仅当这些药物无效、不耐受或禁忌时，才考虑白细胞介素-1 抑制剂（FitzGerald 等，2020）。启动降尿酸治疗时，指南还建议短期内加用抗炎药物，以减少诱发发作（FitzGerald 等，2020）。

各药的限制并不相同。别嘌醇在 HLA-B*5801 阳性个体中与严重皮肤不良反应密切相关（Hung 等，2005）。非布司他在已有心血管疾病的痛风患者中，CARES 试验报告心血管死亡高于别嘌醇（White 等，2018）；其后 FAST 试验未再现这一信号（Mackenzie 等，2020），但对高心血管风险人群，共识仍主张谨慎使用（Borghi 等，2024）。苯溴马隆因肝毒性在多地限制使用（Lee 等，2008）。lesinurad 单药治疗时肾相关不良事件明显增加，不能单独给药（Tausche 等，2017）。控制发作的抗炎药并不降低血尿酸，也不能溶解已沉积的晶体（Dalbeth 等，2021）。MCC950 可选择性抑制 NLRP3（Coll 等，2015），但 II 期试验因肝酶升高终止开发（Tang 等，2024），目前没有这类小分子进入痛风常规治疗。痛风还常与心血管病、代谢病和慢性肾脏病并存（Dalbeth 等，2021；Lin 等，2024）。

全新化学实体从发现到临床可用周期长、淘汰率高；已进入人体研究的化合物通常已有药代和安全性信息，因此常被用作重定位的起点（Pushpakom 等，2019）。从长叶阔蕊苏出发的湿法双作用化学表明，同一分子可以同时触及尿酸盐转运与 NLRP3 炎症（Zhang 等，2025），但其终点是新实体的合成与动物验证。它并不回答：在已有临床小分子空间里，如何在两靶公开数据不对称的条件下，按事先规定的规则筛出值得实验跟进的假说。

计算上，URAT1 与 NLRP3 不能当作同质任务。URAT1 是多构象膜转运体，公开活性数据噪声较大；若回归模型对已知尿酸药的基准回收不足，则不宜再用其分数主排临床库，而应以结构对接提供 URAT1 一侧的排序证据。即便采用对接，也不能默认某一引擎或某一诱饵设定：随机负样本往往夸大表观富集，含实验弱活种子与性质匹配分子的难诱饵集更能检验排序是否稳健。因此，本研究将 Gu 等提出的 TrueDecoy 与 RandomDecoy 双诱饵框架（Gu 等，2025）作为第一阶段：在 inward-open URAT1（PDB 9DKB）上比较开源对接与重打分方案，按事先固定的规则选定生产排序协议，再进入临床库筛选。NLRP3 一侧测定异质但样本量足以支撑监督学习，故以分类模型缩小化学空间；结构分析采用 NACHT 结构（PDB 7ALV，共晶配体为 MCC950 类类似物 NP3-146）（Dekker 等，2021）。机器学习用于缩库；进入跟进短名单时，两靶均须达到结构百分位门控，以免仅靠模型分抬进名单。这与以激酶为数据丰富臂、以固定等权融合商业库命中的不对称虚拟筛选不同：本文对象是临床阶段库与转运体口袋，对接协议先经双诱饵选择，并且把 Pareto 非支配与药物化学提名分成两层决策。

本文据此报告一条事先规定的 URAT1–NLRP3 双节点计算流程：对接协议筛选、NLRP3 机器学习缩库、双靶对接百分位与 Pareto 整合，以及结构警报、类药性与骨架去冗余审计后的假说提名，并对少数提名分子做姿态稳定性分子动力学（URAT1 按膜蛋白体系）。对接分数仅在池内转换为百分位，不解释为物理亲和力。输出是可证伪的计算假说与方法学对照，并不等同于已经证实的双靶抑制剂，也不构成临床用药建议。

---

## 引用清单（定稿时改为期刊格式）

1. Dalbeth N, et al. *Lancet* 2021. doi:10.1016/S0140-6736(21)00569-9  
2. Leask MP, et al. *Nat Rev Rheumatol* 2024. doi:10.1038/s41584-024-01137-1  
3. Martinon F, et al. *Nature* 2006. doi:10.1038/nature04516  
4. Dai Y, Lee CH. *Cell Res* 2024. doi:10.1038/s41422-024-01023-1  
5. Lin X, et al. *Signal Transduct Target Ther* 2024. doi:10.1038/s41392-024-01916-y  
6. Liu Y, Li W, Deng Y. Role of NLRP3 in the pathogenesis and treatment of gout arthritis. *Front Immunol.* 2023;14:1137822. doi:10.3389/fimmu.2023.1137822  
7. FitzGerald JD, Dalbeth N, Mikuls T, et al. 2020 American College of Rheumatology guideline for the management of gout. *Arthritis Rheumatol.* 2020;72(6):879-895. doi:10.1002/art.41247  
8. Fedor JG, Suo Y, Zhang H, et al. *Nat Commun.* 2025;16:5178. doi:10.1038/s41467-025-60480-3  
9. Wu C, Zhang C, Jin S, et al. *Cell Discov.* 2025;11:33. doi:10.1038/s41421-025-00779-z  
10. Hung SI, Chung WH, Liou LB, et al. HLA-B*5801 allele as a genetic marker for severe cutaneous adverse reactions caused by allopurinol. *Proc Natl Acad Sci USA.* 2005;102(11):4134-4139. doi:10.1073/pnas.0409500102  
11. White WB, Saag KG, Becker MA, et al. Cardiovascular safety of febuxostat or allopurinol in patients with gout (CARES). *N Engl J Med.* 2018;378(13):1200-1210. doi:10.1056/NEJMoa1710895  
12. Mackenzie IS, Ford I, Nuki G, et al. Long-term cardiovascular safety of febuxostat compared with allopurinol in patients with gout (FAST). *Lancet.* 2020;396(10264):1745-1757. doi:10.1016/S0140-6736(20)32234-0  
13. Borghi C, Domienik-Karłowicz J, Tykarski A, et al. Expert consensus for the diagnosis and treatment of patients with hyperuricemia and high cardiovascular risk: 2023 update. *Cardiol J.* 2024;31(1):1-14. doi:10.5603/cj.98254  
14. Lee MH, Graham GG, Williams KM, Day RO. A benefit-risk assessment of benzbromarone in the treatment of gout. *Drug Saf.* 2008;31(8):643-665. doi:10.2165/00002018-200831080-00002  
15. Tausche AK, Alten R, Dalbeth N, et al. Lesinurad monotherapy in gout patients intolerant to a xanthine oxidase inhibitor. *Rheumatology (Oxford).* 2017;56(12):2170-2178. doi:10.1093/rheumatology/kex305  
16. Coll RC, Robertson AAB, Chae JJ, et al. A small-molecule inhibitor of the NLRP3 inflammasome for the treatment of inflammatory diseases. *Nat Med.* 2015;21(3):248-255. doi:10.1038/nm.3806  
17. Tang F, Kunder R, Chu T, et al. First-in-human phase 1 trial of NLRP3 inhibitor GDC-2394. *Clin Transl Sci.* 2024;17(1):e13576. doi:10.1111/cts.13576 （文中明确：MCC950/CP-456,773 的 II 期因肝酶升高终止）  
18. Pushpakom S, Iorio F, Eyers PA, et al. Drug repurposing: progress, challenges and recommendations. *Nat Rev Drug Discov.* 2019;18(1):41-58. doi:10.1038/nrd.2018.168  
19. Zhang H, et al. Discovery of multi-target anti-gout agents from *Eurycoma longifolia* Jack. *Nat Commun.* 2025;16:7430. doi:10.1038/s41467-025-62645-6  
20. Gu S, et al. *Nat Mach Intell* 2025. doi:10.1038/s42256-025-00993-0  
21. Dekker A, et al. *J Mol Biol.* 2021;433(20):167189. doi:10.1016/j.jmb.2021.167189  
