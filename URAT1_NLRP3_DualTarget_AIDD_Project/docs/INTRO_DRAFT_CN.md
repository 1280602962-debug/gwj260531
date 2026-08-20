# 引言初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用“作者，年份”，DOI 见文末清单。

痛风是成人最常见的晶体性关节炎之一。当血尿酸长期升高，单钠尿酸盐可在关节及其周围沉积；巨噬细胞识别晶体后激活 NLRP3 炎症小体，经 caspase-1 将前体白细胞介素-1β 加工为活性形式，从而触发急性炎症发作（Dalbeth 等，2021；Leask 等，2024；Martinon 等，2006）。尿酸主要来自肝脏嘌呤代谢。多数患者的高尿酸血症并非生成过多，而是肾脏或肠道排泄不足；近端小管尿酸盐转运体 URAT1（由 *SLC22A12* 编码）承担约九成尿酸重吸收，因而是降尿酸治疗的关键靶点（Dai 与 Lee，2024；Lin 等，2024）。临床管理往往需要同时应对两个相互关联却机制不同的节点：降低血尿酸，以及抑制晶体驱动的炎症（Liu 等，2023）。

现有药物分降尿酸和控制发作两条路径。降尿酸一端，黄嘌呤氧化酶抑制剂减少尿酸生成；促排泄一端则用 URAT1 抑制剂，包括苯溴马隆、lesinurad、dotinurad，以及仍在临床试验中的 verinurad（FitzGerald 等，2020；Dai 与 Lee，2024）。冷冻电镜显示，这四种抑制剂结合在尿酸盐所在的中央口袋，并把转运体稳定在内向开放构象（Wu 等，2025）；苯溴马隆与 lesinurad 的内向开放结合在另一项结构研究中也得到证实（Fedor 等，2025）。急性发作的一线药物是秋水仙碱、非甾体抗炎药或糖皮质激素；仅当这些药物无效、不耐受或禁忌时，才考虑白细胞介素-1 抑制剂（FitzGerald 等，2020）。启动降尿酸治疗时，指南还建议短期内加用抗炎药物，以减少诱发发作（FitzGerald 等，2020）。

各药的限制并不相同。别嘌醇在 HLA-B*5801 阳性个体中与严重皮肤不良反应密切相关（Hung 等，2005）。非布司他在已有心血管疾病的痛风患者中，CARES 试验报告心血管死亡高于别嘌醇（White 等，2018）；其后 FAST 试验未再现这一信号（Mackenzie 等，2020），但对高心血管风险人群，共识仍主张谨慎使用（Borghi 等，2024）。苯溴马隆因肝毒性在多地限制使用（Lee 等，2008）。lesinurad 单药治疗时肾相关不良事件明显增加，不能单独给药（Tausche 等，2017）。控制发作的抗炎药并不降低血尿酸，也不能溶解已沉积的晶体（Dalbeth 等，2021）。MCC950 可选择性抑制 NLRP3（Coll 等，2015），但 II 期试验因肝酶升高终止开发（Tang 等，2024），目前没有这类小分子进入痛风常规治疗。痛风还常与心血管病、代谢病和慢性肾脏病并存（Dalbeth 等，2021；Lin 等，2024）。

把两个节点分开处理的现有药物各自受安全性约束，而新化学实体从发现到临床可用周期长、淘汰率高。痛风双节点因此可以先在已经进入人体研究的小分子中，寻找同时具备 URAT1 结合与 NLRP3 相关活性的计算线索。这类分子通常已有药代和安全性信息，后续推进的成本和周期都低于从零开始的新实体，因而常被用作重定位的起点（Pushpakom 等，2019）。本文据此把筛选对象限定为 ChEMBL 中已进入临床阶段或已经上市的小分子。

不选用尚未进入人体研究的大规模商品库，是因为那一类筛选回答的是还有没有新骨架，而不是现有临床阶段分子里有没有可检验的双节点线索。天然产物衍生那条路则已有湿法先例：从长叶阔蕊苏出发的工作表明，经过结构优化的同一分子可以同时影响尿酸盐转运和 NLRP3 介导的炎症（Zhang 等，2025）；那一项研究回答的是新实体能否做成双作用分子，本文不再沿同一路线重复。要把同样的双节点逻辑用到上述临床阶段集合上，两边的公开证据却并不对等，筛选方法也不能直接照搬。

URAT1 是交替开放的膜转运体，已知抑制剂结合在 inward-open 构象的底物口袋中（Dai 与 Lee，2024；Wu 等，2025；Fedor 等，2025）。该靶点在 ChEMBL 中的活性来自不同测定，相互冲突的记录在清洗中被剔除后训练覆盖变窄，若干已上市尿酸药的骨架也不在其中，因而不能默认回归模型可以外推到这批临床阶段分子，更不宜用模型分去给它们排序，URAT1 一侧的排序应交给结构对接。对接读出本身也依赖于诱饵如何构成：Gu 等表明，同一套对接或重打分方法在实验无活构成的难诱饵集与随机商业库诱饵集上可以给出很不一样的虚拟筛选表现（Gu 等，2025），生产协议便只能先在独立的活性物–诱饵基准上比较，而不能拿这批分子在对接后的排名回头挑选方法。NLRP3 的公开测定条件同样不均一，但化合物数量足以训练分类器，可以先用来缩小待对接集合；分类分数只说明分子在训练测定语境下更像活性物，不足以单独决定跟进名单。NLRP3 结构分析采用 NACHT 结构 7ALV，共晶配体是 MCC950 类类似物 NP3-146，而不是 MCC950 本身（Dekker 等，2021）。

因此本文先在这批临床阶段分子之外完成 URAT1 对接协议的选择。参照 Gu 等的 TrueDecoy 与 RandomDecoy 思路，我们在单靶 URAT1 上分别构建实验弱活加性质匹配的难诱饵、以及与之等规模且无重叠的随机诱饵，于 lesinurad 复合的 inward-open 结构 9DKB（Fedor 等，2025）上比较开源对接与重打分，并按事先写明的规则锁定生产读出。协议选定之后，才用 NLRP3 分类模型压缩该临床阶段集合，将入选分子按同一协议分别对接至 9DKB 与 7ALV，再把结构分转为池内百分位。Pareto 非支配集合只用来查看哪些分子在两轴上不被支配；跟进短名单另要求两靶对接百分位同时过门，并做结构警报、类药性和骨架去冗余。NLRP3 模型分只负责缩库，不单独把分子送进短名单。

本文的目的是在上述证据不对称的条件下，从已进入人体研究的小分子中提出可供实验检验的 URAT1–NLRP3 双节点假说，并说明对接协议如何选定、名单如何审计。对接分数只在同一池内比较高低，不解释为亲和力；文中提名也不等于已经证实的双靶抑制剂，更不构成临床用药建议。

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
20. Gu S, Zhang X, Shen C, et al. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat Mach Intell.* 2025. doi:10.1038/s42256-025-00993-0  
21. Dekker C, Mattes H, Wright M, et al. Crystal structure of NLRP3 NACHT domain with an inhibitor defines mechanism of inflammasome inhibition. *J Mol Biol.* 2021;433(24):167309. doi:10.1016/j.jmb.2021.167309  
