# 引言初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用“作者，年份”，DOI 见文末清单。

痛风是成人最常见的晶体性关节炎之一。当血尿酸长期升高，单钠尿酸盐可在关节及其周围沉积；巨噬细胞识别晶体后激活 NLRP3 炎症小体，经 caspase-1 将前体白细胞介素-1β 加工为活性形式，从而触发急性炎症发作（Dalbeth 等，2021；Leask 等，2024；Martinon 等，2006）。多数高尿酸血症患者的尿酸升高主要与肾脏或肠道排泄不足有关，而非生成过多。肾脏尿酸重吸收由多个转运体共同完成，其中近端小管转运体 URAT1（由 *SLC22A12* 编码）是排泄侧的关键节点之一，GLUT9/SLC2A9 与 ABCG2 等也参与尿酸处理（Dai 与 Lee，2024；Suo 等，2025）。因此，URAT1 与 NLRP3 并非同一分子信号通路中的直接上下游蛋白，而是分别位于尿酸稳态和晶体驱动炎症两个相互衔接的病理节点（Liu 等，2023）。

现有药物也分别作用于这两个节点（FitzGerald 等，2020）。降尿酸一端，黄嘌呤氧化酶抑制剂减少尿酸生成；促排泄一端则用 URAT1 抑制剂，包括苯溴马隆、lesinurad、dotinurad，以及仍在临床试验中的 verinurad（Dai 与 Lee，2024；Wu 等，2025）。近年 URAT1 冷冻电镜结构显示，多种抑制剂结合于其中央配体结合腔体并稳定 inward-open 构象，为该靶点的结构基础筛选提供了实验依据（Wu 等，2025；Suo 等，2025）。急性发作的一线药物是秋水仙碱、非甾体抗炎药或糖皮质激素；仅当这些药物无效、不耐受或禁忌时，才考虑白细胞介素-1 抑制剂（FitzGerald 等，2020）。启动降尿酸治疗时，指南还建议同时给予抗炎预防、疗程至少 3–6 个月，以减少诱发发作（FitzGerald 等，2020）。这一临床实践也反映出，降尿酸与抗炎治疗分别针对病理链条上的不同阶段，临床管理往往需要对二者进行协同而又相对独立的药理干预。

两类治疗的作用范围并不完全重叠。降尿酸治疗能够降低尿酸负荷，但对已经沉积的 MSU 晶体所诱导的急性炎症没有直接、快速的消除作用；抗炎治疗能够控制发作期炎症，却不降低尿酸负荷，也不清除已沉积的晶体（Dalbeth 等，2021）。二者各自还受到药物安全性及适用范围的限制。别嘌醇存在 HLA-B\*5801 相关严重皮肤不良反应风险（Hung 等，2005）；非布司他在高心血管风险患者中的安全性仍需谨慎评估（White 等，2018；Mackenzie 等，2020；Borghi 等，2024）；苯溴马隆和 lesinurad 也分别受到肝脏或肾脏安全性问题限制（Lee 等，2008；Tausche 等，2017）。另一方面，尽管 NLRP3 抑制具有明确的抗炎药理基础（Coll 等，2015），相关小分子的临床转化仍面临安全性和疗效验证等挑战。例如，NLRP3 小分子抑制剂 GDC-2394 在首次人体研究中曾因严重肝损伤事件而终止开发（Tang 等，2023），目前尚无 NLRP3 小分子抑制剂成为痛风的常规治疗药物。痛风还常与心血管病、代谢病和慢性肾脏病并存（Dalbeth 等，2021；Du 等，2024）。

上述分工提示，现行治疗需要通过不同药物分别控制尿酸负荷与晶体炎症；这提示同时作用于病理链条上游尿酸负荷和下游晶体炎症环节的分子，理论上可能具有药理学价值。但这并不等于双节点分子必然优于现有的降尿酸联合抗炎方案，本文也不检验这一优效性问题。真正尚不清楚的是：在已经进入人体研究的小分子中，是否存在同时具备 URAT1 结构结合证据与 NLRP3 相关活性线索的双节点候选分子。

从头设计新化学实体来回答这一问题，周期长、淘汰率高。药物重定位转而在已经进入人体研究的分子中寻找新用途：这类分子可在一定程度上利用已有的药代、安全性或临床开发信息，降低部分早期成药性风险，为后续实验验证提供更明确的转化背景（Pushpakom 等，2019）。本文据此把筛选对象限定为 ChEMBL 中已进入临床阶段或已经上市的小分子（具体纳入口径见 Methods）。与面向新化学空间的大规模虚拟筛选不同，本研究关注一个更受约束的问题：在已有一定人体研究或临床开发信息的这批分子中，能否找到可进一步验证的 URAT1–NLRP3 双节点候选。近期已有研究通过天然产物衍生和结构优化获得同时影响尿酸盐转运与 NLRP3 炎症的双作用分子，证明了这一药理组合具有实验可行性；但该路线聚焦于新化学实体的发现，而本研究关注的是已进入人体研究小分子中的潜在双节点重定位机会（Zhang 等，2025）。

要在这批分子上核验上述线索，URAT1 与 NLRP3 两侧的公开证据并不对等。由于 URAT1 公开活性数据的测定条件和化学空间覆盖存在局限，单纯依赖数据驱动模型可能难以可靠外推至临床阶段候选分子，因此本研究采用结构基础方法作为 URAT1 侧的主要筛选依据。相比之下，NLRP3 公开活性数据规模较大，可用于建立活性优先级模型，但其预测结果仅作为候选缩库证据，而非独立的活性确证。对接排名对诱饵构成和基准设计高度敏感（Gu 等，2025），URAT1 侧协议因而不能仅依据公开基准或单一打分函数预先指定。近年 NLRP3 NACHT 结构解析同样为基于结构的配体识别提供了实验基础（Dekker 等，2021）。

为回答这一问题，本文分别为两个靶点建立独立的证据筛选框架，并在统一候选池上整合 URAT1 结构对接证据与 NLRP3 相关活性线索，提出可供实验检验的 URAT1–NLRP3 双节点候选分子。文中提名仍需通过体外生化、细胞及进一步药理实验验证。

---

## 引用清单（定稿时改为期刊格式）

1. Dalbeth N, et al. *Lancet* 2021. doi:10.1016/S0140-6736(21)00569-9  
2. Leask MP, et al. *Nat Rev Rheumatol* 2024. doi:10.1038/s41584-024-01137-1  
3. Martinon F, et al. *Nature* 2006. doi:10.1038/nature04516  
4. Dai Y, Lee CH. *Cell Res* 2024. doi:10.1038/s41422-024-01023-1  
5. Du L, Zong Y, Li H, et al. Hyperuricemia and its related diseases: mechanisms and advances in therapy. *Signal Transduct Target Ther.* 2024;9:212. doi:10.1038/s41392-024-01916-y  
6. Liu Y, Li W, Deng Y. Role of NLRP3 in the pathogenesis and treatment of gout arthritis. *Front Immunol.* 2023;14:1137822. doi:10.3389/fimmu.2023.1137822  
7. FitzGerald JD, Dalbeth N, Mikuls T, et al. 2020 American College of Rheumatology guideline for the management of gout. *Arthritis Rheumatol.* 2020;72(6):879-895. doi:10.1002/art.41247  
8. Suo Y, Fedor JG, Zhang H, et al. *Nat Commun.* 2025;16:5178. doi:10.1038/s41467-025-60480-3 （Suo 与 Fedor 同等贡献；期刊作 Suo et al.）  
9. Wu C, Zhang C, Jin S, et al. *Cell Discov.* 2025;11:33. doi:10.1038/s41421-025-00779-z  
10. Hung SI, Chung WH, Liou LB, et al. HLA-B*5801 allele as a genetic marker for severe cutaneous adverse reactions caused by allopurinol. *Proc Natl Acad Sci USA.* 2005;102(11):4134-4139. doi:10.1073/pnas.0409500102  
11. White WB, Saag KG, Becker MA, et al. Cardiovascular safety of febuxostat or allopurinol in patients with gout (CARES). *N Engl J Med.* 2018;378(13):1200-1210. doi:10.1056/NEJMoa1710895  
12. Mackenzie IS, Ford I, Nuki G, et al. Long-term cardiovascular safety of febuxostat compared with allopurinol in patients with gout (FAST). *Lancet.* 2020;396(10264):1745-1757. doi:10.1016/S0140-6736(20)32234-0  
13. Borghi C, Domienik-Karłowicz J, Tykarski A, et al. Expert consensus for the diagnosis and treatment of patients with hyperuricemia and high cardiovascular risk: 2023 update. *Cardiol J.* 2024;31(1):1-14. doi:10.5603/cj.98254  
14. Lee MH, Graham GG, Williams KM, Day RO. A benefit-risk assessment of benzbromarone in the treatment of gout. *Drug Saf.* 2008;31(8):643-665. doi:10.2165/00002018-200831080-00002  
15. Tausche AK, Alten R, Dalbeth N, et al. Lesinurad monotherapy in gout patients intolerant to a xanthine oxidase inhibitor. *Rheumatology (Oxford).* 2017;56(12):2170-2178. doi:10.1093/rheumatology/kex305  
16. Coll RC, Robertson AAB, Chae JJ, et al. A small-molecule inhibitor of the NLRP3 inflammasome for the treatment of inflammatory diseases. *Nat Med.* 2015;21(3):248-255. doi:10.1038/nm.3806  
17. Tang F, Kunder R, Chu T, et al. First-in-human phase 1 trial of NLRP3 inhibitor GDC-2394. *Clin Transl Sci.* 2023;16(9):1653-1666. doi:10.1111/cts.13576 （文中明确：MCC950/CP-456,773 的 II 期因肝酶升高终止）  
18. Pushpakom S, Iorio F, Eyers PA, et al. Drug repurposing: progress, challenges and recommendations. *Nat Rev Drug Discov.* 2019;18(1):41-58. doi:10.1038/nrd.2018.168  
19. Zhang Z, Shi X, Wu T, et al. Discovery of multi-target anti-gout agents from *Eurycoma longifolia* Jack. *Nat Commun.* 2025;16:7430. doi:10.1038/s41467-025-62645-6  
20. Gu S, Zhang X, Shen C, et al. Benchmarking AI-powered docking methods from the perspective of virtual screening. *Nat Mach Intell.* 2025. doi:10.1038/s42256-025-00993-0  
21. Dekker C, Mattes H, Wright M, et al. Crystal structure of NLRP3 NACHT domain with an inhibitor defines mechanism of inflammasome inhibition. *J Mol Biol.* 2021;433(24):167309. doi:10.1016/j.jmb.2021.167309  
