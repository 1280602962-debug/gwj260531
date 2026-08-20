# 引言初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)。Results：[`RESULTS_DRAFT_CN.md`](RESULTS_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用“作者，年份”，DOI 见文末清单。

痛风是成人最常见的晶体性关节炎之一，其发生发展与持续性高尿酸血症及单钠尿酸（monosodium urate，MSU）晶体沉积密切相关。多数患者发病并非尿酸生成过多，主要是与尿酸排泄不足有关。当血尿酸水平长期升高并超过尿酸盐的溶解度阈值时，MSU晶体可在关节及其周围组织中形成和沉积（Dalbeth 等，2021）。MSU晶体形成后，可被巨噬细胞等先天免疫细胞识别并诱导NLRP3炎症小体活化，并进一步启动炎症级联反应，最终导致急性痛风性关节炎的发生（Martinon 等，2006；Dalbeth 等，2021；Leask 等，2024）。因此，痛风的病理过程并非单一环节所致，而是呈现出由尿酸稳态失衡和MSU晶体形成向晶体驱动的炎症反应连续发展的病理链条。其中，介导肾近端小管尿酸重吸收的关键转运体URAT1位于上游尿酸稳态调控环节，影响高尿酸血症及后续晶体形成风险（Dai 与 Lee，2024；Suo 等，2025）；NLRP3则位于下游晶体诱导炎症的关键环节，参与介导MSU晶体所触发的炎症反应（Martinon 等，2006）。基于这一连续的病理轴，降低尿酸水平以减少MSU晶体形成与沉积，以及抑制NLRP3介导的炎症反应，分别对应痛风治疗中的两个关键干预节点（Liu 等，2023）。

现有痛风药物治疗主要包括降尿酸治疗和控制急性炎症发作两类策略（FitzGerald 等，2020）。降尿酸一端，黄嘌呤氧化酶抑制剂减少尿酸生成；促排泄一端则用 URAT1 抑制剂，包括苯溴马隆、lesinurad、dotinurad，以及仍在临床试验中的 verinurad（Dai 与 Lee，2024；Wu 等，2025）。近年 URAT1 冷冻电镜结构显示，多种抑制剂结合于其中央配体结合腔体并稳定 inward-open 构象，为该靶点的结构基础筛选提供了实验依据（Wu 等，2025；Suo 等，2025）。急性发作的一线药物是秋水仙碱、非甾体抗炎药或糖皮质激素；仅当这些药物无效、不耐受或禁忌时，才考虑白细胞介素-1 抑制剂（FitzGerald 等，2020）。启动降尿酸治疗时，指南还建议同时给予抗炎预防、疗程至少 3–6 个月，以减少诱发发作（FitzGerald 等，2020）。因此，现行痛风治疗需要通过不同药物分别控制尿酸负荷和晶体驱动的炎症反应，提示同时作用于病理链条上游尿酸负荷和下游炎症环节的分子可能具有潜在的药理学价值。近期已有研究通过天然产物衍生和结构优化获得同时影响尿酸盐转运与 NLRP3 炎症的双作用分子，为这一双节点药理策略提供了实验可行性的先例（Zhang 等，2025）。

然而，这两种治疗方式的作用范围并不完全重叠。降尿酸治疗能够降低尿酸负荷，但对已经沉积的 MSU 晶体所诱导的急性炎症没有直接、快速的消除作用；抗炎治疗能够控制发作期炎症，却不降低尿酸负荷，也不清除已沉积的晶体（Dalbeth 等，2021）。二者各自还受到药物安全性及适用范围的限制。别嘌醇存在 HLA-B\*5801 相关严重皮肤不良反应风险（Hung 等，2005）；非布司他在高心血管风险患者中的安全性仍需谨慎评估（White 等，2018；Mackenzie 等，2020；Borghi 等，2024）；苯溴马隆和 lesinurad 也分别受到肝脏或肾脏安全性问题限制（Lee 等，2008；Tausche 等，2017）。另一方面，尽管 NLRP3 抑制具有明确的抗炎药理基础（Coll 等，2015），相关小分子的临床转化仍面临安全性和疗效验证等挑战。例如，NLRP3 小分子抑制剂 GDC-2394 在首次人体研究中曾因严重肝损伤事件而终止开发（Tang 等，2023），目前尚无 NLRP3 小分子抑制剂成为痛风的常规治疗药物。痛风还常与心血管病、代谢病和慢性肾脏病并存（Dalbeth 等，2021；Du 等，2024）。

近年来，全新的小分子药物从发现到临床可用周期长、淘汰率高（Pushpakom 等，2019）。痛风双节点因此可以先在已经进入人体研究的小分子中，寻找同时具备 URAT1 结合与 NLRP3 相关活性的计算线索。这类分子通常已有一定的人体药理、药代或安全性信息，可在一定程度上降低早期成药性和转化风险，并为后续实验验证提供更明确的开发背景（Pushpakom 等，2019）。综上所述，本研究试图从 ChEMBL 数据库中临床阶段及已上市小分子里筛选出可以分别作用于痛风发病机制轴上下两处的双作用候选分子。

URAT1与NLRP3两侧的公开活性证据在数据规模、测定条件及化学空间覆盖方面并不对等，因此两侧不能简单采用相同的数据驱动筛选策略。URAT1侧公开活性数据的化学空间覆盖有限，单纯依赖数据驱动模型可能难以可靠外推至临床阶段候选分子，因此本研究以结构基础方法作为主要筛选依据；NLRP3侧则利用规模较大的公开活性数据建立优先级模型，其结果仅用于候选缩库，而不作为独立的活性确证。

因此，本研究聚焦于一个具体的转化问题：在已经进入人体研究或临床使用的小分子中，是否存在能够同时提供URAT1结构结合和NLRP3相关活性证据的候选分子？为降低单一预测模型或单一对接评分带来的偏差，本研究分别针对两个靶点建立独立的证据筛选框架，并在统一的候选池中进行结构和活性证据整合，从而提出可供后续实验验证的URAT1–NLRP3双节点候选分子。

---

## 引用清单（定稿时改为期刊格式）

1. Martinon F, et al. Gout-associated uric acid crystals activate the NALP3 inflammasome. *Nature* 2006. doi:10.1038/nature04516  
2. Dalbeth N, et al. Gout. *Lancet* 2021. doi:10.1016/S0140-6736(21)00569-9  
3. Leask MP, et al. The pathogenesis of gout: molecular insights from genetic, epigenomic and transcriptomic studies. *Nat Rev Rheumatol* 2024. doi:10.1038/s41584-024-01137-1  
4. Liu Y, Li W, Deng Y. Role of NLRP3 in the pathogenesis and treatment of gout arthritis. *Front Immunol.* 2023;14:1137822. doi:10.3389/fimmu.2023.1137822  
5. FitzGerald JD, Dalbeth N, Mikuls T, et al. 2020 American College of Rheumatology guideline for the management of gout. *Arthritis Rheumatol.* 2020;72(6):879-895. doi:10.1002/art.41247  
6. Dai Y, Lee CH. Transport mechanism and structural pharmacology of human urate transporter URAT1. *Cell Res* 2024. doi:10.1038/s41422-024-01023-1  
7. Wu C, Zhang C, Jin S, et al. *Cell Discov.* 2025;11:33. doi:10.1038/s41421-025-00779-z  
8. Suo Y, Fedor JG, Zhang H, et al. *Nat Commun.* 2025;16:5178. doi:10.1038/s41467-025-60480-3 （Suo 与 Fedor 同等贡献；期刊作 Suo et al.）  
9. Zhang Z, Shi X, Wu T, et al. Discovery of multi-target anti-gout agents from *Eurycoma longifolia* Jack. *Nat Commun.* 2025;16:7430. doi:10.1038/s41467-025-62645-6  
10. Hung SI, Chung WH, Liou LB, et al. HLA-B\*5801 allele as a genetic marker for severe cutaneous adverse reactions caused by allopurinol. *Proc Natl Acad Sci USA.* 2005;102(11):4134-4139. doi:10.1073/pnas.0409500102  
11. White WB, Saag KG, Becker MA, et al. Cardiovascular safety of febuxostat or allopurinol in patients with gout (CARES). *N Engl J Med.* 2018;378(13):1200-1210. doi:10.1056/NEJMoa1710895  
12. Mackenzie IS, Ford I, Nuki G, et al. Long-term cardiovascular safety of febuxostat compared with allopurinol in patients with gout (FAST). *Lancet.* 2020;396(10264):1745-1757. doi:10.1016/S0140-6736(20)32234-0  
13. Borghi C, Domienik-Karłowicz J, Tykarski A, et al. Expert consensus for the diagnosis and treatment of patients with hyperuricemia and high cardiovascular risk: 2023 update. *Cardiol J.* 2024;31(1):1-14. doi:10.5603/cj.98254  
14. Lee MH, Graham GG, Williams KM, Day RO. A benefit-risk assessment of benzbromarone in the treatment of gout. *Drug Saf.* 2008;31(8):643-665. doi:10.2165/00002018-200831080-00002  
15. Tausche AK, Alten R, Dalbeth N, et al. Lesinurad monotherapy in gout patients intolerant to a xanthine oxidase inhibitor. *Rheumatology (Oxford).* 2017;56(12):2170-2178. doi:10.1093/rheumatology/kex305  
16. Coll RC, Robertson AAB, Chae JJ, et al. A small-molecule inhibitor of the NLRP3 inflammasome for the treatment of inflammatory diseases. *Nat Med.* 2015;21(3):248-255. doi:10.1038/nm.3806  
17. Tang F, Kunder R, Chu T, et al. First-in-human phase 1 trial of NLRP3 inhibitor GDC-2394. *Clin Transl Sci.* 2023;16(9):1653-1666. doi:10.1111/cts.13576  
18. Du L, Zong Y, Li H, et al. Hyperuricemia and its related diseases: mechanisms and advances in therapy. *Signal Transduct Target Ther.* 2024;9:212. doi:10.1038/s41392-024-01916-y  
19. Pushpakom S, Iorio F, Eyers PA, et al. Drug repurposing: progress, challenges and recommendations. *Nat Rev Drug Discov.* 2019;18(1):41-58. doi:10.1038/nrd.2018.168  
