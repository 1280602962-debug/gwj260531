# 引言初稿（中文）

> 投稿正文入口。目标期刊：*Molecular Diversity*（拒稿后可转 *JCAMD*）。  
> Methods：[`METHODS_DRAFT_CN.md`](METHODS_DRAFT_CN.md)。大纲：[`MANUSCRIPT.md`](MANUSCRIPT.md)。  
> 正文引用用“作者，年份”，DOI 见文末清单。

痛风是成人最常见的晶体性关节炎之一。当血尿酸长期升高，单钠尿酸盐可在关节及其周围沉积；巨噬细胞识别晶体后激活 NLRP3 炎症小体，经 caspase-1 将前体白细胞介素-1β 加工为活性形式，从而触发急性炎症发作（Dalbeth 等，2021；Leask 等，2024；Martinon 等，2006）。尿酸主要来自肝脏嘌呤代谢；多数患者的高尿酸血症并非生成过多，而是肾脏或肠道排泄不足。肾脏尿酸重吸收由多个转运体共同完成，其中近端小管转运体 URAT1（由 *SLC22A12* 编码）承担约九成尿酸重吸收，是排泄侧的关键节点之一，GLUT9/SLC2A9 与 ABCG2 等转运体也参与尿酸处理（Dai 与 Lee，2024；Suo 等，2025）。URAT1 与 NLRP3 之间并非直接的蛋白信号通路，而是经尿酸稳态与晶体形成串联起来的两个节点：URAT1 所在的上游尿酸稳态节点，与 NLRP3 所在的下游晶体驱动炎症节点（Liu 等，2023）。

现有药物也分别作用于这两个节点（FitzGerald 等，2020）。降尿酸一端，黄嘌呤氧化酶抑制剂减少尿酸生成；促排泄一端则用 URAT1 抑制剂，包括苯溴马隆、lesinurad、dotinurad，以及仍在临床试验中的 verinurad（Dai 与 Lee，2024；Wu 等，2025）。近年 URAT1 冷冻电镜结构显示，这些抑制剂结合于其尿酸盐结合位点并稳定其抑制性构象，为该靶点后续的结构对接提供了实验基础（Wu 等，2025；Suo 等，2025）。急性发作的一线药物是秋水仙碱、非甾体抗炎药或糖皮质激素；仅当这些药物无效、不耐受或禁忌时，才考虑白细胞介素-1 抑制剂（FitzGerald 等，2020）。启动降尿酸治疗时，指南还建议同时给予抗炎预防、疗程至少 3–6 个月，以减少诱发发作（FitzGerald 等，2020）——这一临床操作本身提示，降尿酸与抗炎需要分别给药、分别管理，二者控制的是病理链条上的不同阶段。

两类治疗的作用范围并不完全重叠。降尿酸治疗能够降低尿酸负荷，但对已经沉积的 MSU 晶体所诱导的急性炎症没有直接、快速的消除作用；抗炎治疗能够控制发作期炎症，却不降低尿酸负荷，也不清除已沉积的晶体（Dalbeth 等，2021）。二者各自还受安全性约束：别嘌醇在 HLA-B\*5801 阳性个体中与严重皮肤不良反应密切相关（Hung 等，2005）；非布司他在已有心血管疾病的痛风患者中，CARES 试验报告心血管死亡高于别嘌醇（White 等，2018），其后 FAST 试验未再现这一信号（Mackenzie 等，2020），但对高心血管风险人群，共识仍主张谨慎使用（Borghi 等，2024）；苯溴马隆因肝毒性在多地限制使用（Lee 等，2008）；lesinurad 单药治疗时肾相关不良事件明显增加，不能单独给药（Tausche 等，2017）。NLRP3 抑制已有明确的抗炎药理基础，MCC950 即是这类选择性抑制剂的代表（Coll 等，2015），但其小分子的临床转化仍面临安全性、选择性和疗效验证方面的挑战：同属该化学系列的 GDC-2394 在健康志愿者首次人体 I 期试验中，因两名受试者出现严重肝损伤而终止（Tang 等，2023），目前没有 NLRP3 小分子抑制剂进入痛风常规治疗。痛风还常与心血管病、代谢病和慢性肾脏病并存（Dalbeth 等，2021；Du 等，2024）。

上述分工提示，现行治疗需要通过不同药物分别控制尿酸负荷与晶体炎症；这提示同时调控上游尿酸稳态节点与下游晶体炎症节点的分子，理论上可能覆盖病理链条上更完整的区段。但这并不等于双节点分子必然优于现有的降尿酸联合抗炎方案，本文也不检验这一优效性问题。真正尚不清楚的是：在已经进入人体研究的小分子中，是否存在同时具备 URAT1 结构结合证据与 NLRP3 相关活性线索的候选分子——这是本文关注的具体转化问题。

从头设计新化学实体来回答这一问题，周期长、淘汰率高。药物重定位转而在已经进入人体研究的分子中寻找新用途：这类分子可在一定程度上利用已有的药代、安全性或临床开发信息，降低部分早期成药性风险，为后续实验验证提供更明确的转化背景（Pushpakom 等，2019）。本文据此把筛选对象限定为 ChEMBL 中已进入临床阶段或已经上市的小分子（具体纳入口径见 Methods）。与面向新化学空间的大规模虚拟筛选不同，本研究关注一个更受约束的问题：在已具有临床药理和人体暴露证据的这批分子中，能否找到可进一步验证的 URAT1–NLRP3 双节点候选。天然产物衍生路线已有湿法先例：从长叶阔蕊苏出发的工作表明，经结构优化后的同一分子可以同时影响尿酸盐转运和 NLRP3 介导的炎症（Zhang 等，2025）；那项研究回答的是新实体能否做成双作用分子，本文回答的则是已进入人体研究的分子中是否存在同样的双节点线索，两者互补而非重复。

要在这批临床阶段分子上核验上述线索，URAT1 与 NLRP3 两侧的公开证据并不对等，筛选方法也不能直接照搬。URAT1 侧，ChEMBL 活性来自不同测定，清洗后训练覆盖变窄，若干已上市尿酸药的骨架也不在其中，回归模型难以外推到这批分子，因而该侧排序应交给结构对接，前述 URAT1 抑制剂的结构证据为此提供了基础。但对接排名对诱饵构成和基准设计高度敏感（Gu 等，2025），生产协议因此不能仅依据公开基准或单一打分函数预先指定，而需要在独立构建的活性物–诱饵基准上完成评价与选择（具体构建见 Methods）。NLRP3 侧，公开测定条件同样不均一，但化合物数量足以训练分类器；该分类分数用于优先级排序与候选缩库，而非对 NLRP3 抑制活性的独立确证。近年 NLRP3 NACHT 结构解析同样为基于结构的配体识别提供了实验基础（Dekker 等，2021）。

因此，本研究聚焦于一个具体的转化问题：在已经进入人体研究或已经上市的小分子中，是否存在能够同时提供 URAT1 结构对接证据与 NLRP3 相关活性线索的候选分子。为降低单一预测模型或单一对接评分带来的偏差，本文分别为两个靶点建立独立的证据筛选框架，并在统一候选池上整合结构与活性证据，提出可供实验检验的 URAT1–NLRP3 双节点假说。对接分数只在同一池内比较高低，不解释为亲和力；文中提名也不等于已经证实的双靶抑制剂，更不构成临床用药建议，双节点分子是否优于现有联合用药方案也不在本文验证范围内。

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
