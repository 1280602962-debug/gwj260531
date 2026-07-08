# NLRP3抑制剂开发领域最新进展综述

**课程作业**

---

## 摘要

NLRP3炎症小体在多种炎症性疾病中发挥重要作用，是近年来药物研发关注较多的靶点之一。本文结合2023年至2026年发表的相关文献，对NLRP3抑制剂的分子机制研究、药物化学进展和临床开发现状进行梳理和讨论。全文共分四个部分，从NLRP3的激活机制与结构研究入手，进而讨论代表性抑制剂的化学特征和临床管线布局，最后分析该领域面临的主要问题和发展趋势。目前尚无直接靶向NLRP3的药物获批上市，但dapansutrile、NT-0796、VTX2735、BGE-102、usnoflast和selnoflast等多个候选化合物已进入临床II期试验。总体来看，该领域在结构生物学和临床管线两方面均有实质性推进，但临床转化的直接证据仍然有限，正处于从实验室走向临床验证的过渡阶段。

**关键词**：NLRP3炎症小体，小分子抑制剂，dapansutrile，临床开发

---

## 1 引言

NLRP3（NOD样受体热蛋白结构域相关蛋白3）炎症小体是细胞内重要的免疫信号复合物，由传感器蛋白NLRP3、接头蛋白ASC（凋亡相关斑点样蛋白）和效应蛋白酶caspase-1组成。细胞受到病原体相关分子模式或损伤相关分子模式刺激后，NLRP3炎症小体被激活，促进IL-1β和IL-18等促炎细胞因子成熟并释放到细胞外，同时触发焦亡，进一步放大炎症反应。NLRP3过度激活或调节失衡与痛风、动脉粥样硬化、2型糖尿病、神经退行性疾病以及多种自身炎症性疾病密切相关，因此被认为是有潜力的治疗靶点。

在药物开发策略上，目前临床上已有canakinumab、anakinra等靶向IL-1通路的生物制剂，在类风湿关节炎、自身炎症性疾病等领域取得了一定疗效。但这些药物作用于NLRP3的下游，只能阻断IL-1β这一条通路，且多为注射给药，长期使用成本较高，患者依从性也面临一定挑战。直接抑制NLRP3本身的小分子药物，理论上可通过阻断炎症小体组装和caspase-1激活，从而减少IL-1β和IL-18的成熟与释放，口服给药也更为方便。Vande Walle和Lamkanfi在2024年发表于*Nature Reviews Drug Discovery*的综述《Drugging the NLRP3 inflammasome: from signalling mechanisms to therapeutic targets》[1]中系统总结了该领域的研究现状，指出多个小分子和生物制剂候选药物正在向临床推进，但截至该文发表时，尚无NLRP3直接抑制剂获得监管批准。

NLRP3是连接基础免疫学和成药开发的重要节点，围绕其展开综述有助于将机制研究和临床转化两条线串起来理解。本文主要依据2023年至2026年发表的结构生物学、药物化学和临床相关文献进行归纳和评述，所引用的参考文献均为正式发表的论文或综述，可在PubMed或期刊官网检索核实，文末附有DOI链接。需要说明的是，很多候选药物的临床数据仍以综述论文和临床试验注册信息为主，正式发表的大样本临床论文还相对有限，因此本文在讨论临床进展时会对此有所交代，并尽量将结论与可核实的文献来源对应。全文共分四个部分，依次介绍NLRP3的激活机制与结构研究、抑制剂的药物化学进展、临床开发现状，以及当前面临的主要问题与展望。

## 2 NLRP3激活机制与结构生物学研究进展

要理解NLRP3抑制剂为什么有效、应该靶向哪个位点，首先需要了解NLRP3是如何被激活的。根据Vande Walle和Lamkanfi（2024，*Nat Rev Drug Discov*）[1]的综述，在经典（canonical）激活途径中，NLRP3的激活通常需要两个信号。第一个信号称为priming，通常由Toll样受体等通路触发，使细胞表达更多的NLRP3蛋白并做好激活准备。第二个信号称为activation，由钾离子外流、溶酶体损伤、活性氧产生、线粒体功能障碍等多种刺激触发，导致NLRP3炎症小体组装。NLRP3被激活后会招募ASC并形成ASC斑点，进而激活caspase-1，最终促进IL-1β和IL-18成熟释放。NLRP3能够响应理化性质差异很大的激活剂，这一特点长期以来是理解其调控机制的主要难点，也是药物设计需要面对的复杂性所在。

近三年冷冻电镜结构的解析，使NLRP3激活过程的分子图景变得清晰了许多。Le等人在2023年发表于*Nature*的论文《Cryo-EM structures of the active NLRP3 inflammasome disc》[3]报道了活性NLRP3炎症小体"盘状"结构的冷冻电镜结果。研究发现，NLRP3在激活后NACHT结构域发生约90°的旋转，从ADP结合的关闭状态转变为ATP结合的开放状态。NEK7（NIMA相关激酶7）结合LRR结构域后，可促进非活性"笼状"寡聚体向具有组装能力的构象转变，进而促进活性盘状结构的形成。在这一过程中，NLRP3特有的FISNA结构域在稳定活性构象和介导盘状寡聚体相互作用方面发挥了关键作用，这是NLRP3区别于多数其他NLR家族成员的结构性特征（图1）。

![图1 活性NLRP3炎症小体盘状结构（Le et al., 2023, Nature）[3]](./NLRP3_review_figures/fig1_nlrp3_disc_nature.png)

*图1 活性NLRP3炎症小体盘状冷冻电镜结构示意图。NLRP3 NACHT结构域在激活过程中发生约85°–90°旋转，FISNA结构域在稳定活性构象中发挥关键作用。改绘自Le et al., *Nature* 2023[3]。*

2024年，Yu等人在*Nature Communications*发表的《Structural basis for the oligomerization-facilitated NLRP3 activation》[2]进一步解析了NLRP3开放八聚体的结构，填补了从关闭笼状结构到活性盘状结构之间的中间态空白。研究发现，NLRP3在转变过程中NACHT结构域同样发生约90°的铰链旋转，形成具有Head-Face和Tail-Tail相互作用界面的开放八聚体。NEK7的结合有助于促进非活性寡聚体向具备组装能力的物种转化，形成NEK7/NLRP3单体或二聚体，这是后续盘状炎症小体组装的前提。研究团队对开放八聚体界面进行突变后，IL-1β信号明显减弱，说明寡聚化协同激活是NLRP3炎症小体组装的核心机制（图2、图3）。

![图2 NLRP3开放八聚体冷冻电镜结构（Yu et al., 2024, Nat Commun）[2]](./NLRP3_review_figures/fig2_open_octamer.png)

*图2 NLRP3开放八聚体（NLRP3ΔPYD/+ATP）的冷冻电镜结构，显示NACHT域铰链旋转及Head-Face、Tail-Tail等寡聚界面。改绘自Yu et al., *Nat Commun* 2024[2]。*

![图3 NLRP3激活机制模型（Yu et al., 2024, Nat Commun）[2]](./NLRP3_review_figures/fig5_activation_model.png)

*图3 NLRP3从关闭笼状结构经开放八聚体向活性盘状结构转变的激活模型。改绘自Yu et al., *Nat Commun* 2024[2]。*

这些结构工作的重要性在于，它们不仅解释了"NLRP3如何被激活"，也为理解现有抑制剂的作用位点提供了直接的结构参照。目前已有多个小分子抑制剂的结合模式被解析，它们大多结合于NACHT结构域的中央疏水口袋，通过稳定关闭构象阻断向活性盘状结构的转变。NLRP3的激活并不是简单的"开"和"关"两种状态，而是经历多个中间步骤的复杂过程，这意味着药物研发者可以选择在不同构象状态或不同组装环节进行干预，而不必都集中在MCC950所结合的那个口袋上。

除经典的离子流变化和溶酶体损伤等激活途径外，近三年研究还发现氧化DNA（oxDNA），特别是氧化线粒体DNA（ox-mtDNA），可直接参与NLRP3的激活，这是机制研究上的一个重要补充。Cabral等人在2023年*Communications Biology*发表的《Differential binding of NLRP3 to non-oxidized and Ox-mtDNA mediates NLRP3 inflammasome activation》[4]表明，NLRP3对氧化与非氧化mtDNA的结合能力存在明显差异，分离的PYD结构域更倾向于结合ox-mtDNA，不过其精确的分子机制目前仍在深入研究之中（图4）。NLRP3的PYD结构域在蛋白质折叠方式上与DNA糖基化酶hOGG1存在结构同源性，这为进一步理解NLRP3如何"感知"氧化应激提供了分子基础。在病理条件下，线粒体功能障碍和活性氧升高可导致ox-mtDNA释放到细胞质，进而触发NLRP3炎症小体组装，这一机制将代谢应激、氧化损伤和先天免疫激活串联起来。

![图4 NLRP3与ox-mtDNA结合机制（Cabral et al., 2023, Commun Biol）[4]](./NLRP3_review_figures/fig3_oxmtdna.png)

*图4 NLRP3 PYD结构域对氧化与非氧化线粒体DNA（ox-mtDNA）的差异化结合。改绘自Cabral et al., *Commun Biol* 2023[4]。*

Lackner等人在2025年*Trends in Biochemical Sciences*发表的综述《How interactions between oxidized DNA and the NLRP3 inflammasome fuel inflammatory disease》[5]系统讨论了oxDNA与NLRP3的相互作用，认为ox-mtDNA是重要的内源性危险信号，可将氧化应激与先天免疫激活联系起来，并在多种慢性炎症和自身炎症性疾病中发挥病理作用。该综述还提到，hOGG1抑制剂TH5487和SU0268可阻断NLRP3与ox-mtDNA的结合，从而抑制炎症小体激活和IL-1β释放。这一发现提示，除了传统的NACHT结构域结合位点之外，PYD结构域和oxDNA通路也可能成为新的药物靶点。不过需要指出的是，目前基于这一机制开发的临床候选药物还很少，大多数进入临床的化合物仍然结合于NACHT结构域的CRID3类似口袋。从药物设计角度看，oxDNA通路的发现说明NLRP3抑制剂的研发不必局限于"封闭NACHT、抑制ATP酶"这一种思路，未来可能有更多差异化作用机制的化合物出现。

## 3 NLRP3抑制剂的药物化学与临床候选化合物

MCC950（又称CRID3或CP-456773）是目前研究中最常用的NLRP3抑制剂工具化合物，也是许多后续药物的结构起点和比较参照。它属于磺酰脲类化合物，结合于NACHT结构域中央疏水口袋，靠近Walker A motif，通过稳定NLRP3的关闭构象并抑制ATP酶活性发挥作用（图5）。Coll等人在2019年发表于*Nature Chemical Biology*的论文《MCC950 directly targets the NLRP3 ATP-hydrolysis motif for inflammasome inhibition》[11]证实，MCC950直接与NLRP3 NACHT域Walker B motif相互作用，阻断ATP水解。Hochheiser等人在2022年*Nature*发表的《Structure of the NLRP3 decamer bound to the cytokine release inhibitor CRID3》[14]进一步以冷冻电镜解析了CRID3/MCC950结合于NACHT与LRR交界处的分子细节，为其磺酰脲基团与Walker A motif的特异性相互作用提供了结构依据（图5）。MCC950在多种动物炎症模型中显示出较好的抗炎效果，是NLRP3临床前研究中最为常见的阳性对照化合物。这也带来一个问题：许多结论是否能在非MCC950类化合物上重复，仍需要更多验证。MCC950最初由Pfizer发现，但后续临床开发因早期临床评价中观察到肝毒性而中止。此后Genentech开发的GDC-2394（JT001）也因在健康志愿者中出现药物性肝损伤而停止临床开发。Kennedy等人在2021年*SLAS Discovery*的研究[15]揭示，MCC950对碳酸酐酶II（CA2）具有非竞争性抑制活性，这一脱靶效应可能与不良反应有关。上述案例说明临床前模型中的高活性并不能保证临床安全性，后续化合物必须在化学结构上做出实质性改进，而不能简单地在MCC950骨架上做微小修饰。Vande Walle和Lamkanfi（2024）[1]也指出，MCC950的肝毒性可能与化合物结构本身以及脱靶效应有关，这推动了业界对新型化学骨架的迫切需求。

![图5 MCC950与NLRP3 NACHT域结合（PDB 7ALW）](./NLRP3_review_figures/fig_mcc950.png)

*图5 MCC950结合于NLRP3 NACHT域中央疏水口袋的冷冻电镜/晶体结构示意（PDB: 7ALW）。*

MCC950的肝毒性问题推动了新型化学骨架的探索。Vande Walle等人在2024年*Life Science Alliance*发表的《Novel chemotype NLRP3 inhibitors that target the CRID3-binding pocket with high potency》[6]通过PubChem高通量筛选，发现了一类吡咯并三嗪乙酰胺类分子（NIC系列）。需要强调的是，NIC系列并非从MCC950结构优化而来，而是通过独立的高通量筛选获得的全新骨架，但其与NLRP3的物理结合经BRET实验证实，结构模拟显示其结合于CRID3/MCC950的同一口袋，却采取明显不同的结合构象（图6）。值得注意的是，NIC-12缺乏MCC950对碳酸酐酶I和II的脱靶抑制活性，在LPS内毒素血症小鼠模型中可选择性降低循环IL-1β水平，在冷吡啉相关周期性综合征（cryopyrin-associated periodic syndromes，CAPS）患者来源的单核细胞中效力约为CRID3的10倍。该文的通讯作者Lamkanfi长期深耕NLRP3炎症小体研究，这项工作是在MCC950之后寻找替代化学骨架的延续。

![图6 NIC系列化合物筛选与活性（Vande Walle et al., 2024, Life Sci Alliance）[6]](./NLRP3_review_figures/fig4_nic_series.jpg)

*图6 NIC系列吡咯并三嗪乙酰胺类NLRP3抑制剂的筛选流程与代表性结构。改绘自Vande Walle et al., *Life Sci Alliance* 2024[6]。*

Matico等人在2025年*EMBO Molecular Medicine*发表的《Navigating from cellular phenotypic screen to clinical candidate: selective targeting of the NLRP3 inflammasome》[7]则从表型高通量筛选出发，发现了一类化学结构不同的吡咯并三嗪乙酰胺类化合物，经构效关系优化后得到化合物C。冷冻电镜、nanoDSF和氢氘交换质谱均证实该类化合物结合于人源NLRP3的NACHT结构域，且对NLRP1b和NLRC4炎症小体无明显交叉抑制。在CAPS相关动物模型中，化合物C显示出较好的体内疗效，作者提出可优先在CAPS患者中开展概念验证。NIC系列和化合物C之所以值得关注，是因为它们都通过冷冻电镜等实验手段明确了与NLRP3的结合模式，而不只是停留在细胞水平的活性筛选。

除上述结合于CRID3口袋的化合物外，药物化学正在向全新结合位点和作用机制拓展。Hartman等人在2024年*Bioorganic & Medicinal Chemistry Letters*发表的《The discovery of novel and potent indazole NLRP3 inhibitors enabled by DNA-encoded library screening》[12]通过DNA编码文库筛选发现了吲唑类先导化合物BAL-0028；Wilhelmsen等人在2024年bioRxiv预印本《Discovery of a Potent and Selective Inhibitor of Human NLRP3 with a Novel Binding Modality and Mechanism of Action》[13]进一步证实，BAL-0028结合于与MCC950不同的NACHT位点，不抑制ATP酶活性，且对部分CAPS致病突变（如D303H、L353P）的效力优于MCC950。BioAge Labs在此基础上开发的临床候选药BGE-102保留了这一新型结合模式，具有脑穿透性，目前已完成I期试验并进入II期开发。Cabral等人在2025年*Trends in Pharmacological Sciences*的综述《Targeting the NLRP3 inflammasome for inflammatory disease therapy》[9]提到，多个药企正在开发结合位点或化学骨架各不相同的NLRP3抑制剂，说明该领域的药物化学正在从"单一磺酰脲骨架"向多元化方向发展。

在诸多候选化合物中，dapansutrile（OLT1177）是目前临床进展较快、公开资料也相对较多的口服NLRP3抑制剂之一。该药由Olatec Therapeutics开发，为β-磺酰基腈类小分子，分子量仅133 Da，是已报道的分子量最小的NLRP3抑制剂之一，化学结构上与MCC950有明显区别，并非简单的磺酰脲衍生物。磺酰基腈类化合物通常具有较小的分子体积和较高的亲脂性，有利于跨膜转运；dapansutrile的低分子量进一步降低了血脑屏障穿透的分子量门槛，使其在神经炎症相关适应症中具备结构上的先天优势。Amo-Aparicio等人在2023年*Journal of Neuroinflammation*发表的《Pharmacologic inhibition of NLRP3 reduces the levels of α-synuclein and protects dopaminergic neurons in a model of Parkinson's disease》[8]表明，dapansutrile在MPTP帕金森病小鼠模型中能够穿过血脑屏障，在脑内达到有效分布，减少脑内α-突触核蛋白水平，并对多巴胺能神经元起到一定保护作用，改善了小鼠的运动功能。该研究为dapansutrile在神经炎症相关疾病中的应用提供了临床前依据，也将"分子量小"与"能入脑"之间的因果关系在实验层面予以印证。从临床开发进度来看，根据Cabral等人2025年的综述[9]，dapansutrile已在痛风、骨关节炎、心力衰竭等多种适应症中开展临床试验，并启动了针对2型糖尿病（NCT06047262）和急性痛风发作（NCT05658575）的II/III期研究。痛风被认为是NLRP3过度激活的代表性疾病，尿酸盐结晶激活NLRP3炎症小体是痛风急性发作的经典机制，因此痛风试验的结果对于判断NLRP3抑制剂能否在人体中起效具有标杆意义。

除dapansutrile之外，目前进入临床II期阶段的候选药物还包括NT-0796、VTX2735、VTX3232、BGE-102、usnoflast和selnoflast等，分属不同公司，化学结构和脑穿透性也各不相同（图7）。NT-0796是NodThera公司开发的脑穿透性口服NLRP3抑制剂，其异丙酯前药在体内转化为活性代谢物NDT-19795。Cabral等人2025年的综述[9]引用了该药Ib/IIa期试验（NCT06129409）的结果，在肥胖合并心血管风险因素且基线hsCRP升高的受试者中，用药组28天后超过75%的患者hsCRP降至2 mg/L以下，而安慰剂组这一比例不足25%。需要指出的是，hsCRP属于替代终点（surrogate endpoint），而非心血管事件等临床硬终点（clinical endpoint），上述结果说明NT-0796在人体中能够降低系统性炎症水平，但其能否最终转化为临床获益仍有待更大规模试验验证。VTX2735是Ventyx公司开发的外周限制性口服NLRP3抑制剂，已在CAPS患者中完成IIa期概念验证试验（NCT05812781），目前正在复发性心包炎患者中开展开放标签II期试验（NCT06836232）。VTX3232是Ventyx公司另一款脑穿透性NLRP3抑制剂，Bultinck等人在2026年*Molecular Metabolism*发表的《NLRP3 inhibition by VTX3232 tempers inflammation resulting in reduced body weight, hyperglycemia, and hepatic steatosis in obese male mice》[10]表明，VTX3232在肥胖小鼠模型中可减轻炎症反应，并改善体重、高血糖和肝脏脂肪变性。Selnoflast（RO7486967）是Roche从Inflazome收购后推进的外周限制性NLRP3抑制剂，在稳定冠心病、有心肌梗死史且hsCRP≥2 mg/L的患者中，28天口服给药后hsCRP和IL-6等炎症替代终点较安慰剂有明显下降，目前该药正在动脉粥样硬化患者中开展IIa期RIVULET试验（NCT07448038）。Usnoflast（ZYIL1）是Zydus公司自主研发的脑穿透性NLRP3抑制剂，2024年完成了ALS患者的IIa期概念验证试验，2025年已启动IIb期UNITE-ALS试验（NCT07023835），并获FDA快速通道资格。BGE-102是BioAge Labs开发的脑穿透性口服NLRP3抑制剂，2026年完成I期试验并启动心血管风险II期试验QUELL-CV（NCT07656727），其结合位点与MCC950不同，代表了一类非ATP酶抑制型候选药物。

![图7 NLRP3抑制剂临床开发现状（Vande Walle & Lamkanfi, 2024, Nat Rev Drug Discov）[1]](./NLRP3_review_figures/fig7_clinical_landscape.png)

*图7 主要NLRP3抑制剂候选药物的临床开发阶段与适应症布局。改绘自Vande Walle & Lamkanfi, *Nat Rev Drug Discov* 2024[1]；BGE-102等2025–2026年新增管线据公司公开信息更新。*

从上述管线的布局可以看出，当前临床开发呈现出几个明显特点。一是适应症选择上，企业普遍优先推进机制明确、炎症标志物可测量的疾病，如CAPS、复发性心包炎、急性痛风、肥胖伴高hsCRP等，而不是一开始就挑战阿尔茨海默病这类终点复杂、试验周期长的慢病。自身炎症性疾病是NLRP3抑制剂临床验证的重要切入点，因为这类疾病往往与NLRP3功能获得性突变直接相关，发病机制较为明确，疗效评价也相对直观。以CAPS为例，患者体内NLRP3蛋白本身存在激活增强的突变，如果抑制剂有效，理论上应能在较短时间内观察到症状和炎症标志物的改善，Cosson等人（2024，*J Exp Med*）[16]对多种CAPS突变的分层分析也为这一策略提供了分子层面的支持。二是脑穿透性与外周限制性的区分越来越受重视。神经系统疾病优先选择NT-0796、VTX3232、usnoflast、BGE-102等脑穿透性化合物，心血管代谢疾病则更多布局VTX2735、selnoflast等外周限制性化合物，以避免不必要的中枢神经系统免疫抑制风险。三是联合用药策略开始受到关注，NT-0796与司美格鲁肽的联合试验反映了同时针对代谢异常和慢性炎症两条通路的研发思路。Bultinck等人（2026，*Mol Metab*）[10]在肥胖小鼠中的VTX3232实验也提示，NLRP3抑制可能与减重药物产生代谢层面的协同效应，但这一发现在人体中能否重复，仍有待联合试验数据验证。

帕金森病和ALS等神经退行性疾病也有多个候选药物布局，Amo-Aparicio等人2023年的动物实验[8]为dapansutrile在这一方向的应用提供了初步支持，但神经退行性疾病的临床试验周期较长，患者功能评估的标准化程度也不如心血管指标那样成熟。usnoflast在ALS中的IIa期试验虽然显示药物耐受性良好且能进入脑脊液，但神经丝轻链等生物标志物的变化尚未达到统计学显著，这说明神经疾病方向的临床验证并不容易，脑内达到有效浓度只是第一步，能否在疾病进展中显示出有临床意义的获益才是更关键的问题。在阿尔茨海默病方面，尸检组织和脑脊液研究均显示NLRP3炎症小体相关蛋白表达升高，Vande Walle和Lamkanfi的综述[1]讨论了NLRP3在神经退行性疾病中的病理作用，认为小胶质细胞中NLRP3的慢性激活可能加剧Aβ和tau病理相关的神经炎症。Mangan等人（2018）[17]亦曾将阿尔茨海默病列为NLRP3抑制剂的重要潜在适应症，但截至目前尚无NLRP3抑制剂在阿尔茨海默病患者中开展临床试验，这一方向仍有较大空白。对于以神经退行性疾病为主要研究方向的课题组而言，NLRP3抑制剂可能仍需要等待痛风、CAPS或心血管代谢等适应症首先完成临床验证，再考虑向阿尔茨海默病拓展。

## 4 问题与展望

### 4.1 尚未解决的开放性问题

尽管该领域进展迅速，仍有若干根本性问题尚无定论，这些问题将决定NLRP3抑制剂能否成为与TNF-α或IL-17抑制剂比肩的成熟药物类别。

第一，**NLRP3的生理功能是否允许长期全面抑制？** NLRP3在宿主免疫防御中具有一定作用，Vande Walle和Lamkanfi（2024）[1]提醒，完全阻断其活性可能增加感染风险。目前I期和IIa期试验的随访时间多为数周至数月，对于需要长期口服给药的慢性病，长期用药对感染易感性、肿瘤监视和免疫系统稳态的影响仍缺乏充分数据。"抗炎"与"免疫抑制"之间的界限如何把握，尚无临床共识。

第二，**临床前模型能否可靠预测人体安全性？** MCC950和GDC-2394在动物实验中表现良好，却在人体中出现肝毒性，说明现有临床前安全性评价体系存在盲区。小鼠与人NLRP3在序列和药理学响应上存在差异——以BAL-0028/BGE-102为例，其对灵长类NLRP3高效但对小鼠NLRP3几乎无效[13]，常规小鼠模型甚至无法用于评价这类化合物。基因敲除动物还可能存在发育代偿，使得表型与药理学抑制不完全等同。如何建立更具预测力的人源化评价体系，仍是转化医学层面的开放问题。

第三，**炎症替代终点能否转化为临床硬终点？** 多个II期试验以hsCRP、IL-6等炎症标志物为主要疗效指标，但这些属于替代终点而非心血管事件、疾病复发率等硬终点。NT-0796的Ib/IIa期数据[9]虽显示hsCRP显著下降，但其能否降低心肌梗死或卒中风险，需要结局试验来回答。Mangan等人在2018年*Nature Reviews Drug Discovery*的综述《Targeting the NLRP3 inflammasome in inflammatory diseases》[17]曾系统梳理NLRP3抑制剂临床开发路径，指出炎症生物标志物虽便于早期概念验证，却难以单独支撑长期慢病用药的获益—风险评价。在尚无获批NLRP3抑制剂可供参照的情况下，监管机构和药企对疗效终点的选择标准也仍在摸索之中。

第四，**CAPS致病突变是否会导致抑制剂耐药？** Cosson等人在2024年*Journal of Experimental Medicine*的研究[16]表明，多数CAPS相关NLRP3功能获得性突变仍对MCC950敏感，但部分位于抑制剂结合口袋附近的突变（如D303H、L353P）可导致显著耐药。这意味着"同一靶点、不同突变背景"的患者可能对药物反应迥异，未来精准分层治疗可能比广谱抗炎更为现实，也对新型非磺酰脲骨架抑制剂（如BAL-0028/BGE-102[13]）提出了更高期待。

### 4.2 未来2—3年值得关注的方向

**临床管线数据读出**将是近期最值得关注的事件。dapansutrile在急性痛风（NCT05658575）和2型糖尿病（NCT06047262）的II/III期试验、VTX2735在复发性心包炎（NCT06836232）的II期试验、selnoflast的RIVULET动脉粥样硬化试验（NCT07448038）、NT-0796与司美格鲁肽联合用药试验（NCT07220629），以及BGE-102的QUELL-CV心血管风险II期试验（NCT07656727），均有望在2026—2027年陆续公布结果。其中，痛风试验因机制链条清晰（尿酸盐结晶→NLRP3→IL-1β）而被视为"能否在人体起效"的试金石；BGE-102的QUELL-CV试验以12周hsCRP变化为主要终点，将首次为新型非磺酰脲、非ATP酶抑制型NLRP3抑制剂的人体药效提供关键证据。这些试验将初步回答"直接抑制NLRP3能否在人体中发挥预期的抗炎效果"这一核心问题。CAPS和痛风因机制明确，是最有可能率先读出阳性数据的适应症；若至少一个适应症取得成功，有望为整个领域提供首个人体概念验证参照，并参照Mangan等人（2018）[17]所描绘的路径，推动后续向更大患者群体的慢性病拓展。

**新作用机制与新结合位点**的化合物值得持续关注。oxDNA/PYD通路方面，Cabral等人（2023，*Commun Biol*）[4]和Lackner等人（2025，*Trends Biochem Sci*）[5]的工作提示，靶向NLRP3与氧化线粒体DNA的相互作用或hOGG1通路，可能开辟NACHT域之外的干预策略，但目前尚无基于该机制的临床候选药物。NACHT域新位点方面，Hartman等人（2024，*Bioorg Med Chem Lett*）[12]和Wilhelmsen等人（2024，bioRxiv）[13]报道的BAL-0028及其临床衍生物BGE-102不抑制ATP酶、结合于FISNA附近区域，对部分CAPS突变有效，代表了一类与MCC950正交的抑制模式；值得注意的是，这类化合物对灵长类NLRP3高效但对小鼠NLRP3几乎无效，提示临床前评价必须采用人源化模型或灵长类细胞体系。未来若其冷冻电镜结构正式发表，将为"非CRID3口袋"药物设计提供直接模板。Matico等人（2025，*EMBO Mol Med*）[7]报道的化合物C和NIC系列[6]则证明，即使在同一口袋内，全新骨架也可规避MCC950的脱靶毒性。

**脑穿透性与适应症匹配的精细化策略**也将进一步深化。随着NT-0796、VTX3232、usnoflast和BGE-102等脑穿透性化合物推进神经疾病试验，而外周限制性化合物主攻心血管代谢疾病，"是否需要药物入脑"将成为适应症选择的核心决策变量，而非简单的活性优先。Vande Walle和Lamkanfi（2024）[1]在综述中亦指出，中枢神经系统炎症与系统性代谢炎症在病理机制和可及靶点方面存在显著差异，同一NLRP3靶点在不同组织中的药理需求并不相同。未来管线分化可能进一步加剧：脑穿透性化合物主攻帕金森病、ALS、糖尿病黄斑水肿等需要中枢或视网膜分布的疾病，而外周限制性化合物则聚焦痛风、心包炎和动脉粥样硬化等以系统性炎症为主的适应症。

### 4.3 小结

NLRP3抑制剂领域的基础研究与药物化学已相当深入，多个候选药物进入II期临床，但尚无获批品种；未来2—3年II期数据读出和新型非磺酰脲骨架的临床验证，将决定该领域能否完成从"有希望的靶点"到"可成药疗法"的关键跨越。

---

## 参考文献

[1] Vande Walle L, Lamkanfi M. Drugging the NLRP3 inflammasome: from signalling mechanisms to therapeutic targets. *Nat Rev Drug Discov*. 2024;23(1):43-66. doi:10.1038/s41573-023-00822-2  
链接：https://doi.org/10.1038/s41573-023-00822-2

[2] Yu X, Matico RE, Miller R, et al. Structural basis for the oligomerization-facilitated NLRP3 activation. *Nat Commun*. 2024;15:1164. doi:10.1038/s41467-024-45396-8  
链接：https://doi.org/10.1038/s41467-024-45396-8

[3] Le X, Magupalli VG, Wu H. Cryo-EM structures of the active NLRP3 inflammasome disc. *Nature*. 2023;613(7944):595-600. doi:10.1038/s41586-022-05570-8  
链接：https://doi.org/10.1038/s41586-022-05570-8

[4] Cabral A, Cabral JE, Wang A, et al. Differential binding of NLRP3 to non-oxidized and Ox-mtDNA mediates NLRP3 inflammasome activation. *Commun Biol*. 2023;6:578. doi:10.1038/s42003-023-04817-y  
链接：https://doi.org/10.1038/s42003-023-04817-y

[5] Lackner A, Leonidas L, Macapagal A, Lee H, McNulty R. How interactions between oxidized DNA and the NLRP3 inflammasome fuel inflammatory disease. *Trends Biochem Sci*. 2025;50(11):931-944. doi:10.1016/j.tibs.2025.07.007  
链接：https://doi.org/10.1016/j.tibs.2025.07.007

[6] Vande Walle L, Said MS, Paerewijck O, et al. Novel chemotype NLRP3 inhibitors that target the CRID3-binding pocket with high potency. *Life Sci Alliance*. 2024;7(6):e202402644. doi:10.26508/lsa.202402644  
链接：https://doi.org/10.26508/lsa.202402644

[7] Matico R, Grauwen K, Chauhan D, et al. Navigating from cellular phenotypic screen to clinical candidate: selective targeting of the NLRP3 inflammasome. *EMBO Mol Med*. 2025;17(1):54-84. doi:10.1038/s44321-024-00181-4  
链接：https://doi.org/10.1038/s44321-024-00181-4

[8] Amo-Aparicio J, Daly J, Højen JF, et al. Pharmacologic inhibition of NLRP3 reduces the levels of α-synuclein and protects dopaminergic neurons in a model of Parkinson's disease. *J Neuroinflammation*. 2023;20:147. doi:10.1186/s12974-023-02830-w  
链接：https://doi.org/10.1186/s12974-023-02830-w

[9] Cabral JE, Wu A, Zhou H, Pham MA, Lin S, McNulty R. Targeting the NLRP3 inflammasome for inflammatory disease therapy. *Trends Pharmacol Sci*. 2025;46(6):503-519. doi:10.1016/j.tips.2025.04.007  
链接：https://doi.org/10.1016/j.tips.2025.04.007

[10] Bultinck J, Yuan S, Cantuti-Castelvetri L, et al. NLRP3 inhibition by VTX3232 tempers inflammation resulting in reduced body weight, hyperglycemia, and hepatic steatosis in obese male mice. *Mol Metab*. 2026;103:102282. doi:10.1016/j.molmet.2025.102282  
链接：https://doi.org/10.1016/j.molmet.2025.102282

[11] Coll RC, Hill JR, Day CJ, et al. MCC950 directly targets the NLRP3 ATP-hydrolysis motif for inflammasome inhibition. *Nat Chem Biol*. 2019;15(6):556-559. doi:10.1038/s41589-019-0277-7  
链接：https://doi.org/10.1038/s41589-019-0277-7

[12] Hartman G, Humphries P, Hughes R, et al. The discovery of novel and potent indazole NLRP3 inhibitors enabled by DNA-encoded library screening. *Bioorg Med Chem Lett*. 2024;104:129454. doi:10.1016/j.bmcl.2024.129454  
链接：https://doi.org/10.1016/j.bmcl.2024.129454

[13] Wilhelmsen K, Deshpande A, Tronnes S, et al. Discovery of a potent and selective inhibitor of human NLRP3 with a novel binding modality and mechanism of action. *bioRxiv*. 2024. doi:10.1101/2024.12.21.629867  
链接：https://doi.org/10.1101/2024.12.21.629867

[14] Hochheiser IV, Pilsl M, Hagelueken G, et al. Structure of the NLRP3 decamer bound to the cytokine release inhibitor CRID3. *Nature*. 2022;604(7904):184-189. doi:10.1038/s41586-022-04467-w  
链接：https://doi.org/10.1038/s41586-022-04467-w

[15] Kennedy CR, Tan YS, Muller M, et al. Unraveling the specificity of MCC950 activity, inhibition of the NLRP3 inflammasome and carbonic anhydrase 2. *SLAS Discov*. 2021;26(10):1264-1272. doi:10.1177/24725552211044351  
链接：https://doi.org/10.1177/24725552211044351

[16] Cosson C, Belot A, Lambotte O, et al. Functional diversity of NLRP3 gain-of-function mutants associated with CAPS autoinflammation. *J Exp Med*. 2024;221(5):e20231200. doi:10.1084/jem.20231200  
链接：https://doi.org/10.1084/jem.20231200

[17] Mangan MSJ, Olhava EJ, Roush WR, et al. Targeting the NLRP3 inflammasome in inflammatory diseases. *Nat Rev Drug Discov*. 2018;17(8):588-606. doi:10.1038/nrd.2018.97  
链接：https://doi.org/10.1038/nrd.2018.97
