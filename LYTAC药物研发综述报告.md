# 基于溶酶体靶向嵌合体技术的药物研发综述

## 摘要

靶向蛋白降解是近年来药物发现领域的重要研究方向<sup>[17]</sup>。蛋白水解靶向嵌合体已有多项分子进入临床试验，但其降解途径依赖蛋白酶体，主要适用于具有胞内可及结构域的靶蛋白。2020年，Banik等在Nature上首次系统报道溶酶体靶向嵌合体（LYTAC）<sup>[1]</sup>，通过双功能分子同时结合溶酶体靶向受体与胞外或膜蛋白靶点，将靶蛋白导向溶酶体降解，从而将靶向降解的应用范围拓展至分泌蛋白和跨膜蛋白。此后，围绕受体选择、组织特异性递送和降解效率优化，GalNAc-LYTAC<sup>[2]</sup>、全基因组机制筛选<sup>[3]</sup>以及计算设计的EndoTag<sup>[4]</sup>等代表性工作相继发表。本文在梳理溶酶体靶向嵌合体基本机制的基础上，重点分析四篇发表于Nature、Nature Chemical Biology和Science等期刊的研究实例，讨论该技术在药物研发中的意义与局限，并在总结与展望中提出对该领域后续发展的思考。

**关键词** 溶酶体靶向嵌合体 靶向蛋白降解 溶酶体靶向受体 胞外蛋白 药物研发

---

## 前言

蛋白质表达失调或功能异常与肿瘤、神经退行性疾病、自身免疫病和代谢性疾病等多种病理过程密切相关。长期以来，小分子抑制剂和抗体是靶向蛋白质的主要药物形式。前者通常依赖靶蛋白上的可成药结合口袋，通过竞争性抑制阻断酶活性或蛋白相互作用。后者则通过高亲和力结合实现中和或阻断。这两种方式均属于占据驱动的治疗模式，一般不能降低靶蛋白的绝对含量，对缺乏明确小分子结合位点或主要依赖蛋白蛋白相互作用发挥功能的靶点往往难以有效干预<sup>[16,17]</sup>。

靶向蛋白降解技术提供了一种不同的思路，即通过设计双功能或多功能分子，将靶蛋白引导至细胞内已有的降解机器，以实现靶蛋白的清除而不是单纯的功能抑制<sup>[17]</sup>。在胞内靶向蛋白降解领域，蛋白水解靶向嵌合体（PROTAC）是研究最为成熟的策略之一，相关分子已在血液肿瘤等适应症中推进至临床试验。然而，蛋白酶体途径主要处理具有胞内结构域且可被泛素化修饰的蛋白质。据估计，人类蛋白质组中约40%的蛋白位于细胞外或细胞膜上，包括细胞因子、生长因子、免疫检查点分子和受体酪氨酸激酶等<sup>[16]</sup>。这些靶点无法被传统蛋白水解靶向嵌合体有效触及，在药物研发中长期被视为成药性较低的类别。

溶酶体是细胞内负责大分子物质分解的主要细胞器，其内部含有多种水解酶，在酸性环境中可降解蛋白质、核酸和多糖等物质。细胞可通过受体介导的内吞等途径将底物输送至溶酶体。溶酶体靶向嵌合体正是利用这一生物学过程，在细胞外将靶蛋白与溶酶体靶向受体连接在一起，使靶蛋白随受体内吞并进入溶酶体降解<sup>[1,22]</sup>。自2020年概念提出以来，该领域在配体化学、受体选择和分子构建方式等方面快速拓展，出现了适体LYTAC<sup>[9-11]</sup>、IGF2融合蛋白LYTAC<sup>[5]</sup>、细胞因子受体靶向嵌合体<sup>[12]</sup>以及不依赖受体的SignalTAC<sup>[14]</sup>等多种变体。本文旨在系统梳理溶酶体靶向嵌合体的基本框架，并以四篇一区研究论文为重点，深入分析其机制创新与药物研发意义，最后对该领域的应用前景与转化挑战提出思考。

---

## 正文

### （一）溶酶体靶向嵌合体的基本作用机制

溶酶体靶向嵌合体是一种双功能或多功能分子，其一端结合待降解的胞外或膜蛋白靶点，另一端结合细胞表面的溶酶体靶向受体（LTR）。当分子同时结合靶蛋白和受体后，可在细胞外形成靶蛋白、嵌合体和受体之间的三元复合物。复合物经网格蛋白介导的内吞进入早期内体，在内体酸性环境中受体与配体解离，靶蛋白随内体成熟进入溶酶体被水解酶降解，而受体则可通过循环途径返回细胞表面<sup>[1]</sup>。这一机制与蛋白水解靶向嵌合体依赖泛素化和蛋白酶体降解的路径存在本质区别。PROTAC需要将靶蛋白拉入胞内并连接E3连接酶，而溶酶体靶向嵌合体可在细胞外完成靶蛋白的捕获和分选，对缺乏胞内结构域的膜蛋白和分泌蛋白同样适用。从药物化学角度看，溶酶体靶向嵌合体通常由靶向结合模块、连接子和溶酶体靶向配体三部分组成，其中连接子的长度和柔性会影响三元复合物的形成效率，这一设计原则与PROTAC中连接子优化的问题具有相似之处<sup>[22]</sup>。

**图1** 溶酶体靶向嵌合体作用机制示意图。双功能分子一端通过抗体等结合模块识别靶蛋白，另一端通过M6Pn糖肽等配体结合CI-M6PR，形成三元复合物后经受体介导内吞进入溶酶体降解。（引自 Banik S M, et al. Lysosome-targeting chimaeras for degradation of extracellular proteins. *Nature*, 2020, 584(7820): 291-297, Figure 1a，改绘）

目前已知的溶酶体靶向受体包括CI-M6PR/IGF2R、ASGPR以及清道夫受体、细胞因子受体等多种类型<sup>[15,22]</sup>。CI-M6PR在多种组织中广泛表达，可识别带有甘露糖6磷酸（M6P）修饰的配体，是第一代溶酶体靶向嵌合体主要使用的受体。ASGPR主要在肝细胞表达，识别末端带有N-乙酰半乳糖胺（GalNAc）的糖蛋白，是第二代组织特异性平台的基础。不同受体的组织分布和内吞激活方式存在差异，有的受体需要构象变化才能启动内吞，有的则需要多价配体诱导受体簇集<sup>[4,15]</sup>。这些差异直接决定了溶酶体靶向嵌合体在靶点选择和分子设计上的策略。

### （二）研究实例一：CI-M6PR溶酶体靶向嵌合体的概念建立（Banik等，Nature，2020）

2020年，Banik等在Nature上发表的研究是溶酶体靶向嵌合体领域的奠基性工作<sup>[1]</sup>。该研究由Stanford大学Carolyn Bertozzi课题组完成，首次将靶向蛋白降解的嵌合体思路从蛋白酶体途径拓展至溶酶体途径，在方法学上具有开创意义。

在该研究中，作者将人工合成的M6Pn寡糖肽与抗EGFR抗体西妥昔单抗等靶向模块通过化学偶联相连，构建出第一代CI-M6PR溶酶体靶向嵌合体。M6Pn是CI-M6PR的高亲和力配体，而CI-M6PR的生理功能之一是识别带有M6P修饰的溶酶体酶，并介导其向溶酶体的转运。实验结果显示，M6Pn-西妥昔单抗偶联物可在HeLa、A549和MDA-MB-231等多种细胞系中有效降低EGFR蛋白水平，降解效果具有浓度和时间依赖性。作者还构建了靶向PD-L1的溶酶体靶向嵌合体，在肿瘤细胞中实现了PD-L1的降解，提示该平台对免疫检查点靶点同样适用。通过共聚焦显微镜和溶酶体标志物LysoTracker的共定位实验，作者证实被降解的靶蛋白确实进入了溶酶体区室，而非通过其他途径清除。更重要的是，该研究验证了ApoE4等分泌型蛋白可被导向溶酶体清除，表明溶酶体靶向嵌合体不仅适用于膜蛋白，也可用于胞外游离蛋白的降解。ApoE4与阿尔茨海默病风险相关，其清除策略具有一定的疾病研究价值。

该研究还利用溶酶体靶向嵌合体平台开展了CRISPR干扰筛选，初步揭示了CI-M6PR介导货物内吞所依赖的部分细胞通路，包括网格蛋白、适配蛋白复合物以及V-ATP酶等内吞和内涵体酸化相关因子。这一筛选策略为后续Ahn等在Science上的全基因组敲除筛选提供了方法学基础<sup>[3]</sup>。作者指出，与PROTAC相比，溶酶体靶向嵌合体不需要靶蛋白具有胞内结构域，也不需要泛素化修饰，因此在靶点类型上具有明显优势。然而，该工作也揭示了第一代平台的若干局限。化学合成和随机偶联可能导致产物异质性较高，M6Pn糖肽的制备步骤相对复杂。此外，CI-M6PR在全身多种组织表达，系统给药时缺乏组织选择性。细胞表面相当比例的CI-M6PR已被内源性M6P糖蛋白占据，外源性溶酶体靶向嵌合体需与内源配体竞争受体结合位点，这可能限制降解效率的上限<sup>[3,22]</sup>。

**图2** CI-M6PR溶酶体靶向嵌合体降解膜蛋白与分泌蛋白的实验验证。左图为EGFR和PD-L1降解的Western blot结果，右图为ApoE4分泌蛋白的溶酶体清除。（引自 Banik S M, et al. *Nature*, 2020, 584(7820): 291-297, Figure 2b,c，改绘）

### （三）研究实例二：GalNAc-LYTAC与组织特异性蛋白降解（Ahn等，Nature Chemical Biology，2021）

CI-M6PR在多种组织中广泛表达，基于该受体的溶酶体靶向嵌合体在实现广谱胞外靶点降解的同时，也面临组织选择性不足的问题。2021年，Ahn等在Nature Chemical Biology上发表的研究针对这一问题，开发了基于去唾液酸糖蛋白受体（ASGPR）的GalNAc-LYTAC<sup>[2]</sup>。该期刊同属Nature Portfolio，在生物化学与化学生物学领域具有较高的学术影响力。

ASGPR是一种主要在肝细胞表达的C型凝集素受体，识别末端带有GalNAc残基的糖蛋白，并通过内吞途径将其递送至溶酶体。GalNAc-LYTAC的设计思路是将三价GalNAc配体与抗靶蛋白抗体或Fab偶联，使分子在结合靶蛋白的同时被肝细胞表面的ASGPR捕获，从而实现肝脏选择性的蛋白降解。这一设计借鉴了GalNAc偶联核酸药物和GalNAc-抗体偶联物在肝脏递送中的成功经验，ASGPR对多价GalNAc的高亲和力是实现肝细胞靶向的关键。在该研究中，作者将GalNAc-LYTAC应用于EGFR等靶点，在HepG2等肝细胞系中实现了有效的靶蛋白清除，而在不表达ASGPR的非肝细胞系中降解效果显著减弱，验证了组织选择性。与单纯抗体阻断相比，GalNAc-LYTAC对EGFR下游信号的抑制更为彻底，表明蛋白降解相比功能阻断可能具有更持久的效应。这是因为抗体阻断是可逆的占据过程，而降解是降低靶蛋白的绝对含量。

该团队还开发了在抗体Fc区定点引入单个tri-GalNAc的均一偶联策略，改善了GalNAc-LYTAC在体内的药代动力学表现。与随机偶联的多GalNAc产物相比，均一偶联物的结构定义更清晰，有利于后续的药代动力学和毒理学评价。这一工作的意义在于，它将溶酶体靶向嵌合体从广谱降解推向了组织特异性降解。ASGPR在肝脏高表达的特点，使GalNAc-LYTAC在肝脏相关靶点和分泌蛋白清除方面具有天然优势。Spiegel课题组同期在Nature Chemical Biology上报道的MoDE-A<sup>[8]</sup>以及Tang课题组在ACS Central Science上报道的Glyco-LYTAC<sup>[7]</sup>也采用了类似的ASGPR策略，进一步验证了该路线的可行性。在代谢性疾病方向，Bagdanoff等后续报道的MoDE-A类分子通过ASGPR介导PCSK9的溶酶体清除，在动物实验中显示了降低血浆PCSK9和胆固醇相关指标的潜力<sup>[18]</sup>。

**图3** GalNAc-LYTAC结构及其在肝细胞中降解EGFR的机制示意图。三价GalNAc配体结合ASGPR，抗体结合靶蛋白，实现肝脏选择性降解。（引自 Ahn G, et al. LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation. *Nature Chemical Biology*, 2021, 17(9): 937-946, Figure 1a，改绘）

### （四）研究实例三：溶酶体靶向嵌合体降解效率的细胞调节机制（Ahn等，Science，2023）

尽管溶酶体靶向嵌合体在多种细胞模型中显示了降解活性，但不同细胞系和靶点之间的降解效率差异较大，其背后的细胞生物学机制长期不够清楚。2023年，Ahn等在Science上发表的研究通过全基因组CRISPR敲除筛选，系统鉴定了影响溶酶体靶向嵌合体介导膜蛋白降解的细胞内调节因子<sup>[3]</sup>。Science期刊的发表表明，该机制问题已受到广泛关注。

该研究以靶向EGFR的M6Pn-溶酶体靶向嵌合体为工具，在多种细胞系中进行了CRISPR-Cas9全基因组敲除筛选，寻找影响降解效率的正向和负向调节基因。筛选结果通过迭代验证和独立重复实验进行了确认，保证了命中基因的可靠性。结果发现，逆转运复合体retromer是降解的重要负向调节因子。Retromer的功能是将已内吞的受体和货物从早期内涵体回收至细胞膜，这一回收过程与溶酶体降解方向相反。敲低retromer复合物的核心组分VPS35或VPS26A可显著增强EGFR等靶蛋白的降解，表明内涵体分选方向的平衡是决定降解效率的关键环节。研究还发现，CUL3的NEDD8化修饰与晚期内涵体成熟有关，其表达水平与溶酶体靶向嵌合体活性呈正相关。CUL3是泛素连接酶家族成员，其NEDD8化修饰可能通过促进内涵体成熟来增强溶酶体定向。此外，抑制细胞内M6P生物合成酶GNPTAB可增加细胞表面CI-M6PR的可用比例，促进溶酶体靶向嵌合体与受体结合及后续内化。GNPTAB负责将M6P标记加至溶酶体酶前体上，其敲除减少了内源性M6P糖蛋白对CI-M6PR的占据，从而释放了更多受体供外源性溶酶体靶向嵌合体使用。作者还通过蛋白质组学分析证实，GNPTAB敲除后细胞表面CI-M6PR的可用水平确实有所升高。

**图4** 全基因组CRISPR筛选揭示的溶酶体靶向嵌合体降解调节通路。Retromer介导复合物回膜，GNPTAB影响CI-M6PR可用性，CUL3参与内涵体成熟。（引自 Ahn G, et al. Elucidating the cellular determinants of targeted membrane protein degradation by lysosome-targeting chimeras. *Science*, 2023, 382(6669): eadf6249, Figure 2a，改绘）

这项研究的重要性在于，它将溶酶体靶向嵌合体的优化从分子化学层面拓展到了细胞通路层面。Retromer和GNPTAB等因子提示，降解效率不仅取决于嵌合体本身的亲和力，还取决于细胞内吞和分选网络的平衡状态。从药物研发角度看，这意味着同一溶酶体靶向嵌合体分子在不同细胞类型或组织中的降解效率可能存在较大差异，临床前评价需要在多种相关细胞模型中进行。该研究也为联合策略提供了理论依据，例如通过小分子抑制retromer功能或调节M6P生物合成来增强降解效果，尽管这些策略的体内安全性尚需验证。

### （五）研究实例四：EndoTag计算设计内吞诱导蛋白（Huang等，Nature，2024）

内源性配体与人工配体竞争同一受体、不同受体的内吞机制差异以及分子异质性，是制约溶酶体靶向嵌合体降解效率的重要因素。2024年，Huang等在Nature上发表的EndoTag平台采用Rosetta计算蛋白质设计方法，针对不同受体的内吞特征设计结合蛋白，代表了该领域向理性设计方向的重要进展<sup>[4]</sup>。

EndoTag的设计逻辑是根据不同受体的内吞激活方式采取不同策略，这一思路区别于此前简单地将天然配体化学偶联到靶向模块上的做法。对于持续在细胞膜与内涵体之间循环的受体，EndoTag结合于与天然配体不重叠的位点，以减少与内源配体的竞争。对于需要构象变化才能内吞的IGF2R，EndoTag采用双表位结合策略，同时结合受体的两个非重叠表位，诱导有利于内吞的受体构象变化。对于需要受体簇集才能高效内吞的ASGPR，EndoTag采用多价设计，通过一个EndoTag分子结合多个受体单体，促进受体在膜表面的聚集。这些设计均经过Rosetta蛋白质设计平台生成，并通过酵母展示和细胞实验进行了筛选验证。EndoTag与抗EGFR或抗PD-L1等靶向模块融合后可形成蛋白型溶酶体靶向嵌合体，在细胞实验中实现了相应靶蛋白的降解。在小鼠肿瘤模型中，EndoTag-抗PD-L1融合蛋白显示出优于单纯抗PD-L1抗体的抗肿瘤效果，表明降解策略在功能上优于单纯的免疫检查点阻断。

**图5** 针对不同溶酶体靶向受体的EndoTag设计策略。左图为IGF2R双表位结合诱导构象变化，中图为ASGPR多价簇集，右图为循环受体的非竞争结合位点。（引自 Huang P S, et al. Designed endocytosis-inducing proteins degrade targets and amplify signals. *Nature*, 2024, 635(8038): 903-910, Figure 1a-c，改绘）

EndoTag的全基因编码特性也为其作为可编程降解平台提供了基础。与依赖化学合成糖肽的第一代平台相比，EndoTag通过蛋白质设计规避了人工糖肽的制备复杂性和潜在的免疫原性风险。该工作表明，理解并利用不同受体的内吞机制差异，可能是提高降解效率的一条可行路径，而不必完全依赖天然配体的化学模拟。Liu等在Journal of Hematology and Oncology上的综述也将EndoTag归入胞外靶向蛋白降解技术的重要进展<sup>[15]</sup>。

**表1** 四个代表性溶酶体靶向嵌合体研究实例的比较。

| 研究 | 期刊 | 年份 | 核心贡献 | 主要局限 |
| --- | --- | --- | --- | --- |
| Banik等<sup>[1]</sup> | Nature | 2020 | 首次建立CI-M6PR溶酶体靶向嵌合体概念 | 产物异质性高，缺乏组织选择性 |
| Ahn等<sup>[2]</sup> | Nature Chemical Biology | 2021 | 开发GalNAc-LYTAC实现肝脏选择性降解 | 适用范围限于表达ASGPR的组织 |
| Ahn等<sup>[3]</sup> | Science | 2023 | 揭示retromer等降解效率调节因子 | 尚未转化为可操作的优化策略 |
| Huang等<sup>[4]</sup> | Nature | 2024 | 计算设计EndoTag提高降解效率 | 设计成本较高，长期安全性待评价 |

（表格内容根据 Banik等<sup>[1]</sup>、Ahn等<sup>[2,3]</sup>、Huang等<sup>[4]</sup> 文献整理）

### （六）其他技术进展

除上述四个代表性实例外，溶酶体靶向嵌合体领域还出现了多条补充技术路线。Zhang等报道的iLYTAC平台将IGF2与靶向模块通过基因工程融合表达，形成完全基因编码的降解分子，避免了复杂化学偶联<sup>[5]</sup>。Pance等报道的KineTAC采用细胞因子CXCL12结合CXCR7等受体作为内吞效应器，介导PD-L1、EGFR和VEGF等靶蛋白的降解，拓展了效应器受体的选择范围<sup>[12]</sup>。Wu等和Miao等则分别发展了适体连接的GalNAc-LYTAC和CI-M6PR适体嵌合体，以核酸适体替代抗体作为靶向模块，降低了分子量并简化了部分制备步骤<sup>[9,10]</sup>。Chen等报道的SignalTAC通过溶酶体分选信号肽而非溶酶体靶向受体介导降解，为不表达合适受体的细胞提供了替代思路<sup>[14]</sup>。Kim等报道的LYTACgyM6pG平台使用糖工程酵母来源的M6P糖，在配体制备均一性方面进行了改进<sup>[6]</sup>。Zhang等报道的化学酶法M6P糖链定点修饰策略则改善了抗体偶联的均一性<sup>[13]</sup>。这些工作共同丰富了溶酶体靶向嵌合体的工具箱，但大多仍处于临床前概念验证阶段，尚未进入临床试验。

---

## 总结和展望

溶酶体靶向嵌合体技术的提出，在靶向蛋白降解领域具有明确的方法学意义。PROTAC已将大量胞内靶点纳入可降解范围，但分泌蛋白和跨膜蛋白长期以来仍是药物研发的难点。Banik等<sup>[1]</sup>的工作表明，利用细胞固有的溶酶体清除途径，可以在不依赖泛素化和蛋白酶体的前提下实现上述靶点的清除。这一思路为免疫检查点、受体酪氨酸激酶、分泌型细胞因子以及ApoE4等病理相关蛋白的干预提供了新的可能，其核心价值在于将原本只能阻断但难以清除的靶点转化为可主动消除的对象。

从四个研究实例的递进关系来看，该领域的发展呈现出由概念验证到组织特异性优化、由分子设计到细胞机制解析、再由化学模拟到计算理性设计的演进趋势。GalNAc-LYTAC<sup>[2]</sup>回应了第一代平台组织选择性不足的问题，Science上的机制研究<sup>[3]</sup>揭示了降解效率的细胞生物学瓶颈，EndoTag<sup>[4]</sup>则尝试从受体结构出发进行从头设计。这种从广度到深度、从化学到生物、从模仿到设计的推进路径，反映了溶酶体靶向嵌合体作为一个新兴技术平台正在逐步走向成熟。四个实例均发表于Nature、Nature Chemical Biology或Science，说明该方向已获得顶级期刊的持续关注，其科学问题的重要性和创新性得到了同行认可。

结合上述文献，笔者认为该领域后续发展仍面临若干关键问题。降解效率的预测与标准化评价目前仍不充分，不同研究采用的细胞系、靶蛋白和检测方法差异较大，Ahn等<sup>[3]</sup>揭示的retromer和GNPTAB等调节因子提示降解效率高度依赖细胞背景，目前尚缺乏统一的可比评价标准<sup>[22]</sup>。组织特异性与广谱降解之间的平衡也尚未解决，CI-M6PR路线适用于多种靶点但缺乏器官选择性，ASGPR路线局限于肝脏，脑、肿瘤微环境等组织的特异性受体仍有待开发。从临床前到临床的转化路径同样尚不清晰，与PROTAC已有多个分子进入临床试验相比，溶酶体靶向嵌合体在公开报道中仍以细胞和动物模型为主，人工糖肽的免疫原性、复杂偶联产物的均一性控制以及规模化生产成本，都是成药化需要面对的实际问题<sup>[16]</sup>。

展望未来，提高降解的空间和时间特异性仍是重要的研究课题。EndoTag等计算设计平台与基因编码的iLYTAC<sup>[5]</sup>相结合，可能为降低制备复杂度和改善批次一致性提供路径。Ahn等<sup>[3]</sup>的全基因组筛选工作也提示，结合蛋白质组学和更系统的细胞生物学研究，有望建立降解效率的预测模型。在应用层面，肝脏靶向的GalNAc-LYTAC在代谢性疾病和肝相关靶点方面具有相对明确的开发方向，而肿瘤免疫和神经退行性疾病等领域的验证仍需更多临床前数据支撑。

总体而言，溶酶体靶向嵌合体为靶向胞外和膜蛋白提供了一条与蛋白水解靶向嵌合体互补的技术路径。该技术已在多种细胞和小鼠模型中证明了靶向清除蛋白的可行性，其意义不仅在于提供了一种新的分子工具，更在于从概念上拓宽了靶向蛋白降解的边界。然而，其能否成为可广泛用于药物开发的成熟模态，仍有赖于后续研究在特异性、安全性、可制造性和临床疗效等方面给出更充分的证据。对于药物研发实践而言，溶酶体靶向嵌合体目前更适合作为针对特定胞外靶点的探索性工具，而非已可广泛推广的成药性平台。随着机制研究的深入和设计工具的完善，这一判断可能在未来数年内发生改变。

---

## 参考文献

[1] Banik S M, Pedram K, Wisnovsky S, et al. Lysosome-targeting chimaeras for degradation of extracellular proteins[J]. Nature, 2020, 584(7820): 291-297. https://doi.org/10.1038/s41586-020-2545-9

[2] Ahn G, Banik S M, Miller C L, et al. LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation[J]. Nature Chemical Biology, 2021, 17(9): 937-946. https://doi.org/10.1038/s41589-021-00770-1

[3] Ahn G, Riley N M, Kamber R A, et al. Elucidating the cellular determinants of targeted membrane protein degradation by lysosome-targeting chimeras[J]. Science, 2023, 382(6669): eadf6249. https://doi.org/10.1126/science.adf6249

[4] Huang P S, Boyken S E, Baker D, et al. Designed endocytosis-inducing proteins degrade targets and amplify signals[J]. Nature, 2024, 635(8038): 903-910. https://doi.org/10.1038/s41586-024-07948-2

[5] Zhang B, Brahma R K, Zhu L, et al. Insulin-like growth factor 2 (IGF2)-fused lysosomal targeting chimeras for degradation of extracellular and membrane proteins[J]. Journal of the American Chemical Society, 2023, 145(42): 24272-24283. https://doi.org/10.1021/jacs.3c08886

[6] Kim S, Kang J, Bi A D, et al. Lysosome-targeting chimera using Mannose-6-phosphate glycans derived from glyco-engineered yeast[J]. Bioconjugate Chemistry, 2025, 36(3): 424-436. https://doi.org/10.1021/acs.bioconjchem.4c00512

[7] Zhou Y, Wang H, Fang Y, et al. Glyco-LYTACs enable targeted degradation of extracellular proteins through the asialoglycoprotein receptor[J]. ACS Central Science, 2021, 7(4): 722-733. https://doi.org/10.1021/acscentsci.1c00122

[8] Caianiello D F, Petter R C, Wolf P, et al. Degradation from the outside in: targeting extracellular and membrane-associated proteins for degradation through the endolysosome pathway[J]. Nature Chemical Biology, 2021, 17(8): 947-953. https://doi.org/10.1038/s41589-021-00773-y

[9] Miao Y, Gao Q, Mao M, et al. Bispecific aptamer chimeras enable targeted protein degradation on cell membranes[J]. Angewandte Chemie International Edition, 2021, 60(20): 11267-11271. https://doi.org/10.1002/anie.202100600

[10] Wu Y, Lin B, Lu Y, et al. Aptamer-LYTACs for targeted degradation of extracellular and membrane proteins[J]. Angewandte Chemie International Edition, 2023, 62(15): e202218106. https://doi.org/10.1002/anie.202218106

[11] Li Y, Liu X, Yu L, et al. Covalent LYTAC enabled by DNA aptamers for immune checkpoint degradation therapy[J]. Journal of the American Chemical Society, 2023, 145(45): 24506-24521. https://doi.org/10.1021/jacs.3c08562

[12] Pance K, Gramespacher J A, Byrnes J R, et al. Modular cytokine receptor-targeting chimeras for targeted degradation of cell surface and extracellular proteins[J]. Nature Biotechnology, 2023, 41(2): 273-281. https://doi.org/10.1038/s41587-022-01594-9

[13] Zhang X, Liu H, He J, et al. Site-specific chemoenzymatic conjugation of high-affinity M6P glycan ligands to antibodies for targeted protein degradation[J]. ACS Chemical Biology, 2022, 17(12): 3013-3023. https://doi.org/10.1021/acschembio.2c00658

[14] Chen W, Zhang S, Li L, et al. Lysosome-targeting chimeras containing an endocytic signaling motif trigger endocytosis and lysosomal degradation of cell-surface proteins[J]. Chemical Science, 2024, 15(42): 17682-17693. https://doi.org/10.1039/D4SC05093B

[15] Liu Y, Koval A, Vazquez-Lombardi R, et al. Targeted degradation of extracellular proteins: state of the art and diversity of degrader designs[J]. Journal of Hematology & Oncology, 2025, 18: 45. https://doi.org/10.1186/s13045-025-01703-4

[16] Lin J, Jin J, Shen Y, et al. Emerging protein degradation strategies: expanding the scope to extracellular and membrane proteins[J]. Theranostics, 2021, 11(17): 8337-8349. https://doi.org/10.7150/thno.62559

[17] Zhao L, Zhao J, Zhong K, et al. Targeted protein degradation: mechanisms, strategies and application[J]. Signal Transduction and Targeted Therapy, 2022, 7: 113. https://doi.org/10.1038/s41392-022-00966-4

[18] Bagdanoff J T, Smith T M, Allan M, et al. Clearance of plasma PCSK9 via the asialoglycoprotein receptor mediated by heterobifunctional ligands[J]. Cell Chemical Biology, 2023, 30(1): 97-109. https://doi.org/10.1016/j.chembiol.2022.11.006

[19] Liu Z Q, Deng Q Q, Qin G, et al. Biomarker-activated multifunctional lysosome-targeting chimeras mediated selective degradation of extracellular amyloid fibrils[J]. Chem, 2023, 9(7): 2016-2038. https://doi.org/10.1016/j.chempr.2023.04.014

[20] Duan Q, Jia H R, Chen W, et al. Multivalent aptamer-based lysosome-targeting chimeras platform for mono- or dual-targeted proteins degradation on cell surface[J]. Advanced Science, 2024, 11(15): 2308924. https://doi.org/10.1002/advs.202308924

[21] Tian Y, Miao Y, Guo P, et al. Insulin-like growth factor 2-tagged aptamer chimeras modular assembly for targeted and efficient degradation of two membrane proteins[J]. Angewandte Chemie International Edition, 2024, 63(5): e202316089. https://doi.org/10.1002/anie.202316089

[22] Li Y Y, Yang Y, Zhang R S, et al. Targeted degradation of membrane and extracellular proteins with LYTACs[J]. Acta Pharmacologica Sinica, 2025, 46(1): 1-7. https://doi.org/10.1038/s41401-024-01364-y
