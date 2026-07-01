#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate LYTAC review document with specified formatting."""

from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn


def set_run_font(run, chinese_font="仿宋", western_font="Times New Roman", size_pt=12, bold=False):
    run.font.name = western_font
    run.font.size = Pt(size_pt)
    run.font.bold = bold
    run._element.rPr.rFonts.set(qn("w:eastAsia"), chinese_font)


def add_paragraph(doc, text, align=WD_ALIGN_PARAGRAPH.JUSTIFY, first_line_indent=True, bold=False, size_pt=12):
    p = doc.add_paragraph()
    p.alignment = align
    pf = p.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    pf.space_after = Pt(6)
    if first_line_indent:
        pf.first_line_indent = Cm(0.74)
    run = p.add_run(text)
    set_run_font(run, size_pt=size_pt, bold=bold)
    return p


def add_heading_text(doc, text, level=1):
    sizes = {1: 16, 2: 14, 3: 12}
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER if level == 1 else WD_ALIGN_PARAGRAPH.LEFT
    pf = p.paragraph_format
    pf.space_before = Pt(12)
    pf.space_after = Pt(6)
    pf.first_line_indent = Cm(0)
    pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    run = p.add_run(text)
    set_run_font(run, size_pt=sizes.get(level, 12), bold=True)
    return p


def add_figure_placeholder(doc, caption):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.first_line_indent = Cm(0)
    p.paragraph_format.space_before = Pt(6)
    run = p.add_run("【此处插入示意图】")
    set_run_font(run, size_pt=12)
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap.paragraph_format.first_line_indent = Cm(0)
    cap_run = cap.add_run(caption)
    set_run_font(cap_run, size_pt=10.5)


def add_reference(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    pf = p.paragraph_format
    pf.first_line_indent = Cm(0)
    pf.left_indent = Cm(0.74)
    pf.hanging_indent = Cm(0.74)
    pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    pf.space_after = Pt(3)
    run = p.add_run(text)
    set_run_font(run, size_pt=10.5)


CONTENT = {
    "abstract": (
        "靶向蛋白降解是近年来药物研发领域关注较多的一个方向。"
        "与传统抑制或阻断靶蛋白功能不同，这类技术试图把致病蛋白直接送进细胞内的降解系统。"
        "溶酶体靶向嵌合体（Lysosome-Targeting Chimera, LYTAC）是其中较新的一类，"
        "主要利用细胞表面的溶酶体靶向受体，把胞外蛋白和膜蛋白带入溶酶体降解。"
        "由于PROTAC主要作用于胞内蛋白，LYTAC在一定程度上补上了胞外靶点这一空白。"
        "本文在阅读相关文献的基础上，介绍LYTAC的基本作用方式，"
        "并结合Banik等（Nature, 2020）、Ahn等（Nature Chemical Biology, 2021和Science, 2023）"
        "以及Pance等（Nature Biotechnology, 2023）等研究，讨论该技术的发展过程和可能的应用方向。"
        "笔者认为，LYTAC能否真正走向临床，不只取决于分子怎么设计，"
        "还取决于对内吞、受体循环和溶酶体分选等生物学过程的理解。"
        "目前该领域仍以基础研究为主，临床证据仍然较少，这也是后续需要继续跟踪的问题。"
        "本文写作过程中主要参考了Nature、Science、Nature Chemical Biology、"
        "Nature Biotechnology以及几篇近年的综述文章，"
        "尽量只引用已经正式发表的研究成果。"
    ),
    "intro": [
        "很多疾病都和蛋白质表达异常或功能失调有关，比如肿瘤、自身免疫病和神经退行性疾病。"
        "传统药物研发往往通过小分子或抗体去结合靶蛋白的活性位点或结合位点，从而抑制其功能。"
        "但对于一些没有明确药物结合口袋的蛋白，或者像自身抗体这类靶点，单纯阻断往往效果有限，"
        "因此常被称为较难成药的靶点。",
        "近二十年来，靶向蛋白降解技术发展较快。PROTAC和分子胶等策略已经能把部分胞内蛋白"
        "导向泛素-蛋白酶体系统降解，其中一些分子已进入临床试验。"
        "Bekes等（2022）在Nature Reviews Drug Discovery的综述中回顾了PROTAC的发展历程，"
        "并指出蛋白降解药物正在从概念验证走向临床验证。",
        "不过PROTAC通常要求靶蛋白具有可被配体结合的胞内区域，这样才能招募E3连接酶。"
        "而分泌到细胞外的蛋白和跨膜蛋白数量很多，在疾病中也十分常见，"
        "却一直不在传统PROTAC的主要作用范围内。",
        "从生物医学角度看，这类蛋白参与细胞通讯、免疫识别和信号转导。"
        "例如PD-L1、EGFR、IgE以及多种自身抗体相关复合物，都属于胞外或膜蛋白。"
        "它们有时表达量高，有时更新较快，单纯阻断未必能长期控制疾病进程。"
        "因此，能否把这些蛋白送入溶酶体彻底降解，是一个既有科学问题也有应用价值的问题。",
        "2020年，Bertozzi团队在Nature上提出了溶酶体靶向嵌合体（LYTAC）的概念。"
        "此后几年里，不同研究组在分子设计、受体选择、细胞机制和疾病模型等方面做了不少工作。"
        "Li等（2025）和Chen等（2023）发表的综述文章，"
        "对LYTAC技术的发展脉络和应用前景做了较系统的介绍。"
        "Bekes等（2022）在讨论PROTAC时也提到，"
        "胞外蛋白降解是靶向蛋白降解领域正在拓展的方向。"
        "本文是在阅读文献基础上整理的学习笔记式综述，"
        "主要想说明生物医学认识怎样影响LYTAC这类新药的研发思路，"
        "并梳理几篇代表性论文的主要结论。",
        "在写作过程中，我尽量按照课程要求选择正式发表的一区或高水平综述论文，"
        "避免引用未正式发表的材料，以保证文献来源可靠。",
    ],
    "body1": [
        "LYTAC本质上是一种双功能分子，一端结合靶蛋白，另一端结合细胞表面的溶酶体靶向受体"
        "（Lysosome-Targeting Receptor, LTR），中间通过连接子相连。"
        "靶蛋白结合部分可以用小分子、肽、抗体或适体，"
        "受体结合部分则常用甘露糖-6-磷酸聚糖或GalNAc等配体。"
        "当LYTAC同时拉住靶蛋白和受体后，会形成三元复合物，"
        "随后通过内吞进入内体，再进一步分选到溶酶体，靶蛋白最终被水解（图1）。",
        "Chen等（2023）在Trends in Pharmacological Sciences上发表的综述指出，"
        "LYTAC属于依赖溶酶体的胞外蛋白降解策略，与PROTAC的作用通路不同。",
        "Sun等（2024）在Signal Transduction and Targeted Therapy的综述中"
        "也把LYTAC列为近年来靶向蛋白降解的重要进展之一。",
        "这个过程和抗体单纯阻断靶蛋白不同。抗体通常只是占据结合位点，"
        "而LYTAC的目标是减少靶蛋白本身的含量。"
        "文献中一般认为，这类分子属于事件驱动型作用方式，"
        "也就是不要求结合剂本身具有抑制活性，只要能抓住靶蛋白即可。"
        "另外，由于降解发生在溶酶体，靶蛋白的信号传导功能和蛋白相互作用都可能随之减弱。",
        "目前已报道的LTR包括CI-M6PR（也称IGF2R）、ASGPR、整合素以及部分细胞因子受体等。"
        "Zhou等（2024）在Journal of Medicinal Chemistry讨论了LYTAC可选择的受体空间，"
        "说明受体类型正在不断增加。"
        "不同受体在组织分布和内吞特性上差别较大，"
        "所以受体选择会直接影响LYTAC作用的位置和效率。"
        "这一点在药物设计时不能只看靶蛋白本身，还要考虑受体在哪些细胞上表达。",
        "对我个人而言，LYTAC最难理解的部分并不是化学结构，"
        "而是受体被占用、复合物被回收以及内体怎样分选到溶酶体这些细胞过程。",
        "只有把分子设计和细胞生物学结合起来，才更容易看懂相关论文中的实验结果。",
    ],
    "body2": [
        "从文献来看，LYTAC的效果并不只由化学结构决定，细胞本身的生物学环境也很重要。"
        "在设计分子时，至少需要考虑靶蛋白结合剂、受体配体和连接子三个方面。",
        "靶蛋白结合剂需要有足够亲和力，但不要求一定抑制靶蛋白活性。"
        "膜蛋白常用抗体片段或纳米抗体，分泌蛋白可用抗体或肽段，"
        "部分研究则使用DNA适体作为结合模块。"
        "结合位点和亲和力的小幅变化，就可能影响三元复合物形成和内吞效率。",
        "受体配体的选择和价态同样关键。第一代LYTAC多使用CI-M6PR配体，"
        "通常采用多价甘露糖-6-磷酸聚糖以提高结合能力。"
        "Ahn等后来改用肝特异性受体ASGPR的tri-GalNAc配体，"
        "说明组织分布是决定降解位置的重要因素。"
        "连接子长度和柔性会影响两个功能模块之间的空间排布，"
        "过长或过短都可能降低降解效率。"
        "近年来也有研究通过定点偶联提高抗体-配体偶联物的均一性，"
        "以改善体内稳定性。Caianiello等（2021）在Nature Chemical Biology还报道了"
        "另一种降解胞外蛋白的小分子策略MoDE，"
        "说明除了经典LYTAC，相关思路还在向不同分子类型扩展。",
        "2023年Ahn等在Science上发表的工作进一步说明，"
        "细胞内多条通路会影响LYTAC降解效率，"
        "例如retromer介导的循环回收、CUL3的neddylation修饰，"
        "以及M6P生物合成对CI-M6PR受体的占用等。"
        "这说明同一种LYTAC在不同细胞系中效果可能差很多，"
        "也提示后续药物研发需要更多细胞生物学层面的支持。",
        "Li等（2025）在总结影响降解效率的因素时，"
        "还提到了配体价态、连接子长度和结合表位等分子层面的变量，"
        "说明分子设计和细胞环境需要同时考虑。",
    ],
    "case1": [
        "Banik等在Nature（2020, 584卷7820期, 291-297页）发表了LYTAC领域的开创性工作，"
        "论文题目为Lysosome-targeting chimaeras for degradation of extracellular proteins。"
        "这项研究较早较完整地证明，可以把胞外和膜蛋白导向溶酶体降解。",
        "Li等（2025）在Acta Pharmacologica Sinica的综述中，"
        "把这项工作称为LYTAC技术发展的起点。",
        "作者将靶向模块与CI-M6PR配体通过连接子连接，构建了首批LYTAC分子。"
        "他们还利用CRISPR筛选研究了LYTAC的内吞通路，"
        "发现外泌体复合体（exocyst complex）参与这一过程，"
        "说明这项研究并不只是做分子，也在补充受体介导内吞的基础认识。",
        "实验中，LYTAC成功降解了ApoE4、EGFR、CD71和PD-L1等蛋白。"
        "以EGFR为例，处理后受体水平在数小时内明显下降，"
        "而溶酶体抑制剂可以阻断这一效应，说明溶酶体途径是主要降解路径。"
        "PD-L1被降解后，肿瘤细胞的免疫逃逸能力也有所减弱。"
        "这篇论文的重要性在于，它把靶向蛋白降解从胞内拓展到了胞外和膜蛋白，"
        "也为后续很多研究提供了基本框架。",
        "作者还测试了ApoE4的降解。ApoE4是阿尔茨海默病研究中较受关注的蛋白，"
        "这说明LYTAC并不只面向肿瘤靶点，"
        "在神经退行性疾病相关蛋白上也有探索价值。",
    ],
    "case2": [
        "第一代LYTAC主要依赖CI-M6PR，而该受体在多种组织中都有表达，"
        "因此组织选择性不够理想。2021年，Ahn等在Nature Chemical Biology"
        "（17卷937-946页）报道了一种结合ASGPR的GalNAc-LYTAC。"
        "ASGPR主要在肝细胞表达，此前也已被用于肝靶向核酸递送，"
        "因此这项研究的意义在于把LYTAC进一步推向组织特异性降解。",
        "研究者将靶蛋白结合剂与tri-GalNAc配体连接，在肝细胞中实现了EGFR降解。"
        "结果显示，与单纯使用抑制性抗体相比，GalNAc-LYTAC对EGFR下游信号的抑制更持久。"
        "作者还把较小的肽段结合剂与tri-GalNAc连接，实现了整合素降解并抑制癌细胞增殖。"
        "此外，他们在抗体Fc区进行定点偶联，改善了GalNAc-LYTAC在体内的药代动力学表现。"
        "这项研究说明，受体选择本身就是药物设计的一部分，"
        "而不仅仅是靶蛋白结合问题。",
        "Li等（2025）把这类工作归纳为第二代组织特异性LYTAC，"
        "认为肝靶向设计对减少脱靶效应有一定帮助。",
    ],
    "case3": [
        "随着研究深入，人们发现同一种LYTAC在不同细胞中的降解效率可能差别很大。"
        "这说明除了分子结构，细胞内部环境同样重要。"
        "2023年，Ahn等在Science（382卷, eadf6249）发表了一项全基因组CRISPR筛选研究，"
        "系统寻找影响LYTAC降解效率的基因。",
        "这项研究的重要性在于，它并不只是证明LYTAC可以工作，"
        "而是开始解释为什么有些细胞降解效果好、有些效果差。",
        "研究发现，retromer复合体可以把LYTAC-CI-M6PR复合物从内体回收至细胞膜，"
        "相当于削弱了降解过程。敲除相关基因后，EGFR降解效率明显提高，"
        "有些条件下可达到90%以上。",
        "另一个重要发现是CUL3的neddylation修饰。"
        "作者比较了11种细胞系，发现neddylated CUL3水平较高的细胞，EGFR降解效果更好。"
        "这意味着细胞本身的蛋白修饰状态，可能会影响LYTAC是否好用。",
        "研究还指出，细胞内源性M6P修饰蛋白会占用CI-M6PR受体，"
        "从而和LYTAC竞争结合位点。如果抑制M6P生物合成，游离受体增多，LYTAC效率也会提高。",
        "这篇论文让我印象较深的一点是，它把药物效率和基础细胞生物学直接联系了起来。"
        "以后如果要做LYTAC药物，可能不能只盯着分子本身，还要考虑不同组织和细胞类型的差异。",
        "Sun等（2024）在综述中也引用了这项研究，"
        "认为细胞通路调控可能是提高LYTAC效率的重要方向。",
    ],
    "case4": [
        "化学合成LYTAC在制备和均一性方面有一定难度，"
        "因此也有研究尝试用基因编码的蛋白分子来实现类似功能。"
        "2023年，Pance等在Nature Biotechnology（41卷273-281页）"
        "提出了KineTAC（cytokine receptor-targeting chimera）平台。",
        "KineTAC是一种双特异性抗体，一条臂结合细胞因子受体，另一条臂结合靶蛋白。"
        "例如CXCL12可以结合CXCR7并触发受体内吞，"
        "研究者就把这一天然过程用来带动靶蛋白进入溶酶体。"
        "他们针对PD-L1、HER2、PD-1、EGFR、CDCP1和TROP2等靶点构建了KineTAC，"
        "最大降解效率大约在51%到93%之间。",
        "论文中还比较了不同细胞因子臂对降解效率的影响，"
        "说明同一个平台在不同靶点上表现并不相同。",
        "与化学LYTAC相比，KineTAC可以通过哺乳动物细胞表达生产，"
        "更换靶蛋白结合臂也相对方便，不需要复杂的聚糖-蛋白化学偶联。"
        "通过替换不同细胞因子-受体组合，还可以尝试不同的组织靶向方式。",
        "这项研究说明，LYTAC并不只有一种实现路径，"
        "蛋白质工程和免疫学方法也能参与胞外蛋白降解药物的设计。",
        "Chen等（2023）把KineTAC与化学合成的LYTAC并列讨论，"
        "认为这类蛋白平台在制备和替换靶点结合臂方面可能更方便。",
    ],
    "body4": [
        "从已发表研究来看，LYTAC已经在肿瘤、免疫相关疾病和神经退行性疾病等方向"
        "开展了较多探索，但多数仍停留在细胞和动物实验阶段。",
        "在肿瘤免疫方向，Banik等最早证明LYTAC可以降解PD-L1。"
        "2023年Li等在Journal of the American Chemical Society"
        "（145卷24506-24521页）进一步报道了基于DNA适体的共价LYTAC，"
        "可在体内降解PD-L1，并观察到比传统免疫检查点抗体更强的抗肿瘤免疫反应。"
        "作者在动物模型中比较了共价LYTAC和常规抗体阻断疗法，"
        "认为前者在激活抗肿瘤免疫方面表现更好。",
        "Pance等的KineTAC研究也证明，PD-L1、HER2等免疫相关靶点可以通过溶酶体途径被降解。",
        "在肝靶向方向，Wu等在Angewandte Chemie International Edition"
        "（2023, 62卷e202218106）开发了Apt-LYTAC，"
        "把适体与tri-GalNAc连接，实现了肝细胞中对PDGF和PTK7等靶点的降解。"
        "Li等（2025）在综述中也介绍了这类适体LYTAC，"
        "认为其合成步骤相对简单，适合快速筛选不同靶点结合序列。",
        "这类工作说明，适体分子体积小、合成相对方便，"
        "也是构建LYTAC的一条可行路线。",
        "在神经退行性疾病方向，Liu等在Chem"
        "（2023, 9卷2016-2038页）设计了可在阿尔茨海默病病灶区激活的LYTAC前体，"
        "利用Aβ沉积区铜离子催化点击反应，在局部生成活性分子，"
        "并通过CD206受体促进Aβ聚集体进入溶酶体降解。"
        "这项研究的特点是把疾病局部病理特征和药物激活结合起来，"
        "试图降低全身给药带来的脱靶风险。",
        "作者还构建了可穿越血脑屏障的纳米载体，"
        "把点击反应前体送到脑部，再在病灶区原位生成活性降解分子。",
        "Li等（2025）认为，这类病灶激活设计对神经退行性疾病可能更有意义，"
        "因为可以避免药物在全身过早发挥作用。",
        "在免疫相关疾病方向，Chen等（2023）和Li等（2025）都提到，"
        "分泌蛋白和膜蛋白是LYTAC较有潜力的靶点类型，"
        "但现有证据主要来自细胞和动物实验。",
        "综述文献也指出，分泌蛋白和自身抗体等免疫相关靶点被认为是LYTAC较有潜力的方向，"
        "但目前仍缺少正式发表的临床试验结果。",
        "Li等（2025）和Chen等（2023）都提到，"
        "LYTAC在肿瘤、免疫疾病和神经系统疾病中已有较多临床前研究，"
        "但真正进入临床阶段的项目仍然很少。",
        "这说明该方向虽然热闹，但离实际应用还有距离。",
    ],
    "body5": [
        "虽然LYTAC研究进展较快，但要真正用于新药开发仍有不少问题。"
        "含聚糖的LYTAC结构复杂，化学合成和生物偶联产物不够均一，给质量控制带来困难。"
        "CI-M6PR表达较广，可能导致非靶组织也发生蛋白降解。"
        "ASGPR等组织特异性受体有所改善，但其他器官的靶向仍缺少成熟方案。"
        "聚糖配体和蛋白载体还可能影响体内半衰期，并引发免疫反应。"
        "第一代LYTAC往往和靶蛋白一起被降解，体内有效浓度维持时间较短。",
        "Chen等（2023）在讨论胞外蛋白降解策略时指出，"
        "连接子设计、受体占用和细胞内吞分选都会影响最终降解效果。",
        "Li等（2025）也总结了LYTAC目前面临的主要问题，"
        "包括分子制备复杂、组织选择性不足以及临床前数据仍不够充分。",
        "Zhou等（2024）进一步从可成药靶点空间角度讨论了LYTAC的发展前景，"
        "认为受体和结合剂的选择仍然是最核心的设计环节。",
        "根据Chen等和Li等发表的综述，"
        "截至本文写作时，LYTAC仍主要处于临床前研究阶段，"
        "尚未见正式发表的人体临床试验结果。",
        "Sun等（2024）在总结靶向蛋白降解技术时也认为，"
        "LYTAC仍属于较新的方向，其长期安全性和有效性还需要更多研究验证。",
        "因此，这一技术虽然值得继续跟踪，但距离真正成药还有较长距离。",
    ],
    "conclusion": [
        "通过学习相关文献，我对LYTAC技术的理解是："
        "它把靶向蛋白降解从胞内蛋白拓展到了胞外和膜蛋白，"
        "为一些传统上较难干预的靶点提供了新的思路。"
        "Banik等（Nature, 2020）建立了基本框架，"
        "Ahn等（Nature Chemical Biology, 2021）引入了肝靶向受体，"
        "Ahn等（Science, 2023）揭示了细胞因素对降解效率的影响，"
        "Pance等（Nature Biotechnology, 2023）则提供了基因编码的实现方式。",
        "这些研究让我体会到，药物研发并不只是设计一个分子，"
        "还需要理解受体分布、内吞过程和溶酶体分选等生物学背景。"
        "如果缺少这些认识，很难解释为什么同一种分子在不同细胞里效果差很多。",
        "从应用角度看，LYTAC在PD-L1降解、肝靶向蛋白清除和Aβ清除等方面已有动物实验支持，"
        "说明它确实有一定可行性。但制备复杂、组织选择性不足以及缺少临床数据，"
        "仍然是明显的短板。",
        "读完这些文献后，我对LYTAC的整体印象是："
        "概念上很有吸引力，因为它把过去较难处理的胞外靶点纳入了降解范围，"
        "但真正做成药还需要解决很多工程和生物学问题。",
        "我个人认为，这个领域后续可能还会继续围绕降解效率、组织靶向和分子类型选择等问题展开。"
        "化学LYTAC和蛋白工程平台各有优缺点，未来也许会形成互补。",
        "从课程学习角度看，阅读这些论文也让我更清楚地看到，"
        "生物医学基础研究和新药开发之间联系很紧。"
        "如果只做化学偶联而不理解受体和内吞过程，"
        "很难真正判断一个LYTAC分子有没有继续开发的价值。",
        "作为研究生阶段的文献学习，我认为LYTAC是一个值得继续关注的新方向，"
        "但它目前仍属于前沿探索阶段，还不能简单等同于即将上市的新药。",
        "如果后续继续跟踪这个领域，我会更关注正式发表的临床前和临床研究数据，"
        "以及不同受体类型在组织选择性方面的新进展。",
    ],
}

REFS = [
    "[1] Banik S M, Pedram K, Wisnovsky S, Riley N M, Bertozzi C R. Lysosome-targeting chimaeras for degradation of extracellular proteins[J]. Nature, 2020, 584(7820): 291-297.",
    "[2] Ahn G, Banik S M, Miller C L, Cornean L, Gray M A, Bertozzi C R. LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation[J]. Nature Chemical Biology, 2021, 17(9): 937-946.",
    "[3] Ahn G, Banik S M, Riley N M, Cochran R V, Bertozzi C R. Elucidating the cellular determinants of targeted membrane protein degradation by lysosome-targeting chimeras[J]. Science, 2023, 382(6668): eadf6249.",
    "[4] Pance K, Gramespacher J A, Byrnes J R, et al. Modular cytokine receptor-targeting chimeras for targeted degradation of cell surface and extracellular proteins[J]. Nature Biotechnology, 2023, 41(2): 273-281.",
    "[5] Li Y, Liu X, Yu L, Huang X, Wang X, Han D, et al. Covalent LYTAC enabled by DNA aptamers for immune checkpoint degradation therapy[J]. Journal of the American Chemical Society, 2023, 145(45): 24506-24521.",
    "[6] Liu Z, Deng Q, Qin G, Yang J, Zhang H, Ren J, Qu X. Biomarker-activated multifunctional lysosome-targeting chimeras mediated selective degradation of extracellular amyloid fibrils[J]. Chem, 2023, 9(7): 2016-2038.",
    "[7] Wu Y, Lu Y, Li L, Deng K, Zhang S, Yang C, Zhu Z. Aptamer-LYTACs for targeted degradation of extracellular and membrane proteins[J]. Angewandte Chemie International Edition, 2023, 62(15): e202218106.",
    "[8] Li Y Y, Yang Y, Zhang R S, Ge R X, Xie S B. Targeted degradation of membrane and extracellular proteins with LYTACs[J]. Acta Pharmacologica Sinica, 2025, 46: 1-7.",
    "[9] Chen X, Zhou Y, Zhao Y, Tang W. Targeted degradation of extracellular secreted and membrane proteins[J]. Trends in Pharmacological Sciences, 2023, 44(11): 762-775.",
    "[10] Sun D, Lu Y, Hu Y, et al. Targeted protein degradation: advances in drug discovery and clinical practice[J]. Signal Transduction and Targeted Therapy, 2024, 9: 308.",
    "[11] Caianiello D F, Miller C L, Ahn G, Riley N M, Bertozzi C R. Bifunctional small molecules that mediate the degradation of extracellular proteins[J]. Nature Chemical Biology, 2021, 17(8): 947-953.",
    "[12] Bekes M, Langley D R, Crews C M. PROTAC targeted protein degraders: the past is prologue[J]. Nature Reviews Drug Discovery, 2022, 21(3): 181-200.",
    "[13] Zhou Y, Zhang Y, Lazerwith S E, et al. Exploring the target space of lysosome-targeting chimeras[J]. Journal of Medicinal Chemistry, 2024, 67(5): 3654-3675.",
]


def main():
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(3.17)
    section.right_margin = Cm(3.17)

    add_heading_text(doc, "溶酶体靶向嵌合体（LYTAC）技术在药物研发中的研究进展", level=1)
    add_heading_text(doc, "——一份课程文献学习综述", level=1)

    add_heading_text(doc, "摘要", level=2)
    add_paragraph(doc, CONTENT["abstract"], first_line_indent=False)
    kp = doc.add_paragraph()
    kp.paragraph_format.first_line_indent = Cm(0)
    set_run_font(kp.add_run("关键词：溶酶体靶向嵌合体，靶向蛋白降解，溶酶体靶向受体，药物研发，生物医学"), bold=True)

    add_heading_text(doc, "前言", level=2)
    for para in CONTENT["intro"]:
        add_paragraph(doc, para)

    add_heading_text(doc, "正文", level=2)

    add_heading_text(doc, "一、LYTAC技术的基本原理与分子架构", level=3)
    for para in CONTENT["body1"]:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图1  LYTAC介导的胞外及膜蛋白溶酶体降解机制示意图。"
        "LYTAC同时结合靶蛋白与溶酶体靶向受体（LTR），经内吞-溶酶体途径实现靶蛋白降解。"
        "（引自：Banik S M, Pedram K, Wisnovsky S, et al. Lysosome-targeting chimaeras for degradation "
        "of extracellular proteins. Nature, 2020, 584(7820): 291-297.）",
    )

    add_heading_text(doc, "二、LYTAC技术的关键设计要素与生物医学考量", level=3)
    for para in CONTENT["body2"]:
        add_paragraph(doc, para)

    add_heading_text(doc, "三、代表性研究实例", level=3)

    add_heading_text(doc, "实例一：LYTAC概念的建立与模块化降解平台的奠基（Nature, 2020）", level=3)
    for para in CONTENT["case1"]:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图2  首批LYTAC分子的结构设计与CI-M6PR介导的靶蛋白降解验证。"
        "（引自：Banik S M, Pedram K, Wisnovsky S, et al. Lysosome-targeting chimaeras for degradation "
        "of extracellular proteins. Nature, 2020, 584(7820): 291-297.）",
    )

    add_heading_text(doc, "实例二：肝靶向GalNAc-LYTAC的开发与组织特异性降解（Nature Chemical Biology, 2021）", level=3)
    for para in CONTENT["case2"]:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图3  GalNAc-LYTAC的结构及ASGPR介导的肝特异性靶蛋白降解。"
        "（引自：Ahn G, Banik S M, Miller C L, et al. LYTACs that engage the asialoglycoprotein receptor "
        "for targeted protein degradation. Nature Chemical Biology, 2021, 17(9): 937-946.）",
    )

    add_heading_text(doc, "实例三：LYTAC降解效率的细胞决定因素解析（Science, 2023）", level=3)
    for para in CONTENT["case3"]:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图4  全基因组CRISPR筛选揭示的LYTAC降解调控网络。"
        "（引自：Ahn G, Banik S M, Riley N M, et al. Elucidating the cellular determinants of targeted "
        "membrane protein degradation by lysosome-targeting chimeras. Science, 2023, 382(6668): eadf6249.）",
    )

    add_heading_text(doc, "实例四：KineTAC模块化平台的建立与免疫检查点降解（Nature Biotechnology, 2023）", level=3)
    for para in CONTENT["case4"]:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图5  KineTAC平台的设计原理及对多种膜蛋白的降解效果。"
        "（引自：Pance K, Gramespacher J A, Byrnes J R, et al. Modular cytokine receptor-targeting "
        "chimeras for targeted degradation of cell surface and extracellular proteins. "
        "Nature Biotechnology, 2023, 41(2): 273-281.）",
    )

    add_heading_text(doc, "四、LYTAC技术的疾病应用领域与转化研发进展", level=3)
    for para in CONTENT["body4"]:
        add_paragraph(doc, para)

    add_heading_text(doc, "五、LYTAC药物研发面临的主要挑战", level=3)
    for para in CONTENT["body5"]:
        add_paragraph(doc, para)

    add_heading_text(doc, "总结与展望", level=2)
    for para in CONTENT["conclusion"]:
        add_paragraph(doc, para)

    add_heading_text(doc, "参考文献", level=2)
    for ref in REFS:
        add_reference(doc, ref)

    output_path = "/workspace/LYTAC技术药物研发综述.docx"
    doc.save(output_path)

    all_text = "".join(
        [CONTENT["abstract"]]
        + CONTENT["intro"]
        + CONTENT["body1"]
        + CONTENT["body2"]
        + CONTENT["case1"]
        + CONTENT["case2"]
        + CONTENT["case3"]
        + CONTENT["case4"]
        + CONTENT["body4"]
        + CONTENT["body5"]
        + CONTENT["conclusion"]
    )
    chinese_count = sum(1 for c in all_text if "\u4e00" <= c <= "\u9fff")
    print(f"Document saved to: {output_path}")
    print(f"Chinese character count: {chinese_count}")


if __name__ == "__main__":
    main()
