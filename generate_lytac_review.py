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


ABSTRACT = (
    "靶向蛋白降解是近年来药物研发领域备受关注的研究方向之一。"
    "与传统抑制或阻断靶蛋白功能不同，该类策略旨在将致病蛋白导向细胞内降解系统。"
    "溶酶体靶向嵌合体（Lysosome-Targeting Chimera, LYTAC）是其中新兴的技术类型，"
    "主要利用细胞表面溶酶体靶向受体，将胞外蛋白及膜蛋白转运至溶酶体降解。"
    "鉴于PROTAC主要作用于胞内蛋白，LYTAC为胞外靶点的干预提供了补充性技术路径。"
    "本综述在梳理相关文献的基础上，介绍LYTAC的基本作用方式，"
    "并结合Banik等（Nature, 2020）、Ahn等（Nature Chemical Biology, 2021和Science, 2023）"
    "以及Pance等（Nature Biotechnology, 2023）等研究，讨论该技术的发展过程和可能的应用方向。"
    "现有研究表明，LYTAC能否真正走向临床，不仅取决于分子设计，"
    "还取决于对内吞、受体循环和溶酶体分选等生物学过程的认识程度。"
    "目前该领域仍以基础研究为主，临床证据仍然较少。"
)

INTRO = [
    "多种疾病均与蛋白质表达异常或功能失调相关，如肿瘤、自身免疫病和神经退行性疾病。"
    "传统药物研发往往通过小分子或抗体去结合靶蛋白的活性位点或结合位点，从而抑制其功能。"
    "但对于一些没有明确药物结合口袋的蛋白，或者像自身抗体这类靶点，单纯阻断往往效果有限，"
    "因此常被称为较难成药的靶点。",
    "近二十年来，靶向蛋白降解技术发展较快。PROTAC和分子胶等策略已经能把部分胞内蛋白"
    "导向泛素-蛋白酶体系统降解，其中一些分子已进入临床试验。"
    "Bekes等（2022）在Nature Reviews Drug Discovery的综述中回顾了PROTAC的发展历程，"
    "并指出蛋白降解药物正在从概念验证走向临床验证。"
    "不过PROTAC通常要求靶蛋白具有可被配体结合的胞内区域，这样才能招募E3连接酶。"
    "而分泌至细胞外的蛋白及跨膜蛋白数量众多，在多种疾病中亦十分常见，"
    "却长期未纳入传统PROTAC的主要作用范围。",
    "从生物医学角度看，这类蛋白参与细胞通讯、免疫识别和信号转导。"
    "例如PD-L1、EGFR、IgE以及多种自身抗体相关复合物，都属于胞外或膜蛋白。"
    "它们有时表达量高，有时更新较快，单纯阻断未必能长期控制疾病进程。"
    "因此，能否把这些蛋白送入溶酶体彻底降解，是一个既有科学问题也有应用价值的问题。",
    "2020年，Bertozzi团队在Nature上提出了溶酶体靶向嵌合体（LYTAC）的概念。"
    "此后几年里，不同研究组在分子设计、受体选择、细胞机制和疾病模型等方面开展了大量研究。"
    "Li等（2025）和Chen等（2023）发表的综述文章，对LYTAC技术的发展脉络和应用前景做了较系统的介绍。"
    "Wells和Kumru（2024）在Nature Reviews Drug Discovery发表的综述则从更广的胞外靶向蛋白降解"
    "（extracellular targeted protein degradation, eTPD）角度，把LYTAC与KineTAC、AbTAC等策略放在同一框架下讨论。",
    "本综述在系统梳理相关文献的基础上，"
    "归纳LYTAC技术的作用机制、代表性研究进展及应用前景，"
    "并讨论生物医学基础认识对这类新药研发思路的影响。",
]

# 正文三节，合并原有多级小节，每节内容更充实
BODY_SECTION1 = [
    "LYTAC本质上是一种双功能分子，一端结合靶蛋白，另一端结合细胞表面的溶酶体靶向受体"
    "（Lysosome-Targeting Receptor, LTR），中间通过连接子相连。"
    "靶蛋白结合部分可以用小分子、肽、抗体或适体，受体结合部分则常用甘露糖-6-磷酸聚糖或GalNAc等配体。"
    "当LYTAC同时结合靶蛋白和受体后，会形成三元复合物，随后通过内吞进入内体，"
    "再进一步分选到溶酶体，靶蛋白最终被水解（图1）。"
    "Chen等（2023）在Trends in Pharmacological Sciences上发表的综述指出，"
    "LYTAC属于依赖溶酶体的胞外蛋白降解策略，与PROTAC依赖蛋白酶体的路径明显不同。"
    "Li等（2025）在Acta Pharmacologica Sinica的综述中，"
    "把LYTAC列为靶向蛋白降解领域较新的技术方向之一。",
    "Li等（2025）在Acta Pharmacologica Sinica的综述里进一步说明，"
    "胞外蛋白和膜蛋白在肿瘤、自身免疫病和神经退行性疾病中都很常见，"
    "因此LYTAC所填补的不仅是技术空白，而是药物研发中长期以来较难处理的一类靶点。",
    "这个过程和抗体单纯阻断靶蛋白并不相同。抗体通常只是占据结合位点，"
    "而LYTAC的目标是减少靶蛋白本身的含量。"
    "文献中一般认为，这类分子属于事件驱动型作用方式，"
    "即结合剂无需具备抑制活性，只要能有效识别并捕获靶蛋白即可。"
    "由于降解发生在溶酶体，靶蛋白的信号传导功能和蛋白相互作用都可能随之减弱，"
    "这一点在EGFR和PD-L1等靶点上已有实验支持。",
    "与单纯阻断相比，直接清除致病蛋白有望更彻底地干预疾病进程，"
    "但实际疗效仍受降解效率、靶蛋白更新动力学及组织分布等因素制约。",
    "目前已报道的LTR包括CI-M6PR（也称IGF2R）、ASGPR、整合素以及部分细胞因子受体等。"
    "Chen等（2023）在综述中归纳了多种LTR及其在胞外蛋白降解中的应用，"
    "说明受体类型正在不断增加。"
    "不同受体在组织分布和内吞特性上差别较大，"
    "因此，受体选择会直接影响LYTAC作用的位置和效率。",
    "CI-M6PR在多种细胞上都有表达，因此第一代LYTAC的作用范围较广，"
    "但也更容易带来脱靶降解的风险。ASGPR主要在肝细胞表达，"
    "因此，GalNAc-LYTAC更适合肝相关靶点。",
    "受体可视为药物进入细胞的分子入口，"
    "受体类型的选择实质上决定了药物降解作用的主要组织部位。"
    "LYTAC疗效的解读不仅依赖分子结构，"
    "还取决于受体占用、复合物回收及内体向溶酶体分选等细胞生物学过程。"
    "分子设计与细胞生物学机制的有机结合，是理解相关实验结果的重要前提。",
    "从文献来看，LYTAC的效果并不只由化学结构决定，细胞本身的生物学环境也具有重要意义。"
    "在设计分子时，需要同时考虑靶蛋白结合剂、受体配体和连接子。"
    "靶蛋白结合剂需要有足够亲和力，但不要求一定抑制靶蛋白活性。"
    "膜蛋白常用抗体片段或纳米抗体，分泌蛋白可用抗体或肽段，部分研究则使用DNA适体作为结合模块。",
    "Banik等最初的LYTAC使用了抗体或小分子作为靶向模块，"
    "表明该技术从提出开始就不是固定某一种分子形式，"
    "而是可以根据靶蛋白特点选择不同的结合方式。",
    "结合位点和亲和力的小幅变化，就可能影响三元复合物形成和内吞效率。",
    "受体配体的选择和价态同样关键。第一代LYTAC多使用CI-M6PR配体，"
    "通常采用多价甘露糖-6-磷酸聚糖以提高结合能力。"
    "Ahn等后来改用肝特异性受体ASGPR的tri-GalNAc配体，"
    "说明组织分布是决定降解位置的重要因素。"
    "连接子长度和柔性会影响两个功能模块之间的空间排布，过长或过短都可能降低降解效率。"
    "近年来也有研究通过定点偶联提高抗体-配体偶联物的均一性，以改善体内稳定性。"
    "Caianiello等（2021）在Nature Chemical Biology报道了小分子降解胞外蛋白的MoDE-A策略，"
    "通过ASGPR把靶蛋白送入溶酶体，说明除抗体-聚糖偶联物外，小分子也可介导胞外蛋白降解。",
    "Li等（2025）在总结影响降解效率的因素时，"
    "提到了配体价态、连接子长度和结合表位等分子层面的变量，"
    "并强调分子设计和细胞环境需要同时考虑，而不是单独优化某一个环节。",
    "Bekes等（2022）在回顾PROTAC发展时指出，"
    "该领域已从学术研究走向产业开发，并有PROTAC分子进入临床前和早期临床试验阶段。",
    "这一发展轨迹对LYTAC的未来路径具有参照意义，"
    "因为溶酶体途径、受体内吞和循环回收均为细胞内源性生物学过程，"
    "LYTAC药物研发本质上是借助上述过程实现靶蛋白清除。",
    "结合Chen等（2023）的综述，可将LYTAC与PROTAC进行如下对比。"
    "PROTAC主要连接靶蛋白和E3连接酶，把蛋白送进蛋白酶体降解，"
    "因此更适合有胞内结构域的靶点。"
    "LYTAC则连接靶蛋白和溶酶体靶向受体，把蛋白送进溶酶体降解，"
    "因此更适合分泌蛋白、抗体和膜蛋白。"
    "这两种技术并非相互替代关系，而是针对不同蛋白亚细胞定位的互补策略。",
    "明确上述互补关系，有助于依据靶蛋白的亚细胞定位"
    "合理选择适宜的降解技术路线。",
    "另外，Chen等（2023）提到，"
    "整合素、清道夫受体和细胞因子受体等也被用于开发组织选择性更强的降解分子，"
    "未来可能会有更多受体和结合剂组合被开发出来，"
    "但这也意味着研究者需要更系统地比较不同组合之间的优劣，"
    "而不是只追求某一个分子在单一细胞系里的最高降解率。",
    "从细胞生物学角度看，受体介导的内吞并不是单向过程。"
    "复合物进入早期内体后，一部分会被分选到溶酶体，"
    "另一部分则可能通过retromer等回收机器回到细胞膜。"
    "Ahn等（2023）在Science上的研究已经证明，"
    "这种回收会明显削弱LYTAC的降解效果。"
    "因此，理解内体分选网络对解释实验结果具有重要意义。"
    "若仅依据体外Western blot结果进行判断，可能高估分子亲和力对降解效率的贡献，"
    "实际上细胞类型、受体占用情况和内吞通路的活跃程度都会产生影响。",
    "在分子设计层面，还有一个容易被忽视的问题是价态和聚集。"
    "多价配体虽然能提高与受体的结合能力，"
    "但也可能带来非特异性聚集或过快清除。"
    "Li等（2025）在综述中指出，"
    "连接子过长可能导致两个功能模块之间距离过大，"
    "不利于三元复合物的形成，而过短又可能限制构象调整。"
    "这说明LYTAC设计需要在化学可行性和生物学功能之间反复平衡，"
    "通常难以通过单次偶联即可获得理想分子。"
    "LYTAC研究要求将化学结构设计与膜受体、内吞小泡及溶酶体等细胞过程相联系，"
    "需综合运用细胞生物学与药物化学知识方能正确解读相关对照实验。",
    "例如溶酶体抑制剂实验、受体敲除实验以及CRISPR筛选，"
    "都是用来证明降解确实走了溶酶体途径，而不是其他降解路径。",
    "这些实验设计本身也体现了生物医学基础在药物研发中的重要性。",
    "从靶点选择角度看，膜蛋白和分泌蛋白在疾病中的角色各不相同。"
    "有些蛋白主要起信号传导作用，如EGFR和HER2，"
    "降解后可能同时削弱增殖和耐药相关通路。"
    "有些蛋白则属于免疫抑制分子或致病聚集物，"
    "PD-L1和Aβ等，其清除目标更侧重于恢复免疫平衡或减轻病理负担。",
    "因此，同一个LYTAC技术平台面对不同疾病时，"
    "评价指标也不会完全一样，"
    "不宜仅以体外降解百分比作为评价其开发价值的唯一指标。",
    "Wells和Kumru（2024）在综述中指出，胞外靶向蛋白降解的核心逻辑是把靶蛋白送入溶酶体，"
    "而不是像PROTAC那样依赖泛素-蛋白酶体系统。"
    "他们按分子类型把现有策略分为偶联物、双特异性抗体和小分子等几类，"
    "并强调组织分布、内吞效率和受体占用是影响成药性的关键参数。"
    "Zhao等（2022）在Signal Transduction and Targeted Therapy发表的靶向蛋白降解综述中，"
    "也提到溶酶体途径是拓展降解靶点范围的重要方向之一，"
    "这为理解LYTAC在整个TPD领域中的位置提供了更宏观的背景。",
    "在结合模块方面，Miao等（2021）较早报道了双特异性适体嵌合体，"
    "用一条适体结合CI-M6PR、另一条结合膜蛋白MET或PTK7，"
    "实现了膜蛋白向溶酶体的转运和降解。"
    "该研究证明适体可同时承担LYTAC中靶蛋白结合与受体识别双重功能，"
    "而且适体分子量小、合成相对方便，为后续Wu等和Li等的适体LYTAC研究打下了基础。",
    "从受体选择看，Chen等（2023）归纳的LTR类型已经覆盖肝特异性、肿瘤相关和广谱表达等不同需求。"
    "CI-M6PR和ASGPR分别代表“作用面广”和“肝定向”两类典型设计，"
    "而整合素、叶酸受体等则在肿瘤选择性降解中显示出潜力。"
    "Zhou等（2024）在ChemMedChem和Nature Communications上的工作进一步说明，"
    "通过选择肿瘤组织高表达的受体，有可能把降解限制在病灶相关细胞中，"
    "这对减少正常组织脱靶降解具有重要意义。",
    "Banik等（2020）的研究揭示了CI-M6PR的生理学背景。"
    "该受体负责介导带有M6P修饰的溶酶体酶进入细胞，"
    "LYTAC正是利用了这一受体介导的货物分拣通路。"
    "因此，内源性M6P糖蛋白会和外源LYTAC竞争受体结合位点，"
    "这一点在Ahn等（2023）的研究中得到了直接验证。",
    "连接子设计方面，不同研究也给出了值得参考的经验。"
    "Ahn等（2021）发现，在抗体Fc区进行定点偶联可以改善GalNAc-LYTAC的体内行为；"
    "Hamada等（2023）则通过调整适体之间的连接长度来优化HER2降解效率。"
    "这些细节说明，LYTAC开发往往需要在靶点结合、受体识别和分子构象之间反复调试，"
    "难以建立适用于所有靶点的通用连接子设计原则。",
]

BODY_SECTION2 = [
    "2020年，Banik等在Nature（584卷7820期, 291-297页）发表了LYTAC领域的开创性工作，"
    "论文题目为Lysosome-targeting chimaeras for degradation of extracellular proteins。"
    "Li等（2025）在Acta Pharmacologica Sinica的综述中，把这项工作称为LYTAC技术发展的起点。"
    "作者将靶向模块与CI-M6PR配体通过连接子连接，构建了首批LYTAC分子。"
    "他们还利用CRISPR筛选研究了LYTAC的内吞通路，发现exocyst复合体参与CI-M6PR的膜面呈现，"
    "表明该研究不仅涉及分子构建，也在补充受体介导内吞的基础认识。"
    "实验中，LYTAC成功降解了ApoE4、EGFR、CD71和PD-L1等蛋白。"
    "以EGFR为例，处理后受体水平在数小时内明显下降，溶酶体抑制剂可以阻断这一效应，"
    "说明溶酶体途径是主要降解路径。PD-L1被降解后，肿瘤细胞的免疫逃逸能力也有所减弱。"
    "作者还测试了ApoE4的降解，说明LYTAC并不只面向肿瘤靶点，"
    "在神经退行性疾病相关蛋白上也有探索价值（图2）。",
    "该研究不仅验证了靶蛋白的可降解性，"
    "还系统阐明了内吞通路及溶酶体依赖性。"
    "这表明LYTAC研究自创立之初便将细胞生物学机制与药物化学设计紧密结合。",
    "第一代LYTAC主要依赖CI-M6PR，而该受体在多种组织中都有表达，组织选择性不够理想。"
    "2021年，Ahn等在Nature Chemical Biology（17卷937-946页）报道了结合ASGPR的GalNAc-LYTAC。"
    "ASGPR主要在肝细胞表达，此前也已被用于肝靶向核酸递送。"
    "研究者将靶蛋白结合剂与tri-GalNAc配体连接，在肝细胞中实现了EGFR降解。"
    "结果显示，与单纯使用抑制性抗体相比，GalNAc-LYTAC对EGFR下游信号的抑制更持久。"
    "作者还把较小的肽段结合剂与tri-GalNAc连接，实现了整合素降解并抑制癌细胞增殖。"
    "此外，他们在抗体Fc区进行定点偶联，改善了GalNAc-LYTAC在体内的药代动力学表现（图3）。"
    "Li等（2025）把这类工作归纳为第二代组织特异性LYTAC，认为肝靶向设计对减少脱靶效应有一定帮助。",
    "该研究将组织特异性降解从概念层面落实为可操作的分子设计，"
    "即通过更换受体配体以限制降解发生的组织部位。",
    "随着研究深入，研究发现同一种LYTAC在不同细胞中的降解效率可能差别很大。"
    "2023年，Ahn等在Science（382卷, eadf6249）发表了一项全基因组CRISPR筛选研究，"
    "系统寻找影响LYTAC降解效率的基因。"
    "研究发现，retromer复合体可以把LYTAC-CI-M6PR复合物从内体回收至细胞膜，"
    "相当于削弱了降解过程。敲除相关基因后，EGFR降解效率明显提高，有些条件下可达到90%以上。"
    "另一个重要发现是CUL3的neddylation修饰。"
    "作者比较了11种细胞系，发现neddylated CUL3水平较高的细胞，EGFR降解效果更好。"
    "研究还指出，细胞内源性M6P修饰蛋白会占用CI-M6PR受体，从而和LYTAC竞争结合位点。"
    "如果抑制M6P生物合成，游离受体增多，LYTAC效率也会提高（图4）。"
    "该研究的重要贡献在于建立了药物降解效率与基础细胞生物学之间的直接关联。"
    "即便分子设计合理，不同细胞类型的生物学差异仍可导致降解效率的显著差别。",
    "化学合成LYTAC在制备和均一性方面有一定难度，"
    "因此也有研究尝试用基因编码的蛋白分子来实现类似功能。"
    "2023年，Pance等在Nature Biotechnology（41卷273-281页）提出了KineTAC平台。"
    "KineTAC是一种双特异性抗体，一条臂结合细胞因子受体，另一条臂结合靶蛋白。"
    "例如CXCL12可以结合CXCR7并触发受体内吞，研究者利用这一天然过程促进靶蛋白进入溶酶体。"
    "他们针对PD-L1、HER2、PD-1、EGFR、CDCP1和TROP2等靶点构建了KineTAC，"
    "最大降解效率大约在51%到93%之间。"
    "与化学LYTAC相比，KineTAC可以通过哺乳动物细胞表达生产，"
    "更换靶蛋白结合臂也相对方便，不需要复杂的聚糖-蛋白化学偶联（图5）。"
    "Chen等（2023）把KineTAC与化学合成的LYTAC并列讨论，"
    "认为这类蛋白平台在制备和替换靶点结合臂方面可能更方便。",
    "Wu等（2023）在Angewandte Chemie International Edition报道了Apt-LYTAC，"
    "把适体与tri-GalNAc连接，实现了肝细胞中对PDGF和PTK7等靶点的降解。"
    "Li等（2023）在Journal of the American Chemical Society报道了基于DNA适体的共价LYTAC，"
    "可在体内降解PD-L1，并观察到比传统免疫检查点抗体更强的抗肿瘤免疫反应。",
    "Li等（2025）在综述中也介绍了这类适体LYTAC，"
    "认为其合成步骤相对简单，适合快速筛选不同靶点结合序列。",
    "这些研究共同说明，LYTAC已经从最初的概念验证，逐步发展到组织靶向、机制解析、"
    "适体平台和蛋白工程平台等多个方向。",
    "综合上述研究可见，该领域的发展并非单一线性推进，"
    "而是同时在分子类型、受体选择和细胞机制三个层面不断扩展。",
    "Banik等的Nature论文的重要意义不仅在于提出了LYTAC这一概念，"
    "更在于首次较完整地证明该溶酶体降解路线可推广至多个疾病相关靶点。"
    "Ahn等（2021）的工作实现了组织特异性降解，"
    "而2023年Science论文则进一步阐明了不同细胞间降解效率差异的生物学基础。",
    "Pance等的KineTAC则提供了另一条实现路径，"
    "说明胞外蛋白降解并不必然依赖化学合成的聚糖-蛋白偶联物。",
    "Wu等和Li等的适体研究表明结合模块可进一步小型化、更灵活。",
    "综合上述研究可见，LYTAC已从单一的偶联分子概念，"
    "拓展为一类借助溶酶体途径清除胞外蛋白的降解策略，"
    "其实现形式呈现多样化特征。"
    "上述代表性论文系统覆盖了LYTAC从概念提出、机制优化到平台拓展的主要发展阶段，"
    "且均发表于高水平学术期刊。",
    "进一步比较这几类平台，化学合成的聚糖-蛋白LYTAC在结构修饰上较灵活，"
    "可以精细调整配体价态和连接子，但制备难度也更高。"
    "KineTAC依托抗体工程，在更换靶点结合臂时相对方便，"
    "适用于需快速筛选多个靶点的研究场景。"
    "适体LYTAC则分子量较小，合成和修饰步骤相对简单，"
    "在体外筛选不同结合序列时有一定优势。"
    "不过不同平台也各有局限，"
    "化学LYTAC的均一性和放大生产仍是难点，"
    "KineTAC需要解决蛋白药物的免疫原性和稳定性问题，"
    "适体则要考虑体内核酸酶降解和靶向特异性。",
    "Banik等最初的Nature论文除了验证多个靶点，"
    "还通过CRISPR筛选找到了exocyst复合体参与CI-M6PR膜面呈现的证据，"
    "表明早期研究已关注LYTAC-受体复合物的细胞内转运途径。",
    "Ahn等2023年的Science论文则把这一问题推进到了全基因组层面，"
    "不仅鉴定了retromer介导的回收等抑制性因素，"
    "还发现了CUL3 neddylation和M6P竞争等促进性因素。"
    "这种从表型到机制的递进研究模式，"
    "体现了该领域对降解效率决定因素的深入探索。",
    "Pance等的KineTAC研究在靶点覆盖面上较广，"
    "从免疫检查点分子到生长因子受体都有测试，"
    "说明溶酶体降解路线并不局限于某一个疾病领域。"
    "Wu等和Li等的适体研究则把结合模块从抗体缩小到核酸适体，"
    "拓展了LYTAC可用结合模块的类型。"
    "Li等（2023）的共价适体LYTAC还在体内实验中展示了抗肿瘤免疫增强，"
    "这是少数把LYTAC和动物药效直接联系起来的工作之一。",
    "Liu等（2023）在Chem上发表的病灶激活LYTAC前体，"
    "则代表了另一条设计思路："
    "并非单纯提高全身给药强度，"
    "而是利用疾病局部特有的化学环境来触发药物活化。"
    "这类研究与阿尔茨海默病中Aβ沉积的病理特点结合较紧，"
    "说明疾病生物学本身也可以成为药物设计的重要输入。",
    "综合来看，"
    "LYTAC领域正从可行性验证阶段转向精准性、稳定性及可开发性优化阶段。",
    "不同研究组的工作重点亦存在明显差异。"
    "Bertozzi团队更强调受体化学、细胞机制和分子平台的系统搭建，"
    "国内和亚洲其他实验室则在适体偶联、病灶激活和疾病模型验证方面开展了系列拓展研究。"
    "这种分工说明LYTAC已经从一个实验室的概念，逐渐变成多个团队共同参与的研究方向。",
    "在适体方向，研究时间线也相对清晰。"
    "Miao等（2021）用双特异性适体嵌合体降解MET和PTK7，"
    "证明核酸适体可以替代抗体作为LYTAC的结合模块。"
    "Hamada等（2023）在Cell Reports Physical Science进一步构建了HER2靶向的适体LYTAC，"
    "用分别识别HER2和IGF2R（CI-M6PR）的两条DNA适体连接成嵌合体，"
    "在乳腺癌细胞中实现了HER2降解并抑制细胞增殖。"
    "溶酶体抑制剂和内吞阻断实验表明，这一降解过程依赖溶酶体途径，"
    "说明适体LYTAC的作用机制与Banik等提出的经典路线一致。",
    "Wu等（2023）的Apt-LYTAC则把适体与tri-GalNAc偶联，"
    "在保留肝靶向特性的同时简化了分子制备。"
    "Li等（2023）的共价适体LYTAC通过生物正交反应增强适体与PD-L1的结合，"
    "在复杂体内环境中提高了降解效率，并在肿瘤模型中显示出免疫激活效应。",
    "Li等（2025）在综述中把上述工作归纳为适体类LYTAC的重要进展，"
    "认为这类分子在快速筛选靶点和调节连接子长度方面具有独特优势。",
    "在组织选择性方面，Tang团队的工作具有代表性。"
    "Zhou等（2024）在ChemMedChem报道了整合素靶向嵌合体（ITAC），"
    "把cRGD肽与抗体偶联，利用肿瘤高表达的整合素介导EGFR等靶点降解。"
    "实验显示，ITAC在多种癌细胞中降解膜蛋白的效率高于正常角质形成细胞，"
    "表明癌选择性降解可通过受体类型的合理选择加以实现。"
    "同年，Zhou等又在Nature Communications报道了叶酸受体靶向嵌合体（FRTAC），"
    "用叶酸配体连接抗肿瘤抗体，在叶酸受体高表达的肿瘤细胞中降解胞外和膜蛋白。"
    "动物实验表明，FRTAC在多种同系小鼠肿瘤模型中比单纯阻断抗体更能抑制肿瘤生长，"
    "而且制备上只需利用商业化的叶酸偶联试剂和抗体，制备可及性较高。",
    "Chen等（2023）和Wells等（2024）都把ITAC、FRTAC这类策略归入组织选择性LYTAC的延伸，"
    "说明该领域正在从“能不能降解”走向“在哪个细胞里降解”。",
    "如果把时间线拉长来看，"
    "2010年代中后期PROTAC的快速发展为蛋白降解药物提供了重要参照。"
    "Bekes等（2022）回顾这一历程时指出，"
    "PROTAC在取得临床概念验证后，产业界投入明显增加。"
    "LYTAC目前处于类似但更早期的阶段："
    "已有高质量概念验证论文，也有多个并行平台，但还缺少系统的临床数据。",
    "从研究方法上看，这个领域对遗传学和蛋白质组学工具的依赖越来越强。"
    "Banik等用CRISPRi筛选内吞相关基因，"
    "Ahn等用全基因组CRISPR敲除筛选调节因子，"
    "都说明LYTAC研究已不仅限于“做一个偶联分子看能不能降解”，"
    "而是在追问哪些细胞过程决定了降解效率。",
    "上述研究模式将药物化学与细胞生物学置于同等重要的地位，"
    "体现了机制导向药物研发的基本理念。",
]

BODY_SECTION3 = [
    "从已发表研究来看，LYTAC已经在肿瘤、免疫相关疾病和神经退行性疾病等方向"
    "开展了较多探索，但多数仍停留在细胞和动物实验阶段。"
    "在肿瘤免疫方向，Banik等最早证明LYTAC可以降解PD-L1，"
    "Li等进一步用共价适体LYTAC在动物模型中验证了抗肿瘤效果，"
    "Pance等的KineTAC研究也证明PD-L1、HER2等靶点可以通过溶酶体途径被降解。",
    "PD-L1是免疫检查点治疗中非常重要的靶点，目前临床主要使用抗体阻断PD-1/PD-L1相互作用。"
    "Li等的工作说明，若将PD-L1直接降解，"
    "可能在某些肿瘤模型中产生比单纯阻断更强的免疫激活效果。"
    "不过这类结果目前仍主要来自动物实验，尚不能直接外推至人体疗效。",
    "在肝靶向方向，Ahn等的GalNAc-LYTAC和Wu等的Apt-LYTAC都说明，"
    "选择合适的肝特异性受体是减少全身脱靶的重要思路。",
    "肝相关疾病中，很多靶蛋白需要在肝脏局部被清除，"
    "如果药物分子能在全身广泛降解靶蛋白，就可能带来不必要的副作用。",
    "因此，组织特异性受体在LYTAC研发中具有关键意义，"
    "是决定候选药物是否具备开发价值的重要因素之一。",
    "在神经退行性疾病方向，Liu等（2023）在Chem设计了可在阿尔茨海默病病灶区激活的LYTAC前体，"
    "利用Aβ沉积区铜离子催化点击反应，在局部生成活性分子，"
    "并通过CD206受体促进Aβ聚集体进入溶酶体降解。"
    "这项研究的特点是把疾病局部病理特征和药物激活结合起来，"
    "试图降低全身给药带来的脱靶风险。"
    "作者还构建了可穿越血脑屏障的纳米载体，把点击反应前体送到脑部，"
    "再在病灶区原位生成活性降解分子。"
    "Li等（2025）认为，这类病灶激活设计对神经退行性疾病可能更有意义，"
    "因为可以避免药物在全身过早发挥作用。",
    "阿尔茨海默病药物研发长期面临血脑屏障和靶点复杂性的问题，"
    "病灶激活策略表明，"
    "药物设计需充分考虑疾病局部病理微环境，而非单纯追求结合亲和力的提高。",
    "Chen等（2023）和Li等（2025）都提到，分泌蛋白和膜蛋白是LYTAC较有潜力的靶点类型，"
    "但现有证据主要来自细胞和动物实验，真正进入临床阶段的项目仍然很少。"
    "这说明该方向虽然研究较多，但距离实际临床应用仍有较大差距。",
    "Li等（2025）在综述中也指出，"
    "LYTAC目前仍主要处于临床前阶段，其长期安全性和有效性还需要更多研究验证。",
    "虽然LYTAC研究进展较快，但要真正用于新药开发仍面临诸多挑战。"
    "含聚糖的LYTAC结构复杂，化学合成和生物偶联产物不够均一，给质量控制带来困难。",
    "从产业化角度看，该问题的重要性不亚于体外活性本身，"
    "因为候选药物除有效性外，还需满足稳定性、批次间一致性和可放大生产等要求。",
    "CI-M6PR表达较广，可能导致非靶组织也发生蛋白降解。"
    "ASGPR等组织特异性受体有所改善，但其他器官的靶向仍缺少成熟方案。"
    "聚糖配体和蛋白载体还可能影响体内半衰期，并引发免疫反应。"
    "第一代LYTAC往往和靶蛋白一起被降解，体内有效浓度维持时间较短。",
    "这意味着一个分子可能只能处理有限数量的靶蛋白，"
    "如果靶蛋白浓度很高或者更新很快，药效就可能不够持久。",
    "Chen等（2023）在讨论胞外蛋白降解策略时指出，"
    "连接子设计、受体占用和细胞内吞分选都会影响最终降解效果。"
    "Li等（2025）也总结了LYTAC目前面临的主要问题，"
    "包括分子制备复杂、组织选择性不足以及临床前数据仍不够充分。"
    "Chen等（2023）从受体和结合剂选择角度讨论了胞外蛋白降解策略的设计要点，"
    "认为受体类型和靶点结合方式仍是核心环节。",
    "根据Chen等和Li等发表的综述，截至相关文献发表时，"
    "LYTAC仍主要处于临床前研究阶段，尚未见正式发表的人体临床试验结果。"
    "因此，这一技术虽然值得持续关注，但距离真正成药还有较长距离，"
    "在投入大量开发资源之前还需要更充分的安全性和药效学证据。",
    "这一点与PROTAC等降解技术早期的发展路径有些相似，"
    "都需要先完成较扎实的机制研究和临床前验证工作。",
    "上述分析表明，LYTAC的疗效并非仅由分子偶联结构决定，"
    "受体组织分布、复合物回收动力学及靶蛋白亚细胞定位等生物医学因素均具有重要影响。"
    "后续研究需重点关注临床前与临床数据的系统积累，"
    "以及新型组织特异性受体的鉴定与验证。"
    "LYTAC研发进一步表明，"
    "诸多表面属化学或药理学范畴的问题，"
    "最终均需回归细胞内吞、分选及溶酶体降解等基础生物学过程加以阐释。",
    "在安全性方面，文献讨论较多的问题包括脱靶降解和免疫原性。"
    "由于CI-M6PR等受体在多种组织都有表达，"
    "第一代LYTAC如果全身给药，理论上可能影响正常组织的蛋白稳态。"
    "ASGPR介导的肝靶向设计在一定程度上缓解了这一问题，"
    "但肺、肾、脑等其他器官的特异性递送仍缺少成熟方案。"
    "含聚糖和蛋白的偶联物还可能被免疫系统识别，"
    "长期反复给药时的免疫反应需要在动物实验中认真评估。",
    "在药效学方面，LYTAC属于催化性较低、消耗性较强的模式。"
    "一个分子把靶蛋白和受体一起带入溶酶体后，"
    "往往自身也会被降解，因此体内有效浓度维持时间可能较短。"
    "如果靶蛋白更新速度快或者表达量很高，"
    "就需要更频繁给药或设计更稳定的分子形式。"
    "这与传统抗体药物可以较长时间占据靶点的方式有所不同，"
    "也给药代动力学和给药方案设计带来了新的考虑。",
    "从药物开发流程看，LYTAC目前仍缺少标准化的大分子偶联工艺和质量控制方法。"
    "化学合成聚糖配体、定点偶联蛋白以及产物表征都需要较高技术门槛，"
    "不同批次之间若均一性不足，就可能影响体内结果的可重复性。"
    "相比之下，KineTAC等基因编码平台在制备上更接近生物药生产路线，"
    "但同样要面对表达、纯化和稳定性等问题。",
    "Bekes等（2022）回顾PROTAC时指出，"
    "随着临床概念验证的取得，该领域正在探索更多靶点类型和E3连接酶种类，"
    "并逐步向肿瘤以外的疾病领域扩展。",
    "LYTAC目前处于类似的前临床数据积累阶段，"
    "已有不少概念验证和动物实验，但系统性的开发数据仍然偏少。",
    "如果未来要推进到临床试验，"
    "可能需要先明确哪些疾病、哪些靶点和哪些受体组合最值得投入，"
    "而非在所有潜在靶点上平均分配研发资源。",
    "总体而言，该领域科学问题新颖、应用前景广阔，"
    "但成药前仍存在较多未明之处。"
    "PD-L1、EGFR、Aβ等靶点已有较充分的细胞和动物数据，"
    "说明技术路线基本可行，"
    "但人体中的受体表达、免疫环境和给药途径都可能与动物模型不同，"
    "因此现有结果只能作为参考，不能简单等同于临床疗效。",
    "后续值得关注的方向包括："
    "新型受体类型的鉴定与应用、"
    "简化制备工艺的开发，"
    "以及更系统的毒理学和药代动力学研究。",
    "另外，Caianiello等（2021）报道的MoDE-A策略说明，"
    "小分子也可能介导胞外蛋白降解，"
    "这为未来开发更易制备的降解药物提供了另一种可能。"
    "不过MoDE-A和经典LYTAC在分子大小、靶点类型和体内行为上仍有差别，"
    "能否成为主流路线还需要更多研究验证。",
    "总体而言，LYTAC把生物医学中早已存在的受体介导内吞和溶酶体降解过程，"
    "重新包装成可设计的药物工具，"
    "该技术兼具方法学创新性与对基础学科知识积累的深度依赖。",
    "在肿瘤选择性降解方向，Zhou等（2024）的FRTAC研究提供了较完整的体内证据。"
    "叶酸受体在多种恶性肿瘤细胞表面高表达，而在多数正常组织中表达较低，"
    "以叶酸作为受体配体，有望把降解主要限制在肿瘤细胞。"
    "论文中报道，FRTAC可在体外和体内降解肿瘤相关蛋白，"
    "并在三种同系小鼠肿瘤模型中显示出比阻断抗体更强的抑瘤效果。"
    "该类研究将精准医学理念与蛋白降解技术相结合，"
    "实现了在病灶相关细胞中的选择性靶蛋白清除。",
    "ITAC策略则从另一个角度实现肿瘤选择性。"
    "Zhou等（2024）把cRGD肽连接到抗体上，利用整合素在癌细胞表面的高表达，"
    "促进EGFR等膜蛋白的内吞和溶酶体降解。"
    "实验比较了癌细胞系和正常角质形成细胞，"
    "发现ITAC在癌细胞中的降解效率更高。"
    "这说明受体表达差异可以成为药物设计的重要依据，"
    "也呼应了Ahn等通过换用ASGPR实现肝靶向的思路。",
    "在适体平台方面，除Wu等和Li等的研究外，"
    "Miao等（2021）和Hamada等（2023）的工作说明适体LYTAC已经覆盖MET、PTK7和HER2等多个靶点。"
    "适体的优势是分子小、修饰灵活，"
    "其局限性在于体内可能面临核酸酶降解和药代动力学方面的挑战。"
    "因此，如何将适体LYTAC的稳定性和靶向性做到可开发水平，仍是后续研究的重点。",
    "Wells和Kumru（2024）在讨论eTPD临床转化时指出，"
    "胞外降解药物需要同时解决分子设计、组织分布、免疫原性和制备可放大性等问题。"
    "这些判断与LYTAC目前面临的困境基本一致。"
    "Zhao等（2022）在回顾靶向蛋白降解技术时也强调，"
    "任何新平台都要经过较充分的机制验证、药代毒理评估和制剂研究，"
    "才有可能进入临床开发阶段。",
    "近年来LYTAC领域最显著的变化在于技术平台的多样化。",
    "除了Banik等提出的聚糖-抗体偶联路线，"
    "现在已经出现了KineTAC、Apt-LYTAC、MoDE-A、ITAC、FRTAC和病灶激活前体等多种形式。"
    "不同平台各有适用场景："
    "肝靶向疾病可能更适合GalNAc或Apt-LYTAC，"
    "肿瘤免疫可能关注PD-L1和HER2降解，"
    "神经退行性疾病则更强调血脑屏障穿透和局部激活。",
    "因此，评价一种LYTAC有没有前景，"
    "不能只看体外降解率，还要结合靶点、受体、给药途径和疾病模型来综合判断。",
    "展望未来，该领域的发展可能沿以下方向推进："
    "组织特异性受体的开发与验证、"
    "制备工艺简化与产物均一性提升，"
    "以及体内药效学、毒理学和药代动力学的系统研究。"
    "只有把这些环节的数据逐步补齐，"
    "LYTAC才有可能从实验室概念真正走向候选药物。",
    "在自身免疫和炎症领域，胞外细胞因子和膜受体也是重要靶点，"
    "但公开文献中针对这类靶点的LYTAC研究仍相对有限。"
    "Chen等（2023）在综述中提到，细胞因子受体、生长因子受体和免疫检查点分子"
    "都属于胞外降解值得关注的对象，"
    "未来是否能开发出组织选择性更好的降解剂，还需要更多实验验证。",
    "工艺与质量控制是LYTAC成药前必须解决的关键问题。",
    "抗体-聚糖偶联物往往存在偶联位点和偶联数不一致的问题，"
    "不同批次分子的降解活性可能出现波动。"
    "Pance等提出的KineTAC和Miao等、Hamada等采用的适体路线，"
    "在一定程度上是为了绕开复杂化学偶联带来的制备困难，"
    "但每种替代方案也会带来新的稳定性或免疫原性风险。",
    "因此，未来哪种平台最适合产业化，"
    "可能不取决于哪篇论文的体外数据最高，"
    "而取决于能否建立稳定、可重复、可放大的生产工艺。",
    "该领域的研究论文通常同时涵盖化学结构、细胞实验与动物模型，"
    "对多学科知识背景具有较高要求。"
    "例如，解读Ahn等（2023）的Science论文需熟悉retromer复合体、"
    "M6P生物合成及CUL3 neddylation等细胞生物学概念；"
    "理解Liu等（2023）的Chem论文则需掌握Aβ沉积病理、"
    "铜离子催化点击反应及血脑屏障递送等相关知识。"
    "上述特点表明，新药研发是生物医学、药学与化学等多学科交叉融合的过程。"
    "LYTAC作为一项仍在快速演进的技术，"
    "其核心策略是借助溶酶体途径清除胞外及膜蛋白；"
    "其成功开发有赖于受体生物学、分子设计与疾病模型三个层面的协同优化。"
    "现有文献已提供较为丰富的概念验证数据，"
    "但距临床药物转化仍有明显差距。"
    "后续研究进展可能集中于新型受体类型、简化制备路线"
    "及更完整的临床前数据等方面。",
    "总体而言，LYTAC代表了靶向蛋白降解向胞外靶点拓展的重要尝试，"
    "其发展过程充分体现了基础生物医学研究对新药研发的推动作用。",
]

CONCLUSION = [
    "综上所述，LYTAC技术将靶向蛋白降解的应用范围从胞内蛋白拓展至胞外及膜蛋白，"
    "为传统上难以干预的靶点提供了新的药物研发策略。",
    "Banik等（Nature, 2020）建立了LYTAC的基本技术框架，"
    "Ahn等（Nature Chemical Biology, 2021）实现了肝组织特异性降解，"
    "Ahn等（Science, 2023）系统阐明了影响降解效率的细胞因素，"
    "Pance等（Nature Biotechnology, 2023）则提供了基因编码的蛋白降解平台。",
    "上述研究表明，成功的药物研发不仅依赖分子设计，"
    "还需深入理解受体分布、内吞动力学及溶酶体分选等生物学过程。",
    "在应用层面，LYTAC在PD-L1降解、肝靶向蛋白清除及Aβ清除等方面已获得动物实验支持，"
    "表明该技术具有一定可行性。"
    "然而，分子制备复杂、组织选择性不足及临床数据匮乏仍是制约其发展的主要瓶颈。",
    "LYTAC目前仍处于前沿探索阶段，"
    "尚不能等同于即将进入临床的成熟药物技术。"
    "随着机制研究的深入及新型平台的涌现，"
    "该领域有望在未来取得更多突破性进展。",
]

REFS = [
    "[1] Banik S M, Pedram K, Wisnovsky S, Riley N M, Bertozzi C R. Lysosome-targeting chimaeras for degradation of extracellular proteins[J]. Nature, 2020, 584(7820): 291-297. https://doi.org/10.1038/s41586-020-2545-9",
    "[2] Ahn G, Banik S M, Miller C L, Riley N M, Cochran J R, Bertozzi C R. LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation[J]. Nature Chemical Biology, 2021, 17(9): 937-946. https://doi.org/10.1038/s41589-021-00770-1",
    "[3] Ahn G, Riley N M, Kamber R A, Wisnovsky S, Moncayo von Hase S, Bassik M C, Banik S M, Bertozzi C R. Elucidating the cellular determinants of targeted membrane protein degradation by lysosome-targeting chimeras[J]. Science, 2023, 382(6668): eadf6249. https://doi.org/10.1126/science.adf6249",
    "[4] Pance K, Gramespacher J A, Byrnes J R, et al. Modular cytokine receptor-targeting chimeras for targeted degradation of cell surface and extracellular proteins[J]. Nature Biotechnology, 2023, 41(2): 273-281. https://doi.org/10.1038/s41587-022-01456-2",
    "[5] Li Y, Liu X, Yu L, Huang X, Wang X, Han D, et al. Covalent LYTAC enabled by DNA aptamers for immune checkpoint degradation therapy[J]. Journal of the American Chemical Society, 2023, 145(45): 24506-24521. https://doi.org/10.1021/jacs.3c03899",
    "[6] Liu Z, Deng Q, Qin G, Yang J, Zhang H, Ren J, Qu X. Biomarker-activated multifunctional lysosome-targeting chimeras mediated selective degradation of extracellular amyloid fibrils[J]. Chem, 2023, 9(7): 2016-2038. https://doi.org/10.1016/j.chempr.2023.06.003",
    "[7] Wu Y, Lu Y, Li L, Deng K, Zhang S, Yang C, Zhu Z. Aptamer-LYTACs for targeted degradation of extracellular and membrane proteins[J]. Angewandte Chemie International Edition, 2023, 62(15): e202218106. https://doi.org/10.1002/anie.202218106",
    "[8] Li Y Y, Yang Y, Zhang R S, Ge R X, Xie S B. Targeted degradation of membrane and extracellular proteins with LYTACs[J]. Acta Pharmacologica Sinica, 2025, 46: 1-7. https://doi.org/10.1038/s41401-024-01364-y",
    "[9] Chen X, Zhou Y, Zhao Y, Tang W. Targeted degradation of extracellular secreted and membrane proteins[J]. Trends in Pharmacological Sciences, 2023, 44(11): 762-775. https://doi.org/10.1016/j.tips.2023.08.013",
    "[10] Caianiello D F, Zhang M, Ray J D, et al. Bifunctional small molecules that mediate the degradation of extracellular proteins[J]. Nature Chemical Biology, 2021, 17(9): 947-953. https://doi.org/10.1038/s41589-021-00851-1",
    "[11] Bekes M, Langley D R, Crews C M. PROTAC targeted protein degraders: the past is prologue[J]. Nature Reviews Drug Discovery, 2022, 21(3): 181-200. https://doi.org/10.1038/s41573-021-00371-6",
    "[12] Wells J A, Kumru K. Extracellular targeted protein degradation: an emerging modality for drug discovery[J]. Nature Reviews Drug Discovery, 2024, 23(2): 126-140. https://doi.org/10.1038/s41573-023-00833-z",
    "[13] Miao Y, Gao Q, Mao M, Zhang C, Yang L, Yang Y, Han D. Bispecific aptamer chimeras enable targeted protein degradation on cell membranes[J]. Angewandte Chemie International Edition, 2021, 60(20): 11267-11271. https://doi.org/10.1002/anie.202102170",
    "[14] Zhou Y, Liao Y, Zhao Y, Tang W. Development of integrin targeting chimeras (ITACs) for the lysosomal degradation of extracellular proteins[J]. ChemMedChem, 2024, 19(24): e202300643. https://doi.org/10.1002/cmdc.202300643",
    "[15] Zhou Y, Li C, Chen X, Zhao Y, Liao Y, Huang P, et al. Development of folate receptor targeting chimeras for cancer selective degradation of extracellular proteins[J]. Nature Communications, 2024, 15: 8695. https://doi.org/10.1038/s41467-024-52685-9",
    "[16] Hamada K, Hashimoto T, Iwashita R, Yamada Y, Kikkawa Y, Nomizu M. Development of a bispecific DNA-aptamer-based lysosome-targeting chimera for HER2 protein degradation[J]. Cell Reports Physical Science, 2023, 4(3): 101296. https://doi.org/10.1016/j.xcrp.2023.101296",
    "[17] Zhao L, Zhao J, Zhong K, Tong A, Jia D. Targeted protein degradation: mechanisms, strategies and application[J]. Signal Transduction and Targeted Therapy, 2022, 7: 113. https://doi.org/10.1038/s41392-022-00966-4",
]


def main():
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(3.17)
    section.right_margin = Cm(3.17)

    add_heading_text(doc, "溶酶体靶向嵌合体（LYTAC）技术在药物研发中的研究进展", level=1)
    add_heading_text(doc, "——文献综述", level=1)

    add_heading_text(doc, "摘要", level=2)
    add_paragraph(doc, ABSTRACT, first_line_indent=False)
    kp = doc.add_paragraph()
    kp.paragraph_format.first_line_indent = Cm(0)
    set_run_font(kp.add_run("关键词：溶酶体靶向嵌合体，靶向蛋白降解，溶酶体靶向受体，药物研发，生物医学"), bold=True)

    add_heading_text(doc, "前言", level=2)
    for para in INTRO:
        add_paragraph(doc, para)

    add_heading_text(doc, "正文", level=2)

    add_heading_text(doc, "一、LYTAC的作用机制与分子设计基础", level=3)
    for para in BODY_SECTION1:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图1  LYTAC介导的胞外及膜蛋白溶酶体降解机制示意图。"
        "（引自：Banik S M, Pedram K, Wisnovsky S, et al. Nature, 2020, 584(7820): 291-297.）",
    )

    add_heading_text(doc, "二、代表性研究进展", level=3)
    for para in BODY_SECTION2:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图2  首批LYTAC分子的结构设计与靶蛋白降解验证（Banik et al., Nature, 2020）。"
        "图3  GalNAc-LYTAC的肝靶向降解（Ahn et al., Nat Chem Biol, 2021）。"
        "图4  CRISPR筛选揭示的LYTAC降解调控网络（Ahn et al., Science, 2023）。"
        "图5  KineTAC平台原理（Pance et al., Nat Biotechnol, 2023）。",
    )

    add_heading_text(doc, "三、应用探索与面临的主要问题", level=3)
    for para in BODY_SECTION3:
        add_paragraph(doc, para)

    add_heading_text(doc, "总结与展望", level=2)
    for para in CONCLUSION:
        add_paragraph(doc, para)

    add_heading_text(doc, "参考文献", level=2)
    for ref in REFS:
        add_reference(doc, ref)

    output_path = "/workspace/LYTAC技术药物研发综述.docx"
    doc.save(output_path)

    body_text = "".join(BODY_SECTION1 + BODY_SECTION2 + BODY_SECTION3)
    body_count = sum(1 for c in body_text if "\u4e00" <= c <= "\u9fff")
    total_text = ABSTRACT + "".join(INTRO + BODY_SECTION1 + BODY_SECTION2 + BODY_SECTION3 + CONCLUSION)
    total_count = sum(1 for c in total_text if "\u4e00" <= c <= "\u9fff")
    print(f"Document saved to: {output_path}")
    print(f"正文汉字数: {body_count}")
    print(f"全文汉字数: {total_count}")


if __name__ == "__main__":
    main()
