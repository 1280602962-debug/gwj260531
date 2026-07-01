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
        "靶向蛋白降解（Targeted Protein Degradation, TPD）已成为当代创新药物研发的重要方向。"
        "与传统占据性抑制策略不同，TPD通过事件驱动机制实现致病蛋白的选择性、彻底性清除，"
        "在克服「不可成药」靶点方面展现出独特优势。溶酶体靶向嵌合体（Lysosome-Targeting Chimera, LYTAC）"
        "作为TPD技术家族的重要成员，利用细胞表面溶酶体靶向受体（Lysosome-Targeting Receptor, LTR）"
        "介导的内吞-溶酶体途径，将降解范围从胞内蛋白拓展至约占人类蛋白质组40%的胞外及膜蛋白，"
        "填补了PROTAC等泛素-蛋白酶体降解技术无法覆盖的生物学空间。本文从生物医学与药物研发的交叉视角出发，"
        "系统梳理LYTAC技术的基本原理、分子设计策略及关键生物学决定因素，"
        "重点剖析Nature、Science、Nature Chemical Biology、Nature Biotechnology等高水平期刊发表的"
        "代表性研究成果，涵盖技术奠基、组织特异性靶向、细胞降解机制解析及模块化平台构建等维度。"
        "同时，本文结合自身免疫病、肿瘤免疫及神经退行性疾病等领域的转化进展，"
        "讨论LYTAC从实验室概念验证走向临床候选药物过程中所面临的配体设计、组织选择性、"
        "药代动力学及免疫原性等挑战。研究表明，深入理解受体占用、内吞分选、"
        "溶酶体成熟等细胞生物学过程，是LYTAC药物理性设计与疗效预测的核心基础；"
        "生物医学研究不仅为LYTAC提供了机制框架，更直接驱动了cataLYTAC等下一代降解平台的诞生。"
        "展望未来，LYTAC技术有望在过敏性疾病、自身免疫病及肿瘤等领域实现突破，"
        "成为连接基础生命科学与临床药物开发的重要桥梁。"
    ),
    "intro": [
        "蛋白质作为生命活动的主要执行者，其表达失调、错误折叠或异常聚集与肿瘤、自身免疫病、"
        "神经退行性疾病及代谢性疾病等多种重大疾病密切相关。现代药物研发长期以「占据性抑制」"
        "（occupancy-based inhibition）为主要逻辑，即通过小分子或抗体与靶蛋白活性位点或配体结合位点结合，"
        "阻断其生物学功能。然而，对于缺乏明确催化活性口袋的支架蛋白、构象多变的多聚体蛋白、"
        "以及致病性自身抗体等靶点，传统抑制策略往往难以实现充分且持久的药效，"
        "这些靶点因而被统称为「不可成药」（undruggable）靶点。",
        "21世纪以来，以泛素-蛋白酶体系统（Ubiquitin-Proteasome System, UPS）为劫持对象的"
        "蛋白降解技术迅速崛起。蛋白水解靶向嵌合体（Proteolysis-Targeting Chimera, PROTAC）"
        "和分子胶（Molecular Glue）等策略已在血液肿瘤等领域进入临床应用或后期临床试验，"
        "证明了「降解优于抑制」（degradation over inhibition）的治疗理念具有坚实的转化价值。"
        "然而，PROTAC的作用机制要求靶蛋白必须具有可被配体结合的胞内结构域，"
        "以便招募E3泛素连接酶并启动蛋白酶体降解，这使得约占蛋白质组40%的胞外分泌蛋白和跨膜蛋白"
        "长期处于TPD技术的「盲区」。",
        "从生物医学角度看，胞外及膜蛋白在细胞间通讯、免疫识别、信号转导和微环境重塑中扮演核心角色。"
        "免疫检查点分子（如PD-L1）、生长因子受体（如EGFR、HER2）、致病性免疫球蛋白（如IgE）"
        "以及多种自身抗原-抗体复合物，均属于这一蛋白类别。它们往往以高丰度、多价态或快速再合成的方式"
        "参与疾病进程，使得单纯「阻断」难以从根本上逆转病理状态。因此，开发能够将这些蛋白"
        "主动运送至溶酶体进行彻底降解的新技术，既是靶向蛋白降解领域的重大科学问题，"
        "也是具有明确临床需求的药物研发方向。",
        "2020年，Bertozzi团队在Nature杂志上首次提出溶酶体靶向嵌合体（LYTAC）概念，"
        "开创了利用溶酶体靶向受体介导胞外及膜蛋白降解的新范式。此后五年间，"
        "全球多个研究团队从分子设计、受体选择、细胞机制、疾病模型及产业化等多个层面持续推进，"
        "使LYTAC迅速成为靶向蛋白降解领域最具活力的分支之一。本文旨在从生物医学在药物研发中"
        "发挥关键作用这一主线出发，系统综述LYTAC技术的发展脉络、代表性研究进展及未来展望，"
        "以期为相关领域的科研人员及药物研发从业者提供参考。",
    ],
    "body1": [
        "LYTAC是一类异双功能分子（heterobifunctional molecule），其设计理念源于「分子胶」"
        "和「分子伴侣」（molecular glue）等诱导邻近效应（induced proximity）策略，"
        "但作用通路指向溶酶体而非蛋白酶体。典型的LYTAC分子由三部分构成："
        "（1）靶蛋白结合模块（Target Binder），可为小分子、肽段、抗体、适体（aptamer）等；"
        "（2）溶酶体靶向受体配体模块（LTR Ligand），用于结合细胞表面的LTR；"
        "（3）连接子（Linker），将上述两个功能模块以适当的空间距离和柔性连接起来。",
        "LYTAC发挥功能的生物学过程可概括为「识别-内吞-分选-降解」四个连续步骤（图1）。"
        "首先，LYTAC在细胞外或细胞表面同时结合靶蛋白和LTR，形成三元复合物"
        "（LTR-LYTAC-靶蛋白）。随后，该复合物通过网格蛋白介导的内吞作用（clathrin-mediated endocytosis）"
        "进入早期内体（early endosome）。在内体成熟过程中，复合物被分选至晚期内体/多泡体（MVB），"
        "最终与溶酶体融合，靶蛋白在溶酶体蛋白酶作用下被彻底水解。与抗体中和相比，"
        "LYTAC的优势在于：其一，事件驱动（event-driven）机制使单个降解分子可循环处理多个靶蛋白分子；"
        "其二，不依赖靶蛋白的催化活性或信号传导功能，仅需结合即可实现降解；"
        "其三，彻底消除靶蛋白可同时阻断其信号传导、蛋白-蛋白相互作用及支架功能。",
        "目前已知的LTR主要包括：阳离子非依赖性甘露糖-6-磷酸受体（CI-M6PR/IGF2R）、"
        "去唾液酸糖蛋白受体（ASGPR）、脱唾液酸糖蛋白受体家族成员、整合素、"
        "以及多种细胞因子受体等。不同LTR在组织分布、内吞速率、循环回收特性等方面存在显著差异，"
        "这为LYTAC的组织靶向设计提供了生物学基础，也使受体选择成为药物研发中的关键决策节点。",
    ],
    "body2": [
        "LYTAC的药物研发并非单纯的化学偶联问题，而是深度依赖对细胞生物学过程的系统理解。"
        "以下从生物医学视角阐述影响LYTAC降解效率的核心设计要素。",
        "（一）靶蛋白结合模块的选择。靶蛋白结合剂需具备足够的亲和力和特异性，"
        "但不必抑制靶蛋白功能——这是LYTAC相对传统药物的重要优势。"
        "对于膜蛋白，抗体Fab片段、单链可变区（scFv）及纳米抗体是常用选择；"
        "对于分泌蛋白，可采用全抗体或高亲和力肽段；对于难以获得抗体的靶点，"
        "DNA/RNA适体及小分子配体提供了更灵活的替代方案。值得注意的是，"
        "结合表位、亲和力的微小差异可显著影响三元复合物的形成及内吞效率，"
        "提示靶蛋白结合模块的优化必须结合结构生物学和细胞功能实验进行迭代。",
        "（二）LTR配体与价态设计。以CI-M6PR为例，其天然配体为甘露糖-6-磷酸（M6P）修饰的溶酶体酶，"
        "LYTAC通常采用多价M6P聚糖（M6Pn，n=3-9）作为配体以提高受体结合亲和力。"
        "Ahn等的研究表明，三价N-乙酰半乳糖胺（tri-GalNAc）与ASGPR的结合效能优于单价配体，"
        "为肝靶向LYTAC设计提供了重要参数。配体的化学合成、糖基化修饰及与蛋白载体的偶联位点"
        "均会影响分子的均一性、免疫原性及体内稳定性。",
        "（三）连接子与分子拓扑。连接子的长度、柔性及化学性质决定两个功能模块的空间取向，"
        "直接影响三元复合物的形成概率。过短的连接子可能导致空间位阻，"
        "而过长的连接子则可能降低有效浓度并增加非特异性结合。"
        "近年来，定点偶联（site-specific conjugation）技术的发展，"
        "使抗体-聚糖偶联产物的均一性显著提高，改善了药代动力学性质。",
        "（四）细胞层面的降解决定因素。2023年Ahn等在Science发表的全基因组CRISPR筛选研究揭示，"
        "LYTAC的降解效率受到多条细胞通路的协同调控，包括retromer复合体介导的循环回收、"
        "CUL3的neddylation修饰以及M6P生物合成通路对CI-M6PR的配体占位等。"
        "这些发现表明，同一LYTAC分子在不同细胞系或组织中的降解效能可能存在显著差异，"
        "为临床疗效预测和患者分层提供了潜在的生物标志物。",
    ],
    "case1": [
        "2020年，斯坦福大学Bertozzi实验室在Nature杂志发表了LYTAC领域的奠基性工作，"
        "标志着靶向蛋白降解技术正式拓展至胞外及膜蛋白领域。该研究由Banik、Pedram、Wisnovsky等完成，"
        "论文题目为《Lysosome-targeting chimaeras for degradation of extracellular proteins》，"
        "发表于Nature第584卷第7820期，页码291-297。",
        "该团队将靶向蛋白降解的「诱导邻近」理念与溶酶体分选途径相结合，设计并合成了首批LYTAC分子。"
        "这些分子由靶向模块（小分子或抗体）与CI-M6PR激动剂——化学合成的甘露糖-6-磷酸聚糖配体（M6Pn）"
        "通过连接子偶联而成。在机制验证方面，研究团队建立了基于CRISPR筛选的功能基因组学平台，"
        "系统鉴定了CI-M6PR依赖的LYTAC内吞通路，并意外发现外泌体复合体（exocyst complex）"
        "是LYTAC介导内吞的关键组分，为理解受体介导的内吞提供了新的生物学视角。",
        "在靶蛋白降解验证中，该研究展示了LYTAC对多种具有重要治疗意义的蛋白的降解能力，"
        "包括载脂蛋白E4（ApoE4，阿尔茨海默病遗传风险因子）、表皮生长因子受体（EGFR）、"
        "转铁蛋白受体（CD71）及程序性死亡配体1（PD-L1）等。以EGFR为例，"
        "LYTAC可在数小时内实现受体蛋白的显著下调，其降解动力学与溶酶体抑制剂的处理结果一致，"
        "证实了溶酶体途径在降解中的核心作用。该研究还证明，LYTAC对PD-L1的降解可有效削弱"
        "肿瘤细胞的免疫逃逸能力，提示其在肿瘤免疫治疗中的潜在应用价值。",
        "从药物研发视角审视，该工作具有三重重要意义：第一，在概念层面证明了「降解优于抑制」"
        "原则可推广至胞外蛋白；第二，在工具层面提供了可模块化替换的分子架构，"
        "为后续多样化LYTAC设计奠定了基础；第三，在生物学层面揭示了LYTAC内吞的分子机制，"
        "体现了基础生物医学研究对技术创新方向的引领作用。该论文发表后迅速成为LYTAC领域的"
        "引用基石，截至2025年已被引用超过千次，并直接催生了Lycia Therapeutics等"
        "专注于LYTAC药物开发的企业。",
    ],
    "case2": [
        "第一代LYTAC采用CI-M6PR作为溶酶体靶向受体，该受体在多种组织广泛表达，"
        "限制了LYTAC的组织选择性。实现组织特异性蛋白降解是LYTAC走向临床应用的关键需求，"
        "尤其在肿瘤和代谢性疾病中，脱靶降解可能带来严重毒副作用。",
        "2021年，Ahn、Banik、Miller和Bertozzi等在Nature Chemical Biology（第17卷，页码937-946）"
        "发表了题为《LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation》"
        "的研究。该工作将LYTAC的LTR模块从CI-M6PR拓展至去唾液酸糖蛋白受体（ASGPR），"
        "后者是主要在肝细胞表面高表达的LTR，已被广泛验证为肝靶向递送的有效锚点"
        "（如GalNAc-siRNA偶联物）。",
        "研究团队将靶蛋白结合剂与三价N-乙酰半乳糖胺（tri-GalNAc）配体偶联，构建了GalNAc-LYTAC平台。"
        "在肝细胞中，GalNAc-LYTAC可有效降解EGFR，且其信号阻断效应优于传统EGFR抑制性抗体，"
        "表明蛋白清除较受体占据可产生更持久的功能抑制。进一步地，研究者将3.4 kDa的肽段结合剂"
        "与tri-GalNAc连接，实现了对整合素的降解及癌细胞增殖的抑制。"
        "该研究还系统优化了抗体骨架上的定点偶联策略：在抗体Fc区特定位点引入单个tri-GalNAc配体，"
        "显著改善了GalNAc-LYTAC在体内的药代动力学性质，延长了半衰期并提高了生物利用度。",
        "这一研究从生物医学与药物化学的交叉视角，解决了LYTAC研发中的两个核心问题："
        "组织选择性与成药性。ASGPR的肝特异性表达使GalNAc-LYTAC成为肝相关疾病"
        "（如肝细胞癌、代谢性肝病）的理性设计选择；而定点偶联技术的引入则表明，"
        "蛋白质工程与糖化学的协同优化是提升LYTAC药物性质的关键路径。"
        "该工作也为后来cataLYTAC等平台的ASGPR模块选择提供了直接依据。",
    ],
    "case3": [
        "随着LYTAC分子设计的日趋成熟，研究者逐渐认识到：相同结构的LYTAC在不同细胞类型中的降解效率"
        "可能存在数量级差异。这一现象提示，细胞内在的生物学特征——而非仅仅是分子结构——"
        "是决定LYTAC疗效的关键因素。阐明这些细胞决定因素，对于LYTAC的理性设计、"
        "临床疗效预测及患者分层具有至关重要的意义。",
        "2023年，Ahn、Banik和Bertozzi团队在Science杂志（第382卷，文章编号eadf6249）"
        "发表了题为《Elucidating the cellular determinants of targeted membrane protein degradation "
        "by lysosome-targeting chimeras》的研究。该团队利用全基因组CRISPR敲除筛选技术，"
        "系统鉴定了调控LYTAC介导膜蛋白降解的正向及负向调节因子，"
        "为LYTAC药物研发提供了迄今最全面的细胞生物学图谱。",
        "该研究的核心发现包括三个方面。第一，retromer复合体介导的循环回收是LYTAC降解的重要负调控通路："
        "retromer将LYTAC-CI-M6PR复合物从早期内体回收至细胞表面，"
        "使LYTAC在发挥降解功能前即脱离内吞途径；敲除retromer基因可显著增强LYTAC的降解效率，"
        "在某些条件下使EGFR降解率提升至90%以上。第二，CUL3的neddylation修饰是LYTAC-靶蛋白复合物"
        "从内体向溶酶体成熟转运的关键步骤；在11种来自8种不同组织的人源细胞系中，"
        "neddylated CUL3的表达水平与EGFR降解效率呈正相关，提示其可作为LYTAC疗效预测的潜在生物标志物。"
        "第三，M6P生物合成通路通过为CI-M6PR提供内源性竞争性配体，"
        "占据受体结合位点从而抑制LYTAC活性；抑制M6P生物合成可增加细胞表面游离CI-M6PR的比例，"
        "增强LYTAC的受体结合和内化效率。",
        "该研究深刻体现了生物医学基础研究的「反哺」价值：它不仅解释了LYTAC在不同细胞中"
        "效能差异的分子基础，更为下一代LYTAC的设计提供了明确方向——"
        "可通过调控retromer活性、联合CUL3 neddylation调节剂或优化M6P配体竞争等方式增强降解效率。"
        "从药物研发管线角度，neddylated CUL3水平的检测有望成为LYTAC疗法的伴随诊断指标，"
        "推动LYTAC从「一刀切」走向精准医疗。",
    ],
    "case4": [
        "化学合成的LYTAC在分子均一性和规模化生产方面面临挑战，"
        "这促使研究者探索基因编码的蛋白降解平台。2023年，"
        "加州大学旧金山分校Wells实验室的Pance、Gramespacher、Byrnes等在"
        "Nature Biotechnology（第41卷第2期，页码273-281）发表了KineTAC（细胞因子受体靶向嵌合体）"
        "平台的研究，为LYTAC技术提供了重要的生物学替代方案。",
        "KineTAC是完全基因编码的双特异性抗体，由两个功能臂组成："
        "细胞因子臂结合内源性细胞因子受体（如CXCR7），靶蛋白结合臂识别目标蛋白。"
        "以CXCL12-CXCR7轴为例，CXCL12与CXCR7结合后可触发受体内吞并导向溶酶体降解，"
        "KineTAC将这一天然生物学过程「劫持」用于靶蛋白的共转运降解。"
        "研究团队针对PD-L1、HER2、PD-1、EGFR、CDCP1和TROP2等8种治疗相关靶蛋白"
        "构建了相应的KineTAC分子，最大降解效率（Dmax）介于51%至93%之间，"
        "展示了平台的广泛适用性。",
        "与化学LYTAC相比，KineTAC具有多方面优势：模块化的基因编码设计使靶蛋白结合臂"
        "可快速替换，无需复杂的化学合成和糖-蛋白偶联；可利用哺乳动物细胞表达系统生产，"
        "保证了产品的均一性和可放大性；通过选择不同细胞因子-受体轴"
        "（如CXCL11-CXCR7、IL-2-IL-2R），可实现对不同组织或细胞类型的靶向。"
        "以PD-L1降解为例，KineTAC不仅消除了免疫检查点信号，"
        "还可能通过改变膜蛋白组构成影响肿瘤免疫微环境，为免疫治疗提供了不同于抗体阻断的新机制。",
        "KineTAC的出现表明，LYTAC药物研发不必局限于化学偶联路径，"
        "融合蛋白质工程、免疫学和细胞生物学的生物降解平台同样具有强大的竞争力。"
        "从药物开发角度，KineTAC更接近传统抗体药物的开发和质控流程，"
        "可能更快实现从实验室到生产的转化，体现了多学科交叉融合对创新药物研发的推动作用。",
    ],
    "body4": [
        "基于上述基础研究积累，LYTAC技术已在多个疾病领域展现出转化潜力，"
        "以下从生物医学与药物研发结合的角度进行梳理。",
        "（一）过敏及免疫性疾病。IgE是I型过敏反应的关键效应分子，"
        "目前临床标准治疗药物奥马珠单抗（omalizumab）通过中和游离IgE发挥作用，"
        "但无法降低总IgE水平，且疗效受患者基线IgE浓度限制。"
        "Lycia Therapeutics基于Bertozzi实验室的技术开发了cataLYTAC平台，"
        "通过稳定化ASGPR配体、pH敏感性靶蛋白结合及FcRn介导的循环回收，"
        "实现了对IgE的催化性（catalytic）降解。临床前研究显示，"
        "cataLYTAC在体外可降解超化学计量水平的IgE，"
        "在非人灵长类动物中单次给药即可实现超过98%的内源性IgE清除并维持约两周，"
        "疗效和持续时间均优于奥马珠单抗。候选药物LCA-0061已进入IND申报阶段，"
        "拟用于食物过敏、过敏性哮喘等疾病，有望成为首个进入临床试验的LYTAC药物。",
        "在自身免疫病领域，致病性自身抗体（如Graves病中的TSH受体抗体、"
        "重症肌无力中的抗MuSK抗体）是理想的LYTAC靶点。"
        "Lycia的LCA-0321和LCA-0391分别靶向TSHR自身抗体和抗MuSK抗体，"
        "旨在实现抗原特异性抗体清除而不引起广泛免疫抑制，"
        "体现了LYTAC在精准免疫调节中的独特优势。",
        "（二）肿瘤免疫治疗。PD-L1/PD-1轴的免疫检查点阻断疗法已改变多种肿瘤的治疗格局，"
        "但客观缓解率仍有限。2023年，Li等在Journal of the American Chemical Society"
        "（第145卷第45期，页码24506-24521）报道了基于DNA适体的共价LYTAC（covalent LYTAC），"
        "通过生物正交共价连接强化PD-L1结合，实现了高效的PD-L1溶酶体降解。"
        "该研究首次证明LYTAC介导的PD-L1降解可诱导肿瘤细胞的免疫原性凋亡，"
        "激活肿瘤特异性免疫应答，在抗肿瘤疗效与炎症损伤之间取得了优于抗体疗法的平衡。"
        "这一工作将适体技术、化学生物学与肿瘤免疫学相结合，"
        "为免疫检查点降解（Immune Checkpoint Degradation, ICD）疗法开辟了新路径。",
        "（三）神经退行性疾病。胞外蛋白聚集体（如β淀粉样蛋白、Tau蛋白）"
        "是阿尔茨海默病等神经退行性疾病的病理标志物。2023年，"
        "中国科学技术大学团队在Cell姊妹刊Chem（第9卷，页码2016-2038）发表了"
        "生物标志物激活的多功能LYTAC（KPLY）研究，设计了可穿越血脑屏障的"
        "聚多巴胺纳米载体，负载点击化学前体，在Aβ-Cu复合物催化下于病灶区域原位生成活性LYTAC，"
        "通过CD206受体介导Aβ多态性聚集体（寡聚体、原纤维等）的溶酶体降解。"
        "该研究巧妙利用了阿尔茨海默病病灶区铜离子异常积累的病理特征，"
        "实现了LYTAC的时空可控激活，降低了全身给药的外周脱靶风险，"
        "展示了生物医学对疾病微环境的深刻理解如何驱动创新型药物设计。",
    ],
    "body5": [
        "尽管LYTAC技术取得了令人瞩目的进展，但从实验室到临床的转化仍面临多重挑战。",
        "第一，分子设计与制备的复杂性。含聚糖链的LYTAC分子结构复杂，"
        "化学合成和生物偶联的产物异质性较高，给质量控制和大规模生产带来困难。"
        "如何在保证降解活性的同时简化分子结构、提高批次间一致性，是CMC（化学、生产和控制）"
        "层面的核心问题。",
        "第二，组织选择性与脱靶效应。广泛表达的CI-M6PR可能导致全身性靶蛋白降解，"
        "引发非预期的组织毒性。虽然ASGPR等组织特异性LTR提供了部分解决方案，"
        "但对于非肝脏靶点，仍需开发更多器官特异性受体或前药激活策略。",
        "第三，药代动力学与免疫原性。聚糖配体可能激活先天免疫系统，"
        "异源蛋白载体可能诱导抗药抗体（ADA）反应。第一代LYTAC与靶蛋白共降解导致"
        "治疗分子无法循环利用，限制了体内暴露和降解容量；cataLYTAC等催化性设计"
        "正是针对这一瓶颈的生物医学工程解决方案。",
        "第四，临床验证的缺乏。截至2025年，尚无LYTAC药物获批上市或完成临床试验，"
        "其人体安全性、有效性和最佳给药方案仍需通过系统的临床研究开发。"
        "细胞层面的降解标志物（如neddylated CUL3）能否在人体中可靠预测疗效，"
        "也需要临床数据的验证。",
    ],
    "conclusion": [
        "溶酶体靶向嵌合体（LYTAC）技术作为靶向蛋白降解领域的重要突破，"
        "将蛋白降解的适用范围从胞内拓展至约占人类蛋白质组40%的胞外及膜蛋白，"
        "为「不可成药」靶点的药物开发提供了全新范式。本文系统综述了LYTAC的基本原理、"
        "设计要素、代表性研究进展及疾病应用前景，重点通过Nature、Science、"
        "Nature Chemical Biology和Nature Biotechnology等期刊的四项高水平研究，"
        "展示了从概念奠基、组织靶向、机制解析到平台模块化的技术演进路径。",
        "贯穿这些研究的主线，是生物医学对药物研发的深度赋能。"
        "没有对内吞-溶酶体途径的系统理解，就无法建立LYTAC的分子设计逻辑；"
        "没有对ASGPR肝特异性表达的认知，就难以实现组织靶向降解；"
        "没有全基因组筛选揭示的细胞调控网络，就无法解释疗效差异并开发伴随诊断；"
        "没有对细胞因子受体循环通路的借鉴，就难以诞生KineTAC等基因编码平台。"
        "这些实例充分说明，创新药物的研发并非化学或药理学的「单兵突进」，"
        "而是需要细胞生物学、免疫学、糖生物学、蛋白质组学等多学科知识的深度融合。",
        "从笔者对领域发展的思考来看，LYTAC技术的未来可能呈现以下趋势：",
        "其一，从「降解」到「催化降解」的范式升级。cataLYTAC通过FcRn循环实现治疗分子的再利用，"
        "突破了第一代LYTAC「一次性消耗」的药代动力学瓶颈，"
        "这一创新本质上是对受体循环生物学和抗体工程学的综合运用，"
        "预示着LYTAC药物将走向更高的催化效率和更低的给药剂量。",
        "其二，从「广谱靶向」到「精准靶向」的深化。随着更多LTR（如GPC3、整合素、CXCR7等）"
        "被纳入LYTAC设计，联合病灶微环境激活策略（如Aβ-Cu触发的原位生成），"
        "LYTAC有望实现器官、细胞类型乃至病理微环境的多层次精准靶向，"
        "最大限度拓展治疗窗。",
        "其三，从「实验室工具」到「临床药物」的跨越。LCA-0061等候选药物进入临床开发阶段，"
        "标志着LYTAC技术正经历从0到1的关键跃迁。未来5-10年，"
        "临床试验数据将回答LYTAC在人体中的安全性、免疫原性及疗效等核心问题，"
        "也将验证细胞生物标志物指导的精准用药策略是否可行。",
        "其四，从「单一技术」到「融合平台」的拓展。LYTAC与PROTAC、AUTAC、分子胶等"
        "其他TPD技术的互补性，使同时靶向胞内、胞外及膜蛋白的联合降解策略成为可能；"
        "与适体、抗体、小分子及基因治疗载体的融合，将进一步丰富LYTAC的分子形式和应用场景。",
        "综上所述，LYTAC技术代表了生物医学研究与药物研发深度融合的典范。"
        "它既拓展了人类对抗疾病靶点的能力边界，也提醒我们："
        "真正具有变革性的药物创新，往往源于对生命科学基本问题的深刻洞察"
        "与临床需求的紧密结合。随着基础研究的持续深入、转化平台的日趋成熟"
        "及监管科学的逐步完善，LYTAC有望成为继PROTAC之后靶向蛋白降解领域的"
        "又一里程碑，为过敏性疾病、自身免疫病、肿瘤及神经退行性疾病患者"
        "带来新的治疗希望。",
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
    "[8] Wang Y, Song Y, Liu Z, et al. Targeted degradation of membrane and extracellular proteins with LYTACs[J]. Acta Pharmacologica Sinica, 2024, 45(7): 1333-1345.",
    "[9] Sun D, Lu Y, Hu Y, et al. Targeted protein degradation: advances in drug discovery and clinical practice[J]. Signal Transduction and Targeted Therapy, 2024, 9: 308.",
    "[10] Caianiello D F, Miller C L, Ahn G, Riley N M, Bertozzi C R. Bifunctional small molecules that mediate the degradation of extracellular proteins[J]. Nature Chemical Biology, 2021, 17(8): 947-953.",
    "[11] Ji C H, Kim H Y, Lee M J. The AUTOTAC chemical biology platform for targeted protein degradation via the autophagy-lysosome system[J]. Nature Communications, 2022, 13(1): 904.",
    "[12] Cotton A D, Nguyen D P, Gramespacher J A, Seiple I B, Wells J A. Development of antibody-based PROTACs for the degradation of the cell-surface immune checkpoint protein PD-L1[J]. Journal of the American Chemical Society, 2021, 143(2): 593-598.",
    "[13] Lycia Therapeutics. LYTAC Platform[EB/OL]. https://lyciatx.com/technology/lytac-platform/, 2025.",
    "[14] Lycia Therapeutics. Pipeline Programs[EB/OL]. https://lyciatx.com/pipeline/programs/, 2025.",
    "[15] Bertozzi C R, et al. Catalytic degradation of circulating targets with FcRn-mediated cycling LYTACs (cataLYTACs)[R]. bioRxiv, 2025. doi: 10.1101/2025.01.12.632472.",
    "[16] Bekes M, Langley D R, Crews C M. PROTAC targeted protein degraders: the past is prologue[J]. Nature Reviews Drug Discovery, 2022, 21(3): 181-200.",
    "[17] Zhou Y, Zhang Y, Lazerwith S E, et al. Exploring the target space of lysosome-targeting chimeras[J]. Journal of Medicinal Chemistry, 2024, 67(5): 3654-3675.",
    "[18] Wang Y, Liu Z, Deng Q, et al. Unlocking the undruggable: current landscape and emerging frontiers in lysosomal receptor-mediated protein degradation[J]. Drug Delivery, 2026, 33(1): 2636323.",
]


def main():
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(3.17)
    section.right_margin = Cm(3.17)

    add_heading_text(doc, "溶酶体靶向嵌合体（LYTAC）技术在药物研发中的进展与展望", level=1)
    add_heading_text(doc, "——生物医学视角下的靶向蛋白降解新范式", level=1)

    add_heading_text(doc, "摘要", level=2)
    add_paragraph(doc, CONTENT["abstract"], first_line_indent=False)
    kp = doc.add_paragraph()
    kp.paragraph_format.first_line_indent = Cm(0)
    set_run_font(kp.add_run("关键词：溶酶体靶向嵌合体；靶向蛋白降解；溶酶体靶向受体；药物研发；生物医学"), bold=True)

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
