#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate LYTAC review document with specified formatting."""

import re

from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn

# 参考文献编号与英文题名（正文标注用，不计入正文字数）
REF_TITLES = {
    1: "Lysosome-targeting chimaeras for degradation of extracellular proteins",
    2: "LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation",
    3: "Modular cytokine receptor-targeting chimeras for targeted degradation of cell surface and extracellular proteins",
    4: "Covalent LYTAC enabled by DNA aptamers for immune checkpoint degradation therapy",
    5: "Biomarker-activated multifunctional lysosome-targeting chimeras mediated selective degradation of extracellular amyloid fibrils",
    6: "Aptamer-LYTACs for targeted degradation of extracellular and membrane proteins",
    7: "Targeted degradation of membrane and extracellular proteins with LYTACs",
    8: "Targeted degradation of extracellular secreted and membrane proteins",
    9: "Bifunctional small molecules that mediate the degradation of extracellular proteins",
    10: "PROTAC targeted protein degraders: the past is prologue",
    11: "Emerging protein degradation strategies: expanding the scope to extracellular and membrane proteins",
    12: "Designed endocytosis-inducing proteins degrade targets and amplify signals",
    13: "Insulin-like Growth Factor 2 (IGF2)-Fused Lysosomal Targeting Chimeras for Degradation of Extracellular and Membrane Proteins",
    14: "Lysosome-Targeting Chimera Using Mannose-6-Phosphate Glycans Derived from Glyco-Engineered Yeast",
    15: "Lysosome-targeting chimeras containing an endocytic signaling motif trigger endocytosis and lysosomal degradation of cell-surface proteins",
    16: "Targeted degradation of extracellular proteins: state of the art and diversity of degrader designs",
    17: "Targeted protein degradation: advances in drug discovery and clinical practice",
    18: "Targeted protein degradation: mechanisms, strategies and application",
}


def cite(*nums):
    """在段落文本中暂存文献编号；实际显示时集中移到段末。"""
    return "".join(f"{{CITE:{n}}}" for n in nums)


def format_citations(text):
    """将分散在句中的引用集中为段末文献名标注，使正文阅读更连贯。"""
    refs = []

    def collect(match):
        ref_id = int(match.group(1))
        if ref_id not in refs:
            refs.append(ref_id)
        return ""

    body = re.sub(r"\{CITE:(\d+)\}", collect, text)
    if not refs:
        return body
    note = "（参考文献：" + "；".join(f"《{REF_TITLES[n]}》[{n}]" for n in refs) + "）"
    return body + note


def count_chinese_chars(text):
    """统计汉字数；排除文献名标注块"""
    cleaned = re.sub(r"\{CITE:\d+\}", "", text)
    cleaned = re.sub(r"（参考文献：[^）]*）", "", cleaned)
    return sum(1 for c in cleaned if "\u4e00" <= c <= "\u9fff")


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
    run = p.add_run(format_citations(text))
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
    "溶酶体靶向嵌合体（LYTAC）通过偶联靶蛋白与溶酶体靶向受体，"
    "将胞外及膜蛋白导入溶酶体降解，是面向胞外靶点的重要降解策略。"
    "本综述依受体生物学基础、模块化分子设计、一区研究实例、"
    "疾病应用验证和转化瓶颈五个层面展开论述，"
    "分析受体选择对降解效率与组织定位的决定作用，"
    "归纳靶向臂、受体臂和连接子三模块的平台属性，"
    "并以Banik等（Nature, 2020）、Ahn等（Nat Chem Biol, 2021）、"
    "Pance等（Nat Biotechnol, 2023）和Li等（JACS, 2023）四篇一区论文为实例，"
    "评估降解策略相对单纯阻断的增益及成药前主要障碍。"
    + cite(1, 2, 3, 4, 8, 16)
)

INTRO = [
    "靶向蛋白降解为传统上难以药化的蛋白提供了新路径。"
    "PROTAC通过招募胞内E3连接酶将靶蛋白导向蛋白酶体降解，"
    "已有候选分子进入临床评估"
    + cite(10)
    + "；"
    "然而该策略对分泌蛋白和跨膜蛋白等胞外靶点覆盖有限。",
    "2020年Banik等提出溶酶体靶向嵌合体（LYTAC），"
    "通过同时结合靶蛋白与细胞表面溶酶体靶向受体，"
    "经内吞-溶酶体途径实现蛋白清除"
    + cite(1)
    + "。"
    "本综述在简述上述背景后，"
    "依次讨论受体选择决策、三模块分子工程、"
    "一区代表性研究实例、疾病应用证据及转化瓶颈，"
    "以期为LYTAC药物研发提供逻辑清晰的分析框架"
    + cite(8, 16)
    + "。",
]

BODY_SECTION1 = [
    "LYTAC的降解效率首先取决于所选用LTR的生物学特性，而非仅由化学结构决定。"
    "三元复合物形成后，经网格蛋白介导内吞进入早期内体；"
    "此后复合物或被分选至溶酶体完成蛋白水解，"
    "或通过retromer等途径回收至细胞膜，两条路径的竞争直接决定有效降解通量（图1）。"
    "Lin等和Mamun等的综述均将这一分选-回收平衡视为解释不同细胞系降解效率差异的核心变量"
    + cite(11, 16)
    + "。",
    "与抗体阻断不同，LYTAC属于事件驱动型策略："
    "结合剂无需抑制靶蛋白活性，只要有效捕获即可启动降解程序。"
    "降解发生在溶酶体后，靶蛋白的膜信号传导和胞外相互作用可同时终止，"
    "这一机制差异是后续“降解是否优于阻断”比较的生物学基础"
    + cite(7, 8)
    + "。",
    "验证溶酶体途径的标准实验体系包括：溶酶体抑制剂阻断降解、"
    "受体敲除或配体竞争削弱效应，以及遗传学筛选鉴定内吞分选相关基因。"
    "Banik等在初代研究中用溶酶体抑制剂确认EGFR降解依赖溶酶体，"
    "并以CRISPR干扰筛选鉴定CI-M6PR内吞的细胞机器"
    + cite(1)
    + "。",
    "从受体组织分布出发，可建立LYTAC的受体选择决策框架。"
    "第一类是以CI-M6PR（IGF2R）为代表的广谱受体："
    "其在多种组织均有表达，适合需要广泛降解的靶点，"
    "但脱靶风险较高，且内源性M6P修饰蛋白可占用受体结合位点，"
    "削弱外源LYTAC的降解效率"
    + cite(1, 12)
    + "。"
    "Banik等通过CRISPR筛选还发现，exocyst复合体影响CI-M6PR的膜面呈现，"
    "提示受体可及性本身即是降解效率的上游决定因素"
    + cite(1)
    + "。",
    "第二类是以ASGPR为代表的组织定向受体："
    "其表达集中于肝细胞，tri-GalNAc配体可将降解限制在肝脏微环境，"
    "Ahn等的GalNAc-LYTAC和Wu等的Apt-LYTAC均利用这一策略实现肝靶向降解"
    + cite(2, 6)
    + "。"
    "当治疗目标为肝内靶蛋白或需降低全身脱靶时，应优先评估ASGPR通路；"
    "当靶蛋白分布广泛或需作用于非肝组织时，则需转向其他LTR或蛋白工程策略。",
    "第三类是以细胞因子受体和人工内吞模块为代表的可编程受体臂："
    "Pance等的KineTAC利用CXCL12-CXCR7等天然内吞过程触发靶蛋白降解，"
    "Huang等的EndoTag和Zhang等的IGF2融合蛋白（iLYTAC）"
    "则分别通过计算设计和IGF2配体结合CI-M6PR实现溶酶体导向"
    + cite(3, 12, 13)
    + "。"
    "不同EndoTag因结合受体组织分布不同，可实现降解的组织定向；"
    "iLYTAC则通过基因编码方式绕开复杂聚糖化学合成。",
    "Huang等的研究进一步揭示，"
    "内源性M6P糖蛋白与LYTAC竞争CI-M6PR结合位点，"
    "是广谱受体路线中不可忽视的生物学约束"
    + cite(12)
    + "。"
    "该发现意味着：即便分子亲和力很高，"
    "细胞内源性受体占用仍可能限制降解通量；"
    "解读体外降解数据时需同步考虑受体竞争状态。",
    "第四类是绕过外源LTR配体、直接利用内吞分选信号的工程策略："
    "Fang等的SignalTAC将CI-M6PR来源的酪氨酸基内吞信号肽P3与靶蛋白结合模块融合，"
    "在不依赖外源聚糖配体的情况下诱导靶蛋白溶酶体降解"
    + cite(15)
    + "。"
    "该路线将“受体臂”的功能转化为内源性分选信号，"
    "拓展了受体选择的逻辑边界。",
    "综合上述四类路径，受体选择可遵循以下决策逻辑："
    "若追求最大靶蛋白覆盖且可接受脱靶风险，优先考虑CI-M6PR/M6P配体体系；"
    "若需肝定向降解，选择ASGPR/GalNAc或肝靶向适体偶联物；"
    "若需组织特异性且希望基因编码生产，评估EndoTag或iLYTAC；"
    "若化学合成聚糖受限，可考察酵母来源M6P聚糖或内吞信号肽融合策略"
    + cite(14, 15)
    + "。"
    "Chen等的综述从效应器类型角度进一步归纳了LTR、整合素及细胞因子受体等可选方案"
    + cite(8)
    + "，"
    "为上述决策提供了靶点-受体配对参考。",
    "Zhao等将溶酶体途径定位为拓展靶向蛋白降解（TPD）靶点范围的重要方向，"
    "指出胞外及膜蛋白约占蛋白质组的相当比例，"
    "是LYTAC存在的结构性需求"
    + cite(18)
    + "。"
    "在实际决策中，建议按以下顺序评估："
    "（1）靶蛋白的细胞外定位与更新速率；"
    "（2）靶标组织的LTR表达谱；"
    "（3）可接受的脱靶降解范围；"
    "（4）制备平台与受体臂类型的匹配性。"
    "只有四项指标对齐后，方可进入靶向臂和连接子优化。",
    "Lin等的综述进一步指出，"
    "整合素、清道夫受体等亦被用于组织选择性胞外蛋白降解，"
    "提示LTR的选择空间仍在持续扩展"
    + cite(11)
    + "。"
    "对于膜蛋白靶点，还需关注其与LTR在同一细胞表面的空间距离和共定位关系，"
    "因为这将影响三元复合物形成的几何可行性。",
    "从疾病相关靶点分布看，"
    "PD-L1、EGFR和HER2等膜蛋白在肿瘤中表达异质性高，"
    "CI-M6PR广谱降解可能带来肿瘤外组织脱靶；"
    "ApoE4和Aβ等靶点则涉及分泌蛋白或聚集体的清除，"
    "对受体通量和局部给药策略提出不同要求"
    + cite(1, 5, 8)
    + "。"
    "因此，受体决策框架必须与具体疾病场景联用，"
    "不能脱离靶点生物学单独优化分子结构。",
    "在评价降解效率时，"
    "还需区分“受体结合亲和力”与“有效降解通量”两个层次。"
    "高亲和力仅保证三元复合物形成，"
    "而通量取决于内吞速率、溶酶体分选效率和受体回收动力学的综合平衡。"
    "同一LYTAC在不同细胞系中降解效率的差异，"
    "往往反映的是后者而非前者的问题"
    + cite(11, 16)
    + "。",
    "本节小结：LYTAC疗效的关键瓶颈不在于靶蛋白能否被结合，"
    "而在于所选LTR的组织分布、内吞通量及内体分选倾向是否与疾病场景匹配。"
    "受体选择应作为药物设计的首要决策节点，而非分子构建完成后的附带调整。",
    "Caianiello等的MoDE-A研究从另一角度丰富了受体决策框架："
    "小分子亦可同时承担靶向臂和受体结合功能，"
    "以ASGPR为受体降解整合素等靶点，"
    "为分子量敏感的应用场景提供了补充选项"
    + cite(9)
    + "。"
    "Kim等则证明，"
    "即使沿用CI-M6PR受体臂，"
    "聚糖配体的来源和制备方式也可显著影响分子的可开发性，"
    "酵母发酵来源的gyM6pG聚糖在制备可及性上优于化学合成多价M6P聚糖"
    + cite(14)
    + "。",
]

BODY_SECTION2 = [
    "在明确受体决策框架后，LYTAC分子可拆解为三个可独立替换的功能模块："
    "靶向臂（结合靶蛋白）、受体臂（触发内吞与溶酶体分选）和连接子（调控模块间距与偶联均一性）。"
    "三模块的解耦设计使LYTAC具备平台化属性——"
    "同一受体臂可搭配不同靶向臂以切换靶点，"
    "同一靶向臂亦可更换受体臂以改变组织分布"
    + cite(8, 16)
    + "。",
    "靶向臂的可选范围已从抗体扩展至小分子、肽段、适体和基因编码结合蛋白。"
    "Banik等在初代LYTAC中分别使用抗体和小分子结合EGFR等靶点，"
    "证明靶向臂化学类型并非固定"
    + cite(1)
    + "。"
    "Caianiello等的MoDE-A以双功能小分子同时结合靶蛋白和ASGPR，"
    "实现了整合素等靶点的溶酶体降解，分子量显著低于抗体偶联物"
    + cite(9)
    + "。"
    "Ahn等以肽段结合剂降解整合素，表明靶向臂不必局限于完整抗体"
    + cite(2)
    + "。"
    "Wu等和Li等分别以适体连接tri-GalNAc或共价结合PD-L1，"
    "在简化靶向臂制备的同时实现了肝细胞或肿瘤微环境中的靶蛋白降解"
    + cite(6, 4)
    + "。"
    "Zhang等的iLYTAC以affibody、纳米抗体或IgG结合Z结构域作为靶向臂，"
    "与IGF2受体臂融合表达，覆盖了EGFR、PD-L1、CD20和α-突触核蛋白等靶点"
    + cite(13)
    + "。",
    "受体臂的模块化同样显著。"
    "化学合成路线以多价M6P聚糖（CI-M6PR）和tri-GalNAc（ASGPR）为代表"
    + cite(1, 2)
    + "；"
    "Kim等以糖工程酵母发酵制备人源兼容的gyM6pG聚糖，"
    "通过无铜点击化学与纳米抗体偶联，降低了聚糖受体臂的合成门槛"
    + cite(14)
    + "。"
    "基因编码路线中，Pance等的KineTAC以细胞因子受体介导内吞，"
    "Huang等的EndoTag以人工设计蛋白结合IGF2R或ASGPR，"
    "Zhang等的iLYTAC以IGF2多肽结合CI-M6PR，"
    "三者均避免了聚糖化学修饰"
    + cite(3, 12, 13)
    + "。"
    "Fang等的SignalTAC则将受体臂功能转化为CI-M6PR来源的内吞信号肽，"
    "在膜蛋白靶向降解中显示出优于亲本抗体的活性"
    + cite(15)
    + "。",
    "连接子模块决定两个功能臂的空间排布与偶联均一性。"
    "Ahn等在抗体Fc区定点偶联tri-GalNAc，改善了GalNAc-LYTAC的体内稳定性"
    + cite(2)
    + "。"
    "Li等在综述中指出，连接子过长不利于三元复合物形成，过短则限制构象调整；"
    "配体价态过高亦可能引发非特异性聚集"
    + cite(7)
    + "。"
    "在模块化平台中，连接子优化应服务于“靶向臂-受体臂”组合的功能匹配，"
    "而非孤立追求某一理化参数。",
    "从平台可替换性看，三模块组合呈现出清晰的“即插即用”特征。"
    "Pance等的KineTAC可通过更换靶蛋白结合臂降解PD-L1、HER2、TROP2等多个靶点，"
    "而保留同一细胞因子受体臂"
    + cite(3)
    + "；"
    "Huang等的EndoTag融合不同靶蛋白结合蛋白后，"
    "可定向至IGF2R或ASGPR等不同受体"
    + cite(12)
    + "；"
    "Zhang等的iLYTAC型-II以IgG结合Z结构域为通用适配器，"
    "可与临床已有抗体偶联实现降解"
    + cite(13)
    + "。"
    "这些设计表明，LYTAC研发正从“逐个分子定制”转向“平台化模块组装”。",
    "Mamun等的综述按效应器类型将eTPD降解剂分为LYTAC、双特异性抗体和小分子等类别，"
    "与本节的三模块框架在逻辑上相互印证"
    + cite(16)
    + "。",
    "按三模块组合逻辑，现有LYTAC变体可归纳为若干典型平台："
    "聚糖-抗体平台（Banik、Kim）、GalNAc-肝靶向平台（Ahn、Wu）、"
    "适体-共价平台（Li）、双特异性蛋白平台（Pance、Huang、Zhang）"
    "和小分子-受体平台（Caianiello）"
    + cite(1, 2, 4, 6, 9, 12, 13, 14)
    + "。"
    "各平台在制备难度、靶点切换灵活性和组织定向能力上各有权衡，"
    "但均遵循同一底层逻辑：通过可替换模块实现靶蛋白向溶酶体的定向转运。",
    "制备路线可按受体臂类型分为三条主线："
    "化学偶联聚糖臂（Banik、Kim）、化学合成小分子臂（Caianiello）"
    "和基因编码蛋白臂（Pance、Huang、Zhang、Fang）。"
    "前者的优势在于结构可精细调控，后者的优势在于批次均一性和靶点切换速度。"
    "连接子优化在三条主线上均为共性挑战，"
    "尤其化学偶联路线中偶联位点和偶联数不一致会直接影响降解活性重现性"
    + cite(1, 9, 12, 14)
    + "。",
    "本节小结：LYTAC并非单一分子类型的线性迭代，"
    "而是以靶向臂、受体臂和连接子三模块为基础的可组合平台。"
    "药物设计应先在受体决策框架内确定受体臂类型，"
    "再匹配适宜的靶向臂和连接子策略。",
    "需要指出的是，模块化并不意味着任意组合均有效。"
    "靶向臂亲和力不足将导致三元复合物形成效率低下；"
    "受体臂与靶组织不匹配将造成降解位置错误；"
    "连接子长度不当则影响两个功能臂的协同构象。"
    "因此，平台化的前提是建立针对特定“靶点-组织-制备”组合的标准优化流程，"
    "而非简单拼接已有模块。",
    "在化学合成与基因编码两条制备路线之间，"
    "选择依据应包括靶点切换频率、所需生产规模和分子均一性要求。"
    "需要频繁筛选靶点的发现阶段更适合KineTAC或iLYTAC等基因编码平台；"
    "需要精细调节配体价态和连接子长度的优化阶段"
    "则更适合化学偶联平台"
    + cite(3, 13, 14)
    + "。"
    "Mamun等指出，"
    "eTPD降解剂的多样性要求建立统一的体外活性评价和体内验证标准，"
    "以促进不同平台之间的可比性"
    + cite(16)
    + "。",
    "从靶向臂角度，"
    "抗体靶向臂亲和力强但分子量大，"
    "小分子和肽段靶向臂渗透性较好但亲和力有限，"
    "适体靶向臂则兼具分子量适中和可化学修饰的优势。"
    "Wu等的Apt-LYTAC和Li等的共价适体LYTAC分别从非共价和共价两条路径"
    "探索了适体靶向臂的应用边界"
    + cite(4, 6)
    + "。"
    "Zhang等的iLYTAC和Huang等的EndoTag则证明，"
    "基因编码结合蛋白可作为靶向臂与受体臂的通用连接单元，"
    "显著缩短靶点切换的周期"
    + cite(12, 13)
    + "。",
    "连接子模块的优化需与受体臂类型联用："
    "化学偶联路线中，"
    "连接子长度影响靶向臂与受体臂的空间几何关系，"
    "偶联位点影响分子构象和体内稳定性；"
    "基因编码路线中，"
    "连接子通常为柔性肽链，"
    "其长度和氨基酸组成影响融合蛋白的表达和折叠。"
    "Fang等在SignalTAC中通过优化靶蛋白结合模块与内吞信号肽之间的肽链接头，"
    "实现了膜蛋白的高效降解，"
    "说明连接子优化在基因编码平台中同样不可或缺"
    + cite(15)
    + "。",
]


BODY_SECTION_CASES = [
    "为将上述理论框架与具体研究对接，"
    "本节选取四篇发表于Nature、Nature Chemical Biology、Nature Biotechnology"
    "和Journal of the American Chemical Society的一区研究论文进行实例剖析。"
    "所选实例分别对应LYTAC概念建立、组织定向优化、蛋白平台模块化"
    "和体内药效验证四个关键节点，"
    "构成理解该领域演进逻辑的主线样本。",
    "实例一：Banik等（Nature, 2020）建立LYTAC概念并验证多靶点可降解性。"
    "该研究要解决的核心问题是：能否将受体介导的溶酶体转运机制"
    "转化为可设计的胞外蛋白清除工具。"
    "在分子设计上，作者将靶向模块（抗体或小分子）"
    "与CI-M6PR配体（多价M6P聚糖）通过连接子偶联，"
    "构成典型的“靶向臂-受体臂-连接子”三模块结构。"
    "实验结果表明，LYTAC可在数小时内显著降低EGFR、CD71、PD-L1及ApoE4等蛋白水平，"
    "溶酶体抑制剂可阻断该效应，确证溶酶体途径依赖性。"
    "尤为重要的是，PD-L1被降解后肿瘤细胞免疫逃逸能力减弱，"
    "首次将LYTAC与肿瘤免疫调控联系起来。"
    "在机制层面，作者利用CRISPR干扰筛选发现exocyst复合体参与CI-M6PR膜面呈现，"
    "说明LYTAC研究不能脱离对受体可及性和内吞机器的生物学解析。"
    "该工作的学术价值在于："
    "其一，证明溶酶体降解路线可推广至多个疾病相关靶点；"
    "其二，将药物化学设计与细胞生物学机制研究紧密结合；"
    "其三，为后续组织定向和平台化研究提供了概念起点（图2）"
    + cite(1)
    + "。",
    "实例二：Ahn等（Nature Chemical Biology, 2021）实现肝组织特异性降解。"
    "第一代CI-M6PR依赖LYTAC虽靶点覆盖广，"
    "但受体广泛表达带来全身脱靶风险。"
    "Ahn等将受体臂替换为肝特异性ASGPR的tri-GalNAc配体，"
    "构建GalNAc-LYTAC，在肝细胞中实现EGFR高效降解。"
    "与抑制性抗体相比，GalNAc-LYTAC对EGFR下游信号抑制更持久，"
    "表明降解策略在持续时间上可能优于单纯阻断。"
    "研究还将肽段作为靶向臂降解整合素，证明靶向模块不必局限于完整抗体；"
    "通过在抗体Fc区定点偶联tri-GalNAc，"
    "改善了分子的体内药代动力学表现。"
    "该实例的启示在于："
    "受体臂的更换可直接改变降解发生的组织部位，"
    "是降低脱靶效应的首要工程手段；"
    "连接子位点和偶联化学对体内行为具有决定性影响（图3）"
    + cite(2)
    + "。",
    "实例三：Pance等（Nature Biotechnology, 2023）提出模块化KineTAC蛋白平台。"
    "化学合成聚糖-蛋白LYTAC在制备均一性和放大生产方面面临瓶颈。"
    "Pance等以双特异性抗体为骨架，"
    "一条臂结合细胞因子受体（如CXCR7），另一条臂结合靶蛋白，"
    "利用CXCL12等配体触发受体内吞，"
    "将靶蛋白导向溶酶体降解。"
    "该平台针对PD-L1、HER2、PD-1、EGFR、CDCP1和TROP2等靶点构建KineTAC，"
    "降解效率介于51%至93%之间，"
    "且可通过哺乳动物细胞表达生产，"
    "更换靶蛋白结合臂无需重新进行聚糖化学合成。"
    "从三模块视角看，KineTAC将受体臂功能嵌入细胞因子受体介导的内吞过程，"
    "将靶向臂设计为可替换的抗体片段，"
    "充分体现了LYTAC的平台化属性。"
    "该研究的贡献在于为化学LYTAC提供了基因编码替代路线，"
    "并证明受体臂不必局限于CI-M6PR或ASGPR等传统LTR配体（图4）"
    + cite(3)
    + "。",
    "实例四：Li等（Journal of the American Chemical Society, 2023）"
    "实现共价适体LYTAC的体内抗肿瘤验证。"
    "PD-L1是免疫检查点治疗的重要靶点，"
    "临床以抗体阻断为主，"
    "但高表达靶蛋白可持续补充膜面水平。"
    "Li等设计DNA适体靶向PD-L1，"
    "通过生物正交反应增强适体与靶蛋白的共价结合，"
    "再偶联溶酶体靶向模块实现PD-L1降解。"
    "在肿瘤动物模型中，"
    "共价适体LYTAC不仅显著降低PD-L1水平，"
    "还诱导了强于抗PD-L1抗体的T细胞浸润和细胞因子释放，"
    "直接提供了“降解优于阻断”的体内证据。"
    "该实例表明："
    "适体可作为轻量化的靶向臂，"
    "在简化分子构建的同时实现强效降解；"
    "体内药效比较是评价LYTAC临床转化潜力的关键环节，"
    "不能仅依赖体外Western blot数据。"
    "这一工作将LYTAC研究从细胞水平推进至动物药效验证，"
    "为肿瘤免疫方向的后续开发提供了重要参照"
    + cite(4)
    + "。",
    "从实验设计角度看，"
    "Banik等建立了溶酶体抑制剂、受体配体竞争和遗传学筛选三位一体的机制验证范式，"
    "这一范式已成为后续LYTAC研究的参照标准。"
    "Ahn等则首次在分子设计中引入组织定向思维，"
    "将“受体臂选择”从化学优化问题上升为药物定位策略问题。"
    "Pance等和Li等分别解决了可放大生产和体内药效验证两个产业化前置环节，"
    "使LYTAC研究从“能否降解”推进至“能否成药”"
    + cite(2, 3, 4)
    + "。",
    "综合四个实例可见，"
    "LYTAC领域的发展并非简单的分子结构微调，"
    "而是围绕“受体选择—模块组合—体内验证”三条主线递进："
    "Banik等解决“能否降解”，"
    "Ahn等解决“在哪里降解”，"
    "Pance等解决“如何模块化生产”，"
    "Li等解决“降解是否优于阻断”。"
    "四者共同构成了评价新LYTAC研究价值的基本坐标系。",
    "值得注意的是，"
    "四个实例分别采用了不同的受体臂策略："
    "Banik等使用CI-M6PR/M6P聚糖，"
    "Ahn等使用ASGPR/GalNAc，"
    "Pance等使用细胞因子受体内吞触发，"
    "Li等则通过适体-受体偶联物实现溶酶体导向。"
    "这一差异本身说明，"
    "不存在适用于所有疾病场景的“最优LYTAC”，"
    "只存在与特定靶点-组织组合最匹配的平台选择。",
    "此外，"
    "Li等实例的独特价值还在于将生物正交化学引入适体LYTAC设计，"
    "通过共价捕获增强靶蛋白结合在复杂体内环境中的稳定性，"
    "这一策略对于更新快、表达波动大的免疫检查点蛋白尤为重要。"
    "该工作提示，"
    "靶向臂工程化（如共价结合、亲和力成熟）"
    "与受体臂选择同等重要，"
    "二者共同决定体内降解效率的下限"
    + cite(4)
    + "。",
    "本节小结："
    "一区代表性研究的核心贡献分别体现在概念验证、组织定向、平台模块化和体内药效四个方面。"
    "后续新分子研究应明确自身定位于上述坐标系的哪一环节，"
    "避免在已充分验证的方向上重复低增量工作。",
    "上述四个一区实例的共同特征是："
    "均将分子构建与严格的机制验证相结合，"
    "均围绕明确的科学问题而非单纯追求更多靶点数据，"
    "均为后续研究提供了可复制的实验范式。"
    "这一经验对LYTAC领域的健康发展具有方法论意义。",
    "上述实例所建立的实验范式，"
    "包括溶酶体途径确证、组织定向验证和体内药效比较，"
    "可作为评价后续LYTAC研究质量的参照基准。",
]

BODY_SECTION3 = [
    "在四个一区实例奠定方法学基础之后，"
    "本节从疾病场景出发，"
    "系统评估降解策略相对单纯阻断的增益证据及其局限。",
    "评价的核心指标包括："
    "靶蛋白清除幅度、下游信号抑制持续时间、"
    "免疫功能激活程度及全身脱靶风险。"
    "不同疾病场景中上述指标的权重并不相同，"
    "需分别加以分析"
    + cite(7, 8)
    + "。",
    "在肿瘤免疫方向，PD-L1是最具代表性的比较靶点。"
    "Banik等首先证明LYTAC可降解PD-L1并减弱肿瘤细胞免疫逃逸"
    + cite(1)
    + "。"
    "Li等的共价适体LYTAC在肿瘤模型中不仅降解PD-L1，"
    "还诱导了强于抗PD-L1抗体的T细胞浸润和细胞因子释放，"
    "直接回答了“降解优于阻断”的问题"
    + cite(4)
    + "。"
    "Huang等的EndoTag融合PD-L1抗体在小鼠肿瘤模型中的疗效亦显著优于抗体单药"
    + cite(12)
    + "；"
    "Kim等的LYTACgyM6pG降解PD-L1后，增强T细胞杀伤肿瘤细胞的效率高于纳米抗体单用"
    + cite(14)
    + "。"
    "上述结果支持：在PD-L1靶点上，降解策略在激活抗肿瘤免疫方面可能具有超越阻断的增益；"
    "但该结论目前仍限于动物模型，人体验证尚缺。",
    "在肝靶向方向，Ahn等比较了GalNAc-LYTAC与抑制性抗体对EGFR的效应："
    "降解策略对EGFR下游信号的抑制更持久，"
    "且以肽段为靶向臂亦可降解整合素并抑制癌细胞增殖"
    + cite(2)
    + "。"
    "Wu等的Apt-LYTAC在肝细胞中实现了PDGF和PTK7的降解，"
    "证明适体靶向臂与肝定向受体臂的组合具有可行性"
    + cite(6)
    + "。"
    "肝靶向场景下，降解的增益主要体现在局部靶蛋白清除和全身脱靶风险降低，"
    "而非与阻断剂的头对头药效比较。",
    "在神经退行性疾病方向，Liu等设计了病灶激活LYTAC前体："
    "利用Aβ沉积区铜离子催化点击反应原位生成活性分子，"
    "经CD206受体介导Aβ纤维溶酶体降解"
    + cite(5)
    + "。"
    "Banik等亦验证了ApoE4等神经相关分泌蛋白的可降解性"
    + cite(1)
    + "。"
    "该领域尚缺乏降解与阻断的直接对照实验，"
    "但病灶激活策略本身旨在解决全身给药脱靶这一阻断策略难以回避的问题。",
    "跨疾病比较可见，降解优于阻断的证据在PD-L1肿瘤免疫中最为充分；"
    "肝靶向场景的优势主要体现在组织选择性而非药效幅度；"
    "神经退行性疾病则尚处于概念验证阶段。"
    "Pance等的KineTAC在PD-L1、HER2等多个靶点上实现了51%—93%的降解效率，"
    "但尚未系统比较其与阻断抗体的功能差异"
    + cite(3)
    + "。",
    "在生长因子受体方向，Banik等证明LYTAC可降解EGFR并削弱下游增殖信号，"
    "溶酶体抑制剂可逆转该效应，确证了降解途径的特异性"
    + cite(1)
    + "。"
    "Ahn等进一步表明肝靶向GalNAc-LYTAC对EGFR下游信号的抑制较抗体更持久"
    + cite(2)
    + "。"
    "然而，EGFR降解是否能在临床耐药场景中优于酪氨酸激酶抑制剂或单抗，"
    "仍需在更接近临床的模型中验证。",
    "Zhang等的iLYTAC在异种移植瘤模型中验证了抗肿瘤活性，"
    "表明基因编码平台兼具体外降解和体内药效验证的潜力"
    + cite(13)
    + "。"
    "对于更新速率较快的膜蛋白靶点，"
    "即使降解效率较高，"
    "仍须评估靶蛋白重新合成是否会抵消降解效应——"
    "这是降解策略相对阻断策略的潜在劣势之一，"
    "需在长期给药实验中加以考察"
    + cite(7)
    + "。",
    "在自身免疫与炎症方向，"
    "胞外细胞因子和免疫球蛋白等靶点理论上亦适合LYTAC降解，"
    "但公开文献中系统的体内验证仍较少。"
    "Chen等指出细胞因子受体和生长因子受体属于胞外降解值得关注的对象，"
    "该方向可能是LYTAC拓展适应症的重要突破口"
    + cite(8)
    + "。",
    "在药效学评价方法上，"
    "建议将“蛋白清除动力学”与“功能终点”联合考察："
    "前者反映降解机制是否有效，"
    "后者反映降解是否产生预期的生物学后果。"
    "仅报告某一时间点的蛋白残留百分比，"
    "不足以支持降解优于阻断的判断。",
    "Fang等的SignalTAC在膜蛋白降解中显示出优于亲本抗体的抗肿瘤活性，"
    "为“降解优于阻断”提供了膜蛋白靶点方向的补充证据"
    + cite(15)
    + "。"
    "该研究与Li等的PD-L1实例形成互补："
    "前者通过内吞信号肽绕过外源受体配体，"
    "后者通过共价适体增强靶向臂结合，"
    "从不同工程路径达成了相同的药理学目标——"
    "即彻底清除靶蛋白而非仅阻断其功能。",
    "本节小结：现有动物实验支持PD-L1降解在免疫激活上优于单纯阻断；"
    "肝靶向和神经退行性疾病方向的证据侧重于可行性与局部选择性，"
    "尚不足以普遍外推“降解优于阻断”的结论。"
    "后续研究需在各疾病场景中建立降解剂与阻断剂的标准化头对头比较。",
]

BODY_SECTION4 = [
    "LYTAC从概念验证走向临床转化，面临制备、选择性、药代和监管证据四重瓶颈。",
    "制备方面，化学合成多价M6P聚糖偶联物存在步骤繁琐、批次均一性差和放大困难等问题。"
    "Kim等的酵母聚糖路线和Zhang等、Huang等的基因编码平台"
    "分别代表了简化聚糖制备和绕开化学偶联的两种解决思路"
    + cite(14, 13, 12)
    + "，"
    "但各路线仍需建立与生物药或化药相对应的GMP质控标准。",
    "组织选择性方面，CI-M6PR广谱表达带来的脱靶降解风险尚未完全解决；"
    "ASGPR和EndoTag等策略可改善肝或特定组织定向，"
    "但肺、肾、脑等器官的精准递送仍缺乏成熟方案"
    + cite(2, 12, 16)
    + "。",
    "药代动力学方面，LYTAC与靶蛋白及受体形成三元复合物后常被一同送入溶酶体降解，"
    "体内有效浓度维持时间较短；"
    "高更新速率或高表达靶点可能需要更频繁给药"
    + cite(7, 16)
    + "。"
    "含聚糖或外源蛋白的偶联物还存在潜在免疫原性，需在反复给药动物实验中系统评估。",
    "监管与临床证据方面，Zhong等和Mamun等的综述均指出，"
    "LYTAC及相关eTPD技术仍主要处于临床前阶段，"
    "尚无正式发表的人体临床试验结果"
    + cite(16, 17)
    + "。"
    "Zhao等强调，新兴降解平台需完成充分的机制验证、毒理评估和制剂研究方可进入临床"
    + cite(18)
    + "。"
    "与PROTAC相比，LYTAC的临床转化路径更缺少可参照的审评先例"
    + cite(10, 17)
    + "。",
    "适体靶向臂虽具合成简便优势，"
    "但体内可能面临核酸酶降解和肾脏快速清除，"
    "需要通过化学修饰或递送系统改善药代性质"
    + cite(4, 6)
    + "。"
    "Fang等的SignalTAC和Liu等的病灶激活前体表明，"
    "疾病病理特征本身可作为药物设计输入，"
    "为克服组织选择性不足提供了非常规思路"
    + cite(5, 15)
    + "。",
    "与PROTAC产业化路径相比，"
    "Bekes等指出降解药物从概念验证到临床开发需经历严格的靶点拓展和工艺固化过程"
    + cite(10)
    + "。"
    "LYTAC目前缺乏类似PROTAC的临床概念验证案例，"
    "其在监管审评中可能面临更大的机制与安全性解释负担。",
    "面向后续研发，优先级较高的方向包括："
    "建立受体-靶点-组织三联匹配的标准决策流程；"
    "开发均一性更好的模块化制备工艺；"
    "在关键疾病模型中完成降解剂对阻断剂的头对头比较；"
    "以及开展系统的体内毒理和药代研究。",
    "从更宏观的TPD格局看，"
    "Zhong等指出降解药物研发需同步解决理性设计、联合用药和耐药克服等问题"
    + cite(17)
    + "；"
    "Zhao等则强调溶酶体途径与蛋白酶体途径在靶点覆盖面上互补，"
    "共同构成TPD技术版图"
    + cite(18)
    + "。"
    "LYTAC作为溶酶体途径的代表，"
    "其后续发展取决于能否在工艺、选择性、药效学和安全性四条线上同步推进，"
    "而非继续在单一细胞系中追求更高的体外降解百分比。",
    "从监管科学角度，"
    "LYTAC作为新型大分子偶联物或基因编码蛋白，"
    "其质量属性定义、杂质控制策略和免疫原性评价方法"
    "均尚未形成行业共识"
    + cite(16, 17)
    + "。"
    "在缺乏临床先例的情况下，"
    "临床前研究的完整性和可重复性将成为监管审评的关键依据。",
    "此外，"
    "Huang等揭示的内源性配体竞争、"
    "Fang等开发的内吞信号肽策略以及Liu等的病灶激活前体，"
    "分别从生物学约束和工程应对两个维度"
    "为克服现有瓶颈提供了新思路"
    + cite(5, 12, 15)
    + "。",
    "在成药性评价维度，"
    "还需关注制剂稳定性、给药途径选择和联合用药策略。"
    "Zhong等指出，"
    "TPD药物与化疗、免疫治疗或靶向药物的联合应用"
    "可能产生协同效应，但也增加了安全性评价的复杂性"
    + cite(17)
    + "。"
    "LYTAC若进入临床开发，"
    "其联合用药方案的设计需充分考虑脱靶降解对正常组织蛋白稳态的影响。",
    "Zhao等从TPD技术版图角度指出，"
    "溶酶体途径与蛋白酶体途径在靶点类型上互补，"
    "LYTAC的独特定位在于覆盖分泌蛋白、膜蛋白和蛋白聚集体等"
    "传统小分子和PROTAC均难以充分干预的靶点类别"
    + cite(18)
    + "。"
    "因此，"
    "评价LYTAC的研究价值不应简单类比PROTAC的临床进度，"
    "而应基于其是否解决了胞外靶点降解这一特定科学问题。",
    "Lin等的综述亦强调，"
    "胞外蛋白降解技术的成功开发需要药物化学、"
    "细胞生物学和临床药理学等多学科团队的持续协作，"
    "单一学科的优化难以解决LYTAC成药的全链条问题"
    + cite(11)
    + "。",
    "本节小结：LYTAC的核心科学问题已基本明确，"
    "成药瓶颈主要集中在工艺可放大性、组织选择性和临床前证据体系三方面。"
    "在上述瓶颈未解决前，不宜将LYTAC视为接近临床的成熟技术平台。",
]

CONCLUSION = [
    "LYTAC技术将靶向蛋白降解从胞内蛋白酶体途径拓展至胞外及膜蛋白，"
    "其重要意义在于为约占蛋白质组相当比例的胞外靶点提供了可操作的干预策略，"
    "有望突破传统抗体和小分子抑制剂在“难以彻底清除致病蛋白”方面的局限"
    + cite(1, 8, 18)
    + "。",
    "本综述的梳理表明，"
    "LYTAC药物研发并非单纯的化学偶联问题，"
    "而是受体生物学、模块化分子工程与疾病药理学三者耦合的系统工程。"
    "受体选择决定降解的组织定位与通量，"
    "三模块平台设计决定靶点切换效率与制备可行性，"
    "而体内降解与阻断的头对头比较才是评价临床转化潜力的最终标尺。"
    "现有证据支持PD-L1降解在动物模型中可产生优于抗体阻断的免疫激活效应，"
    "但肝靶向和神经退行性疾病方向仍缺乏充分的对比数据"
    + cite(2, 4, 5)
    + "。",
    "从研究发展规律看，"
    "LYTAC正经历PROTAC早期所经历的“概念验证—平台分化—临床前聚焦”阶段。"
    "Bekes等和Zhong等关于PROTAC临床转化的回顾提示，"
    "降解技术从实验室走向临床需同步解决机制确证、工艺放大和监管路径三大问题"
    + cite(10, 17)
    + "。"
    "LYTAC目前尚无任何人体临床试验报道，"
    "在可放大制备、组织选择性、药代动力学和长期安全性方面仍存在明显短板"
    + cite(16)
    + "。",
    "基于上述分析，本综述认为后续研究应优先聚焦以下方向："
    "第一，建立“靶点-受体-组织”三联匹配的标准决策流程，"
    "避免盲目追求体外降解率；"
    "第二，推动化学聚糖路线与基因编码路线的工艺标准化，"
    "解决批次均一性和放大生产问题；"
    "第三，在PD-L1、EGFR等关键靶点上系统开展降解剂对阻断剂的头对头动物实验；"
    "第四，探索病灶激活、内吞信号肽融合等非常规策略，"
    "以克服全身给药脱靶和受体竞争等生物学约束"
    + cite(5, 12, 14, 15)
    + "。",
    "总体而言，LYTAC代表了靶向蛋白降解领域的重要拓展方向。"
    "其长期价值不仅在于提供新的分子工具，"
    "更在于推动药物研发从“功能抑制”向“蛋白清除”的理念转变。"
    "随着受体生物学认识的深入和模块化制备技术的成熟，"
    "LYTAC有望在肿瘤免疫、肝靶向代谢病和神经退行性疾病等领域"
    "成为具有独特竞争优势的候选药物策略，"
    "但其成药前景仍有赖于临床前数据的系统积累与严格验证。",
    "对本领域的深入研读亦提示："
    "LYTAC研究正在从“证明新概念”转向“解决可开发性问题”，"
    "这一转变要求研究者同时具备细胞生物学、药物化学和疾病药理学的交叉视野。"
    "只有将受体决策、模块工程和体内验证三条主线贯通，"
    "才能推动LYTAC从学术热点成长为具有临床竞争力的药物技术。",
]
REFS = [
    "[1] Banik S M, Pedram K, Wisnovsky S, Riley N M, Bertozzi C R. Lysosome-targeting chimaeras for degradation of extracellular proteins[J]. Nature, 2020, 584(7820): 291-297. https://doi.org/10.1038/s41586-020-2545-9",
    "[2] Ahn G, Banik S M, Miller C L, Riley N M, Cochran J R, Bertozzi C R. LYTACs that engage the asialoglycoprotein receptor for targeted protein degradation[J]. Nature Chemical Biology, 2021, 17(9): 937-946. https://doi.org/10.1038/s41589-021-00770-1",
    "[3] Pance K, Gramespacher J A, Byrnes J R, et al. Modular cytokine receptor-targeting chimeras for targeted degradation of cell surface and extracellular proteins[J]. Nature Biotechnology, 2023, 41(2): 273-281. https://doi.org/10.1038/s41587-022-01456-2",
    "[4] Li Y, Liu X, Yu L, Huang X, Wang X, Han D, et al. Covalent LYTAC enabled by DNA aptamers for immune checkpoint degradation therapy[J]. Journal of the American Chemical Society, 2023, 145(45): 24506-24521. https://doi.org/10.1021/jacs.3c03899",
    "[5] Liu Z, Deng Q, Qin G, Yang J, Zhang H, Ren J, Qu X. Biomarker-activated multifunctional lysosome-targeting chimeras mediated selective degradation of extracellular amyloid fibrils[J]. Chem, 2023, 9(7): 2016-2038. https://doi.org/10.1016/j.chempr.2023.06.003",
    "[6] Wu Y, Lu Y, Li L, Deng K, Zhang S, Yang C, Zhu Z. Aptamer-LYTACs for targeted degradation of extracellular and membrane proteins[J]. Angewandte Chemie International Edition, 2023, 62(15): e202218106. https://doi.org/10.1002/anie.202218106",
    "[7] Li Y Y, Yang Y, Zhang R S, Ge R X, Xie S B. Targeted degradation of membrane and extracellular proteins with LYTACs[J]. Acta Pharmacologica Sinica, 2025, 46: 1-7. https://doi.org/10.1038/s41401-024-01364-y",
    "[8] Chen X, Zhou Y, Zhao Y, Tang W. Targeted degradation of extracellular secreted and membrane proteins[J]. Trends in Pharmacological Sciences, 2023, 44(11): 762-775. https://doi.org/10.1016/j.tips.2023.08.013",
    "[9] Caianiello D F, Zhang M, Ray J D, et al. Bifunctional small molecules that mediate the degradation of extracellular proteins[J]. Nature Chemical Biology, 2021, 17(9): 947-953. https://doi.org/10.1038/s41589-021-00851-1",
    "[10] Bekes M, Langley D R, Crews C M. PROTAC targeted protein degraders: the past is prologue[J]. Nature Reviews Drug Discovery, 2022, 21(3): 181-200. https://doi.org/10.1038/s41573-021-00371-6",
    "[11] Lin J, Jin J, Shen Y, Zhang L, Gong G, Bian H, Chen H, Nagle D G, Wu Y, Zhang W, Luan X. Emerging protein degradation strategies: expanding the scope to extracellular and membrane proteins[J]. Theranostics, 2021, 11(17): 8337-8349. https://doi.org/10.7150/thno.62686",
    "[12] Huang B, Abedi M, Coventry B, et al. Designed endocytosis-inducing proteins degrade targets and amplify signals[J]. Nature, 2024, 632(8024): 191-200. https://doi.org/10.1038/s41586-024-07948-2",
    "[13] Zhang B, Brahma R K, Zhu L, Feng J, Hu S, Qian L, et al. Insulin-like Growth Factor 2 (IGF2)-Fused Lysosomal Targeting Chimeras for Degradation of Extracellular and Membrane Proteins[J]. Journal of the American Chemical Society, 2023, 145(44): 24272-24283. https://doi.org/10.1021/jacs.3c08886",
    "[14] Kim S, Kang J Y, Bi A D, Seo J, Oh D B. Lysosome-Targeting Chimera Using Mannose-6-Phosphate Glycans Derived from Glyco-Engineered Yeast[J]. Bioconjugate Chemistry, 2025, 36(3): 424-436. https://doi.org/10.1021/acs.bioconjchem.4c00512",
    "[15] Fang T, Zheng Z, Li N, Zhang Y, Ma J, Yun C, Cai X. Lysosome-targeting chimeras containing an endocytic signaling motif trigger endocytosis and lysosomal degradation of cell-surface proteins[J]. Chemical Science, 2024, 15(42): 17652-17662. https://doi.org/10.1039/d4sc05093b",
    "[16] Mamun A A, Uzunparmak B, Crews C M. Targeted degradation of extracellular proteins: state of the art and diversity of degrader designs[J]. Journal of Hematology & Oncology, 2025, 18: 23. https://doi.org/10.1186/s13045-025-01703-4",
    "[17] Zhong G, Chang X, Xie W, Zhou X. Targeted protein degradation: advances in drug discovery and clinical practice[J]. Signal Transduction and Targeted Therapy, 2024, 9: 308. https://doi.org/10.1038/s41392-024-02004-x",
    "[18] Zhao L, Zhao J, Zhong K, Tong A, Jia D. Targeted protein degradation: mechanisms, strategies and application[J]. Signal Transduction and Targeted Therapy, 2022, 7: 113. https://doi.org/10.1038/s41392-022-00966-4",
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

    add_heading_text(doc, "一、溶酶体靶向受体的生物学基础与受体选择决策", level=3)
    for para in BODY_SECTION1:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图1  内体分选决定LYTAC-受体复合物进入溶酶体降解或回收至细胞膜，"
        "是降解效率的核心分叉点。"
        "（引自：Banik S M, et al. Nature, 2020, 584(7820): 291-297 [1]）",
    )

    add_heading_text(doc, "二、LYTAC的模块化分子工程设计", level=3)
    for para in BODY_SECTION2:
        add_paragraph(doc, para)

    add_heading_text(doc, "三、一区代表性研究实例剖析", level=3)
    for para in BODY_SECTION_CASES:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图2  Banik等证实CI-M6PR依赖LYTAC可同步降解EGFR、PD-L1等多种靶蛋白，"
        "奠定LYTAC概念框架。"
        "（引自：Banik S M, et al. Nature, 2020, 584(7820): 291-297 [1]）"
        "图3  Ahn等证明GalNAc-LYTAC可实现肝组织定向降解，"
        "且对EGFR下游信号抑制优于阻断性抗体。"
        "（引自：Ahn G, et al. Nat Chem Biol, 2021, 17(9): 937-946 [2]）"
        "图4  Pance等证实KineTAC可通过更换靶向臂实现多靶点降解，"
        "体现受体臂与靶向臂的模块化可替换性。"
        "（引自：Pance K, et al. Nat Biotechnol, 2023, 41(2): 273-281 [3]）"
        "图5  Li等证实共价适体LYTAC在动物模型中降解PD-L1并激活抗肿瘤免疫，"
        "效应强于抗PD-L1抗体。"
        "（引自：Li Y, et al. J Am Chem Soc, 2023, 145(45): 24506-24521 [4]）",
    )

    add_heading_text(doc, "四、疾病模型中的应用验证：降解优于阻断的证据与局限", level=3)
    for para in BODY_SECTION3:
        add_paragraph(doc, para)

    add_heading_text(doc, "五、转化瓶颈与发展展望", level=3)
    for para in BODY_SECTION4:
        add_paragraph(doc, para)

    add_heading_text(doc, "总结与展望", level=2)
    for para in CONCLUSION:
        add_paragraph(doc, para)

    add_heading_text(doc, "参考文献", level=2)
    for ref in REFS:
        add_reference(doc, ref)

    output_path = "/workspace/LYTAC技术药物研发综述.docx"
    doc.save(output_path)

    body_text = "".join(BODY_SECTION1 + BODY_SECTION2 + BODY_SECTION_CASES + BODY_SECTION3 + BODY_SECTION4)
    body_count = count_chinese_chars(body_text)
    total_text = ABSTRACT + "".join(INTRO + BODY_SECTION1 + BODY_SECTION2 + BODY_SECTION_CASES + BODY_SECTION3 + BODY_SECTION4 + CONCLUSION)
    total_count = count_chinese_chars(total_text)
    print(f"Document saved to: {output_path}")
    print(f"正文汉字数（不含文献名标注）: {body_count}")
    print(f"全文汉字数（不含文献名标注）: {total_count}")


if __name__ == "__main__":
    main()
