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
    "溶酶体靶向嵌合体（LYTAC）通过偶联靶蛋白与溶酶体靶向受体（LTR），"
    "将胞外及膜蛋白导入溶酶体降解，是PROTAC之外面向胞外靶点的重要降解策略。"
    "本综述围绕受体生物学基础、模块化分子设计、疾病应用验证和转化瓶颈四个层面，"
    "分析受体选择如何决定降解效率与组织定位，"
    "讨论靶向臂、受体臂和连接子三大模块的可替换平台属性，"
    "并评估降解策略相对单纯阻断的增益证据及成药前主要障碍。"
    + cite(1, 8, 16)
)

INTRO = [
    "靶向蛋白降解为传统上难以药化的蛋白提供了新路径。"
    "PROTAC通过招募胞内E3连接酶将靶蛋白导向蛋白酶体降解，"
    "已有候选分子进入临床评估"
    + cite(10)
    + "。"
    "然而，该策略要求靶蛋白具备可被配体结合的胞内结构域，"
    "对分泌蛋白、抗体和跨膜蛋白等胞外靶点覆盖有限。",
    "2020年，Banik等提出溶酶体靶向嵌合体（LYTAC），"
    "其核心设计是同时结合靶蛋白与细胞表面LTR，"
    "经受体介导内吞将三元复合物分选至溶酶体完成蛋白清除"
    + cite(1)
    + "。"
    "本综述依“受体生物学基础—模块化分子设计—疾病应用验证—转化瓶颈展望”"
    "四段式框架，讨论LYTAC药物研发的逻辑主线"
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
    "本节小结：LYTAC疗效的关键瓶颈不在于靶蛋白能否被结合，"
    "而在于所选LTR的组织分布、内吞通量及内体分选倾向是否与疾病场景匹配。"
    "受体选择应作为药物设计的首要决策节点，而非分子构建完成后的附带调整。",
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
]

BODY_SECTION3 = [
    "疾病应用评估的核心问题，是蛋白降解是否优于单纯阻断。"
    "抗体等阻断剂通过占据结合位点抑制靶蛋白功能，"
    "但高表达或快速更新的靶蛋白可持续补充膜面水平；"
    "降解策略则旨在降低靶蛋白总量，"
    "理论上可同时削弱其信号传导和蛋白相互作用网络"
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
    "Fang等的SignalTAC在膜蛋白降解中显示出优于亲本抗体的抗肿瘤活性，"
    "为“降解优于阻断”提供了另一类型的实验证据"
    + cite(15)
    + "。",
    "在神经退行性疾病方向，Liu等构建了可穿越血脑屏障的纳米载体，"
    "将点击反应前体递送至脑部，"
    "在Aβ沉积区原位激活LYTAC并介导纤维清除"
    + cite(5)
    + "。"
    "该设计的关键增益不在于与阻断剂比较药效，"
    "而在于利用病理微环境实现局部激活，"
    "从给药策略层面规避全身脱靶。",
    "本节小结：现有动物实验支持PD-L1降解在免疫激活上优于单纯阻断；"
    "肝靶向和神经退行性疾病方向的证据侧重于可行性与局部选择性，"
    "尚不足以普遍外推“降解优于阻断”的结论。"
    "后续研究需在各疾病场景中建立降解剂与阻断剂的标准化头对头比较。",
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
    "本节小结：LYTAC的核心科学问题已基本明确，"
    "成药瓶颈主要集中在工艺可放大性、组织选择性和临床前证据体系三方面。"
    "在上述瓶颈未解决前，不宜将LYTAC视为接近临床的成熟技术平台。",
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
]

CONCLUSION = [
    "本综述按受体生物学基础、模块化分子设计、疾病应用验证和转化瓶颈四个层面，"
    "梳理了LYTAC药物研发的逻辑主线。"
    "受体选择决定降解的组织定位与通量，"
    "三模块平台设计决定靶点切换与制备可行性，"
    "疾病模型中的增益证据目前主要集中于PD-L1降解优于阻断，"
    "而成药前仍需在工艺、选择性和临床前数据方面取得系统性突破"
    + cite(1, 8, 16)
    + "。",
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
        "是降解效率的核心分叉点。（引自[1]）",
    )

    add_heading_text(doc, "二、LYTAC的模块化分子工程设计", level=3)
    for para in BODY_SECTION2:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图2  CI-M6PR依赖LYTAC可同时降解多种膜蛋白与分泌蛋白，验证了广谱降解路线的可行性。（引自[1]）"
        "图3  ASGPR介导GalNAc-LYTAC实现肝组织定向降解，对EGFR下游信号的抑制优于阻断性抗体。（引自[2]）"
        "图4  KineTAC以细胞因子受体为内吞触发模块，证明受体臂可脱离聚糖化学合成独立替换。（引自[3]）"
        "图5  酵母来源M6P聚糖受体臂可简化LYTAC制备并实现PD-L1降解及免疫激活。（引自[14]）",
    )

    add_heading_text(doc, "三、疾病模型中的应用验证：降解优于阻断的证据与局限", level=3)
    for para in BODY_SECTION3:
        add_paragraph(doc, para)

    add_heading_text(doc, "四、转化瓶颈与发展展望", level=3)
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

    body_text = "".join(BODY_SECTION1 + BODY_SECTION2 + BODY_SECTION3 + BODY_SECTION4)
    body_count = count_chinese_chars(body_text)
    total_text = ABSTRACT + "".join(INTRO + BODY_SECTION1 + BODY_SECTION2 + BODY_SECTION3 + BODY_SECTION4 + CONCLUSION)
    total_count = count_chinese_chars(total_text)
    print(f"Document saved to: {output_path}")
    print(f"正文汉字数（不含文献名标注）: {body_count}")
    print(f"全文汉字数（不含文献名标注）: {total_count}")


if __name__ == "__main__":
    main()
