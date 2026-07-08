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
    summary_prefix = "【本节小结】"
    if text.startswith(summary_prefix):
        label_run = p.add_run(summary_prefix)
        set_run_font(label_run, size_pt=size_pt, bold=True)
        body_text = format_citations(text[len(summary_prefix):])
        body_run = p.add_run(body_text)
        set_run_font(body_run, size_pt=size_pt, bold=bold)
    else:
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
    "本综述依受体生物学基础、模块化分子设计、"
    "疾病应用验证和转化瓶颈四个层面展开论述，"
    "分析受体选择对降解效率与组织定位的决定作用，"
    "归纳靶向臂、受体臂和连接子三模块的平台属性，"
    "并在分子工程论述中嵌入Banik、Ahn、Pance和Li等代表性研究的实证节点，"
    "评估降解策略相对单纯阻断的增益及成药前主要障碍。"
    + cite(1, 2, 3, 4)
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
    "Chen等和Mamun等的综述已对胞外蛋白降解的技术版图与分子类型做过宏观梳理"
    + cite(8, 16)
    + "。"
    "本综述在简述上述背景后，"
    "依次讨论受体选择决策、三模块分子工程、疾病应用证据及转化瓶颈，"
    "以期为LYTAC药物研发提供逻辑清晰的分析框架。",
]

BODY_SECTION1 = [
    "LYTAC的降解效率首先取决于所选用溶酶体靶向受体（LTR）的生物学特性，"
    "而非仅由化学结构决定。"
    "三元复合物形成后，经网格蛋白介导内吞进入早期内体；"
    "此后复合物或被分选至溶酶体完成蛋白水解，"
    "或通过retromer等途径回收至细胞膜，两条路径的竞争直接决定有效降解通量（图1）。"
    "Lin等的早期综述已将这一分选-回收平衡视为解释不同细胞系降解效率差异的核心变量"
    + cite(11)
    + "。",
    "与抗体阻断不同，LYTAC属于事件驱动型策略："
    "结合剂无需抑制靶蛋白活性，只要有效捕获即可启动降解程序。"
    "降解发生在溶酶体后，靶蛋白的膜信号传导和胞外相互作用可同时终止，"
    "这一机制差异是后续“降解是否优于阻断”比较的生物学基础"
    + cite(7)
    + "。",
    "验证溶酶体途径通常需联合溶酶体抑制剂、受体敲除或配体竞争，"
    "以及遗传学筛选鉴定内吞分选相关基因。"
    "Banik等通过CRISPR干扰筛选发现，"
    "exocyst复合体参与调控CI-M6PR的膜面呈现，"
    "提示受体可及性本身即是降解效率的上游决定因素"
    + cite(1)
    + "。",
    "依组织分布与功能特征，"
    "可将LTR归纳为三类并据此建立受体选择决策框架。",
    "广谱型受体的代表是CI-M6PR（IGF2R）："
    "其在多种组织均有表达，内吞通量较高，"
    "适合不需要严格组织定向、且可接受一定脱靶风险的场景。"
    "Huang等通过内源性配体竞争实验揭示，"
    "细胞内M6P修饰糖蛋白可与外源配体竞争CI-M6PR结合位点，"
    "是广谱受体路线中不可忽视的生物学约束"
    + cite(12)
    + "。"
    "该发现意味着：即便分子亲和力很高，"
    "内源性受体占用仍可能限制降解通量；"
    "解读体外降解数据时需同步考虑受体竞争状态。",
    "组织定向型受体的代表是ASGPR："
    "其表达集中于肝细胞，"
    "可将内吞-降解过程限制在肝脏微环境，"
    "降解位置可控、全身脱靶风险相对较低。"
    "当治疗目标为肝内靶蛋白或需降低全身脱靶时，应优先评估ASGPR通路；"
    "当靶蛋白分布广泛或需作用于非肝组织时，则需转向其他LTR。",
    "可拓展型受体包括整合素、清道夫受体等："
    "其组织表达谱各异，"
    "相关降解策略尚处探索阶段，"
    "但为非常规组织靶向提供了潜在选项"
    + cite(11)
    + "。"
    "选择可拓展型受体时，"
    "需重点评估其在靶组织的表达丰度、内吞通量及与靶蛋白的共定位条件，"
    "避免因受体数据不完整而导致降解位置偏离预期。",
    "综合上述三类受体，决策逻辑可概括为："
    "若追求最大靶蛋白覆盖且可接受脱靶风险，优先考虑CI-M6PR体系；"
    "若需肝定向降解，选择ASGPR通路；"
    "若现有LTR难以匹配靶组织表达谱，可评估整合素、清道夫受体等可拓展型受体。",
    "Zhao等将溶酶体途径定位为拓展靶向蛋白降解（TPD）靶点范围的重要方向，"
    "指出胞外及膜蛋白约占蛋白质组的相当比例，"
    "是LYTAC存在的结构性需求"
    + cite(18)
    + "。"
    "在实际决策中，建议按以下顺序评估："
    "（1）靶蛋白的细胞外定位与更新速率；"
    "（2）靶标组织的LTR表达谱；"
    "（3）可接受的脱靶降解范围；"
    "（4）内源性配体竞争与受体回收动力学。",
    "对于膜蛋白靶点，还需关注其与LTR在同一细胞表面的空间距离和共定位关系，"
    "因为这将影响三元复合物形成的几何可行性。",
    "从疾病相关靶点分布看，"
    "PD-L1、EGFR和HER2等膜蛋白在肿瘤中表达异质性高，"
    "CI-M6PR广谱降解可能带来肿瘤外组织脱靶；"
    "Aβ等靶点则涉及分泌蛋白或聚集体的清除，"
    "对受体通量和局部给药策略提出不同要求"
    + cite(5)
    + "。"
    "因此，受体决策框架必须与具体疾病场景联用，"
    "不能脱离靶点生物学单独优化分子结构。",
    "在评价降解效率时，"
    "还需区分“受体结合亲和力”与“有效降解通量”两个层次。"
    "高亲和力仅保证三元复合物形成，"
    "而通量取决于内吞速率、溶酶体分选效率和受体回收动力学的综合平衡。"
    "同一LYTAC在不同细胞系中降解效率的差异，"
    "往往反映的是后者而非前者的问题。",
    "受体表达密度与膜周转速率亦构成重要的生物学变量："
    "在受体表达水平较低的细胞系中，"
    "即便配体亲和力很高，"
    "有效三元复合物形成的概率仍可能受限；"
    "在膜周转较快的细胞中，"
    "内吞通量较高但受体回收亦更活跃，"
    "需通过时间分辨的降解动力学实验加以区分。",
    "对于分泌蛋白靶点，"
    "胞外蛋白池的动态平衡（合成、分泌、清除）"
    "决定了降解剂是否能在稳态下维持靶蛋白低水平，"
    "这一考量与膜蛋白靶点存在本质差异"
    + cite(7)
    + "。",
    "【本节小结】LYTAC疗效的关键瓶颈不在于靶蛋白能否被结合，"
    "而在于所选LTR的组织分布、内吞通量、内体分选倾向及内源性配体竞争"
    "是否与疾病场景匹配。"
    "受体选择应作为药物设计的首要决策节点。"
    "第二节将分别介绍对应上述受体类型的分子工程实现路径。",
    "在跨物种转化时，"
    "还需注意LTR表达谱和溶酶体分选机器在人源细胞与常用细胞系之间的差异，"
    "避免将体外筛选得到的“最优受体”直接等同于临床适用方案。",
    "这一原则同样适用于解读不同细胞系间的降解效率差异。",
    "受体回收动力学与溶酶体分选效率的联立评估，"
    "是区分“结合成功”与“降解成功”的关键实验环节之一。",
]

BODY_SECTION2 = [
    "在第一节受体决策框架确定后，"
    "LYTAC分子可拆解为三个可独立替换的功能模块："
    "靶向臂（结合靶蛋白）、受体臂（触发内吞与溶酶体分选）和连接子（调控模块间距与偶联均一性）。"
    "三模块的解耦设计使LYTAC具备平台化属性——"
    "同一受体臂可搭配不同靶向臂以切换靶点，"
    "同一靶向臂亦可更换受体臂以改变组织分布"
    + cite(7)
    + "。",
    "（一）靶向臂：按结合分子类型分组",
    "靶向臂的可选范围已从抗体扩展至小分子、肽段、适体和基因编码结合蛋白。"
    "按制备与偶联方式，可归纳为四条主线。",
    "抗体靶向臂亲和力强、可及性成熟，是多数早期平台的默认选择。"
    "Banik等在初代LYTAC中分别使用抗体和小分子结合EGFR等靶点，"
    "数小时内显著降低EGFR、CD71、PD-L1及ApoE4等蛋白水平，"
    "溶酶体抑制剂可阻断该效应，"
    "首次将“靶向臂-受体臂-连接子”三模块概念与多靶点可降解性一并验证"
    + cite(1)
    + "。",
    "小分子与肽段靶向臂分子量较低、组织渗透性较好。"
    "Caianiello等的MoDE-A以双功能小分子同时结合靶蛋白和ASGPR，"
    "实现了整合素等靶点的溶酶体降解，分子量显著低于抗体偶联物"
    + cite(9)
    + "。"
    "Ahn等以肽段结合剂降解整合素，"
    "表明靶向臂不必局限于完整抗体"
    + cite(2)
    + "。",
    "适体靶向臂兼具合成简便和可化学修饰的优势。"
    "Wu等以适体连接tri-GalNAc实现肝细胞中PDGF和PTK7的降解"
    + cite(6)
    + "；"
    "Li等设计DNA适体靶向PD-L1，"
    "通过生物正交反应增强适体与靶蛋白的共价结合，"
    "再偶联溶酶体靶向模块，"
    "在肿瘤动物模型中诱导了强于抗PD-L1抗体的T细胞浸润和细胞因子释放，"
    "为“降解优于阻断”提供了体内证据"
    + cite(4)
    + "。",
    "基因编码结合蛋白（affibody、纳米抗体、IgG结合Z结构域等）"
    "可与受体臂融合表达，缩短靶点切换周期。"
    "Zhang等的iLYTAC以多种结合蛋白为靶向臂、IGF2为受体臂，"
    "覆盖EGFR、PD-L1、CD20和α-突触核蛋白等靶点"
    + cite(13)
    + "；"
    "Huang等的EndoTag融合不同靶蛋白结合蛋白后，"
    "可定向至IGF2R或ASGPR等不同受体"
    + cite(12)
    + "。",
    "（二）受体臂：按制备路线分组",
    "受体臂的工程实现可按“化学偶联”与“基因编码”两条主线组织，"
    "并以前述第一节的受体生物学依据为选型前提。",
    "基于CI-M6PR广谱受体（第一节已述其组织分布与内源性配体竞争约束），"
    "化学偶联路线以多价M6P聚糖为典型受体臂。"
    "Banik等采用化学合成M6P聚糖与靶向模块偶联，"
    "奠定了聚糖-蛋白LYTAC的基本架构"
    + cite(1)
    + "。"
    "Kim等以糖工程酵母发酵制备人源兼容的gyM6pG聚糖，"
    "通过无铜点击化学与纳米抗体偶联，"
    "在沿用CI-M6PR受体臂的同时降低了聚糖合成门槛"
    + cite(14)
    + "。",
    "基于ASGPR肝靶向受体（第一节已述其组织定向优势），"
    "化学偶联路线主要采用tri-GalNAc配体与靶向臂偶联。"
    "Ahn等构建GalNAc-LYTAC，在肝细胞中实现EGFR高效降解；"
    "与抑制性抗体相比，对EGFR下游信号抑制更持久，"
    "表明受体臂更换可直接改变降解发生的组织部位"
    + cite(2)
    + "。"
    "Wu等的Apt-LYTAC将适体靶向臂与tri-GalNAc受体臂组合，"
    "进一步验证了肝定向化学偶联平台的模块可替换性"
    + cite(6)
    + "。",
    "基因编码受体臂绕开复杂聚糖化学，"
    "适合需要频繁切换靶点的发现阶段。"
    "其中一条路线是利用细胞因子受体等天然高周转膜蛋白触发靶蛋白共内化："
    "Pance等以双特异性蛋白为骨架，"
    "一条臂结合细胞因子受体（如CXCR7），另一条臂结合靶蛋白，"
    "利用CXCL12等配体触发受体内吞，"
    "针对PD-L1、HER2、PD-1、EGFR、CDCP1和TROP2等靶点构建KineTAC，"
    "降解效率介于51%至93%之间，"
    "且可通过哺乳动物细胞表达生产，"
    "更换靶蛋白结合臂无需重新进行聚糖化学合成"
    + cite(3)
    + "。"
    "该策略将受体臂功能嵌入天然内吞过程，"
    "但需评估配体占用对受体生理信号传导的干扰。",
    "另一条基因编码路线仍以经典LTR为终点："
    "Zhang等的iLYTAC以IGF2多肽结合CI-M6PR，"
    "Huang等的EndoTag以人工设计蛋白结合IGF2R或ASGPR，"
    "二者均避免了聚糖化学修饰"
    + cite(12, 13)
    + "。",
    "此外，还存在绕过外源LTR配体、直接利用内吞分选信号的策略："
    "Fang等的SignalTAC将CI-M6PR来源的酪氨酸基内吞信号肽P3"
    "与靶蛋白结合模块融合，"
    "在不依赖外源聚糖配体的情况下诱导靶蛋白溶酶体降解，"
    "在膜蛋白靶向降解中显示出优于亲本抗体的活性"
    + cite(15)
    + "。"
    "该路线将受体臂功能转化为内源性分选信号，"
    "组织定向性取决于内吞信号来源与靶细胞分选机器的差异。",
    "（三）连接子与偶联化学",
    "连接子模块决定两个功能臂的空间排布与偶联均一性。"
    "Ahn等在抗体Fc区定点偶联tri-GalNAc，改善了GalNAc-LYTAC的体内稳定性"
    + cite(2)
    + "。"
    "Li等在综述中指出，连接子过长不利于三元复合物形成，过短则限制构象调整；"
    "配体价态过高亦可能引发非特异性聚集"
    + cite(7)
    + "。"
    "Fang等在SignalTAC中通过优化靶蛋白结合模块与内吞信号肽之间的肽链接头，"
    "实现了膜蛋白的高效降解，"
    "说明连接子优化在基因编码平台中同样不可或缺"
    + cite(15)
    + "。"
    "化学偶联路线中，偶联位点和偶联数不一致会直接影响降解活性重现性；"
    "基因编码路线中，柔性肽链的长度和氨基酸组成影响融合蛋白的表达与折叠。",
    "（四）平台组合与制备路线选择",
    "按三模块组合逻辑，现有LYTAC变体可归纳为若干典型平台："
    "聚糖-抗体平台（Banik、Kim）、GalNAc-肝靶向平台（Ahn、Wu）、"
    "适体-共价平台（Li）、双特异性蛋白平台（Pance、Huang、Zhang）"
    "和小分子-受体平台（Caianiello）"
    + cite(1, 2, 4, 6, 9, 12, 13, 14)
    + "。"
    "制备路线可按受体臂类型分为三条主线："
    "化学偶联聚糖臂（Banik、Kim）、化学合成小分子臂（Caianiello）"
    "和基因编码蛋白臂（Pance、Huang、Zhang、Fang）。"
    "前者的优势在于结构可精细调控，后者的优势在于批次均一性和靶点切换速度。",
    "各平台在制备难度、靶点切换灵活性和组织定向能力上各有权衡，"
    "选择时应与第一节确定的受体类型和疾病场景相匹配，"
    "而非单纯追求体外降解率最高。",
    "从平台可替换性看，"
    "Pance等的KineTAC可通过更换靶蛋白结合臂降解多个靶点而保留同一受体臂"
    + cite(3)
    + "；"
    "Zhang等的iLYTAC型-II以IgG结合Z结构域为通用适配器，"
    "可与临床已有抗体偶联实现降解"
    + cite(13)
    + "。"
    "这些设计表明，LYTAC研发正从“逐个分子定制”转向“平台化模块组装”。",
    "选择化学合成与基因编码路线时，"
    "应综合靶点切换频率、生产规模和分子均一性要求："
    "发现阶段更适合KineTAC或iLYTAC等基因编码平台；"
    "需要精细调节配体价态和连接子长度的优化阶段"
    "则更适合化学偶联平台"
    + cite(3, 13, 14)
    + "。",
    "从代表性研究的演进逻辑看，"
    "Banik等、Ahn等、Pance等和Li等分别对应"
    "“能否降解”“在哪里降解”“如何模块化生产”和“降解是否优于阻断”"
    "四个工程节点，构成评价新LYTAC研究价值的基本坐标系。",
    "在活性评价层面，"
    "Banik等建立的溶酶体抑制剂确证范式，"
    "已成为化学偶联与基因编码平台的共同参照"
    + cite(1)
    + "。",
    "需要强调的是，"
    "受体-模块匹配并非一次性静态决策："
    "同一靶点在不同疾病分期、"
    "不同给药途径或联合用药背景下，"
    "最优受体臂与制备路线可能发生变化。"
    "例如，肿瘤免疫治疗早期可能优先考虑广谱降解以快速清除免疫检查点，"
    "而维持治疗阶段则可能更需要组织定向策略以降低脱靶免疫毒性。",
    "在实验设计层面，"
    "比较不同工程策略时应尽量固定靶向臂和连接子条件，"
    "仅改变受体臂类型，"
    "以避免将模块差异误判为受体生物学差异；"
    "评估内源性配体竞争时，"
    "应在多种细胞背景下重复验证，"
    "因为受体表达水平和内源性配体丰度均存在细胞系特异性。",
    "抗体靶向臂与小分子、适体靶向臂的选择，"
    "本质上是在亲和力、分子量、可修饰性和制备成本之间权衡："
    "较小分子量的靶向臂有利于组织渗透，"
    "高亲和力抗体或共价适体则更有利于维持高更新速率靶点的降解通量"
    + cite(4, 9)
    + "。",
    "在质量属性层面，"
    "化学偶联平台需关注偶联位点、偶联数和聚集体控制；"
    "基因编码平台需关注融合蛋白表达量、折叠正确性和宿主细胞糖基化差异。"
    "Kim等的酵母聚糖路线在降低受体臂合成难度的同时，"
    "仍需解决聚糖结构表征和批次间一致性评价问题"
    + cite(14)
    + "。"
    "无论采用哪条制备主线，"
    "均需在临床前阶段明确关键质量属性与降解活性的关联，"
    "否则难以支撑工艺放大和监管申报。",
    "【本节小结】LYTAC并非单一分子类型的线性迭代，"
    "而是以靶向臂、受体臂和连接子三模块为基础的可组合平台。"
    "药物设计应先在第一节受体决策框架内确定受体臂类型，"
    "再匹配适宜的靶向臂和连接子策略；"
    "平台化的前提是建立针对特定“靶点-组织-制备”组合的标准优化流程。",
]

BODY_SECTION3 = [
    "本节从疾病场景出发，"
    "系统评估降解策略相对单纯阻断的增益证据及其局限。"
    "评价的核心指标包括："
    "靶蛋白清除幅度、下游信号抑制持续时间、"
    "免疫功能激活程度及全身脱靶风险。",
    "在肿瘤免疫方向，PD-L1是最具代表性的比较靶点。"
    "Banik等首次证明PD-L1可被LYTAC降解"
    + cite(1)
    + "；"
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
    "在肿瘤膜蛋白靶点拓展方面，"
    "Pance等的KineTAC在HER2、TROP2等靶点上实现了51%—93%的降解效率"
    + cite(3)
    + "；"
    "Zhang等的iLYTAC在异种移植瘤模型中验证了抗肿瘤活性"
    + cite(13)
    + "。"
    "对于高更新速率的膜蛋白靶点，"
    "仍须评估靶蛋白重新合成是否会抵消降解效应，"
    "需在长期给药实验中加以考察"
    + cite(7)
    + "。",
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
    "在药效学评价方法上，"
    "建议将“蛋白清除动力学”与“功能终点”联合考察："
    "前者反映降解机制是否有效，"
    "后者反映降解是否产生预期的生物学后果。"
    "仅报告某一时间点的蛋白残留百分比，"
    "不足以支持降解优于阻断的判断。"
    "Li等和Fang等的研究分别从免疫检查点和膜蛋白靶点提供了"
    "“功能终点优于阻断对照”的体内证据"
    + cite(4, 15)
    + "，"
    "但其外推范围仍受动物模型和给药方案所限。",
    "在神经退行性疾病方向，Liu等设计了病灶激活LYTAC前体："
    "利用Aβ沉积区铜离子催化点击反应原位生成活性分子，"
    "经CD206受体介导Aβ纤维溶酶体降解"
    + cite(5)
    + "。"
    "该领域尚缺乏降解与阻断的直接对照实验，"
    "但病灶激活策略旨在解决全身给药脱靶这一阻断策略难以回避的问题。",
    "在自身免疫与炎症方向，"
    "胞外细胞因子和免疫球蛋白等靶点理论上亦适合LYTAC降解，"
    "但公开文献中系统的体内验证仍较少。",
    "对于分泌蛋白和聚集体靶点，"
    "评价逻辑需从“膜蛋白清除”转向“胞外池耗竭”，"
    "病理沉积物的清除效率是核心终点而非与传统抑制剂的活性比较。",
    "Fang等的SignalTAC在膜蛋白降解中显示出优于亲本抗体的抗肿瘤活性，"
    "与Li等的PD-L1研究形成互补，"
    "表明不同工程路径均可达成彻底清除靶蛋白的药理学目标"
    + cite(15)
    + "。",
    "跨疾病比较可见，降解优于阻断的证据在PD-L1肿瘤免疫中最为充分；"
    "肝靶向场景的优势主要体现在组织选择性而非药效幅度；"
    "神经退行性疾病则尚处于概念验证阶段。",
    "LYTAC的价值主张应与其所针对的靶点更新动力学相匹配："
    "对于更新快、持续补充的免疫检查点蛋白，"
    "彻底清除可能比单纯阻断更具优势；"
    "对于需要局部干预的肝代谢靶点，"
    "组织定向降解的价值主要体现在安全性；"
    "对于聚集体和纤维状靶点，"
    "核心指标是病理沉积物的清除效率。",
    "【本节小结】现有动物实验支持PD-L1降解在免疫激活上优于单纯阻断；"
    "肝靶向和神经退行性疾病方向的证据侧重于可行性与局部选择性，"
    "尚不足以普遍外推“降解优于阻断”的结论。"
    "后续研究需在各疾病场景中建立降解剂与阻断剂的标准化头对头比较。",
    "建立疾病场景特异性的评价标准，"
    "有助于避免将肿瘤免疫领域的结论不加区分地外推至其他适应症领域，"
    "从而提高临床前研究的判断准确性。",
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
    + cite(2, 12)
    + "。",
    "药代动力学方面，LYTAC与靶蛋白及受体形成三元复合物后常被一同送入溶酶体降解，"
    "体内有效浓度维持时间较短；"
    "高更新速率或高表达靶点可能需要更频繁给药"
    + cite(7)
    + "。",
    "监管与临床证据方面，"
    "LYTAC及相关eTPD技术仍主要处于临床前阶段，"
    "尚无正式发表的人体临床试验结果"
    + cite(17)
    + "。"
    "Zhao等强调，新兴降解平台需完成充分的机制验证、毒理评估和制剂研究方可进入临床"
    + cite(18)
    + "。"
    "与PROTAC相比，LYTAC的临床转化路径更缺少可参照的审评先例"
    + cite(10, 17)
    + "。",
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
    "其中，头对头比较实验应纳入功能终点而不仅是蛋白清除百分比，"
    "以明确降解策略相对阻断的真实获益边界。",
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
    "适体靶向臂虽具合成简便优势，"
    "但体内可能面临核酸酶降解和肾脏快速清除，"
    "需要通过化学修饰或递送系统改善药代性质；"
    "含聚糖或外源蛋白的偶联物还存在潜在免疫原性，"
    "需在反复给药动物实验中系统评估"
    + cite(4, 6)
    + "。",
    "从监管科学角度，"
    "LYTAC作为新型大分子偶联物或基因编码蛋白，"
    "其质量属性定义、杂质控制策略和免疫原性评价方法"
    "均尚未形成行业共识"
    + cite(17)
    + "。"
    "在缺乏临床先例的情况下，"
    "临床前研究的完整性和可重复性将成为监管审评的关键依据。",
    "此外，"
    "Huang等揭示的内源性配体竞争是CI-M6PR路线的生物学约束之一"
    + cite(12)
    + "；"
    "Fang等开发的内吞信号肽策略以及Liu等的病灶激活前体，"
    "则为克服组织选择性不足提供了工程应对思路"
    + cite(5, 15)
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
    "胞外蛋白降解技术的成药开发本质上需要药物化学、"
    "细胞生物学和临床药理学的交叉协作："
    "前者解决模块偶联与均一性，"
    "后者解释受体分选与组织定向，"
    "临床药理学则定义降解相对阻断的获益边界。"
    "单一学科视角的优化难以覆盖从分子设计到体内验证的全链条问题。",
    "在工艺放大层面，"
    "化学偶联路线需解决多步合成、纯化收率和结构确证问题；"
    "基因编码路线需解决宿主细胞株稳定性、"
    "下游纯化工艺和长期储存条件下的蛋白活性保持。"
    "Pance等的KineTAC和Zhang等的iLYTAC虽在实验室规模展现出良好的模块切换能力，"
    "但其向GMP生产体系迁移时，"
    "仍需重新定义关键工艺参数和放行标准"
    + cite(3, 13)
    + "。",
    "在安全性评价层面，"
    "脱靶降解正常组织蛋白是LYTAC特有的风险类型，"
    "不同于传统抗体阻断的“占用受体”风险。"
    "广谱CI-M6PR路线可能放大脱靶效应的不可预测性；"
    "肝定向ASGPR路线虽可缩小降解地理范围，"
    "但需评估肝细胞非靶蛋白的意外清除"
    + cite(2)
    + "。"
    "系统的组织病理学和蛋白组学监测，"
    "应成为LYTAC临床前安全性评价的常规组成部分。",
    "【本节小结】LYTAC的核心科学问题已基本明确，"
    "成药瓶颈主要集中在工艺可放大性、组织选择性和临床前证据体系三方面。"
    "在上述瓶颈未解决前，不宜将LYTAC视为接近临床的成熟技术平台。",
]

CONCLUSION = [
    "LYTAC技术将靶向蛋白降解从胞内蛋白酶体途径拓展至胞外及膜蛋白，"
    "其重要意义在于为约占蛋白质组相当比例的胞外靶点提供了可操作的干预策略，"
    "有望突破传统抗体和小分子抑制剂在“难以彻底清除致病蛋白”方面的局限"
    + cite(1, 18)
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
    + cite(17)
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
    "本综述的结构安排亦体现了上述逻辑："
    "第一节聚焦“选哪个受体、为什么”，"
    "第二节聚焦“如何装配分子模块”，"
    "第三、四节分别评估疾病证据与转化障碍。"
    "这一分层有助于避免受体生物学与分子工程论述的交叉重复，"
    "使读者在每一节获得清晰、不重叠的信息增量。",
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
        "说明降解效率取决于分选-回收通量竞争而非单纯的受体结合亲和力。"
        "（引自：Banik S M, et al. Nature, 2020, 584(7820): 291-297 [1]）",
    )

    add_heading_text(doc, "二、LYTAC的模块化分子工程设计", level=3)
    for para in BODY_SECTION2:
        add_paragraph(doc, para)
    add_figure_placeholder(
        doc,
        "图2  Banik等设计的LYTAC分子结构（左）及其对EGFR和PD-L1的降解效果（右），"
        "表明三模块架构可适用于多种膜蛋白靶点。"
        "（引自：Banik S M, et al. Nature, 2020, 584(7820): 291-297 [1]）"
        "图3  Ahn等GalNAc-LYTAC的肝定向降解示意及EGFR下游信号抑制时程，"
        "表明受体臂更换可将降解限制在肝脏并延长药效持续时间。"
        "（引自：Ahn G, et al. Nat Chem Biol, 2021, 17(9): 937-946 [2]）"
        "图4  Pance等KineTAC的模块化架构及多靶点降解效率，"
        "表明更换靶向臂即可在保留受体臂的前提下切换靶蛋白。"
        "（引自：Pance K, et al. Nat Biotechnol, 2023, 41(2): 273-281 [3]）"
        "图5  Li等共价适体LYTAC在肿瘤模型中的PD-L1清除及T细胞浸润，"
        "表明降解策略可在体内产生优于抗体阻断的免疫激活效应。"
        "（引自：Li Y, et al. J Am Chem Soc, 2023, 145(45): 24506-24521 [4]）",
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
