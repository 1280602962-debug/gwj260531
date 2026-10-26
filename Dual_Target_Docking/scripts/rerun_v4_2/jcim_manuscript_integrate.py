"""Build current manuscript tables, review PDFs and a source-data package.

Presentation only: no AUROC, bootstrap, model fitting or ranking is recomputed.
Workbook authoring uses @oai/artifact-tool through jcim_manuscript_tables.mjs.
"""
from __future__ import annotations
import csv
import hashlib
import json
import os
import sys
import subprocess
import tempfile
import zipfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parent
MAN = ROOT / 'manuscript'
SCRIPTS = Path(__file__).parent
ORDER = ['EGFR/HER2','JAK1/JAK2','JAK1/TYK2','PIK3CA/mTOR','AChE/BChE','F2/F10','PPARG/PPARA','PPARA/PPARD']
Q = ROOT / 'extension_jcim_validation/pik3ca_mtor_q110_v42'
PROVENANCE = ROOT / 'results/current/provenance_audit_20261023'
CAP5 = ROOT / 'results/current/chemical_diversity_cap5_20261023'
EXTRA_SHEETS = [
    ('S19_Source_Coverage', PROVENANCE / 'CURRENT_867_SOURCE_COVERAGE_SUMMARY.csv'),
    ('S20_Source_Crosswalk', PROVENANCE / 'CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv'),
    ('S21_Q110_Activity_Rows', PROVENANCE / 'Q110_CHEMBL37_ACTIVITY_ROWS_LONG.csv'),
    ('S22_Selection_Provenance', PROVENANCE / 'SEVEN_PAIR_SELECTION_PROVENANCE.csv'),
    ('S23_CAP5_Class_Counts', CAP5 / 'CAP5_CLASS_COUNTS.csv'),
    ('S24_CAP5_Excluded', CAP5 / 'CAP5_EXCLUDED_LIGANDS.csv'),
    ('S25_CAP5_PRIMARY', CAP5 / 'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv'),
    ('S26_CAP5_Method_Deltas', CAP5 / 'CAP5_PRIMARY_METHOD_DELTAS.csv'),
    ('S27_CAP5_Negative_Class', CAP5 / 'CAP5_MODULE_B_NEGATIVE_CLASS.csv'),
    ('S28_CAP5_Pocket', CAP5 / 'CAP5_MODULE_C_POCKET_CORRESPONDENCE.csv'),
    ('S29_CAP5_Chemistry', CAP5 / 'CAP5_MODULE_D_CHEMISTRY_UNITS.csv'),
    ('S30_CAP5_Topk', CAP5 / 'CAP5_MODULE_I_TOPK.csv'),
    ('S31_CAP5_Vs_Primary', CAP5 / 'CAP5_PRIMARY_VS_ORIGINAL.csv'),
]



def read_csv(path):
    with path.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames),list(reader)


def table(name,path):
    fields,rows=read_csv(path)
    return {'name':name,'matrix':[fields]+[[r.get(k,'') for k in fields] for r in rows]}


def build_reviews():
    import fitz
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,KeepTogether
    from xml.sax.saxutils import escape
    with tempfile.TemporaryDirectory() as tmp:
        font=Path(tmp)/'cjk.ttf'
        font.write_bytes(fitz.Font('china-s').buffer)
        pdfmetrics.registerFont(TTFont('ReviewCJK',str(font)))
        styles={
            'body':ParagraphStyle('body',fontName='ReviewCJK',fontSize=9.5,leading=15,spaceAfter=7,wordWrap='CJK'),
            'h1':ParagraphStyle('h1',fontName='ReviewCJK',fontSize=15,leading=21,spaceAfter=13),
            'h2':ParagraphStyle('h2',fontName='ReviewCJK',fontSize=11,leading=17,spaceBefore=10,spaceAfter=7),
        }
        def paragraphs(text):
            result=[]
            for block in text.split('\n\n'):
                block=block.strip()
                if not block:continue
                if block.startswith('## R8.'):
                    result.append(PageBreak())
                level='h1' if block.startswith('# ') else 'h2' if block.startswith('## ') else 'body'
                block=block.lstrip('# ')
                result.append(Paragraph(escape(block).replace('\n','<br/>'),styles[level]))
            return result
        main_figs=['Scheme01_Study_Design','Fig01_Primary_Directional_and_Deltas','Fig02_Negative_Class_Contrast','Fig03_Chemistry_and_Docking_Increment','Fig04_Pocket_Contrast','Fig05_Seed_and_Receptor_Sensitivity','Fig06_Directional_and_Topk_Changes']
        si_figs=['FigS01_Descriptors_D_vs_A','FigS02_Descriptors_D_vs_B','FigS03_Cross_Class_Similarity','FigS04_All_Method_AUROC','FigS05_Five_Seeds','FigS06_Alternative_Receptors','FigS07_Missingness','FigS08_Topk_and_Concordance']
        for output,source,figures in [('Results_Review.pdf','Results_CN.md',main_figs),('Supporting_Information_Review.pdf','Supporting_Information_CN.md',si_figs)]:
            story=paragraphs((MAN/source).read_text())
            story += [Paragraph('完整表格见对应 Excel 文件，原始精度见 CSV 数据包。',styles['body'])]
            for stem in figures:
                story.append(PageBreak())
                story.append(Paragraph(stem.replace('_',' '),styles['h2']))
                image=Image(str(ROOT/'results/jcim_stage2_figures'/f'{stem}.png'))
                scale=min(500/image.imageWidth,690/image.imageHeight)
                image.drawWidth=image.imageWidth*scale
                image.drawHeight=image.imageHeight*scale
                story.append(image)
            story.append(PageBreak())
            story.extend(paragraphs((MAN/'Figure_Captions_CN.md').read_text()))
            def footer(canvas,doc):
                canvas.setFont('ReviewCJK',8)
                canvas.drawRightString(A4[0]-40,24,str(doc.page))
            SimpleDocTemplate(str(MAN/output),pagesize=A4,rightMargin=40,leftMargin=40,topMargin=35,bottomMargin=38).build(story,onFirstPage=footer,onLaterPages=footer)


def main():
    MAN.mkdir(exist_ok=True)
    fields,pop=read_csv(ROOT/'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    eligible=[r for r in pop if r['activity_eligible']=='1' and r['class']]
    t1=[['index','pair','target_A','pdb_A','target_B','pdb_B','dual','A_only','B_only','neither']]
    for i,pair in enumerate(ORDER,1):
        rows=[r for r in eligible if r['pair']==pair]
        pdbs={(r['pdb_A'],r['pdb_B']) for r in rows}
        assert len(pdbs)==1
        pa,pb=next(iter(pdbs));a,b=pair.split('/');c=Counter(r['class'] for r in rows)
        t1.append([i,pair,a,pa,b,pb,*[c[k] for k in ['dual','A_only','B_only','neither']]])
    desc=read_csv(ROOT/'results/jcim_stage2/module_I_descriptive.csv')[1]
    by={(r['pair'],r['method'],r['population_kind']):r for r in desc}
    header=['index','pair']
    for method in ['M1','M1b','M2','M3']:
        tag=method+'_supporting' if method=='M1b' else method
        header.extend([tag+'_N_common',tag+'_k',tag+'_dual_M0',tag+'_dual_new'])
    t2=[header]
    for i,pair in enumerate(ORDER,1):
        values=[i,pair]
        for method in ['M1','M1b','M2','M3']:
            old=by[pair,'M0','pairwise_common_vs_'+method];new=by[pair,method,'pairwise_common_vs_M0']
            assert (old['n_total'],old['k'])==(new['n_total'],new['k'])
            values.extend(int(x) for x in [old['n_total'],old['k'],old['n_dual_top'],new['n_dual_top']])
        t2.append(values)
    sheets=[('S1_Mapping',ROOT/'results/current/pair_ligand_mapping.csv'),('S2_Population_Layers',ROOT/'analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv'),('S3_PRIMARY_Metrics',ROOT/'results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv'),('S4_PRIMARY_Deltas',ROOT/'results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv')]
    for name,file in [('S5_Negative_Class','module_B_negative_class.csv'),('S6_Wrong_Pocket','module_C_wrong_pocket.csv'),('S7_Chemistry_Units','module_D_chemistry_units.csv'),('S8_Composition','module_D1_composition.csv'),('S9_Chemistry_OOF','module_D_chemistry_oof_ligands.csv'),('S10_Five_Seeds','module_E_five_seed.csv'),('S11_Seed_Summaries','module_E_five_seed_summary.csv'),('S12_Alt_Receptors','module_F_alt_receptor.csv'),('S13_Missingness','module_H_missingness.csv'),('S14_Ranking_Descriptive','module_I_descriptive.csv'),('S15_Ranking_Comparison','module_I_method_compare.csv'),('S16_Ligand_Ranks','module_I_ligand_ranks.csv')]:
        sheets.append((name,ROOT/'results/jcim_stage2'/file))
    sheets.extend([('S17_Concordance_Counts',ROOT/'results/jcim_stage2_figures/module_I_concordance_counts.csv'),('S18_Q110_Strict_Subset',Q/'Q110_STRICT_SUBSET_METRICS.csv')])
    core_sheets = list(sheets)
    sheets.extend(EXTRA_SHEETS)
    supporting_only = '--supporting-update' in sys.argv
    workbooks = [] if supporting_only else [{'file':'Main_Tables.xlsx','sheets':[{'name':'Table1_Label_Counts','matrix':t1},{'name':'Table2_Common_Topk','matrix':t2}]}]
    selected = EXTRA_SHEETS if supporting_only else sheets
    workbooks.append({'file':'Supplementary_Tables.xlsx','sheets':[table(name,path) for name,path in selected]})
    payload={'root':str(MAN),'workbooks':workbooks}
    with tempfile.TemporaryDirectory() as tmp:
        source=Path(tmp)/'tables.json';source.write_text(json.dumps(payload))
        node=os.environ.get('CODEX_PRIMARY_RUNTIME_NODE','node')
        subprocess.run([node,str(SCRIPTS/'jcim_manuscript_tables.mjs'),str(source)],check=True)
    files=[p for _,p in sheets]
    files.extend([ROOT/'ANALYSIS_BASELINE_README.md',ROOT/'ANALYSIS_BASELINE_MANIFEST.csv',ROOT/'CURRENT_RESULT_AUTHORITIES.json',ROOT/'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv',ROOT/'results/current/RESULT_SYNC_REPORT.json',ROOT/'results/canonical/model_fold_assignments.csv',ROOT/'results/jcim_stage2/RUN_MANIFEST.json',ROOT/'reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/parent_alias_registry.csv',ROOT/'reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/FORMAL_AUTHORITY_PATHS.yaml',ROOT/'reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml',ROOT/'reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml',Q/'PIK3CA_MTOR_Q110_MEMBERSHIP.csv',Q/'Q110_POPULATION.csv',Q/'Q110_FOLD_ASSIGNMENT_RULE.json',Q/'Q110_CURRENT_RESULT_AUTHORITIES.json',MAN/'Source_Register.csv'])
    # Keep original provenance bytes; the current crosswalk defines which
    # original observations are used in the 867-observation main population.
    for directory in (PROVENANCE, CAP5, ROOT / 'activity_label_source_audit_20261008'):
        files.extend(sorted(p for p in directory.iterdir() if p.is_file() and p.suffix in {'.csv','.json','.md','.py'}))
    files = list(dict.fromkeys(files))
    manifest=MAN/'source_views/SHA256SUMS.csv';manifest.parent.mkdir(exist_ok=True)
    with manifest.open('w',newline='') as out:
        w=csv.writer(out);w.writerow(['relative_path','sha256','bytes'])
        for p in files:w.writerow([p.relative_to(REPO).as_posix(),hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_size])
    with zipfile.ZipFile(MAN/'Supplementary_Source_Data.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,p.relative_to(REPO).as_posix())
        z.write(manifest,'SHA256SUMS.csv')
    build_reviews()
    print('CURRENT_MANUSCRIPT_BUILT',len(pop),len(eligible),len(sheets))


if __name__=='__main__':main()
