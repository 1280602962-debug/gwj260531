"""Editable publication tables and a paginated PDF proof from frozen tables.

Only display transformations are performed. All original table cells are kept
in Tables_Cell_Provenance.json, and source CSVs are copied without alteration.
"""
from __future__ import annotations
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import re
import textwrap
import zipfile
import sys
from decimal import Decimal,ROUND_HALF_UP

import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
from openpyxl import Workbook,load_workbook
from openpyxl.styles import Alignment,Font,Border,Side
from openpyxl.worksheet.pagebreak import Break

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'submission'
OUT=SOURCE/'Molecular_Diversity'
TABLES=OUT/'Tables'
CHECKS=OUT/'Checks'
ORDER=['EGFR/HER2','JAK1/JAK2','JAK1/TYK2','PIK3CA/mTOR','AChE/BChE','F2/F10','PPARG/PPARA','PPARA/PPARD']
MEASURE=fitz.Font(fontfile=font_manager.findfont('Arial'))

def wrap_measured(text,width):
    if '\n' in text:
        return '\n'.join(wrap_measured(part,width) for part in text.split('\n'))
    lines=[];line=''
    for word in text.replace('\n',' ').split():
        candidate=(line+' '+word).strip()
        if MEASURE.text_length(candidate,fontsize=8)<=width:
            line=candidate;continue
        if line:lines.append(line);line=''
        # Preserve normal words; wrap unusually long stored identifiers.
        while MEASURE.text_length(word,fontsize=8)>width:
            n=1
            while n<len(word) and MEASURE.text_length(word[:n+1],fontsize=8)<=width:n+=1
            lines.append(word[:n]);word=word[n:]
        line=word
    if line:lines.append(line)
    return '\n'.join(lines)
captions={}
for block in (SOURCE/'TABLE_CAPTIONS.md').read_text(encoding='utf8').split('## ')[1:]:
    lines=block.strip().splitlines()
    head=lines[0]
    en=[x for x in lines[1:] if x.strip() and x.isascii()]
    captions[head.split('.')[0]]=(head.split('. ',1)[1], ' '.join(en))

def parse(path):
    panels=[];heading='';rows=[]
    for line in path.read_text(encoding='utf8').splitlines():
        if line.startswith('## '):
            if rows:panels.append((heading,rows));rows=[]
            heading=line[3:]
        if line.startswith('|'):
            cells=[c.strip() for c in line.strip().strip('|').split('|')]
            if all(re.fullmatch(r':?-+:?',c) for c in cells):continue
            rows.append(cells)
        elif rows:
            panels.append((heading,rows));rows=[]
    if rows:panels.append((heading,rows))
    return panels

def clean(s):
    s=s.replace('Summary minuend','Comparison minimum AUROC').replace('Summary subtrahend','Reference minimum AUROC').replace('Summary difference','Minimum AUROC difference').replace('summary_min','minimum AUROC').replace('auc_B','AUROC B').replace('auc_A','AUROC A').replace('Physicochemical AUROC','Descriptor AUROC').replace('Alternative PDB','Alt. PDB').replace('Replaced side','Side')
    return s.replace('Summary min','Minimum AUROC').replace('summary min','minimum AUROC').replace('D_vs_A','Dual vs A-only').replace('D_vs_B','Dual vs B-only').replace('AUROC B','AUROC, target B').replace('AUROC A','AUROC, target A').replace(' (A)',' (Å)').replace('<2 A','<2 Å').replace('CONCORDANT_IMPROVE','Concordant improvement').replace('CONCORDANT_WORSEN','Concordant worsening').replace('NO_CHANGE','No change').replace('DISCORDANT','Discordant').replace('DOCUMENTED_NOT_REPLAYED','Documented; not replayed').replace('+A_only:', ' / A-only: ').replace('+B_only:', ' / B-only: ').replace('+neither:', ' / neither: ').replace('dual/A_only/B_only/neither','dual / A-only / B-only / neither').replace('Threshold-discordant repeats','Threshold crossing (n)').replace('Observations','n').replace('5 A ','5 Å ').replace('20 A.','20 Å.')

def typed(s):
    if re.fullmatch(r'-?\d+',s):return int(s)
    if re.fullmatch(r'-?\d+\.\d+(e[-+]?\d+)?',s,re.I):return float(s)
    return s

def combine(headers,rows,new,cols,fmt):
    """Merge display columns without changing numeric source values."""
    if not all(c in headers for c in cols):return headers,rows
    pos=min(headers.index(c) for c in cols)
    idx=[headers.index(c) for c in cols]
    keep=[i for i in range(len(headers)) if i not in idx]
    nh=[headers[i] for i in keep];nh.insert(pos,new)
    nr=[]
    for r in rows:
        rr=[r[i] for i in keep];rr.insert(pos,fmt([r[i] for i in idx]));nr.append(rr)
    return nh,nr

def transform(key,panel,raw):
    headers=raw[0].copy();rows=[r.copy() for r in raw[1:]]
    # Place confidence limits next to their own estimate, rather than scattering
    # lower/upper bounds across separate columns.
    for estimate,lo,hi in [('AUROC B','B CI lower','B CI upper'),('AUROC A','A CI lower','A CI upper'),('Summary min','Summary CI lower','Summary CI upper'),('Difference','Difference CI lower','Difference CI upper'),('Summary difference','Difference CI lower','Difference CI upper'),('Primary summary min','Primary CI lower','Primary CI upper'),('CAP5 summary min','CAP5 CI lower','CAP5 CI upper')]:
        headers,rows=combine(headers,rows,estimate+' (95% CI)',[estimate,lo,hi],lambda v:f'{v[0]} [{v[1]}, {v[2]}]')
    for cols in [['Dual n','A-only n','B-only n'],['Common dual n','Common A-only n','Common B-only n']]:
        headers,rows=combine(headers,rows,'n (dual / A-only / B-only)',cols,lambda v:' / '.join(v))
    if key=='S6':
        headers,rows=combine(headers,rows,'n (dual / single / neither)',['Dual n','Single-target n','Neither n'],lambda v:' / '.join(v))
        headers,rows=combine(headers,rows,'n (dual / negative)',['Dual n','Common negative n'],lambda v:' / '.join(v))
        headers=[h.replace('AUROC corresponding','Corresp. AUROC').replace('AUROC noncorresponding','Noncorresp. AUROC') for h in headers]
    if key=='S3':
        for prefix,label in [('Top-1','Top-1 RMSD median [min, max] (Å)'),('Best-saved','Best-saved RMSD median [min, max] (Å)')]:
            headers,rows=combine(headers,rows,label,[prefix+' median RMSD (A)',prefix+' min RMSD (A)',prefix+' max RMSD (A)'],lambda v:f'{v[0]} [{v[1]}, {v[2]}]')
        short={'Cognate CCD':'CCD','Seeds tested':'Seeds','Successful jobs':'Jobs','Modes returned':'Modes','Top-1 RMSD <2 A (n)':'Top-1 <2 Å (n)','Best-saved RMSD <2 A (n)':'Best <2 Å (n)'}
        headers=[short.get(h,h) for h in headers]
    if key=='S10' and panel.startswith('C.'):
        headers=[h.replace('Concordant improvement','Concordant\nimprovement').replace('Concordant worsening','Concordant\nworsening') for h in headers]
    headers=[clean(h) for h in headers];rows=[[clean(c) for c in r] for r in rows]
    if key=='S5':
        headers=[h.replace('n (dual / A-only / B-only)','n (dual/A/B)').replace('AUROC, target B (95% CI)','AUROC B [95% CI]').replace('AUROC, target A (95% CI)','AUROC A [95% CI]').replace('Minimum AUROC (95% CI)','Min. AUROC [95% CI]') for h in headers]
    if key=='S8' and panel.startswith('A.'):
        j=headers.index('Metric')
        for r in rows:r[j]=r[j].replace('AUROC, target B','AUROC B').replace('AUROC, target A','AUROC A').replace('minimum AUROC','Minimum AUROC')
    if key=='S2':
        # Identification and box geometry serve different reading tasks.
        keep=[i for i,h in enumerate(headers) if not h.startswith(('Center ','Size '))]
        geo=[headers.index(x) for x in ['Role','Target','PDB']]
        center=[i for i,h in enumerate(headers) if h.startswith('Center ')]
        size=[i for i,h in enumerate(headers) if h.startswith('Size ')]
        return [('A. Receptor identification',[headers[i] for i in keep],[[r[i] for i in keep] for r in rows]),
                ('B. Docking-box coordinates',['Role','Target','PDB','Center (x, y, z), Å','Size (x, y, z), Å'],[[r[i] for i in geo]+[', '.join(r[i] for i in center),', '.join(r[i] for i in size)] for r in rows])]
    if key=='S10' and panel.startswith('A.'):
        return [('A1. Method-specific top-10% counts',headers[:8],[r[:8] for r in rows]),
                ('A2. Method-specific top-10% fractions',headers[:2]+headers[8:],[r[:2]+r[8:] for r in rows])]
    return [(clean(panel),headers,rows)]

def main_table2():
    with (SOURCE/'01_Main_Tables/Table_2.csv').open(encoding='utf8') as h:raw=list(csv.DictReader(h))
    h=['Target pair','Scheme','N','k','Dual selected\nS1 → scheme','Dual recall (%)\nS1 → scheme']
    rows=[]
    for r in raw:
        pct=lambda x:str((Decimal(x)*100).quantize(Decimal('0.01'),rounding=ROUND_HALF_UP))
        scheme=r['Scheme vs S1']+('ᵃ' if r['Scheme vs S1']=='S4' else '')
        rows.append([r['Target pair'],scheme,r['Common N'],r['k'],r['Dual selected S1']+' → '+r['Dual selected new'],pct(r['Dual recall S1'])+' → '+pct(r['Dual recall new'])])
    return [('',h,rows)],raw

def md_table(h,rows):
    return '\n'.join(['| '+' | '.join(x.replace('\n',' ') for x in h)+' |','| '+' | '.join(['---']*len(h))+' |']+['| '+' | '.join(x.replace('\n',' ') for x in r)+' |' for r in rows])

def col_weights(h,rows):
    weights=[]
    for i,header in enumerate(h):
        vals=[r[i] for r in rows]
        longest=max([len(header)]+[len(x) for x in vals])
        if i==0 and header in ['Target pair','Target']:w=1.8
        elif header in ['Scheme','PDB','Role','Replaced side','k','N']:w=.9
        elif any(len(x)>85 for x in vals):w=4.3
        elif '95% CI' in header:w=2.8
        elif '→' in header:w=2.4
        elif 'x, y, z' in header:w=3.3
        elif longest<8:w=.9
        else:w=min(2.3,max(1.3,longest/13))
        weights.append(w)
    total=sum(weights);ratios=[x/total for x in weights]
    physical_pt=174/25.4*72*.97
    minimum={}
    for i,header in enumerate(h):
        if header in ['Target pair','Target','PDB','Scheme','Role','Contrast','Metric','Seeds','Modes','Side','Jobs','n (dual/A/B)']:
            value=max([MEASURE.text_length(header,fontsize=8)]+[MEASURE.text_length(r[i],fontsize=8) for r in rows])+16
            minimum[i]=value/physical_pt
    fixed={}
    for _ in range(len(h)):
        free=[i for i in range(len(h)) if i not in fixed]
        rem=1-sum(fixed.values());sw=sum(weights[i] for i in free)
        ratios=[fixed[i] if i in fixed else rem*weights[i]/sw for i in range(len(h))]
        bad=[i for i in free if ratios[i]<minimum.get(i,0)]
        if not bad:break
        for i in bad:fixed[i]=minimum[i]
    return ratios

def render_table(pdf,label,title,panel,headers,rows,footnote,checks):
    width=174/25.4;left=.015;right=.985
    ratios=col_weights(headers,rows)
    # Fixed 8 pt at final size. Row height grows with wrapped content.
    physical_pt=174/25.4*72*.97
    wrapwidth=[w*physical_pt-9 for w in ratios]
    wrapped_head=[wrap_measured(h,n) for h,n in zip(headers,wrapwidth)]
    wrapped_rows=[[wrap_measured(x,n) for x,n in zip(r,wrapwidth)] for r in rows]
    headpt=max(x.count('\n')+1 for x in wrapped_head)*9+10
    rowpts=[max(c.count('\n')+1 for c in r)*9+5 for r in wrapped_rows]
    # Split only between records. Repeated header and continued caption on each page.
    footerlines=[line for paragraph in footnote.splitlines() for line in textwrap.wrap(' '.join(paragraph.split()),115,break_long_words=False)]
    footpt=9*len(footerlines)+9
    bodylimit=min(450,663-headpt-footpt-82)
    groups=[]
    for i,r in enumerate(rows):
        if groups and rows[groups[-1][-1]][0]==r[0]:groups[-1].append(i)
        else:groups.append([i])
    pages=[];batch=[];pts=0
    for group in groups:
        rp=sum(rowpts[i] for i in group)
        if rp>bodylimit:raise ValueError(label+' contains an oversized row group')
        if batch and pts+rp>bodylimit:
            pages.append(batch);batch=[];pts=0
        batch.extend(group);pts+=rp
    if batch:pages.append(batch)
    for pi,idxs in enumerate(pages):
        footerlines=[line for paragraph in footnote.splitlines() for line in textwrap.wrap(' '.join(paragraph.split()),115,break_long_words=False)]
        footpt=9*len(footerlines)+9
        tablept=headpt+sum(rowpts[i] for i in idxs)
        heightpt=max(220,tablept+footpt+82)
        if heightpt>663:raise ValueError(label+' table exceeds final page height')
        fig=plt.figure(figsize=(width,heightpt/72));ax=fig.add_axes([left,.01,right-left,.98]);ax.axis('off')
        caption=label+'  '+title+(' (continued)' if pi else '')
        caplines=textwrap.wrap(caption,100,break_long_words=False)
        top=heightpt-8
        for line in caplines:fig.text(.018,top/heightpt,line,fontsize=9,fontweight='bold',va='top');top-=12
        if panel:fig.text(.018,(top-3)/heightpt,panel,fontsize=8.5,fontweight='bold',va='top');top-=18
        top-=8
        celltext=[wrapped_rows[i] for i in idxs]
        table=ax.table(cellText=celltext,colLabels=wrapped_head,colWidths=ratios,
                       bbox=[0,(top-tablept)/(.98*heightpt)-.01,1,tablept/(.98*heightpt)],cellLoc='center')
        table.auto_set_font_size(False);table.set_fontsize(8)
        for (r,c),cell in table.get_celld().items():
            cell.set_edgecolor('#222222');cell.set_linewidth(.45);cell.visible_edges=''
            cell.set_height((headpt if r==0 else rowpts[idxs[r-1]])/tablept)
            if r==0:cell.visible_edges='TB';cell.set_text_props(weight='bold')
            elif r==len(celltext):cell.visible_edges='B'
            if r and c==0:cell.set_text_props(ha='left')
            if r>1 and rows[idxs[r-1]][0]!=rows[idxs[r-2]][0]:cell.visible_edges+='T'
        footer_y=(top-tablept-12)/heightpt
        fig.text(.018,footer_y,'\n'.join(footerlines),fontsize=8,va='top',linespacing=1.2)
        fig.canvas.draw()
        renderer=fig.canvas.get_renderer();fw,fh=fig.canvas.get_width_height()
        overflow=[]
        for (r,c),cell in table.get_celld().items():
            t=cell.get_text();bb=t.get_window_extent(renderer);cb=cell.get_window_extent(renderer)
            if bb.width>cb.width+1 or bb.height>cb.height+1:overflow.append([r,c,t.get_text()])
        if overflow:raise ValueError(label+' clipped cells: '+str(overflow)[:450])
        pdf.savefig(fig);checks.append({'table':label,'panel':panel,'page_part':pi+1,'rows':len(idxs),'width_mm':174,'height_mm':round(heightpt/72*25.4,2),'font_pt':8,'cell_overflow':overflow})
        plt.close(fig)

def add_sheet(wb,name,title,headers,rows,note):
    ws=wb.create_sheet(name[:31]);n=len(headers)
    ws.sheet_view.showGridLines=False
    ws.merge_cells(start_row=1,start_column=1,end_row=1,end_column=n)
    ws.cell(1,1,title).font=Font(name='Arial',size=11,bold=True,color='000000')
    ws.row_dimensions[1].height=32
    for j,h in enumerate(headers,1):ws.cell(3,j,h.replace('\n',' '))
    for row in rows:ws.append([typed(x) for x in row])
    thin=Side(style='thin',color='000000')
    for row in ws.iter_rows(min_row=3,max_row=3+len(rows),max_col=n):
        for c in row:
            c.font=Font(name='Arial',size=10,color='000000',bold=c.row==3)
            c.alignment=Alignment(horizontal='left' if c.column==1 else 'center',vertical='center',wrap_text=True)
            if c.row==3:c.border=Border(top=thin,bottom=thin)
            elif c.row==3+len(rows):c.border=Border(bottom=thin)
            if isinstance(c.value,float):c.number_format='0.0000'
    weights=col_weights(headers,rows)
    for j,w in enumerate(weights,1):ws.column_dimensions[ws.cell(3,j).column_letter].width=max(9,w*145)
    ws.row_dimensions[3].height=48
    for i,r in enumerate(rows,4):
        lines=max(math.ceil(len(str(v))/max(8,ws.column_dimensions[ws.cell(i,j).column_letter].width-2)) for j,v in enumerate(r,1))
        ws.row_dimensions[i].height=max(24,lines*14+7)
    end=5+len(rows);ws.merge_cells(start_row=end,start_column=1,end_row=end+3,end_column=n)
    c=ws.cell(end,1,note);c.font=Font(name='Arial',size=9,color='000000');c.alignment=Alignment(wrap_text=True,vertical='top')
    ws.row_dimensions[end].height=max(45,math.ceil(len(note)/(sum(ws.column_dimensions[ws.cell(3,j).column_letter].width for j in range(1,n+1))*.85))*12+12)
    ws.freeze_panes='B4';ws.print_title_rows='1:3';ws.print_options.horizontalCentered=True
    ws.sheet_properties.pageSetUpPr.fitToPage=True
    ws.page_setup.orientation='landscape' if n>7 else 'portrait'
    ws.page_setup.paperSize=ws.PAPERSIZE_A4;ws.page_setup.fitToWidth=1;ws.page_setup.fitToHeight=0
    ws.print_area=f'A1:{ws.cell(end+3,n).coordinate}'
    return ws

def main():
    update_workbooks='--no-workbooks' not in sys.argv
    for p in [TABLES,CHECKS]:p.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42})
    mainwb=Workbook();mainwb.remove(mainwb.active)
    siwb=Workbook();siwb.remove(siwb.active)
    allpanels=[];provenance=[];checks=[];source_checks=[]
    archive=zipfile.ZipFile(SOURCE/'03_Supplementary_Data/Supplementary_Data.zip')
    rawdest=OUT/'Source_Data/Tables';rawdest.mkdir(parents=True,exist_ok=True)
    for name in archive.namelist():
        if name.startswith('Numbered_tables/') and name.endswith('.csv'):
            (rawdest/Path(name).name).write_bytes(archive.read(name))
    paths=[('1',SOURCE/'01_Main_Tables/Table_1.md'),('2',SOURCE/'01_Main_Tables/Table_2.md')]+[(f'S{i}',SOURCE/f'02_SI_Tables/Table_S{i}.md') for i in range(1,11)]
    for key,path in paths:
        label='Table '+key
        title,note=captions[label]
        if key=='2':
            panels,raw=main_table2();note='Each row uses its own common two-target complete-case population (N); k = ceil(0.10 N). Arrows show S1 → comparison scheme. Dual recall (%) is the fraction of all dual ligands recovered. Counts and recall retain comparison-specific S1 baselines. Full class composition: Table S10.\nᵃ Supporting comparison against S1.'
            provenance.append({'table':key,'source':str(path.relative_to(ROOT)),'source_csv':raw})
        else:
            parsed=parse(path);panels=[]
            for panel,raw in parsed:
                source_name='Numbered_tables/Table_'+key+('_'+panel.split('.')[0] if panel else '')+'.csv'
                original=list(csv.reader(archive.read(source_name).decode('utf-8-sig').splitlines()))
                if len(original)!=len(raw) or original[0]!=raw[0]:raise ValueError('Source table layout mismatch '+source_name)
                for ri,(a,z) in enumerate(zip(raw[1:],original[1:]),1):
                    for ci,(actual,full) in enumerate(zip(a,z)):
                        try:
                            af=float(actual);zf=float(full)
                            good=abs(af-zf)<=.0000500001
                        except ValueError:good=actual==full
                        if not good:raise ValueError(f'{source_name} cell {ri},{ci}: {actual} != {full}')
                        source_checks.append({'table':key,'panel':panel,'row':ri,'column':ci,'pass':True})
                provenance.append({'table':key,'panel':panel,'source':str(path.relative_to(ROOT)),'source_cells':raw})
                panels.extend(transform(key,panel,raw))
        if key=='1':note='Counts are target-pair–ligand observations; the 867 observations represent 845 distinct ligand entities. Target A and target B follow the order in each pair name. Activity is defined by pChEMBL ≥ 6.0. PDB identifiers refer to the receptors used in the primary analysis'
        if key=='S6':note+=' Corresp., corresponding pocket; noncorresp., noncorresponding pocket. Grouped counts follow the class order in the column heading'
        if key=='S7':note+=' OOF refers to the learned chemistry models; S1 is a direct docking score'
        if key=='S8':note+=' AUROC B compares dual with A-only; AUROC A compares dual with B-only. Minimum AUROC is the smaller directional AUROC'
        if key=='S3':note+=' CCD denotes the cognate ligand Chemical Component Dictionary identifier; jobs denotes successful docking jobs; modes denotes saved poses per job'
        if key=='S5':note='AUROC B: dual versus A-only; AUROC A: dual versus B-only. Counts: dual/A-only/B-only. S1, Vina; S2, fixed-pose GNINA; S3, fixed-pose RTMScore; S4, maximum saved-pose CNNscore; S5, independent GNINA. Higher scores are favorable. Paired contrasts use their own shared members and need not equal differences between estimates in panel A. CI: 95% class-stratified percentile bootstrap (10,000 resamples; seed 271828). S4-S2: saved-pose reranking. Supporting S4-S1 in Table S10 has no added CI.'
        note=clean(note).rstrip().rstrip('.')+'.'
        md=['# '+label+' '+title,'']
        for j,(panel,h,rows) in enumerate(panels):
            if panel:md+=['## '+panel,'']
            md += [md_table(h,rows),'']
            allpanels.append((label,title,panel,h,rows,note))
            wb=mainwb if key in ['1','2'] else siwb
            if update_workbooks:add_sheet(wb,key+('_'+str(j+1) if len(panels)>1 else ''),label+' '+title+((' — '+panel) if panel else ''),h,rows,note)
        md += [note,'']
        (TABLES/('Table'+key+'.md')).write_text('\n'.join(md),encoding='utf8')
    # Preserve machine-readable main data exactly, separately from formatting.
    for f in (SOURCE/'01_Main_Tables').glob('*.csv'):(TABLES/f.name).write_bytes(f.read_bytes())
    if update_workbooks:
        mainwb.save(TABLES/'Main_Tables.xlsx');siwb.save(TABLES/'Supplementary_Tables.xlsx')
    for name,original in ([('Main_Tables.xlsx',mainwb),('Supplementary_Tables.xlsx',siwb)] if update_workbooks else []):
        reloaded=load_workbook(TABLES/name,data_only=False)
        assert reloaded.sheetnames==original.sheetnames
        for a,bb in zip(original,reloaded):
            assert [[c.value for c in r] for r in a]==[[c.value for c in r] for r in bb]
    with PdfPages(OUT/'Tables_Review.pdf') as pdf:
        for args in allpanels:render_table(pdf,*args,checks)
    proof=fitz.open(OUT/'Tables_Review.pdf')
    for i,p in enumerate(proof):p.get_pixmap(matrix=fitz.Matrix(1.7,1.7)).save(CHECKS/f'TablePage{i+1:02d}.png')
    proof.close()
    (CHECKS/'Tables_Cell_Provenance.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    (CHECKS/'Tables_QA.json').write_text(json.dumps({'scientific_recomputation':False,'source_tables':12,'display_panels':len(allpanels),'source_cell_checks_passed':len(source_checks),'source_cell_checks_failed':0,'review_pages':len(checks),'workbook_roundtrip':'all values identical','full_precision_source_csvs':'Source_Data/Tables','pages':checks},indent=2)+'\n')
    # Native HTML tables can be pasted into a manuscript editor as editable cells.
    body=[]
    for label,title,panel,h,rows,note in allpanels:
        body.append('<section><h2>'+html.escape(label+' '+title)+'</h2>')
        if panel:body.append('<h3>'+html.escape(panel)+'</h3>')
        body.append('<table><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in h)+'</tr></thead><tbody>')
        for row in rows:body.append('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in row)+'</tr>')
        body.append('</tbody></table><p class="note">'+html.escape(note).replace('\n','<br>')+'</p></section>')
    page='<!doctype html><html lang="en"><meta charset="utf-8"><title>Molecular Diversity tables</title><style>body{font:10pt Arial;color:#000;max-width:1100px;margin:30px auto}h2{font-size:12pt}h3{font-size:10pt}table{border-collapse:collapse;width:100%;font-size:9pt}th{border-top:1.5px solid #000;border-bottom:1px solid #000;padding:7px}td{padding:6px;text-align:center}td:first-child{text-align:left}tbody tr:last-child td{border-bottom:1.5px solid #000}.note{font-size:9pt;line-height:1.4}section{margin-bottom:45px}thead{display:table-header-group}@media print{section{break-before:page}}</style>'+''.join(body)+'</html>'
    (TABLES/'Editable_Tables.html').write_text(page,encoding='utf8')
    print(json.dumps({'tables':12,'panels':len(allpanels),'pages':len(checks),'workbooks':2}))

if __name__=='__main__':main()
