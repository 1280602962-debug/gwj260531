"""Molecular Diversity artwork from frozen result tables; no scientific reruns.

Run from any directory. The legacy builder supplies figure-specific source joins.
This module changes presentation, audits the plotted values and writes journal
exports separately from the frozen source results.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import importlib.util

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgb
import numpy as np
import fitz
from PIL import Image

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('legacy', HERE/'build_submission_figures.py')
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)
ROOT = b.ROOT
OUT = ROOT/'submission/Molecular_Diversity'
MAIN = OUT/'Main_Figures'
SI = OUT/'Supplementary_Figures'
QA = OUT/'Checks'
DATA = OUT/'Source_Data'
BASE_COMMIT = 'ce81b011bab449c6cf5e2b9664626ed808be5857'
WIDTH = 174/25.4
audit = []
figures = []

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def style():
    family = b.style()
    plt.rcParams.update({'font.size':8.5,'axes.labelsize':8.5,
        'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,
        'legend.fontsize':8,'svg.fonttype':'path', 'hatch.linewidth':0.45,
        'pdf.fonttype':42, 'ps.fonttype':42})
    return family

def title(ax, letter, text):
    # Panel identifiers and short condition labels; full captions are external.
    short = {'Evaluation-set composition':'Experimental classes',
        'Activity-source coverage':'Activity records',
        'Directional Vina discrimination':'Directional AUROC',
        'Minimum directional AUROC':'Minimum directional AUROC',
        'Same-member AUROC':'AUROC', 'Paired AUROC difference':'AUROC difference',
        'Chemistry models and Vina':'Matched-population AUROC',
        'Addition of the Vina score':'ECFP4 + S1 minus ECFP4',
        'Five search seeds':'Search seeds',
        'One-target receptor substitutions':'Receptor substitutions',
        'CAP5: Vina summary_min':'CAP5', 'Strict-activity subset':'Strict subset'}
    ax.set_title(letter+'  '+short.get(text,text),loc='left',pad=10,
                 fontsize=9,fontweight='bold')

def contrast_text(color):
    rgb=np.asarray(to_rgb(color))
    linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
    luminance=linear@np.array([.2126,.7152,.0722])
    return 'black' if (luminance+.05)/.05>=4.5 else 'white'

def check_records(records, sources):
    """Every exported row must match its frozen source, by content, not position."""
    source_rows=[r for s in sources if s.endswith('.csv') for r in b.read(s)]
    matches=0
    # Some graphics derive counts and fractions from ligand-level records.
    if records and source_rows and not all('panel'in r and 'n'in r for r in records):
        for r in records:
            clean={k:str(v) for k,v in r.items() if k not in ('panel','population')}
            candidates=[s for s in source_rows if all(k in s and str(s[k])==v for k,v in clean.items())]
            if not candidates:
                raise ValueError('Figure row not found verbatim in authoritative inputs: '+str(clean)[:250])
            matches+=1
    return matches

def finish(fig,number,name,section,table,sources,records,caption,caption_cn,si=False):
    folder=SI if si else MAIN
    stem='SchemeS1' if number=='Scheme_S1' else 'Fig'+number
    # Keep journal width fixed rather than scaling exported artwork afterwards.
    fig.set_size_inches(WIDTH,fig.get_figheight(),forward=True)
    for a in fig.axes:
        a.tick_params(axis='both',labelsize=8)
        for which in ['x','y']:
            labels=(a.get_xticklabels() if which=='x' else a.get_yticklabels())
            if any('AUROC_' in t.get_text() for t in labels):
                vals=a.get_xticks() if which=='x' else a.get_yticks()
                clean=[t.get_text().replace('AUROC_B','AUROC B').replace('AUROC_A','AUROC A') for t in labels]
                (a.set_xticks if which=='x' else a.set_yticks)(vals,clean)
        for im in a.images:
            # Preserve native heat-map cells in vector exports; resampling at
            # 600 dpi would inflate EPS files into tens of megabytes.
            im.set_interpolation('none')
            im.set_resample(False)
            data=np.asarray(im.get_array())
            if data.ndim!=2:continue
            for t in a.texts:
                x,y=t.get_position();i,j=round(y),round(x)
                if abs(y-i)<1e-8 and abs(x-j)<1e-8 and 0<=i<data.shape[0] and 0<=j<data.shape[1]:
                    t.set_color(contrast_text(im.cmap(im.norm(data[i,j]))[:3]))
    if number=='S2':
        fig.subplots_adjust(bottom=.19)
        cb=fig.colorbar(fig.axes[0].images[0],cax=fig.add_axes([.31,.11,.50,.012]),orientation='horizontal')
        maximum=getattr(fig,'_missing_count_scale',None)
        cb.set_ticks(range(maximum+1) if maximum is not None else [0,.5,1]);cb.ax.tick_params(labelsize=8)
        cb.set_label('Unavailable scores (n)' if maximum is not None else 'Valid / expected score fraction',fontsize=8)
    for t in fig.findobj(Text):
        if t.get_text().strip():
            t.set_text(t.get_text().replace('summary_min','Minimum AUROC').replace('AUROC_B','AUROC B').replace('AUROC_A','AUROC A').replace('square angstrom','Å²').replace('(angstrom)','(Å)'))
            t.set_fontsize(max(8,t.get_fontsize()))
            if t.get_text()=='Scoring and pose-selection workflows':t.set_visible(False)
    # EPS does not support alpha. Composite point colors against white before
    # writing any format so all exports represent the same artwork.
    for c in fig.findobj(PathCollection):
        alpha=c.get_alpha()
        if alpha is not None and alpha<1:
            for getter,setter in [(c.get_facecolors,c.set_facecolors),(c.get_edgecolors,c.set_edgecolors)]:
                colors=getter().copy()
                if len(colors):
                    colors[:,:3]=colors[:,:3]*alpha+1-alpha
                    colors[:,3]=1
                    setter(colors)
            c.set_alpha(1)
    hook=globals().get('publication_layout_hook')
    if hook is not None:hook(fig,number)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    fw,fh=fig.canvas.get_width_height()
    outside=[]; sizes=[]
    off_axis_ticks=set()
    for a in fig.axes:
        xmin,xmax=sorted(a.get_xlim());ymin,ymax=sorted(a.get_ylim())
        off_axis_ticks.update(id(t) for v,t in zip(a.get_xticks(),a.get_xticklabels()) if not xmin-1e-10<=v<=xmax+1e-10)
        off_axis_ticks.update(id(t) for v,t in zip(a.get_yticks(),a.get_yticklabels()) if not ymin-1e-10<=v<=ymax+1e-10)
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip() or id(t) in off_axis_ticks:continue
        box=t.get_window_extent(renderer)
        if box.width<=0 or box.height<=0:continue
        # Ignore invisible tick labels that Matplotlib retains beyond limits.
        if t.axes is not None and t in t.axes.get_xticklabels():
            pass
        sizes.append(t.get_fontsize())
        if box.x0< -1 or box.y0< -1 or box.x1>fw+1 or box.y1>fh+1:
            outside.append(t.get_text())
    if outside:raise ValueError(stem+' has text outside canvas: '+str(outside))
    if fig.get_figheight()*25.4>234:raise ValueError(stem+' exceeds page height')
    for ext in ['pdf','eps','svg','png']:
        fig.savefig(folder/(stem+'.'+ext),format=ext,dpi=600)
    png_path=folder/(stem+'.png')
    with Image.open(png_path) as png:rgb=png.convert('RGB')
    rgb.save(png_path,dpi=(600,600))
    pd=fitz.open(folder/(stem+'.pdf'))
    pd[0].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(QA/(stem+'.png'))
    font_list=pd[0].get_fonts(full=True)
    if not font_list or any(not row[3].split('+')[0] for row in font_list):
        raise ValueError('Missing embedded font '+stem)
    font_details=[]
    for row in font_list:
        extracted=pd.extract_font(row[0])
        if not extracted[-1]:raise ValueError('Unembedded PDF font '+str(row))
        font_details.append({'font':row[3],'type':row[2],'embedded_bytes':len(extracted[-1])})
    pd.close()
    nmatch=check_records(records,sources)
    if records:
        keys=list(dict.fromkeys(k for row in records for k in row))
        with (DATA/(stem+'.csv')).open('w',newline='',encoding='utf8') as h:
            w=csv.DictWriter(h,keys);w.writeheader();w.writerows(records)
    # Clean publisher style and remove working-script wording from captions.
    caption=caption.replace('stored ','').replace('recorded ','')
    caption=caption.replace('summary_min','minimum AUROC').replace('paired directional minimum AUROC change','paired change in minimum AUROC').replace('AUROC_B','AUROC B').replace('AUROC_A','AUROC A')
    caption=caption.replace('No new confidence interval is added. ','')
    caption=caption.replace('no additional member is lost','no additional member is excluded')
    caption=caption.replace('(b) minimum AUROC','(b) Minimum AUROC').rstrip('.')
    figures.append({'number':number,'stem':stem,'group':'SI' if si else 'Main',
        'section':section,'table':table,'sources':sources,'caption':caption})
    audit.append({'figure':number,'width_mm':174,'height_mm':round(fig.get_figheight()*25.4,2),
        'min_font_pt':min(sizes),'outside_canvas':outside,'matched_source_rows':nmatch,
        'embedded_fonts':font_details,'scientific_recomputation':False})
    plt.close(fig)

def fig1():
    pop=[r for r in b.read('results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv') if r['activity_eligible']=='1' and r['class']]
    cross=b.read('results/current/provenance_audit_20261023/CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv')
    fig,ax=plt.subplots(1,2,figsize=(WIDTH,4.7),gridspec_kw={'width_ratios':[1.05,1]})
    fig.subplots_adjust(left=.19,right=.985,bottom=.27,top=.90,wspace=.22)
    y=np.arange(8)[::-1];left=np.zeros(8);records=[]
    for cls,label,color,hatch in [('dual','dual',b.BLUE,''),('A_only','A-only',b.ORANGE,'///'),('B_only','B-only',b.GREEN,'\\\\'),('neither','neither','#B7BEC7','..')]:
        vals=[sum(r['pair']==p and r['class']==cls for r in pop) for p in b.ORDER]
        ax[0].barh(y,vals,left=left,color=color,edgecolor='white',linewidth=.6,hatch=hatch,label=label,height=.68)
        for yy,v,bb,p in zip(y,vals,left,b.ORDER):
            ax[0].text(bb+v/2,yy,str(v),ha='center',va='center',fontsize=8,
                color='black',bbox={'facecolor':'white','edgecolor':'none','pad':.45})
            records.append({'panel':'a','pair':p,'class':cls,'n':v})
        left+=vals
    ax[0].set_yticks(y,b.ORDER);ax[0].set_xlim(0,115)
    ax[0].set_xticks([0,25,50,75,100]);ax[0].set_xlabel('Pair-ligand observations (n)')
    b.axes_style(ax[0]);title(ax[0],'a','Evaluation-set composition')
    ax[0].legend(loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2,frameon=False,
                 columnspacing=1.2,handlelength=1.5)
    fields=['both_ends_raw_activity','literature_complete','cross_source_or_endpoint_mix','replicate_threshold_discordance']
    nn=np.array([sum(r['pair']==p for r in cross) for p in b.ORDER])
    counts=np.array([[sum(r['pair']==p and r[k]=='YES' for r in cross) for k in fields] for p in b.ORDER])
    frac=counts/nn[:,None]
    im=ax[1].imshow(frac,vmin=0,vmax=1,cmap='cividis',aspect='auto')
    for i,p in enumerate(b.ORDER):
        for j,k in enumerate(fields):
            ax[1].text(j,i,str(counts[i,j]),ha='center',va='center',fontsize=8,
                          color=contrast_text(im.cmap(im.norm(frac[i,j]))[:3]))
            records.append({'panel':'b','pair':p,'flag':k,'n':int(counts[i,j]),'denominator':int(nn[i]),'fraction':float(frac[i,j])})
    ax[1].set_yticks(range(8),[])
    ax[1].set_xticks(range(4),['Both-end\nrecords','Complete\nIDs','Source /\nendpoint\nmix','Threshold\ncrossing'])
    ax[1].tick_params(axis='x',length=0,pad=7)
    title(ax[1],'b','Activity-source coverage')
    cb=fig.colorbar(im,cax=fig.add_axes([.60,.105,.35,.025]),orientation='horizontal')
    cb.set_ticks([0,.5,1]);cb.set_label('Fraction within target pair',fontsize=8)
    assert len(pop)==867 and sum(left)==867 and np.array_equal(nn,left)
    finish(fig,'1','Evaluation_Set','Results 3.1','Tables 1, S1B, S4',
      ['results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv','results/current/provenance_audit_20261023/CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv'],records,
      'Experimental classes and activity-source coverage. (a) Numbers of dual, A-only, B-only and neither observations in each target pair. (b) Counts with activity records on both targets, complete document identifiers, mixed sources or endpoints, and threshold-discordant repeated measurements. Cell color gives the fraction within each pair; labels give counts. Coverage flags overlap. Complete identifiers require DOI, PubMed or patent identifiers for all tied activity maxima at both endpoints. The 867 observations represent 845 distinct ligand entities','')

def task_rows(ax):
    y=np.arange(16)[::-1]
    ax.set_yticks(y,['A-only','B-only']*8)
    ax.tick_params(axis='y',length=0,pad=7)
    for i,p in enumerate(b.ORDER):
        ax.text(-.30,14.5-2*i,p,transform=ax.get_yaxis_transform(),
                ha='right',va='center',fontsize=8)
    for a in [ax]:a.set_ylim(-.7,15.7)

def control(number,pocket=False):
    source=b.ST+('module_C_wrong_pocket.csv' if pocket else 'module_B_negative_class.csv')
    data=b.read(source);arm='arm' if pocket else 'direction'
    by={(r['pair'],r[arm]):r for r in data};rr=[by[u] for u in b.units()]
    one,two,delta,lo,hi=('auc_correct','auc_wrong','delta_pocket','delta_pocket_ci_lo','delta_pocket_ci_hi') if pocket else ('auc_single','auc_neither','delta_negative','delta_negative_ci_lo','delta_negative_ci_hi')
    fig,ax=plt.subplots(1,2,figsize=(WIDTH,6.6),sharey=True,gridspec_kw={'width_ratios':[1.06,1]})
    fig.subplots_adjust(left=.34,right=.985,bottom=.16,top=.92,wspace=.22)
    y=np.arange(16)[::-1]
    for yy,r in zip(y,rr):ax[0].plot([b.f(r,one),b.f(r,two)],[yy,yy],color='#AAB1B8',lw=.8)
    ax[0].scatter([b.f(r,one) for r in rr],y,s=26,color=b.BLUE,marker='o',label='Corresponding pocket' if pocket else 'Single-target negatives',zorder=3)
    ax[0].scatter([b.f(r,two) for r in rr],y,s=26,color=b.ORANGE,marker='s',label='Noncorresponding pocket' if pocket else 'Neither negatives',zorder=3)
    b.interval(ax[1],[b.f(r,delta) for r in rr],[b.f(r,lo) for r in rr],[b.f(r,hi) for r in rr],y)
    task_rows(ax[0]);ax[0].set_xlim(0,1);ax[0].set_xticks([0,.5,1]);ax[0].set_xlabel('AUROC')
    for a in ax:b.band_pairs(a);a.set_ylim(-.7,15.7)
    ax[1].set_xlim((-.32,.45) if pocket else (-.25,.75))
    ax[1].set_xlabel('Corresponding minus\nnoncorresponding' if pocket else 'Neither minus\nsingle-target')
    b.axes_style(ax[0],reference=.5);b.axes_style(ax[1],reference=0)
    title(ax[0],'a','AUROC');title(ax[1],'b','AUROC difference')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.65,.015),frameon=False,ncol=1)
    fig.text(.04,.955,'Target pair',fontsize=8,fontweight='bold')
    fig.text(.23,.955,'Negative class',fontsize=8,fontweight='bold')
    caption=('Corresponding and noncorresponding pocket controls. (a) AUROC on the same members using corresponding and noncorresponding target scores. The corresponding target is B for dual versus A-only and A for dual versus B-only. (b) Corresponding minus noncorresponding AUROC with 95% paired bootstrap confidence intervals. All scores are Vina S1. Dual is the positive class in every row' if pocket else 'Effect of negative-class composition. (a) AUROC using single-target or neither negatives; the dual members and score target are held fixed within each comparison. (b) Neither minus single-target AUROC with 95% bootstrap confidence intervals. The two negative classes are resampled separately with shared dual resamples. All scores are Vina S1. Row labels identify the single-target negative class: A-only uses target B scores; B-only uses target A scores')
    finish(fig,str(number),'','Results 3.5' if pocket else 'Results 3.3','Table S6B' if pocket else 'Table S6A',[source],data,caption,'')

def fig5():
    source=b.ST+'module_D_chemistry_units.csv';data=b.read(source)
    by={(r['pair'],r['arm']):r for r in data};rr=[by[u] for u in b.units()];y=np.arange(16)[::-1]
    fig,ax=plt.subplots(1,2,figsize=(WIDTH,6.6),sharey=True,gridspec_kw={'width_ratios':[1.15,1]})
    fig.subplots_adjust(left=.34,right=.985,bottom=.18,top=.92,wspace=.23)
    for key,label,c,m,off in [('auroc_physchem','Descriptors',b.BLUE,'o',.18),('auroc_ecfp4','ECFP4',b.ORANGE,'s',.06),('auroc_ecfp4_m0','ECFP4 + S1',b.GREEN,'^',-.06),('auroc_SUPPORTING_CHEMISTRY_MATCHED_M0','Vina S1',b.PURPLE,'D',-.18)]:
        ax[0].scatter([b.f(r,key) for r in rr],y+off,s=22,color=c,marker=m,label=label,zorder=3)
    ax[1].scatter([b.f(r,'delta_ecfp4_plus_m0') for r in rr],y,s=28,color=b.GREEN,marker='^')
    task_rows(ax[0]);ax[0].set_xlim(0,1);ax[0].set_xticks([0,.5,1]);ax[0].set_xlabel('AUROC')
    ax[1].set_xlim(-.025,.045);ax[1].set_xticks([-.02,0,.02,.04]);ax[1].set_xlabel('AUROC difference')
    for a in ax:b.band_pairs(a);a.set_ylim(-.7,15.7)
    b.axes_style(ax[0],reference=.5);b.axes_style(ax[1],reference=0)
    title(ax[0],'a','AUROC');title(ax[1],'b','Addition of the Vina score')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.65,.02),ncol=2,frameon=False)
    fig.text(.04,.955,'Target pair',fontsize=8,fontweight='bold');fig.text(.23,.955,'Negative class',fontsize=8,fontweight='bold')
    finish(fig,'5','','Results 3.4','Table S7',[source],data,
      'Ligand chemistry and the incremental contribution of Vina. (a) Descriptor, ECFP4 and ECFP4-plus-S1 models use five-fold out-of-fold predictions. Vina S1 is a direct score on identical members within each task. (b) AUROC(ECFP4 plus S1) minus AUROC(ECFP4). Points show estimates. Confidence intervals were not calculated for these chemistry comparisons. Dual is the positive class; A-only rows use target B scores and B-only rows use target A scores','')

style_original=b.style
b.style=style
b.title=title
b.finish=finish
b.fig1=fig1
b.control=control
b.fig5=fig5

def main():
    for folder in [MAIN,SI,QA,DATA]:folder.mkdir(parents=True,exist_ok=True)
    # Hash scientific inputs before any presentation output is written.
    protected={str(p.relative_to(ROOT)):sha(p) for parent in ['results','extension_jcim_validation','reruns'] for p in (ROOT/parent).rglob('*') if p.is_file()}
    # Avoid recursive style replacement.
    b.style=style_original
    family=style()
    b.fig1();b.fig2();b.fig3();b.control(4);b.fig5();b.control(6,True)
    b.fig7();b.fig8();b.fig9();b.si1();b.si2();b.si3()
    b.descriptor_figure('S4','D_vs_A','A_only');b.descriptor_figure('S5','D_vs_B','B_only')
    b.si6();b.si7();b.si8();b.si9();b.scheme()
    for p,digest in protected.items():
        if sha(ROOT/p)!=digest:raise ValueError('Scientific input modified: '+p)
    (QA/'Protected_Inputs.json').write_text(json.dumps(protected,indent=2)+'\n')
    (QA/'Artwork_QA.json').write_text(json.dumps({'base_commit':BASE_COMMIT,'font':family,
        'scientific_recomputation':False,'protected_files_unchanged':len(protected),
        'figures':audit},indent=2)+'\n')
    with (OUT/'Figure_Index.csv').open('w',newline='',encoding='utf8') as h:
        w=csv.writer(h);w.writerow(['figure','group','file','citation','tables','sources'])
        for r in figures:w.writerow([r['number'],r['group'],('Main_Figures/' if r['group']=='Main' else 'Supplementary_Figures/')+r['stem'],r['section'],r['table'],';'.join(r['sources'])])
    caps=[]
    for r in figures:
        label='Scheme S1' if r['number']=='Scheme_S1' else 'Fig. '+r['number']
        caps+=['**'+label+'** '+r['caption'],'']
    (OUT/'Figure_Captions_EN.md').write_text('\n'.join(caps),encoding='utf8')
    review=fitz.open()
    for r in figures:
        source=fitz.open((SI if r['group']=='SI' else MAIN)/(r['stem']+'.pdf'))
        review.insert_pdf(source)
    review.set_toc([[1,'Scheme S1' if r['number']=='Scheme_S1' else 'Fig. '+r['number'],i+1] for i,r in enumerate(figures)])
    review.save(OUT/'Figures_Review.pdf');review.close()
    print(json.dumps({'figures':len(figures),'min_font':min(x['min_font_pt'] for x in audit),'protected_files':len(protected)}))

if __name__=='__main__':main()
