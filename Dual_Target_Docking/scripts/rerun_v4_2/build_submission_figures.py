"""Draw submission figures from accepted results. No docking/model/statistic reruns.

Main figures 1-9 follow Results 3.1-3.7. Earlier working graphics are retained
as provenance; the submission figure index is the publication numbering source.
"""
from __future__ import annotations
import csv
import json
import io
import math
import zlib
from pathlib import Path
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch,Patch
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'submission'
MAIN=OUT/'05_Main_Figures'; SI=OUT/'06_SI_Figures'; DATA=OUT/'07_Figure_Source_Data'
ORDER=['EGFR/HER2','JAK1/JAK2','JAK1/TYK2','PIK3CA/mTOR','AChE/BChE','F2/F10','PPARG/PPARA','PPARA/PPARD']
METHODS=['M0','M1','M2','M1b','M3']; LABEL=dict(zip(METHODS,['S1','S2','S3','S4','S5']))
ARMS=['D_vs_A','D_vs_B']; CONTRAST=['dual_vs_A_only','dual_vs_B_only','summary_min']
BLUE='#0072B2'; ORANGE='#D55E00'; GREEN='#009E73'; PURPLE='#9B59B6'; GREY='#6B7280'
MC={'M0':GREY,'M1':ORANGE,'M2':GREEN,'M1b':PURPLE,'M3':BLUE}
MK={'M0':'o','M1':'s','M2':'^','M1b':'D','M3':'v'}
STATE={'CONCORDANT_IMPROVE':(BLUE,'o','Concordant improvement'),'CONCORDANT_WORSEN':(ORANGE,'s','Concordant worsening'),'DISCORDANT':(PURPLE,'D','Discordant'),'NO_CHANGE':(GREY,'^','No change')}
ST='results/jcim_stage2/'; CAP='results/current/chemical_diversity_cap5_20261023/'
PR='results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv'; DELTA='results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv'
FIGURES=[]; QUALITY=[]; SOURCE={}

def read(path):
    if path not in SOURCE:
        with (ROOT/path).open(encoding='utf-8-sig',newline='')as f: SOURCE[path]=list(csv.DictReader(f))
    return SOURCE[path]
def f(row,key):return float(row[key])
def pub_delta(s):return s.replace('M1b','S4').replace('M0','S1').replace('M1','S2').replace('M2','S3').replace('M3','S5')
def style():
    available={x.name for x in font_manager.fontManager.ttflist}
    family=next((s for s in ['Arial','Liberation Sans','DejaVu Sans']if s in available),'DejaVu Sans')
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':[family],'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.5,'axes.linewidth':0.7,'lines.linewidth':0.8,'patch.linewidth':0.7,'xtick.major.width':0.6,'ytick.major.width':0.6,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
    return family
def axes_style(ax,axis='x',reference=None):
    ax.set_axisbelow(True);ax.grid(axis=axis,color='#E5E7EB',linewidth=0.5)
    if reference is not None:
        (ax.axvline if axis=='x' else ax.axhline)(reference,color='#777777',ls=(0,(3,3)),lw=0.7,zorder=0)
def title(ax,letter,text):ax.set_title(letter+'  '+text,loc='left',pad=9,fontsize=9,fontweight='bold')
def labels_units():return [p+'  '+('dual vs A-only' if a=='D_vs_A' else 'dual vs B-only')for p in ORDER for a in ARMS]
def units():return [(p,a)for p in ORDER for a in ARMS]
def band_pairs(ax,n=16):
    for i in range(0,n,4):ax.axhspan(i-0.5,min(i+1.5,n-0.5),color='#F2F4F6',zorder=-2)
def interval(ax,points,lo,hi,y,color=BLUE,marker='o',label=None,open_=False):
    ax.hlines(y,lo,hi,color=color,lw=0.9,zorder=2)
    ax.scatter(points,y,s=24,marker=marker,edgecolors=color,facecolors='white'if open_ else color,linewidths=0.8,zorder=3,label=label)
def export_data(stem,records):
    if not records:return
    keys=list(dict.fromkeys(k for r in records for k in r))
    with(DATA/(stem+'.csv')).open('w',encoding='utf-8',newline='')as h:
        w=csv.DictWriter(h,fieldnames=keys,lineterminator='\n');w.writeheader();w.writerows(records)
def finish(fig,number,name,section,table,sources,records,caption,caption_cn,si=False):
    stem=(number if number.startswith('Scheme') else 'Figure_'+number.zfill(2)if not si else 'Figure_'+number)+'_'+name
    folder=SI if si else MAIN
    fig.canvas.draw();renderer=fig.canvas.get_renderer();fw,fh=fig.canvas.get_width_height()
    outside=[];sizes=[];off_axis_ticks=set()
    for a in fig.axes:
        xmin,xmax=sorted(a.get_xlim());ymin,ymax=sorted(a.get_ylim())
        for loc,tick in zip(a.get_xticks(),a.get_xticklabels()):
            if not xmin-1e-10<=loc<=xmax+1e-10:off_axis_ticks.add(id(tick))
        for loc,tick in zip(a.get_yticks(),a.get_yticklabels()):
            if not ymin-1e-10<=loc<=ymax+1e-10:off_axis_ticks.add(id(tick))
    for t in fig.findobj(Text):
        if not t.get_visible()or not t.get_text().strip()or id(t)in off_axis_ticks:continue
        bb=t.get_window_extent(renderer)
        if bb.width>0 and bb.height>0:
            sizes.append(t.get_fontsize())
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1>fw+1 or bb.y1>fh+1:outside.append(t.get_text())
    for extension in ['pdf','svg','png']:
        buffer=io.BytesIO();fig.savefig(buffer,format=extension,dpi=600);content=buffer.getvalue()
        if extension=='pdf'and not content.rstrip().endswith(b'%%EOF'):raise ValueError('Incomplete PDF export: '+stem)
        destination=folder/(stem+'.'+extension);temporary=destination.with_suffix(destination.suffix+'.tmp')
        temporary.write_bytes(content);temporary.replace(destination)
    export_data(stem,records)
    QUALITY.append({'figure':number,'stem':stem,'width_in':float(fig.get_figwidth()),'height_in':float(fig.get_figheight()),'png_dpi':600,'smallest_visible_font_pt':min(sizes)if sizes else 0,'outside_canvas_text':outside})
    FIGURES.append({'number':number,'stem':stem,'group':'SI'if si else'Main','first_citation':section,'table':table,'sources':sources,'caption':caption,'caption_cn':caption_cn})
    plt.close(fig)

def fig1():
    pop=[r for r in read('results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')if r['activity_eligible']=='1'and r['class']]
    cross=read('results/current/provenance_audit_20261023/CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv')
    fig,ax=plt.subplots(1,2,figsize=(7,4.5),gridspec_kw={'width_ratios':[1.1,1]});fig.subplots_adjust(left=.18,right=.96,bottom=.26,top=.86,wspace=.20)
    yy=np.arange(8)[::-1];bottom=np.zeros(8);records=[]
    for cls,label,color,hatch in [('dual','dual',BLUE,''),('A_only','A-only',ORANGE,'///'),('B_only','B-only',GREEN,'\\\\'),('neither','neither','#B7BEC7','..')]:
        vals=[sum(r['pair']==p and r['class']==cls for r in pop)for p in ORDER]
        ax[0].barh(yy,vals,left=bottom,color=color,edgecolor='white',linewidth=.6,hatch=hatch,label=label,height=.70)
        for y,v,b,p in zip(yy,vals,bottom,ORDER):
            ax[0].text(b+v/2,y,str(v),ha='center',va='center',color='#111111',fontsize=7.5,bbox={'facecolor':'white','edgecolor':'none','pad':.45})
            records.append({'panel':'a','pair':p,'class':cls,'n':v})
        bottom+=vals
    ax[0].set_yticks(yy,ORDER);ax[0].set_xlim(0,115);ax[0].set_xlabel('Pair-ligand observations (n)');axes_style(ax[0]);title(ax[0],'a','Evaluation-set composition')
    ax[0].legend(loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2,frameon=False,handlelength=1.6)
    fields=['both_ends_raw_activity','literature_complete','cross_source_or_endpoint_mix','replicate_threshold_discordance']
    counts=np.array([[sum(r['pair']==p and r[k]=='YES'for r in cross)for k in fields]for p in ORDER]);nn=np.array([sum(r['pair']==p for r in cross)for p in ORDER]);fractions=counts/nn[:,None]
    im=ax[1].imshow(fractions,vmin=0,vmax=1,cmap='cividis',aspect='auto',extent=(-.5,3.5,-.5,7.5))
    for i,p in enumerate(ORDER):
        for j,k in enumerate(fields):
            ax[1].text(j,7-i,str(counts[i,j]),ha='center',va='center',fontsize=8,color='white'if fractions[i,j]<.45 else'black');records.append({'panel':'b','pair':p,'flag':k,'n':int(counts[i,j]),'denominator':int(nn[i]),'fraction':float(fractions[i,j])})
    ax[1].set_yticks(yy,[]);ax[1].set_xticks(range(4),['Both-end\nrecords','Complete\nIDs','Source /\nendpoint mix','Threshold\ncrossing'],rotation=35,ha='right',fontsize=7.5);title(ax[1],'b','Activity-source coverage')
    cax=fig.add_axes([.58,.12,.34,.025]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Fraction within target pair',fontsize=7.5);cb.ax.tick_params(labelsize=7.5,width=.6)
    fig.text(.18,.96,'867 pair-ligand observations; 845 distinct ligand entities',fontsize=9)
    finish(fig,'1','Evaluation_Set','Results 3.1','Tables 1, S1B, S4',['results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv','results/current/provenance_audit_20261023/CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv'],records,'Experimental classes and activity-source coverage. (a) Numbers of dual, A-only, B-only and neither observations in each target pair. (b) Counts of observations with raw records on both targets, complete specified document identifiers, mixed sources or endpoints, and threshold-discordant repeated measurements. Cell color represents the fraction within each pair; text gives counts. Audit flags overlap. Complete identifiers require DOI, PubMed or patent identifiers for all tied activity maxima at both endpoints. The 867 observations represent 845 distinct global ligand entities.','评价集的四类组成及活性来源覆盖。颜色表示各靶点对内比例，格内数字为观察数；来源标记可以重叠。')

def fig2():
    data=[r for r in read(PR)if r['method']=='M0'];by={(r['pair'],r['contrast']):r for r in data}
    fig,ax=plt.subplots(1,2,figsize=(7,4.7),sharey=True,gridspec_kw={'width_ratios':[1.25,1]});fig.subplots_adjust(left=.19,right=.97,bottom=.17,top=.86,wspace=.18);y=np.arange(8)[::-1]
    for con,color,marker,off,label in [('dual_vs_A_only',BLUE,'o',.12,'dual vs A-only (target B)'),('dual_vs_B_only',ORANGE,'s',-.12,'dual vs B-only (target A)')]:
        rr=[by[p,con]for p in ORDER];interval(ax[0],[f(r,'auroc')for r in rr],[f(r,'ci_lo')for r in rr],[f(r,'ci_hi')for r in rr],y+off,color,marker,label)
    rr=[by[p,'summary_min']for p in ORDER];interval(ax[1],[f(r,'auroc')for r in rr],[f(r,'ci_lo')for r in rr],[f(r,'ci_hi')for r in rr],y,GREEN,'D')
    ax[0].set_yticks(y,ORDER)
    for a in ax:a.set_xlim(0,1);a.set_ylim(-.6,7.6);a.set_xlabel('AUROC');axes_style(a,reference=.5)
    title(ax[0],'a','Directional Vina discrimination');title(ax[1],'b','Minimum directional AUROC')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.55,.01),ncol=1,frameon=False)
    finish(fig,'2','Vina_Directional','Results 3.2, first three paragraphs','Table S5A',[PR],data,'Vina directional discrimination and summary_min. Points and horizontal bars show estimates and stored 95% class-stratified percentile bootstrap intervals. (a) Dual versus A-only uses target B; dual versus B-only uses target A. (b) summary_min is the smaller directional AUROC, calculated within each resample. The vertical line marks AUROC 0.5 as a reference. The summary interval is not a simultaneous test that both directions exceed or fall below 0.5.','Vina 两个方向的 AUROC 和 summary_min。点与横线分别为点估计及已有的 95% bootstrap 区间，虚线为 0.5 参考值。')

def fig3():
    data=read(DELTA); fig,axs=plt.subplots(2,2,figsize=(7,6.3),sharex=True,sharey=True);fig.subplots_adjust(left=.18,right=.97,bottom=.12,top=.92,wspace=.18,hspace=.33);yy=np.arange(8)[::-1]
    for a,d,letter in zip(axs.ravel(),['M1-M0','M2-M0','M1b-M1','M3-M0'],'abcd'):
        rr=[next(r for r in data if r['pair']==p and r['delta']==d)for p in ORDER];m=d.split('-')[0]
        interval(a,[f(r,'delta_summary_min')for r in rr],[f(r,'ci_lo')for r in rr],[f(r,'ci_hi')for r in rr],yy,MC[m],MK[m]);a.set_yticks(yy,ORDER);a.set_xlim(-.55,.55);a.set_ylim(-.6,7.6);axes_style(a,reference=0);title(a,letter,pub_delta(d))
    for a in axs[1]:a.set_xlabel('Change in summary_min')
    finish(fig,'3','Paired_Scoring_Comparisons','Results 3.2, scoring workflow paragraphs','Table S5B',[DELTA],data,'Paired changes in summary_min between scoring workflows. (a) S2-S1, fixed-pose GNINA versus Vina. (b) S3-S1, fixed-pose RTMScore versus Vina. (c) S4-S2, selection of the highest CNNscore among saved Vina poses versus the original mode-1 pose. (d) S5-S1, independent GNINA docking versus Vina. Every contrast uses its own complete-case paired members. Bars are stored 95% paired bootstrap intervals; no across-pair pooled estimate is calculated.','四类评分工作流的配对差值。每个比较使用自身的共同有效成员，横线为已有的 95% 配对 bootstrap 区间。')

def control(number,pocket=False):
    source=ST+('module_C_wrong_pocket.csv'if pocket else'module_B_negative_class.csv');data=read(source);arm='arm'if pocket else'direction';by={(r['pair'],r[arm]):r for r in data};rr=[by[u]for u in units()]
    one,two,delta,low,high=('auc_correct','auc_wrong','delta_pocket','delta_pocket_ci_lo','delta_pocket_ci_hi')if pocket else('auc_single','auc_neither','delta_negative','delta_negative_ci_lo','delta_negative_ci_hi')
    fig,ax=plt.subplots(1,2,figsize=(7,6.6),sharey=True,gridspec_kw={'width_ratios':[1.12,1]});fig.subplots_adjust(left=.33,right=.97,bottom=.15,top=.92,wspace=.20);y=np.arange(16)[::-1]
    for yy,r in zip(y,rr):ax[0].plot([f(r,one),f(r,two)],[yy,yy],color='#B6BDC5',lw=.7,zorder=1)
    ax[0].scatter([f(r,one)for r in rr],y,s=23,color=BLUE,marker='o',label='Corresponding pocket'if pocket else'Single-target negative',zorder=3)
    ax[0].scatter([f(r,two)for r in rr],y,s=23,color=ORANGE,marker='s',label='Noncorresponding pocket'if pocket else'Neither negative',zorder=3)
    interval(ax[1],[f(r,delta)for r in rr],[f(r,low)for r in rr],[f(r,high)for r in rr],y)
    ax[0].set_yticks(y,labels_units(),fontsize=7.5);ax[0].set_xlim(0,1);ax[0].set_xlabel('AUROC');axes_style(ax[0],reference=.5)
    for a in ax:a.set_ylim(-.6,15.6);band_pairs(a)
    ax[1].set_xlim((-.32,.45)if pocket else(-.25,.75));axes_style(ax[1],reference=0);ax[1].set_xlabel('Corresponding minus\nnoncorresponding AUROC'if pocket else'Neither - single-target')
    title(ax[0],'a','Same-member AUROC');title(ax[1],'b','Paired AUROC difference')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.60,.015),frameon=False,ncol=1)
    caption=('Corresponding and noncorresponding pocket controls. (a) AUROC values on the same members using corresponding and noncorresponding target scores. The corresponding target is B for dual versus A-only and A for dual versus B-only. (b) Corresponding minus noncorresponding AUROC with stored 95% paired bootstrap intervals. All scores are Vina S1. No additional member is lost relative to the corresponding directional S1 population.'if pocket else'Effect of negative-class composition. (a) AUROC using single-target-selective or neither negatives, with the dual set and score target held fixed within each directional comparison. (b) Neither minus single-target AUROC with stored 95% bootstrap intervals. All scores are Vina S1. Directional counts are listed in Table S6A; neither and single-target negatives are sampled separately while sharing dual resamples.')
    finish(fig,str(number),'Pocket_Correspondence'if pocket else'Negative_Class', 'Results 3.5'if pocket else'Results 3.3','Table S6B'if pocket else'Table S6A',[source],data,caption,'对应与非对应口袋的同成员对照及配对差值。'if pocket else'单靶选择性与 neither 负类的 AUROC 对照及差值，固定每个比较中的 dual 成员与评分靶点。')

def fig5():
    source=ST+'module_D_chemistry_units.csv';data=read(source);by={(r['pair'],r['arm']):r for r in data};rr=[by[u]for u in units()];y=np.arange(16)[::-1]
    fig,ax=plt.subplots(1,2,figsize=(7,6.6),sharey=True,gridspec_kw={'width_ratios':[1.15,1]});fig.subplots_adjust(left=.33,right=.97,bottom=.17,top=.92,wspace=.20)
    for key,label,color,marker,offset in [('auroc_physchem','Descriptors',BLUE,'o',.18),('auroc_ecfp4','ECFP4',ORANGE,'s',.06),('auroc_ecfp4_m0','ECFP4 + S1',GREEN,'^',-.06),('auroc_SUPPORTING_CHEMISTRY_MATCHED_M0','Vina S1',PURPLE,'D',-.18)]:ax[0].scatter([f(r,key)for r in rr],y+offset,s=20,marker=marker,color=color,label=label,zorder=3)
    ax[1].scatter([f(r,'delta_ecfp4_plus_m0')for r in rr],y,s=25,color=GREEN,marker='^',zorder=3)
    ax[0].set_yticks(y,labels_units(),fontsize=7.5);ax[0].set_xlim(0,1);ax[0].set_xlabel('AUROC on common\nchemistry members');ax[1].set_xlim(-.025,.045);ax[1].set_xticks([-.02,0,.02,.04]);ax[1].set_xlabel('Change after adding S1\nECFP4 + S1 minus ECFP4')
    for a in ax:a.set_ylim(-.6,15.6);band_pairs(a)
    axes_style(ax[0],reference=.5);axes_style(ax[1],reference=0);title(ax[0],'a','Chemistry models and Vina');title(ax[1],'b','Addition of the Vina score')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.60,.005),frameon=False,ncol=2)
    finish(fig,'5','Chemistry_and_Vina_Increment','Results 3.4','Table S7',[source],data,'Ligand chemistry and the incremental contribution of Vina. (a) Descriptor, ECFP4 and ECFP4-plus-S1 models use existing out-of-fold predictions. Vina S1 is a direct score evaluated on the same members and is not an out-of-fold learned model. (b) AUROC(ECFP4 plus S1) minus AUROC(ECFP4). All 16 directional tasks use matched chemistry members. Points are estimates; no confidence interval for the chemistry-model increment was calculated.','配体化学模型与 Vina 的共同成员比较。三个化学模型使用折外预测，单独 Vina 为直接评分；图中增量为点估计，没有新增置信区间。')

def fig7():
    seed=read(ST+'module_E_five_seed.csv');summ=read(ST+'module_E_five_seed_summary.csv');alt=read(ST+'module_F_alt_receptor.csv');by={(r['pair'],int(r['seed'])):r for r in seed};ss={r['pair']:r for r in summ if r['metric']=='summary_min'}
    fig,ax=plt.subplots(2,1,figsize=(7,7.1),gridspec_kw={'height_ratios':[1.15,1]});fig.subplots_adjust(left=.29,right=.97,bottom=.13,top=.91,hspace=.54)
    yy=np.arange(8)[::-1];colors=[BLUE,ORANGE,GREEN,PURPLE,'#555555'];markers=['o','s','^','D','v']
    for seedno,off,color,m in zip([17,29,42,71,101],np.linspace(-.16,.16,5),colors,markers):ax[0].scatter([f(by[p,seedno],'summary_min')for p in ORDER],yy+off,s=22,color=color,marker=m,label=str(seedno),zorder=3)
    ax[0].hlines(yy,[f(ss[p],'min')for p in ORDER],[f(ss[p],'max')for p in ORDER],color='#222222',lw=.9,zorder=2);ax[0].scatter([f(ss[p],'median')for p in ORDER],yy,s=80,color='black',marker='|',linewidths=1,zorder=4,label='Median')
    ax[0].set_yticks(yy,ORDER);ax[0].set_xlim(0,1);ax[0].set_ylim(-.6,7.6);ax[0].set_xlabel('summary_min');axes_style(ax[0],reference=.5);title(ax[0],'a','Vina random seeds')
    fig.legend(*ax[0].get_legend_handles_labels(),frameon=False,ncol=6,loc='upper center',bbox_to_anchor=(.63,.99),handletextpad=.3,columnspacing=.8)
    rr=sorted(alt,key=lambda r:(ORDER.index(r['pair']),r['replaced_side']));yy=np.arange(6)[::-1]
    for field,label,color,marker,off in [('delta_AUC_affected','Affected direction',BLUE,'o',.13),('delta_summary_min','summary_min',PURPLE,'D',-.13)]:interval(ax[1],[f(r,field)for r in rr],[f(r,field+'_ci_lo')for r in rr],[f(r,field+'_ci_hi')for r in rr],yy+off,color,marker,label)
    ax[1].set_yticks(yy,[r['pair']+'  '+r['alt_id'].replace('ALT_','')+' ('+r['replaced_side']+')'for r in rr],fontsize=7.5);ax[1].set_xlim(-.24,.20);ax[1].set_ylim(-.6,5.6);ax[1].set_xlabel('Alternative - original receptor');axes_style(ax[1],reference=0);title(ax[1],'b','One-target receptor substitutions');ax[1].legend(frameon=False,loc='lower center',bbox_to_anchor=(.5,-.40),ncol=2)
    records=[{'panel':'a',**r}for r in seed]+[{'panel':'a_summary',**r}for r in summ]+[{'panel':'b',**r}for r in rr]
    finish(fig,'7','Seed_and_Receptor_Sensitivity','Results 3.6, first two paragraphs','Table S8A,B',[ST+'module_E_five_seed.csv',ST+'module_E_five_seed_summary.csv',ST+'module_F_alt_receptor.csv'],records,'Sensitivity to search seeds and receptor structures. (a) Five seed-specific Vina summary_min estimates. Short horizontal segments show the stored minimum-to-maximum range; black ticks show medians. These segments are not confidence intervals. Each seed uses its valid members; AChE/BChE seed 17 loses one A-only ligand relative to seed 42. (b) Paired changes in affected-direction AUROC and summary_min for six single-target receptor replacements. Parentheses identify replaced side A or B; the other receptor is held fixed. Bars in (b) are stored 95% paired bootstrap intervals.','随机种子与单端受体替换敏感性。a 中短线表示五种子的范围，非置信区间；b 中横线为受体比较已有的 95% 配对 bootstrap 区间。')

def fig8():
    primary=read(PR);cap=read(CAP+'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv');strict=read('extension_jcim_validation/pik3ca_mtor_q110_v42/Q110_STRICT_SUBSET_METRICS.csv');counts=read(CAP+'CAP5_CLASS_COUNTS.csv')
    fig,ax=plt.subplots(1,2,figsize=(7,4.8),gridspec_kw={'width_ratios':[1.2,1]});fig.subplots_adjust(left=.23,right=.97,bottom=.23,top=.86,wspace=.35);y=np.arange(8)[::-1]
    for rows,color,marker,offset,label,open_ in [(primary,GREY,'o',.12,'Primary',True),(cap,PURPLE,'D',-.12,'CAP5',False)]:
        by={r['pair']:r for r in rows if r['method']=='M0'and r['contrast']=='summary_min'};rr=[by[p]for p in ORDER];interval(ax[0],[f(r,'auroc')for r in rr],[f(r,'ci_lo')for r in rr],[f(r,'ci_hi')for r in rr],y+offset,color,marker,label,open_)
    nn={p:(sum(int(r['n_original'])for r in counts if r['pair']==p),sum(int(r['n_cap5'])for r in counts if r['pair']==p))for p in ORDER}
    ax[0].set_yticks(y,[p+'\n'+str(nn[p][0])+' to '+str(nn[p][1])for p in ORDER],fontsize=7.5);ax[0].set_ylim(-.7,7.7);title(ax[0],'a','CAP5: Vina summary_min')
    y=np.arange(5)[::-1]
    for rows,color,marker,off,label,open_ in [(primary,GREY,'o',.10,'Q110 (n=110)',True),(strict,GREEN,'s',-.10,'Strict subset (n=56)',False)]:
        by={r['method']:r for r in rows if r['pair']=='PIK3CA/mTOR'and r['contrast']=='summary_min'};rr=[by[m]for m in METHODS];interval(ax[1],[f(r,'auroc')for r in rr],[f(r,'ci_lo')for r in rr],[f(r,'ci_hi')for r in rr],y+off,color,marker,label,open_)
    ax[1].set_yticks(y,[LABEL[m]for m in METHODS]);ax[1].set_ylim(-.5,4.5);title(ax[1],'b','Strict-activity subset')
    for a in ax:a.set_xlim(0,1);axes_style(a,reference=.5);a.set_xlabel('summary_min');a.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=1)
    records=[{'population':'Primary',**r}for r in primary if r['method']=='M0'or r['pair']=='PIK3CA/mTOR']+[{'population':'CAP5',**r}for r in cap if r['method']=='M0']+[{'population':'Strict subset',**r}for r in strict]
    finish(fig,'8','Membership_Sensitivity','Results 3.6, CAP5 and strict-subset paragraphs','Table S9A-D',[PR,CAP+'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv','extension_jcim_validation/pik3ca_mtor_q110_v42/Q110_STRICT_SUBSET_METRICS.csv'],records,'Sensitivity to evaluation-set membership. (a) Vina summary_min for the primary and CAP5 populations; labels give original and retained pair sizes. CAP5 retains 861 of 867 observations. (b) summary_min for the five scoring schemes in the 110-member PIK3CA/mTOR set and its 56-member strict-activity subset. Population-size labels include all four experimental classes; directional sample sizes are reported in Tables S5A and S9C. Bars are stored 95% bootstrap intervals for each population; no paired interval for the between-population change is shown. Full directional CAP5 comparisons, including the PPARA/PPARD S5 dual-versus-A-only interval-boundary change, are shown in Figure S8. Chemistry sensitivity is tabulated in Table S9D.','评价成员组成敏感性。a 为 CAP5 前后的 Vina 指标，b 为 PIK3CA/mTOR 主评价集与严格子集的五种评分。不同人群的原有区间分别展示，未新增人群差值区间。')

def place_labels(ax,points):
    ax.figure.canvas.draw();renderer=ax.figure.canvas.get_renderer();used=[]
    from matplotlib.transforms import Bbox
    point_boxes=[]
    for _,x,y in points:
        px,py=ax.transData.transform((x,y));point_boxes.append(Bbox.from_bounds(px-3,py-3,6,6))
    candidates=[(7,7),(-8,7),(7,-10),(-8,-10),(13,1),(-14,1),(1,15),(1,-17),(18,14),(-20,-14)]
    for label,x,y in points:
        selected=None
        for dx,dy in candidates:
            t=ax.annotate(str(label),(x,y),xytext=(dx,dy),textcoords='offset points',fontsize=7.5,ha='center',va='center',zorder=5)
            bb=t.get_window_extent(renderer).expanded(1.12,1.15)
            collision=any(bb.overlaps(b)for b in used+point_boxes);inside=ax.bbox.contains(bb.x0,bb.y0)and ax.bbox.contains(bb.x1,bb.y1)
            if not collision and inside:selected=(t,bb,dx,dy);break
            t.remove()
        if selected is None:
            t=ax.annotate(str(label),(x,y),xytext=(7,7),textcoords='offset points',fontsize=7.5,ha='center',va='center',zorder=5);bb=t.get_window_extent(renderer);selected=(t,bb,7,7)
        used.append(selected[1])
        if abs(selected[2])+abs(selected[3])>20:ax.annotate('',(x,y),xytext=(selected[2]*.75,selected[3]*.75),textcoords='offset points',arrowprops={'arrowstyle':'-','color':'#777777','lw':.5},zorder=2)

def fig9():
    source=ST+'module_I_method_compare.csv';data=read(source);by={(r['pair'],r['method']):r for r in data}
    fig,axs=plt.subplots(2,2,figsize=(7,7.2),sharex=True,sharey=True);fig.subplots_adjust(left=.13,right=.97,bottom=.25,top=.91,wspace=.25,hspace=.42)
    for a,m,letter in zip(axs.ravel(),METHODS[1:],'abcd'):
        labels=[]
        for i,p in enumerate(ORDER,1):
            r=by[p,m];x=f(r,'delta_summary_min');y=100*f(r,'delta_top10_single_target_fraction');color,marker,_=STATE[r['concordance']]
            a.scatter(x,y,s=36,marker=marker,edgecolors=color,facecolors='white'if m=='M1b'else color,linewidths=.9,zorder=3);labels.append((i,x,y))
        a.axhline(0,color='#777777',lw=.7);a.axvline(0,color='#777777',lw=.7);a.set_xlim(-.39,.34);a.set_ylim(-39,58);a.set_xticks([-.3,0,.3]);a.set_yticks([-25,0,25,50]);axes_style(a);title(a,letter,LABEL[m]+' vs S1'+(' (supporting)'if m=='M1b'else''));place_labels(a,labels)
    for a in axs[1]:a.set_xlabel('Change in summary_min')
    for a in axs[:,0]:a.set_ylabel('Single-target fraction change\n(percentage points)')
    handles=[Line2D([],[],color=c,marker=mk,linestyle='none',label=label,markersize=5)for c,mk,label in STATE.values()];fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.55,.12),ncol=2,frameon=False)
    fig.text(.12,.065,'1 EGFR/HER2   2 JAK1/JAK2   3 JAK1/TYK2   4 PIK3CA/mTOR',fontsize=7.5)
    fig.text(.12,.035,'5 AChE/BChE   6 F2/F10   7 PPARG/PPARA   8 PPARA/PPARD',fontsize=7.5)
    finish(fig,'9','Directional_and_Joint_Ranking','Results 3.7','Tables 2, S10B,C',[source],data,'Relationship between directional discrimination and candidate-list composition. Each point is one target pair. Both axes show the new scheme minus S1. The horizontal axis is the paired directional summary_min change; the vertical axis is the change in single-target-selective fraction among the top 10% of joint ranks, expressed in percentage points. Ranking and directional changes use their respective complete-case populations, as defined in Methods. Colors and symbols identify the four change categories. Open symbols identify the supporting S4-S1 comparison. No new confidence interval is added. Matched top-list counts are shown in Table 2 and Figure S9.','方向性指标变化与联合名单组成变化的关系。纵轴为百分点，不是相对百分比；每个方法比较分别使用自己的共同成员。空心点表示 S4 对 S1 的支持性比较。')

def pair_key(fig,bottom=.015):
    fig.text(.12,bottom+.025,'1 EGFR/HER2   2 JAK1/JAK2   3 JAK1/TYK2   4 PIK3CA/mTOR',fontsize=7.5)
    fig.text(.12,bottom,'5 AChE/BChE   6 F2/F10   7 PPARG/PPARA   8 PPARA/PPARD',fontsize=7.5)

def si1():
    source='reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/redocking_master.csv';data=read(source);reg=read('reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/receptor_registry.csv');reg=list({r['pdb_id']:r for r in reg}.values());names=[r['pdb_id']for r in reg]
    fig,ax=plt.subplots(1,2,figsize=(7,6.0),sharey=True);fig.subplots_adjust(left=.23,right=.97,bottom=.16,top=.86,wspace=.18);yy=np.arange(14)[::-1]
    colors=[BLUE,ORANGE,GREEN,PURPLE,GREY];markers=['o','s','^','D','v']
    for a,key,letter,label in zip(ax,['top1_rmsd_A','min_saved_pose_rmsd_A'],'ab',['Top-ranked pose','Best saved pose']):
        for seed,off,c,m in zip([17,29,42,71,101],np.linspace(-.16,.16,5),colors,markers):
            rr=[next(r for r in data if r['pdb_id']==p and int(r['seed'])==seed)for p in names]
            a.scatter([f(r,key)for r in rr],yy+off,s=20,color=c,marker=m,label=str(seed))
        a.set_xlim((0,11.8)if key=='top1_rmsd_A'else(0,3));a.set_ylim(-.6,13.6);a.set_xlabel('Heavy-atom RMSD (angstrom)');axes_style(a,reference=2);title(a,letter,label)
    ax[0].set_yticks(yy,[r['target']+' '+r['pdb_id']for r in reg]);fig.legend(*ax[0].get_legend_handles_labels(),ncol=5,frameon=False,loc='lower center',bbox_to_anchor=(.58,.025))
    finish(fig,'S1','Cocrystal_Redocking','Methods 2.3','Table S3',[source],data,'Cocrystal-ligand redocking for 14 primary receptors and five search seeds. (a) Top-ranked-pose RMSD. (b) Minimum RMSD among all saved poses. RMSDs are symmetry-aware heavy-atom distances computed in the receptor frame. The dashed line marks 2 angstrom. Each point is a recorded seed result; no coordinate or RMSD is recomputed for this figure. Axes differ between panels to preserve visibility of large top-ranked-pose errors.','14 个正式受体的五种子重对接，分别展示首位姿态及保存姿态中的最佳 RMSD。两个面板的横轴范围不同。',si=True)

def si2():
    source=ST+'module_H_missingness.csv';data=read(source);by={(r['pair'],r['method'],r['class'],r['target']):r for r in data}
    fig,axs=plt.subplots(4,2,figsize=(7,7.4));fig.subplots_adjust(left=.08,right=.94,bottom=.15,top=.94,wspace=.18,hspace=.48)
    for a,(cls,target),letter in zip(axs.ravel(),[(c,t)for c in ['dual','A_only','B_only','neither']for t in ['A','B']],'abcdefgh'):
        grid=np.array([[f(by[p,m,cls,target],'valid')/f(by[p,m,cls,target],'expected')for m in METHODS]for p in ORDER]);a.imshow(grid,vmin=0,vmax=1,cmap='cividis',aspect='auto')
        for i,p in enumerate(ORDER):
            for j,m in enumerate(METHODS):r=by[p,m,cls,target];a.text(j,i,r['valid']+'/'+r['expected'],ha='center',va='center',fontsize=7.5,color='black'if grid[i,j]>.45 else'white')
        a.set_xticks(range(5),[LABEL[m]for m in METHODS]);a.set_yticks(range(8),[str(i)for i in range(1,9)]);title(a,letter,cls.replace('_','-')+', target '+target)
    pair_key(fig)
    finish(fig,'S2','Score_Availability','Methods 2.5; Results 3.1','Table S4B',[source],data,'Availability of scores by experimental class, target and scoring scheme. Each cell gives valid/expected counts; its color represents that fraction. Pair numbers are identified in the figure. Missing and invalid scores retain their original status, and are not imputed. The S1-S5 order follows the Methods scoring definitions.','按活性类别、评分靶点与评分方案展示有效评分人数；每格为有效数/应有数。',si=True)

def si3():
    data=read(PR);by={(r['pair'],r['method'],r['contrast']):r for r in data}
    fig,axs=plt.subplots(1,3,figsize=(7,4.5),sharey=True);fig.subplots_adjust(left=.19,right=.97,bottom=.22,top=.86,wspace=.15)
    for a,con,letter,text in zip(axs,CONTRAST,'abc',['dual vs A-only','dual vs B-only','summary_min']):
        grid=np.array([[f(by[p,m,con],'auroc')for m in METHODS]for p in ORDER]);im=a.imshow(grid,vmin=0,vmax=1,cmap='cividis',aspect='auto')
        for i in range(8):
            for j in range(5):a.text(j,i,f'{grid[i,j]:.2f}',ha='center',va='center',fontsize=7.5,color='white'if grid[i,j]<.45 else'black')
        a.set_xticks(range(5),[LABEL[m]for m in METHODS]);a.set_yticks(range(8),ORDER);title(a,letter,text)
    cb=fig.colorbar(im,cax=fig.add_axes([.31,.08,.5,.025]),orientation='horizontal');cb.set_label('AUROC',fontsize=8)
    finish(fig,'S3','All_Scoring_Schemes','Results 3.2','Table S5A',[PR],data,'Directional AUROC and summary_min for all five scoring schemes. Cell annotations show estimates rounded to two decimals; colors use a common 0-1 scale. Intervals and sample sizes are given in Table S5A. The heat map does not replace paired method comparisons, which use their own shared members.','五种评分方案的两个方向 AUROC 与 summary_min；图中仅为点估计，完整区间见 Table S5A。',si=True)

DESCRIPTORS=[('MolWt','Molecular weight (Da)'),('MolLogP','Calculated LogP'),('TPSA','TPSA (square angstrom)'),('NumHDonors','H-bond donors'),('NumHAcceptors','H-bond acceptors'),('NumRotatableBonds','Rotatable bonds'),('RingCount','Ring count'),('FormalCharge','Formal charge')]
def jitter(label):return ((zlib.crc32(label.encode())%10001)/10000-.5)*.22
def descriptive_box(ax,values,pos,horizontal=False):
    if not values:return
    a=np.array(values);q1,med,q3=np.percentile(a,[25,50,75],method='linear');iqr=q3-q1;inside=a[(a>=q1-1.5*iqr)&(a<=q3+1.5*iqr)];lo=float(inside.min())if len(inside)else q1;hi=float(inside.max())if len(inside)else q3
    stats=[{'med':med,'q1':q1,'q3':q3,'whislo':lo,'whishi':hi,'fliers':[]}]
    ax.bxp(stats,positions=[pos],widths=.30,vert=not horizontal,showfliers=False,manage_ticks=False,boxprops={'linewidth':.7},whiskerprops={'linewidth':.7},capprops={'linewidth':.7},medianprops={'linewidth':.9,'color':'#111111'},zorder=1)
def descriptor_figure(number,arm,negative):
    source=ST+'module_D1_composition.csv';allrows=read(source);data=[r for r in allrows if r['arm']==arm]
    fig,axs=plt.subplots(4,2,figsize=(7,7.4));fig.subplots_adjust(left=.10,right=.97,bottom=.15,top=.89,wspace=.22,hspace=.58)
    for a,(key,label),letter in zip(axs.ravel(),DESCRIPTORS,'abcdefgh'):
        for i,p in enumerate(ORDER):
            for cls,off,color,marker in [('dual',-.16,BLUE,'o'),(negative,.16,ORANGE,'^')]:
                rr=[r for r in data if r['pair']==p and r['class']==cls];v=[f(r,key)for r in rr];descriptive_box(a,v,i+off)
                a.scatter([i+off+jitter(r['canonical_ligand_id'])for r in rr],v,s=7,marker=marker,color=color,alpha=.55,linewidths=0,zorder=2)
        a.set_xticks(range(8),[str(i)for i in range(1,9)]);a.set_xlim(-.6,7.6);title(a,letter,label);a.tick_params(labelsize=7.5);axes_style(a,axis='y')
    handles=[Line2D([],[],color=BLUE,marker='o',linestyle='none',label='dual'),Line2D([],[],color=ORANGE,marker='^',linestyle='none',label=negative.replace('_','-'))];fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.97),ncol=2,frameon=False);pair_key(fig)
    finish(fig,number,'Descriptors_'+arm,'Results 3.4','Table S7; Supplementary Data',[source],data,'Physicochemical descriptor distributions for dual and '+negative.replace('_','-')+' ligands. Points are labeled-population observations before score and fold filtering. Boxes show the median and interquartile range, with whiskers to the most extreme values within 1.5 IQR. All observations are retained as points, including those beyond whiskers. Pair numbers are identified in the figure. These distributions are descriptive and do not imply independently validated classification.','dual 与 '+negative.replace('_','-')+' 配体的八项理化描述符分布。展示标签层全部观察；点为配体，箱线为中位数、四分位区间及 1.5 IQR 范围。',si=True)

def si6():
    source=ST+'module_D1_composition.csv';data=read(source);fig,axs=plt.subplots(1,2,figsize=(7,4.7),sharey=True);fig.subplots_adjust(left=.19,right=.97,bottom=.18,top=.86,wspace=.18)
    for a,arm,neg,letter in zip(axs,ARMS,['A_only','B_only'],'ab'):
        for i,p in enumerate(ORDER):
            rr=[r for r in data if r['pair']==p and r['arm']==arm];descriptive_box(a,[f(r,'max_ecfp4_tanimoto_other_class')for r in rr],7-i,horizontal=True)
            for cls,c,m in [('dual',BLUE,'o'),(neg,ORANGE,'^')]:
                rows=[r for r in rr if r['class']==cls];a.scatter([f(r,'max_ecfp4_tanimoto_other_class')for r in rows],[7-i+jitter(r['canonical_ligand_id'])for r in rows],s=7,color=c,marker=m,alpha=.55,linewidths=0,label=cls.replace('_','-')if i==0 else None,zorder=2)
        a.set_yticks(range(7,-1,-1),ORDER);a.set_xlim(0,1);a.set_ylim(-.6,7.6);a.set_xlabel('Maximum cross-class ECFP4 Tanimoto');title(a,letter,'dual vs '+neg.replace('_','-'));axes_style(a);a.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2)
    finish(fig,'S6','Cross_Class_Similarity','Results 3.4','Supplementary Data',[source],data,'Structural similarity between the two experimental classes within each directional task. Each point is the maximum ECFP4 Tanimoto similarity to any ligand in the other class of the same pair and direction. Fingerprints use radius 2, 2048 bits, useChirality=False and useFeatures=False. Pooled boxes summarize the two classes together; class identity is shown by point color and shape.','每个配体与同一靶点对、同一方向另一活性类别配体的最大 ECFP4 Tanimoto 相似性。箱线汇总两类配体，点的颜色与形状区分类别。',si=True)

def si7():
    data=read(ST+'module_E_five_seed.csv');summary=read(ST+'module_E_five_seed_summary.csv');by={(r['pair'],int(r['seed'])):r for r in data}
    fig,axs=plt.subplots(3,1,figsize=(7,7.4),sharex=True);fig.subplots_adjust(left=.20,right=.97,bottom=.08,top=.92,hspace=.45)
    for a,metric,letter,text in zip(axs,['auc_B','auc_A','summary_min'],'abc',['dual vs A-only (target B)','dual vs B-only (target A)','summary_min']):
        y=np.arange(8)[::-1]
        for seed,off,c,m in zip([17,29,42,71,101],np.linspace(-.16,.16,5),[BLUE,ORANGE,GREEN,PURPLE,GREY],['o','s','^','D','v']):a.scatter([f(by[p,seed],metric)for p in ORDER],y+off,s=16,color=c,marker=m,label=str(seed))
        ss={r['pair']:r for r in summary if r['metric']==metric};a.hlines(y,[f(ss[p],'min')for p in ORDER],[f(ss[p],'max')for p in ORDER],color='black',lw=.8);a.scatter([f(ss[p],'median')for p in ORDER],y,color='black',marker='|',s=55,linewidths=1)
        a.set_yticks(y,ORDER,fontsize=7.5);a.set_ylim(-.6,7.6);a.set_xlim(0,1);axes_style(a,reference=.5);title(a,letter,text)
    axs[-1].set_xlabel('AUROC');fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.60,.995),ncol=5,frameon=False)
    finish(fig,'S7','Five_Seeds_All_Metrics','Results 3.6','Table S8A',[ST+'module_E_five_seed.csv',ST+'module_E_five_seed_summary.csv'],data,'Full seed sensitivity for directional AUROC and summary_min. Each panel shows five recorded seed-specific estimates, stored minimum-to-maximum ranges and medians. Range segments are not confidence intervals. Effective populations follow each seed\'s valid scores; seed 17 loses one A-only AChE/BChE observation relative to seed 42.','五种子下两个方向 AUROC 和 summary_min 的完整变化。短线为种子间范围，非置信区间。',si=True)

def si8():
    source=CAP+'CAP5_PRIMARY_VS_ORIGINAL.csv';data=[r for r in read(source)if r['layer']=='PRIMARY'and r['pair']in ['JAK1/JAK2','PPARA/PPARD']];by={(r['pair'],r['method'],r['contrast']):r for r in data}
    fig,ax=plt.subplots(1,2,figsize=(7,6.6),sharex=True);fig.subplots_adjust(left=.17,right=.97,bottom=.17,top=.88,wspace=.38);keys=[(m,c)for m in METHODS for c in CONTRAST];y=np.arange(15)[::-1]
    for a,p,letter in zip(ax,['JAK1/JAK2','PPARA/PPARD'],'ab'):
        rr=[by[p,m,c]for m,c in keys]
        for prefix,color,m,off,label,open_ in [('original',GREY,'o',.12,'Primary',True),('cap5',PURPLE,'D',-.12,'CAP5',False)]:interval(a,[f(r,prefix+'_value')for r in rr],[f(r,prefix+'_ci_lo')for r in rr],[f(r,prefix+'_ci_hi')for r in rr],y+off,color,m,label,open_)
        a.set_yticks(y,[LABEL[m]+' '+{'dual_vs_A_only':'AUROC_B','dual_vs_B_only':'AUROC_A','summary_min':'min'}[c]+('*'if p=='PPARA/PPARD'and m=='M3'and c=='dual_vs_A_only'else'')for m,c in keys],fontsize=7.5);a.set_xlim(0,1);a.set_ylim(-.6,14.6);a.set_xlabel('AUROC');axes_style(a,reference=.5);title(a,letter,p)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.58,.07),ncol=2,frameon=False)
    fig.text(.17,.035,'* PPARA/PPARD S5 AUROC_B: lower 95% CI 0.4922 to 0.5030',fontsize=7.5)
    finish(fig,'S8','CAP5_All_Scoring_Workflows','Results 3.6, CAP5 paragraph','Tables S5A, S9; Supplementary Data',[source],data,'Full directional CAP5 sensitivity in the two pairs with excluded members. Each scheme has AUROC_B (dual versus A-only), AUROC_A (dual versus B-only) and summary_min. Primary and CAP5 estimates retain their respective stored 95% intervals; no interval for the difference is added. The asterisk identifies PPARA/PPARD S5 AUROC_B, whose interval changes from [0.4921875, 0.767578125] to [0.5030241935, 0.7772177419]. This boundary change is reported explicitly without a multiplicity-adjusted significance claim. The other six pairs have unchanged members.','CAP5 实际改变成员的两个靶点对中，五种评分的两个方向及 summary_min 的完整比较；标明 PPARA/PPARD S5 的区间边界变化。',si=True)

def si9():
    source=ST+'module_I_descriptive.csv';allrows=read(source);by={(r['pair'],r['method'],r['population_kind']):r for r in allrows};fig,axs=plt.subplots(2,2,figsize=(7,7.4),sharex=True,sharey=True);fig.subplots_adjust(left=.19,right=.97,bottom=.14,top=.92,wspace=.15,hspace=.25);records=[]
    for a,m,letter in zip(axs.ravel(),METHODS[1:],'abcd'):
        for i,p in enumerate(ORDER):
            for rr,off,hatch in [(by[p,'M0','pairwise_common_vs_'+m],.18,'///'),(by[p,m,'pairwise_common_vs_M0'],-.18,'')]:
                bottom=0;records.append(rr)
                for field,c in [('n_dual_top',BLUE),('n_A_only_top',ORANGE),('n_B_only_top',GREEN)]:
                    v=int(rr[field]);a.barh(7-i+off,v,left=bottom,height=.32,color=c,hatch=hatch,edgecolor='white',linewidth=.5)
                    if v:a.text(bottom+v/2,7-i+off,str(v),ha='center',va='center',color='black',fontsize=7.5,bbox={'facecolor':'white','edgecolor':'none','pad':.15})
                    bottom+=v
        a.set_yticks(range(7,-1,-1),ORDER,fontsize=7.5);a.set_ylim(-.6,7.6);a.set_xlim(0,10.5);a.set_xticks([0,5,10]);axes_style(a);title(a,letter,LABEL[m]+' vs S1'+(' (supporting)'if m=='M1b'else''))
    for a in axs[1]:a.set_xlabel('Ligands selected in the top 10% (n)')
    handles=[Patch(facecolor=c,label=l)for c,l in [(BLUE,'dual'),(ORANGE,'A-only'),(GREEN,'B-only')]]+[Patch(facecolor='white',edgecolor='#444444',hatch='///',label='S1 (upper bar)')];fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.57,.015),ncol=2,frameon=False)
    finish(fig,'S9','Matched_Candidate_Composition','Results 3.7','Tables 2, S10',[source],records,'Class composition of matched top-10% joint-ranking lists. Each target pair has an upper hatched S1 bar and a lower solid new-scheme bar, evaluated on that comparison\'s common two-target complete-case members. Segment labels give counts. Population sizes and k may differ between method comparisons; the same comparison has identical N and k for its two bars. S4-S1 is supporting. No method-specific S1 list is substituted for these paired baselines.','各方法与 Vina 在各自共同成员上的前 10% 候选名单。上方斜线为 S1，下方实心为补充评分；每段数字为类别人数。',si=True)

def scheme():
    fig,ax=plt.subplots(figsize=(7,3.6));fig.subplots_adjust(left=.03,right=.97,bottom=.06,top=.94);ax.set_xlim(0,10);ax.set_ylim(0,5);ax.axis('off')
    def box(x,y,w,h,text,fill):
        p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.03',facecolor=fill,edgecolor='#455365',lw=.7);ax.add_patch(p);ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=8)
    def arrow(x,y,xx,yy):ax.annotate('',(xx,yy),xytext=(x,y),arrowprops={'arrowstyle':'->','lw':.8,'color':'#455365'})
    box(.15,3.4,2,1.1,'Vina search\nseed 42','#E7EEF5');box(2.8,3.4,2,1.1,'Vina mode 1\nS1: -affinity','#E7EEF5');arrow(2.2,3.95,2.75,3.95)
    box(5.6,3.5,3.8,.9,'S2: mode-1 CNNscore','#FFF0E7');arrow(4.85,4.08,5.55,4.0)
    box(5.6,2.25,3.8,.9,'S3: mode-1 RTMScore','#E5F2EB');arrow(4.85,3.6,5.55,2.7)
    box(2.8,1.8,2,1,'All saved\nVina poses','#F0E8F7');arrow(1.2,3.35,2.75,2.3);box(5.6,1.0,3.8,.9,'S4: maximum saved CNNscore','#F0E8F7');arrow(4.85,2.2,5.55,1.45)
    box(.15,.1,4.65,.95,'Same prepared receptor and ligand inputs\nIndependent GNINA search\nS5: top CNNscore-ranked pose','#E5F2FA')
    ax.text(.15,4.85,'Scoring and pose-selection workflows',fontsize=10,fontweight='bold')
    finish(fig,'Scheme_S1','Scoring_Workflows','Methods 2.4','Scoring definitions',['reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml'],[], 'Scoring workflows S1-S5. S2 and S3 rescore the same Vina mode-1 pose used by S1. S4 chooses the highest CNNscore among all poses saved from that Vina search. S5 uses an independent GNINA search and CNNscore pose ordering. S4-S2 tests reranking, whereas S5-S1 compares complete docking workflows.','S1–S5 的评分和姿态选择关系。S4−S2 为保存姿态重新排序比较，S5−S1 为独立工作流比较。',si=True)

def main():
    for p in [MAIN,SI,DATA]:p.mkdir(parents=True,exist_ok=True)
    family=style()
    fig1();fig2();fig3();control(4);fig5();control(6,True);fig7();fig8();fig9()
    si1();si2();si3();descriptor_figure('S4','D_vs_A','A_only');descriptor_figure('S5','D_vs_B','B_only');si6();si7();si8();si9();scheme()
    with(OUT/'FIGURE_INDEX.csv').open('w',encoding='utf-8',newline='')as h:
        w=csv.writer(h,lineterminator='\n');w.writerow(['figure','group','first_citation','file_stem','tables','source_paths'])
        for r in FIGURES:w.writerow([r['number'],r['group'],r['first_citation'],('06_SI_Figures/'if r['group']=='SI'else'05_Main_Figures/')+r['stem'],r['table'],';'.join(r['sources'])])
    captions=[]
    for r in FIGURES:
        name=('Scheme S1'if r['number']=='Scheme_S1'else'Figure '+r['number']);captions += ['## '+name,'',r['caption'],'',r['caption_cn'],'']
    (OUT/'FIGURE_CAPTIONS.md').write_text('\n'.join(captions),encoding='utf-8')
    renumber=[['Fig01_Primary_Directional_and_Deltas','Figure 2; Figure 3','Split directional estimates from method differences'],['Fig02_Negative_Class_Contrast','Figure 4','Redrawn'],['Fig03_Chemistry_and_Docking_Increment','Figure 5','Direct Vina score is no longer called OOF'],['Fig04_Pocket_Contrast','Figure 6','Corresponding/noncorresponding terminology'],['Fig05_Seed_and_Receptor_Sensitivity','Figure 7','Readable labels; seed range distinct from CI'],['Fig06_Directional_and_Topk_Changes','Figure 9','S labels; categories and pair key added'],['FigS01_Descriptors_D_vs_A','Figure S4','Redrawn'],['FigS02_Descriptors_D_vs_B','Figure S5','Redrawn'],['FigS03_Cross_Class_Similarity','Figure S6','Redrawn'],['FigS04_All_Method_AUROC','Figure S3','S1-S5 order'],['FigS05_Five_Seeds','Figure S7','Redrawn'],['FigS06_Alternative_Receptors','Figure 7b; Table S8B','Unreplaced zero retained in table; redundant figure removed'],['FigS07_Missingness','Figure S2','Redrawn'],['FigS08_Topk_and_Concordance','Figure S9; Table S10C','Matched ranking populations replace method-specific composition view'],['Scheme01_Study_Design','Scheme S1; Methods text','Workflow redraw; no extra main-text figure'],['NEW','Figure 1','Evaluation-set composition and source audit'],['NEW','Figure 8','CAP5 and strict-subset sensitivity'],['NEW','Figure S1','Recorded redocking outcomes'],['NEW','Figure S8','Full CAP5 directional comparisons in affected pairs']]
    with(OUT/'FIGURE_RENUMBERING.csv').open('w',encoding='utf-8',newline='')as h:w=csv.writer(h,lineterminator='\n');w.writerow(['previous_graphic','submission_graphic','change']);w.writerows(renumber)
    from pypdf import PdfWriter
    writer=PdfWriter()
    for r in FIGURES:
        p=(SI if r['group']=='SI'else MAIN)/(r['stem']+'.pdf');writer.append(str(p),outline_item=('Scheme S1'if r['number']=='Scheme_S1'else'Figure '+r['number']))
    with(OUT/'Figure_Review.pdf').open('wb')as h:writer.write(h)
    report={'source_commit':'a7b9d7521dff56ece18ec10a55734392f417e684','font':family,'main_figures':9,'supplementary_figures':9,'supplementary_schemes':1,'scientific_recomputation':False,'formats':['pdf','svg','png'],'dpi':600,'figures':QUALITY}
    (OUT/'04_Editorial_Checks/FIGURE_TECHNICAL_QA.json').write_text(json.dumps(report,indent=2)+'\n')
    outside=[r['figure']for r in QUALITY if r['outside_canvas_text']]
    print(json.dumps({'main':9,'SI':9,'scheme':1,'font':family,'outside_canvas_figures':outside,'min_font_pt':min(r['smallest_visible_font_pt']for r in QUALITY)}))

if __name__=='__main__':main()
