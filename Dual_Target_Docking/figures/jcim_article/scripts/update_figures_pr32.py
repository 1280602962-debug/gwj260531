"""Rebuild local JCIM artwork from the pinned PR32 inputs, preserving v3 style.

First run: --source-root PATH_TO_PR32/Dual_Target_Docking
Subsequent runs use the small input_snapshot stored alongside the figures.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
from PIL import Image
import jcim_figure_style as style
import plot_jcim_article_figures_v3 as v

OUT = Path(__file__).resolve().parents[1]
SNAP = OUT / 'input_snapshot'
SHA = 'abb61a20a04eb6a085ad526876624eadb518c4cc'
C = style.C
PAIRS = style.PRIMARY_PAIRS
SMIN = r'summary$_{\mathrm{min}}$'
GENERATED = []
READS = {}
P = {}
SOURCE = None

def snapshot_path(rel):
    # Short flat names avoid Windows MAX_PATH for the nested census table paths.
    return SNAP / (hashlib.sha256(rel.encode()).hexdigest()[:16]+Path(rel).suffix)

def source_path(rel):
    return snapshot_path(rel) if SOURCE==SNAP.resolve() else SOURCE/rel

def read(rel):
    path = source_path(rel)
    READS[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def save(fig, stem, toc=False):
    # Preserve established 7-inch, 300-dpi RGB triplets and embedded PDF fonts.
    if stem=='Fig3_ligand_chemistry':
        # Lift 3A title+legend above the EGFR/HER2 row; park 3B legend inside the axes.
        fig.set_size_inches(7,6.9)
        fig.subplots_adjust(hspace=.50,bottom=.09,top=.82,left=.16,right=.98)
        old=fig.axes[0].get_legend()
        handles=old.legend_handles; labels=[t.get_text() for t in old.get_texts()]
        old.remove()
        fig.axes[0].legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,1.04),ncol=4,fontsize=6,
                           columnspacing=.9,handletextpad=.35,borderaxespad=0)
        old=fig.axes[1].get_legend()
        handles=old.legend_handles; labels=[t.get_text() for t in old.get_texts()]
        old.remove()
        fig.axes[1].legend(handles,labels,loc='lower right',fontsize=6,ncol=1,frameon=True,
                           fancybox=False,edgecolor='none',facecolor='white',framealpha=.92,
                           borderpad=.35,handletextpad=.35)
    if stem=='Fig5_mismatched_pocket':
        for ax in fig.axes[:2]:ax.set_xlabel('Δ'+SMIN+'\n(matched − mismatched)',fontsize=7)
    GENERATED.append(stem)
    style.save_all(fig, stem, toc=toc)
    plt.close(fig)

def label(ax, letter):
    style.panel_label(ax, letter, x=-0.12, y=1.04)

def dotci(ax, x, lo, hi, y, color, marker='o'):
    ax.plot([lo, hi], [y, y], color=color, lw=1.1)
    ax.plot(x, y, marker=marker, color=color, ms=4.1, ls='none')

def fig1(D):
    census = next(r for r in read('data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv') if r['slice']=='all')
    fig = plt.figure(figsize=(7,4.9))
    gs = fig.add_gridspec(2,2,height_ratios=[1,1.25],wspace=.36,hspace=.40)
    ax = fig.add_subplot(gs[0,0]); label(ax,'A')
    ax.set(xlim=(0,2),ylim=(0,2)); ax.set_aspect('equal')
    for x,y,name,col in [(0,1,'dual',C['dual']),(1,1,'A-only',C['a_only']),(0,0,'B-only',C['b_only']),(1,0,'neither',C['neither'])]:
        ax.add_patch(plt.Rectangle((x,y),1,1,fc=col,alpha=.23,ec='white',lw=2))
        ax.text(x+.5,y+.5,name,ha='center',va='center',fontsize=8)
    ax.set_xticks([.5,1.5],['B active','B low activity'],fontsize=7)
    ax.set_yticks([.5,1.5],['A low activity','A active'],fontsize=7)
    ax.tick_params(length=0); ax.set_title('Experimental activity states',fontsize=8)
    for s in ax.spines.values(): s.set_visible(False)
    ax = fig.add_subplot(gs[0,1]);label(ax,'B');ax.axis('off')
    for y,txt,col in [(.69,'Dual vs A-only  →  target B score',C['vina']),(.36,'Dual vs B-only  →  target A score',C['a_only'])]:
        ax.add_patch(FancyBboxPatch((.02,y-.11),.96,.22,boxstyle='round,pad=.02',fc='white',ec=col,lw=1))
        ax.text(.5,y,txt,ha='center',va='center',fontsize=7.2)
    ax.text(.5,.05,'Two directional AUROCs; lower value = '+SMIN,ha='center',fontsize=6.8)
    ax.set_title('Directional evaluation',fontsize=8)
    ax = fig.add_subplot(gs[1,:]);style.panel_label(ax,'C',x=0,y=1.04);ax.axis('off')
    ax.set_title('ChEMBL target-pair census',fontsize=8,pad=7)
    keys=['n_pairs_n_both_ge_1','n_pairs_n_both_ge_10','n_directional_n10','n_strict_thick']
    names=['≥1 ligand measured at both targets','≥10 ligands measured at both targets','Dual, A-only and B-only each ≥10 (θ=6.0)','Strict bidirectional supply criterion (6.5/5.5)']
    for i,(k,txt) in enumerate(zip(keys,names)):
        yy=.91-i*.205
        ax.add_patch(FancyBboxPatch((.03,yy-.073),.94,.145,boxstyle='round,pad=.007',fc='#F4F7FA',ec='#D5DDE4',lw=.7))
        ax.text(.075,yy,format(int(census[k]),','),va='center',fontsize=9,fontweight='bold',color=C['vina'])
        ax.text(.32,yy,txt,va='center',fontsize=7.5)
    ax.text(.5,.04,'8 primary target pairs after panel construction and structural assessment',ha='center',fontsize=7.4)
    P['fig1C']=census
    fig.subplots_adjust(left=.16,right=.97,top=.92,bottom=.055)
    save(fig,'Fig1_four_state_and_supply')

def fig1_c_data(D):
    """Standalone Python-rendered data panel for submission as a separate figure."""
    census = next(r for r in read('data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv') if r['slice']=='all')
    fig, ax = plt.subplots(figsize=(7, 3.25))
    style.panel_label(ax, 'C', x=0.0, y=1.08)
    ax.axis('off')
    ax.set_title('ChEMBL target-pair census', fontsize=9, pad=8)
    keys=['n_pairs_n_both_ge_1','n_pairs_n_both_ge_10','n_directional_n10','n_strict_thick']
    names=['≥1 ligand measured at both targets','≥10 ligands measured at both targets','Dual, A-only and B-only each ≥10 (θ=6.0)','Strict bidirectional supply criterion (6.5/5.5)']
    for i,(k,txt) in enumerate(zip(keys,names)):
        yy=.88-i*.22
        ax.add_patch(FancyBboxPatch((.03,yy-.075),.94,.15,boxstyle='round,pad=.007',fc='#F4F7FA',ec='#D5DDE4',lw=.75))
        ax.text(.075,yy,format(int(census[k]),','),va='center',fontsize=10,fontweight='bold',color=C['vina'])
        ax.text(.32,yy,txt,va='center',fontsize=8.1)
    ax.text(.5,.035,'8 primary target pairs after panel construction and structural assessment',ha='center',fontsize=7.8)
    P['fig1C_standalone']=census
    fig.subplots_adjust(left=.02,right=.98,top=.86,bottom=.05)
    save(fig,'Fig1_C_chEMBL_supply')

def fig2(D):
    fig,axs=plt.subplots(3,1,figsize=(7,7.6),gridspec_kw={'height_ratios':[1,1.15,1.15]})
    for ax,letter in zip(axs,'ABC'):
        label(ax,letter);ax.set_yticks(range(8),PAIRS,fontsize=7);ax.set_ylim(7.65,-.65)
    for i,p in enumerate(PAIRS):
        r=v.primary_row(D,p)
        axs[0].plot([r['da'],r['db']],[i,i],color='#B8B8B8',lw=.8)
        axs[0].plot(r['da'],i,'o',color=C['vina'],ms=4.2)
        axs[0].plot(r['db'],i,'s',color=C['a_only'],ms=4)
        dotci(axs[1],r['smin'],r['lo'],r['hi'],i-.16,C['vina'])
        dotci(axs[1],r['nei'],r['nei_lo'],r['nei_hi'],i+.16,C['desc'],'D' if r['n_neg']<10 else 's')
        src=D['equal'] if p in v.ORIGINAL_THREE else D['five_s34']
        for j,(contrast,col,m) in enumerate([('D_vs_B_or_neither_pocketA',C['vina'],'o'),('D_vs_A_or_neither_pocketB',C['a_only'],'s')]):
            rr=src[(p,contrast)]
            dotci(axs[2],float(rr['delta_neither_minus_selective']),float(rr['delta_ci_lo']),float(rr['delta_ci_hi']),i+(-.16 if j==0 else .16),col,'D' if rr['underpowered_neither']=='1' else m)
    for ax in axs[:2]:ax.axvline(.5,color=C['chance'],ls='--',lw=.8);ax.set(xlim=(0,1.02),xlabel='AUROC')
    axs[2].axvline(0,color=C['chance'],ls='--',lw=.8);axs[2].set(xlim=(-.58,.72),xlabel='ΔAUROC (Dual vs neither − Dual vs selective)')
    titles=['Pocket-matched directional AUROC','Directional summary and mean-score Dual vs neither','Negative-class contrast with the score channel held fixed']
    labs=[['Dual vs A-only, target B','Dual vs B-only, target A'],['Directional '+SMIN,'Dual vs neither, mean score'],['Target A score: B-only → neither','Target B score: A-only → neither']]
    for n,(ax,title,ll) in enumerate(zip(axs,titles,labs)):
        ax.set_title(title,fontsize=8,pad=5)
        handles=[Line2D([],[],color=C['vina'],marker='o',ls='none',ms=4,label=ll[0]),
                 Line2D([],[],color=C['desc'] if n==1 else C['a_only'],marker='s',ls='none',ms=4,label=ll[1])]
        if n==1:
            handles.append(Line2D([],[],color=C['desc'],marker='D',ls='none',ms=4,label='neither n=4'))
        if n==2:
            handles.append(Line2D([],[],color=C['desc'],marker='D',ls='none',ms=4,label='underpowered neither'))
        ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.24),ncol=len(handles),fontsize=6.3)
    P['fig2']={p:{'primary':v.primary_row(D,p),'fixed_A':D['equal' if p in v.ORIGINAL_THREE else 'five_s34'][(p,'D_vs_B_or_neither_pocketA')],'fixed_B':D['equal' if p in v.ORIGINAL_THREE else 'five_s34'][(p,'D_vs_A_or_neither_pocketB')]} for p in PAIRS}
    fig.subplots_adjust(left=.19,right=.97,top=.945,bottom=.085,hspace=.76)
    save(fig,'Fig2_negative_class_formulation')

def e8_value():
    panel={r['panel_id']:r for r in read('data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv')}
    scores=read('data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv')
    groups={k:[] for k in ['dual','A_only','B_only']}
    for r in scores:
        cls=panel[r['ligand']]['class']
        if cls in groups and r['4L23_affinity_E8'] and r['4JT6_affinity_E8']:
            groups[cls].append([-float(r['4L23_affinity_E8']),-float(r['4JT6_affinity_E8'])])
    def auc(a,b):
        dif=np.asarray(a)[:,None]-np.asarray(b)[None,:]
        return float(np.mean((dif>0)+.5*(dif==0)))
    return min(auc([x[1] for x in groups['dual']],[x[1] for x in groups['A_only']]),auc([x[0] for x in groups['dual']],[x[0] for x in groups['B_only']]))

def counts_heatmap(ax,rows,columns,title,gate):
    mat=np.array([[int(r[k]) for k in columns] for r in rows])
    # Threshold-relative color prevents the largest dual counts masking shortages.
    ax.imshow(np.minimum(mat/gate,1),cmap='Blues',vmin=0,vmax=1,aspect='auto')
    for i in range(len(rows)):
        for j in range(3):ax.text(j,i,str(mat[i,j]),ha='center',va='center',fontsize=7,color='white' if mat[i,j]>=gate*.65 else C['ink'])
    ax.set_yticks(range(len(rows)),[r['pair'] for r in rows],fontsize=6.5)
    ax.set_xticks(range(3),['Dual','A-only','B-only'],fontsize=7)
    ax.set_title(title,fontsize=8,pad=6)
    ax.tick_params(length=0)
    return mat.tolist()

def fig6(D):
    # Stack B/C/D in the right column so the two-point traces match panel A in scale.
    fig=plt.figure(figsize=(7,6.85))
    gs=fig.add_gridspec(3,2,height_ratios=[.90,.90,1.42],width_ratios=[1.18,1.00],wspace=.40,hspace=.48)
    ax=fig.add_subplot(gs[:,0]);label(ax,'A')
    rules=['theta_5.5','theta_6.0','theta_6.5','strict_6.5_5.5']
    recs=[[v.theta_grid_record(D,p,r) for r in rules] for p in PAIRS]
    mat=np.array([[r['value'] for r in row] for row in recs])
    im=ax.imshow(mat,cmap='RdBu',vmin=.25,vmax=.75,aspect='auto')
    for i in range(8):
        for j in range(4):ax.text(j,i,f'{mat[i,j]:.3f}'+('†' if recs[i][j]['under'] else ''),ha='center',va='center',fontsize=6,color='white' if mat[i,j]<.35 or mat[i,j]>.67 else C['ink'])
    ax.set_yticks(range(8),PAIRS,fontsize=6.7);ax.set_xticks(range(4),['θ=5.5','θ=6.0','θ=6.5','strict'],fontsize=6.5);ax.set_title('Activity definitions',fontsize=8)
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.035,pad=.065);cb.set_label(SMIN,fontsize=7)
    P['fig6A']=mat.tolist()
    for letter,row,title,vals,names in [('B',0,'Panel size',[float(D['pm110'][('PM48','vina')]['summary_min']),float(D['pm110'][('PM110','vina')]['summary_min'])],['PM48','PM110']),('C',1,'Exhaustiveness',[v.primary_row(D,'PIK3CA/mTOR')['smin'],e8_value()],['16','8'])]:
        ax=fig.add_subplot(gs[row,1]);label(ax,letter)
        ax.plot(range(2),vals,'-o',color=C['vina'],ms=5.5,lw=1.0);ax.axhline(.5,color=C['chance'],ls='--',lw=.7)
        ax.set_xticks(range(2),names,fontsize=7);ax.set(xlim=(-.35,1.35),ylim=(.45,.82));ax.set_title(title,fontsize=8)
        ax.set_ylabel(SMIN,fontsize=7)
        P['fig6'+letter]=dict(zip(names,vals))
    ax=fig.add_subplot(gs[2,1]);label(ax,'D')
    rows=[next(r for r in D['native'] if r['pair']==p) for p in PAIRS]
    P['fig6D']=counts_heatmap(ax,rows,['n_dual','n_A_only','n_B_only'],'BindingDB after independence filters',20)
    ax.text(.5,-.16,'Cell color saturates at n=20\n0/8 pairs meet all external criteria',transform=ax.transAxes,ha='center',fontsize=6.5)
    fig.subplots_adjust(left=.16,right=.98,top=.94,bottom=.10)
    save(fig,'Fig6_evidence_boundary')

def supplements(D):
    v.fig_s4_forest(D);v.fig_s5_holdout(D)
    sim=[r for r in read('data/jcim_novelty_v0/tables/detectable_effect_simulation_v1.csv') if r['contrast']=='summary_min']
    pp=[p for p in PAIRS if any(r['pair']==p for r in sim)]
    grid=sorted({float(r['true_auroc']) for r in sim})
    vals=[[float(next(r for r in sim if r['pair']==p and float(r['true_auroc'])==a)['p_ci_excludes_0p5']) for a in grid] for p in pp]
    fig,ax=plt.subplots(figsize=(7,2.9));im=ax.imshow(vals,cmap='YlGnBu',vmin=0,vmax=1,aspect='auto')
    for i,row in enumerate(vals):
        for j,z in enumerate(row):ax.text(j,i,f'{z:.3f}',ha='center',va='center',fontsize=7,color='white' if z>.55 else C['ink'])
    ax.set_yticks(range(len(pp)),pp);ax.set_xticks(range(len(grid)),[f'{g:.2f}' for g in grid]);ax.set_xlabel('True AUROC on both directional arms')
    fig.colorbar(im,ax=ax).set_label('P(95% CI excludes 0.5)')
    ax.set_title('Detectable-effect simulation: available three-pair results',fontsize=8)
    fig.subplots_adjust(left=.18,right=.95,bottom=.23,top=.85);P['figS6']=sim;save(fig,'FigS6_detectable_effect')
    rows=[next(r for r in D['native'] if r['pair']==p) for p in PAIRS]
    fig,axs=plt.subplots(1,2,figsize=(7,4.2));
    for ax,cols,title,gate,letter in [(axs[0],['n_dual','n_A_only','n_B_only'],'Compounds (criterion: n≥20 per class)',20,'A'),(axs[1],['n_sources_dual','n_sources_A_only','n_sources_B_only'],'Sources (criterion: ≥3 per class)',3,'B')]:
        label(ax,letter);counts_heatmap(ax,rows,cols,title,gate)
    fig.subplots_adjust(left=.16,right=.98,wspace=.55,bottom=.14,top=.87);P['figS8']=rows;save(fig,'FigS8_bindingdb_native_slice')
    clusters=read('data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv')
    fig,axs=plt.subplots(1,2,figsize=(7,3),sharex=True)
    for ax,p,letter in zip(axs,['EGFR/HER2','JAK1/TYK2'],'AB'):
        label(ax,letter);r=v.s34_row(D,p);dotci(ax,r['delta'],r['lo'],r['hi'],0,C['vina'])
        for i,key in enumerate(['scaffold_cluster','document_cluster'],1):
            r=next(r for r in clusters if r['pair']==p and r['estimator']==key)
            dotci(ax,float(r['delta_point']),float(r['delta_ci_lo']),float(r['delta_ci_hi']),i,C['desc'] if i==1 else C['a_only'])
        ax.axvline(0,color=C['chance'],ls='--',lw=.7);ax.set(xlim=(-.12,.75),ylim=(2.5,-.5),xlabel='ΔAUROC, target A score')
        ax.set_yticks(range(3),['Ligand','Scaffold cluster','Document cluster'],fontsize=6.5);ax.set_title(p,fontsize=8)
    fig.subplots_adjust(left=.18,right=.97,wspace=.68,top=.83,bottom=.23);P['figS11']=clusters;save(fig,'FigS11_cluster_uncertainty')
    top=next(r for r in read('data/jcim_novelty_v0/tables/mixed_library_enrichment_v1.csv') if r['pair']=='EGFR/HER2' and r['score']=='vina_mean' and r['cutoff']=='Top10')
    filt=next(r for r in read('data/jcim_novelty_v0/tables/and_filter_operating_point_v1.csv') if r['pair']=='EGFR/HER2' and r['score']=='vina_worst' and r['dual_percentile']=='50')
    fig,axs=plt.subplots(1,2,figsize=(7,2.9))
    for ax,vals,title,letter in [(axs[0],[int(top[k]) for k in ['n_dual_top','n_A_only_top','n_B_only_top','n_neither_top']],'Top-10 by mean Vina score (library n=110)','A'),(axs[1],[int(filt[k]) for k in ['n_dual_pass','n_A_only_pass','n_B_only_pass']],'Dual-median AND filter (library n=98)','B')]:
        label(ax,letter); names=['Dual','A-only','B-only','Neither'][:len(vals)]
        ax.bar(range(len(vals)),vals,color=[C['dual'],C['a_only'],C['b_only'],C['neither']][:len(vals)],width=.6)
        for i,z in enumerate(vals):ax.text(i,z+.35,str(z),ha='center',fontsize=7)
        ax.set_xticks(range(len(vals)),names,fontsize=7);ax.set_ylabel('Retained compounds');ax.set_ylim(0,max(vals)*1.22+1);ax.set_title(title,fontsize=7.5)
    fig.subplots_adjust(left=.1,right=.97,wspace=.35,top=.8,bottom=.2);P['figS7']={'top10':top,'and_filter':filt};save(fig,'FigS7_posthoc_diagnostics')

def main():
    global SOURCE
    parser=argparse.ArgumentParser();parser.add_argument('--source-root',type=Path,default=SNAP);args=parser.parse_args();SOURCE=args.source_root.resolve()
    style.OUT=OUT;style.apply_style();v.ROOT=SOURCE;v.DATA=SOURCE/'data';v.OUT=OUT;v.save_all=save
    v._read=lambda path:read(str(path.relative_to(SOURCE)).replace('\\','/'))
    D=v.load()
    # Canonical Table-2/3 inputs are loaded before any artwork is overwritten.
    assert set(D['theta6'])==set(v.ORIGINAL_THREE)
    assert set(D['native'][i]['pair'] for i in range(len(D['native'])))==set(PAIRS)
    assert all(float(r['packaged_as_external_evaluation'])==0 for r in D['native'])
    fig1(D);fig1_c_data(D);fig2(D);v.fig3_chemistry(D);v.fig4_realization(D);v.fig5_mismatched(D);fig6(D);supplements(D);v.toc_graphic()
    P.update(v.PROVENANCE['plotted'])
    P['primary']={p:v.primary_row(D,p) for p in PAIRS}
    assert abs(P['fig3B_max_abs']-.023)<.001
    # Keep exactly the table inputs read by this generator, pinned to the PR head.
    for rel in READS:
        dest=snapshot_path(rel);dest.parent.mkdir(parents=True,exist_ok=True)
        if SOURCE!=SNAP.resolve():shutil.copyfile(SOURCE/rel,dest)
    for rel in ['docs/MANUSCRIPT_JCIM_ZH.md','docs/MANUSCRIPT_JCIM_EN.md','docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md']:
        path=source_path(rel);READS[rel]=hashlib.sha256(path.read_bytes()).hexdigest();dest=snapshot_path(rel);dest.parent.mkdir(parents=True,exist_ok=True)
        if SOURCE!=SNAP.resolve():shutil.copyfile(path,dest)
    audit={'commit':SHA,'source':'https://github.com/1280602962-debug/gwj260531/pull/32','inputs_sha256':READS,'input_files':{rel:str(snapshot_path(rel).relative_to(OUT)) for rel in READS},'generated':GENERATED,'plotted':P}
    (OUT/'plotted_values.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    for stem in GENERATED:
        for ext in ['png','tif']:
            with Image.open(OUT/(stem+'.'+ext)) as im:
                assert im.mode=='RGB';assert abs(im.info['dpi'][0]-300)<1
                assert im.width<=2101 and im.height<=2751
    print(f'PASS: {len(GENERATED)} figures; {len(READS)} pinned inputs; RGB / 300 dpi / size checks')

if __name__=='__main__':main()
