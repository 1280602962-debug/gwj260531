#!/usr/bin/env python3
"""Static manuscript figures sourced exclusively from manuscript tables."""
from pathlib import Path
import hashlib
import json
import os
import tempfile
import pandas as pd
_cache=tempfile.TemporaryDirectory(prefix='manuscript-mpl-')
os.environ.setdefault('MPLCONFIGDIR',_cache.name)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'data/manuscript/figures'


def sha(path):
    b=path.read_bytes()
    if path.suffix in {'.csv','.py','.json','.svg'}:b=b.replace(b'\r\n',b'\n')
    return hashlib.sha256(b).hexdigest()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'svg.hashsalt':'urat1-nlrp3-manuscript'})
    inputs=[ROOT/'data/manuscript/tables/screening_funnel.csv',ROOT/'data/manuscript/tables/acid_descriptors.csv',
            ROOT/'data/manuscript/tables/candidate_seed_metrics.csv',Path(__file__)]
    funnel=pd.read_csv(inputs[0]).set_index('stage')
    n=lambda key:int(funnel.loc[key,'output_n'])
    fig,ax=plt.subplots(figsize=(11,7.6));ax.set_xlim(0,11);ax.set_ylim(0,8);ax.axis('off')
    def box(x,y,w,label,color='#e9eff5'):
        ax.add_patch(FancyBboxPatch((x,y),w,.66,boxstyle='round,pad=0.05',linewidth=.7,
                                   facecolor=color,edgecolor='#597084'))
        ax.text(x+w/2,y+.33,label,ha='center',va='center',fontsize=10)
    def arrow(x,y,xx,yy):
        ax.annotate('',xy=(xx,yy),xytext=(x,y),arrowprops={'arrowstyle':'->','color':'#597084'})
    x,w=3.35,4.3
    stages=[('clinical_library','Clinical library'),('nlrp3_activity_pool','NLRP3-related activity pool'),
            ('acid_equivalents','Acid-equivalent pool'),('chemistry_soft','Chemistry-soft docking pool')]
    ys=[6.85,5.65,4.45,3.25]
    for i,((key,label),y) in enumerate(zip(stages,ys)):
        box(x,y,w,f'{label}: {n(key):,}')
        if i<3:arrow(5.5,y,5.5,ys[i+1]+.66)
    box(.65,1.85,4.2,f'Dual structural >=2/3 seeds: {n("dual_structural_ge2of3")}')
    # Forty is the historical chemistry set, not the output of a sequential 54->40 gate.
    ledger=pd.read_csv(ROOT/'data/manuscript/tables/screening_decision_ledger.csv')
    inputs.append(ROOT/'data/manuscript/tables/screening_decision_ledger.csv')
    forty=int(ledger.stored_chemistry_eligible.sum())
    box(6.15,1.85,4.2,f'Historical chemistry eligible: {forty}')
    arrow(4.8,3.25,2.8,2.51);arrow(6.2,3.25,8.2,2.51)
    ax.text(5.5,1.85,'intersection',ha='center',va='center',fontsize=9)
    arrow(2.8,1.85,4.65,1.42);arrow(8.2,1.85,6.35,1.42)
    ax.text(5.5,1.3,f'{n("intersection_with_historical_chemistry_40")} before final exclusions',ha='center')
    box(1.25,.28,3.7,f'Primary: {n("primary")} (13 - 1 control)','#dbece7')
    box(6.05,.28,3.7,f'Reserve: {n("reserve_parallel_branch")} (24 - 3 beta-lactams)','#f4e8d9')
    arrow(5.1,1.22,3.1,.94);arrow(5.9,1.22,7.9,.94)
    fig.suptitle('Clinical-library screening and frozen nomination',fontsize=16,x=.5,y=.98)
    fig.text(.5,.015,'Geometry-supported hypotheses; historical chemistry exclusions disclosed; reserve is not a potency rank.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.95))
    files=[]
    for ext in ['png','svg']:
        p=OUT/f'figure_1_screening_funnel.{ext}';fig.savefig(p,dpi=220,metadata={'Date':None} if ext=='svg' else {});files.append(p)
    plt.close(fig)
    d=pd.read_csv(inputs[1]);c=pd.read_csv(inputs[2]);primary=set(c.loc[c.nomination_status.eq('primary'),'ligand_id'])
    soft=d[d.chemistry_soft_pass];pr=d[d.ligand_id.isin(primary)]
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for ax,y,label in [(axes[0],'tpsa',r'TPSA ($\AA^2$)'),(axes[1],'clogp','Calculated logP')]:
        ax.scatter(d.mw_recomputed,d[y],c='#cbd3dc',s=18,label='Acid pool (303)',alpha=.75)
        ax.scatter(soft.mw_recomputed,soft[y],c='#3e7599',s=24,label='Docking pool (156)',alpha=.8)
        ax.scatter(pr.mw_recomputed,pr[y],facecolor='#d27432',edgecolor='white',linewidth=.5,s=65,label='Primary (12)',zorder=4)
        pf=pr[pr.ligand_id.eq('REP_07580')].iloc[0]
        ax.annotate('PF-03882845',(pf.mw_recomputed,pf[y]),xytext=(.36,.80),textcoords='axes fraction',
                    arrowprops={'arrowstyle':'-','color':'#65452f'},fontsize=8,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':3})
        ax.set(xlabel='Molecular weight (Da)',ylabel=label)
        ax.grid(alpha=.16);ax.set_axisbelow(True)
    axes[0].legend(frameon=False,fontsize=8,loc='upper left')
    fig.suptitle('Chemical space of deposited acid-equivalent SMILES',fontsize=15)
    fig.text(.5,.01,'Descriptors recomputed for annotation. Groups are nested; logP was not a hard selection threshold.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.94))
    for ext in ['png','svg']:
        p=OUT/f'figure_2_chemical_space.{ext}';fig.savefig(p,dpi=220,metadata={'Date':None} if ext=='svg' else {});files.append(p)
    plt.close(fig)
    (OUT/'provenance.json').write_text(json.dumps({'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in inputs+files]},indent=2)+'\n',encoding='utf-8')
    print(f'Generated {len(files)} figure files from manuscript tables.')


if __name__=='__main__':main()
