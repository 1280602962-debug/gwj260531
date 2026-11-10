"""Independent numeric verification of the frozen values used in artwork.

No docking, learned model, bootstrap or study analysis is rerun. Exact AUROC
and top-list checks validate already accepted estimates from deposited scores.
"""
from __future__ import annotations
import csv
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'submission/Molecular_Diversity/Checks'
METHODS=['M0','M1','M2','M1b','M3']
checks=[]

def read(p):
    with (ROOT/p).open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))

def finite(r,m,t):
    try:return r[f'{m}_{t}_valid']=='1' and math.isfinite(float(r[f'{m}_score_{t}']))
    except (KeyError,ValueError):return False

def auc(pos,neg):
    x=np.asarray(pos)[:,None];y=np.asarray(neg)[None,:]
    return float(((x>y)+.5*(x==y)).mean())

def verify(label,observed,expected,tol=1e-9):
    ok=abs(float(observed)-float(expected))<=tol
    checks.append({'check':label,'observed':observed,'expected':expected,'tolerance':tol,'pass':ok})
    if not ok:raise ValueError(label+': '+str(observed)+' != '+str(expected))

def ranks_desc(values):
    a=np.asarray(values);order=np.argsort(-a,kind='stable');rank=np.zeros(len(a))
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and a[order[j]]==a[order[i]]:j+=1
        rank[order[i:j]]=(i+1+j)/2;i=j
    return rank

def main():
    pop=[r for r in read('results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv') if r['activity_eligible']=='1' and r['class']]
    verify('eligible observations',len(pop),867,0)
    verify('distinct ligand entities',len({r['global_ligand_entity_id'] for r in pop}),845,0)
    for cls,n in [('dual',248),('A_only',256),('B_only',252),('neither',111)]:verify('class '+cls,sum(r['class']==cls for r in pop),n,0)
    metric=read('results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv')
    for r in metric:
        pair,m=r['pair'],r['method'];rows=[x for x in pop if x['pair']==pair]
        dual=[x for x in rows if x['class']=='dual' and all(finite(x,m,t) for t in 'AB')]
        values={}
        for contrast,cls,t in [('dual_vs_A_only','A_only','B'),('dual_vs_B_only','B_only','A')]:
            neg=[x for x in rows if x['class']==cls and finite(x,m,t)]
            values[contrast]=auc([float(x[f'{m}_score_{t}']) for x in dual],[float(x[f'{m}_score_{t}']) for x in neg])
        observed=min(values.values()) if r['contrast']=='summary_min' else values[r['contrast']]
        verify(pair+' '+m+' '+r['contrast'],observed,r['auroc'])
        verify(pair+' '+m+' dual n',len(dual),r['n_dual'],0)
        if not 0<=float(r['ci_lo'])<=float(r['ci_hi'])<=1:raise ValueError('Invalid interval '+str(r))
        verify('bootstrap B',r['boot_B'],10000,0);verify('bootstrap seed',r['boot_seed'],271828,0)
    for r in read('results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv'):
        x,y=r['delta'].split('-');rows=[a for a in pop if a['pair']==r['pair']]
        dual=[a for a in rows if a['class']=='dual' and all(finite(a,m,t) for m in [x,y] for t in 'AB')]
        mins=[]
        for m in [x,y]:
            vals=[]
            for cls,t in [('A_only','B'),('B_only','A')]:
                neg=[a for a in rows if a['class']==cls and all(finite(a,mm,t) for mm in [x,y])]
                vals.append(auc([float(a[f'{m}_score_{t}']) for a in dual],[float(a[f'{m}_score_{t}']) for a in neg]))
            mins.append(min(vals))
        verify(r['pair']+' paired '+r['delta'],mins[0]-mins[1],r['delta_summary_min'])
        verify(r['pair']+' shared dual '+r['delta'],len(dual),r['n_dual_four_sided'],0)
    for r in read('results/jcim_stage2/module_B_negative_class.csv'):
        rows=[a for a in pop if a['pair']==r['pair']];t='B' if r['direction']=='D_vs_A' else 'A';cls='A_only' if t=='B' else 'B_only'
        dual=[a for a in rows if a['class']=='dual' and all(finite(a,'M0',z) for z in 'AB')]
        vals=[]
        for c in [cls,'neither']:
            neg=[a for a in rows if a['class']==c and finite(a,'M0',t)]
            vals.append(auc([float(a[f'M0_score_{t}']) for a in dual],[float(a[f'M0_score_{t}']) for a in neg]))
        for value,key in zip([*vals,vals[1]-vals[0]],['auc_single','auc_neither','delta_negative']):verify(r['pair']+' '+r['direction']+' '+key,value,r[key])
    for r in read('results/jcim_stage2/module_C_wrong_pocket.csv'):
        rows=[a for a in pop if a['pair']==r['pair']];t='B' if r['arm']=='D_vs_A' else 'A';other='A' if t=='B' else 'B';cls='A_only' if t=='B' else 'B_only'
        dual=[a for a in rows if a['class']=='dual' and all(finite(a,'M0',z) for z in 'AB')]
        neg=[a for a in rows if a['class']==cls and all(finite(a,'M0',z) for z in 'AB')]
        vals=[auc([float(a[f'M0_score_{z}']) for a in dual],[float(a[f'M0_score_{z}']) for a in neg]) for z in [t,other]]
        for value,key in zip([*vals,vals[0]-vals[1]],['auc_correct','auc_wrong','delta_pocket']):verify(r['pair']+' '+r['arm']+' '+key,value,r[key])
    descriptive=read('results/jcim_stage2/module_I_descriptive.csv')
    for r in descriptive:
        m=r['method'];others=[m]
        if r['population_kind'].startswith('pairwise_common_vs_'):others.append(r['population_kind'].replace('pairwise_common_vs_',''))
        rows=[a for a in pop if a['pair']==r['pair'] and a['class']!='neither' and all(finite(a,mm,z) for mm in others for z in 'AB')]
        a=ranks_desc([float(x[f'{m}_score_A']) for x in rows]);bb=ranks_desc([float(x[f'{m}_score_B']) for x in rows])
        ids=sorted(range(len(rows)),key=lambda i:(max(a[i],bb[i]),a[i]+bb[i],rows[i]['global_ligand_entity_id']))
        k=math.ceil(.1*len(rows));top=[rows[i] for i in ids[:k]];nd=sum(x['class']=='dual' for x in rows)
        label=r['pair']+' '+m+' '+r['population_kind']
        for value,key in [(len(rows),'n_total'),(k,'k'),(nd,'n_dual')]+[(sum(x['class']==cls for x in top),key) for cls,key in [('dual','n_dual_top'),('A_only','n_A_only_top'),('B_only','n_B_only_top')]]:
            verify(label+' '+key,value,r[key],0)
        verify(label+' dual recall',sum(x['class']=='dual' for x in top)/nd,r['dual_recall'])
    # Validate the exact main-table baselines against their paired ranking rows.
    for r in read('submission/01_Main_Tables/Table_2.csv'):
        m={'S2':'M1','S3':'M2','S4':'M1b','S5':'M3'}[r['Scheme vs S1']]
        old=next(a for a in descriptive if a['pair']==r['Target pair'] and a['method']=='M0' and a['population_kind']=='pairwise_common_vs_'+m)
        new=next(a for a in descriptive if a['pair']==r['Target pair'] and a['method']==m and a['population_kind']=='pairwise_common_vs_M0')
        for key,row,field in [('Common N',old,'n_total'),('k',old,'k'),('Dual selected S1',old,'n_dual_top'),('Dual selected new',new,'n_dual_top'),('Dual recall S1',old,'dual_recall'),('Dual recall new',new,'dual_recall')]:verify('Table 2 '+r['Target pair']+' '+m+' '+key,r[key],row[field])
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'Independent_Numeric_Audit.json').write_text(json.dumps({'base_commit':'ce81b011bab449c6cf5e2b9664626ed808be5857','checks_passed':len(checks),'checks_failed':0,'new_scientific_analysis':False,'numeric_verification_from_deposited_scores':True,'checks':checks},indent=2)+'\n')
    print(json.dumps({'numeric_checks':len(checks),'failed':0}))

if __name__=='__main__':main()
