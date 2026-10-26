"""Check the supplied Results draft against accepted CSVs without rerunning science."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'submission/04_Editorial_Checks'
PR='results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv'
DD='results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv'
ST='results/jcim_stage2/'
CAP='results/current/chemical_diversity_cap5_20261023/'
SOURCE_COMMIT='a7b9d7521dff56ece18ec10a55734392f417e684'

def read(p): return list(csv.DictReader((ROOT/p).open(encoding='utf-8-sig',newline='')))

def main():
    OUT.mkdir(parents=True,exist_ok=True); audit=[]
    def check(section,claim,expected,actual,source,digits=None):
        equal=round(float(actual),digits)==float(expected) if digits is not None else actual==expected
        audit.append({'section':section,'claim':claim,'draft_value':expected,'source_value':actual,'display_decimals':digits,'status':'PASS' if equal else 'FAIL','source_path':source})
    def numrow(section,label,source,row,field,expected): check(section,label,expected,float(row[field]),source,4)
    pop=read('results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv'); pop=[r for r in pop if r['activity_eligible']=='1' and r['class']]
    check('3.1','Eligible observations',867,len(pop),'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    check('3.1','Distinct global ligand entities',845,len({r['global_ligand_entity_id'] for r in pop}),'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    for cls,n in [('dual',248),('A_only',256),('B_only',252),('neither',111)]: check('3.1',cls+' observations',n,sum(r['class']==cls for r in pop),'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    for pair,n in [('EGFR/HER2',109),('AChE/BChE',98),('JAK1/JAK2',110),('JAK1/TYK2',110),('PIK3CA/mTOR',110),('F2/F10',110),('PPARG/PPARA',110),('PPARA/PPARD',110)]: check('3.1',pair+' size',n,sum(r['pair']==pair for r in pop),'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    source=read('results/current/provenance_audit_20261023/CURRENT_867_SOURCE_COVERAGE_SUMMARY.csv'); sd={r['item']:int(r['n']) for r in source}
    for field,n in [('both_ends_raw_activity',867),('literature_complete',855),('document_identifier_missing',12),('cross_source_aggregation',81),('endpoint_heterogeneity',123),('value_or_label_mismatch',0)]:check('3.1',field,n,sd[field],'results/current/provenance_audit_20261023/CURRENT_867_SOURCE_COVERAGE_SUMMARY.csv')
    layers=read('analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv')
    for field,n in [('n_labeled_population',1004),('n_m0_directional_population',1000),('n_final_chemistry_oof',994),('n_m0_missing',4),('n_fold_unassigned_after_m0_filter',6)]:check('3.1',field,n,sum(int(r[field])for r in layers),'analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv')
    pr=read(PR);m0=[r for r in pr if r['method']=='M0']; by={(r['pair'],r['contrast']):r for r in m0}
    directional=[r for r in m0 if r['contrast']!='summary_min']; summaries=[r for r in m0 if r['contrast']=='summary_min']
    check('3.2','Directional estimates',16,len(directional),PR);check('3.2','Directions below 0.5',4,sum(float(r['auroc'])<0.5 for r in directional),PR)
    check('3.2','Directional range lower',0.2467,min(float(r['auroc'])for r in directional),PR,4);check('3.2','Directional range upper',0.7485,max(float(r['auroc'])for r in directional),PR,4)
    check('3.2','Summary range upper',0.6341,max(float(r['auroc'])for r in summaries),PR,4)
    check('3.2','Summary CIs containing 0.5',6,sum(float(r['ci_lo'])<=0.5<=float(r['ci_hi'])for r in summaries),PR)
    points={'EGFR/HER2':(0.6274,0.2467,0.2467),'JAK1/TYK2':(0.6475,0.4219,0.4219),'PPARA/PPARD':(0.6787,0.5107,0.5107),'F2/F10':(0.4082,0.3623,0.3623),'PIK3CA/mTOR':(0.7363,0.5801,0.5801),'AChE/BChE':(0.6341,0.6442,0.6341)}
    for pair,values in points.items():
        for contrast,value in zip(['dual_vs_A_only','dual_vs_B_only','summary_min'],values):numrow('3.2',pair+' '+contrast,PR,by[pair,contrast],'auroc',value)
    for pair,lo,hi in [('EGFR/HER2',0.1328,0.3783),('F2/F10',0.2305,0.4805),('PIK3CA/mTOR',0.4346,0.7129),('AChE/BChE',0.4607,0.7540)]:
        numrow('3.2',pair+' summary CI lower',PR,by[pair,'summary_min'],'ci_lo',lo);numrow('3.2',pair+' summary CI upper',PR,by[pair,'summary_min'],'ci_hi',hi)
    dd=read(DD);db={(r['pair'],r['delta']):r for r in dd}
    for pair,delta,pt,lo,hi in [('F2/F10','M2-M0',0.2471,0.0654,0.3896),('PPARA/PPARD','M2-M0',0.1582,0.0059,0.2979),('PPARG/PPARA','M2-M0',-0.2915,-0.4380,-0.1104),('PPARG/PPARA','M1-M0',-0.3089,-0.4688,-0.1313),('JAK1/TYK2','M3-M0',0.2471,0.0312,0.4082),('PIK3CA/mTOR','M1-M0',-0.0049,-0.1924,0.1680),('JAK1/JAK2','M1-M0',-0.1523,-0.3672,0.0605)]:
        for key,v in [('delta_summary_min',pt),('ci_lo',lo),('ci_hi',hi)]:numrow('3.2/3.7',pair+' '+delta+' '+key,DD,db[pair,delta],key,v)
    numrow('3.2','PIK3CA/mTOR S3-S1',DD,db['PIK3CA/mTOR','M2-M0'],'delta_summary_min',0.1240)
    for delta,n in [('M1b-M1',8),('M3-M0',7)]:check('3.2',delta+' CIs containing zero',n,sum(float(r['ci_lo'])<=0<=float(r['ci_hi'])for r in dd if r['delta']==delta),DD)
    neg=read(ST+'module_B_negative_class.csv');nb={(r['pair'],r['direction']):r for r in neg}
    for condition,n in [('increase',12),('decrease',4),('CI excludes zero',2)]:check('3.3',condition,n,sum((float(r['delta_negative'])>0 if condition=='increase' else float(r['delta_negative'])<0 if condition=='decrease' else not(float(r['delta_negative_ci_lo'])<=0<=float(r['delta_negative_ci_hi'])))for r in neg),ST+'module_B_negative_class.csv')
    for pair,arm,a,b,d,lo,hi in [('EGFR/HER2','D_vs_B',0.2467,0.7470,0.5004,0.3025,0.6819),('JAK1/TYK2','D_vs_B',0.4219,0.8058,0.3839,0.1953,0.5693),('AChE/BChE','D_vs_B',0.6442,0.5679,-0.0763,-0.2019,0.0469)]:
        for key,v in [('auc_single',a),('auc_neither',b),('delta_negative',d),('delta_negative_ci_lo',lo),('delta_negative_ci_hi',hi)]:numrow('3.3',pair+' '+arm+' '+key,ST+'module_B_negative_class.csv',nb[pair,arm],key,v)
    for pair,arm,d in [('PPARA/PPARD','D_vs_A',-0.0269),('PPARA/PPARD','D_vs_B',-0.0241),('PIK3CA/mTOR','D_vs_A',0.0315),('PIK3CA/mTOR','D_vs_B',0.1164)]:numrow('3.3',pair+' '+arm+' change',ST+'module_B_negative_class.csv',nb[pair,arm],'delta_negative',d)
    chem=read(ST+'module_D_chemistry_units.csv');cb={(r['pair'],r['arm']):r for r in chem}
    for key,lo,hi in [('auroc_physchem',0.2939,0.9014),('auroc_ecfp4',0.6633,0.9678),('delta_ecfp4_plus_m0',-0.0164,0.0363)]:
        check('3.4',key+' minimum',lo,min(float(r[key])for r in chem),ST+'module_D_chemistry_units.csv',4);check('3.4',key+' maximum',hi,max(float(r[key])for r in chem),ST+'module_D_chemistry_units.csv',4)
    check('3.4','ECFP4 exceeds descriptors',15,sum(float(r['auroc_ecfp4'])>float(r['auroc_physchem'])for r in chem),ST+'module_D_chemistry_units.csv')
    check('3.4','ECFP4 exceeds matched Vina',16,sum(float(r['auroc_ecfp4'])>float(r['auroc_SUPPORTING_CHEMISTRY_MATCHED_M0'])for r in chem),ST+'module_D_chemistry_units.csv')
    check('3.4','Positive Vina increments',5,sum(float(r['delta_ecfp4_plus_m0'])>0 for r in chem),ST+'module_D_chemistry_units.csv')
    check('3.4','Negative Vina increments',11,sum(float(r['delta_ecfp4_plus_m0'])<0 for r in chem),ST+'module_D_chemistry_units.csv')
    for pair,arm,field,value in [('EGFR/HER2','D_vs_B','auroc_ecfp4',0.8790),('EGFR/HER2','D_vs_B','auroc_SUPPORTING_CHEMISTRY_MATCHED_M0',0.2408),('JAK1/JAK2','D_vs_A','auroc_ecfp4',0.9185),('JAK1/JAK2','D_vs_B','auroc_ecfp4',0.9678),('PIK3CA/mTOR','D_vs_A','auroc_ecfp4',0.8213),('PIK3CA/mTOR','D_vs_B','auroc_ecfp4',0.7285),('F2/F10','D_vs_B','auroc_ecfp4',0.6633),('F2/F10','D_vs_B','auroc_ecfp4_m0',0.6996),('PIK3CA/mTOR','D_vs_A','delta_ecfp4_plus_m0',0.0215),('PIK3CA/mTOR','D_vs_B','delta_ecfp4_plus_m0',-0.0029),('EGFR/HER2','D_vs_A','delta_ecfp4_plus_m0',-0.0164)]:numrow('3.4',pair+' '+arm+' '+field,ST+'module_D_chemistry_units.csv',cb[pair,arm],field,value)
    pocket=read(ST+'module_C_wrong_pocket.csv');pb={(r['pair'],r['arm']):r for r in pocket}
    for label,n in [('positive',10),('negative',5),('equal',1),('CI contains zero',14)]:check('3.5',label,n,sum((float(r['delta_pocket'])>0 if label=='positive' else float(r['delta_pocket'])<0 if label=='negative' else float(r['delta_pocket'])==0 if label=='equal' else float(r['delta_pocket_ci_lo'])<=0<=float(r['delta_pocket_ci_hi']))for r in pocket),ST+'module_C_wrong_pocket.csv')
    check('3.5','Additional pocket exclusions',0,sum(int(r['extra_drop_vs_directional_m0'])for r in pocket),ST+'module_C_wrong_pocket.csv')
    for pair,arm,a,b,d,lo,hi in [('AChE/BChE','D_vs_A',0.6341,0.4681,0.1659,0.0548,0.2830),('PPARA/PPARD','D_vs_A',0.6787,0.4434,0.2354,0.0693,0.4043),('PPARG/PPARA','D_vs_B',0.6021,0.6909,-0.0889,-0.2563,0.0869)]:
        for key,v in [('auc_correct',a),('auc_wrong',b),('delta_pocket',d),('delta_pocket_ci_lo',lo),('delta_pocket_ci_hi',hi)]:numrow('3.5',pair+' '+arm+' '+key,ST+'module_C_wrong_pocket.csv',pb[pair,arm],key,v)
    numrow('3.5','PIK3CA/mTOR noncorresponding AUROC',ST+'module_C_wrong_pocket.csv',pb['PIK3CA/mTOR','D_vs_A'],'auc_wrong',0.6519)
    seed=read(ST+'module_E_five_seed_summary.csv');sb={r['pair']:r for r in seed if r['metric']=='summary_min'}
    for pair,lo,hi in [('EGFR/HER2',0.2305,0.2467),('F2/F10',0.3594,0.3682),('PIK3CA/mTOR',0.5596,0.5869),('PPARA/PPARD',0.4800,0.5151)]:
        numrow('3.6',pair+' seed minimum',ST+'module_E_five_seed_summary.csv',sb[pair],'min',lo);numrow('3.6',pair+' seed maximum',ST+'module_E_five_seed_summary.csv',sb[pair],'max',hi)
    alt=read(ST+'module_F_alt_receptor.csv');check('3.6','Alternative conditions',6,len(alt),ST+'module_F_alt_receptor.csv');check('3.6','Alternative CIs contain zero',6,sum(float(r['delta_summary_min_ci_lo'])<=0<=float(r['delta_summary_min_ci_hi'])for r in alt),ST+'module_F_alt_receptor.csv')
    for pdb,before,after,d,lo,hi in [('ALT_4EY6',0.6184,0.5670,-0.0514,-0.1852,0),('ALT_3SHC',0.3828,0.4333,0.0505,0,0.1441),('ALT_6KAX',0.5107,0.4438,-0.0669,-0.1484,0.0044)]:
        r=next(r for r in alt if r['alt_id']==pdb)
        for key,v in [('summary_min_primary',before),('summary_min_alt',after),('delta_summary_min',d),('delta_summary_min_ci_lo',lo),('delta_summary_min_ci_hi',hi)]:numrow('3.6',pdb+' '+key,ST+'module_F_alt_receptor.csv',r,key,v)
    cc=read(CAP+'CAP5_CLASS_COUNTS.csv');check('3.6','CAP5 retained',861,sum(int(r['n_cap5'])for r in cc),CAP+'CAP5_CLASS_COUNTS.csv');check('3.6','CAP5 excluded',6,sum(int(r['n_excluded'])for r in cc),CAP+'CAP5_CLASS_COUNTS.csv')
    cp=read(CAP+'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv');cpb={(r['pair'],r['method'],r['contrast']):r for r in cp}
    for pair,pt,lo,hi in [('JAK1/JAK2',0.6315,0.4871,0.7597),('PPARA/PPARD',0.4969,0.3427,0.6354)]:
        for key,v in [('auroc',pt),('ci_lo',lo),('ci_hi',hi)]:numrow('3.6',pair+' CAP5 '+key,CAP+'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv',cpb[pair,'M0','summary_min'],key,v)
    cchem=read(CAP+'CAP5_MODULE_D_CHEMISTRY_UNITS.csv');check('3.6','All ECFP4 changes below 0.02',16,sum(abs(float(r['auroc_ecfp4'])-float(cb[r['pair'],r['arm']]['auroc_ecfp4']))<0.02 for r in cchem),CAP+'CAP5_MODULE_D_CHEMISTRY_UNITS.csv')
    strict=read('extension_jcim_validation/pik3ca_mtor_q110_v42/Q110_STRICT_SUBSET_METRICS.csv'); ss={r['contrast']:r for r in strict if r['method']=='M0'}
    for contrast,pt in [('dual_vs_A_only',0.7900),('dual_vs_B_only',0.4571),('summary_min',0.4571)]:numrow('3.6','Strict S1 '+contrast,'extension_jcim_validation/pik3ca_mtor_q110_v42/Q110_STRICT_SUBSET_METRICS.csv',ss[contrast],'auroc',pt)
    for key,v in [('ci_lo',0.2743),('ci_hi',0.6429)]:numrow('3.6','Strict S1 summary '+key,'extension_jcim_validation/pik3ca_mtor_q110_v42/Q110_STRICT_SUBSET_METRICS.csv',ss['summary_min'],key,v)
    rank=read(ST+'module_I_method_compare.csv');check('3.7','Method comparisons',32,len(rank),ST+'module_I_method_compare.csv');check('3.7','Common population minimum',79,min(int(r['n_rank_common'])for r in rank),ST+'module_I_method_compare.csv');check('3.7','Common population maximum',97,max(int(r['n_rank_common'])for r in rank),ST+'module_I_method_compare.csv')
    counts=Counter(r['concordance']for r in rank)
    for label,n in [('CONCORDANT_IMPROVE',14),('CONCORDANT_WORSEN',10),('DISCORDANT',6),('NO_CHANGE',2)]:check('3.7',label,n,counts[label],ST+'module_I_method_compare.csv')
    rt=Counter(r['concordance']for r in rank if r['method']=='M2');check('3.7','RTMScore concordant improvement',5,rt['CONCORDANT_IMPROVE'],ST+'module_I_method_compare.csv');check('3.7','RTMScore concordant worsening',3,rt['CONCORDANT_WORSEN'],ST+'module_I_method_compare.csv')
    desc=read(ST+'module_I_descriptive.csv');db={(r['pair'],r['method'],r['population_kind']):r for r in desc}
    for pair,m,before,after in [('F2/F10','M2',4,7),('PPARA/PPARD','M2',3,6),('PIK3CA/mTOR','M2',4,7),('PPARG/PPARA','M2',6,4),('PIK3CA/mTOR','M1',4,6),('JAK1/JAK2','M1',6,7)]:
        check('3.7',pair+' '+m+' baseline dual count',before,int(db[pair,'M0','pairwise_common_vs_'+m]['n_dual_top']),ST+'module_I_descriptive.csv');check('3.7',pair+' '+m+' new dual count',after,int(db[pair,m,'pairwise_common_vs_M0']['n_dual_top']),ST+'module_I_descriptive.csv')
    with (OUT/'RESULTS_NUMERIC_AUDIT_20261027.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
    summary={'source_commit':SOURCE_COMMIT,'checks':len(audit),'passed':sum(r['status']=='PASS'for r in audit),'failed':sum(r['status']=='FAIL'for r in audit),'scientific_recomputation':False}
    (OUT/'RESULTS_NUMERIC_AUDIT_20261027.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary));print('FAILURES',[r for r in audit if r['status']=='FAIL'])

if __name__=='__main__':main()
