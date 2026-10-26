"""Create numbered submission tables from accepted results; no scientific reruns.

Usage: python build_submission_tables.py --project-root Dual_Target_Docking
The accompanying JS builder requires @oai/artifact-tool in the Codex runtime.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import zipfile
from collections import Counter

PAIRS = ['EGFR/HER2','JAK1/JAK2','JAK1/TYK2','PIK3CA/mTOR','AChE/BChE','F2/F10','PPARG/PPARA','PPARA/PPARD']
METHODS = ['M0','M1','M2','M1b','M3']
SCHEMES = dict(zip(METHODS, ['S1','S2','S3','S4','S5']))
NAMES = ['Vina mode 1','GNINA fixed-pose CNNscore','RTMScore fixed pose','GNINA best saved CNNscore','GNINA independent docking']
BASE_COMMIT = '1090ed86d8cb553ee365e52bc51c34e48c0e5f6c'
RP = 'reruns/UNIFORM_RERUN_V4_2_20260921/'
PR = 'results/current/provenance_audit_20261023/'
CAP = 'results/current/chemical_diversity_cap5_20261023/'
Q = 'extension_jcim_validation/pik3ca_mtor_q110_v42/'

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def scalar(value):
    if value is None or value == '': return None
    if isinstance(value, (int,float)): return value
    try:
        f = float(value)
        if math.isfinite(f): return int(f) if f.is_integer() else f
    except (ValueError,TypeError): pass
    return value

def delta_name(value):
    return value.replace('M1b','S4').replace('M0','S1').replace('M1','S2').replace('M2','S3').replace('M3','S5')

def csv_bytes(headers, rows):
    s = io.StringIO(newline='')
    w = csv.writer(s, lineterminator='\n'); w.writerow(headers); w.writerows(rows)
    return s.getvalue().encode('utf-8')

def display(value):
    if value is None: return 'NA'
    if isinstance(value,float): return f'{value:.4f}'
    return str(value).replace('|','\\|').replace('\n',' ')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[2])
    args=parser.parse_args(); root=args.project_root.resolve(); out=root/'submission'
    for d in ['01_Main_Tables','02_SI_Tables','03_Supplementary_Data','04_Editorial_Checks']:
        (out/d).mkdir(parents=True,exist_ok=True)
    sources=set(); checks=[]; raw_extra={}; tables=[]
    def rows(p):
        sources.add(p); return read_csv(root/p)
    def check(name, condition, detail=''):
        checks.append({'check':name,'status':'PASS' if condition else 'FAIL','detail':detail})
        if not condition: raise ValueError(name+': '+str(detail))
    def sortrows(data, field='method'):
        return sorted(data,key=lambda r:(PAIRS.index(r['pair']), METHODS.index(r[field]) if r.get(field) in METHODS else 0, r.get('arm',r.get('direction',r.get('contrast','')))))
    def panel(label,title,headers,data,note=''):
        return {'label':label,'title':title,'headers':headers,'rows':[[scalar(v) for v in r] for r in data],'note':note}
    def table(number,title,cn,first,panels,note=''):
        tables.append({'number':number,'title':title,'title_cn':cn,'first_citation':first,'panels':panels,'note':note})
    pop=rows('results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    eligible=[r for r in pop if r['activity_eligible']=='1' and r['class']]
    check('formal eligible observations',len(eligible)==867)
    check('formal class totals',Counter(r['class'] for r in eligible)==Counter(dual=248,A_only=256,B_only=252,neither=111))
    primary=rows('results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv')
    deltas=rows('results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv')
    desc=rows('results/jcim_stage2/module_I_descriptive.csv')
    comp=rows('results/jcim_stage2/module_I_method_compare.csv')
    check('ranking population counts',all(int(r['n_total'])==sum(int(r[k]) for k in ['n_dual','n_A_only','n_B_only']) for r in desc))
    check('ranking top-list sizes',all(int(r['k'])==math.ceil(0.10*int(r['n_total']))==sum(int(r[k]) for k in ['n_dual_top','n_A_only_top','n_B_only_top']) for r in desc))
    check('ranking recall denominators',all(math.isclose(float(r['dual_recall']),int(r['n_dual_top'])/int(r['n_dual']),abs_tol=1e-12) for r in desc))
    mdict={(r['pair'],r['method'],r['contrast']):r for r in primary}
    ddict={(r['pair'],r['method'],r['population_kind']):r for r in desc}
    t1=[]
    for pair in PAIRS:
        ps=[r for r in eligible if r['pair']==pair]; c=Counter(r['class'] for r in ps)
        pa,pb={(r['pdb_A'],r['pdb_B']) for r in ps}.pop()
        t1.append([pair,pa,pb,*[c[k] for k in ['dual','A_only','B_only','neither']],len(ps)])
    table('1','Target pairs and experimental activity classes','靶点对与实验活性类别','Methods 2.2',[
        panel('','', ['Target pair','PDB A','PDB B','Dual (n)','A-only (n)','B-only (n)','Neither (n)','Total (n)'],t1)],
        'Counts are target-pair observations. A molecule may occur in more than one pair. Final activity labels use pChEMBL >= 6.0. Two ineligible observations and the AB_040 alias are excluded.')
    t2=[]
    for pair in PAIRS:
        for method in METHODS[1:]:
            a=ddict[pair,'M0','pairwise_common_vs_'+method]; b=ddict[pair,method,'pairwise_common_vs_M0']
            check('matched ranking population '+pair+' '+method,(a['n_total'],a['k'])==(b['n_total'],b['k']))
            t2.append([pair,SCHEMES[method],a['n_total'],a['k'],a['n_dual_top'],b['n_dual_top'],a['dual_recall'],b['dual_recall']])
    table('2','Dual-ligand recovery in matched top-10% candidate lists','共同成员上前 10% 候选名单的双靶配体检出','Results: joint ranking',[
        panel('','',['Target pair','Scheme vs S1','Common N','k','Dual selected S1','Dual selected new','Dual recall S1','Dual recall new'],t2)],
        'Every comparison uses its own common population with valid scores on both targets. k = ceil(0.10 N). S4 versus S1 is a supporting comparison. Class composition and directional changes are given in Table S10.')

    selection=rows(PR+'SEVEN_PAIR_SELECTION_PROVENANCE.csv'); sel={r['pair']:r for r in selection}
    selection_data=[]
    for pair in PAIRS:
        n=sum(r['pair']==pair for r in eligible)
        if pair=='PIK3CA/mTOR':
            selection_data.append([pair,'6.0','Fixed 110-molecule membership','32 / 32 / 32 / 14',n,'Not independently replayed in this compilation'])
        else:
            r=sel[pair]
            rule='6.5 / 5.5' if 'strict' in r['activity_rule_at_sampling'] else '6.0'
            chem='Murcko cap 5 within class' if pair=='EGFR/HER2' else ('ChEMBL ID-prefix diversity key in deposited script' if pair=='AChE/BChE' else 'Small-molecule / size filters; no Murcko cap')
            selection_data.append([pair,rule,chem,r['quota'],n,'DOCUMENTED_NOT_REPLAYED'])
    cross=rows(PR+'CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv')
    coverage=[]
    for pair in PAIRS:
        rr=[r for r in cross if r['pair']==pair]
        flags=['both_ends_raw_activity','literature_complete','document_identifier_missing','cross_source_or_endpoint_mix','replicate_threshold_discordance','value_or_label_mismatch']
        coverage.append([pair,len(rr),*[sum(r[f]=='YES' for r in rr) for f in flags]])
    check('provenance crosswalk count',len(cross)==867)
    check('provenance members and classes', {(r['pair'],r['panel_id'],r['class']) for r in cross}=={(r['pair'],r['panel_id'],r['class']) for r in eligible})
    table('S1','Evaluation-set construction and activity provenance','评价集构建与实验活性来源','Methods 2.2',[
        panel('A','Selection conditions',['Target pair','Selection threshold','Chemical selection record','Deposited quota / composition','Eligible n','Selection replay status'],selection_data),
        panel('B','Activity-source coverage',['Target pair','Observations','Both ends traced','Complete identifiers','Missing identifiers','Source / endpoint mix','Threshold-discordant repeats','Value / label mismatch'],coverage)],
        'All final labels use pChEMBL >= 6.0. Provenance flags overlap. Complete identifiers require DOI, PubMed or patent identifiers for all tied maxima at both endpoints. Across all pairs, cross-document aggregation is present in 81 observations and endpoint heterogeneity in 123 (union 168). The full source crosswalk, raw activity rows, selection evidence and exceptions are included in Supplementary_Data.zip. CAP5 is a post hoc sensitivity analysis and does not redefine initial selection.')

    registry=rows(RP+'00_protocol/receptor_registry.csv')
    prep={r['pdb_id']:r for r in rows(RP+'13_qa/receptor_preparation_manifest.csv')}
    altprep={r['pdb_id']:r for r in rows(RP+'13_qa/alt_receptor_preparation_manifest.csv')}
    altboxes={r['pdb_id']:r for r in rows(RP+'13_qa/alt_box_manifest.csv')}
    alt=rows('results/jcim_stage2/module_F_alt_receptor.csv')
    alt=[r for r in alt if r.get('pair')]
    check('current alternative conditions',len(alt)==6 and {r['pair'] for r in alt}=={'AChE/BChE','F2/F10','PPARA/PPARD'})
    contexts={'4L23':'p85 removed','4JT6':'Second mTOR copy and mLST8 removed','4EY7':'AChE chains A and B retained','4UDW':'Heavy chain, light chain and hirudin retained','2JKH':'Heavy and light chains retained','9V8H':'PPARG and PG08-NL peptide retained','5U3Q':'PPARD chain B retained','6N7A':'Chains A and B retained; PTR converted to TYR','8BXH':'PTR converted to TYR','3LXP':'PTR converted to TYR'}
    structures=[]; seen=set()
    for r in registry:
        pdb=r['pdb_id']
        if pdb in seen: continue
        seen.add(pdb); p=RP+'03_boxes/'+pdb+'_box.json'; sources.add(p); b=json.loads((root/p).read_text())
        check('receptor prepared '+pdb,prep[pdb]['status']=='PREPARED')
        structures.append(['Primary',r['target'],pdb,r['protein_auth_asym_id'],r['cognate_ccd'],r['box_defining_instance'],contexts.get(pdb,'Selected protein chain retained'),*[b[k] for k in ['center_x','center_y','center_z','size_x','size_y','size_z']]])
    for r in sorted(alt,key=lambda x:(PAIRS.index(x['pair']),x['replaced_side'])):
        pdb=r['alt_id'].replace('ALT_',''); b=altboxes[pdb]; ap=altprep[pdb]
        target=r['pair'].split('/')[0 if r['replaced_side']=='A' else 1]
        context={'4EY6':'AChE chains A and B retained','3SHC':'Heavy chain, light chain and hirudin retained','2Y5F':'Heavy and light chains retained'}.get(pdb,'Selected protein chain retained')
        check('alternative prepared '+pdb,ap['status']=='PREPARED')
        structures.append(['Alternative',target,pdb,ap['protein_auth'],b['ccd'],b['reference_instance'],context,*[b[k] for k in ['center_x','center_y','center_z','size_x','size_y','size_z']]])
    check('14 primary plus 6 current alternative receptors',len(seen)==14 and len(structures)==20)
    table('S2','Receptor structures and docking boxes','受体结构与对接盒参数','Methods 2.3',[
        panel('','',['Role','Target','PDB','Auth chains','Cognate CCD','Cognate instance','Retained / removed components','Center x (A)','Center y (A)','Center z (A)','Size x (A)','Size y (A)','Size z (A)'],structures)],
        'Coordinates and dimensions are in angstrom. Waters and cocrystal ligands are removed from the receptor. Box centers are bounding-box midpoints; padding is 5 A with minimum edge 20 A. Chain identifiers are author chain identifiers. Formal prepared files and preparation exceptions, including observed-fragment templates, are documented in the accompanying source files. Alternative structures cover only the six conditions used in the current analysis; no Q110 alternative-receptor result is included.')
    redock=rows(RP+'13_qa/redocking_master.csv')
    qc=[]
    for r in registry:
        pdb=r['pdb_id']
        if any(x[1]==pdb for x in qc): continue
        rr=[x for x in redock if x['pdb_id']==pdb]; valid=[x for x in rr if x['status']=='OK']
        top=[float(x['top1_rmsd_A']) for x in valid if x['top1_rmsd_A']]; best=[float(x['min_saved_pose_rmsd_A']) for x in valid if x['min_saved_pose_rmsd_A']]
        qc.append([r['target'],pdb,r['cognate_ccd'],len(rr),len(valid),statistics.median(top) if top else None,min(top) if top else None,max(top) if top else None,statistics.median(best) if best else None,min(best) if best else None,max(best) if best else None,sum(x<2 for x in top),sum(x<2 for x in best),'/'.join(str(x) for x in sorted({int(x['n_modes_returned']) for x in valid}))])
    check('cocrystal redocking jobs',len(redock)==70 and all(len([r for r in redock if r['pdb_id']==p])==5 for p in seen))
    table('S3','Cocrystal-ligand redocking across five random seeds','五随机种子的共晶配体重对接','Methods 2.3',[
        panel('','',['Target','PDB','Cognate CCD','Seeds tested','Successful jobs','Top-1 median RMSD (A)','Top-1 min RMSD (A)','Top-1 max RMSD (A)','Best-saved median RMSD (A)','Best-saved min RMSD (A)','Best-saved max RMSD (A)','Top-1 RMSD <2 A (n)','Best-saved RMSD <2 A (n)','Modes returned'],qc)],
        'Seeds: 17, 29, 42, 71 and 101. Symmetry-aware heavy-atom RMSD is calculated in the receptor coordinate frame. Top-ranked-pose success and success among saved poses are separate outcomes. The summary uses recorded RMSDs without recomputing coordinates; every seed-level record is in Supplementary_Data.zip.')

    layers=rows('analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv'); missing=rows('results/jcim_stage2/module_H_missingness.csv')
    population_data=[[r[k] for k in ['pair','arm','n_labeled_population','n_m0_directional_population','n_m0_missing','n_fold_unassigned_after_m0_filter','n_final_chemistry_oof']] for r in sortrows(layers)]
    miss_data=[]
    for pair in PAIRS:
        for m in METHODS:
            rr=[r for r in missing if r['pair']==pair and r['method']==m]
            get=lambda target,k:sum(int(r[k]) for r in rr if r['target']==target)
            miss_data.append([pair,SCHEMES[m],get('A','expected'),get('A','valid'),get('A','missing'),get('B','expected'),get('B','valid'),get('B','missing')])
    check('population-layer totals',[sum(int(r[k]) for r in layers) for k in ['n_labeled_population','n_m0_directional_population','n_final_chemistry_oof']]==[1004,1000,994])
    check('missingness reconciliation',all(int(r['expected'])==int(r['valid'])+int(r['missing']) for r in missing))
    table('S4','Analysis populations and score availability','分析人群与评分可用性','Methods 2.5',[
        panel('A','Directional and chemistry populations',['Target pair','Direction','Labeled n','S1 directional n','S1 missing n','Unassigned fold after S1','Chemistry OOF n'],population_data),
        panel('B','Score availability across activity classes',['Target pair','Scheme','Expected A','Valid A','Missing A','Expected B','Valid B','Missing B'],miss_data)],
        'Panel A counts ligand-direction records. Panel B aggregates the four activity classes; class-specific missingness, exclusions and ligand identifiers remain in the complete source files. D_vs_A means dual versus A-only using target B; D_vs_B uses target A. Missing fold assignments are excluded from chemistry evaluation. Paired method comparisons use their own complete-case populations, reported in Table S5B.')
    metric_headers=['Target pair','Scheme','Dual n','A-only n','B-only n','AUROC B','B CI lower','B CI upper','AUROC A','A CI lower','A CI upper','Summary min','Summary CI lower','Summary CI upper']
    def metric_matrix(data):
        by={(r['pair'],r['method'],r['contrast']):r for r in data}; matrix=[]
        for pair in PAIRS:
            if not any(r['pair']==pair for r in data): continue
            for method in METHODS:
                a=by[pair,method,'dual_vs_A_only']; b=by[pair,method,'dual_vs_B_only']; s=by[pair,method,'summary_min']
                matrix.append([pair,SCHEMES[method],a['n_dual'],a['n_negative'],b['n_negative'],*[a[k] for k in ['auroc','ci_lo','ci_hi']],*[b[k] for k in ['auroc','ci_lo','ci_hi']],*[s[k] for k in ['auroc','ci_lo','ci_hi']]])
        return matrix
    delta_order=['M1-M0','M2-M0','M1b-M1','M3-M0']
    delta_matrix=[]
    for pair in PAIRS:
        for d in delta_order:
            r=next(r for r in deltas if r['pair']==pair and r['delta']==d)
            delta_matrix.append([pair,delta_name(d),*[r[k] for k in ['n_dual_four_sided','n_A_only_both_score_B','n_B_only_both_score_A','summary_min_x','summary_min_y','delta_summary_min','ci_lo','ci_hi']]])
    check('primary metrics and paired deltas coverage',len(primary)==120 and len(deltas)==32)
    table('S5','Directional discrimination and paired scoring-workflow comparisons','方向性判别及评分工作流配对比较','Methods 2.5',[
        panel('A','Method-specific estimates',metric_headers,metric_matrix(primary)),
        panel('B','Paired differences in summary min',['Target pair','Contrast','Common dual n','Common A-only n','Common B-only n','Summary minuend','Summary subtrahend','Difference','Difference CI lower','Difference CI upper'],delta_matrix)],
        'AUROC B compares dual with A-only; AUROC A compares dual with B-only. S1 Vina; S2 fixed-pose GNINA; S3 fixed-pose RTMScore; S4 maximum saved-pose CNNscore; S5 independent GNINA. All scores have higher-better direction. Each paired contrast is recalculated on its own shared members; it need not equal subtraction of the method-specific estimates in panel A. CI: 95% class-stratified percentile bootstrap, 10,000 resamples, seed 271828. S4-S2 isolates saved-pose reranking; the supporting S4-S1 link in Table S10 has no added confidence interval.')
    neg=rows('results/jcim_stage2/module_B_negative_class.csv'); pocket=rows('results/jcim_stage2/module_C_wrong_pocket.csv')
    table('S6','Negative-class and pocket-correspondence controls','负类组成与口袋对应对照','Methods 2.6.1; reused in 2.6.3',[
        panel('A','Negative-class replacement',['Target pair','Direction','Dual n','Single-target n','Neither n','AUROC single-target','AUROC neither','Difference','Difference CI lower','Difference CI upper'],[[r[k] for k in ['pair','direction','n_dual','n_single','n_neither','auc_single','auc_neither','delta_negative','delta_negative_ci_lo','delta_negative_ci_hi']] for r in sortrows(neg)]),
        panel('B','Corresponding and noncorresponding pockets',['Target pair','Direction','Dual n','Common negative n','AUROC corresponding','AUROC noncorresponding','Difference','Difference CI lower','Difference CI upper','Additional exclusions'],[[r[k] for k in ['pair','arm','n_dual','n_negative_both_finite','auc_correct','auc_wrong','delta_pocket','delta_pocket_ci_lo','delta_pocket_ci_hi','extra_drop_vs_directional_m0']] for r in sortrows(pocket)])],
        'Panel A difference = neither minus single-target negative class. Panel B difference = corresponding minus noncorresponding pocket. All estimates use S1. Individual AUROC confidence intervals and ligand identifiers are retained in the complete source CSVs.')
    chem=rows('results/jcim_stage2/module_D_chemistry_units.csv')
    table('S7','Out-of-fold performance of ligand-chemistry models','配体化学模型的折外判别表现','Methods 2.6.2',[
        panel('','',['Target pair','Direction','Common OOF n','Physicochemical AUROC','ECFP4 AUROC','ECFP4 + S1 AUROC','Matched S1 AUROC','ECFP4 + S1 increment'],[[r[k] for k in ['pair','arm','n_final_chemistry_oof','auroc_physchem','auroc_ecfp4','auroc_ecfp4_m0','auroc_SUPPORTING_CHEMISTRY_MATCHED_M0','delta_ecfp4_plus_m0']] for r in sortrows(chem)])],
        'All models and S1 are evaluated on identical members within each direction. The increment is AUROC(ECFP4 + S1) minus AUROC(ECFP4). Existing frozen five-fold assignments are retained; Q110 uses sorted canonical Murcko scaffolds assigned cyclically to folds 0-4. Standardization is fitted in training folds. Full descriptors, scaffolds, similarities, folds and OOF predictions are supplied as data files; no new model is fitted during table preparation.')
    seeds=rows('results/jcim_stage2/module_E_five_seed_summary.csv')
    table('S8','Random-seed and receptor-structure sensitivity','随机种子与受体结构敏感性','Methods 2.7',[
        panel('A','Five-seed distributions',['Target pair','Metric','Seed 42','Median','Minimum','Maximum','IQR','Finite seeds'],[[r[k] for k in ['pair','metric','seed42','median','min','max','iqr','n_finite_seeds']] for r in sorted(seeds,key=lambda r:(PAIRS.index(r['pair']),['auc_B','auc_A','summary_min'].index(r['metric'])))]),
        panel('B','One-receptor substitutions',['Target pair','Alternative PDB','Replaced side','Common dual n','Common A-only n','Common B-only n','Original summary min','Alternative summary min','Summary difference','Difference CI lower','Difference CI upper','Affected AUROC difference','Unreplaced AUROC difference'],[[r['pair'],r['alt_id'].replace('ALT_',''),r['replaced_side'],*[r[k] for k in ['n_dual_common','n_A_only_common','n_B_only_common','summary_min_primary','summary_min_alt','delta_summary_min','delta_summary_min_ci_lo','delta_summary_min_ci_hi','delta_AUC_affected','delta_AUC_unreplaced']]] for r in sortrows(alt)])],
        'Five Vina seeds: 17, 29, 42, 71 and 101. Each seed uses its own valid population. Receptor comparisons use paired common members and replace only one target. Full direction-specific estimates, confidence intervals and members are in the source data. PIK3CA/mTOR Q110 has no alternative-receptor calculation.')
    capcount=rows(CAP+'CAP5_CLASS_COUNTS.csv'); capmetrics=rows(CAP+'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv'); capcompare=rows(CAP+'CAP5_PRIMARY_VS_ORIGINAL.csv')
    capdict={(r['pair'],r['method'],r['contrast']):r for r in capmetrics}
    counts=[]; capsummary=[]
    for pair in PAIRS:
        rr=[r for r in capcount if r['pair']==pair]
        counts.append([pair,sum(int(r['n_original']) for r in rr),sum(int(r['n_cap5']) for r in rr),sum(int(r['n_excluded']) for r in rr),*[sum(int(r['n_excluded']) for r in rr if r['class']==c) for c in ['dual','A_only','B_only','neither']]])
        a=mdict[pair,'M0','summary_min']; b=capdict[pair,'M0','summary_min']
        capsummary.append([pair,*[a[k] for k in ['auroc','ci_lo','ci_hi']],*[b[k] for k in ['auroc','ci_lo','ci_hi']],float(b['auroc'])-float(a['auroc'])])
    check('CAP5 membership totals',sum(int(r['n_cap5']) for r in capcount)==861 and sum(int(r['n_excluded']) for r in capcount)==6)
    strict=rows(Q+'Q110_STRICT_SUBSET_METRICS.csv')
    capchem=rows(CAP+'CAP5_MODULE_D_CHEMISTRY_UNITS.csv')
    original_chem={(r['pair'],r['arm']):r for r in chem}
    chem_sensitivity=[]
    for r in sortrows(capchem):
        a=original_chem[r['pair'],r['arm']]
        chem_sensitivity.append([r['pair'],r['arm'],a['n_final_chemistry_oof'],r['n_final_chemistry_oof'],a['auroc_ecfp4'],r['auroc_ecfp4'],float(r['auroc_ecfp4'])-float(a['auroc_ecfp4']),a['delta_ecfp4_plus_m0'],r['delta_ecfp4_plus_m0']])
    table('S9','Chemical-diversity and activity-threshold sensitivity','化学多样性与活性阈值敏感性','Methods 2.7',[
        panel('A','CAP5 membership',['Target pair','Original n','CAP5 n','Excluded n','Excluded dual','Excluded A-only','Excluded B-only','Excluded neither'],counts),
        panel('B','Vina summary min before and after CAP5',['Target pair','Primary summary min','Primary CI lower','Primary CI upper','CAP5 summary min','CAP5 CI lower','CAP5 CI upper','Point change'],capsummary),
        panel('C','PIK3CA/mTOR strict 6.5/5.5 subset',metric_headers,metric_matrix(strict)),
        panel('D','Chemistry-model sensitivity to CAP5',['Target pair','Direction','Original OOF n','CAP5 OOF n','Original ECFP4 AUROC','CAP5 ECFP4 AUROC','ECFP4 AUROC change','Original Vina increment','CAP5 Vina increment'],chem_sensitivity)],
        'CAP5 is post hoc supporting sensitivity: at most five ligands per nonempty Murcko scaffold within pair and activity class; acyclic ligands are exempt. Six ligands are excluded without replacement. Complete CAP5 method comparisons, chemical models, controls and ranking outputs are retained in Supplementary_Data.zip. The PPARA/PPARD S5 dual-versus-A-only CI changes from [0.4921875, 0.767578125] to [0.5030241935, 0.7772177419]; this boundary change is retained, without a general significance claim. The strict Q110 subset has 56 members (25 dual, 12 A-only, 14 B-only, 5 neither). Main and subset populations differ; no new paired CI for their difference is calculated.')
    concord=rows('results/jcim_stage2_figures/module_I_concordance_counts.csv')
    descriptive=[r for r in desc if r['population_kind']=='method_specific']
    table('S10','Joint-ranking composition and directional concordance','联合排序候选组成与方向性变化关系','Methods 2.8',[
        panel('A','Method-specific top-10% lists',['Target pair','Scheme','N','k','Dual in population','Dual selected','A-only selected','B-only selected','Dual fraction','Single-target fraction','Dual recall','Population dual fraction'],[[r['pair'],SCHEMES[r['method']],*[r[k] for k in ['n_total','k','n_dual','n_dual_top','n_A_only_top','n_B_only_top','top10_dual_fraction','top10_single_target_fraction','dual_recall','expected_dual_fraction']]] for r in sortrows(descriptive)]),
        panel('B','Paired changes versus Vina',['Target pair','Scheme vs S1','Ranking common n','Directional common dual n','Summary min difference','Single-target fraction difference','Concordance','Supporting comparison'],[[r['pair'],SCHEMES[r['method']],*[r[k] for k in ['n_rank_common','n_directional_common_dual','delta_summary_min','delta_top10_single_target_fraction','concordance','supporting_only']]] for r in sortrows(comp)]),
        panel('C','Concordance counts',['Scheme vs S1','Pairs','Concordant improvement','Concordant worsening','Discordant','No change'],[[SCHEMES[r['method']],*[r[k] for k in ['n_pairs','CONCORDANT_IMPROVE','CONCORDANT_WORSEN','DISCORDANT','NO_CHANGE']]] for r in sorted(concord,key=lambda r:METHODS.index(r['method']))])],
        'Panel A uses each method\'s own members and is descriptive. Table 2 and panel B use separate matched populations for every comparison. Joint rank = max(rank A, rank B), with average single-target ranks; ties are resolved by rank sum then global molecular ID. Both differences are new scheme minus S1. Opposite signs define concordant improvement or worsening; equal nonzero signs are discordant; either exact zero defines no change. S4-S1 is supporting, using verified identical directional members and no additional CI. Complete paired composition, per-ligand ranks and all CAP5 ranking results are supplied as data files.')

    check('main table numbering',[t['number'] for t in tables[:2]]==['1','2'])
    check('supplement table numbering',[t['number'] for t in tables[2:]]==['S'+str(i) for i in range(1,11)])
    check('S3 and S4 scoring definitions',SCHEMES['M2']=='S3' and SCHEMES['M1b']=='S4')
    check('concordance totals',all(sum(int(r[k]) for k in ['CONCORDANT_IMPROVE','CONCORDANT_WORSEN','DISCORDANT','NO_CHANGE'])==8 for r in concord))
    # Domain-specific verification of the displayed CAP5 interval boundary.
    boundary=next(r for r in capcompare if r['layer']=='PRIMARY' and r['pair']=='PPARA/PPARD' and r['method']=='M3' and r['contrast']=='dual_vs_A_only')
    check('CAP5 CI boundary retained',float(boundary['original_ci_lo'])<0.5<float(boundary['cap5_ci_lo']))

    index=[]; panel_files={}
    captions=[]
    for t in tables:
        stem='Table_'+t['number']; folder='01_Main_Tables' if not t['number'].startswith('S') else '02_SI_Tables'
        text=[f"# Table {t['number']}. {t['title']}",'',t['title_cn'],'']
        for p in t['panels']:
            p_stem=stem+('_'+p['label'] if p['label'] else '')
            data=csv_bytes(p['headers'],p['rows']); panel_files['Numbered_tables/'+p_stem+'.csv']=data
            if folder=='01_Main_Tables': (out/folder/(p_stem+'.csv')).write_bytes(data)
            if p['label']: text.extend([f"## {p['label']}. {p['title']}",''])
            text.extend(['| '+' | '.join(p['headers'])+' |','| '+' | '.join(['---']*len(p['headers']))+' |'])
            text.extend('| '+' | '.join(display(v) for v in row)+' |' for row in p['rows'])
            text.extend(['',p['note'],''] if p['note'] else [''])
            index.append([t['number'],p['label'],t['title'],t['first_citation'],folder+'/'+('Main_Tables.xlsx' if folder=='01_Main_Tables' else 'Supplementary_Tables.xlsx'),p_stem,len(p['rows']),'Supplementary_Data.zip/Numbered_tables/'+p_stem+'.csv'])
        text.extend(['',t['note'],''])
        (out/folder/(stem+'.md')).write_text('\n'.join(text),encoding='utf-8')
        captions.extend([f"## Table {t['number']}. {t['title']}",'',t['title_cn'],'',t['note'],''])
    (out/'TABLE_CAPTIONS.md').write_text('\n'.join(captions),encoding='utf-8')
    (out/'TABLE_INDEX.csv').write_bytes(csv_bytes(['table','panel','title','first_citation','workbook','data_table','rows','csv_in_data_package'],index))
    source_map={
        '1':['results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv'],
        '2':['results/jcim_stage2/module_I_descriptive.csv'],
        'S1A':[PR+'SEVEN_PAIR_SELECTION_PROVENANCE.csv',Q+'PIK3CA_MTOR_Q110_MEMBERSHIP.csv'],
        'S1B':[PR+'CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv',PR+'CURRENT_867_SOURCE_COVERAGE_SUMMARY.csv'],
        'S2':[RP+'00_protocol/receptor_registry.csv',RP+'13_qa/receptor_preparation_manifest.csv',RP+'13_qa/alt_receptor_preparation_manifest.csv',RP+'13_qa/alt_box_manifest.csv',RP+'03_boxes/'],
        'S3':[RP+'13_qa/redocking_master.csv',RP+'00_protocol/receptor_registry.csv'],
        'S4A':['analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv'],
        'S4B':['results/jcim_stage2/module_H_missingness.csv'],
        'S5A':['results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv'],
        'S5B':['results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv'],
        'S6A':['results/jcim_stage2/module_B_negative_class.csv'],
        'S6B':['results/jcim_stage2/module_C_wrong_pocket.csv'],
        'S7':['results/jcim_stage2/module_D_chemistry_units.csv'],
        'S8A':['results/jcim_stage2/module_E_five_seed_summary.csv'],
        'S8B':['results/jcim_stage2/module_F_alt_receptor.csv'],
        'S9A':[CAP+'CAP5_CLASS_COUNTS.csv'],
        'S9B':['results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv',CAP+'CAP5_PRIMARY_DIRECTIONAL_METRICS.csv'],
        'S9C':[Q+'Q110_STRICT_SUBSET_METRICS.csv'],
        'S9D':['results/jcim_stage2/module_D_chemistry_units.csv',CAP+'CAP5_MODULE_D_CHEMISTRY_UNITS.csv'],
        'S10A':['results/jcim_stage2/module_I_descriptive.csv'],
        'S10B':['results/jcim_stage2/module_I_method_compare.csv'],
        'S10C':['results/jcim_stage2_figures/module_I_concordance_counts.csv'],
    }
    source_register=csv_bytes(['table_panel','source_path','source_commit'],[[key,'Dual_Target_Docking/'+p,BASE_COMMIT] for key,paths in source_map.items() for p in paths])
    (out/'TABLE_SOURCES.csv').write_bytes(source_register)
    scoring=[['Publication scheme','Source method','Definition']]+[[SCHEMES[m],m,n] for m,n in zip(METHODS,NAMES)]
    (out/'SCORING_SCHEMES.csv').write_bytes(csv_bytes(scoring[0],scoring[1:]))
    # Current scores are copied from the accepted population, with published names.
    score_headers=['pair','panel_id','global_ligand_entity_id','class','pdb_A','pdb_B']
    for m in METHODS:
        for side in ['A','B']: score_headers += [SCHEMES[m]+'_valid_'+side,SCHEMES[m]+'_score_'+side,SCHEMES[m]+'_missing_'+side]
    score_rows=[]
    for r in eligible:
        vals=[r[k] for k in score_headers[:6]]
        for m in METHODS:
            for side in ['A','B']: vals += [r[m+'_'+side+'_valid'],r[m+'_score_'+side],r[m+'_'+side+'_missing']]
        score_rows.append(vals)
    raw_extra['Current_scores/CURRENT_SCORES_S1_S5.csv']=csv_bytes(score_headers,score_rows)
    raw_extra['Current_scores/SCORING_SCHEMES.csv']=csv_bytes(scoring[0],scoring[1:])
    # Existing accepted source package remains a source archive. Add missing
    # receptor/redocking evidence and table views without rewriting original bytes.
    oldzip=root/'manuscript/Supplementary_Source_Data.zip'
    with zipfile.ZipFile(oldzip) as z: archive={n:z.read(n) for n in z.namelist() if n!='SHA256SUMS.csv'}
    extra_paths=sorted(sources | {RP+'00_protocol/software_versions.txt',RP+'00_protocol/receptor_component_policy.csv',RP+'00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml',RP+'00_protocol/GNINA_RESCORE_REPRESENTATION_ADDENDUM.yaml',RP+'00_protocol/M3_MACROCYCLE_INPUT_REPAIR_ADDENDUM.yaml','scripts/rerun_v4_2/phase5_cognate_redock.py','results/methods_provenance_recovery_20261011/08_ligand_preparation_report_CN.md'})
    for p in extra_paths:
        if (root/p).exists(): archive['Dual_Target_Docking/'+p]=(root/p).read_bytes()
    archive.update(raw_extra); archive.update(panel_files)
    archive['TABLE_SOURCES.csv']=source_register
    archive['README_CN.md']=(
        '# 补充数据\n\n正式主分析为 867 条靶点对观察，CAP5 为 861 条支持性分析观察。\n'
        'Current_scores/CURRENT_SCORES_S1_S5.csv 是从当前正式人群表复制的完整双端评分；较大数值更有利。\n'
        'S1=M0，S2=M1，S3=M2，S4=M1b，S5=M3。原始来源 CSV 保留 M 编号。\n'
        'Numbered_tables/ 保存正文 2 张表及 SI 10 张表各面板的数值，表中显示精度不改变来源文件精度。\n'
        'Dual_Target_Docking/results/current/chemical_diversity_cap5_20261023/ 含完整 CAP5 评分方法比较、化学 OOF 及联合排序，包括单独的 CAP5_MODULE_I_METHOD_COMPARISON.csv。\n'
        '活性来源目录中的 805 行历史核查附件仅通过当前 867 行 crosswalk 使用其匹配成员，不是另一套当前评价集。\n'
        'receptor_registry.csv 的旧 run_id 不决定结构来源；结构和盒参数对应 V4.2 正式资产。\n'
        '冻结文件中的旧规模或 PM48 条件仅保留为方法来源；当前成员与有效统计以 results/current 和 Q110 当前人群为准。\n'
        '本包不包含完整 ChEMBL 数据库、全部对接姿态或分子结构文件；这些文件的位置见分支中既有路径说明。\n'
    ).encode('utf-8')
    data_index=[]
    for name,content in sorted(archive.items()):
        count=''
        if name.endswith('.csv'):
            try: count=max(0,len(list(csv.reader(io.StringIO(content.decode('utf-8-sig')))))-1)
            except UnicodeDecodeError: pass
        data_index.append([name,len(content),count])
    idx=csv_bytes(['path_in_archive','bytes','data_rows'],data_index)
    (out/'03_Supplementary_Data/Data_Index.csv').write_bytes(idx); archive['Data_Index.csv']=idx
    with zipfile.ZipFile(out/'03_Supplementary_Data/Supplementary_Data.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,content in sorted(archive.items()):
            info=zipfile.ZipInfo(name,date_time=(2026,10,26,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED; z.writestr(info,content)
    old_to_new={1:('Data','pair_ligand_mapping.csv'),2:('S4A','Population layers'),3:('S5A','Directional estimates'),4:('S5B','Paired differences'),5:('S6A','Negative classes'),6:('S6B','Pocket correspondence'),7:('S7','Chemistry OOF summary'),8:('Data','module_D1_composition.csv'),9:('Data','module_D_chemistry_oof_ligands.csv'),10:('Data; S8A summary','module_E_five_seed.csv'),11:('S8A','Seed summaries'),12:('S8B','Six current alternative conditions'),13:('S4B; Data class details','module_H_missingness.csv'),14:('S10A; Table 2; Data','Own-method and matched ranking populations'),15:('S10B','Directional and top-list change'),16:('Data','module_I_ligand_ranks.csv'),17:('S10C','Concordance counts'),18:('S9C','Strict Q110 subset'),19:('S1B; Data','Source coverage'),20:('Data; S1B summary','CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv'),21:('Data','Q110_CHEMBL37_ACTIVITY_ROWS_LONG.csv'),22:('S1A; Data','Selection evidence plus Q110 membership'),23:('S9A; Data class details','CAP5_CLASS_COUNTS.csv'),24:('Data','CAP5_EXCLUDED_LIGANDS.csv'),25:('S9B; Data complete methods','CAP5_PRIMARY_DIRECTIONAL_METRICS.csv'),26:('Data','CAP5_PRIMARY_METHOD_DELTAS.csv'),27:('Data','CAP5_MODULE_B_NEGATIVE_CLASS.csv'),28:('Data','CAP5_MODULE_C_POCKET_CORRESPONDENCE.csv'),29:('S9D; Data','CAP5_MODULE_D_CHEMISTRY_UNITS.csv'),30:('Data','CAP5_MODULE_I_TOPK.csv'),31:('S9B summary; Data complete','CAP5_PRIMARY_VS_ORIGINAL.csv')}
    (out/'RENUMBERING.csv').write_bytes(csv_bytes(['previous_workbook_table','submission_destination','content'],[['S'+str(i),*old_to_new[i]] for i in range(1,32)]))
    with tempfile.TemporaryDirectory(prefix='submission_tables_') as tmp:
        runtime_modules=os.environ.get('CODEX_PRIMARY_RUNTIME_NODE_MODULES')
        if runtime_modules: os.symlink(runtime_modules,Path(tmp)/'node_modules',target_is_directory=True)
        payload=Path(tmp)/'layout.json'; payload.write_text(json.dumps({'output':str(out),'preview_dir':str(Path(tmp)/'previews'),'tables':tables}),encoding='utf-8')
        node=os.environ.get('CODEX_PRIMARY_RUNTIME_NODE','node')
        subprocess.run([node,str(Path(__file__).with_name('build_submission_workbooks.mjs')),str(payload)],check=True,cwd=tmp)
        for p in out.rglob('*.xlsx.inspect.ndjson'): p.unlink()
        # Kept beside the external validation output, outside the submission data.
        preview_dest=root.parent/'submission_previews'; preview_dest.mkdir(exist_ok=True)
        for p in (Path(tmp)/'previews').glob('*.png'): (preview_dest/p.name).write_bytes(p.read_bytes())
    # Validate every exported panel against the matrix used to write its CSV.
    import openpyxl  # read-only verification; all workbook authoring is JS.
    for kind in ['main','si']:
        chosen=tables[:2] if kind=='main' else tables[2:]
        path=out/('01_Main_Tables/Main_Tables.xlsx' if kind=='main' else '02_SI_Tables/Supplementary_Tables.xlsx')
        wb=openpyxl.load_workbook(path,read_only=True,data_only=False)
        check(kind+' workbook sheets',wb.sheetnames==['Table_'+t['number'] for t in chosen])
        for t in chosen:
            ws=wb['Table_'+t['number']]; rr=list(ws.values); line=4
            for p in t['panels']:
                if p['label']: line+=1
                check('header '+t['number']+p['label'],list(rr[line-1][:len(p['headers'])])==p['headers'])
                for i,want in enumerate(p['rows']):
                    actual=list(rr[line+i][:len(want)])
                    valid=all((a==b) if not isinstance(b,float) else (isinstance(a,(int,float)) and math.isclose(a,b,rel_tol=1e-13,abs_tol=1e-13)) for a,b in zip(actual,want))
                    check('cells '+t['number']+p['label']+' row '+str(i+1),valid)
                line+=len(p['rows'])+3
            check('no Excel error '+t['number'],not any(c.data_type=='e' for row in ws for c in row))
        wb.close()
    with zipfile.ZipFile(out/'03_Supplementary_Data/Supplementary_Data.zip') as z:
        check('data archive integrity',z.testzip() is None)
        check('CAP5 paired ranking evidence present','Dual_Target_Docking/'+CAP+'CAP5_MODULE_I_METHOD_COMPARISON.csv' in z.namelist())
        check('current scores exported',len(list(csv.DictReader(io.StringIO(z.read('Current_scores/CURRENT_SCORES_S1_S5.csv').decode()))))==867)
        for p in sources:
            check('source bytes '+p,z.read('Dual_Target_Docking/'+p)==(root/p).read_bytes())
    report={'source_commit':BASE_COMMIT,'presentation_only':True,'scientific_recomputation':False,'main_tables':2,'supplementary_tables':10,'formal_observations':867,'cap5_observations':861,'data_files':len(archive),'passed':len(checks),'failed':0,'checks':checks}
    (out/'04_Editorial_Checks/Validation_Report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'}))

if __name__=='__main__': main()
