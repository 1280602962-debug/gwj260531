"""Verify current CSV/workbook/package consistency against immutable source rows."""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import subprocess
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parent
PAIR='PIK3CA/mTOR'
BASELINE='dc028b3350e804cbf8f40acea08261ed4f89c4fa'
Q=ROOT/'extension_jcim_validation/pik3ca_mtor_q110_v42'


def read(path):
    with path.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-ref',default=BASELINE)
    args=parser.parse_args()
    checks=[]
    def check(name,condition,detail=''):
        checks.append({'check':name,'status':'PASS' if condition else 'FAIL','detail':detail})
    current=ROOT/'results/current'
    pop=read(current/'PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv')
    eligible=[r for r in pop if r['activity_eligible']=='1' and r['class']]
    check('population',len(pop)==869 and len(eligible)==867)
    check('classes',collections.Counter(r['class'] for r in eligible)=={'dual':248,'A_only':256,'B_only':252,'neither':111})
    qp=[r for r in pop if r['pair']==PAIR]
    check('Q110_population_exact',qp==read(Q/'Q110_POPULATION.csv'))
    check('Q110_finite_scores',all(r[f'{m}_{s}_valid']=='1' and math.isfinite(float(r[f'{m}_score_{s}'])) for r in qp for m in ['M0','M1','M1b','M2','M3'] for s in ['A','B']))
    for path in [*sorted((ROOT/'results/formal_metrics').glob('*.csv')),*sorted((ROOT/'results/jcim_stage2').glob('*.csv'))]:
        rel=path.relative_to(REPO).as_posix()
        old=list(csv.DictReader(io.StringIO(subprocess.check_output(['git','show',f'{args.baseline_ref}:{rel}'],cwd=REPO,text=True))))
        new=read(path)
        check('seven_pairs_unchanged:'+path.name,[r for r in old if r['pair']!=PAIR]==[r for r in new if r['pair']!=PAIR])
    for filename,qfile,keys in [('PRIMARY_DIRECTIONAL_METRICS.csv','Q110_PRIMARY_METRICS.csv',['method','contrast']),('PRIMARY_METHOD_DELTA_METRICS.csv','Q110_PRIMARY_DELTAS.csv',['delta'])]:
        now={tuple(r[k] for k in keys):r for r in read(ROOT/'results/formal_metrics'/filename) if r['pair']==PAIR}
        qrows=read(Q/qfile)
        check('Q110_fields_exact:'+filename,all(all(now[tuple(qr[k] for k in keys)][k]==v for k,v in qr.items() if k in now[tuple(qr[k] for k in keys)]) for qr in qrows))
    metrics=read(ROOT/'results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv')
    # Independent pairwise Mann-Whitney evaluation, solely as verification.
    for row in [r for r in metrics if r['pair']==PAIR and r['contrast']!='summary_min']:
        side,negative=('B','A_only') if row['contrast']=='dual_vs_A_only' else ('A','B_only')
        pos=[float(r[f'{row["method"]}_score_{side}']) for r in qp if r['class']=='dual']
        neg=[float(r[f'{row["method"]}_score_{side}']) for r in qp if r['class']==negative]
        value=sum((x>y)+.5*(x==y) for x in pos for y in neg)/(len(pos)*len(neg))
        check('Q110_AUROC:'+row['method']+':'+side,abs(value-float(row['auroc']))<1e-9)
    layers=read(ROOT/'analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv')
    check('population_layers',[sum(int(r[k]) for r in layers) for k in ['n_labeled_population','n_m0_directional_population','n_final_chemistry_oof']]==[1004,1000,994])
    alt=read(ROOT/'results/jcim_stage2/module_F_alt_receptor.csv')
    check('Q110_no_alt_rows',len(alt)==6 and all(r['pair']!=PAIR for r in alt))
    comparison=read(ROOT/'results/jcim_stage2/module_I_method_compare.csv')
    counts=read(ROOT/'results/jcim_stage2_figures/module_I_concordance_counts.csv')
    states=['CONCORDANT_IMPROVE','CONCORDANT_WORSEN','DISCORDANT','NO_CHANGE']
    check('concordance_counts',all(int(r[s])==sum(x['method']==r['method'] and x['concordance']==s for x in comparison) for r in counts for s in states))
    check('rank_count',len(read(ROOT/'results/jcim_stage2/module_I_ligand_ranks.csv'))==9734)
    # Read workbooks independently; no workbook authoring occurs here.
    from openpyxl import load_workbook
    main=load_workbook(ROOT/'manuscript/Main_Tables.xlsx',read_only=True,data_only=True)
    table=list(main['Table1_Label_Counts'].values)
    pm=next(r for r in table[1:] if r[1]==PAIR)
    check('Table1_Q110',tuple(map(int,pm[6:10]))==(32,32,32,14))
    sup=load_workbook(ROOT/'manuscript/Supplementary_Tables.xlsx',read_only=True,data_only=True)
    tabs={'S3_PRIMARY_Metrics':ROOT/'results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv','S4_PRIMARY_Deltas':ROOT/'results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv','S12_Alt_Receptors':ROOT/'results/jcim_stage2/module_F_alt_receptor.csv','S16_Ligand_Ranks':ROOT/'results/jcim_stage2/module_I_ligand_ranks.csv','S17_Concordance_Counts':ROOT/'results/jcim_stage2_figures/module_I_concordance_counts.csv','S18_Q110_Strict_Subset':Q/'Q110_STRICT_SUBSET_METRICS.csv'}
    for name,path in tabs.items():
        source=read(path);fields=list(source[0]);expected=[fields]+[[r.get(k,'') for k in fields] for r in source]
        actual=[[str(v) if v is not None else '' for v in row[:len(fields)]] for row in sup[name].values if any(v is not None and v != '' for v in row)]
        check('workbook_matches_csv:'+name,actual==expected)
    main.close();sup.close()
    with zipfile.ZipFile(ROOT/'manuscript/Supplementary_Source_Data.zip') as z:
        manifest=list(csv.DictReader(io.StringIO(z.read('SHA256SUMS.csv').decode())))
        check('package_sha256',all(hashlib.sha256(z.read(r['relative_path'])).hexdigest()==r['sha256'] for r in manifest))
        check('package_current_files',all(z.read(r['relative_path'])==(REPO/r['relative_path']).read_bytes() for r in manifest))
        check('package_no_PM48_regression_outputs',not any('regression_pm48' in p for p in z.namelist()))
    for path in ['Results_Review.pdf','Supporting_Information_Review.pdf']:
        import fitz
        doc=fitz.open(ROOT/'manuscript'/path)
        check('PDF_nonempty:'+path,len(doc)>0 and all(page.get_text().strip() for page in doc));doc.close()
    report={'source_commit':BASELINE,'n_checks':len(checks),'n_pass':sum(r['status']=='PASS' for r in checks),'n_fail':sum(r['status']=='FAIL' for r in checks),'checks':checks}
    (current/'CURRENT_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'}))
    for r in checks:
        if r['status']=='FAIL':print(r)
    return bool(report['n_fail'])


if __name__=='__main__':raise SystemExit(main())
