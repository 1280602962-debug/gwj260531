"""Independent table/manuscript/figure checks and contact sheets for visual QA."""
from pathlib import Path
import csv
import hashlib
import json
import re
from PIL import Image, ImageOps, ImageDraw

OUT=Path(__file__).resolve().parents[1]
audit=json.loads((OUT/'plotted_values.json').read_text(encoding='utf-8'))
checks=[]
def check(ok,what):
    checks.append((bool(ok),what))

def content(rel):
    return (OUT/audit['input_files'][rel]).read_text(encoding='utf-8-sig')

def rows(rel):return list(csv.DictReader(content(rel).splitlines()))
def nums(s):return [float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',s)]
def close(a,b):return len(a)==len(b) and all(abs(x-y)<=.00051 for x,y in zip(a,b))

for rel,digest in audit['inputs_sha256'].items():
    check(hashlib.sha256((OUT/audit['input_files'][rel]).read_bytes()).hexdigest()==digest,'SHA256 '+rel)

primary=audit['plotted']['primary']
order=['EGFR/HER2','JAK1/JAK2','JAK1/TYK2','PIK3CA/mTOR','AChE/BChE','F2/F10','PPARG/PPARA','PPARA/PPARD']
check(audit['plotted'].get('pair_order')==order,'scientific display order on figures')
check(list(primary)==order,'primary plotted dict follows scientific order')
fa=float(audit['plotted']['fig2']['EGFR/HER2']['fixed_A']['delta_neither_minus_selective'])
fb=float(audit['plotted']['fig2']['JAK1/TYK2']['fixed_A']['delta_neither_minus_selective'])
check(abs(fa-0.378)<0.001,'Figure 2A EGFR/HER2 fixed-score Δ is 0.378')
check(abs(fb-0.444)<0.001,'Figure 2A JAK1/TYK2 fixed-score Δ is 0.444')
check(audit.get('data_snapshot_commit')=='abb61a20a04eb6a085ad526876624eadb518c4cc','data snapshot commit is pinned separately from artwork')
check(bool(audit.get('artwork_git_head')),'artwork git HEAD recorded at generation')
captions=(OUT/'MANUSCRIPT_FIGURE_CAPTIONS.md').read_text(encoding='utf-8')
for phrase in ['initially evaluated','added after the census','original three','census five','later five','Horizontal gray rules','underpowered','flagship','data-collection sequence','were then selected','python3','PowerPoint','historical original-set','Regenerate']:
    check(phrase.lower() not in captions.lower(),'caption has no '+phrase)
check('The census summarizes the availability' in captions,'Fig1C caption separates census from the eight-pair evaluation')
check('## Figure 4. Matched- versus mismatched-pocket scoring controls.' in captions,'Figure 4 caption is pocket correspondence')
check('## Figure 5. Computational realization.' in captions,'Figure 5 caption is computational realization')
check('## Figure S1. Post-hoc formulation and screening diagnostics.' in captions,'SI figures start at Top-10/AND diagnostics')
check('Primary rank-based Vina AUROC' in captions,'Fig3A caption distinguishes Vina rank from ECFP4 GroupKFold')
check('FigS3_protocol_sensitivity' in audit['generated'] and 'FigS4_cognate_rmsd' in audit['generated'],'SI protocol and cognate RMSD figures generated')
check('Fig4_mismatched_pocket' in audit['generated'] and 'Fig5_computational_realization' in audit['generated'],'main Figures 4/5 are pocket then computational')
check('FigS1_posthoc_diagnostics' in audit['generated'] and 'FigS2_pocket_matched_forest' in audit['generated'],'SI figures start with Top-10 then descriptor forest')
for lang in ['ZH','EN']:
    manuscript=content('docs/MANUSCRIPT_JCIM_'+lang+'.md')
    for pair,r in primary.items():
        lines=[line for line in manuscript.splitlines() if line.startswith('| '+pair+' |')]
        t2=next(line for line in lines if len(line.split('|'))==7 and re.fullmatch(r'\s*\d+ / \d+ / \d+\s*',line.split('|')[2]) and re.fullmatch(r'\s*0?\.\d+\s*',line.split('|')[3]))
        cells=[x.strip() for x in t2.split('|')[1:-1]]
        # Table 2 cells: n_scored, D/A, D/B, summary_min [lo, hi].
        observed=[float(cells[2]),float(cells[3])]+nums(cells[4])
        check(close(observed,[r['da'],r['db'],r['smin'],r['lo'],r['hi']]),lang+' Table 2 plotted AUROCs/CI '+pair)
        t3=next(line for line in lines if len(line.split('|'))==6 and '[' in line)
        cells=[x.strip() for x in t3.split('|')[1:-1]]
        # Table 3 cells: directional summary [lo, hi], neither, n_neither.
        observed=nums(cells[1])+nums(cells[2])+[float(cells[3])]
        check(close(observed,[r['smin'],r['lo'],r['hi'],r['nei'],r['nei_lo'],r['nei_hi'],r['n_neg']]),lang+' Table 3 plotted AUROCs/CI/n_neither '+pair)

external=rows('data/jcim_novelty_v0/tables/external_slice_summary_v1.csv')
check(len(external)==8 and all(r['packaged_as_external_evaluation']=='0' for r in external),'BindingDB: eight pairs and zero external evaluations')
check(all(min(int(r[k]) for k in ['n_dual','n_A_only','n_B_only'])<20 or min(int(r[k]) for k in ['n_sources_dual','n_sources_A_only','n_sources_B_only'])<3 for r in external),'BindingDB failures follow class-count/source-count criteria')
check('PIK3CA/PIK3CB' not in primary,'Withdrawn PIK3CB absent from primary figures')
check(abs(audit['plotted']['fig3B_max_abs']-.0231)<.0006,'ECFP4 incremental maximum agrees with manuscript 0.023')
sim=audit['plotted']['figS6']
check(len({r['pair'] for r in sim})==3 and {float(r['true_auroc']) for r in sim}=={.50,.55,.60,.65,.70,.75},'Simulation: three pairs, complete six-point grid')
off=audit['plotted'].get('figS12_offscale') or []
check(any(abs(float(r['rmsd'])-9.505)<.002 and r.get('series')=='top1' for r in off),'EGFR 3POZ top-1 is plotted off-scale')
check(all('top3' not in r for r in audit['plotted']['figS12']),'Figure S4 plots only top-1 and best-of-9')
cl=audit['plotted']['figS11']
jdoc=next(r for r in cl if r['pair']=='JAK1/TYK2' and r['estimator']=='document_cluster')
check(float(jdoc['delta_ci_lo'])<0<float(jdoc['delta_ci_hi']),'JAK1/TYK2 document-cluster interval crosses zero')
top=audit['plotted']['figS7']['top10'];filt=audit['plotted']['figS7']['and_filter']
check([int(top[k]) for k in ['n_dual_top','n_A_only_top','n_B_only_top','n_neither_top']]==[1,5,4,0],'Top-10 class counts: 1/5/4/0')
check([int(filt[k]) for k in ['n_dual_pass','n_A_only_pass','n_B_only_pass']]==[14,9,24],'AND-filter class counts: 14/9/24')
files={}
for stem in audit['generated']:
    for ext in (['png','tif'] if stem=='TOC_graphic' else ['png','tif','pdf']):
        path=OUT/(stem+'.'+ext);files[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        if ext!='pdf':
            with Image.open(path) as im:check(im.mode=='RGB' and abs(im.info['dpi'][0]-300)<1 and im.width<=2101 and im.height<=2751,'RGB/300dpi/size '+path.name)
(OUT/'FIGURE_SHA256.json').write_text(json.dumps(files,indent=2),encoding='utf-8')
qa=OUT/'qa';qa.mkdir(exist_ok=True)
for group,stems in [('main',[s for s in audit['generated'] if re.match(r'Fig[1-6]_',s)]),('supplement',[s for s in audit['generated'] if s.startswith('FigS')]+['TOC_graphic'])]:
    canvas=Image.new('RGB',(1400,550*((len(stems)+1)//2)),'#eeeeee');draw=ImageDraw.Draw(canvas)
    for i,s in enumerate(stems):
        with Image.open(OUT/(s+'.png')) as im:thumb=ImageOps.contain(im,(680,510))
        xx=(i%2)*700;yy=(i//2)*550;canvas.paste(thumb,(xx+(700-thumb.width)//2,yy+25));draw.text((xx+12,yy+5),s,fill='black')
    canvas.save(qa/(group+'_contact_sheet.jpg'),quality=92)
report=['# PR32 figure numerical audit','',
        f"Pinned numerical-data snapshot: `{audit.get('data_snapshot_commit', audit.get('commit'))}`",
        f"Artwork generated from git HEAD: `{audit.get('artwork_git_head','')}`",
        '',f"{sum(x for x,_ in checks)} PASS / {sum(not x for x,_ in checks)} FAIL",'', '| Status | Check |','|---|---|']
report += [f'| {"PASS" if ok else "FAIL"} | {what} |' for ok,what in checks]
(OUT/'FIGURE_AUDIT_PR32.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
print(next(line for line in report if 'PASS /' in line))
if any(not ok for ok,_ in checks):raise SystemExit('Figure audit failed')
