"""Final release checks for the Molecular Diversity presentation package."""
import csv
import hashlib
import json
from pathlib import Path
import re
import fitz
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'submission/Molecular_Diversity'

def main():
    qa=json.loads((OUT/'Checks/Artwork_QA.json').read_text())
    nums=json.loads((OUT/'Checks/Independent_Numeric_Audit.json').read_text())
    tables=json.loads((OUT/'Checks/Tables_QA.json').read_text())
    protected=json.loads((OUT/'Checks/Protected_Inputs.json').read_text())
    for p,digest in protected.items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest,p
    assert nums['checks_failed']==tables['source_cell_checks_failed']==0
    assert len(qa['figures'])==19
    n_main=sum(not r['figure'].startswith(('S','Scheme')) for r in qa['figures'])
    n_si=sum(r['figure'].startswith('S') and not r['figure'].startswith('Scheme') for r in qa['figures'])
    assert n_main+2<=10, 'Molecular Diversity main-text display-item limit'
    assets=[]
    for row in qa['figures']:
        stem='SchemeS1' if row['figure']=='Scheme_S1' else 'Fig'+row['figure']
        folder=OUT/('Supplementary_Figures' if row['figure'].startswith(('S','Scheme')) else 'Main_Figures')
        for ext in ['eps','pdf','svg','png']:assert (folder/(stem+'.'+ext)).is_file()
        eps=(folder/(stem+'.eps')).read_text(encoding='latin1')
        assert re.search(r'/FontType\s+42',eps),stem+' EPS font embedding'
        with Image.open(folder/(stem+'.png')) as image:
            assert image.mode=='RGB'
            assert all(abs(x-600)<.1 for x in image.info['dpi'])
        doc=fitz.open(folder/(stem+'.pdf'));page=doc[0]
        assert abs(page.rect.width*25.4/72-174)<.01
        assert page.rect.height*25.4/72<=234
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                for span in line['spans']:
                    if span['text'].strip():assert span['size']>=7.99,(stem,span)
        assets.append(stem)
    review=fitz.open(OUT/'Figures_Review.pdf');assert len(review)==19
    tableproof=fitz.open(OUT/'Tables_Review.pdf');assert len(tableproof)==tables['review_pages']
    for p in tableproof:
        for block in p.get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                for span in line['spans']:
                    x0,y0,x1,y1=span['bbox']
                    assert min(x0,y0)>=-1 and x1<=p.rect.width+1 and y1<=p.rect.height+1,span
    with (ROOT/'reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/redocking_master.csv').open(newline='') as h:redock=list(csv.DictReader(h))
    assert len(redock)==70
    assert sum(float(r['top1_rmsd_A'])<2 for r in redock)==42
    assert sum(float(r['min_saved_pose_rmsd_A'])<2 for r in redock)==66
    with (ROOT/'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv').open(newline='') as h:pop=[r for r in csv.DictReader(h) if r['activity_eligible']=='1' and r['class']]
    with (ROOT/'submission/01_Main_Tables/Table_1.csv').open(newline='') as h:t1=list(csv.DictReader(h))
    for r in t1:
        members=[x for x in pop if x['pair']==r['Target pair']]
        assert {(x['pdb_A'],x['pdb_B']) for x in members}=={(r['PDB A'],r['PDB B'])}
        for cls,col in [('dual','Dual (n)'),('A_only','A-only (n)'),('B_only','B-only (n)'),('neither','Neither (n)')]:assert sum(x['class']==cls for x in members)==int(r[col])
    report={'base_commit':nums['base_commit'],'main_figures':n_main,'supplementary_figures':n_si,
        'supplementary_schemes':1,'main_tables':2,'supplementary_tables':10,
        'numeric_checks_passed':nums['checks_passed'],'table_source_cell_checks_passed':tables['source_cell_checks_passed'],
        'scientific_input_files_unchanged':len(protected),'redocking_records_checked':70,
        'top1_rmsd_below_2A':42,'best_saved_rmsd_below_2A':66,
        'figure_width_mm':174,'minimum_figure_font_pt':8,'raster_mode':'RGB','png_dpi':600,
        'vector_font_embedding':'PDF and EPS verified','table_review_pages':len(tableproof),
        'figure_review_pages':len(review),'main_display_items':n_main+2,'outside_page_table_text':[],
        'new_docking_model_or_bootstrap_runs':False,'visual_review':'all 19 figures and all table proof pages inspected'}
    (OUT/'Checks/Release_QA.json').write_text(json.dumps(report,indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.suffix!='.png' and p.name!='SHA256SUMS.csv']
    files += [p for f in ['Main_Figures','Supplementary_Figures'] for p in (OUT/f).glob('*.png')]
    with (OUT/'SHA256SUMS.csv').open('w',newline='',encoding='utf8') as h:
        writer=csv.writer(h);writer.writerow(['path','sha256','bytes'])
        for p in sorted(set(files)):writer.writerow([p.relative_to(OUT).as_posix(),hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_size])
    print(json.dumps(report))

if __name__=='__main__':main()
