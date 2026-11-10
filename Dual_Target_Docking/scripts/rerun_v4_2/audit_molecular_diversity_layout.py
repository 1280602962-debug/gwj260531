"""Verify exported PDF geometry and assemble a separate captioned review proof."""
from pathlib import Path
import csv
import json
import re
import fitz
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'submission/Molecular_Diversity'
MM=72/25.4


def spans(page):
    return [s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text'].strip()]


def main():
    index=list(csv.DictReader((OUT/'Figure_Index.csv').open()))
    report={'figure_margins':[],'caption_gaps_mm':[],'figure_caption_count':0,'figure_caption_end_periods':True}
    physical=json.loads((OUT/'Checks/Publication_Layout_QA.json').read_text())
    physical_by={r['figure']:r for r in physical['figures']}
    contrast=[]
    def luminance(rgb):
        return sum(w*(v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4) for w,v in zip([.2126,.7152,.0722],rgb))
    for r in index:
        source=OUT/(r['file']+'.pdf');doc=fitz.open(source);p=doc[0];ss=spans(p)
        margins=[min(s['bbox'][0] for s in ss),min(s['bbox'][1] for s in ss),p.rect.width-max(s['bbox'][2] for s in ss),p.rect.height-max(s['bbox'][3] for s in ss)]
        margins=[v/MM for v in margins];assert min(margins)>=4,(source,margins)
        report['figure_margins'].append({'figure':r['figure'],'left_top_right_bottom_mm':margins})
        labels=physical_by[r['figure']]['heatmap_labels']
        if labels:
            pix=p.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
            for item in labels:
                x,y=item['background_sample_pdf_pt'];bg=pix.pixel(round(x*2),round(y*2))[:3]
                cx,cy=item['text_center_pdf_pt'];choices=[s for s in ss if s['text']==item['text']]
                assert choices,(r['figure'],item['text'])
                text=min(choices,key=lambda s:abs((s['bbox'][0]+s['bbox'][2])/2-cx)+abs((s['bbox'][1]+s['bbox'][3])/2-cy))
                n=text['color'];fg=[(n>>16)&255,(n>>8)&255,n&255]
                a=luminance([v/255 for v in bg]);b=luminance([v/255 for v in fg]);ratio=(max(a,b)+.05)/(min(a,b)+.05)
                assert ratio>=4.5,(r['figure'],item['text'],ratio,bg,fg)
                contrast.append({'figure':r['figure'],'text':item['text'],'ratio':ratio,'background_rgb':bg,'foreground_rgb':fg})
    review_text=(OUT/'Figure_Captions_EN.md').read_text(encoding='utf8')
    captions=[p.strip() for p in review_text.split('\n\n') if p.strip()]
    assert len(captions)==len(index)==21
    externalized=['6','8','S2','S4','S5']
    for row,caption in zip(index,captions):
        if row['figure'] in externalized:
            assert 'Target-pair codes: 1, EGFR/HER2;' in caption
            assert '8, PPARA/PPARD' in caption
            page=fitz.open(OUT/(row['file']+'.pdf'))[0]
            assert '1 EGFR/HER2' not in page.get_text()
            assert '5 AChE/BChE' not in page.get_text()
    report['target_pair_keys_in_external_captions']=externalized
    assert all(p.endswith('.') for p in captions)
    submission=(OUT/'Figure_Captions_EN_Submission.md').read_text(encoding='utf8')
    assert all(not p.strip().endswith('.') for p in submission.split('\n\n') if p.strip())
    assert '(b) minimum AUROC' not in review_text
    report['figure_caption_count']=len(captions)
    proof=fitz.open();font_path=font_manager.findfont('Arial')
    for r,caption in zip(index,captions):
        src=fitz.open(OUT/(r['file']+'.pdf'));page=proof.new_page(width=210*MM,height=297*MM)
        left=(210*MM-src[0].rect.width)/2;top=12*MM
        box=fitz.Rect(left,top,left+src[0].rect.width,top+src[0].rect.height)
        page.show_pdf_page(box,src,0)
        page.insert_font(fontname='ArialReview',fontfile=font_path)
        plain=caption.replace('**','')
        capbox=fitz.Rect(left,box.y1+5*MM,page.rect.width-left,page.rect.height-12*MM)
        free=page.insert_textbox(capbox,plain,fontname='ArialReview',fontsize=9,lineheight=1.2)
        assert free>=0,(r['figure'],'caption overflow')
        first=min(s['bbox'][1] for s in spans(page) if s['bbox'][1]>box.y1)
        gap=(first-box.y1)/MM
        assert 4<=gap<=6,(r['figure'],gap)
        report['caption_gaps_mm'].append({'figure':r['figure'],'gap':gap})
    proof.set_toc([[1,'Scheme S1' if r['figure']=='Scheme_S1' else 'Fig. '+r['figure'],i+1] for i,r in enumerate(index)])
    proof.save(OUT/'Figures_With_Captions_Review.pdf');proof.close()
    # Check footnotes and the two dense tables without changing any data rows.
    table_qa=json.loads((OUT/'Checks/Tables_QA.json').read_text())
    assert table_qa['review_pages']==29,table_qa['review_pages']
    for p in (OUT/'Tables').glob('Table*.md'):
        text=p.read_text(encoding='utf8').rstrip();assert text.endswith('.'),p
    dense=[r for r in table_qa['pages'] if r['table']=='Table 2' or (r['table']=='Table S5' and r['panel'].startswith('A.'))]
    assert next(r for r in dense if r['table']=='Table 2')['rows']==32
    assert all(r['height_mm']<=213 for r in dense),dense
    report['dense_table_pages']=dense
    report['heatmap_exported_pdf_contrast_checks']=len(contrast)
    report['minimum_heatmap_pdf_contrast']=min(r['ratio'] for r in contrast)
    report['heatmap_contrast_samples']=contrast
    report['minimum_pdf_margin_mm']=min(min(r['left_top_right_bottom_mm']) for r in report['figure_margins'])
    report['source_data_recomputed']=False
    (OUT/'Checks/Exported_Layout_QA.json').write_text(json.dumps(report,indent=2)+'\n')
    # Render representative proof pages for inspection; they are not submission assets.
    doc=fitz.open(OUT/'Figures_With_Captions_Review.pdf')
    for i in [1,4,5,7,12,20]:doc[i].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(OUT/'Checks'/f'CaptionProof{i+1:02d}.png')
    print(json.dumps({'pdf_figures_checked':len(index),'minimum_pdf_margin_mm':report['minimum_pdf_margin_mm'],
                      'captioned_pages':len(captions),'caption_gap_range_mm':[min(r['gap'] for r in report['caption_gaps_mm']),max(r['gap'] for r in report['caption_gaps_mm'])],
                      'table_pages':table_qa['review_pages']}))


if __name__=='__main__':main()
