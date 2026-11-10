"""Render the Table S11 addendum without regenerating accepted Tables S1-S10."""
from pathlib import Path
import csv
import shutil
import sys

import fitz
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

from supplement_single_target_controls import NOTE

root = Path(__file__).resolve().parents[2]
data = root/'analysis/single_target_controls_20261110'
dest = root/'submission/Molecular_Diversity'
rows = list(csv.DictReader((data/'Table_S11_metrics.csv').open()))
pdfmetrics.registerFont(TTFont('EmbeddedSans','/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'))
pdfmetrics.registerFont(TTFont('EmbeddedSansBold','/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf'))
pdfmetrics.registerFontFamily('EmbeddedSans',normal='EmbeddedSans',bold='EmbeddedSansBold')
style = ParagraphStyle('body', fontName='Helvetica', fontSize=9, leading=12, textColor=colors.HexColor('#222222'))
title = ParagraphStyle('title', parent=style, fontName='Helvetica-Bold', fontSize=12, leading=16)
subtitle = ParagraphStyle('subtitle', parent=style, fontName='Helvetica-Bold', fontSize=10, leading=14)
header = ParagraphStyle('header', parent=style, fontName='Helvetica-Bold', fontSize=8.5, leading=11)
style.fontName='EmbeddedSans'; title.fontName='EmbeddedSansBold'
subtitle.fontName='EmbeddedSansBold'; header.fontName='EmbeddedSansBold'
pdf = dest/'Tables/Table_S11.pdf'
story = []
for panel, heading in [('A','Single-target-active versus neither'), ('B','Single-target-active versus opposite single-target-active')]:
    if panel == 'B':
        story.append(PageBreak())
    story += [Paragraph('Table S11. Exploratory single-target controls', title), Spacer(1,8),
              Paragraph(panel+'. '+heading, subtitle), Spacer(1,10)]
    values = [[Paragraph(x, header) for x in ['Target pair','Score target','Positive class','n positive','n negative','AUROC (95% CI)']]]
    for r in rows:
        if r['panel'] != panel:
            continue
        values.append([r['pair'], r['score_target'], r['positive_class'].replace('_','-'), r['n_positive'], r['n_negative'],
                       f"{float(r['auroc']):.4f} [{float(r['ci_lower']):.4f}, {float(r['ci_upper']):.4f}]"])
    table = Table(values, colWidths=[93,62,75,53,53,187], repeatRows=1, hAlign='LEFT')
    table.setStyle(TableStyle([
        ('FONTNAME',(0,0),(-1,-1),'EmbeddedSans'), ('FONTSIZE',(0,1),(-1,-1),9),
        ('VALIGN',(0,0),(-1,-1),'MIDDLE'), ('ALIGN',(3,1),(-1,-1),'RIGHT'),
        ('LEFTPADDING',(0,0),(-1,-1),5), ('RIGHTPADDING',(0,0),(-1,-1),5),
        ('TOPPADDING',(0,0),(-1,-1),6), ('BOTTOMPADDING',(0,0),(-1,-1),6),
        ('LINEABOVE',(0,0),(-1,0),.8,colors.black), ('LINEBELOW',(0,0),(-1,0),.6,colors.black),
        ('LINEBELOW',(0,-1),(-1,-1),.8,colors.black),
        *[('LINEABOVE',(0,i),(-1,i),.3,colors.HexColor('#BBBBBB')) for i in range(3,17,2)]
    ]))
    story += [table, Spacer(1,14), Paragraph(NOTE,style)]

SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=36,rightMargin=36,topMargin=32,bottomMargin=32,
                  title='Table S11. Exploratory single-target controls',author='').build(story)
proof = fitz.open(pdf)
assert len(proof)==2
assert all('PPARA/PPARD' in page.get_text() for page in proof)
proof.close()
shutil.copy2(data/'TableS11.md',dest/'Tables/TableS11.md')
shutil.copy2(data/'Table_S11_metrics.csv',dest/'Source_Data/Tables/Table_S11.csv')

# Idempotent append: retained S1-S10 proof pages are copied, not rerendered.
combined_path = dest/'Tables_Review.pdf'
old = fitz.open(combined_path)
new = fitz.open()
existing = next((i for i,p in enumerate(old) if 'Table S11.' in p.get_text()),len(old))
new.insert_pdf(old,from_page=0,to_page=existing-1)
add = fitz.open(pdf)
new.insert_pdf(add)
tmp = combined_path.with_suffix('.tmp.pdf')
new.save(tmp)
new.close(); old.close(); add.close()
tmp.replace(combined_path)
html_path=dest/'Tables/Editable_Tables.html'
old_html=html_path.read_text()
marker='<section id="table-s11-a">'
if marker in old_html:
    old_html=old_html.split(marker)[0]+'</html>'
new_html=old_html.replace('</html>',(data/'Table_S11_sections.html').read_text()+'</html>')
html_path.write_text(new_html)
print('Rendered two Table S11 pages; appended to retained S1-S10 proof and HTML.')
