"""Refine visual encodings of frozen estimates; do not fit or resample data."""
from pathlib import Path
import csv
import importlib.util
import json
import shutil
import xml.etree.ElementTree as ET

import fitz
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.collections import PathCollection

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('md_artwork', HERE/'build_molecular_diversity_artwork.py')
md = importlib.util.module_from_spec(spec)
spec.loader.exec_module(md)
b = md.b
OUT = md.OUT
PREVIOUS = OUT/'Comparison_Baseline'
original_finish = md.finish


def finish(fig, number, name, section, table, sources, records, caption, caption_cn, si=False):
    if number=='S2':
        by={(r['pair'],r['method'],r['class'],r['target']):r for r in records}
        maximum=max(int(r['expected'])-int(r['valid']) for r in records)
        assert all(int(r['missing'])==int(r['expected'])-int(r['valid']) for r in records)
        for ax,(cls,target) in zip(fig.axes,[(c,t) for c in ['dual','A_only','B_only','neither'] for t in ['A','B']]):
            missing=np.array([[int(by[p,m,cls,target]['missing']) for m in b.METHODS] for p in b.ORDER])
            ax.images[0].set_data(missing);ax.images[0].set_cmap('Blues');ax.images[0].set_clim(0,maximum)
        fig._missing_count_scale=maximum
        caption=caption.replace('its color represents that fraction',
          'its color gives the number of unavailable or invalid scores on a common 0-'+str(maximum)+' scale, with white indicating complete availability')
    for ax in fig.axes:
        xmin,xmax=sorted(ax.get_xlim());ymin,ymax=sorted(ax.get_ylim())
        for c in ax.collections:
            if isinstance(c,PathCollection):
                offsets=np.asarray(c.get_offsets())
                if not offsets.size:continue
                assert np.all((offsets[:,0]>=xmin-1e-8)&(offsets[:,0]<=xmax+1e-8)),(number,'clipped x coordinate')
                assert np.all((offsets[:,1]>=ymin-1e-8)&(offsets[:,1]<=ymax+1e-8)),(number,'clipped y coordinate')
    # Main figures plus tables must not exceed ten display items.
    if number == '8':
        number, si = 'S10', True
    elif number == '9':
        number = '8'
    plt.rcParams['svg.fonttype'] = 'none'
    original_finish(fig, number, name, section, table, sources, records, caption, caption_cn, si)


def chemistry():
    source = b.ST+'module_D_chemistry_units.csv'
    data = b.read(source)
    by = {(r['pair'],r['arm']):r for r in data}
    rr = [by[u] for u in b.units()]
    keys = ['auroc_physchem','auroc_ecfp4','auroc_ecfp4_m0','auroc_SUPPORTING_CHEMISTRY_MATCHED_M0']
    grid = np.array([[b.f(r,k) for k in keys] for r in rr])
    fig, ax = plt.subplots(1,2,figsize=(md.WIDTH,6.6),gridspec_kw={'width_ratios':[1.20,1]})
    fig.subplots_adjust(left=.34,right=.975,bottom=.21,top=.90,wspace=.25)
    im = ax[0].imshow(grid,aspect='auto',cmap='cividis',vmin=0,vmax=1,interpolation='none')
    for i in range(16):
        for j in range(4):
            ax[0].text(j,i,f'{grid[i,j]:.2f}',ha='center',va='center',fontsize=8,
                          color=md.contrast_text(im.cmap(im.norm(grid[i,j]))[:3]))
    ax[0].set_xticks(range(4),['Physchem','ECFP4','ECFP4\n+ S1','Vina S1'])
    ax[0].tick_params(axis='both',length=0)
    ax[0].set_yticks(range(16),['A-only','B-only']*8)
    for i,p in enumerate(b.ORDER):
        ax[0].text(-.30,2*i+.5,p,transform=ax[0].get_yaxis_transform(),ha='right',va='center',fontsize=8)
    for i in range(1,8):ax[0].axhline(2*i-.5,color='white',lw=1.4)
    yy=np.arange(16)
    delta=np.array([b.f(r,'delta_ecfp4_plus_m0') for r in rr])
    ax[1].hlines(yy,0,delta,color='#8F9AA4',lw=1)
    ax[1].scatter(delta,yy,s=25,color=b.GREEN,marker='^',zorder=3)
    for y,d in zip(yy,delta):
        ax[1].annotate(f'{d:+.3f}',(d,y),xytext=(5 if d>=0 else -5,0),
                       textcoords='offset points',ha='left' if d>=0 else 'right',va='center',fontsize=8)
    ax[1].set_xlim(-.036,.057);ax[1].set_xticks([-.02,0,.02,.04])
    ax[1].set_yticks(yy,[]);ax[1].set_ylim(15.5,-.5)
    ax[1].axvline(0,color='#777777',lw=.7,ls=(0,(3,3)))
    ax[1].set_xlabel('AUROC difference')
    for i in range(0,8,2):ax[1].axhspan(2*i-.5,2*i+1.5,color='#F2F4F6',zorder=-2)
    md.title(ax[0],'a','Matched-population AUROC')
    md.title(ax[1],'b','ECFP4 + S1 minus ECFP4')
    cb=fig.colorbar(im,cax=fig.add_axes([.37,.065,.28,.018]),orientation='horizontal')
    cb.set_ticks([0,.5,1]);cb.set_label('AUROC',fontsize=8)
    fig.text(.04,.955,'Target pair',fontsize=8,fontweight='bold')
    fig.text(.23,.955,'Negative class',fontsize=8,fontweight='bold')
    finish(fig,'5','','Results 3.4','Table S7',[source],data,
      'Ligand chemistry and the incremental contribution of Vina. (a) Matched-population AUROC for physicochemical descriptor (Physchem), ECFP4, ECFP4-plus-S1 models and direct Vina S1 scores. Cells show two-decimal estimates on a common 0-1 color scale. The three learned models use five-fold out-of-fold predictions. (b) AUROC(ECFP4 plus S1) minus AUROC(ECFP4); stems start at zero and labels show signed three-decimal differences. The two panels use different numeric scales. Confidence intervals were not calculated for these chemistry comparisons. Dual is the positive class; A-only and B-only rows use target B and A scores, respectively','')


def seeds_receptors():
    seed=b.read(b.ST+'module_E_five_seed.csv')
    summ=b.read(b.ST+'module_E_five_seed_summary.csv')
    alt=b.read(b.ST+'module_F_alt_receptor.csv')
    by={(r['pair'],int(r['seed'])):r for r in seed}
    ss={r['pair']:r for r in summ if r['metric']=='summary_min'}
    ref={p:b.f(by[p,42],'summary_min') for p in b.ORDER}
    fig,ax=plt.subplots(2,1,figsize=(md.WIDTH,7.0),gridspec_kw={'height_ratios':[1.15,1]})
    fig.subplots_adjust(left=.35,right=.965,bottom=.15,top=.88,hspace=.58)
    y=np.arange(8)[::-1]
    for sno,off,c,m in zip([17,29,42,71,101],np.linspace(-.17,.17,5),
                            [b.BLUE,b.ORANGE,b.GREEN,b.PURPLE,'#555555'],['o','s','^','D','v']):
        ax[0].scatter([b.f(by[p,sno],'summary_min')-ref[p] for p in b.ORDER],y+off,
                        s=24,color=c,marker=m,label=str(sno),zorder=3)
    ax[0].hlines(y,[b.f(ss[p],'min')-ref[p] for p in b.ORDER],
                    [b.f(ss[p],'max')-ref[p] for p in b.ORDER],color='#333333',lw=.8)
    ax[0].scatter([b.f(ss[p],'median')-ref[p] for p in b.ORDER],y,marker='|',color='black',s=65,label='Median',zorder=4)
    ax[0].set_yticks(y,[p+f'  [{ref[p]:.3f}]' for p in b.ORDER])
    changes=[b.f(by[p,s],'summary_min')-ref[p] for p in b.ORDER for s in [17,29,42,71,101]]
    bound=max(abs(min(changes)),abs(max(changes)))+.006
    bound=max(.04,np.ceil(bound*100)/100)
    ax[0].set_ylim(-.6,7.6);ax[0].set_xlim(-bound,bound)
    ax[0].set_xticks([-.04,-.02,0,.02,.04]);ax[0].axvline(0,color='#777777',ls=(0,(3,3)),lw=.7)
    ax[0].set_xlabel('Change in minimum AUROC relative to seed 42')
    md.title(ax[0],'a','Search seeds; brackets give seed-42 AUROC')
    fig.legend(*ax[0].get_legend_handles_labels(),frameon=False,ncol=6,
               loc='upper center',bbox_to_anchor=(.63,.98),handletextpad=.3,columnspacing=.8)
    rr=sorted(alt,key=lambda r:(b.ORDER.index(r['pair']),r['replaced_side']))
    yy=np.arange(6)[::-1]
    for field,label,c,m,off in [('delta_AUC_affected','Affected direction',b.BLUE,'o',.13),
                                ('delta_summary_min','Minimum AUROC',b.PURPLE,'D',-.13)]:
        b.interval(ax[1],[b.f(r,field) for r in rr],[b.f(r,field+'_ci_lo') for r in rr],
                      [b.f(r,field+'_ci_hi') for r in rr],yy+off,c,m,label)
    ax[1].set_yticks(yy,[r['pair']+'  '+r['alt_id'].replace('ALT_','')+' ('+r['replaced_side']+')' for r in rr])
    ax[1].set_xlim(-.24,.20);ax[1].set_ylim(-.6,5.6)
    ax[1].set_xlabel('Alternative minus original receptor AUROC')
    ax[1].axvline(0,color='#777777',ls=(0,(3,3)),lw=.7)
    md.title(ax[1],'b','Receptor substitutions')
    ax[1].legend(frameon=False,loc='lower center',bbox_to_anchor=(.5,-.43),ncol=2)
    records=[{'panel':'a',**r} for r in seed]+[{'panel':'a_summary',**r} for r in summ]+[{'panel':'b',**r} for r in rr]
    finish(fig,'7','','Results 3.6','Table S8A,B',
      [b.ST+'module_E_five_seed.csv',b.ST+'module_E_five_seed_summary.csv',b.ST+'module_F_alt_receptor.csv'],records,
      'Sensitivity to search seeds and receptor structures. (a) Seed-specific Vina minimum AUROC minus the corresponding seed-42 value. Brackets in row labels give the absolute seed-42 AUROC. Segments show the observed minimum-to-maximum seed range after subtracting seed 42; black ticks show the median change. These ranges are not confidence intervals. All five seed results are retained; each uses its valid members. AChE/BChE seed 17 has one fewer A-only ligand than seed 42, so the displayed difference includes that membership change. Absolute seed-specific estimates are shown in Figure S7. (b) Paired changes in affected-direction AUROC and minimum AUROC for six single-target receptor replacements, with 95% paired bootstrap intervals. Parentheses identify replaced side A or B; the other receptor is held fixed','')


def similarity():
    source=b.ST+'module_D1_composition.csv';data=b.read(source)
    fig,axs=plt.subplots(1,2,figsize=(md.WIDTH,6.1),sharey=True)
    fig.subplots_adjust(left=.19,right=.975,bottom=.17,top=.90,wspace=.20)
    for a,arm,neg,letter in zip(axs,b.ARMS,['A_only','B_only'],'ab'):
        for i,p in enumerate(b.ORDER):
            for cls,off,c,m in [('dual',.19,b.BLUE,'o'),(neg,-.19,b.ORANGE,'^')]:
                rows=[r for r in data if r['pair']==p and r['arm']==arm and r['class']==cls]
                vals=[b.f(r,'max_ecfp4_tanimoto_other_class') for r in rows]
                pos=7-i+off
                # Each class gets its own summary; all observations remain visible.
                a.boxplot([vals],positions=[pos],vert=False,widths=.27,showfliers=False,
                          manage_ticks=False,patch_artist=True,
                          boxprops={'facecolor':'white','edgecolor':c,'linewidth':.8},
                          medianprops={'color':'black','linewidth':1},
                          whiskerprops={'color':c,'linewidth':.7},capprops={'color':c,'linewidth':.7},zorder=2)
                a.scatter(vals,[pos+.65*b.jitter(r['canonical_ligand_id']) for r in rows],
                          s=8,color=c,marker=m,alpha=.65,linewidths=0,zorder=3,
                          label=cls.replace('_','-') if i==0 else None)
        a.set_yticks(range(7,-1,-1),b.ORDER);a.set_xlim(0,1);a.set_ylim(-.6,7.6)
        a.set_xlabel('Maximum cross-class ECFP4 Tanimoto')
        md.title(a,letter,'dual vs '+neg.replace('_','-'))
        a.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.14),ncol=2)
    finish(fig,'S6','','Results 3.4','Supplementary Data',[source],data,
      'Structural similarity between the two experimental classes within each directional task. Each point is the maximum ECFP4 Tanimoto similarity to any ligand in the other class of the same target pair and direction. Fingerprints use radius 2, 2048 bits, useChirality=False and useFeatures=False. Dual and single-target classes have separate boxes and vertically separated points. Boxes show the median and interquartile range; whiskers extend to the most extreme values within 1.5 IQR. All observations, including those beyond the whiskers, are shown. These are descriptive class-conditioned nearest-neighbor similarities, not independent replicate measurements or a cross-validation leakage test','',si=True)


def postprocess():
    # Remove superseded exports only in this explicitly scoped artwork folder.
    for ext in ['eps','pdf','png','svg']:
        p=OUT/'Main_Figures'/('Fig9.'+ext)
        if p.is_file():p.unlink()
    p=OUT/'Source_Data/Fig9.csv'
    if p.is_file():p.unlink()
    md.figures.sort(key=lambda r:(r['group']!='Main',r['number']=='Scheme_S1',int(r['number'].lstrip('S')) if r['number']!='Scheme_S1' else 99))
    review=fitz.open()
    for r in md.figures:
        doc=fitz.open((md.MAIN if r['group']=='Main' else md.SI)/(r['stem']+'.pdf'))
        review.insert_pdf(doc)
    review.set_toc([[1,'Scheme S1' if r['number']=='Scheme_S1' else 'Fig. '+r['number'],i+1] for i,r in enumerate(md.figures)])
    review.save(OUT/'Figures_Review.pdf');review.close()
    with (OUT/'Figure_Index.csv').open('w',newline='',encoding='utf8') as h:
        w=csv.writer(h);w.writerow(['figure','group','file','citation','tables','sources'])
        for r in md.figures:w.writerow([r['number'],r['group'],('Main_Figures/' if r['group']=='Main' else 'Supplementary_Figures/')+r['stem'],r['section'],r['table'],';'.join(r['sources'])])
    caps=[]
    for r in md.figures:
        label='Scheme S1' if r['number']=='Scheme_S1' else 'Fig. '+r['number']
        caps+=['**'+label+'** '+r['caption'],'']
    (OUT/'Figure_Captions_EN.md').write_text('\n'.join(caps),encoding='utf8')
    # Figure-specific machine data must remain the same apart from renumbering.
    encoding_audit=[]
    for name in ['Fig5','Fig7','FigS6','FigS2']:
        p=OUT/'Source_Data'/(name+'.csv');q=PREVIOUS/(name+'.csv')
        assert p.read_bytes()==q.read_bytes(),name+' source rows changed'
        encoding_audit.append({'figure':name,'source_rows_byte_identical':True})
    vectors=[]
    for r in md.figures:
        p=(md.MAIN if r['group']=='Main' else md.SI)/(r['stem']+'.svg')
        tree=ET.parse(p)
        text_nodes=tree.findall('.//{http://www.w3.org/2000/svg}text')
        assert text_nodes,p
        vectors.append({'figure':r['number'],'svg_text_nodes':len(text_nodes)})
    qa={'main_figures':8,'main_tables':2,'main_display_items':10,'supplementary_figures':10,
        'supplementary_schemes':1,'retained_scientific_estimates':True,'scatter_centers_within_axes':True,
        'data_checks':encoding_audit,'editable_svg':vectors,
        'renumbering':{'previous Figure 8':'Figure S10','previous Figure 9':'Figure 8'},
        'scientific_recomputation':False}
    (md.QA/'Visualization_Encoding_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf8')
    with (OUT/'Figure_Renumbering.csv').open('w',newline='',encoding='utf8') as h:
        w=csv.writer(h);w.writerow(['previous','current','reason'])
        w.writerows([['Figure 8','Figure S10','Membership sensitivity supports the main analysis'],
                     ['Figure 9','Figure 8','Consecutive main-text numbering; eight figures plus two tables']])
    comparison=fitz.open()
    for name in ['Fig5','Fig7','FigS6','FigS2']:
        for version,p in [('Previous',PREVIOUS/(name+'.pdf')),('Revised',(md.SI if 'S' in name else md.MAIN)/(name+'.pdf'))]:
            doc=fitz.open(p);page=comparison.new_page(width=doc[0].rect.width,height=doc[0].rect.height+24)
            page.insert_text((10,16),name+' - '+version,fontsize=10)
            page.show_pdf_page(fitz.Rect(0,24,page.rect.width,page.rect.height),doc,0)
    comparison.save(OUT/'Visualization_Comparison.pdf');comparison.close()
    print(json.dumps(qa))


def main():
    PREVIOUS.mkdir(parents=True,exist_ok=True)
    for name in ['Fig5','Fig7','FigS6','FigS2']:
        for ext,folder in [('pdf',md.SI if 'S' in name else md.MAIN),('csv',md.DATA)]:
            target=PREVIOUS/(name+'.'+ext)
            if not target.exists():shutil.copyfile(folder/(name+'.'+ext),target)
    md.finish=finish;b.finish=finish
    b.fig5=chemistry;b.fig7=seeds_receptors;b.si6=similarity
    md.main()
    postprocess()


if __name__=='__main__':main()
