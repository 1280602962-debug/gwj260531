"""Physical layout and contrast checks at final 174-mm publication width."""
import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.text import Annotation, Text
from matplotlib.transforms import Bbox
import numpy as np

REPORT=[]


def luminance(rgb):
    rgb=np.asarray(rgb)
    lin=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
    return float(lin@np.array([.2126,.7152,.0722]))


def readable_color(rgb):
    lum=luminance(rgb)
    black=(lum+.05)/.05;white=1.05/(lum+.05)
    return ('black',black) if black>=white else ('white',white)


def visible_text(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();skip=set()
    for ax in fig.axes:
        for vals,labels,lim in [(ax.get_xticks(),ax.get_xticklabels(),ax.get_xlim()),
                                (ax.get_yticks(),ax.get_yticklabels(),ax.get_ylim())]:
            lo,hi=sorted(lim)
            skip.update(id(t) for v,t in zip(vals,labels) if not lo-1e-8<=v<=hi+1e-8)
    return [(t,t.get_window_extent(renderer)) for t in fig.findobj(Text)
            if t.get_visible() and t.get_text().strip() and id(t) not in skip]


def margins(fig):
    items=visible_text(fig);w,h=fig.canvas.get_width_height();u=25.4/fig.dpi
    return dict(left=min(b.x0 for _,b in items)*u,right=(w-max(b.x1 for _,b in items))*u,
                bottom=min(b.y0 for _,b in items)*u,top=(h-max(b.y1 for _,b in items))*u)


def pad_canvas(fig,minimum=4.6):
    # Move axes, figure text and figure legends together; font sizes stay fixed.
    for _ in range(5):
        m=margins(fig)
        if min(m.values())>=minimum-.02:return
        w,h=fig.get_size_inches()*25.4
        l=max(0,minimum-m['left'])/w;r=max(0,minimum-m['right'])/w
        b=max(0,minimum-m['bottom'])/h;t=max(0,minimum-m['top'])/h
        for ax in fig.axes:
            p=ax.get_position();ax.set_position([l+p.x0*(1-l-r),b+p.y0*(1-b-t),p.width*(1-l-r),p.height*(1-b-t)])
        for text in fig.texts:
            x,y=text.get_position();text.set_position((l+x*(1-l-r),b+y*(1-b-t)))
        for leg in fig.legends:
            bb=leg.get_bbox_to_anchor().transformed(fig.transFigure.inverted())
            leg.set_bbox_to_anchor((l+bb.x0*(1-l-r),b+bb.y0*(1-b-t)),transform=fig.transFigure)
    assert min(margins(fig).values())>=4.4,margins(fig)


def legend_gap(fig,axes):
    fig.canvas.draw();rr=fig.canvas.get_renderer()
    legend_top=max(l.get_window_extent(rr).y1 for l in fig.legends)
    xlabel_bottom=min(a.xaxis.label.get_window_extent(rr).y0 for a in axes if a.get_xlabel())
    return (xlabel_bottom-legend_top)*25.4/fig.dpi


def ensure_legend_gap(fig,axes,target=3.8):
    # Raise each affected plotting bottom while retaining the plotting top.
    gap=legend_gap(fig,axes)
    if gap<target:
        d=(target-gap)/(fig.get_figheight()*25.4)
        for ax in axes:
            p=ax.get_position();ax.set_position([p.x0,p.y0+d,p.width,p.height-d])


def relabel_numbers(ax,minimum_gap_mm=2.0):
    points=[]
    for text in list(ax.texts):
        if isinstance(text,Annotation) and text.get_text().isdigit():
            points.append((int(text.get_text()),*text.xy));text.remove()
        elif isinstance(text,Annotation) and not text.get_text() and text.arrow_patch is not None:text.remove()
    if not points:return None
    ax.figure.canvas.draw();rr=ax.figure.canvas.get_renderer()
    gap=minimum_gap_mm*ax.figure.dpi/25.4
    point_boxes=[];occupied=[]
    for _,x,y in points:
        px,py=ax.transData.transform((x,y));point_boxes.append(Bbox.from_bounds(px-5-gap/2,py-5-gap/2,10+gap,10+gap))
    candidates=[(12,12),(-12,12),(12,-14),(-12,-14),(22,8),(-22,8),(8,24),(-8,-24),
                (29,18),(-29,18),(29,-20),(-29,-20),(0,34),(0,-34),(40,0),(-40,0),
                (42,26),(-42,26),(42,-28),(-42,-28)]
    for label,x,y in points:
        for dx,dy in candidates:
            t=ax.annotate(str(label),(x,y),xytext=(dx,dy),textcoords='offset points',
                          fontsize=8,ha='center',va='center',zorder=5)
            box=t.get_window_extent(rr).padded(gap/2)
            if not any(box.overlaps(z) for z in occupied+point_boxes) and ax.bbox.contains(box.x0,box.y0) and ax.bbox.contains(box.x1,box.y1):
                occupied.append(box)
                ax.annotate('',(x,y),xytext=(dx*.78,dy*.78),textcoords='offset points',
                            arrowprops={'arrowstyle':'-','color':'#777777','lw':.5},zorder=2)
                break
            t.remove()
        else:raise ValueError('Cannot place separated numeric label '+str(label))
    return {'labels':len(points),'minimum_label_gap_mm':minimum_gap_mm}


def polish(fig,number):
    before=margins(fig);extra={};heat=[]
    if number in ['4','5','S12']:
        for t in fig.texts:
            if t.get_text()=='Negative class':t.set_x(.21)
        for ax in fig.axes[:2]:
            p=ax.get_position();ax.set_position([p.x0,p.y0,p.width,p.height-.025])
    if number=='S3':
        ax=fig.axes[-1];p=ax.get_position();ax.set_position([p.x0,p.y0+.04,p.width,p.height])
    if number=='5':
        fig.canvas.draw();rr=fig.canvas.get_renderer();cb=fig.axes[-1]
        bottom=min(t.get_window_extent(rr).y0 for t in fig.axes[0].get_xticklabels() if t.get_visible())
        p=cb.get_position();new_top=(bottom-3.8*fig.dpi/25.4)/(fig.get_figheight()*fig.dpi)
        cb.set_position([p.x0,new_top-p.height,p.width,p.height])
    if number=='6':
        fig.canvas.draw();rr=fig.canvas.get_renderer()
        key=[t for t in fig.texts if t.get_text().startswith(('1 EGFR','5 AChE'))]
        top=max(t.get_window_extent(rr).y1 for t in key)
        bottom=min(l.get_window_extent(rr).y0 for l in fig.legends)
        d=(bottom-top-2.8*fig.dpi/25.4)/(fig.get_figheight()*fig.dpi)
        for t in key:x,y=t.get_position();t.set_position((x,y+d))
    if number=='8':
        leg=fig.legends[0];handles=list(leg.legend_handles);labels=[t.get_text() for t in leg.get_texts()];leg.remove()
        handles.append(Line2D([],[],marker='o',linestyle='none',color='#555555',markerfacecolor='white',
                              label='Open symbols: S4–S1 (supporting)'))
        labels.append('Open symbols: S4–S1 (supporting)')
        fig.legend(handles=handles,labels=labels,loc='lower center',bbox_to_anchor=(.55,.12),ncol=2,frameon=False)
        # Preserve equal panel heights when reserving room for the extra key.
        fig.subplots_adjust(bottom=.31)
    if number=='S8':
        for t in fig.texts:
            t.set_text(t.get_text().replace('lower 95% CI 0.4922 to 0.5030',
                         'lower bound of the 95% CI changes from 0.4922 to 0.5030'))
    if number=='S1':
        fig.set_figheight(7.2)
        for ax in fig.axes:
            collections=[c for c in ax.collections if isinstance(c,PathCollection)]
            for i,c in enumerate(collections):
                v=np.asarray(c.get_offsets()).copy();v[:,1]+=np.linspace(-.28,.28,5)[i]-np.linspace(-.16,.16,5)[i]
                c.set_offsets(v);c.set_sizes([12])
        extra['seed_points_per_panel']=70
        extra['seed_display_offsets']=[-.28,-.14,0,.14,.28]
    for ax in fig.axes:
        for im in ax.images:
            if im.cmap.name=='cividis':
                colors=im.cmap(np.linspace(0,1,256));colors[:,:3]=.8*colors[:,:3]+.2
                im.set_cmap(ListedColormap(colors,name='cividis_lighter'))
            grid=np.asarray(im.get_array())
            if grid.ndim!=2:continue
            for t in ax.texts:
                x,y=t.get_position();i,j=round(y),round(x)
                if abs(y-i)<1e-8 and abs(x-j)<1e-8 and 0<=i<grid.shape[0] and 0<=j<grid.shape[1]:
                    rgb=im.cmap(im.norm(grid[i,j]))[:3];color,ratio=readable_color(rgb);t.set_color(color)
                    heat.append({'text':t.get_text(),'row':i,'column':j,'axes_index':fig.axes.index(ax),'foreground':color,'contrast':ratio})
    if number in ['2','4','S12'] and fig.legends:ensure_legend_gap(fig,fig.axes[:2])
    pad_canvas(fig)
    if number in ['6','8']:
        extra['numeric_labels']=[relabel_numbers(ax) for ax in fig.axes]
    after=margins(fig)
    for item in heat:
        ax=fig.axes[item['axes_index']];x,y=ax.transData.transform((item['column']+.35,item['row']+.35))
        cx,cy=ax.transData.transform((item['column'],item['row']))
        item['background_sample_pdf_pt']=[x*72/fig.dpi,(fig.get_figheight()*fig.dpi-y)*72/fig.dpi]
        item['text_center_pdf_pt']=[cx*72/fig.dpi,(fig.get_figheight()*fig.dpi-cy)*72/fig.dpi]
    assert min(after.values())>=4.4,(number,after)
    if heat:assert min(r['contrast'] for r in heat)>=4.5,(number,'low contrast')
    if number in ['2','4','S12']:
        extra['axis_legend_gap_mm']=legend_gap(fig,fig.axes[:2])
        assert extra['axis_legend_gap_mm']>=3.4,(number,extra)
    if number in ['4','5','S12']:
        fig.canvas.draw();rr=fig.canvas.get_renderer()
        header=next(t for t in fig.texts if t.get_text()=='Negative class').get_window_extent(rr)
        title=fig.axes[0].title.get_window_extent(rr) if fig.axes[0].title.get_text() else fig.axes[0]._left_title.get_window_extent(rr)
        assert not header.overlaps(title),(number,'header/title overlap')
    if number=='6':
        fig.canvas.draw();rr=fig.canvas.get_renderer()
        key=[t for t in fig.texts if t.get_text().startswith(('1 EGFR','5 AChE'))]
        extra['legend_key_gap_mm']=(min(l.get_window_extent(rr).y0 for l in fig.legends)-max(t.get_window_extent(rr).y1 for t in key))*25.4/fig.dpi
        assert 2<=extra['legend_key_gap_mm']<=3.2,extra
    REPORT.append({'figure':number,'before_margin_mm':before,'after_margin_mm':after,
                   'heatmap_labels':heat,**extra})


def save_report(out):
    (Path(out)/'Checks/Publication_Layout_QA.json').write_text(json.dumps({'minimum_margin_target_mm':4.6,
        'minimum_heatmap_contrast':4.5,'figures':REPORT},indent=2)+'\n')
