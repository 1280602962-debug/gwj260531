"""Rebuild JCIM artwork from results/canonical and explicit frozen inputs.

Figure contract (evaluation paper, five main figures + five SI figures):
  1 setup/supply → 2 candidate-ranking consequence → 3 ligand chemistry
  → 4 pocket correspondence → 5 computational robustness
  S1 chemistry detail → S2 PIK3CA protocol → S3 cognate RMSD
  → S4 label/source robustness → S5 external-data eligibility.

Analysis-derived AUROC/CI/n/EF/Top-10/RMSD values are read from results/canonical.
Frozen experimental/eligibility inputs are used only when they cannot be recomputed.
"""
from pathlib import Path
import argparse
import csv
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
from PIL import Image
import jcim_figure_style as style
import plot_jcim_article_figures_v3 as v

OUT = Path(__file__).resolve().parents[1]
C = style.C
PAIRS = style.PRIMARY_PAIRS
SMIN = r'summary$_{\mathrm{min}}$'
GENERATED = []
READS = {}
P = {}
SOURCE = None
FLAGSHIP = {'EGFR/HER2', 'JAK1/TYK2'}


def source_path(rel):
    return SOURCE / rel


def read(rel):
    path = source_path(rel)
    READS[rel] = str(path)
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def save(fig, stem, toc=False):
    GENERATED.append(stem)
    style.save_all(fig, stem, toc=toc)
    plt.close(fig)


def label(ax, letter, x=-0.12, y=1.04):
    style.panel_label(ax, letter, x=x, y=y)


def L(color, marker, text, ms=4.2, mfc=None, mec=None, ls='none', lw=1.0):
    """Legend handle matching the plotted marker exactly."""
    edge = mec or color
    face = color if mfc is None else mfc
    return Line2D([], [], color=edge, marker=marker, ls=ls, lw=lw, ms=ms,
                  markerfacecolor=face, markeredgecolor=edge, label=text)


def dotci(ax, x, lo, hi, y, color, marker='o', ms=4.1, mfc=None, mec=None, z=4):
    edge = mec or color
    face = color if mfc is None else mfc
    ax.plot([lo, hi], [y, y], color=edge, lw=1.1, zorder=z - 1)
    ax.plot(x, y, marker=marker, color=edge, mfc=face, mec=edge, ms=ms, ls='none', zorder=z)


def pair_yticks(ax, fontsize=7, egfr_note=False):
    labels = [p + ('†' if egfr_note and p == 'EGFR/HER2' else '') for p in PAIRS]
    ax.set_yticks(range(len(PAIRS)), labels, fontsize=fontsize)
    ax.set_ylim(len(PAIRS) - 0.35, -0.65)


def equal_src(D, pair):
    return D['fixed']


def ranking_row(D, pair):
    return D['ranking'][pair]


def operating_point_row(D, pair):
    return D['and_filter'][pair]


def draw_census_and_primary(ax, census):
    """Compact supply drop; eight-pair box is not drawn as 86→8."""
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    ax.text(.50, .98, 'Paired data decrease under four-state requirements',
            ha='center', va='top', fontsize=7.6, fontweight='bold')
    keys = ['n_pairs_n_both_ge_1', 'n_pairs_n_both_ge_10', 'n_directional_n10', 'n_strict_thick']
    names = ['≥1 both-end\nligand', '≥10 both-end\nligands', 'four-state\nn≥10 (θ=6.0)', 'strict\n6.5/5.5']
    xs = [.125, .375, .625, .875]
    w, h = .20, .36
    y0 = .46
    for i, (k, txt, xc) in enumerate(zip(keys, names, xs)):
        ax.add_patch(FancyBboxPatch((xc - w / 2, y0), w, h, boxstyle='round,pad=.006',
                                    fc='#F4F7FA', ec='#D5DDE4', lw=.7))
        ax.text(xc, y0 + .24, format(int(census[k]), ','), ha='center', va='center',
                fontsize=7.3, fontweight='bold', color=C['vina'])
        ax.text(xc, y0 + .08, txt, ha='center', va='center', fontsize=6.0, color='#555555',
                linespacing=1.05)
        if i < 3:
            ax.text((xc + xs[i + 1]) / 2, y0 + h / 2, '→', ha='center', va='center',
                    fontsize=9, color='#888888')
    ax.add_patch(FancyBboxPatch((.03, .05), .94, .34, boxstyle='round,pad=.008',
                                fc='#FFF8F0', ec=C['a_only'], lw=1.05))
    ax.text(.50, .27, 'Primary evaluation: 8 target pairs', ha='center', va='center',
            fontsize=7.6, fontweight='bold')
    ax.text(.50, .14, 'selected using panel-construction and structural criteria',
            ha='center', va='center', fontsize=6.3, color='#555555')


def build_fig1_setup(D):
    census = next(r for r in D['census'] if r['slice'] == 'all')
    fig = plt.figure(figsize=(7, 4.55))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.70, 0.90], wspace=.34, hspace=.30)
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'A')
    ax.set(xlim=(0, 2), ylim=(0, 2)); ax.set_aspect('equal')
    for x, y, name, col in [(0, 1, 'dual', C['dual']), (1, 1, 'A-only', C['a_only']),
                            (0, 0, 'B-only', C['b_only']), (1, 0, 'neither', C['neither'])]:
        ax.add_patch(plt.Rectangle((x, y), 1, 1, fc=col, alpha=.23, ec='white', lw=2))
        ax.text(x + .5, y + .5, name, ha='center', va='center', fontsize=8)
    ax.set_xticks([.5, 1.5], [r'B $\geq\theta$', r'B $<\theta$'], fontsize=7)
    ax.set_yticks([.5, 1.5], [r'A $<\theta$', r'A $\geq\theta$'], fontsize=7)
    ax.tick_params(length=0); ax.set_title('Experimental activity states', fontsize=8)
    for s in ax.spines.values():
        s.set_visible(False)

    ax = fig.add_subplot(gs[0, 1]); label(ax, 'B'); ax.axis('off')
    for y, left, right, col in [
        (.70, 'dual — A-only', 'score B', C['vina']),
        (.42, 'dual — B-only', 'score A', C['a_only']),
    ]:
        ax.add_patch(FancyBboxPatch((.04, y - .10), .92, .20, boxstyle='round,pad=.02', fc='white', ec=col, lw=1))
        ax.text(.10, y, left, ha='left', va='center', fontsize=7.2, fontweight='bold')
        ax.text(.90, y, right, ha='right', va='center', fontsize=7.2, color=col)
    ax.text(.5, .14, r'$\mathrm{summary}_{\mathrm{min}}=\min[\mathrm{AUROC}_{D/A}(B),\,\mathrm{AUROC}_{D/B}(A)]$',
            ha='center', fontsize=6.4)
    ax.set_title('Directional evaluation', fontsize=8)

    ax = fig.add_subplot(gs[1, :]); style.panel_label(ax, 'C', x=0, y=1.02)
    draw_census_and_primary(ax, census)
    P['fig1C'] = census
    fig.subplots_adjust(left=.16, right=.97, top=.93, bottom=.04)
    save(fig, 'Fig1_four_state_and_supply')


def build_fig2_ranking_consequence(D):
    fig, axs = plt.subplots(4, 1, figsize=(7, 8.50),
                            gridspec_kw={'height_ratios': [1.12, 1.00, 1.10, 0.88]})
    for ax, letter in zip(axs[:3], 'ABC'):
        label(ax, letter)
        pair_yticks(ax, fontsize=7)
    label(axs[3], 'D')

    # A — fixed-score ΔAUROC forest. Encoding matches panel B:
    # blue circle = A-only comparison / target B score
    # orange square = B-only comparison / target A score
    for i, p in enumerate(PAIRS):
        src = equal_src(D, p)
        for j, (contrast, col, m) in enumerate([
            ('D_vs_A_or_neither_pocketB', C['vina'], 'o'),
            ('D_vs_B_or_neither_pocketA', C['a_only'], 's'),
        ]):
            rr = src[(p, contrast)]
            y = i + (-.16 if j == 0 else .16)
            dlt = float(rr['delta_neither_minus_selective'])
            lo, hi = float(rr['delta_ci_lo']), float(rr['delta_ci_hi'])
            dotci(axs[0], dlt, lo, hi, y, col, m)
            if p in FLAGSHIP and contrast.endswith('pocketA'):
                axs[0].text(hi + 0.018, y, f'{dlt:.3f}', va='center', fontsize=6.4, color=col)
    axs[0].axvline(0, color=C['chance'], ls='--', lw=.8)
    axs[0].set(xlim=(-.58, .78),
               xlabel=r'$\Delta$AUROC (dual–neither $-$ dual–single-target-active)')
    axs[0].legend(handles=[
        L(C['vina'], 'o', 'B score, A-only → neither', ms=4.4),
        L(C['a_only'], 's', 'A score, B-only → neither', ms=4.2),
    ], loc='upper center', bbox_to_anchor=(.5, -.22), ncol=2, fontsize=6.2)

    # B — directional dumbbell (was panel A)
    for i, p in enumerate(PAIRS):
        r = v.primary_row(D, p)
        axs[1].plot([r['da'], r['db']], [i, i], color='#B8B8B8', lw=.85, zorder=2)
        axs[1].plot(r['da'], i, 'o', color=C['vina'], ms=4.3, zorder=4)
        axs[1].plot(r['db'], i, 's', color=C['a_only'], ms=4.1, zorder=4)
    axs[1].axvline(.5, color=C['chance'], ls='--', lw=.8)
    axs[1].set(xlim=(0.18, 1.02), xlabel='AUROC')
    axs[1].legend(handles=[
        L(C['vina'], 'o', 'Dual vs A-only, target B', ms=4.4),
        L(C['a_only'], 's', 'Dual vs B-only, target A', ms=4.2),
    ], loc='upper center', bbox_to_anchor=(.5, -.22), ncol=2, fontsize=6.2)

    # C — descriptive summary_min vs dual–neither (was panel B)
    for i, p in enumerate(PAIRS):
        r = v.primary_row(D, p)
        dotci(axs[2], r['smin'], r['lo'], r['hi'], i - .16, C['vina'], 'o')
        nei_m = 'D' if r['n_neg'] < 10 else 's'
        dotci(axs[2], r['nei'], r['nei_lo'], r['nei_hi'], i + .16, C['desc'], nei_m)
        if r['n_neg'] < 10:
            axs[2].text(r['nei'] + 0.025, i + .16, f"n={r['n_neg']}",
                        va='center', fontsize=6.0, color=C['desc'])
    axs[2].axvline(.5, color=C['chance'], ls='--', lw=.8)
    axs[2].set(xlim=(0.08, 1.05), xlabel='AUROC')
    n_neither = {p: int(v.primary_row(D, p)['n_neg']) for p in PAIRS}
    small_n = sorted({n for n in n_neither.values() if n < 10})
    legend_c = [
        L(C['vina'], 'o', SMIN, ms=4.4),
        L(C['desc'], 's', 'Dual vs neither, mean Vina', ms=4.2),
    ]
    if small_n:
        legend_c.append(L(C['desc'], 'D', f'neither n={small_n[0]}', ms=4.2))
    axs[2].legend(handles=legend_c, loc='upper center', bbox_to_anchor=(.5, -.22),
                  ncol=len(legend_c), fontsize=6.1)

    axd = axs[3]
    ranks = [ranking_row(D, p) for p in PAIRS]
    keys = ['top_dual', 'top_A_only', 'top_B_only', 'top_neither']
    cols = [C['dual'], C['a_only'], C['b_only'], C['neither']]
    for i, op in enumerate(ranks):
        vals = [int(op[k]) for k in keys]
        k = max(sum(vals), 1)
        left = 0.0
        for val, col in zip(vals, cols):
            axd.barh(i, val / k, left=left, color=col, height=0.62, linewidth=0)
            left += val / k
    axd.set_xlim(0, 1.02)
    axd.set_xlabel('Fraction of top 10%')
    axd.set_yticks(range(len(PAIRS)))
    axd.set_yticklabels(
        [f"{p}  k={int(ranking_row(D, p)['top_k'])}" for p in PAIRS],
        fontsize=6.4,
    )
    axd.invert_yaxis()
    axd.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc=C['dual'], label='dual'),
        plt.Rectangle((0, 0), 1, 1, fc=C['a_only'], label='A-only'),
        plt.Rectangle((0, 0), 1, 1, fc=C['b_only'], label='B-only'),
        plt.Rectangle((0, 0), 1, 1, fc=C['neither'], label='neither'),
    ], loc='upper center', bbox_to_anchor=(.5, -0.28), ncol=4, fontsize=6.1, frameon=False)

    P['fig2'] = {
        p: {
            'primary': v.primary_row(D, p),
            'fixed_A': equal_src(D, p)[(p, 'D_vs_B_or_neither_pocketA')],
            'fixed_B': equal_src(D, p)[(p, 'D_vs_A_or_neither_pocketB')],
        } for p in PAIRS
    }
    P['fig2D'] = {p: ranking_row(D, p) for p in PAIRS}
    P['fig2C_neither_n'] = n_neither
    fig.subplots_adjust(left=.22, right=.97, top=.97, bottom=.10, hspace=.86)
    save(fig, 'Fig2_negative_class_formulation')


def build_fig3_ligand_chemistry(D):
    fig = plt.figure(figsize=(7, 6.40))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.18, 1.00], hspace=.38)
    ax = fig.add_subplot(gs[0, 0])
    label(ax, 'A', x=-0.08, y=1.08)
    off = {'vina_da': .30, 'vina_db': .10, 'ecfp_da': -.10, 'ecfp_db': -.30}
    cols = {'vina_da': C['vina'], 'vina_db': C['vina'], 'ecfp_da': C['desc'], 'ecfp_db': C['desc']}
    marks = {'vina_da': 'o', 'vina_db': 's', 'ecfp_da': 'o', 'ecfp_db': 's'}
    filled = {'vina_da': C['vina'], 'vina_db': C['vina'], 'ecfp_da': 'white', 'ecfp_db': 'white'}
    plotted = {k: [] for k in off}
    for i, p in enumerate(PAIRS):
        a = v.ecfp_row(D, p, 'D_vs_A')
        b = v.ecfp_row(D, p, 'D_vs_B')
        vals = {'vina_da': a['vina'], 'vina_db': b['vina'], 'ecfp_da': a['ecfp'], 'ecfp_db': b['ecfp']}
        for k, val in vals.items():
            plotted[k].append(val)
            ax.plot(val, i + off[k], marks[k], color=cols[k], mfc=filled[k], mec=cols[k], ms=5.0, zorder=4)
    ax.axvline(.5, color=C['chance'], ls='--', lw=.85)
    ax.set_yticks(range(len(PAIRS)), PAIRS, fontsize=6.5)
    ax.set_ylim(7.55, -1.85)
    ax.set_xlabel('AUROC'); ax.set_xlim(.20, 1.05)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Vina D/A', ms=5.2),
        L(C['vina'], 's', 'Vina D/B', ms=4.8),
        L(C['desc'], 'o', 'ECFP4 D/A', ms=5.2, mfc='white', mec=C['desc']),
        L(C['desc'], 's', 'ECFP4 D/B', ms=4.8, mfc='white', mec=C['desc']),
    ], loc='upper center', ncol=4, fontsize=6.5, frameon=True, fancybox=False,
       edgecolor='none', facecolor='white', framealpha=.92,
       columnspacing=.9, handletextpad=.35)
    P['fig3A'] = plotted

    ax = fig.add_subplot(gs[1, 0])
    label(ax, 'B', x=-0.08, y=1.06)
    deltas = []
    ink = C['ink']
    for i, p in enumerate(PAIRS):
        for contrast, yoff, m in (('D_vs_A', .12, 'o'), ('D_vs_B', -.12, 's')):
            r = v.ecfp_row(D, p, contrast)
            dlt = r['ecfp_plus'] - r['ecfp_base']
            deltas.append(dlt)
            ax.plot([0, dlt], [i + yoff, i + yoff], color='#888888', lw=.9, zorder=2)
            ax.plot(dlt, i + yoff, m, color=ink, ms=4.2, zorder=4)
    ax.axvline(0, color=C['ink'], lw=.8)
    ax.set_yticks(range(len(PAIRS)), PAIRS, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel(r'$\Delta$AUROC (ECFP4+Vina − ECFP4)')
    ax.set_xlim(-0.03, 0.03)
    ax.legend(handles=[
        L(ink, 'o', 'D/A', ms=4.4),
        L(ink, 's', 'D/B', ms=4.2),
    ], loc='lower right', bbox_to_anchor=(.98, .06), fontsize=6.5, frameon=True,
       fancybox=False, edgecolor='none', facecolor='white', framealpha=.92)
    P['fig3B_deltas'] = deltas
    canon_max = max(
        abs(float(r['delta_ECFP4_plus_docking_minus_ECFP4']))
        for r in D['canon_ecfp'].values()
    )
    P['fig3B_max_abs'] = canon_max
    if abs(float(max(abs(d) for d in deltas)) - canon_max) > 1e-12:
        print('FAIL: fig3B_max_abs does not match results/canonical/ecfp4_incremental_information.csv',
              file=sys.stderr)
        raise SystemExit(1)
    fig.subplots_adjust(left=.16, right=.98, top=.92, bottom=.09)
    save(fig, 'Fig3_ligand_chemistry')


def tpsa_panel(ax, D):
    rng = np.random.default_rng(20260729)
    data = [D['tpsa']['dual'], D['tpsa']['A_only'], D['tpsa']['B_only']]
    colors = [C['dual'], C['a_only'], C['b_only']]
    ns = [len(d) for d in data]
    for i, (vals, col) in enumerate(zip(data, colors), start=1):
        arr = np.asarray(vals, float)
        jitter = rng.uniform(-0.12, 0.12, size=len(arr))
        ax.scatter(np.full(len(arr), i) + jitter, arr, s=10, color=col, alpha=.75,
                   edgecolors=C['ink'], linewidths=.25, zorder=3)
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        ax.plot([i - .22, i + .22], [med, med], color=C['ink'], lw=1.3, zorder=4)
        ax.plot([i, i], [q1, q3], color=C['ink'], lw=1.0, zorder=4)
    ax.set_xticks([1, 2, 3], [f'dual\nn={ns[0]}', f'A-only\nn={ns[1]}', f'B-only\nn={ns[2]}'], fontsize=6.2)
    ax.set_ylabel(r'TPSA ($\mathrm{\AA}^2$)')
    return {'n': ns, 'mean': [float(np.mean(d)) for d in data], 'median': [float(np.median(d)) for d in data]}


def build_fig5_computational_robustness(D):
    fig = plt.figure(figsize=(7, 6.20))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.22, 1.18], hspace=.52, wspace=.34)
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'A', x=-0.18, y=1.06)
    vina_smin, vina_nei, g_smin, g_nei = [], [], [], []
    yy = np.arange(len(v.GNINA_INDEP_PAIRS))
    for i, p in enumerate(v.GNINA_INDEP_PAIRS):
        pr, g = v.primary_row(D, p), v.gnina_indep(D, p)
        vina_smin.append(pr['smin']); vina_nei.append(pr['nei'])
        g_smin.append(g['smin']); g_nei.append(g['nei'])
        # Gray connectors emphasize the task gap, not engine rank.
        ax.plot([pr['smin'], pr['nei']], [i - .14, i - .14], color='#B8B8B8', lw=.9, zorder=2)
        ax.plot([g['smin'], g['nei']], [i + .14, i + .14], color='#B8B8B8', lw=.9, zorder=2)
        ax.plot(pr['smin'], i - .14, 'o', color=C['vina'], ms=5.0, zorder=4)
        ax.plot(pr['nei'], i - .14, 'o', color=C['desc'], ms=5.0, zorder=4)
        ax.plot(g['smin'], i + .14, 's', color=C['vina'], mfc='white', mec=C['vina'], ms=4.8, zorder=4)
        ax.plot(g['nei'], i + .14, 's', color=C['desc'], mfc='white', mec=C['desc'], ms=4.8, zorder=4)
    ax.axvline(.5, color=C['chance'], ls='--', lw=.85)
    ax.set_yticks(yy, v.GNINA_INDEP_PAIRS, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel('AUROC'); ax.set_xlim(.12, .95)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Vina', ms=5.0),
        L(C['vina'], 's', 'GNINA', ms=4.8, mfc='white', mec=C['vina']),
        L(C['vina'], 'o', SMIN, ms=5.0),
        L(C['desc'], 'o', 'Dual vs neither', ms=5.0),
    ], loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=6.0, frameon=False,
       columnspacing=.8, labelspacing=.28, handletextpad=.35)
    P['fig4A'] = {
        'pairs': list(v.GNINA_INDEP_PAIRS),
        'vina_smin': vina_smin, 'gnina_smin': g_smin,
        'vina_neither': vina_nei, 'gnina_neither': g_nei,
    }

    ax = fig.add_subplot(gs[0, 1]); label(ax, 'B', x=-0.18, y=1.06)
    items = [
        ('4L23 / 4JT6', v.primary_row(D, 'PIK3CA/mTOR')['smin'],
         v.primary_row(D, 'PIK3CA/mTOR')['lo'], v.primary_row(D, 'PIK3CA/mTOR')['hi'], C['vina'], 'primary'),
        ('4JPS', float(D['jps']['summary_min']), float(D['jps']['summary_min_ci_lo']),
         float(D['jps']['summary_min_ci_hi']), C['a_only'], 'PIK3CA swapped'),
        ('5DXT', float(D['dxt']['summary_min']), float(D['dxt']['summary_min_ci_lo']),
         float(D['dxt']['summary_min_ci_hi']), C['a_only'], 'PIK3CA swapped'),
        ('4JSX', float(D['jsx']['summary_min']), float(D['jsx']['summary_min_ci_lo']),
         float(D['jsx']['summary_min_ci_hi']), C['holdout'], 'mTOR swapped'),
    ]
    fig4b = []
    for i, (name, y, lo, hi, col, grp) in enumerate(items):
        ax.errorbar(i, y, yerr=[[y - lo], [hi - y]], fmt='o', color=col, ecolor=col,
                    elinewidth=1.2, capsize=2.0, markersize=6.0, zorder=4)
        fig4b.append({'label': name, 'y': y, 'lo': lo, 'hi': hi, 'group': grp})
    ax.axhline(.5, color=C['chance'], ls='--', lw=.85)
    ax.set_xticks(range(4), [n for n, *_ in items], fontsize=6.5)
    ax.set_ylabel(SMIN); ax.set_ylim(.12, 1.02); ax.set_xlim(-.55, 3.45)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Primary', ms=5.2),
        L(C['a_only'], 'o', 'PIK3CA substituted', ms=5.2),
        L(C['holdout'], 'o', 'mTOR substituted', ms=5.2),
    ], loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=6.0, frameon=False,
       columnspacing=.8, labelspacing=.28, handletextpad=.35)
    P['fig4B'] = fig4b

    ax = fig.add_subplot(gs[1, :]); label(ax, 'C', x=-0.08, y=1.04)
    plotted_s = {}
    for i, p in enumerate(PAIRS):
        r = v.five_seed_range(D, p)
        ax.plot([r['min'], r['max']], [i, i], color=C['vina'], lw=1.4, zorder=3)
        ax.plot(r['median'], i, 'o', color=C['vina'], ms=5.2, zorder=4)
        ax.plot(r['primary'], i, 'D', color=C['desc'], ms=4.4, zorder=5)
        plotted_s[p] = r
    ax.axvline(.5, color=C['chance'], ls='--', lw=.85)
    pair_yticks(ax, fontsize=6.5)
    ax.set_xlabel(SMIN + ' across five Vina seeds')
    ax.set_xlim(.22, .82)
    ax.legend(handles=[
        Line2D([], [], color=C['vina'], ls='-', lw=1.4, marker='', label='range'),
        L(C['vina'], 'o', 'median', ms=5.0),
        L(C['desc'], 'D', 'primary seed', ms=4.6),
    ], loc='upper center', bbox_to_anchor=(.5, -.18), ncol=3, fontsize=6.2)
    P['fig4C'] = plotted_s
    fig.subplots_adjust(left=.16, right=.98, top=.93, bottom=.13)
    save(fig, 'Fig5_computational_realization')


def build_fig4_pocket_correspondence(D):
    fig, axs = plt.subplots(2, 1, figsize=(7, 6.80))
    ax = axs[0]; label(ax, 'A'); pair_yticks(ax, fontsize=7)
    recs_a = []
    for i, p in enumerate(PAIRS):
        m = v.wp_main(D, p)
        recs_a.append({'pair': p, 'y': m['delta'], 'lo': m['lo'], 'hi': m['hi'], 'excl': m['excl']})
        dotci(ax, m['delta'], m['lo'], m['hi'], i, C['main'], 'o')
    ax.axvline(0, color=C['chance'], ls='--', lw=.85)
    ax.set_xlim(-.28, .40)
    ax.set_xlabel(r'$\Delta$' + SMIN + ' (matched − mismatched), main panel')
    ax.legend(handles=[L(C['main'], 'o', 'main panel', ms=4.6)],
              loc='upper center', bbox_to_anchor=(.5, -.18), ncol=1, fontsize=6.3)
    P['fig5A'] = recs_a

    ax = axs[1]; label(ax, 'B'); pair_yticks(ax, fontsize=7, egfr_note=True)
    recs_b = []
    for i, p in enumerate(PAIRS):
        if p == 'EGFR/HER2':
            recs_b.append({'pair': p, 'missing': True})
            continue
        h = v.wp_hold(D, p)
        recs_b.append({'pair': p, 'y': h['delta'], 'lo': h['lo'], 'hi': h['hi'], 'excl': h['excl']})
        dotci(ax, h['delta'], h['lo'], h['hi'], i, C['holdout'], 's')
    ax.axvline(0, color=C['chance'], ls='--', lw=.85)
    ax.set_xlim(-.28, .40)
    ax.set_xlabel(r'$\Delta$' + SMIN + ' (matched − mismatched), unused-pool internal holdout')
    ax.legend(handles=[L(C['holdout'], 's', 'unused-pool internal holdout', ms=4.4)],
              loc='upper center', bbox_to_anchor=(.5, -.18), ncol=1, fontsize=6.3)
    P['fig5B'] = recs_b
    fig.text(.57, .018, r'$\dagger$ EGFR/HER2 has no unused-pool internal holdout', ha='center', fontsize=6.4, color='#555555')
    fig.subplots_adjust(left=.20, right=.97, top=.96, bottom=.14, hspace=.50)
    save(fig, 'Fig4_mismatched_pocket')


def counts_heatmap(ax, rows, columns, title, gate):
    mat = np.array([[int(r[k]) for k in columns] for r in rows])
    ax.imshow(np.minimum(mat / gate, 1), cmap='Blues', vmin=0, vmax=1, aspect='auto')
    for i in range(len(rows)):
        for j in range(3):
            ax.text(j, i, str(mat[i, j]), ha='center', va='center', fontsize=7,
                    color='white' if mat[i, j] >= gate * .65 else C['ink'])
    ax.set_yticks(range(len(rows)), [r['pair'] for r in rows], fontsize=6.5)
    ax.set_xticks(range(3), ['Dual', 'A-only', 'B-only'], fontsize=7)
    ax.set_title(title, fontsize=8, pad=6)
    ax.tick_params(length=0)
    return mat.tolist()


def build_figS4_label_source(D):
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.80))
    ax = axs[0]; label(ax, 'A', x=-0.18, y=1.04)
    rules = ['theta_5.5', 'theta_6.0', 'theta_6.5', 'strict_6.5_5.5']
    recs = [[v.theta_grid_record(D, p, r) for r in rules] for p in PAIRS]
    mat = np.array([[r['value'] for r in row] for row in recs])
    im = ax.imshow(mat, cmap='RdBu', vmin=.25, vmax=.75, aspect='auto')
    for i in range(8):
        for j in range(4):
            ax.text(j, i, f"{mat[i, j]:.3f}" + ('†' if recs[i][j]['under'] else ''),
                    ha='center', va='center', fontsize=6,
                    color='white' if mat[i, j] < .35 or mat[i, j] > .67 else C['ink'])
    ax.set_yticks(range(8), PAIRS, fontsize=6.5)
    ax.set_xticks(range(4), ['θ=5.5', 'θ=6.0', 'θ=6.5', 'strict\n6.5/5.5'], fontsize=6.1)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=.045, pad=.10)
    cb.set_label(SMIN, fontsize=7)
    ax.text(0.0, 8.55, r'$\dagger$ min($n_{\mathrm{dual}}$, $n_{\mathrm{A}}$ , $n_{\mathrm{B}}$) < 10',
            fontsize=6.0, color='#555555', clip_on=False)
    P['fig6A'] = mat.tolist()

    ax = axs[1]; label(ax, 'B', x=-0.28, y=1.04)
    clusters = D['cluster']
    ylabels, plotted = [], {}
    y = 0
    styles = {
        'ligand_stratified': (C['vina'], 'o', 'ligand'),
        'scaffold_cluster': (C['desc'], 's', 'scaffold'),
        'document_cluster': (C['a_only'], 'D', 'document'),
    }
    for pair in ('EGFR/HER2', 'JAK1/TYK2'):
        ylabels.append(pair)
        y += 1
        for key in ('ligand_stratified', 'scaffold_cluster', 'document_cluster'):
            rr = next((z for z in clusters if z['pair'] == pair and z['estimator'] == key), None)
            if rr is None or rr.get('delta_ci_lo') in ('', None):
                continue
            col, m, lab = styles[key]
            ylabels.append(lab)
            dotci(ax, float(rr['delta_point']), float(rr['delta_ci_lo']), float(rr['delta_ci_hi']), y, col, m)
            plotted[(pair, key)] = rr
            y += 1
    ax.axvline(0, color=C['chance'], ls='--', lw=.7)
    ax.set_yticks(range(len(ylabels)), ylabels, fontsize=6.5)
    for tick, lab in zip(ax.get_yticklabels(), ylabels):
        if lab in ('EGFR/HER2', 'JAK1/TYK2'):
            tick.set_fontweight('bold')
            tick.set_fontsize(7.0)
    ax.set_ylim(len(ylabels) - .4, -.6)
    ax.set(xlim=(-.12, .78), xlabel=r'$\Delta$AUROC, target A score')
    P['fig6B'] = {f'{p}|{k}': rec for (p, k), rec in plotted.items()}
    P['figS11'] = clusters
    fig.subplots_adjust(left=.16, right=.98, top=.90, bottom=.16, wspace=.42)
    save(fig, 'FigS4_label_source_robustness')


def build_figS5_external_eligibility(D):
    rows = [next(r for r in D['native'] if r['pair'] == p) for p in PAIRS]
    fig2, axs = plt.subplots(1, 2, figsize=(7, 4.2))
    for ax, cols, title, gate, letter in [
        (axs[0], ['n_dual', 'n_A_only', 'n_B_only'], 'Compounds (criterion: n≥20 per class)', 20, 'A'),
        (axs[1], ['n_sources_dual', 'n_sources_A_only', 'n_sources_B_only'], 'Sources (criterion: ≥3 per class)', 3, 'B'),
    ]:
        label(ax, letter); counts_heatmap(ax, rows, cols, title, gate)
    fig2.text(.57, .02, '0/8 pairs met the full external docking gate. Eligibility screen, not external validation.',
              ha='center', fontsize=6.5)
    fig2.subplots_adjust(left=.16, right=.98, wspace=.55, bottom=.14, top=.87)
    P['figS5'] = rows
    save(fig2, 'FigS5_external_eligibility')


def build_figS2_protocol(D):
    fig, axs = plt.subplots(1, 2, figsize=(7, 2.85))
    y48 = float(D['protocol'][('PM48', 'E16_primary')]['summary_min'])
    y110 = float(D['protocol'][('PM110', 'E16')]['summary_min'])
    e16 = float(D['protocol'][('PM48', 'E16_primary')]['summary_min'])
    e8 = float(D['protocol'][('PM48', 'E8')]['summary_min'])
    ax = axs[0]; label(ax, 'A', x=-0.18, y=1.08)
    ax.plot([0, 1], [y48, y110], '-', color=C['vina'], lw=1.0)
    ax.plot(0, y48, 'o', color=C['vina'], ms=6)
    ax.plot(1, y110, 's', color=C['holdout'], ms=5.6)
    ax.axhline(.5, color=C['chance'], ls='--', lw=.7)
    ax.set_xticks([0, 1], ['PM48\n(primary)', 'PM110'], fontsize=6.4)
    ax.set(xlim=(-.35, 1.35), ylim=(.45, .82), ylabel=SMIN)
    ax = axs[1]; label(ax, 'B', x=-0.18, y=1.08)
    ax.plot([0, 1], [e16, e8], '-', color=C['vina'], lw=1.0)
    ax.plot(0, e16, 'o', color=C['vina'], ms=6)
    ax.plot(1, e8, 's', color=C['gnina'], ms=5.6)
    ax.axhline(.5, color=C['chance'], ls='--', lw=.7)
    ax.set_xticks([0, 1], ['E = 16', 'E = 8'], fontsize=7)
    ax.set(xlim=(-.35, 1.35), ylim=(.45, .82), ylabel=SMIN)
    P['figS2'] = {'PM48': y48, 'PM110': y110, 'E16': e16, 'E8': e8}
    fig.subplots_adjust(left=.10, right=.98, wspace=.38, top=.88, bottom=.18)
    save(fig, 'FigS2_protocol_sensitivity')


def cognate_rmsd_rows(D):
    """Assemble 14 main-receptor RMSD markers from the canonical CalcRMS table."""
    unified = D['rmsd']
    by_pdb = {}
    n_by_pdb = {}
    for r in unified:
        by_pdb[r['pdb']] = {
            'protein': r['protein'], 'pdb': r['pdb'],
            'top1': float(r['calcrrms_top1_A']),
            'best': float(r['calcrrms_best_A']),
            'source': 'results/canonical/cognate_rmsd.csv',
            'pose_status': r['pose_status'],
        }
        n_by_pdb[r['pdb']] = int(r['n_modes'])
    order = [
        ('EGFR', '3POZ'), ('HER2', '3RCD'), ('JAK1', '6N7A'), ('JAK2', '8BXH'),
        ('TYK2', '3LXP'), ('PIK3CA', '4L23'), ('mTOR', '4JT6'), ('AChE', '4EY7'),
        ('BChE', '4BDS'), ('F2', '4UDW'), ('F10', '2JKH'), ('PPARG', '9V8H'),
        ('PPARA', '6LXA'), ('PPARD', '5U3Q'),
    ]
    rows = []
    for protein, pdb in order:
        rec = dict(by_pdb[pdb])
        rec['protein'] = protein
        rec['n_modes'] = n_by_pdb.get(pdb, '')
        rec['label'] = f'{protein} {pdb}'
        rec['best9'] = rec['best']
        rows.append(rec)
    return rows


def build_figS3_cognate(D):
    rows = cognate_rmsd_rows(D)
    xmax = 3.5
    fig, ax = plt.subplots(figsize=(7, 6.35))
    offscale = []
    for i, r in enumerate(rows):
        pts = [('top1', r['top1'], i - .14, 'o', C['vina']),
               ('best', r['best'], i + .14, 'D', C['a_only'])]
        for name, val, yy, m, col in pts:
            if val <= xmax:
                ax.plot(val, yy, m, color=col, ms=4.5, zorder=4)
            else:
                offscale.append({'label': r['label'], 'series': name, 'rmsd': val})
                ax.annotate(
                    '',
                    xy=(xmax + 0.22, yy),
                    xytext=(xmax - 0.05, yy),
                    arrowprops=dict(arrowstyle='-|>', color=col, lw=0.9),
                    annotation_clip=False,
                )
                ax.text(xmax + 0.26, yy, f'{val:.3f}', va='center', ha='left', fontsize=6.5, color=col)
    ax.axvline(2.0, color=C['chance'], ls='--', lw=.9)
    ax.set_yticks(range(len(rows)), [r['label'] for r in rows], fontsize=7)
    ax.set_ylim(len(rows) - .35, -.65)
    ax.set_xlabel(r'Heavy-atom RMSD ($\mathrm{\AA}$)')
    ax.set_xlim(-0.12, 4.35)
    ax.legend(handles=[
        L(C['vina'], 'o', 'top-1', ms=4.8),
        L(C['a_only'], 'D', 'lowest saved-pose RMSD', ms=4.4),
    ], loc='upper center', bbox_to_anchor=(.5, -.11), ncol=2, fontsize=6.5)
    P['figS12'] = rows
    P['figS12_offscale'] = offscale
    fig.subplots_adjust(left=.18, right=.90, top=.93, bottom=.14)
    save(fig, 'FigS3_cognate_rmsd')


def toc_graphic():
    fig = plt.figure(figsize=(3.25, 1.75), dpi=300)
    ax = fig.add_axes([0.02, 0.06, 0.96, 0.88])
    ax.set_xlim(0, 32); ax.set_ylim(0, 16); ax.axis('off')
    ax.text(5.2, 14.7, 'Four states', ha='center', fontsize=8.0, fontweight='bold')
    cells = [(1.15, 7.6, 'dual', C['dual']), (5.55, 7.6, 'A-only', C['a_only']),
             (1.15, 3.15, 'B-only', C['b_only']), (5.55, 3.15, 'neither', C['neither'])]
    for x, y, name, col in cells:
        ax.add_patch(Rectangle((x, y), 4.15, 3.7, fc=col, alpha=.22, ec=col, lw=.5))
        ax.text(x + 2.07, y + 1.85, name, ha='center', va='center', fontsize=7.0)
    ax.text(5.2, 1.35, r'B $\geq\theta$          B $<\theta$', ha='center', fontsize=6.5, color='#555555')
    ax.text(11.0, 8.4, '→', ha='center', va='center', fontsize=11, color='#888888')
    ax.text(17.4, 14.7, 'Same docking', ha='center', fontsize=8.0, fontweight='bold')
    ax.add_patch(FancyBboxPatch((11.7, 8.55), 11.4, 4.55, boxstyle='round,pad=0.06,rounding_size=0.18',
                                fc='#F7F7F7', ec=C['vina'], lw=.7))
    ax.text(17.4, 11.85, 'dual — A-only  |  score B', ha='center', fontsize=6.5)
    ax.text(17.4, 10.05, 'dual — B-only  |  score A', ha='center', fontsize=6.5)
    ax.add_patch(FancyBboxPatch((11.7, 3.15), 11.4, 4.55, boxstyle='round,pad=0.06,rounding_size=0.18',
                                fc='#F7F7F7', ec=C['desc'], lw=.7))
    ax.text(17.4, 5.45, 'dual — neither', ha='center', fontsize=6.5)
    ax.text(23.2, 8.4, '→', ha='center', va='center', fontsize=11, color='#888888')
    ax.text(27.7, 14.7, 'Controls', ha='center', fontsize=8.0, fontweight='bold')
    ax.add_patch(FancyBboxPatch((23.7, 3.15), 7.9, 9.95, boxstyle='round,pad=0.06,rounding_size=0.18',
                                fc='#F4F7FA', ec='#CCCCCC', lw=.7))
    ax.text(27.65, 10.15, 'ligand-only', ha='center', fontsize=6.5, fontweight='bold')
    ax.text(27.65, 6.05, 'pocket', ha='center', fontsize=6.5, fontweight='bold')
    ax.text(27.65, 4.55, 'correspondence', ha='center', fontsize=6.5)
    save(fig, 'TOC_graphic', toc=True)


def build_figS1_chemistry(D):
    fig, axs = plt.subplots(1, 2, figsize=(7, 4.40), gridspec_kw={'width_ratios': [1.35, 0.85]})
    ax = axs[0]; label(ax, 'A', x=-0.18, y=1.04)
    plotted = {}
    for i, p in enumerate(PAIRS):
        r = v.primary_row(D, p)
        name, dval = v.best_desc(D, p)
        plotted[p] = {'vina_smin': r['smin'], 'desc_name': name, 'desc': dval}
        ax.plot([r['lo'], r['hi']], [i + 0.12, i + 0.12], color=C['vina'], lw=1.5, zorder=3)
        ax.plot(r['smin'], i + 0.12, 'o', color=C['vina'], markersize=5.2, zorder=4)
        ax.plot(dval, i - 0.14, 's', color=C['desc'], markersize=4.6, zorder=4)
    ax.axvline(0.5, color=C['chance'], ls='--', lw=0.85)
    ax.set_yticks(range(len(PAIRS)), PAIRS, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel(SMIN + ' AUROC')
    ax.set_xlim(0.12, 1.02)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Vina', ms=5.2),
        L(C['desc'], 's', 'best single descriptor', ms=4.6),
    ], loc='lower right', fontsize=6.2, frameon=False)
    P['figS1A'] = plotted
    ax = axs[1]; label(ax, 'B', x=-0.22, y=1.04)
    P['figS1B'] = tpsa_panel(ax, D)
    fig.subplots_adjust(left=.16, right=.98, wspace=.38, top=.90, bottom=.14)
    save(fig, 'FigS1_ligand_chemistry_detail')


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    if isinstance(obj, Path):
        return str(obj)
    return obj


def emit_plotted_values_postfix(audit):
    """Authoritative plotted-value lock. One file only."""
    rows = []

    def add(figure, panel, pair, metric, source, source_row, raw, displayed=None):
        if raw is None or raw == '':
            return
        try:
            raw_f = float(raw)
            disp = displayed if displayed is not None else f'{raw_f:.3f}'
            raw_out = raw_f
        except (TypeError, ValueError):
            raw_out = raw
            disp = displayed if displayed is not None else str(raw)
        rows.append({
            'figure': figure,
            'panel': panel,
            'pair': pair,
            'metric': metric,
            'source_file': source,
            'source_row': source_row,
            'raw_value': raw_out,
            'displayed_value': disp,
        })

    census = P.get('fig1C') or {}
    for key in ['n_pairs_n_both_ge_1', 'n_pairs_n_both_ge_10', 'n_directional_n10', 'n_strict_thick']:
        if key in census:
            add('Figure 1', 'C', 'census', key,
                'data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv',
                'slice=all', census[key], format(int(census[key]), ','))

    for p, rec in (P.get('fig2') or {}).items():
        src = 'results/canonical/fixed_score_negative_class_delta.csv'
        fa = rec.get('fixed_A') or {}
        fb = rec.get('fixed_B') or {}
        add('Figure 2', 'A', p, 'fixed_delta_pocketA', src, p,
            fa.get('delta_neither_minus_selective'))
        add('Figure 2', 'A', p, 'fixed_delta_pocketB', src, p,
            fb.get('delta_neither_minus_selective'))
        pr = rec.get('primary') or {}
        add('Figure 2', 'B', p, 'AUROC_D_vs_A_pocketB',
            'results/canonical/primary_directional_auroc.csv', p, pr.get('da'))
        add('Figure 2', 'B', p, 'AUROC_D_vs_B_pocketA',
            'results/canonical/primary_directional_auroc.csv', p, pr.get('db'))
        add('Figure 2', 'C', p, 'summary_min',
            'results/canonical/primary_summary_min.csv', p, pr.get('smin'))
        add('Figure 2', 'C', p, 'D_vs_neither_vina_mean',
            'results/canonical/two_pocket_mean_ranking.csv', p, pr.get('nei'))
        add('Figure 2', 'C', p, 'n_neither',
            'results/canonical/class_counts.csv', p, pr.get('n_neg'),
            str(pr.get('n_neg')))

    for p, op in (P.get('fig2D') or {}).items():
        src = 'results/canonical/top10_operating_points.csv'
        for col in ['top_dual', 'top_A_only', 'top_B_only', 'top_neither', 'top_k']:
            add('Figure 2', 'D', p, col, src, p, op.get(col), str(op.get(col)))

    fig3_src = 'results/canonical/ecfp4_incremental_information.csv'
    fig3a = P.get('fig3A') or {}
    for i, p in enumerate(PAIRS):
        add('Figure 3', 'A', p, 'vina_D_vs_A', fig3_src, p,
            (fig3a.get('vina_da') or [None] * len(PAIRS))[i])
        add('Figure 3', 'A', p, 'vina_D_vs_B', fig3_src, p,
            (fig3a.get('vina_db') or [None] * len(PAIRS))[i])
        add('Figure 3', 'A', p, 'ECFP4_D_vs_A', fig3_src, p,
            (fig3a.get('ecfp_da') or [None] * len(PAIRS))[i])
        add('Figure 3', 'A', p, 'ECFP4_D_vs_B', fig3_src, p,
            (fig3a.get('ecfp_db') or [None] * len(PAIRS))[i])
    deltas = P.get('fig3B_deltas') or []
    k = 0
    for p in PAIRS:
        for contrast in ('D_vs_A', 'D_vs_B'):
            if k < len(deltas):
                add('Figure 3', 'B', p, f'delta_{contrast}', fig3_src, p, deltas[k])
                k += 1
    add('Figure 3', 'B', 'all', 'fig3B_max_abs', fig3_src, 'max_|delta|',
        P.get('fig3B_max_abs'))

    for p, rec in (P.get('figS1A') or {}).items():
        add('Figure S1', 'A', p, 'vina_summary_min',
            'results/canonical/primary_summary_min.csv', p, rec.get('vina_smin'))
        add('Figure S1', 'A', p, 'best_descriptor_auroc',
            'results/canonical/descriptor_baselines.csv', p, rec.get('desc'))
        add('Figure S1', 'A', p, 'best_descriptor_name',
            'results/canonical/descriptor_baselines.csv', p, rec.get('desc_name'),
            rec.get('desc_name'))

    tpsa = P.get('figS1B') or {}
    if tpsa:
        add('Figure S1', 'B', 'AChE/BChE', 'n_dual_A_B',
            'results/canonical/current_score_master.csv', 'AChE/BChE TPSA from SMILES',
            ','.join(str(x) for x in tpsa.get('n', [])),
            '/'.join(str(x) for x in tpsa.get('n', [])))

    for rec in P.get('figS12') or []:
        add('Figure S3', 'A', rec.get('label'), 'top1_rmsd_A',
            'results/canonical/cognate_rmsd.csv',
            rec.get('pdb'), rec.get('top1'))
        add('Figure S3', 'A', rec.get('label'), 'lowest_saved_rmsd_A',
            'results/canonical/cognate_rmsd.csv',
            rec.get('pdb'), rec.get('best'))

    for rec in P.get('fig5A') or []:
        add('Figure 4', 'A', rec.get('pair'), 'matched_minus_mismatched_delta',
            'results/canonical/matched_mismatched_pocket.csv',
            rec.get('pair'), rec.get('y'))
    for rec in P.get('fig5B') or []:
        if rec.get('missing'):
            add('Figure 4', 'B', rec.get('pair'), 'holdout_status',
                'results/canonical/holdout_metrics.csv',
                rec.get('pair'), 'no_holdout', 'no unused-pool holdout')
            continue
        add('Figure 4', 'B', rec.get('pair'), 'holdout_matched_minus_mismatched_delta',
            'results/canonical/holdout_metrics.csv',
            rec.get('pair'), rec.get('y'))

    g = P.get('fig4A') or {}
    for i, p in enumerate(g.get('pairs') or []):
        add('Figure 5', 'A', p, 'vina_summary_min',
            'results/canonical/primary_summary_min.csv', p,
            g['vina_smin'][i])
        add('Figure 5', 'A', p, 'gnina_summary_min',
            'results/canonical/computational_robustness.csv', p,
            g['gnina_smin'][i])
        add('Figure 5', 'A', p, 'vina_D_vs_neither',
            'results/canonical/two_pocket_mean_ranking.csv', p,
            g['vina_neither'][i])
        add('Figure 5', 'A', p, 'gnina_D_vs_neither',
            'results/canonical/computational_robustness.csv', p,
            g['gnina_neither'][i])
    for rec in P.get('fig4B') or []:
        add('Figure 5', 'B', 'PIK3CA/mTOR', rec.get('label'),
            'results/canonical/receptor_substitution.csv', rec.get('label'), rec.get('y'))
    for p, rec in (P.get('fig4C') or {}).items():
        add('Figure 5', 'C', p, 'five_seed_min',
            'results/canonical/five_seed_summary_min.csv', p, rec.get('min'))
        add('Figure 5', 'C', p, 'five_seed_max',
            'results/canonical/five_seed_summary_min.csv', p, rec.get('max'))
        add('Figure 5', 'C', p, 'five_seed_median',
            'results/canonical/five_seed_summary_min.csv', p, rec.get('median'))
        add('Figure 5', 'C', p, 'five_seed_primary',
            'results/canonical/primary_summary_min.csv', p, rec.get('primary'))

    proto = P.get('figS2') or {}
    for metric, val in proto.items():
        add('Figure S2', 'A' if metric.startswith('PM') else 'B', 'PIK3CA/mTOR', metric,
            'results/canonical/protocol_sensitivity.csv', metric, val)

    mat = P.get('fig6A') or []
    rules = ['theta_5.5', 'theta_6.0', 'theta_6.5', 'strict_6.5_5.5']
    for i, p in enumerate(PAIRS):
        if i < len(mat):
            for j, rule in enumerate(rules):
                add('Figure S4', 'A', p, rule,
                    'results/canonical/label_aggregation_sensitivity.csv', p, mat[i][j])
    for key, rec in (P.get('fig6B') or {}).items():
        add('Figure S4', 'B', key, 'cluster_delta',
            'results/canonical/cluster_bootstrap_sensitivity.csv', key,
            rec.get('delta_point') if isinstance(rec, dict) else None)
    for rec in P.get('figS5') or []:
        add('Figure S5', 'A', rec.get('pair'), 'n_dual',
            'results/canonical/external_eligibility.csv', rec.get('pair'), rec.get('n_dual'),
            str(rec.get('n_dual')))
        add('Figure S5', 'B', rec.get('pair'), 'n_sources_dual',
            'results/canonical/external_eligibility.csv', rec.get('pair'), rec.get('n_sources_dual'),
            str(rec.get('n_sources_dual')))

    payload = {
        'schema': 'plotted_values_postfix.v4',
        'lock': 'docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md',
        'generated': audit.get('generated'),
        'n_records': len(rows),
        'records': rows,
        'nested_plotted': _jsonable(P),
        'input_files': audit.get('input_files'),
    }
    (OUT / 'plotted_values_postfix.json').write_text(
        json.dumps(payload, indent=2, default=str), encoding='utf-8')
    plotted_json = OUT / 'plotted_values.json'
    if plotted_json.is_file():
        plotted_json.unlink()


def main():
    global SOURCE
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--dry-run', action='store_true', help='validate inputs without writing canonical figures')
    args = parser.parse_args()
    SOURCE = args.source_root.resolve()
    style.OUT = OUT
    style.apply_style()
    if style.FONT != 'Arial' or style.ARIAL_PATH is None:
        print('FAIL: Arial is required for canonical artwork', file=sys.stderr)
        raise SystemExit(1)
    v.ROOT = SOURCE
    v.DATA = SOURCE / 'data'
    v.CANON = SOURCE / 'results' / 'canonical'
    v.OUT = OUT
    v.save_all = save
    D = v.load()
    if set(D['smin']) != set(PAIRS):
        print('FAIL: canonical pair set mismatch', file=sys.stderr)
        raise SystemExit(1)
    native_pairs = set(r['pair'] for r in D['native'])
    if native_pairs != set(PAIRS):
        print('FAIL: native pair set mismatch', file=sys.stderr)
        raise SystemExit(1)
    if not all(float(r['packaged_as_external_evaluation']) == 0 for r in D['native']):
        print('FAIL: holdout packaged as external evaluation', file=sys.stderr)
        raise SystemExit(1)
    if args.dry_run:
        print('dry-run: figure inputs validated; canonical artwork not written')
        return
    P['pair_order'] = list(PAIRS)
    build_fig1_setup(D)
    build_fig2_ranking_consequence(D)
    build_fig3_ligand_chemistry(D)
    build_fig4_pocket_correspondence(D)
    build_fig5_computational_robustness(D)
    build_figS1_chemistry(D)
    build_figS2_protocol(D)
    build_figS3_cognate(D)
    build_figS4_label_source(D)
    build_figS5_external_eligibility(D)
    toc_graphic()
    P.update({k: val for k, val in v.PROVENANCE['plotted'].items() if k not in P})
    P['primary'] = {p: v.primary_row(D, p) for p in PAIRS}
    egfr_delta = float(D['fixed'][('EGFR/HER2', 'D_vs_B_or_neither_pocketA')]['delta_neither_minus_selective'])
    jak_delta = float(D['fixed'][('JAK1/TYK2', 'D_vs_B_or_neither_pocketA')]['delta_neither_minus_selective'])
    P['flagship_fixed_delta'] = {'EGFR/HER2': egfr_delta, 'JAK1/TYK2': jak_delta}
    P['and_filter_table_s10'] = D['and_filter']['JAK1/TYK2']
    audit = {
        'source': 'results/canonical',
        'input_files': {rel: rel for rel in READS},
        'generated': GENERATED,
        'plotted': P,
    }
    emit_plotted_values_postfix(audit)
    for stem in GENERATED:
        for ext in ['png', 'tif']:
            with Image.open(OUT / (stem + '.' + ext)) as im:
                if im.mode != 'RGB':
                    print(f'FAIL: {stem}.{ext} mode {im.mode} != RGB', file=sys.stderr)
                    raise SystemExit(1)
                if abs(im.info['dpi'][0] - 300) >= 1:
                    print(f'FAIL: {stem}.{ext} dpi {im.info["dpi"]}', file=sys.stderr)
                    raise SystemExit(1)
                if im.width > 2101 or im.height > 2751:
                    print(f'FAIL: {stem}.{ext} size {im.width}x{im.height}', file=sys.stderr)
                    raise SystemExit(1)
    print(f'PASS: {len(GENERATED)} figures; {len(READS)} extra frozen reads; RGB / 300 dpi / size checks')


if __name__ == '__main__':
    main()
