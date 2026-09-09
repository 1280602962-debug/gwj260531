"""Rebuild JCIM artwork from pinned PR32 inputs (v3 style).

Figure contract (evaluation paper, six main figures):
  1 setup/supply → 2 experimental-state comparison (hero) → 3 ligand chemistry
  → 4 computational realization → 5 pocket correspondence → 6 evidence boundary.

Display order is protein-system grouped. ORIGINAL_THREE / CENSUS_FIVE are CSV
routing keys only and must not appear on figures or captions.

First run: --source-root PATH_TO_PR32/Dual_Target_Docking
Subsequent runs use the small input_snapshot stored alongside the figures.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
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
SNAP = OUT / 'input_snapshot'
SHA = 'abb61a20a04eb6a085ad526876624eadb518c4cc'
C = style.C
PAIRS = style.PRIMARY_PAIRS
SMIN = r'summary$_{\mathrm{min}}$'
GENERATED = []
READS = {}
P = {}
SOURCE = None
FLAGSHIP = {'EGFR/HER2', 'JAK1/TYK2'}


def snapshot_path(rel):
    # Short flat names avoid Windows MAX_PATH for the nested census table paths.
    return SNAP / (hashlib.sha256(rel.encode()).hexdigest()[:16] + Path(rel).suffix)


def source_path(rel):
    return snapshot_path(rel) if SOURCE == SNAP.resolve() else SOURCE / rel


def read(rel):
    path = source_path(rel)
    READS[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
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
    return D['equal'] if pair in v.ORIGINAL_THREE else D['five_s34']


def draw_census_and_primary(ax, census):
    """ChEMBL supply census with an independent 8-pair evaluation box (no 86→8 arrow)."""
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    ax.text(.50, .975, 'ChEMBL data-supply census', ha='center', va='top', fontsize=8, fontweight='bold')
    keys = ['n_pairs_n_both_ge_1', 'n_pairs_n_both_ge_10', 'n_directional_n10', 'n_strict_thick']
    names = ['≥1 ligand measured at both targets',
             '≥10 ligands measured at both targets',
             'Dual, A-only and B-only each ≥10 (θ=6.0)',
             'Strict bidirectional supply criterion (6.5/5.5)']
    tops = [.855, .695, .535, .375]
    h = .112
    for i, (k, txt, yc) in enumerate(zip(keys, names, tops)):
        ax.add_patch(FancyBboxPatch((.03, yc - h / 2), .94, h, boxstyle='round,pad=.006',
                                    fc='#F4F7FA', ec='#D5DDE4', lw=.7))
        ax.text(.07, yc, format(int(census[k]), ','), va='center', fontsize=8.6, fontweight='bold', color=C['vina'])
        ax.text(.30, yc, txt, va='center', fontsize=7.1)
        if i < 3:
            gap = (tops[i] - h / 2 + tops[i + 1] + h / 2) / 2
            ax.text(.50, gap, '↓', ha='center', va='center', fontsize=8, color='#888888')
    ax.add_patch(FancyBboxPatch((.03, .035), .94, .215, boxstyle='round,pad=.008',
                                fc='#FFF8F0', ec=C['a_only'], lw=1.05))
    ax.text(.50, .175, 'Primary evaluation: 8 target pairs', ha='center', va='center',
            fontsize=7.6, fontweight='bold')
    ax.text(.50, .085, 'data availability  ·  structural eligibility  ·  pair-specific panel criteria',
            ha='center', va='center', fontsize=6.3, color='#555555')


def fig1(D):
    census = next(r for r in read('data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv') if r['slice'] == 'all')
    fig = plt.figure(figsize=(7, 5.20))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.48], wspace=.34, hspace=.38)
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
    ax.text(.5, .12, r'summary$_{\mathrm{min}}=\min(\mathrm{AUC}_{D/A},\mathrm{AUC}_{D/B})$',
            ha='center', fontsize=7.0)
    ax.set_title('Directional evaluation', fontsize=8)

    ax = fig.add_subplot(gs[1, :]); style.panel_label(ax, 'C', x=0, y=1.02)
    draw_census_and_primary(ax, census)
    P['fig1C'] = census
    fig.subplots_adjust(left=.16, right=.97, top=.93, bottom=.04)
    save(fig, 'Fig1_four_state_and_supply')


def fig1_c_data(D):
    census = next(r for r in read('data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv') if r['slice'] == 'all')
    fig, ax = plt.subplots(figsize=(7, 3.55))
    style.panel_label(ax, 'C', x=0.0, y=1.06)
    draw_census_and_primary(ax, census)
    P['fig1C_standalone'] = census
    fig.subplots_adjust(left=.02, right=.98, top=.90, bottom=.04)
    save(fig, 'Fig1_C_chEMBL_supply')


def fig2(D):
    fig, axs = plt.subplots(3, 1, figsize=(7, 7.55), gridspec_kw={'height_ratios': [1.12, 1.00, 1.10]})
    for ax, letter in zip(axs, 'ABC'):
        label(ax, letter)
        pair_yticks(ax, fontsize=7)

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
    axs[0].set(xlim=(-.58, .78), xlabel=r'$\Delta$AUROC (Dual–neither $-$ Dual–selective)')
    axs[0].set_title('Fixed score channel; experimental-state comparison', fontsize=8, pad=4)
    axs[0].legend(handles=[
        L(C['vina'], 'o', 'Target B score: A-only → neither', ms=4.4),
        L(C['a_only'], 's', 'Target A score: B-only → neither', ms=4.2),
    ], loc='upper center', bbox_to_anchor=(.5, -.22), ncol=2, fontsize=6.2)

    # B — directional dumbbell (was panel A)
    for i, p in enumerate(PAIRS):
        r = v.primary_row(D, p)
        axs[1].plot([r['da'], r['db']], [i, i], color='#B8B8B8', lw=.85, zorder=2)
        axs[1].plot(r['da'], i, 'o', color=C['vina'], ms=4.3, zorder=4)
        axs[1].plot(r['db'], i, 's', color=C['a_only'], ms=4.1, zorder=4)
    axs[1].axvline(.5, color=C['chance'], ls='--', lw=.8)
    axs[1].set(xlim=(0.18, 1.02), xlabel='AUROC')
    axs[1].set_title('Pocket-matched directional tasks', fontsize=8, pad=4)
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
    axs[2].axvline(.5, color=C['chance'], ls='--', lw=.8)
    axs[2].set(xlim=(0.08, 1.05), xlabel='AUROC')
    axs[2].set_title('Descriptive two-pocket comparison', fontsize=8, pad=4)
    axs[2].legend(handles=[
        L(C['vina'], 'o', 'Directional ' + SMIN, ms=4.4),
        L(C['desc'], 's', 'Dual vs neither, mean score', ms=4.2),
        L(C['desc'], 'D', 'neither n=4', ms=4.2),
    ], loc='upper center', bbox_to_anchor=(.5, -.22), ncol=3, fontsize=6.1)

    P['fig2'] = {
        p: {
            'primary': v.primary_row(D, p),
            'fixed_A': equal_src(D, p)[(p, 'D_vs_B_or_neither_pocketA')],
            'fixed_B': equal_src(D, p)[(p, 'D_vs_A_or_neither_pocketB')],
        } for p in PAIRS
    }
    fig.subplots_adjust(left=.19, right=.97, top=.955, bottom=.088, hspace=.78)
    save(fig, 'Fig2_negative_class_formulation')


def fig3(D):
    fig = plt.figure(figsize=(7, 6.55))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 1.00], hspace=.42, wspace=.32)
    ax = fig.add_subplot(gs[0, :])
    label(ax, 'A', x=-0.08, y=1.08)
    # Color = method (Vina blue, ECFP4 orange); shape = direction (circle D/A, square D/B).
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
    ax.set_title('Vina rank AUROC and ECFP4 scaffold GroupKFold', fontsize=8, pad=6)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Vina D/A', ms=5.2),
        L(C['vina'], 's', 'Vina D/B', ms=4.8),
        L(C['desc'], 'o', 'ECFP4 D/A', ms=5.2, mfc='white', mec=C['desc']),
        L(C['desc'], 's', 'ECFP4 D/B', ms=4.8, mfc='white', mec=C['desc']),
    ], loc='upper center', ncol=4, fontsize=6.0, frameon=True, fancybox=False,
       edgecolor='none', facecolor='white', framealpha=.92,
       columnspacing=.9, handletextpad=.35)
    P['fig3A'] = plotted

    ax = fig.add_subplot(gs[1, 0])
    label(ax, 'B', x=-0.22, y=1.06)
    deltas = []
    for i, p in enumerate(PAIRS):
        for contrast, yoff, col, m in (('D_vs_A', .12, C['vina'], 'o'), ('D_vs_B', -.12, C['a_only'], 's')):
            r = v.ecfp_row(D, p, contrast)
            dlt = r['ecfp_plus'] - r['ecfp_base']
            deltas.append(dlt)
            ax.plot([0, dlt], [i + yoff, i + yoff], color=col, lw=.9, zorder=2)
            ax.plot(dlt, i + yoff, m, color=col, ms=4.2, zorder=4)
    ax.axvline(0, color=C['ink'], lw=.8)
    ax.set_yticks(range(len(PAIRS)), PAIRS, fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel(r'$\Delta$AUROC (ECFP4+Vina − ECFP4)')
    ax.set_xlim(-0.03, 0.03)
    ax.set_title('Increment from adding Vina', fontsize=8, pad=3)
    ax.legend(handles=[
        L(C['vina'], 'o', 'D/A', ms=4.4),
        L(C['a_only'], 's', 'D/B', ms=4.2),
    ], loc='lower right', bbox_to_anchor=(.98, .06), fontsize=6.0, frameon=True,
       fancybox=False, edgecolor='none', facecolor='white', framealpha=.92)
    P['fig3B_deltas'] = deltas
    P['fig3B_max_abs'] = float(max(abs(d) for d in deltas))

    ax = fig.add_subplot(gs[1, 1])
    label(ax, 'C', x=-0.22, y=1.06)
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
    ax.set_title('Illustrative TPSA separation in AChE/BChE', fontsize=8, pad=3)
    P['fig3C'] = {'n': ns, 'mean': [float(np.mean(d)) for d in data], 'median': [float(np.median(d)) for d in data]}
    fig.subplots_adjust(left=.16, right=.98, top=.90, bottom=.09)
    save(fig, 'Fig3_ligand_chemistry')


def fig4(D):
    fig = plt.figure(figsize=(7, 6.15))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.18], hspace=.42, wspace=.34)
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
    ax.set_yticks(yy, v.GNINA_INDEP_PAIRS, fontsize=6.3)
    ax.invert_yaxis()
    ax.set_xlabel('AUROC'); ax.set_xlim(.12, .95)
    ax.set_title('Independent pose generation', fontsize=8, pad=3)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Vina directional', ms=5.0),
        L(C['desc'], 'o', 'Vina Dual vs neither', ms=5.0),
        L(C['vina'], 's', 'GNINA directional', ms=4.8, mfc='white', mec=C['vina']),
        L(C['desc'], 's', 'GNINA Dual vs neither', ms=4.8, mfc='white', mec=C['desc']),
    ], loc='lower left', fontsize=5.5, frameon=True, fancybox=False, edgecolor='none',
       facecolor='white', framealpha=.92, ncol=1, labelspacing=.25)
    P['fig4A'] = {
        'pairs': list(v.GNINA_INDEP_PAIRS),
        'vina_smin': vina_smin, 'gnina_smin': g_smin,
        'vina_neither': vina_nei, 'gnina_neither': g_nei,
    }

    ax = fig.add_subplot(gs[0, 1]); label(ax, 'B', x=-0.18, y=1.06)
    items = [
        ('Primary\n4L23 / 4JT6', v.primary_row(D, 'PIK3CA/mTOR')['smin'],
         v.primary_row(D, 'PIK3CA/mTOR')['lo'], v.primary_row(D, 'PIK3CA/mTOR')['hi'], C['vina'], 'primary'),
        ('PIK3CA\n4JPS', float(D['jps']['summary_min']), float(D['jps']['summary_min_ci_lo']),
         float(D['jps']['summary_min_ci_hi']), C['a_only'], 'PIK3CA swapped'),
        ('PIK3CA\n5DXT', float(D['dxt']['summary_min']), float(D['dxt']['summary_min_ci_lo']),
         float(D['dxt']['summary_min_ci_hi']), C['a_only'], 'PIK3CA swapped'),
        ('mTOR\n4JSX', float(D['jsx']['summary_min']), float(D['jsx']['summary_min_ci_lo']),
         float(D['jsx']['summary_min_ci_hi']), C['holdout'], 'mTOR swapped'),
    ]
    fig4b = []
    for i, (name, y, lo, hi, col, grp) in enumerate(items):
        ax.errorbar(i, y, yerr=[[y - lo], [hi - y]], fmt='o', color=col, ecolor=col,
                    elinewidth=1.2, capsize=2.0, markersize=6.0, zorder=4)
        fig4b.append({'label': name, 'y': y, 'lo': lo, 'hi': hi, 'group': grp})
    ax.axhline(.5, color=C['chance'], ls='--', lw=.85)
    ax.set_xticks(range(4), [n for n, *_ in items], fontsize=6.1)
    ax.set_ylabel(SMIN); ax.set_ylim(.12, 1.02); ax.set_xlim(-.55, 3.45)
    ax.set_title('PIK3CA/mTOR receptor structures', fontsize=8, pad=3)
    ax.legend(handles=[
        L(C['vina'], 'o', 'Primary (4L23 / 4JT6)', ms=5.2),
        L(C['a_only'], 'o', 'PIK3CA substituted', ms=5.2),
        L(C['holdout'], 'o', 'mTOR substituted', ms=5.2),
    ], loc='lower left', fontsize=5.6, frameon=True, fancybox=False, edgecolor='none',
       facecolor='white', framealpha=.92, labelspacing=.25)
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
    ax.set_title('Five Vina seeds (8 pairs × 5 seeds)', fontsize=8, pad=3)
    ax.legend(handles=[
        Line2D([], [], color=C['vina'], ls='-', lw=1.4, marker='', label='five-seed range'),
        L(C['vina'], 'o', 'median of five seeds', ms=5.0),
        L(C['desc'], 'D', 'primary seed', ms=4.6),
    ], loc='upper center', bbox_to_anchor=(.5, -.18), ncol=3, fontsize=6.2)
    P['fig4C'] = plotted_s
    fig.subplots_adjust(left=.16, right=.98, top=.93, bottom=.11)
    save(fig, 'Fig4_computational_realization')


def fig5(D):
    fig, axs = plt.subplots(2, 1, figsize=(7, 6.70))
    ax = axs[0]; label(ax, 'A'); pair_yticks(ax, fontsize=7, egfr_note=True)
    recs_a, recs_b = [], []
    for i, p in enumerate(PAIRS):
        m = v.wp_main(D, p)
        recs_a.append({'pair': p, 'y': m['delta'], 'lo': m['lo'], 'hi': m['hi'], 'excl': m['excl']})
        dotci(ax, m['delta'], m['lo'], m['hi'], i - .16, C['main'], 'o')
        if p == 'EGFR/HER2':
            recs_b.append({'pair': p, 'missing': True})
            continue
        h = v.wp_hold(D, p)
        recs_b.append({'pair': p, 'y': h['delta'], 'lo': h['lo'], 'hi': h['hi'], 'excl': h['excl']})
        dotci(ax, h['delta'], h['lo'], h['hi'], i + .16, C['holdout'], 's')
    ax.axvline(0, color=C['chance'], ls='--', lw=.85)
    ax.set_xlim(-.28, .40)
    ax.set_xlabel(r'$\Delta$' + SMIN + ' (matched − mismatched)')
    ax.set_title('Pocket correspondence on main panels and holdouts', fontsize=8, pad=4)
    ax.legend(handles=[
        L(C['main'], 'o', 'main', ms=4.6),
        L(C['holdout'], 's', 'holdout', ms=4.4),
    ], loc='upper center', bbox_to_anchor=(.5, -.20), ncol=2, fontsize=6.3)
    P['fig5A'] = recs_a
    P['fig5B'] = recs_b  # holdout Δ retained under former key for lock continuity

    ax = axs[1]; label(ax, 'B'); pair_yticks(ax, fontsize=7, egfr_note=True)
    recs_c = []
    for i, p in enumerate(PAIRS):
        main = v.primary_row(D, p)
        dotci(ax, main['smin'], main['lo'], main['hi'], i - .16, C['main'], 'o')
        if p == 'EGFR/HER2':
            recs_c.append({'pair': p, 'missing': True})
            continue
        h = v.holdout_smin(D, p)
        recs_c.append({'pair': p, 'main': main['smin'], 'main_lo': main['lo'], 'main_hi': main['hi'],
                       'hold': h['y'], 'hold_lo': h['lo'], 'hold_hi': h['hi']})
        dotci(ax, h['y'], h['lo'], h['hi'], i + .16, C['holdout'], 's')
    ax.axvline(.5, color=C['chance'], ls='--', lw=.85)
    ax.set_xlim(.10, 1.02)
    ax.set_xlabel(SMIN)
    ax.set_title('Main versus unused-pool holdout ' + SMIN, fontsize=8, pad=4)
    ax.legend(handles=[
        L(C['main'], 'o', 'main', ms=4.6),
        L(C['holdout'], 's', 'holdout', ms=4.4),
    ], loc='upper center', bbox_to_anchor=(.5, -.20), ncol=2, fontsize=6.3)
    P['fig5C'] = recs_c
    fig.text(.57, .015, r'$\dagger$ no unused-pool holdout available', ha='center', fontsize=6.4, color='#555555')
    fig.subplots_adjust(left=.20, right=.97, top=.95, bottom=.08, hspace=.50)
    save(fig, 'Fig5_mismatched_pocket')


def e8_value():
    panel = {r['panel_id']: r for r in read('data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv')}
    scores = read('data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv')
    groups = {k: [] for k in ['dual', 'A_only', 'B_only']}
    for r in scores:
        cls = panel[r['ligand']]['class']
        if cls in groups and r['4L23_affinity_E8'] and r['4JT6_affinity_E8']:
            groups[cls].append([-float(r['4L23_affinity_E8']), -float(r['4JT6_affinity_E8'])])

    def auc(a, b):
        dif = np.asarray(a)[:, None] - np.asarray(b)[None, :]
        return float(np.mean((dif > 0) + .5 * (dif == 0)))

    return min(auc([x[1] for x in groups['dual']], [x[1] for x in groups['A_only']]),
               auc([x[0] for x in groups['dual']], [x[0] for x in groups['B_only']]))


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


def fig6(D):
    fig = plt.figure(figsize=(7, 7.45))
    gs = fig.add_gridspec(2, 2, hspace=.38, wspace=.36)
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'A', x=-0.18, y=1.04)
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
    ax.set_title('Activity definitions', fontsize=8)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=.045, pad=.10)
    cb.set_label(SMIN, fontsize=7)
    P['fig6A'] = mat.tolist()

    ax = fig.add_subplot(gs[0, 1]); label(ax, 'B', x=-0.28, y=1.04)
    clusters = read('data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv')
    ylabels, plotted = [], {}
    y = 0
    for p in ['EGFR/HER2', 'JAK1/TYK2']:
        ylabels.append(p)
        y += 1
        r = v.s34_row(D, p)
        ylabels.append('ligand')
        dotci(ax, r['delta'], r['lo'], r['hi'], y, C['vina'], 'o')
        plotted[(p, 'ligand')] = r
        y += 1
        for key, col, m, lab in [
            ('scaffold_cluster', C['desc'], 's', 'scaffold'),
            ('document_cluster', C['a_only'], 'D', 'document'),
        ]:
            rr = next(z for z in clusters if z['pair'] == p and z['estimator'] == key)
            ylabels.append(lab)
            dotci(ax, float(rr['delta_point']), float(rr['delta_ci_lo']), float(rr['delta_ci_hi']), y, col, m)
            plotted[(p, key)] = rr
            y += 1
    ax.axvline(0, color=C['chance'], ls='--', lw=.7)
    ax.set_yticks(range(len(ylabels)), ylabels, fontsize=5.8)
    for tick, lab in zip(ax.get_yticklabels(), ylabels):
        if lab in ('EGFR/HER2', 'JAK1/TYK2'):
            tick.set_fontweight('bold')
            tick.set_fontsize(6.3)
    ax.set_ylim(len(ylabels) - .4, -.6)
    ax.set(xlim=(-.12, .78), xlabel=r'$\Delta$AUROC, target A score')
    ax.set_title('Cluster-resampling sensitivity\nof fixed-score differences', fontsize=7.4)
    ax.legend(handles=[
        L(C['vina'], 'o', 'ligand', ms=4.4),
        L(C['desc'], 's', 'scaffold cluster', ms=4.2),
        L(C['a_only'], 'D', 'document cluster', ms=4.2),
    ], loc='upper left', fontsize=5.6, frameon=True, fancybox=False, edgecolor='none',
       facecolor='white', framealpha=.92, labelspacing=.25)
    P['fig6B'] = {f'{p}|{k}': rec for (p, k), rec in plotted.items()}
    P['figS11'] = clusters  # same table; SI figure still generated

    rows = [next(r for r in D['native'] if r['pair'] == p) for p in PAIRS]
    ax = fig.add_subplot(gs[1, 0]); label(ax, 'C', x=-0.18, y=1.04)
    P['fig6C'] = counts_heatmap(ax, rows, ['n_dual', 'n_A_only', 'n_B_only'],
                                'BindingDB compounds (gate: n ≥ 20 / class)', 20)
    ax = fig.add_subplot(gs[1, 1]); label(ax, 'D', x=-0.18, y=1.04)
    P['fig6D'] = counts_heatmap(ax, rows, ['n_sources_dual', 'n_sources_A_only', 'n_sources_B_only'],
                                'BindingDB sources (gate: ≥ 3 / class)', 3)
    fig.text(.57, .02, '0/8 pairs met both compound-count and independent-source criteria.',
             ha='center', fontsize=6.5)
    fig.subplots_adjust(left=.16, right=.98, top=.94, bottom=.08)
    save(fig, 'Fig6_evidence_boundary')


def fig_s1_protocol(D):
    fig, axs = plt.subplots(1, 2, figsize=(7, 2.85))
    y48 = float(D['pm110'][('PM48', 'vina')]['summary_min'])
    y110 = float(D['pm110'][('PM110', 'vina')]['summary_min'])
    e16 = v.primary_row(D, 'PIK3CA/mTOR')['smin']
    e8 = e8_value()
    ax = axs[0]; label(ax, 'A', x=-0.18, y=1.08)
    ax.plot([0, 1], [y48, y110], '-', color=C['vina'], lw=1.0)
    ax.plot(0, y48, 'o', color=C['vina'], ms=6)
    ax.plot(1, y110, 's', color=C['holdout'], ms=5.6)
    ax.axhline(.5, color=C['chance'], ls='--', lw=.7)
    ax.set_xticks([0, 1], ['PM48', 'PM110'], fontsize=7)
    ax.set(xlim=(-.35, 1.35), ylim=(.45, .82), ylabel=SMIN)
    ax.set_title('PIK3CA/mTOR panel size', fontsize=8)
    ax.legend(handles=[L(C['vina'], 'o', 'PM48', ms=5.4), L(C['holdout'], 's', 'PM110', ms=5.2)],
              loc='lower right', fontsize=6.0)
    ax = axs[1]; label(ax, 'B', x=-0.18, y=1.08)
    ax.plot([0, 1], [e16, e8], '-', color=C['vina'], lw=1.0)
    ax.plot(0, e16, 'o', color=C['vina'], ms=6)
    ax.plot(1, e8, 's', color=C['gnina'], ms=5.6)
    ax.axhline(.5, color=C['chance'], ls='--', lw=.7)
    ax.set_xticks([0, 1], ['E = 16', 'E = 8'], fontsize=7)
    ax.set(xlim=(-.35, 1.35), ylim=(.45, .82), ylabel=SMIN)
    ax.set_title('Search exhaustiveness', fontsize=8)
    ax.legend(handles=[L(C['vina'], 'o', 'E = 16', ms=5.4), L(C['gnina'], 's', 'E = 8', ms=5.2)],
              loc='lower right', fontsize=6.0)
    P['figS1'] = {'PM48': y48, 'PM110': y110, 'E16': e16, 'E8': e8}
    # Former main-text keys retained for any downstream numeric lock.
    P['fig6B_protocol'] = {'PM48': y48, 'PM110': y110}
    fig.subplots_adjust(left=.10, right=.98, wspace=.38, top=.82, bottom=.18)
    save(fig, 'FigS1_protocol_sensitivity')


def cognate_rmsd_rows():
    """Assemble 14 main-receptor RMSD markers from frozen CSVs only."""
    rank = read('data/jcim_novelty_v0/tables/cognate_rank_rmsd_reaudit_v1.csv')
    layer = read('data/jcim_chembl_universe_v0/local_track_b_v0/tables/layer3_cognate_rmsd_v1.csv')
    pm = read('data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/tables/pm48_01_rmsd_E16.csv')
    by_pdb = {}
    for r in rank:
        pdb = r['pdb']
        if pdb in by_pdb:
            continue
        by_pdb[pdb] = {
            'protein': r['target'], 'pdb': pdb,
            'top1': float(r['best_top1_A']),
            'top3': float(r['best_top3_A']),
            'best9': float(r['best_all_deposited_A']),
        }
    for r in layer:
        mode = int(float(r['best_mode']))
        rec = {
            'protein': r['protein'], 'pdb': r['pdb'],
            'top1': float(r['top1_rmsd']),
            'best9': float(r['best_of_9_rmsd']),
            'top3': float(r['best_of_9_rmsd']) if mode <= 3 else None,
        }
        by_pdb[r['pdb']] = rec
    seed_row = {'4L23': None, '4JT6': None}
    for r in pm:
        if r['seed'] != '20260727':
            continue
        pdb = r['target']
        mode = int(float(r['best_of_9_mode']))
        seed_row[pdb] = {
            'protein': 'PIK3CA' if pdb == '4L23' else 'mTOR', 'pdb': pdb,
            'top1': float(r['rmsd_mode1']),
            'best9': float(r['rmsd_best_of_9']),
            'top3': float(r['rmsd_best_of_9']) if mode <= 3 else None,
        }
    by_pdb.update({k: v for k, v in seed_row.items() if v})
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
        rec['label'] = f'{protein} {pdb}'
        rows.append(rec)
    return rows


def fig_s12_cognate():
    rows = cognate_rmsd_rows()
    xmax = 3.5
    fig, ax = plt.subplots(figsize=(7, 6.35))
    offscale = []
    for i, r in enumerate(rows):
        pts = [('top1', r['top1'], i - .18, 'o', C['vina']),
               ('top3', r['top3'], i, 's', C['desc']),
               ('best9', r['best9'], i + .18, 'D', C['a_only'])]
        for name, val, yy, m, col in pts:
            if val is None:
                continue
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
                ax.text(xmax + 0.26, yy, f'{val:.3f}', va='center', ha='left', fontsize=6.0, color=col)
    ax.axvline(2.0, color=C['chance'], ls='--', lw=.9)
    ax.set_yticks(range(len(rows)), [r['label'] for r in rows], fontsize=7)
    ax.set_ylim(len(rows) - .35, -.65)
    ax.set_xlabel(r'Heavy-atom RMSD ($\mathrm{\AA}$)')
    ax.set_xlim(-0.12, 4.35)
    ax.set_title('Cognate redocking of the 14 primary receptors', fontsize=8, pad=4)
    ax.legend(handles=[
        L(C['vina'], 'o', 'top-1', ms=4.8),
        L(C['desc'], 's', 'top-3 minimum', ms=4.6),
        L(C['a_only'], 'D', 'best-of-9 minimum', ms=4.4),
    ], loc='upper center', bbox_to_anchor=(.5, -.11), ncol=3, fontsize=6.4)
    P['figS12'] = rows
    P['figS12_offscale'] = offscale
    fig.subplots_adjust(left=.18, right=.90, top=.93, bottom=.14)
    save(fig, 'FigS12_cognate_rmsd')


def toc_graphic():
    fig = plt.figure(figsize=(3.25, 1.75), dpi=300)
    ax = fig.add_axes([0.02, 0.06, 0.96, 0.88])
    ax.set_xlim(0, 32); ax.set_ylim(0, 16); ax.axis('off')
    ax.text(5.2, 14.7, 'Four states', ha='center', fontsize=6.4, fontweight='bold')
    cells = [(1.15, 7.6, 'dual', C['dual']), (5.55, 7.6, 'A-only', C['a_only']),
             (1.15, 3.15, 'B-only', C['b_only']), (5.55, 3.15, 'neither', C['neither'])]
    for x, y, name, col in cells:
        ax.add_patch(Rectangle((x, y), 4.15, 3.7, fc=col, alpha=.22, ec=col, lw=.5))
        ax.text(x + 2.07, y + 1.85, name, ha='center', va='center', fontsize=5.8)
    ax.text(5.2, 1.35, r'B $\geq\theta$          B $<\theta$', ha='center', fontsize=5.2, color='#555555')
    ax.text(11.0, 8.4, '→', ha='center', va='center', fontsize=11, color='#888888')
    ax.text(17.4, 14.7, 'Same docking results', ha='center', fontsize=6.4, fontweight='bold')
    ax.add_patch(FancyBboxPatch((11.7, 8.55), 11.4, 4.55, boxstyle='round,pad=0.06,rounding_size=0.18',
                                fc='#F7F7F7', ec=C['vina'], lw=.7))
    ax.text(17.4, 11.85, 'dual — A-only  |  score B', ha='center', fontsize=5.8)
    ax.text(17.4, 10.05, 'dual — B-only  |  score A', ha='center', fontsize=5.8)
    ax.add_patch(FancyBboxPatch((11.7, 3.15), 11.4, 4.55, boxstyle='round,pad=0.06,rounding_size=0.18',
                                fc='#F7F7F7', ec=C['desc'], lw=.7))
    ax.text(17.4, 5.45, 'dual — neither', ha='center', fontsize=5.8)
    ax.text(23.2, 8.4, '→', ha='center', va='center', fontsize=11, color='#888888')
    ax.text(27.7, 14.7, 'Interpretation', ha='center', fontsize=6.4, fontweight='bold')
    ax.add_patch(FancyBboxPatch((23.7, 3.15), 7.9, 9.95, boxstyle='round,pad=0.06,rounding_size=0.18',
                                fc='#F4F7FA', ec='#CCCCCC', lw=.7))
    ax.text(27.65, 10.55, 'ligand-only', ha='center', fontsize=6.0, fontweight='bold')
    ax.text(27.65, 8.55, 'control', ha='center', fontsize=6.0)
    ax.text(27.65, 6.15, 'pocket', ha='center', fontsize=6.0, fontweight='bold')
    ax.text(27.65, 4.25, 'correspondence', ha='center', fontsize=6.0)
    save(fig, 'TOC_graphic', toc=True)


def supplements(D):
    fig_s1_protocol(D)
    v.fig_s4_forest(D); v.fig_s5_holdout(D)
    sim = [r for r in read('data/jcim_novelty_v0/tables/detectable_effect_simulation_v1.csv') if r['contrast'] == 'summary_min']
    pp = [p for p in PAIRS if any(r['pair'] == p for r in sim)]
    grid = sorted({float(r['true_auroc']) for r in sim})
    vals = [[float(next(r for r in sim if r['pair'] == p and float(r['true_auroc']) == a)['p_ci_excludes_0p5']) for a in grid] for p in pp]
    fig, ax = plt.subplots(figsize=(7, 2.9)); im = ax.imshow(vals, cmap='YlGnBu', vmin=0, vmax=1, aspect='auto')
    for i, row in enumerate(vals):
        for j, z in enumerate(row):
            ax.text(j, i, f'{z:.3f}', ha='center', va='center', fontsize=7, color='white' if z > .55 else C['ink'])
    ax.set_yticks(range(len(pp)), pp)
    ax.set_xticks(range(len(grid)), [f'{g:.2f}' for g in grid])
    ax.set_xlabel('True AUROC on both directional arms')
    fig.colorbar(im, ax=ax).set_label('P(95% CI excludes 0.5)')
    ax.set_title('Detectable-effect simulation: available three-pair results', fontsize=8)
    fig.subplots_adjust(left=.18, right=.95, bottom=.23, top=.85)
    P['figS6'] = sim
    save(fig, 'FigS6_detectable_effect')
    rows = [next(r for r in D['native'] if r['pair'] == p) for p in PAIRS]
    fig, axs = plt.subplots(1, 2, figsize=(7, 4.2))
    for ax, cols, title, gate, letter in [
        (axs[0], ['n_dual', 'n_A_only', 'n_B_only'], 'Compounds (criterion: n≥20 per class)', 20, 'A'),
        (axs[1], ['n_sources_dual', 'n_sources_A_only', 'n_sources_B_only'], 'Sources (criterion: ≥3 per class)', 3, 'B'),
    ]:
        label(ax, letter); counts_heatmap(ax, rows, cols, title, gate)
    fig.subplots_adjust(left=.16, right=.98, wspace=.55, bottom=.14, top=.87)
    P['figS8'] = rows
    save(fig, 'FigS8_bindingdb_native_slice')
    clusters = read('data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv')
    fig, axs = plt.subplots(1, 2, figsize=(7, 3), sharex=True)
    for ax, p, letter in zip(axs, ['EGFR/HER2', 'JAK1/TYK2'], 'AB'):
        label(ax, letter)
        r = v.s34_row(D, p); dotci(ax, r['delta'], r['lo'], r['hi'], 0, C['vina'], 'o')
        for i, (key, col, m) in enumerate([('scaffold_cluster', C['desc'], 's'), ('document_cluster', C['a_only'], 'D')], 1):
            rr = next(z for z in clusters if z['pair'] == p and z['estimator'] == key)
            dotci(ax, float(rr['delta_point']), float(rr['delta_ci_lo']), float(rr['delta_ci_hi']), i, col, m)
        ax.axvline(0, color=C['chance'], ls='--', lw=.7)
        ax.set(xlim=(-.12, .75), ylim=(2.5, -.5), xlabel=r'$\Delta$AUROC, target A score')
        ax.set_yticks(range(3), ['Ligand', 'Scaffold cluster', 'Document cluster'], fontsize=6.5)
        ax.set_title(p, fontsize=8)
        ax.legend(handles=[
            L(C['vina'], 'o', 'ligand', ms=4.4),
            L(C['desc'], 's', 'scaffold', ms=4.2),
            L(C['a_only'], 'D', 'document', ms=4.2),
        ], loc='upper center', bbox_to_anchor=(.5, -.28), ncol=3, fontsize=5.6)
    fig.subplots_adjust(left=.18, right=.97, wspace=.55, top=.83, bottom=.28)
    P['figS11'] = clusters
    save(fig, 'FigS11_cluster_uncertainty')
    top = next(r for r in read('data/jcim_novelty_v0/tables/mixed_library_enrichment_v1.csv')
               if r['pair'] == 'EGFR/HER2' and r['score'] == 'vina_mean' and r['cutoff'] == 'Top10')
    filt = next(r for r in read('data/jcim_novelty_v0/tables/and_filter_operating_point_v1.csv')
                if r['pair'] == 'EGFR/HER2' and r['score'] == 'vina_worst' and r['dual_percentile'] == '50')
    fig, axs = plt.subplots(1, 2, figsize=(7, 2.9))
    for ax, vals, title, letter in [
        (axs[0], [int(top[k]) for k in ['n_dual_top', 'n_A_only_top', 'n_B_only_top', 'n_neither_top']],
         'Top-10 by mean Vina score (library n=110)', 'A'),
        (axs[1], [int(filt[k]) for k in ['n_dual_pass', 'n_A_only_pass', 'n_B_only_pass']],
         'Dual-median AND filter (library n=98)', 'B'),
    ]:
        label(ax, letter)
        names = ['Dual', 'A-only', 'B-only', 'Neither'][:len(vals)]
        ax.bar(range(len(vals)), vals, color=[C['dual'], C['a_only'], C['b_only'], C['neither']][:len(vals)], width=.6)
        for i, z in enumerate(vals):
            ax.text(i, z + .35, str(z), ha='center', fontsize=7)
        ax.set_xticks(range(len(vals)), names, fontsize=7)
        ax.set_ylabel('Retained compounds')
        ax.set_ylim(0, max(vals) * 1.22 + 1)
        ax.set_title(title, fontsize=7.5)
    fig.subplots_adjust(left=.10, right=.97, wspace=.35, top=.80, bottom=.20)
    P['figS7'] = {'top10': top, 'and_filter': filt}
    save(fig, 'FigS7_posthoc_diagnostics')
    fig_s12_cognate()


def main():
    global SOURCE
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', type=Path, default=SNAP)
    args = parser.parse_args()
    SOURCE = args.source_root.resolve()
    style.OUT = OUT
    style.apply_style()
    v.ROOT = SOURCE
    v.DATA = SOURCE / 'data'
    v.OUT = OUT
    v.save_all = save
    v._read = lambda path: read(str(path.relative_to(SOURCE)).replace('\\', '/'))
    D = v.load()
    assert set(D['theta6']) == set(v.ORIGINAL_THREE)
    assert set(D['native'][i]['pair'] for i in range(len(D['native']))) == set(PAIRS)
    assert all(float(r['packaged_as_external_evaluation']) == 0 for r in D['native'])
    P['pair_order'] = list(PAIRS)
    fig1(D); fig1_c_data(D); fig2(D); fig3(D); fig4(D); fig5(D); fig6(D); supplements(D); toc_graphic()
    P.update({k: val for k, val in v.PROVENANCE['plotted'].items() if k not in P})
    P['primary'] = {p: v.primary_row(D, p) for p in PAIRS}
    assert abs(P['fig3B_max_abs'] - .023) < .001
    for rel in READS:
        dest = snapshot_path(rel)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if SOURCE != SNAP.resolve():
            shutil.copyfile(SOURCE / rel, dest)
    for rel in ['docs/MANUSCRIPT_JCIM_ZH.md', 'docs/MANUSCRIPT_JCIM_EN.md', 'docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md']:
        path = source_path(rel)
        READS[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
        dest = snapshot_path(rel)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if SOURCE != SNAP.resolve():
            shutil.copyfile(path, dest)
    audit = {
        'commit': SHA, 'source': 'https://github.com/1280602962-debug/gwj260531/pull/32',
        'inputs_sha256': READS,
        'input_files': {rel: str(snapshot_path(rel).relative_to(OUT)) for rel in READS},
        'generated': GENERATED, 'plotted': P,
    }
    (OUT / 'plotted_values.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    for stem in GENERATED:
        for ext in ['png', 'tif']:
            with Image.open(OUT / (stem + '.' + ext)) as im:
                assert im.mode == 'RGB'
                assert abs(im.info['dpi'][0] - 300) < 1
                assert im.width <= 2101 and im.height <= 2751
    print(f'PASS: {len(GENERATED)} figures; {len(READS)} pinned inputs; RGB / 300 dpi / size checks')


if __name__ == '__main__':
    main()
