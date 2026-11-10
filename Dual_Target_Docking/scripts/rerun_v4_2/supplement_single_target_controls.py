"""Reproduce the exploratory single-target controls in manuscript v6 (Table S11).

This is a new archived reproduction, not recovery of an original bootstrap file.
Only frozen S1 scores are read. No docking, model training, or primary metric
file is changed. Run this script from any directory; --root is optional.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import platform
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

PAIRS = ['EGFR/HER2', 'JAK1/JAK2', 'JAK1/TYK2', 'PIK3CA/mTOR',
         'AChE/BChE', 'F2/F10', 'PPARG/PPARA', 'PPARA/PPARD']
B, SEED = 10000, 271828


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def auc(pos, neg):
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + .5 * (diff == 0).sum()) / diff.size)


PROTOCOL = ('S1 is the frozen higher-is-better Vina score (negative docking energy; '
        'search seed 42). Positive class: the single-target class active at the '
        'scored target. Panel A uses neither negatives; panel B uses the opposite '
        'single-target class. Dual ligands are excluded. Counts include only '
        'activity-eligible observations with a valid S1 score at the scored target. '
        'CI: 95% class-stratified percentile bootstrap, 10,000 resamples, '
        'NumPy default_rng seed 271828 reset for each pair and score target; '
        'linear 2.5th and 97.5th percentiles. Within each replicate, positive, '
        'opposite-single-target, then neither indices are drawn; the same positive '
        'resample is shared between panels. Input row order is retained. '
        'Ties count as 0.5. These exploratory intervals are unadjusted for '
        'multiple comparisons and are not simultaneous intervals. They were '
        'reproduced from frozen scores on 10 November 2026; an original archived '
        'bootstrap output for this comparison was not recovered.')

NOTE = ('S1: negative Vina docking energy, with higher values favorable (search seed 42). '
        'The positive class is the single-target class active at the scored target. '
        'Neither ligands are negative in panel A; the opposite single-target class '
        'is negative in panel B. Dual ligands are excluded. Counts include eligible '
        'observations with a valid score at that target. CI: 95% class-stratified '
        'percentile bootstrap, 10,000 resamples, seed 271828; linear 2.5th and 97.5th '
        'percentiles. Positive resamples are shared between the two comparisons. '
        'Score ties count as 0.5. Intervals are exploratory and unadjusted for multiplicity.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    source = root / 'results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv'
    before = sha(source)
    out = root / 'analysis/single_target_controls_20261110'
    out.mkdir(parents=True, exist_ok=True)
    original = list(csv.DictReader(source.open(newline='', encoding='utf-8')))
    rows = [r for r in original if r['activity_eligible'] == '1'
            and r['analysis_set'] in ('main', 'q110_primary')]
    assert len(original) == 869 and len(rows) == 867
    metrics, members, replicas, checks = [], [], {}, []
    for pair in PAIRS:
        for side, positive, opposite in [('A', 'A_only', 'B_only'), ('B', 'B_only', 'A_only')]:
            target = pair.split('/')[0 if side == 'A' else 1]
            key = f'M0_score_{side}'
            groups = []
            for cls in (positive, opposite, 'neither'):
                eligible = [r for r in rows if r['pair'] == pair and r['class'] == cls]
                kept = [r for r in eligible if r[f'M0_{side}_valid'] == '1']
                assert len({r['canonical_ligand_id'] for r in kept}) == len(kept)
                values = np.array([float(r[key]) for r in kept])
                assert len(values) and np.isfinite(values).all()
                groups.append(values)
                order = 0
                for r in eligible:
                    valid = r in kept
                    if valid:
                        order += 1
                    members.append(dict(pair=pair, score_side=side, score_target=target,
                        positive_class=positive, experimental_class=cls,
                        role='positive' if cls == positive else 'negative',
                        canonical_ligand_id=r['canonical_ligand_id'],
                        global_ligand_entity_id=r['global_ligand_entity_id'],
                        analysis_set=r['analysis_set'], included=int(valid),
                        bootstrap_group_order=order if valid else '',
                        score_higher_better=r[key] if valid else '',
                        exclusion_reason='' if valid else r[f'M0_{side}_missing']))
            pos, single, neither = groups
            rng = np.random.default_rng(SEED)
            reps = np.empty((B, 2), dtype=float)
            for rep in range(B):
                pp = pos[rng.integers(len(pos), size=len(pos))]
                ss = single[rng.integers(len(single), size=len(single))]
                nn = neither[rng.integers(len(neither), size=len(neither))]
                reps[rep] = auc(pp, ss), auc(pp, nn)
            replicas[f'{PAIRS.index(pair)+1:02d}_{side}'] = reps
            for panel, negcls, neg, col in [('A', 'neither', neither, 1), ('B', opposite, single, 0)]:
                point = auc(pos, neg)
                y = [1] * len(pos) + [0] * len(neg)
                sklearn_point = roc_auc_score(y, np.r_[pos, neg])
                assert abs(point - sklearn_point) < 1e-14
                assert abs(point + auc(neg, pos) - 1) < 1e-14
                lo, hi = np.quantile(reps[:, col], [.025, .975], method='linear')
                assert np.isfinite(reps[:, col]).all() and 0 <= lo <= hi <= 1
                metrics.append(dict(panel=panel, pair=pair, score_side=side,
                    score_target=target, positive_class=positive, negative_class=negcls,
                    n_positive=len(pos), n_negative=len(neg), auroc=point,
                    ci_lower=float(lo), ci_upper=float(hi), bootstrap_B=B,
                    bootstrap_seed=SEED, ci_method='class_stratified_percentile_linear',
                    ci_relative_to_0_5='above' if lo > .5 else 'below' if hi < .5 else 'includes',
                    analysis_status='exploratory_reproduced_20261110'))
                checks.append(dict(pair=pair, score_side=side, panel=panel,
                    auroc_sklearn_difference=point-sklearn_point, reverse_sum=point+auc(neg, pos),
                    finite_bootstrap_replicates=int(np.isfinite(reps[:, col]).sum())))
    metrics.sort(key=lambda r: (r['panel'], PAIRS.index(r['pair']), r['score_side']))
    a, b = ([r for r in metrics if r['panel'] == p] for p in ('A', 'B'))
    assert sum(r['auroc'] > .5 for r in a) == 16
    assert sum(r['ci_lower'] > .5 for r in a) == 4
    assert sum(r['auroc'] > .5 for r in b) == 11
    assert [(r['pair'], r['score_target']) for r in b if r['ci_upper'] < .5] == [
        ('EGFR/HER2', 'EGFR'), ('JAK1/TYK2', 'JAK1')]
    assert abs(b[0]['auroc'] - .18074324324324326) < 1e-14
    assert abs(b[4]['auroc'] - .3125) < 1e-14
    assert sha(source) == before
    write_csv(out/'Table_S11_metrics.csv', metrics)
    write_csv(out/'Table_S11_membership.csv', members)
    write_csv(out/'Table_S11_QA.csv', checks)
    np.savez_compressed(out/'Table_S11_bootstrap_replicates.npz', **replicas)
    config = dict(date='2026-11-10', source_commit='2017a5d7f9815715f5203a7b4e413c9719e9e032',
        source_path=str(source.relative_to(root)), source_sha256_before=before,
        source_sha256_after=sha(source), source_rows=869, activity_eligible_rows=867,
        Python=platform.python_version(), NumPy=np.__version__,
        bootstrap_B=B, bootstrap_seed=SEED, score_column='M0_score_A or M0_score_B',
        scoring_scheme='S1 Vina higher-better, docking seed 42',
        replication_status='new_reproduction_of_exploratory_control_from_frozen_scores',
        original_bootstrap_output_recovered=False, protocol=PROTOCOL,
        npz_columns=['opposite_single_target_AUROC', 'neither_AUROC'],
        npz_keys='01..08 in fixed pair order, suffix A/B = score side',
        summary=dict(panel_A_point_above_0_5=16, panel_A_ci_above_0_5=4,
                     panel_B_point_above_0_5=11, panel_B_ci_below_0_5=2),
        script_sha256=sha(Path(__file__)))
    (out/'Table_S11_config.json').write_text(json.dumps(config, indent=2, ensure_ascii=False)+'\n')
    title = 'Table S11 Exploratory single-target controls using frozen Vina scores'
    md, sections = ['# '+title, ''], []
    headers = ['Target pair', 'Score target', 'Positive class', 'n positive', 'n negative', 'AUROC (95% CI)']
    for panel, rr, subtitle in [('A', a, 'Single-target-active versus neither'),
                               ('B', b, 'Single-target-active versus opposite single-target-active')]:
        display = [[r['pair'], r['score_target'], r['positive_class'].replace('_','-'),
                    str(r['n_positive']), str(r['n_negative']),
                    f"{r['auroc']:.4f} [{r['ci_lower']:.4f}, {r['ci_upper']:.4f}]"] for r in rr]
        md += ['## '+panel+'. '+subtitle, '', '| '+' | '.join(headers)+' |',
               '| '+' | '.join(['---']*len(headers))+' |']
        md += ['| '+' | '.join(row)+' |' for row in display]
        md += ['']
        sections.append('<section id="table-s11-'+panel.lower()+'"><h2>'+title+'</h2><h3>'+panel+'. '+subtitle+'</h3><table><thead><tr>'+
            ''.join('<th>'+html.escape(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+
            ''.join('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in row)+'</tr>' for row in display)+
            '</tbody></table><p class="note">'+html.escape(NOTE)+'</p></section>')
    md += [NOTE, '']
    (out/'TableS11.md').write_text('\n'.join(md), encoding='utf-8')
    (out/'Table_S11_sections.html').write_text(''.join(sections), encoding='utf-8')
    print(json.dumps(config['summary']))
    for r in b:
        if r['ci_upper'] < .5:
            print(r['pair'],r['score_target'],r['auroc'],r['ci_lower'],r['ci_upper'])


if __name__ == '__main__':
    main()
