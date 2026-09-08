#!/usr/bin/env python3
"""Build the manuscript evidence layer from immutable deposited results.

No docking, model fitting, nomination or MD is performed. Recomputed descriptors
are annotations; historical membership is checked rather than silently changed.
Run from any directory. Inputs/outputs and hashes are recorded for verification.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import re
from importlib.metadata import version
import subprocess
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
C1 = 'data/campaigns/c1/'
C5 = 'data/campaigns/c5/'
OUT = ROOT / 'data/manuscript'
INPUTS: dict[str, dict] = {}
OUTPUTS: set[str] = set()
CHECKS: list[dict] = []


def digest(path):
    # Git on Windows checks text out as CRLF. Hash text in canonical LF form so
    # the same committed evidence verifies on Windows, Linux and GitHub.
    if Path(path).suffix.lower() in {'.csv','.tsv','.json','.yaml','.yml','.md','.py','.sh','.txt','.sdf','.pdbqt','.pdb','.cif','.smi','.log'}:
        return hashlib.sha256(Path(path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def source(rel):
    p = ROOT / rel
    INPUTS[rel] = {'path': rel, 'sha256': digest(p), 'bytes': p.stat().st_size}
    return p


def frame(rel):
    return pd.read_csv(source(rel))


def obj(rel):
    return json.loads(source(rel).read_text(encoding='utf-8-sig'))


def write(rel, value):
    p = ROOT / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(value, pd.DataFrame):
        value.to_csv(p, index=False, lineterminator='\n')
    elif p.suffix == '.yaml':
        p.write_text(yaml.safe_dump(value, allow_unicode=True, sort_keys=False), encoding='utf-8')
    elif isinstance(value, (list, dict)):
        p.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    else:
        p.write_text(value, encoding='utf-8')
    OUTPUTS.add(rel)


def check(name, passed, detail=''):
    CHECKS.append({'check': name, 'passed': bool(passed), 'detail': detail})


def truth(v):
    return str(v).strip().lower() in {'true', '1', '1.0'}


def ids(df, key='ligand_id'):
    return set(df[key].astype(str))


def metrics(tp, fn, fp, tn):
    sens, spec = tp/(tp+fn), tn/(fp+tn)
    def wilson(k, n):
        z = 1.959963984540054
        p = k/n
        c = (p+z*z/(2*n))/(1+z*z/n)
        r = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
        return c-r, c+r
    sl, sh = wilson(tp, tp+fn)
    tl, th = wilson(tn, fp+tn)
    odds, p = fisher_exact([[tp, fn], [fp, tn]])
    return dict(tp=tp, fn=fn, fp=fp, tn=tn, sensitivity=sens,
                specificity=spec, lr_plus=sens/(1-spec) if spec < 1 else None,
                sensitivity_wilson95_low=sl, sensitivity_wilson95_high=sh,
                specificity_wilson95_low=tl, specificity_wilson95_high=th,
                odds_ratio_mle=float(odds) if math.isfinite(odds) else None,
                odds_ratio_haldane=(tp+.5)*(tn+.5)/((fn+.5)*(fp+.5)),
                fisher_exact_p=float(p))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inventory', action='store_true', help='Hash every file in the audited base commit')
    args = ap.parse_args()
    policy = yaml.safe_load(source('config/campaign_final.yaml').read_text(encoding='utf-8'))
    base = policy['audited_base_commit']
    records = frame('data/processed/nlrp3_records.csv')
    library = frame('data/repurposing/repurposing_manifest.csv')
    pool = frame('data/repurposing/screening/docking_pool_p05.csv')
    acids = frame(C1+'07_clinical_dock/acid_pool/acid_clinical_pool.csv')
    soft = frame(C1+'07_clinical_dock/acid_pool/acid_clinical_pool_chemistry_pass.csv')
    prim = frame(C5+'04_shortlist_frozen/primary_candidates.csv')
    backup = frame(C5+'04_shortlist_frozen/backup_candidates.csv')
    summary = obj(C5+'04_shortlist_frozen/shortlist_freeze_summary.json')
    train = obj('data/models/training_report.json')
    screen = obj('data/repurposing/screening/nlrp3_screening_summary_clinical_all.json')
    sens = obj('data/si/nlrp3_threshold_sensitivity/summary.json')
    for label, df, n, key in [('library', library, 8319, 'repurposing_id'),
                             ('pool', pool, 1588, 'repurposing_id'), ('acids', acids, 303, 'repurposing_id'),
                             ('soft', soft, 156, 'repurposing_id'), ('primary', prim, 12, 'ligand_id'),
                             ('reserve', backup, 21, 'ligand_id')]:
        check(label+'_count', len(df)==n, str(len(df)))
        check(label+'_unique_ids', df[key].is_unique)
    check('pool_nested_in_library', ids(pool, 'repurposing_id') <= ids(library, 'repurposing_id'))
    check('acids_nested_in_pool', ids(acids, 'repurposing_id') <= ids(pool, 'repurposing_id'))
    check('soft_is_exact_soft_gate_subset', ids(soft, 'repurposing_id') == ids(acids[acids.pass_chemistry_soft.map(truth)], 'repurposing_id'))
    check('pool_scores_ge_0.5', (pool.p_active_nlrp3 >= .5).all())
    check('primary_reserve_disjoint', not ids(prim) & ids(backup))
    check('candidate_ids_in_docked_pool', (ids(prim)|ids(backup)) <= ids(soft, 'repurposing_id'))
    check('training_records', len(records)==609)
    check('training_molecules', records.canonical_smiles.nunique()==513)
    check('training_assays', records['Assay ChEMBL ID'].nunique()==39)

    # Explicit provenance, preserving heterogeneous source descriptions/metadata.
    assay_rows = []
    def joincol(g, c):
        return ' | '.join(sorted(set(g[c].dropna().astype(str)))) if c in g else ''
    for aid, g in records.groupby('Assay ChEMBL ID', sort=True):
        description = joincol(g, 'Assay Description')
        text_lower=description.lower()
        inferred_cell='not stated'
        for term,label in [('whole blood','whole blood'),('thp','THP-1 / differentiated derivative'),
                           ('pbmc','PBMC'),('bone marrow','bone-marrow-derived macrophages'),
                           ('monocyte derived','monocyte-derived macrophages'),('glial','primary glial cells')]:
            if term in text_lower:
                inferred_cell=label
                break
        assay_rows.append(dict(assay_id=aid, assay_description=description,
            system=joincol(g, 'BAO Label'), organism=joincol(g, 'Assay Organism'),
            cell_type=joincol(g, 'Assay Cell Type'), inferred_cell_type_source=joincol(g, 'assay_cell_type'),
            cell_type_from_description=inferred_cell,
            system_from_description='ex_vivo_whole_blood' if 'whole blood' in text_lower else 'cell_based_functional',
            tissue=joincol(g, 'Assay Tissue Name'), assay_type=joincol(g, 'Assay Type'),
            endpoint='IL-1beta-related readout; see verbatim assay description',
            standard_type=joincol(g, 'Standard Type'), n_records=len(g),
            n_molecules=g.canonical_smiles.nunique(), evidence_class='NLRP3-related activity; not direct binding',
            classification_basis='curated endpoint and source description; ChEMBL assay type B alone is not binding proof',
            missing_system_or_cell=(not joincol(g, 'BAO Label') or not joincol(g, 'Assay Cell Type')),
            source_file='data/processed/nlrp3_records.csv', source_url=f'https://www.ebi.ac.uk/chembl/explore/assay/{aid}'))
    write('data/manuscript/si/nlrp3_assay_provenance.csv', pd.DataFrame(assay_rows))
    write('data/manuscript/si/nlrp3_threshold_sensitivity.json', sens)
    write('data/manuscript/tables/model_validation.json', {
        'nlrp3_activity_classifier': train['nlrp3']['cv_metrics'],
        'n_unique_assays_in_records':39, 'n_named_one_hot_assay_categories':25,'n_one_hot_dimensions_including_other':26,
        'n_production_prediction_contexts':len(screen['ensemble_assays']),
        'production_assay_contexts':screen['ensemble_assays'],
        'urat1_regression_SI_only':train['urat1']['cv_metrics'],
        'production_pool_is_frozen':True, 'stored_model_generation_link':'not established by a contemporaneous hash',
        'retraining_is_sensitivity_not_production_replay':True})
    model = source('data/models/nlrp3_model.joblib')
    write('data/frozen/model_manifest.json', {
        'path':'data/models/nlrp3_model.joblib', 'sha256':digest(model),
        'git_blob':'819d2b4ea9bab0197748381f4c5dc0bd82637d78',
        'first_commit':'72b132282e941a308dd7541d8e0d5711fc07d82e',
        'first_commit_utc':'2026-07-01T08:45:18Z',
        'unchanged_since_first_commit':True,
        'role':'deposited model artifact from the same commit as the 1588 pool; bit-identical regeneration of those scores is not independently demonstrated',
        'unpickled_in_this_audit':False,
        'pickle_strings':{'sklearn_version':'1.9.0','n_jobs':-1,'xgboost_package_version':None},
        'production_membership':'data/repurposing/screening/docking_pool_p05.csv',
        'pool_sha256':digest(source('data/repurposing/screening/docking_pool_p05.csv')),
        'env_search':'data/manuscript/audit/nlrp3_scoring_env_search.json'})

    # Reconstruct the actual set operations, retaining historical chemistry decisions.
    a1 = frame(C1+'07_clinical_dock/acid_dual_a1_frozen/acid_dual_keep_seed42.csv')
    a1ids = ids(a1[a1.keep_dual_acid_geometry.map(truth)])
    structural_count, loose_count = Counter(), Counter()
    seed_metrics = {}
    loose42 = set()
    for seed in [42,43,44]:
        st = frame(C1+f'07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{seed}.csv')
        structural_count.update(ids(st[st.keep_dual_acid_structural.map(truth)]))
        loose = frame(C1+f'07_clinical_dock/acid_dual_a2/acid_dual_keep_a2_seed{seed}.csv')
        loose_count.update(ids(loose[loose.keep_dual_acid_geometry.map(truth)]))
        if seed == 42:
            loose42 = ids(loose[loose.keep_dual_acid_geometry.map(truth)])
        seed_metrics[seed] = frame(C1+f'07_clinical_dock/acid_dual_a2/nlrp3_structural_metrics_seed{seed}.csv').set_index('ligand_id')
    structural = {k for k,v in structural_count.items() if v>=2}
    audited = frame(C1+'08_nomination/acid_a2_eligible_audited.csv')
    eligible = ids(audited) & structural
    t1, t2 = frame(C5+'03_tiering/tier1_candidates.csv'), frame(C5+'03_tiering/tier2_candidates.csv')
    check('tier1_reconstructed', ids(t1)==(eligible&a1ids))
    check('tier2_reconstructed', ids(t2)==(eligible-a1ids))
    for df, frozen, label in [(t1,prim,'primary'),(t2,backup,'reserve')]:
        accepted=df[~df.beta_lactam_flag.map(truth)&~df.is_structural_control_not_candidate.map(truth)]
        check(label+'_freeze_membership', ids(accepted)==ids(frozen))
    chemistry = yaml.safe_load(source('config/chemistry_final.yaml').read_text(encoding='utf-8'))
    filt = frame('data/repurposing/p2/filters_pool.csv')
    # Name joins were inherited; fail on conflicting annotations instead of picking one.
    conflict = filt.groupby('name')[['pains_any','brenk','nih']].nunique().max().max()
    check('legacy_alert_name_join_unambiguous', conflict<=1)
    filt = filt.drop_duplicates('name').set_index('name')
    patt = Chem.MolFromSmarts(chemistry['beta_lactam_smarts'])
    ledger=[]
    for _, r in soft.iterrows():
        lid, name = r.repurposing_id, r['name']
        f = filt.loc[name] if name in filt.index else None
        alerts_known=f is not None
        clean=alerts_known and not truth(f.pains_any) and not truth(f.brenk)
        matches=[x for x in chemistry['historical_name_exclusion_substrings'] if x in name.upper()]
        beta=Chem.MolFromSmiles(r.canonical_smiles).HasSubstructMatch(patt)
        old_eligible=(lid in loose42 and loose_count[lid]>=2 and clean and not matches)
        ledger.append(dict(ligand_id=lid,name=name,a1_seed42_dual_pass=lid in a1ids,
            urat1_pose_atom_matcher_in_scope=truth(r.has_carboxylate),
            structural_seed_passes=structural_count[lid],loose_seed_passes=loose_count[lid],
            stored_chemistry_eligible=lid in ids(audited),reconstructed_historical_eligible=old_eligible,
            loose_seed42_pass=lid in loose42,alerts_known=alerts_known,
            pains_any=truth(f.pains_any) if alerts_known else None,brenk=truth(f.brenk) if alerts_known else None,
            nih=truth(f.nih) if alerts_known else None,beta_lactam_flag=beta,
            historical_name_exclusion=';'.join(matches),structural_control=lid=='REP_07907',
            structural_only_counterfactual_eligible=(lid in structural and lid in loose42 and clean and not beta and lid!='REP_07907'),
            nomination_status='primary' if lid in ids(prim) else 'reserve_admissible' if lid in ids(backup) else 'not_nominated'))
    ledger=pd.DataFrame(ledger)
    check('historical_chemistry_40_reconstructed',
          ids(ledger[ledger.reconstructed_historical_eligible])==ids(audited))
    write('data/manuscript/tables/screening_decision_ledger.csv',ledger)
    delta=ledger[ledger.structural_only_counterfactual_eligible & ~ledger.ligand_id.isin(ids(prim)|ids(backup))]
    write('data/manuscript/si/structure_only_chemistry_counterfactual.csv',delta)

    # Descriptors are recomputed on deposited SMILES, never fed back into nomination.
    descriptors=[]
    acid_patterns={k:Chem.MolFromSmarts(v) for k,v in chemistry['acid_smarts'].items()}
    for _,r in acids.iterrows():
        m=Chem.MolFromSmiles(r.canonical_smiles)
        if m is None:
            raise ValueError(f'Invalid acid SMILES: {r.repurposing_id}')
        descriptors.append(dict(ligand_id=r.repurposing_id,name=r['name'],canonical_smiles=r.canonical_smiles,
            inchi_key_recomputed=Chem.MolToInchiKey(m),murcko_scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=m),
            mw_recomputed=Descriptors.MolWt(m),clogp=Crippen.MolLogP(m),tpsa=rdMolDescriptors.CalcTPSA(m),
            hbd=Lipinski.NumHDonors(m),hba=Lipinski.NumHAcceptors(m),rotatable_bonds=Lipinski.NumRotatableBonds(m),
            qed_recomputed=QED.qed(m),acid_types=';'.join(k for k,p in acid_patterns.items() if m.HasSubstructMatch(p)),
            chemistry_soft_pass=truth(r.pass_chemistry_soft),rdkit_version=rdBase.rdkitVersion))
    desc=pd.DataFrame(descriptors)
    write('data/manuscript/tables/acid_descriptors.csv',desc)
    chem_summary=[]
    for label, subset in [('acid_303',desc),('docked_156',desc[desc.chemistry_soft_pass]),('primary_12',desc[desc.ligand_id.isin(ids(prim))])]:
        row=dict(population=label,n=len(subset),n_murcko_scaffolds=subset.murcko_scaffold.nunique(),
                 n_acyclic=int(subset.murcko_scaffold.fillna('').eq('').sum()))
        for typ in acid_patterns:
            row['n_'+typ]=int(subset.acid_types.str.contains(typ,regex=False).sum())
        for col in ['mw_recomputed','clogp','tpsa','hbd','hba','rotatable_bonds','qed_recomputed']:
            for stat,value in [('min',subset[col].min()),('median',subset[col].median()),('max',subset[col].max())]:
                row[col+'_'+stat]=float(value)
        chem_summary.append(row)
    write('data/manuscript/tables/chemical_space_summary.csv',pd.DataFrame(chem_summary))
    prep=frame(C1+'01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv')
    write('data/manuscript/si/ligand_microstates.csv',prep)
    assay_ids=ids(records.rename(columns={'Molecule ChEMBL ID':'ligand_id'}))
    candidates=[]
    for role, df in [('primary',prim),('reserve_admissible',backup)]:
        c=df.drop(columns=['nlrp3_percentile'],errors='ignore').copy()
        c=c.merge(desc.drop(columns=['name','canonical_smiles']),on='ligand_id',validate='one_to_one')
        c=c.merge(library[['repurposing_id','inchi_key']].rename(columns={'repurposing_id':'ligand_id','inchi_key':'inchi_key_source'}),on='ligand_id',validate='one_to_one')
        c['nomination_status']=role
        c['structural_seed_passes']=c.ligand_id.map(structural_count)
        c['historical_n_seed_pass_meaning']='A2 dual loose geometry, not A1 or W2 consistency'
        for seed, dfm in seed_metrics.items():
            for col in ['keep_nlrp3_structural','selected_mode','pocket_overlap_frac','n_key_contacts','ifp_jaccard_vs_np3146']:
                c[f'nlrp3_{col}_seed{seed}']=c.ligand_id.map(dfm[col])
        check(role+'_not_exact_chembl_training_id', not set(c.chembl_id.dropna())&assay_ids)
        write('data/frozen/'+('primary_candidates_frozen.csv' if role=='primary' else 'backup_candidates_frozen.csv'),c)
        candidates.append(c)
    write('data/manuscript/tables/candidate_seed_metrics.csv',pd.concat(candidates,ignore_index=True))

    # Auditable 2x2 statistics; do not pool repeated seeds as independent ligands.
    w2=frame(C5+'02_urat1_ifp/w2_ifp_gate_per_mol.csv')
    old_ci=frame(C5+'05_presubmission_audit/gate_interval_metrics.csv')
    gates=[]
    for gate,col,oldname in [('A1','keep_a1','A1_arg_locked_7.7027'),('A2','keep_a2','A2_geometry_then_cnn'),('W2','keep_urat1_ifp','W2_IFP_cnn_top1')]:
        positive=w2.label.eq(1); passed=w2[col].map(truth)
        values=metrics(int((positive&passed).sum()),int((positive&~passed).sum()),int((~positive&passed).sum()),int((~positive&~passed).sum()))
        old=old_ci[old_ci.gate.eq(oldname)].iloc[0]
        check(gate+'_2x2_matches_deposit',all(values[k]==int(old[k]) for k in ['tp','fn','fp','tn']))
        values.update(gate=gate,seed=42,ci_method='Wilson binomial; deposited stratified bootstrap for LR/OR',
                      lr_boot95_low=float(old.lr_boot_lo),lr_boot95_high=float(old.lr_boot_hi),
                      or_haldane_boot95_low=float(old.or_boot_lo),or_haldane_boot95_high=float(old.or_boot_hi))
        gates.append(values)
    check('W2_subset_of_A1',not (w2.keep_urat1_ifp.map(truth)&~w2.keep_a1.map(truth)).any())
    write('data/manuscript/tables/urat1_gate_validation.csv',pd.DataFrame(gates))
    w4=frame(C5+'02_nlrp3_panel/w4_panel_metrics_all_seeds.csv')
    w4_rows=[]
    for seed in [42,43,44]:
        sub=w4[w4.seed.eq(seed)]
        pos=sub[sub.role.eq('positive')]
        dec=sub[sub.role.eq('decoy')]
        failed=dec.error.fillna('').ne('') if 'error' in dec else pd.Series(False,index=dec.index)
        for panel,bg in [('intent_to_score',dec),('complete_case',dec[~failed]),('clinical_background',sub[sub.role.eq('clinical_acid_background')])]:
            for gate,col in [('loose','keep_nlrp3_pose'),('structural','keep_nlrp3_structural')]:
                tp=int(pos[col].map(truth).sum()); fp=int(bg[col].map(truth).sum())
                w4_rows.append(dict(seed=seed,panel=panel,gate=gate,failed_decoys=int(failed.sum()),
                    **metrics(tp,len(pos)-tp,fp,len(bg)-fp)))
    w4_table=pd.DataFrame(w4_rows)
    write('data/manuscript/tables/nlrp3_gate_validation.csv',w4_table)
    w4source=obj(C5+'02_nlrp3_panel/w4_structural_gate_summary.json')
    s42=w4_table[(w4_table.seed==42)&(w4_table.panel=='intent_to_score')&(w4_table.gate=='structural')].iloc[0]
    check('W4_seed42_2x2_matches_deposit',all(int(s42[k])==w4source['primary']['structural'][k] for k in ['tp','fn','fp','tn']))

    l2=frame(C1+'02_selfdock/l2_selfdock_metrics.csv')
    redock=l2[['ligand_id','seed','target','rmsd_cnnscore_selected','rmsd_best_of_n']].rename(columns={
        'ligand_id':'ligand','rmsd_cnnscore_selected':'top1_pose_rmsd_A','rmsd_best_of_n':'best_of_9_pose_rmsd_A'})
    redock['source_file']=C1+'02_selfdock/l2_selfdock_metrics.csv'
    bz=frame(C5+'01_crossdock/gate_audit_summary_dual.csv')
    br=pd.DataFrame({'ligand':'benzbromarone','seed':bz.seed,'target':'urat1_9dka',
        'top1_pose_rmsd_A':bz.cnn_top1_rmsd_pose_rmsd,'best_of_9_pose_rmsd_A':bz.best9_pose_rmsd,
        'source_file':C5+'01_crossdock/gate_audit_summary_dual.csv'})
    redock=pd.concat([redock,br],ignore_index=True)
    redock['top1_le_2A']=redock.top1_pose_rmsd_A<=2
    redock['rmsd_definition']='symmetry-aware heavy-atom pose RMSD in receptor frame; no ligand superposition'
    write('data/manuscript/tables/urat1_redocking_summary.csv',redock[redock.target.str.startswith('urat1')])
    write('data/manuscript/tables/nlrp3_selfdocking_summary.csv',redock[redock.target.str.startswith('nlrp3')])
    maps=[]
    for target,path in [('URAT1',C5+'02_urat1_ifp/urat1_key_residues.json'),('NLRP3',C1+'01_ligand_prep/selfdock_refs/nlrp3_key_residues.json')]:
        j=obj(path)
        for label,r in j['residues'].items():
            maps.append(dict(target=target,pdb=j['pdb'],literature_label=label,
                prepared_chain=r['prep_chain'],prepared_residue=r['prep_resi'],prepared_resname=r['prep_resn'],status='mapped',source=path))
        for label,r in j.get('unmatched',{}).items():
            maps.append(dict(target=target,pdb=j['pdb'],literature_label=label,status='unmapped_sequence_mismatch',source=path))
    write('data/manuscript/tables/residue_mapping.csv',pd.DataFrame(maps))
    receptor_rows=[]
    for pdb,target,ligand,role in [('9DK9','URAT1','apo','optional_crossdock_not_completed'),
        ('9DKA','URAT1','benzbromarone_R75','selfdock_and_crystal_anchor_and_rigid_transfer'),
        ('9DKB','URAT1','lesinurad_A1AIL','production'),
        ('9DKC','URAT1','TD-3_A1A45','crystal_anchor_and_rigid_transfer'),
        ('9B1I','URAT1','verinurad','rigid_transfer_only'),
        ('7ALV','NLRP3','NP3-146_RM5','production'),('8ETR','NLRP3','GDC-2394_WTN','rigid_transfer_only')]:
        rel=f'data/structures/pdb/{pdb}.'+('pdb' if pdb=='7ALV' else 'cif')
        text=source(rel).read_text(encoding='utf-8')
        match=re.search(r'_em_3d_reconstruction.resolution\s+([\d.]+)',text)
        if not match:
            match=re.search(r'REMARK   2 RESOLUTION.\s+([\d.]+)',text)
        resolution=float(match.group(1)) if match else None
        receptor_rows.append(dict(pdb=pdb,target=target,ligand=ligand,resolution_A=resolution,
            experiment='X-RAY DIFFRACTION' if pdb=='7ALV' else 'ELECTRON MICROSCOPY',
            structural_state='apo' if pdb=='9DK9' else 'ligand_bound',role=role,
            enters_production_library_docking=role=='production',source=rel,source_sha256=INPUTS[rel]['sha256'],
            rcsb_url=f'https://www.rcsb.org/structure/{pdb}'))
    write('data/manuscript/tables/receptor_manifest.csv',pd.DataFrame(receptor_rows))

    # Export supportive audits without altering their historical outputs.
    for name in ['novelty_leakage_primary.csv','w4_positive_scaffolds.csv','arg_threshold_sensitivity_primary.csv',
                 'arg_threshold_sensitivity_benchmark.csv','receptor_sensitivity_urat1.csv','receptor_sensitivity_nlrp3.csv']:
        write('data/manuscript/si/'+name,frame(C5+'05_presubmission_audit/'+name))
    funnel=[('clinical_library',None,len(library),'clinical_all deposited library','data/repurposing/repurposing_manifest.csv','scripts/build_repurposing_library.py'),
        ('nlrp3_activity_pool',len(library),len(pool),'frozen model score >=0.5','data/repurposing/screening/docking_pool_p05.csv','scripts/screen_repurposing_library.py'),
        ('acid_equivalents',len(pool),len(acids),'deposited acid SMARTS',C1+'07_clinical_dock/acid_pool/acid_clinical_pool.csv','scripts/build_c1_acid_clinical_pool.py'),
        ('chemistry_soft',len(acids),len(soft),'MW 200-550; TPSA<=140; rotatable<=10; HBD<=5; HBA<=10',C1+'07_clinical_dock/acid_pool/acid_clinical_pool_chemistry_pass.csv','scripts/build_c1_acid_clinical_pool.py'),
        ('dual_structural_ge2of3',len(soft),len(structural),'URAT1 A2 AND NLRP3 structural in >=2/3 seeds',C1+'07_clinical_dock/acid_dual_a2/','scripts/build_c5_tier_assignment.py'),
        ('intersection_with_historical_chemistry_40',len(structural),len(eligible),'intersection with historical eligible set (not sequential 54->40)',C1+'08_nomination/acid_a2_eligible_audited.csv','scripts/build_c5_tier_assignment.py'),
        ('primary',len(eligible),len(prim),'A1 seed42; remove structural control',C5+'04_shortlist_frozen/primary_candidates.csv','scripts/freeze_c5_shortlist.py'),
        ('reserve_parallel_branch',len(eligible),len(backup),'not A1; beta-lactam exclusion; not second-best potency',C5+'04_shortlist_frozen/backup_candidates.csv','scripts/freeze_c5_shortlist.py')]
    write('data/manuscript/tables/screening_funnel.csv',pd.DataFrame(funnel,columns=['stage','input_n','output_n','rule','source','historical_script']))
    facts={'schema_version':1,'audited_base_commit':base,'status':'evidence_layer_audited_limitations_accepted_md_not_run',
        'counts':{'clinical_library':len(library),'nlrp3_training_molecules':records.canonical_smiles.nunique(),
        'nlrp3_training_records':len(records),'nlrp3_unique_assays':records['Assay ChEMBL ID'].nunique(),
        'nlrp3_one_hot_assays':25,'prediction_assay_contexts':5,'production_pool':len(pool),'acid_equivalents':len(acids),
        'chemistry_soft':len(soft),'primary':len(prim),'reserve':len(backup),'historical_chemistry_eligible':len(audited),
        'dual_structural_ge2of3':len(structural),'intersection_before_final_exclusions':len(eligible)},
        'nlrp3_cv':train['nlrp3']['cv_metrics'], 'urat1_a1_arg_max_A':7.7027,
        'urat1_a1_lr_plus':gates[0]['lr_plus'],'nlrp3_w4_seed42':{k:int(s42[k]) for k in ['tp','fn','fp','tn']},
        'production_receptors':['9DKB','7ALV'],'docking_seeds':[42,43,44],
        'a1_nomination_seed':42,'w2_role':'annotation; 8 of 12 primary pass >=2/3; not membership criterion',
        'pf_followup':'PF-03882845','pf_ligand_id':'REP_07580','md_authorized':False,'md_results_available':False,
        'structure_only_counterfactual_additions':delta.ligand_id.tolist(),
        'limitations':['historical chemistry name exclusions remain part of frozen cohort provenance',
        'acid atom matcher recognizes COOH/carboxylate only, not all pool acid equivalents',
        'items 2 chemistry rules accepted as limitations; 12/21 not replaced',
        'model hash identifies deposited artifact, not proven historical generator of frozen pool',
        '1377 is SI retrain only; multithreading is not an established cause',
        'W4 includes scaffold-concentrated positives and five all-seed failed decoys',
        'A1 has no explicit clash or IFP gate; W2 annotations cannot be assigned to all 12',
        'C1 production acid docks have SDFs without GNINA version logs',
        'GNINA 1.3.1 is documented in old text but not observed in deposited logs',
        'PF-03882845 MR papers do not establish URAT1 apical exposure or a new target',
        'MD protocol is locked but unauthorized and unexecuted; no trajectory numbers']}
    write('MANUSCRIPT_FACTS.yaml',facts)
    write('data/manuscript/audit/checks.json',CHECKS)
    write('data/manuscript/audit/environment.json',{'python':platform.python_version(),'pandas':pd.__version__,
        'rdkit':rdBase.rdkitVersion,
        'packages':{p:version(p) for p in ['pandas','numpy','scipy','PyYAML','rdkit']},
        'audit_kind':'deposited rows/set membership/statistics/descriptors; no new docking or ML fit'})
    if args.inventory:
        inventory=[]
        entries=subprocess.check_output(['git','ls-tree','-r','--long',base],cwd=ROOT.parent, text=True).splitlines()
        batch=subprocess.Popen(['git','cat-file','--batch'],cwd=ROOT.parent,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        for line in entries:
            meta,path=line.split('\t',1); fields=meta.split()
            p=ROOT.parent/path
            relative=path.removeprefix(ROOT.name+'/')
            if not path.startswith(ROOT.name+'/'):
                role='other_project_out_of_edit_scope'
            elif relative.startswith('data/raw/'):
                role='raw_immutable'
            elif relative.startswith(('data/campaigns/','docking_export_','data/repurposing/p2/')):
                role='historical_evidence_retained_in_place'
            elif relative.startswith('docs/'):
                role='development_document_archived_or_superseded'
            else:
                role='retained_source_or_implementation'
            batch.stdin.write((fields[2]+'\n').encode());batch.stdin.flush()
            header=batch.stdout.readline().decode().split()
            body=batch.stdout.read(int(header[2]));batch.stdout.read(1)
            inventory.append(dict(path=path,git_blob=fields[2],bytes=fields[3],classification=role,
                sha256=hashlib.sha256(body).hexdigest(),base_commit=base))
        batch.stdin.close();batch.stdout.close();batch.wait()
        write('data/manuscript/audit/base_file_inventory.csv',pd.DataFrame(inventory))
    # Source files are stable; include manuscript generators/config as provenance too.
    source('scripts/manuscript_pipeline/build_evidence.py')
    for rel in ['config/docking_final.yaml','config/urat1_gate_definition.yaml','config/nlrp3_gate_definition.yaml','config/md_protocol.yaml',
                'data/manuscript/tables/gnina_execution_log_summary.csv',
                'data/manuscript/si/literature_primary_sources.md','data/manuscript/si/chemistry_limitation_acceptance.md',
                'data/manuscript/si/published_dual_structures.csv','docs/FOLLOWUP_EXPERIMENTS.md']:
        source(rel)
    for rel in ['data/manuscript/audit/base_file_inventory.csv','data/manuscript/audit/geometry_recheck.json',
                'data/manuscript/audit/clinical_pose_recheck.csv']:
        if (ROOT/rel).exists():
            OUTPUTS.add(rel)
    write('data/manuscript/provenance.json',{'base_commit':base,
        'hash_policy':'SHA256; text extensions normalized CRLF to LF; binary exact bytes; base inventory hashes Git blob bytes',
        'inputs':list(INPUTS.values()),
        'outputs':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(OUTPUTS)],
        'raw_and_legacy_outputs_modified':False})
    failures=[x for x in CHECKS if not x['passed']]
    print(json.dumps({'checks':len(CHECKS),'failed':failures,'counts':facts['counts'],
                      'counterfactual_additions':facts['structure_only_counterfactual_additions']},ensure_ascii=False,indent=2))
    if failures:
        raise SystemExit(1)


if __name__=='__main__':
    main()
