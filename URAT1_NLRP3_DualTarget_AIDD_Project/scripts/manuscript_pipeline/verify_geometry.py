#!/usr/bin/env python3
"""Re-read all 936 clinical SDFs and audit gates without writing legacy results.

Uses deposited metric definitions; tests reproducibility of implementation, not
external scientific validity. W2 is checked for all frozen candidate/reserve poses.
"""
from pathlib import Path
import json
import sys
import math
import hashlib
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from parse_c1_sdf_readouts import load_ref_centroid,load_poses
from c1_acid_pose_selection import evaluate_urat1_acid_sdf
from c1_nlrp3_pose_metrics import evaluate_nlrp3_structural,load_key_map,load_receptor_heavy,crystal_reference_ifp
from run_c5_w2_urat1_ifp_gate import build_crystal_anchors,evaluate_sdf


def sha(p):
    return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    c1=ROOT/'data/campaigns/c1'
    ref=c1/'01_ligand_prep/selfdock_refs'
    arg=ref/'arg477_coords.json'
    # Use the actual reference path recorded in the deposited metadata.
    les=ref/'lesinurad_A1AIL_crystal_ref.sdf'
    if not les.exists():
        candidates=list(ref.glob('*lesinurad*ref.sdf'))
        if len(candidates)!=1: raise RuntimeError('Ambiguous lesinurad reference')
        les=candidates[0]
    com=load_ref_centroid(les)
    nk=load_key_map(); nr=load_receptor_heavy(); nifp,nh,ncom=crystal_reference_ifp()
    uk=json.loads((ROOT/'data/campaigns/c5/02_urat1_ifp/urat1_key_residues.json').read_text())['residues']
    ur=load_receptor_heavy(ROOT/'data/structures/prepared/9DKB_receptor.pdbqt')
    atoms=json.loads(arg.read_text())['atoms']; anchors=build_crystal_anchors(uk,ur,atoms)
    soft=pd.read_csv(c1/'07_clinical_dock/acid_pool/acid_clinical_pool_chemistry_pass.csv')
    membership=pd.read_csv(ROOT/'data/manuscript/tables/candidate_seed_metrics.csv')
    chosen=set(membership.ligand_id)
    rows=[];discrepancies=[]
    def compare(lid,seed,target,key,actual,stored,tol=1e-5):
        equal=(actual==stored) if isinstance(actual,bool) else (pd.isna(actual) and pd.isna(stored)) or (actual is not None and pd.notna(stored) and abs(float(actual)-float(stored))<=tol)
        if not equal:
            discrepancies.append(dict(ligand_id=lid,seed=seed,target=target,metric=key,actual=str(actual),stored=str(stored)))
    for seed in [42,43,44]:
        folder=c1/'07_clinical_dock'/('acid_dual' if seed==42 else 'acid_dual_a2')
        ast=pd.read_csv(c1/f'07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{seed}.csv').set_index('ligand_id')
        a1=pd.read_csv(c1/'07_clinical_dock/acid_dual_a1_frozen/acid_dual_keep_seed42.csv').set_index('ligand_id')
        for lid in soft.repurposing_id:
            us=folder/'urat1_9dkb'/f'seed{seed}'/f'{lid}_out.sdf'
            ns=folder/'nlrp3_7alv'/f'seed{seed}'/f'{lid}_out.sdf'
            if not ns.exists():
                override=ROOT/'data/campaigns/c5/02_nlrp3_panel/background'/f'seed{seed}'/f'{lid}_out.sdf'
                if override.exists():
                    ns=override
            u=evaluate_urat1_acid_sdf(us,arg,com,lid,seed,rule='a2')
            n=evaluate_nlrp3_structural(ns,lid,seed,key_map=nk,ref_heavy=nh,ref_com=ncom,receptor_heavy=nr,ref_ifp=nifp)
            s=ast.loc[lid]
            compare(lid,seed,'URAT1','keep_urat1_acid',bool(u['keep_urat1_acid']),bool(s.keep_urat1_acid))
            compare(lid,seed,'NLRP3','keep_nlrp3_structural',bool(n['keep_nlrp3_structural']),bool(s.keep_nlrp3_structural))
            if seed==42:
                ua=evaluate_urat1_acid_sdf(us,arg,com,lid,seed,rule='a1')
                compare(lid,seed,'URAT1','A1_dual',bool(ua['keep_urat1_acid'] and n['keep_nlrp3_pose']),bool(a1.loc[lid].keep_dual_acid_geometry))
            if lid in chosen:
                w=evaluate_sdf(us,lid,uk,ur,anchors['ref_heavy'],anchors['ref_ifp'],atoms,anchors['thresholds'])
                old=pd.read_csv(ROOT/f'data/campaigns/c5/04_shortlist_frozen/per_seed/{lid}_w2.csv')
                old=old[old.seed==seed].iloc[0]
                compare(lid,seed,'URAT1','W2',bool(w['keep_urat1_ifp']),bool(old.keep_urat1_ifp))
                for col in ['acid_arg477_min_A','ifp_jaccard_vs_crystal_union','n_key_contacts']:
                    compare(lid,seed,'URAT1',col,w.get(col),old[col])
            rows.append(dict(ligand_id=lid,seed=seed,urat1_sdf=us.relative_to(ROOT).as_posix(),
                nlrp3_sdf=ns.relative_to(ROOT).as_posix(),urat1_sha256=sha(us),nlrp3_sha256=sha(ns),
                urat1_nposes=u['n_poses'],nlrp3_nposes=n.get('n_poses',0),
                urat1_a2_pass=u['keep_urat1_acid'],nlrp3_structural_pass=n['keep_nlrp3_structural']))
        print(f'Checked seed {seed}: {len(rows)} paired jobs',flush=True)
    out=ROOT/'data/manuscript/audit'
    pd.DataFrame(rows).to_csv(out/'clinical_pose_recheck.csv',index=False)
    inputs=[Path(__file__),ROOT/'scripts/c1_acid_pose_selection.py',ROOT/'scripts/c1_nlrp3_pose_metrics.py',
            ROOT/'scripts/parse_c1_sdf_readouts.py',ROOT/'scripts/run_c5_w2_urat1_ifp_gate.py',arg,les,
            ROOT/'data/structures/prepared/9DKB_receptor.pdbqt',ROOT/'data/structures/prepared/7ALV_receptor.pdbqt',
            ROOT/'data/campaigns/c5/02_urat1_ifp/urat1_key_residues.json',ref/'nlrp3_key_residues.json']
    result={'paired_jobs':len(rows),'sdf_files':2*len(rows),'gate_mismatches':discrepancies,
            'inputs':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in inputs],
            'hash_policy':'SHA256 text normalized CRLF to LF',
            'method':'re-evaluation of deposited definitions; no docking and no frozen membership changes'}
    (out/'geometry_recheck.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'sdf_files':2*len(rows),'mismatches':len(discrepancies)},indent=2))
    if discrepancies: raise SystemExit(1)


if __name__=='__main__':main()
