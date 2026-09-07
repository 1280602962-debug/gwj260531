#!/usr/bin/env python3
"""One-time mechanical document retirement; do not run during evidence rebuild.

Original Markdown remains byte-identical in docs/archive or docs/internal. Small
redirects at original paths preserve historical links. Original scripts/configs
remain at their import paths with explicit roles in machine-readable registries.
"""
from pathlib import Path
import csv
import hashlib
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
BASE='cf5fd75833baac94be603ebaf6188022de4d0275'
INTERNAL={'C5_CRITIQUE_RESPONSE.md','C5_JOURNAL_COMPARISON_MD_AFTER.md','MOL_DIVERS_REVISION_PLAN.md',
          'DIFFERENTIATION_VS_PLK1_NLRP3.md','IF_STILL_WANT_CANDIDATES.md','LOCAL_AGENT_TASKS.md'}
CURRENT_IMPL={'00_prepare_data.py','02_train_asymmetric_models.py','screen_repurposing_library.py',
 'build_repurposing_library.py','build_c1_acid_clinical_pool.py','prepare_ligands_c1.py',
 'prepare_receptor_vina.py','run_gnina_batch.py','gnina_dock.py','run_c1_acid_dual_dock.py',
 'run_c1_acid_dual_a2.py','c1_acid_pose_selection.py','c1_nlrp3_pose_metrics.py',
 'parse_c1_sdf_readouts.py','rescore_nlrp3_structural.py','run_acid_gate_benchmark.py',
 'run_c5_w2_urat1_ifp_gate.py','score_c5_w4_nlrp3_panel.py','run_c5_w4_nlrp3_panel.py',
 'sample_c5_w4_decoys.py','build_c5_tier_assignment.py','freeze_c5_shortlist.py',
 'run_c5_presubmission_audit.py','utils_ml.py','dock_score_utils.py'}


def table(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    moves=[]
    old_docs=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE,'docs'],cwd=ROOT,text=True).splitlines()
    for rel in old_docs:
        p=ROOT/rel
        if p.suffix!='.md':
            continue
        if p.read_text(encoding='utf-8-sig').startswith('<!-- manuscript-redirect -->'):
            continue
        area='internal' if p.name in INTERNAL or 'paper_spine_ars_analysis' in rel else 'archive'
        target=ROOT/'docs'/area/'development'/p.relative_to(ROOT/'docs')
        if target.exists():
            raise RuntimeError(f'Archive already exists: {target}')
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,target)
        original_hash=hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
        # Redirect uses relative paths from the old document's parent.
        import os
        arch=os.path.relpath(target,p.parent).replace('\\','/')
        master=os.path.relpath(ROOT/'docs/MANUSCRIPT_MASTER.md',p.parent).replace('\\','/')
        p.write_text('<!-- manuscript-redirect -->\n# Archived development document\n\n'
            '此文件已退出论文主线；保留原路径作为兼容入口。历史数字、候选和路线不能作为现行结论。\n\n'
            f'- [现行论文入口]({master})\n- [原文完整存档]({arch})\n\n'
            f'原文快照：`{BASE}`。存档文件内的相对路径按原位置解释，原始位置见归档索引。\n',encoding='utf-8')
        moves.append(dict(original_path=rel,archive_path=target.relative_to(ROOT).as_posix(),
            original_sha256=original_hash,role=area,original_commit=BASE))
    if moves:
        table(ROOT/'docs/archive/document_registry.csv',moves)
    registry=ROOT/'docs/archive/document_registry.csv'
    if registry.exists():
        with registry.open(encoding='utf-8',newline='') as f:
            registry_rows=list(csv.DictReader(f))
        for r in registry_rows:
            r['original_sha256']=hashlib.sha256((ROOT/r['archive_path']).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            r['hash_policy']='SHA256 after CRLF to LF normalization'
        table(registry,registry_rows)
    # Retain project README before replacing it in a reviewed patch.
    old_readme=ROOT/'docs/archive/README_PROJECT_before_consolidation.md'
    if not old_readme.exists():
        shutil.copyfile(ROOT/'README.md',old_readme)
    configs=[]
    for p in sorted((ROOT/'config').glob('*.yaml')):
        if p.name in {'campaign_final.yaml','docking_final.yaml','chemistry_final.yaml','urat1_gate_definition.yaml',
                      'nlrp3_gate_definition.yaml','md_protocol.yaml'}:
            continue
        target=ROOT/'config/archive'/p.name
        if not target.exists():
            target.parent.mkdir(exist_ok=True)
            shutil.copyfile(p,target)
        text=p.read_text(encoding='utf-8-sig')
        if not text.startswith('# RETAINED HISTORICAL CONFIG'):
            p.write_text('# RETAINED HISTORICAL CONFIG — manuscript policy: config/campaign_final.yaml\n'
                '# Kept at original path for script compatibility; original snapshot in config/archive/.\n'+text,encoding='utf-8')
        configs.append(dict(path=p.relative_to(ROOT).as_posix(),role='historical_configuration_retained_for_compatibility',
                            archive=target.relative_to(ROOT).as_posix()))
    table(ROOT/'config/archive/config_registry.csv',configs)
    scripts=[]
    for p in sorted((ROOT/'scripts').iterdir()):
        if p.suffix not in {'.py','.sh'}:
            continue
        role='retained_implementation_dependency' if p.name in CURRENT_IMPL else 'legacy_or_supporting_only'
        if p.name=='build_c1_acid_shortlist_a2.py':
            role='retired_hardcoded_nomination_do_not_use_for_current_names'
        scripts.append(dict(path=p.relative_to(ROOT).as_posix(),role=role,
            current_entry='scripts/manuscript_pipeline/build_evidence.py',
            note='Original path preserved for imports and audit; not a default execution workflow'))
    table(ROOT/'scripts/archive/script_registry.csv',scripts)
    print(f'Archived {len(moves)} documents; registered {len(configs)} configs and {len(scripts)} scripts.')


if __name__=='__main__':
    main()
