#!/usr/bin/env python3
"""Read-only checks of the manuscript export; exits nonzero on drift."""
from pathlib import Path
import csv
import hashlib
import json
import re
import sys
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[2]


def sha(path):
    if path.suffix.lower() in {'.csv','.tsv','.json','.yaml','.yml','.md','.py','.sh','.txt','.sdf','.pdbqt','.pdb','.cif','.smi','.log','.svg'}:
        return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):
            h.update(b)
    return h.hexdigest()


def main():
    p=json.loads((ROOT/'data/manuscript/provenance.json').read_text(encoding='utf-8'))
    failures=[]
    for row in p['inputs']+p['outputs']:
        path=ROOT/row['path']
        if not path.exists() or sha(path)!=row['sha256']:
            failures.append('hash drift: '+row['path'])
    checks=json.loads((ROOT/'data/manuscript/audit/checks.json').read_text(encoding='utf-8'))
    failures.extend('failed assertion: '+r['check'] for r in checks if not r['passed'])
    # Internal historical links remain interpreted at the original path; only
    # active docs and redirects are checked as live navigation.
    docs=[ROOT/'README.md',ROOT.parent/'README.md']
    for directory in ['docs','data/manuscript','data/frozen','scripts/manuscript_pipeline','config','data/archive','scripts/archive']:
        docs.extend((ROOT/directory).glob('*.md'))
    docs.extend((ROOT/'docs/archive').glob('README.md'))
    for path in docs:
        text=path.read_text(encoding='utf-8-sig')
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',text):
            if link.startswith(('http:','https:','#','mailto:')):
                continue
            target=unquote(link.split('#',1)[0].strip('<>'))
            if target and not (path.parent/target).exists():
                failures.append(f'broken link: {path.relative_to(ROOT.parent)} -> {target}')
    with (ROOT/'docs/archive/document_registry.csv').open(encoding='utf-8',newline='') as f:
        for r in csv.DictReader(f):
            if sha(ROOT/r['archive_path'])!=r['original_sha256']:
                failures.append('archive changed: '+r['archive_path'])
    for name,n in [('primary_candidates_frozen.csv',12),('backup_candidates_frozen.csv',21)]:
        with (ROOT/'data/frozen'/name).open(encoding='utf-8',newline='') as f:
            rows=list(csv.DictReader(f))
        if len(rows)!=n or len({r['ligand_id'] for r in rows})!=n:
            failures.append('candidate count/ID drift: '+name)
    geometry=ROOT/'data/manuscript/audit/geometry_recheck.json'
    if geometry.exists():
        g=json.loads(geometry.read_text(encoding='utf-8'))
        if g['gate_mismatches'] or g['sdf_files']!=936:
            failures.append('clinical geometry recheck failed')
        for r in g.get('inputs',[]):
            if sha(ROOT/r['path'])!=r['sha256']:
                failures.append('geometry input drift: '+r['path'])
        with (ROOT/'data/manuscript/audit/clinical_pose_recheck.csv').open(encoding='utf-8',newline='') as f:
            for r in csv.DictReader(f):
                for target in ['urat1','nlrp3']:
                    if sha(ROOT/r[target+'_sdf'])!=r[target+'_sha256']:
                        failures.append('pose file drift: '+r[target+'_sdf'])
    fig_manifest=ROOT/'data/manuscript/figures/provenance.json'
    if fig_manifest.exists():
        f=json.loads(fig_manifest.read_text(encoding='utf-8'))
        for r in f['files']:
            if sha(ROOT/r['path'])!=r['sha256']:
                failures.append('figure source/output drift: '+r['path'])
    print(json.dumps({'checked_hashes':len(p['inputs'])+len(p['outputs']),
                      'scientific_checks':len(checks),'failures':failures},ensure_ascii=False,indent=2))
    return 1 if failures else 0


if __name__=='__main__':
    sys.exit(main())
