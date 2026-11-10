#!/usr/bin/env python3
"""Publish new runtime identities while retaining previous snapshots and inputs."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/linux/native_pair'))
from snapshot import digest, encoded, load_pair
from runtime import runtime_binding

START='4387e56c6b218411d755d27148f654aea5bc345a'
BASE=Path('docs/validation/results-first-bnv-pilot')
OLD=ROOT/BASE/'native-runtime-retry/snapshots'
TOLERANCES={'O1':{'relative':1e-12,'absolute':[1e-17,1e-23,1e-23]},
            'O2':{'relative':1e-13,'absolute':[1e-18,5e-27,5e-27]}}
CHANGED={'CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp','tests/native/stage.cpp'}


def verify_publish(run):
    run=Path(run).resolve()
    old=[json.loads((OLD/(r+'.json')).read_text()) for r in ('source','control')]
    ids=json.loads((OLD/'identities.json').read_text())
    for role in ('source','control'):
        assert digest(OLD/(role+'.json'))==ids[role]
    changes=[]
    for item in old[0]['files']:
        actual=digest(item['path'])
        if actual!=item['sha256']:
            relative=str(Path(item['path']).relative_to(ROOT))
            assert relative in CHANGED,('unauthorized snapshot dependency change',relative)
            changes.append(dict(path=relative,before=item['sha256'],after=actual))
    assert {x['path'] for x in changes}==CHANGED
    name='CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp'
    original=subprocess.check_output(['git','show',f'{START}:{name}'],cwd=ROOT)
    expected=original.replace(b'{{1.0e-18,1.0e-24,1.0e-24}}',b'{{1.0e-18,5.0e-27,5.0e-27}}')
    assert (ROOT/name).read_bytes()==expected,'production change exceeds two approved tolerance literals'
    assert '100% tests passed' in (run/'contracts.log').read_text()
    s=copy.deepcopy(old[0])
    verify=run/'build/tests/native_p2_verify'
    stage=run/'build/tests/native_p2_stage'
    args=[str(verify),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],
          s['archive'],s['archive_sha256'],str(run/'static-work'),str(run/'static-rhs.tsv')]
    with (run/'verify.log').open('w') as f:
        subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
    old_static=Path(s['archive']).parents[1]/'static-rhs.tsv'
    assert (run/'static-rhs.tsv').read_bytes()==old_static.read_bytes(),'static RHS changed'
    replacements={s['verify_binary']:str(verify),s['stage_binary']:str(stage)}
    s['verify_binary'],s['stage_binary']=str(verify),str(stage)
    s['runtime_bindings']={str(b):runtime_binding(b) for b in (verify,stage)}
    old_bindings=list(old[0]['runtime_bindings'].values())
    assert all(b==old_bindings[0] for b in s['runtime_bindings'].values()),'runtime linkage changed'
    paths={Path(replacements.get(item['path'],item['path'])) for item in s['files']}
    paths.add(ROOT/'tests/native/o2_configuration.cpp')
    paths.add(ROOT/'tests/diagnostics/checkpoint_evidence.hpp')
    paths.update(ROOT/'tools/linux/native_qualification'/f for f in
                 ('recovery_snapshots.py','recovery_build.sbatch','recovery_stage.sbatch'))
    paths.add(ROOT/BASE/'stage-a-recovery/PREDECLARATION.md')
    s['files']=[dict(path=str(p),sha256=digest(p)) for p in sorted(paths)]
    s['assembly_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    s['checkpoint_configuration']=TOLERANCES
    common_before=copy.deepcopy(old[0]);common_after=copy.deepcopy(s)
    for key in ('files','runtime_bindings','assembly_commit','verify_binary','stage_binary'):
        del common_before[key];del common_after[key]
    del common_after['checkpoint_configuration']
    assert common_before==common_after,'physical/common configuration drift'
    out=run/'snapshots';out.mkdir(exist_ok=False)
    new_ids={}
    for role,previous in zip(('source','control'),old):
        s['source']=previous['source']
        p=out/(role+'.json');p.write_bytes(encoded(s));new_ids[role]=digest(p)
    (out/'identities.json').write_bytes(encoded(new_ids))
    load_pair(out)
    (run/'snapshot-diff.json').write_bytes(encoded(dict(status='PASS',old_identities=ids,
        new_identities=new_ids,changed_bound_sources=changes,binary_replacements=replacements,
        added_checkpoint_configuration=TOLERANCES,physical_inputs_unchanged=True,
        static_RHS_identical=True,archive_sha256=s['archive_sha256'],runtime_linkage_identical=True)))
    print('PASS fresh static owners, unchanged twelve RHS states, seven contracts and matched snapshots')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['verify-publish'])
    p.add_argument('run',type=Path)
    a=p.parse_args();verify_publish(a.run)
