#!/usr/bin/env python3
"""Authenticate, regress retained outputs and publish a scheduling-only snapshot.

No source trajectory, control trajectory, stellar model or new reference solve.
"""
import copy
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'docs/validation/results-first-bnv-pilot'
sys.path.insert(0,str(ROOT/'tools/linux/native_pair'))
from snapshot import digest, encoded, load_pair
from runtime import runtime_binding, verify_runtime

ALLOWED = {'tests/native/stage.cpp','tests/CMakeLists.txt'}
BASE_COMMIT = '88da80156ae2be37bd8ad57a50034f133eaff8ab'

def authenticate():
    schema=json.loads((BASE/'parallel-preparation/TRANSPORT_SCHEMA.json').read_text())
    for name,value in schema['production_header_sha256'].items():
        assert digest(ROOT/name)==value,('IPC schema requires explicit review for changed declaration',name)
    directory = BASE/'restart-preparation/snapshots'
    ids = json.loads((directory/'identities.json').read_text())
    pair = []
    changes = []
    for role in ('source','control'):
        path=directory/(role+'.json')
        assert digest(path)==ids[role]
        snapshot=json.loads(path.read_text());pair.append(snapshot)
        for item in snapshot['files']:
            p=Path(item['path']);actual=digest(p)
            if actual!=item['sha256']:
                name=str(p.relative_to(ROOT))
                assert name in ALLOWED,('unapproved bound source change',name)
                committed=subprocess.check_output(['git','show',BASE_COMMIT+':'+name],cwd=ROOT)
                import hashlib
                assert hashlib.sha256(committed).hexdigest()==item['sha256']
                if role=='source':changes.append(dict(path=name,before=item['sha256'],after=actual))
        verify_runtime(snapshot['runtime_bindings'],{i['path']:i['sha256'] for i in snapshot['files']})
    a,b=copy.deepcopy(pair)
    del a['source'];del b['source'];assert a==b
    assert set(i['path'] for i in changes)==ALLOWED
    assert not subprocess.check_output(['git','diff',BASE_COMMIT,'--name-only','--','CompactStar'],cwd=ROOT)
    # Authenticate every retained input/output, including the original failed run.
    for rel in ('stage-a-source-retry/SHA256SUMS','stage-a-source-retry/evidence/external-SHA256SUMS'):
        for line in (BASE/rel).read_text().splitlines():
            expected,name=line.split('  ',1);path=Path(name)
            if not path.is_absolute():path=ROOT/path
            if path.is_relative_to(ROOT) and str(path.relative_to(ROOT)) in ALLOWED:
                continue # only the authorized stage.cpp dependency changed
            assert digest(path)==expected,('retained evidence authentication failed',str(path))
    return pair,ids,changes


def run_checked(args,log):
    with log.open('w') as f:subprocess.run(list(map(str,args)),cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)


def prepare(run):
    pair,ids,changes=authenticate();s=copy.deepcopy(pair[0])
    assert '100% tests passed' in (run/'contracts.log').read_text()
    assert 'PASS 24 distinct' in (run/'pool-tests.log').read_text()
    assert 'PASS all IPC fields' in (run/'transport-tests.log').read_text()
    original=ROOT.parents[1]/'stage-a-restart-retry-20261109-r1/source-stage/run/trajectory/source/tier-0-cap-0'
    regression=run/'build/tests/native_parallel_regression'
    run_checked([regression,s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],
        s['archive'],s['archive_sha256'],run/'regression-assembly',original,run/'regression'],run/'regression.log')
    comparisons=[]
    for workers in (1,24):
        actual=run/'regression'/('workers-'+str(workers))
        names=['diagnostics-240.tsv','uncertainties-240.tsv']
        names += [f'checkpoints-240/point-{i}/{name}' for i in range(241)
                  for name in ('checkpoint.json','state.tsv','diagnostics.tsv','failure.txt')]
        for name in names:
            assert (original/name).read_bytes()==(actual/name).read_bytes(),('retained byte mismatch',workers,name)
        old_rows=(original/'R20-grids.tsv').read_text().splitlines()
        assert (actual/'R20-grids.tsv').read_text().splitlines()==old_rows[:2],('R20 bytes changed',workers)
        comparisons.append(dict(workers=workers,identical_files=len(names)+1,
            grid_intervals=240,packets=241,grid_wall_seconds=float((actual/'timing.txt').read_text()),
            R20_byte_identical=True))
    receipt=dict(status='PASS',retained_original=str(original),comparisons=comparisons,
        speedup=comparisons[0]['grid_wall_seconds']/comparisons[1]['grid_wall_seconds'],
        main_integrations=0,stellar_solves=0,new_reference_trajectories=0,
        scope='Scheduling regression only. Does not qualify source accuracy or discharge actual-output HP/prefix/roundoff gates.')
    (run/'regression.json').write_bytes(encoded(receipt))
    # Reuse authenticated static verifier; only the harness ELF needs replacement.
    verifier=Path(s['verify_binary']);stage=run/'build/tests/native_p2_stage'
    run_checked([verifier,s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],
        s['archive'],s['archive_sha256'],run/'static-work',run/'static-rhs.tsv'],run/'verify.log')
    assert (run/'static-rhs.tsv').read_bytes()==(BASE/'stage-a-recovery/evidence/static-rhs.tsv').read_bytes()
    smoke=subprocess.run([str(stage)],cwd=ROOT,capture_output=True,text=True)
    (run/'stage-smoke.log').write_text(smoke.stdout+smoke.stderr)
    assert smoke.returncode==1 and 'STOP profile thermal entry' in smoke.stderr
    old_stage=s['stage_binary'];s['stage_binary']=str(stage)
    s['runtime_bindings']={str(p):runtime_binding(p) for p in (verifier,stage)}
    assert s['runtime_bindings'][str(stage)]==pair[0]['runtime_bindings'][old_stage]
    paths={Path(str(stage) if x['path']==old_stage else x['path']) for x in s['files']}
    paths.update(ROOT/'tests/native'/p for p in ('process_pool.hpp','checkpoint_transport.hpp',
        'process_pool_test.cpp','checkpoint_transport_test.cpp','parallel_regression.cpp'))
    paths.update(ROOT/'tools/linux/native_qualification'/p for p in
        ('parallel_preparation.py','parallel_preparation.sbatch','parallel_stage.sbatch'))
    paths.add(BASE/'parallel-preparation/OWNER_APPROVAL.md')
    paths.add(BASE/'parallel-preparation/TRANSPORT_SCHEMA.json')
    s['files']=[dict(path=str(p),sha256=digest(p)) for p in sorted(paths)]
    s['assembly_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    s['checkpoint_execution']=dict(method='24 isolated pinned Linux processes; immutable accepted history',
        reduction='coordinator calls unchanged ComputeCheckpointR20 in ascending observation order',
        workers=24,parallel_scope='source Stage A only',ipc='same-build exact binary64 complete CheckpointOutput',
        failure='stop dispatch, preserve evidence, terminate and reap all workers; no retry')
    before,after=copy.deepcopy(pair[0]),copy.deepcopy(s)
    for k in ('files','runtime_bindings','assembly_commit','stage_binary'):
        del before[k];del after[k]
    del after['checkpoint_execution'];assert before==after,'physical/numerical setting changed'
    out=run/'snapshots';out.mkdir(exist_ok=False);new_ids={}
    for role,prior in zip(('source','control'),pair):
        s['source']=prior['source'];p=out/(role+'.json');p.write_bytes(encoded(s));new_ids[role]=digest(p)
    (out/'identities.json').write_bytes(encoded(new_ids));load_pair(out)
    (run/'snapshot-diff.json').write_bytes(encoded(dict(status='PASS',old_identities=ids,new_identities=new_ids,
        changed_bound_sources=changes,old_stage=old_stage,new_stage=str(stage),archive_sha256=s['archive_sha256'],
        unchanged_physics_and_numerics=True,static_rhs_identical=True,runtime_linkage_identical=True,
        source_qualified=False,control_executed=False)))
    print('PASS serial/24-worker retained-grid identity, ordered R20, private workers, authenticated replacement snapshots')

if __name__=='__main__':
    if sys.argv[1]=='authenticate':
        pair,ids,changes=authenticate()
        print(json.dumps(dict(status='PASS',old_identities=ids,changed_bound_sources=changes),indent=2))
    else:prepare(Path(sys.argv[1]).resolve())
