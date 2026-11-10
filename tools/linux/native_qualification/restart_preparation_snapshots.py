#!/usr/bin/env python3
"""Publish the approved restart-harness identity; no trajectory invocation."""
import copy
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'docs/validation/results-first-bnv-pilot'
sys.path.insert(0, str(ROOT / 'tools/linux/native_pair'))
from snapshot import digest, encoded, load_pair
from runtime import runtime_binding


def publish(run):
    old_dir = BASE / 'stage-a-recovery/snapshots'
    ids = json.loads((old_dir / 'identities.json').read_text())
    old = [json.loads((old_dir / (role+'.json')).read_text()) for role in ('source','control')]
    for role in ('source','control'):
        assert digest(old_dir/(role+'.json')) == ids[role]
    expected = json.loads((BASE/'restart-coverage/evidence/result.json').read_text())['candidate_hashes']
    for name, value in expected.items():
        assert digest(ROOT/name) == value, ('not the exact approved patch', name)
    assert not subprocess.check_output(['git','diff','3557da9','--name-only','--','CompactStar'], cwd=ROOT)
    changes = []
    for item in old[0]['files']:
        actual = digest(item['path'])
        if actual != item['sha256']:
            relative = str(Path(item['path']).relative_to(ROOT))
            assert relative == 'tests/native/stage.cpp', ('unexpected dependency change', relative)
            changes.append(dict(path=relative,before=item['sha256'],after=actual))
    assert len(changes) == 1
    assert '100% tests passed' in (run/'contracts.log').read_text()
    assert 'exact state roundtrip' in (run/'selector-roundtrip.log').read_text()
    s = copy.deepcopy(old[0])
    # The verifier did not change: reuse and authenticate its original ELF.
    verifier = Path(s['verify_binary'])
    stage = run/'build/tests/native_p2_stage'
    args = [str(verifier),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],
            s['archive'],s['archive_sha256'],str(run/'static-work'),str(run/'static-rhs.tsv')]
    with (run/'verify.log').open('w') as log:
        subprocess.run(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    reference = BASE/'stage-a-recovery/evidence/static-rhs.tsv'
    assert (run/'static-rhs.tsv').read_bytes() == reference.read_bytes(), 'static RHS changed'
    smoke = subprocess.run([str(stage)],cwd=ROOT,capture_output=True,text=True)
    (run/'stage-linkage-smoke.log').write_text(smoke.stdout+smoke.stderr)
    assert smoke.returncode == 1 and 'STOP profile thermal entry frozen coefficients archive sha fresh-work years output' in smoke.stderr
    old_stage = s['stage_binary']
    s['stage_binary'] = str(stage)
    s['runtime_bindings'] = {str(p):runtime_binding(p) for p in (verifier,stage)}
    assert s['runtime_bindings'][str(stage)] == old[0]['runtime_bindings'][old_stage], 'stage runtime linkage changed'
    assert s['runtime_bindings'][str(verifier)] == old[0]['runtime_bindings'][str(verifier)]
    paths = {Path(str(stage) if item['path']==old_stage else item['path']) for item in s['files']}
    paths.update(ROOT/p for p in expected)
    paths.update(ROOT/'tools/linux/native_qualification'/p for p in
                 ('restart_preparation_snapshots.py','restart_preparation.sbatch'))
    paths.add(BASE/'restart-preparation/OWNER_APPROVAL.md')
    paths.add(BASE/'restart-coverage/proposed-selector.patch')
    s['files'] = [dict(path=str(p),sha256=digest(p)) for p in sorted(paths)]
    s['assembly_commit'] = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    s['restart_test_configuration'] = dict(stage='A',horizon_s=31557600,
        selection='latest accepted noninitial endpoint at or before half-time',
        missing_endpoint='fail explicitly; no passive or initial-state substitution',
        minimum_tail_s=15778800,stage_B_rule_unchanged=True,
        approval='docs/validation/results-first-bnv-pilot/restart-preparation/OWNER_APPROVAL.md')
    before,after = copy.deepcopy(old[0]),copy.deepcopy(s)
    for key in ('files','runtime_bindings','assembly_commit','stage_binary'):
        del before[key];del after[key]
    del after['restart_test_configuration']
    assert before == after, 'physical/numerical configuration drift'
    out = run/'snapshots';out.mkdir(exist_ok=False)
    new_ids = {}
    for role,prior in zip(('source','control'),old):
        s['source'] = prior['source']
        p = out/(role+'.json');p.write_bytes(encoded(s));new_ids[role] = digest(p)
    (out/'identities.json').write_bytes(encoded(new_ids))
    load_pair(out)
    (run/'snapshot-diff.json').write_bytes(encoded(dict(status='PASS',old_identities=ids,
        new_identities=new_ids,changed_bound_sources=changes,
        stage_binary_replacement=dict(before=old_stage,after=str(stage)),
        verifier_binary_retained=str(verifier),physical_inputs_unchanged=True,
        production_sources_unchanged=True,numerical_settings_unchanged=True,
        exact_approved_patch=True,static_twelve_RHS_states_identical=True,
        runtime_linkage_identical=True,archive_sha256=s['archive_sha256'],
        trajectory_calls=0,prepared_assembly_commit=s['assembly_commit'])))
    print('PASS exact amendment, static RHS identity, contracts, roundtrip and authenticated matched snapshots; zero trajectories')


if __name__ == '__main__':
    publish(Path(sys.argv[1]).resolve())
