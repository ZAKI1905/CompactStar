#!/usr/bin/env python3
"""Serial candidate gate; fresh characterized profiles are a CTest fixture.

Expected background and old-route mathematics remain in the independent
planning producer. Every candidate local/global/reduction/spin operation is C++.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('executable', type=Path)
    p.add_argument('characterization_root', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--repeat', action='store_true')
    args = p.parse_args()
    evidence = list(args.characterization_root.glob('run-*/result.json'))
    assert evidence, 'fresh characterization fixture missing'
    characterized = max(evidence, key=lambda path: path.stat().st_mtime_ns)
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix='run-', dir=args.output))
    runner = Path(__file__).with_name('chemical_production_evidence.py')
    # Fail the actual certificate builder if its robust neutron enclosure is
    # replaced by the bare mathematical onset. No observed-answer padding.
    import chemical_production_evidence as certificate_builder
    fixed=json.loads((Path(__file__).resolve().parents[2]/'docs/validation/phase5c2_preproduction_evidence.json').read_text())
    character=json.loads(characterized.read_text())
    original=certificate_builder.Model.neutron_source_upper
    certificate_builder.Model.neutron_source_upper=lambda model:(model.onsets[0],0.,0.)
    try:
        try:
            certificate_builder.certificates(characterized.parent/'t8192-r80000',character,fixed['goals'])
        except AssertionError:
            pass
        else:
            raise RuntimeError('N7 bare-onset mutant was not refused')
    finally:
        certificate_builder.Model.neutron_source_upper=original
    artifacts = []
    for fixture in ('t4096-r40000', 't8192-r40000', 't8192-r80000'):
        command = [sys.executable, str(runner), str(args.executable),
                   str(characterized), str(root/fixture), '--fixture', fixture]
        if fixture != 't8192-r80000':
            command.append('--global-only')
        with (root/(fixture+'.log')).open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        print(f'{fixture} raw rc={result.returncode}', flush=True)
        assert result.returncode == 0, f'production failure: {root/fixture}'
        artifacts.append(json.loads((root/fixture/'candidate.json').read_text()))
    gs = [np.array(a['values']['G']) for a in artifacts]
    spread = 2*(abs(gs[1]-gs[0])+abs(gs[2]-gs[1]))
    assert np.all(spread <= np.array(artifacts[-1]['values']['E_background']))
    first = root/'t8192-r80000/candidate.json'
    hashes = [hashlib.sha256(first.read_bytes()).hexdigest()]
    if args.repeat:
        command = [sys.executable, str(runner), str(args.executable),
                   str(characterized), str(root/'independent-repeat')]
        with (root/'independent-repeat.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        print(f'independent repeat raw rc={result.returncode}', flush=True)
        assert result.returncode == 0
        second = root/'independent-repeat/candidate.json'
        hashes.append(hashlib.sha256(second.read_bytes()).hexdigest())
        assert first.read_bytes() == second.read_bytes()
    report = {'PASS': True, 'raw_return_codes': [0]*len(artifacts),
              'production_table_profile_ladder': [a['values']['G'] for a in artifacts],
              'background_componentwise_spread': spread.tolist(),
              'candidate_sha256': hashes, 'repeat_requested': args.repeat,
              'candidate_only': True, 'governed_baseline': False}
    (root/'validation.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(root, flush=True)


if __name__ == '__main__':
    main()
