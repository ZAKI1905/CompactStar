#!/usr/bin/env python3
"""Run the real TaskManager and compare every deterministic artifact exactly."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--build-mode', choices=['Debug', 'Release'], required=True)
    p.add_argument('--mode', choices=['T1', 'T2', 'replay'], required=True)
    p.add_argument('--eos', type=Path)
    p.add_argument('--observe', action='store_true')
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    fixtures = Path(__file__).with_name('fixtures')
    parent = Path(tempfile.mkdtemp(prefix='csm-', dir='/private/tmp'))
    root = parent / 'run'
    try:
        start = time.monotonic()
        if a.mode == 'replay':
            expected = fixtures / f'taskmanager-{a.build_mode}.oracle.tsv'
            actual = parent / 'replay.tsv'
            subprocess.run([a.exe, fixtures / f'taskmanager-{a.build_mode}.bits', actual], check=True)
            differences = [] if actual.read_bytes() == expected.read_bytes() else ['replay.tsv']
            count = len(actual.read_text().splitlines())
        else:
            record = json.loads((fixtures / f'{a.build_mode}-{a.mode}.json').read_text())
            (root / 'EOS').mkdir(parents=True)
            if hashlib.sha256(a.eos.read_bytes()).hexdigest() != '5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd':
                raise RuntimeError('EOS identity mismatch')
            shutil.copyfile(a.eos, root / 'EOS/DS(CMF)-1_with_crust.eos')
            if a.mode == 'T1':
                for rel in ['NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv', 'EOS/Fermi_Gas_0.8mn.eos']:
                    dest = root / rel
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    dest.write_bytes(gzip.decompress((fixtures / a.build_mode / (rel + '.gz')).read_bytes()))
            with (parent / 'run.log').open('w') as log:
                subprocess.run([a.exe, root, a.mode], stdout=log, stderr=log, check=True,
                               env={**os.environ, 'MPLBACKEND': 'Agg',
                                    'CSM_OBSERVATIONS': str(parent / 'observations.tsv')})
            actual = {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest()
                      for f in root.rglob('*') if f.is_file() and f.suffix in ('.tsv', '.eos')
                      and f.name != 'DS(CMF)-1_with_crust.eos'
                      and (a.mode == 'T2' or str(f.relative_to(root)) not in (
                          'NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv', 'EOS/Fermi_Gas_0.8mn.eos'))}
            differences = [k for k in sorted(actual.keys() | record['files'].keys())
                           if actual.get(k) != record['files'].get(k)]
            count = len(actual)
            if a.observe:
                expected_observations = gzip.decompress((fixtures / f'{a.build_mode}-observations.tsv.gz').read_bytes())
                if (parent / 'observations.tsv').read_bytes() != expected_observations:
                    differences.append('observations.tsv')
        result = dict(mode=a.mode, build_mode=a.build_mode, count=count,
                      differing_files=differences, wall_seconds=time.monotonic() - start,
                      retained_run=str(parent) if differences else None)
        result['executable_sha256'] = hashlib.sha256(a.exe.read_bytes()).hexdigest()
        if a.mode != 'replay':
            result['artifact_sha256'] = actual
        if a.mode == 'T2':
            lifetime_root = 'NStar/Dark_Core/0.8/B_conts/2.01/BNV_tau/'
            result['bnv_lifetime_sha256'] = {
                species: actual[lifetime_root + f'BNV_tau_{species}.tsv']
                for species in ('neutron', 'lambda', 'sigmam')
            }
        if a.observe:
            observed = (parent / 'observations.tsv').read_bytes()
            result['observation_records'] = len(observed.splitlines())
            result['observations_sha256'] = hashlib.sha256(observed).hexdigest()
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result))
        if differences:
            raise RuntimeError('Exact historical comparison failed')
    except Exception:
        print(f'Failure artifacts retained at {parent}')
        raise
    else:
        shutil.rmtree(parent)


if __name__ == '__main__':
    main()
