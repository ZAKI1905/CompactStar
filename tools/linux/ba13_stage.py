#!/usr/bin/env python3
"""Isolated BA13 recovery stages. Unpublished *.partial directories retain failures.

The recovered producers are invoked unchanged. Only successful, structurally complete
stages receive a hash manifest and are atomically renamed. No trajectory is invoked.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
sys.dont_write_bytecode = True


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(p, count):
    with p.open() as f:
        data = list(csv.DictReader(f, delimiter='\t'))
    if len(data) != count:
        raise RuntimeError(f'{p}: expected {count} rows, got {len(data)}')
    for row in data:
        if any(value is None or not math.isfinite(float(value)) for value in row.values()):
            raise RuntimeError(f'{p}: incomplete or nonfinite numeric row')
    return data


def verify(stage):
    manifest = json.loads((stage/'COMPLETE.json').read_text())
    actual = {str(p.relative_to(stage)) for p in stage.rglob('*') if p.is_file() and p.name != 'COMPLETE.json'}
    if actual != set(manifest['files']):
        raise RuntimeError('completed stage file inventory mismatch')
    for name, sha in manifest['files'].items():
        if digest(stage/name) != sha:
            raise RuntimeError(f'changed completed stage input: {stage/name}')


def main():
    if sys.flags.optimize:
        raise RuntimeError('BA13 requires active upstream Python assertions')
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['stars', 'adapters', 'coefficient'])
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--characterization', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--input', type=Path)
    p.add_argument('--index', type=int)
    p.add_argument('--native-tail', action='store_true')
    a = p.parse_args()
    source, exe, char, dest = [x.resolve() for x in (a.source,a.executable,a.characterization,a.output)]
    profile = char.parent/'t8192-r80000'
    inherited = json.loads(char.read_text())
    for name, expected in inherited['files'].items():
        if digest(char.parent/name) != expected:
            raise RuntimeError('characterization input hash mismatch: '+name)
    if dest.exists() or dest.with_suffix('.partial').exists():
        raise RuntimeError('fresh output root required')
    partial = dest.with_suffix('.partial')
    partial.parent.mkdir(parents=True, exist_ok=True)
    if a.stage == 'stars':
        cmd = [str(exe), 'stars', str(profile), str(partial)]
    else:
        if a.input is None:
            raise RuntimeError('completed input stage required')
        inp = a.input.resolve()
        verify(inp)
        partial.mkdir()
        for name in ('target_stars.tsv', 'tangent.tsv'):
            shutil.copyfile(inp/name, partial/name)
        if a.stage == 'adapters':
            shutil.copytree(inp/'chemical-profiles', partial/'chemical-profiles')
            cmd = [sys.executable, '-B', str(source/'tests/bnv'/('produce_native_ba13_inputs.py' if a.native_tail else 'produce_frozen_certificates.py')),
                   '--source-root', str(source), '--profiles', str(partial/'chemical-profiles'),
                   '--phase5-characterization', str(char)]
        else:
            if a.index is None or not 0 <= a.index < 21:
                raise RuntimeError('coefficient index must be 0..20')
            name = f'target-{a.index}'
            shutil.copytree(inp/'chemical-profiles'/name, partial/'chemical-profiles'/name)
            cmd = [str(exe), 'coefficient', str(profile), str(partial), str(a.index)]
    log = dest.with_suffix('.log')
    start = time.monotonic()
    with log.open('x') as f:
        result = subprocess.run(cmd, cwd=source, stdout=f, stderr=subprocess.STDOUT,
                                env={**os.environ, 'PYTHONDONTWRITEBYTECODE':'1', 'PYTHONOPTIMIZE':'0'})
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt = dict(command=cmd, cwd=str(source), exit_code=result.returncode,
        wall_seconds=time.monotonic()-start, user_seconds=usage.ru_utime,
        system_seconds=usage.ru_stime, max_rss_kib=usage.ru_maxrss,
        slurm_job=os.environ.get('SLURM_JOB_ID'), node=os.uname().nodename,
        executable_sha256=digest(exe), characterization_sha256=digest(char))
    dest.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    if result.returncode:
        raise RuntimeError(f'producer refused; retained {partial} and {log}')
    targets = rows(partial/'target_stars.tsv', 21)
    tangent = rows(partial/'tangent.tsv', 21)
    if [int(r['index']) for r in targets] != list(range(21)):
        raise RuntimeError('target order mismatch')
    if a.stage == 'stars':
        for i in range(21):
            for name in ('profile.tsv','model.txt','windows.tsv','freegas.tsv'):
                if not (partial/'chemical-profiles'/f'target-{i}'/name).stat().st_size:
                    raise RuntimeError('empty star input')
    elif a.stage == 'adapters':
        for i in range(21):
            for name in ('certificate.txt','certificate.json','thermal/eos.t','thermal/eos.nb','thermal/eos.yq','thermal/eos.thermo'):
                if not (partial/'chemical-profiles'/f'target-{i}'/name).stat().st_size:
                    raise RuntimeError('empty adapter input')
    else:
        rows(partial/f'coefficients-shard-{a.index}.tsv', 1)
        rows(partial/f'thermal_surface-shard-{a.index}.tsv', 13)
        rows(partial/f'profile_metrics-shard-{a.index}.tsv', 2049)
    files = {str(x.relative_to(partial)):digest(x) for x in sorted(partial.rglob('*')) if x.is_file()}
    # Close, flush and sync every output before committing the completion marker.
    for path in partial.rglob('*'):
        if path.is_file():
            with path.open('rb') as f: os.fsync(f.fileno())
    marker = partial/'COMPLETE.json'
    with marker.open('x') as f:
        f.write(json.dumps(dict(stage=a.stage,files=files),indent=2)+'\n')
        f.flush(); os.fsync(f.fileno())
    fd=os.open(partial,os.O_RDONLY); os.fsync(fd); os.close(fd)
    partial.rename(dest)
    fd=os.open(dest.parent,os.O_RDONLY); os.fsync(fd); os.close(fd)
    print(f'COMPLETE {a.stage} {dest}', flush=True)

if __name__ == '__main__':
    main()
