#!/usr/bin/env python3
"""Run existing bounded Phase-6 executable with authenticated relocated inputs.

No input rewriting, Mac platform impersonation, or full-campaign mode.
"""
import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tests/bnv'))
from run_breakpoint_qualification import digest, tree, validate_output

p = argparse.ArgumentParser()
p.add_argument('--handoff', type=Path, required=True)
p.add_argument('--executable', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--mode', choices=['source', 'control', 'bounded', 'events', 'crossings', 'probe'], required=True)
p.add_argument('--timeout', type=int, default=1800)
a = p.parse_args()
a.handoff = a.handoff.resolve()
a.output = a.output.resolve()
a.executable = a.executable.resolve()
# Pin the transferred mapping in addition to the original six-input authority.
map_path = a.handoff/'scientific/input-map.json'
EXPECTED_MAP = 'e88e96fa95661352739d2c1232a62bb96d7c8a7481b91db9627e742693b0675f'
if digest(map_path) != EXPECTED_MAP:
    raise SystemExit('Transferred input map hash mismatch')
inputs = json.loads(map_path.read_text())
paths = {}
for item in inputs:
    path = (a.handoff/'scientific'/item['relative_path']).resolve(strict=True)
    path.relative_to(a.handoff/'scientific')
    if tree(path) != item['sha256']:
        raise SystemExit('Scientific identity mismatch: '+item['key'])
    paths[item['key']] = path
original = json.loads((ROOT/'docs/validation/evidence/obs32-forensics/authentication.json').read_text())
actual = {key: tree(paths[key]) for key in original['input_hashes']}
if actual != original['input_hashes']:
    raise SystemExit('Historical six-input authority mismatch')
a.output.mkdir(parents=True, exist_ok=False)
cmd = [str(a.executable)] + [str(paths[k]) for k in ('profile_tree', 'certificate', 'thermal')]
cmd += [str(a.output/'work')] + [str(paths[k]) for k in ('entry', 'frozen', 'coefficients')]
cmd += [str(a.handoff/'scientific/obs32/inputs'), str(a.output/'output.partial'), a.mode]
record = {'status': 'CANDIDATE_UNQUALIFIED', 'mode': a.mode, 'command': cmd,
          'source_head': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
          'executable_sha256': digest(a.executable), 'inputs': actual,
          'relocated_paths': {k: str(v) for k, v in paths.items()}, 'published': False}
start = time.monotonic()
with (a.output/'console.log').open('w') as log:
    try:
        result = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=a.timeout)
        record['returncode'] = result.returncode
    except subprocess.TimeoutExpired:
        record['returncode'] = 124
        record['timeout_seconds'] = a.timeout
record['seconds'] = time.monotonic() - start
record['maxrss_kib'] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
if record['returncode'] == 0:
    try:
        validate_output(a.output/'output.partial', a.mode)
        (a.output/'output.partial').rename(a.output/'output')
        record['published'] = True
        record['outputs'] = {f.name: digest(f) for f in (a.output/'output').iterdir()}
    except Exception as error:
        record['validation_error'] = str(error)
        record['returncode'] = 2
(a.output/'execution.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record, indent=2))
raise SystemExit(record['returncode'])
