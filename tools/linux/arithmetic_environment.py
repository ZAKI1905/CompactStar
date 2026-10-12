#!/usr/bin/env python3
"""Runtime qualification receipt, supplementing the unchanged package lock."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import os

p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
probe=a.build/'tests/linux_arithmetic_probe'
raw=subprocess.check_output([str(probe)],text=True)
loaded=sorted({r.split('\t',1)[1] for r in raw.splitlines() if r.startswith('loaded\t')})
info=Path('/proc/cpuinfo').read_text().splitlines()
cpus={key:sorted({r.split(':',1)[1].strip() for r in info if r.startswith(key+'\t') or r.startswith(key+' ')}) for key in ('vendor_id','model name','flags')}
commands=json.loads((a.build/'compile_commands.json').read_text())
assert all('-ffp-contract=off' in c['command'] for c in commands)
assert not any('-ffast-math' in c['command'] or '-Ofast' in c['command'] for c in commands)
lock=a.build/'linux-dependency-authentication.json'
tc=Path(os.environ['TC']);dep=tc/'dependency-lock-linux.json'
d=json.loads(dep.read_text())
result=dict(status='CANDIDATE_RUNTIME_RECEIPT',host=platform.node(),kernel=platform.release(),cpu=cpus,
    types=next(r.split('\t')[1:] for r in raw.splitlines() if r.startswith('types\t')),
    loaded_files={f:sha(f) for f in loaded},probe_sha256=sha(probe),raw_probe=raw,
    dependency_lock_sha256=sha(dep),dependency_closure=d['closure'],link_inputs=d['link_inputs'],
    compiler=d['toolchain'] if 'toolchain' in d else d.get('compiler'),
    compile_commands_sha256=sha(a.build/'compile_commands.json'),
    effective_commands=[c for c in commands if any(n in c['file'] for n in ('LocalThermodynamics.cpp','TrackRFreeGasThermodynamics.cpp','ChemicalResponse.cpp','GlobalUrcaChannelCoefficient.cpp','arithmetic_probe.cpp'))],
    runtime_dispatch_note='CPU feature class and loaded bytes observed; arithmetic results compared across allocated nodes, no inference solely from CPU labels')
a.output.write_text(json.dumps(result,indent=2)+'\n')
