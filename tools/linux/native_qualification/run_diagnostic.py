#!/usr/bin/env python3
"""One source-only diagnostic replay, authenticated separately from original snapshots."""
import json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'tools/linux/native_pair'))
from snapshot import load_pair,digest
from runtime import runtime_binding
snapshots=root/'docs/validation/results-first-bnv-pilot/native-runtime-retry/snapshots'
build=pathlib.Path(sys.argv[1]);retained=pathlib.Path(sys.argv[2]);work=pathlib.Path(sys.argv[3])
m=json.loads((build/'manifest.json').read_text())
for item in m['files']:
 if digest(item['path'])!=item['sha256']:raise ValueError('diagnostic build dependency mismatch: '+item['path'])
for role in ['source','control']:
 if digest(snapshots/(role+'.json'))!=m['snapshots'][role]:raise ValueError('original snapshot identity mismatch')
s=load_pair(snapshots)[0]
if runtime_binding(build/'checkpoint_replay')!=m['runtime_binding']:raise ValueError('diagnostic runtime mismatch')
# Authenticate the retained scientific files before they can be used as references.
checks=root/'docs/validation/results-first-bnv-pilot/native-runtime-retry/evidence/external-SHA256SUMS'
for line in checks.read_text().splitlines():
 sha,name=line.split(maxsplit=1)
 if str(pathlib.Path(name)).startswith(str(retained)+'/') and digest(name)!=sha:raise ValueError('retained evidence mismatch: '+name)
work.mkdir(parents=True,exist_ok=False)
args=[str(build/'checkpoint_replay'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(work/'assembly'),str(retained/'run/trajectory/source'),str(work/'output')]
(work/'invocation.json').write_text(json.dumps(dict(args=args,snapshots=m['snapshots'],build_manifest_sha256=digest(build/'manifest.json')),indent=2)+'\n')
sys.exit(subprocess.run(args,cwd=root).returncode)
