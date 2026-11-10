#!/usr/bin/env python3
"""Authenticate on the producer, then dispatch unchanged native binaries via Bash.

Avoids Python startup stalls on compute nodes. Native archive/currentness guards
and all physics/numerical settings remain unchanged. No already-started task is
eligible for this dispatcher.
"""
import json,pathlib,shlex,sys
ROOT=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/linux/preliminary_pair'))
from run import SNAP,load_pair,digest,save
out=pathlib.Path(sys.argv[1]);pair=load_pair(SNAP);build=out/'build'
manifest=json.loads((build/'build.json').read_text())
for p,h in manifest['files'].items():assert digest(p)==h,p
base=out/'direct';base.mkdir(exist_ok=False)
critical={str(build/'driver'):digest(build/'driver')}
for p,h in manifest['files'].items():
 if p.endswith(('.so','.a')):critical[p]=h
for s in pair:
 for key in ('archive','coefficients','entry','frozen'):
  p=s[key];critical[p]=digest(p)
 p=str(pathlib.Path(s['coefficients']).parent/'tangent.tsv');critical[p]=digest(p)
(base/'inputs.sha256').write_text(''.join(f'{h}  {p}\n' for p,h in critical.items()))
for index in range(2,6):
 role='source' if index<3 else 'control';i=index%3;s=pair[role=='control']
 assert not (out/f'{role}-tier-{i}').exists(),'a task already entered native dispatch'
 work=out/f'{role}-tier-{i}-direct';work.mkdir(exist_ok=False)
 args=[str(build/'driver'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(work/'assembly'),str(out/'native'/role),role,'-1',str(i)]
 script=base/f'task-{index}.sh'
 script.write_text('#!/bin/bash\nset -u\nexport RF_STAGE_WORKERS=24\n'+
  'sha256sum --check '+shlex.quote(str(base/'inputs.sha256'))+' > '+shlex.quote(str(work/'runtime-hashes.log'))+' || exit 90\n'+
  'start=$(date +%s)\n'+shlex.join(args)+' > '+shlex.quote(str(work/'run.log'))+' 2>&1\n'+
  'code=$?\nend=$(date +%s)\n'+
  'printf \'{"exit_code":%s,"wall_seconds":%s,"job":"%s"}\\n\' "$code" "$((end-start))" "$SLURM_JOB_ID" > '+shlex.quote(str(work/'completion.json'))+'\nexit "$code"\n')
 save(work/'invocation.json',dict(args=args,snapshot=digest(SNAP/(role+'.json')),build=digest(build/'build.json'),script=str(script),script_sha256=digest(script),producer_authentication='load_pair and every compiled/input binding checked before dispatch; essential runtime bytes checked again on compute node; native archive and currentness guards unchanged.'))
save(base/'producer-authentication.json',dict(status='PASS',snapshot_ids=json.loads((SNAP/'identities.json').read_text()),build_sha256=digest(build/'build.json'),critical_runtime_files=critical,executed_files={str(p):digest(p) for p in base.glob('task-*.sh')}))
