#!/usr/bin/env python3
"""Preserve completed trajectories and migrate interrupted diagnostic collection.

Run only after the old writers have stopped. Partial diagnostic grids remain in
r2, never promoted to complete or overwritten. New grid collections use the same
completed main histories; no original-t0 trajectory is repeated.
"""
import json,pathlib,shutil,sys
from piecewise import ROOT
from run import digest,save
old=pathlib.Path(sys.argv[1]);new=pathlib.Path(sys.argv[2]);new.mkdir(exist_ok=False)
for name in ('preflight.json','hp-source-1000000.json','hp-control-1000000.json'):
 shutil.copy2(old/name,new/name)
shutil.copytree(old/'build',new/'build')
files={}
for role in ('source','control'):
 target=new/'native'/role;target.mkdir(parents=True)
 for p in (old/'native'/role).glob('reference-*'):
  shutil.copy2(p,target/p.name);files[str(target/p.name)]=digest(target/p.name)
 for tier in range(3):
  if role=='source' and tier<2:continue
  src=old/'native'/role/f'tier-{tier}-cap-0';dst=target/src.name;dst.mkdir()
  for folder in (pathlib.Path('.'),pathlib.Path('endpoint-schedule-repeat')):
   (dst/folder).mkdir(exist_ok=True)
   for name in ('main.accepted_states.tsv','events.tsv','main.internal_steps.tsv','main.audit.tsv','main.steps','main.observations.tsv','schedule.tsv'):
    p=src/folder/name
    if p.exists():
     q=dst/folder/name;shutil.copy2(p,q);assert digest(p)==digest(q);files[str(q)]=digest(q)
  assert (dst/'main.accepted_states.tsv').exists()
save(new/'retained-inputs.json',dict(files=files,source=str(old),reason='NFS-stalled diagnostic workers migrated after cancellation; immutable completed main histories retained; partial diagnostics remain untouched in source directory.'))
